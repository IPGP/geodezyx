#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Author: Samuel Nahmani
# Description: Compute VMF1 mapping values (hydro & wet) for a station/epoch/elevations.
# - Simple mode (default): 1 epoch + list of elevations → prints to stdout
# - Batch mode: --pairs <2-col file epoch_iso,elev_deg> → vectorized grids read; use --out to save CSV
#
# Notes:
#  - VMF1 grids directory default: Mapping_Fcn/vmf1 ; missing grids are downloaded automatically.
#  - --pairs file: separator auto-detected (comma, semicolon, tab, pipe, spaces). Header optional. Lines starting with # ignored.

import os
import re
import csv
import argparse
from math import floor
from datetime import datetime, timezone, timedelta
from pathlib import Path
from typing import Iterable, List, Tuple, Dict, Optional

import numpy as np

# --- Local project modules ---
import download_VMF as dl
import read_vmf1_grid as rvg
import vmf1_ht

# Robustify path handling inside read_vmf1_grid if it builds Windows-like paths
import builtins
rvg.open = lambda path, *args, **kwargs: builtins.open(path.replace("\\", os.sep), *args, **kwargs)

# ---------------- time utilities ----------------
def parse_iso_utc(iso_utc: str) -> datetime:
    """Parse ISO 8601 string and return timezone-aware datetime in UTC."""
    s = iso_utc.strip()
    if s.endswith("Z"):
        s = s[:-1]
    try:
        dt = datetime.fromisoformat(s)
    except ValueError:
        # tolerate space instead of 'T'
        dt = datetime.fromisoformat(s.replace(" ", "T"))
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    else:
        dt = dt.astimezone(timezone.utc)
    return dt

def iso_to_mjd(iso_utc: str) -> float:
    """Convert ISO UTC string to MJD (UTC)."""
    dt = parse_iso_utc(iso_utc)
    Y, M, D = dt.year, dt.month, dt.day
    frac_day = (dt.hour + dt.minute/60 + dt.second/3600 + dt.microsecond/3.6e9)/24.0
    if M <= 2:
        Yp, Mp = Y - 1, M + 12
    else:
        Yp, Mp = Y, M
    A = floor(Yp/100)
    B = 2 - A + floor(A/4)
    JD = floor(365.25*(Yp+4716)) + floor(30.6001*(Mp+1)) + (D + frac_day) + B - 1524.5
    return JD - 2400000.5

def six_hour_bracketing_dates(dt_utc: datetime):
    """Return the calendar dates potentially involved in 6-hour interpolation around dt."""
    base = dt_utc.replace(minute=0, second=0, microsecond=0)
    floor6 = base.replace(hour=(base.hour//6)*6)
    ceil6  = floor6 + timedelta(hours=6)
    return {
        (floor6.year, floor6.month, floor6.day),
        (ceil6.year,  ceil6.month,  ceil6.day)
    }

# ---------------- VMF1 grids handling ----------------
def ensure_vmf1_files_for_dates(dates_set: set, vmf_dir: str):
    """
    Ensure needed VMF1 grids exist in vmf_dir for all (Y,M,D) in dates_set.
    If missing, call download_VMF.download_VMF once per date.
    """
    os.makedirs(vmf_dir, exist_ok=True)
    for (y, m, d) in sorted(dates_set):
        wanted = [f"VMFG_{y}{m:02d}{d:02d}.H{hh}" for hh in ("00", "06", "12", "18")]
        have = all(os.path.exists(os.path.join(vmf_dir, w)) for w in wanted)
        if not have:
            print(f"[info] Grilles VMF1 manquantes pour {y}-{m:02d}-{d:02d}. Téléchargement…")
            cwd = os.getcwd()
            try:
                os.chdir(vmf_dir)
                dl.download_VMF(year=str(y), month=f"{m:02d}", day=d)
            finally:
                os.chdir(cwd)

def ensure_vmf1_files_for_epoch(dt_utc: datetime, vmf_dir: str):
    """Backward-compatible helper: ensure for one epoch."""
    ensure_vmf1_files_for_dates(six_hour_bracketing_dates(dt_utc), vmf_dir)

# ---------------- robust pairs reader ----------------
def _try_parse_epoch(s: str) -> Optional[datetime]:
    """Return UTC datetime if parseable, else None."""
    s = s.strip()
    if not s:
        return None
    if s.endswith("Z"):
        s2 = s[:-1]
    else:
        s2 = s
    try:
        dt = datetime.fromisoformat(s2)
    except ValueError:
        # tolerate space instead of 'T'
        try:
            dt = datetime.fromisoformat(s2.replace(" ", "T"))
        except ValueError:
            return None
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    else:
        dt = dt.astimezone(timezone.utc)
    return dt

def read_pairs_file(path: Path, forced_delim: Optional[str] = None) -> List[Tuple[datetime, float]]:
    """
    Read 2-column file: epoch_iso, elev_deg.
    - Auto-detects delimiter (csv.Sniffer -> ',', ';', '\\t', '|', ' ').
    - Falls back to regex split on spaces.
    - Header optional (first line skipped if not a date).
    - Lines starting with '#' are ignored.
    """
    text = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    # drop comments and empty
    lines = [ln.strip("\ufeff").strip() for ln in text if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        return []

    # 1) forced delimiter wins
    if forced_delim:
        candidates = [forced_delim]
    else:
        # 2) try csv.Sniffer
        try:
            dialect = csv.Sniffer().sniff("\n".join(lines[:10]), delimiters=[",",";","\t","|"," "])
            candidates = [dialect.delimiter]
        except Exception:
            candidates = [",",";","\t","|"," "]

    def parse_with_delim(delim: str) -> List[Tuple[datetime, float]]:
        out: List[Tuple[datetime, float]] = []
        if delim == " ":
            for i, ln in enumerate(lines):
                parts = re.split(r"\s+", ln)
                if len(parts) < 2:
                    continue
                dt = _try_parse_epoch(parts[0])
                if dt is None:
                    if i == 0:  # header
                        continue
                    else:
                        continue
                try:
                    elev = float(parts[1])
                except ValueError:
                    continue
                out.append((dt, elev))
            return out

        rdr = csv.reader(lines, delimiter=delim)
        first = True
        for row in rdr:
            if not row or len(row) < 2:
                continue
            if first:
                dt = _try_parse_epoch(row[0])
                if dt is None:
                    first = False
                    continue  # header
                first = False
            else:
                dt = _try_parse_epoch(row[0])
                if dt is None:
                    continue
            try:
                elev = float(row[1])
            except ValueError:
                continue
            out.append((dt, elev))
        return out

    for d in candidates:
        parsed = parse_with_delim(d)
        if parsed:
            return parsed

    # last resort: spaces regex
    fallback: List[Tuple[datetime, float]] = []
    for i, ln in enumerate(lines):
        parts = re.split(r"\s+", ln)
        if len(parts) < 2:
            continue
        dt = _try_parse_epoch(parts[0])
        if dt is None:
            if i == 0:
                continue
            else:
                continue
        try:
            elev = float(parts[1])
        except ValueError:
            continue
        fallback.append((dt, elev))
    return fallback

# ---------------- presentation ----------------
def print_header(station: Dict, epoch_iso: str, ah: float, aw: float):
    mjd = iso_to_mjd(epoch_iso)
    print(f"# Station   : {station['name']}")
    print(f"# Lat, Lon  : {station['lat_deg']:.6f}°, {station['lon_deg']:.6f}°")
    print(f"# Height    : {station['h_m']:.3f} m")
    print(f"# Epoch UTC : {epoch_iso}  (MJD {mjd:.6f})\n")
    print(f"ah (grid) = {ah:.6e}")
    print(f"aw (grid) = {aw:.6e}\n")
    print("Elevation(deg),  VMF1_hydro,  VMF1_wet")

def print_table(rows: List[Tuple[float, float, float]]):
    for el_deg, vmf1h, vmf1w in rows:
        print(f"{el_deg:>6.1f}          {vmf1h:>10.6f}   {vmf1w:>10.6f}")

def write_csv(out_path: Path, station: Dict, records: List[Tuple[str, float, float, float]]):
    """
    Write CSV with columns:
    station,lat_deg,lon_deg,height_m,epoch_iso,mjd,elev_deg,vmf1_hydro,vmf1_wet
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["station","lat_deg","lon_deg","height_m","epoch_iso","mjd","elev_deg","vmf1_hydro","vmf1_wet"])
        for epoch_iso, elev, vmf1h, vmf1w in records:
            mjd = iso_to_mjd(epoch_iso)
            w.writerow([
                station["name"],
                f"{station['lat_deg']:.6f}",
                f"{station['lon_deg']:.6f}",
                f"{station['h_m']:.3f}",
                epoch_iso,
                f"{mjd:.6f}",
                f"{elev:.1f}",
                f"{vmf1h:.6f}",
                f"{vmf1w:.6f}",
            ])

# ---------------- core computation ----------------
def compute_for_epoch(station: Dict, epoch_iso: str, elevations_deg: Iterable[float],
                      vmf_dir: Path) -> Tuple[float, float, List[Tuple[float, float, float]]]:
    """
    Return (ah, aw, rows) where rows = [(elev_deg, vmf1_h, vmf1_w), ...]
    """
    # ensure grids present for this epoch
    dt = parse_iso_utc(epoch_iso)
    ensure_vmf1_files_for_epoch(dt, str(vmf_dir))

    # read (ah,aw) from grid
    mjd = iso_to_mjd(epoch_iso)
    ah_vec, aw_vec = rvg.read_vmf1_grid(mjd=[mjd],
                                        ell=[station["lat_deg"], station["lon_deg"], station["h_m"]])
    ah = float(np.atleast_1d(ah_vec).ravel()[0])
    aw = float(np.atleast_1d(aw_vec).ravel()[0])

    # compute mapping via vmf1_ht (zd = 90° - elev)
    dlat_rad = np.deg2rad(station["lat_deg"])
    ht_m = station["h_m"]

    rows: List[Tuple[float, float, float]] = []
    for el_deg in elevations_deg:
        zd_rad = np.deg2rad(90.0 - float(el_deg))
        vmf1h, vmf1w = vmf1_ht.vmf1_ht(ah=ah, aw=aw, dmjd=mjd, dlat=dlat_rad, ht=ht_m, zd=zd_rad)
        rows.append((float(el_deg), float(vmf1h), float(vmf1w)))

    return ah, aw, rows

# ---------------- argparse ----------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Compute VMF1 mapping values (hydro & wet) for a station and epoch.\n\n"
            "Two modes:\n"
            "  - Simple mode (default): --epoch + --elev ... → prints results to stdout.\n"
            "  - Batch mode: --pairs <CSV/TSV epoch_iso,elev_deg> → computes many rows; "
            "use --out to save a CSV file."
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    # station & coordinates
    p.add_argument("--station", default="DUNQ00FRA", help="Station label (for traceability)")
    p.add_argument("--lat", type=float, default=51.048100, help="Latitude in degrees (WGS84)")
    p.add_argument("--lon", type=float, default=2.366699, help="Longitude in degrees (WGS84)")
    p.add_argument("--height", type=float, default=51.111900, help="Ellipsoidal height (meters)")
    # time & elevations
    p.add_argument("--epoch", default="2021-07-10T01:00:18",
                   help="Epoch in ISO 8601 (UTC), e.g. 2021-07-10T01:00:18")
    p.add_argument("--elev", type=float, nargs="*", default=[5,10,15,20,30,45,60],
                   help="Elevation angles in degrees (space-separated). Default: 5 10 15 20 30 45 60")
    # batch mode
    p.add_argument(
        "--pairs", type=Path,
        help=("CSV/TSV with 2 columns: epoch_iso,elev_deg (header optional). "
              "Separator auto-detected among comma, semicolon, tab, pipe, or spaces. "
              "Lines starting with '#' are ignored. "
              "Batch mode; use --out to save results.")
    )
    # (optional) force delimiter if needed
    p.add_argument(
        "--pairs_delim",
        choices=[",",";","tab","|","space"],
        help=("Force delimiter for --pairs (default: auto-detect). "
              "Use 'tab' for \\t and 'space' for one/multiple spaces.")
    )
    # CSV output
    p.add_argument(
        "--out", type=Path,
        help=("Optional CSV output file. "
              "⚠️ Only active when --pairs is used (batch mode). "
              "In simple mode results are printed to stdout.")
    )
    # VMF grids directory
    p.add_argument(
        "--vmf_dir", type=Path, default=Path("Mapping_Fcn")/"vmf1",
        help=("Directory where VMF1 grids are stored/downloaded (default: Mapping_Fcn/vmf1). "
              "If missing, required grids are downloaded automatically.")
    )
    return p

# ---------------- main ----------------
def main():
    args = build_parser().parse_args()

    station = {
        "name": args.station,
        "lon_deg": float(args.lon),
        "lat_deg": float(args.lat),
        "h_m": float(args.height),
    }

    vmf_dir: Path = args.vmf_dir

    # Batch mode: file with (epoch_iso, elev_deg)
    if args.pairs:
        # 0) read pairs with robust autodetection
        forced = None  # type: Optional[str]
        if args.pairs_delim == "tab":
            forced = "\t"
        elif args.pairs_delim == "space":
            forced = " "
        elif args.pairs_delim in {",",";","|"}:
            forced = args.pairs_delim

        pairs = read_pairs_file(args.pairs, forced_delim=forced)
        if not pairs:
            print(f"[warn] No valid pairs found in {args.pairs}")
            return

        # group by epoch (UTC) → deduplicate epochs, reuse ah/aw per epoch
        from collections import defaultdict
        elevs_by_epoch_dt: Dict[datetime, List[float]] = defaultdict(list)
        for dt, ev in pairs:
            elevs_by_epoch_dt[dt].append(float(ev))

        # 1) build the set of calendar days needed for ALL epochs (for 6-hour interpolation)
        needed_calendar_days = set()
        for dt in elevs_by_epoch_dt.keys():
            needed_calendar_days |= six_hour_bracketing_dates(dt)

        # 2) ensure grids once
        ensure_vmf1_files_for_dates(needed_calendar_days, str(vmf_dir))

        # 3) vectorized read of ah/aw for all unique epochs
        unique_epochs_sorted = sorted(elevs_by_epoch_dt.keys())
        mjd_list = [iso_to_mjd(dt.isoformat()) for dt in unique_epochs_sorted]
        ah_vec, aw_vec = rvg.read_vmf1_grid(mjd=mjd_list,
                                            ell=[station["lat_deg"], station["lon_deg"], station["h_m"]])

        # 4) compute rows per epoch (reuse ah/aw)
        all_csv_rows: List[Tuple[str, float, float, float]] = []
        dlat_rad = np.deg2rad(station["lat_deg"])
        ht_m = station["h_m"]

        for idx, dt in enumerate(unique_epochs_sorted):
            epoch_iso = dt.isoformat()
            ah = float(np.atleast_1d(ah_vec)[idx])
            aw = float(np.atleast_1d(aw_vec)[idx])

            rows: List[Tuple[float, float, float]] = []
            for el_deg in elevs_by_epoch_dt[dt]:
                zd_rad = np.deg2rad(90.0 - float(el_deg))
                vmf1h, vmf1w = vmf1_ht.vmf1_ht(ah=ah, aw=aw,
                                               dmjd=mjd_list[idx], dlat=dlat_rad, ht=ht_m, zd=zd_rad)
                rows.append((float(el_deg), float(vmf1h), float(vmf1w)))

            print_header(station, epoch_iso, ah, aw)
            print_table(rows)
            print()

            for el_deg, vmf1h, vmf1w in rows:
                all_csv_rows.append((epoch_iso, el_deg, vmf1h, vmf1w))

        # 5) write CSV once
        if args.out:
            write_csv(args.out, station, all_csv_rows)
            print(f"[info] CSV written to {args.out}")
        return

    # Simple mode: single epoch + multiple elevations
    epoch_iso = args.epoch
    ah, aw, rows = compute_for_epoch(station, epoch_iso, args.elev, vmf_dir)
    print_header(station, epoch_iso, ah, aw)
    print_table(rows)

if __name__ == "__main__":
    main()

