#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Author: Samuel Nahmani
# Description: Compute VMF1 mapping values (hydro & wet) for a station/epoch/elevations.
# Defaults reproduce a known example. Supports batch CSV of (epoch,elev).

import os
import csv
import argparse
from math import floor
from datetime import datetime, timezone, timedelta
from pathlib import Path
from typing import Iterable, List, Tuple, Dict

import numpy as np

# --- Local project modules (present in your repo) ---
import download_VMF as dl
import read_vmf1_grid as rvg
import vmf1_ht

# Robustify path handling inside read_vmf1_grid if it builds Windows-like paths
import builtins
rvg.open = lambda path, *args, **kwargs: builtins.open(path.replace("\\", os.sep), *args, **kwargs)

# ---------------- time utilities ----------------
def parse_iso_utc(iso_utc: str) -> datetime:
    """Parse ISO 8601 string and return timezone-aware datetime in UTC."""
    dt = datetime.fromisoformat(iso_utc)
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
def ensure_vmf1_files_for_epoch(dt_utc: datetime, vmf_dir: str):
    """
    Ensure needed VMF1 grids exist in vmf_dir.
    If missing, call download_VMF.download_VMF(year, month, day) for each required date.
    """
    os.makedirs(vmf_dir, exist_ok=True)
    needed_dates = six_hour_bracketing_dates(dt_utc)  # at most two dates
    for (y, m, d) in sorted(needed_dates):
        wanted = [f"VMFG_{y}{str(m).zfill(2)}{str(d).zfill(2)}.H{hh}" for hh in ("00", "06", "12", "18")]
        have = all(os.path.exists(os.path.join(vmf_dir, w)) for w in wanted)
        if not have:
            print(f"[info] Grilles VMF1 manquantes pour {y}-{str(m).zfill(2)}-{str(d).zfill(2)}. Téléchargement…")
            cwd = os.getcwd()
            try:
                os.chdir(vmf_dir)
                dl.download_VMF(year=str(y), month=str(m).zfill(2), day=d)
            finally:
                os.chdir(cwd)

# ---------------- pairs file ----------------
def read_pairs_file(path: Path) -> List[Tuple[datetime, float]]:
    """
    Read CSV/TSV with two columns: epoch_iso,elev_deg
    Auto-detect comma or tab. Header optional.
    """
    text = path.read_text(encoding="utf-8").strip().splitlines()
    if not text:
        return []
    delim = "," if ("," in text[0]) else "\t"
    out: List[Tuple[datetime, float]] = []
    for row in csv.reader(text, delimiter=delim):
        if not row or len(row) < 2:
            continue
        if row[0].strip().lower().startswith("epoch"):
            continue
        dt = parse_iso_utc(row[0].strip())
        ev = float(row[1].strip())
        out.append((dt, ev))
    return out

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
    # ensure grids present
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
              "Batch mode; use --out to save results.")
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
        help=("Directory where VMF1 grids are stored/downloaded "
              "(default: Mapping_Fcn/vmf1). "
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
        pairs = read_pairs_file(args.pairs)
        if not pairs:
            print(f"[warn] No valid pairs found in {args.pairs}")
            return

        # group by epoch to reuse ah/aw
        from collections import defaultdict
        elevs_by_epoch: Dict[str, List[float]] = defaultdict(list)
        for dt, ev in pairs:
            elevs_by_epoch[dt.isoformat()].append(float(ev))

        all_csv_rows: List[Tuple[str, float, float, float]] = []

        for epoch_iso in sorted(elevs_by_epoch.keys()):
            ah, aw, rows = compute_for_epoch(station, epoch_iso, elevs_by_epoch[epoch_iso], vmf_dir)
            print_header(station, epoch_iso, ah, aw)
            print_table(rows)
            print()
            for el_deg, vmf1h, vmf1w in rows:
                all_csv_rows.append((epoch_iso, el_deg, vmf1h, vmf1w))

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

