#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# compute_vmf1_values-TER.py
#
# Default behavior (no options): SIMPLE mode with historical defaults,
# prints table to stdout, auto-downloads VMF grids if missing.
#
# Adds:
#  - zhd/zwd read from VMF1 grids (when available)
#  - SlantHydro = vmf1h * zhd, SlantWet = vmf1w * zwd, SlantTotal
#  - BATCH mode with --pairs <2-col file epoch_iso,elev_deg> and optional --out CSV
#  - Optional ZHD altitude correction (OFF by default to preserve original behavior)
#

import os
import re
import csv
import argparse
from math import floor
from datetime import datetime, timezone, timedelta
from pathlib import Path
from typing import Iterable, List, Tuple, Dict, Optional

import numpy as np

# --- Project-local modules (expected in PYTHONPATH or same dir) ---
import download_VMF as dl
import read_vmf1_grid as rvg   # must provide read_vmf1_grid_ex (ah,aw,zhd,zwd)
import vmf1_ht

# -------------- Time / MJD utilities --------------
def parse_iso_utc(iso_utc: str) -> datetime:
    s = iso_utc.strip()
    if s.endswith("Z"):
        s = s[:-1]
    try:
        dt = datetime.fromisoformat(s)
    except ValueError:
        dt = datetime.fromisoformat(s.replace(" ", "T"))
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    else:
        dt = dt.astimezone(timezone.utc)
    return dt

def iso_to_mjd(iso_utc: str) -> float:
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
    base = dt_utc.replace(minute=0, second=0, microsecond=0)
    floor6 = base.replace(hour=(base.hour//6)*6)
    ceil6  = floor6 + timedelta(hours=6)
    return {
        (floor6.year, floor6.month, floor6.day),
        (ceil6.year,  ceil6.month,  ceil6.day)
    }

# -------------- VMF1 grids handling --------------
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

# -------------- Robust pairs reader (batch) --------------
def _try_parse_epoch(s: str) -> Optional[datetime]:
    s = s.strip()
    if not s:
        return None
    s2 = s[:-1] if s.endswith("Z") else s
    try:
        dt = datetime.fromisoformat(s2)
    except ValueError:
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
    lines = [ln.strip("\ufeff").strip() for ln in text if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        return []

    if forced_delim:
        candidates = [forced_delim]
    else:
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
                    if i == 0:
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

    # Fallback: spaces
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

# -------------- Presentation --------------
def print_header(station: Dict, epoch_iso: str, ah: float, aw: float, zhd: float, zwd: float):
    mjd = iso_to_mjd(epoch_iso)
    print(f"# Station   : {station['name']}")
    print(f"# Lat, Lon  : {station['lat_deg']:.6f}°, {station['lon_deg']:.6f}°")
    print(f"# Height    : {station['h_m']:.3f} m")
    print(f"# Epoch UTC : {epoch_iso}  (MJD {mjd:.6f})\n")
    print(f"ah (grid) = {ah:.6e}")
    print(f"aw (grid) = {aw:.6e}")
    if np.isfinite(zhd):
        print(f"zhd (grid) = {zhd:.6e}  [m]")
    else:
        print("zhd (grid) = NaN        [m]")
    if np.isfinite(zwd):
        print(f"zwd (grid) = {zwd:.6e}  [m]\n")
    else:
        print("zwd (grid) = NaN        [m]\n")
    print("Elevation(deg),  VMF1_hydro,  VMF1_wet,   SlantHydro[m],  SlantWet[m],  SlantTotal[m]")

def print_table(rows: List[Tuple[float, float, float, float, float, float]]):
    for el_deg, vmf1h, vmf1w, sh, sw, stot in rows:
        print(f"{el_deg:>6.1f}          {vmf1h:>10.6f}   {vmf1w:>10.6f}   {sh:>12.6f}   {sw:>10.6f}   {stot:>12.6f}")

def write_csv(out_path: Path, station: Dict, records: List[Tuple[str, float, float, float, float, float, float, float, float]]):
    """
    CSV columns:
    station,lat_deg,lon_deg,height_m,epoch_iso,mjd,elev_deg,vmf1_hydro,vmf1_wet,ah,aw,zhd_m,zwd_m,slant_h_m,slant_w_m,slant_total_m
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["station","lat_deg","lon_deg","height_m","epoch_iso","mjd",
                    "elev_deg","vmf1_hydro","vmf1_wet","ah","aw","zhd_m","zwd_m","slant_h_m","slant_w_m","slant_total_m"])
        for (epoch_iso, elev, vmf1h, vmf1w, ah, aw, zhd, zwd, sh, sw, stot) in records:
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
                f"{ah:.6e}",
                f"{aw:.6e}",
                f"{zhd:.6e}",
                f"{zwd:.6e}",
                f"{sh:.6e}",
                f"{sw:.6e}",
                f"{stot:.6e}",
            ])

# -------------- ZHD altitude correction helpers --------------
import math
def saastamoinen_zhd_from_pressure(P_hPa, lat_deg, h_m):
    h_km = h_m / 1000.0
    denom = 1.0 - 0.00266*math.cos(math.radians(2.0*lat_deg)) - 0.00028*h_km
    return 0.0022768 * (P_hPa) / denom  # m

def pressure_from_zhd(zhd_m, lat_deg, h_m):
    h_km = h_m / 1000.0
    denom = 1.0 - 0.00266*math.cos(math.radians(2.0*lat_deg)) - 0.00028*h_km
    return zhd_m * denom / 0.0022768  # hPa

def barometric_adjust(P_ref_hPa, dh_m, H_km=7.8):
    # P(h+dh) = P(h) * exp(-dh/H); H in km
    return P_ref_hPa * math.exp(-(dh_m/1000.0)/H_km)

# -------------- Core computation --------------
def compute_for_epoch(station: Dict, epoch_iso: str, elevations_deg: Iterable[float],
                      vmf_dir: Path) -> Tuple[float, float, float, float, List[Tuple[float, float, float, float, float, float]]]:
    """
    Return (ah, aw, zhd, zwd, rows) where
      rows = [(elev_deg, vmf1_h, vmf1_w, slant_h, slant_w, slant_total), ...]
    """
    # ensure grids present for this epoch
    dt = parse_iso_utc(epoch_iso)
    ensure_vmf1_files_for_epoch(dt, str(vmf_dir))

    # read (ah,aw,zhd,zwd) from grid (time+space interpolation inside)
    mjd = iso_to_mjd(epoch_iso)
    ah_vec, aw_vec, zhd_vec, zwd_vec = rvg.read_vmf1_grid_ex(
        mjd=[mjd],
        ell=[station["lat_deg"], station["lon_deg"], station["h_m"]]
    )
    ah = float(np.atleast_1d(ah_vec).ravel()[0])
    aw = float(np.atleast_1d(aw_vec).ravel()[0])
    zhd = float(np.atleast_1d(zhd_vec).ravel()[0]) if np.atleast_1d(zhd_vec).size else float("nan")
    zwd = float(np.atleast_1d(zwd_vec).ravel()[0]) if np.atleast_1d(zwd_vec).size else float("nan")

    # compute mapping via vmf1_ht (zd = 90° - elev)
    dlat_rad = np.deg2rad(station["lat_deg"])
    ht_m = station["h_m"]

    rows: List[Tuple[float, float, float, float, float, float]] = []
    for el_deg in elevations_deg:
        zd_rad = np.deg2rad(90.0 - float(el_deg))
        vmf1h, vmf1w = vmf1_ht.vmf1_ht(ah=ah, aw=aw, dmjd=mjd, dlat=dlat_rad, ht=ht_m, zd=zd_rad)
        # Slant delays along LOS if zenith delays available
        slant_h = vmf1h * zhd if np.isfinite(zhd) else float("nan")
        slant_w = vmf1w * zwd if np.isfinite(zwd) else float("nan")
        slant_total = (slant_h + slant_w) if (np.isfinite(slant_h) and np.isfinite(slant_w)) else float("nan")
        rows.append((float(el_deg), float(vmf1h), float(vmf1w), float(slant_h), float(slant_w), float(slant_total)))

    return ah, aw, zhd, zwd, rows

# -------------- argparse (keeps original defaults) --------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Compute VMF1 mapping values (hydro & wet) for a station and epoch, "
            "and also read zhd/zwd from the grids (when available).\n\n"
            "Two modes:\n"
            "  - Simple mode (default): --epoch + --elev ... → prints results to stdout.\n"
            "  - Batch mode: --pairs <CSV/TSV epoch_iso,elev_deg> → computes many rows; "
            "use --out to save a CSV file."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=True,
    )

    # Defaults identical to original example
    p.add_argument("--station", default="DUNQ00FRA", help="Station label (for traceability)")
    p.add_argument("--lat", type=float, default=51.048100, help="Latitude in degrees (WGS84)")
    p.add_argument("--lon", type=float, default=2.366699, help="Longitude in degrees (WGS84)")
    p.add_argument("--height", type=float, default=51.111900, help="Ellipsoidal height (meters)")

    p.add_argument("--epoch", default="2021-07-10T01:00:18",
                   help="Epoch in ISO 8601 (UTC), e.g. 2021-07-10T01:00:18")
    p.add_argument("--elev", type=float, nargs="*", default=[5,10,15,20,30,45,60],
                   help="Elevation angles in degrees (space-separated). Default: 5 10 15 20 30 45 60")

    # Batch mode
    p.add_argument("--pairs", type=Path,
                   help=("CSV/TSV with 2 columns: epoch_iso,elev_deg (header optional). "
                         "Separator auto-detected among comma, semicolon, tab, pipe, or spaces. "
                         "Lines starting with '#' are ignored."))
    p.add_argument("--pairs_delim", choices=[",",";","tab","|","space"],
                   help=("Force delimiter for --pairs (default: auto-detect). "
                         "Use 'tab' for \\t and 'space' for one/multiple spaces."))
    p.add_argument("--out", type=Path, help="Optional CSV output file (batch mode only).")

    # VMF grids directory
    p.add_argument("--vmf_dir", type=Path, default=Path("Mapping_Fcn")/"vmf1",
                   help="Directory where VMF1 grids are stored/downloaded (default: Mapping_Fcn/vmf1).")

    # ZHD altitude correction (OFF by default to preserve legacy outputs)
    p.add_argument("--zhd_alt_correction", choices=["on","off"], default="off",
                   help="Apply altitude correction to ZHD (barometric law + Saastamoinen). Default: off")
    p.add_argument(
    "--Hkm", type=float, default=7.8,
    help=("Barometric scale height in km (used only if --zhd_alt_correction on). "
          "Typical range 7.0–8.5 km; default = 7.8 km (US Standard Atmosphere).")
)


    return p

# -------------- main --------------
def main():
    args = build_parser().parse_args()

    station = {
        "name": args.station,
        "lon_deg": float(args.lon),
        "lat_deg": float(args.lat),
        "h_m": float(args.height),
    }
    vmf_dir: Path = args.vmf_dir

    # ----- BATCH MODE -----
    if args.pairs is not None:
        forced = None
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

        # Group elevations by epoch
        from collections import defaultdict
        elevs_by_epoch_dt: Dict[datetime, List[float]] = defaultdict(list)
        for dt, ev in pairs:
            elevs_by_epoch_dt[dt].append(float(ev))

        # Ensure grids
        needed_calendar_days = set()
        for dt in elevs_by_epoch_dt.keys():
            needed_calendar_days |= six_hour_bracketing_dates(dt)
        ensure_vmf1_files_for_dates(needed_calendar_days, str(vmf_dir))

        # Vectorized (ah,aw,zhd,zwd) for all epochs
        unique_epochs_sorted = sorted(elevs_by_epoch_dt.keys())
        mjd_list = [iso_to_mjd(dt.isoformat()) for dt in unique_epochs_sorted]
        ah_vec, aw_vec, zhd_vec, zwd_vec = rvg.read_vmf1_grid_ex(
            mjd=mjd_list,
            ell=[station["lat_deg"], station["lon_deg"], station["h_m"]]
        )

        all_csv_rows: List[Tuple[str, float, float, float, float, float, float, float, float]] = []
        dlat_rad = np.deg2rad(station["lat_deg"])
        ht_m = station["h_m"]

        def zhd_alt_adjust(zhd_raw):
            if args.zhd_alt_correction != "on" or not np.isfinite(zhd_raw):
                return float(zhd_raw)
            P_g = pressure_from_zhd(zhd_raw, station["lat_deg"], 0.0)
            P_s = barometric_adjust(P_g, ht_m - 0.0, H_km=args.Hkm)
            return saastamoinen_zhd_from_pressure(P_s, station["lat_deg"], ht_m)

        for idx, dt in enumerate(unique_epochs_sorted):
            epoch_iso = dt.isoformat()
            ah = float(np.atleast_1d(ah_vec)[idx])
            aw = float(np.atleast_1d(aw_vec)[idx])
            zhd_raw = float(np.atleast_1d(zhd_vec)[idx]) if np.atleast_1d(zhd_vec).size else float("nan")
            zwd = float(np.atleast_1d(zwd_vec)[idx]) if np.atleast_1d(zwd_vec).size else float("nan")
            zhd = zhd_alt_adjust(zhd_raw)

            rows: List[Tuple[float, float, float, float, float, float]] = []
            for el_deg in elevs_by_epoch_dt[dt]:
                zd_rad = np.deg2rad(90.0 - float(el_deg))
                vmf1h, vmf1w = vmf1_ht.vmf1_ht(ah=ah, aw=aw, dmjd=mjd_list[idx],
                                               dlat=dlat_rad, ht=ht_m, zd=zd_rad)
                slant_h = vmf1h * zhd if np.isfinite(zhd) else float("nan")
                slant_w = vmf1w * zwd if np.isfinite(zwd) else float("nan")
                slant_total = (slant_h + slant_w) if (np.isfinite(slant_h) and np.isfinite(slant_w)) else float("nan")
                rows.append((float(el_deg), float(vmf1h), float(vmf1w), float(slant_h), float(slant_w), float(slant_total)))

            print_header(station, epoch_iso, ah, aw, zhd, zwd)
            print_table(rows)
            print()

            for el_deg, vmf1h, vmf1w, sh, sw, stot in rows:
                all_csv_rows.append((epoch_iso, el_deg, vmf1h, vmf1w, ah, aw, zhd, zwd, sh, sw, stot))

        if args.out:
            write_csv(args.out, station, all_csv_rows)
            print(f"[info] CSV written to {args.out}")
        return

    # ----- SIMPLE MODE (default) -----
    epoch_iso = args.epoch
    ensure_vmf1_files_for_epoch(parse_iso_utc(epoch_iso), str(vmf_dir))

    ah, aw, zhd, zwd, rows = compute_for_epoch(station, epoch_iso, args.elev, vmf_dir)

    # Optional ZHD altitude correction (OFF by default)
    if args.zhd_alt_correction == "on" and np.isfinite(zhd):
        P_g = pressure_from_zhd(zhd, station["lat_deg"], 0.0)
        P_s = barometric_adjust(P_g, station["h_m"] - 0.0, H_km=args.Hkm)
        zhd = saastamoinen_zhd_from_pressure(P_s, station["lat_deg"], station["h_m"])
        # Recompute slant hydro with corrected ZHD
        dlat_rad = np.deg2rad(station["lat_deg"])
        mjd = iso_to_mjd(epoch_iso)
        new_rows = []
        for el_deg, vmf1h, vmf1w, _sh, sw, _st in rows:
            zd_rad = np.deg2rad(90.0 - float(el_deg))
            vmf1h2, vmf1w2 = vmf1_ht.vmf1_ht(ah=ah, aw=aw, dmjd=mjd, dlat=dlat_rad, ht=station["h_m"], zd=zd_rad)
            sh2 = vmf1h2 * zhd
            st2 = (sh2 + sw) if np.isfinite(sw) else float("nan")
            new_rows.append((el_deg, vmf1h2, vmf1w2, sh2, sw, st2))
        rows = new_rows

    print_header(station, epoch_iso, ah, aw, zhd, zwd)
    print_table(rows)

# -------------- Entrypoint --------------
if __name__ == "__main__":
    main()

