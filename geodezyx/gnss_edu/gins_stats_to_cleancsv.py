#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GINS statistics (GNSS type=11) → CSV propre, entête façon .tropo.

- Parse *fixe* des colonnes GNSS (i3, e16.10, e9.3, f15.9, a8, f8.3, f8.3, i8, f8.3, f8.3, f12.3, f9.3, f9.3).
- Conversion de la date "e16.10" en ISO UTC (suffixe 'Z', sans µs), MJD, et SoD/JCN selon le cas.
- Détection automatique SoD / JCN(JUL50) / MJD / JD.
- "Hint JCN" récupéré depuis le nom du fichier: ..._<YYYY>_<DDD>_<JCN>__...
- Métadonnées station lues dans SPOTGINS_<STA>.tropo (Latitude/Longitude/Height).
- Délimiteur configurable (',' par défaut ; possible ';' ou ' ').

Usage :
  python3 gins_stats_to_cleancsv_fixed.py \
    --stats SPOTGINS_DUNQ00FRA_2021_191_26123__.yml.250819_135206.1 \
    --tropo SPOTGINS_DUNQ00FRA.tropo \
    --out DUNQ00FRA_stats.clean.csv \
    --delimiter ' '
"""

import re
import csv
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, Tuple, Dict
from datetime import datetime, timezone, timedelta
from math import floor

# ---------- Temps ----------
def iso_z(dt: datetime) -> str:
    """ISO8601 UTC sans microsecondes avec suffixe 'Z'."""
    return dt.astimezone(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')

def parse_year_doy_jcn_from_path(path: str) -> Tuple[int, int, Optional[float]]:
    """
    Extrait YYYY, DOY et éventuellement JCN (JUL50) du nom :
    ..._<YYYY>_<DDD>_<JCN>__...
    Ex. SPOTGINS_DUNQ00FRA_2021_191_26123__.yml...
    """
    m = re.search(r'[_\-.](20\d{2}|19\d{2})[_\-.](\d{3})(?:[_\-.](\d{5})(?:[_\-.]|$))?', path)
    if not m:
        raise ValueError("Impossible d'extraire YYYY/DDD/JCN du nom du fichier statistiques.")
    year = int(m.group(1)); doy = int(m.group(2))
    jcn_hint = float(m.group(3)) if m.group(3) is not None else None
    if not (1 <= doy <= 366):
        raise ValueError(f"DOY invalide: {doy}")
    return year, doy, jcn_hint

def doy_to_date(year: int, doy: int) -> datetime:
    return datetime(year, 1, 1, tzinfo=timezone.utc) + timedelta(days=doy-1)

def mjd_to_datetime(mjd: float) -> datetime:
    jd = mjd + 2400000.5
    Z = int(jd + 0.5)
    F = (jd + 0.5) - Z
    if Z < 2299161:
        A = Z
    else:
        alpha = int((Z - 1867216.25)/36524.25)
        A = Z + 1 + alpha - int(alpha/4)
    B = A + 1524
    C = int((B - 122.1)/365.25)
    D = int(365.25*C)
    E = int((B - D)/30.6001)
    day = B - D - int(30.6001*E) + F
    month = (E-1) if E < 14 else (E-13)
    year  = (C-4716) if month > 2 else (C-4715)
    day_whole = int(day); frac = day - day_whole
    sec = frac*86400.0
    hh = int(sec//3600); mm = int((sec%3600)//60); ss = int(sec%60)
    us = int(round((sec - (hh*3600+mm*60+ss))*1e6))
    return datetime(year, month, day_whole, hh, mm, ss, us, tzinfo=timezone.utc)

def datetime_to_mjd(dt: datetime) -> float:
    if dt.tzinfo is None: dt = dt.replace(tzinfo=timezone.utc)
    else: dt = dt.astimezone(timezone.utc)
    Y,M,D = dt.year, dt.month, dt.day
    frac = (dt.hour + dt.minute/60 + dt.second/3600 + dt.microsecond/3.6e9)/24.0
    if M <= 2: Yp, Mp = Y-1, M+12
    else: Yp, Mp = Y, M
    A = floor(Yp/100); B = 2 - A + floor(A/4)
    JD = floor(365.25*(Yp+4716)) + floor(30.6001*(Mp+1)) + (D + frac) + B - 1524.5
    return JD - 2400000.5

def detect_and_to_datetime(x: float,
                           year_hint: Optional[int],
                           doy_hint: Optional[int],
                           jcn_hint: Optional[float] = None,
                           mode: str = "auto"):
    """
    Renvoie (dt, sod, jcn). Modes: 'auto'|'sod'|'jcn'|'mjd'|'jd'.
    Règle 'auto':
      - si jcn_hint présent ET 18000<=x<=40000 → JCN (JUL50)
      - sinon si 0<=x<=86400 → SoD
      - sinon MJD/JD comme d'habitude
    """
    if mode == "jcn":
        base = datetime(1950, 1, 1, tzinfo=timezone.utc)
        return base + timedelta(days=x), None, float(x)
    if mode == "sod":
        if year_hint is None or doy_hint is None:
            raise ValueError("SoD forcé mais YYYY/DOY manquants.")
        base = doy_to_date(year_hint, doy_hint)
        return base + timedelta(seconds=x), float(x), None
    if mode == "mjd":
        return mjd_to_datetime(x), None, None
    if mode == "jd":
        return mjd_to_datetime(x - 2400000.5), None, None

    # AUTO
    if jcn_hint is not None and 18000.0 <= x <= 40000.0:
        base = datetime(1950, 1, 1, tzinfo=timezone.utc)
        return base + timedelta(days=x), None, float(x)

    if 0.0 <= x <= 86400.0:
        if year_hint is None or doy_hint is None:
            raise ValueError("SoD détecté mais YYYY/DOY manquants.")
        base = doy_to_date(year_hint, doy_hint)
        return base + timedelta(seconds=x), float(x), None

    if 50000.0 <= x <= 80000.0:
        return mjd_to_datetime(x), None, None

    if 2_400_000.0 <= x <= 2_600_000.0:
        return mjd_to_datetime(x - 2400000.5), None, None

    # dernier recours : si DOY dispo et x pas énorme → SoD
    if (year_hint is not None and doy_hint is not None) and (x < 200000.0):
        base = doy_to_date(year_hint, doy_hint)
        return base + timedelta(seconds=x), float(x), None

    raise ValueError(f"Date inconnue: {x}")

# ---------- Métadonnées station ----------
@dataclass
class StationMeta:
    name: str; lat: float; lon: float; h: float

def read_station_meta_from_tropo(tropo_path: Path) -> StationMeta:
    name = None; lat = None; lon = None; h = None
    with tropo_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if s.startswith("# STATION"):
                m = re.search(r':\s*(\S+)', s); name = m.group(1) if m else name
            elif s.startswith("# Latitude"):
                m = re.search(r':\s*([+-]?\d+(?:\.\d+)?)', s); lat = float(m.group(1)) if m else lat
            elif s.startswith("# Longitude"):
                m = re.search(r':\s*([+-]?\d+(?:\.\d+)?)', s); lon = float(m.group(1)) if m else lon
            elif s.startswith("# Height"):
                m = re.search(r':\s*([+-]?\d+(?:\.\d+)?)', s); h = float(m.group(1)) if m else h
            if name and (lat is not None) and (lon is not None) and (h is not None): break
    if name is None or lat is None or lon is None or h is None:
        raise ValueError(f"Métadonnées station manquantes dans {tropo_path}")
    return StationMeta(name=name, lat=lat, lon=lon, h=h)

# ---------- Parse fixe GNSS type=11 ----------
# i3,1x,e16.10,1x,e9.3,1x,f15.9,1x,a8,1x,f8.3,1x,f8.3,1x,i8,1x,f8.3,1x,f8.3,1x,f12.3,1x,f9.3,1x,f9.3
SPECS = [
    ("type", 3, "int"),
    ("sp",   1, "sp"),
    ("date", 16, "float"),
    ("sp",   1, "sp"),
    ("weight", 9, "float"),
    ("sp",   1, "sp"),
    ("residual", 15, "float"),
    ("sp",   1, "sp"),
    ("sat", 8, "str"),
    ("sp",   1, "sp"),
    ("sat_lat", 8, "float"),
    ("sp",   1, "sp"),
    ("sat_lon", 8, "float"),
    ("sp",   1, "sp"),
    ("station_code", 8, "int"),
    ("sp",   1, "sp"),
    ("az_deg", 8, "float"),
    ("sp",   1, "sp"),
    ("el_deg", 8, "float"),
    ("sp",   1, "sp"),
    ("pressure", 12, "float"),
    ("sp",   1, "sp"),
    ("temperature", 9, "float"),
    ("sp",   1, "sp"),
    ("humidity", 9, "float"),
]

def parse_fixed_gnss_line(line: str) -> Optional[Dict]:
    if not line.strip() or line.lstrip().startswith("#"):
        return None
    minlen = sum(w for _, w, _ in SPECS)
    if len(line) < minlen:
        return None
    pos = 0; rec: Dict = {}
    try:
        for key, width, kind in SPECS:
            field = line[pos:pos+width]; pos += width
            if kind == "sp":
                continue
            if kind == "int":
                rec[key] = int(field.strip())
            elif kind == "float":
                rec[key] = float(field.strip())
            elif kind == "str":
                rec[key] = field.strip()
        rec["tail"] = line[pos:].rstrip("\n")
    except Exception:
        return None
    if rec.get("type") != 11:
        return None
    return rec

# ---------- Conversion principale ----------
def convert_stats_to_cleancsv(stats_path: Path, tropo_path: Path, out_path: Path, delimiter: str = ","):
    year, doy, jcn_hint = parse_year_doy_jcn_from_path(stats_path.name)
    meta = read_station_meta_from_tropo(tropo_path)

    rows = []
    n_in = 0; n_out = 0
    with stats_path.open("r", encoding="utf-8", errors="ignore") as fin:
        for ln in fin:
            n_in += 1
            rec = parse_fixed_gnss_line(ln)
            if not rec:
                continue
            # Date → ISO/MJD/SoD/JCN
            dt, sod, jcn = detect_and_to_datetime(rec["date"], year, doy, jcn_hint=jcn_hint, mode="auto")
            mjd = datetime_to_mjd(dt)
            rows.append({
                "Date": iso_z(dt),
                "Type": rec["type"],
                "Sat": rec["sat"],
                "Az_deg": f"{rec['az_deg']:.3f}",
                "El_deg": f"{rec['el_deg']:.3f}",
                "Residual_m": f"{rec['residual']:.9f}",
                "Weight": f"{rec['weight']:.3f}",
                "Pressure": f"{rec['pressure']:.3f}",
                "Temperature": f"{rec['temperature']:.3f}",
                "Humidity": f"{rec['humidity']:.3f}",
                "MJD": f"{mjd:.6f}",
                "SoD": f"{sod:.3f}" if sod is not None else "",
                "JCN": f"{jcn:.6f}" if jcn is not None else "",
                "Tail": rec["tail"].strip()
            })
            n_out += 1

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        # en-tête façon .tropo (lisible)
        f.write("# SPOTGINS STATISTICS [GNSS TYPE=11] CLEAN EXPORT\n")
        f.write("# Creation (UTC)      : " + iso_z(datetime.utcnow()) + "\n")
        f.write("#--------------------------------------\n")
        f.write(f"# STATION             : {meta.name}\n")
        f.write("# ANALYSIS_CENTRE     : SPOTGINS - IPGP\n")
        f.write("# UNITS               : angles in degrees, residual in meters\n")
        f.write("#--------------------------------------\n")
        f.write(f"# Latitude            : {meta.lat:.6f}\n")
        f.write(f"# Longitude           : {meta.lon:.6f}\n")
        f.write(f"# Height              : {meta.h:.6f}\n")
        f.write("#--------------------------------------\n")
        f.write("# Columns:\n")
        f.write("#   Date, Type, Sat, Az_deg, El_deg, Residual_m, Weight, Pressure, Temperature, Humidity, MJD, SoD, JCN, Tail\n")
        f.write("#--------------------------------------\n")
        # ligne d’en-tête + données
        fieldnames = ["Date","Type","Sat","Az_deg","El_deg","Residual_m","Weight","Pressure","Temperature","Humidity","MJD","SoD","JCN","Tail"]
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter=delimiter, lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"[info] Station   : {meta.name} (lat={meta.lat:.6f}, lon={meta.lon:.6f}, h={meta.h:.3f})")
    print(f"[info] Day hint  : {year}-DOY{doy}  |  JCN hint: {jcn_hint if jcn_hint is not None else 'n/a'}")
    print(f"[info] Written   : {n_out} / {n_in} lignes → {out_path}")

# ---------- CLI ----------
if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="GINS stats GNSS (type=11) → CSV propre (parse *fixe*, dates ISO propres).")
    p.add_argument("--stats", required=True, type=Path, help="Fichier statistiques GINS (SPOTGINS_<STA>_YYYY_DDD_...JCN__...).yml...")
    p.add_argument("--tropo", required=True, type=Path, help="Fichier tropo SPOTGINS_<STA>.tropo (métadonnées station)")
    p.add_argument("--out", type=Path, help="Chemin de sortie (défaut: <stats>.clean.csv)")
    p.add_argument("--delimiter", default=",", help="Séparateur (',' par défaut). Exemple: ';' ou ' '")
    args = p.parse_args()

    out = args.out or args.stats.with_suffix(args.stats.suffix + ".clean.csv")
    convert_stats_to_cleancsv(args.stats, args.tropo, out, delimiter=args.delimiter)

