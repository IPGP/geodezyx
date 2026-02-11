#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Convertit un fichier statistiques GINS en CSV propre avec dates ISO.

Entrées :
  - Fichier statistiques, ex.:
      SPOTGINS_DUNQ00FRA_2021_191_26123__.yml.250819_135206.1
    (le nom encode station + année + DOY + (souvent) JUL50)
  - Fichier tropo correspondant, ex.:
      SPOTGINS_DUNQ00FRA.tropo
    (l’entête contient Latitude / Longitude / Height)

Sortie :
  - CSV avec en-tête : station, lat_deg, lon_deg, height_m, datetime_utc, mjd, sod, jcn, type, sat, az_deg, el_deg, residual, weight, extra_fields

Remarque :
  - Le parseur “type 11” (GNSS non différencié) suit le format standard montré :
      0:type 1:date 2:poids 3:résidu 4:sat 5:lat_sat 6:lon_sat 7:code_sta 8:az 9:el ...
    Les autres types sont passés en "best-effort" (si az/el détectables).
"""

import re
import csv
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, Tuple, List
from datetime import datetime, timezone, timedelta
from math import floor

# ----------------- Utilitaires temps -----------------
def parse_year_doy_from_path(path: str) -> Tuple[int, int]:
    """
    Extrait (YYYY, DOY) du nom, ex: _2021_191_ .
    """
    m = re.search(r'[_\-.](20\d{2}|19\d{2})[_\-.](\d{3})(?:[_\-.]|$)', path)
    if not m:
        raise ValueError("Impossible d'extraire YYYY_DDD du nom du fichier statistiques.")
    year = int(m.group(1)); doy = int(m.group(2))
    if not (1 <= doy <= 366):
        raise ValueError(f"DOY invalide: {doy}")
    return year, doy

def doy_to_date(year: int, doy: int) -> datetime:
    jan1 = datetime(year, 1, 1, tzinfo=timezone.utc)
    return jan1 + timedelta(days=doy - 1)

def mjd_to_datetime(mjd: float) -> datetime:
    jd = mjd + 2400000.5
    Z = int(jd + 0.5)
    F = (jd + 0.5) - Z
    if Z < 2299161:
        A = Z
    else:
        alpha = int((Z - 1867216.25) / 36524.25)
        A = Z + 1 + alpha - int(alpha / 4)
    B = A + 1524
    C = int((B - 122.1) / 365.25)
    D = int(365.25 * C)
    E = int((B - D) / 30.6001)
    day = B - D - int(30.6001 * E) + F
    month = (E - 1) if E < 14 else (E - 13)
    year = (C - 4716) if month > 2 else (C - 4715)
    day_whole = int(day)
    frac = day - day_whole
    seconds = frac * 86400.0
    hh = int(seconds // 3600)
    mm = int((seconds % 3600) // 60)
    ss = int(seconds % 60)
    us = int(round((seconds - (hh * 3600 + mm * 60 + ss)) * 1e6))
    return datetime(year, month, day_whole, hh, mm, ss, us, tzinfo=timezone.utc)

def datetime_to_mjd(dt: datetime) -> float:
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    else:
        dt = dt.astimezone(timezone.utc)
    Y, M, D = dt.year, dt.month, dt.day
    frac = (dt.hour + dt.minute/60 + dt.second/3600 + dt.microsecond/3.6e9) / 24.0
    if M <= 2:
        Yp, Mp = Y - 1, M + 12
    else:
        Yp, Mp = Y, M
    A = floor(Yp/100)
    B = 2 - A + floor(A/4)
    JD = floor(365.25*(Yp+4716)) + floor(30.6001*(Mp+1)) + (D + frac) + B - 1524.5
    return JD - 2400000.5

def detect_and_to_datetime(x: float, year_hint: Optional[int], doy_hint: Optional[int]) -> Tuple[datetime, Optional[float], Optional[float]]:
    """
    Détecte SoD / JCN(JUL50) / MJD / JD et renvoie (datetime, sod, jcn)
    - SoD : [0..86400]
    - JCN : ~[18000..40000] (jours depuis 1950-01-01 00:00:00)
    - MJD : [50000..80000]
    - JD  : [2.4e6..2.6e6]
    """
    # SoD ?
    if 0.0 <= x <= 86400.0:
        if year_hint is None or doy_hint is None:
            raise ValueError("SoD détecté mais YYYY/DOY inconnus (nécessaires).")
        base = doy_to_date(year_hint, doy_hint)
        dt = base + timedelta(seconds=x)
        return dt, float(x), None

    # JCN (JUL50) ?
    if 18000.0 <= x <= 40000.0:
        base = datetime(1950, 1, 1, tzinfo=timezone.utc)
        dt = base + timedelta(days=x)
        return dt, None, float(x)

    # MJD ?
    if 50000.0 <= x <= 80000.0:
        dt = mjd_to_datetime(x)
        return dt, None, None

    # JD ?
    if 2_400_000.0 <= x <= 2_600_000.0:
        dt = mjd_to_datetime(x - 2400000.5)
        return dt, None, None

    # Dernier recours : si on a YYYY/DOY, essayer SoD (formats exotiques)
    if (year_hint is not None and doy_hint is not None) and (x < 200000.0):
        base = doy_to_date(year_hint, doy_hint)
        dt = base + timedelta(seconds=x)
        return dt, float(x), None

    raise ValueError(f"Impossible d'interpréter le champ Date = {x}")

# ----------------- Métadonnées station (tropo) -----------------
@dataclass
class StationMeta:
    name: str
    lat: float
    lon: float
    h: float

def read_station_meta_from_tropo(tropo_path: Path) -> StationMeta:
    name = None; lat = None; lon = None; h = None
    with tropo_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if s.startswith("# STATION"):
                # "# STATION          : DUNQ00FRA"
                m = re.search(r':\s*(\S+)', s)
                if m: name = m.group(1)
            elif s.startswith("# Latitude"):
                m = re.search(r':\s*([+-]?\d+(?:\.\d+)?)', s)
                if m: lat = float(m.group(1))
            elif s.startswith("# Longitude"):
                m = re.search(r':\s*([+-]?\d+(?:\.\d+)?)', s)
                if m: lon = float(m.group(1))
            elif s.startswith("# Height"):
                m = re.search(r':\s*([+-]?\d+(?:\.\d+)?)', s)
                if m: h = float(m.group(1))
            if name and (lat is not None) and (lon is not None) and (h is not None):
                break
    if name is None or lat is None or lon is None or h is None:
        raise ValueError(f"Impossible de lire les métadonnées station dans {tropo_path}")
    return StationMeta(name=name, lat=lat, lon=lon, h=h)

# ----------------- Parsing d'une ligne stats -----------------
@dataclass
class ParsedLine:
    dt_iso: str
    mjd: float
    sod: Optional[float]
    jcn: Optional[float]
    typ: int
    sat: Optional[str]
    az: Optional[float]
    el: Optional[float]
    residual: Optional[float]
    weight: Optional[float]
    extra: str  # tout le reste, pour traçabilité

def parse_stats_line(line: str, year_hint: int, doy_hint: int) -> Optional[ParsedLine]:
    s = line.strip()
    if not s or s.startswith("#"):  # ignorer commentaires/vides
        return None
    parts = s.split()
    # type
    try:
        typ = int(parts[0])
    except Exception:
        return None

    # date (champ e16.10)
    try:
        date_field = float(parts[1])
    except Exception:
        return None
    dt, sod, jcn = detect_and_to_datetime(date_field, year_hint, doy_hint)
    mjd = datetime_to_mjd(dt)

    # poids, résidu si dispo
    weight = None; residual = None
    sat = None; az = None; el = None

    # best-effort selon format standard vu pour type 11
    # indices attendus: 0:type 1:date 2:poids 3:résidu 4:sat 5:lat_sat 6:lon_sat 7:code_sta 8:az 9:el ...
    if len(parts) >= 10:
        # poids
        try:
            weight = float(parts[2])
        except Exception:
            weight = None
        # résidu
        try:
            residual = float(parts[3])
        except Exception:
            residual = None
        # sat
        sat = parts[4]
        # az/el
        try:
            az_cand = float(parts[8])
            el_cand = float(parts[9])
            if 0.0 <= az_cand <= 360.0 and 0.0 <= el_cand <= 90.0:
                az, el = az_cand, el_cand
        except Exception:
            az = el = None

    # extra (trace brute hors colonnes standard)
    extra = " ".join(parts[10:]) if len(parts) > 10 else ""

    return ParsedLine(
        dt_iso=dt.isoformat(),
        mjd=mjd,
        sod=sod,
        jcn=jcn,
        typ=typ,
        sat=sat,
        az=az,
        el=el,
        residual=residual,
        weight=weight,
        extra=extra
    )

# ----------------- Programme principal -----------------
def convert_stats_to_csv(stats_path: Path, tropo_path: Path, out_csv: Path):
    year, doy = parse_year_doy_from_path(stats_path.name)
    meta = read_station_meta_from_tropo(tropo_path)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        # en-tête riche et stable
        w.writerow([
            "station","lat_deg","lon_deg","height_m",
            "datetime_utc","mjd","sod","jcn",
            "type","sat","az_deg","el_deg",
            "residual","weight","extra_fields"
        ])
        # parse lignes
        n_in = 0; n_out = 0
        with stats_path.open("r", encoding="utf-8", errors="ignore") as fin:
            for ln in fin:
                n_in += 1
                rec = parse_stats_line(ln, year, doy)
                if rec is None:
                    continue
                w.writerow([
                    meta.name, f"{meta.lat:.6f}", f"{meta.lon:.6f}", f"{meta.h:.3f}",
                    rec.dt_iso, f"{rec.mjd:.6f}",
                    (f"{rec.sod:.3f}" if rec.sod is not None else ""),
                    (f"{rec.jcn:.6f}" if rec.jcn is not None else ""),
                    rec.typ, (rec.sat or ""),
                    (f"{rec.az:.3f}" if rec.az is not None else ""),
                    (f"{rec.el:.3f}" if rec.el is not None else ""),
                    (f"{rec.residual:.9f}" if rec.residual is not None else ""),
                    (f"{rec.weight:.3f}" if rec.weight is not None else ""),
                    rec.extra
                ])
                n_out += 1

    print(f"[info] Station   : {meta.name}  (lat={meta.lat:.6f}, lon={meta.lon:.6f}, h={meta.h:.3f} m)")
    print(f"[info] Day hint  : {year}-DOY{doy}")
    print(f"[info] Parsed    : {n_out} / {n_in} lignes → {out_csv}")

# -------------- CLI --------------
if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(
        description="Convertit un fichier statistiques GINS en CSV propre (dates ISO, MJD, SoD/JCN, colonnes principales)."
    )
    p.add_argument("--stats", required=True, type=Path, help="Chemin du fichier statistiques GINS")
    p.add_argument("--tropo", required=True, type=Path, help="Chemin du fichier SPOTGINS_<STA>.tropo (pour métadonnées station)")
    p.add_argument("--out", type=Path, help="CSV de sortie (défaut: <stats>.csv)")
    args = p.parse_args()

    out = args.out or args.stats.with_suffix(args.stats.suffix + ".csv")
    convert_stats_to_csv(args.stats, args.tropo, out)

