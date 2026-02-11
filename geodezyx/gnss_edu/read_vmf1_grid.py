#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
read_vmf1_grid.py — autonome (sans jdutil)
Lecture + interpolation des coefficients/grilles VMF1 aux coordonnées station
et à des dates MJD arbitraires, à partir des fichiers VMFG_YYYYMMDD.HH.

- Interpolation spatiale:
    * tente scipy.interpolate.griddata (linear)
    * sinon repli nearest neighbor (NumPy) si SciPy indisponible
- Interpolation temporelle: linéaire entre pas 6h (t0, t1=t0+6h)
- Chemins: ./Mapping_Fcn/vmf1/VMFG_YYYYMMDD.HH

Format attendu (après lignes commençant par '!'):
    lat lon ah aw [zhd zwd]    ← zhd/zwd optionnels
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Iterable, List, Tuple
import os
import math

# SciPy optionnelle
try:
    import numpy as np
    _HAS_NUMPY = True
except Exception:
    _HAS_NUMPY = False

try:
    from scipy.interpolate import griddata  # type: ignore
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False


# ----------- Conversions de dates (autonomes) -----------

def calendar_to_jd(year: int, month: int, day: int, hour: float = 0.0) -> float:
    """Calendrier → Julian Day (double)."""
    Y, M = int(year), int(month)
    D = float(day)
    frac_day = float(hour) / 24.0
    if M <= 2:
        Yp, Mp = Y - 1, M + 12
    else:
        Yp, Mp = Y, M
    A = math.floor(Yp / 100)
    B = 2 - A + math.floor(A / 4)
    JD = math.floor(365.25 * (Yp + 4716)) + math.floor(30.6001 * (Mp + 1)) + (D + frac_day) + B - 1524.5
    return float(JD)


def jd_to_mjd(jd: float) -> float:
    return float(jd) - 2400000.5


def mjd_to_jd(mjd: float) -> float:
    return float(mjd) + 2400000.5


def jd_to_calendar(jd: float) -> Tuple[int, int, int, float]:
    """Julian Day → (year, month, day, hour_float)."""
    jd = float(jd)
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
    month = E - 1 if E < 14 else E - 13
    year = C - 4716 if month > 2 else C - 4715

    day_int = int(day)
    frac = day - day_int
    hour = frac * 24.0
    return int(year), int(month), int(day_int), float(hour)


def datetime_to_mjd(dt) -> float:
    """datetime naïf (UTC) → MJD."""
    y, m, d = dt.year, dt.month, dt.day
    h = dt.hour + dt.minute / 60.0 + (dt.second + dt.microsecond / 1e6) / 3600.0
    return jd_to_mjd(calendar_to_jd(y, m, d, h))


def mjd_to_datetime(mjd: float):
    """MJD → datetime naïf (UTC)."""
    from datetime import datetime
    y, m, d, h = jd_to_calendar(mjd_to_jd(mjd))
    hour = int(h)
    minute = int((h - hour) * 60.0)
    sec = int(round((((h - hour) * 60.0) - minute) * 60.0))
    if sec == 60:
        sec = 59
    return datetime(int(y), int(m), int(d), hour, minute, sec)


def floor_to_6h(dt):
    """Aligne un datetime à l'heure 0/6/12/18 la plus proche ≤ dt."""
    from datetime import datetime
    h = (dt.hour // 6) * 6
    return datetime(dt.year, dt.month, dt.day, h, 0, 0)


def add_hours(dt, hrs):
    from datetime import timedelta
    return dt + timedelta(hours=int(hrs))


# ----------- Lecture grilles VMF -----------

@dataclass
class VMFPoint:
    lat: float
    lon: float
    ah: float
    aw: float
    zhd: float | None = None
    zwd: float | None = None


def _vmf_dir() -> str:
    """Répertoire où se trouvent les VMFG_YYYYMMDD.HH."""
    return os.path.join(".", "Mapping_Fcn", "vmf1")


def _vmf_filename(dt_6h) -> str:
    """Nom de fichier VMFG_YYYYMMDD.HH pour un datetime aligné 6h."""
    y = dt_6h.year
    m = f"{dt_6h.month:02d}"
    d = f"{dt_6h.day:02d}"
    h = f"{dt_6h.hour:02d}"
    return f"VMFG_{y}{m}{d}.H{h}"


def _read_vmf_file(filepath: str) -> List[VMFPoint]:
    """Lit un VMFG_* et retourne une liste de VMFPoint (ah, aw, zhd?, zwd?)."""
    pts: List[VMFPoint] = []
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"VMF grid file not found: {filepath}")

    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("!"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                lat = float(parts[0]); lon = float(parts[1])
                ah  = float(parts[2]); aw  = float(parts[3])
                zhd = float(parts[4]) if len(parts) >= 5 else float("nan")
                zwd = float(parts[5]) if len(parts) >= 6 else float("nan")
                pts.append(VMFPoint(lat=lat, lon=lon, ah=ah, aw=aw, zhd=zhd, zwd=zwd))
            except Exception:
                continue
    if not pts:
        raise ValueError(f"No valid data rows in {os.path.basename(filepath)}")
    return pts


def _interp_spatial(points: List[VMFPoint], lat0: float, lon0: float) -> tuple[float, float, float, float]:
    """Interpolation spatiale (ah, aw, zhd, zwd) aux coordonnées demandées."""
    if not _HAS_NUMPY:
        raise RuntimeError("NumPy requis pour l'interpolation.")

    import numpy as np
    lats = np.array([p.lat for p in points], dtype=float)
    lons = np.array([p.lon for p in points], dtype=float)
    ahs  = np.array([p.ah  for p in points], dtype=float)
    aws  = np.array([p.aw  for p in points], dtype=float)
    zhds = np.array([p.zhd for p in points], dtype=float)
    zwds = np.array([p.zwd for p in points], dtype=float)

    def interp(arr):
        if _HAS_SCIPY:
            pts = np.column_stack([lats, lons])
            val = griddata(pts, arr, (lat0, lon0), method="linear")  # type: ignore
            if (val is None) or np.isnan(val):
                val = griddata(pts, arr, (lat0, lon0), method="nearest")  # type: ignore
            return float(val)
        # repli: plus proche voisin
        d2 = (lats - lat0)**2 + (lons - lon0)**2
        return float(arr[int(np.argmin(d2))])

    ah_i  = interp(ahs)
    aw_i  = interp(aws)
    zhd_i = interp(zhds) if np.isfinite(zhds).any() else float("nan")
    zwd_i = interp(zwds) if np.isfinite(zwds).any() else float("nan")
    return ah_i, aw_i, zhd_i, zwd_i


def _coeffs_from_file_for_station(dt_6h, lat_deg: float, lon_deg: float) -> tuple[float, float, float, float]:
    """Charge le VMFG de t=dt_6h et interpole (ah, aw, zhd, zwd) à (lat, lon)."""
    path = os.path.join(_vmf_dir(), _vmf_filename(dt_6h))
    points = _read_vmf_file(path)
    return _interp_spatial(points, lat_deg, lon_deg)


# ----------- API publique -----------

def read_vmf1_grid_ex(mjd: Iterable[float], ell: Iterable[float]):
    """
    Version étendue: renvoie (ah, aw, zhd, zwd), interpolés en temps et espace.
    - ell = (lat_deg, lon_deg, h_m)  # h_m non utilisé ici
    - Interpolation temporelle: linéaire entre t0=floor_6h(dt), t1=t0+6h
    """
    if not _HAS_NUMPY:
        raise RuntimeError("NumPy requis (installé automatiquement avec SciPy).")

    import numpy as np
    mjd_arr = np.atleast_1d(np.array(list(mjd), dtype=float))
    lat_deg, lon_deg, _h_m = float(ell[0]), float(ell[1]), float(ell[2])

    # Recenser les pas 6h nécessaires
    from datetime import datetime
    needed_slots = set()
    dt_by_slot = {}
    for m in mjd_arr:
        dt = mjd_to_datetime(m)
        t0 = floor_to_6h(dt)
        t1 = add_hours(t0, 6)
        for t in (t0, t1):
            key = (t.year, t.month, t.day, t.hour)
            needed_slots.add(key)
            dt_by_slot[key] = t

    # Cache pas 6h → (ah, aw, zhd, zwd)
    slot_cache: dict[tuple[int,int,int,int], tuple[float,float,float,float]] = {}

    def get_slot(y, mo, d, h):
        key = (y, mo, d, h)
        if key not in slot_cache:
            slot_cache[key] = _coeffs_from_file_for_station(dt_by_slot[key], lat_deg, lon_deg)
        return slot_cache[key]

    ah_out = np.empty_like(mjd_arr)
    aw_out = np.empty_like(mjd_arr)
    zhd_out = np.empty_like(mjd_arr)
    zwd_out = np.empty_like(mjd_arr)

    for i, m in enumerate(mjd_arr):
        dt = mjd_to_datetime(m)
        t0 = floor_to_6h(dt)
        t1 = add_hours(t0, 6)
        m0 = datetime_to_mjd(t0); m1 = datetime_to_mjd(t1)

        ah0, aw0, zhd0, zwd0 = get_slot(t0.year, t0.month, t0.day, t0.hour)
        ah1, aw1, zhd1, zwd1 = get_slot(t1.year, t1.month, t1.day, t1.hour)

        w = 0.0 if m1 == m0 else min(1.0, max(0.0, (m - m0) / (m1 - m0)))

        ah_out[i]  = (1 - w) * ah0  + w * ah1
        aw_out[i]  = (1 - w) * aw0  + w * aw1

        # zhd/zwd: si une des extrémités est NaN on renvoie NaN
        if _np_all_finite([zhd0, zhd1]):
            zhd_out[i] = (1 - w) * zhd0 + w * zhd1
        else:
            zhd_out[i] = float("nan")

        if _np_all_finite([zwd0, zwd1]):
            zwd_out[i] = (1 - w) * zwd0 + w * zwd1
        else:
            zwd_out[i] = float("nan")

    return ah_out, aw_out, zhd_out, zwd_out


def read_vmf1_grid(mjd: Iterable[float], ell: Iterable[float]):
    """
    API historique: renvoie uniquement (ah, aw).
    Reste 100% compatible avec l’existant.
    """
    ah, aw, _zhd, _zwd = read_vmf1_grid_ex(mjd, ell)
    return ah, aw


# ----------- utilitaires internes -----------

def _np_all_finite(x) -> bool:
    if not _HAS_NUMPY:
        return all((not (isinstance(v, float) and (math.isnan(v) or math.isinf(v)))) for v in x)
    import numpy as np
    return np.isfinite(np.array(x, dtype=float)).all()


# ----------- test manuel -----------

if __name__ == "__main__":
    # Petit test (nécessite que les fichiers VMF soient présents):
    # Station: COCO (approx.)
    ell = [-12.192, 96.834, 10.0]
    from datetime import datetime, timezone
    dt = datetime(2020, 9, 24, 12, 0, tzinfo=timezone.utc).replace(tzinfo=None)
    mjd_val = datetime_to_mjd(dt)
    try:
        ah, aw, zhd, zwd = read_vmf1_grid_ex([mjd_val], ell)
        print("ah =", ah, "\naw =", aw, "\nzhd =", zhd, "\nzwd =", zwd)
        # Exemple de slant delay si élévation e connue: SD ≈ ah*zhd + aw*zwd
    except Exception as e:
        print("Test échoué (fichiers manquants ?):", e)

