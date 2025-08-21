#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
read_vmf1_grid.py — autonome (sans jdutil)
Lecture + interpolation des coefficients VMF1 (ah, aw) aux coordonnées station
et à des dates MJD arbitraires, à partir des fichiers VMFG_YYYYMMDD.HHH.

- Interpolation spatiale:
    * tente scipy.interpolate.griddata (linear)
    * sinon repli nearest neighbor (NumPy) si SciPy indisponible
- Interpolation temporelle: linéaire entre pas 6h
- Chemins: ./Mapping_Fcn/vmf1/VMFG_YYYYMMDD.HHH

Format attendu des fichiers (après les lignes commençant par '!'):
    lat lon ah aw zhd zwd   (au minimum: lat, lon, ah, aw)

Auteur: adapté pour être autonome (remplacement de jdutil)
"""

from __future__ import annotations
import os
from math import floor
from dataclasses import dataclass
from typing import Iterable, Tuple, List
import numpy as np

# ---------------------------------------------------------------------
# Optionnel: SciPy pour une interpolation spatiale "propre".
# Sinon on bascule en nearest neighbor.
# ---------------------------------------------------------------------
try:
    from scipy.interpolate import griddata  # type: ignore
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False


# =======================
#  UTILITAIRES TEMPS
# =======================
def datetime_to_mjd(dt) -> float:
    """Convertit un datetime UTC naïf/aware en MJD (UTC)."""
    from datetime import timezone
    if dt.tzinfo is not None:
        dt = dt.astimezone(timezone.utc).replace(tzinfo=None)
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


def mjd_to_jd(mjd: float) -> float:
    return float(mjd) + 2400000.5


def jd_to_calendar(jd: float) -> Tuple[int, int, int, float]:
    """
    Convertit un Julian Day (float) en (year, month, day, hour_float).
    Implémentation autonome (remplace jdutil.jd_to_date).
    """
    jd = float(jd)
    Z = int(jd + 0.5)
    F = (jd + 0.5) - Z
    if Z < 2299161:
        A = Z
    else:
        alpha = int((Z - 1867216.25)/36524.25)
        A = Z + 1 + alpha - int(alpha/4)
    B = A + 1524
    C = int((B - 122.1)/365.25)
    D = int(365.25 * C)
    E = int((B - D)/30.6001)
    day = B - D - int(30.6001 * E) + F
    if E < 14:
        month = E - 1
    else:
        month = E - 13
    if month > 2:
        year = C - 4716
    else:
        year = C - 4715
    # heure fractionnaire (0–24)
    hour = (day % 1.0) * 24.0
    day_int = int(day)
    return year, month, day_int, hour


def mjd_to_datetime(mjd: float):
    """Conversion MJD -> datetime UTC naïf (sans tzinfo)."""
    from datetime import datetime
    y, m, d, hour = jd_to_calendar(mjd_to_jd(mjd))
    hh = int(hour)
    mmf = (hour - hh) * 60.0
    mm = int(mmf)
    ss = (mmf - mm) * 60.0
    sec = int(ss)
    usec = int(round((ss - sec) * 1e6))
    return datetime(y, m, d, hh, mm, sec, usec)


def floor_to_6h(dt):
    """Ramène dt à l'heure multiple de 6h immédiatement inférieure (UTC)."""
    from datetime import datetime
    return dt.replace(minute=0, second=0, microsecond=0, hour=(dt.hour // 6) * 6)


def add_hours(dt, hours: int):
    from datetime import timedelta
    return dt + timedelta(hours=hours)


# =======================
#  DONNÉES & FICHIERS
# =======================
@dataclass
class VMFPoint:
    lat: float
    lon: float
    ah: float
    aw: float


def _vmf_dir() -> str:
    """Répertoire où se trouvent les VMFG_YYYYMMDD.HHH."""
    return os.path.join(".", "Mapping_Fcn", "vmf1")


def _vmf_filename(dt_6h) -> str:
    """Nom de fichier VMFG_YYYYMMDD.HHH pour un datetime aligné 6h."""
    y = dt_6h.year
    m = f"{dt_6h.month:02d}"
    d = f"{dt_6h.day:02d}"
    h = f"{dt_6h.hour:02d}"
    return f"VMFG_{y}{m}{d}.H{h}"


def _read_vmf_file(filepath: str) -> List[VMFPoint]:
    """
    Lit un fichier VMFG_* et retourne une liste de VMFPoint.
    Ignore les lignes commençant par '!'.
    Attend au minimum: lat lon ah aw ...
    """
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
                lat = float(parts[0])
                lon = float(parts[1])
                ah = float(parts[2])
                aw = float(parts[3])
                pts.append(VMFPoint(lat=lat, lon=lon, ah=ah, aw=aw))
            except Exception:
                # ligne malformée -> ignorer
                continue
    if not pts:
        raise ValueError(f"No valid data rows in {os.path.basename(filepath)}")
    return pts


def _interp_spatial(points: List[VMFPoint], lat0: float, lon0: float) -> Tuple[float, float]:
    """
    Interpolation spatiale (ah, aw) au point (lat0, lon0).
    - Si SciPy dispo: griddata 'linear' sur les points (lat,lon).
    - Sinon: nearest neighbor (NumPy pur).
    """
    lats = np.array([p.lat for p in points], dtype=float)
    lons = np.array([p.lon for p in points], dtype=float)
    ahs = np.array([p.ah for p in points], dtype=float)
    aws = np.array([p.aw for p in points], dtype=float)

    if _HAS_SCIPY:
        pts = np.column_stack([lats, lons])
        ah_i = griddata(pts, ahs, (lat0, lon0), method="linear")
        aw_i = griddata(pts, aws, (lat0, lon0), method="linear")
        # Si en dehors de la convex hull -> fallback nearest
        if (ah_i is None) or np.isnan(ah_i):
            ah_i = griddata(pts, ahs, (lat0, lon0), method="nearest")
        if (aw_i is None) or np.isnan(aw_i):
            aw_i = griddata(pts, aws, (lat0, lon0), method="nearest")
        return float(ah_i), float(aw_i)
    else:
        # Fallback: nearest neighbor simple
        d2 = (lats - lat0)**2 + (lons - lon0)**2
        idx = int(np.argmin(d2))
        return float(ahs[idx]), float(aws[idx])


def _coeffs_from_file_for_station(dt_6h, lat_deg: float, lon_deg: float) -> Tuple[float, float]:
    """
    Lit le fichier VMFG correspondant à dt_6h et interpole à (lat, lon).
    """
    path = os.path.join(_vmf_dir(), _vmf_filename(dt_6h))
    points = _read_vmf_file(path)
    return _interp_spatial(points, lat_deg, lon_deg)


# =======================
#  API PRINCIPALE
# =======================
def read_vmf1_grid(mjd: Iterable[float], ell: Iterable[float]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Retourne (ah_vec, aw_vec) pour chaque MJD donné et une station ell=[lat_deg, lon_deg, h_m].

    Paramètres
    ----------
    mjd : Iterable[float]
        Liste/array de dates en Modified Julian Day (UTC).
    ell : [lat_deg, lon_deg, h_m]
        Latitude (deg), Longitude (deg), Hauteur ellipsoïdale (m). (La hauteur n'est
        pas utilisée ici pour ah/aw, mais on garde la signature d’interface.)

    Renvoie
    -------
    ah_vec : np.ndarray shape (N,)
    aw_vec : np.ndarray shape (N,)

    Notes
    -----
    - Les fichiers requis sont:
        ./Mapping_Fcn/vmf1/VMFG_YYYYMMDD.H{00,06,12,18}
      pour les 6h entourant chaque MJD.
    - Interpolation temporelle: linéaire entre t0=⌊MJD⌋_6h et t1=t0+6h.
    """
    mjd_arr = np.atleast_1d(np.array(list(mjd), dtype=float))
    lat_deg, lon_deg, _h_m = float(ell[0]), float(ell[1]), float(ell[2])

    # Pré-chargement des fichiers nécessaires (évite relecture multiple)
    # On boucle d’abord pour recenser tous les pas 6h requis.
    from datetime import datetime
    needed_slots = set()
    dt_by_slot = {}

    for m in mjd_arr:
        dt = mjd_to_datetime(m)
        t0 = floor_to_6h(dt)
        t1 = add_hours(t0, 6)
        needed_slots.add((t0.year, t0.month, t0.day, t0.hour))
        needed_slots.add((t1.year, t1.month, t1.day, t1.hour))
        dt_by_slot[(t0.year, t0.month, t0.day, t0.hour)] = t0
        dt_by_slot[(t1.year, t1.month, t1.day, t1.hour)] = t1

    # Cache des (ah,aw) spatialisés aux pas 6h
    slot_cache = {}

    def get_slot_ahaw(y, mo, d, h):
        key = (y, mo, d, h)
        if key in slot_cache:
            return slot_cache[key]
        dt6 = dt_by_slot[key]
        ah, aw = _coeffs_from_file_for_station(dt6, lat_deg, lon_deg)
        slot_cache[key] = (ah, aw)
        return ah, aw

    # Interpolation temporelle
    ah_out = np.empty_like(mjd_arr, dtype=float)
    aw_out = np.empty_like(mjd_arr, dtype=float)

    for i, m in enumerate(mjd_arr):
        dt = mjd_to_datetime(m)
        t0 = floor_to_6h(dt)
        t1 = add_hours(t0, 6)
        m0 = datetime_to_mjd(t0)
        m1 = datetime_to_mjd(t1)

        ah0, aw0 = get_slot_ahaw(t0.year, t0.month, t0.day, t0.hour)
        ah1, aw1 = get_slot_ahaw(t1.year, t1.month, t1.day, t1.hour)

        # coefficient d'interpolation (linéaire)
        if m1 == m0:
            w = 0.0
        else:
            w = (m - m0) / (m1 - m0)
            if w < 0.0: w = 0.0
            if w > 1.0: w = 1.0

        ah_out[i] = (1.0 - w) * ah0 + w * ah1
        aw_out[i] = (1.0 - w) * aw0 + w * aw1

    return ah_out, aw_out


# =======================
#  TEST RAPIDE (optionnel)
# =======================
if __name__ == "__main__":
    # Petit test sur une date fictive (nécessite fichiers présents)
    # Station: Paris (approx.)
    ell = [48.85, 2.35, 100.0]
    # Exemple MJD: 2021-07-10T01:00:18
    # ->  Converti ici via nos utilitaires pour test
    from datetime import datetime, timezone
    dt = datetime(2021, 7, 10, 1, 0, 18, tzinfo=timezone.utc).replace(tzinfo=None)
    mjd_val = datetime_to_mjd(dt)
    try:
        ah, aw = read_vmf1_grid([mjd_val], ell)
        print("ah =", ah, "aw =", aw)
    except Exception as e:
        print("Test lecture échoué (fichiers manquants ?):", e)

