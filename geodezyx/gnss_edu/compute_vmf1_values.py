#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Télécharge (si besoin) les grilles VMF1 puis calcule les mapping functions VMF1
pour DUNQ00FRA au 2021-07-10T01:00:18, en utilisant:
  - download_VMF.download_VMF (ton script de téléchargement)
  - read_vmf1_grid.read_vmf1_grid
  - vmf1_ht.vmf1_ht
Exécuter depuis la racine du projet.
"""

import os
from math import floor
from datetime import datetime, timezone, timedelta
import numpy as np

# --- briques du projet ---
import download_VMF as dl
import read_vmf1_grid as rvg
import vmf1_ht

# --- Patch de portabilité: read_vmf1_grid construit des chemins avec '\\'.
#     On monkey-patch "open" dans le module pour convertir en séparateur OS.
import builtins
rvg.open = lambda path, *args, **kwargs: builtins.open(path.replace("\\", os.sep), *args, **kwargs)

# ---------------- outils temps ----------------
def iso_to_mjd(iso_utc: str) -> float:
    """Convertit une date ISO UTC en MJD (UTC)."""
    dt = datetime.fromisoformat(iso_utc)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    else:
        dt = dt.astimezone(timezone.utc)
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
    """Retourne l'ensemble des dates (YYYY,MM,DD) impliquées par l'interpolation 6h autour de dt."""
    base = dt_utc.replace(minute=0, second=0, microsecond=0)
    # borne basse (multiple de 6h)
    floor6 = base.replace(hour=(base.hour//6)*6)
    ceil6  = floor6 + timedelta(hours=6)
    return {
        (floor6.year, floor6.month, floor6.day),
        (ceil6.year,  ceil6.month,  ceil6.day)
    }

# ---------------- téléchargement si besoin ----------------
def ensure_vmf1_files_for_epoch(dt_utc: datetime, vmf_dir: str):
    """
    Vérifie la présence des grilles nécessaires dans vmf_dir.
    Si manquantes, appelle download_VMF.download_VMF(year, month, day) pour chaque date utile.
    """
    os.makedirs(vmf_dir, exist_ok=True)

    needed_dates = six_hour_bracketing_dates(dt_utc)  # au pire 2 dates (avant/après)
    # On considère que la fonction de lecture peut avoir besoin des H00/H06/H12/H18 de ces dates.
    # On télécharge donc la journée complète pour chaque date impliquée (c'est ainsi que marche ton downloader).
    for (y, m, d) in sorted(needed_dates):
        # fichiers attendus
        wanted = [
            f"VMFG_{y}{str(m).zfill(2)}{str(d).zfill(2)}.H{hh}"
            for hh in ("00", "06", "12", "18")
        ]
        have = all(os.path.exists(os.path.join(vmf_dir, w)) for w in wanted)
        if not have:
            print(f"[info] Grilles VMF1 manquantes pour {y}-{str(m).zfill(2)}-{str(d).zfill(2)}. "
                  f"Téléchargement…")
            # Le downloader télécharge dans le répertoire courant: on se place donc dans vmf_dir le temps du DL.
            cwd = os.getcwd()
            try:
                os.chdir(vmf_dir)
                # Très important: ton download_VMF attend un month *tel quel* dans le nom du fichier.
                # On lui passe une chaîne zero-padded pour garantir 'YYYYMMDD'.
                dl.download_VMF(year=str(y), month=str(m).zfill(2), day=d)
            finally:
                os.chdir(cwd)
        else:
            print(f"[ok] Grilles VMF1 déjà présentes pour {y}-{str(m).zfill(2)}-{str(d).zfill(2)}.")

# ---------------- paramètres station & epoch ----------------
station = {
    "name": "DUNQ00FRA",
    "lon_deg": 2.366699,
    "lat_deg": 51.048100,
    "h_m": 51.111900,
}
epoch_iso_utc = "2021-07-10T01:00:18"  # UTC
elevations_deg = [5, 10, 15, 20, 30, 45, 60]  # à adapter si besoin

def main():
    # chemins
    project_root = os.getcwd()
    vmf_dir = os.path.join(project_root, "Mapping_Fcn", "vmf1")

    # 1) S'assurer que les grilles nécessaires sont présentes (sinon les télécharger)
    dt = datetime.fromisoformat(epoch_iso_utc).replace(tzinfo=timezone.utc)
    ensure_vmf1_files_for_epoch(dt, vmf_dir)

    # 2) Temps → MJD
    mjd = iso_to_mjd(epoch_iso_utc)

    # 3) Lire coefficients a_h, a_w à la station/datetime
    ah_vec, aw_vec = rvg.read_vmf1_grid(
        mjd=[mjd],
        ell=[station["lat_deg"], station["lon_deg"], station["h_m"]]
    )
    ah = float(np.atleast_1d(ah_vec).ravel()[0])
    aw = float(np.atleast_1d(aw_vec).ravel()[0])

    # 4) Calculer VMF1 (hydro & humide) pour une liste d’élévations
    dlat_rad = np.deg2rad(station["lat_deg"])
    ht_m = station["h_m"]

    print(f"# Station   : {station['name']}")
    print(f"# Lat, Lon  : {station['lat_deg']:.6f}°, {station['lon_deg']:.6f}°")
    print(f"# Height    : {ht_m:.3f} m")
    print(f"# Epoch UTC : {epoch_iso_utc}  (MJD {mjd:.6f})\n")
    print(f"ah (grid) = {ah:.6e}")
    print(f"aw (grid) = {aw:.6e}\n")

    print("Elevation(deg),  VMF1_hydro,  VMF1_wet")
    for el_deg in elevations_deg:
        zd_rad = np.deg2rad(90.0 - el_deg)  # distance zénithale
        vmf1h, vmf1w = vmf1_ht.vmf1_ht(
            ah=ah, aw=aw, dmjd=mjd, dlat=dlat_rad, ht=ht_m, zd=zd_rad
        )
        print(f"{el_deg:>6.1f}         {vmf1h:>10.6f}  {vmf1w:>10.6f}")

if __name__ == "__main__":
    main()

