#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filière ING3 - PPMD - Traitement de la mesure de phase par ligne de base (GNSS)

Created on Thu Feb  1 16:44:03 2024

@author: Samuel Nahmani (1,2)
Web: https://www.ipgp.fr/annuaire/nahmani/
Contact: nahmani@ipgp.fr | samuel.nahmani@ign.fr

(1) Université Paris Cité, Institut de physique du globe de Paris (IPGP), CNRS, IGN, F-75005 Paris, France
(2) Univ Gustave Eiffel, Géodata Paris, IGN, F-75238 Paris, France

Version: 1.0

Dependencies (Python packages)
-----------------------------
- geodezyx
- pandas
- numpy

Standard library
----------------
- datetime
- pathlib
- os
- gc
"""

# %%
# GeodeZYX toolbox reference
# Sakic, P., Mansur, G., Chaiyaporn, K., & Ballu, V. (2019).
# The geodeZYX toolbox: a versatile Python 3 toolbox for geodetic-oriented purposes (v4.0).
# GFZ Data Services. https://doi.org/10.5880/GFZ.1.1.2019.002
#
# Documentation: https://geodezyx.github.io/geodezyx-toolbox/
# Installation:  pip install git+https://github.com/GeodeZYX/geodezyx-toolbox

# %%
# Standard library
import datetime as dt
import os
import gc
from pathlib import Path

# Third-party
import numpy as np
import pandas as pd

# GeodeZYX
from geodezyx import files_rw     # read/write module
from geodezyx import conv         # conversion module
from geodezyx import operational  # download rinex module
from geodezyx import gnss_edu     # learning module
from geodezyx import reffram      # reference frame / higher geodesy module

# %%
# création du dossier gnss_edu_data qui va contenir les données et les résultats du TP
# see gnss_edu_phase_TP01.py

my_directory = os.environ["HOME"] + "/gnss_edu_data/"

# Chemin avec expansion du ~ vers le home
folder = Path(my_directory).expanduser()

# Création du dossier s'il n'existe pas
folder.mkdir(parents=True, exist_ok=True)

# %%
# Téléchargement automatique des données RINEX des stations de SMNE et MLVL distance d'une
# dizaine de kilomètres dans la région de Paris, France sur le serveur IGN (France)
# données pour 1 jour (2019-176) à 30s
# see gnss_edu_phase_TP01.py

# Création d'un datetime pour gérer le jour de traitement et ne pas à avoir à gérer les doy, les jjul etc !
my_date_to_process = dt.datetime(2019,6,25)

dwl_output_station = operational.download_gnss_rinex(statdico={"rgp" : ["SMNE","MLVL"]},
                                output_dir=my_directory,
                                startdate= my_date_to_process ,
                                enddate= my_date_to_process ,
                                parallel_download = 1) 


# Téléchargement automatique des données ORBIT et CLOCK pour le jour du traitement
# ici  (2019-176) qui correspond à la semaine GPS 2059
dwl_output_satellite = operational.download_gnss_products(archive_dir= my_directory,
                                   startdate= my_date_to_process,
                                   enddate= my_date_to_process,
                                   archtype= 'year/doy',
                                   AC_names=("IGS",),
                                   repro=0,
                                   archive_center="ign",
                                   parallel_download = 1,
                                   ) 
# %%
###############################################################################
# RINEX paths returned by GeodeZYX download utility
#
# GeodeZYX returns a list of tuples: (local_path, success_flag).
# Here we explicitly assign:
#   - BASE  station: SMNE
#   - ROVER station: MLVL
###############################################################################

rinex_base_path, ok_base = dwl_output_station[0]
rinex_rover_path, ok_rover = dwl_output_station[1]

if not ok_base:
    raise RuntimeError(f"Download failed for BASE RINEX: {rinex_base_path}")
if not ok_rover:
    raise RuntimeError(f"Download failed for ROVER RINEX: {rinex_rover_path}")

print("BASE  (SMNE):", rinex_base_path)
print("ROVER (MLVL):", rinex_rover_path)

# Read RINEX observation files
df_base = files_rw.read_rinex_obs(rinex_base_path)
df_rover = files_rw.read_rinex_obs(rinex_rover_path)

print(f"BASE rows : {len(df_base)}")
print(f"ROVER rows: {len(df_rover)}")


# %%
###############################################################################
# Minimal RINEX cleaning (KISS principle)
#
# Educational objective
# ---------------------
# Keep observations as close as possible to raw measurements.
# No frequency combination.
# No unit conversion.
# Carrier phase remains expressed in cycles.
###############################################################################

def clean_rinex_kiss(df):
    """
    Minimal RINEX cleaning following the KISS principle.

    Philosophy
    ----------
    - Keep observations as close as possible to raw RINEX data
    - Do NOT convert carrier phase units
    - Do NOT select observables yet
    - Only remove structurally unusable data

    Operations performed
    --------------------
    1. Convert observable columns to numeric when possible
    2. Remove completely empty columns
    3. Keep GPS satellites only
    4. Remove rows without epoch/prn identifiers
    5. Reset indexing and create a stable row index

    Returns
    -------
    Cleaned pandas DataFrame
    """

    df = df.copy()

    # ------------------------------------------------------------------
    # Convert observable columns to numeric (future-proof pandas)
    # ------------------------------------------------------------------
    protected_cols = {"epoch", "sys", "prn"}

    for col in df.columns:
        if col not in protected_cols:
            try:
                df[col] = pd.to_numeric(df[col])
            except Exception:
                pass

    # ------------------------------------------------------------------
    # Remove columns containing only NaN values
    # ------------------------------------------------------------------
    df = df.loc[:, df.notna().any(axis=0)]

    # ------------------------------------------------------------------
    # Keep GPS satellites only
    # ------------------------------------------------------------------
    if "sys" in df.columns:
        df = df[df["sys"] == "G"]
    elif "prn" in df.columns:
        df = df[df["prn"].astype(str).str.startswith("G")]

    # ------------------------------------------------------------------
    # Remove rows without essential identifiers
    # ------------------------------------------------------------------
    df = df.dropna(subset=["epoch", "prn"])

    # ------------------------------------------------------------------
    # Reset indexing (important before MultiIndex construction)
    # ------------------------------------------------------------------
    df = df.reset_index(drop=True)

    # Stable row index useful for matrix construction later
    df["ind_ligne"] = np.arange(len(df))

    # ------------------------------------------------------------------
    # Diagnostic print (pedagogical)
    # ------------------------------------------------------------------
    empty_cols = df.columns[df.isna().all()]
    if len(empty_cols) > 0:
        print("Warning: remaining empty columns:", list(empty_cols))

    print("KISS cleaning applied.")
    print(f"Remaining observations: {len(df)} rows")

    return df


df_base = clean_rinex_kiss(df_base)
df_rover = clean_rinex_kiss(df_rover)

print("KISS cleaning applied.")
print("Observables kept:", df_base.columns)

# %%
###############################################################################
# GNSS canonical indexing (same structure as TP02)
#
# Each observation is uniquely identified by:
#     (epoch, satellite)
###############################################################################

df_base = df_base.set_index(["epoch", "prn"]).sort_index()
df_rover = df_rover.set_index(["epoch", "prn"]).sort_index()

print("BASE index:", df_base.index.names)
print("ROVER index:", df_rover.index.names)

#%%
# Workspace cleanup
#
# We remove intermediate variables that are no longer needed in order to:
#   - keep the workspace readable,
#   - avoid accidental reuse of temporary data,
#   - reduce memory usage when working with large GNSS datasets.
#
# Note:
# Python normally manages memory automatically. The explicit cleanup below
# simply helps when working interactively with large SP3 and RINEX files.
###############################################################################

for var in [
    "dwl_output_station",
    "dwl_output_satellite",
    "ok_base", "ok_rover",
]:
    if var in globals():
        del globals()[var]

del var
gc.collect()
