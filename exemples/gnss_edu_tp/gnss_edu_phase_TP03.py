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
# Chargement du fichier sp3 dans un dataframe
# SP3 (load once)
sp3_path = dwl_output_satellite[0][0] if isinstance(dwl_output_satellite[0], tuple) else dwl_output_satellite[0]
print(f"Using SP3 file: {sp3_path}")
df_sp3 = files_rw.read_sp3(sp3_path, returns_pandas=True, new_col_names=True)

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


# %%
###############################################################################
# Synchronization of BASE and ROVER observations
#
# Educational objective
# ---------------------
# Differential GNSS processing requires both stations to observe
# the SAME satellite at the SAME epoch.
#
# We therefore keep only common (epoch, satellite) observations.
###############################################################################

print("***** Synchronizing BASE and ROVER observations *****")

# Intersection of MultiIndex
common_index = df_base.index.intersection(df_rover.index)

print(f"Common observations: {len(common_index)}")

# Keep synchronized observations
df_base_sync  = df_base.loc[common_index].copy()
df_rover_sync = df_rover.loc[common_index].copy()

print("BASE synchronized rows :", len(df_base_sync))
print("ROVER synchronized rows:", len(df_rover_sync))


# %%
###############################################################################
# Satellite states computed independently for BASE and ROVER (TP02 reuse)
#
# Educational objective
# ---------------------
# Compute satellite positions at emission time separately for BASE and ROVER.
# This allows us to quantify the (small) difference in satellite position due
# to the slightly different signal travel times.
###############################################################################

def add_satellite_state_from_sp3(df_rnx_in, df_sp3, label=""):
    """
    Reuse TP02 logic: compute satellite position at emission time + clk + relativistic term.
    Inputs
    ------
    df_rnx_in : DataFrame with columns at least ['epoch','prn','C1']
    df_sp3    : SP3 DataFrame with columns at least ['epoch','prn','x','y','z','clk'] (x,y,z in km)
    label     : used only for prints

    Output
    ------
    df_out : same DF with columns added:
        X_sat, Y_sat, Z_sat (m), dte_sat (s), dRelat (s)
    """
    df_rnx = df_rnx_in.copy()

    # Ensure we have columns (not MultiIndex)
    if isinstance(df_rnx.index, pd.MultiIndex) and set(df_rnx.index.names) >= {"epoch", "prn"}:
        df_rnx = df_rnx.reset_index()

    needed = {"epoch", "prn", "C1"}
    missing = needed - set(df_rnx.columns)
    if missing:
        raise ValueError(f"{label} df_rnx is missing columns: {missing}")

    # Ensure SP3 is in expected "column" form for the interpolation routine
    df_sp3_use = df_sp3.copy()
    if isinstance(df_sp3_use.index, pd.MultiIndex) and set(df_sp3_use.index.names) >= {"epoch", "prn"}:
        df_sp3_use = df_sp3_use.reset_index()

    # Init output columns
    df_rnx["X_sat"] = np.nan
    df_rnx["Y_sat"] = np.nan
    df_rnx["Z_sat"] = np.nan
    df_rnx["dte_sat"] = np.nan
    df_rnx["dRelat"] = np.nan

    # Loop per PRN (same as TP02)
    for prn in df_rnx["prn"].unique():
        df_rnx_prn = df_rnx[df_rnx["prn"] == prn]
        if prn not in df_sp3_use["prn"].unique():
            print(f"{label} - PRN {prn} not in SP3: skipped")
            continue

        df_sp3_prn = df_sp3_use[df_sp3_use["prn"] == prn]

        # Signal flight time τ ≈ C1 / c
        fly_time = pd.to_timedelta(df_rnx_prn["C1"] / conv.SPEED_OF_LIGHT, unit="s")

        t_rec = df_rnx_prn["epoch"]
        t_emi_approx = t_rec - fly_time

        # Approx satellite state at emission time
        orb_df_approx = reffram.orb_df_lagrange_interpolate(df_sp3_prn, t_emi_approx.values)
        orb_df_approx[["x", "y", "z"]] = orb_df_approx[["x", "y", "z"]] * 1e3  # km -> m

        # Relativistic correction (finite-difference velocity)
        delta_t = pd.to_timedelta(1e-3, unit="s")
        orb_fwd = reffram.orb_df_lagrange_interpolate(df_sp3_prn, t_emi_approx.values + delta_t)
        orb_bak = reffram.orb_df_lagrange_interpolate(df_sp3_prn, t_emi_approx.values - delta_t)
        orb_fwd[["x", "y", "z"]] = orb_fwd[["x", "y", "z"]] * 1e3
        orb_bak[["x", "y", "z"]] = orb_bak[["x", "y", "z"]] * 1e3

        v_xyz = (orb_fwd[["x", "y", "z"]] - orb_bak[["x", "y", "z"]]) / (2 * delta_t.total_seconds())
        r_xyz = orb_df_approx[["x", "y", "z"]]
        dRelat_v = -2.0 * (v_xyz * r_xyz).sum(axis=1) / (conv.SPEED_OF_LIGHT ** 2)

        # Refined emission time (apply satellite clock, in microseconds)
        t_emi_ok = t_emi_approx - pd.to_timedelta(orb_df_approx["clk"].values, unit="us")

        orb_ok = reffram.orb_df_lagrange_interpolate(df_sp3_prn, t_emi_ok.values)
        orb_ok[["x", "y", "z"]] = orb_ok[["x", "y", "z"]] * 1e3  # km -> m

        # Write back
        df_rnx.loc[df_rnx_prn.index, "X_sat"] = orb_ok["x"].values
        df_rnx.loc[df_rnx_prn.index, "Y_sat"] = orb_ok["y"].values
        df_rnx.loc[df_rnx_prn.index, "Z_sat"] = orb_ok["z"].values
        df_rnx.loc[df_rnx_prn.index, "dte_sat"] = orb_ok["clk"].values * 1e-6  # us -> s
        df_rnx.loc[df_rnx_prn.index, "dRelat"] = dRelat_v.values

    # Back to canonical GNSS index (epoch, prn)
    df_rnx = df_rnx.set_index(["epoch", "prn"]).sort_index()
    return df_rnx

print("***** Satellite states at emission time (BASE vs ROVER) *****")

# Compute sat state for each station (same routine, different C1 -> different τ -> different t_emit)
df_base_sat  = add_satellite_state_from_sp3(df_base_sync.reset_index(),  df_sp3, label="BASE")
df_rover_sat = add_satellite_state_from_sp3(df_rover_sync.reset_index(), df_sp3, label="ROVER")


# %%
###############################################################################
# Compare satellite positions (ROVER - BASE)
#
# Educational objective
# ---------------------
# Quantify how much the satellite ECEF position differs when computed using
# the receiver-specific emission time (via C1).
# Expected magnitude: centimetric to decimetric for ~10 km baseline.
###############################################################################

df_diff = pd.DataFrame(index=df_base_sat.index)
df_diff["dX"] = df_rover_sat["X_sat"] - df_base_sat["X_sat"]
df_diff["dY"] = df_rover_sat["Y_sat"] - df_base_sat["Y_sat"]
df_diff["dZ"] = df_rover_sat["Z_sat"] - df_base_sat["Z_sat"]
df_diff["dr_norm"] = np.sqrt(df_diff["dX"]**2 + df_diff["dY"]**2 + df_diff["dZ"]**2)

print(df_diff["dr_norm"].describe())
print("\nMax ||Δr|| (m):", df_diff["dr_norm"].max())

print("\nPer-PRN max ||Δr|| (m):")
print(df_diff["dr_norm"].groupby(level="prn").max().sort_values(ascending=False).head(10))


# %%
###############################################################################
# Receiver-dependent satellite geometry
#
# Educational objective
# ---------------------
# Demonstrate that satellite positions depend on the receiver location.
#
# Physical explanation
# --------------------
# Even if BASE and ROVER observe the same satellite at the same reception
# epoch, the signal travel time differs slightly because the geometric
# distance satellite–receiver is not identical.
#
# Consequently:
#
#       t_emit(BASE) ≠ t_emit(ROVER)
#
# and satellites are interpolated at slightly different emission epochs.
#
# This leads to centimetric differences in satellite coordinates between
# receivers, typically at the 5–15 cm level for short GNSS baselines.
#
# Importance for GNSS processing
# ------------------------------
# This effect becomes significant when working with carrier-phase
# observations and motivates the differential GNSS formulation used
# in baseline processing.
###############################################################################

print("\n------------------------------------------------------------")
print("Satellite position depends on receiver location.")
print("BASE and ROVER do NOT see satellites at the same emission time.")
print("Typical satellite position difference:", 
      f"{df_diff['dr_norm'].mean():.3f} m")
print("This effect becomes critical for carrier-phase positioning.")
print("------------------------------------------------------------\n")






# %%
###############################################################################
# Sagnac correction (vectorized Z-rotation) applied to BASE and ROVER
#
# Educational objective
# ---------------------
# Correct satellite ECEF coordinates for Earth rotation during signal travel time.
#
# Physical idea
# -------------
# During the signal flight time τ ≈ C1 / c, the Earth rotates by:
#     dθ = ω_E * τ
#
# To express the satellite coordinates in the Earth-fixed frame at reception time,
# we "undo" the Earth rotation during propagation. This is implemented as a
# backward rotation around the Z-axis:
#     α = -dθ
#
# Practical note (baseline)
# -------------------------
# τ depends on the receiver-satellite geometry, so BASE and ROVER generally have
# slightly different τ and therefore slightly different Sagnac rotations.
###############################################################################

print("***** Sagnac correction (BASE and ROVER) *****")

def add_sagnac_rotated_coords(df_sat, C1_series, suffix="_sagnac"):
    """
    Add Sagnac-rotated satellite coordinates to a satellite-state dataframe.

    Parameters
    ----------
    df_sat : DataFrame (indexed by ['epoch','prn'])
        Must contain X_sat, Y_sat, Z_sat in meters (ECEF).
    C1_series : Series aligned with df_sat index
        Pseudorange in meters, used to approximate flight time τ = C1/c.
    suffix : str
        Suffix used for new columns.

    Returns
    -------
    df_out : DataFrame
        Copy of df_sat with added columns:
        X_sat{suffix}, Y_sat{suffix}, Z_sat{suffix}
    """
    df_out = df_sat.copy()

    # Signal flight time (seconds): τ ≈ C1 / c
    tau_s = C1_series.to_numpy(dtype=float) / conv.SPEED_OF_LIGHT

    # Physical Earth rotation angle during signal travel time
    dtheta = tau_s * conv.EARTH_ROTATION_MEAN_ANGULAR_VELOCITY

    # Backward rotation to express satellite coordinates in reception frame
    alpha = -dtheta

    cos_a = np.cos(alpha)
    sin_a = np.sin(alpha)

    x = df_out["X_sat"].to_numpy(dtype=float)
    y = df_out["Y_sat"].to_numpy(dtype=float)
    z = df_out["Z_sat"].to_numpy(dtype=float)

    x_rot = cos_a * x - sin_a * y
    y_rot = sin_a * x + cos_a * y
    z_rot = z

    df_out[f"X_sat{suffix}"] = x_rot
    df_out[f"Y_sat{suffix}"] = y_rot
    df_out[f"Z_sat{suffix}"] = z_rot

    return df_out


# Apply to BASE and ROVER (C1 differs slightly -> rotation differs slightly)
df_base_sat  = add_sagnac_rotated_coords(df_base_sat,  df_base_sync["C1"])
df_rover_sat = add_sagnac_rotated_coords(df_rover_sat, df_rover_sync["C1"])

# -------------------------------------------------------------------------
# Diagnostics (pedagogical)
# -------------------------------------------------------------------------
# 1) Magnitude of Sagnac correction (per station): ||r_sagnac - r||
def sagnac_magnitude(df_sat):
    dx = df_sat["X_sat_sagnac"] - df_sat["X_sat"]
    dy = df_sat["Y_sat_sagnac"] - df_sat["Y_sat"]
    dz = df_sat["Z_sat_sagnac"] - df_sat["Z_sat"]
    return np.sqrt(dx*dx + dy*dy + dz*dz)

mag_base  = sagnac_magnitude(df_base_sat)
mag_rover = sagnac_magnitude(df_rover_sat)

print("\nSagnac correction magnitude (meters):")
print("BASE :", pd.Series(mag_base, index=df_base_sat.index).describe())
print("ROVER:", pd.Series(mag_rover, index=df_rover_sat.index).describe())

# 2) Difference BASE vs ROVER on Sagnac-corrected satellite positions
df_sagnac_diff = pd.DataFrame(index=df_base_sat.index)
df_sagnac_diff["dX_sagnac"] = df_rover_sat["X_sat_sagnac"] - df_base_sat["X_sat_sagnac"]
df_sagnac_diff["dY_sagnac"] = df_rover_sat["Y_sat_sagnac"] - df_base_sat["Y_sat_sagnac"]
df_sagnac_diff["dZ_sagnac"] = df_rover_sat["Z_sat_sagnac"] - df_base_sat["Z_sat_sagnac"]
df_sagnac_diff["dr_sagnac_norm"] = np.sqrt(
    df_sagnac_diff["dX_sagnac"]**2 +
    df_sagnac_diff["dY_sagnac"]**2 +
    df_sagnac_diff["dZ_sagnac"]**2
)

print("\nBASE vs ROVER satellite position difference AFTER Sagnac (meters):")
print(df_sagnac_diff["dr_sagnac_norm"].describe())



