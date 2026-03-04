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
- matplotlib

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
import matplotlib.pyplot as plt

# GeodeZYX
from geodezyx import files_rw     # read/write module
from geodezyx import conv         # conversion module
from geodezyx import operational  # download rinex module
from geodezyx import gnss_edu     # learning module
from geodezyx import reffram      # reference frame / higher geodesy module

# %%

###############################################################################
# ⚠ Development note
#
# If the module `gnss_edu` was modified during the session, Python keeps the
# old version in memory. The following lines force Python to reload it.
#
# Tip: in Spyder, press Ctrl+G on `gnss_edu` to open the module and explore
# the available functions.
###############################################################################

from importlib import reload
from geodezyx import gnss_edu

gnss_edu = reload(gnss_edu)

# Minimal check: confirm the function exists
print("plot_gnss_sd_by_prn available:", hasattr(gnss_edu, "plot_gnss_sd_by_prn"))

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

def clean_rinex_kiss(df: pd.DataFrame) -> pd.DataFrame:
    """
    Minimal RINEX cleaning (KISS) for baseline TP.

    Key points
    ----------
    - Keep observables as in RINEX (phases stay in cycles)
    - Do not pre-select observables
    - Robustly remove columns that are effectively empty
      (all missing values, including blank strings)

    Steps
    -----
    1) Normalize "missing" string patterns to real NaN
    2) Convert non-identifier columns to numeric when possible
    3) Keep GPS only
    4) Drop rows without (epoch, prn)
    5) Drop columns that are all NaN (do it AFTER filtering too)
    6) Reset index and add a stable row id
    """
    df = df.copy()

    # ------------------------------------------------------------------
    # 1) Normalize common "missing" encodings to NaN
    # ------------------------------------------------------------------
    # If some columns are object dtype, blank strings can survive and prevent
    # "all-NaN" detection. This makes missingness explicit.
    df.replace(
        to_replace=[r"^\s*$", "nan", "NaN", "NA", "N/A", "null", "None"],
        value=np.nan,
        regex=True,
        inplace=True,
    )

    # ------------------------------------------------------------------
    # 2) Convert observable columns to numeric when possible
    # ------------------------------------------------------------------
    protected_cols = {"epoch", "sys", "prn"}
    for col in df.columns:
        if col in protected_cols:
            continue
        # Only attempt conversion for object-like columns (cheap + safe)
        if df[col].dtype == object:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # ------------------------------------------------------------------
    # 3) Keep GPS only
    # ------------------------------------------------------------------
    if "sys" in df.columns:
        df = df[df["sys"] == "G"]
    elif "prn" in df.columns:
        df = df[df["prn"].astype(str).str.startswith("G")]

    # ------------------------------------------------------------------
    # 4) Drop rows without essential identifiers
    # ------------------------------------------------------------------
    df = df.dropna(subset=["epoch", "prn"])

    # ------------------------------------------------------------------
    # 5) Drop fully empty columns (do it after filtering!)
    # ------------------------------------------------------------------
    df = df.loc[:, df.notna().any(axis=0)]

    # ------------------------------------------------------------------
    # 6) Reset indexing and add stable row id
    # ------------------------------------------------------------------
    df = df.reset_index(drop=True)
    df["ind_ligne"] = np.arange(len(df))

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
df_base_sat  = add_sagnac_rotated_coords(df_base_sat,  df_base_sat["C1"])
df_rover_sat = add_sagnac_rotated_coords(df_rover_sat, df_rover_sat["C1"])

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


# %%
###############################################################################
# Synchronization check between BASE and ROVER
###############################################################################

same_index = df_base_sat.index.equals(df_rover_sat.index)

print("BASE / ROVER synchronized:", same_index)

if not same_index:
    raise RuntimeError(
        "BASE and ROVER datasets are not perfectly synchronized."
    )
    
print(df_base_sat.index[:5])
print(df_rover_sat.index[:5])

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
    "df_base","df_base_sync",
    "df_rover","df_rover_sync",
    "df_diff", "df_sp3",
    "common_index", "df_sagnac_diff", "mag_base", "mag_rover", "same_index"
]:
    if var in globals():
        del globals()[var]

del var
gc.collect()

# %%
###############################################################################
# Sanity checks before differential processing (SD / DD)
#
# We now work ONLY with:
#   - df_base_sat
#   - df_rover_sat
#
# These dataframes must be:
#   1) perfectly synchronized (same MultiIndex)
#   2) containing satellite coordinates at emission time
#   3) containing Sagnac-corrected coordinates for later geometry (optional)
###############################################################################

required_cols = [
    "C1",
    "X_sat", "Y_sat", "Z_sat",
    "X_sat_sagnac", "Y_sat_sagnac", "Z_sat_sagnac",
]

if not df_base_sat.index.equals(df_rover_sat.index):
    raise RuntimeError("df_base_sat and df_rover_sat are not synchronized (index mismatch).")

missing_base = [c for c in required_cols if c not in df_base_sat.columns]
missing_rover = [c for c in required_cols if c not in df_rover_sat.columns]

if missing_base:
    raise RuntimeError(f"Missing columns in df_base_sat: {missing_base}")
if missing_rover:
    raise RuntimeError(f"Missing columns in df_rover_sat: {missing_rover}")

print("OK: df_base_sat and df_rover_sat are synchronized and contain Sagnac-corrected satellite coordinates.")
print("Index example:", df_base_sat.index[:3].tolist())

for var in [    "required_cols","missing_base", "missing_rover","df_rover_sync"]:
    if var in globals():
        del globals()[var]
del var
gc.collect()


# %%
###############################################################################
# Single Differences (SD): ROVER − BASE on meaningful observables (C*, P*, L*)
#
# Educational objective
# ---------------------
# Build single differences for observables that carry geometric information:
#   - Code         : C* and P*  (meters in RINEX)
#   - Carrier phase: L*         (cycles, KISS choice)
#
# We intentionally exclude:
#   - Doppler (D*) and SNR (S*) : diagnostics only, not used for baseline SD/DD
#   - Flags (LLI/SSI)           : indicators, not measurements
#
# Inputs
# ------
# df_base_sat and df_rover_sat must be synchronized:
#   - same MultiIndex (epoch, prn)
#
# Output
# ------
# df_SD: DataFrame indexed by (epoch, prn) with columns prefixed by "SD_"
###############################################################################

print("***** Single Differences (ROVER − BASE) on C*, P*, L* *****")


def compute_single_differences(
    df_base_sat: pd.DataFrame,
    df_rover_sat: pd.DataFrame,
    prefixes=("C", "P", "L"),
    add_prefix="SD_",
) -> pd.DataFrame:
    """
    Compute GNSS single differences (ROVER − BASE) for meaningful observables.

    Meaningful observables (default)
    -------------------------------
    - Code  : C*, P*
    - Phase : L*        (kept in cycles here, KISS choice)

    Excluded
    --------
    - Doppler (D*), SNR (S*)
    - Indicators (*_LLI, *_SSI)

    Inputs
    ------
    df_base_sat and df_rover_sat must be synchronized and indexed by (epoch, prn).

    Returns
    -------
    df_SD : DataFrame indexed by (epoch, prn) with columns SD_*
    """

    # -------------------------------------------------------------------------
    # Safety check: perfect synchronization
    # -------------------------------------------------------------------------
    if not df_base_sat.index.equals(df_rover_sat.index):
        raise RuntimeError("df_base_sat and df_rover_sat are not synchronized (index mismatch).")

    # -------------------------------------------------------------------------
    # Select common columns (present in both dataframes)
    # -------------------------------------------------------------------------
    common_cols = [c for c in df_base_sat.columns if c in df_rover_sat.columns]

    def is_meaningful_obs(col: str) -> bool:
        # Exclude non-measurement indicators
        if col.endswith(("_LLI", "_SSI")):
            return False
        # Keep only code and phase observables
        return col.startswith(prefixes)

    cand_cols = [c for c in common_cols if is_meaningful_obs(c)]

    # Keep only columns producing at least one valid SD value
    usable_cols = []
    for c in cand_cols:
        sd_tmp = df_rover_sat[c] - df_base_sat[c]
        if sd_tmp.notna().any():
            usable_cols.append(c)

    # -------------------------------------------------------------------------
    # Compute SD (vectorized)
    # -------------------------------------------------------------------------
    df_SD = df_rover_sat[usable_cols].subtract(df_base_sat[usable_cols]).add_prefix(add_prefix)

    return df_SD


# -------------------------------------------------------------------------
# Compute SD using the dedicated function
# -------------------------------------------------------------------------
df_SD = compute_single_differences(df_base_sat, df_rover_sat)

# -------------------------------------------------------------------------
# Quick diagnostics for students
# -------------------------------------------------------------------------
print(f"SD computed for {df_SD.shape[1]} observables.")
print("SD columns:", list(df_SD.columns))

print("\nValid SD counts (top 15):")
print(df_SD.notna().sum().sort_values(ascending=False).head(15))

# Optional: split by type (useful for next steps)
sd_phase_cols = [c for c in df_SD.columns if c.startswith("SD_L")]
sd_code_cols  = [c for c in df_SD.columns if c.startswith(("SD_C", "SD_P"))]

print(f"\nPhase SD columns: {len(sd_phase_cols)}")
print(f"Code  SD columns: {len(sd_code_cols)}")

print("\nExample SD values:")
print(df_SD.head())

#%%
###############################################################################
# Visual diagnostics on simple differences (SD)
#
# Educational objective
# ---------------------
# 1) Visualize SD on carrier phase (SD_L1) and code (SD_C1) by satellite (PRN).
# 2) Highlight discontinuities (arcs) using "gap handling".
# 3) Use time-derivative to reveal sharp events (cycle slips / bad code at low elevation).
# 4) Automatically detect "holes" (missing epochs inside arcs), which often imply
#    that a NEW ambiguity parameter must be introduced for phase processing.
#
# Reminder
# --------
# - SD_L1 is in cycles  (KISS choice: no wavelength conversion here)
# - SD_C1 is in meters
###############################################################################

print("\n===== SD diagnostics (by satellite) =====")

# -------------------------------------------------------------------------
# 1) Carrier phase SD_L1 (cycles): should show smooth arcs + possible jumps
# -------------------------------------------------------------------------
print("\n[1/4] Plot SD_L1 (phase, cycles): arcs + potential phase jumps")
gnss_edu.plot_gnss_sd_by_prn(
    df_SD,
    observable="SD_L1",
    label_arcs=True,
    show_legend=False
)

# -------------------------------------------------------------------------
# 2) Code SD_C1 (meters): typically noisier and more sensitive at rise/set
# -------------------------------------------------------------------------
print("\n[2/4] Plot SD_C1 (code, meters): usually noisier (multipath, low elevation)")
gnss_edu.plot_gnss_sd_by_prn(
    df_SD,
    observable="SD_C1",
    label_arcs=True,
    show_legend=False
)

# -------------------------------------------------------------------------
# 3) Time derivative: reveals sharp events (cycle slips / bad code segments)
# -------------------------------------------------------------------------
print("\n[3/4] Plot d/dt(SD_L1): sharp peaks often indicate cycle slips or resets")

gnss_edu.plot_sd_derivative_by_prn(
    df_SD,
    obs="SD_C1",
    normalize_by_dt=False,              # keep raw ΔSD between epochs
    gap=pd.Timedelta(seconds=30)        # expected sampling (30 s)
    
)

gnss_edu.plot_sd_derivative_by_prn(
    df_SD,
    obs="SD_L1",
    normalize_by_dt=False,              # keep raw ΔSD between epochs
    gap=pd.Timedelta(seconds=30)        # expected sampling (30 s)
)

gnss_edu.plot_sd_derivative_by_prn(
    df_SD,
    obs="SD_L2",
    normalize_by_dt=False,              # keep raw ΔSD between epochs
    gap=pd.Timedelta(seconds=30)        # expected sampling (30 s)
)

# -------------------------------------------------------------------------
# 4) Detect holes in raw RINEX observation timeline (per PRN)
#    Holes inside an arc are important: they often force a new ambiguity.
# -------------------------------------------------------------------------
print("\n[4/4] Detect intra-arc holes (missing epochs) in ROVER observations")
holes = gnss_edu.detect_intra_arc_holes(
    df_rover_sat,
    sampling=pd.Timedelta(seconds=30),  # expected sampling
    gap=pd.Timedelta(minutes=30)        # arc-segmentation threshold (same as plots)
)

print(f"Holes found: {len(holes)}")
display(holes.sort_values("dt", ascending=False).head(20))



