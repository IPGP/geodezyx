#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 03 - Differential GNSS processing on a short baseline with carrier phase

Educational scope
-----------------
This third step introduces differential GNSS processing on a short baseline.
The pedagogical objective is to move from synchronized observations at two
receivers to single differences and double differences, with special attention
to pivot selection, arc continuity, and carrier-phase behavior.

Learning goals
--------------
By the end of this script, students should be able to:
1. synchronize BASE and ROVER observations on common (epoch, satellite) pairs,
2. compute receiver-dependent satellite states at emission time,
3. quantify the effect of satellite geometry, Earth rotation, and data gaps,
4. build single differences on code and phase observables,
5. define a stable pivot strategy for double-difference processing,
6. construct and interpret a rich double-difference table.

Pedagogical position in the curriculum
--------------------------------------
Step 02 focused on code-based point positioning with progressive model
enrichment for a single receiver.

Step 03 shifts to relative GNSS processing. The main pedagogical idea is that
short-baseline carrier-phase positioning is built from synchronized
observations, receiver-dependent satellite geometry, and a carefully managed
pivot strategy.

Important note
--------------
This script deliberately keeps some modeling choices simple in order to expose
the structure of differential GNSS processing. The emphasis is on
understanding the workflow and the observables rather than delivering a final
survey-grade baseline solution.

Author: Samuel Nahmani (1,2)
Web: https://www.ipgp.fr/annuaire/nahmani/
Contact: nahmani@ipgp.fr | samuel.nahmani@ign.fr

(1) Université Paris Cité, Institut de physique du globe de Paris (IPGP),
    CNRS, IGN, F-75005 Paris, France
(2) Univ Gustave Eiffel, Géodata Paris, IGN, F-75238 Paris, France

Version: 1.0

Dependencies (Python packages)
------------------------------
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
###############################################################################
# Reference
###############################################################################

# Sakic, P., Mansur, G., Chaiyaporn, K., & Ballu, V. (2019).
# The geodeZYX toolbox: a versatile Python 3 toolbox for geodetic-oriented
# purposes (v4.0). GFZ Data Services.
# https://doi.org/10.5880/GFZ.1.1.2019.002

# %%
###############################################################################
# Imports
###############################################################################

import datetime as dt
import gc
import os
from pathlib import Path

from importlib import reload

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from IPython.display import display

from geodezyx import conv
from geodezyx import files_rw
from geodezyx import gnss_edu
from geodezyx import operational
from geodezyx import reffram

# %%
###############################################################################
# Configuration
###############################################################################

PROCESSING_DATE = dt.datetime(2019, 6, 25)
WORK_DIR = Path(os.environ["HOME"]).expanduser() / "gnss_edu_data"

BASE_STATION = "SMNE"
ROVER_STATION = "MLVL"
STATION_DICT = {"rgp": [BASE_STATION, ROVER_STATION]}

SATELLITE_SYSTEM = "G"
PRIMARY_CODE_OBSERVABLE = "C1"
SNR_COLUMN = "S1"
PIVOT_ELEVATION_COLUMN = "elevation_deg"

RINEX_SAMPLING = pd.Timedelta(seconds=30)
ARC_GAP = pd.Timedelta(minutes=30)
PIVOT_ELEVATION_CUTOFF_DEG = 15.0
PIVOT_SNR_MIN = 48.0
PIVOT_SWITCH_MARGIN_DB = 2.0
# Examples:
#   ("L1",)       -> mono-frequency pivot compatibility
#   ("L1", "L5")  -> keep only pivots compatible with iono-free L1/L5 work
PIVOT_REQUIRED_OBSERVABLES = ("L1",)
PIVOT_MIN_DURATION = pd.Timedelta(hours=1)
PIVOT_MIN_DWELL = pd.Timedelta(minutes=45)

print("Configuration loaded.")
print("Processing date          :", PROCESSING_DATE)
print("Working directory        :", WORK_DIR)
print("Base station             :", BASE_STATION)
print("Rover station            :", ROVER_STATION)
print("Satellite system         :", SATELLITE_SYSTEM)
print("Primary code observable  :", PRIMARY_CODE_OBSERVABLE)
print("SNR observable           :", SNR_COLUMN)
print("Elevation cutoff [deg]   :", PIVOT_ELEVATION_CUTOFF_DEG)
print("Required pivot observables:", PIVOT_REQUIRED_OBSERVABLES)


# %%
###############################################################################
# Pedagogical roadmap
###############################################################################

print(
    "Step 03 roadmap\n"
    "----------------\n"
    "M0: synchronize BASE and ROVER observations\n"
    "M1: compute receiver-dependent satellite states at emission time\n"
    "M2: quantify Sagnac and geometry differences between receivers\n"
    "M3: build single differences on code and carrier phase\n"
    "M4: inspect arcs, gaps, and derivative diagnostics\n"
    "M5: define a stable pivot schedule from SNR\n"
    "M6: build rich double differences with pivot traceability"
)

###############################################################################
# Interactive development note
###############################################################################

gnss_edu = reload(gnss_edu)

print(
    "Interactive helper check: plot_gnss_sd_by_prn available:",
    hasattr(gnss_edu, "plot_gnss_sd_by_prn"),
)

# %%
###############################################################################
# Helper functions
###############################################################################

def extract_download_path(entry):
    """Return the local file path from a geodezyx download entry."""
    if isinstance(entry, tuple):
        return entry[0]
    return entry


def prepare_work_directory(work_dir: Path) -> Path:
    """Create the working directory used for downloads and outputs."""
    work_dir.mkdir(parents=True, exist_ok=True)
    return work_dir


def download_station_and_product_data(
    work_dir: Path,
    processing_date: dt.datetime,
    station_dict: dict[str, list[str]],
):
    """Download RINEX observations and precise products for one processing day."""
    station_downloads = operational.download_gnss_rinex(
        statdico=station_dict,
        output_dir=str(work_dir),
        startdate=processing_date,
        enddate=processing_date,
        parallel_download=1,
    )

    product_downloads = operational.download_gnss_products(
        archive_dir=str(work_dir),
        startdate=processing_date,
        enddate=processing_date,
        archtype="year/doy",
        AC_names=("IGS",),
        repro=0,
        archive_center="ign",
        parallel_download=1,
    )

    return station_downloads, product_downloads


def assign_base_and_rover_paths(
    station_downloads,
    base_station: str,
    rover_station: str,
) -> tuple[str, str]:
    """Map downloaded RINEX files to the requested BASE and ROVER stations."""
    station_paths = {}

    for path, ok in station_downloads:
        if not ok:
            raise RuntimeError(f"Download failed for station entry: {path}")
        station_paths[Path(path).name[:4].upper()] = path

    try:
        return station_paths[base_station], station_paths[rover_station]
    except KeyError as exc:
        available = sorted(station_paths)
        raise RuntimeError(
            f"Missing requested station in downloaded files: {exc}. "
            f"Available stations: {available}"
        ) from exc


def read_sp3_product(product_downloads: list) -> tuple[str, pd.DataFrame]:
    """Read the first downloaded SP3 product as a pandas DataFrame."""
    sp3_path = extract_download_path(product_downloads[0])
    df_sp3 = files_rw.read_sp3(sp3_path, returns_pandas=True, new_col_names=True)
    return sp3_path, df_sp3


def extract_approx_position_from_rinex_header(rinex_header: list[str]) -> np.ndarray:
    """Extract APPROX POSITION XYZ from RINEX header lines."""
    approx_position_lines = [
        line for line in rinex_header if "APPROX POSITION XYZ" in line
    ]
    if len(approx_position_lines) == 0:
        raise RuntimeError("Approximate receiver position not found in RINEX header.")
    return np.array(approx_position_lines[0].split()[:3], dtype=float)


def cleanup_variables(*variable_names: str) -> None:
    """Delete intermediate globals when running the script interactively."""
    for variable_name in variable_names:
        if variable_name in globals():
            del globals()[variable_name]
    gc.collect()

# %%
###############################################################################
# Working directory
###############################################################################

WORK_DIR = prepare_work_directory(WORK_DIR)
print("Working directory ready:", WORK_DIR)


# %%
###############################################################################
# Download RINEX observations and precise products
###############################################################################

station_downloads, product_downloads = download_station_and_product_data(
    work_dir=WORK_DIR,
    processing_date=PROCESSING_DATE,
    station_dict=STATION_DICT,
)


# %%
###############################################################################
# Precise orbit product
###############################################################################

sp3_path, df_sp3 = read_sp3_product(product_downloads)
print(f"Using SP3 file: {sp3_path}")

# %%
###############################################################################
# RINEX observation files
###############################################################################

rinex_base_path, rinex_rover_path = assign_base_and_rover_paths(
    station_downloads,
    base_station=BASE_STATION,
    rover_station=ROVER_STATION,
)

print(f"BASE  ({BASE_STATION}):", rinex_base_path)
print(f"ROVER ({ROVER_STATION}):", rinex_rover_path)

# Read RINEX observation files
df_base, rinex_base_header = files_rw.read_rinex_obs(
    rinex_base_path,
    return_header=True,
)
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
# GNSS canonical indexing (same structure as step02.py)
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

cleanup_variables("station_downloads", "product_downloads")


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



# %%
###############################################################################
# Data gaps in SD time series (per PRN) — why we care before SD plots
#
# Key idea
# --------
# In baseline processing, the SD dataset (ROVER − BASE) exists only when
# BOTH stations tracked the SAME satellite at the SAME epoch.
#
# Therefore, every "white" zone in the timeline below means:
#   - no SD can be formed (at least one station missed the observation)
#
# GNSS gaps have two different meanings:
#   (A) Long gaps   -> satellite not visible (rise/set, masking) => new arc
#   (B) Short gaps  -> missing epochs inside an arc             => tracking loss
#
# Practical consequence
# ---------------------
# Carrier-phase ambiguities are at least "one per arc".
# Short gaps inside an arc often force a NEW ambiguity parameter as well.
###############################################################################

print("\n===== SD gap analysis (base + rover common tracking) =====")

# Expected RINEX sampling (here: 30 s)
sampling = RINEX_SAMPLING

# -------------------------------------------------------------------------
# Step 0 — Visual inspection (fast intuition)
# -------------------------------------------------------------------------
print("\n[0] SD tracking timeline (black = present SD, white = missing)")
print("Look for: (i) long white blocks (visibility breaks), (ii) small holes inside black arcs.")
gnss_edu.plot_tracking_timeline(df_SD, sampling=sampling)

# -------------------------------------------------------------------------
# Step 1 — Detect ALL gaps (no arc segmentation)
#   We intentionally use a huge 'gap' so each PRN is treated as one long arc.
#   This reveals the full distribution of gap durations.
# -------------------------------------------------------------------------
gap_all = pd.Timedelta(hours=30)   # huge on purpose

holes_all = gnss_edu.detect_intra_arc_holes(
    df_SD,
    sampling=sampling,
    gap=gap_all
)

print("\n[1/3] ALL gaps detected (no arc segmentation)")
print(f"holes_all: {len(holes_all)} gaps found (gap={gap_all}, sampling={sampling})")
print("Largest gaps (often rise/set or masking):")
print(holes_all.sort_values("dt", ascending=False).head(20))

# -------------------------------------------------------------------------
# Step 2 — Gap-duration distribution
#   Goal: separate the two regimes visually:
#     - hours-scale gaps   -> visibility breaks (new arcs)
#     - minutes/seconds    -> true missing epochs inside arcs
# -------------------------------------------------------------------------
print("\n[2/3] Gap-duration distribution (minutes)")

dt_minutes = holes_all["dt"].dt.total_seconds() / 60.0

fig, ax = plt.subplots(figsize=(10, 4))
ax.hist(dt_minutes, bins=80)
ax.set_title("Distribution of gap durations (minutes)")
ax.set_xlabel("dt (minutes)")
ax.set_ylabel("count")
plt.tight_layout()
plt.show()

print("Interpretation:")
print("  - hours-scale gaps   -> visibility breaks (rise/set, masking) => new arc")
print("  - minutes/seconds    -> true missing epochs inside arcs       => tracking loss")

# -------------------------------------------------------------------------
# Step 3 — Detect holes INSIDE arcs (pedagogical arc segmentation)
#   Now we set a realistic arc boundary (e.g., 30 min).
#   Any missing epochs smaller than that occur *inside* a tracking arc.
# -------------------------------------------------------------------------
gap_arc = ARC_GAP

holes_arc = gnss_edu.detect_intra_arc_holes(
    df_SD,
    sampling=sampling,
    gap=gap_arc
)

print("\n[3/3] Holes INSIDE arcs (after arc segmentation)")
print(f"holes_arc: {len(holes_arc)} gaps found (gap={gap_arc}, sampling={sampling})")
print("Largest intra-arc holes (these often trigger new ambiguity parameters):")
print(holes_arc.sort_values("dt", ascending=False).head(20))

# -------------------------------------------------------------------------
# Optional — Arc count proxy per PRN (very useful discussion point)
#   Each visibility break (dt >= gap_arc) splits tracking into a new arc.
#   So: n_arcs ≈ 1 + number_of_visibility_breaks
# -------------------------------------------------------------------------
visibility_breaks = holes_all[holes_all["dt"] >= gap_arc]
arcs_per_prn = (visibility_breaks.groupby(level="prn").size() + 1).rename("n_arcs_proxy")

print("\nOptional: arc count proxy per PRN (n_arcs ≈ 1 + #visibility breaks)")
print(arcs_per_prn.sort_index().to_frame().head(30))




#%%
###############################################################################
# Visual diagnostics on Single Differences (SD)
#
# Pedagogical goal
# ---------------
# Now that we have quantified data gaps, we can choose plotting parameters
# from the data instead of arbitrary values.
#
# What we want to show
# --------------------
# (1) SD_L1 (phase, cycles): smooth arcs, possible discontinuities (cycle slips)
# (2) SD_C1 (code, meters): more affected by noise/multipath, especially at rise/set
# (3) Epoch-to-epoch increments ΔSD: spikes reveal abnormal events
# (4) Time-normalized derivative d(SD)/dt: same idea, but expressed per second
#
# Notes
# -----
# - SD_L* is in cycles (KISS choice here)
# - For increments/derivatives, we convert phase cycles -> meters to make the
#   noise level comparable to code.
###############################################################################

print("\n===== SD diagnostics (by satellite PRN) =====")

# -------------------------------------------------------------------------
# Parameters justified by the previous "hole detection" cell
# -------------------------------------------------------------------------
sampling = RINEX_SAMPLING                # nominal RINEX sampling for this dataset
gap_arc = ARC_GAP                        # arc segmentation threshold (visibility breaks)

# For increments/derivatives, we want a much tighter gap threshold:
# accept up to 3 missed epochs (<= 90 s) before splitting the segment.
gap_inc  = 3 * sampling + pd.Timedelta(seconds=1)   # 91 s safe bound

print("\nPlot parameters derived from the gap analysis:")
print(f"  sampling = {sampling}")
print(f"  gap_arc  = {gap_arc}   (used to split visibility arcs)")
print(f"  gap_inc  = {gap_inc}   (used for increments/derivatives; do not bridge missing epochs)")

# %%
# -------------------------------------------------------------------------
# 1) Carrier phase SD_L1 (cycles): arcs + potential phase jumps
# -------------------------------------------------------------------------
print("\n[1/4] SD_L1 (phase, cycles): arcs + potential discontinuities")
gnss_edu.plot_gnss_sd_by_prn(
    df_SD,
    observable="SD_L1",
    gap=gap_arc,
    label_arcs=True,
    show_legend=False
)

# %%
# -------------------------------------------------------------------------
# 2) Code SD_C1 (meters): noisier, especially near rise/set
# -------------------------------------------------------------------------
print("\n[2/4] SD_C1 (code, meters): typically noisier near rise/set (multipath)")
gnss_edu.plot_gnss_sd_by_prn(
    df_SD,
    observable="SD_C1",
    gap=gap_arc,
    label_arcs=True,
    show_legend=False
)

# %%
# -------------------------------------------------------------------------
# 3) Epoch-to-epoch increments ΔSD: highlight sharp events (raw jumps)
#    - normalize_by_dt=False -> raw jump between consecutive epochs
#    - gap_inc is tight      -> do NOT connect points across missing epochs
# -------------------------------------------------------------------------
print("\n[3/4] Epoch-to-epoch increments ΔSD (raw jumps)")
print("      (phase increments are converted to meters for comparability)")

gnss_edu.plot_sd_derivative_by_prn(
    df_SD,
    obs="SD_C1",
    normalize_by_dt=False,
    gap=gap_inc,
    title="ΔSD_C1 between consecutive epochs (meters)"
)

gnss_edu.plot_sd_derivative_by_prn(
    df_SD,
    obs="SD_L1",
    phase_to_meters=True,
    normalize_by_dt=False,
    gap=gap_inc,
    title="ΔSD_L1 between consecutive epochs (meters, phase converted)"
)

if "SD_L2" in df_SD.columns:
    gnss_edu.plot_sd_derivative_by_prn(
        df_SD,
        obs="SD_L2",
        phase_to_meters=True,
        normalize_by_dt=False,
        gap=gap_inc,
        title="ΔSD_L2 between consecutive epochs (meters, phase converted)"
    )
else:
    print("Note: SD_L2 not available in df_SD (skipped).")

# %%
# -------------------------------------------------------------------------
# 4) Time-normalized derivative d(SD)/dt: units per second
#    Useful when dt is not perfectly constant, or to compare datasets.
# -------------------------------------------------------------------------
print("\n[4/4] Time-normalized derivative d(SD)/dt (units per second)")
print("      (more 'physical', but spikes may look smaller if dt is larger)")

gnss_edu.plot_sd_derivative_by_prn(
    df_SD,
    obs="SD_C1",
    normalize_by_dt=True,
    gap=gap_inc,
    title="d(SD_C1)/dt (m/s)"
)

gnss_edu.plot_sd_derivative_by_prn(
    df_SD,
    obs="SD_L1",
    phase_to_meters=True,
    normalize_by_dt=True,
    gap=gap_inc,
    title="d(SD_L1)/dt (m/s, phase converted)"
)

if "SD_L2" in df_SD.columns:
    gnss_edu.plot_sd_derivative_by_prn(
        df_SD,
        obs="SD_L2",
        phase_to_meters=True,
        normalize_by_dt=True,
        gap=gap_inc,
        title="d(SD_L2)/dt (m/s, phase converted)"
    )
    


# %%

#            satellite s
# BASE  -------------------------
# ROVER -------------------------

# Single difference:   SD_s = obs_rover − obs_base       <------ OK

# Choose a pivot p                                       <------ NOW

# Double difference:   DD_s = SD_s − SD_p

###############################################################################
def add_satellite_elevation_columns(
    df_sat: pd.DataFrame,
    receiver_xyz: np.ndarray,
    x_col: str = "X_sat_sagnac",
    y_col: str = "Y_sat_sagnac",
    z_col: str = "Z_sat_sagnac",
    elev_col: str = PIVOT_ELEVATION_COLUMN,
) -> pd.DataFrame:
    """Add azimuth/elevation columns following the Step 02 M5 logic."""
    df_out = df_sat.copy()

    x0, y0, z0 = receiver_xyz
    sat_xyz = df_out[[x_col, y_col, z_col]].to_numpy(dtype=float)

    azimuth_rad_list = []
    elevation_rad_list = []
    slant_range_m_list = []

    for sat_xyz_i in sat_xyz:
        azi_rad_i, ele_rad_i, slant_range_i = conv.xyz2azi_ele(
            sat_xyz_i[0], sat_xyz_i[1], sat_xyz_i[2],
            x0, y0, z0,
            outdeg=False,
        )
        azimuth_rad_list.append(azi_rad_i)
        elevation_rad_list.append(ele_rad_i)
        slant_range_m_list.append(slant_range_i)

    df_out["azimuth_rad"] = np.array(azimuth_rad_list, dtype=float)
    df_out["elevation_rad"] = np.array(elevation_rad_list, dtype=float)
    df_out["azimuth_deg"] = np.degrees(df_out["azimuth_rad"])
    df_out[elev_col] = np.degrees(df_out["elevation_rad"])
    df_out["slant_range_m"] = np.array(slant_range_m_list, dtype=float)
    return df_out


def plot_tracking_timeline_from_mask(
    df: pd.DataFrame,
    mask_col: str,
    sampling: pd.Timedelta,
    title: str,
):
    """Plot a black/white tracking timeline from a boolean availability mask."""
    if mask_col not in df.columns:
        raise ValueError(f"mask_col='{mask_col}' not found in df.")

    epochs_obs = df.index.get_level_values("epoch")
    prns = np.array(sorted(df.index.get_level_values("prn").unique()))

    t0, t1 = epochs_obs.min(), epochs_obs.max()
    expected_epochs = pd.date_range(start=t0, end=t1, freq=sampling)
    epoch_to_col = {t: i for i, t in enumerate(expected_epochs)}

    avail = np.zeros((len(prns), len(expected_epochs)), dtype=np.uint8)

    for i, prn in enumerate(prns):
        d = df.xs(prn, level="prn").sort_index()
        d = d[d[mask_col].fillna(False).astype(bool)]

        cols = [epoch_to_col.get(t) for t in d.index]
        cols = [c for c in cols if c is not None]
        if cols:
            avail[i, cols] = 1

    fig, ax = plt.subplots(figsize=(14, max(6, 0.25 * len(prns))))
    cmap = plt.matplotlib.colors.ListedColormap(["white", "black"])
    ax.imshow(avail, aspect="auto", interpolation="nearest", cmap=cmap)

    ax.set_yticks(np.arange(len(prns)))
    ax.set_yticklabels(prns)
    ax.set_ylabel("PRN")

    nticks = 8
    xt = np.linspace(0, len(expected_epochs) - 1, nticks).astype(int)
    ax.set_xticks(xt)
    ax.set_xticklabels([expected_epochs[j].strftime("%m-%d %H:%M") for j in xt])
    ax.set_xlabel("Time (epoch)")
    ax.set_title(title)

    plt.tight_layout()
    plt.show()
    return fig, ax


def plot_tracking_timeline_colored(
    df: pd.DataFrame,
    value_col: str,
    sampling: pd.Timedelta,
    title: str,
    cmap: str = "viridis",
    vmin: float | None = None,
    vmax: float | None = None,
):
    """Plot satellite presence with color driven by a per-epoch value."""
    if value_col not in df.columns:
        raise ValueError(f"value_col='{value_col}' not found in df.")

    epochs_obs = df.index.get_level_values("epoch")
    prns = np.array(sorted(df.index.get_level_values("prn").unique()))
    prn_to_row = {prn: i for i, prn in enumerate(prns)}

    t0, t1 = epochs_obs.min(), epochs_obs.max()
    expected_epochs = pd.date_range(start=t0, end=t1, freq=sampling)
    epoch_to_col = {t: i for i, t in enumerate(expected_epochs)}

    values = df[value_col].dropna()
    rows = np.array([prn_to_row[prn] for prn in values.index.get_level_values("prn")])
    cols = np.array([epoch_to_col[t] for t in values.index.get_level_values("epoch")])

    fig, ax = plt.subplots(figsize=(14, max(6, 0.25 * len(prns))))
    sc = ax.scatter(
        cols,
        rows,
        c=values.to_numpy(dtype=float),
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        marker="s",
        s=18,
        linewidths=0,
    )

    ax.set_yticks(np.arange(len(prns)))
    ax.set_yticklabels(prns)
    ax.set_ylabel("PRN")

    nticks = 8
    xt = np.linspace(0, len(expected_epochs) - 1, nticks).astype(int)
    ax.set_xticks(xt)
    ax.set_xticklabels([expected_epochs[j].strftime("%m-%d %H:%M") for j in xt])
    ax.set_xlabel("Time (epoch)")
    ax.set_title(title)

    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label(value_col)

    plt.tight_layout()
    plt.show()
    return fig, ax


def build_observable_compatibility_mask(
    df: pd.DataFrame,
    required_observables: tuple[str, ...],
) -> pd.Series:
    """Return True where all required observables are available."""
    if len(required_observables) == 0:
        return pd.Series(True, index=df.index)

    missing_cols = [obs for obs in required_observables if obs not in df.columns]
    if missing_cols:
        raise ValueError(
            f"Required observables are missing from dataframe: {missing_cols}"
        )

    return df.loc[:, list(required_observables)].notna().all(axis=1)


def build_pivot_support_dataframe(
    df: pd.DataFrame,
    elev_col: str,
    elev_cutoff_deg: float,
    snr_col: str,
    snr_min: float,
    required_observables: tuple[str, ...],
    support_col: str = "pivot_support_score",
) -> pd.DataFrame:
    """
    Build cumulative pivot-availability flags from elevation, SNR and observables.

    The support column keeps the original SNR only when all filters are satisfied.
    It can therefore be passed to the existing set-cover and scheduling utilities.
    """
    df_out = df.copy()

    df_out["pivot_elev_ok"] = df_out[elev_col] >= elev_cutoff_deg
    df_out["pivot_snr_ok"] = df_out[snr_col].notna() & (df_out[snr_col] >= snr_min)
    df_out["pivot_obs_ok"] = build_observable_compatibility_mask(
        df_out,
        required_observables=required_observables,
    )
    df_out["pivot_elev_snr_ok"] = df_out["pivot_elev_ok"] & df_out["pivot_snr_ok"]
    df_out["pivot_usable"] = df_out["pivot_elev_snr_ok"] & df_out["pivot_obs_ok"]
    df_out[support_col] = df_out[snr_col].where(df_out["pivot_usable"])

    return df_out


###############################################################################
# Pivot satellites for Double Differences (DD)
# --------------------------------------------
# Progressive pedagogical strategy:
#   0) inspect raw availability,
#   1) apply an elevation mask,
#   2) refine with SNR,
#   3) keep only satellites compatible with required observables,
#   4) build a minimal and stable pivot schedule over the full session.
#
# Pedagogical assumption
# ----------------------
# In this short-baseline educational setup, pivot satellites are selected from
# the BASE receiver perspective. This choice is justified because:
#   - the base station is fixed,
#   - it is assumed to offer the most stable observing conditions,
#   - the rover is assumed to remain close enough to the base that satellite
#     visibility is very similar.
#
# This approximation is reasonable for short baselines, but it should be
# questioned when the rover moves farther away from the base.
###############################################################################

snr_col = SNR_COLUMN
elev_col = PIVOT_ELEVATION_COLUMN
sampling = RINEX_SAMPLING
gap_arc = ARC_GAP
elev_cutoff_deg = PIVOT_ELEVATION_CUTOFF_DEG
snr_min = PIVOT_SNR_MIN
required_observables = PIVOT_REQUIRED_OBSERVABLES
switch_margin_db = PIVOT_SWITCH_MARGIN_DB
support_col = "pivot_support_score"

base_receiver_xyz = extract_approx_position_from_rinex_header(rinex_base_header)

# Keep the same pedagogical logic as Step 02 / M5:
# define explicitly the receiver position used as geometry reference.
receiver_xyz_geom_ref = base_receiver_xyz.copy()
print("Using approximate RINEX-header BASE position as geometry reference.")

df_base_sat = add_satellite_elevation_columns(df_base_sat, receiver_xyz_geom_ref)
print("Elevation statistics before cutoff [deg]:")
print(df_base_sat[elev_col].describe())

df_base_pivot = build_pivot_support_dataframe(
    df=df_base_sat,
    elev_col=elev_col,
    elev_cutoff_deg=elev_cutoff_deg,
    snr_col=snr_col,
    snr_min=snr_min,
    required_observables=required_observables,
    support_col=support_col,
)

print("\n=== Pivot strategy parameters ===")
print(f"Elevation column       : {elev_col}")
print(f"Elevation cutoff [deg] : {elev_cutoff_deg}")
print(f"SNR observable         : {snr_col}")
print(f"SNR threshold          : {snr_min}")
print(f"Required observables   : {required_observables}")
print(f"Sampling               : {sampling}")
print(f"Hysteresis             : {switch_margin_db} SNR units")

# --------------------------------------------------
# Level 0 — Raw satellite availability
# --------------------------------------------------
print("\n[0/4] Raw satellite availability at BASE")

gnss_edu.plot_tracking_timeline(
    df_base_pivot,
    sampling=sampling,
    title="Raw tracking timeline at BASE (present / missing)",
)

# --------------------------------------------------
# Level 1 — Elevation-driven availability
# --------------------------------------------------
print("\n[1/4] Elevation-based availability")
print("Color = elevation angle. The cutoff removes the low-elevation arc edges.")

plot_tracking_timeline_colored(
    df_base_pivot,
    value_col=elev_col,
    sampling=sampling,
    title=f"Satellite tracking colored by elevation ({elev_col})",
    cmap="viridis",
    vmin=0.0,
    vmax=90.0,
)

plot_tracking_timeline_from_mask(
    df_base_pivot,
    mask_col="pivot_elev_ok",
    sampling=sampling,
    title=f"Tracking above elevation cutoff ({elev_cutoff_deg:.0f} deg)",
)

# --------------------------------------------------
# Level 2 — SNR refinement after elevation mask
# --------------------------------------------------
print("\n[2/4] SNR refinement after elevation mask")

gnss_edu.plot_gnss_timeseries_by_prn(
    df_base_pivot,
    y=snr_col,
    gap=gap_arc,
    label_arcs=True,
    show_legend=False,
    title=f"{snr_col} time series by satellite PRN (BASE)",
)

plot_tracking_timeline_from_mask(
    df_base_pivot,
    mask_col="pivot_elev_snr_ok",
    sampling=sampling,
    title=(
        f"Tracking after elevation + SNR filters "
        f"({elev_cutoff_deg:.0f} deg, {snr_col} >= {snr_min:g})"
    ),
)

# --------------------------------------------------
# Level 3 — Observable compatibility
# --------------------------------------------------
print("\n[3/4] Observable compatibility for the target treatment")
print(
    "At the BASE station, a pivot must carry all observables required "
    "by the DD or combination to be formed."
)
print(
    "For this short-baseline setup, we assume that a pivot compatible at BASE "
    "is also a relevant candidate for the nearby rover."
)

plot_tracking_timeline_from_mask(
    df_base_pivot,
    mask_col="pivot_usable",
    sampling=sampling,
    title=(
        "Tracking after elevation + SNR + observable compatibility "
        f"{required_observables}"
    ),
)

obs_compat_per_prn = (
    df_base_pivot.groupby(level="prn")["pivot_usable"]
    .sum()
    .sort_values(ascending=False)
    .rename("usable_epochs")
    .to_frame()
)
print("\nUsable epochs per PRN after all cumulative filters:")
display(obs_compat_per_prn)

# --------------------------------------------------
# Level 4 — Minimal coverage + stable pivot schedule
# --------------------------------------------------
print("\n[4/4] Minimal pivot set cover and stable active schedule")

selected_prns, diag = gnss_edu.greedy_pivot_set_cover(
    df=df_base_pivot,
    snr_col=support_col,
    snr_min=0.0,
    sampling=sampling,
    require_full_coverage=True,
    return_diagnostics=True,
)

print("\nCandidate pivot satellites:", selected_prns)
print("Number of pivots selected :", diag["n_selected"])
print("Coverage ratio            :", diag["coverage_ratio"])
print("Uncovered epochs          :", diag["n_uncovered"])

active_pivot = gnss_edu.build_active_pivot_schedule(
    df=df_base_pivot,
    selected_prns=selected_prns,
    snr_col=support_col,
    snr_min=0.0,
    sampling=sampling,
    switch_margin_db=switch_margin_db,
)

cov = gnss_edu.check_full_coverage_from_active_pivot(active_pivot)
print("Coverage check before stability post-processing:", cov)
assert cov["n_none"] == 0, (
    "Active pivot has gaps after the cumulative filters. "
    "Relax one of the filters or revise the observable requirement."
)

active_pivot = gnss_edu.fix_short_pivot_segments(
    df=df_base_pivot,
    active_pivot=active_pivot,
    snr_col=support_col,
    snr_min=0.0,
    sampling=sampling,
    min_duration=PIVOT_MIN_DURATION,
    pool="selected_only",
    selected_prns=selected_prns,
)

active_pivot = gnss_edu.enforce_min_dwell(
    df=df_base_pivot,
    active_pivot=active_pivot,
    snr_col=support_col,
    snr_min=0.0,
    sampling=sampling,
    min_dwell=PIVOT_MIN_DWELL,
)

cov = gnss_edu.check_full_coverage_from_active_pivot(active_pivot)
print("Coverage check after stability post-processing:", cov)
assert cov["n_none"] == 0, (
    "Coverage was broken by pivot post-processing. "
    "Relax the stability constraints or the candidate filters."
)

fig, ax, info = gnss_edu.plot_tracking_timeline_with_pivots(
    df=df_base_pivot,
    selected_prns=selected_prns,
    sampling=sampling,
    snr_col=support_col,
    snr_min=0.0,
    active_pivot=active_pivot,
    title=(
        "Final pivot strategy after elevation + SNR + observable compatibility "
        "+ stable scheduling"
    ),
)


# %%
###############################################################################
# Pivot schedule -> segments -> diagnostics
###############################################################################

segments = gnss_edu.pivot_schedule_to_segments(
    active_pivot,
    sampling=sampling,
    drop_none=True
)

print("\nNumber of pivot segments:", len(segments))
display(segments)

diag_cov = gnss_edu.check_full_coverage_from_active_pivot(active_pivot)
print("\nCoverage diagnostics:")
print(diag_cov)

segments2 = segments.copy()
segments2["is_short"] = segments2["duration"] < PIVOT_MIN_DWELL

pivot_seg_stats = (
    segments2.groupby("prn")
    .agg(
        n_segments=("prn", "size"),
        min_duration=("duration", "min"),
        median_duration=("duration", "median"),
        mean_duration=("duration", "mean"),
        max_duration=("duration", "max"),
        n_short=("is_short", "sum"),
    )
)

pivot_seg_stats["short_fraction"] = (
    pivot_seg_stats["n_short"] / pivot_seg_stats["n_segments"]
)

print(f"\nContinuous pivot segment statistics (short < {PIVOT_MIN_DWELL}):")
display(
    pivot_seg_stats.sort_values(
        ["short_fraction", "n_short", "median_duration", "n_segments"],
        ascending=[False, False, True, False]
    )
)

short_segments = (
    segments2.loc[segments2["is_short"]]
    .sort_values("duration")
    .assign(duration_min=lambda df: df["duration"].dt.total_seconds() / 60.0)
)

print("\nShortest pivot segments (most suspicious pivot choices):")
display(short_segments[["prn", "start", "end", "duration", "duration_min", "n_epochs"]])

print("\nPedagogical discussion point:")
print(
    "Up to what rover-base distance does a pivot strategy defined from the "
    "BASE perspective remain a reasonable approximation?"
)


# %%

#            satellite s
# BASE  -------------------------
# ROVER -------------------------

# Single difference:   SD_s = obs_rover − obs_base       <------ OK

# Choose a pivot p                                       <------ OK

# Double difference:   DD_s = SD_s − SD_p                <------ NOW

def compute_double_differences_from_sd_rich(
    df_SD: pd.DataFrame,
    active_pivot: pd.Series,
    sd_cols: list[str] | None = None,
    keep_pivot_rows: bool = True,
    add_pivot_sd: bool = True,
    pivot_sd_prefix: str = "PIV_",
    dd_prefix: str = "DD_",
) -> pd.DataFrame:
    """
    Build a *rich* Double Differences (DD) table from Single Differences (SD)
    using a time-varying pivot schedule.

    Outputs include:
    - pivot_prn : pivot used at each epoch
    - is_pivot  : True if the row corresponds to the pivot PRN
    - DD_*      : DD observables for all requested SD columns
    - PIV_*     : (optional) pivot SD values replicated on each row of the epoch

    Parameters
    ----------
    df_SD : DataFrame
        MultiIndex (epoch, prn). Columns contain SD observables (e.g., SD_L1, SD_C1, ...).
    active_pivot : Series
        Indexed by epoch, values are pivot PRNs (e.g., "G05") or None/NaN.
        Must cover the epochs to be processed (no None for full DD computation).
    sd_cols : list[str] | None
        Which SD columns to difference. If None, uses all columns in df_SD.
    keep_pivot_rows : bool
        If True, keep pivot PRN rows (DD will be 0 by construction).
        If False, drop pivot PRN rows.
    add_pivot_sd : bool
        If True, add columns with pivot SD values replicated per epoch (PIV_*).
    pivot_sd_prefix : str
        Prefix for pivot SD columns (default "PIV_").
    dd_prefix : str
        Prefix for DD columns (default "DD_").

    Returns
    -------
    df_DD : DataFrame
        MultiIndex (epoch, prn) with metadata columns + DD columns (+ optional PIV columns).
    """
    # ----------------------- checks -----------------------
    if not isinstance(df_SD.index, pd.MultiIndex):
        raise ValueError("df_SD must have MultiIndex ('epoch','prn').")
    if "epoch" not in df_SD.index.names or "prn" not in df_SD.index.names:
        raise ValueError("df_SD index must have levels ('epoch','prn').")

    if sd_cols is None:
        sd_cols = list(df_SD.columns)

    missing = [c for c in sd_cols if c not in df_SD.columns]
    if missing:
        raise ValueError(f"Some requested SD columns are missing in df_SD: {missing}")

    # Epoch grid to process = epochs present in df_SD
    epochs = pd.DatetimeIndex(sorted(df_SD.index.get_level_values("epoch").unique()))
    piv = active_pivot.reindex(epochs)

    if piv.isna().any():
        n = int(piv.isna().sum())
        first_bad = piv[piv.isna()].index[0]
        raise RuntimeError(
            f"active_pivot has {n} None/NaN epochs after alignment to df_SD. "
            f"First missing pivot at epoch={first_bad}. "
            "DD requires a pivot at each processed epoch."
        )

    # -------------------- main join -----------------------
    df = df_SD[sd_cols].copy()
    df = df.join(piv.rename("pivot_prn").to_frame(), on="epoch")

    # Reset for easier pivot-row extraction
    dfr = df.reset_index()  # epoch, prn, SD..., pivot_prn

    # Identify pivot rows
    dfr["is_pivot"] = dfr["prn"] == dfr["pivot_prn"]

    pivot_rows = dfr[dfr["is_pivot"]]
    if pivot_rows.empty:
        raise RuntimeError(
            "No pivot rows found in df_SD: pivot PRNs are not present for the corresponding epochs."
        )

    # Pivot SD values indexed by epoch
    pivot_vals = pivot_rows.set_index("epoch")[sd_cols]

    # Join pivot SD back to every row of the same epoch
    dfr = dfr.join(pivot_vals, on="epoch", rsuffix="_PIV")

    # --------------------- compute DD ---------------------
    out = {}
    for c in sd_cols:
        # Use a readable DD name: DD_L1 from SD_L1, DD_C1 from SD_C1, etc.
        name = dd_prefix + c.removeprefix("SD_")
        out[name] = dfr[c] - dfr[c + "_PIV"]

        if add_pivot_sd:
            out[pivot_sd_prefix + c.removeprefix("SD_")] = dfr[c + "_PIV"]

    # Metadata
    out["pivot_prn"] = dfr["pivot_prn"]
    out["is_pivot"] = dfr["is_pivot"]

    df_DD = pd.DataFrame(out)
    df_DD["epoch"] = dfr["epoch"].values
    df_DD["prn"] = dfr["prn"].values
    df_DD = df_DD.set_index(["epoch", "prn"]).sort_index()

    if not keep_pivot_rows:
        df_DD = df_DD.loc[~df_DD["is_pivot"]]

    return df_DD



def add_pivot_change_flag(df_DD: pd.DataFrame) -> pd.DataFrame:
    """
    Add a boolean flag 'pivot_changed' True at epochs where pivot differs from previous epoch.
    """
    if "pivot_prn" not in df_DD.columns:
        raise ValueError("df_DD must contain 'pivot_prn'.")

    out = df_DD.copy()
    piv = out["pivot_prn"].groupby(level="epoch").first()
    changed = piv.ne(piv.shift(1)).rename("pivot_changed")
    out = out.join(changed.to_frame(), on="epoch")
    return out



# %%
###############################################################################
# What is df_DD?
#
# The dataframe `df_DD` contains the GNSS Double Differences (DD) computed
# from the Single Differences (SD) using the pivot satellite selected at
# each epoch.
#
# Index
# -----
# df_DD uses a MultiIndex (epoch, prn):
#
#   epoch : observation time
#   prn   : satellite identifier
#
# Each row therefore corresponds to one satellite observed at one epoch.
#
# Pivot information
# -----------------
# Several columns describe how the DD were formed:
#
#   pivot_prn     : pivot satellite used at that epoch
#   is_pivot      : True if the row corresponds to the pivot satellite itself
#   pivot_changed : True when the pivot satellite changes compared to the
#                   previous epoch
#
# When `is_pivot=True`, the DD values are zero by construction because:
#
#       DD_pivot = SD_pivot − SD_pivot = 0
#
# Observables
# -----------
# For each SD observable available in df_SD, the dataframe contains the
# corresponding double difference:
#
#       DD_obs = SD_obs(satellite) − SD_obs(pivot)
#
# Example:
#
#       DD_L1 = SD_L1(satellite) − SD_L1(pivot)
#
# Why track pivot changes?
# ------------------------
# When the pivot satellite changes, the definition of the double differences
# also changes. In carrier phase processing this typically introduces a new
# set of ambiguity parameters. The column `pivot_changed` therefore marks
# epochs where the DD geometry changes.
#
###############################################################################

# %%
###############################################################################
# Double Differences (DD) from SD using a time-varying pivot schedule
# (rich output: pivot info + all DD + optional pivot SD)
###############################################################################

print("\n===== Build rich DD table =====")

df_DD = compute_double_differences_from_sd_rich(
    df_SD=df_SD,
    active_pivot=active_pivot,
    sd_cols=None,              # <= all SD columns available in df_SD
    keep_pivot_rows=True,      # keep pivot rows (DD=0), useful for teaching
    add_pivot_sd=True,         # add PIV_* columns for traceability
)

print("df_DD columns (first 20):")
print(list(df_DD.columns)[:20])
print("... total columns:", len(df_DD.columns))

print("\nSanity checks:")
print("Unique pivots used:", df_DD["pivot_prn"].unique()[:20])
print("Pivot rows:", int(df_DD["is_pivot"].sum()))
print("Non-pivot rows:", int((~df_DD["is_pivot"]).sum()))

print("\nExample rows:")
print(df_DD.head(10))

df_DD = add_pivot_change_flag(df_DD)
print("Number of pivot change epochs:", int(df_DD["pivot_changed"].groupby(level="epoch").first().sum()))
