#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step 02 - Code-based GNSS positioning with progressive model enrichment

Educational scope
-----------------
This second step introduces code-based GNSS positioning from real observations.
The pedagogical objective is to start from a deliberately simplified model and
progressively enrich it with the main physical effects that govern GNSS code
measurements.

Learning goals
--------------
By the end of this script, students should be able to:
1. load GNSS observation data,
2. extract a code observable suitable for positioning,
3. retrieve the approximate receiver position from the RINEX header,
4. build a clean working DataFrame for code-based positioning,
5. enrich the observation table with precise satellite information,
6. observe how successive model refinements improve the estimated position.

Pedagogical position in the curriculum
--------------------------------------
Step 01 focused on reading RINEX files and structuring GNSS observations.

Step 02 focuses on transforming GNSS code observations into a positioning
model. The key pedagogical idea is that the estimated position should
progressively approach the approximate receiver position given in the RINEX
header as the model becomes more realistic.

Important note
--------------
In this script, the approximate receiver position stored in the RINEX header is
used as a practical pedagogical reference. It is not treated as a rigorous
ground-truth product, but as a convenient benchmark to visualize how the
solution improves as additional physical effects are introduced.

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
- os
- pathlib
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

import numpy as np
import pandas as pd

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

# Two permanent GNSS stations about 10 km apart in the Paris region
STATION_DICT = {"rgp": ["SMNE", "MLVL"]}

# First code observable to try
CODE_OBSERVABLE = "C1"
CODE_OBSERVABLE_FALLBACKS = ["P1", "C2", "P2"]

# Keep only GPS for this first pedagogical version
SATELLITE_SYSTEM = "G"

# Convergence threshold for the iterative least-squares adjustment
POSITION_CONVERGENCE_THRESHOLD_M = 1.0

print("Configuration loaded.")
print("Processing date            :", PROCESSING_DATE)
print("Working directory          :", WORK_DIR)
print("Satellite system           :", SATELLITE_SYSTEM)
print("Preferred code observable  :", CODE_OBSERVABLE)


# %%
###############################################################################
# Pedagogical roadmap
###############################################################################

print("Step 02 roadmap")
print("----------------")
print("M0 : naive code-based positioning")
print("M1 : add satellite clock and relativistic corrections")
print("M2 : estimate one receiver clock parameter per epoch")
print("M3 : add Earth-rotation (Sagnac) correction")
print("Later extensions: ionosphere, troposphere, antenna effects")

# This dictionary stores the successive positioning solutions.
# It is the main pedagogical backbone of Step 02.
solutions = {}


# %%
###############################################################################
# Helper functions
###############################################################################

def extract_download_path(entry):
    """Return the local file path from a geodezyx download entry."""
    if isinstance(entry, tuple):
        return entry[0]
    return entry


def build_receiver_clock_block(index: pd.MultiIndex) -> tuple[np.ndarray, np.ndarray]:
    """
    Build the receiver-clock indicator block: one clock parameter per epoch.

    Returns
    -------
    block_dt_r : np.ndarray
        Indicator matrix of shape (n_obs, n_epochs).
    epoch_unique : np.ndarray
        Ordered unique epochs.
    """
    epoch_codes, epoch_unique = pd.factorize(index.get_level_values("epoch"))
    block_dt_r = np.zeros((len(index), len(epoch_unique)))
    block_dt_r[np.arange(len(index)), epoch_codes] = 1.0
    return block_dt_r, epoch_unique


def store_solution(
    key: str,
    description: str,
    P_est: np.ndarray,
    P_ref: np.ndarray,
    residual_vector: np.ndarray,
    n_iterations: int,
    active_effects: list[str],
) -> None:
    """
    Store one positioning solution in the pedagogical solutions dictionary.
    """
    E, N, U = conv.xyz2enu(
        P_ref[0], P_ref[1], P_ref[2],
        P_est[0], P_est[1], P_est[2]
    )
    distance_to_header = float(np.linalg.norm(P_est - P_ref))
    residual_rms = float(np.sqrt(np.mean(residual_vector ** 2)))

    solutions[key] = {
        "description": description,
        "receiver_xyz_m": P_est.copy(),
        "enu_error_m": np.array([E, N, U], dtype=float),
        "distance_to_rinex_header_m": distance_to_header,
        "residual_rms_m": residual_rms,
        "n_iterations": int(n_iterations),
        "active_effects": active_effects.copy(),
    }

    print(f"{key} summary")
    print("-" * (len(key) + 8))
    print("Description                    :", description)
    print("Distance to RINEX header [m]   :", distance_to_header)
    print("ENU error [m]                  :", np.array([E, N, U]))
    print("Residual RMS [m]               :", residual_rms)
    print("Iterations                     :", n_iterations)
    print()


def summarize_solutions(solutions_dict: dict) -> pd.DataFrame:
    """
    Convert the pedagogical solution dictionary into a summary DataFrame.
    """
    rows = []
    for model_name, result in solutions_dict.items():
        rows.append({
            "model": model_name,
            "description": result["description"],
            "distance_to_rinex_header_m": result["distance_to_rinex_header_m"],
            "residual_rms_m": result["residual_rms_m"],
            "east_m": result["enu_error_m"][0],
            "north_m": result["enu_error_m"][1],
            "up_m": result["enu_error_m"][2],
            "n_iterations": result["n_iterations"],
        })
    return pd.DataFrame(rows)


# %%
###############################################################################
# Create the working directory
###############################################################################

WORK_DIR.mkdir(parents=True, exist_ok=True)

print("Working directory is ready:")
print(WORK_DIR)


# %%
###############################################################################
# Download sample RINEX observation files
###############################################################################

download_output = operational.download_gnss_rinex(
    statdico=STATION_DICT,
    output_dir=str(WORK_DIR),
    startdate=PROCESSING_DATE,
    enddate=PROCESSING_DATE,
    parallel_download=1,
)

print("Download output:")
print(download_output)


# %%
###############################################################################
# Identify the BASE and ROVER files
###############################################################################

base_rinex_file = extract_download_path(download_output[0])
rover_rinex_file = extract_download_path(download_output[1])

print("BASE RINEX file :")
print(base_rinex_file)
print()
print("ROVER RINEX file:")
print(rover_rinex_file)


# %%
###############################################################################
# Read the ROVER RINEX observation file
#
# Educational objective
# ---------------------
# In this step, we position one receiver from its code observations.
###############################################################################

df_obs, rinex_header = files_rw.read_rinex_obs(
    rover_rinex_file,
    return_header=True,
    set_index=["epoch", "prn"],
)

print("Observation DataFrame loaded.")
print("Shape      :", df_obs.shape)
print("Index type :", type(df_obs.index).__name__)
print("Index names:", df_obs.index.names)
print("First columns:", list(df_obs.columns[:12]))


# %%
###############################################################################
# Extract the approximate receiver position from the RINEX header
###############################################################################

approx_position_lines = [
    line for line in rinex_header if "APPROX POSITION XYZ" in line
]

if len(approx_position_lines) == 0:
    raise RuntimeError("Approximate receiver position not found in RINEX header.")

approx_receiver_xyz = np.array(
    approx_position_lines[0].split()[:3],
    dtype=float,
)

approx_receiver_flh = conv.XYZ2GEO(*approx_receiver_xyz)
lat_deg, lon_deg, h_m = approx_receiver_flh

print("Approximate receiver position from RINEX header")
print("-----------------------------------------------")
print("XYZ [m]           :", approx_receiver_xyz)
print("Latitude  [deg]   :", lat_deg)
print("Longitude [deg]   :", lon_deg)
print("Height    [m]     :", h_m)


# %%
###############################################################################
# Keep only GPS observations
###############################################################################

prn_values = df_obs.index.get_level_values("prn").astype(str)
gps_mask = prn_values.str.startswith(SATELLITE_SYSTEM)
df_obs = df_obs.loc[gps_mask].copy()

print("GPS-only observation table:")
print("Shape:", df_obs.shape)
print("PRNs :", df_obs.index.get_level_values("prn").unique().tolist())


# %%
###############################################################################
# Normalize missing values and convert observables to numeric
###############################################################################

df_obs = df_obs.replace(
    to_replace=[r"^\s*$", "nan", "NaN", "NA", "N/A", "null", "None"],
    value=np.nan,
    regex=True,
)

for col in df_obs.columns:
    if col not in {"sys", "prn", "prni", "row_id"}:
        df_obs[col] = pd.to_numeric(df_obs[col], errors="coerce")

empty_cols = df_obs.columns[~df_obs.notna().any(axis=0)].tolist()
df_obs = df_obs.drop(columns=empty_cols)

print("Columns removed because they were fully empty:")
print(empty_cols if empty_cols else "None")


# %%
###############################################################################
# Select the code observable used for the first positioning model
###############################################################################

candidate_code_observables = [CODE_OBSERVABLE] + CODE_OBSERVABLE_FALLBACKS
available_code_observables = [
    obs for obs in candidate_code_observables if obs in df_obs.columns
]

if len(available_code_observables) == 0:
    raise RuntimeError(
        f"No suitable code observable found. Tried: {candidate_code_observables}"
    )

code_obs_used = available_code_observables[0]

print("Code observable selected for the first model:", code_obs_used)


# %%
###############################################################################
# Build the working DataFrame for code-based positioning
#
# Educational objective
# ---------------------
# Keep only the information needed for the first code-based model.
#
# Pedagogical note
# ----------------
# From this point on, df_code is the main table of the TP. Intermediate flat
# tables are used only temporarily and are deleted afterward to keep the
# Spyder Variable Explorer readable.
###############################################################################

df_code = df_obs[[code_obs_used]].copy()
df_code = df_code.rename(columns={code_obs_used: "code_m"})
df_code = df_code.dropna(subset=["code_m"]).copy()
df_code["row_id"] = np.arange(len(df_code))

print("Code-based working DataFrame")
print("----------------------------")
print("Observable used :", code_obs_used)
print("Shape           :", df_code.shape)
print("Columns         :", df_code.columns.tolist())
print()
print(df_code.head())


# %%
###############################################################################
# Inspect the code measurements
###############################################################################

print("Number of epochs      :", df_code.index.get_level_values("epoch").nunique())
print("Number of satellites  :", df_code.index.get_level_values("prn").nunique())
print("Number of code rows   :", len(df_code))

obs_per_epoch = df_code.groupby(level="epoch").size()

print("\nObservation count per epoch")
print("---------------------------")
print(obs_per_epoch.describe())


# %%
###############################################################################
# Prepare the initial receiver coordinates
###############################################################################

receiver_xyz_approx = approx_receiver_xyz.copy()

print("Initial receiver coordinates used for the iteration [m]:")
print(receiver_xyz_approx)


# %%
###############################################################################
# Download precise satellite products
###############################################################################

download_products = operational.download_gnss_products(
    archive_dir=str(WORK_DIR),
    startdate=PROCESSING_DATE,
    enddate=PROCESSING_DATE,
    archtype="year/doy",
    AC_names=("IGS",),
    repro=0,
    archive_center="ign",
    parallel_download=1,
)

print("Precise product download output:")
print(download_products)


# %%
###############################################################################
# Load the SP3 precise orbit/clock product
###############################################################################

sp3_path = extract_download_path(download_products[0])

print("Using SP3 file:")
print(sp3_path)

df_sp3 = files_rw.read_sp3(
    sp3_path,
    returns_pandas=True,
    new_col_names=True,
)

print("SP3 DataFrame loaded.")
print("Shape      :", df_sp3.shape)
print("Columns    :", list(df_sp3.columns[:12]))


# %%
###############################################################################
# Enrich df_code with satellite state at signal emission time
#
# Educational objective
# ---------------------
# Attach one satellite state to each code observation.
#
# Method
# ------
# For each PRN:
#   1. estimate the signal flight time from the pseudorange,
#   2. derive an approximate emission epoch,
#   3. interpolate the SP3 orbit at that epoch,
#   4. compute the relativistic correction,
#   5. refine the emission epoch using the satellite clock,
#   6. store the final satellite state in the MultiIndex working table.
###############################################################################

df_tmp = df_code.reset_index().copy()

df_tmp["X_sat"] = np.nan
df_tmp["Y_sat"] = np.nan
df_tmp["Z_sat"] = np.nan
df_tmp["dte_sat"] = np.nan
df_tmp["dRelat"] = np.nan

print("***** Satellite-state computation at signal emission time *****")

for prn in df_tmp["prn"].unique():

    prn_mask = df_tmp["prn"] == prn
    df_tmp_prn = df_tmp.loc[prn_mask].copy()
    df_sp3_prn = df_sp3[df_sp3["prn"] == prn].copy()

    if df_sp3_prn.empty:
        print(f"PRN {prn} not found in SP3 product: skipped")
        continue

    fly_time = pd.to_timedelta(
        df_tmp_prn["code_m"] / conv.SPEED_OF_LIGHT,
        unit="s",
    )

    t_rec = df_tmp_prn["epoch"]
    t_emit_approx = t_rec - fly_time

    orb_df_approx = reffram.orb_df_lagrange_interpolate(
        df_sp3_prn,
        t_emit_approx.values,
    ).copy()
    orb_df_approx[["x", "y", "z"]] = orb_df_approx[["x", "y", "z"]] * 1e3

    delta_t = pd.to_timedelta(1e-3, unit="s")

    orb_df_fwd = reffram.orb_df_lagrange_interpolate(
        df_sp3_prn,
        t_emit_approx.values + delta_t,
    ).copy()
    orb_df_bwd = reffram.orb_df_lagrange_interpolate(
        df_sp3_prn,
        t_emit_approx.values - delta_t,
    ).copy()

    orb_df_fwd[["x", "y", "z"]] = orb_df_fwd[["x", "y", "z"]] * 1e3
    orb_df_bwd[["x", "y", "z"]] = orb_df_bwd[["x", "y", "z"]] * 1e3

    v_xyz = (orb_df_fwd[["x", "y", "z"]] - orb_df_bwd[["x", "y", "z"]]) / (
        2 * delta_t.total_seconds()
    )
    r_xyz = orb_df_approx[["x", "y", "z"]]

    dRelat_prn = -2.0 * (v_xyz * r_xyz).sum(axis=1) / (conv.SPEED_OF_LIGHT ** 2)

    t_emit_refined = t_emit_approx - pd.to_timedelta(
        orb_df_approx["clk"].values,
        unit="us",
    )

    orb_df_refined = reffram.orb_df_lagrange_interpolate(
        df_sp3_prn,
        t_emit_refined.values,
    ).copy()
    orb_df_refined[["x", "y", "z"]] = orb_df_refined[["x", "y", "z"]] * 1e3

    df_tmp.loc[prn_mask, "X_sat"] = orb_df_refined["x"].values
    df_tmp.loc[prn_mask, "Y_sat"] = orb_df_refined["y"].values
    df_tmp.loc[prn_mask, "Z_sat"] = orb_df_refined["z"].values
    df_tmp.loc[prn_mask, "dte_sat"] = orb_df_refined["clk"].values * 1e-6
    df_tmp.loc[prn_mask, "dRelat"] = dRelat_prn.values

rows_before = len(df_tmp)
df_tmp = df_tmp.dropna(subset=["X_sat", "Y_sat", "Z_sat", "dte_sat", "dRelat"]).copy()
rows_after = len(df_tmp)

df_code = df_tmp.set_index(["epoch", "prn"]).sort_index()

print("Satellite-state filtering summary")
print("---------------------------------")
print("Rows before filtering :", rows_before)
print("Rows after filtering  :", rows_after)
print("Rows removed          :", rows_before - rows_after)
print()
print(df_code.head())

# Keep only the main MultiIndex working DataFrame visible in the workspace.
del df_tmp, df_tmp_prn, df_sp3_prn
del fly_time, t_rec, t_emit_approx, t_emit_refined
del orb_df_approx, orb_df_fwd, orb_df_bwd, orb_df_refined
del v_xyz, r_xyz, dRelat_prn, delta_t, prn_mask, prn
gc.collect()


# %%
###############################################################################
# Status before the first positioning model
###############################################################################

print("df_code is now ready for the first code-based positioning model.")
print("Columns available:")
print(df_code.columns.tolist())


# %%
###############################################################################
# Model M0 - Naive code-based positioning
#
# Educational objective
# ---------------------
# Start with a deliberately simplified geometric model:
#   code_m ≈ geometric_range
#
# This model ignores satellite clock, receiver clock, relativistic, and
# Earth-rotation effects. It is intentionally incomplete.
###############################################################################

print("***** M0 - Naive code-based positioning *****")

P_app = np.array([0.0, 0.0, 0.0])
dP_est = np.array([100.0, 100.0, 100.0])
n_iter = 0

while np.linalg.norm(dP_est) > POSITION_CONVERGENCE_THRESHOLD_M:
    distances = np.sqrt(
        (df_code["X_sat"].values - P_app[0]) ** 2 +
        (df_code["Y_sat"].values - P_app[1]) ** 2 +
        (df_code["Z_sat"].values - P_app[2]) ** 2
    )

    B = df_code["code_m"].values - distances

    dX = (P_app[0] - df_code["X_sat"].values) / distances
    dY = (P_app[1] - df_code["Y_sat"].values) / distances
    dZ = (P_app[2] - df_code["Z_sat"].values) / distances
    A = np.column_stack((dX, dY, dZ))

    dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
    P_est = P_app + dP_est

    n_iter += 1
    print(f"Iteration {n_iter}: X={P_est[0]:.3f}, Y={P_est[1]:.3f}, Z={P_est[2]:.3f}")
    P_app = P_est

gnss_edu.plot_residual_analysis(
    A,
    B,
    dP_est,
    figure_title="M0 - Naive code-based positioning",
    save_path=WORK_DIR / "01_M0_naive_code.png",
    P_est=P_est,
    P_rnx_header=approx_receiver_xyz,
)

store_solution(
    key="M0_naive_code",
    description="Code-only positioning with no clock or Sagnac correction",
    P_est=P_est,
    P_ref=approx_receiver_xyz,
    residual_vector=B - A @ dP_est,
    n_iterations=n_iter,
    active_effects=["satellite positions at transmission time"],
)

del P_app, dP_est, P_est, n_iter, distances, dX, dY, dZ, A, B
gc.collect()


# %%
###############################################################################
# Numerical aside - Least-squares solvers
#
# Educational objective
# ---------------------
# Compare different numerical ways to solve the same least-squares problem.
#
# Pedagogical note
# ----------------
# The explicit inverse is useful as a lecture formula, but it is usually not
# the preferred numerical implementation in scientific software.
###############################################################################

print("***** Numerical aside - inv vs solve vs lstsq *****")

P_app = np.array([0.0, 0.0, 0.0])

distances = np.sqrt(
    (df_code["X_sat"].values - P_app[0]) ** 2 +
    (df_code["Y_sat"].values - P_app[1]) ** 2 +
    (df_code["Z_sat"].values - P_app[2]) ** 2
)

B = df_code["code_m"].values - distances
dX = (P_app[0] - df_code["X_sat"].values) / distances
dY = (P_app[1] - df_code["Y_sat"].values) / distances
dZ = (P_app[2] - df_code["Z_sat"].values) / distances
A = np.column_stack((dX, dY, dZ))

x_inv = np.linalg.inv(A.T @ A) @ (A.T @ B)
x_solve = np.linalg.solve(A.T @ A, A.T @ B)
x_lstsq, _, _, _ = np.linalg.lstsq(A, B, rcond=None)

print("Norm difference inv  - solve :", np.linalg.norm(x_inv - x_solve))
print("Norm difference inv  - lstsq :", np.linalg.norm(x_inv - x_lstsq))
print("Norm difference solve- lstsq :", np.linalg.norm(x_solve - x_lstsq))

del P_app, distances, B, dX, dY, dZ, A, x_inv, x_solve, x_lstsq
gc.collect()


# %%
###############################################################################
# Model M1 - Add satellite clock and relativistic corrections
#
# Observation model
# -----------------
#   code_m ≈ geometric_range - c (dt_sat + dRelat)
###############################################################################

print("***** M1 - Satellite clock and relativistic corrections *****")

P_app = np.array([0.0, 0.0, 0.0])
dP_est = np.array([100.0, 100.0, 100.0])
n_iter = 0

while np.linalg.norm(dP_est) > POSITION_CONVERGENCE_THRESHOLD_M:
    distances = np.sqrt(
        (df_code["X_sat"].values - P_app[0]) ** 2 +
        (df_code["Y_sat"].values - P_app[1]) ** 2 +
        (df_code["Z_sat"].values - P_app[2]) ** 2
    )

    B = (
        df_code["code_m"].values
        - distances
        + conv.SPEED_OF_LIGHT * (df_code["dte_sat"].values + df_code["dRelat"].values)
    )

    dX = (P_app[0] - df_code["X_sat"].values) / distances
    dY = (P_app[1] - df_code["Y_sat"].values) / distances
    dZ = (P_app[2] - df_code["Z_sat"].values) / distances
    A = np.column_stack((dX, dY, dZ))

    dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
    P_est = P_app + dP_est

    n_iter += 1
    print(f"Iteration {n_iter}: X={P_est[0]:.3f}, Y={P_est[1]:.3f}, Z={P_est[2]:.3f}")
    P_app = P_est

gnss_edu.plot_residual_analysis(
    A,
    B,
    dP_est,
    figure_title="M1 - Satellite clock and relativistic corrections",
    save_path=WORK_DIR / "02_M1_satellite_clock.png",
    P_est=P_est,
    P_rnx_header=approx_receiver_xyz,
)

store_solution(
    key="M1_satellite_clock",
    description="Code positioning with satellite clock and relativistic corrections",
    P_est=P_est,
    P_ref=approx_receiver_xyz,
    residual_vector=B - A @ dP_est,
    n_iterations=n_iter,
    active_effects=[
        "satellite positions at transmission time",
        "satellite clock correction",
        "relativistic correction",
    ],
)

del P_app, dP_est, P_est, n_iter, distances, dX, dY, dZ, A, B
gc.collect()


# %%
###############################################################################
# Model M2 - Add one receiver clock parameter per epoch
#
# Observation model
# -----------------
#   code_m ≈ geometric_range - c(dt_sat + dRelat) + dt_r_epoch
#
# Educational objective
# ---------------------
# Show that the receiver clock is not a negligible nuisance parameter.
###############################################################################

print("***** M2 - Add receiver clock estimation *****")

block_dt_r, epoch_unique = build_receiver_clock_block(df_code.index)

P_app = np.array([0.0, 0.0, 0.0])
dP_est = np.array([100.0, 100.0, 100.0])
n_iter = 0

while np.linalg.norm(dP_est[0:3]) > POSITION_CONVERGENCE_THRESHOLD_M:
    distances = np.sqrt(
        (df_code["X_sat"].values - P_app[0]) ** 2 +
        (df_code["Y_sat"].values - P_app[1]) ** 2 +
        (df_code["Z_sat"].values - P_app[2]) ** 2
    )

    B = (
        df_code["code_m"].values
        - distances
        + conv.SPEED_OF_LIGHT * (df_code["dte_sat"].values + df_code["dRelat"].values)
    )

    dX = (P_app[0] - df_code["X_sat"].values) / distances
    dY = (P_app[1] - df_code["Y_sat"].values) / distances
    dZ = (P_app[2] - df_code["Z_sat"].values) / distances
    A = np.column_stack((dX, dY, dZ, block_dt_r))

    dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
    P_est = P_app + dP_est[0:3]

    n_iter += 1
    print(f"Iteration {n_iter}: X={P_est[0]:.3f}, Y={P_est[1]:.3f}, Z={P_est[2]:.3f}")
    P_app = P_est

gnss_edu.plot_residual_analysis(
    A,
    B,
    dP_est,
    figure_title="M2 - Satellite and receiver clock corrections",
    save_path=WORK_DIR / "03_M2_satellite_receiver_clock.png",
    P_est=P_est,
    P_rnx_header=approx_receiver_xyz,
)

store_solution(
    key="M2_receiver_clock",
    description="Code positioning with satellite clock, relativistic correction, and one receiver clock per epoch",
    P_est=P_est,
    P_ref=approx_receiver_xyz,
    residual_vector=B - A @ dP_est,
    n_iterations=n_iter,
    active_effects=[
        "satellite positions at transmission time",
        "satellite clock correction",
        "relativistic correction",
        "receiver clock estimation",
    ],
)

del P_app, dP_est, P_est, n_iter, distances, dX, dY, dZ, A, B
del block_dt_r, epoch_unique
gc.collect()


# %%
###############################################################################
# Model M3 - Add Earth-rotation (Sagnac) correction
#
# Educational objective
# ---------------------
# Express the satellite coordinates in the Earth-fixed frame at reception time.
#
# Physical idea
# -------------
# During the signal travel time, the Earth rotates. To compare the receiver and
# the satellite in the same frame, the satellite coordinates must be rotated
# back by the Earth rotation accumulated during signal propagation.
###############################################################################

print("***** M3 - Add Earth-rotation (Sagnac) correction *****")

tau_s = df_code["code_m"].to_numpy(dtype=float) / conv.SPEED_OF_LIGHT
dtheta = tau_s * conv.EARTH_ROTATION_MEAN_ANGULAR_VELOCITY
alpha = -dtheta

cos_a = np.cos(alpha)
sin_a = np.sin(alpha)

x = df_code["X_sat"].to_numpy(dtype=float)
y = df_code["Y_sat"].to_numpy(dtype=float)
z = df_code["Z_sat"].to_numpy(dtype=float)

x_rot = cos_a * x - sin_a * y
y_rot = sin_a * x + cos_a * y
z_rot = z

df_code["X_sat_sagnac"] = x_rot
df_code["Y_sat_sagnac"] = y_rot
df_code["Z_sat_sagnac"] = z_rot

block_dt_r, epoch_unique = build_receiver_clock_block(df_code.index)

P_app = np.array([0.0, 0.0, 0.0])
dP_est = np.array([100.0, 100.0, 100.0])
n_iter = 0

while np.linalg.norm(dP_est[0:3]) > POSITION_CONVERGENCE_THRESHOLD_M:
    distances = np.sqrt(
        (df_code["X_sat_sagnac"].values - P_app[0]) ** 2 +
        (df_code["Y_sat_sagnac"].values - P_app[1]) ** 2 +
        (df_code["Z_sat_sagnac"].values - P_app[2]) ** 2
    )

    B = (
        df_code["code_m"].values
        - distances
        + conv.SPEED_OF_LIGHT * (df_code["dte_sat"].values + df_code["dRelat"].values)
    )

    dX = (P_app[0] - df_code["X_sat_sagnac"].values) / distances
    dY = (P_app[1] - df_code["Y_sat_sagnac"].values) / distances
    dZ = (P_app[2] - df_code["Z_sat_sagnac"].values) / distances
    A = np.column_stack((dX, dY, dZ, block_dt_r))

    dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
    P_est = P_app + dP_est[0:3]

    n_iter += 1
    print(f"Iteration {n_iter}: X={P_est[0]:.3f}, Y={P_est[1]:.3f}, Z={P_est[2]:.3f}")
    P_app = P_est

gnss_edu.plot_residual_analysis(
    A,
    B,
    dP_est,
    figure_title="M3 - Satellite and receiver clocks + Sagnac",
    save_path=WORK_DIR / "04_M3_satellite_receiver_clock_sagnac.png",
    P_est=P_est,
    P_rnx_header=approx_receiver_xyz,
)

store_solution(
    key="M3_sagnac",
    description="Code positioning with satellite clock, relativistic correction, receiver clock, and Sagnac correction",
    P_est=P_est,
    P_ref=approx_receiver_xyz,
    residual_vector=B - A @ dP_est,
    n_iterations=n_iter,
    active_effects=[
        "satellite positions at transmission time",
        "satellite clock correction",
        "relativistic correction",
        "receiver clock estimation",
        "Earth-rotation (Sagnac) correction",
    ],
)

del P_app, dP_est, P_est, n_iter, distances, dX, dY, dZ, A, B
del block_dt_r, epoch_unique
del tau_s, dtheta, alpha, cos_a, sin_a, x, y, z, x_rot, y_rot, z_rot
gc.collect()


# %%
###############################################################################
# Compare the successive positioning models
###############################################################################

df_summary = summarize_solutions(solutions)

print("Model comparison summary")
print("------------------------")
print(df_summary)


# %%
###############################################################################
# Interpretation of the remaining discrepancy
#
# Educational message
# -------------------
# After applying precise satellite positions, satellite clock correction,
# relativistic correction, receiver clock estimation, and Sagnac correction,
# the remaining discrepancy with the RINEX-header position is usually dominated
# by propagation effects, especially in the vertical component.
#
# This naturally motivates the next pedagogical extensions:
#   - ionosphere-free combinations,
#   - tropospheric modeling,
#   - antenna phase-center effects.
###############################################################################

print("Interpretation")
print("--------------")
print(
    "The remaining discrepancy with the RINEX-header position is expected to "
    "be dominated mainly by propagation effects, especially in the vertical "
    "component. This motivates the next model enrichments: ionosphere, "
    "troposphere, and antenna modeling."
)


# %%
###############################################################################
# Step 02 conclusion
###############################################################################

print("Step 02 completed.")
print("Main workspace objects to keep in mind:")
print("  - df_code    : the main MultiIndex working DataFrame")
print("  - solutions  : the successive positioning solutions")
print("  - df_summary : the final model comparison table")






# %%

# %%
###############################################################################
# Export a compact LaTeX table with Model, short Description, and ENU components
# ENU coordinates are rounded to the nearest millimeter.
###############################################################################

import numpy as np
import pandas as pd


def solutions_to_latex_table_enu(
    solutions: dict,
    caption: str = "Progressive evolution of the estimated position in the local ENU frame.",
    label: str = "tab:step02_enu_progression",
    filepath: str | None = None,
) -> tuple[pd.DataFrame, str]:
    """
    Convert the pedagogical `solutions` dictionary into a compact LaTeX table
    showing the model name, a short description, and the ENU components of the
    successive positioning solutions.
    """

    model_name_map = {
        "M0_naive_code": "M0",
        "M1_satellite_clock": "M1",
        "M2_receiver_clock": "M2",
        "M3_sagnac": "M3",
    }

    model_description_map = {
    "M0_naive_code": "Naive code model (Geometry only)",
    "M1_satellite_clock": "+ sat. clock + relativity",
    "M2_receiver_clock": "+ receiver clock",
    "M3_sagnac": "+ Earth rotation (Sagnac effect)",
}

    rows = []

    for model_name, result in solutions.items():
        enu = result.get("enu_error_m", [np.nan, np.nan, np.nan])

        rows.append({
            "Model": model_name_map.get(model_name, model_name),
            "Description": model_description_map.get(
                model_name,
                result.get("description", "")
            ),
            "E [m]": float(enu[0]),
            "N [m]": float(enu[1]),
            "U [m]": float(enu[2]),
        })

    df_summary = pd.DataFrame(rows)

    desired_order = ["M0", "M1", "M2", "M3"]
    df_summary["Model"] = pd.Categorical(
        df_summary["Model"],
        categories=desired_order,
        ordered=True,
    )
    df_summary = df_summary.sort_values("Model").reset_index(drop=True)

    # Round to the nearest millimeter
    df_summary[["E [m]", "N [m]", "U [m]"]] = (
        df_summary[["E [m]", "N [m]", "U [m]"]].round(3)
    )

    latex_table = df_summary.to_latex(
        index=False,
        escape=True,
        caption=caption,
        label=label,
        column_format="llrrr",
        float_format="%.3f",
    )

    if filepath is not None:
        with open(filepath, "w", encoding="utf-8") as f:
            f.write(latex_table)

    return df_summary, latex_table


# Example of use
df_enu, latex_enu = solutions_to_latex_table_enu(
    solutions,
    caption=(
        "Progressive evolution of the estimated position in the local ENU frame "
        "as the code-based GNSS model is enriched."
    ),
    label="tab:step02_enu_progression",
    filepath=str(WORK_DIR / "step02_enu_progression.tex"),
)

print(df_enu)
print()
print(latex_enu)