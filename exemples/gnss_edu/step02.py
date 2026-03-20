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

from _io_guard import ensure_matplotlib_cache
from _io_guard import extract_download_path
from _io_guard import resolve_rinex_files
from _io_guard import resolve_sp3_product

ensure_matplotlib_cache()
import matplotlib.pyplot as plt
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
FIGURES_DIR = WORK_DIR / "figures"
SHOW_FIGURES = False

# Two permanent GNSS stations about 10 km apart in the Paris region
BASE_STATION = "SMNE"
ROVER_STATION = "MLVL"
STATION_DICT = {"rgp": [BASE_STATION, ROVER_STATION]}

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
print("Figures directory          :", FIGURES_DIR)
print("Show figures               :", SHOW_FIGURES)
print("Satellite system           :", SATELLITE_SYSTEM)
print("Preferred code observable  :", CODE_OBSERVABLE)

WORK_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

if SHOW_FIGURES:
    plt.ion()
else:
    plt.ioff()
    plt.show = lambda *args, **kwargs: None


# %%
###############################################################################
# Pedagogical roadmap
###############################################################################

print(
    "Step 02 roadmap\n"
    "----------------\n"
    "M0: basic code-based positioning\n"
    "M1: add satellite clock and relativistic corrections\n"
    "M2: estimate one receiver clock offset per epoch\n"
    "M3: add Earth rotation (Sagnac) correction\n"
    "M4: switch to the ionosphere-free code combination\n"
    "M5: add simple tropospheric modeling\n"
    "M6: apply carrier smoothing to the ionosphere-free code\n"
    "Next steps: antenna phase-center effects"
)

# This dictionary stores the successive positioning solutions.
# It is the main pedagogical backbone of Step 02.
solutions = {}


# %%
###############################################################################
# Helper functions
###############################################################################


def build_station_file_map(download_entries) -> dict[str, str]:
    """Map 4-character station codes to downloaded local RINEX paths."""
    station_map = {}

    for entry in download_entries:
        path = extract_download_path(entry)
        station_code = Path(path).name[:4].upper()
        station_map[station_code] = path

    return station_map


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


def finalize_figure(fig, output_name, show=SHOW_FIGURES):
    """Save a figure to disk and optionally display it."""
    output_path = FIGURES_DIR / output_name
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Figure saved: {output_path}")

    if show:
        fig.show()

    plt.close(fig)


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

download_output = resolve_rinex_files(
    statdico=STATION_DICT,
    date=PROCESSING_DATE,
    work_dir=WORK_DIR,
)


# %%
###############################################################################
# Resolve the downloaded files for the configured BASE and ROVER stations
###############################################################################

station_file_map = build_station_file_map(download_output)

missing_stations = [
    station for station in [BASE_STATION, ROVER_STATION]
    if station not in station_file_map
]
if missing_stations:
    raise RuntimeError(
        "Missing downloaded RINEX files for requested stations: "
        f"{missing_stations}"
    )

base_rinex_file = station_file_map[BASE_STATION]
rover_rinex_file = station_file_map[ROVER_STATION]

print("BASE station     :", BASE_STATION)
print("BASE RINEX file  :", base_rinex_file)
print()
print("ROVER station    :", ROVER_STATION)
print("ROVER RINEX file :", rover_rinex_file)


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

sp3_path = resolve_sp3_product(
    work_dir=WORK_DIR,
    processing_date=PROCESSING_DATE,
)

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

# Reattach residuals to the observations used in the adjustment.
# This makes the residual analysis physically meaningful, since each residual
# remains linked to its observation context (epoch, PRN, elevation, etc.).
df_residuals = gnss_edu.build_residual_dataframe(
    df_obs_used=df_code,
    A=A,
    B=B,
    dP_est=dP_est,
)

# Choose the residual plot shown in the top-left panel:
# - "timeseries"
# - "by_prn"
residual_plot_mode = "by_prn"
gnss_edu.plot_residual_analysis(
    df_residuals=df_residuals,
    figure_title="M0 - Naive code-based positioning",
    save_path=WORK_DIR / f"01_M0_naive_code_{residual_plot_mode}.png",
    P_est=P_est[:3],
    P_rnx_header=approx_receiver_xyz,
    top_left_plot=residual_plot_mode,
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


# Reattach residuals to the observations used in the adjustment.
# This makes the residual analysis physically meaningful, since each residual
# remains linked to its observation context (epoch, PRN, elevation, etc.).
df_residuals = gnss_edu.build_residual_dataframe(
    df_obs_used=df_code,
    A=A,
    B=B,
    dP_est=dP_est,
)

# Choose the residual plot shown in the top-left panel:
# - "timeseries"
# - "by_prn"
residual_plot_mode = "by_prn"
gnss_edu.plot_residual_analysis(
    df_residuals=df_residuals,
    figure_title="M1 - Satellite clock and relativistic corrections",
    save_path=WORK_DIR / f"02_M1_satellite_clock_{residual_plot_mode}.png",
    P_est=P_est[:3],
    P_rnx_header=approx_receiver_xyz,
    top_left_plot=residual_plot_mode,
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

# Reattach residuals to the observations used in the adjustment.
# This makes the residual analysis physically meaningful, since each residual
# remains linked to its observation context (epoch, PRN, elevation, etc.).
df_residuals = gnss_edu.build_residual_dataframe(
    df_obs_used=df_code,
    A=A,
    B=B,
    dP_est=dP_est,
)

# Choose the residual plot shown in the top-left panel:
# - "timeseries"
# - "by_prn"
residual_plot_mode = "by_prn"
gnss_edu.plot_residual_analysis(
    df_residuals=df_residuals,
    figure_title="M2 - Satellite and receiver clock corrections",
    save_path=WORK_DIR / f"03_M2_satellite_receiver_clock_{residual_plot_mode}.png",
    P_est=P_est[:3],
    P_rnx_header=approx_receiver_xyz,
    top_left_plot=residual_plot_mode,
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

# Reattach residuals to the observations used in the adjustment.
# This makes the residual analysis physically meaningful, since each residual
# remains linked to its observation context (epoch, PRN, elevation, etc.).
df_residuals = gnss_edu.build_residual_dataframe(
    df_obs_used=df_code,
    A=A,
    B=B,
    dP_est=dP_est,
)

# Choose the residual plot shown in the top-left panel:
# - "timeseries"
# - "by_prn"
residual_plot_mode = "by_prn"
gnss_edu.plot_residual_analysis(
    df_residuals=df_residuals,
    figure_title="M3 - Satellite and receiver clocks + Sagnac",
    save_path=WORK_DIR / f"04_M3_satellite_receiver_clock_sagnac_{residual_plot_mode}.png",
    P_est=P_est[:3],
    P_rnx_header=approx_receiver_xyz,
    top_left_plot=residual_plot_mode,
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
# Model M4 - Ionosphere-free code combination
#
# Educational objective
# ---------------------
# Build a dual-frequency ionosphere-free code observable and use it in the
# same positioning framework as in M3.
#
# Physical idea
# -------------
# The first-order ionospheric delay depends on the inverse square of the
# carrier frequency. By combining two code observables measured at different
# frequencies, we can remove this first-order effect.
#
# Pedagogical note
# ----------------
# In this step, we keep the satellite-state columns already computed in df_code.
# This keeps the workflow simple and lets students focus on the effect of the
# new observable itself.
###############################################################################

print("***** M4 - Ionosphere-free code combination *****")

required_if_obs = ["C1", "P2"]
missing_if_obs = [obs for obs in required_if_obs if obs not in df_obs.columns]

if missing_if_obs:
    raise RuntimeError(
        "The ionosphere-free step is skipped. "
        f"Missing observables: {', '.join(missing_if_obs)}"
    )


# -------------------------------------------------------------------------
# Build the first-order ionosphere-free code combination
#
# Principle
# ---------
# The first-order ionospheric delay scales as 1 / f^2. By combining two code
# observables measured at different carrier frequencies, we eliminate this
# dominant ionospheric term.
#
# Important consequence
# ---------------------
# The ionosphere-free code is less biased by the ionosphere, but it is noisier
# than each original code observable taken separately.
# -------------------------------------------------------------------------
f1 = conv.L1_CARRIER_FREQUENCY
f2 = conv.L2_CARRIER_FREQUENCY
f1_sq = f1 ** 2
f2_sq = f2 ** 2

df_if = df_obs[required_if_obs].copy()
df_if = df_if.dropna(subset=required_if_obs).copy()

# P_IF = (f1^2 * P1 - f2^2 * P2) / (f1^2 - f2^2)
# Here, C1 is used as the L1 code observable and P2 as the L2 code observable.
df_if["code_if_m"] = (
    f1_sq * df_if["C1"] - f2_sq * df_if["P2"]
) / (f1_sq - f2_sq)



# -------------------------------------------------------------------------
# Reuse the satellite-state columns already available in df_code
# -------------------------------------------------------------------------
sat_cols_needed = [
    "X_sat", "Y_sat", "Z_sat",
    "dte_sat", "dRelat",
    "X_sat_sagnac", "Y_sat_sagnac", "Z_sat_sagnac",
]

missing_sat_cols = [col for col in sat_cols_needed if col not in df_code.columns]
if missing_sat_cols:
    raise RuntimeError(
        f"Missing satellite-state columns in df_code: {missing_sat_cols}"
    )

df_if = df_if.join(df_code[sat_cols_needed], how="inner")
df_if = df_if.dropna(
    subset=["code_if_m"] + sat_cols_needed
).copy()

df_if["row_id"] = np.arange(len(df_if))

print("Ionosphere-free working DataFrame")
print("---------------------------------")
print("Shape   :", df_if.shape)
print("Columns :", df_if.columns.tolist())
print()
print(df_if.head())

# -------------------------------------------------------------------------
# Receiver clock block
# -------------------------------------------------------------------------
block_dt_r, epoch_unique = build_receiver_clock_block(df_if.index)

# -------------------------------------------------------------------------
# Solve the ionosphere-free positioning model
# -------------------------------------------------------------------------
P_app = np.array([0.0, 0.0, 0.0])
dP_est = np.array([100.0, 100.0, 100.0])
n_iter = 0

while np.linalg.norm(dP_est[0:3]) > POSITION_CONVERGENCE_THRESHOLD_M:
    distances = np.sqrt(
        (df_if["X_sat_sagnac"].values - P_app[0]) ** 2 +
        (df_if["Y_sat_sagnac"].values - P_app[1]) ** 2 +
        (df_if["Z_sat_sagnac"].values - P_app[2]) ** 2
    )

    B = (
        df_if["code_if_m"].values
        - distances
        + conv.SPEED_OF_LIGHT * (df_if["dte_sat"].values + df_if["dRelat"].values)
    )

    dX = (P_app[0] - df_if["X_sat_sagnac"].values) / distances
    dY = (P_app[1] - df_if["Y_sat_sagnac"].values) / distances
    dZ = (P_app[2] - df_if["Z_sat_sagnac"].values) / distances
    A = np.column_stack((dX, dY, dZ, block_dt_r))

    dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)

    P_est = P_app.copy()
    P_est[0:3] = P_app[0:3] + dP_est[0:3]

    n_iter += 1
    print(
        f"Iteration {n_iter}: "
        f"X={P_est[0]:.3f}, Y={P_est[1]:.3f}, Z={P_est[2]:.3f}"
    )

    P_app = P_est.copy()

# Reattach residuals to the observations used in the adjustment.
# This makes the residual analysis physically meaningful, since each residual
# remains linked to its observation context (epoch, PRN, elevation, etc.).
df_residuals = gnss_edu.build_residual_dataframe(
    df_obs_used=df_if,
    A=A,
    B=B,
    dP_est=dP_est,
)

# Choose the residual plot shown in the top-left panel:
# - "timeseries"
# - "by_prn"
residual_plot_mode = "by_prn"
gnss_edu.plot_residual_analysis(
    df_residuals=df_residuals,
    figure_title="M4 - Ionosphere-free code combination",
    save_path=WORK_DIR / f"05_M4_iono_free_{residual_plot_mode}.png",
    P_est=P_est[:3],
    P_rnx_header=approx_receiver_xyz,
    top_left_plot=residual_plot_mode,
)


store_solution(
    key="M4_iono_free",
    description="Iono-free code combination",
    P_est=P_est[:3],
    P_ref=approx_receiver_xyz,
    residual_vector=B - A @ dP_est,
    n_iterations=n_iter,
    active_effects=[
        "ionosphere-free code combination (C1, P2)",
        "satellite clock correction",
        "relativistic correction",
        "receiver clock estimation",
        "Earth-rotation (Sagnac) correction",
    ],
)

# -------------------------------------------------------------------------
# Cleanup
# -------------------------------------------------------------------------
del f1, f2, f1_sq, f2_sq
del block_dt_r, epoch_unique
del P_app, dP_est, P_est, n_iter
del distances, dX, dY, dZ, A, B
gc.collect()

# %%
###############################################################################
# Model M5 - Simple tropospheric modeling with piecewise-constant ZTD
#
# Educational objective
# ---------------------
# Introduce a first explicit tropospheric parameterization by estimating one
# Zenith Total Delay (ZTD) parameter per fixed time interval.
#
# Physical idea
# -------------
# In this simplified model, the slant tropospheric delay is written as:
#
#     STD(elevation) = ZTD / sin(elevation)
#
# where:
#   - STD is the slant tropospheric delay,
#   - ZTD is the zenith total delay,
#   - elevation is the satellite elevation angle.
#
# Pedagogical note
# ----------------
# This is a deliberately simple mapping function. The purpose is not yet to
# introduce a sophisticated tropospheric model, but to show how an atmospheric
# parameter can be estimated directly in the least-squares system.
###############################################################################

print("***** M5 - Simple tropospheric modeling *****")

# User-defined duration (in hours) of each [a,b) ZTD interval
TROPO_ZTD_INTERVAL_HOURS = 1.0

# Elevation cutoff
TROPO_ELEVATION_CUTOFF_DEG = 10.0
TROPO_ELEVATION_CUTOFF_RAD = np.radians(TROPO_ELEVATION_CUTOFF_DEG)

if "df_if" not in locals():
    raise RuntimeError(
        "The simple tropospheric step cannot run because df_if is not available. "
        "Run M4_iono_free first."
    )

# -------------------------------------------------------------------------
# Choose the receiver position used to compute azimuth/elevation
# -------------------------------------------------------------------------
if "M4_iono_free" in solutions:
    receiver_xyz_geom_ref = solutions["M4_iono_free"]["receiver_xyz_m"].copy()
    print("Using M4_iono_free receiver position as geometry reference.")
else:
    receiver_xyz_geom_ref = approx_receiver_xyz.copy()
    print("Using approximate RINEX-header position as geometry reference.")

x0, y0, z0 = receiver_xyz_geom_ref

# -------------------------------------------------------------------------
# Compute azimuth and elevation from the Sagnac-corrected satellite positions
# -------------------------------------------------------------------------
sat_xyz = df_if[["X_sat_sagnac", "Y_sat_sagnac", "Z_sat_sagnac"]].to_numpy(dtype=float)

azimuth_rad_list = []
elevation_rad_list = []

for sat_xyz_i in sat_xyz:
    azi_rad_i, ele_rad_i, _ = conv.xyz2azi_ele(
        sat_xyz_i[0], sat_xyz_i[1], sat_xyz_i[2],
        x0, y0, z0,
        outdeg=False,
    )
    azimuth_rad_list.append(azi_rad_i)
    elevation_rad_list.append(ele_rad_i)

df_if["azimuth_rad"] = np.array(azimuth_rad_list, dtype=float)
df_if["elevation_rad"] = np.array(elevation_rad_list, dtype=float)
df_if["azimuth_deg"] = np.degrees(df_if["azimuth_rad"])
df_if["elevation_deg"] = np.degrees(df_if["elevation_rad"])

print("Elevation statistics before cutoff [deg]:")
print(df_if["elevation_deg"].describe())

# -------------------------------------------------------------------------
# Apply the elevation cutoff
# -------------------------------------------------------------------------
rows_before_cutoff = len(df_if)
df_if = df_if[df_if["elevation_rad"] >= TROPO_ELEVATION_CUTOFF_RAD].copy()
rows_after_cutoff = len(df_if)

print()
print("Elevation cutoff applied")
print("------------------------")
print("Cutoff [deg]       :", TROPO_ELEVATION_CUTOFF_DEG)
print("Rows before cutoff :", rows_before_cutoff)
print("Rows after cutoff  :", rows_after_cutoff)
print("Rows removed       :", rows_before_cutoff - rows_after_cutoff)

# -------------------------------------------------------------------------
# Simple tropospheric mapping function: 1 / sin(elevation)
# -------------------------------------------------------------------------
df_if["tropo_map"] = 1.0 / np.sin(df_if["elevation_rad"].to_numpy(dtype=float))

print()
print("Elevation statistics after cutoff [deg]:")
print(df_if["elevation_deg"].describe())
print()
print("Tropospheric mapping function statistics:")
print(df_if["tropo_map"].describe())

# -------------------------------------------------------------------------
# Build [a,b) ZTD intervals covering all remaining observations
# -------------------------------------------------------------------------
epoch_values = df_if.index.get_level_values("epoch")
epoch_start = epoch_values.min()
epoch_end = epoch_values.max()

interval_length = pd.Timedelta(hours=TROPO_ZTD_INTERVAL_HOURS)
interval_id = ((epoch_values - epoch_start) // interval_length).astype(int)

df_if["ztd_interval_id"] = interval_id

interval_unique = np.sort(df_if["ztd_interval_id"].unique())
n_intervals = len(interval_unique)

print()
print("ZTD interval definition")
print("-----------------------")
print("First epoch         :", epoch_start)
print("Last epoch          :", epoch_end)
print("Data span           :", epoch_end - epoch_start)
print("Interval length [h] :", TROPO_ZTD_INTERVAL_HOURS)
print("Unique interval ids :", interval_unique)
print("Number of intervals :", n_intervals)

for interval_k in interval_unique:
    a_k = epoch_start + interval_k * interval_length
    b_k = a_k + interval_length
    print(f"Interval {interval_k}: [{a_k}, {b_k})")

# -------------------------------------------------------------------------
# Build the receiver-clock block and the ZTD block
# -------------------------------------------------------------------------
block_dt_r, epoch_unique = build_receiver_clock_block(df_if.index)

block_ztd = np.zeros((len(df_if), n_intervals))
interval_to_col = {interval_k: j for j, interval_k in enumerate(interval_unique)}

ztd_interval_array = df_if["ztd_interval_id"].to_numpy()
tropo_map_array = df_if["tropo_map"].to_numpy(dtype=float)

for i_obs in range(len(df_if)):
    j = interval_to_col[ztd_interval_array[i_obs]]
    block_ztd[i_obs, j] = tropo_map_array[i_obs]

print()
print("Receiver clock block shape :", block_dt_r.shape)
print("ZTD block shape            :", block_ztd.shape)
print("First rows of block_ztd:")
print(block_ztd[:5, :])

# -------------------------------------------------------------------------
# Solve the ionosphere-free + simple troposphere positioning model
# -------------------------------------------------------------------------
P_app = np.array([0.0, 0.0, 0.0])
dP_est = np.array([100.0, 100.0, 100.0])
n_iter = 0

while np.linalg.norm(dP_est[0:3]) > POSITION_CONVERGENCE_THRESHOLD_M:
    distances = np.sqrt(
        (df_if["X_sat_sagnac"].values - P_app[0]) ** 2 +
        (df_if["Y_sat_sagnac"].values - P_app[1]) ** 2 +
        (df_if["Z_sat_sagnac"].values - P_app[2]) ** 2
    )

    B = (
        df_if["code_if_m"].values
        - distances
        + conv.SPEED_OF_LIGHT * (df_if["dte_sat"].values + df_if["dRelat"].values)
    )

    dX = (P_app[0] - df_if["X_sat_sagnac"].values) / distances
    dY = (P_app[1] - df_if["Y_sat_sagnac"].values) / distances
    dZ = (P_app[2] - df_if["Z_sat_sagnac"].values) / distances

    A = np.column_stack((dX, dY, dZ, block_dt_r, block_ztd))

    dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)

    P_est = P_app.copy()
    P_est[0:3] = P_app[0:3] + dP_est[0:3]

    n_iter += 1
    print(
        f"Iteration {n_iter}: "
        f"X={P_est[0]:.3f}, Y={P_est[1]:.3f}, Z={P_est[2]:.3f}"
    )

    P_app = P_est.copy()
    
# Reattach residuals to the observations used in the adjustment.
# This makes the residual analysis physically meaningful, since each residual
# remains linked to its observation context (epoch, PRN, elevation, etc.).
df_residuals = gnss_edu.build_residual_dataframe(
    df_obs_used=df_if,
    A=A,
    B=B,
    dP_est=dP_est,
)

# Choose the residual plot shown in the top-left panel:
# - "timeseries"
# - "by_prn"
residual_plot_mode = "by_prn"
gnss_edu.plot_residual_analysis(
    df_residuals=df_residuals,
    figure_title=("M5 - Iono-free code positioning with simple tropospheric modeling\n"
                  f"Iono-free C1/P2, cutoff = {TROPO_ELEVATION_CUTOFF_DEG:.1f} deg,\n"
                  f"STD = ZTD/sin(e), piecewise-constant ZTD every {TROPO_ZTD_INTERVAL_HOURS:.1f} h"
    ),
    save_path=WORK_DIR / f"06_M5_tropo_simple_{residual_plot_mode}.png",
    P_est=P_est[:3],
    P_rnx_header=approx_receiver_xyz,
    top_left_plot=residual_plot_mode,
)


store_solution(
    key="M5_tropo_simple",
    description="Iono-free code + simple troposphere",
    P_est=P_est[:3],
    P_ref=approx_receiver_xyz,
    residual_vector=B - A @ dP_est,
    n_iterations=n_iter,
    active_effects=[
        "ionosphere-free code combination (C1, P2)",
        "satellite clock correction",
        "relativistic correction",
        "receiver clock estimation",
        "Earth-rotation (Sagnac) correction",
        f"observations below {TROPO_ELEVATION_CUTOFF_DEG:.1f} deg excluded",
        f"piecewise-constant ZTD every {TROPO_ZTD_INTERVAL_HOURS:.1f} h",
        "simple mapping function 1/sin(elevation)",
    ],
)

# -------------------------------------------------------------------------
# Store the estimated ZTD parameters for later inspection
# -------------------------------------------------------------------------
n_clock_params = block_dt_r.shape[1]
ztd_estimates_m = dP_est[3 + n_clock_params : 3 + n_clock_params + n_intervals]

solutions["M5_tropo_simple"]["elevation_cutoff_deg"] = TROPO_ELEVATION_CUTOFF_DEG
solutions["M5_tropo_simple"]["ztd_interval_hours"] = TROPO_ZTD_INTERVAL_HOURS
solutions["M5_tropo_simple"]["ztd_interval_ids"] = interval_unique.copy()
solutions["M5_tropo_simple"]["ztd_estimates_m"] = ztd_estimates_m.copy()

print()
print("Estimated ZTD parameters [m]:")
for interval_k, ztd_k in zip(interval_unique, ztd_estimates_m):
    a_k = epoch_start + interval_k * interval_length
    b_k = a_k + interval_length
    print(f"Interval {interval_k} [{a_k}, {b_k}) : ZTD = {ztd_k:.4f} m")

# -------------------------------------------------------------------------
# Cleanup
# -------------------------------------------------------------------------
del receiver_xyz_geom_ref, x0, y0, z0
del sat_xyz, azimuth_rad_list, elevation_rad_list
del rows_before_cutoff, rows_after_cutoff
del epoch_values, epoch_start, epoch_end, interval_length
del interval_id, interval_unique, n_intervals
del block_dt_r, epoch_unique, block_ztd
del interval_to_col, ztd_interval_array, tropo_map_array
del P_app, dP_est, P_est, n_iter
del distances, dX, dY, dZ, A, B
del n_clock_params, ztd_estimates_m
gc.collect()

# %%# %%
###############################################################################
# Plot the estimated ZTD parameters as a piecewise-constant function
###############################################################################

if "M5_tropo_simple" not in solutions:
    print("No M5_tropo_simple solution available.")
else:
    ztd_interval_hours = solutions["M5_tropo_simple"]["ztd_interval_hours"]
    ztd_interval_ids = solutions["M5_tropo_simple"]["ztd_interval_ids"]
    ztd_estimates_m = solutions["M5_tropo_simple"]["ztd_estimates_m"]

    epoch_values = df_if.index.get_level_values("epoch")
    epoch_start = epoch_values.min()
    interval_length = pd.Timedelta(hours=ztd_interval_hours)

    interval_starts = [
        epoch_start + int(k) * interval_length
        for k in ztd_interval_ids
    ]
    interval_ends = [t + interval_length for t in interval_starts]

    # Build step coordinates
    x_step = []
    y_step = []

    for t0, t1, ztd in zip(interval_starts, interval_ends, ztd_estimates_m):
        x_step.extend([t0, t1])
        y_step.extend([ztd, ztd])

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(x_step, y_step)
    ax.set_title("Estimated ZTD parameters (piecewise-constant)")
    ax.set_xlabel("Time")
    ax.set_ylabel("ZTD [m]")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    finalize_figure(fig, "step02_ztd_piecewise_constant.png")


# %%
###############################################################################
# Model M6 - Carrier-smoothed ionosphere-free code
#
# Educational objective
# ---------------------
# Keep the same physical model as in M5, but improve the code observable by
# smoothing the ionosphere-free code with the ionosphere-free carrier phase.
#
# Principle
# ---------
# The code pseudorange provides an absolute range-like observable, but it is
# relatively noisy.
#
# The carrier phase is much less noisy and describes very well how the range
# evolves from one epoch to the next, as long as it remains continuous
# (no cycle slip).
#
# The idea is therefore:
#    - keep the absolute information from the code,
#    - use the carrier phase to reduce the code noise.
#
# Hatch smoothing
# ---------------
# Let:
#    P_tilde(k) = smoothed code at epoch k
#    P(k)       = raw code at epoch k
#    L(k)       = carrier phase in meters at epoch k
#
# The Hatch recursive formula is:
#
#     P_tilde(k) = (1/n) * P(k)
#                + ((n - 1)/n) * [P_tilde(k-1) + L(k) - L(k-1)]
#
# Interpretation:
#    - the phase increment L(k) - L(k-1) provides a precise short-term
#      evolution of the range,
#    - the code P(k) keeps the observable tied to the absolute range level,
#    - the recursion progressively reduces the code noise.
#
# Why ionosphere-free?
# --------------------
# Since M5 uses the ionosphere-free code observable, we must smooth it with
# the corresponding ionosphere-free carrier phase:
#
#     P_IF = ionosphere-free code
#     L_IF = ionosphere-free carrier phase
#
# This keeps the smoothing physically consistent.
#
# What changes compared with M5?
# ------------------------------
# Only the observable changes:
#    M5 : raw ionosphere-free code
#    M6 : carrier-smoothed ionosphere-free code
#
# The positioning model itself remains the same.
#
# Important note
# --------------
# This is still a code-based positioning approach.
# The carrier phase is only used here to improve the code observable.
#
# Quick summary
# -------------
# Raw IF code      : noisy absolute observable
# IF carrier phase : precise relative evolution
# Smoothed IF code : absolute observable with reduced noise
###############################################################################

print("***** M6 - Carrier-smoothed ionosphere-free code *****")

M6_SMOOTHING_METHOD = "forward_backward"   # "forward" or "forward_backward"
M6_SMOOTHING_WINDOW = 100
M6_SAMPLING_SECONDS = 30.0
M6_MAX_GAP_FACTOR = 1.5
M6_SLIP_THRESHOLD_M = 10.0

required_phase_obs = ["L1", "L2"]
missing_phase_obs = [obs for obs in required_phase_obs if obs not in df_obs.columns]

if missing_phase_obs:
    print("The carrier-smoothed ionosphere-free step is skipped.")
    print("Missing observables:", missing_phase_obs)

else:
    # -------------------------------------------------------------------------
    # Select smoothing strategy once, outside the main processing flow
    # -------------------------------------------------------------------------
    if M6_SMOOTHING_METHOD == "forward":
        smoothing_callable = gnss_edu.hatch_carrier_smoothing_by_satellite
        smoothed_col_name = "code_if_smooth_m"
        smoothing_description = "forward Hatch smoothing"
        figure_method_label = "forward"

    elif M6_SMOOTHING_METHOD == "forward_backward":
        smoothing_callable = gnss_edu.hatch_carrier_smoothing_by_satellite_fb
        smoothed_col_name = "code_if_m_smooth_fb"
        smoothing_description = "forward-backward Hatch smoothing"
        figure_method_label = "forward_backward"

    else:
        raise ValueError(
            "Unsupported M6 smoothing method. Use 'forward' or "
            "'forward_backward'."
        )

    # -------------------------------------------------------------------------
    # Build ionosphere-free code and carrier phase
    # -------------------------------------------------------------------------
    f1 = conv.L1_CARRIER_FREQUENCY
    f2 = conv.L2_CARRIER_FREQUENCY
    f1_sq = f1 ** 2
    f2_sq = f2 ** 2

    lambda_1 = conv.SPEED_OF_LIGHT / f1
    lambda_2 = conv.SPEED_OF_LIGHT / f2

    df_phase_if = df_obs[["C1", "P2", "L1", "L2"]].copy()
    df_phase_if = df_phase_if.dropna(subset=["C1", "P2", "L1", "L2"]).copy()

    df_phase_if["code_if_m"] = (
        f1_sq * df_phase_if["C1"] - f2_sq * df_phase_if["P2"]
    ) / (f1_sq - f2_sq)

    df_phase_if["phase_if_m"] = (
        f1_sq * (df_phase_if["L1"] * lambda_1)
        - f2_sq * (df_phase_if["L2"] * lambda_2)
    ) / (f1_sq - f2_sq)

    # -------------------------------------------------------------------------
    # Apply smoothing
    # -------------------------------------------------------------------------
    df_smooth = smoothing_callable(
        df=df_phase_if,
        code_col="code_if_m",
        phase_col="phase_if_m",
        window=M6_SMOOTHING_WINDOW,
        sampling_seconds=M6_SAMPLING_SECONDS,
        max_gap_factor=M6_MAX_GAP_FACTOR,
        slip_threshold_m=M6_SLIP_THRESHOLD_M,
    )

    if isinstance(df_smooth, pd.Series):
        df_phase_if["code_if_smooth_m"] = df_smooth
    else:
        df_phase_if = df_phase_if.join(df_smooth)
        df_phase_if["code_if_smooth_m"] = df_phase_if[smoothed_col_name]

    # -------------------------------------------------------------------------
    # Reuse the M5 observation table and attach the smoothed IF code
    # -------------------------------------------------------------------------
    df_m6 = df_if.copy()

    join_cols = ["phase_if_m", "code_if_smooth_m"]
    for col in [
        "code_if_m_smooth_fwd",
        "code_if_m_smooth_bwd",
        "code_if_m_smooth_fb",
    ]:
        if col in df_phase_if.columns:
            join_cols.append(col)

    df_m6 = df_m6.join(df_phase_if[join_cols], how="left")

    rows_before_dropna = len(df_m6)
    df_m6 = df_m6.dropna(subset=["code_if_smooth_m"]).copy()
    rows_after_dropna = len(df_m6)

    print("Carrier-smoothed IF DataFrame")
    print("-----------------------------")
    print("Rows before dropna :", rows_before_dropna)
    print("Rows after dropna  :", rows_after_dropna)
    print("Rows removed       :", rows_before_dropna - rows_after_dropna)
    print()
    print(df_m6[["code_if_m", "code_if_smooth_m", "phase_if_m"]].head())

    # -------------------------------------------------------------------------
    # Rebuild the receiver-clock block and the ZTD block
    # -------------------------------------------------------------------------
    block_dt_r_m6, epoch_unique_m6 = build_receiver_clock_block(df_m6.index)

    interval_unique_m6 = np.sort(df_m6["ztd_interval_id"].unique())
    n_intervals_m6 = len(interval_unique_m6)

    block_ztd_m6 = np.zeros((len(df_m6), n_intervals_m6))
    interval_to_col_m6 = {
        interval_k: j for j, interval_k in enumerate(interval_unique_m6)
    }

    ztd_interval_array_m6 = df_m6["ztd_interval_id"].to_numpy()
    tropo_map_array_m6 = df_m6["tropo_map"].to_numpy(dtype=float)

    for i_obs in range(len(df_m6)):
        j = interval_to_col_m6[ztd_interval_array_m6[i_obs]]
        block_ztd_m6[i_obs, j] = tropo_map_array_m6[i_obs]

    print()
    print("Receiver clock block shape (M6) :", block_dt_r_m6.shape)
    print("ZTD block shape (M6)            :", block_ztd_m6.shape)

    # -------------------------------------------------------------------------
    # Solve the same model as M5, but with carrier-smoothed IF code
    # -------------------------------------------------------------------------
    P_app = np.array([0.0, 0.0, 0.0])
    dP_est = np.array([100.0, 100.0, 100.0])
    n_iter = 0

    while np.linalg.norm(dP_est[0:3]) > POSITION_CONVERGENCE_THRESHOLD_M:
        distances = np.sqrt(
            (df_m6["X_sat_sagnac"].values - P_app[0]) ** 2
            + (df_m6["Y_sat_sagnac"].values - P_app[1]) ** 2
            + (df_m6["Z_sat_sagnac"].values - P_app[2]) ** 2
        )

        B = (
            df_m6["code_if_smooth_m"].values
            - distances
            + conv.SPEED_OF_LIGHT
            * (df_m6["dte_sat"].values + df_m6["dRelat"].values)
        )

        dX = (P_app[0] - df_m6["X_sat_sagnac"].values) / distances
        dY = (P_app[1] - df_m6["Y_sat_sagnac"].values) / distances
        dZ = (P_app[2] - df_m6["Z_sat_sagnac"].values) / distances

        A = np.column_stack((dX, dY, dZ, block_dt_r_m6, block_ztd_m6))

        dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)

        P_est = P_app.copy()
        P_est[0:3] = P_app[0:3] + dP_est[0:3]

        n_iter += 1
        print(
            f"Iteration {n_iter}: "
            f"X={P_est[0]:.3f}, Y={P_est[1]:.3f}, Z={P_est[2]:.3f}"
        )

        P_app = P_est.copy()
        
        
    # Reattach residuals to the observations used in the adjustment.
    # This makes the residual analysis physically meaningful, since each residual
    # remains linked to its observation context (epoch, PRN, elevation, etc.).
    df_residuals = gnss_edu.build_residual_dataframe(
        df_obs_used=df_m6,
        A=A,
        B=B,
        dP_est=dP_est,
    )
    
    # Choose the residual plot shown in the top-left panel:
    # - "timeseries"
    # - "by_prn"
    residual_plot_mode = "by_prn"
    gnss_edu.plot_residual_analysis(
        df_residuals=df_residuals,
        figure_title=(
            "M6 - Carrier-smoothed iono-free code positioning\n"
            f"Method = {figure_method_label}, "
            f"C1/P2 iono-free smoothed by L1/L2, "
            f"cutoff = {TROPO_ELEVATION_CUTOFF_DEG:.1f} deg,\n "
            f"STD = ZTD/sin(e), "
            f"constant ZTD every {TROPO_ZTD_INTERVAL_HOURS:.1f} h"
        ),
        save_path=WORK_DIR / f"07_M6_code_phase_smoothing_{figure_method_label}_{residual_plot_mode}.png",
        P_est=P_est[:3],
        P_rnx_header=approx_receiver_xyz,
        top_left_plot=residual_plot_mode,
    )


    solution_key = f"M6_code_phase_smoothing_{figure_method_label}"

    store_solution(
        key=solution_key,
        description=(
            "Carrier-smoothed iono-free code + simple troposphere "
            f"({smoothing_description})"
        ),
        P_est=P_est[:3],
        P_ref=approx_receiver_xyz,
        residual_vector=B - A @ dP_est,
        n_iterations=n_iter,
        active_effects=[
            "ionosphere-free code combination (C1, P2)",
            "ionosphere-free carrier phase combination (L1, L2)",
            smoothing_description,
            f"smoothing window = {M6_SMOOTHING_WINDOW:d} epochs",
            f"nominal sampling = {M6_SAMPLING_SECONDS:.1f} s",
            f"max gap factor = {M6_MAX_GAP_FACTOR:.1f}",
            (
                f"slip threshold = {M6_SLIP_THRESHOLD_M:.1f} m"
                if M6_SLIP_THRESHOLD_M is not None
                else "no slip-threshold consistency test"
            ),
            "satellite clock correction",
            "relativistic correction",
            "receiver clock estimation",
            "Earth-rotation (Sagnac) correction",
            f"observations below {TROPO_ELEVATION_CUTOFF_DEG:.1f} deg excluded",
            f"piecewise-constant ZTD every {TROPO_ZTD_INTERVAL_HOURS:.1f} h",
            "simple mapping function 1/sin(elevation)",
        ],
    )

    # -------------------------------------------------------------------------
    # Store useful M6 quantities for later inspection
    # -------------------------------------------------------------------------
    n_clock_params_m6 = block_dt_r_m6.shape[1]
    ztd_estimates_m6 = dP_est[
        3 + n_clock_params_m6 : 3 + n_clock_params_m6 + n_intervals_m6
    ]

    solutions[solution_key]["elevation_cutoff_deg"] = TROPO_ELEVATION_CUTOFF_DEG
    solutions[solution_key]["ztd_interval_hours"] = TROPO_ZTD_INTERVAL_HOURS
    solutions[solution_key]["ztd_interval_ids"] = interval_unique_m6.copy()
    solutions[solution_key]["ztd_estimates_m"] = ztd_estimates_m6.copy()
    solutions[solution_key]["smoothing_method"] = M6_SMOOTHING_METHOD
    solutions[solution_key]["smoothing_window_epochs"] = M6_SMOOTHING_WINDOW
    solutions[solution_key]["sampling_seconds"] = M6_SAMPLING_SECONDS
    solutions[solution_key]["max_gap_factor"] = M6_MAX_GAP_FACTOR
    solutions[solution_key]["slip_threshold_m"] = M6_SLIP_THRESHOLD_M
    solutions[solution_key]["smoothed_code_used_m"] = df_m6["code_if_smooth_m"].copy()

    if "code_if_m_smooth_fwd" in df_m6.columns:
        solutions[solution_key]["smoothed_code_forward_m"] = (
            df_m6["code_if_m_smooth_fwd"].copy()
        )

    if "code_if_m_smooth_bwd" in df_m6.columns:
        solutions[solution_key]["smoothed_code_backward_m"] = (
            df_m6["code_if_m_smooth_bwd"].copy()
        )

    if "code_if_m_smooth_fb" in df_m6.columns:
        solutions[solution_key]["smoothed_code_forward_backward_m"] = (
            df_m6["code_if_m_smooth_fb"].copy()
        )

    # -------------------------------------------------------------------------
    # Cleanup
    # -------------------------------------------------------------------------
    del f1, f2, f1_sq, f2_sq, lambda_1, lambda_2
    del block_dt_r_m6, epoch_unique_m6
    del interval_unique_m6, n_intervals_m6, interval_to_col_m6
    del ztd_interval_array_m6, tropo_map_array_m6, block_ztd_m6
    del rows_before_dropna, rows_after_dropna
    del P_app, dP_est, P_est, n_iter
    del distances, dX, dY, dZ, A, B
    del n_clock_params_m6, ztd_estimates_m6
    del smoothing_callable, smoothed_col_name, smoothing_description
    del figure_method_label, solution_key, df_smooth
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
# After introducing satellite clock correction, relativistic correction,
# receiver clock estimation, Earth-rotation correction, ionosphere-free code,
# simple tropospheric modeling, and carrier smoothing, the remaining mismatch
# with the RINEX-header position should be interpreted with caution.
#
# At this stage, the residual discrepancy is no longer driven only by the main
# first-order propagation effects. It can also reflect the limits of the code-
# based model itself, residual atmospheric effects, antenna phase-center
# effects, measurement noise, and the fact that the RINEX-header position is
# only an approximate reference.
###############################################################################

print("Interpretation")
print("--------------")
print(
    "After the successive model enrichments up to M6, the remaining \n"
    "discrepancy with the RINEX-header position is expected to reflect both \n"
    "residual propagation effects and the intrinsic limits of a code-based \n"
    "positioning model. Remaining contributors may include residual \n"
    "tropospheric effects, antenna phase-center effects, measurement noise, \n"
    "and the approximate nature of the RINEX-header coordinates themselves."
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
###############################################################################
# Export a compact LaTeX table with Model, short Description, and ENU components
# ENU coordinates are rounded to the nearest millimeter.
###############################################################################


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
        "M4_iono_free": "M4",
        "M5_tropo_simple": "M5",
        "M6_code_phase_smoothing_forward": "M6",
        "M6_code_phase_smoothing_forward_backward": "M6",
    }

    model_description_map = {
        "M0_naive_code": "Naive code model (Geometry only)",
        "M1_satellite_clock": "+ sat. clock + relativity",
        "M2_receiver_clock": "+ receiver clock",
        "M3_sagnac": "+ Earth rotation (Sagnac effect)",
        "M4_iono_free": "+ iono-free code",
        "M5_tropo_simple": "+ simple troposphere",
        "M6_code_phase_smoothing_forward": "+ carrier-smoothed iono-free code (forward Hatch)",
        "M6_code_phase_smoothing_forward_backward": "+ carrier-smoothed iono-free code (forward-backward Hatch)",
    }

    rows = []

    for key, sol in solutions.items():
        if "enu_error_m" not in sol:
            continue

        enu = np.asarray(sol["enu_error_m"], dtype=float).ravel()
        if enu.size != 3:
            raise ValueError(
                f"Solution '{key}' has an invalid 'enu_error_m'. Expected 3 values."
            )

        model = model_name_map.get(key, key)
        description = model_description_map.get(key, sol.get("description", ""))

        rows.append(
            {
                "Model": model,
                "Description": description,
                "E [m]": enu[0],
                "N [m]": enu[1],
                "U [m]": enu[2],
            }
        )

    df_table = pd.DataFrame(rows)

    if df_table.empty:
        raise ValueError(
            "No valid ENU solutions were found in the input dictionary. "
            "Expected key 'enu_error_m' in each solution."
        )

    desired_order = ["M0", "M1", "M2", "M3", "M4", "M5", "M6"]
    df_table["Model"] = pd.Categorical(
        df_table["Model"],
        categories=desired_order,
        ordered=True,
    )
    df_table = df_table.sort_values("Model").reset_index(drop=True)

    latex_table = df_table.to_latex(
        index=False,
        caption=caption,
        label=label,
        float_format="%.3f",
        column_format="llrrr",
        escape=False,
    )

    if filepath is not None:
        with open(filepath, "w", encoding="utf-8") as f:
            f.write(latex_table)

    return df_table, latex_table


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
