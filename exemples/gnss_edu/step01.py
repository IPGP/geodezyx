#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step 01 - From RINEX observations to structured GNSS data tables

Educational scope
-----------------
This first step introduces the data structures that will be reused in the
following GNSS teaching scripts.

Learning goals
--------------
By the end of this script, students should be able to:
1. download sample GNSS RINEX observation files with geodezyx,
2. read a RINEX observation file into pandas DataFrames,
3. compare a flat table and a MultiIndex table,
4. extract observations for one satellite and one epoch or time window,
5. visualize one observable over time,
6. identify incomplete observations before moving to GNSS modeling.

Pedagogical position in the curriculum
--------------------------------------
This script corresponds to the "from data to model" stage:
students manipulate GNSS observations with pandas before introducing
least-squares models in Step 02.

@author: Samuel Nahmani (1,2)
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
"""

# %%
###############################################################################
# Reference
#
# Sakic, P., Mansur, G., Chaiyaporn, K., & Ballu, V. (2019).
# The geodeZYX toolbox: a versatile Python 3 toolbox for geodetic-oriented
# purposes (v4.0). GFZ Data Services.
# https://doi.org/10.5880/GFZ.1.1.2019.002
###############################################################################


# %%
###############################################################################
# Imports
###############################################################################

import datetime as dt
import os
from itertools import count
from pathlib import Path

from _io_guard import ensure_matplotlib_cache
from _io_guard import resolve_rinex_files

ensure_matplotlib_cache()
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from geodezyx import files_rw

# Robust import for the pedagogical plotting helper
try:
    from geodezyx.gnss_edu import plot_gnss_timeseries_by_prn
except Exception:
    try:
        from geodezyx import gnss_edu
        plot_gnss_timeseries_by_prn = gnss_edu.plot_gnss_timeseries_by_prn
    except Exception as exc:
        raise ImportError(
            "Could not import 'plot_gnss_timeseries_by_prn' from gnss_edu. "
            "Please check where this function is exposed in your local geodezyx branch."
        ) from exc


# %%
###############################################################################
# Configuration
###############################################################################

PROCESSING_DATE = dt.datetime(2019, 6, 25)
WORK_DIR = Path(os.environ["HOME"]).expanduser() / "gnss_edu_data"
FIGURES_DIR = WORK_DIR / "figures"
SHOW_FIGURES = False

# Two permanent GNSS stations about 10 km apart in the Paris region
STATION_DICT = {"rgp": ["SMNE", "MLVL"]}

# Examples used throughout the script
EXAMPLE_PRN = "G05"
EXAMPLE_OBSERVABLE = "L1"
EXAMPLE_START_TIME = PROCESSING_DATE + dt.timedelta(hours=2, minutes=40)
EXAMPLE_END_TIME = PROCESSING_DATE + dt.timedelta(hours=2, minutes=45, seconds=30)

# Time gap threshold used to split continuous GNSS arcs
GAP_THRESHOLD = pd.Timedelta(minutes=30)

print("Configuration loaded.")
print("Processing date :", PROCESSING_DATE)
print("Working directory:", WORK_DIR)
print("Figures directory:", FIGURES_DIR)
print("Show figures     :", SHOW_FIGURES)
print("Example PRN     :", EXAMPLE_PRN)
print("Example observable:", EXAMPLE_OBSERVABLE)


# %%
###############################################################################
# Create the working directory
#
# Educational objective
# ---------------------
# Prepare a clean local folder to store downloaded data and generated outputs.
###############################################################################

WORK_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

print("Working directory is ready:")
print(WORK_DIR)
print("Figures directory is ready:")
print(FIGURES_DIR)

if SHOW_FIGURES:
    plt.ion()
else:
    plt.ioff()
    plt.show = lambda *args, **kwargs: None

FIGURE_COUNTER = count(1)


def finalize_figure(fig, output_name, show=SHOW_FIGURES):
    """Save a figure to disk and optionally display it."""
    figure_index = next(FIGURE_COUNTER)
    output_name_clean = output_name
    if output_name_clean.startswith("step01_"):
        output_name_clean = output_name_clean[len("step01_") :]
    numbered_output_name = f"step01_{figure_index:02d}_{output_name_clean}"
    output_path = FIGURES_DIR / numbered_output_name
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Figure saved: {output_path}")

    if show:
        fig.show()

    plt.close(fig)

# %%
###############################################################################
# Download sample RINEX observation files
#
# Educational objective
# ---------------------
# Obtain real GNSS observation files that will serve as the raw input data.
#
# Why this matters
# ----------------
# A teaching workflow is much more meaningful when students manipulate
# real observation files rather than artificial toy arrays.
###############################################################################

download_output = resolve_rinex_files(
    statdico=STATION_DICT,
    date=PROCESSING_DATE,
    work_dir=WORK_DIR,
)


# %%
###############################################################################
# Identify the BASE and ROVER files
#
# Note
# ----
# geodezyx usually returns tuples of the form:
#   (local_file_path, success_flag)
###############################################################################

def _extract_download_path(entry):
    """Return the local file path from a geodezyx download entry."""
    if isinstance(entry, tuple):
        if len(entry) >= 2 and not entry[1]:
            raise RuntimeError(f"Download reported as failed for entry: {entry}")
        return entry[0]
    return entry

base_rinex_file = _extract_download_path(download_output[0])
rover_rinex_file = _extract_download_path(download_output[1])

print("BASE RINEX file :")
print(base_rinex_file)
print()
print("ROVER RINEX file:")
print(rover_rinex_file)


# %%
###############################################################################
# Read one RINEX file in two different pandas formats
#
# Educational objective
# ---------------------
# Understand that the same GNSS observations can be organized in different
# tabular structures depending on the intended use.
#
# Flat table
# ----------
# Convenient for explicit boolean filtering and row-by-row inspection.
#
# MultiIndex table
# ----------------
# Convenient when GNSS observations are naturally indexed by (epoch, prn).
###############################################################################

df_flat, rinex_header = files_rw.read_rinex_obs(
    base_rinex_file,
    return_header=True,
)

df_index = files_rw.read_rinex_obs(
    base_rinex_file,
    set_index=["epoch", "prn"],
)

print("Flat DataFrame")
print("--------------")
print("Shape      :", df_flat.shape)
print("Index type :", type(df_flat.index).__name__)
print("Columns    :", list(df_flat.columns[:12]))

print("\nMultiIndex DataFrame")
print("--------------------")
print("Shape      :", df_index.shape)
print("Index type :", type(df_index.index).__name__)
print("Index names:", df_index.index.names)
print("Columns    :", list(df_index.columns[:12]))


# %%
###############################################################################
# Optional inspection in Spyder
#
# Open the Variable Explorer and compare:
#   - df_flat
#   - df_index
#
# Question for students
# ---------------------
# What are the practical differences between these two DataFrames?
###############################################################################

df_flat
# df_index


# %%
###############################################################################
# Extract the approximate receiver position from the RINEX header
#
# Educational objective
# ---------------------
# Show that useful metadata are stored in the RINEX header and can already
# provide a first reference position.
###############################################################################

approx_position_lines = [
    line for line in rinex_header if "APPROX POSITION XYZ" in line
]

if len(approx_position_lines) > 0:
    approx_position_xyz = np.array(
        approx_position_lines[0].split()[:3],
        dtype=float,
    )
else:
    approx_position_xyz = np.array([0.0, 0.0, 0.0])

print("Approximate receiver XYZ position from RINEX header [m]:")
print(approx_position_xyz)


# %%
###############################################################################
# Flat-table selection: one PRN at one epoch
#
# Educational objective
# ---------------------
# Learn explicit selection using boolean masks.
#
# Why this matters
# ----------------
# This is the simplest way to show students how GNSS observations are stored
# before introducing more advanced indexing strategies.
###############################################################################

example_epoch = pd.Timestamp(EXAMPLE_START_TIME)

bool_prn = df_flat["prn"] == EXAMPLE_PRN
bool_epoch = df_flat["epoch"] == example_epoch

flat_single_epoch = df_flat.loc[
    bool_prn & bool_epoch,
    EXAMPLE_OBSERVABLE,
]

print("Flat-table selection")
print("--------------------")
print(f"Satellite  : {EXAMPLE_PRN}")
print(f"Epoch      : {example_epoch}")
print(f"Observable : {EXAMPLE_OBSERVABLE}")
print()
print(flat_single_epoch)


# %%
###############################################################################
# Flat-table selection: one PRN over a time window
#
# Educational objective
# ---------------------
# Extend the previous selection from a single epoch to a small time interval.
###############################################################################

bool_prn = df_flat["prn"] == EXAMPLE_PRN
bool_time_window = (
    (df_flat["epoch"] >= pd.Timestamp(EXAMPLE_START_TIME))
    & (df_flat["epoch"] <= pd.Timestamp(EXAMPLE_END_TIME))
)

flat_time_window = df_flat.loc[
    bool_prn & bool_time_window,
    [EXAMPLE_OBSERVABLE, "epoch", "prn"],
]

print("Flat-table time-window selection")
print("--------------------------------")
print(f"Satellite  : {EXAMPLE_PRN}")
print(f"Start time : {EXAMPLE_START_TIME}")
print(f"End time   : {EXAMPLE_END_TIME}")
print()
print(flat_time_window)


# %%
###############################################################################
# MultiIndex selection: one PRN at one epoch
#
# Educational objective
# ---------------------
# Show how hierarchical indexing makes GNSS queries more natural when
# observations are identified by (epoch, prn).
###############################################################################

index_single_epoch = df_index.loc[
    (example_epoch, EXAMPLE_PRN),
    EXAMPLE_OBSERVABLE,
]

print("MultiIndex selection")
print("--------------------")
print(f"Satellite  : {EXAMPLE_PRN}")
print(f"Epoch      : {example_epoch}")
print(f"Observable : {EXAMPLE_OBSERVABLE}")
print()
print(index_single_epoch)


# %%
###############################################################################
# MultiIndex selection: one PRN over a time window
#
# Educational objective
# ---------------------
# Use slice-based selection in a structure that is closer to GNSS logic.
###############################################################################

index_time_window = df_index.loc[
    (
        slice(pd.Timestamp(EXAMPLE_START_TIME), pd.Timestamp(EXAMPLE_END_TIME)),
        EXAMPLE_PRN,
    ),
    [EXAMPLE_OBSERVABLE],
]

print("MultiIndex time-window selection")
print("--------------------------------")
print(index_time_window)


# %%
###############################################################################
# Add an explicit row identifier
#
# Educational objective
# ---------------------
# Make the future link between data rows and least-squares equations visible.
#
# Why this matters
# ----------------
# In Step 02, each selected observation will progressively contribute to a row
# of the design matrix. Adding a row identifier now helps students see this
# connection explicitly.
###############################################################################

df_index = df_index.copy()
df_index["row_id"] = np.arange(len(df_index))

row_demo = df_index.loc[
    (
        slice(pd.Timestamp(EXAMPLE_START_TIME), pd.Timestamp(EXAMPLE_END_TIME)),
        EXAMPLE_PRN,
    ),
    [EXAMPLE_OBSERVABLE, "row_id"],
]

print("Observation rows with explicit row identifiers")
print("----------------------------------------------")
print(row_demo)


# %%
###############################################################################
# Visualize one observable by satellite PRN
#
# Educational objective
# ---------------------
# Inspect the temporal evolution of one observable and identify tracking arcs.
#
# Why this matters
# ----------------
# GNSS observations are not always continuous in time for each satellite.
# Gap-aware plotting helps students see that the observable is naturally split
# into arcs.
###############################################################################

my_obs_to_extract = EXAMPLE_OBSERVABLE

fig, ax = plot_gnss_timeseries_by_prn(
    df=df_index,
    y=my_obs_to_extract,
    gap=GAP_THRESHOLD,
    label_arcs=True,
    arc_label_fontsize=8,
    show_legend=False,
    title=f"Time series of {my_obs_to_extract} by satellite PRN",
    xlabel="Time (epoch)",
    ylabel=my_obs_to_extract,
    figsize=(10, 6),
)

ax.grid(True, alpha=0.3)
finalize_figure(fig, "step01_timeseries_by_prn_arcs.png")


# %%
###############################################################################
# Optional variant: the same plot with a legend
#
# Note
# ----
# For dense GNSS plots, direct arc labels are often more readable than
# a large legend. This cell is kept as an optional comparison.
###############################################################################

fig, ax = plot_gnss_timeseries_by_prn(
    df=df_index,
    y=my_obs_to_extract,
    gap=GAP_THRESHOLD,
    label_arcs=False,
    show_legend=True,
    legend_outside=True,
    legend_ncol=2,
    title=f"Time series of {my_obs_to_extract} by satellite PRN",
    xlabel="Time (epoch)",
    ylabel=my_obs_to_extract,
    figsize=(10, 6),
)

ax.grid(True, alpha=0.3)
finalize_figure(fig, "step01_timeseries_by_prn_legend.png")


# %%
###############################################################################
# Clean the observation table before later modeling
#
# Educational objective
# ---------------------
# Identify missing values robustly and keep only rows that are complete for
# selected observables.
#
# Why this matters
# ----------------
# In GNSS tables, missing values are not always stored as true NaN values.
# They may appear as blank strings or textual placeholders, which must be
# normalized before column and row cleaning.
###############################################################################

required_observables = ["L1", "L2"]

df_clean = df_index.copy()

# Normalize common textual encodings of missing values
df_clean = df_clean.replace(
    to_replace=[r"^\s*$", "nan", "NaN", "NA", "N/A", "null", "None"],
    value=np.nan,
    regex=True,
)

protected_cols = {"sys", "prn", "prni", "row_id"}

# -------------------------------------------------------------------------
# First pass: remove columns already fully empty
# -------------------------------------------------------------------------
for col in df_clean.columns:
    if col not in protected_cols:
        df_clean[col] = pd.to_numeric(df_clean[col], errors="coerce")

empty_cols_before = df_clean.columns[~df_clean.notna().any(axis=0)].tolist()

print("Columns removed before row cleaning:")
print(empty_cols_before if empty_cols_before else "None")

df_clean = df_clean.drop(columns=empty_cols_before)

# -------------------------------------------------------------------------
# Row cleaning based on required observables
# -------------------------------------------------------------------------
required_observables_available = [
    obs for obs in required_observables if obs in df_clean.columns
]

if not required_observables_available:
    print("\nNo required observable is available. No row-based cleaning performed.")
    df_removed = df_clean.iloc[0:0].copy()
else:
    rows_with_missing = df_clean[required_observables_available].isna().any(axis=1)
    df_removed = df_clean.loc[rows_with_missing].copy()
    df_clean = df_clean.dropna(subset=required_observables_available)

# -------------------------------------------------------------------------
# Second pass: remove columns that became empty after row cleaning
# -------------------------------------------------------------------------
empty_cols_after = df_clean.columns[~df_clean.notna().any(axis=0)].tolist()

print("\nColumns removed after row cleaning:")
print(empty_cols_after if empty_cols_after else "None")

df_clean = df_clean.drop(columns=empty_cols_after)

# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------
print("\nCleaning summary")
print("----------------")
print("Original shape            :", df_index.shape)
print("Final clean shape         :", df_clean.shape)
print("Removed rows              :", len(df_removed))
print("Required observables used :", required_observables_available)

if len(df_removed) > 0:
    print("\nFirst removed rows:")
    print(df_removed.head())

# %%
###############################################################################
# Suggested student questions
#
# 1. Why is the MultiIndex structure convenient for GNSS observations?
# 2. Why can a row identifier be useful before introducing least squares?
# 3. Why should we avoid plotting across large temporal gaps?
# 4. Why do dual-frequency combinations require complete observations?
###############################################################################


# %%
###############################################################################
# Step 01 conclusion
#
# At this stage, students should be able to:
#   - read a RINEX observation file with geodezyx,
#   - compare flat and MultiIndex pandas tables,
#   - extract one satellite at one epoch or over a time window,
#   - visualize one observable by PRN,
#   - identify missing observations before modeling.
#
# Bridge to Step 02
# -----------------
# In Step 02, this structured observation table will be reused to build a
# simple GNSS positioning model. The key pedagogical idea is that each
# observation row will progressively contribute to a least-squares system.
###############################################################################

print("Step 01 completed.")
print("The observation table is now ready for the data-to-model transition in Step 02.")
