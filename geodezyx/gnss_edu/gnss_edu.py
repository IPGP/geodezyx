#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 16:39:34 2025

@author: snahmani
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import seaborn as sns
import statsmodels.api as sm
from scipy import stats

import plotly.graph_objects as go
import plotly.io as pio

# GeodeZYX Toolbox’s - [Sakic et al., 2019]
from geodezyx import files_rw     # Import the read/write module
from geodezyx import conv         # Import the conversion module


import datetime as dt
from .klobuchar import *


import warnings


# Pour chercher des chaînes de caractère dans des fichiers
import re



def detect_intra_arc_holes(df, sampling=pd.Timedelta(seconds=30),
                          gap=pd.Timedelta(minutes=30)):
    """
    Detect missing observations inside satellite arcs.

    A hole is defined as:
        sampling*2 < dt < gap

    Parameters
    ----------
    df : pandas.DataFrame
        GNSS dataframe indexed by (epoch, prn)
    sampling : Timedelta
        Expected RINEX sampling interval
    gap : Timedelta
        Arc segmentation threshold

    Returns
    -------
    DataFrame with columns:
        epoch_prev : last valid observation
        epoch      : first observation after the hole
        dt         : time difference
        n_missing  : number of missing epochs
    """

    thr_min = 2 * sampling
    thr_max = gap

    rows = []

    for prn in df.index.get_level_values("prn").unique():

        t = df.xs(prn, level="prn").index.to_series().sort_values()
        dt = t.diff()

        mask = (dt > thr_min) & (dt < thr_max)

        if mask.any():

            idx = dt.index[mask]
            dt_sel = dt.loc[idx]
            prev = t.shift(1).loc[idx]

            n_missing = (dt_sel / sampling).round().astype(int) - 1

            rows.append(pd.DataFrame({
                "prn": prn,
                "epoch_prev": prev.values,
                "epoch": idx.values,
                "dt": dt_sel.values,
                "n_missing": n_missing.values
            }))

    if rows:

        holes = (
            pd.concat(rows, ignore_index=True)
            .set_index(["prn", "epoch"])
            .sort_values(["n_missing", "dt"], ascending=False)
        )

    else:

        holes = pd.DataFrame(
            columns=["epoch_prev", "dt", "n_missing"],
            index=pd.MultiIndex.from_arrays([[], []], names=["prn", "epoch"])
        )

    print(f"Intra-arc holes detected: {len(holes)}")

    return holes





def grep_file(pattern, filename):
    """Recherche un motif dans un fichier et retourne les lignes correspondantes."""
    results = []
    with open(filename, "r", encoding="utf-8") as file:
        for line in file:
            if re.search(pattern, line):
                results.append(line.strip())  # Supprime les espaces et les sauts de ligne

    return results  # Retourne les lignes trouvées

def get_approx_position(fichier):
    """
    Recherche la ligne contenant 'APPROX POSITION XYZ' dans le fichier
    et retourne les trois premières valeurs sous forme de numpy.array.
    """
    lignes = grep_file(r"APPROX POSITION XYZ", fichier)
    if lignes:
        # Récupérer la première ligne correspondante et les trois premières valeurs
        valeurs = lignes[0].split()[:3]
        return np.array(valeurs, dtype=float)
    else:
        return None
    
def hatch_carrier_smoothing_by_satellite(
    df: pd.DataFrame,
    code_col: str,
    phase_col: str,
    window: int = 100,
    sampling_seconds: float = 30.0,
    max_gap_factor: float = 1.5,
    slip_threshold_m: float | None = 10.0,
) -> pd.Series:
    """
    Apply Hatch carrier smoothing to a code observable independently for each
    satellite.

    This function implements a simple code-phase smoothing strategy in which
    the code observable provides the absolute range level, while the carrier
    phase provides a much less noisy short-term range evolution.

    The smoothing is performed separately for each satellite track, assuming
    that the input DataFrame uses a MultiIndex containing at least:
        - 'epoch' : observation epoch
        - 'prn'   : satellite identifier

    Parameters
    ----------
    df : pd.DataFrame
        Input observation DataFrame indexed by a MultiIndex containing
        at least ('epoch', 'prn').
    code_col : str
        Name of the code observable column, in meters.
    phase_col : str
        Name of the carrier-phase observable column, in meters.
        Important: this must already be expressed in meters, not in cycles.
    window : int, default=100
        Maximum Hatch smoothing length, in number of epochs.

        The effective smoothing duration is approximately:

            window * sampling_seconds

        With the default values used in this TP:
            100 * 30 s = 3000 s ≈ 50 min

        A larger window gives stronger smoothing but requires longer
        continuity of the phase observable.
    sampling_seconds : float, default=30.0
        Nominal time interval, in seconds, between two consecutive epochs
        for a given satellite.

        This is used to detect abnormal time gaps. For example:
            - 30.0 for a 30-second RINEX
            - 1.0  for a 1-second RINEX
    max_gap_factor : float, default=1.5
        Time-gap tolerance factor used to decide whether the smoothing filter
        must be reset.

        The filter is reset when the time interval between two consecutive
        epochs exceeds:

            max_gap_factor * sampling_seconds

        With the default values:
            1.5 * 30 s = 45 s

        Therefore, the filter is reset if a time gap larger than 45 seconds
        is detected.
    slip_threshold_m : float | None, default=10.0
        Optional crude consistency threshold, in meters, used as a simple
        protection against cycle slips or other discontinuities.

        At each epoch, the function compares:
            - the code increment   : P(k) - P(k-1)
            - the phase increment  : L(k) - L(k-1)

        If:

            abs( [P(k)-P(k-1)] - [L(k)-L(k-1)] ) > slip_threshold_m

        the filter is reset.

        This is only a simple pedagogical safeguard, not a rigorous
        cycle-slip detection method.

        If set to None, this additional test is disabled.

    Returns
    -------
    pd.Series
        Smoothed code observable, in meters, aligned with df.index.

    Notes
    -----
    The recursive Hatch formula is:

        P_tilde(k) = (1/n) * P(k)
                   + ((n - 1)/n) * [P_tilde(k-1) + L(k) - L(k-1)]

    where:
        - P_tilde(k) is the smoothed code,
        - P(k) is the raw code,
        - L(k) is the carrier phase in meters,
        - n is the current smoothing length, increasing progressively up to
          the maximum value 'window'.

    The filter is reset when:
        - the code or phase is missing,
        - the time gap is too large,
        - the optional code-phase increment consistency test fails.

    In this TP, the default values are chosen for a simple and robust
    introductory use case based on 30-second GNSS data.
    """
    # -------------------------------------------------------------------------
    # Input validation
    # -------------------------------------------------------------------------
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")

    if not isinstance(df.index, pd.MultiIndex):
        raise ValueError(
            "The input DataFrame must use a MultiIndex containing at least "
            "'epoch' and 'prn'."
        )

    required_index_names = {"epoch", "prn"}
    if not required_index_names.issubset(df.index.names):
        raise ValueError(
            "The DataFrame MultiIndex must contain at least the levels "
            "'epoch' and 'prn'."
        )

    if code_col not in df.columns:
        raise KeyError(f"Column '{code_col}' was not found in the DataFrame.")

    if phase_col not in df.columns:
        raise KeyError(f"Column '{phase_col}' was not found in the DataFrame.")

    if not isinstance(window, int):
        raise TypeError("window must be an integer.")

    if window < 1:
        raise ValueError("window must be >= 1.")

    if sampling_seconds <= 0:
        raise ValueError("sampling_seconds must be > 0.")

    if max_gap_factor <= 0:
        raise ValueError("max_gap_factor must be > 0.")

    if slip_threshold_m is not None and slip_threshold_m <= 0:
        raise ValueError("slip_threshold_m must be > 0 when provided.")

    # -------------------------------------------------------------------------
    # Prepare output
    # -------------------------------------------------------------------------
    smoothed = pd.Series(
        index=df.index,
        dtype=float,
        name=f"{code_col}_smoothed",
    )

    epoch_level = df.index.names.index("epoch")

    # -------------------------------------------------------------------------
    # Process each satellite independently
    # -------------------------------------------------------------------------
    for _, df_sat in df.groupby(level="prn", sort=False):
        df_sat = df_sat.sort_index(level="epoch")

        prev_epoch = None
        prev_code = None
        prev_phase = None
        prev_smoothed = None
        n_valid = 0

        for idx, row in df_sat.iterrows():
            epoch = idx[epoch_level]
            code = row[code_col]
            phase = row[phase_col]

            reset_filter = False

            # Missing data break the smoothing continuity
            if pd.isna(code) or pd.isna(phase):
                reset_filter = True

            # Large time gaps also break the continuity
            if not reset_filter and prev_epoch is not None:
                dt_seconds = (epoch - prev_epoch).total_seconds()
                if dt_seconds > max_gap_factor * sampling_seconds:
                    reset_filter = True

            # Optional crude code-phase consistency check
            if (
                not reset_filter
                and slip_threshold_m is not None
                and prev_code is not None
                and prev_phase is not None
            ):
                code_increment = code - prev_code
                phase_increment = phase - prev_phase

                if abs(code_increment - phase_increment) > slip_threshold_m:
                    reset_filter = True

            # Initialize or continue the Hatch recursion
            if reset_filter or prev_smoothed is None:
                smoothed.loc[idx] = code
                n_valid = 1
            else:
                n_valid = min(n_valid + 1, window)
                smoothed.loc[idx] = (
                    (1.0 / n_valid) * code
                    + ((n_valid - 1.0) / n_valid)
                    * (prev_smoothed + (phase - prev_phase))
                )

            prev_epoch = epoch
            prev_code = code
            prev_phase = phase
            prev_smoothed = smoothed.loc[idx]

    return smoothed

def hatch_carrier_smoothing_by_satellite_fb(
    df: pd.DataFrame,
    code_col: str,
    phase_col: str,
    window: int = 100,
    sampling_seconds: float = 30.0,
    max_gap_factor: float = 1.5,
    slip_threshold_m: float | None = 10.0,
) -> pd.DataFrame:
    """
    Apply forward-backward Hatch carrier smoothing independently for each
    satellite.

    The function returns three aligned series:
        - forward  Hatch smoothing
        - backward Hatch smoothing
        - forward-backward average

    Parameters
    ----------
    df : pd.DataFrame
        Input observation DataFrame indexed by a MultiIndex containing at
        least ('epoch', 'prn').
    code_col : str
        Name of the code observable column, in meters.
    phase_col : str
        Name of the carrier-phase observable column, in meters.
        Important: this phase must already be expressed in meters.
    window : int, default=100
        Maximum Hatch smoothing length, in number of epochs.
    sampling_seconds : float, default=30.0
        Nominal sampling interval in seconds.
    max_gap_factor : float, default=1.5
        Reset the filter when the time gap exceeds
        max_gap_factor * sampling_seconds.
    slip_threshold_m : float | None, default=10.0
        Optional crude consistency threshold, in meters, applied to the
        difference between code and phase increments. If exceeded, the filter
        is reset. If None, no such test is applied.

    Returns
    -------
    pd.DataFrame
        DataFrame with three columns:
            - f"{code_col}_smooth_fwd"
            - f"{code_col}_smooth_bwd"
            - f"{code_col}_smooth_fb"

    Notes
    -----
    The forward-backward average is a simple way to reduce the edge effects
    of a purely forward causal smoother. It remains pedagogical and easy to
    inspect, but it is not an optimal state-space smoother.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")

    if not isinstance(df.index, pd.MultiIndex):
        raise ValueError(
            "The input DataFrame must use a MultiIndex containing at least "
            "'epoch' and 'prn'."
        )

    required_index_names = {"epoch", "prn"}
    if not required_index_names.issubset(df.index.names):
        raise ValueError(
            "The DataFrame MultiIndex must contain at least the levels "
            "'epoch' and 'prn'."
        )

    if code_col not in df.columns:
        raise KeyError(f"Column '{code_col}' was not found in the DataFrame.")

    if phase_col not in df.columns:
        raise KeyError(f"Column '{phase_col}' was not found in the DataFrame.")

    if not isinstance(window, int):
        raise TypeError("window must be an integer.")

    if window < 1:
        raise ValueError("window must be >= 1.")

    if sampling_seconds <= 0:
        raise ValueError("sampling_seconds must be > 0.")

    if max_gap_factor <= 0:
        raise ValueError("max_gap_factor must be > 0.")

    if slip_threshold_m is not None and slip_threshold_m <= 0:
        raise ValueError("slip_threshold_m must be > 0 when provided.")

    def _hatch_one_direction(
        df_in: pd.DataFrame,
        reverse: bool = False,
    ) -> pd.Series:
        """
        Apply one-direction Hatch smoothing satellite by satellite.

        If reverse=False, smoothing is applied forward in time.
        If reverse=True, smoothing is applied backward in time.
        """
        smoothed = pd.Series(
            index=df_in.index,
            dtype=float,
            name=f"{code_col}_smoothed",
        )

        epoch_level = df_in.index.names.index("epoch")

        for _, df_sat in df_in.groupby(level="prn", sort=False):
            df_sat = df_sat.sort_index(level="epoch", ascending=not reverse)

            prev_epoch = None
            prev_code = None
            prev_phase = None
            prev_smoothed = None
            n_valid = 0

            for idx, row in df_sat.iterrows():
                epoch = idx[epoch_level]
                code = row[code_col]
                phase = row[phase_col]

                reset_filter = False

                if pd.isna(code) or pd.isna(phase):
                    reset_filter = True

                if not reset_filter and prev_epoch is not None:
                    dt_seconds = abs((epoch - prev_epoch).total_seconds())
                    if dt_seconds > max_gap_factor * sampling_seconds:
                        reset_filter = True

                if (
                    not reset_filter
                    and slip_threshold_m is not None
                    and prev_code is not None
                    and prev_phase is not None
                ):
                    code_increment = code - prev_code
                    phase_increment = phase - prev_phase

                    if abs(code_increment - phase_increment) > slip_threshold_m:
                        reset_filter = True

                if reset_filter or prev_smoothed is None:
                    smoothed.loc[idx] = code
                    n_valid = 1
                else:
                    n_valid = min(n_valid + 1, window)
                    smoothed.loc[idx] = (
                        (1.0 / n_valid) * code
                        + ((n_valid - 1.0) / n_valid)
                        * (prev_smoothed + (phase - prev_phase))
                    )

                prev_epoch = epoch
                prev_code = code
                prev_phase = phase
                prev_smoothed = smoothed.loc[idx]

        return smoothed

    smooth_fwd = _hatch_one_direction(df, reverse=False)
    smooth_bwd = _hatch_one_direction(df, reverse=True)
    smooth_fb = 0.5 * (smooth_fwd + smooth_bwd)

    return pd.DataFrame(
        {
            f"{code_col}_smooth_fwd": smooth_fwd,
            f"{code_col}_smooth_bwd": smooth_bwd,
            f"{code_col}_smooth_fb": smooth_fb,
        },
        index=df.index,
    )



def load_and_clean_rinex(path):
    """
    Charge et nettoie un fichier RINEX d'observation.

    Paramètres :
      - path : chemin vers le fichier RINEX.

    Retourne :
      - df : DataFrame nettoyé, indexé par ['epoch', 'prn'],
             et contenant une colonne 'ind_ligne' indiquant le numéro de ligne.
    """
    # Chargement du fichier dans un DataFrame avec l'index ['epoch', 'prn']
    df = files_rw.read_rinex_obs(path, set_index=['epoch', 'prn'])

    # Suppression des colonnes entièrement vides
    df = df.dropna(axis=1, how='all')

    # Suppression des lignes manquant des mesures essentielles
    df = df.dropna(subset=['C1', 'L1', 'L2'])

    # Filtrage pour ne garder que les mesures des satellites GPS ('G')
    df = df[df['sys'].str.contains('G')]

    # On refait un dropna au cas où après filtrage certaines colonnes seraient entièrement vides
    df = df.dropna(axis=1, how='all')

    # Ajout d'une colonne d'indice de ligne basée sur la longueur du DataFrame
    df['ind_ligne'] = range(len(df))

    return df




def Sagnac_rotate_around_z(row):
    """
    Calcule la correction Sagnac pour les coordonnées satellites.
    Retourne un pd.Series avec les nouvelles colonnes X_sat_corr, Y_sat_corr, Z_sat_corr.
    """
    # Calcul de l'angle de rotation en radians
    alpha_rad = row['C1'] / conv.SPEED_OF_LIGHT * gnss_const.Omega_e
    # Construction de la matrice de rotation autour de l'axe Z
    Rz = np.array([[np.cos(alpha_rad), -np.sin(alpha_rad), 0],
                   [np.sin(alpha_rad),  np.cos(alpha_rad), 0],
                   [0, 0, 1]])
    # Vecteur de position original du satellite
    original_vector = np.array([row['X_sat'], row['Y_sat'], row['Z_sat']])
    # Calcul du vecteur corrigé
    rotated_vector = Rz.dot(original_vector)
    # Retourne les nouvelles coordonnées dans un pd.Series
    return pd.Series({
        'X_sat_sagnac': rotated_vector[0],
        'Y_sat_sagnac': rotated_vector[1],
        'Z_sat_sagnac': rotated_vector[2]
    })





def plot_tracking_timeline(
    df: pd.DataFrame,
    sampling: pd.Timedelta = pd.Timedelta(seconds=30),
    snr_col: str | None = None,
    snr_min: float | None = None,
    title: str | None = None,
):
    """
    Timeline heatmap of tracking availability by PRN.

    Black = epoch exists for this PRN (and SNR >= snr_min if provided)
    White = missing epoch (or SNR < snr_min)

    Parameters
    ----------
    df : DataFrame
        Must be indexed by MultiIndex ('epoch','prn').
    sampling : Timedelta
        Expected observation interval (e.g. 30 s).
    snr_col : str or None
        Column name for SNR (e.g., 'S1'). If None, only presence/absence is shown.
    snr_min : float or None
        If provided with snr_col, only epochs with SNR >= snr_min are considered "present".
    """
    if not isinstance(df.index, pd.MultiIndex):
        raise ValueError("df must have MultiIndex ('epoch','prn').")
    if "epoch" not in df.index.names or "prn" not in df.index.names:
        raise ValueError("df index must have levels ('epoch','prn').")

    if snr_col is not None and snr_col not in df.columns:
        raise ValueError(f"snr_col='{snr_col}' not found in df columns.")

    epochs_obs = df.index.get_level_values("epoch")
    prns = np.array(sorted(df.index.get_level_values("prn").unique()))

    t0, t1 = epochs_obs.min(), epochs_obs.max()
    expected_epochs = pd.date_range(start=t0, end=t1, freq=sampling)

    epoch_to_col = {t: i for i, t in enumerate(expected_epochs)}
    avail = np.zeros((len(prns), len(expected_epochs)), dtype=np.uint8)

    for i, prn in enumerate(prns):
        d = df.xs(prn, level="prn").sort_index()

        if snr_col is not None and snr_min is not None:
            d = d.dropna(subset=[snr_col])
            d = d[d[snr_col] >= snr_min]

        cols = [epoch_to_col.get(t) for t in d.index]
        cols = [c for c in cols if c is not None]
        if cols:
            avail[i, cols] = 1

    fig, ax = plt.subplots(figsize=(14, 0.28 * len(prns) + 2))
    cmap = ListedColormap(["white", "black"])
    ax.imshow(avail, aspect="auto", interpolation="nearest", cmap=cmap)

    ax.set_yticks(np.arange(len(prns)))
    ax.set_yticklabels(prns)
    ax.set_ylabel("PRN")

    nticks = 8
    xt = np.linspace(0, len(expected_epochs) - 1, nticks).astype(int)
    ax.set_xticks(xt)
    ax.set_xticklabels([expected_epochs[j].strftime("%m-%d %H:%M") for j in xt], rotation=0)
    ax.set_xlabel("Time (epoch)")

    if title is None:
        if snr_col is None or snr_min is None:
            title = "Tracking timeline (black = present, white = missing)"
        else:
            title = f"Tracking timeline with SNR filter (black = {snr_col} ≥ {snr_min})"
    ax.set_title(title)

    plt.tight_layout()
    plt.show()
    return fig, ax

def plot_tracking_timeline_with_pivots(
    df: pd.DataFrame,
    selected_prns: list[str],
    sampling: pd.Timedelta = pd.Timedelta(seconds=30),
    snr_col: str = "S1",
    snr_min: float = 40.0,
    active_pivot: pd.Series | None = None,
    title: str | None = None,
):
    """
    Tracking timeline with SNR threshold + pivot satellites highlighted.

    Color meaning
    -------------
    0 = white : missing epoch OR SNR < snr_min
    1 = black : usable (SNR >= snr_min) but NOT selected as pivot
    2 = green : usable (SNR >= snr_min) AND selected as pivot
    3 = dark green : active pivot at this epoch (provided by `active_pivot`)

    Parameters
    ----------
    df : DataFrame
        MultiIndex (epoch, prn). Must contain snr_col.
    selected_prns : list[str]
        PRNs selected by the greedy algorithm (candidate pivots).
    sampling : Timedelta
        Expected sampling interval (e.g. 30 s).
    snr_col : str
        SNR column name (e.g., "S1").
    snr_min : float
        Minimum SNR threshold to consider a satellite usable.
    active_pivot : pd.Series | None
        Series indexed by expected_epochs, values are PRN (string) or None/NaN.
        If provided, it is used to mark state=3 (dark green).
    title : str | None
        Plot title override.

    Returns
    -------
    fig, ax, info
    """

    # ----------------------- checks -----------------------
    if not isinstance(df.index, pd.MultiIndex):
        raise ValueError("df must have MultiIndex ('epoch','prn').")
    if "epoch" not in df.index.names or "prn" not in df.index.names:
        raise ValueError("df index must have levels ('epoch','prn').")
    if snr_col not in df.columns:
        raise ValueError(f"snr_col='{snr_col}' not found in df.")

    # ---------------------- setup ------------------------
    epochs_obs = df.index.get_level_values("epoch")
    t0, t1 = epochs_obs.min(), epochs_obs.max()
    expected_epochs = pd.date_range(start=t0, end=t1, freq=sampling)

    prns = pd.Index(sorted(df.index.get_level_values("prn").unique()))
    prn_to_row = {p: i for i, p in enumerate(prns)}
    epoch_to_col = {t: j for j, t in enumerate(expected_epochs)}

    # 0 white, 1 black, 2 green, 3 dark green
    mat = np.zeros((len(prns), len(expected_epochs)), dtype=np.uint8)
    selected_set = set(selected_prns)

    # ------------------- fill usable states --------------
    for prn in prns:
        d = df.xs(prn, level="prn").sort_index()

        d = d.dropna(subset=[snr_col])
        if d.empty:
            continue

        d = d.loc[d.index.intersection(expected_epochs)]
        if d.empty:
            continue

        good = d[d[snr_col] >= snr_min]
        if good.empty:
            continue

        r = prn_to_row[prn]
        cols = [epoch_to_col[t] for t in good.index]
        mat[r, cols] = 2 if prn in selected_set else 1

    # ------------------ overlay active pivot -------------
    if active_pivot is not None:
        # Align to expected grid
        ap = active_pivot.reindex(expected_epochs)

        # Mark dark green only where:
        # - an active pivot is defined
        # - the pivot PRN is actually usable at this epoch (mat state >= 1)
        for t, prn in ap.dropna().items():
            r = prn_to_row.get(prn)
            c = epoch_to_col.get(t)
            if r is None or c is None:
                continue
            if mat[r, c] >= 1:
                mat[r, c] = 3  # dark green

    # ------------------------- plot -----------------------
    fig, ax = plt.subplots(figsize=(14, max(6, 0.25 * len(prns))))

    cmap = ListedColormap(["white", "black", "#2ca02c", "#006400"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)

    ax.imshow(mat, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)

    ax.set_yticks(np.arange(len(prns)))
    ax.set_yticklabels(prns)
    ax.set_ylabel("PRN")

    nticks = 8
    xt = np.linspace(0, len(expected_epochs) - 1, nticks).astype(int)
    ax.set_xticks(xt)
    ax.set_xticklabels([expected_epochs[j].strftime("%m-%d %H:%M") for j in xt])
    ax.set_xlabel("Time (epoch)")

    if title is None:
        title = (
            f"Tracking timeline (snr_min={snr_min:g}) "
            f"+ candidate pivots (green) + active pivot (dark green)"
        )
    ax.set_title(title)

    plt.tight_layout()
    plt.show()

    info = {"expected_epochs": expected_epochs, "prns": prns}
    return fig, ax, info



def plot_gnss_timeseries_by_prn(
    df: pd.DataFrame,
    y: str,
    gap: pd.Timedelta = pd.Timedelta(minutes=30),
    label_arcs: bool = True,
    arc_label_fontsize: int = 8,
    arc_label_xytext: tuple[int, int] = (3, 3),  # offset in points
    show_legend: bool = False,
    legend_outside: bool = True,
    legend_ncol: int = 4,
    title: str | None = None,
    xlabel: str = "Time (epoch)",
    ylabel: str | None = None,
    figsize: tuple[int, int] = (10, 6),
):
    """
    Plot a GNSS time series by satellite PRN with gap handling.

    This function is generic: it works for any DataFrame that:
      - is indexed by a MultiIndex with levels ('epoch', 'prn')
      - contains a numeric column `y`

    Parameters
    ----------
    df : pandas.DataFrame
        MultiIndex (epoch, prn).
    y : str
        Column name to plot.
    gap : pandas.Timedelta
        Time gap threshold used to split continuous segments ("arcs").
    label_arcs : bool
        If True, write PRN text at the start of each arc.
    arc_label_fontsize : int
        Font size for PRN labels.
    arc_label_xytext : tuple
        (dx, dy) offset in points for the arc labels.
    show_legend : bool
        If True, show legend (often unnecessary if label_arcs=True).
    legend_outside : bool
        If True, place legend outside the axes (right side).
    legend_ncol : int
        Legend columns.
    title : str or None
        Plot title. If None, a default title is generated.
    xlabel : str
        X axis label.
    ylabel : str or None
        Y axis label. If None, uses `y`.
    figsize : tuple
        Figure size.

    Returns
    -------
    fig, ax
    """

    # --- checks
    if not isinstance(df.index, pd.MultiIndex):
        raise ValueError("df must have a MultiIndex (epoch, prn).")

    if "epoch" not in df.index.names or "prn" not in df.index.names:
        raise ValueError("df index must have levels ('epoch','prn').")

    if y not in df.columns:
        raise ValueError(f"Column '{y}' not found in df.")

    prns = df.index.get_level_values("prn").unique()

    fig, ax = plt.subplots(figsize=figsize)

    for prn in prns:
        # Extract one PRN time series (index becomes epochs)
        data = df.xs(prn, level="prn")[[y]].copy().sort_index()

        # Skip if nothing to plot
        if not data[y].notna().any():
            continue

        # Drop NaNs for segmentation + plotting
        data = data.dropna(subset=[y])
        if data.empty:
            continue

        # Segment id based on gaps
        dt = data.index.to_series().diff()
        seg_id = (dt > gap).cumsum()

        color = None

        for _, seg in data.groupby(seg_id):
            if seg.empty:
                continue

            line, = ax.plot(
                seg.index,
                seg[y].to_numpy(),
                color=color,
                label=prn if (color is None and show_legend) else None,
            )

            if color is None:
                color = line.get_color()

            # Label at the beginning of each arc (segment)
            if label_arcs:
                x0 = seg.index[0]
                y0 = float(seg[y].iloc[0])
                ax.annotate(
                    prn,
                    xy=(x0, y0),
                    xytext=arc_label_xytext,
                    textcoords="offset points",
                    fontsize=arc_label_fontsize,
                    color=color,
                    ha="left",
                    va="bottom",
                )

    # labels
    if title is None:
        title = f"Time series by satellite PRN ({y})"
    ax.set_title(title)

    ax.set_xlabel(xlabel)

    if ylabel is None:
        ylabel = y
    ax.set_ylabel(ylabel)

    # legend
    if show_legend:
        if legend_outside:
            ax.legend(
                title="PRN",
                bbox_to_anchor=(1.02, 1),
                loc="upper left",
                borderaxespad=0,
                ncol=legend_ncol,
            )
            plt.tight_layout(rect=[0, 0, 0.85, 1])
        else:
            ax.legend(title="PRN", ncol=legend_ncol)
            plt.tight_layout()
    else:
        plt.tight_layout()

    plt.show()
    return fig, ax

# -------------------------------------------------------------------------
# Plot Single Differences (SD) by satellite PRN (with gap handling)
#
# Input expected:
#   df_SD : DataFrame indexed by (epoch, prn) with columns like "SD_L1", "SD_C1", ...
# -------------------------------------------------------------------------

def plot_gnss_sd_by_prn(df_SD: pd.DataFrame, observable: str = "SD_L1", **kwargs):
    """
    Convenience wrapper for SD plots (keeps the TP API stable).
    """
    return plot_gnss_timeseries_by_prn(
        df=df_SD,
        y=observable,
        title=f"Single differences time series by satellite PRN ({observable})",
        ylabel="Value (cycles for SD_L*, meters for SD_C*/SD_P*)",
        **kwargs,
    )





def plot_sd_derivative_by_prn(
    df_SD: pd.DataFrame,
    obs: str,
    gap: pd.Timedelta = pd.Timedelta(minutes=30),
    label_arcs: bool = True,
    normalize_by_dt: bool = True,
    elev_col: str | None = None,
    cutoff_deg: float | None = None,
    title: str | None = None,
    phase_to_meters: bool = False,
    wavelength_m: float | None = None,
):
    """
    Plot temporal derivative of a single-difference observable by PRN.

    Notes (pedagogical)
    -------------------
    - If obs is a carrier-phase SD_L* (stored in cycles), you may set
      phase_to_meters=True to convert cycles -> meters using the wavelength.
      This is ONLY for comparison with code (meters); it does not remove
      ambiguities.

    Parameters
    ----------
    df_SD : DataFrame
        MultiIndex (epoch, prn). Must contain column `obs`.
    obs : str
        Column to differentiate (e.g., "SD_L1" cycles, "SD_C1" meters).
    phase_to_meters : bool
        If True and obs starts with "SD_L", convert cycles to meters before
        differencing.
    wavelength_m : float or None
        Override wavelength (meters) used for cycles->meters conversion.
        If None and obs is SD_L1/SD_L2/SD_L5, use conv.L*_WAVELENGTH.
    """

    if not isinstance(df_SD.index, pd.MultiIndex):
        raise ValueError("df_SD must have a MultiIndex (epoch, prn).")

    if obs not in df_SD.columns:
        raise ValueError(f"Column '{obs}' not found in df_SD.")

    if elev_col is not None and elev_col not in df_SD.columns:
        raise ValueError(f"elev_col='{elev_col}' not found in df_SD.")

    # ------------------------------------------------------------------
    # Decide units / optional conversion for phase observables
    # ------------------------------------------------------------------
    is_phase = obs.startswith("SD_L")
    use_phase_to_m = bool(phase_to_meters and is_phase)

    if use_phase_to_m:
        # Auto wavelength from observable name unless overridden
        if wavelength_m is None:
            if obs == "SD_L1":
                wavelength_m = conv.L1_WAVELENGTH
            elif obs == "SD_L2":
                wavelength_m = conv.L2_WAVELENGTH
            elif obs == "SD_L5":
                wavelength_m = conv.L5_WAVELENGTH
            else:
                raise ValueError(
                    f"Cannot infer wavelength for '{obs}'. "
                    "Provide wavelength_m explicitly."
                )

        y_unit = "m"
        obs_label = f"{obs} (cycles→m)"
    else:
        # Keep original units
        y_unit = "cycles" if is_phase else "m"
        obs_label = obs

    prns = df_SD.index.get_level_values("prn").unique()
    fig, ax = plt.subplots(figsize=(12, 6))

    for prn in prns:
        # Extract one PRN series
        s = df_SD.xs(prn, level="prn")[[obs]].copy()
        if elev_col is not None:
            s[elev_col] = df_SD.xs(prn, level="prn")[elev_col]

        s = s.dropna(subset=[obs])
        if len(s) < 2:
            continue

        # Optional cutoff filter
        if elev_col is not None and cutoff_deg is not None:
            s = s[s[elev_col] >= cutoff_deg]
            if len(s) < 2:
                continue

        s = s.sort_index()

        # Optional conversion cycles -> meters (local only)
        if use_phase_to_m:
            s[obs] = s[obs] * float(wavelength_m)

        # Segment arcs (gaps)
        dt_epoch = s.index.to_series().diff()
        seg_id = (dt_epoch > gap).cumsum()

        color = None
        for _, seg in s.groupby(seg_id):
            if len(seg) < 2:
                continue

            # Derivative
            y = seg[obs].diff()

            if normalize_by_dt:
                dt_s = seg.index.to_series().diff().dt.total_seconds()
                y = y / dt_s  # units per second

            seg_plot = pd.DataFrame({"y": y}).dropna()
            if seg_plot.empty:
                continue

            line, = ax.plot(seg_plot.index, seg_plot["y"], color=color)
            if color is None:
                color = line.get_color()

            if label_arcs:
                x0 = seg_plot.index[0]
                y0 = seg_plot["y"].iloc[0]
                ax.text(x0, y0, prn, color=color, fontsize=8, ha="left", va="bottom")

    ax.axhline(0, color="black", linewidth=0.8)

    # Y label depending on normalization
    rate = "per second" if normalize_by_dt else "per epoch-step"
    ax.set_ylabel(f"Δ({obs_label}) [{y_unit}] {rate}")

    if title is None:
        cutoff_txt = ""
        if elev_col is not None and cutoff_deg is not None:
            cutoff_txt = f" (cutoff {cutoff_deg:.0f}°)"
        title = f"Temporal derivative of {obs_label}{cutoff_txt}"

    ax.set_title(title)
    ax.set_xlabel("Time (epoch)")
    plt.tight_layout()
    plt.show()
    return fig, ax



def plot_series(df, col1, col2=None, coeff1=1.0, coeff2=1.0, seuil=3600, renderer="browser"):
    """
    Affiche les séries temporelles pour chaque satellite.
    Si col2 est fourni, affiche la série : coeff1 * col1 - coeff2 * col2.
    Sinon, affiche la série de la colonne col1 directement.
    La série est découpée en segments lorsqu'un "trou" (écart > seuil) est détecté.

    Paramètres
    ----------
    df : DataFrame
        DataFrame avec un index multi-niveaux contenant au moins le niveau 'prn'
        et les colonnes col1 (et éventuellement col2).
    col1 : str
        Nom de la première colonne.
    col2 : str, optional
        Nom de la seconde colonne (optionnel). Si None, on affiche col1.
    coeff1 : float
        Coefficient multiplicateur pour la première colonne (défaut 1.0).
    coeff2 : float
        Coefficient multiplicateur pour la seconde colonne (défaut 1.0).
    seuil : float
        Seuil en secondes pour considérer un "trou" dans la série (défaut 3600).
    renderer : str
        Renderer Plotly (ex : "browser" ou "iframe").

    Retourne
    --------
    fig : Figure
        Figure Plotly contenant les courbes tracées.
    """
    # Configuration du renderer
    pio.renderers.default = renderer

    # Extraction et tri de la liste des satellites (niveau 'prn')
    satellites = sorted(df.index.get_level_values('prn').unique())

    fig = go.Figure()

    for sat in satellites:
        # Sélection des données pour le satellite
        sub_df = df.xs(sat, level='prn')

        # Calcul de la série
        if col2 is not None:
            ts = (sub_df[col1] * coeff1 - sub_df[col2] * coeff2).dropna()
            y_label = f"{coeff1}*{col1} - {coeff2}*{col2}"
        else:
            ts = (sub_df[col1] * coeff1).dropna()
            y_label = f"{coeff1}*{col1}"

        dates = ts.index

        # Calcul des écarts entre dates successives (en secondes)
        ecarts = dates.to_series().diff().dt.total_seconds()

        # Identification des positions où l'écart dépasse le seuil
        breaks = np.where(ecarts > seuil)[0]

        # Découpage de la série en segments
        segments_dates = np.split(dates, breaks)
        segments_values = np.split(ts.values, breaks)

        # Reconstitution de la série en insérant des valeurs None pour créer des "trous"
        combined_dates = []
        combined_values = []
        for seg_dates, seg_values in zip(segments_dates, segments_values):
            if len(seg_dates) > 0:
                combined_dates.extend(seg_dates)
                combined_values.extend(seg_values)
                # Insertion d'un "trou" (None) pour forcer une cassure
                combined_dates.append(None)
                combined_values.append(None)
        # Supprimer le dernier trou ajouté si présent
        if combined_dates and combined_dates[-1] is None:
            combined_dates = combined_dates[:-1]
            combined_values = combined_values[:-1]

        # Ajout d'une trace unique par satellite
        fig.add_trace(go.Scatter(
            x=combined_dates,
            y=combined_values,
            mode='lines',
            name=f"Sat {sat}",
            legendgroup=f"sat_{sat}"
        ))

    # Boutons pour sélectionner/désélectionner les traces
    buttons = [
        dict(
            label="Tout sélectionner",
            method="update",
            args=[{"visible": [True] * len(fig.data)}]
        ),
        dict(
            label="Tout désélectionner",
            method="update",
            args=[{"visible": ["legendonly"] * len(fig.data)}]
        )
    ]

    # Mise à jour de la mise en page du graphique
    fig.update_layout(
        title=dict(
            text=y_label+" par satellite",
            x=0.5,
            xanchor="center"
        ),
        xaxis_title="Date",
        yaxis_title=y_label,
        hovermode="x unified",
        updatemenus=[dict(
            type="buttons",
            direction="left",
            buttons=buttons,
            pad={"r": 10, "t": 10},
            showactive=True,
            x=0.0,
            xanchor="left",
            y=1.1,   # Position des boutons en dehors de la zone de tracé, en haut
            yanchor="top"
        )],
        margin=dict(t=120)  # Augmentation de la marge supérieure pour que le titre et les boutons ne se chevauchent pas
    )

    fig.show()
    return fig



def build_residual_dataframe(
    df_obs_used: pd.DataFrame,
    A: np.ndarray,
    B: np.ndarray,
    dP_est: np.ndarray,
    predicted_col: str = "predicted_m",
    residual_col: str = "residual_m",
) -> pd.DataFrame:
    """
    Build a residual DataFrame aligned with the observation DataFrame actually
    used in the least-squares adjustment.

    Parameters
    ----------
    df_obs_used : pd.DataFrame
        Observation DataFrame used to build A and B.
        Its row order must match the row order of A and B.
    A : np.ndarray
        Design matrix.
    B : np.ndarray
        Observation vector.
    dP_est : np.ndarray
        Estimated parameter vector.
    predicted_col : str, default="predicted_m"
        Name of the column storing the predicted values A @ dP_est.
    residual_col : str, default="residual_m"
        Name of the column storing the residuals B - A @ dP_est.

    Returns
    -------
    pd.DataFrame
        Copy of df_obs_used enriched with predicted values and residuals.
    """
    if not isinstance(df_obs_used, pd.DataFrame):
        raise TypeError("df_obs_used must be a pandas DataFrame.")

    A = np.asarray(A)
    B = np.asarray(B).ravel()
    dP_est = np.asarray(dP_est).ravel()

    if A.ndim != 2:
        raise ValueError("A must be a 2D array.")

    n_obs = len(df_obs_used)

    if A.shape[0] != n_obs:
        raise ValueError(
            "Row mismatch: A and df_obs_used must contain the same number "
            "of observations."
        )

    if len(B) != n_obs:
        raise ValueError(
            "Length mismatch: B and df_obs_used must contain the same number "
            "of observations."
        )

    if A.shape[1] != len(dP_est):
        raise ValueError(
            "Dimension mismatch: the number of columns in A must match the "
            "length of dP_est."
        )

    predicted = A @ dP_est
    residuals = B - predicted

    df_res = df_obs_used.copy()
    df_res[predicted_col] = predicted
    df_res[residual_col] = residuals
    df_res["abs_residual_m"] = np.abs(residuals)

    return df_res


def _plot_residual_analysis_from_df(
    df_residuals: pd.DataFrame,
    figure_title: str | None = None,
    save_path=None,
    P_est=None,
    P_rnx_header=None,
    residual_col: str = "residual_m",
    predicted_col: str = "predicted_m",
    top_left_plot: str = "timeseries",
    prn_gap: pd.Timedelta = pd.Timedelta(minutes=30),
    label_arcs: bool = True,
):
    """
    New implementation: create a residual diagnostic figure from a residual
    DataFrame.
    """
    if not isinstance(df_residuals, pd.DataFrame):
        raise TypeError("df_residuals must be a pandas DataFrame.")

    if residual_col not in df_residuals.columns:
        raise ValueError(f"Column '{residual_col}' not found in df_residuals.")

    if predicted_col not in df_residuals.columns:
        raise ValueError(f"Column '{predicted_col}' not found in df_residuals.")

    v_est = df_residuals[residual_col].to_numpy(dtype=float)
    b_est = df_residuals[predicted_col].to_numpy(dtype=float)

    mean_val = np.mean(v_est)
    variance = np.var(v_est)
    std_dev = np.std(v_est)
    rms_val = np.sqrt(np.mean(v_est ** 2))
    median_val = np.median(v_est)
    skewness = stats.skew(v_est)
    kurtosis = stats.kurtosis(v_est)

    nrows = 4 if (P_est is not None and P_rnx_header is not None) else 3
    height_ratios = [1, 1, 0.5, 0.5] if nrows == 4 else [1, 1, 0.5]

    fig = plt.figure(figsize=(15, 15))
    gs = fig.add_gridspec(nrows, 2, height_ratios=height_ratios)

    # ------------------------------------------------------------------
    # 1. Top-left panel
    # ------------------------------------------------------------------
    ax_top_left = fig.add_subplot(gs[0, 0])

    if top_left_plot == "timeseries":
        ax_top_left.scatter(np.arange(len(v_est)), v_est, color="green", s=12)
        ax_top_left.set_title("Time Series of Residuals")
        ax_top_left.set_xlabel("Observation index")
        ax_top_left.set_ylabel("Residuals (m)")

    elif top_left_plot == "by_prn":
        if not isinstance(df_residuals.index, pd.MultiIndex):
            raise ValueError(
                "For top_left_plot='by_prn', df_residuals must use a MultiIndex "
                "with levels ('epoch', 'prn')."
            )

        if "epoch" not in df_residuals.index.names or "prn" not in df_residuals.index.names:
            raise ValueError(
                "For top_left_plot='by_prn', df_residuals index must contain "
                "levels ('epoch', 'prn')."
            )

        prns = df_residuals.index.get_level_values("prn").unique()

        for prn in prns:
            data = df_residuals.xs(prn, level="prn")[[residual_col]].copy()
            data = data.sort_index().dropna(subset=[residual_col])

            if data.empty:
                continue

            dt = data.index.to_series().diff()
            seg_id = (dt > prn_gap).cumsum()

            color = None

            for _, seg in data.groupby(seg_id):
                if seg.empty:
                    continue

                (line,) = ax_top_left.plot(
                    seg.index,
                    seg[residual_col].to_numpy(),
                    color=color,
                    linewidth=1.0,
                )

                if color is None:
                    color = line.get_color()

                if label_arcs:
                    x0 = seg.index[0]
                    y0 = float(seg[residual_col].iloc[0])
                    ax_top_left.annotate(
                        prn,
                        xy=(x0, y0),
                        xytext=(3, 3),
                        textcoords="offset points",
                        fontsize=8,
                        color=color,
                        ha="left",
                        va="bottom",
                    )

        ax_top_left.set_title("Residuals by Satellite PRN")
        ax_top_left.set_xlabel("Epoch")
        ax_top_left.set_ylabel("Residuals (m)")

    else:
        raise ValueError(
            "Unsupported top_left_plot. Use 'timeseries' or 'by_prn'."
        )

    # ------------------------------------------------------------------
    # 2. Histogram
    # ------------------------------------------------------------------
    ax_hist = fig.add_subplot(gs[0, 1])
    ax_hist.hist(v_est, bins=30, edgecolor="black")
    ax_hist.set_title("Histogram of Residuals")
    ax_hist.set_xlabel("Residuals (m)")
    ax_hist.set_ylabel("Number of Observations")

    # ------------------------------------------------------------------
    # 3. Q-Q plot
    # ------------------------------------------------------------------
    ax_qq = fig.add_subplot(gs[1, 0])
    sm.qqplot(v_est, line="s", ax=ax_qq)
    ax_qq.set_title("Q-Q Plot of Residuals")

    # ------------------------------------------------------------------
    # 4. Residuals vs predicted values
    # ------------------------------------------------------------------
    ax_scatter = fig.add_subplot(gs[1, 1])
    ax_scatter.scatter(b_est, v_est, alpha=0.7, color="darkorange", s=12)
    ax_scatter.axhline(0, color="red", linestyle="--")
    ax_scatter.set_xlabel("Predicted values")
    ax_scatter.set_ylabel("Residuals (m)")
    ax_scatter.set_title("Residuals vs. Predicted Values")

    # ------------------------------------------------------------------
    # 5. Statistics text box
    # ------------------------------------------------------------------
    ax_stats = fig.add_subplot(gs[2, :])
    ax_stats.axis("off")

    stats_text = (
        f"Mean       : {mean_val:.4f}\n"
        f"Median     : {median_val:.4f}\n"
        f"Variance   : {variance:.4f}\n"
        f"Std Dev    : {std_dev:.4f}\n"
        f"RMS        : {rms_val:.4f}\n"
        f"Skewness   : {skewness:.4f}\n"
        f"Kurtosis   : {kurtosis:.4f}"
    )

    ax_stats.text(
        0.5,
        0.5,
        stats_text,
        transform=ax_stats.transAxes,
        fontsize=14,
        verticalalignment="center",
        horizontalalignment="center",
        bbox=dict(facecolor="wheat", edgecolor="black", boxstyle="round,pad=1"),
    )
    ax_stats.set_title("Residual Statistics", fontsize=16)

    # ------------------------------------------------------------------
    # 6. Optional position information
    # ------------------------------------------------------------------
    if nrows == 4:
        dist = np.sqrt(np.sum((P_est - P_rnx_header) ** 2))
        E, N, U = conv.xyz2enu(
            P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],
            P_est[0], P_est[1], P_est[2]
        )

        extra_text = (
            "Distance between estimated position and initial RINEX header "
            f"position: {dist:.4f} m\n"
            "Local ENU coordinates:\n"
            f"  East (E)  : {E[0]:.4f} m\n"
            f"  North (N) : {N[0]:.4f} m\n"
            f"  Up (U)    : {U[0]:.4f} m"
        )

        ax_extra = fig.add_subplot(gs[3, :])
        ax_extra.axis("off")
        ax_extra.text(
            0.5,
            0.5,
            extra_text,
            transform=ax_extra.transAxes,
            fontsize=14,
            verticalalignment="center",
            horizontalalignment="center",
            bbox=dict(facecolor="lightcyan", edgecolor="black", boxstyle="round,pad=1"),
        )
        ax_extra.set_title("Position Information", fontsize=16)

    if figure_title is not None:
        fig.suptitle(figure_title, fontsize=20)
        plt.subplots_adjust(top=0.92)

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    if save_path is not None:
        fig.savefig(save_path)
        print(f"Figure saved to: {save_path}")

    plt.show()
    return fig


def _plot_residual_analysis_legacy(
    A,
    B,
    dP_est,
    figure_title=None,
    save_path=None,
    P_est=None,
    P_rnx_header=None,
):
    """
    Legacy implementation preserved for backward compatibility.

    Behavior is intentionally kept as close as possible to the published
    version, with only non-breaking micro-corrections.
    """
    # Compute residuals and predicted values
    v_est = B - A @ dP_est
    b_est = A @ dP_est

    # Compute residual statistics
    mean_val = np.mean(v_est)
    variance = np.var(v_est)
    std_dev = np.std(v_est)
    skewness = stats.skew(v_est)
    kurtosis = stats.kurtosis(v_est)

    # Create figure with GridSpec
    nrows = 4 if (P_est is not None and P_rnx_header is not None) else 3
    height_ratios = [1, 1, 0.5, 0.5] if nrows == 4 else [1, 1, 0.5]

    fig = plt.figure(figsize=(15, 15))
    gs = fig.add_gridspec(nrows, 2, height_ratios=height_ratios)

    # 1. Time series of residuals (points only) – Top left
    ax_time = fig.add_subplot(gs[0, 0])
    ax_time.scatter(np.arange(len(v_est)), v_est, color="green")
    ax_time.set_title("Time Series of Residuals")
    ax_time.set_xlabel("Time / Index")
    ax_time.set_ylabel("Residuals (m)")

    # 2. Histogram of residuals – Top right
    ax_hist = fig.add_subplot(gs[0, 1])
    ax_hist.hist(v_est, bins=30, edgecolor="black")
    ax_hist.set_title("Histogram of Residuals")
    ax_hist.set_xlabel("Residuals")
    ax_hist.set_ylabel("Number of Observations")

    # 3. Q-Q Plot of residuals
    ax_qq = fig.add_subplot(gs[1, 0])
    sm.qqplot(v_est, line="s", ax=ax_qq)
    ax_qq.set_title("Q-Q Plot of Residuals")

    # 4. Residuals vs. predicted values plot
    ax_scatter = fig.add_subplot(gs[1, 1])
    ax_scatter.scatter(b_est, v_est, alpha=0.7, color="darkorange")
    ax_scatter.axhline(0, color="red", linestyle="--")
    ax_scatter.set_xlabel("Predicted Values")
    ax_scatter.set_ylabel("Residuals (m)")
    ax_scatter.set_title("Residuals vs. Predicted Values")

    # 5. Statistics text box
    ax_stats = fig.add_subplot(gs[2, :])
    ax_stats.axis("off")
    stats_text = (
        f"Mean       : {mean_val:.4f}\n"
        f"Variance   : {variance:.4f}\n"
        f"Std Dev    : {std_dev:.4f}\n"
        f"Skewness   : {skewness:.4f}\n"
        f"Kurtosis   : {kurtosis:.4f}"
    )
    ax_stats.text(
        0.5,
        0.5,
        stats_text,
        transform=ax_stats.transAxes,
        fontsize=14,
        verticalalignment="center",
        horizontalalignment="center",
        bbox=dict(facecolor="wheat", edgecolor="black", boxstyle="round,pad=1"),
    )
    ax_stats.set_title("Residual Statistics", fontsize=16)

    # 6. Position information text box – optional
    if nrows == 4:
        dist = np.sqrt(np.sum((P_est - P_rnx_header) ** 2))
        E, N, U = conv.xyz2enu(
            P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],
            P_est[0], P_est[1], P_est[2]
        )
        extra_text = (
            f"Distance between estimated position and initial RINEX header position: {dist:.4f} m\n"
            f"Local ENU coordinates:\n"
            f"  East (E)  : {E[0]:.4f} m\n"
            f"  North (N) : {N[0]:.4f} m\n"
            f"  Up (U)    : {U[0]:.4f} m"
        )
        ax_extra = fig.add_subplot(gs[3, :])
        ax_extra.axis("off")
        ax_extra.text(
            0.5,
            0.5,
            extra_text,
            transform=ax_extra.transAxes,
            fontsize=14,
            verticalalignment="center",
            horizontalalignment="center",
            bbox=dict(facecolor="lightcyan", edgecolor="black", boxstyle="round,pad=1"),
        )
        ax_extra.set_title("Position Information", fontsize=16)

    if figure_title is not None:
        fig.suptitle(figure_title, fontsize=20)
        plt.subplots_adjust(top=0.92)

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    if save_path is not None:
        fig.savefig(save_path)
        print(f"Figure saved to: {save_path}")

    plt.show()
    return fig


def plot_residual_analysis(
    *args,
    df_residuals: pd.DataFrame | None = None,
    figure_title: str | None = None,
    save_path=None,
    P_est=None,
    P_rnx_header=None,
    residual_col: str = "residual_m",
    predicted_col: str = "predicted_m",
    top_left_plot: str = "timeseries",
    prn_gap: pd.Timedelta = pd.Timedelta(minutes=30),
    label_arcs: bool = True,
    **kwargs,
):
    """
    Public residual-analysis function with backward-compatible routing.

    Preferred API
    -------------
    plot_residual_analysis(
        df_residuals=df_residuals,
        figure_title=...,
        save_path=...,
        P_est=...,
        P_rnx_header=...,
        top_left_plot="timeseries" or "by_prn",
    )

    Legacy API (deprecated)
    -----------------------
    plot_residual_analysis(
        A, B, dP_est,
        figure_title=...,
        save_path=...,
        P_est=...,
        P_rnx_header=...,
    )
    """
    if kwargs:
        unexpected = ", ".join(kwargs.keys())
        raise TypeError(f"Unexpected keyword argument(s): {unexpected}")

    # --------------------------------------------------------------
    # New API
    # --------------------------------------------------------------
    if df_residuals is not None:
        if len(args) != 0:
            raise TypeError(
                "When 'df_residuals' is provided, do not pass positional "
                "arguments A, B, dP_est."
            )

        return _plot_residual_analysis_from_df(
            df_residuals=df_residuals,
            figure_title=figure_title,
            save_path=save_path,
            P_est=P_est,
            P_rnx_header=P_rnx_header,
            residual_col=residual_col,
            predicted_col=predicted_col,
            top_left_plot=top_left_plot,
            prn_gap=prn_gap,
            label_arcs=label_arcs,
        )

    # --------------------------------------------------------------
    # Legacy API
    # --------------------------------------------------------------
    if len(args) >= 3:
        A, B, dP_est = args[:3]

        warnings.warn(
            "Calling plot_residual_analysis(A, B, dP_est, ...) is deprecated "
            "and will be removed in a future release. "
            "Please use build_residual_dataframe(...) and then call "
            "plot_residual_analysis(df_residuals=..., ...).",
            DeprecationWarning,
            stacklevel=2,
        )

        return _plot_residual_analysis_legacy(
            A=A,
            B=B,
            dP_est=dP_est,
            figure_title=figure_title,
            save_path=save_path,
            P_est=P_est,
            P_rnx_header=P_rnx_header,
        )

    raise TypeError(
        "Invalid call. Use either:\n"
        "  plot_residual_analysis(df_residuals=..., ...)\n"
        "or the deprecated legacy form:\n"
        "  plot_residual_analysis(A, B, dP_est, ...)"
    )



# =============================================================================
# Helpers (epoch grid + per-PRN SNR extraction)
# =============================================================================

def expected_epochs(t0: pd.Timestamp, t1: pd.Timestamp, sampling: pd.Timedelta) -> pd.DatetimeIndex:
    """Return the expected epoch grid (inclusive endpoints) for a given sampling."""
    return pd.date_range(start=t0, end=t1, freq=sampling)


def prn_snr_on_grid(
    df: pd.DataFrame,
    prn: str,
    snr_col: str,
    epochs: pd.DatetimeIndex,
) -> pd.Series:
    """
    Return SNR series for one PRN reindexed on `epochs`.
    Missing epochs -> NaN.
    """
    try:
        s = df.xs(prn, level="prn")[snr_col]
    except KeyError:
        return pd.Series(index=epochs, dtype="float64")
    return s.reindex(epochs)


# =============================================================================
# Coverage predicates (the only one you really need)
# =============================================================================

def prn_covers_interval(
    df: pd.DataFrame,
    prn: str,
    t0: pd.Timestamp,
    t1: pd.Timestamp,
    snr_col: str,
    snr_min: float,
    sampling: pd.Timedelta,
) -> bool:
    """
    True if PRN has non-NaN SNR >= snr_min at ALL expected epochs in [t0, t1].
    """
    epochs = expected_epochs(t0, t1, sampling)
    s = prn_snr_on_grid(df, prn, snr_col, epochs)
    ok = s.notna() & (s >= snr_min)
    return bool(ok.all())


def best_prn_for_interval(
    df: pd.DataFrame,
    t0: pd.Timestamp,
    t1: pd.Timestamp,
    snr_col: str,
    snr_min: float,
    sampling: pd.Timedelta,
    pool_prns: list[str],
) -> str | None:
    """
    Pick a PRN from pool_prns that fully covers [t0, t1] with SNR>=snr_min.
    Tie-breaker: highest mean SNR over the interval.

    Returns None if no PRN can fully cover the interval.
    """
    epochs = expected_epochs(t0, t1, sampling)
    best_prn = None
    best_score = -np.inf

    for prn in pool_prns:
        s = prn_snr_on_grid(df, prn, snr_col, epochs)
        ok = s.notna() & (s >= snr_min)
        if not bool(ok.all()):
            continue
        score = float(s.mean())
        if score > best_score:
            best_prn = prn
            best_score = score

    return best_prn


# =============================================================================
# 1) Greedy set cover (candidate pivots)
# =============================================================================

def greedy_pivot_set_cover(
    df: pd.DataFrame,
    snr_col: str,
    snr_min: float,
    sampling: pd.Timedelta = pd.Timedelta(seconds=30),
    require_full_coverage: bool = True,
    return_diagnostics: bool = True,
):
    """
    Greedy set-cover selection of PRNs to cover ALL expected epochs with SNR >= snr_min.

    Guarantee:
    - If require_full_coverage=True: raise RuntimeError if full coverage is impossible.

    Returns
    -------
    selected_prns : list[str]
    diagnostics : dict (optional)
    """
    if not isinstance(df.index, pd.MultiIndex):
        raise ValueError("df must have MultiIndex ('epoch','prn').")
    if snr_col not in df.columns:
        raise ValueError(f"snr_col='{snr_col}' not found.")
    if "epoch" not in df.index.names or "prn" not in df.index.names:
        raise ValueError("df index must have ('epoch','prn').")

    epochs_obs = df.index.get_level_values("epoch")
    t0, t1 = epochs_obs.min(), epochs_obs.max()
    grid = expected_epochs(t0, t1, sampling)
    universe = set(grid)

    prns = sorted(df.index.get_level_values("prn").unique())
    cover: dict[str, set[pd.Timestamp]] = {}
    mean_snr: dict[str, float] = {}

    for prn in prns:
        d = df.xs(prn, level="prn").sort_index()
        d = d.dropna(subset=[snr_col])
        good = d[d[snr_col] >= snr_min]
        good_idx = good.index.intersection(grid)
        s = set(good_idx)
        if s:
            cover[prn] = s
            mean_snr[prn] = float(good.loc[good_idx, snr_col].mean())

    union_all = set().union(*cover.values()) if cover else set()
    impossible = sorted(universe - union_all)
    if require_full_coverage and impossible:
        raise RuntimeError(
            "Full coverage is IMPOSSIBLE with the current SNR threshold.\n"
            f"- snr_col={snr_col}, snr_min={snr_min}\n"
            f"- sampling={sampling}\n"
            f"- missing epochs (first 10): {impossible[:10]}\n"
            "Lower snr_min (or change quality rule) if you need guaranteed coverage."
        )

    uncovered = set(universe)
    selected: list[str] = []

    while uncovered and cover:
        best_prn = None
        best_gain = -1
        best_snr = -np.inf

        for prn, s in cover.items():
            gain = len(s & uncovered)
            ms = mean_snr.get(prn, -np.inf)
            if gain > best_gain or (gain == best_gain and ms > best_snr):
                best_prn, best_gain, best_snr = prn, gain, ms

        if best_prn is None or best_gain <= 0:
            break

        selected.append(best_prn)
        uncovered -= cover[best_prn]
        del cover[best_prn]

    uncovered_sorted = pd.DatetimeIndex(sorted(uncovered))
    coverage_ratio = 1.0 - (len(uncovered_sorted) / len(grid))

    if require_full_coverage and len(uncovered_sorted) > 0:
        raise RuntimeError(
            "Greedy selection did not achieve full coverage (unexpected if feasibility passed).\n"
            f"- coverage_ratio={coverage_ratio:.6f}\n"
            f"- uncovered epochs={len(uncovered_sorted)} (first 10: {list(uncovered_sorted[:10])})\n"
        )

    if not return_diagnostics:
        return selected, coverage_ratio

    diagnostics = {
        "coverage_ratio": coverage_ratio,
        "uncovered_epochs": uncovered_sorted,
        "n_uncovered": int(len(uncovered_sorted)),
        "snr_min": float(snr_min),
        "sampling": sampling,
        "t0": t0,
        "t1": t1,
        "n_selected": len(selected),
    }
    return selected, diagnostics


# =============================================================================
# 2) Stable active pivot schedule (hysteresis)
# =============================================================================

def build_active_pivot_schedule(
    df: pd.DataFrame,
    selected_prns: list[str],
    snr_col: str = "S1",
    snr_min: float = 40.0,
    sampling: pd.Timedelta = pd.Timedelta(seconds=30),
    switch_margin_db: float = 2.0,
) -> pd.Series:
    """
    Build a stable active pivot schedule from a pre-selected PRN set.
    Switch only if the best candidate is better than current by switch_margin_db.

    Returns
    -------
    active_pivot : pd.Series indexed by expected epochs, values are PRN or None
    """
    if not isinstance(df.index, pd.MultiIndex):
        raise ValueError("df must have MultiIndex ('epoch','prn').")
    if snr_col not in df.columns:
        raise ValueError(f"snr_col='{snr_col}' not found.")
    if "epoch" not in df.index.names or "prn" not in df.index.names:
        raise ValueError("df index must have ('epoch','prn').")

    epochs_obs = df.index.get_level_values("epoch")
    t0, t1 = epochs_obs.min(), epochs_obs.max()
    grid = expected_epochs(t0, t1, sampling)

    active = None
    out = {}

    for epoch in grid:
        try:
            row = df.xs(epoch, level="epoch")
        except KeyError:
            out[epoch] = None
            active = None
            continue

        cand = row.loc[row.index.isin(selected_prns)]
        cand = cand.dropna(subset=[snr_col])
        cand = cand[cand[snr_col] >= snr_min]

        if cand.empty:
            out[epoch] = None
            active = None
            continue

        best_prn = cand[snr_col].idxmax()
        best_snr = float(cand.loc[best_prn, snr_col])

        if active is None:
            active = best_prn
        else:
            if active not in cand.index:
                active = best_prn
            else:
                current_snr = float(cand.loc[active, snr_col])
                if best_snr > current_snr + switch_margin_db:
                    active = best_prn

        out[epoch] = active

    return pd.Series(out, name="active_pivot")


# =============================================================================
# 3) Schedule -> segments + coverage check
# =============================================================================

def pivot_schedule_to_segments(
    active_pivot: pd.Series,
    sampling: pd.Timedelta,
    drop_none: bool = True,
) -> pd.DataFrame:
    """
    Convert an active pivot schedule to segments:
    columns = prn, start, end, duration, n_epochs.
    """
    s = active_pivot.copy()
    if drop_none:
        s = s.dropna()

    if s.empty:
        return pd.DataFrame(columns=["prn", "start", "end", "duration", "n_epochs"])

    # A new segment starts when PRN changes or when epochs are not consecutive.
    dt = s.index.to_series().diff()
    new_seg = (s != s.shift(1)) | (dt != sampling)
    seg_id = new_seg.cumsum()

    rows = []
    for _, g in s.groupby(seg_id):
        prn = str(g.iloc[0])
        start = g.index[0]
        end = g.index[-1]
        n = int(len(g))
        duration = (n - 1) * sampling  # span covered by consecutive epochs
        rows.append((prn, start, end, duration, n))

    return pd.DataFrame(rows, columns=["prn", "start", "end", "duration", "n_epochs"])


def check_full_coverage_from_active_pivot(active_pivot: pd.Series) -> dict:
    """
    Simple proof helper:
    - counts None epochs
    - returns first/last None epoch if any
    """
    n_none = int(active_pivot.isna().sum())
    none_epochs = active_pivot.index[active_pivot.isna()]
    return {
        "n_epochs": int(len(active_pivot)),
        "n_none": n_none,
        "first_none": none_epochs[0] if n_none else None,
        "last_none": none_epochs[-1] if n_none else None,
    }


# =============================================================================
# 4) Post-processing: remove too-short segments (merge prev/next, else fallback)
# =============================================================================

def fix_short_pivot_segments(
    df: pd.DataFrame,
    active_pivot: pd.Series,
    snr_col: str,
    snr_min: float,
    sampling: pd.Timedelta,
    min_duration: pd.Timedelta,
    pool: str = "selected_only",          # "selected_only" or "all"
    selected_prns: list[str] | None = None,
) -> pd.Series:
    """
    Remove segments shorter than min_duration.

    Priority
    --------
    1. merge into previous pivot if it covers the short interval
    2. merge into next pivot if it covers the short interval
    3. fallback: pick best PRN that covers the short interval (optional pool)

    Returns
    -------
    active_pivot_fixed : pd.Series
    """
    if pool not in ("selected_only", "all"):
        raise ValueError("pool must be 'selected_only' or 'all'.")
    if pool == "selected_only" and not selected_prns:
        raise ValueError("selected_prns must be provided when pool='selected_only'.")

    seg = pivot_schedule_to_segments(active_pivot, sampling=sampling, drop_none=True)
    if seg.empty:
        return active_pivot

    pool_prns = list(selected_prns) if pool == "selected_only" else \
        sorted(df.index.get_level_values("prn").unique())

    ap = active_pivot.copy()

    for i in range(len(seg)):
        if seg.loc[i, "duration"] >= min_duration:
            continue

        t0 = seg.loc[i, "start"]
        t1 = seg.loc[i, "end"]

        merged = False

        # Try previous
        if i > 0:
            prev_prn = seg.loc[i - 1, "prn"]
            if prn_covers_interval(df, prev_prn, t0, t1, snr_col, snr_min, sampling):
                ap.loc[t0:t1] = prev_prn
                merged = True

        # Try next
        if (not merged) and (i < len(seg) - 1):
            next_prn = seg.loc[i + 1, "prn"]
            if prn_covers_interval(df, next_prn, t0, t1, snr_col, snr_min, sampling):
                ap.loc[t0:t1] = next_prn
                merged = True

        # Fallback: best PRN in pool
        if not merged:
            fb = best_prn_for_interval(df, t0, t1, snr_col, snr_min, sampling, pool_prns)
            if fb is not None:
                ap.loc[t0:t1] = fb

    return ap


# =============================================================================
# 5) Small elegant stability improvement: minimum dwell time
# =============================================================================

def enforce_min_dwell(
    df: pd.DataFrame,
    active_pivot: pd.Series,
    snr_col: str,
    snr_min: float,
    sampling: pd.Timedelta,
    min_dwell: pd.Timedelta,
) -> pd.Series:
    """
    Stability post-filter:
    If a segment is shorter than min_dwell, try to absorb it into the previous pivot
    if the previous pivot can cover the short interval.

    This reduces last-minute tiny segments (e.g. 25 min) without harming coverage.
    """
    seg = pivot_schedule_to_segments(active_pivot, sampling=sampling, drop_none=True)
    if seg.empty:
        return active_pivot

    ap = active_pivot.copy()

    for i in range(len(seg)):
        if seg.loc[i, "duration"] >= min_dwell:
            continue
        if i == 0:
            continue

        t0 = seg.loc[i, "start"]
        t1 = seg.loc[i, "end"]
        prev_prn = seg.loc[i - 1, "prn"]

        if prn_covers_interval(df, prev_prn, t0, t1, snr_col, snr_min, sampling):
            ap.loc[t0:t1] = prev_prn

    return ap



