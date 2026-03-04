#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 16:39:34 2025

@author: snahmani
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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


# def enrich_df_with_sat_positions(df, mysp3):
#     """
#     Pour chaque ligne du DataFrame df (avec un index composé de (epoch, prn)),
#     calcule la position du satellite en prenant en compte le retard, l'effet
#     relativiste et met à jour le DataFrame avec les colonnes X_sat, Y_sat, Z_sat,
#     dte_sat et dRelat.

#     Parameters:
#       - df : pandas DataFrame dont l'index est (epoch, prn) et qui contient la colonne 'C1'
#       - mysp3 : instance de l'objet d'orbite contenant la méthode calcSatCoord
#       - t : instance de gpst.gpsdatetime() utilisée pour la conversion des temps

#     Returns:
#       - df enrichi avec les colonnes calculées.
#     """

#     # Listes pour stocker les résultats
#     x_sat  = []
#     y_sat  = []
#     z_sat  = []
#     dte_sat = []
#     d_relat = []

#     # Création de l'objet temps (pour la conversion)
#     t = gpst.gpsdatetime()

#     # Pour chaque observation du DataFrame
#     for (time_i, prn_i) in df.index:
#         # Conversion du temps d'observation au format attendu
#         # Par exemple, on formate time_i en chaîne de caractère
#         # t.rinex_t(time_i.to_pydatetime().strftime('%y %m %d %H %M %S.%f'))

#         # Calcul du temps d'émission initial (mjd)
#         t_emission_mjd = t.mjd - df.loc[(time_i, prn_i), 'C1'] / conv.SPEED_OF_LIGHT / 86400.0

#         # Calcul initial de la position du satellite
#         (x_sat_v, y_sat_v, z_sat_v, dte_sat_v) = mysp3.calcSatCoord(prn_i[0], int(prn_i[1:]), t_emission_mjd)

#         # Calcul de l'effet relativiste à partir d'une dérivée numérique
#         delta_t = 1e-3  # écart de temps en secondes
#         (xs1, ys1, zs1, _) = mysp3.calcSatCoord(prn_i[0], int(prn_i[1:]), t_emission_mjd - delta_t/86400.0)
#         (xs2, ys2, zs2, _) = mysp3.calcSatCoord(prn_i[0], int(prn_i[1:]), t_emission_mjd + delta_t/86400.0)

#         # Estimation de la vitesse par différence centrée
#         vx  = np.array([xs2 - xs1, ys2 - ys1, zs2 - zs1]) / (2.0 * delta_t)
#         vx0 = np.array([x_sat_v, y_sat_v, z_sat_v])

#         d_relat_v = -2.0 * vx0.T @ vx / (conv.SPEED_OF_LIGHT ** 2)

#         # Correction du temps d'émission tenant compte du retard d'horloge et de l'effet relativiste
#         t_emission_corr = t_emission_mjd - dte_sat_v / 86400.0 - d_relat_v / 86400.0

#         # Recalcul de la position au temps corrigé
#         (x_sat_v, y_sat_v, z_sat_v, dte_sat_v) = mysp3.calcSatCoord(prn_i[0], int(prn_i[1:]), t_emission_corr)

#         # Stockage des résultats
#         x_sat.append(x_sat_v)
#         y_sat.append(y_sat_v)
#         z_sat.append(z_sat_v)
#         dte_sat.append(dte_sat_v)
#         d_relat.append(d_relat_v)

#     # Ajout des nouvelles colonnes au DataFrame
#     df['X_sat']   = x_sat
#     df['Y_sat']   = y_sat
#     df['Z_sat']   = z_sat
#     df['dte_sat'] = dte_sat
#     df['dRelat']  = d_relat

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


# def add_az_el_iono_columns(df, P_rnx_header, mynav):
#     """
#     Ajoute dans le DataFrame les colonnes Az (azimut en degrés), Ele (élévation en degrés)
#     et d_ion1 (correction ionosphérique selon le modèle de Klobuchar).

#     Paramètres :
#       - df : DataFrame contenant les colonnes 'X_sat', 'Y_sat', 'Z_sat' et 'ind_ligne'.
#              L'index doit être un MultiIndex (time, prn).
#       - P_rnx_header : tableau ou liste contenant les coordonnées de la station (X, Y, Z)
#       - mynav : objet de navigation contenant les attributs ion_alpha_gps et ion_beta_gps.
#       - t : objet gpsdatetime (issu de gpst.gpsdatetime) pour la gestion des temps.
#       - tools : module (par ex. gnsstoolbox.gnsstools) fournissant les fonctions toolAzEle et toolCartGeoGRS80.
#       - klobuchar : module contenant la fonction klobuchar pour le calcul de la correction iono.

#     Retourne :
#       - df : le DataFrame enrichi avec les colonnes 'Az', 'Ele' et 'd_ion1'.
#     """

#     # Création de l'objet temps (pour la conversion)
#     t = gpst.gpsdatetime()

#     # Conversion des coordonnées de la station en géographiques
#     lon, lat, h = conv.xyz2geo(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2])
#     rad2deg = 180 / np.pi
#     lon_d = lon * rad2deg
#     lat_d = lat * rad2deg

#     # Calcul de l'azimut et de l'élévation en radians à partir des positions satellites
#     # az_rad, ele_rad = tools.toolAzEle(P_rnx_header[0],
#     #                                    P_rnx_header[1],
#     #                                    P_rnx_header[2],
#     #                                    df['X_sat'],
#     #                                    df['Y_sat'],
#     #                                    df['Z_sat'])

#     az_rad , ele_rad = conv.xyz2azi_ele(df['X_sat'], df['Y_sat'], df['Z_sat'],
#                                         P_rnx_header[0], P_rnx_header[1], P_rnx_header[2])

#     # Conversion en degrés
#     az_deg = az_rad * rad2deg
#     ele_deg = ele_rad * rad2deg

#     # Récupération des coefficients ionosphériques depuis le fichier de navigation
#     alpha = mynav.ion_alpha_gps
#     beta  = mynav.ion_beta_gps

#     # Calcul de la correction iono pour chaque observation
#     d_ion1 = []
#     for (time_i, prn_i) in df.index:
#         # Mise à jour de l'objet temps pour cette observation (conversion via rinex_t)
#         #t.rinex_t(time_i.to_pydatetime().strftime('%y %m %d %H %M %S.%f'))
#         #wsec_v = t.wsec  # secondes dans la semaine GPS

#         _, wsec_v = conv.dt2gpstime(time_i, secinweek=True)

#         # Récupération de l'indice de ligne correspondant (doit être un entier valide)
#         i = df.loc[(time_i, prn_i), 'ind_ligne']

#         # Calcul de la correction ionosphérique (la fonction klobuchar doit accepter ces paramètres)
#         d_ion1_v = klobuchar(lat_d, lon_d, ele_deg[i], az_deg[i], wsec_v, alpha, beta)
#         d_ion1.append(d_ion1_v)

#     # Ajout des nouvelles colonnes dans le DataFrame
#     df['Az']    = az_deg
#     df['Ele']   = ele_deg
#     df['d_ion1'] = d_ion1

#     return df



#%% 
# -------------------------------------------------------------------------
# Plot Single Differences (SD) by satellite PRN (with gap handling)
#
# Input expected:
#   df_SD : DataFrame indexed by (epoch, prn) with columns like "SD_L1", "SD_C1", ...
# -------------------------------------------------------------------------

def plot_gnss_sd_by_prn(
    df_SD,
    observable="SD_L1",
    gap=pd.Timedelta(minutes=30),
    label_arcs=True,
    arc_label_fontsize=8,
    arc_label_xytext=(3, 3),   # offset in points
    show_legend=False,
    legend_outside=True,
    legend_ncol=4
):
    """
    Plot GNSS single differences time series by satellite PRN, with gap handling.
    Optionally write PRN labels directly on the plot at the beginning of each arc.

    Parameters
    ----------
    df_SD : pandas.DataFrame
        MultiIndex (epoch, prn). Columns include SD observables (e.g., "SD_L1").
    observable : str
        Column name to plot.
    gap : pandas.Timedelta
        Time gap threshold to break line segments.
    label_arcs : bool
        If True, place PRN text at the start of each continuous segment ("arc").
    arc_label_fontsize : int
        Font size for PRN labels.
    arc_label_xytext : tuple
        (dx, dy) offset in points for the text label, to avoid overlapping the line.
    show_legend : bool
        If True, show legend. Often unnecessary if label_arcs=True.
    legend_outside : bool
        If True, place legend outside the axes (right side).
    legend_ncol : int
        Number of columns in legend.

    Returns
    -------
    fig, ax
    """

    if not isinstance(df_SD.index, pd.MultiIndex):
        raise ValueError("df_SD must have a MultiIndex (epoch, prn).")

    if observable not in df_SD.columns:
        raise ValueError(f"{observable} not found in df_SD.")

    if "epoch" not in df_SD.index.names or "prn" not in df_SD.index.names:
        raise ValueError("df_SD index must have levels ('epoch','prn').")

    prns = df_SD.index.get_level_values("prn").unique()

    fig, ax = plt.subplots(figsize=(10, 6))

    for prn in prns:
        # Extract one PRN time series (index becomes epochs)
        data = df_SD.xs(prn, level="prn").copy().sort_index()

        # Skip if nothing to plot
        if not data[observable].notna().any():
            continue

        # Remove NaNs for segmentation + plotting
        data = data.dropna(subset=[observable])
        if data.empty:
            continue

        # Segment id based on gaps
        dt = data.index.to_series().diff()
        seg_id = (dt > gap).cumsum()

        color = None

        for _, seg in data.groupby(seg_id):
            if seg.empty:
                continue

            # Plot this segment
            line, = ax.plot(
                seg.index,
                seg[observable].to_numpy(),
                color=color,
                label=prn if (color is None and show_legend) else None
            )

            if color is None:
                color = line.get_color()

            # Label at the beginning of each arc (segment)
            if label_arcs:
                x0 = seg.index[0]
                y0 = float(seg[observable].iloc[0])

                ax.annotate(
                    prn,
                    xy=(x0, y0),
                    xytext=arc_label_xytext,
                    textcoords="offset points",
                    fontsize=arc_label_fontsize,
                    color=color,
                    ha="left",
                    va="bottom"
                )

    ax.set_title(f"Single differences time series ({observable}) by satellite PRN")
    ax.set_xlabel("Time (epoch)")
    ax.set_ylabel("Value (cycles for SD_L*, meters for SD_C*/SD_P*)")

    if show_legend:
        if legend_outside:
            ax.legend(
                title="PRN",
                bbox_to_anchor=(1.02, 1),
                loc="upper left",
                borderaxespad=0,
                ncol=legend_ncol
            )
            plt.tight_layout(rect=[0, 0, 0.85, 1])
        else:
            ax.legend(title="PRN", ncol=legend_ncol)
            plt.tight_layout()
    else:
        plt.tight_layout()

    plt.show()
    return fig, ax






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

    Paramètres :
        df      : DataFrame avec un index multi-niveaux contenant au moins le niveau 'prn'
                  et les colonnes col1 (et éventuellement col2).
        col1    : Nom de la première colonne.
        col2    : Nom de la seconde colonne (optionnel). Si None, on affiche col1.
        coeff1  : Coefficient multiplicateur pour la première colonne (défaut 1.0).
        coeff2  : Coefficient multiplicateur pour la seconde colonne (défaut 1.0).
        seuil   : Seuil en secondes pour considérer un "trou" dans la série (défaut 3600).
        renderer: Renderer Plotly (ex : "browser" ou "iframe").

    Retourne :
        fig     : Figure Plotly contenant les courbes tracées.
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


def plot_residual_analysis(A, B, dP_est, figure_title=None, save_path=None,
                           P_est=None, P_rnx_header=None):
    """
    Computes residuals (v_est = B - A @ dP_est) and creates a figure containing:
      1. Time series of residuals (displayed as points)
      2. Histogram of residuals (number of observations per bin)
      3. Q-Q Plot of residuals
      4. Scatter plot of residuals vs predicted values
      5. Text box displaying statistics (mean, variance, std dev, skewness, kurtosis)
      6. (Optional) Text box with additional position information:
         - Distance between estimated position and initial RINEX header position
         - Local ENU coordinates computed via tools.toolCartLocGRS80

    Parameters:
      - A : array-like, design matrix.
      - B : array-like, observation vector.
      - dP_est : array-like, estimated parameter vector.
      - figure_title (optional) : str, overall figure title.
      - save_path (optional) : str, full path (name + extension) to save the figure.
      - P_est (optional) : array-like, estimated position (for additional info computation).
      - P_rnx_header (optional) : array-like, initial position from RINEX header.


    Returns:
      - fig : matplotlib Figure object containing all subplots.
    """

    # Compute residuals and predicted values
    v_est = B - A @ dP_est
    b_est = A @ dP_est

    # Compute residual statistics
    mean_val   = np.mean(v_est)
    variance   = np.var(v_est)
    std_dev    = np.std(v_est)
    skewness   = stats.skew(v_est)
    kurtosis   = stats.kurtosis(v_est)

    # Create figure with GridSpec
    # Using 4 rows and 2 columns:
    # - Row 0: Time series (col 0) and Histogram (col 1)
    # - Row 1: Q-Q Plot (col 0) and Residuals vs. Predicted values (col 1)
    # - Row 2: Statistics text box (spanning 2 columns)
    # - Row 3: (Optional) Position information text box (spanning 2 columns)
    nrows = 4 if (P_est is not None and P_rnx_header is not None) else 3
    height_ratios = [1, 1, 0.5, 0.5] if nrows == 4 else [1, 1, 0.5]

    fig = plt.figure(figsize=(15, 15))
    gs = fig.add_gridspec(nrows, 2, height_ratios=height_ratios)

    # 1. Time series of residuals (points only) – Top left
    ax_time = fig.add_subplot(gs[0, 0])
    ax_time.scatter(np.arange(len(v_est)), v_est, color='green')
    ax_time.set_title("Time Series of Residuals")
    ax_time.set_xlabel("Time / Index")
    ax_time.set_ylabel("Residuals (m)")

    # 2. Histogram of residuals (raw observation count) – Top right
    ax_hist = fig.add_subplot(gs[0, 1])
    sns.histplot(v_est, bins=30, stat="count", color='skyblue', edgecolor='black', ax=ax_hist)
    ax_hist.set_title("Histogram of Residuals")
    ax_hist.set_xlabel("Residuals")
    ax_hist.set_ylabel("Number of Observations")

    # 3. Q-Q Plot of residuals – Row 1, Column 0
    ax_qq = fig.add_subplot(gs[1, 0])
    sm.qqplot(v_est, line='s', ax=ax_qq)
    ax_qq.set_title("Q-Q Plot of Residuals")

    # 4. Residuals vs. predicted values plot – Row 1, Column 1
    ax_scatter = fig.add_subplot(gs[1, 1])
    ax_scatter.scatter(b_est, v_est, alpha=0.7, color='darkorange')
    ax_scatter.axhline(0, color='red', linestyle='--')
    ax_scatter.set_xlabel("Predicted Values")
    ax_scatter.set_ylabel("Residuals (m)")
    ax_scatter.set_title("Residuals vs. Predicted Values")

    # 5. Statistics text box – Next row (spans 2 columns)
    ax_stats = fig.add_subplot(gs[2, :])
    ax_stats.axis('off')  # Hide axes
    stats_text = (
        f"Mean       : {mean_val:.4f}\n"
        f"Variance   : {variance:.4f}\n"
        f"Std Dev    : {std_dev:.4f}\n"
        f"Skewness   : {skewness:.4f}\n"
        f"Kurtosis   : {kurtosis:.4f}"
    )
    ax_stats.text(0.5, 0.5, stats_text, transform=ax_stats.transAxes,
                  fontsize=14, verticalalignment='center', horizontalalignment='center',
                  bbox=dict(facecolor='wheat', edgecolor='black', boxstyle='round,pad=1'))
    ax_stats.set_title("Residual Statistics", fontsize=16)

    # 6. Position information text box – (optional)
    if nrows == 4:
        # Compute distance between estimated position and RINEX header position
        dist = np.sqrt(np.sum((P_est - P_rnx_header)**2))
        # Compute local ENU coordinates via toolCartLocGRS80 function
        E, N, U = conv.xyz2enu(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],
                                          P_est[0], P_est[1], P_est[2])
        extra_text = (
            f"Distance between estimated position and initial RINEX header position: {dist:.4f} m\n"
            f"Local ENU coordinates:\n"
            f"  East (E)  : {E[0]:.4f} m"
            f"  North (N) : {N[0]:.4f} m"
            f"  Up (U)    : {U[0]:.4f} m"
        )
        ax_extra = fig.add_subplot(gs[3, :])
        ax_extra.axis('off')
        ax_extra.text(0.5, 0.5, extra_text, transform=ax_extra.transAxes,
                      fontsize=14, verticalalignment='center', horizontalalignment='center',
                      bbox=dict(facecolor='lightcyan', edgecolor='black', boxstyle='round,pad=1'))
        ax_extra.set_title("Position Information", fontsize=16)

    # Overall figure title if provided
    if figure_title is not None:
        fig.suptitle(figure_title, fontsize=20)
        plt.subplots_adjust(top=0.92)

    # Adjust layout
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Save figure if path is provided
    if save_path is not None:
        fig.savefig(save_path)
        print(f"Figure saved to: {save_path}")

    plt.show()
    return fig



