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
import geodezyx.files_rw as files_rw  # Import the file reading module
import geodezyx.conv as conv                  # Import the conversion module
import datetime as dt


import gnsstoolbox.orbits as orb
import gpsdatetime as gpst
import gnsstoolbox.gnss_const as gnss_const

from .klobuchar import *
import gnsstoolbox.gnsstools as tools

# Pour chercher des chaînes de caractère dans des fichiers
import re


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


def enrich_df_with_sat_positions(df, mysp3):
    """
    Pour chaque ligne du DataFrame df (avec un index composé de (epoch, prn)),
    calcule la position du satellite en prenant en compte le retard, l'effet
    relativiste et met à jour le DataFrame avec les colonnes X_sat, Y_sat, Z_sat,
    dte_sat et dRelat.

    Parameters:
      - df : pandas DataFrame dont l'index est (epoch, prn) et qui contient la colonne 'C1'
      - mysp3 : instance de l'objet d'orbite contenant la méthode calcSatCoord
      - t : instance de gpst.gpsdatetime() utilisée pour la conversion des temps

    Returns:
      - df enrichi avec les colonnes calculées.
    """

    # Listes pour stocker les résultats
    x_sat  = []
    y_sat  = []
    z_sat  = []
    dte_sat = []
    d_relat = []

    # Création de l'objet temps (pour la conversion)
    t = gpst.gpsdatetime()

    # Pour chaque observation du DataFrame
    for (time_i, prn_i) in df.index:
        # Conversion du temps d'observation au format attendu
        # Par exemple, on formate time_i en chaîne de caractère
        # t.rinex_t(time_i.to_pydatetime().strftime('%y %m %d %H %M %S.%f'))

        # Calcul du temps d'émission initial (mjd)
        t_emission_mjd = t.mjd - df.loc[(time_i, prn_i), 'C1'] / conv.SPEED_OF_LIGHT / 86400.0

        # Calcul initial de la position du satellite
        (x_sat_v, y_sat_v, z_sat_v, dte_sat_v) = mysp3.calcSatCoord(prn_i[0], int(prn_i[1:]), t_emission_mjd)

        # Calcul de l'effet relativiste à partir d'une dérivée numérique
        delta_t = 1e-3  # écart de temps en secondes
        (xs1, ys1, zs1, _) = mysp3.calcSatCoord(prn_i[0], int(prn_i[1:]), t_emission_mjd - delta_t/86400.0)
        (xs2, ys2, zs2, _) = mysp3.calcSatCoord(prn_i[0], int(prn_i[1:]), t_emission_mjd + delta_t/86400.0)

        # Estimation de la vitesse par différence centrée
        vx  = np.array([xs2 - xs1, ys2 - ys1, zs2 - zs1]) / (2.0 * delta_t)
        vx0 = np.array([x_sat_v, y_sat_v, z_sat_v])

        d_relat_v = -2.0 * vx0.T @ vx / (conv.SPEED_OF_LIGHT ** 2)

        # Correction du temps d'émission tenant compte du retard d'horloge et de l'effet relativiste
        t_emission_corr = t_emission_mjd - dte_sat_v / 86400.0 - d_relat_v / 86400.0

        # Recalcul de la position au temps corrigé
        (x_sat_v, y_sat_v, z_sat_v, dte_sat_v) = mysp3.calcSatCoord(prn_i[0], int(prn_i[1:]), t_emission_corr)

        # Stockage des résultats
        x_sat.append(x_sat_v)
        y_sat.append(y_sat_v)
        z_sat.append(z_sat_v)
        dte_sat.append(dte_sat_v)
        d_relat.append(d_relat_v)

    # Ajout des nouvelles colonnes au DataFrame
    df['X_sat']   = x_sat
    df['Y_sat']   = y_sat
    df['Z_sat']   = z_sat
    df['dte_sat'] = dte_sat
    df['dRelat']  = d_relat

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


def add_az_el_iono_columns(df, P_rnx_header, mynav):
    """
    Ajoute dans le DataFrame les colonnes Az (azimut en degrés), Ele (élévation en degrés)
    et d_ion1 (correction ionosphérique selon le modèle de Klobuchar).

    Paramètres :
      - df : DataFrame contenant les colonnes 'X_sat', 'Y_sat', 'Z_sat' et 'ind_ligne'.
             L'index doit être un MultiIndex (time, prn).
      - P_rnx_header : tableau ou liste contenant les coordonnées de la station (X, Y, Z)
      - mynav : objet de navigation contenant les attributs ion_alpha_gps et ion_beta_gps.
      - t : objet gpsdatetime (issu de gpst.gpsdatetime) pour la gestion des temps.
      - tools : module (par ex. gnsstoolbox.gnsstools) fournissant les fonctions toolAzEle et toolCartGeoGRS80.
      - klobuchar : module contenant la fonction klobuchar pour le calcul de la correction iono.

    Retourne :
      - df : le DataFrame enrichi avec les colonnes 'Az', 'Ele' et 'd_ion1'.
    """

    # Création de l'objet temps (pour la conversion)
    t = gpst.gpsdatetime()

    # Conversion des coordonnées de la station en géographiques
    lon, lat, h = conv.xyz2geo(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2])
    rad2deg = 180 / np.pi
    lon_d = lon * rad2deg
    lat_d = lat * rad2deg

    # Calcul de l'azimut et de l'élévation en radians à partir des positions satellites
    # az_rad, ele_rad = tools.toolAzEle(P_rnx_header[0],
    #                                    P_rnx_header[1],
    #                                    P_rnx_header[2],
    #                                    df['X_sat'],
    #                                    df['Y_sat'],
    #                                    df['Z_sat'])

    az_rad , ele_rad = conv.xyz2azi_ele(df['X_sat'], df['Y_sat'], df['Z_sat'],
                                        P_rnx_header[0], P_rnx_header[1], P_rnx_header[2])

    # Conversion en degrés
    az_deg = az_rad * rad2deg
    ele_deg = ele_rad * rad2deg

    # Récupération des coefficients ionosphériques depuis le fichier de navigation
    alpha = mynav.ion_alpha_gps
    beta  = mynav.ion_beta_gps

    # Calcul de la correction iono pour chaque observation
    d_ion1 = []
    for (time_i, prn_i) in df.index:
        # Mise à jour de l'objet temps pour cette observation (conversion via rinex_t)
        #t.rinex_t(time_i.to_pydatetime().strftime('%y %m %d %H %M %S.%f'))
        #wsec_v = t.wsec  # secondes dans la semaine GPS

        _, wsec_v = conv.dt2gpstime(time_i, secinweek=True)

        # Récupération de l'indice de ligne correspondant (doit être un entier valide)
        i = df.loc[(time_i, prn_i), 'ind_ligne']

        # Calcul de la correction ionosphérique (la fonction klobuchar doit accepter ces paramètres)
        d_ion1_v = klobuchar.klobuchar(lat_d, lon_d, ele_deg[i], az_deg[i], wsec_v, alpha, beta)
        d_ion1.append(d_ion1_v)

    # Ajout des nouvelles colonnes dans le DataFrame
    df['Az']    = az_deg
    df['Ele']   = ele_deg
    df['d_ion1'] = d_ion1

    return df


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
                           P_est=None, P_rnx_header=None, tools=None):
    """
    Calcule les résidus (v_est = B - A @ dP_est) et trace une figure contenant :
      1. La série temporelle des résidus (affichée en points)
      2. L'histogramme des résidus (nombre d'observations par bin)
      3. Le Q-Q Plot des résidus
      4. Un scatter plot des résidus en fonction des valeurs prédites
      5. Un encadré affichant quelques statistiques (moyenne, variance, écart-type, skewness, kurtosis)
      6. (Optionnel) Un encadré avec des informations supplémentaires sur la position :
         - Distance entre la position estimée et la position initiale du header RINEX
         - Coordonnées ENU locales calculées via tools.toolCartLocGRS80
         
    Paramètres :
      - A : array-like, matrice des variables explicatives.
      - B : array-like, vecteur des observations.
      - dP_est : array-like, vecteur des paramètres estimés.
      - figure_title (optionnel) : str, titre global de la figure.
      - save_path (optionnel) : str, chemin complet (nom + extension) pour sauvegarder la figure.
      - P_est (optionnel) : array-like, position estimée (pour le calcul des informations supplémentaires).
      - P_rnx_header (optionnel) : array-like, position initiale extraite du header RINEX.
      - tools (optionnel) : module ou objet possédant la fonction toolCartLocGRS80.
      
    Renvoie :
      - fig : l'objet Figure de matplotlib contenant l'ensemble des graphiques.
    """
    
    # Calcul des résidus et des valeurs prédites
    v_est = B - A @ dP_est
    b_est = A @ dP_est
    
    # Calcul des statistiques sur les résidus
    moyenne    = np.mean(v_est)
    variance   = np.var(v_est)
    ecart_type = np.std(v_est)
    skewness   = stats.skew(v_est)
    kurtosis   = stats.kurtosis(v_est)
    
    # Création de la figure avec GridSpec.
    # On utilise 4 lignes et 2 colonnes :
    # - Ligne 0 : Série temporelle (col 0) et Histogramme (col 1)
    # - Ligne 1 : Q-Q Plot (col 0) et Résidus vs. Valeurs prédites (col 1)
    # - Ligne 2 : Encadré des statistiques (s'étendant sur 2 colonnes)
    # - Ligne 3 : (Optionnel) Encadré des informations de position (s'étendant sur 2 colonnes)
    nrows = 4 if (P_est is not None and P_rnx_header is not None and tools is not None) else 3
    height_ratios = [1, 1, 0.5, 0.5] if nrows == 4 else [1, 1, 0.5]
    
    fig = plt.figure(figsize=(15, 15))
    gs = fig.add_gridspec(nrows, 2, height_ratios=height_ratios)
    
    # 1. Série temporelle des résidus (affichage uniquement des points) – Top gauche
    ax_time = fig.add_subplot(gs[0, 0])
    ax_time.scatter(np.arange(len(v_est)), v_est, color='green')
    ax_time.set_title("Série temporelle des résidus")
    ax_time.set_xlabel("Temps / Index")
    ax_time.set_ylabel("Résidus")
    
    # 2. Histogramme des résidus (nombre brut d'observations) – Top droite
    ax_hist = fig.add_subplot(gs[0, 1])
    sns.histplot(v_est, bins=30, stat="count", color='skyblue', edgecolor='black', ax=ax_hist)
    ax_hist.set_title("Histogramme des résidus")
    ax_hist.set_xlabel("Résidus")
    ax_hist.set_ylabel("Nombre d'observations")
    
    # 3. Q-Q Plot des résidus – Ligne 1, Colonne 0
    ax_qq = fig.add_subplot(gs[1, 0])
    sm.qqplot(v_est, line='s', ax=ax_qq)
    ax_qq.set_title("Q-Q Plot des résidus")
    
    # 4. Graphique des résidus vs. valeurs prédites – Ligne 1, Colonne 1
    ax_scatter = fig.add_subplot(gs[1, 1])
    ax_scatter.scatter(b_est, v_est, alpha=0.7, color='darkorange')
    ax_scatter.axhline(0, color='red', linestyle='--')
    ax_scatter.set_xlabel("Valeurs prédites")
    ax_scatter.set_ylabel("Résidus")
    ax_scatter.set_title("Résidus vs. Valeurs prédites")
    
    # 5. Encadré avec les statistiques – Ligne suivante (s'étend sur 2 colonnes)
    ax_stats = fig.add_subplot(gs[2, :])
    ax_stats.axis('off')  # Masquer les axes
    texte_stats = (
        f"Moyenne    : {moyenne:.4f}\n"
        f"Variance   : {variance:.4f}\n"
        f"Écart-type : {ecart_type:.4f}\n"
        f"Skewness   : {skewness:.4f}\n"
        f"Kurtosis   : {kurtosis:.4f}"
    )
    ax_stats.text(0.5, 0.5, texte_stats, transform=ax_stats.transAxes,
                  fontsize=14, verticalalignment='center', horizontalalignment='center',
                  bbox=dict(facecolor='wheat', edgecolor='black', boxstyle='round,pad=1'))
    ax_stats.set_title("Statistiques des résidus", fontsize=16)
    
    # 6. Encadré avec les informations de position – (optionnel)
    if nrows == 4:
        # Calcul de la distance entre la position estimée et la position du header RINEX
        dist = np.sqrt(np.sum((P_est - P_rnx_header)**2))
        # Calcul des coordonnées ENU locales via la fonction toolCartLocGRS80
        E, N, U = tools.toolCartLocGRS80(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],
                                          P_est[0], P_est[1], P_est[2])
        texte_extra = (
            f"Distance entre la position estimée et la position initiale du header RINEX : {dist:.4f}\n"
            f"Coordonnées ENU locales :\n"
            f"  Est (E) : {E:.4f}\n"
            f"  Nord (N): {N:.4f}\n"
            f"  Haut (U) : {U:.4f}"
        )
        ax_extra = fig.add_subplot(gs[3, :])
        ax_extra.axis('off')
        ax_extra.text(0.5, 0.5, texte_extra, transform=ax_extra.transAxes,
                      fontsize=14, verticalalignment='center', horizontalalignment='center',
                      bbox=dict(facecolor='lightcyan', edgecolor='black', boxstyle='round,pad=1'))
        ax_extra.set_title("Informations de position", fontsize=16)
    
    # Titre global de la figure si fourni
    if figure_title is not None:
        fig.suptitle(figure_title, fontsize=20)
        plt.subplots_adjust(top=0.92)
    
    # Ajuster la mise en page
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    # Sauvegarder la figure si un chemin est fourni
    if save_path is not None:
        fig.savefig(save_path)
        print(f"Figure sauvegardée dans : {save_path}")
    
    plt.show()
    return fig
