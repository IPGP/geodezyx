#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filière ING3 - PPMD - Traitement de la mesure de code (GNSS)

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
# Chargement des fichiers RINEX d'observation
fichier_rnx = dwl_output_station[0][0]

# Chargement des données RINEX d'observation dans un pandas dataframe via  GeodeZYX
df_rnx_orig, l_rnx_head = files_rw.read_rinex_obs(fichier_rnx,
                                           return_header=True)
df_rnx = df_rnx_orig.copy()


# Position approchée lue dans le header du fichier RINEX
columns = [l for l in l_rnx_head if r"APPROX POSITION XYZ" in l]
if columns:
    # Séparer la ligne en supprimant le texte "APPROX POSITION XYZ"
    valeurs = columns[0].split()[:3]  # Récupérer uniquement les 3 premières valeurs
    P_rnx_header = np.array(valeurs, dtype=float)  # Convertir en numpy array
    print("Coordonnées XYZ :", P_rnx_header)
else:
    P_rnx_header = np.array([0,0,0])
    print("Aucune correspondance trouvée.")
del columns, valeurs


# il faut consulter la doc et vous constaterez que les unités de L1 et L2 sont en cycle
# on multiplie les valeurs de L1 et L2 par les longueurs d'onde correspondantes pour avoir
# des distances (ambigues mais en mètre).

df_rnx['L1'] = df_rnx['L1']*conv.L1_WAVELENGTH
df_rnx['L2'] = df_rnx['L2']*conv.L2_WAVELENGTH
df_rnx['L5'] = df_rnx['L5']*conv.L5_WAVELENGTH


# nettoyage cf fin du TP01
df_rnx = df_rnx.dropna(axis=1, how='all')
rows_with_nan = df_rnx[['C1','L1', 'L2']].isna().any(axis=1)
df_removed = df_rnx[rows_with_nan]
df_rnx = df_rnx.dropna(subset=['C1','L1', 'L2'])
del rows_with_nan

# on ne garde que les mesures des satellites GPS
# Filtrer pour garder uniquement les lignes où la colonne 'sys' contient 'G'
df_filtre = df_rnx[df_rnx['sys'].str.contains('G')]

df_rnx = df_filtre
del df_filtre

df_rnx = df_rnx.dropna(axis=1, how='all')

# ajout de l'indice de ligne
df_rnx['ind_ligne'] = range(len(df_rnx)) 


# %%
# Chargement du fichier sp3 dans un dataframe
fichier_sp3 = dwl_output_satellite[0]

print(f"Using SP3 file: {fichier_sp3}")

df_sp3 = files_rw.read_sp3(fichier_sp3)

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
    "df_removed",
    "df_rnx_orig",
    "dwl_output_station",
    "dwl_output_satellite",
    "dwl_output_navigation",
    "fichier_sp3",
    "l_rnx_head",
]:
    if var in globals():
        del globals()[var]

del var
gc.collect()

# %%
# Satellite position computation at signal emission time (GNSS observation model)
#
# This block computes satellite positions and clock corrections at the true
# signal emission epoch for each GNSS observation contained in the RINEX
# dataframe.
#
# Methodology
# -----------
# 1. For each satellite (PRN), the signal travel time is estimated from the
#    pseudorange observation (C1):
#
#           τ ≈ ρ / c
#
#    where ρ is the pseudorange and c the speed of light.
#
# 2. An approximate emission time is derived:
#
#           t_emit ≈ t_receive − τ
#
# 3. Satellite positions are interpolated from precise SP3 orbits using
#    Lagrange interpolation (GeodeZYX orbital interpolation tools).
#
# 4. The relativistic correction is computed using:
#
#           dRel = -2 (r · v) / c²
#
#    where satellite velocity is numerically estimated by centered finite
#    differences around the emission epoch.
#
# 5. Satellite clock correction contained in the SP3 file is applied to refine
#    the emission epoch.
#
# 6. Final satellite coordinates (ECEF, meters), clock offsets (seconds),
#    and relativistic corrections are stored both:
#       - in dedicated output arrays
#       - directly inside the observation dataframe.
#
# Notes
# -----
# - SP3 positions are provided in kilometers and converted here to meters.
# - Clock corrections are converted from microseconds to seconds.
# - Computations are performed independently for each satellite to ensure
#   consistency of orbital interpolation.
#
# Output columns added to df_rnx
# ------------------------------
#   X_sat   : satellite X coordinate (m, ECEF)
#   Y_sat   : satellite Y coordinate (m, ECEF)
#   Z_sat   : satellite Z coordinate (m, ECEF)
#   dte_sat : satellite clock correction (s)
#   dRelat  : relativistic correction (s)
#
###############################################################################

X_sat = []
Y_sat = []
Z_sat = []
dte_sat = []
dRelat = []

#### NEW GEODEZYX STYLE : on utilise les fonctions d'interpolation de GeodeZYX pour calculer la position des satellites à chaque temps d'émission
for prn in df_rnx['prn'].unique():
    df_rnx_prn = df_rnx[df_rnx['prn'] == prn]
    if not prn in df_sp3['prn'].unique():
        print(f"Satellite {prn} non présent dans le fichier SP3, on ne peut pas calculer sa position à chaque temps d'émission")
        continue
    df_sp3_prn = df_sp3[df_sp3['prn'] == prn]
    fly_time = df_rnx_prn['C1'] / conv.SPEED_OF_LIGHT
    fly_time = pd.to_timedelta(fly_time, unit="s")
    t_rec = df_rnx_prn['epoch']
    t_emi_approx = t_rec - fly_time
    orb_df_approx = reffram.orb_df_lagrange_interpolate(df_sp3_prn, t_emi_approx.values)
    # les positions dans les fichiers SP3 sont en km, on les convertit en mètre pour être cohérent avec les unités de la mesure de pseudodistance
    orb_df_approx[["x","y","z"]] = orb_df_approx[["x","y","z"]] * 10**3 # km -> m

    # calcul de l'effet relativiste
    delta_t = pd.to_timedelta(1e-3, unit="s") # écart de temps en +/- pour calculer la dérivée
    orb_df_rel_fwd = reffram.orb_df_lagrange_interpolate(df_sp3_prn, t_emi_approx.values + delta_t)
    orb_df_rel_fwd[["x","y","z"]] = orb_df_rel_fwd[["x","y","z"]] * 10**3 # km -> m
    orb_df_rel_bak = reffram.orb_df_lagrange_interpolate(df_sp3_prn, t_emi_approx.values - delta_t)
    orb_df_rel_bak[["x","y","z"]] = orb_df_rel_bak[["x","y","z"]] * 10**3 # km -> m

    xyz_diff = (orb_df_rel_fwd[["x","y","z"]] - orb_df_rel_bak[["x","y","z"]]) / (2 * delta_t.total_seconds())
    xyz = orb_df_approx[["x","y","z"]]
    dRelat_v = -2.0 * (xyz_diff * xyz).sum(axis=1) / (conv.SPEED_OF_LIGHT **2)

    t_emi_ok = t_emi_approx -  pd.to_timedelta(orb_df_approx["clk"].values, unit="us") # - pd.to_timedelta(dRelat_v, unit="s")
    orb_df_ok = reffram.orb_df_lagrange_interpolate(df_sp3_prn, t_emi_ok.values)
    orb_df_ok[["x","y","z"]] = orb_df_ok[["x","y","z"]] * 10**3 # km -> m

    X_sat.extend(orb_df_ok["x"].values)
    Y_sat.extend(orb_df_ok["y"].values)
    Z_sat.extend(orb_df_ok["z"].values)
    dte_sat.extend(orb_df_ok["clk"].values * 1e-6) # microsecondes -> secondes
    dRelat.extend(dRelat_v.values)

    df_rnx.loc[df_rnx_prn.index, 'X_sat'] = orb_df_ok["x"].values
    df_rnx.loc[df_rnx_prn.index, 'Y_sat'] = orb_df_ok["y"].values
    df_rnx.loc[df_rnx_prn.index, 'Z_sat'] = orb_df_ok["z"].values
    df_rnx.loc[df_rnx_prn.index, 'dte_sat'] = orb_df_ok["clk"].values * 1e-6 # microsecondes -> secondes
    df_rnx.loc[df_rnx_prn.index, 'dRelat'] = dRelat_v.values

df_rinex_flat = df_rnx.copy()
df_rnx = df_rnx.set_index(['epoch', 'prn'], drop=True)

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
    "X_sat",
    "Y_sat",
    "Z_sat",
    "xyz",
    "xyz_diff",
    "t_emi_approx",
    "delta_t",
    "df_rinex_flat",
    "df_rnx_prn",
    "df_sp3_prn",
    "df_sp3",
    "dte_sat",
    "orb_df_approx",
    "orb_df_ok",
    "orb_df_rel_bak",
    "orb_df_rel_fwd",
    "dRelat",
    "dRelat_v",
    "fly_time",
    "prn",
    "t_emi_ok",
    "t_rec",
]:
    if var in globals():
        del globals()[var]

del var
gc.collect()


# %%
# Traitement classique sur le code pour la position d'un récepteur GNSS
# cas trivial

print('*****  Calcul trivial *****')


# Initialisation des coordonnées approximatives du récepteur
P_app = np.array([0, 0, 0])

# Itération pour l'affinement de la position du récepteur
dP_est=np.array([100, 100, 100]) # initialisation : juste pour entrer dans la boucle while ci-après
i=1
while np.linalg.norm(dP_est)>1:
    # Calcul des distances approximatives satellite-récepteur
    distances = np.sqrt((df_rnx['X_sat'].values - P_app[0])**2 +
                        (df_rnx['Y_sat'].values - P_app[1])**2 +
                        (df_rnx['Z_sat'].values - P_app[2])**2)

    # Vecteur des observations corrigées des différents modèles
    B = df_rnx['C1'].values - distances

    # Construction de la matrice des dérivées partielles
    df_dX = (P_app[0] - df_rnx['X_sat'].values) / distances
    df_dY = (P_app[1] - df_rnx['Y_sat'].values) / distances
    df_dZ = (P_app[2] - df_rnx['Z_sat'].values) / distances
    A = np.column_stack((df_dX, df_dY, df_dZ))

    # Résolution par moindres carrés pour estimer le déplacement
    # formule vue en cours mais qui est plus lente 
    #dP_est = np.linalg.inv(A.T@A)@A.T@B
    dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
    
    # Mise à jour de la position estimée
    P_est = P_app + dP_est
    
    # Affichage de la position estimée à chaque itération
    print(f"Iteration {i}: Position estimée - X: {P_est[0]}, Y: {P_est[1]}, Z: {P_est[2]}")
    P_app = P_est  # Mise à jour de la position approximative pour la prochaine itération
    i+=1

    # Calcul de la distance finale entre la position estimée et la position initiale du header RINEX
    dist_P_est_P_rnx_header = np.sqrt(np.sum((P_est - P_rnx_header)**2))
    print("\n")
    print("Distance entre la position estimée et la position initiale du header RINEX:", dist_P_est_P_rnx_header)
    
    # Calculer les coordonnées ENU locales
    E, N, U = conv.xyz2enu(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
    print("Est (E):", E)
    print("North (N):", N)
    print("Up (U):", U)
    print("\n")


fig = gnss_edu.plot_residual_analysis(A, B, dP_est, figure_title="Trivial GPS data processing", 
                                      save_path = folder / "01_trivial.png",
                                      P_est=P_est, P_rnx_header=P_rnx_header)

# -------------------------------------------------------------------------
# Store results for later reuse
# -------------------------------------------------------------------------
S01_trivial_solution = {
    "P_est": P_est.copy(),
    "ENU": np.array([E, N, U]),
    "distance_to_header": dist_P_est_P_rnx_header,
    "iterations": i-1
}

#%%
# Numerical comparison of least-squares solution methods
#
# Educational objective
# ---------------------
# Compare different numerical implementations of the least-squares solution:
#
#       dP = (AᵀA)⁻¹ AᵀB
#
# Although this analytical formula is commonly introduced in lectures,
# it is rarely implemented directly in scientific software.
#
# Three approaches are compared:
#
# 1) Explicit matrix inverse (theoretical formulation)
#       -> simple but numerically inefficient and discouraged
#
# 2) Normal equation solved with np.linalg.solve
#       -> faster and numerically safer
#
# 3) np.linalg.lstsq (QR/SVD decomposition)
#       -> most robust method, used in professional GNSS processing
#
# The computation is repeated many times to make execution time differences
# clearly visible.
###############################################################################

import time

ATA = A.T @ A
ATB = A.T @ B

N_REPEAT = 5000

# -------------------------------------------------------------------------
# 1) Explicit inverse (shown for pedagogical purposes only)
# -------------------------------------------------------------------------
t0 = time.perf_counter()
for _ in range(N_REPEAT):
    x1 = np.linalg.inv(ATA) @ ATB
t1 = time.perf_counter()

# -------------------------------------------------------------------------
# 2) Solve normal equations (recommended fast approach)
# -------------------------------------------------------------------------
t2 = time.perf_counter()
for _ in range(N_REPEAT):
    x2 = np.linalg.solve(ATA, ATB)
t3 = time.perf_counter()

# -------------------------------------------------------------------------
# 3) Least-squares solver (robust reference method)
# -------------------------------------------------------------------------
t4 = time.perf_counter()
for _ in range(N_REPEAT):
    x3, *_ = np.linalg.lstsq(A, B, rcond=None)
t5 = time.perf_counter()

# -------------------------------------------------------------------------
# Results
# -------------------------------------------------------------------------
print("===== Computational cost comparison =====")
print(f"Inverse method : {t1 - t0:.3f} s")
print(f"Solve method   : {t3 - t2:.3f} s")
print(f"lstsq method   : {t5 - t4:.3f} s")

print("\n===== Numerical differences =====")
print("max|x_inv - x_solve| =", np.max(np.abs(x1 - x2)))
print("max|x_solve - x_lstsq| =", np.max(np.abs(x2 - x3)))

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
    "E",
    "N",
    "U",
    "ATA",
    "ATB",
    "P_est",
    "dP_est",
    "P_app",
    "A",
    "B",
    "df_dX",
    "df_dY",
    "df_dZ",
    "distances",
    "dist_P_est_P_rnx_header",
    "fig",
    "i","x1","x2","x3","t0","t1","t2","t3","t4","t5","N_REPEAT"
]:
    if var in globals():
        del globals()[var]

del var
gc.collect()



# %%
# GNSS receiver positioning using pseudorange observations
#
# This block estimates the receiver position from GNSS code measurements (C1)
# using an iterative least-squares adjustment.
#
# Processing strategy
# -------------------
# Satellite clock errors and relativistic corrections are assumed known
# and are directly applied to the observation vector.
#
# At each iteration:
#
#   1. Compute approximate geometric distances between receiver and satellites
#   2. Form corrected observation equation:
#
#          ρ_obs = ρ_geom + c (dt_sat + dRel)
#
#   3. Linearize the observation model around an approximate receiver position
#   4. Build the design matrix A (line-of-sight unit vectors)
#   5. Estimate receiver position update using least squares:
#
#          dP = (AᵀA)⁻¹ Aᵀ B
#
#   6. Update receiver coordinates until convergence
#
# Convergence criterion
# ---------------------
# Iterations stop when the position update norm becomes smaller than 1 meter.
#
# Outputs
# -------
# P_est : estimated receiver position in ECEF coordinates (meters)
# ENU   : local East-North-Up offsets with respect to RINEX header position
#
# Educational note
# ----------------
# This implementation corresponds to the classical GNSS positioning model
# used in Single Point Positioning (SPP), simplified here by assuming that
# receiver clock bias has already been corrected or neglected.
###############################################################################


print('*****  Prise en compte des erreurs d''horloge satellites *****')

# Initialisation des coordonnées approximatives du récepteur
P_app = np.array([0, 0, 0])
dP_est=100
i=1
# Itération pour l'affinement de la position du récepteur
while np.linalg.norm(dP_est)>1:
    # Calcul des distances approximatives satellite-récepteur
    distances = np.sqrt((df_rnx['X_sat'].values - P_app[0])**2 +
                        (df_rnx['Y_sat'].values - P_app[1])**2 +
                        (df_rnx['Z_sat'].values - P_app[2])**2)

    # Vecteur des observations corrigées des différents modèles
    B = df_rnx['C1'].values - distances + conv.SPEED_OF_LIGHT * (df_rnx['dte_sat'].values + df_rnx['dRelat'].values )

    # Construction de la matrice des dérivées partielles
    df_dX = (P_app[0] - df_rnx['X_sat'].values) / distances
    df_dY = (P_app[1] - df_rnx['Y_sat'].values) / distances
    df_dZ = (P_app[2] - df_rnx['Z_sat'].values) / distances
    A = np.column_stack((df_dX, df_dY, df_dZ))

    # Résolution par moindres carrés pour estimer le déplacement
    dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
    
    # Mise à jour de la position estimée
    P_est = P_app + dP_est
    
    # Affichage de la position estimée à chaque itération
    print(f"Iteration {i}: Position estimée - X: {P_est[0]}, Y: {P_est[1]}, Z: {P_est[2]}")
    P_app = P_est  # Mise à jour de la position approximative pour la prochaine itération
    i+=1

# Calcul de la distance finale entre la position estimée et la position initiale du header RINEX
dist_P_est_P_rnx_header = np.sqrt(np.sum((P_est - P_rnx_header)**2))
print("\n")
print("Distance entre la position estimée et la position initiale du header RINEX:", dist_P_est_P_rnx_header)

# Calculer les coordonnées ENU locales
E, N, U = conv.xyz2enu(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
print("Est (E):", E)
print("Nord (N):", N)
print("Up (U):", U)
print("\n")


gnss_edu.plot_residual_analysis(A, B, dP_est, figure_title="Satellite clocks correction", 
                                save_path = folder / "02_sat_clocks_only.png",
                                P_est=P_est, P_rnx_header=P_rnx_header)

# -------------------------------------------------------------------------
# Store results for later reuse
# -------------------------------------------------------------------------
S02_Clk_Sat_solution = {
    "P_est": P_est.copy(),
    "ENU": np.array([E, N, U]),
    "distance_to_header": dist_P_est_P_rnx_header,
    "iterations": i-1
}


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
    "E",
    "N",
    "U",
    "P_est",
    "dP_est",
    "P_app",
    "A",
    "B",
    "df_dX",
    "df_dY",
    "df_dZ",
    "distances",
    "dist_P_est_P_rnx_header",
    "fig",
    "i",
]:
    if var in globals():
        del globals()[var]

del var
gc.collect()

# %%
#
# GNSS receiver positioning using code observations
# Satellite and receiver clock corrections
#
# Educational objective
# ---------------------
# Estimate the receiver position from GNSS pseudorange observations while:
#   - correcting satellite clock errors (known from SP3 products),
#   - estimating one receiver clock bias per epoch.
#
# Observation model (simplified SPP formulation)
# ----------------------------------------------
#   ρ = R + c(dt_r - dt_s) + c·dRel
#
# with:
#   R    : geometric range satellite–receiver
#   dt_s : satellite clock correction (known from SP3)
#   dt_r : receiver clock bias (estimated here per epoch)
#   dRel : relativistic correction
#
# Implementation note
# -------------------
# The receiver clock is represented by a block of columns appended to the
# design matrix A: one column per epoch (indicator matrix).
# A fast vectorized construction is preferred over a Python loop.
###############################################################################

print("*****  Prise en compte des erreurs d'horloge satellites et récepteur *****")

# -------------------------------------------------------------------------
# Receiver clock block (one parameter per epoch)
# -------------------------------------------------------------------------
# We build an indicator matrix block_dt_r of size (n_obs, n_epochs) such that:
#   block_dt_r[k, i] = 1 if observation k belongs to epoch i, else 0.

# --- (Reference / naive approach) ---
# This approach is explicit but slower and relies on an auxiliary 'ind_ligne'.
#
# epoch_uniques = df_rnx.index.get_level_values("epoch").unique()
# nb_epochs = len(epoch_uniques)
# block_dt_r = np.zeros((len(df_rnx), nb_epochs))
# for i, epoch in enumerate(epoch_uniques):
#     block_dt_r[df_rnx.loc[epoch, "ind_ligne"], i] = 1

# --- Optimized vectorized approach (recommended) ----------------------------
#
# Goal:
# Build the receiver clock block appended to the design matrix A.
#
# Each observation must be associated with ONE receiver clock parameter
# corresponding to its observation epoch.
#
# Principle:
# pd.factorize assigns an integer identifier to each unique epoch while
# preserving the order of observations in df_rnx.
#
# Example:
#   epochs → [t1, t1, t1, t2, t2, t3]
#   codes  → [0 , 0 , 0 , 1 , 1 , 2 ]
#
# This allows direct construction of an indicator matrix without looping
# over epochs (vectorized implementation).
#
# Result:
# block_dt_r[k, i] = 1
#     if observation k belongs to epoch i
#     else 0
# ---------------------------------------------------------------------------

# Assign integer code to each observation epoch
epoch_codes, epoch_uniques = pd.factorize(
    df_rnx.index.get_level_values("epoch")
)

# Number of receiver clock parameters (one per epoch)
nb_epochs = len(epoch_uniques)

# Initialize receiver clock design matrix block
block_dt_r = np.zeros((len(df_rnx), nb_epochs))

# Vectorized filling of the indicator matrix
block_dt_r[np.arange(len(df_rnx)), epoch_codes] = 1.0


# -------------------------------------------------------------------------
# Initial receiver position guess
# -------------------------------------------------------------------------
P_app = np.array([0.0, 0.0, 0.0])
dP_est = np.array([100.0, 100.0, 100.0])  # dummy init to enter loop
i = 1


# -------------------------------------------------------------------------
# Iterative least-squares adjustment
# -------------------------------------------------------------------------
while np.linalg.norm(dP_est[0:3]) > 1.0:

    # Approximate geometric distances
    distances = np.sqrt(
        (df_rnx["X_sat"].values - P_app[0]) ** 2 +
        (df_rnx["Y_sat"].values - P_app[1]) ** 2 +
        (df_rnx["Z_sat"].values - P_app[2]) ** 2
    )

    # Observation vector (meters)
    # C1 (m) - R (m) + c*(dt_s + dRel) (m)
    B = (
        df_rnx["C1"].values
        - distances
        + conv.SPEED_OF_LIGHT * (df_rnx["dte_sat"].values + df_rnx["dRelat"].values)
    )

    # Design matrix:
    # - first three columns: line-of-sight unit vectors (geometry)
    # - appended block: receiver clock parameters (one per epoch)
    df_dX = (P_app[0] - df_rnx["X_sat"].values) / distances
    df_dY = (P_app[1] - df_rnx["Y_sat"].values) / distances
    df_dZ = (P_app[2] - df_rnx["Z_sat"].values) / distances

    A = np.column_stack((df_dX, df_dY, df_dZ, block_dt_r))

    # Least-squares solution
    dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)

    # Update receiver coordinates
    P_est = P_app + dP_est[0:3]

    print(f"Iteration {i}: X={P_est[0]:.3f}, Y={P_est[1]:.3f}, Z={P_est[2]:.3f}")

    P_app = P_est
    i += 1


# -------------------------------------------------------------------------
# Final diagnostics
# -------------------------------------------------------------------------
dist_P_est_P_rnx_header = np.sqrt(np.sum((P_est - P_rnx_header) ** 2))
print("\nDistance entre la position estimée et le header RINEX:",
      dist_P_est_P_rnx_header)

E, N, U = conv.xyz2enu(
    P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],
    P_est[0], P_est[1], P_est[2]
)
print(f"Est (E): {E}")
print(f"Nord (N): {N}")
print(f"Haut (U): {U}\n")


# -------------------------------------------------------------------------
# Residual analysis
# -------------------------------------------------------------------------
gnss_edu.plot_residual_analysis(
    A, B, dP_est,
    figure_title="Satellite clocks and receiver correction",
    save_path=folder / "03_sat_rec_clocks.png",
    P_est=P_est,
    P_rnx_header=P_rnx_header
)


# -------------------------------------------------------------------------
# Store results for later reuse
# -------------------------------------------------------------------------
S03_Clk_Sat_Rec_solution = {
    "P_est": P_est.copy(),
    "ENU": np.array([E, N, U]),
    "distance_to_header": dist_P_est_P_rnx_header,
    "iterations": i - 1
}


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
    "E",
    "N",
    "U",
    "P_est",
    "dP_est",
    "P_app",
    "A",
    "B",
    "df_dX",
    "df_dY",
    "df_dZ",
    "distances",
    "dist_P_est_P_rnx_header",
    "fig",
    "i","nb_epochs","epoch","block_dt_r","epoch_uniques"
]:
    if var in globals():
        del globals()[var]

del var
gc.collect()


# %%
# GNSS receiver positioning: satellite + receiver clocks, plus Sagnac effect
#
# Educational objective
# ---------------------
# Introduce the Earth rotation (Sagnac) correction in the computation of
# satellite positions used in the pseudorange model.
#
# Physical idea
# -------------
# During the signal travel time τ ≈ C1 / c, the Earth rotates by an angle:
#
#       Δθ = ω_E * τ
#
# where ω_E is the Earth's mean angular rotation rate.
#
# Reference-frame consistency
# ---------------------------
# Satellite coordinates from SP3 are expressed in the Earth-fixed frame at
# signal emission time, whereas the receiver is naturally expressed in the
# Earth-fixed frame at reception time.
#
# To compare satellite and receiver in the SAME frame (reception frame),
# we compensate the Earth's rotation that occurred during propagation:
#
#       X_sat(reception) = Rz(-Δθ) · X_sat(emission)
#
# Hence the rotation applied to satellite coordinates uses the negative angle.
#
# Implementation note (important)
# -------------------------------
# Using df.apply(axis=1) is convenient but slow because it loops in Python.
# Here we implement a fully vectorized version to keep the TP fast and to
# illustrate good scientific Python practices.
###############################################################################

print("*****  Prise en compte des erreurs d'horloge satellites et recepteur *****")
print("*****  + Sagnac *****")

# -------------------------------------------------------------------------
# Vectorized Sagnac rotation of satellite coordinates
# -------------------------------------------------------------------------
# Inputs expected in df_rnx:
#   - C1     : pseudorange (m)
#   - X_sat, Y_sat, Z_sat : satellite coordinates in ECEF (m)
#
# Output:
#   df_Sagnac with rotated satellite coordinates (m), expressed in the
#   Earth-fixed frame at reception time.
# -------------------------------------------------------------------------

# Signal flight time (seconds): τ ≈ C1 / c
tau_s = df_rnx["C1"].to_numpy(dtype=float) / conv.SPEED_OF_LIGHT

# Physical Earth rotation angle during signal travel time (radians)
dtheta = tau_s * conv.EARTH_ROTATION_MEAN_ANGULAR_VELOCITY

# We rotate satellite coordinates backward to express them in the reception frame
alpha = -dtheta

cos_a = np.cos(alpha)
sin_a = np.sin(alpha)

x = df_rnx["X_sat"].to_numpy(dtype=float)
y = df_rnx["Y_sat"].to_numpy(dtype=float)
z = df_rnx["Z_sat"].to_numpy(dtype=float)

# Rotation around Z-axis (ECEF):
# [x'] = [ cos -sin  0 ] [x]
# [y']   [ sin  cos  0 ] [y]
# [z']   [  0    0   1 ] [z]
x_rot = cos_a * x - sin_a * y
y_rot = sin_a * x + cos_a * y
z_rot = z

# Store results in a DataFrame aligned with df_rnx index
df_Sagnac = pd.DataFrame(
    {"X_sat": x_rot, "Y_sat": y_rot, "Z_sat": z_rot},
    index=df_rnx.index
)

# Option (recommended for clarity in later steps):
# keep both versions in df_rnx to compare models
df_rnx["X_sat_sagnac"] = df_Sagnac["X_sat"]
df_rnx["Y_sat_sagnac"] = df_Sagnac["Y_sat"]
df_rnx["Z_sat_sagnac"] = df_Sagnac["Z_sat"]


# Assign integer code to each observation epoch
epoch_codes, epoch_uniques = pd.factorize(
    df_rnx.index.get_level_values("epoch")
)

# Number of receiver clock parameters (one per epoch)
nb_epochs = len(epoch_uniques)

# Initialize receiver clock design matrix block
block_dt_r = np.zeros((len(df_rnx), nb_epochs))

# Vectorized filling of the indicator matrix
block_dt_r[np.arange(len(df_rnx)), epoch_codes] = 1.0


# -------------------------------------------------------------------------
# Initial receiver position guess
# -------------------------------------------------------------------------
P_app = np.array([0.0, 0.0, 0.0])
dP_est = np.array([100.0, 100.0, 100.0])  # dummy init to enter loop
i = 1


# -------------------------------------------------------------------------
# Iterative least-squares adjustment
# -------------------------------------------------------------------------
while np.linalg.norm(dP_est[0:3]) > 1.0:

    # Approximate geometric distances
    distances = np.sqrt(
        (df_rnx["X_sat_sagnac"].values - P_app[0]) ** 2 +
        (df_rnx["Y_sat_sagnac"].values - P_app[1]) ** 2 +
        (df_rnx["Z_sat_sagnac"].values - P_app[2]) ** 2
    )

    # Observation vector (meters)
    # C1 (m) - R (m) + c*(dt_s + dRel) (m)
    B = (
        df_rnx["C1"].values
        - distances
        + conv.SPEED_OF_LIGHT * (df_rnx["dte_sat"].values + df_rnx["dRelat"].values)
    )

    # Design matrix:
    # - first three columns: line-of-sight unit vectors (geometry)
    # - appended block: receiver clock parameters (one per epoch)
    df_dX = (P_app[0] - df_rnx["X_sat_sagnac"].values) / distances
    df_dY = (P_app[1] - df_rnx["Y_sat_sagnac"].values) / distances
    df_dZ = (P_app[2] - df_rnx["Z_sat_sagnac"].values) / distances

    A = np.column_stack((df_dX, df_dY, df_dZ, block_dt_r))

    # Least-squares solution
    dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)

    # Update receiver coordinates
    P_est = P_app + dP_est[0:3]

    print(f"Iteration {i}: X={P_est[0]:.3f}, Y={P_est[1]:.3f}, Z={P_est[2]:.3f}")

    P_app = P_est
    i += 1


# -------------------------------------------------------------------------
# Final diagnostics
# -------------------------------------------------------------------------
dist_P_est_P_rnx_header = np.sqrt(np.sum((P_est - P_rnx_header) ** 2))
print("\nDistance entre la position estimée et le header RINEX:",
      dist_P_est_P_rnx_header)

E, N, U = conv.xyz2enu(
    P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],
    P_est[0], P_est[1], P_est[2]
)
print(f"Est (E): {E}")
print(f"Nord (N): {N}")
print(f"Haut (U): {U}\n")


# -------------------------------------------------------------------------
# Residual analysis
# -------------------------------------------------------------------------
gnss_edu.plot_residual_analysis(
    A, B, dP_est,
    figure_title="Satellite clocks and receiver correction + Sagnac",
    save_path=folder / "04_sat_rec_clocks_sagnac.png",
    P_est=P_est,
    P_rnx_header=P_rnx_header
)


# -------------------------------------------------------------------------
# Store results for later reuse
# -------------------------------------------------------------------------
S04_Clk_Sat_Rec_Sagnac_solution = {
    "P_est": P_est.copy(),
    "ENU": np.array([E, N, U]),
    "distance_to_header": dist_P_est_P_rnx_header,
    "iterations": i - 1
}



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
    "E", "N", "U",
    "P_est", "dP_est", "P_app",
    "A", "B", "df_dX", "df_dY", "df_dZ",
    "distances", "dist_P_est_P_rnx_header",
    "fig",
    "i","nb_epochs","epoch","block_dt_r","epoch_uniques","dtheta","epoch_codes",
    "alpha","cos_a","sin_a", "tau_s","x","y","z","x_rot","y_rot","z_rot","df_Sagnac",
]:
    if var in globals():
        del globals()[var]

del var
gc.collect()



# %%
# -------------------------------------------------------------------------
# Interpretation of remaining positioning error
# -------------------------------------------------------------------------
#
# After applying precise orbit information, satellite and receiver clock
# corrections, relativistic effects and Sagnac correction, the remaining
# discrepancy between:
#
#   - the receiver position estimated from GNSS observations, and
#   - the approximate receiver coordinates provided in the RINEX header
#
# is now largely dominated by signal propagation effects.
#
# Observation:
# The residual error appears mainly on the vertical component (Up).
#
# Interpretation:
# Atmospheric propagation delays (ionosphere and troposphere) act primarily
# along the satellite line-of-sight. Due to GNSS satellite geometry, these
# errors are weakly constrained horizontally and are therefore mostly
# absorbed by the vertical component of the estimated position.
#
# This explains why vertical positioning is generally less accurate than
# horizontal positioning and motivates the need for atmospheric modeling.
# -------------------------------------------------------------------------

print(
    "Note: The remaining discrepancy with the RINEX header position is "
    "mainly vertical, illustrating the impact of atmospheric propagation "
    "errors (ionosphere and troposphere)."
)

# %%
# -------------------------------------------------------------------------
# Validity of local satellite geometry computation
# -------------------------------------------------------------------------
#
# The estimated receiver position now differs from the approximate RINEX
# header coordinates by about a few tens of meters.
#
# Although this accuracy is insufficient for precise positioning, it is
# largely adequate to define the local reference frame centered on the
# receiver antenna.
#
# In particular, an approximate receiver position with an error of a few
# tens of meters introduces only negligible errors in:
#   - satellite azimuth,
#   - satellite elevation angles.
#
# Therefore, the current receiver position estimate can safely be used to
# compute satellite directions in the local ENU frame.
#
# This step is essential before introducing atmospheric models, since both
# ionospheric and tropospheric delays strongly depend on satellite elevation.
# -------------------------------------------------------------------------

print(
    "\nReceiver position accuracy (~30 m) is sufficient "
    "to compute satellite azimuth and elevation angles.\n"
    "We can now analyze signal propagation effects "
    "(ionosphere and troposphere)."
)

