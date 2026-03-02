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
df_rnx_orig, l_rnx_head = files_rw.read_rinex2_obs(fichier_rnx,
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
# Traitement classique sur le code pour la position d'un récepteur GNSS
# corrections des erreurs d'horloges satellites
# -> elles sont connues -> correction directe du vecteur B

print('*****  Prise en compte des erreurs d''horloge satellites et récepteur *****')

# Obtention des époques uniques
epoch_uniques = df_rnx.index.get_level_values('epoch').unique()
nb_epochs = len(epoch_uniques)
# Squelette du bloc à concaténer à la matrice modèle A
block_dt_r = np.zeros((len(df_rnx), nb_epochs))
# Remplissage du bloc correspondant à l'estimation des erreurs d'horloge recepteur

for i, epoch in enumerate(epoch_uniques):
    block_dt_r[df_rnx.loc[epoch, 'ind_ligne'],i]=1

# Initialisation des coordonnées approximatives du récepteur
P_app = np.array([0, 0, 0])
dP_est= np.array([100, 100, 100])
i=1
# Itération pour l'affinement de la position du récepteur
while np.linalg.norm(dP_est[0:3])>1:
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
    A = np.column_stack((df_dX, df_dY, df_dZ, block_dt_r))

    # Résolution par moindres carrés pour estimer le déplacement
    dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
    
    # Mise à jour de la position estimée
    P_est = P_app + dP_est[0:3]
    
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
print("Haut (U):", U)
print("\n")

gnss_edu.plot_residual_analysis(A, B, dP_est, figure_title="Satellite clocks and receiver correction",
                                save_path = folder / "03_sat_rec_clocks.png",
                                P_est=P_est, P_rnx_header=P_rnx_header)


# -------------------------------------------------------------------------
# Store results for later reuse
# -------------------------------------------------------------------------
S03_Clk_Sat_Rec_solution = {
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
    "i","nb_epochs","epoch","block_dt_r","epoch_uniques"
]:
    if var in globals():
        del globals()[var]

del var
gc.collect()


# %%
# Traitement classique sur le code pour la position d'un récepteur GNSS
# corrections des erreurs d'horloges satellites
# -> elles sont connues -> correction directe du vecteur B

print('*****  Prise en compte des erreurs d''horloge satellites *****')
print('*****  + Sagnac *****')

# introduire une visualisation du résultat avec folium
def rotate_around_z(row):
    # C1 / c -> temps de vol
    # temps de vol * vitesse angulaire rotation terrestre -> alpha_rad
    alpha_rad = row['C1'] / conv.SPEED_OF_LIGHT * conv.EARTH_ROTATION_MEAN_ANGULAR_VELOCITY
    # Matrice de rotation autour de l'axe Z
    rz = np.array([[np.cos(alpha_rad), -np.sin(alpha_rad), 0],
                   [np.sin(alpha_rad), np.cos(alpha_rad), 0],
                   [0, 0, 1]])
    # Vecteur de position original
    original_vector = np.array([row['X_sat'], row['Y_sat'], row['Z_sat']])
    # Calculer le vecteur de position rotatif
    rotated_vector = rz.dot(original_vector)
    return pd.Series(rotated_vector, index=['X_sat', 'Y_sat', 'Z_sat'])

# Appliquer la rotation à chaque ligne et créer un nouveau dataframe avec les résultats
df_Sagnac = df_rnx.apply(rotate_around_z, axis=1)


# Initialisation des coordonnées approximatives du récepteur
P_app = np.array([0, 0, 0])
dP_est=100
i=1
# Itération pour l'affinement de la position du récepteur
while np.linalg.norm(dP_est)>1:
    # Calcul des distances approximatives satellite-récepteur
    distances = np.sqrt((df_Sagnac['X_sat'].values - P_app[0])**2 +
                        (df_Sagnac['Y_sat'].values - P_app[1])**2 +
                        (df_Sagnac['Z_sat'].values - P_app[2])**2)

    # Vecteur des observations corrigées des différents modèles
    B = df_rnx['C1'].values - distances + conv.SPEED_OF_LIGHT * (df_rnx['dte_sat'].values + df_rnx['dRelat'].values )

    # Construction de la matrice des dérivées partielles
    df_dX = (P_app[0] - df_Sagnac['X_sat'].values) / distances
    df_dY = (P_app[1] - df_Sagnac['Y_sat'].values) / distances
    df_dZ = (P_app[2] - df_Sagnac['Z_sat'].values) / distances
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
print("Haut (U):", U)
print("\n")

del E, N, U, P_est, dP_est, P_app, A, B, df_dX, df_dY, df_dZ, distances, i


# %%
# Traitement classique sur le code pour la position d'un récepteur GNSS
# corrections des erreurs d'horloges satellites
# -> elles sont connues -> correction directe du vecteur B
# estimation des erreurs d'horloges récepteur
# -> elles sont inconnues et doivent être estimées à chaque époque


print('*****  Prise en compte des erreurs d''horloge satellites *****')
print('*****  + Sagnac *****')
print('*****  + horloges récepteur *****')

# Obtention des époques uniques
epoch_uniques = df_rnx.index.get_level_values('epoch').unique()
nb_epochs = len(epoch_uniques)
# Squelette du bloc à concaténer à la matrice modèle A
block_dt_r = np.zeros((len(df_rnx), nb_epochs))
# Remplissage du bloc correspondant à l'estimation des erreurs d'horloge recepteur

for i, epoch in enumerate(epoch_uniques):
    block_dt_r[df_rnx.loc[epoch, 'ind_ligne'],i]=1


# Initialisation des coordonnées approximatives du récepteur
P_app = np.array([0, 0, 0])

dP_est=np.array([100, 100, 100])
i=1
# Itération pour l'affinement de la position du récepteur
while np.linalg.norm(dP_est[0:3])>1:
    
    # Calcul des distances approximatives satellite-récepteur
    distances = np.sqrt((df_Sagnac['X_sat'].values - P_app[0])**2 +
                        (df_Sagnac['Y_sat'].values - P_app[1])**2 +
                        (df_Sagnac['Z_sat'].values - P_app[2])**2)

    # Vecteur des observations corrigées des différents modèles
    B = df_rnx['C1'].values - distances + conv.SPEED_OF_LIGHT * (df_rnx['dte_sat'].values + df_rnx['dRelat'].values )

    # Construction de la matrice des dérivées partielles
    df_dX = (P_app[0] - df_Sagnac['X_sat'].values) / distances
    df_dY = (P_app[1] - df_Sagnac['Y_sat'].values) / distances
    df_dZ = (P_app[2] - df_Sagnac['Z_sat'].values) / distances
    A = np.column_stack((df_dX, df_dY, df_dZ, block_dt_r))

    # Résolution par moindres carrés pour estimer le déplacement
    dP_est = np.linalg.inv(A.T@A)@A.T@B     
    
    # Mise à jour de la position estimée
    # Attention : on a estimé nb_epochs paramètres d'horloge récepteur
    P_est = np.zeros(len(dP_est))
    
    P_est[0:3] = P_app[0:3]+dP_est[0:3]
    P_est[3:]  = dP_est[3:]
    
    
    # Affichage de la position estimée à chaque itération
    print(f"Iteration {i}: Position estimée - X: {P_est[0]}, Y: {P_est[1]}, Z: {P_est[2]}")
    P_app = P_est  # Mise à jour de la position approximative pour la prochaine itération
    i+=1


# Calcul de la distance finale entre la position estimée et la position initiale du header RINEX
dist_P_est_P_rnx_header = np.sqrt(np.sum((P_est[:3] - P_rnx_header)**2))
print("\n")
print("Distance entre la position estimée et la position initiale du header RINEX:", dist_P_est_P_rnx_header)

# Calculer les coordonnées ENU locales
E, N, U = conv.xyz2enu(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
print("Est (E):", E)
print("Nord (N):", N)
print("Haut (U):", U)
print("\n")

gnss_edu.plot_residual_analysis(A, B, dP_est, figure_title="Calcul corr Sat Rec et Sagnac", save_path="./corr_sat_clock_sagnac.png",
                           P_est=P_est[:3], P_rnx_header=P_rnx_header)


del E, N, U, P_est, dP_est, P_app, A, block_dt_r, B, df_dX, df_dY, df_dZ, distances, i, epoch, epoch_uniques, nb_epochs

# introduire une visualisation du résultat avec folium


# %%
# Traitement classique sur le code pour la position d'un récepteur GNSS
# corrections des erreurs d'horloges satellites
# -> elles sont connues -> correction directe du vecteur B
# estimation des erreurs d'horloges récepteur
# -> elles sont inconnues et doivent être estimées à chaque époque
# correction des effets ionosphériques via le modèle de Klobuchar
# Questions : quel signe à apporter à la correction +/- ?
#             correction sur le code ou sur la phase    ?
# Réponses :
#



# Pour Klobchar -> besoin du fichier de navigation
# Téléchargement automatique des données NAVIGATION pour le jour du traitement
# ici  (2019-176) qui correspond à la semaine GPS 2059
dwl_output_navigation = operational.download_gnss_rinex(statdico={"nav" : ["brdc"]},
                                output_dir=my_directory,
                                startdate= my_date_to_process ,
                                enddate= my_date_to_process ,
                                parallel_download = 1)

fichier_nav = dwl_output_navigation[0][0]
# Chargement du fichier de navigation RINEX via GeodeZYX
# comme Dusa et al. on retraité tous les fichiers de navigation RINEX2 -> RINEX3
# on en profite :-) !
iono_cor_dic, time_sys_corr_dic = files_rw.read_rinex_nav_v3_header(fichier_nav)


print('*****  Prise en compte des erreurs d''horloge satellites *****')
print('*****  + Sagnac *****')
print('*****  + horloges récepteur *****')
print('*****  + iono Klobuchar *****')

# Calcul des Az et des Ele de chaque mesure
# toolAzEle : azimut and elevation (radians) for one or several satellites Xs,Ys,Zs (scalar or vector) seen from a point with X,Y,Z coordinates.

Az_rad, Ele_rad = conv.xyz2azi_ele(P_rnx_header[0],P_rnx_header[1],P_rnx_header[2],df_rnx.X_sat,df_rnx.Y_sat,df_rnx.Z_sat)

# xyz2geo : cartesian to geographic coordinates conversion. All angles are given in radians.
lon,lat,h = conv.xyz2geo(P_rnx_header[0],P_rnx_header[1],P_rnx_header[2])

# radian to degree
rad2deg = 180 / np.pi
lon_d = lon * rad2deg
lat_d = lat * rad2deg

Az_deg = Az_rad * rad2deg
Ele_deg =  Ele_rad * rad2deg

# coefficients alpha et beta extrait du fichier de navigation
alpha = iono_cor_dic["GPSA"]
beta = iono_cor_dic["GPSB"]

# t = gpst.gpsdatetime() # création d'un objet temps de la gnsstime
dIon1 = []
for (time_i,prn_i) in df_rnx.index:
    #t.rinex_t(time_i.to_pydatetime().strftime('%y %m %d %H %M %S.%f'))
    #wsec_v = t.wsec       # seconds in GPS week
    t_dt = time_i.to_pydatetime()
    _ , wsec_v = conv.dt2gpstime(t_dt, secinweek=True)
    i = df_rnx.loc[(time_i,prn_i), 'ind_ligne'] 
    dIon1_v = gnss_edu.klobuchar(lat_d, lon_d, Ele_deg[i], Az_deg[i] , wsec_v, alpha, beta )
    dIon1.append(dIon1_v)

df_rnx['Az']      = Az_deg
df_rnx['Ele']     = Ele_deg
df_rnx['dIon1']   = dIon1 

del Az_rad, Ele_rad, lon, lat, lat_d, lon_d, h, rad2deg, Az_deg, Ele_deg, alpha, beta, dIon1, wsec_v, i, dIon1_v, time_i, prn_i


# Obtention des époques uniques
epoch_uniques = df_rnx.index.get_level_values('epoch').unique()
nb_epochs = len(epoch_uniques)
# Squelette du bloc à concaténer à la matrice modèle A
block_dt_r = np.zeros((len(df_rnx), nb_epochs))
# Remplissage du bloc correspondant à l'estimation des erreurs d'horloge recepteur

for i, epoch in enumerate(epoch_uniques):
    block_dt_r[df_rnx.loc[epoch, 'ind_ligne'],i]=1


# Initialisation des coordonnées approximatives du récepteur
P_app = np.array([0, 0, 0])

dP_est=np.array([100, 100, 100])
i=1
# Itération pour l'affinement de la position du récepteur
while np.linalg.norm(dP_est[0:3])>1:

    # Calcul des distances approximatives satellite-récepteur
    distances = np.sqrt((df_Sagnac['X_sat'].values - P_app[0])**2 +
                        (df_Sagnac['Y_sat'].values - P_app[1])**2 +
                        (df_Sagnac['Z_sat'].values - P_app[2])**2)

    # Vecteur des observations corrigées des différents modèles
    B = df_rnx['C1'].values - distances + conv.SPEED_OF_LIGHT * (df_rnx['dte_sat'].values + df_rnx['dRelat'].values ) - df_rnx['dIon1'].values

    # Construction de la matrice des dérivées partielles
    df_dX = (P_app[0] - df_Sagnac['X_sat'].values) / distances
    df_dY = (P_app[1] - df_Sagnac['Y_sat'].values) / distances
    df_dZ = (P_app[2] - df_Sagnac['Z_sat'].values) / distances
    A = np.column_stack((df_dX, df_dY, df_dZ, block_dt_r))

    # Résolution par moindres carrés pour estimer le déplacement
    dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)

    
    # Mise à jour de la position estimée
    # Attention : on a estimé nb_epochs paramètres d'horloge récepteur
    P_est = np.zeros(len(dP_est))
    
    P_est[0:3] = P_app[0:3]+dP_est[0:3]
    P_est[3:]  = dP_est[3:]
    

    # Affichage de la position estimée à chaque itération
    print(f"Iteration {i}: Position estimée - X: {P_est[0]}, Y: {P_est[1]}, Z: {P_est[2]}")
    P_app = P_est  # Mise à jour de la position approximative pour la prochaine itération
    i+=1


# Calcul de la distance finale entre la position estimée et la position initiale du header RINEX
dist_P_est_P_rnx_header = np.sqrt(np.sum((P_est[:3] - P_rnx_header)**2))
print("\n")
print("Distance entre la position estimée et la position initiale du header RINEX:", dist_P_est_P_rnx_header)

# Calculer les coordonnées ENU locales
E, N, U = conv.xyz2enu(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
print("Est (E):", E)
print("Nord (N):", N)
print("Haut (U):", U)
print("\n")

del E, N, U, P_est, dP_est, P_app, A, block_dt_r, B, df_dX, df_dY, df_dZ, distances, i, epoch, epoch_uniques, nb_epochs

# introduire une visualisation du résultat avec folium


# %%
# Traitement classique sur le code pour la position d'un récepteur GNSS
# corrections des erreurs d'horloges satellites
# -> elles sont connues -> correction directe du vecteur B
# estimation des erreurs d'horloges récepteur
# -> elles sont inconnues et doivent être estimées à chaque époque
# correction des effets ionosphériques en utilisant une CL des observations sur
# le code
# -> correction directement du vecteur B

print('*****  Prise en compte des erreurs d''horloge satellites *****')
print('*****  + Sagnac *****')
print('*****  + horloges récepteur *****')
print('*****  + iono CL(C1,P2) *****')

f1 = conv.L1_CARRIER_FREQUENCY
f2 = conv.L2_CARRIER_FREQUENCY

l1 = conv.L1_WAVELENGTH
l2 = conv.L2_WAVELENGTH

# On applique la formule de la combinaison linéaire
# l3 = (f1**2*l1- f2**2*l2)/(f1**2-f2**2);
# qui se simplifie par:
l3=conv.SPEED_OF_LIGHT/(f1+f2)


# attention : on se rend compte que les mesures L1 et L2 sont en cycle alors que l'on
# a besoin d'une mesure de distance (ambigue mais distance quand même)

df_rnx['L3']  = (f1**2*df_rnx['L1'] - f2**2*df_rnx['L2'])/(f1**2-f2**2)
df_rnx['P3']  = (f1**2*df_rnx['C1'] - f2**2*df_rnx['P2'])/(f1**2-f2**2)


# %%

del f1, f2, l1, l2

# On peut s'amuser à recalculer les positions des satellites avec les nouveaux
# temps de vol obtenus via P3
# impact inférieur à 2 mm sur chacune des composantes X_sat, Y_sat et Z_sat

# %%
# # Il faut calculer la position de chaque satellite GNSS à chaque temps d'émission
# t = gpst.gpsdatetime()

# X_sat = []
# Y_sat = []
# Z_sat = []
# dte_sat = []
# for (time_i,prn_i) in df_rnx.index:

#     t.rinex_t(time_i.to_pydatetime().strftime('%y %m %d %H %M %S.%f'))

#     t_emission_mjd  = t.mjd - df_rnx.loc[(time_i,prn_i), 'P3'] / conv.SPEED_OF_LIGHT / 86400.0

#     (X_sat_v,Y_sat_v,Z_sat_v,dte_sat_v)	 = mysp3.calcSatCoord(prn_i[0], int(prn_i[1:]),t_emission_mjd)


#     X_sat.append(X_sat_v)
#     Y_sat.append(Y_sat_v)
#     Z_sat.append(Z_sat_v)
#     dte_sat.append(dte_sat_v)


# df_rnx['X_sat_P3']   = X_sat
# df_rnx['Y_sat_P3']   = Y_sat
# df_rnx['Z_sat_P3']   = Z_sat
# df_rnx['dte_sat_P3'] = dte_sat


# del X_sat, Y_sat, Z_sat, dte_sat, time_i, prn_i, t_emission_mjd, t, X_sat_v, Y_sat_v, Z_sat_v, dte_sat_v

# #%%
# prns = df_rnx.index.get_level_values('prn').unique()

# # Créer une figure et un axe pour le plot
# fig, ax = plt.subplots(figsize=(10, 6))

# # Boucler sur chaque PRN et tracer sa série temporelle
# for prn in prns:
#     # Sélectionner les données pour le PRN actuel
#     data = df_rnx.xs(prn, level='prn')
#     # Tracer les données
#     ax.plot(data.index.get_level_values('epoch'), data['Z_sat_P3']-data['Z_sat'], label=prn)

# # Configurer le graphique
# ax.set_title('Séries temporelles Z_sat_P3 - Z_sat par PRN de satellite')
# ax.set_xlabel('Temps')
# ax.set_ylabel('Valeur')
# ax.legend(title='PRN')

# # Afficher le graphique
# plt.show()


# %%

# Obtention des époques uniques
epoch_uniques = df_rnx.index.get_level_values('epoch').unique()
nb_epochs = len(epoch_uniques)
# Squelette du bloc à concaténer à la matrice modèle A
block_dt_r = np.zeros((len(df_rnx), nb_epochs))
# Remplissage du bloc correspondant à l'estimation des erreurs d'horloge recepteur

for i, epoch in enumerate(epoch_uniques):
    block_dt_r[df_rnx.loc[epoch, 'ind_ligne'],i]=1


# Initialisation des coordonnées approximatives du récepteur
P_app = np.array([0, 0, 0])

dP_est=np.array([100, 100, 100])
i=1
# Itération pour l'affinement de la position du récepteur
while np.linalg.norm(dP_est[0:3])>1:

    # Calcul des distances approximatives satellite-récepteur
    distances = np.sqrt((df_Sagnac['X_sat'].values - P_app[0])**2 +
                        (df_Sagnac['Y_sat'].values - P_app[1])**2 +
                        (df_Sagnac['Z_sat'].values - P_app[2])**2)

    # Vecteur des observations corrigées des différents modèles
    B = df_rnx['P3'].values - distances + conv.SPEED_OF_LIGHT * (df_rnx['dte_sat'].values + df_rnx['dRelat'].values )

    # Construction de la matrice des dérivées partielles
    df_dX = (P_app[0] - df_Sagnac['X_sat'].values) / distances
    df_dY = (P_app[1] - df_Sagnac['Y_sat'].values) / distances
    df_dZ = (P_app[2] - df_Sagnac['Z_sat'].values) / distances
    A = np.column_stack((df_dX, df_dY, df_dZ, block_dt_r))

    # Résolution par moindres carrés pour estimer le déplacement
    dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
    dP_est = np.linalg.inv(A.T@A)@A.T@B
    # Mise à jour de la position estimée
    # Attention : on a estimé nb_epochs paramètres d'horloge récepteur
    P_est = np.zeros(len(dP_est))
    
    
    P_est[0:3] = P_app[0:3]+dP_est[0:3]
    P_est[3:]  = dP_est[3:]
    

    
    # Affichage de la position estimée à chaque itération
    print(f"Iteration {i}: Position estimée - X: {P_est[0]}, Y: {P_est[1]}, Z: {P_est[2]}")
    P_app = P_est  # Mise à jour de la position approximative pour la prochaine itération
    i+=1


# Calcul de la distance finale entre la position estimée et la position initiale du header RINEX
dist_P_est_P_rnx_header = np.sqrt(np.sum((P_est[:3] - P_rnx_header)**2))
print("\n")
print("Distance entre la position estimée et la position initiale du header RINEX:", dist_P_est_P_rnx_header)

# Calculer les coordonnées ENU locales
E, N, U = conv.xyz2enu(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
print("Est (E):", E)
print("Nord (N):", N)
print("Haut (U):", U)
print("\n")


gnss_edu.plot_residual_analysis(A, B, dP_est, figure_title="Calcul corr Sat Rec et Sagnac et iono", save_path="./corr_sat_clock_sagnac_iono.png",
                           P_est=P_est[:3], P_rnx_header=P_rnx_header)

del E, N, U, P_est, dP_est, P_app, A, block_dt_r, B, df_dX, df_dY, df_dZ, distances, i, epoch, epoch_uniques, nb_epochs

# introduire une visualisation du résultat avec folium


# %%
# Traitement classique sur le code pour la position d'un récepteur GNSS
# corrections des erreurs d'horloges satellites
# -> elles sont connues -> correction directe du vecteur B
# estimation des erreurs d'horloges récepteur
# -> elles sont inconnues et doivent être estimées à chaque époque
# correction des effets ionosphériques en utilisant une CL des observations sur
# le code
# -> correction directement du vecteur B

# correction de l'effet Sagnac
# c'est lui qui se voit comme le nez au milieu de la figure alors autant ne plus
# tourner autour du pot

print('*****  On doit obtenir la même chose que précédemment *****')

def rotate_around_z(row):
    # C1 / c -> temps de vol
    # temps de vol * vitesse angulaire rotation terrestre -> alpha_rad
    alpha_rad = row['C1'] / conv.SPEED_OF_LIGHT * conv.EARTH_ROTATION_MEAN_ANGULAR_VELOCITY
    # Matrice de rotation autour de l'axe Z
    Rz = np.array([[np.cos(alpha_rad), -np.sin(alpha_rad), 0],
                   [np.sin(alpha_rad), np.cos(alpha_rad), 0],
                   [0, 0, 1]])
    # Vecteur de position original
    original_vector = np.array([row['X_sat'], row['Y_sat'], row['Z_sat']])
    # Calculer le vecteur de position rotatif
    rotated_vector = Rz.dot(original_vector)
    return pd.Series(rotated_vector, index=['X_sat', 'Y_sat', 'Z_sat'])

# Appliquer la rotation à chaque ligne et créer un nouveau dataframe avec les résultats
df_Sagnac = df_rnx.apply(rotate_around_z, axis=1)


# Obtention des époques uniques
epoch_uniques = df_rnx.index.get_level_values('epoch').unique()
nb_epochs = len(epoch_uniques)
# Squelette du bloc à concaténer à la matrice modèle A
block_dt_r = np.zeros((len(df_rnx), nb_epochs))
# Remplissage du bloc correspondant à l'estimation des erreurs d'horloge recepteur

for i, epoch in enumerate(epoch_uniques):
    block_dt_r[df_rnx.loc[epoch, 'ind_ligne'],i]=1


# Initialisation des coordonnées approximatives du récepteur
P_app = np.array([0, 0, 0])

dP_est=np.array([100, 100, 100])
i=1
# Itération pour l'affinement de la position du récepteur
while np.linalg.norm(dP_est[0:3])>1:

    # Calcul des distances approximatives satellite-récepteur
    distances = np.sqrt((df_Sagnac['X_sat'].values - P_app[0])**2 +
                        (df_Sagnac['Y_sat'].values - P_app[1])**2 +
                        (df_Sagnac['Z_sat'].values - P_app[2])**2)

    # Vecteur des observations corrigées des différents modèles
    B = df_rnx['P3'].values - distances + conv.SPEED_OF_LIGHT * (df_rnx['dte_sat'].values + df_rnx['dRelat'].values )

    # Construction de la matrice des dérivées partielles
    df_dX = (P_app[0] - df_Sagnac['X_sat'].values) / distances
    df_dY = (P_app[1] - df_Sagnac['Y_sat'].values) / distances
    df_dZ = (P_app[2] - df_Sagnac['Z_sat'].values) / distances
    A = np.column_stack((df_dX, df_dY, df_dZ, block_dt_r))

    # Résolution par moindres carrés pour estimer le déplacement
    dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
    
    # Mise à jour de la position estimée
    # Attention : on a estimé nb_epochs paramètres d'horloge récepteur
    P_est = np.zeros(len(dP_est))
    
    
    P_est[0:3] = P_app[0:3]+dP_est[0:3]
    P_est[3:]  = dP_est[3:]
    

    
    # Affichage de la position estimée à chaque itération
    print(f"Iteration {i}: Position estimée - X: {P_est[0]}, Y: {P_est[1]}, Z: {P_est[2]}")
    P_app = P_est  # Mise à jour de la position approximative pour la prochaine itération
    i+=1


# Calcul de la distance finale entre la position estimée et la position initiale du header RINEX
dist_P_est_P_rnx_header = np.sqrt(np.sum((P_est[:3] - P_rnx_header)**2))
print("\n")
print("Distance entre la position estimée et la position initiale du header RINEX:", dist_P_est_P_rnx_header)

# Calculer les coordonnées ENU locales
E, N, U = conv.xyz2enu(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
print("Est (E):", E)
print("Nord (N):", N)
print("Haut (U):", U)
print("\n")

del E, N, U, P_est, dP_est, P_app, A, block_dt_r, B, df_dX, df_dY, df_dZ, distances, i, epoch, epoch_uniques, nb_epochs

# introduire une visualisation du résultat avec folium


# %%
# Traitement classique sur le code pour la position d'un récepteur GNSS
# corrections des erreurs d'horloges satellites
# -> elles sont connues -> correction directe du vecteur B
# estimation des erreurs d'horloges récepteur
# -> elles sont inconnues et doivent être estimées à chaque époque
# correction des effets ionosphériques en utilisant une CL des observations sur
# le code
# -> correction directement du vecteur B

# correction de l'effet Sagnac
# c'est lui qui se voit comme le nez au milieu de la figure alors autant ne plus
# tourner autour du pot

# correction des effets troposphériques

print('*****  Prise en compte des erreurs d''horloge satellites *****')
print('*****  + Sagnac *****')
print('*****  + horloges récepteur *****')
print('*****  + iono CL(C1,P2) *****')
print('*****  + TROPO     *****')
print('*****  + angle de coupure     *****')


# xyz2geo : cartesian to geographic coordinates conversion. All angles are given in radians.
lon, lat, h = conv.xyz2geo(P_rnx_header[0],P_rnx_header[1],P_rnx_header[2])

# grille gpt3
gpt3_5 = gnss_edu.gpt3_5_fast_readGrid(filename ='gpt3_5.grd')


# Il faut calculer les corrections et les estimations pour chaque satellite et
# pour chaque époque
# t = gpst.gpsdatetime()

ZHD = []
mfh = []
mfw = []

for (time_i,prn_i) in df_rnx.index:

    #t.rinex_t(time_i.to_pydatetime().strftime('%y %m %d %H %M %S.%f'))
    t_mjd = conv.dt2mjd(time_i.to_pydatetime())
    #gmfh, gmfw = gpt3.gmf(dmjd=t.mjd, dlat=lat, dlon=lon, dhgt=h, zd =np.pi/2-df_rnx.loc[(time_i,prn_i), 'Ele']*np.pi/180)
    gmfh = 1/np.sin(df_rnx.loc[(time_i,prn_i), 'Ele']*np.pi/180)
    gmfw = gmfh
    T = gnss_edu.gpt3_5_fast(mjd=t_mjd, lat=np.array([lat]), lon=np.array([lon]), h_ell=np.array([h]), it=0, grid=gpt3_5)
    gm = 1 - 0.00265 * np.cos(2*lat) - 0.000285 * h*1e-3
    pression = T[0][0][0] #hPa
    ZHD_v = 0.0022768 * pression / gm
    
    
    mfh.append(gmfh)
    mfw.append(gmfw)
    ZHD.append(ZHD_v)


df_rnx['ZHD']   = ZHD 
df_rnx['mfh']   = mfh 
df_rnx['mfw']   = mfw


del gm, pression, ZHD, mfh, mfw, ZHD_v, gmfh, gmfw, gpt3_5, lat, lon, h


# %% Angle de coupure

Az_rad, Ele_rad = conv.xyz2azi_ele(P_rnx_header[0],P_rnx_header[1],P_rnx_header[2],df_rnx.X_sat,df_rnx.Y_sat,df_rnx.Z_sat)

# radian to degree
rad2deg = 180 / np.pi

Ele_deg =  Ele_rad * rad2deg


df_rnx_new = df_rnx[df_rnx['Ele']>7].copy()
df_Sagnac_new = df_Sagnac[df_rnx['Ele']>7].copy()

# ajout de l'indice de ligne
df_rnx_new['ind_ligne'] = range(len(df_rnx_new)) 

# %%
# Obtention des époques uniques
epoch_uniques = df_rnx_new.index.get_level_values('epoch').unique()
nb_epochs = len(epoch_uniques)
# Squelette du bloc à concaténer à la matrice modèle A
block_dt_r = np.zeros((len(df_rnx_new), nb_epochs))

# Remplissage du bloc correspondant à l'estimation des erreurs d'horloge recepteur
for i, epoch in enumerate(epoch_uniques):
    block_dt_r[df_rnx_new.loc[epoch, 'ind_ligne'],i]=1

# Remplissage du bloc correspondant à l'estimation des retards troposphériques humides (estimation horaire)
start = df_rnx_new.index[0][0]
delta_T = pd.Timedelta(days=0, hours=2, minutes=0)
delta_sec = pd.Timedelta(seconds=1)
end   = start + delta_T - delta_sec


import math

delta = df_rnx_new.index[-1][0] - df_rnx_new.index[0][0]
nb_par_zwd = math.ceil(delta / pd.Timedelta(hours=2))


# Squelette du bloc à concaténer à la matrice modèle A
block_zwd = np.zeros((len(df_rnx_new), nb_par_zwd))
ind_c=0

while start <= df_rnx_new.index[-1][0]:
    
    end   = start + delta_T - delta_sec
    extract2c = df_rnx_new.loc[(slice(start, end)), ('ind_ligne','mfw')]   
    block_zwd[extract2c['ind_ligne'],ind_c]=extract2c['mfw']
    
    start = start + delta_T
    ind_c = ind_c + 1

# Initialisation des coordonnées approximatives du récepteur
P_app = np.array([0, 0, 0])

dP_est=np.array([100, 100, 100])
i=1
# Itération pour l'affinement de la position du récepteur
while np.linalg.norm(dP_est[0:3])>1:

    # Calcul des distances approximatives satellite-récepteur
    distances = np.sqrt((df_Sagnac_new['X_sat'].values - P_app[0])**2 +
                        (df_Sagnac_new['Y_sat'].values - P_app[1])**2 +
                        (df_Sagnac_new['Z_sat'].values - P_app[2])**2)

    # Vecteur des observations corrigées des différents modèles
    B = df_rnx_new['P3'].values - distances + conv.SPEED_OF_LIGHT * (df_rnx_new['dte_sat'].values
    + df_rnx_new['dRelat'].values )
    - df_rnx_new['ZHD'].values* df_rnx_new['mfh'].values

    # Construction de la matrice des dérivées partielles
    df_dX = (P_app[0] - df_Sagnac_new['X_sat'].values) / distances
    df_dY = (P_app[1] - df_Sagnac_new['Y_sat'].values) / distances
    df_dZ = (P_app[2] - df_Sagnac_new['Z_sat'].values) / distances
    A = np.column_stack((df_dX, df_dY, df_dZ, block_dt_r, block_zwd))

    # Résolution par moindres carrés pour estimer le déplacement
    dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
    
    # Mise à jour de la position estimée
    # Attention : on a estimé nb_epochs paramètres d'horloge récepteur
    P_est = np.zeros(len(dP_est))
    
    
    P_est[0:3] = P_app[0:3]+dP_est[0:3]
    P_est[3:]  = dP_est[3:]
    

    
    # Affichage de la position estimée à chaque itération
    print(f"Iteration {i}: Position estimée - X: {P_est[0]}, Y: {P_est[1]}, Z: {P_est[2]}")
    P_app = P_est  # Mise à jour de la position approximative pour la prochaine itération
    i+=1


# Calcul de la distance finale entre la position estimée et la position initiale du header RINEX
dist_P_est_P_rnx_header = np.sqrt(np.sum((P_est[:3] - P_rnx_header)**2))
print("\n")
print("Distance entre la position estimée et la position initiale du header RINEX:", dist_P_est_P_rnx_header)

# Calculer les coordonnées ENU locales
E, N, U = conv.xyz2enu(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
print("Est (E):", E)
print("Nord (N):", N)
print("Haut (U):", U)
print("\n")

# del E, N, U, P_est, dP_est, P_app, A, block_dt_r, B, df_dX, df_dY, df_dZ, distances, i, epoch, epoch_uniques, nb_epochs

# introduire une visualisation du résultat avec folium

gnss_edu.plot_residual_analysis(A, B, dP_est, figure_title="Calcul corr Sat Rec et Sagnac et iono et tropo", save_path="./corr_sat_clock_sagnac_iono_tropo.png",
                           P_est=P_est[:3], P_rnx_header=P_rnx_header)

# %%


# %%


# %%
# Recherche des observations continues par satellite
df_reset = df_rnx_new.reset_index()
df_reset['Time_Diff'] = df_reset.groupby('prn')['epoch'].diff().dt.total_seconds()

# génération d'un dataframe listant les ruptures d'observation pour les satellites concernés
# attention, certains satellites n'ont qu'une période d'observation et il n'y a qu'une seule ambiguïté
# à introduire

filtered_df = df_reset[df_reset['Time_Diff'] > 30]


prns_unique_all = sorted(df_reset['prn'].unique())
prns_unique_multi_amb = sorted(filtered_df['prn'].unique())  
prns_unique_one_amb = [x for x in prns_unique_all if x not in prns_unique_multi_amb]

# Ajout à la matrice A d'un bloc correspond aux ambiguïtés prns_unique_one_amb
block_one_amb = np.zeros((len(df_rnx_new), len(prns_unique_one_amb)))

ind_c = 0
for prn_i in prns_unique_one_amb:
    ind_ligne = np.where(df_rnx_new.index.get_level_values('prn') == prn_i )
    block_one_amb[ ind_ligne, ind_c] = 1 
    ind_c = ind_c + 1

# %%
# Création de blocs spécifiques pour les satellites qui se lèvent/couchent
for prn_i in prns_unique_multi_amb:
    print(prn_i)
    nb_amb_prn = sorted(filtered_df['prn']).count(prn_i)+1
    
    print(nb_amb_prn) # je dois rajouter nb_amb_prn colonnes dans la matrice block_multi_amb
    
    df_filtered_var = df_reset[np.isnan(df_reset['Time_Diff']) & (df_reset['prn'] == prn_i)]
    first_epoch_value = df_filtered_var['epoch'].reset_index(drop=True)
    extract_df = filtered_df[filtered_df['prn']==prn_i]['epoch'].reset_index(drop=True)
    
    extract_df = pd.concat([first_epoch_value, extract_df], ignore_index=True)
    
    # Extraction de la dernière ligne qui satisfait la condition
    derniere_ligne = df_reset[df_reset['prn'] == prn_i]['epoch'].iloc[-1]
    
    # Assurez-vous que 'extract_df' est un DataFrame avec une colonne nommée 'epoch'
    extract_df = pd.DataFrame(extract_df, columns=['epoch'])

    # Assurez-vous que 'derniere_ligne_df' est également un DataFrame avec une colonne 'epoch'
    derniere_ligne_df = pd.DataFrame([derniere_ligne], columns=['epoch'])
    
    # Convertir derniere_ligne_df à datetime64[us] pour correspondre à extract_df
    derniere_ligne_df['epoch'] = derniere_ligne_df['epoch'].astype('datetime64[us]')

    # Concaténation le long de l'axe 0 (ajout de lignes)
    extract_df = pd.concat([extract_df, derniere_ligne_df], ignore_index=True, axis=0)

    print(extract_df)
    
    add_block_amb=np.zeros((len(df_rnx_new), nb_amb_prn))
    
    for  amb in range(nb_amb_prn):
        print(amb)  
        print(extract_df.iloc[amb]['epoch'])


# %%
# Traitement classique sur la phase pour la position d'un récepteur GNSS
# corrections des erreurs d'horloges satellites
# -> elles sont connues -> correction directe du vecteur B
# estimation des erreurs d'horloges récepteur
# -> elles sont inconnues et doivent être estimées à chaque époque
# correction des effets ionosphériques en utilisant une CL des observations sur
# le code
# -> correction directement du vecteur B

# correction de l'effet Sagnac
# c'est lui qui se voit comme le nez au milieu de la figure alors autant ne plus
# tourner autour du pot

# correction des effets troposphériques

print('*****  Prise en compte des erreurs d''horloge satellites *****')
print('*****  + Sagnac *****')
print('*****  + horloges récepteur *****')
print('*****  + iono CL(L1,L2) *****')
print('*****  + TROPO     *****')
print('*****  + angle de coupure     *****')
print('*****  + estimation amb.     *****')

print('*****  ET NORMALEMENT ... ARG :) on est dans PLOUF  *****')


# %%
# Obtention des époques uniques
epoch_uniques = df_rnx_new.index.get_level_values('epoch').unique()
nb_epochs = len(epoch_uniques)
# Squelette du bloc à concaténer à la matrice modèle A
block_dt_r = np.zeros((len(df_rnx_new), nb_epochs))

# Remplissage du bloc correspondant à l'estimation des erreurs d'horloge recepteur
for i, epoch in enumerate(epoch_uniques):
    block_dt_r[df_rnx_new.loc[epoch, 'ind_ligne'],i]=1

# Remplissage du bloc correspondant à l'estimation des retards troposphériques humides (estimation horaire)
start = df_rnx_new.index[0][0]
delta_T = pd.Timedelta(days=0, hours=2, minutes=0)
delta_sec = pd.Timedelta(seconds=1)
end   = start + delta_T - delta_sec

import math

delta = df_rnx_new.index[-1][0] - df_rnx_new.index[0][0]
nb_par_zwd = math.ceil(delta / pd.Timedelta(hours=2))


# Squelette du bloc à concaténer à la matrice modèle A
block_zwd = np.zeros((len(df_rnx_new), nb_par_zwd))
ind_c=0

while start <= df_rnx_new.index[-1][0]:
    
    end   = start + delta_T - delta_sec
    extract2c = df_rnx_new.loc[(slice(start, end)), ('ind_ligne','mfw')]   
    block_zwd[extract2c['ind_ligne'],ind_c]=extract2c['mfw']
    
    start = start + delta_T
    ind_c = ind_c + 1


# Initialisation des coordonnées approximatives du récepteur
P_app = np.array([0, 0, 0])

dP_est=np.array([100, 100, 100])
i=1
# Itération pour l'affinement de la position du récepteur
while np.linalg.norm(dP_est[0:3])>1:

    # Calcul des distances approximatives satellite-récepteur
    distances = np.sqrt((df_Sagnac_new['X_sat'].values - P_app[0])**2 +
                        (df_Sagnac_new['Y_sat'].values - P_app[1])**2 +
                        (df_Sagnac_new['Z_sat'].values - P_app[2])**2)

    # Vecteur des observations corrigées des différents modèles
    B = df_rnx_new['L3'].values - distances + conv.SPEED_OF_LIGHT * (df_rnx_new['dte_sat'].values \
    + df_rnx_new['dRelat'].values ) \
    - df_rnx_new['ZHD'].values* df_rnx_new['mfh'].values

    # Construction de la matrice des dérivées partielles
    df_dX = (P_app[0] - df_Sagnac_new['X_sat'].values) / distances
    df_dY = (P_app[1] - df_Sagnac_new['Y_sat'].values) / distances
    df_dZ = (P_app[2] - df_Sagnac_new['Z_sat'].values) / distances
    A = np.column_stack((df_dX, df_dY, df_dZ, block_dt_r, block_zwd, block_one_amb))

    # Résolution par moindres carrés pour estimer le déplacement
    dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
    
    # Mise à jour de la position estimée
    # Attention : on a estimé nb_epochs paramètres d'horloge récepteur
    P_est = np.zeros(len(dP_est))
    
    
    P_est[0:3] = P_app[0:3]+dP_est[0:3]
    P_est[3:]  = dP_est[3:]
    

    
    # Affichage de la position estimée à chaque itération
    print(f"Iteration {i}: Position estimée - X: {P_est[0]}, Y: {P_est[1]}, Z: {P_est[2]}")
    P_app = P_est  # Mise à jour de la position approximative pour la prochaine itération
    i+=1


# Calcul de la distance finale entre la position estimée et la position initiale du header RINEX
dist_P_est_P_rnx_header = np.sqrt(np.sum((P_est[:3] - P_rnx_header)**2))
print("\n")
print("Distance entre la position estimée et la position initiale du header RINEX:", dist_P_est_P_rnx_header)

# Calculer les coordonnées ENU locales
E, N, U = conv.xyz2enu(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
print("Est (E):", E)
print("Nord (N):", N)
print("Haut (U):", U)
print("\n")

# del E, N, U, P_est, dP_est, P_app, A, block_dt_r, B, df_dX, df_dY, df_dZ, distances, i, epoch, epoch_uniques, nb_epochs

# introduire une visualisation du résultat avec folium

# plot_residual_analysis(A, B, dP_est, figure_title="Calcul corr Sat Rec et Sagnac et iono et tropo", save_path="./corr_sat_clock_sagnac_iono_tropo.png",
#                           P_est=P_est[:3], P_rnx_header=P_rnx_header)


# %%
f1 = conv.L1_CARRIER_FREQUENCY
f2 = conv.L2_CARRIER_FREQUENCY

l1 = conv.L1_WAVELENGTH
l2 = conv.L2_WAVELENGTH

# On applique la formule de la combinaison linéaire
# l3 = (f1**2*l1- f2**2*l2)/(f1**2-f2**2);
# qui se simplifie par:
l3=conv.SPEED_OF_LIGHT/(f1+f2)


# attention : on se rend compte que les mesures L1 et L2 sont en cycle alors que l'on
# a besoin d'une mesure de distance (ambigue mais distance quand même)

df_rnx_new['L3']  = (f1**2*df_rnx_new['L1'] - f2**2*df_rnx_new['L2'])/(f1**2-f2**2)
df_rnx_new['P3']  = (f1**2*df_rnx_new['C1'] - f2**2*df_rnx_new['P2'])/(f1**2-f2**2)


gnss_edu.plot_series(df=df_rnx_new, col1='L3', col2='P3' , coeff1=1.0, coeff2=1.0, seuil=3600, renderer="browser")

# %%
# calcul d'autres combinaisons linéaires

df_rnx_new['Lw'] = l2/(l2-l1)*df_rnx_new['L1'] - l1/(l2-l1)*df_rnx_new['L2']
df_rnx_new['Pw'] = l2/(l2-l1)*df_rnx_new['C1'] - l1/(l2-l1)*df_rnx_new['P2']

lw = l1*l2/(l2-l1)


df_rnx_new['Ln'] = l2/(l2+l1)*df_rnx_new['L1'] + l1/(l2+l1)*df_rnx_new['L2']
df_rnx_new['Pn'] = l2/(l2+l1)*df_rnx_new['C1'] + l1/(l2+l1)*df_rnx_new['P2']

ln = l1*l2/(l2+l1)


df_rnx_new['Lmw'] = lw*(df_rnx_new['L1']/l1 - df_rnx_new['L2']/l2) - (f1*df_rnx_new['C1']+f2*df_rnx_new['P2'])/(f1+f2)

lmw = 0.86
df_rnx_new['Lnew'] =  df_rnx_new['L3'] - (f2/(f1+f2)) * df_rnx_new['Lmw']
# %%


fig = gnss_edu.plot_series(df=df_rnx_new, col1='Lmw', col2=None , coeff1=1.0 , coeff2=1.0, seuil=3600, renderer="browser")
