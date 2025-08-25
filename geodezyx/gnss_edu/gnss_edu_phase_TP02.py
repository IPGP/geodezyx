#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 16:44:03 2024

Filière ING3 - PPMD - Traitement de la mesure de phase

@author: Samuel Nahmani (1,2)
https://www.ipgp.fr/annuaire/nahmani/)
contact : nahmani@ipgp.fr ou samuel.nahmani@ign.fr
(1) Université Paris Cité, Institut de physique du globe de Paris, CNRS, IGN, F-75005 Paris, France.
(2) Univ Gustave Eiffel, ENSG, IGN, F-77455 Marne-la-Vallée, France. 

Version: 1.0
Dépendances: pandas, numpy, geodezyx, datetime, gpsdatetime, gnsstoolbox

"""

# %%
# GeodeZYX Toolbox’s
# [Sakic et al., 2019]
# Sakic, Pierre; Mansur, Gustavo; Chaiyaporn, Kitpracha; Ballu, Valérie (2019):
# The geodeZYX toolbox: a versatile Python 3 toolbox for geodetic-oriented purposes. v. 4.0.
# GFZ Data Services. http://doi.org/10.5880/GFZ.1.1.2019.002
#
# Documentation
# https://geodezyx.github.io/geodezyx-toolbox/
#
# Installation
# pip install git+https://github.com/GeodeZYX/geodezyx-toolbox
# pip uninstall geodezyx

# %%
# gpsdatetime
# Python GPS date/time management package
# Copyright (C) 2014-2023, Jacques Beilin / ENSG-Geomatique
# Distributed under terms of the CECILL-C licence.
# %%
# GnssToolbox - Python package for GNSS learning
# Copyright (C) 2014-2023, Jacques Beilin / ENSG-Geomatique
# Distributed under terms of the CECILL-C licence.

# %%
# GeodeZYX Toolbox’s - [Sakic et al., 2019]
import geodezyx
import geodezyx.conv as conv                  # Import the conversion module
import datetime as dt
import geodezyx.files_rw as files_rw         # Import the file reading/writing module
import geodezyx.reffram as reffram         # Import the operational module
from geodezyx import utils


import gpsdatetime as gpst

import gnsstoolbox.orbits as orb
import gnsstoolbox.gnss_const as gnss_const
import gnsstoolbox.gnsstools as tools
import gnsstoolbox.gnss_process as proc
import gnsstoolbox.gnss_corr as gnss_corr

import pandas as pd
import numpy as np

# pour visualiser les données
import matplotlib.pyplot as plt

# mon module
from gnss_edu import *

# %%
# Chargement des fichiers RINEX d'observation
# fichier_rnx='./data/data-2019/mlvl176z.18o'
# fichier_rnx='./data/data-2019/mlvl1760.18o'

fichier_rnx='./data/data-2019/mlvl1760.18d.Z'

# Chargement des données RINEX d'observation dans un pandas dataframe via  GeodeZYX
df_rnx_flat, l_rnx_head = geodezyx.files_rw.read_rinex2_obs(fichier_rnx,
                                           return_header=True)

df_rnx = df_rnx_flat.set_index(['epoch', 'prn'], drop=True)


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

l1 = gnss_const.lambda1
l2 = gnss_const.lambda2

# il faut consulter la doc et vous constaterez que les unités de L1 et L2 sont en cycle
df_rnx['L1'] = df_rnx['L1']*l1 
df_rnx['L2'] = df_rnx['L2']*l2

# attention : on se rend compte que les mesures L1 et L2 sont en cycle alors que l'on
# a besoin d'une mesure de distance (ambigue mais distance quand même)


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


# # Chargement des données RINEX d'observation via  GnssToolbox
# # Trop de lourdeur pour pas grand chose ... autant passer par un équivalent de grep
# import gnsstoolbox.rinex_o as rx
# my_rnx =  rx.rinex_o()
# my_rnx.loadRinexO(fichier_rnx)

# # Position initiale du récepteur (à partir du header RINEX)
# P_rnx_header = np.array([my_rnx.headers[0].X, my_rnx.headers[0].Y, my_rnx.headers[0].Z])
# del my_rnx


# %%
# Chargement des fichiers d'orbites
fichier_sp3  = ['./data/data-2019/igs20071.sp3.Z', './data/data-2019/igs20072.sp3.Z']
fichier_brdc = './data/data-2019/mlvl176z.18n'

mysp3 = orb.orbit()
mysp3.loadSp3(fichier_sp3)

dforb = pd.concat([files_rw.read_sp3(f) for f in fichier_sp3])

mynav = orb.orbit()
mynav.loadRinexN('./data/data-2019/mlvl176z.18n')

# Il faut calculer la position de chaque satellite GNSS à chaque temps d'émission
t = gpst.gpsdatetime()

X_sat = []
Y_sat = []
Z_sat = []
dte_sat = []
dRelat = []


epochs_rec = df_rnx_flat["epoch"].unique()

from geodezyx import  reffram
dforb_itrp = reffram.orb_df_lagrange_interpolate(dforb, epochs)
dforb_itrp_aft = reffram.orb_df_lagrange_interpolate(dforb, epochs + dt.timedelta(milliseconds=1))
dforb_itrp_bef = reffram.orb_df_lagrange_interpolate(dforb, epochs - dt.timedelta(milliseconds=1))

dforb_itrp_aft[["x","y","z"]] - dforb_itrp_bef[["x","y","z"]] / 2.0 / 1e3


for (time_i,prn_i) in df_rnx.index:

    t.rinex_t(time_i.to_pydatetime().strftime('%y %m %d %H %M %S.%f'))

    t_emission_mjd  = t.mjd - df_rnx.loc[(time_i,prn_i), 'C1'] / gnss_const.c / 86400.0

    (X_sat_v,Y_sat_v,Z_sat_v,dte_sat_v)	 = mysp3.calcSatCoord(prn_i[0], int(prn_i[1:]),t_emission_mjd)

    # calcul de l'effet relativiste
    delta_t = 1e-3 # écart de temps en +/- pour calculer la dérivée 
    (Xs1,Ys1,Zs1,clocks1) = mysp3.calcSatCoord(prn_i[0], int(prn_i[1:]),t_emission_mjd - delta_t / 86400.0)    
    (Xs2,Ys2,Zs2,clocks2) = mysp3.calcSatCoord(prn_i[0], int(prn_i[1:]),t_emission_mjd + delta_t / 86400.0)  

    VX      = (np.array([Xs2-Xs1, Ys2-Ys1, Zs2-Zs1]))/2.0/delta_t
    VX0     = np.array([X_sat_v,Y_sat_v,Z_sat_v])

    dRelat_v  = -2.0 * VX0.T@ VX /(gnss_const.c **2)

    # temps d'emission du signal GNSS en temps GNSS (mjd)
    t_emission_mjd = t_emission_mjd - dte_sat_v / 86400.0 - dRelat_v / 86400.0

    # Recalcul de la position du satellite au temps d'emission (temps GNSS en mjd)
    (X_sat_v,Y_sat_v,Z_sat_v,dte_sat_v)	 = mysp3.calcSatCoord(prn_i[0], int(prn_i[1:]),t_emission_mjd)

    X_sat.append(X_sat_v)
    Y_sat.append(Y_sat_v)
    Z_sat.append(Z_sat_v)
    dte_sat.append(dte_sat_v)
    dRelat.append(dRelat_v)


df_rnx['X_sat']   = X_sat 
df_rnx['Y_sat']   = Y_sat 
df_rnx['Z_sat']   = Z_sat
df_rnx['dte_sat'] = dte_sat
df_rnx['dRelat']  = dRelat


del X_sat, Y_sat, Z_sat, dte_sat, time_i, prn_i, t_emission_mjd, t, X_sat_v, Y_sat_v, Z_sat_v, dte_sat_v, dRelat_v, dRelat


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
    dP_est = np.linalg.inv(A.T@A)@A.T@B
    #dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
    
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
    E, N, U = tools.toolCartLocGRS80(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
    print("Est (E):", E)
    print("Nord (N):", N)
    print("Haut (U):", U)
    print("\n")


fig = plot_residual_analysis(A, B, dP_est, figure_title="Calcul Trivial", save_path="./trivial_bis.png",
                           P_est=P_est, P_rnx_header=P_rnx_header, tools=tools);


# del E, N, U, P_est, dP_est, P_app, A, B, df_dX, df_dY, df_dZ, distances, i

# %%
# Traitement classique sur le code pour la position d'un récepteur GNSS
# corrections des erreurs d'horloges satellites
# -> elles sont connues -> correction directe du vecteur B

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
    B = df_rnx['C1'].values - distances + gnss_const.c * (df_rnx['dte_sat'].values + df_rnx['dRelat'].values )

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
E, N, U = tools.toolCartLocGRS80(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
print("Est (E):", E)
print("Nord (N):", N)
print("Haut (U):", U)
print("\n")


plot_residual_analysis(A, B, dP_est, figure_title="Prise en compte des horloges satellites", save_path="./sat_clocks_only.png",
                           P_est=P_est, P_rnx_header=P_rnx_header, tools=tools)


# del E, N, U, P_est, dP_est, P_app, A, B, df_dX, df_dY, df_dZ, distances, i

# %%
# Traitement classique sur le code pour la position d'un récepteur GNSS
# corrections des erreurs d'horloges satellites
# -> elles sont connues -> correction directe du vecteur B

print('*****  Prise en compte des erreurs d''horloge satellites et récepteur *****')

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
    B = df_rnx['C1'].values - distances + gnss_const.c * (df_rnx['dte_sat'].values + df_rnx['dRelat'].values )

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
E, N, U = tools.toolCartLocGRS80(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
print("Est (E):", E)
print("Nord (N):", N)
print("Haut (U):", U)
print("\n")


plot_residual_analysis(A, B, dP_est, figure_title="Calcul corr Sat", save_path="./corr_sat_clock.png",
                           P_est=P_est, P_rnx_header=P_rnx_header, tools=tools)


del E, N, U, P_est, dP_est, P_app, A, B, df_dX, df_dY, df_dZ, distances, i


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
    alpha_rad = row['C1'] / gnss_const.c * gnss_const.Omega_e
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
    B = df_rnx['C1'].values - distances + gnss_const.c * (df_rnx['dte_sat'].values + df_rnx['dRelat'].values )

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
E, N, U = tools.toolCartLocGRS80(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
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
    B = df_rnx['C1'].values - distances + gnss_const.c * (df_rnx['dte_sat'].values + df_rnx['dRelat'].values )

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
E, N, U = tools.toolCartLocGRS80(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
print("Est (E):", E)
print("Nord (N):", N)
print("Haut (U):", U)
print("\n")

plot_residual_analysis(A, B, dP_est, figure_title="Calcul corr Sat Rec et Sagnac", save_path="./corr_sat_clock_sagnac.png",
                           P_est=P_est[:3], P_rnx_header=P_rnx_header, tools=tools)


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

print('*****  Prise en compte des erreurs d''horloge satellites *****')
print('*****  + Sagnac *****')
print('*****  + horloges récepteur *****')
print('*****  + iono Klobuchar *****')

import klobuchar
# Calcul des Az et des Ele de chaque mesure
# import gnsstoolbox.gnsstools as tools
# toolAzEle : azimut and elevation (radians) for one or several satellites Xs,Ys,Zs (scalar or vector) seen from a point with X,Y,Z coordinates.

Az_rad, Ele_rad = tools.toolAzEle(P_rnx_header[0],P_rnx_header[1],P_rnx_header[2],df_rnx.X_sat,df_rnx.Y_sat,df_rnx.Z_sat)

# toolCartGeoGRS80 : cartesian to geographic coordinates conversion. All angles are given in radians.
lon,lat,h = tools.toolCartGeoGRS80(P_rnx_header[0],P_rnx_header[1],P_rnx_header[2])

# radian to degree
rad2deg = 180 / np.pi
lon_d = lon * rad2deg
lat_d = lat * rad2deg

Az_deg = Az_rad * rad2deg
Ele_deg =  Ele_rad * rad2deg

# coefficients alpha et beta extrait du fichier de navigation
alpha = mynav.ion_alpha_gps
beta = mynav.ion_beta_gps

t = gpst.gpsdatetime() # création d'un objet temps de la gnsstime
dIon1 = []
for (time_i,prn_i) in df_rnx.index:
    t.rinex_t(time_i.to_pydatetime().strftime('%y %m %d %H %M %S.%f'))
    
    wsec_v = t.wsec       # seconds in GPS week
    i = df_rnx.loc[(time_i,prn_i), 'ind_ligne'] 
    dIon1_v = klobuchar.klobuchar(lat_d, lon_d, Ele_deg[i], Az_deg[i] , wsec_v, alpha, beta )
    dIon1.append(dIon1_v)

df_rnx['Az']      = Az_deg
df_rnx['Ele']     = Ele_deg
df_rnx['dIon1']   = dIon1 

del Az_rad, Ele_rad, lon, lat, lat_d, lon_d, h, rad2deg, Az_deg, Ele_deg, alpha, beta, t, dIon1, wsec_v, i, dIon1_v, time_i, prn_i


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
    B = df_rnx['C1'].values - distances + gnss_const.c * (df_rnx['dte_sat'].values + df_rnx['dRelat'].values ) - df_rnx['dIon1'].values

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
E, N, U = tools.toolCartLocGRS80(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
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

f1 = gnss_const.f1  
f2 = gnss_const.f2

l1 = gnss_const.lambda1
l2 = gnss_const.lambda2

# On applique la formule de la combinaison linéaire
# l3 = (f1**2*l1- f2**2*l2)/(f1**2-f2**2);
# qui se simplifie par:
l3=gnss_const.c/(f1+f2)


# attention : on se rend compte que les mesures L1 et L2 sont en cycle alors que l'on
# a besoin d'une mesure de distance (ambigue mais distance quand même)

df_rnx['L3']  = (f1**2*df_rnx['L1'] - f2**2*df_rnx['L2'])/(f1**2-f2**2);
df_rnx['P3']  = (f1**2*df_rnx['C1'] - f2**2*df_rnx['P2'])/(f1**2-f2**2);


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

#     t_emission_mjd  = t.mjd - df_rnx.loc[(time_i,prn_i), 'P3'] / gnss_const.c / 86400.0

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
    B = df_rnx['P3'].values - distances + gnss_const.c * (df_rnx['dte_sat'].values + df_rnx['dRelat'].values )

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
E, N, U = tools.toolCartLocGRS80(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
print("Est (E):", E)
print("Nord (N):", N)
print("Haut (U):", U)
print("\n")


plot_residual_analysis(A, B, dP_est, figure_title="Calcul corr Sat Rec et Sagnac et iono", save_path="./corr_sat_clock_sagnac_iono.png",
                           P_est=P_est[:3], P_rnx_header=P_rnx_header, tools=tools)

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
    alpha_rad = row['C1'] / gnss_const.c * gnss_const.Omega_e
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
    B = df_rnx['P3'].values - distances + gnss_const.c * (df_rnx['dte_sat'].values + df_rnx['dRelat'].values )

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
E, N, U = tools.toolCartLocGRS80(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
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


import gpt3 as gpt3

# toolCartGeoGRS80 : cartesian to geographic coordinates conversion. All angles are given in radians.
lon, lat, h = tools.toolCartGeoGRS80(P_rnx_header[0],P_rnx_header[1],P_rnx_header[2])

# grille gpt3
gpt3_5 = gpt3.gpt3_5_fast_readGrid(filename ='gpt3_5.grd')


# Il faut calculer les corrections et les estimations pour chaque satellite et
# pour chaque époque
t = gpst.gpsdatetime()

ZHD = []
mfh = []
mfw = []

for (time_i,prn_i) in df_rnx.index:

    t.rinex_t(time_i.to_pydatetime().strftime('%y %m %d %H %M %S.%f'))
    #gmfh, gmfw = gpt3.gmf(dmjd=t.mjd, dlat=lat, dlon=lon, dhgt=h, zd =np.pi/2-df_rnx.loc[(time_i,prn_i), 'Ele']*np.pi/180)
    gmfh = 1/np.sin(df_rnx.loc[(time_i,prn_i), 'Ele']*np.pi/180)
    gmfw = gmfh
    T = gpt3.gpt3_5_fast(mjd=t.mjd, lat=np.array([lat]), lon=np.array([lon]), h_ell=np.array([h]), it=0, grid=gpt3_5)
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

Az_rad, Ele_rad = tools.toolAzEle(P_rnx_header[0],P_rnx_header[1],P_rnx_header[2],df_rnx.X_sat,df_rnx.Y_sat,df_rnx.Z_sat)

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
    B = df_rnx_new['P3'].values - distances + gnss_const.c * (df_rnx_new['dte_sat'].values \
    + df_rnx_new['dRelat'].values ) \
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
E, N, U = tools.toolCartLocGRS80(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
print("Est (E):", E)
print("Nord (N):", N)
print("Haut (U):", U)
print("\n")

# del E, N, U, P_est, dP_est, P_app, A, block_dt_r, B, df_dX, df_dY, df_dZ, distances, i, epoch, epoch_uniques, nb_epochs

# introduire une visualisation du résultat avec folium

plot_residual_analysis(A, B, dP_est, figure_title="Calcul corr Sat Rec et Sagnac et iono et tropo", save_path="./corr_sat_clock_sagnac_iono_tropo.png",
                           P_est=P_est[:3], P_rnx_header=P_rnx_header, tools=tools);

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
    B = df_rnx_new['L3'].values - distances + gnss_const.c * (df_rnx_new['dte_sat'].values \
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
E, N, U = tools.toolCartLocGRS80(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
print("Est (E):", E)
print("Nord (N):", N)
print("Haut (U):", U)
print("\n")

# del E, N, U, P_est, dP_est, P_app, A, block_dt_r, B, df_dX, df_dY, df_dZ, distances, i, epoch, epoch_uniques, nb_epochs

# introduire une visualisation du résultat avec folium

# plot_residual_analysis(A, B, dP_est, figure_title="Calcul corr Sat Rec et Sagnac et iono et tropo", save_path="./corr_sat_clock_sagnac_iono_tropo.png",
#                           P_est=P_est[:3], P_rnx_header=P_rnx_header, tools=tools)


# %%


f1 = gnss_const.f1  
f2 = gnss_const.f2

l1 = gnss_const.lambda1
l2 = gnss_const.lambda2

# On applique la formule de la combinaison linéaire
# l3 = (f1**2*l1- f2**2*l2)/(f1**2-f2**2);
# qui se simplifie par:
l3=gnss_const.c/(f1+f2)


# attention : on se rend compte que les mesures L1 et L2 sont en cycle alors que l'on
# a besoin d'une mesure de distance (ambigue mais distance quand même)

df_rnx_new['L3']  = (f1**2*df_rnx_new['L1'] - f2**2*df_rnx_new['L2'])/(f1**2-f2**2);
df_rnx_new['P3']  = (f1**2*df_rnx_new['C1'] - f2**2*df_rnx_new['P2'])/(f1**2-f2**2);


plot_series(df=df_rnx_new, col1='L3', col2='P3' , coeff1=1.0, coeff2=1.0, seuil=3600, renderer="browser")

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


fig = plot_series(df=df_rnx_new, col1='Lmw', col2=None , coeff1=1.0 , coeff2=1.0, seuil=3600, renderer="browser")
