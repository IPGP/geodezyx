#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  1 23:16:42 2025

@author: Samuel Nahmani (1,2)
https://www.ipgp.fr/annuaire/nahmani/)
contact : nahmani@ipgp.fr ou samuel.nahmani@ign.fr
(1) Université Paris Cité, Institut de physique du globe de Paris, CNRS, IGN, F-75005 Paris, France.
(2) Univ Gustave Eiffel, ENSG, IGN, F-77455 Marne-la-Vallée, France. 

Version: 1.0
Dépendances: pandas, numpy, geodezyx, datetime, gpsdatetime, gnsstoolbox

"""

#%%
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

#%%
# gpsdatetime
# Python GPS date/time management package
# Copyright (C) 2014-2023, Jacques Beilin / ENSG-Geomatique
# Distributed under terms of the CECILL-C licence.
#%%
# GnssToolbox - Python package for GNSS learning
# Copyright (C) 2014-2023, Jacques Beilin / ENSG-Geomatique
# Distributed under terms of the CECILL-C licence.

#%%
# GeodeZYX Toolbox’s - [Sakic et al., 2019]
import geodezyx
import geodezyx.conv as conv                  # Import the conversion module
import datetime as dt

#
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
from PPMD_mon_module import *

# pour chercher dans des dossiers
from pathlib import Path
import matplotlib.dates as mdates

#%%
# Chargement des fichiers RINEX d'observation
# <10 km de l'Etna : EMFN ESPC EINT ESLN ECPN EPLU ECNE DPDN
# <10-20 km du volcan : ECRI EMFN EBAG ESPC ESLN EMGL EMAL EDAM EMCN ECOR
# + EDEN station qui est loin de l'Etna


repertoire = Path("/home/snahmani/Traitements/IAG_WG4.3.9_IAVCEI_ETNA/ETNA_2018/GNSS_2018/rinex_2018")

# Filtrer et trier les fichiers
fichiers_sta = sorted(
    [f for f in repertoire.rglob("*ebag3610*.??o") if f.is_file() and f.stat().st_size > 0],
    key=lambda f: f.name
)

fichiers_sta = [str(f) for f in fichiers_sta]  # convertir en strings

# Lire les fichiers un par un
dfs = []
for fichier in fichiers_sta:
    print(f"Lecture de : {fichier}")
    try:
        df = geodezyx.files_rw.read_rinex2_obs(fichier, set_index=['epoch', 'prn'])  # <== [fichier] ici
        dfs.append(df)
    except Exception as e:
        print(f"Erreur lors de la lecture de {fichier} : {e}")
        
# Position approchée lue dans le header du fichier RINEX
columns =  grep_file(r"APPROX POSITION XYZ", fichier)
if columns:
    # Séparer la ligne en supprimant le texte "APPROX POSITION XYZ"
    valeurs = columns[0].split()[:3]  # Récupérer uniquement les 3 premières valeurs
    P_rnx_header = np.array(valeurs, dtype=float)  # Convertir en numpy array
    print("Coordonnées XYZ :", P_rnx_header)
else:
    print("Aucune correspondance trouvée.")

# Concaténer tous les DataFrames
if dfs:
    df_total = pd.concat(dfs)
    print("✅ Données combinées avec succès")
else:
    print("❌ Aucun fichier valide chargé")

#%%
# nettoyage
df_total = df_total.dropna(axis=1, how='all')

#%%
#
from datetime import datetime, timedelta
import os
import locale

# Définir la locale en anglais pour le formatage des dates
locale.setlocale(locale.LC_TIME, 'en_US.UTF-8')
# Extraire uniquement le nom du fichier (sans le chemin)
filename = os.path.basename(fichier)  # Résultat : 'kokb1430.01o'

try:
    # Extraction des informations DOY et année à partir du nom du fichier
    doy = int(filename[4:7])  # DOY (exemple : 143)
    year_suffix = filename[9:11]  # Suffixe d'année (exemple : '01')

    # Vérification que le suffixe d'année est bien numérique
    if not year_suffix.isdigit():
        raise ValueError(f"Invalid year suffix: {year_suffix}")

    # Conversion du suffixe de l'année en année complète
    if int(year_suffix) < 80:
        year = 2000 + int(year_suffix)
    else:
        year = 1900 + int(year_suffix)

    # Conversion du DOY en date complète
    date = datetime(year, 1, 1) + timedelta(days=doy - 1)

    # Formatage en date américaine
    date_american = date.strftime("%B %d, %Y")  # Exemple : "May 23, 2001"
    
    station = filename[0:4]

    # Créer le titre
    titre_generique = f"{station.upper()} - {date_american}"
    print(titre_generique)

except (ValueError, IndexError) as e:
    print(f"Error processing file '{filename}': {e}")


#%% Chargement des fichiers de navigation

repertoire = Path("/home/snahmani/Traitements/IAG_WG4.3.9_IAVCEI_ETNA/ETNA_2018/GNSS_2018/orb_2018")


# Filtrer et trier les fichiers
fichiers_sp3= sorted(
    [f for f in repertoire.rglob("igs*.sp3") if f.is_file() and f.stat().st_size > 0],
    key=lambda f: f.name
)

fichiers_sp3 = [str(f) for f in fichiers_sp3]  # convertir en strings

mysp3 = orb.orbit()
mysp3.loadSp3(fichiers_sp3)

#%%
#Il faut calculer la position de chaque satellite GNSS à chaque temps d'émission
t = gpst.gpsdatetime()

X_sat = []
Y_sat = []
Z_sat = []
dte_sat = []
dRelat = []

for (time_i,prn_i) in df_total.index:

    t.rinex_t(time_i.to_pydatetime().strftime('%y %m %d %H %M %S.%f'))
    
    t_emission_mjd  = t.mjd - df_total.loc[(time_i,prn_i), 'C1'] / gnss_const.c / 86400.0
    
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
    
    
df_total['X_sat']   = X_sat 
df_total['Y_sat']   = Y_sat 
df_total['Z_sat']   = Z_sat
df_total['dte_sat'] = dte_sat
df_total['dRelat']  = dRelat


del X_sat, Y_sat, Z_sat, dte_sat, time_i, prn_i, t_emission_mjd, t, X_sat_v, Y_sat_v, Z_sat_v, dte_sat_v, dRelat_v, dRelat

#%%


Az_rad, Ele_rad = tools.toolAzEle(P_rnx_header[0],P_rnx_header[1],P_rnx_header[2],df_total.X_sat,df_total.Y_sat,df_total.Z_sat)
# radian to degree
rad2deg = 180 / np.pi

Ele_deg =  Ele_rad * rad2deg
Az_deg  =  Az_rad * rad2deg

# Ajout des colonnes au DataFrame
df_total['Az_deg'] = Az_deg
df_total['Ele_deg'] = Ele_deg

#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Paramètres
gap_threshold = pd.Timedelta(minutes=15)

# Liste des PRNs uniques
prns = df_total.index.get_level_values('prn').unique()

# Création du plot
fig, ax = plt.subplots(figsize=(12, 6))

for prn in prns:
    # Extraire les données pour ce PRN
    data = df_total.xs(prn, level='prn').sort_index()
    
    # Récupérer les timestamps (epochs)
    times = data.index
    values = data['S2'].values

    # Calculer les deltas de temps
    time_deltas = times.to_series().diff()

    # Identifier les ruptures > 15 min
    rupture_idx = time_deltas > gap_threshold

    # On insère NaN juste après les ruptures
    times_with_nan = []
    values_with_nan = []

    for i in range(len(times)):
        times_with_nan.append(times[i])
        values_with_nan.append(values[i])
        
        # Si rupture après ce point → insérer NaN
        if i + 1 < len(times) and rupture_idx.iloc[i + 1]:
            times_with_nan.append(times[i] + pd.Timedelta(seconds=30))  # ajouter un point "vide"
            values_with_nan.append(np.nan)

    # Tracer la courbe
    ax.plot(times_with_nan, values_with_nan, label=prn)

# Mise en forme
ax.set_title("Séries temporelles S2 (L2) par PRN – avec coupures > 15 min")
ax.set_xlabel("Temps")
ax.set_ylabel("S2")
ax.legend(title="PRN")
#%%

import matplotlib.pyplot as plt
import numpy as np

# Extraction pour 'G10' sur une période spécifique avec df_index
start_period = pd.Timestamp('2018-12-27 00:00:00')
end_period = pd.Timestamp('2018-12-27 23:59:30')
extract2b = df_total.loc[(slice(start_period, end_period)), ['Az_deg', 'Ele_deg', 'S2']]

'Az_deg', 'Ele_deg',

# ⚠️ Filtrer les données valides
df_valid = extract2b

# Conversion en radians pour le tracé polaire
az = np.deg2rad(df_valid['Az_deg'].values)
ele = df_valid['Ele_deg'].values

# On convertit l’élévation pour que le zénith soit au centre
r = 90 - ele  # rayon = 0 au zénith, 90 à l’horizon

# Valeur à afficher (par exemple S2)
val = df_valid['S2'].values

# Création du skyplot
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, polar=True)

# Scatter polar avec couleurs selon la valeur
sc = ax.scatter(az, r, c=val, cmap='viridis', s=10, alpha=0.8)

# Personnalisation du plot
ax.set_theta_zero_location("N")  # 0° en haut
ax.set_theta_direction(-1)       # azimut dans le sens horaire
ax.set_rlim(0, 90)               # rayon de 0 (zénith) à 90 (horizon)
ax.set_title("Skyplot des satellites – coloration par S2", va='bottom')

# Barre de couleur
cbar = plt.colorbar(sc, ax=ax, orientation='vertical', pad=0.1)
cbar.set_label("S2")

plt.tight_layout()
plt.show()




#%%
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, polar=True)

sc = ax.scatter(az, r, c=val, cmap='viridis', s=10, alpha=0.8)

ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.set_rlim(0, 90)
ax.set_title("Skyplot des satellites – coloration par S2", va='bottom', fontsize=12)

cbar = plt.colorbar(sc, ax=ax, orientation='vertical', pad=0.1)
cbar.set_label("S2")

# Ajouter une légende descriptive avec station + période
start_period = pd.Timestamp('2018-12-24 09:00:00')
end_period = pd.Timestamp('2018-12-24 17:00:00')
info_text = f"Station : EBAG\nPériode : {start_period.strftime('%Y-%m-%d %H:%M')} → {end_period.strftime('%H:%M')}"
fig.text(0.02, 0.98, info_text, fontsize=10, ha='left', va='top', transform=fig.transFigure)

plt.tight_layout()
plt.show()


#%%

import matplotlib.pyplot as plt
import numpy as np

# Extraction et préparation des données
df_reset = df_valid.reset_index()
df_reset_sorted = df_reset.sort_values('epoch')
first_points = df_reset_sorted.groupby('prn').first()

# Skyplot - coordonnées polaires
az = np.radians(df_reset['Az_deg'].values)
r = 90 - df_reset['Ele_deg'].values
val = df_reset['S2'].values

# Tracé du skyplot
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, polar=True)

sc = ax.scatter(az, r, c=val, cmap='viridis', s=10, alpha=0.8)

# Orientation classique : Nord en haut, rotation horaire
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.set_rlim(0, 90)
ax.set_title("Skyplot des satellites – coloration par S2", va='bottom', fontsize=12)

# Affichage du nom des satellites au début de leur arc
for prn, row in first_points.iterrows():
    az_rad = np.radians(row['Az_deg'])
    r_val = 90 - row['Ele_deg']
    ax.text(az_rad, r_val, prn, fontsize=8, fontweight='bold',
            ha='center', va='center', color='black', clip_on=True)

# Barre de couleur
cbar = plt.colorbar(sc, ax=ax, orientation='vertical', pad=0.1)
cbar.set_label("S2")

# Informations station/période
start_period = df_reset['epoch'].min()
end_period = df_reset['epoch'].max()
info_text = f"Station : EBAG\nPériode : {start_period.strftime('%Y-%m-%d %H:%M')} → {end_period.strftime('%H:%M')}"
fig.text(0.02, 0.98, info_text, fontsize=10, ha='left', va='top', transform=fig.transFigure)

plt.tight_layout()
plt.show()


#%%


import matplotlib.pyplot as plt
import numpy as np

# Extraction pour 'G10' sur une période spécifique avec df_index
start_period = pd.Timestamp('2018-12-27 16:00:00')
end_period = pd.Timestamp('2018-12-27 23:59:30')
extract2b = df_total.loc[(slice(start_period, end_period)), ['Az_deg', 'Ele_deg', 'S2']]

'Az_deg', 'Ele_deg',

# ⚠️ Filtrer les données valides
df_valid = extract2b

# Réinitialiser et trier le DataFrame
df_reset = df_valid.reset_index()
df_sorted = df_reset.sort_values('epoch')

# Récupérer le premier et le dernier point pour chaque PRN
first_points = df_sorted.groupby('prn').first()
last_points = df_sorted.groupby('prn').last()

# Skyplot - coordonnées polaires
az = np.radians(df_reset['Az_deg'].values)
r = 90 - df_reset['Ele_deg'].values
val = df_reset['S2'].values

# Tracé du skyplot
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, polar=True)

sc = ax.scatter(az, r, c=val, cmap='viridis', s=10, alpha=0.8)

# Orientation classique
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.set_rlim(0, 90)
ax.set_title("Skyplot des satellites – coloration par S2", va='bottom', fontsize=12)

# ➕ Étiquette au début de l’arc
for prn, row in first_points.iterrows():
    az_rad = np.radians(row['Az_deg'])
    r_val = 90 - row['Ele_deg']
    ax.text(az_rad, r_val, prn, fontsize=8, fontweight='bold',
            ha='center', va='center', color='black', clip_on=True)

# ➕ Étiquette à la fin de l’arc
for prn, row in last_points.iterrows():
    az_rad = np.radians(row['Az_deg'])
    r_val = 90 - row['Ele_deg']
    ax.text(az_rad, r_val, prn, fontsize=8, fontweight='bold',
            ha='center', va='center', color='black', clip_on=True)

# Barre de couleur
cbar = plt.colorbar(sc, ax=ax, orientation='vertical', pad=0.1)
cbar.set_label("S2")

# Informations station/période
start_period = df_reset['epoch'].min()
end_period = df_reset['epoch'].max()
info_text = f"Station : EBAG\nPériode : {start_period.strftime('%Y-%m-%d %H:%M')} → {end_period.strftime('%H:%M')}"
fig.text(0.02, 0.98, info_text, fontsize=10, ha='left', va='top', transform=fig.transFigure)

plt.tight_layout()
plt.show()


