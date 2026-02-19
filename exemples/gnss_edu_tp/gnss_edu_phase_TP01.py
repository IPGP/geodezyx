#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 10:12:46 2024
Filière ING3 - PPMD - Traitement de la mesure de phase
@author: Samuel Nahmani (1,2)
https://www.ipgp.fr/annuaire/nahmani/)
contact : nahmani@ipgp.fr ou samuel.nahmani@ign.fr
(1) Université Paris Cité, Institut de physique du globe de Paris, CNRS, IGN, F-75005 Paris, France.
(2) Univ Gustave Eiffel, ENSG, IGN, F-77455 Marne-la-Vallée, France. 

Version: 1.0
Dépendances: pandas, numpy, geodezyx, datetime

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
# 
from geodezyx import files_rw     # Import the read/write module
from geodezyx import conv         # Import the conversion module
from geodezyx import operational  # Import the download rinex module

                
import datetime as dt
import pandas as pd
import numpy as np

# pour visualiser les données
import matplotlib.pyplot as plt

from pathlib import Path
import os

#%%
# création du dossier gnss_edu_data qui va contenir les données et les résultats du TP
my_directory = os.environ["HOME"] + "/gnss_edu_data/"

# Chemin avec expansion du ~ vers le home
folder = Path(my_directory).expanduser()

# Création du dossier s'il n'existe pas
folder.mkdir(parents=True, exist_ok=True)


#%%
# Téléchargement automatique des données RINEX des stations de SMNE et MLVL distance d'une
# dizaine de kilomètres dans la région de Paris, France sur le serveur IGN (France)
# données pour 1 jour (2019-176) à 30s

# Création d'un datetime pour gérer le jour de traitement et ne pas à avoir à gérer les doy, les jjul etc !
my_date_to_process = dt.datetime(2019,6,25)

dwl_output = operational.download_gnss_rinex(statdico={"rgp" : ["SMNE","MLVL"]},
                                output_dir=my_directory,
                                startdate= my_date_to_process ,
                                enddate= my_date_to_process ,
                                parallel_download = 1) 

#%%
#
fichier_base     =   dwl_output[0][0]
fichier_mobile   =   dwl_output[1][0]


#%% Préambule
# chargement des pandas DataFrame et utilisation
# Lecture des observations RINEX en deux formats de DataFrame différents

df_flat = files_rw.read_rinex_obs(fichier_base)

df_index = files_rw.read_rinex_obs(fichier_base, set_index=['epoch', 'prn'])

# Quelles différences observez-vous entre les deux dataframes df_flat et df_index ?

#df_index['ind_ligne'] = range(len(df_index)) # à décommenter quand c'est compris
# Que s'est il passé sur le dataframe df_index ? Indication, ouvrez le et examiner les colonnes...

#%%
# Selection "FLAT" - Utilisation d'un DataFrame sans index personnalisé
# Avantage : Accès direct aux numéros de lignes pour faciliter la construction 
# de la matrice modèle

# Filtrage des observations pour le satellite 'G10' à un moment précis
my_prn_for_extraction = 'G05'
my_date_for_extraction =  my_date_to_process + dt.timedelta(hours=2, minutes=40, seconds=0)
my_obs_to_extract = 'L1'

print(my_date_for_extraction)

bool_prn = df_flat['prn'] == my_prn_for_extraction
bool_epoch = df_flat['epoch'] == pd.Timestamp(my_date_for_extraction)
extract1a = df_flat.loc[bool_prn & bool_epoch, my_obs_to_extract]
print(extract1a)


# Extraction des observations pour 'G10' sur une période spécifique
my_prn_for_extraction = 'G05'
my_date_for_extraction_start =  my_date_to_process + dt.timedelta(hours=2, minutes=40, seconds=0)
my_date_for_extraction_end =  my_date_to_process + dt.timedelta(hours=2, minutes=45, seconds=30)
my_obs_to_extract = 'L1'

bool_prn = df_flat['prn'] == my_prn_for_extraction
bool_epoch = (df_flat['epoch'] >= pd.Timestamp(my_date_for_extraction_start)) & (df_flat['epoch'] <= pd.Timestamp(my_date_for_extraction_end))
extract1b = df_flat.loc[bool_prn & bool_epoch, my_obs_to_extract]
print(extract1b)


# Extraction des observations pour 'G10' sur une période définie par une date de début et un delta de temps
my_prn_for_extraction = 'G05'
my_date_for_extraction =  my_date_to_process + dt.timedelta(hours=2, minutes=40, seconds=0)
my_obs_to_extract = 'L1'

bool_prn = df_flat['prn'] ==  my_prn_for_extraction
start = pd.Timestamp(my_date_for_extraction)
delta_T = pd.Timedelta(days=0, hours=1, minutes=15) # on peut gérér le delta T par pandas
end = start + delta_T
bool_epoch = (df_flat['epoch'] >= start) & (df_flat['epoch'] <= end)
extract1c = df_flat.loc[bool_prn & bool_epoch, my_obs_to_extract]
print(extract1c)

# Remarquez que l'utilisation des booléens nous permet d'accéder directement aux numéros
# des lignes concernées. Si on veut les numéros de ligne, il suffit de faire : 
serie_bool = bool_prn & bool_epoch ; 
serie_chiffre = serie_bool.astype(int) ; 



#%%
#### selection (multi)-index
# Je fais des filtrages avec .loc 
# Extraction pour le satellite 'G10' à un moment précis avec df_index

my_prn_for_extraction = 'G10'
my_date_for_extraction =  my_date_to_process + dt.timedelta(hours=3, minutes=27, seconds=30)
my_obs_to_extract = 'L1'

extract2a = df_index.loc[(pd.Timestamp(my_date_for_extraction), my_prn_for_extraction), my_obs_to_extract]
print(extract2a)

# Extraction pour 'G10' sur une période spécifique avec df_index
my_prn_for_extraction = 'G10'
my_date_for_extraction_start =  my_date_to_process + dt.timedelta(hours=3, minutes=27, seconds=0)
my_date_for_extraction_end =  my_date_to_process + dt.timedelta(hours=3, minutes=32, seconds=30)
my_obs_to_extract = 'L1'

start_period = pd.Timestamp(my_date_for_extraction_start)
end_period = pd.Timestamp(my_date_for_extraction_end )
extract2b = df_index.loc[(slice(start_period, end_period), my_prn_for_extraction), my_obs_to_extract]
print(extract2b)

# Extraction pour 'G10' sur une période définie par une date de début et un delta de temps avec df_index
my_prn_for_extraction = 'G10'
my_date_for_extraction =  my_date_to_process + dt.timedelta(hours=3, minutes=27, seconds=30)
my_obs_to_extract = 'L1'

start = pd.Timestamp(my_date_for_extraction)
delta_T = pd.Timedelta(days=0, hours=1, minutes=15)
end = start + delta_T
extract2c = df_index.loc[(slice(start, end), my_prn_for_extraction), my_obs_to_extract]
print(extract2c)

#%%
# bloc à commenter une fois que c'est compris.
# Dans ce cas, si je veux récuperer le numéro des lignes concernées par une condition
# je peux utiliser la commande numpy where
# exemple :   
    
my_prn_for_extraction = 'G10'
my_date_for_extraction_start =  my_date_to_process + dt.timedelta(hours=3, minutes=27, seconds=0)
my_date_for_extraction_end =  my_date_to_process + dt.timedelta(hours=3, minutes=32, seconds=30)

condition1 = df_index.index.get_level_values('epoch') >= pd.Timestamp(my_date_for_extraction_start )
condition2 = df_index.index.get_level_values('epoch') <= pd.Timestamp(my_date_for_extraction_end)
condition3 = df_index.index.get_level_values('prn') == my_prn_for_extraction

np.where(condition1 & condition2 & condition3)

# une méthode plus rapide est directement d'ajouter une colonne de numéros de ligne au DataFrame original
# Ajouter une colonne de numéros de ligne au DataFrame original
df_index['ind_ligne'] = range(len(df_index))

# on a accès aux données et à leurs indices ...
my_prn_for_extraction = 'G10'
my_date_for_extraction_start =  my_date_to_process + dt.timedelta(hours=3, minutes=27, seconds=0)
my_date_for_extraction_end =  my_date_to_process + dt.timedelta(hours=3, minutes=42, seconds=30)

start = pd.Timestamp(my_date_for_extraction_start)
end   = pd.Timestamp(my_date_for_extraction_end)
extract2c = df_index.loc[(slice(start, end), my_prn_for_extraction), ('L1','ind_ligne')]

#%%
# A partir d'ici, vous êtes armés pour charger des fichiers RINEX d'observation
# et accéder facilement aux données.
# Obtenir une liste des PRNs uniques

my_obs_to_extract = 'L1'

prns = df_index.index.get_level_values('prn').unique()

# Créer une figure et un axe pour le plot
fig, ax = plt.subplots(figsize=(10, 6))

# Boucler sur chaque PRN et tracer sa série temporelle
for prn in prns:
    # Sélectionner les données pour le PRN actuel
    data = df_index.xs(prn, level='prn')
    # Tracer les données
    ax.plot(data.index.get_level_values('epoch'), data[my_obs_to_extract], label=prn)

# Configurer le graphique
ax.set_title('Séries temporelles '+my_obs_to_extract+' par PRN de satellite')
ax.set_xlabel('Temps')
ax.set_ylabel('Valeur')
ax.legend(title='PRN')

# Afficher le graphique
plt.show()


#%%

my_obs_to_extract = "L1"
gap = pd.Timedelta(minutes=30)

prns = df_index.index.get_level_values("prn").unique()

fig, ax = plt.subplots(figsize=(10, 6))

for prn in prns:
    # 1) données du PRN
    data = df_index.xs(prn, level="prn").copy()

    # 2) récupérer l'index temps (epoch) et trier
    t = data.index.get_level_values("epoch")
    data = data.set_index(t)
    data = data.sort_index()

    # 3) calcul des gaps et création d'un id de segment
    dt = data.index.to_series().diff()
    seg_id = (dt > gap).cumsum()

    # 4) couleur unique pour ce PRN (on la fixe via le 1er plot)
    color = None
    for _, seg in data.groupby(seg_id):
        line, = ax.plot(seg.index, seg[my_obs_to_extract], color=color, label=prn if color is None else None)
        if color is None:
            color = line.get_color()  # récupérer la couleur auto attribuée et la réutiliser

ax.set_title(f"Séries temporelles {my_obs_to_extract} par PRN de satellite")
ax.set_xlabel("Temps")
ax.set_ylabel("Valeur")
ax.legend(title="PRN")
plt.show()


#%%
# Il n'est pas nécessaire de garder dans le dataframe des colonnes inutilisées:
# Supprimer les colonnes où toutes les valeurs sont NaN
df_index = df_index.dropna(axis=1, how='all')

# On remarque que certaines mesures n'ont pas été réalisées sur L1 ou L2, ce qui
# pourrait être problématique lors de la formation des CL :

# Identifier les lignes avec NaN dans 'L1' ou 'L2'
rows_with_nan = df_index[['L1', 'L2']].isna().any(axis=1)

# Créer un DataFrame avec les lignes à supprimer
# potentiellement utile de savoir quand a eu lieu un éventuel problème
df_removed = df_index[rows_with_nan]


# Supprimer les lignes contenant NaN dans 'L1' ou 'L2' du DataFrame original
df_index = df_index.dropna(subset=['L1', 'L2'])

#%%
# Reste à s'occuper des satellites et obtenir leurs positions et corrections aux 
# dates d'intérêt

