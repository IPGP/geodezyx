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
import geodezyx.files_rw

import pandas as pd
import numpy as np

# pour visualiser les données
import matplotlib.pyplot as plt

#%%
#
dir_data = "/home/snahmani/Bureau/FRS/data/data/data-2019/"
dir_data = "./data/data-2019/"

fichier_base= dir_data + '/mlvl1760.18d.Z'
fichier_mobile= dir_data + '/smne1760.18d.Z'

#%% Préambule
# chargement des pandas DataFrame et utilisation
# Lecture des observations RINEX en deux formats de DataFrame différents

df_flat = geodezyx.files_rw.read_rinex2_obs(fichier_base)
df_index = df_flat.set_index(['epoch', 'prn'])
# df_index['ind_ligne'] = range(len(df_index)) # à décommenter quand c'est compris

#%%
# Selection "FLAT" - Utilisation d'un DataFrame sans index personnalisé
# Avantage : Accès direct aux numéros de lignes pour faciliter la construction 
# de la matrice modèle

# Filtrage des observations pour le satellite 'G10' à un moment précis
bool_prn = df_flat['prn'] == 'G10'
bool_epoch = df_flat['epoch'] == pd.Timestamp('2018-06-25 13:20:00')
extract1a = df_flat.loc[bool_prn & bool_epoch, 'L1']
print(extract1a)

# Extraction des observations pour 'G10' sur une période spécifique
bool_prn = df_flat['prn'] == 'G10'
bool_epoch = (df_flat['epoch'] >= pd.Timestamp('2018-06-25 13:20:00')) & (df_flat['epoch'] <= pd.Timestamp('2018-06-25 13:30:00'))
extract1b = df_flat.loc[bool_prn & bool_epoch, 'L1']
print(extract1b)

# Extraction des observations pour 'G10' sur une période définie par une date de début et un delta de temps
bool_prn = df_flat['prn'] == 'G10'
start = pd.Timestamp('2018-06-25 13:20:00')
delta_T = pd.Timedelta(days=0, hours=1, minutes=15)
end = start + delta_T
bool_epoch = (df_flat['epoch'] >= start) & (df_flat['epoch'] <= end)
extract1c = df_flat.loc[bool_prn & bool_epoch, 'L1']
print(extract1c)

# Remarquez que l'utilisation des booléens nous permet d'accéder directement aux numéros
# des lignes concernées. Si on veut les numéros de ligne, il suffit de faire : 
serie_bool = bool_prn & bool_epoch ; 
serie_chiffre = serie_bool.astype(int) ; 

#%%
#### selection (multi)-index
# Je fais des filtrages avec .loc 
# Extraction pour le satellite 'G10' à un moment précis avec df_index
extract2a = df_index.loc[(pd.Timestamp('2018-06-25 13:20:00'), 'G10'), 'L1']
print(extract2a)

# Extraction pour 'G10' sur une période spécifique avec df_index
start_period = pd.Timestamp('2018-06-25 13:20:00')
end_period = pd.Timestamp('2018-06-25 13:30:00')
extract2b = df_index.loc[(slice(start_period, end_period), 'G10'), 'L1']
print(extract2b)

# Extraction pour 'G10' sur une période définie par une date de début et un delta de temps avec df_index
start = pd.Timestamp('2018-06-25 13:20:00')
delta_T = pd.Timedelta(days=0, hours=1, minutes=15)
end = start + delta_T
extract2c = df_index.loc[(slice(start, end), 'G10'), 'L1']
print(extract2c)

#%%
# bloc à commenter une fois que c'est compris.
# Dans ce cas, si je veux récuperer le numéro des lignes concernées par une condition
# je peux utiliser la commande numpy where
# exemple :    
condition1 = df_index.index.get_level_values('epoch') >= pd.Timestamp('2018-06-25 13:20:00')
condition2 = df_index.index.get_level_values('epoch') <= pd.Timestamp('2018-06-25 13:30:00')
condition3 = df_index.index.get_level_values('prn') == 'G10'

np.where(condition1 & condition2 & condition3)

# une méthode plus rapide est directement d'ajouter une colonne de numéros de ligne au DataFrame original
# Ajouter une colonne de numéros de ligne au DataFrame original
df_index['ind_ligne'] = range(len(df_index))

# on a accès aux données et à leurs indices ...
extract2c = df_index.loc[(slice(start, end), 'G10'), ('L1','ind_ligne')]

#%%
# A partir d'ici, vous êtes armés pour charger des fichiers RINEX d'observation
# et accéder facilement aux données.
# Obtenir une liste des PRNs uniques
prns = df_index.index.get_level_values('prn').unique()

# Créer une figure et un axe pour le plot
fig, ax = plt.subplots(figsize=(10, 6))

# Boucler sur chaque PRN et tracer sa série temporelle
for prn in prns:
    # Sélectionner les données pour le PRN actuel
    data = df_index.xs(prn, level='prn')
    # Tracer les données
    ax.plot(data.index.get_level_values('epoch'), data['L1'], label=prn)

# Configurer le graphique
ax.set_title('Séries temporelles L1 par PRN de satellite')
ax.set_xlabel('Temps')
ax.set_ylabel('Valeur')
ax.legend(title='PRN')

# Afficher le graphique
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

