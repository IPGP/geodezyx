#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 10:12:46 2024
Engineering Program ING3 - PPMD - Phase Measurement Processing
@author: Samuel Nahmani (1,2)
https://www.ipgp.fr/annuaire/nahmani/)
contact : nahmani@ipgp.fr ou samuel.nahmani@ign.fr
(1) Université Paris Cité, Institut de physique du globe de Paris, CNRS, IGN, F-75005 Paris, France.
(2) Univ Gustave Eiffel, ENSG, IGN, F-77455 Marne-la-Vallée, France. 

Version: 1.0
Dependencies: pandas, numpy, geodezyx, datetime

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

# To visualize the data
import matplotlib.pyplot as plt

from pathlib import Path
import os

#%%
# Create the gnss_edu_data folder which will contain the data and TP results
my_directory = os.environ["HOME"] + "/gnss_edu_data/"

# Path with expansion of ~ to home
folder = Path(my_directory).expanduser()

# Create the folder if it does not exist
folder.mkdir(parents=True, exist_ok=True)


#%%
# Automatic download of RINEX data from SMNE and MLVL stations approximately 10 km apart
# in the Paris region, France from the IGN (France) server
# data for 1 day (2019-176) at 30s

# Create a datetime to manage the processing day without having to manage doy, jjul, etc!
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


#%% Preamble
# Loading pandas DataFrames and their usage
# Reading RINEX observations in two different DataFrame formats

df_flat = files_rw.read_rinex_obs(fichier_base)

df_index = files_rw.read_rinex_obs(fichier_base, set_index=['epoch', 'prn'])

# What differences do you observe between the two dataframes df_flat and df_index?

#df_index['ind_ligne'] = range(len(df_index)) # uncomment when understood
# What happened to the df_index dataframe? Hint: open it and examine the columns...

#%%
# "FLAT" Selection - Using a DataFrame without custom index
# Advantage: Direct access to line numbers to facilitate the construction
# of the model matrix

# Filter observations for satellite 'G05' at a specific time
my_prn_for_extraction = 'G05'
my_date_for_extraction =  my_date_to_process + dt.timedelta(hours=2, minutes=40, seconds=0)
my_obs_to_extract = 'L1'

print(my_date_for_extraction)

bool_prn = df_flat['prn'] == my_prn_for_extraction
bool_epoch = df_flat['epoch'] == pd.Timestamp(my_date_for_extraction)
extract1a = df_flat.loc[bool_prn & bool_epoch, my_obs_to_extract]
print(extract1a)


# Extract observations for 'G05' over a specific period
my_prn_for_extraction = 'G05'
my_date_for_extraction_start =  my_date_to_process + dt.timedelta(hours=2, minutes=40, seconds=0)
my_date_for_extraction_end =  my_date_to_process + dt.timedelta(hours=2, minutes=45, seconds=30)
my_obs_to_extract = 'L1'

bool_prn = df_flat['prn'] == my_prn_for_extraction
bool_epoch = (df_flat['epoch'] >= pd.Timestamp(my_date_for_extraction_start)) & (df_flat['epoch'] <= pd.Timestamp(my_date_for_extraction_end))
extract1b = df_flat.loc[bool_prn & bool_epoch, my_obs_to_extract]
print(extract1b)


# Extract observations for 'G05' over a period defined by a start date and a time delta
my_prn_for_extraction = 'G05'
my_date_for_extraction =  my_date_to_process + dt.timedelta(hours=2, minutes=40, seconds=0)
my_obs_to_extract = 'L1'

bool_prn = df_flat['prn'] ==  my_prn_for_extraction
start = pd.Timestamp(my_date_for_extraction)
delta_T = pd.Timedelta(days=0, hours=1, minutes=15) # We can manage the delta T with pandas
end = start + delta_T
bool_epoch = (df_flat['epoch'] >= start) & (df_flat['epoch'] <= end)
extract1c = df_flat.loc[bool_prn & bool_epoch, my_obs_to_extract]
print(extract1c)

# Note that using booleans allows us to directly access the line numbers of interest.
# If we want the line numbers, we just need to do:
serie_bool = bool_prn & bool_epoch ;
serie_chiffre = serie_bool.astype(int) ; 



#%%
#### (Multi)-index selection
# I do filtering with .loc
# Extract for satellite 'G10' at a specific moment with df_index

my_prn_for_extraction = 'G10'
my_date_for_extraction =  my_date_to_process + dt.timedelta(hours=3, minutes=27, seconds=30)
my_obs_to_extract = 'L1'

extract2a = df_index.loc[(pd.Timestamp(my_date_for_extraction), my_prn_for_extraction), my_obs_to_extract]
print(extract2a)

# Extract for 'G10' over a specific period with df_index
my_prn_for_extraction = 'G10'
my_date_for_extraction_start =  my_date_to_process + dt.timedelta(hours=3, minutes=27, seconds=0)
my_date_for_extraction_end =  my_date_to_process + dt.timedelta(hours=3, minutes=32, seconds=30)
my_obs_to_extract = 'L1'

start_period = pd.Timestamp(my_date_for_extraction_start)
end_period = pd.Timestamp(my_date_for_extraction_end )
extract2b = df_index.loc[(slice(start_period, end_period), my_prn_for_extraction), my_obs_to_extract]
print(extract2b)

# Extract for 'G10' over a period defined by a start date and a time delta with df_index
my_prn_for_extraction = 'G10'
my_date_for_extraction =  my_date_to_process + dt.timedelta(hours=3, minutes=27, seconds=30)
my_obs_to_extract = 'L1'

start = pd.Timestamp(my_date_for_extraction)
delta_T = pd.Timedelta(days=0, hours=1, minutes=15)
end = start + delta_T
extract2c = df_index.loc[(slice(start, end), my_prn_for_extraction), my_obs_to_extract]
print(extract2c)

#%%
# Section to comment out once this is understood.
# In this case, if I want to retrieve the line numbers affected by a condition
# I can use the numpy where command
# example:

my_prn_for_extraction = 'G10'
my_date_for_extraction_start =  my_date_to_process + dt.timedelta(hours=3, minutes=27, seconds=0)
my_date_for_extraction_end =  my_date_to_process + dt.timedelta(hours=3, minutes=32, seconds=30)

condition1 = df_index.index.get_level_values('epoch') >= pd.Timestamp(my_date_for_extraction_start )
condition2 = df_index.index.get_level_values('epoch') <= pd.Timestamp(my_date_for_extraction_end)
condition3 = df_index.index.get_level_values('prn') == my_prn_for_extraction

np.where(condition1 & condition2 & condition3)

# A faster method is to directly add a column of line numbers to the original DataFrame
# Add a column of line numbers to the original DataFrame
df_index['ind_ligne'] = range(len(df_index))

# We have access to data and their indices ...
my_prn_for_extraction = 'G10'
my_date_for_extraction_start =  my_date_to_process + dt.timedelta(hours=3, minutes=27, seconds=0)
my_date_for_extraction_end =  my_date_to_process + dt.timedelta(hours=3, minutes=42, seconds=30)

start = pd.Timestamp(my_date_for_extraction_start)
end   = pd.Timestamp(my_date_for_extraction_end)
extract2c = df_index.loc[(slice(start, end), my_prn_for_extraction), ('L1','ind_ligne')]

#%%
# From here on, you are equipped to load RINEX observation files
# and easily access the data.
# Get a list of unique PRNs

my_obs_to_extract = 'L1'

prns = df_index.index.get_level_values('prn').unique()

# Create a figure and an axis for the plot
fig, ax = plt.subplots(figsize=(10, 6))

# Loop over each PRN and plot its time series
for prn in prns:
    # Select the data for the current PRN
    data = df_index.xs(prn, level='prn')
    # Plot the data
    ax.plot(data.index.get_level_values('epoch'), data[my_obs_to_extract], label=prn)

# Configure the graph
ax.set_title('Time series '+my_obs_to_extract+' by satellite PRN')
ax.set_xlabel('Time')
ax.set_ylabel('Value')
ax.legend(title='PRN')

# Display the graph
plt.show()


#%%

my_obs_to_extract = "L1"
gap = pd.Timedelta(minutes=30)

prns = df_index.index.get_level_values("prn").unique()

fig, ax = plt.subplots(figsize=(10, 6))

for prn in prns:
    # 1) data of the PRN
    data = df_index.xs(prn, level="prn").copy()

    # 2) get the time index (epoch) and sort
    t = data.index.get_level_values("epoch")
    data = data.set_index(t)
    data = data.sort_index()

    # 3) calculate gaps and create a segment id
    dt = data.index.to_series().diff()
    seg_id = (dt > gap).cumsum()

    # 4) unique color for this PRN (we set it via the 1st plot)
    color = None
    for _, seg in data.groupby(seg_id):
        line, = ax.plot(seg.index, seg[my_obs_to_extract], color=color, label=prn if color is None else None)
        if color is None:
            color = line.get_color()  # retrieve the auto-assigned color and reuse it

ax.set_title(f"Time series {my_obs_to_extract} by satellite PRN")
ax.set_xlabel("Time")
ax.set_ylabel("Value")
ax.legend(title="PRN")
plt.show()



#%%
# It is not necessary to keep unused columns in the dataframe:
# Remove columns where all values are NaN
df_index = df_index.dropna(axis=1, how='all')

# We notice that some measurements were not made on L1 or L2, which
# could be problematic when forming the combinations:

# Identify rows with NaN in 'L1' or 'L2'
rows_with_nan = df_index[['L1', 'L2']].isna().any(axis=1)

# Create a DataFrame with the rows to be removed
# potentially useful to know when a possible problem occurred
df_removed = df_index[rows_with_nan]


# Remove rows containing NaN in 'L1' or 'L2' from the original DataFrame
df_index = df_index.dropna(subset=['L1', 'L2'])

#%%
# Still need to handle the satellites and get their positions and corrections at
# the dates of interest

