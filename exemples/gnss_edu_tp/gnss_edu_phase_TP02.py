#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 16:44:03 2024

ING3 - PPMD - Phase measurement processing

@author: Samuel Nahmani (1,2)
https://www.ipgp.fr/annuaire/nahmani/)
contact : nahmani@ipgp.fr ou samuel.nahmani@ign.fr
(1) Université Paris Cité, Institut de physique du globe de Paris, CNRS, IGN, F-75005 Paris, France.
(2) Univ Gustave Eiffel, ENSG, IGN, F-77455 Marne-la-Vallée, France. 

Version: 1.0
Dependencies: pandas, numpy, geodezyx, datetime, gpsdatetime, gnsstoolbox

"""
# %%
# GeodeZYX Toolbox's
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
# GeodeZYX Toolbox's - [Sakic et al., 2019]
from geodezyx import files_rw     # Import the read/write module
from geodezyx import conv         # Import the conversion module
from geodezyx import operational  # Import the download rinex module
from geodezyx import gnss_edu     # Import the learning module
from geodezyx import reffram      # Import the reference frame/higher geodesy module

import datetime as dt
#
import gpsdatetime as gpst
import gnsstoolbox.orbits as orb


import pandas as pd
import numpy as np

# to visualize data
import matplotlib.pyplot as plt

from pathlib import Path
import os


# %%
# Create the gnss_edu_data folder that will contain the data and results of the practical work
# see gnss_edu_phase_TP01.py

my_directory = os.environ["HOME"] + "/gnss_edu_data/"

# Path with ~ expansion to home directory
folder = Path(my_directory).expanduser()

# Create the folder if it does not exist
folder.mkdir(parents=True, exist_ok=True)

# %%
# Automatic download of RINEX data for stations SMNE and MLVL, about ten
# kilometers apart in the Paris region, France, from the IGN server (France)
# data for 1 day (2019-176) at 30s
# see gnss_edu_phase_TP01.py

# Create a datetime to manage the processing day without having to deal with doy, julian days, etc.
my_date_to_process = dt.datetime(2019,6,25)

dwl_output_station = operational.download_gnss_rinex(statdico={"rgp" : ["SMNE","MLVL"]},
                                output_dir=my_directory,
                                startdate= my_date_to_process ,
                                enddate= my_date_to_process ,
                                parallel_download = 1)

dwl_output_navigation = operational.download_gnss_rinex(statdico={"nav" : ["brdc"]},
                                output_dir=my_directory,
                                startdate= my_date_to_process ,
                                enddate= my_date_to_process ,
                                parallel_download = 1)

dwl_output_satellite = operational.download_gnss_products(archive_dir= my_directory,
                                   startdate= my_date_to_process,
                                   enddate= my_date_to_process,
                                   AC_names=("IGS",),
                                   repro=0,
                                   archive_center="ign",
                                   parallel_download = 1,
                                   )
# %%
# Load RINEX observation files
obs_rnx_file = dwl_output_station[0][0]
nav_file = dwl_output_navigation[0][0]

# Load the RINEX navigation file via GeodeZYX
iono_cor_dic, time_sys_corr_dic = files_rw.read_rinex_nav_v3_header(nav_file)

# Load RINEX observation data into a pandas dataframe via GeodeZYX
df_rnx_orig, l_rnx_head = files_rw.read_rinex2_obs(obs_rnx_file, return_header=True)
df_rnx = df_rnx_orig.copy()


# Approximate position read from the RINEX file header
columns = [l for l in l_rnx_head if r"APPROX POSITION XYZ" in l]
if columns:
    # Split the line by removing the text "APPROX POSITION XYZ"
    values = columns[0].split()[:3]  # Retrieve only the first 3 values
    P_rnx_header = np.array(values, dtype=float)  # Convert to numpy array
    print("XYZ Coordinates:", P_rnx_header)
else:
    P_rnx_header = np.array([0,0,0])
    print("No match found.")
del columns, values


# You should check the documentation and you will notice that the units of L1 and L2 are in cycles.
# We multiply the L1 and L2 values by the corresponding wavelengths to obtain
# distances (ambiguous but in meters).

df_rnx['L1'] = df_rnx['L1']*conv.L1_WAVELENGTH
df_rnx['L2'] = df_rnx['L2']*conv.L2_WAVELENGTH
df_rnx['L5'] = df_rnx['L5']*conv.L5_WAVELENGTH


# cleaning cf end of TP01
df_rnx = df_rnx.dropna(axis=1, how='all')
rows_with_nan = df_rnx[['C1','L1', 'L2']].isna().any(axis=1)
df_removed = df_rnx[rows_with_nan]
df_rnx = df_rnx.dropna(subset=['C1','L1', 'L2'])
del rows_with_nan

# Keep only measurements from GPS satellites
# Filter to keep only rows where the 'sys' column contains 'G'
df_filtered = df_rnx[df_rnx['sys'].str.contains('G')]

df_rnx = df_filtered
del df_filtered

df_rnx = df_rnx.dropna(axis=1, how='all')

# add a row index
df_rnx['row_ind'] = range(len(df_rnx))


# %%
sp3_file = dwl_output_satellite[0]
df_sp3 = files_rw.read_sp3(sp3_file)

# %%

# Load orbit files
brdc_file = './data/data-2019/mlvl176z.18n'

X_sat = []
Y_sat = []
Z_sat = []
dte_sat = []
dRelat = []

#### NEW GEODEZYX STYLE: use GeodeZYX interpolation functions to compute satellite positions at each emission time
for prn in df_rnx['prn'].unique():
    df_rnx_prn = df_rnx[df_rnx['prn'] == prn]
    if not prn in df_sp3['prn'].unique():
        print(f"Satellite {prn} not found in the SP3 file, cannot compute its position at each emission time")
        continue
    df_sp3_prn = df_sp3[df_sp3['prn'] == prn]
    fly_time = df_rnx_prn['C1'] / conv.SPEED_OF_LIGHT
    fly_time = pd.to_timedelta(fly_time, unit="s")
    t_rec = df_rnx_prn['epoch']
    t_emi_gross = t_rec - fly_time
    orb_df_gross = reffram.orb_df_lagrange_interpolate(df_sp3_prn, t_emi_gross.values)
    # positions in SP3 files are in km, convert to meters to be consistent with pseudorange units
    orb_df_gross[["x","y","z"]] = orb_df_gross[["x","y","z"]] * 10**3 # km -> m

    # relativistic effect computation
    delta_t = pd.to_timedelta(1e-3, unit="s") # time offset +/- for computing the derivative
    orb_df_rel_fwd = reffram.orb_df_lagrange_interpolate(df_sp3_prn, t_emi_gross.values + delta_t)
    orb_df_rel_fwd[["x","y","z"]] = orb_df_rel_fwd[["x","y","z"]] * 10**3 # km -> m
    orb_df_rel_bak = reffram.orb_df_lagrange_interpolate(df_sp3_prn, t_emi_gross.values - delta_t)
    orb_df_rel_bak[["x","y","z"]] = orb_df_rel_bak[["x","y","z"]] * 10**3 # km -> m

    xyz_diff = (orb_df_rel_fwd[["x","y","z"]] - orb_df_rel_bak[["x","y","z"]]) / (2 * delta_t.total_seconds())
    xyz = orb_df_gross[["x","y","z"]]

    dRelat_v = -2.0 * (xyz_diff * xyz).sum(axis=1) / (conv.SPEED_OF_LIGHT **2)

    t_emi_ok = t_emi_gross -  pd.to_timedelta(orb_df_gross["clk"].values, unit="us") # - pd.to_timedelta(dRelat_v, unit="s")
    orb_df_ok = reffram.orb_df_lagrange_interpolate(df_sp3_prn, t_emi_ok.values)
    orb_df_ok[["x","y","z"]] = orb_df_ok[["x","y","z"]] * 10**3 # km -> m

    X_sat.extend(orb_df_ok["x"].values)
    Y_sat.extend(orb_df_ok["y"].values)
    Z_sat.extend(orb_df_ok["z"].values)
    dte_sat.extend(orb_df_ok["clk"].values * 1e-6) # microseconds -> seconds
    dRelat.extend(dRelat_v.values)

    df_rnx.loc[df_rnx_prn.index, 'X_sat'] = orb_df_ok["x"].values
    df_rnx.loc[df_rnx_prn.index, 'Y_sat'] = orb_df_ok["y"].values
    df_rnx.loc[df_rnx_prn.index, 'Z_sat'] = orb_df_ok["z"].values
    df_rnx.loc[df_rnx_prn.index, 'dte_sat'] = orb_df_ok["clk"].values * 1e-6 # microseconds -> seconds
    df_rnx.loc[df_rnx_prn.index, 'dRelat'] = dRelat_v.values

df_rinex_flat = df_rnx.copy()
df_rnx = df_rnx.set_index(['epoch', 'prn'], drop=True)


### OLD STYLE ###################################################################

df_rnx_old = df_rnx.copy()

mysp3 = orb.orbit()
mysp3.loadSp3(sp3_file)

mynav = orb.orbit()
mynav.loadRinexN(nav_file)

# Compute the position of each GNSS satellite at each emission time
t = gpst.gpsdatetime()
X_sat_old = []
Y_sat_old = []
Z_sat_old = []
dte_sat_old = []
dRelat_old = []

for (time_i,prn_i) in df_rnx_old.index:
    t.rinex_t(time_i.to_pydatetime().strftime('%y %m %d %H %M %S.%f'))
    t_emission_mjd  = t.mjd - df_rnx_old.loc[(time_i,prn_i), 'C1'] / conv.SPEED_OF_LIGHT / 86400.0
    (X_sat_v,Y_sat_v,Z_sat_v,dte_sat_v)	 = mysp3.calcSatCoord(prn_i[0], int(prn_i[1:]),t_emission_mjd)

    # relativistic effect computation
    delta_t = 1e-3 # time offset +/- for computing the derivative
    (Xs1,Ys1,Zs1,clocks1) = mysp3.calcSatCoord(prn_i[0], int(prn_i[1:]),t_emission_mjd - delta_t / 86400.0)
    (Xs2,Ys2,Zs2,clocks2) = mysp3.calcSatCoord(prn_i[0], int(prn_i[1:]),t_emission_mjd + delta_t / 86400.0)  

    VX      = (np.array([Xs2-Xs1, Ys2-Ys1, Zs2-Zs1]))/2.0/delta_t
    VX0     = np.array([X_sat_v,Y_sat_v,Z_sat_v])

    dRelat_v1  = -2.0 * VX0.T @ VX /(conv.SPEED_OF_LIGHT **2)
    dRelat_v2 = -2.0 * np.dot(VX0, VX) / (conv.SPEED_OF_LIGHT**2)

    # GNSS signal emission time in GPS time (mjd)
    t_emission_mjd = t_emission_mjd - dte_sat_v / 86400.0 - dRelat_v1 / 86400.0

    # Recompute satellite position at emission time (GPS time in mjd)
    (X_sat_v,Y_sat_v,Z_sat_v,dte_sat_v)	 = mysp3.calcSatCoord(prn_i[0], int(prn_i[1:]),t_emission_mjd)

    X_sat_old.append(X_sat_v)
    Y_sat_old.append(Y_sat_v)
    Z_sat_old.append(Z_sat_v)
    dte_sat_old.append(dte_sat_v)
    dRelat_old.append(dRelat_v1)

df_rnx_old["X_sat"] = X_sat_old
df_rnx_old["Y_sat"] = Y_sat_old
df_rnx_old["Z_sat"] = Z_sat_old
df_rnx_old["dte_sat"] = dte_sat_old
df_rnx_old["dRelat"] = dRelat_old

df_rnx_old[["X_sat","Y_sat","Z_sat","dte_sat","dRelat"]]
df_rnx_old[["X_sat","Y_sat","Z_sat","dte_sat","dRelat"]]

# del X_sat, Y_sat, Z_sat, dte_sat, time_i, prn_i, t_emission_mjd, t, X_sat_v, Y_sat_v, Z_sat_v, dte_sat_v, dRelat_v, dRelat

############## OLD STYLE END ###################################################################

# %%
# Classical code-based processing for GNSS receiver positioning
# trivial case

print('*****  Trivial computation *****')


# Initialization of approximate receiver coordinates
P_app = np.array([0, 0, 0])

# Iteration for receiver position refinement
dP_est=np.array([100, 100, 100]) # initialization: just to enter the while loop below
i=1
while np.linalg.norm(dP_est)>1:
    # Compute approximate satellite-receiver distances
    distances = np.sqrt((df_rnx['X_sat'].values - P_app[0])**2 +
                        (df_rnx['Y_sat'].values - P_app[1])**2 +
                        (df_rnx['Z_sat'].values - P_app[2])**2)

    # Observation vector corrected by the various models
    B = df_rnx['C1'].values - distances

    # Build the partial derivatives matrix
    df_dX = (P_app[0] - df_rnx['X_sat'].values) / distances
    df_dY = (P_app[1] - df_rnx['Y_sat'].values) / distances
    df_dZ = (P_app[2] - df_rnx['Z_sat'].values) / distances
    A = np.column_stack((df_dX, df_dY, df_dZ))

    # Least-squares resolution to estimate the displacement
    dP_est = np.linalg.inv(A.T@A)@A.T@B
    #dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
    
    # Update the estimated position
    P_est = P_app + dP_est
    
    # Print the estimated position at each iteration
    print(f"Iteration {i}: Estimated position - X: {P_est[0]}, Y: {P_est[1]}, Z: {P_est[2]}")
    P_app = P_est  # Update the approximate position for the next iteration
    i+=1

    # Compute the final distance between the estimated position and the initial RINEX header position
    dist_P_est_P_rnx_header = np.sqrt(np.sum((P_est - P_rnx_header)**2))
    print("\n")
    print("Distance between estimated position and RINEX header position:", dist_P_est_P_rnx_header)

    # Compute local ENU coordinates
    E, N, U = conv.xyz2enu(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
    print("East (E):", E)
    print("North (N):", N)
    print("Up (U):", U)
    print("\n")


fig = gnss_edu.plot_residual_analysis(A, B, dP_est, figure_title="Trivial calcul.", 
                                      save_path="./trivial_bis.pdf",
                           P_est=P_est, P_rnx_header=P_rnx_header)


# del E, N, U, P_est, dP_est, P_app, A, B, df_dX, df_dY, df_dZ, distances, i

# %%
# Classical code-based processing for GNSS receiver positioning
# satellite clock error corrections
# -> they are known -> direct correction of vector B

print('*****  Satellite clock error correction *****')

# Initialization of approximate receiver coordinates
P_app = np.array([0, 0, 0])
dP_est=100
i=1
# Iteration for receiver position refinement
while np.linalg.norm(dP_est)>1:
    # Compute approximate satellite-receiver distances
    distances = np.sqrt((df_rnx['X_sat'].values - P_app[0])**2 +
                        (df_rnx['Y_sat'].values - P_app[1])**2 +
                        (df_rnx['Z_sat'].values - P_app[2])**2)

    # Observation vector corrected by the various models
    B = df_rnx['C1'].values - distances + conv.SPEED_OF_LIGHT * (df_rnx['dte_sat'].values + df_rnx['dRelat'].values )

    # Build the partial derivatives matrix
    df_dX = (P_app[0] - df_rnx['X_sat'].values) / distances
    df_dY = (P_app[1] - df_rnx['Y_sat'].values) / distances
    df_dZ = (P_app[2] - df_rnx['Z_sat'].values) / distances
    A = np.column_stack((df_dX, df_dY, df_dZ))

    # Least-squares resolution to estimate the displacement
    dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
    
    # Update the estimated position
    P_est = P_app + dP_est
    
    # Print the estimated position at each iteration
    print(f"Iteration {i}: Estimated position - X: {P_est[0]}, Y: {P_est[1]}, Z: {P_est[2]}")
    P_app = P_est  # Update the approximate position for the next iteration
    i+=1

# Compute the final distance between the estimated position and the initial RINEX header position
dist_P_est_P_rnx_header = np.sqrt(np.sum((P_est - P_rnx_header)**2))
print("\n")
print("Distance between estimated position and RINEX header position:", dist_P_est_P_rnx_header)

# Compute local ENU coordinates
E, N, U = conv.xyz2enu(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
print("East (E):", E)
print("North (N):", N)
print("Up (U):", U)
print("\n")


gnss_edu.plot_residual_analysis(A, B, dP_est, figure_title="Satellite clocks correction", save_path="./sat_clocks_only.png",
                           P_est=P_est, P_rnx_header=P_rnx_header)


# del E, N, U, P_est, dP_est, P_app, A, B, df_dX, df_dY, df_dZ, distances, i

# %%
# Classical code-based processing for GNSS receiver positioning
# satellite clock error corrections
# -> they are known -> direct correction of vector B

print('*****  Satellite and receiver clock error correction *****')

# Initialization of approximate receiver coordinates
P_app = np.array([0, 0, 0])
dP_est=100
i=1
# Iteration for receiver position refinement
while np.linalg.norm(dP_est)>1:
    # Compute approximate satellite-receiver distances
    distances = np.sqrt((df_rnx['X_sat'].values - P_app[0])**2 +
                        (df_rnx['Y_sat'].values - P_app[1])**2 +
                        (df_rnx['Z_sat'].values - P_app[2])**2)

    # Observation vector corrected by the various models
    B = df_rnx['C1'].values - distances + conv.SPEED_OF_LIGHT * (df_rnx['dte_sat'].values + df_rnx['dRelat'].values )

    # Build the partial derivatives matrix
    df_dX = (P_app[0] - df_rnx['X_sat'].values) / distances
    df_dY = (P_app[1] - df_rnx['Y_sat'].values) / distances
    df_dZ = (P_app[2] - df_rnx['Z_sat'].values) / distances
    A = np.column_stack((df_dX, df_dY, df_dZ))

    # Least-squares resolution to estimate the displacement
    dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
    
    # Update the estimated position
    P_est = P_app + dP_est
    
    # Print the estimated position at each iteration
    print(f"Iteration {i}: Estimated position - X: {P_est[0]}, Y: {P_est[1]}, Z: {P_est[2]}")
    P_app = P_est  # Update the approximate position for the next iteration
    i+=1

# Compute the final distance between the estimated position and the initial RINEX header position
dist_P_est_P_rnx_header = np.sqrt(np.sum((P_est - P_rnx_header)**2))
print("\n")
print("Distance between estimated position and RINEX header position:", dist_P_est_P_rnx_header)

# Compute local ENU coordinates
E, N, U = conv.xyz2enu(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
print("East (E):", E)
print("North (N):", N)
print("Up (U):", U)
print("\n")


gnss_edu.plot_residual_analysis(A, B, dP_est, figure_title="Calcul corr Sat", save_path="./corr_sat_clock.png",
                           P_est=P_est, P_rnx_header=P_rnx_header)


del E, N, U, P_est, dP_est, P_app, A, B, df_dX, df_dY, df_dZ, distances, i


# %%
# Classical code-based processing for GNSS receiver positioning
# satellite clock error corrections
# -> they are known -> direct correction of vector B

print('*****  Satellite clock error correction *****')
print('*****  + Sagnac *****')

# introduce a visualization of the result with folium
def rotate_around_z(row):
    # C1 / c -> flight time
    # flight time * Earth rotation angular velocity -> alpha_rad
    alpha_rad = row['C1'] / conv.SPEED_OF_LIGHT * conv.EARTH_ROTATION_MEAN_ANGULAR_VELOCITY
    # Rotation matrix around Z axis
    Rz = np.array([[np.cos(alpha_rad), -np.sin(alpha_rad), 0],
                   [np.sin(alpha_rad), np.cos(alpha_rad), 0],
                   [0, 0, 1]])
    # Original position vector
    original_vector = np.array([row['X_sat'], row['Y_sat'], row['Z_sat']])
    # Compute the rotated position vector
    rotated_vector = Rz.dot(original_vector)
    return pd.Series(rotated_vector, index=['X_sat', 'Y_sat', 'Z_sat'])

# Apply rotation to each row and create a new dataframe with the results
df_Sagnac = df_rnx.apply(rotate_around_z, axis=1)


# Initialization of approximate receiver coordinates
P_app = np.array([0, 0, 0])
dP_est=100
i=1
# Iteration for receiver position refinement
while np.linalg.norm(dP_est)>1:
    # Compute approximate satellite-receiver distances
    distances = np.sqrt((df_Sagnac['X_sat'].values - P_app[0])**2 +
                        (df_Sagnac['Y_sat'].values - P_app[1])**2 +
                        (df_Sagnac['Z_sat'].values - P_app[2])**2)

    # Observation vector corrected by the various models
    B = df_rnx['C1'].values - distances + conv.SPEED_OF_LIGHT * (df_rnx['dte_sat'].values + df_rnx['dRelat'].values )

    # Build the partial derivatives matrix
    df_dX = (P_app[0] - df_Sagnac['X_sat'].values) / distances
    df_dY = (P_app[1] - df_Sagnac['Y_sat'].values) / distances
    df_dZ = (P_app[2] - df_Sagnac['Z_sat'].values) / distances
    A = np.column_stack((df_dX, df_dY, df_dZ))

    # Least-squares resolution to estimate the displacement
    dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
    
    # Update the estimated position
    P_est = P_app + dP_est
    
    # Print the estimated position at each iteration
    print(f"Iteration {i}: Estimated position - X: {P_est[0]}, Y: {P_est[1]}, Z: {P_est[2]}")
    P_app = P_est  # Update the approximate position for the next iteration
    i+=1

# Compute the final distance between the estimated position and the initial RINEX header position
dist_P_est_P_rnx_header = np.sqrt(np.sum((P_est - P_rnx_header)**2))
print("\n")
print("Distance between estimated position and RINEX header position:", dist_P_est_P_rnx_header)

# Compute local ENU coordinates
E, N, U = conv.xyz2enu(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
print("East (E):", E)
print("North (N):", N)
print("Up (U):", U)
print("\n")

del E, N, U, P_est, dP_est, P_app, A, B, df_dX, df_dY, df_dZ, distances, i


# %%
# Classical code-based processing for GNSS receiver positioning
# satellite clock error corrections
# -> they are known -> direct correction of vector B
# receiver clock error estimation
# -> they are unknown and must be estimated at each epoch
# ionospheric effect correction using a linear combination of code observations
# -> direct correction of vector B

# Sagnac effect correction

# tropospheric effect correction

print('*****  Satellite clock error correction *****')
print('*****  + Sagnac *****')
print('*****  + receiver clocks *****')
print('*****  + iono LC(C1,P2) *****')
print('*****  + TROPO     *****')
print('*****  + cutoff angle     *****')
print('*****  + ambiguity estimation     *****')

print('*****  AND NORMALLY... ARG :) we are in PLOUF  *****')


# %%
# Get unique epochs
epoch_uniques = df_rnx_new.index.get_level_values('epoch').unique()
nb_epochs = len(epoch_uniques)
# Skeleton block to concatenate to the design matrix A
block_dt_r = np.zeros((len(df_rnx_new), nb_epochs))

# Fill the block corresponding to the estimation of receiver clock errors
for i, epoch in enumerate(epoch_uniques):
    block_dt_r[df_rnx_new.loc[epoch, 'row_ind'],i]=1

# Fill the block corresponding to the estimation of wet tropospheric delays (estimated every 2 hours)
start = df_rnx_new.index[0][0]
delta_T = pd.Timedelta(days=0, hours=2, minutes=0)
delta_sec = pd.Timedelta(seconds=1)
end   = start + delta_T - delta_sec

import math

delta = df_rnx_new.index[-1][0] - df_rnx_new.index[0][0]
nb_par_zwd = math.ceil(delta / pd.Timedelta(hours=2))


# Skeleton block to concatenate to the design matrix A
block_zwd = np.zeros((len(df_rnx_new), nb_par_zwd))
ind_c=0

while start <= df_rnx_new.index[-1][0]:
    
    end   = start + delta_T - delta_sec
    extract2c = df_rnx_new.loc[(slice(start, end)), ('row_ind','mfw')]
    block_zwd[extract2c['row_ind'],ind_c]=extract2c['mfw']

    start = start + delta_T
    ind_c = ind_c + 1


# Initialization of approximate receiver coordinates
P_app = np.array([0, 0, 0])

dP_est=np.array([100, 100, 100])
i=1
# Iteration for receiver position refinement
while np.linalg.norm(dP_est[0:3])>1:

    # Compute approximate satellite-receiver distances
    distances = np.sqrt((df_Sagnac_new['X_sat'].values - P_app[0])**2 +
                        (df_Sagnac_new['Y_sat'].values - P_app[1])**2 +
                        (df_Sagnac_new['Z_sat'].values - P_app[2])**2)

    # Observation vector corrected by the various models
    B = df_rnx_new['L3'].values - distances + conv.SPEED_OF_LIGHT * (df_rnx_new['dte_sat'].values \
    + df_rnx_new['dRelat'].values ) \
    - df_rnx_new['ZHD'].values* df_rnx_new['mfh'].values

    # Build the partial derivatives matrix
    df_dX = (P_app[0] - df_Sagnac_new['X_sat'].values) / distances
    df_dY = (P_app[1] - df_Sagnac_new['Y_sat'].values) / distances
    df_dZ = (P_app[2] - df_Sagnac_new['Z_sat'].values) / distances
    A = np.column_stack((df_dX, df_dY, df_dZ, block_dt_r, block_zwd, block_one_amb))

    # Least-squares resolution to estimate the displacement
    dP_est, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
    
    # Update the estimated position
    # Note: nb_epochs receiver clock parameters have been estimated
    P_est = np.zeros(len(dP_est))
    
    
    P_est[0:3] = P_app[0:3]+dP_est[0:3]
    P_est[3:]  = dP_est[3:]
    

    
    # Print the estimated position at each iteration
    print(f"Iteration {i}: Estimated position - X: {P_est[0]}, Y: {P_est[1]}, Z: {P_est[2]}")
    P_app = P_est  # Update the approximate position for the next iteration
    i+=1


# Compute the final distance between the estimated position and the initial RINEX header position
dist_P_est_P_rnx_header = np.sqrt(np.sum((P_est[:3] - P_rnx_header)**2))
print("\n")
print("Distance between estimated position and RINEX header position:", dist_P_est_P_rnx_header)

# Compute local ENU coordinates
E, N, U = conv.xyz2enu(P_rnx_header[0], P_rnx_header[1], P_rnx_header[2],P_est[0], P_est[1], P_est[2])
print("East (E):", E)
print("North (N):", N)
print("Up (U):", U)
print("\n")

# del E, N, U, P_est, dP_est, P_app, A, block_dt_r, B, df_dX, df_dY, df_dZ, distances, i, epoch, epoch_uniques, nb_epochs

# introduce a visualization of the result with folium

# plot_residual_analysis(A, B, dP_est, figure_title="Calcul corr Sat Rec et Sagnac et iono et tropo", save_path="./corr_sat_clock_sagnac_iono_tropo.png",
#                           P_est=P_est[:3], P_rnx_header=P_rnx_header)


# %%
f1 = conv.L1_CARRIER_FREQUENCY
f2 = conv.L2_CARRIER_FREQUENCY

l1 = conv.L1_WAVELENGTH
l2 = conv.L2_WAVELENGTH

# Apply the linear combination formula
# l3 = (f1**2*l1- f2**2*l2)/(f1**2-f2**2);
# which simplifies to:
l3=conv.SPEED_OF_LIGHT/(f1+f2)


# note: L1 and L2 measurements are in cycles whereas we need
# a distance measurement (ambiguous but still a distance)

df_rnx_new['L3']  = (f1**2*df_rnx_new['L1'] - f2**2*df_rnx_new['L2'])/(f1**2-f2**2);
df_rnx_new['P3']  = (f1**2*df_rnx_new['C1'] - f2**2*df_rnx_new['P2'])/(f1**2-f2**2);


gnss_edu.plot_series(df=df_rnx_new, col1='L3', col2='P3' , coeff1=1.0, coeff2=1.0, seuil=3600, renderer="browser")

# %%
# compute other linear combinations

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


