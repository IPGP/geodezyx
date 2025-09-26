#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 16 22:48:12 2023

@author: psakicki
"""

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names

import netCDF4
import numpy as np






P="/home/psakicki/GFZ_WORK/IPGP_WORK/REVOSIMA/0000_Pressure_Mayotte/010_from_Treden/RawData/transfer_3167556_files_8a99b7f8/204657_20210409_1130_data.txt"
DF = pd.read_table(P,sep=",")
#DF = DF.infer_objects()

#%%

p = "/home/psakicki/aaa_FOURBI/testmk1.nc"
ncfile = netCDF4.Dataset(p,mode='w',format='NETCDF4_CLASSIC')
print(ncfile)

size = len(DF)

time_dim = ncfile.createDimension('time', None) # unlimited axis (can be appended to).
instru_dim = ncfile.createDimension('instrument', size) 
temp_dim = ncfile.createDimension('temperature', size) 
pres_dim = ncfile.createDimension('pressure', size) 

ncfile.title='A0A OBP'
print(ncfile.title)

ncfile.subtitle="My model data subtitle"
ncfile.anything="write anything"
print(ncfile.subtitle)
print(ncfile)
print(ncfile.anything)

time = ncfile.createVariable('time', np.float64, ('time',))
time.units = 'hours since 1800-01-01'
time.long_name = 'time'

pres_s1 = ncfile.createVariable('pres_s1_raw',np.float64,('pressure',)) # note: unlimited dimension is leftmost
pres_s1.units = 'Pa'
pres_s1.standard_name = 'pressure raw sensor 1' # this is a CF standard name

pres_s2 = ncfile.createVariable('pres_s2_raw',np.float64,('pressure',)) # note: unlimited dimension is leftmost
pres_s2.units = 'Pa'
pres_s2.standard_name = 'pressure raw sensor 2' # this is a CF standard name


temp_s1 = ncfile.createVariable('temp_s1',np.float64,('temperature',)) # note: unlimited dimension is leftmost
temp_s1.units = 'C' # degrees Celcius
temp_s1.standard_name = 'temperature sensor 1' # this is a CF standard name


temp_s2 = ncfile.createVariable('temp_s2',np.float64,('temperature',)) # note: unlimited dimension is leftmost
temp_s2.units = 'C' # degrees Celcius
temp_s2.standard_name = 'temperature sensor 2' # this is a CF standard name


pres_baro = ncfile.createVariable('pres_baro',np.float64,('pressure',)) # note: unlimited dimension is leftmost
pres_baro.units = 'Pa'
pres_baro.standard_name = 'pressure sensor 2' # this is a CF standard name

temp_baro = ncfile.createVariable('temp_baro',np.float64,('temperature',)) # note: unlimited dimension is leftmost
temp_baro.units = 'C' # degrees Celcius
temp_baro.standard_name = 'temperature sensor 2' # this is a CF standard name

temp_ext = ncfile.createVariable('temp_extn',np.float64,('temperature',)) # note: unlimited dimension is leftmost
temp_ext.units = 'C' # degrees Celcius
temp_ext.standard_name = 'temperature extern' # this is a CF standard name


drift_pres_s1 = ncfile.createVariable('drift_pres_s1',np.float64,('pressure',)) # note: unlimited dimension is leftmost
drift_pres_s1.units = 'Pa'
drift_pres_s1.standard_name = 'drift pressure sensor 1' # this is a CF standard name

drift_pres_s2 = ncfile.createVariable('drift_pres_s2',np.float64,('pressure',)) # note: unlimited dimension is leftmost
drift_pres_s2.units = 'Pa'
drift_pres_s2.standard_name = 'pressure sensor 2' # this is a CF standard name



pres_s1_corr = ncfile.createVariable('pres_s1_corr',np.float64,('pressure',)) # note: unlimited dimension is leftmost
pres_s1_corr.units = 'Pa'
pres_s1_corr.standard_name = 'pressure raw sensor 1' # this is a CF standard name

pres_s2_corr = ncfile.createVariable('pres_s2_corr',np.float64,('pressure',)) # note: unlimited dimension is leftmost
pres_s2_corr.units = 'Pa'
pres_s2_corr.standard_name = 'pressure raw sensor 2' # this is a CF standard name




pres_s1[:] = DF['BPR pressure'].values
pres_s2[:] = DF['BPR pressure.1'].values

temp_s1[:] = DF['BPR temperature'].values
temp_s2[:] = DF['BPR temperature.1'].values

pres_baro[:] = DF['Barometer pressure'].values
temp_baro[:] = DF['Barometer temperature'].values

temp_ext[:] = DF['Temperature'].values

drift_pres_s1[:] = DF['Pressure drift'].values
drift_pres_s2[:] = DF['Pressure drift.1'].values

pres_s1_corr[:] = DF['BPR corrected pressure'].values
pres_s2_corr[:] = DF['BPR corrected pressure.1'].values



#ncfile.close()
