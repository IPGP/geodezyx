#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 15:54:53 2023

@author: psakicki
"""

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names

import netCDF4
import numpy as np


## https://www.earthinversion.com/utilities/Writing-NetCDF4-Data-using-Python/

ncfile = netCDF4.Dataset('../../data/new.nc',mode='w',format='NETCDF4_CLASSIC')
print(ncfile)

lat_dim = ncfile.createDimension('lat', 73) # latitude axis
lon_dim = ncfile.createDimension('lon', 144) # longitude axis
time_dim = ncfile.createDimension('time', None) # unlimited axis (can be appended to).
for dim in ncfile.dimensions.items():
  print(dim)


ncfile.title='My model data'
print(ncfile.title)


ncfile.subtitle="My model data subtitle"
ncfile.anything="write anything"
print(ncfile.subtitle)
print(ncfile)
print(ncfile.anything)

lat = ncfile.createVariable('lat', np.float32, ('lat',))
lat.units = 'degrees_north'
lat.long_name = 'latitude'
lon = ncfile.createVariable('lon', np.float32, ('lon',))
lon.units = 'degrees_east'
lon.long_name = 'longitude'
time = ncfile.createVariable('time', np.float64, ('time',))
time.units = 'hours since 1800-01-01'
time.long_name = 'time'
temp = ncfile.createVariable('temp',np.float64,('time','lat','lon')) # note: unlimited dimension is leftmost
temp.units = 'K' # degrees Kelvin
temp.standard_name = 'air_temperature' # this is a CF standard name
print(temp)

print("-- Some pre-defined attributes for variable temp:")
print("temp.dimensions:", temp.dimensions)
print("temp.shape:", temp.shape)
print("temp.dtype:", temp.dtype)
print("temp.ndim:", temp.ndim)


nlats = len(lat_dim); nlons = len(lon_dim); ntimes = 3
lat[:] = -90. + (180./nlats)*np.arange(nlats) # south pole to north pole
lon[:] = (180./nlats)*np.arange(nlons) # Greenwich meridian eastward
data_arr = np.random.uniform(low=280,high=330,size=(ntimes,nlats,nlons))
temp[:,:,:] = data_arr # Appends data along unlimited dimension
print("-- Wrote data, temp.shape is now ", temp.shape)
print("-- Min/Max values:", temp[:,:,:].min(), temp[:,:,:].max())


data_slice = np.random.uniform(low=280,high=330,size=(nlats,nlons))
temp[3,:,:] = data_slice
print("-- Wrote more data, temp.shape is now ", temp.shape)

print(time)
times_arr = time[:]
print(type(times_arr),times_arr)

import datetime as dt
from netCDF4 import date2num,num2date
dates = [dt.datetime(2014,10,1,0),dt.datetime(2014,10,2,0),dt.datetime(2014,10,3,0),dt.datetime(2014,10,4,0)]
print(dates)

times = date2num(dates, time.units)
print(times, time.units) # numeric values

# first print the Dataset object to see what we've got
print(ncfile)
# close the Dataset.
ncfile.close(); print('Dataset is closed!')


