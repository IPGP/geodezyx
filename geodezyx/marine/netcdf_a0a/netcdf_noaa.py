#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 23:13:15 2023

@author: psakicki
"""

import netCDF4
import numpy as np
import matplotlib.pyplot as plt

p = "/home/psakicki/Downloads/104577.nc"
ncfile_ifremer = netCDF4.Dataset(p,mode='r',format='NETCDF4_CLASSIC')
print(ncfile_ifremer)
depth_kw = "depth" # NOAA
depth_kw = "PRES"
time = np.array(ncfile_ifremer["TIME"][:])
pres = np.array(ncfile_ifremer["PRES"][:])
pres = pres.squeeze()
time = time.squeeze()
plt.plot(time,pres)
ncfile_ifremer["LATITUDE"][:]



p = "/media/psakicki/F6FF-BB45/DOWNLOAD_SD/OLD_2307_2411/wc32_19910624to19920604.nc"
ncfile_noaa = netCDF4.Dataset(p,mode='r',format='NETCDF4_CLASSIC')
print(ncfile_noaa)

ncfile_noaa["timeSeries"]
ncfile_noaa["pressure_sensor_info"]


p_homemade = "/home/psakicki/aaa_FOURBI/out.nc"
ncfile_homemade = netCDF4.Dataset(p_homemade,mode='r',format='NETCDF4_CLASSIC')
print(ncfile_homemade)
