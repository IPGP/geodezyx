#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# https://chatgpt.com/c/683ec5ef-d6d0-8012-8163-dc4e8a6297b3
"""
Created on Tue Jun  3 11:52:28 2025

@author: psakicki
"""
from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta, UTC

# Create NetCDF file
ncfile = Dataset('insitu_cmems_example.nc', mode='w', format='NETCDF4')

# Define global attributes (partial list per CMEMS/SeaDataNet convention)
ncfile.title = "CMEMS In Situ TAC Example"
ncfile.institution = "Your Institution Name"
ncfile.source = "In situ observations"
ncfile.history = f"Created on {datetime.now(UTC).isoformat()} UTC"
ncfile.references = "http://marine.copernicus.eu"
ncfile.comment = "Prototype NetCDF file for CMEMS In Situ TAC"
ncfile.data_type = "OceanSITES time-series data"
ncfile.Conventions = "CF-1.6, OceanSITES-1.3, SeaDataNet"

# Define dimensions
time_dim = ncfile.createDimension('TIME', None)  # Unlimited
depth_dim = ncfile.createDimension('DEPTH', 1)
lat_dim = ncfile.createDimension('LATITUDE', 1)
lon_dim = ncfile.createDimension('LONGITUDE', 1)

# Define coordinate variables
times = ncfile.createVariable('TIME', 'f8', ('TIME',))
times.units = 'days since 1950-01-01 00:00:00 UTC'
times.standard_name = 'time'
times.long_name = 'Time of measurement'
times.axis = 'T'

depth = ncfile.createVariable('DEPTH', 'f4', ('DEPTH',))
depth.units = 'm'
depth.positive = 'down'
depth.standard_name = 'depth'
depth.long_name = 'Depth of measurement'
depth.axis = 'Z'

lat = ncfile.createVariable('LATITUDE', 'f8', ('LATITUDE',))
lat.units = 'degrees_north'
lat.standard_name = 'latitude'
lat.long_name = 'Latitude'
lat.axis = 'Y'

lon = ncfile.createVariable('LONGITUDE', 'f8', ('LONGITUDE',))
lon.units = 'degrees_east'
lon.standard_name = 'longitude'
lon.long_name = 'Longitude'
lon.axis = 'X'

# Define a sample data variable (e.g., temperature)
temp = ncfile.createVariable('TEMP', 'f4', 
                             ('TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE',),
                             fill_value=-9999.0)

temp.units = 'degree_Celsius'
temp.standard_name = 'sea_water_temperature'
temp.long_name = 'Sea Water Temperature'
temp.coordinates = "TIME DEPTH LATITUDE LONGITUDE"

# Fill data
n_points = 10
base_time = datetime(2023, 1, 1)
times[:] = [(base_time + timedelta(days=i) - datetime(1950, 1, 1)).days for i in range(n_points)]
depth[:] = [10.0]
lat[:] = [42.0]
lon[:] = [3.0]
temp[:, 0, 0, 0] = np.linspace(14.5, 15.5, n_points)  # Mock temperature data

# Close the file
ncfile.close()
print("CMEMS In Situ TAC NetCDF prototype created.")
