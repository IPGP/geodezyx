#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# https://chatgpt.com/c/683e9336-2858-8012-b4a9-3c9308c0bb46

"""
Created on Tue Jun  3 08:23:57 2025

@author: psakicki

"""


from netCDF4 import Dataset, date2num
import numpy as np
import datetime as dt

# Create NetCDF file
ncfile = Dataset('example_oceansites.nc', mode='w', format='NETCDF4_CLASSIC')

# Global attributes (OceanSITES standard)
ncfile.setncatts({
    'title': 'Example OceanSITES Time Series',
    'institution': 'Your Organization',
    'source': 'Mooring',
    'history': 'Created ' + dt.datetime.now(dt.UTC).strftime('%Y-%m-%d %H:%M:%S UTC'),
    'references': 'http://www.oceansites.org',
    'Conventions': 'CF-1.6 OceanSITES-1.3',
    'netcdf_version': 'NETCDF4_CLASSIC',
    'deployment_code': 'example01',
    'id': 'EXAMPLE_001',
    'data_type': 'OceanSITES time-series data',
    'format_version': '1.3',
    'site_code': 'EXAMPLE',
    'platform_code': 'EX01',
    'instrument_number': '001',
    'data_mode': 'R'
})

# Define dimensions
ncfile.createDimension('TIME', None)   # unlimited
ncfile.createDimension('DEPTH', 1)
ncfile.createDimension('LATITUDE', 1)
ncfile.createDimension('LONGITUDE', 1)

# Create variables
times = ncfile.createVariable('TIME', 'f8', ('TIME',))
times.units = 'days since 1950-01-01 00:00:00 UTC'
times.long_name = 'Time'
times.standard_name = 'time'
times.axis = 'T'

depths = ncfile.createVariable('DEPTH', 'f4', ('DEPTH',))
depths.units = 'm'
depths.positive = 'down'
depths.axis = 'Z'
depths.long_name = 'Depth'
depths.standard_name = 'depth'

lat = ncfile.createVariable('LATITUDE', 'f4', ('LATITUDE',))
lat.units = 'degrees_north'
lat.long_name = 'Latitude'
lat.standard_name = 'latitude'
lat.axis = 'Y'

lon = ncfile.createVariable('LONGITUDE', 'f4', ('LONGITUDE',))
lon.units = 'degrees_east'
lon.long_name = 'Longitude'
lon.standard_name = 'longitude'
lon.axis = 'X'

# Example data variable (e.g., temperature)
temp = ncfile.createVariable('TEMP', 'f4', 
                             ('TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE'),
                             fill_value=-9999.0)
temp.units = 'degree_Celsius'
temp.long_name = 'Sea Water Temperature'
temp.standard_name = 'sea_water_temperature'
temp.coordinates = 'TIME DEPTH LATITUDE LONGITUDE'

# Fill dimensions with example data
num_points = 10
start_time = dt.datetime(2023, 1, 1)
time_vals = [start_time + dt.timedelta(hours=6*i) for i in range(num_points)]
times[:] = date2num(time_vals, units=times.units)

depths[:] = [10.0]
lat[:] = [45.0]
lon[:] = [-30.0]

# Add mock temperature data
temp[:, 0, 0, 0] = np.random.uniform(10, 15, size=num_points)

# Close the NetCDF file
ncfile.close()
print("OceanSITES NetCDF file created successfully.")


