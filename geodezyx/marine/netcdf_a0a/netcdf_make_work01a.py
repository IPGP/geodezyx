#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 22:07:01 2025

@author: psakic
"""


# https://chat.deepseek.com/a/chat/s/b979bb2d-15d3-41e7-880c-736e686e6161

from geodezyx import marine
import netCDF4
import numpy as np
import datetime as dt

ptxt = r"/home/psakicki/GFZ_WORK/IPGP_WORK/REVOSIMA/0110_Pressure_Mayotte/0110_RawData/010_A0A_RBR/204657_20210919_1458textexport/204657_20210919_1458textexport\204657_20210919_1458textexport_data.txt"
df = marine.read_rbr_txt_data(ptxt)
df.reset_index(inplace=True)

# Create sample data
time = df["Time"].values
num_obs = len(time)

latitude_scal = 0.
longitude_scal = 0.
depth_scal = 3000.

pres_1 = df['BPR pressure 1'].values
temp_1 = df['BPR temperature 1'].values
pres_2 = df['BPR pressure 2'].values
temp_2 = df['BPR temperature 2'].values
pres_baro = df['Barometer pressure'].values
temp_baro = df['Barometer temperature'].values

temp_tup = (temp_1, temp_2)
pres_tup = (pres_1, pres_2)

latitude = np.array([latitude_scal])
longitude = np.array([longitude_scal])
depth = np.array([depth_scal])

nsensor = len(pres_tup)

# ***************************************************************************
# Create the NetCDF file
output_file = "./output_netcdf/out.nc"
ncdfobj = netCDF4.Dataset(output_file, 'w', format='NETCDF4')

now = dt.datetime.now(dt.UTC)

# ======================
# Global Attributes (Metadata)
# ======================
# Attributes  must follow this document:
# Copernicus Marine in situ NetCDF Attributes list
# https://doi.org/10.13155/95044


ncdfobj.title = "OceanSITES data file"
ncdfobj.naming_authority = "OceanSITES"
ncdfobj.id = "SITE_CODE-PLATFORM_CODE-FILE_CODE"
ncdfobj.data_type = "OceanSITES time-series data"
ncdfobj.format_version = "1.4"
ncdfobj.Conventions = "CF-1.6, OceanSITES-1.3"
ncdfobj.date_created = now.strftime("%Y-%m-%dT%H:%M:%SZ")
ncdfobj.history = "Created " + now.strftime("%Y-%m-%d")
# ncdfobj.net_cdf_version

# Data source information
ncdfobj.source = "mooring"
ncdfobj.source_platform_category_code = "mooring"
ncdfobj.institution = "Your Institution"
ncdfobj.institution_references = "https://www.your-institution.org"
ncdfobj.contact = "user@your-institution.org"
ncdfobj.references = "http://www.oceansites.org"

# Data quality information
ncdfobj.QC_indicator = "mixed"
ncdfobj.QC_procedures = "OceanSITES QC manual v1.2"

# Position information
ncdfobj.cdm_data_type = "void"
ncdfobj.geospatial_lat_min = np.min(latitude)
ncdfobj.geospatial_lat_max = np.max(latitude)
ncdfobj.geospatial_lon_min = np.min(longitude)
ncdfobj.geospatial_lon_max = np.max(longitude)
ncdfobj.geospatial_vertical_min = 0.
ncdfobj.geospatial_vertical_max = np.max(depth)

# Time information
ncdfobj.time_coverage_start = str(np.min(time)) + "Z"
ncdfobj.time_coverage_end = str(np.max(time)) + "Z"

# Publication
ncdfobj.citation = "void"
ncdfobj.data_assembly_center = "void"
ncdfobj.references = "void"
ncdfobj.update_interval = "void"


# ======================
# Dimensions
# ======================
time_dim = ncdfobj.createDimension('TIME', len(time))
depth_dim = ncdfobj.createDimension('DEPTH', len(depth))
time_dim = ncdfobj.createDimension('LATITUDE', len(latitude))
depth_dim = ncdfobj.createDimension('LONGITUDE', len(longitude))
sensor_dim = ncdfobj.createDimension('SENSOR', nsensor)

# ======================
# Coordinate Variables
# ======================

# TIME variable - note fill_value is set during creation
time_var = ncdfobj.createVariable('TIME', 'f8', ('TIME',), fill_value=-9999.0)
time_var.standard_name = "time"
time_var.long_name = "Time"
time_var.units = "days since 1950-01-01 00:00:00 UTC"
time_var.calendar = "gregorian"
time_var.axis = "T"
time_var.valid_min = 0.0
time_var.valid_max = 90000.0

#DEPH variable - note fill_value is set during creation
deph_var = ncdfobj.createVariable('DEPH', 'f4', ('DEPTH',), fill_value=-9999.0)
deph_var.standard_name = "depth"
deph_var.long_name = "Depth below sea level"
deph_var.units = "meters"
deph_var.axis = "Z"
deph_var.positive = "down"
deph_var.valid_min = 0.0
deph_var.valid_max = 12000.0

# LATITUDE variable (scalar) - note fill_value is set during creation
lat_var = ncdfobj.createVariable('LATITUDE', 'f4', fill_value=-9999.0)
lat_var.standard_name = "latitude"
lat_var.long_name = "Latitude of the measurement"
lat_var.units = "degrees_north"
lat_var.valid_min = -90.0
lat_var.valid_max = 90.0

# LONGITUDE variable (scalar) - note fill_value is set during creation
lon_var = ncdfobj.createVariable('LONGITUDE', 'f4', fill_value=-9999.0)
lon_var.standard_name = "longitude"
lon_var.long_name = "Longitude of the measurement"
lon_var.units = "degrees_east"
lon_var.valid_min = -180.0
lon_var.valid_max = 180.0

# ======================
# Data Variables
# ======================

# Variables must follow this document:
# Copernicus Marine In Situ TAC Parameters list
# https://doi.org/10.13155/53381

dim_tup = ('TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE', 'SENSOR')
dim_tup_str = " ".join(dim_tup)

# Temperature - fill_value set during creation
temp_var = ncdfobj.createVariable('TEMP' , 'f4',
                                  dim_tup,
                                  fill_value=-9999.0)
temp_var.standard_name = "sea_water_temperature"
temp_var.long_name = "Sea water temperature"
temp_var.units = "degrees_C"
temp_var.valid_min = -2.5
temp_var.valid_max = 40.0
temp_var.coordinates = dim_tup_str
temp_var.QC_indicator = "1"
temp_var.QC_procedure = "OceanSITES QC manual v1.2"
temp_var[:] = np.column_stack(temp_tup)


# Pressure - fill_value set during creation
pres_var = ncdfobj.createVariable('PRES' , 'f4',
                                  dim_tup,
                                  fill_value=-9999.0)
pres_var.standard_name = "sea_water_pressure_at_sea_floor"
pres_var.long_name = "Sea water pressure at sea floor"
pres_var.units = "dbar"
pres_var.valid_min = 0.0
pres_var.valid_max = 41.0
pres_var.coordinates = dim_tup_str
pres_var.QC_indicator = "1"
pres_var.QC_procedure = "OceanSITES QC manual v1.2"
pres_var[:] = np.column_stack(pres_tup)

# ======================
# Fill with example data
# ======================

# Time data
time_var[:] = time

# Depth data
deph_var[:] = depth

# Position data
lat_var[:] = latitude
lon_var[:] = longitude

# Create temperature and salinity data with some random variation

# ======================
# Close the file
# ======================
ncdfobj.close()
