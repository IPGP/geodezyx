# https://chat.deepseek.com/a/chat/s/4ac9a75c-2fe7-4a36-b206-a2c9c6fb6caa

import numpy as np
from netCDF4 import Dataset
import datetime as dt

def create_oceansites_netcdf(output_file):
    """
    Create an OceanSITES-compliant NetCDF file
    
    Parameters:
    output_file (str): Path to the output NetCDF file
    """
    
    # Create the NetCDF file
    rootgrp = Dataset(output_file, 'w', format='NETCDF4')
    
    now = dt.datetime.now(dt.UTC)
    
    # ======================
    # Global Attributes (Metadata)
    # ======================
    rootgrp.title = "OceanSITES data file"
    rootgrp.naming_authority = "OceanSITES"
    rootgrp.id = "SITE_CODE-PLATFORM_CODE-FILE_CODE"
    rootgrp.data_type = "OceanSITES time-series data"
    rootgrp.format_version = "1.4"
    rootgrp.Conventions = "CF-1.6, OceanSITES-1.3"
    rootgrp.date_created = now.strftime("%Y-%m-%dT%H:%M:%SZ")
    rootgrp.history = "Created " + now.strftime("%Y-%m-%d")
    
    # Data source information
    rootgrp.source = "mooring"
    rootgrp.institution = "Your Institution"
    rootgrp.institution_references = "https://www.your-institution.org"
    rootgrp.contact = "user@your-institution.org"
    rootgrp.references = "http://www.oceansites.org"
    
    # Data quality information
    rootgrp.QC_indicator = "mixed"
    rootgrp.QC_procedures = "OceanSITES QC manual v1.2"
    
    # Position information
    rootgrp.geospatial_lat_min = -30.5
    rootgrp.geospatial_lat_max = -30.5
    rootgrp.geospatial_lon_min = 15.2
    rootgrp.geospatial_lon_max = 15.2
    rootgrp.geospatial_vertical_min = 5.0
    rootgrp.geospatial_vertical_max = 500.0
    
    # Time information
    rootgrp.time_coverage_start = "2020-01-01T00:00:00Z"
    rootgrp.time_coverage_end = "2020-01-31T23:59:59Z"
    
    # ======================
    # Dimensions
    # ======================
    time_dim = rootgrp.createDimension('TIME', None)  # Unlimited dimension
    depth_dim = rootgrp.createDimension('DEPTH', 10)  # 10 depth levels
    
    # ======================
    # Coordinate Variables
    # ======================
    
    # TIME variable - note fill_value is set during creation
    time_var = rootgrp.createVariable('TIME', 'f8', ('TIME',), fill_value=-9999.0)
    time_var.standard_name = "time"
    time_var.long_name = "Time"
    time_var.units = "days since 1950-01-01 00:00:00 UTC"
    time_var.calendar = "gregorian"
    time_var.axis = "T"
    time_var.valid_min = 0.0
    time_var.valid_max = 90000.0
    
    # DEPTH variable - note fill_value is set during creation
    depth_var = rootgrp.createVariable('DEPTH', 'f4', ('DEPTH',), fill_value=-9999.0)
    depth_var.standard_name = "depth"
    depth_var.long_name = "Depth below sea level"
    depth_var.units = "meters"
    depth_var.axis = "Z"
    depth_var.positive = "down"
    depth_var.valid_min = 0.0
    depth_var.valid_max = 12000.0
    
    # LATITUDE variable (scalar) - note fill_value is set during creation
    lat_var = rootgrp.createVariable('LATITUDE', 'f4', fill_value=-9999.0)
    lat_var.standard_name = "latitude"
    lat_var.long_name = "Latitude of the measurement"
    lat_var.units = "degrees_north"
    lat_var.valid_min = -90.0
    lat_var.valid_max = 90.0
    
    # LONGITUDE variable (scalar) - note fill_value is set during creation
    lon_var = rootgrp.createVariable('LONGITUDE', 'f4', fill_value=-9999.0)
    lon_var.standard_name = "longitude"
    lon_var.long_name = "Longitude of the measurement"
    lon_var.units = "degrees_east"
    lon_var.valid_min = -180.0
    lon_var.valid_max = 180.0
    
    # ======================
    # Data Variables
    # ======================
    
    # Temperature - fill_value set during creation
    temp_var = rootgrp.createVariable('TEMP', 'f4', ('TIME', 'DEPTH'), fill_value=-9999.0)
    temp_var.standard_name = "sea_water_temperature"
    temp_var.long_name = "Sea temperature"
    temp_var.units = "degree_Celsius"
    temp_var.valid_min = -2.5
    temp_var.valid_max = 40.0
    temp_var.coordinates = "TIME DEPTH LATITUDE LONGITUDE"
    temp_var.QC_indicator = "1"
    temp_var.QC_procedure = "OceanSITES QC manual v1.2"
    
    # Salinity - fill_value set during creation
    psal_var = rootgrp.createVariable('PSAL', 'f4', ('TIME', 'DEPTH'), fill_value=-9999.0)
    psal_var.standard_name = "sea_water_practical_salinity"
    psal_var.long_name = "Practical salinity"
    psal_var.units = "1"
    psal_var.valid_min = 0.0
    psal_var.valid_max = 41.0
    psal_var.coordinates = "TIME DEPTH LATITUDE LONGITUDE"
    psal_var.QC_indicator = "1"
    psal_var.QC_procedure = "OceanSITES QC manual v1.2"
    
    # ======================
    # Fill with example data
    # ======================
    
    # Time data (30 days)
    time_data = np.arange(0, 30, 1)
    time_var[:] = time_data
    
    # Depth data (10 levels from 5 to 500m)
    depth_data = np.linspace(5, 500, 10)
    depth_var[:] = depth_data
    
    # Position data
    lat_var[:] = -30.5
    lon_var[:] = 15.2
    
    # Create temperature and salinity data with some random variation
    temp_data = 15 + 5 * np.random.randn(30, 10)
    psal_data = 35 + 0.5 * np.random.randn(30, 10)
    
    temp_var[:, :] = temp_data
    psal_var[:, :] = psal_data
    
    # ======================
    # Close the file
    # ======================
    rootgrp.close()
    print(f"OceanSITES NetCDF file created: {output_file}")

# Example usage
if __name__ == "__main__":
    create_oceansites_netcdf("/home/psakicki/aaa_FOURBI/out.nc")
    
    
