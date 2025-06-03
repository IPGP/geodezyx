# https://chat.deepseek.com/a/chat/s/b979bb2d-15d3-41e7-880c-736e686e6161

from geodezyx import *              # Import the GeodeZYX modules
from geodezyx.externlib import *    # Import the external modules

import netCDF4


#def create_insitu_tac_netcdf(output_path):
# """
# Create a prototype Copernicus Marine In Situ TAC NetCDF file.

# Parameters:
# - output_path: str, path where to save the NetCDF file
# """


ptxt = r"/home/psakicki/GFZ_WORK/IPGP_WORK/REVOSIMA/0110_Pressure_Mayotte/0110_RawData/010_A0A_RBR/204657_20210919_1458textexport/204657_20210919_1458textexport\204657_20210919_1458textexport_data.txt_head"
df = marine.read_rbr_txt_data(ptxt)
df.reset_index(inplace=True)


# Create sample data
#time = pd.date_range(start="2023-01-01", periods=num_obs, freq="H")
time = df["Time"].values
num_obs = len(time)

latitude = [0.]
longitude =  [0.]
depth = [3000.]

pres_1 = df['BPR pressure 1'].values
temp_1 = df['BPR temperature 1'].values
pres_2 = df['BPR pressure 2'].values
temp_2 = df['BPR temperature 2'].values
pres_3 = df['Barometer pressure'].values
temp_3 = df['Barometer temperature'].values



### 


# Create the NetCDF file
output_file = "./output_netcdf/out.nc"
ncdfobj = netCDF4.Dataset(output_file, 'w', format='NETCDF4')

now = dt.datetime.now(dt.UTC)

# ======================
# Global Attributes (Metadata)
# ======================
ncdfobj.title = "OceanSITES data file"
ncdfobj.naming_authority = "OceanSITES"
ncdfobj.id = "SITE_CODE-PLATFORM_CODE-FILE_CODE"
ncdfobj.data_type = "OceanSITES time-series data"
ncdfobj.format_version = "1.4"
ncdfobj.Conventions = "CF-1.6, OceanSITES-1.3"
ncdfobj.date_created = now.strftime("%Y-%m-%dT%H:%M:%SZ")
ncdfobj.history = "Created " + now.strftime("%Y-%m-%d")

# Data source information
ncdfobj.source = "mooring"
ncdfobj.institution = "Your Institution"
ncdfobj.institution_references = "https://www.your-institution.org"
ncdfobj.contact = "user@your-institution.org"
ncdfobj.references = "http://www.oceansites.org"

# Data quality information
ncdfobj.QC_indicator = "mixed"
ncdfobj.QC_procedures = "OceanSITES QC manual v1.2"

# Position information
ncdfobj.geospatial_lat_min = -30.5
ncdfobj.geospatial_lat_max = -30.5
ncdfobj.geospatial_lon_min = 15.2
ncdfobj.geospatial_lon_max = 15.2
ncdfobj.geospatial_vertical_min = 5.0
ncdfobj.geospatial_vertical_max = 500.0

# Time information
ncdfobj.time_coverage_start = "2020-01-01T00:00:00Z"
ncdfobj.time_coverage_end = "2020-01-31T23:59:59Z"

# ======================
# Dimensions
# ======================
time_dim = ncdfobj.createDimension('TIME', len(time))  # Unlimited dimension
depth_dim = ncdfobj.createDimension('DEPTH', len(depth))  # 10 depth levels

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

# DEPTH variable - note fill_value is set during creation
depth_var = ncdfobj.createVariable('DEPTH', 'f4', ('DEPTH',), fill_value=-9999.0)
depth_var.standard_name = "depth"
depth_var.long_name = "Depth below sea level"
depth_var.units = "meters"
depth_var.axis = "Z"
depth_var.positive = "down"
depth_var.valid_min = 0.0
depth_var.valid_max = 12000.0

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

# Temperature - fill_value set during creation
temp_var = ncdfobj.createVariable('TEMP', 'f4', ('TIME', 'DEPTH'), fill_value=-9999.0)
temp_var.standard_name = "sea_water_temperature"
temp_var.long_name = "Sea temperature"
temp_var.units = "degree_Celsius"
temp_var.valid_min = -2.5
temp_var.valid_max = 40.0
temp_var.coordinates = "TIME DEPTH LATITUDE LONGITUDE"
temp_var.QC_indicator = "1"
temp_var.QC_procedure = "OceanSITES QC manual v1.2"

# Pressure - fill_value set during creation
pres_var = ncdfobj.createVariable('PRES', 'f4', ('TIME', 'DEPTH'), fill_value=-9999.0)
pres_var.standard_name = "change_me"
pres_var.long_name = "change_me"
pres_var.units = "1"
pres_var.valid_min = 0.0
pres_var.valid_max = 41.0
pres_var.coordinates = "TIME DEPTH LATITUDE LONGITUDE"
pres_var.QC_indicator = "1"
pres_var.QC_procedure = "OceanSITES QC manual v1.2"

# ======================
# Fill with example data
# ======================

# Time data (30 days)
time_var[:] = time

# Depth data (10 levels from 5 to 500m)
depth_var[:] = depth

# Position data
lat_var[:] = -30.5
lon_var[:] =  15.2

# Create temperature and salinity data with some random variation
temp_var[:, :] = temp_1
pres_var[:, :] = pres_1

# ======================
# Close the file
# ======================
ncdfobj.close()

