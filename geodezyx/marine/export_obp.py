#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 06/12/2025 21:28:55

@author: psakic
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  2 18:51:28 2025
Refactored on Dec 6, 2025

@author: psakicki

Create NetCDF files following CF-1.6 Convention
This script generates ocean bottom pressure and temperature time series data
with support for different instrument types and configurations.

The script supports:
- Simple instruments: 1 pressure sensor, 1 seawater temperature, 1 internal temperature
- Advanced instruments: 2 pressure sensors, 2 internal temperatures, 1 seawater temperature,
                        1 barometer with temperature
"""

import numpy as np
from netCDF4 import Dataset
from datetime import datetime, timezone
import pandas as pd
from typing import Dict, Optional, List, Tuple

from geodezyx import conv, utils

# =============================================================================
# NETCDF EXPORT FUNCTION
# =============================================================================

def export_obp_to_netcdf(
        df_obp: pd.DataFrame,
        column_mapping: Dict[str, Optional[List[str]]],
        output_dir: str,
        station_config: Dict,
        conversion_factors: Optional[Dict[str, float]] = None,
        keep_sensor_dimension: bool = True,
        quality_flags: Optional[np.ndarray] = None
) -> str:
    """
    Export ocean bottom pressure data to CF-1.6 compliant NetCDF file.

    Parameters
    ----------
    df_obp : pd.DataFrame
        dataframe with time column and measurement columns
    column_mapping : dict
        Dictionary mapping standard variable names to dataframe column names:
        - "time": time column name
        - "pressure_seafloor": pressure column name(s) (str or list of str)
        - "temperature_sensor": internal temperature column name(s) (str or list of str)
        - "temperature_seawater": seawater temperature column name(s) (str or list of str)
        - "pressure_barometer": barometer pressure column name(s) (str or list of str, optional)
        - "temperature_barometer": barometer temperature column name(s) (str or list of str , optional)
    output_dir : str
        Output directory path
    station_config : dict
        Station metadata configuration containing:
        - station_id: str
        - latitude: float
        - longitude: float
        - depth: float (meters)
        - institution: str
        - source: str
        - references: str
        - station_name: str
        - comment: str
        - project: str
        - creator_name: str
        - creator_email: str
        - creator_url: str
        - processing_level: str
        - summary: str
    conversion_factors : dict, optional
        Conversion factors for each variable type:
        - "pressure": pressure conversion factor (e.g., 0.01 for hPa to dbar)
        - "temperature_seawater": temperature conversion factor
        - "pressure_barometer": barometer pressure conversion factor
    keep_sensor_dimension : bool
        Whether to keep sensor dimension even for single sensor
    quality_flags : np.ndarray, optional
        Quality control flags array (0=good, 1=questionable, 2=bad)

    Returns
    -------
    output_path : str
        Path to the created NetCDF file
    """
    # Set default conversion factors
    if conversion_factors is None:
        conversion_factors = {
            "pressure": 0.01,  # hPa to dbar
            "temperature_seawater": 1.0,
            "pressure_barometer": 0.01
        }

    # Create output directory
    utils.create_dir(output_dir)

    # Extract time information
    colnam_time = column_mapping["time"]
    start_date = df_obp[colnam_time].min()
    end_date = df_obp[colnam_time].max()

    # Calculate sampling interval
    time_diffs = df_obp[colnam_time].diff().dropna()
    if hasattr(time_diffs.iloc[0], 'seconds'):
        sampling_interval_seconds = time_diffs.value_counts().index[0].seconds
    else:
        # If already in numeric format
        sampling_interval_seconds = int(time_diffs.mode()[0])

    num_samples = len(df_obp)
    num_days = (end_date - start_date).days

    # Convert time to POSIX timestamps
    val_time = conv.dt2posix(df_obp[colnam_time]).values

    # Extract data arrays
    colnam_pres = column_mapping["pressure_seafloor"]
    colnam_temp_sens = column_mapping["temperature_sensor"]
    colnam_temp_seaw = column_mapping["temperature_seawater"]
    colnam_pres_baro = column_mapping.get("pressure_barometer")
    colnam_temp_baro = column_mapping.get("temperature_barometer")

    # Prepare data arrays
    val_pres = df_obp[colnam_pres].values * conversion_factors.get("pressure", 1)
    val_temp_sens = df_obp[colnam_temp_sens].values * conversion_factors.get("temperature_sensor", 1)
    val_temp_seaw = df_obp[colnam_temp_seaw].values * conversion_factors.get("temperature_seawater", 1)

    # Barometer data (optional)
    if colnam_pres_baro is not None and len(colnam_pres_baro) > 0:
        val_pres_baro = df_obp[colnam_pres_baro].values * conversion_factors.get("pressure_barometer", 1)
    else:
        val_pres_baro = None

    if colnam_temp_baro is not None and len(colnam_temp_baro) > 0:
        val_temp_baro = df_obp[colnam_temp_baro].values * conversion_factors.get("temperature_barometer", 1)
    else:
        val_temp_baro = None

    num_sensors = len(colnam_pres)

    # Quality flags (0 = good, 1 = questionable, 2 = bad)
    if quality_flags is None:
        quality_flags = np.zeros(num_samples, dtype=np.int8)

    # Create NetCDF file
    srt_str = conv.dt2str(start_date, "%Y%m%d%H%M%S")
    end_str = conv.dt2str(end_date, "%Y%m%d%H%M%S")
    output_file = f"{station_config['station_id']}_{srt_str}to{end_str}_{sampling_interval_seconds}s.nc"
    output_path = output_dir + "/" + output_file
    nc = Dataset(output_path, 'w', format='NETCDF4')

    # Create dimensions
    time_dim = nc.createDimension('time', num_samples)
    if num_sensors > 1 or keep_sensor_dimension:
        sensor_dim = nc.createDimension('sensor', num_sensors)
        dim_sens_tup = ('time', 'sensor')
    else:
        dim_sens_tup = ('time',)

    string_dim = nc.createDimension('string_length', 50)

    # Create coordinate variables
    var_time = nc.createVariable('time', 'f8', ('time',))
    var_time.units = 'seconds since 1970-01-01 00:00:00'
    var_time.long_name = 'Time'
    var_time.standard_name = 'time'
    var_time.calendar = 'gregorian'
    var_time.axis = 'T'
    var_time[:] = val_time

    # Create data variables - Pressure (bottom pressure in dbar)
    var_pres = nc.createVariable('pressure_seafloor', 'f4', dim_sens_tup, fill_value=-9999.0)
    var_pres.units = 'dbar'
    var_pres.long_name = 'Bottom Pressure'
    var_pres.standard_name = 'sea_water_pressure_at_sea_floor'
    var_pres.positive = 'down'
    var_pres.valid_min = np.float32(val_pres.min())
    var_pres.valid_max = np.float32(val_pres.max())
    var_pres.comment = 'Sea water pressure at seafloor'
    var_pres[...] = val_pres.astype(np.float32)

    # Seawater Temperature
    var_temp_seaw = nc.createVariable('temperature_seawater', 'f4', dim_sens_tup, fill_value=-9999.0)
    var_temp_seaw.units = 'degrees_Celsius'
    var_temp_seaw.long_name = 'Seawater Temperature'
    var_temp_seaw.standard_name = 'sea_water_temperature'
    var_temp_seaw.valid_min = np.float32(val_temp_seaw.min())
    var_temp_seaw.valid_max = np.float32(val_temp_seaw.max())
    var_temp_seaw.comment = 'Sea water temperature at seafloor'
    var_temp_seaw[...] = val_temp_seaw.astype(np.float32)

    # Sensor Internal Temperature
    var_temp_sens = nc.createVariable('temperature_sensor', 'f4', dim_sens_tup, fill_value=-9999.0)
    var_temp_sens.units = 'degrees_Celsius'
    var_temp_sens.long_name = 'Sensor Internal Temperature'
    var_temp_sens.valid_min = np.float32(val_temp_sens.min())
    var_temp_sens.valid_max = np.float32(val_temp_sens.max())
    var_temp_sens.comment = 'Internal temperature of the pressure sensor'
    var_temp_sens[...] = val_temp_sens.astype(np.float32)

    # Barometer pressure (optional)
    if val_pres_baro is not None:
        var_pres_baro = nc.createVariable('pressure_barometer', 'f4',
                                          ('time',), fill_value=-9999.0, zlib=True, complevel=4)
        var_pres_baro.units = 'hPa'
        var_pres_baro.long_name = 'Internal Barometer Pressure'
        var_pres_baro.standard_name = 'air_pressure'
        var_pres_baro.valid_min = np.float32(val_pres_baro.min())
        var_pres_baro.valid_max = np.float32(val_pres_baro.max())
        var_pres_baro.comment = 'Internal atmospheric pressure for leak detection and quality control'
        var_pres_baro.ancillary_variables = 'temperature_barometer'
        var_pres_baro[:] = val_pres_baro.astype(np.float32)

    # Barometer temperature (optional)
    if val_temp_baro is not None:
        var_temp_baro = nc.createVariable('temperature_barometer', 'f4',
                                          ('time',), fill_value=-9999.0, zlib=True, complevel=4)
        var_temp_baro.units = 'degrees_Celsius'
        var_temp_baro.long_name = 'Internal Barometer Temperature'
        var_temp_baro.valid_min = np.float32(val_temp_baro.min())
        var_temp_baro.valid_max = np.float32(val_temp_baro.max())
        var_temp_baro.comment = 'Internal temperature of the barometer'
        var_temp_baro[:] = val_temp_baro.astype(np.float32)

    # Quality flags
    var_qc = nc.createVariable('quality_flag', 'i1', ('time',))
    var_qc.long_name = 'Quality Control Flag'
    var_qc.flag_values = np.array([0, 1, 2], dtype=np.int8)
    var_qc.flag_meanings = 'good questionable bad'
    var_qc.comment = '0=good, 1=questionable, 2=bad'
    var_qc[:] = quality_flags

    # Add global attributes
    nc.Conventions = 'CF-1.6'
    nc.title = station_config.get('title', f'Ocean Bottom Pressure Data - Station {station_config["station_id"]}')
    nc.institution = station_config['institution']
    nc.source = station_config['source']
    nc.history = f'{datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")} - Created using Python netCDF4'
    nc.references = station_config['references']
    nc.station_id = station_config['station_id']
    nc.station_name = station_config['station_name']
    nc.comment = station_config['comment']
    nc.geospatial_lat_min = station_config['latitude']
    nc.geospatial_lat_max = station_config['latitude']
    nc.geospatial_lon_min = station_config['longitude']
    nc.geospatial_lon_max = station_config['longitude']
    nc.geospatial_vertical_min = station_config['depth']
    nc.geospatial_vertical_max = station_config['depth']
    nc.geospatial_vertical_units = 'm'
    nc.geospatial_vertical_positive = 'down'
    nc.time_coverage_start = start_date.strftime('%Y-%m-%dT%H:%M:%SZ')
    nc.time_coverage_end = end_date.strftime('%Y-%m-%dT%H:%M:%SZ')
    nc.time_coverage_duration = f'P{num_days}D'
    nc.time_coverage_resolution = f'PT{sampling_interval_seconds}S'
    nc.creator_name = station_config['creator_name']
    nc.creator_email = station_config['creator_email']
    nc.creator_url = station_config['creator_url']
    nc.project = station_config['project']
    nc.processing_level = station_config['processing_level']
    nc.keywords = 'bottom pressure, ocean pressure, BPR, OBP, temperature'
    nc.summary = station_config['summary']
    nc.date_created = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')

    # Add station location variables
    var_lat = nc.createVariable('latitude', 'f8')
    var_lat.units = 'degrees_north'
    var_lat.long_name = 'Station Latitude'
    var_lat.standard_name = 'latitude'
    var_lat.valid_min = -90.0
    var_lat.valid_max = 90.0
    var_lat[:] = station_config['latitude']

    var_lon = nc.createVariable('longitude', 'f8')
    var_lon.units = 'degrees_east'
    var_lon.long_name = 'Station Longitude'
    var_lon.standard_name = 'longitude'
    var_lon.valid_min = -180.0
    var_lon.valid_max = 180.0
    var_lon[:] = station_config['longitude']

    var_depth = nc.createVariable('depth', 'f8')
    var_depth.units = 'm'
    var_depth.long_name = 'Station Depth'
    var_depth.standard_name = 'depth'
    var_depth.positive = 'down'
    var_depth[:] = station_config['depth']

    # Close the file
    nc.close()

    print(f"NetCDF file created: {output_file}")
    print(f"Number of samples: {num_samples}")
    print(f"Time range: {start_date} to {end_date}")
    print(f"Sampling interval: {sampling_interval_seconds} seconds")
    print(
        f"Station location: {station_config['latitude']}°N, {station_config['longitude']}°E at {station_config['depth']}m depth")

    return output_path

# =============================================================================
# DATA READER FUNCTIONS
# =============================================================================

def read_halios_first_deploy(
        path_pressure: str,
        path_temp_sensor: str,
        path_temp_seawater: str,
        time_col: str = "t",
        pressure_col: str = "pres",
        temp_sensor_col: str = "temp_sns",
        temp_seawater_col: str = "temp_sea"
) -> Tuple[pd.DataFrame, Dict[str, Optional[List[str]]]]:
    """
    Read data from a simple HALIOS instrument (1 pressure + 2 temperature sensors).

    Parameters
    ----------
    path_pressure : str
        Path to pressure data pickle file
    path_temp_sensor : str
        Path to sensor internal temperature pickle file
    path_temp_seawater : str
        Path to seawater temperature pickle file
    time_col : str
        Name of the time column in the raw data
    pressure_col : str
        Name of the pressure column in the raw data
    temp_sensor_col : str
        Name of the sensor temperature column in the raw data
    temp_seawater_col : str
        Name of the seawater temperature column in the raw data

    Returns
    -------
    df_out : pd.DataFrame
        Clean dataframe with time and data columns
    column_mapping_out : dict
        Dictionary mapping standard names to actual column names in the dataframe
    """
    # Read pressure data
    df_pres = pd.read_pickle(path_pressure)
    df_pres = df_pres.rename({"val": pressure_col, "cnt": "cnt_pres"}, axis=1)
    df_pres.set_index(time_col, inplace=True)
    df_pres[pressure_col] = df_pres[pressure_col].astype(float)

    # Read sensor internal temperature
    df_temp_sns = pd.read_pickle(path_temp_sensor)
    df_temp_sns = df_temp_sns.rename({"val": temp_sensor_col, "cnt": "cnt_temp_sns"}, axis=1)
    df_temp_sns.set_index(time_col, inplace=True)
    df_temp_sns[temp_sensor_col] = df_temp_sns[temp_sensor_col].astype(float)

    # Read seawater temperature
    df_temp_sea = pd.read_pickle(path_temp_seawater)
    df_temp_sea = df_temp_sea.rename({"val": temp_seawater_col, "cnt": "cnt_temp_sea"}, axis=1)
    df_temp_sea.set_index(time_col, inplace=True)
    df_temp_sea[temp_seawater_col] = pd.to_numeric(df_temp_sea[temp_seawater_col], errors='coerce')
    df_temp_sea[temp_seawater_col] = df_temp_sea[temp_seawater_col].astype(float)

    # Merge all dataframes
    df_out = pd.concat((df_pres, df_temp_sns, df_temp_sea), axis=1)
    df_out.reset_index(inplace=True)
    df_out = df_out.loc[df_out[pressure_col].dropna().index]

    # Create column mapping
    column_mapping_out = {
        "time": time_col,
        "pressure_seafloor": [pressure_col],
        "temperature_sensor": [temp_sensor_col],
        "temperature_seawater": [temp_seawater_col],
        "pressure_barometer": None,
        "temperature_barometer": None
    }

    return df_out, column_mapping_out
