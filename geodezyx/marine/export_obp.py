#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 06/12/2025 21:28:55

@author: psakic

Create NetCDF files following CF-1.6 Convention
This script generates ocean bottom pressure and temperature time series data
with support for different instrument types and configurations.

The script supports:
- Simple instruments: 1 pressure sensor, 1 seawater temperature, 1 internal temperature
- Advanced instruments: 2 pressure sensors, 2 internal temperatures, 1 seawater temperature,
                        1 barometer with temperature
"""

import logging
import os

import numpy as np
from netCDF4 import Dataset
from datetime import datetime, timezone
import pandas as pd
from typing import Dict, Optional, List, Tuple
import warnings

from geodezyx import conv, utils

# =============================================================================
# NETCDF EXPORT FUNCTION
# =============================================================================


def export_obp_to_netcdf(
    df_obp: pd.DataFrame,
    column_mapping: Dict[str, Optional[List[str]]],
    output_dir: str,
    metadata_dict: Dict,
    conversion_factors: Optional[Dict[str, float]] = None,
    keep_sensor_dimension: bool = True,
    quality_flags: Optional[np.ndarray] = None,
    force: bool = False,
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
    metadata_dict : dict
        Station metadata configuration containing:
        - station_id: str (MANDATORY)
        - latitude: float (MANDATORY)
        - longitude: float (MANDATORY)
        - depth: float (MANDATORY, meters)
        - institution: str (MANDATORY)
        - source: str (MANDATORY)
        - references: str (MANDATORY)
        - station_name: str (MANDATORY)
        - comment: str (MANDATORY)
        - project: str (MANDATORY)
        - creator_name: str (MANDATORY)
        - creator_email: str (MANDATORY)
        - creator_url: str (MANDATORY)
        - processing_level: str (MANDATORY)
        - summary: str (MANDATORY)
        - title: str (optional)
    conversion_factors : dict, optional
        Conversion factors for each variable type:
        - "pressure": pressure conversion factor (e.g., 0.01 for hPa to dbar)
        - "temperature_seawater": temperature conversion factor
        - "pressure_barometer": barometer pressure conversion factor
    keep_sensor_dimension : bool
        Whether to keep sensor dimension even for single sensor
    quality_flags : np.ndarray, optional
        Quality control flags array (0=good, 1=questionable, 2=bad)
    force : bool
        Whether to overwrite existing file if it already exists

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
            "pressure_barometer": 0.01,
        }

    # Create output directory
    utils.create_dir(output_dir)

    # Extract time information
    colnam_time = column_mapping["time"]
    start_date = df_obp[colnam_time].min()
    end_date = df_obp[colnam_time].max()

    # Calculate sampling interval
    time_diffs = df_obp[colnam_time].diff().dropna()
    if hasattr(time_diffs.iloc[0], "seconds"):
        sampling_interval_seconds = time_diffs.value_counts().index[0].seconds
    else:
        # If already in numeric format
        sampling_interval_seconds = int(time_diffs.mode()[0])

    num_samples = len(df_obp)
    num_days = (end_date - start_date).days

    # Convert time to POSIX timestamps
    val_time = conv.dt2posix(df_obp[colnam_time].dt.to_pydatetime()).values

    #####  Extract data arrays
    def _stdvarnam2val(std_var_name, optional=False):
        """
        Extract and process data values for a standard variable from the dataframe.

        This helper function retrieves column names from the column_mapping dictionary,
        applies conversion factors, and transposes the data to match NetCDF dimension ordering.

        Parameters
        ----------
        std_var_name : str
            Standard variable name key to lookup in column_mapping
            (e.g., "pressure_seafloor", "temperature_sensor")
        optional : bool, optional
            If True, use dict.get() for safe lookup (returns None if key missing).
            If False, use direct dict access (raises KeyError if key missing).
            Default is False.

        Returns
        -------
        colnam_out : list of str or None
            List of column names extracted from the dataframe.
            Returns None if no column names are found.
        val_out : np.ndarray or None
            Processed data values with conversion factor applied and transposed to
            (sensor, time) dimension ordering for CF-1.6 compliance.
            Returns None if input column name is None.

        Notes
        -----
        - Uses `utils.listify()` to convert column name(s) to a list format
        - Applies conversion factors from the `conversion_factors` dict (or 1.0 if not found)
        - Transposes array to place sensor dimension before time dimension
        - Requires access to `df_obp`, `column_mapping`, and `conversion_factors` from outer scope
        """
        # Retrieve column name(s) - use safe get() if marked optional
        if optional:
            colnam = column_mapping.get(std_var_name)
        else:
            colnam = column_mapping[std_var_name]

        # Get conversion factor for this variable, default to 1.0 (no conversion)
        conv_fact = conversion_factors.get(std_var_name, 1)

        # Handle case where no column is defined for this variable
        if colnam is None:
            val_out = None
            colnam_out = None
        else:
            # Convert column name(s) to list format
            colnam_out = utils.listify(colnam)
            # Extract values from dataframe
            val_out = df_obp[colnam_out].values
            # Apply conversion factor (e.g., hPa to dbar)
            val_out = val_out * conv_fact
            # Transpose to (sensor, time) ordering for CF-1.6 compliance
            val_out = np.transpose(val_out)

        return colnam_out, val_out

    colnam_pres, val_pres = _stdvarnam2val("pressure_seafloor")
    colnam_temp_seaw, val_temp_sens = _stdvarnam2val("temperature_sensor")
    colnam_temp_sens, val_temp_seaw = _stdvarnam2val("temperature_seawater")

    # Barometer data (optional)
    colnam_pres_baro, val_pres_baro = _stdvarnam2val("pressure_barometer", True)
    colnam_temp_baro, val_temp_baro = _stdvarnam2val("temperature_barometer", True)

    num_sensors = len(colnam_pres) if utils.is_iterable(colnam_pres) else 1

    # Quality flags (0 = good, 1 = questionable, 2 = bad)
    if quality_flags is None:
        quality_flags = np.zeros(num_samples, dtype=np.int8)

    # Create NetCDF file
    srt_str = conv.dt2str(start_date, "%Y%m%d%H%M%S")
    end_str = conv.dt2str(end_date, "%Y%m%d%H%M%S")
    output_file = f"{metadata_dict['station_id']}_{srt_str}_{end_str}_{sampling_interval_seconds}s.nc"
    output_path = output_dir + "/" + output_file

    if os.path.isfile(output_path):
        if force:
            msg = f"Output file already exists and will be overwritten: {output_path}"
            logging.warn(msg)
            os.remove(output_path)
        else:
            msg = f"Output file already exists. Use force=True to overwrite: {output_path}"
            logging.info(msg)
            return output_path

    nc = Dataset(output_path, "w", format="NETCDF4")

    # Create dimensions
    time_dim = nc.createDimension("time", num_samples)
    if num_sensors > 1 or keep_sensor_dimension:
        sensor_dim = nc.createDimension("sensor", num_sensors)
        dim_sens_tup = ("sensor", "time")
    else:
        dim_sens_tup = ("time",)

    string_dim = nc.createDimension("string_length", 50)

    # Create coordinate variables
    var_time = nc.createVariable("time", "f8", ("time",))
    var_time.units = "seconds since 1970-01-01 00:00:00"
    var_time.long_name = "Time"
    var_time.standard_name = "time"
    var_time.calendar = "standard"
    var_time.axis = "T"
    var_time[:] = val_time

    # Create data variables - Pressure (bottom pressure in dbar)
    var_pres = nc.createVariable(
        "pressure_seafloor", "f4", dim_sens_tup, fill_value=-9999.0
    )
    var_pres.units = "dbar"
    var_pres.long_name = "Bottom Pressure"
    var_pres.standard_name = "sea_water_pressure_at_sea_floor"
    var_pres.positive = "down"
    var_pres.valid_min = np.float32(val_pres.min())
    var_pres.valid_max = np.float32(val_pres.max())
    var_pres.comment = "Sea water pressure at seafloor"
    var_pres[...] = val_pres.astype(np.float32)

    # Seawater Temperature
    var_temp_seaw = nc.createVariable(
        "temperature_seawater", "f4", ("time",), fill_value=-9999.0
    )
    var_temp_seaw.units = "degrees_Celsius"
    var_temp_seaw.long_name = "Seawater Temperature"
    var_temp_seaw.standard_name = "sea_water_temperature"
    var_temp_seaw.valid_min = np.float32(val_temp_seaw.min())
    var_temp_seaw.valid_max = np.float32(val_temp_seaw.max())
    var_temp_seaw.comment = "Sea water temperature at seafloor"
    var_temp_seaw[...] = val_temp_seaw.astype(np.float32)

    # Sensor Internal Temperature
    if val_temp_sens is not None:
        var_temp_sens = nc.createVariable(
            "temperature_sensor", "f4", dim_sens_tup, fill_value=-9999.0
        )
        var_temp_sens.units = "degrees_Celsius"
        var_temp_sens.long_name = "Sensor Internal Temperature"
        var_temp_sens.valid_min = np.float32(val_temp_sens.min())
        var_temp_sens.valid_max = np.float32(val_temp_sens.max())
        var_temp_sens.comment = "Internal temperature of the pressure sensor"
        var_temp_sens[...] = val_temp_sens.astype(np.float32)

    # Barometer pressure (optional)
    if val_pres_baro is not None:
        var_pres_baro = nc.createVariable(
            "pressure_barometer",
            "f4",
            ("time",),
            fill_value=-9999.0,
            zlib=True,
            complevel=4,
        )
        var_pres_baro.units = "hPa"
        var_pres_baro.long_name = "Internal Barometer Pressure"
        var_pres_baro.standard_name = "air_pressure"
        var_pres_baro.valid_min = np.float32(val_pres_baro.min())
        var_pres_baro.valid_max = np.float32(val_pres_baro.max())
        var_pres_baro.comment = (
            "Internal atmospheric pressure for leak detection and quality control"
        )
        var_pres_baro.ancillary_variables = "temperature_barometer"
        var_pres_baro[:] = val_pres_baro.astype(np.float32)

    # Barometer temperature (optional)
    if val_temp_baro is not None:
        var_temp_baro = nc.createVariable(
            "temperature_barometer",
            "f4",
            ("time",),
            fill_value=-9999.0,
            zlib=True,
            complevel=4,
        )
        var_temp_baro.units = "degrees_Celsius"
        var_temp_baro.long_name = "Internal Barometer Temperature"
        var_temp_baro.valid_min = np.float32(val_temp_baro.min())
        var_temp_baro.valid_max = np.float32(val_temp_baro.max())
        var_temp_baro.comment = "Internal temperature of the barometer"
        var_temp_baro[:] = val_temp_baro.astype(np.float32)

    # Quality flags
    var_qc = nc.createVariable("quality_flag", "i1", ("time",))
    var_qc.long_name = "Quality Control Flag"
    var_qc.flag_values = np.array([0, 1, 2], dtype=np.int8)
    var_qc.flag_meanings = "good questionable bad"
    var_qc.comment = "0=good, 1=questionable, 2=bad"
    var_qc[:] = quality_flags

    # Add global attributes
    nc.Conventions = "CF-1.6"
    nc.title = metadata_dict.get(
        "title", f'Ocean Bottom Pressure Data - Station {metadata_dict["station_id"]}'
    )
    nc.history = f'{datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")} - Created using Python netCDF4'
    nc.keywords = "bottom pressure, ocean pressure, BPR, OBP, temperature"
    nc.date_created = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    # Add time coverage attributes
    nc.time_coverage_start = start_date.strftime("%Y-%m-%dT%H:%M:%SZ")
    nc.time_coverage_end = end_date.strftime("%Y-%m-%dT%H:%M:%SZ")
    nc.time_coverage_duration = f"P{num_days}D"
    nc.time_coverage_resolution = f"PT{sampling_interval_seconds}S"

    # Add geospatial attributes
    nc.geospatial_lat_min = metadata_dict.get("latitude", -9999.0)
    nc.geospatial_lat_max = metadata_dict.get("latitude", -9999.0)
    nc.geospatial_lon_min = metadata_dict.get("longitude", -9999.0)
    nc.geospatial_lon_max = metadata_dict.get("longitude", -9999.0)
    nc.geospatial_vertical_min = metadata_dict.get("depth", -9999.0)
    nc.geospatial_vertical_max = metadata_dict.get("depth", -9999.0)
    nc.geospatial_vertical_units = "m"
    nc.geospatial_vertical_positive = "down"

    # Add custom attributes from metadata_dict using a loop
    # Define mandatory attributes
    mandatory_attrs = [
        "station_id",
        "latitude",
        "longitude",
        "depth",
        "institution",
        "source",
        "references",
        "station_name",
        "comment",
        "project",
        "creator_name",
        "creator_email",
        "creator_url",
        "processing_level",
        "summary",
    ]

    # Check for missing mandatory attributes
    for attr in mandatory_attrs:
        if attr not in metadata_dict.keys():
            metadata_dict[attr] = "NOT_PROVIDED"
            warnings.warn(
                f"Missing mandatory attributes in metadata_dict: {attr}. ", UserWarning
            )

    for config_key, attr_val in metadata_dict.items():
        nc.setncattr(config_key, attr_val)

    # Add station location variables
    var_lat = nc.createVariable("latitude", "f8")
    var_lat.units = "degrees_north"
    var_lat.long_name = "Station Latitude"
    var_lat.standard_name = "latitude"
    var_lat.valid_min = -90.0
    var_lat.valid_max = 90.0
    var_lat[:] = metadata_dict.get("latitude", -9999.0)

    var_lon = nc.createVariable("longitude", "f8")
    var_lon.units = "degrees_east"
    var_lon.long_name = "Station Longitude"
    var_lon.standard_name = "longitude"
    var_lon.valid_min = -180.0
    var_lon.valid_max = 180.0
    var_lon[:] = metadata_dict.get("longitude", -9999.0)

    var_depth = nc.createVariable("depth", "f8")
    var_depth.units = "m"
    var_depth.long_name = "Station Depth"
    var_depth.standard_name = "depth"
    var_depth.positive = "down"
    var_depth[:] = metadata_dict.get("depth", -9999.0)

    # Close the file
    nc.close()

    print(f"NetCDF file created: {output_file}")
    print(f"Number of samples: {num_samples}")
    print(f"Time range: {start_date} to {end_date}")
    print(f"Sampling interval: {sampling_interval_seconds} seconds")
    print(
        f"Station location: {metadata_dict.get('latitude', 'N/A')}°N, {metadata_dict.get('longitude', 'N/A')}°E at {metadata_dict.get('depth', 'N/A')}m depth"
    )

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
    temp_seawater_col: str = "temp_sea",
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
    df_temp_sns = df_temp_sns.rename(
        {"val": temp_sensor_col, "cnt": "cnt_temp_sns"}, axis=1
    )
    df_temp_sns.set_index(time_col, inplace=True)
    df_temp_sns[temp_sensor_col] = df_temp_sns[temp_sensor_col].astype(float)

    # Read seawater temperature
    df_temp_sea = pd.read_pickle(path_temp_seawater)
    df_temp_sea = df_temp_sea.rename(
        {"val": temp_seawater_col, "cnt": "cnt_temp_sea"}, axis=1
    )
    df_temp_sea.set_index(time_col, inplace=True)
    df_temp_sea[temp_seawater_col] = pd.to_numeric(
        df_temp_sea[temp_seawater_col], errors="coerce"
    )
    df_temp_sea[temp_seawater_col] = df_temp_sea[temp_seawater_col].astype(float)

    # Merge all dataframes
    df_out = pd.concat((df_pres, df_temp_sns, df_temp_sea), axis=1)
    df_out.reset_index(inplace=True)
    df_out = df_out.loc[df_out[pressure_col].dropna().index]

    # Create column mapping
    column_mapping_out = {
        "time": time_col,
        "pressure_seafloor": pressure_col,
        "temperature_sensor": temp_sensor_col,
        "temperature_seawater": temp_seawater_col,
        "pressure_barometer": None,
        "temperature_barometer": None,
    }

    return df_out, column_mapping_out
