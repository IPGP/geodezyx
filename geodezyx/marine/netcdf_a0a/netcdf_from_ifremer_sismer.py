#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 26/09/2025 15:00:44

@author: psakic
"""


import pandas as pd
import numpy as np
import xarray as xr

time = pd.date_range("2023-01-01", periods=5, freq="D")

# Pressure measurements for each sensor
pressure1 = np.random.uniform(1000, 1020, size=(5,))
pressure2 = np.random.uniform(1000, 1020, size=(5,))

# Combine into a DataArray: (sensor, time)
pressure_data = np.stack([pressure1, pressure2])

# Create xarray Dataset
ds = xr.Dataset(
    {
        "pressure": (("sensor", "time"), np.stack([pressure1, pressure2])),
    },
    coords={
        "sensor": ["sensor_1", "sensor_2"],
        "time": time,
        "lat": (("sensor",), [45.0, 45.0]),  # same location
        "lon": (("sensor",), [-10.0, -10.0])
    },
    attrs={
        "featureType": "timeSeries",
        "title": "Pressure time series from two sensors at the same location",
        "Conventions": "CF-1.8",
    }
)
ds["pressure"].attrs = {
    "standard_name": "sea_water_pressure",
    "units": "dbar",
    "coordinates": "time depth lat lon"
}
ds["time"].attrs = {"standard_name": "time"}
ds["lat"].attrs = {
    "standard_name": "latitude",
    "units": "degrees_north"
}
ds["lon"].attrs = {
    "standard_name": "longitude",
    "units": "degrees_east"
}
ds.to_netcdf("timeseries_multiple_sensors.nc")
