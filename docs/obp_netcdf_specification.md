# NetCDF File Specification for Ocean Bottom Pressure Data

**Version:** 1.0  
**Date:** December 2025  
**Conventions:** CF-1.6

---

## 1. Overview

This document describes the structure and content of NetCDF files containing Ocean Bottom Pressure (OBP) data. These files are designed to store time-series measurements from seafloor pressure sensors, including bottom pressure, seawater temperature, sensor temperature, and optional barometer data for quality control.

**Purpose:** Reference specification for data users, archive managers, and downstream processing tools.

---

## 2. File Naming Convention

```
{STATION_ID}_{START_TIME}_to_{END_TIME}_{INTERVAL}s.nc
```

**Components:**
- `STATION_ID`: Station identifier (e.g., `OCFC`, `A0Ax`)
- `START_TIME`: UTC start timestamp in format `YYYYMMDDhhmmss`
- `END_TIME`: UTC end timestamp in format `YYYYMMDDhhmmss`
- `INTERVAL`: Nominal sampling interval in seconds

**Example:**
```
OCFC_20210315120000_to_20210415235959_1s.nc
```

---

## 3. Data Model

### 3.1 Conventions and Standards

- **Conventions:** CF-1.6 (Climate and Forecast Metadata Conventions)
- **Time reference:** POSIX time (seconds since 1970-01-01 00:00:00 UTC)
- **Calendar:** Gregorian
- **Missing data:** Fill value `-9999.0` for floating-point variables

### 3.2 Dimensions

| Dimension | Type | Description |
|-----------|------|-------------|
| `time` | UNLIMITED | Number of time samples in the dataset |
| `sensor` | Fixed (optional) | Number of sensors (present only when multiple sensors exist or explicitly requested) |
| `string_length` | Fixed (50) | Maximum string length for text variables |

**Note:** The `sensor` dimension appears when:
- Multiple sensor columns are provided in source data, OR
- `keep_sensor_dimension=True` is specified

---

## 4. Coordinate Variables

### 4.1 time

**Description:** Time coordinate for all observations

| Property | Value |
|----------|-------|
| **Variable name** | `time` |
| **Data type** | `float64` |
| **Dimensions** | `(time)` |
| **Units** | `seconds since 1970-01-01 00:00:00` |
| **Calendar** | `gregorian` |
| **Standard name** | `time` |
| **Axis** | `T` |
| **Long name** | `Time` |

**Example values:**
```
[1615809600.0, 1615809601.0, 1615809602.0, ...]
```

---

## 5. Data Variables

### 5.1 pressure_seafloor

**Description:** Sea water pressure at the seafloor (bottom pressure)

| Property | Value |
|----------|-------|
| **Variable name** | `pressure_seafloor` |
| **Data type** | `float32` |
| **Dimensions** | `(time)` or `(time, sensor)` |
| **Units** | `dbar` (decibar) |
| **Standard name** | `sea_water_pressure_at_sea_floor` |
| **Long name** | `Bottom Pressure` |
| **Fill value** | `-9999.0` |
| **Positive** | `down` |
| **Valid range** | Set dynamically from data (`valid_min`, `valid_max` attributes) |
| **Comment** | `Sea water pressure at seafloor` |

**Typical range:** 1000–6000 dbar depending on deployment depth

---

### 5.2 temperature_seawater

**Description:** Ambient seawater temperature at the seafloor

| Property | Value |
|----------|-------|
| **Variable name** | `temperature_seawater` |
| **Data type** | `float32` |
| **Dimensions** | `(time)` or `(time, sensor)` |
| **Units** | `degrees_Celsius` |
| **Standard name** | `sea_water_temperature` |
| **Long name** | `Seawater Temperature` |
| **Fill value** | `-9999.0` |
| **Valid range** | Set dynamically from data (`valid_min`, `valid_max` attributes) |
| **Comment** | `Sea water temperature at seafloor` |

**Typical range:** 2–30 °C depending on location and depth

---

### 5.3 temperature_sensor

**Description:** Internal temperature of the pressure sensor housing

| Property | Value |
|----------|-------|
| **Variable name** | `temperature_sensor` |
| **Data type** | `float32` |
| **Dimensions** | `(time)` or `(time, sensor)` |
| **Units** | `degrees_Celsius` |
| **Long name** | `Sensor Internal Temperature` |
| **Fill value** | `-9999.0` |
| **Valid range** | Set dynamically from data (`valid_min`, `valid_max` attributes) |
| **Comment** | `Internal temperature of the pressure sensor` |

**Purpose:** Sensor diagnostics and drift correction

---

### 5.4 pressure_barometer (optional)

**Description:** Internal barometer pressure for leak detection

| Property | Value |
|----------|-------|
| **Variable name** | `pressure_barometer` |
| **Data type** | `float32` |
| **Dimensions** | `(time)` |
| **Units** | `hPa` (hectopascal) |
| **Standard name** | `air_pressure` |
| **Long name** | `Internal Barometer Pressure` |
| **Fill value** | `-9999.0` |
| **Valid range** | Set dynamically from data |
| **Ancillary variables** | `temperature_barometer` |
| **Comment** | `Internal atmospheric pressure for leak detection and quality control` |
| **Compression** | zlib level 4 |

**Typical range:** 900–1100 hPa (nominal atmospheric pressure)  
**Quality indicator:** Values significantly different from atmospheric pressure indicate potential instrument leak

---

### 5.5 temperature_barometer (optional)

**Description:** Internal barometer temperature

| Property | Value |
|----------|-------|
| **Variable name** | `temperature_barometer` |
| **Data type** | `float32` |
| **Dimensions** | `(time)` |
| **Units** | `degrees_Celsius` |
| **Long name** | `Internal Barometer Temperature` |
| **Fill value** | `-9999.0` |
| **Valid range** | Set dynamically from data |
| **Comment** | `Internal temperature of the barometer` |
| **Compression** | zlib level 4 |

---

### 5.6 quality_flag (optional)

**Description:** Quality control flags for data points

| Property | Value |
|----------|-------|
| **Variable name** | `quality_flag` |
| **Data type** | `int8` |
| **Dimensions** | `(time)` |
| **Long name** | `Quality Control Flag` |
| **Flag values** | `[0, 1, 2]` |
| **Flag meanings** | `good questionable bad` |
| **Comment** | `0=good, 1=questionable, 2=bad` |
| **Default** | `0` (all good if not specified) |

**Flag definitions:**
- `0`: Good data, passed all quality checks
- `1`: Questionable data, may require additional scrutiny
- `2`: Bad data, should not be used for scientific analysis

---

## 6. Station/Geospatial Variables

### 6.1 latitude

| Property | Value |
|----------|-------|
| **Variable name** | `latitude` |
| **Data type** | `float64` |
| **Dimensions** | Scalar |
| **Units** | `degrees_north` |
| **Standard name** | `latitude` |
| **Long name** | `Station Latitude` |
| **Valid range** | -90.0 to 90.0 |

### 6.2 longitude

| Property | Value |
|----------|-------|
| **Variable name** | `longitude` |
| **Data type** | `float64` |
| **Dimensions** | Scalar |
| **Units** | `degrees_east` |
| **Standard name** | `longitude` |
| **Long name** | `Station Longitude` |
| **Valid range** | -180.0 to 180.0 |

### 6.3 depth

| Property | Value |
|----------|-------|
| **Variable name** | `depth` |
| **Data type** | `float64` |
| **Dimensions** | Scalar |
| **Units** | `m` (meters) |
| **Standard name** | `depth` |
| **Long name** | `Station Depth` |
| **Positive** | `down` |

**Note:** Depth is measured from sea surface to seafloor (positive downward)

---

## 7. Global Attributes

### 7.1 Mandatory Attributes

| Attribute | Type | Description | Example |
|-----------|------|-------------|---------|
| `Conventions` | string | Metadata convention | `CF-1.6` |
| `title` | string | Dataset title | `Ocean Bottom Pressure Data - Station OCFC` |
| `station_id` | string | Station identifier | `OCFC` |
| `station_name` | string | Full station name | `HALIOS OBP Station "Fer à Cheval" OCFC` |
| `latitude` | float | Station latitude | `-12.83448` |
| `longitude` | float | Station longitude | `45.35826` |
| `depth` | float | Station depth (m) | `1481.0` |
| `institution` | string | Data producer institution | `Institut de physique du globe de Paris` |
| `source` | string | Data source/network | `REVOSIMA` |
| `references` | string | Reference URLs or DOIs | `https://www.ipgp.fr/...` |
| `comment` | string | Additional information | `Ocean Bottom Pressure Data - MAYOBS29 Deployment` |
| `project` | string | Project name | `REVOSIMA` |
| `creator_name` | string | Data creator name | `IPGP/REVOSIMA` |
| `creator_email` | string | Contact email | `gnss-ovs-ipgp@services.cnrs.fr` |
| `creator_url` | string | Creator website | `https://www.ipgp.fr/` |
| `processing_level` | string | Processing status | `Quality Controlled` or `Raw Data` |
| `summary` | string | Dataset summary | `Ocean Bottom Pressure and Temperature Data...` |

### 7.2 Temporal Coverage Attributes

| Attribute | Format | Description | Example |
|-----------|--------|-------------|---------|
| `time_coverage_start` | ISO8601 | Start time (UTC) | `2021-03-15T12:00:00Z` |
| `time_coverage_end` | ISO8601 | End time (UTC) | `2021-04-15T23:59:59Z` |
| `time_coverage_duration` | ISO8601 duration | Total duration | `P31D` (31 days) |
| `time_coverage_resolution` | ISO8601 duration | Sampling interval | `PT1S` (1 second) |

### 7.3 Geospatial Coverage Attributes

| Attribute | Description | Example |
|-----------|-------------|---------|
| `geospatial_lat_min` | Minimum latitude | `-12.83448` |
| `geospatial_lat_max` | Maximum latitude | `-12.83448` |
| `geospatial_lon_min` | Minimum longitude | `45.35826` |
| `geospatial_lon_max` | Maximum longitude | `45.35826` |
| `geospatial_vertical_min` | Minimum depth (m) | `1481.0` |
| `geospatial_vertical_max` | Maximum depth (m) | `1481.0` |
| `geospatial_vertical_units` | Depth units | `m` |
| `geospatial_vertical_positive` | Vertical direction | `down` |

### 7.4 Provenance Attributes

| Attribute | Description | Example |
|-----------|-------------|---------|
| `history` | Processing history | `2025-12-06T21:30:00Z: Created by export_obp_to_netcdf v1.0` |
| `date_created` | File creation date | `2025-12-06T21:30:00Z` |
| `keywords` | Search keywords | `bottom pressure, ocean pressure, BPR, OBP, temperature` |

**Note:** Missing mandatory attributes are set to `NOT_PROVIDED` and generate warnings

---

## 8. Sensor Dimension Behavior

### 8.1 Single Sensor Configuration

When data contains measurements from a single sensor:

```
Dimensions:
  time: UNLIMITED (e.g., 86400)
  
Variables:
  pressure_seafloor(time)
  temperature_seawater(time)
  temperature_sensor(time)
```

### 8.2 Multi-Sensor Configuration

When data contains measurements from multiple sensors:

```
Dimensions:
  time: UNLIMITED (e.g., 86400)
  sensor: 2
  
Variables:
  pressure_seafloor(time, sensor)
  temperature_sensor(time, sensor)
  temperature_seawater(time)  # if single sensor for seawater
```

**Trigger conditions for multi-sensor mode:**
- Multiple columns provided for a variable in source data
- `keep_sensor_dimension=True` flag set explicitly

---

## 9. Data Quality and Valid Ranges

### 9.1 Valid Range Attributes

Each data variable includes `valid_min` and `valid_max` attributes computed from the actual data range:

```
pressure_seafloor:valid_min = 1450.2
pressure_seafloor:valid_max = 1452.8
```

**Note:** These are descriptive (data range) not prescriptive (QC thresholds)

### 9.2 Fill Values

All floating-point variables use fill value `-9999.0` to represent missing or invalid data.

**Usage in analysis:**
```python
import xarray as xr
ds = xr.open_dataset('OCFC_20210315120000_to_20210415235959_1s.nc')
# Masked arrays automatically handle fill values
pressure = ds['pressure_seafloor'].values  # NaN where fill value occurs
```

### 9.3 Compression

Barometer variables use lossless compression (zlib level 4) to reduce file size while maintaining precision.

---

## 10. Usage Examples

### 10.1 Reading with xarray (Python)

```python
import xarray as xr

# Open dataset
ds = xr.open_dataset('OCFC_20210315120000_to_20210415235959_1s.nc')

# View structure
print(ds)

# Access variables
pressure = ds['pressure_seafloor']
time = ds['time']

# Convert to pandas DataFrame
df = ds.to_dataframe()

# Select time range
subset = ds.sel(time=slice('2021-03-20', '2021-03-25'))

# Basic statistics
print(ds['pressure_seafloor'].mean())
print(ds['temperature_seawater'].std())
```

### 10.2 Reading with netCDF4 (Python)

```python
from netCDF4 import Dataset
import numpy as np

# Open file
nc = Dataset('OCFC_20210315120000_to_20210415235959_1s.nc', 'r')

# Read variables
time = nc.variables['time'][:]
pressure = nc.variables['pressure_seafloor'][:]
lat = nc.variables['latitude'][:]

# Read attributes
station_id = nc.getncattr('station_id')
institution = nc.getncattr('institution')

# Close file
nc.close()
```

### 10.3 Reading with R

```r
library(ncdf4)

# Open file
nc <- nc_open("OCFC_20210315120000_to_20210415235959_1s.nc")

# Read variables
time <- ncvar_get(nc, "time")
pressure <- ncvar_get(nc, "pressure_seafloor")

# Read attributes
station_id <- ncatt_get(nc, 0, "station_id")$value

# Close file
nc_close(nc)
```

---

## 11. CF Compliance Notes

### 11.1 Standard Names

This specification uses CF standard names where applicable:
- `time`
- `latitude`, `longitude`, `depth`
- `sea_water_pressure_at_sea_floor`
- `sea_water_temperature`
- `air_pressure`

### 11.2 Units

All units follow UDUNITS-2 conventions:
- Time: `seconds since 1970-01-01 00:00:00`
- Pressure: `dbar`, `hPa`
- Temperature: `degrees_Celsius`
- Position: `degrees_north`, `degrees_east`
- Depth: `m`

### 11.3 Coordinates

The `time` variable serves as the coordinate variable for all data variables through explicit dimension specification.

---

## 12. Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2025-12-06 | Initial specification document |

---

## 13. Contact and Support

For questions about this file specification:
- **Technical issues:** Contact data creator (see `creator_email` attribute)
- **Format questions:** Refer to CF Conventions documentation at http://cfconventions.org/
- **Software bugs:** Report to geodezyx package maintainers

---

## 14. References

- CF Conventions: http://cfconventions.org/cf-conventions/cf-conventions.html
- NetCDF User Guide: https://www.unidata.ucar.edu/software/netcdf/docs/
- UDUNITS-2: https://www.unidata.ucar.edu/software/udunits/

---

**Document maintained by:** geodezyx development team  
**Last updated:** December 2025

