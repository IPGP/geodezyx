# https://chat.deepseek.com/a/chat/s/b979bb2d-15d3-41e7-880c-736e686e6161

import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime
import uuid

def create_insitu_tac_netcdf(output_path):
    """
    Create a prototype Copernicus Marine In Situ TAC NetCDF file.
    
    Parameters:
    - output_path: str, path where to save the NetCDF file
    """
    
    # Create sample data
    num_obs = 100
    time = pd.date_range(start="2023-01-01", periods=num_obs, freq="H")
    latitude = np.random.uniform(-90, 90, num_obs)
    longitude = np.random.uniform(-180, 180, num_obs)
    depth = np.random.uniform(0, 5000, num_obs)
    temperature = np.random.uniform(-2, 30, num_obs)
    salinity = np.random.uniform(30, 38, num_obs)
    
    # Create platform IDs and types
    platforms = [f"platform_{i:03d}" for i in range(num_obs)]
    platform_types = ["buoy"] * num_obs
    
    # Create QC flags (0=no QC, 1=good, 2=probably good, etc.)
    qc_flags = np.random.choice([0, 1, 2, 3, 4], size=num_obs)
    
    # Create xarray Dataset
    ds = xr.Dataset(
        {
            "LATITUDE": (["OBS"], latitude),
            "LONGITUDE": (["OBS"], longitude),
            "DEPH": (["OBS"], depth),
            "TEMP": (["OBS"], temperature),
            "PSAL": (["OBS"], salinity),
            "POSITION_QC": (["OBS"], qc_flags),
            "TEMP_QC": (["OBS"], qc_flags),
            "PSAL_QC": (["OBS"], qc_flags),
            "PLATFORM_NUMBER": (["OBS"], platforms),
            "PLATFORM_TYPE": (["OBS"], platform_types),
        },
        coords={
            "OBS": np.arange(num_obs),
            "TIME": ("OBS", time),
        }
    )
    
    # Add global attributes (following Copernicus conventions)
    ds.attrs = {
        "title": "Copernicus Marine In Situ TAC prototype",
        "institution": "Your Institution",
        "source": "In situ observations",
        "history": f"{datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')} - File created",
        "references": "http://marine.copernicus.eu",
        "Conventions": "CF-1.8, Copernicus-InSituTAC-2.0",
        "featureType": "trajectory",
        "id": str(uuid.uuid4()),
        "data_mode": "D",  # D for delayed mode, R for real-time
        "format_version": "2.0",
        "netcdf_version": "4",
        "naming_authority": "Copernicus",
        "cdm_data_type": "Trajectory",
        "comment": "Prototype file - not for operational use",
        "geospatial_lat_min": str(latitude.min()),
        "geospatial_lat_max": str(latitude.max()),
        "geospatial_lon_min": str(longitude.min()),
        "geospatial_lon_max": str(longitude.max()),
        "time_coverage_start": str(time[0]),
        "time_coverage_end": str(time[-1]),
        "date_created": datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
    
    # Add variable attributes (modified to avoid calendar conflict)
    ds["TIME"].attrs = {
        "standard_name": "time",
        "long_name": "Time",
        "axis": "T",
    }
    
    ds["LATITUDE"].attrs = {
        "standard_name": "latitude",
        "long_name": "Latitude",
        "units": "degrees_north",
        "valid_min": "-90",
        "valid_max": "90",
        "axis": "Y",
    }
    
    ds["LONGITUDE"].attrs = {
        "standard_name": "longitude",
        "long_name": "Longitude",
        "units": "degrees_east",
        "valid_min": "-180",
        "valid_max": "180",
        "axis": "X",
    }
    
    ds["DEPH"].attrs = {
        "standard_name": "depth",
        "long_name": "Depth",
        "units": "meters",
        "positive": "down",
        "axis": "Z",
    }
    
    ds["TEMP"].attrs = {
        "standard_name": "sea_water_temperature",
        "long_name": "Sea temperature",
        "units": "degree_Celsius",
        "valid_min": "-5",
        "valid_max": "40",
    }
    
    ds["PSAL"].attrs = {
        "standard_name": "sea_water_salinity",
        "long_name": "Practical salinity",
        "units": "1",
        "valid_min": "0",
        "valid_max": "50",
    }
    
    # QC flags attributes
    for var in ["POSITION_QC", "TEMP_QC", "PSAL_QC"]:
        ds[var].attrs = {
            "long_name": f"Quality flag for {var.split('_')[0]}",
            "conventions": "Copernicus Marine In Situ reference table 2",
            "flag_values": "0,1,2,3,4,5,6,7,8,9",
            "flag_meanings": (
                "no_qc_performed good_data probably_good_data probably_bad_data bad_data "
                "changed_value value_below_detection value_in_excess interpolated_value missing_value"
            ),
        }
    
    ds["PLATFORM_NUMBER"].attrs = {
        "long_name": "Platform unique identifier",
        "conventions": "WMO float identifier : A9IIIII",
    }
    
    ds["PLATFORM_TYPE"].attrs = {
        "long_name": "Platform type",
        "conventions": "Copernicus Marine In Situ reference table 25",
    }
    
    # Save to NetCDF with proper encoding
    encoding = {
        "TIME": {"dtype": "float64", "_FillValue": None, "calendar": "gregorian"},
        "LATITUDE": {"dtype": "float32", "_FillValue": None},
        "LONGITUDE": {"dtype": "float32", "_FillValue": None},
        "DEPH": {"dtype": "float32", "_FillValue": None},
        "TEMP": {"dtype": "float32", "_FillValue": -9999.0},
        "PSAL": {"dtype": "float32", "_FillValue": -9999.0},
    }
    
    ds.to_netcdf(output_path, encoding=encoding, format="NETCDF4")
    print(f"File saved to {output_path}")

# Example usage
if __name__ == "__main__":
    create_insitu_tac_netcdf("copernicus_insitu_tac_prototype.nc")