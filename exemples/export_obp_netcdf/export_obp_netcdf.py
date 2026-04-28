#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 06/12/2025 21:31:54

@author: psakic
"""

from geodezyx import marine, conv
import pandas as pd

# =============================================================================
# EXAMPLE USAGE / MAIN EXECUTION
# =============================================================================

# Choose which example to run
run_halios = False
run_a0a_rbr = True
run_emso_azores = False

# =========================================================================
# Example 1: HALIOS (IPGP)
# =========================================================================
if run_halios:
    # Define file paths
    step = "1"  # "1","10","100"
    p_pres = f"/home/sakic/IPGP_WORK/REVOSIMA/2510_OBSCOM_OBP/2510_paros_cat/paros_p_{step}s_cat.pkl"
    p_temp_sns = "/home/sakic/IPGP_WORK/IPGP_WORK/REVOSIMA/2510_OBSCOM_OBP/2510_paros_cat/paros_p_temp_cat.pkl"
    p_temp_sea = "/home/sakic/IPGP_WORK/IPGP_WORK/REVOSIMA/2510_OBSCOM_OBP/2510_paros_cat/paros_p2t_1s_cat.pkl"

    # Read data using dedicated function
    df, column_mapping = marine.read_halios_first_deploy(
        path_pressure=p_pres,
        path_temp_sensor=p_temp_sns,
        path_temp_seawater=p_temp_sea,
        time_col="t",
        pressure_col="pres",
        temp_sensor_col="temp_sns",
        temp_seawater_col="temp_sea"
    )

    # Define station configuration
    metadata_dict = {
        'station_id': 'OCFC',
        'station_name': 'HALIOS OBP Station "Fer à Cheval" OCFC',
        'latitude': -12.83448,
        'longitude': 45.35826,
        'depth': 1481.0,
        'institution': 'Institut de physique du globe de Paris',
        'source': 'Réseau de surveillance volcanologique et sismologique de Mayotte (REVOSIMA)',
        'references': 'https://www.ipgp.fr/observation/infrastructures-nationales-hebergees/revosima/',
        'comment': 'Ocean Bottom Pressure Data - MAYOBS29 Deployment',
        'project': 'Réseau de surveillance volcanologique et sismologique de Mayotte (REVOSIMA)',
        'creator_name': 'IPGP/REVOSIMA',
        'creator_email': 'gnss-ovs-ipgp@services.cnrs.fr',
        'creator_url': 'https://www.ipgp.fr/',
        'processing_level': 'Quality Controlled',
        'summary': 'Ocean Bottom Pressure and Temperature Data from REVOSIMA'
    }

    # Define conversion factors
    conversion_factors = {
        'pressure': 0.01,  # hPa to dbar
        'temperature_seawater': 0.001,  # milli-degrees to degrees
        'pressure_barometer': 0.01
    }

    # Output directory
    output_dir = "/home/psakicki/GFZ_WORK/IPGP_WORK/REVOSIMA/2510_OBSCOM_OBP/2510_paros_cat/netcdf_cf"

    # Export to NetCDF
    output_path = marine.export_obp_to_netcdf(
        df_obp=df,
        column_mapping=column_mapping,
        output_dir=output_dir,
        metadata_dict=metadata_dict,
        conversion_factors=conversion_factors,
        keep_sensor_dimension=True
    )

# =========================================================================
# Example 2: A0A RBR (LIENSs)
# =========================================================================
if run_a0a_rbr:
    # Define file path
    p_a0a = "/home/sakic/IPGP_WORK/REVOSIMA/0110_Pressure_Mayotte/0100_from_Treden/RawData/transfer_3167556_files_8a99b7f8/204657_20210409_1130_data.txt"
    p_a0a = "/home/sakic/Downloads/204657_20210919_1458_data.txt"
    df_a0a = pd.read_csv(p_a0a, sep=",")
    df_a0a["Time"] = pd.to_datetime(df_a0a["Time"])

    # Create column mapping
    column_mapping = {
        "time": "Time",
        "pressure_seafloor": ["BPR pressure", "BPR pressure.1"],
        "temperature_sensor": ["BPR temperature", "BPR temperature.1"],
        "temperature_seawater": ["Temperature"],
        "pressure_barometer": ["Barometer pressure"],
        "temperature_barometer": ["Barometer temperature"]
    }

    # Define station configuration
    metadata_dict = {
        'station_id': 'A0Ax',
        'station_name': 'Ocean Bottom Pressure Station A0Ax',
        'latitude': 0.000,
        'longitude': 0.000,
        'depth': 3000.0,
        'institution': 'Institut de physique du globe de Paris',
        'source': 'Réseau de surveillance volcanologique et sismologique de Mayotte (REVOSIMA)',
        'references': 'https://www.ipgp.fr/observation/infrastructures-nationales-hebergees/revosima/',
        'comment': 'Ocean Bottom Pressure Data',
        'project': 'Réseau de surveillance volcanologique et sismologique de Mayotte (REVOSIMA)',
        'creator_name': 'LIENSs/REVOSIMA',
        'creator_email': 'xxxxx@xxxxxxx.fr',
        'creator_url': 'https://lienss.univ-larochelle.fr/',
        'processing_level': 'Raw Data',
        'summary': 'Ocean Bottom Pressure and Temperature Data from A0A OBP Station'
    }

    # Define conversion factors
    conversion_factors = {
        "pressure": 0.01,  # hPa to dbar
        "temperature_seawater": 1.0,  # No conversion needed
        "temperature_sensor": 1.0,  # No conversion needed
        "pressure_barometer": 0.01,
    }

    # Output directory
    output_dir = "/home/sakic/IPGP_WORK/REVOSIMA/2510_OBSCOM_OBP/2510_paros_cat/netcdf_cf"

    # Export to NetCDF
    output_path = marine.export_obp_to_netcdf(
        df_obp=df_a0a,
        column_mapping=column_mapping,
        output_dir=output_dir,
        metadata_dict=metadata_dict,
        conversion_factors=conversion_factors,
    )

# =========================================================================
# Example 3: EMSO Azores
# =========================================================================

if run_emso_azores:
    ### step 3.1:  read the data as a Pandas DataFrame
    p = "/home/sakic/IPGP_WORK/REVOSIMA/0110_Pressure_Mayotte/0110_RawData/090_MOMAR_for_exemple/2017_JPPW_SBE27_2016_2017.csv"
    df = pd.read_csv(p, sep="\t")
    df.columns = df.columns.str.strip()
    df["datetime"] = pd.to_datetime(df["datetime"])

    ### step 3.2: Create column mapping
    column_mapping = {
        "time": "datetime",
        "pressure_seafloor": "pressure",
        "temperature_seawater": "temperature",
        "temperature_sensor": "temperature"
    }

    ### step 3.3: Define conversion factors
    conversion_factors = {
        "pressure": 0.01,  # hPa to dbar
        "temperature_seawater": 1.0,  # No conversion needed
        }


    # step 3.5:  Define station configuration
    metadata_dict = {
        'station_id': 'JPPW',
        'station_name': 'JPP West site',
        'latitude': conv.dms2degdec_num(-37, 17.559)  ,
        'longitude': 32.281416,
        'depth': 1729.0,
        'institution': 'CNRS/LIENSs',
        'source': 'EMSO Azores',
        'references': 'xxxx',
        'comment': 'Ocean Bottom Pressure Data',
        'project': 'EMSO Azores',
        'creator_name': 'CNRS/LIENSs',
        'creator_email': 'xxxxx@xxxxxxx.fr',
        'creator_url': 'https://lienss.univ-larochelle.fr/',
        'processing_level': 'Raw Data',
        'summary': 'Ocean Bottom Pressure and Temperature Data from SBE53 pressure sensor'
    }


    ### step 3.6: Export to NetCDF

    # Output directory
    output_dir = "/home/sakic/IPGP_WORK/REVOSIMA/0110_Pressure_Mayotte/0110_RawData/090_MOMAR_for_exemple/output"

    output_path = marine.export_obp_to_netcdf(
        df_obp=df,
        column_mapping=column_mapping,
        output_dir=output_dir,
        metadata_dict=metadata_dict,
        conversion_factors=conversion_factors,
    )


# =========================================================================
# Final test: Load and print the exported NetCDF file for verification
# =========================================================================

test_netcdf = False

if test_netcdf:
    import xarray as xr
    netcdf = xr.open_dataset(output_path)
    df_netcdf = netcdf.to_dataframe()

    print(df_netcdf.to_string())
