#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 06/12/2025 21:31:54

@author: psakic
"""

from geodezyx import marine
import pandas as pd

# =============================================================================
# EXAMPLE USAGE / MAIN EXECUTION
# =============================================================================

# Choose which example to run
run_halios = False
run_a0a_rbr = True

# =========================================================================
# Example 1: HALIOS (IPGP)
# =========================================================================
if run_halios:
    # Define file paths
    step = "100"  # "1","10","100"
    p_pres = f"/home/psakicki/GFZ_WORK/IPGP_WORK/REVOSIMA/2510_OBSCOM_OBP/2510_paros_cat/paros_p_{step}s_cat.pkl"
    p_temp_sns = "/home/psakicki/GFZ_WORK/IPGP_WORK/REVOSIMA/2510_OBSCOM_OBP/2510_paros_cat/paros_p_temp_cat.pkl"
    p_temp_sea = "/home/psakicki/GFZ_WORK/IPGP_WORK/REVOSIMA/2510_OBSCOM_OBP/2510_paros_cat/paros_p2t_1s_cat.pkl"

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
    station_config = {
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
        station_config=station_config,
        conversion_factors=conversion_factors,
        keep_sensor_dimension=True
    )

# =========================================================================
# Example 2: A0A RBR (LIENSs)
# =========================================================================
if run_a0a_rbr:
    # Define file path
    p_a0a = "/home/psakicki/GFZ_WORK/IPGP_WORK/REVOSIMA/0110_Pressure_Mayotte/0100_from_Treden/RawData/transfer_3167556_files_8a99b7f8/204657_20210409_1130_data_head.txt"
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
    station_config = {
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
    output_dir = "/home/psakicki/GFZ_WORK/IPGP_WORK/REVOSIMA/2510_OBSCOM_OBP/2510_paros_cat/netcdf_cf"

    # Export to NetCDF
    output_path = marine.export_obp_to_netcdf(
        df_obp=df_a0a,
        column_mapping=column_mapping,
        output_dir=output_dir,
        station_config=station_config,
        conversion_factors=conversion_factors,
        keep_sensor_dimension=True
    )
