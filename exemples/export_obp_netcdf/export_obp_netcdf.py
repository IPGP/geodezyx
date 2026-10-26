#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 06/12/2025 21:31:54

@author: psakic, rjoyard
"""

from geodezyx import marine, conv
import pandas as pd
import export_obp
import os

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
    p_temp_sns = "/home/sakic/IPGP_WORK/REVOSIMA/2510_OBSCOM_OBP/2510_paros_cat/paros_p_temp_cat.pkl"
    p_temp_sea = "/home/sakic/IPGP_WORK/REVOSIMA/2510_OBSCOM_OBP/2510_paros_cat/paros_p2t_1s_cat.pkl"

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
        'pressure_seafloor': 0.01,  # hPa to dbar
        'temperature_seawater': 0.001,  # milli-degrees to degrees
        'pressure_barometer': 0.01
    }

    # Output directory
    output_dir = "/home/sakic/IPGP_WORK/REVOSIMA/2510_OBSCOM_OBP/2510_paros_cat/netcdf_cf"

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
    p_a0a = "/Users/rjoyar01/Desktop/mayotte/A0A_MAY30_R/208295_20240921_1633_data.txt"
    p_metadata = "/Users/rjoyar01/Desktop/mayotte/A0A_MAY30_R/208295_20240921_1633_metadata.txt"
    df_a0a = pd.read_csv(p_a0a, sep=",")
    df_a0a["Time"] = pd.to_datetime(df_a0a["Time"])
    (
        sensor_ids,
        instrument_serial,
        instrument_model,
        firmware_version,
        ruskin_version,
        sampling_mode,
    ) = export_obp.get_metadata_info(p_metadata)

    # Create column mapping
    column_mapping = {
        "time": "Time",
        "pressure_sensor": ["BPR pressure", "BPR pressure.1"],
        "temperature_sensor": ["BPR temperature", "BPR temperature.1"],
        "temperature_seawater": ["Temperature"],
        "pressure_barometer": ["Barometer pressure"],
        "temperature_barometer": ["Barometer temperature"]
    }

    # Define station configuration
    metadata_dict = {
        'ruskin_version': ruskin_version,
        'source_data': os.path.basename(p_a0a),
        'source_metadata': os.path.basename(p_metadata),
        'campaign_deployment': 'MAYOBS25',
        'doi_deployment': '10.17600/18003404',
        'campaign_recovery': 'MAYOBS30',
        'doi_recovery': '10.17600/18003592',
        'station_id': 'A0AR',
        'station_name': 'Ocean Bottom Pressure Station A0A Reference',
        'sensor_ids': {
            'BPR1': sensor_ids[0],
            'BPR2': sensor_ids[1],
            'Barometer': sensor_ids[2]
        },
        'instrument_model': instrument_model,
        'rbr_instrument_serial': instrument_serial,
        'firmware_version': firmware_version,
        'sampling_mode': sampling_mode,
        'latitude': -12.695,
        'longitude': 45.426,
        'depth': 1730.0,
        'institution': 'LIttoral ENvironnement et Sociétés, UMR7266 CNRS, La Rochelle Université',
        'project': 'Réseau de surveillance volcanologique et sismologique de Mayotte (REVOSIMA)',
        'references': 'https://www.ipgp.fr/observation/infrastructures-nationales-hebergees/revosima/',
        'comment': 'No quality control or scientific processing has been applied.',
        'creator_name': 'LIENSs',
        'creator_email': 'contact_a0a_pin@univ-lr.fr',
        'creator_url': 'https://lienss.univ-larochelle.fr/',
        'processing_level': 'Raw Data',
        'summary': 'Ocean Bottom Pressure and Temperature Data from A0A OBP Station'
    }

    # Define conversion factors
    conversion_factors = {
        "pressure_sensor": 1.0, # 0.01,       # hPa -> dbar
        "temperature_seawater": 1.0,
        "temperature_sensor": 1.0,
        "pressure_barometer": 1.0,       # hPa -> hPa
    }

    # Output directory
    output_dir = "/Users/rjoyar01/Desktop/mayotte/A0A_MAY30_R/test_netcdf"

    # Export to NetCDF
    output_path = export_obp.export_obp_to_netcdf(
        df_obp=df_a0a,
        column_mapping=column_mapping,
        output_dir=output_dir,
        metadata_dict=metadata_dict,
        conversion_factors=conversion_factors,
        force=True
    )

# =========================================================================
# Example 3: EMSO Azores
# =========================================================================

if run_emso_azores:
    ### step 3.1:  read the data as a Pandas DataFrame
    p = "/Users/rjoyar01/Downloads/2023_JPPW_SBE26_2022_2023.tid"
    df = export_obp.read_emso_tid(p)

    ### step 3.2: Create column mapping
    column_mapping = {
        "time": "datetime",
        "pressure_sensor": "pressure",
        "temperature_seawater": "temperature",
        "temperature_sensor": "temperature"
    }

    conversion_factors = {
        "pressure_sensor": 1.0,
        "temperature_seawater": 1.0,
        "temperature_sensor": 1.0,
    }

    # step 3.5:  Define station configuration
    metadata_dict = {
        'source_data': os.path.basename(p),
        'station_id': 'JPPW',
        'station_name': 'JPP West site',
        'latitude': conv.dms2degdec_num(-37, 17.559)  ,
        'longitude': 32.281416,
        'depth': 1729.0,
        'institution': 'LIttoral ENvironnement et Sociétés, UMR7266 CNRS, La Rochelle Université',
        'project': 'EMSO Azores',
        'references': 'xxxx',
        'comment': 'No quality control or scientific processing has been applied',
        'creator_name': 'LIENSs',
        'creator_email': 'valerie.ballu@univ-lr.fr',
        'creator_url': 'https://lienss.univ-larochelle.fr/',
        'processing_level': 'Raw Data',
        'summary': 'Ocean Bottom Pressure and Temperature Data from SBE53 pressure sensor'
    }

    ### step 3.6: Export to NetCDF

    # Output directory
    output_dir = "/Users/rjoyar01/Desktop/mayotte/A0A_MAY30_R/test_netcdf"

    output_path = export_obp.export_obp_to_netcdf(
        df_obp=df,
        column_mapping=column_mapping,
        output_dir=output_dir,
        metadata_dict=metadata_dict,
        conversion_factors=conversion_factors,
        force=True
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