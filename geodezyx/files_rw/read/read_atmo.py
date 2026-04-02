#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 10:52:09 2021

@author: psakic
"""

import re
import os

########## BEGIN IMPORT ##########
#### External modules
import numpy as np
import pandas as pd

#### geodeZYX modules
from geodezyx import conv
from geodezyx import utils

##########  END IMPORT  ##########


############### Reading Tropospheric ################################################

#  _______                              _                     ______ _ _
# |__   __|                            | |                   |  ____(_) |
#    | |_ __ ___  _ __   ___  ___ _ __ | |__   ___ _ __ ___  | |__   _| | ___  ___
#    | | '__/ _ \| '_ \ / _ \/ __| '_ \| '_ \ / _ \ '__/ _ \ |  __| | | |/ _ \/ __|
#    | | | | (_) | |_) | (_) \__ \ |_) | | | |  __/ | |  __/ | |    | | |  __/\__ \
#    |_|_|  \___/| .__/ \___/|___/ .__/|_| |_|\___|_|  \___| |_|    |_|_|\___||___/
#                | |             | |
#                |_|             |_|


def read_snx_trop(snxfile, dataframe_output=True, version=2):
    """
    Read troposphere solutions from Troposphere SINEX file.

    Parses SINEX format tropospheric solutions and returns station, epoch,
    and troposphere parameters (ZTD, horizontal gradients, and their uncertainties).

    Parameters
    ----------
    snxfile : str
        Path to SINEX troposphere solution file.
    dataframe_output : bool, optional
        If True, return results as pandas DataFrame. If False, return tuple.
        Default is True.
    version : int, optional
        SINEX version (1 or 2). Affects date parsing. Default is 2.

    Returns
    -------
    df : pandas.DataFrame or tuple
        If dataframe_output=True: DataFrame with columns STAT, epoc, tro, stro,
        tgn, stgn, tge, stge. If False: tuple with sorted Station, Epoch,
        ZTD, ZTD_sigma, TGN, TGN_sigma, TGE, TGE_sigma.
    """

    STAT, epoc = [], []
    tro, stro, tgn, stgn, tge, stge = [], [], [], [], [], []

    flagtrop = False

    for line in open(snxfile, "r", encoding="ISO-8859-1"):

        if re.compile("TROP/SOLUTION").search(line):
            flagtrop = not flagtrop
            continue

        if line[0] == " ":
            fields = line.split()
        else:
            continue

        if flagtrop:

            STAT.append(fields[0].upper())

            if not ":" in fields[1]:
                epoc.append(conv.year_decimal2dt(fields[1]))
            else:
                date_elts_lis = fields[1].split(":")
                if version == 2:
                    yy = int(date_elts_lis[0])
                else:
                    yy = int(date_elts_lis[0]) + 2000

                doy = int(date_elts_lis[1])
                sec = int(date_elts_lis[2])
                epoc.append(conv.doy2dt(yy, doy, seconds=sec))

            if len(fields) == 8:
                tro.append(np.nan if "*" in fields[2] else fields[2])
                stro.append(np.nan if "*" in fields[3] else fields[3])
                tgn.append(np.nan if "*" in fields[4] else fields[4])
                stgn.append(np.nan if "*" in fields[5] else fields[5])
                tge.append(np.nan if "*" in fields[6] else fields[6])
                stge.append(np.nan if "*" in fields[7] else fields[7])

            elif len(fields) == 4:
                tro.append(np.nan if "*" in fields[2] else fields[2])
                stro.append(np.nan if "*" in fields[3] else fields[3])
                tgn.append(np.nan)
                stgn.append(np.nan)
                tge.append(np.nan)
                stge.append(np.nan)

            else:
                tro.append(np.nan)
                stro.append(np.nan)
                tgn.append(np.nan)
                stgn.append(np.nan)
                tge.append(np.nan)
                stge.append(np.nan)

    outtuple = tuple(zip(*sorted(zip(STAT, epoc, tro, stro, tgn, stgn, tge, stge))))

    if dataframe_output:
        return troposinex2df(outtuple)
    else:
        return outtuple


def read_gfz_trop(trpfile):
    """
    Read GFZ troposphere SINEX solution into pandas DataFrame.

    Parses GFZ SINEX format and extracts station name, epoch, and ZTD with
    associated uncertainties plus horizontal gradients.

    Parameters
    ----------
    trpfile : str
        Path to GFZ troposphere SINEX file.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame with columns: STAT, epoc, year, doy, secofday, ztd_est,
        ztd_est_std, num_sat, tgn_est, tgn_est_std, tge_est, tge_est_std.
        Numeric columns are converted to float type.
    """
    fields = []
    flagtrop = False

    for line in open(trpfile, "r", encoding="ISO-8859-1"):
        if re.compile("TROP/SOLUTION").search(line):
            flagtrop = not flagtrop
            continue

        if flagtrop == True and line[0] == " ":
            field = line.split()
            fields.append(field)
        else:
            continue

    DF = pd.DataFrame(fields)
    DF.drop(columns=[0, 8, 9, 10, 11, 12, 18, 19, 20], inplace=True)
    DF.columns = [
        "STAT",
        "epoc",
        "year",
        "doy",
        "secofday",
        "ztd_est",
        "ztd_est_std",
        "num_sat",
        "tgn_est",
        "tgn_est_std",
        "tge_est",
        "tge_est_std",
    ]
    cols_numeric = [
        "epoc",
        "ztd_est",
        "ztd_est_std",
        "num_sat",
        "tgn_est",
        "tgn_est_std",
        "tge_est",
        "tge_est_std",
    ]
    DF[cols_numeric] = DF[cols_numeric].apply(pd.to_numeric, errors="coerce")
    DF["epoc"] = conv.mjd2dt(DF["epoc"].values)
    DF["epoc"] = DF["epoc"].dt.floor("H")
    return DF


def troposinex2df(read_sinex_result):
    """
    Convert SINEX troposphere data to pandas DataFrame.

    Transforms output from read_snx_trop function into a structured DataFrame
    with proper column names and numeric type conversion.

    Parameters
    ----------
    read_sinex_result : tuple
        Tuple of lists from read_snx_trop function containing station names,
        epochs, and troposphere parameters (ZTD, gradients, uncertainties).

    Returns
    -------
    df_sinex : pandas.DataFrame
        DataFrame with columns: STAT, epoc (datetime), tro, stro, tgn, stgn,
        tge, stge. Numeric columns are converted to float type.
    """
    DF_Sinex = pd.DataFrame.from_records(list(read_sinex_result)).transpose()
    colnam = ["STAT", "epoc", "tro", "stro", "tgn", "stgn", "tge", "stge"]
    DF_Sinex.columns = colnam
    cols_numeric = ["tro", "stro", "tgn", "stgn", "tge", "stge"]
    DF_Sinex[cols_numeric] = DF_Sinex[cols_numeric].apply(
        pd.to_numeric, errors="coerce"
    )

    return DF_Sinex


def read_bernese_trp(trpfile):
    """
    Read tropospheric solution in TRP format from Bernese GNSS software.

    Parses Bernese TRP format and returns troposphere parameters including
    ZTD, north-south and east-west gradients with associated uncertainties.

    Parameters
    ----------
    trpfile : str
        Path to TRP file from Bernese GNSS software.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame with columns: STAT, year, month, day, hour, minute, second,
        MOD_U, CORR_U, SIGMA_U, TOTAL_U, CORR_N, SIGMA_N, CORR_E, SIGMA_E,
        and datetime column 'dt'. Numeric columns are converted to float type.

    Notes
    -----
    Written by Chaiyaporn Kitpracha.
    """
    flagtrop = False
    field = []
    for line in open(trpfile, "r", encoding="ISO-8859-1"):
        if re.compile("STATION NAME").search(line):
            headers = line.split()
            headers.remove("YYYY")
            headers.remove("MM")
            headers.remove("DD")
            headers.remove("HH")
            headers.remove("MM")
            headers.remove("SS")
            headers[3] = "year"
            headers[4] = "month"
            headers[5] = "day"
            headers[6] = "hour"
            headers[7] = "minute"
            headers[8] = "second"
            flagtrop = True
            continue

        if flagtrop and not line == "\n":
            fields = line.split()
            field.append(fields)
        else:
            continue

    DF = pd.DataFrame(field, columns=headers)
    DF["dt"] = pd.to_datetime(DF[["year", "month", "day", "hour", "minute", "second"]])
    cols_num = [
        "MOD_U",
        "CORR_U",
        "SIGMA_U",
        "TOTAL_U",
        "CORR_N",
        "SIGMA_N",
        "CORR_E",
        "SIGMA_E",
    ]
    DF[cols_num] = DF[cols_num].apply(pd.to_numeric, errors="coerce")
    DF.drop(["year", "month", "day", "hour", "minute", "second"], axis=1, inplace=True)
    return DF


def read_rinex_met(metfile):
    """
    Read RINEX meteorological files and convert to pandas DataFrame.

    Handles single file or multiple files (as list or iterable). Concatenates
    multiple files into a single DataFrame with proper time indexing.

    Parameters
    ----------
    metfile : str or list of str
        Path(s) to RINEX meteorological file(s). Can be a single filename string
        or a list/iterable of filenames (e.g., from glob).

    Returns
    -------
    df : pandas.DataFrame
        Meteorological data from RINEX file(s). Index is set to epoch (datetime).
        Includes columns for temperature, pressure, humidity, and associated
        uncertainties if available.

    Notes
    -----
    Written by Chaiyaporn Kitpracha.
    """
    if utils.is_iterable(metfile):
        merge_df = pd.DataFrame()
        for metfile_m in metfile:
            met_df = read_rinex_met_2(str(metfile_m))
            merge_df = pd.concat([merge_df, met_df])
        return merge_df
    else:
        met_df = read_rinex_met_2(metfile)
        return met_df


def read_rinex_met_2(metfile):
    """
    Read a single RINEX meteorological file (internal function).

    Worker function for read_rinex_met. Parses header for sensor metadata and
    converts meteorological observations to DataFrame with proper types and
    datetime indexing.

    Parameters
    ----------
    metfile : str
        Path to a single RINEX meteorological file.

    Returns
    -------
    df : pandas.DataFrame
        Meteorological data with datetime index. Includes temperature, pressure,
        humidity observations and sensor uncertainties. Station code (STA column)
        is added from RINEX header.

    Notes
    -----
    Written by Chaiyaporn Kitpracha.
    """
    ln = 0
    for line in open(metfile, "r", encoding="ISO-8859-1"):
        if re.compile("MARKER NAME").search(line):
            marker = line.split()[0]
            marker = marker.upper()
        if re.compile("# / TYPES OF OBSERV").search(line):
            tmp = line.split()
            headers = tmp[1 : int(tmp[0]) + 1]
        if re.compile("TD SENSOR MOD/TYPE/ACC").search(line):
            tmp = line.split()
            temp_unc = float(tmp[-4])
        if re.compile("PR SENSOR MOD/TYPE/ACC").search(line):
            tmp = line.split()
            press_unc = float(tmp[-4])
        if re.compile("HR SENSOR MOD/TYPE/ACC").search(line):
            tmp = line.split()
            humrel_unc = float(tmp[-4])
        if re.compile("END OF HEADER").search(line):
            break
        ln = ln + 1

    df = pd.read_csv(
        metfile,
        skiprows=range(0, ln + 1),
        delim_whitespace=True,
        names=["year", "month", "day", "hour", "minute", "second"] + headers,
    )
    df["year"] = (
        df["year"] + 2000 if df["year"].any() <= 79 else df["year"].any() + 1900
    )
    df["STA"] = marker
    df["epoch"] = pd.to_datetime(
        df[["year", "month", "day", "hour", "minute", "second"]], errors="coerce"
    )
    df.drop(["year", "month", "day", "hour", "minute", "second"], axis=1, inplace=True)
    if press_unc is not None:
        df["PR_std"] = press_unc
    if temp_unc is not None:
        df["TD_std"] = temp_unc
    if humrel_unc is not None:
        df["HR_std"] = humrel_unc
    df.set_index("epoch", inplace=True)
    return df


def read_spotgins_tropo(filepath):
    """
    Read a SPOTGINS tropospheric time series file (ZTD or GRAD format).

    Parses SPOTGINS ZTD (zenith total delay) or GRAD (horizontal gradients)
    format and returns tropospheric parameters with associated metadata.
    Column names are extracted from file headers and converted to lowercase.
    Header metadata is dynamically inferred from the file (no hardcoded fields).

    Parameters
    ----------
    filepath : str
        Path to a SPOTGINS ``.ztd`` or ``.grad`` file (can be gzip compressed).

    Returns
    -------
    df : pandas.DataFrame
        DataFrame with one row per epoch. Columns are lowercase with the
        datetime column named 'epoch'. No index is set.
    meta : dict
        Header metadata extracted from comment block. Keys are converted to
        lowercase with spaces replaced by underscores. Includes station
        information, constellation, coordinates, and other metadata.
    """
    meta = {}
    data_lines = []
    cols = None
    datetime_col = None

    with open(filepath, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                # The last comment line is the column header
                cols = line.lstrip("#").strip().split()
                # Find the datetime column (contains timestamp format)
                datetime_col = next((col for col in cols if any(x in col.lower() for x in ["hh", "mm", "ss"])), None)
                
                # Extract metadata from key: value lines
                m = re.match(r"^#\s*([^:]+?)\s*:\s*(.+)", line)
                if m:
                    key = m.group(1).strip()
                    value = m.group(2).strip()
                    # Convert key to lowercase with underscores for spaces
                    meta_key = key.lower().replace(" ", "_")
                    try:
                        meta[meta_key] = float(value)
                    except ValueError:
                        meta[meta_key] = value
            else:
                stripped = line.strip()
                if stripped:
                    data_lines.append(stripped)

    if cols is None:
        raise ValueError("Could not find column header in file")

    if not data_lines:
        return pd.DataFrame(columns=cols).rename(str.lower, axis=1), meta

    from io import StringIO

    df = pd.read_csv(
        StringIO("\n".join(data_lines)),
        sep=r"\s+",
        header=None,
        names=cols,
    )

    # Convert to lowercase column names
    df.columns = df.columns.str.lower()

    # Rename datetime column to epoch
    if datetime_col:
        datetime_col_lower = datetime_col.lower()
        if datetime_col_lower in df.columns:
            df = df.rename(columns={datetime_col_lower: "epoch"})

    # Convert numeric columns to float (exclude string/metadata columns)
    # String columns: const, dateofexe, ginsversion, prairieversion, epoch
    string_cols = {"epoch", "const", "dateofexe", "ginsversion", "prairieversion"}
    float_cols = [c for c in df.columns if c not in string_cols]
    df[float_cols] = df[float_cols].apply(pd.to_numeric, errors="coerce")
    df["epoch"] = pd.to_datetime(df["epoch"], errors="coerce")

    # Add station from metadata as first column named 'site'
    if "station" in meta:
        df.insert(0, "site", meta["station"])

    return df, meta
