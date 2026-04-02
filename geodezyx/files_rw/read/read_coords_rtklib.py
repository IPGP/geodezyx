#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.files_rw contains reading functions to
import files containing geodetic time series.

it can be imported directly with:
from geodezyx import files_rw

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""

########## BEGIN IMPORT ##########
#### External modules
import datetime as dt
import glob
import gzip

#### Import the logger
import logging
import os
import re

import numpy as np
import pandas as pd
import scipy

#### geodeZYX modules
from geodezyx import conv
from geodezyx import files_rw
from geodezyx import time_series
from geodezyx import utils

log = logging.getLogger("geodezyx")


##########  END IMPORT  ##########


def _rtklib_header(filein):
    """Parse RTKLIB header to extract metadata.

    Returns coordinate type (initype) and UTC-GPS offset (d_utcgps).
    """
    initype = "FLH"
    d_utcgps = 0

    re_e_baseline = re.compile(f"e-baseline")
    re_x_ecef = re.compile(f"x-ecef")
    re_utc = re.compile(f"%  UTC")
    re_gpst = re.compile(f"%  GPST")
    re_inp_file =re.compile(f"% inp file")

    inp_file_lis = []

    with open(filein) as f:
        for line in f:
            if line[0] != "%":
                break

            if re_e_baseline.search(line):
                initype = "ENU"
            elif re_x_ecef.search(line):
                initype = "XYZ"
            elif "latitude(deg)" in line:
                initype = "FLH"

            if re_utc.search(line):
                d_utcgps = 17
                log.warning("HARDCODED (MAYBE?) WRONG LEAP SECOND!!!!")
            elif re_gpst.search(line):
                d_utcgps = 0

            if re_inp_file.search(line):
                inp_file_lis.append(line)

    return initype, d_utcgps, inp_file_lis


def _rtklib_load_df(filein, d_utcgps=0, inpfilis=None):
    """Load RTKLIB file into DataFrame."""
    try:
        df = pd.read_csv(filein, comment="%", sep=r"\s+", header=None, engine="python")
    except Exception as e:
        log.error("Error reading RTKLIB file with pandas: %s", e)
        df = pd.DataFrame()  # Return empty DataFrame on error

    if df.empty:
        log.warning("File %s is empty", os.path.basename(filein))
        return df

    # Assign column names
    df.columns = [
        "date",
        "time",
        "a",
        "b",
        "c",
        "Q",
        "ns",
        "s_a",
        "s_b",
        "s_c",
        "sdAB",
        "sdBC",
        "sdAC",
        "age",
        "ratio",
    ]
    df["date"] = (
        pd.to_datetime(df["date"])
        + pd.to_timedelta(df["time"])
        + pd.Timedelta(seconds=d_utcgps)
    )

    # get rover and base names from header if possible, else set to "XXXX00XXX"
    if inpfilis is not None:
        rover, base = _rtklib_inpfiles2sites(inpfilis)
    else:
        rover, base = "XXXX00XXX", "XXXX00XXX"

    df["rover"] = rover
    df["base"] = base

    # remove time columns
    df.drop(["time"], axis=1, inplace=True)
    # reorder column 'date' to 'epoch'
    df.rename(columns={"date": "epoch"}, inplace=True)

    return df


# def _rtklib_parse_dt(date_str, time_str, d_utcgps):
#     """Parse date and time strings from RTKLIB file."""
#     re_date = re.compile(r"[\w']+")
#     date_full = date_str + ":" + time_str
#     date_parts_matched = re_date.findall(date_full)
#
#     date_parts = [int(d) for d in date_parts_matched[:-1]]
#     date_parts.append(int(date_parts_matched[-1][:6]))
#
#     t = dt.datetime(
#         date_parts[0],
#         date_parts[1],
#         date_parts[2],
#         date_parts[3],
#         date_parts[4],
#         date_parts[5],
#         date_parts[6],
#     ) + dt.timedelta(seconds=d_utcgps)
#     return t


def _rtklib_ts_point(row, initype, t):
    """Create a Point object from a DataFrame row."""
    a = float(row["a"])
    b = float(row["b"])
    c = float(row["c"])

    s_a = float(row["s_a"])
    s_b = float(row["s_b"])
    s_c = float(row["s_c"])

    if initype == "FLH":
        s_a, s_b, s_c = conv.sigma_enu2geo(a, b, c, s_a, s_b, s_c)

    point = time_series.Point(a, b, c, t, initype, s_a, s_b, s_c)
    point.anex["sdAB"] = float(row["sdAB"])
    point.anex["sdBC"] = float(row["sdBC"])
    point.anex["sdAC"] = float(row["sdAC"])
    point.anex["Q"] = float(row["Q"])

    return point

def _rtklib_inpfiles2sites(inpfilis):
    rov_bn = os.path.basename(inpfilis[0].split()[-1])
    bas_bn = os.path.basename(inpfilis[1].split()[-1])

    if conv.rinex_regex_search_tester(rov_bn,short_name=False):
        rover = rov_bn[0:4].upper()
    else:
        rover = rov_bn[0:4].upper()

    if conv.rinex_regex_search_tester(bas_bn, short_name=False):
        base = bas_bn[0:4].upper()
    else:
        base = bas_bn[0:4].upper()
    return rover, base

def _rtklib_ts_meta(tsout, filein, initype, inpfilis):
    """Add metadata to TimeSeriePoint from RTKLIB file."""
    tsout.meta_set(filein)

    if initype == "ENU":
        tsout.boolENU = True

    try:
        rov, bas = _rtklib_inpfiles2sites(inpfilis)
        tsout.anex["rover"] = rov
        tsout.anex["base"] = bas
    except:
        pass


def _rtklib_df_col_names(df, initype):
    """Rename DataFrame columns based on coordinate type."""

    col0 = ["epoch"]
    col2 = ["Q", "ns"]
    col4 = ["age", "ratio", "rover", "base"]

    if initype == "ENU":
        col1 = ["e", "n", "u"]
        col3 = ["sde", "sdn", "sdu", "sden", "sdnu", "sdeu"]
    elif initype == "XYZ":
        col1 = ["x", "y", "z"]
        col3 = ["sdx", "sdy", "sdz", "sdxy", "sdyz", "sdxz"]
    elif initype == "FLH":
        col1 = ["lat", "lon", "h"]
        col3 = ["sdlat", "sdlon", "sdh", "sdlatlon", "sdlonh", "sdlath"]
    else:
        log.error("Unknown coordinate type: %s", initype)
        col1 = ["a", "b", "c"]
        col3 = ["s_a", "s_b", "s_c", "sdAB", "sdBC", "sdAC"]

    col = col0 + col1 + col2 + col3 + col4
    df.columns = col
    return df


def read_rtklib(filein, return_df=False):
    """Read a RTKLIB file.

    Parameters
    ----------
    filein : str
        Input file path
    return_df : bool, optional
        If True, return DataFrame instead of TimeSeriePoint (default: False)

    Returns
    -------
    TimeSeriePoint or DataFrame
        TimeSeriePoint object or pandas DataFrame depending on return_df
    """
    # Step 1: Parse header for metadata
    initype, d_utcgps, inp_file_lis = _rtklib_header(filein)

    # Step 2: Load data into DataFrame
    df = _rtklib_load_df(filein, d_utcgps, inp_file_lis)

    # Handle empty file
    if df.empty:
        if return_df:
            return df
        else:
            tsout = time_series.TimeSeriePoint()
            tsout.meta_set(filein)
            return tsout

    # Step 3: If return_df requested, rename columns and return DataFrame
    if return_df:
        df = _rtklib_df_col_names(df, initype)
        return df

    # Step 4: Convert DataFrame to TimeSeriePoint
    tsout = time_series.TimeSeriePoint()

    for idx, row in df.iterrows():
        # Parse datetime
        t = conv.pandas_timestamp2dt(row["epoch"])

        # Create and add point
        point = _rtklib_ts_point(row, initype, t)
        tsout.add_point(point)

    # Step 5: Add metadata
    _rtklib_ts_meta(tsout, filein, initype, inp_file_lis)

    return tsout


def _read_rtklib_legacy(filein):
    """
    Legacy line-by-line implementation of read_rtklib.
    Used as fallback if pandas parsing fails.
    """
    tsout = time_series.TimeSeriePoint()

    initype = "FLH"
    d_utcgps = 0

    # Precompile regex patterns outside the loop for better performance
    re_e_baseline = re.compile("e-baseline")
    re_x_ecef = re.compile("x-ecef")
    re_utc = re.compile("%  UTC")
    re_gpst = re.compile("%  GPST")
    re_date = re.compile(r"[\w']+")

    with open(filein) as f:
        for line in f:
            # Check header patterns to determine coordinate type
            if re_e_baseline.search(line):
                initype = "ENU"
            elif re_x_ecef.search(line):
                initype = "XYZ"
            elif "latitude(deg)" in line:
                initype = "FLH"

            if re_utc.search(line):
                d_utcgps = 16
                log.warning("HARDCODED (MAYBE?) WRONG LEAP SECOND!!!!")
            elif re_gpst.search(line):
                d_utcgps = 0

            # Skip comment lines
            if line[0] == "%":
                continue

            fields = line.split()

            # Optimized date parsing - avoid repeated string concatenation
            date_str = fields[0] + ":" + fields[1]
            date1 = re_date.findall(date_str)
            date_parts = [int(d) for d in date1[:-1]]
            date_parts.append(int(date1[-1][:6]))

            t = dt.datetime(
                date_parts[0],
                date_parts[1],
                date_parts[2],
                date_parts[3],
                date_parts[4],
                date_parts[5],
                date_parts[6],
            ) + dt.timedelta(seconds=d_utcgps)

            # Parse coordinates
            a = float(fields[2])
            b = float(fields[3])
            c = float(fields[4])

            # Parse standard deviations - same fields for both XYZ and FLH
            s_a = float(fields[7])
            s_b = float(fields[8])
            s_c = float(fields[9])

            if initype == "FLH":
                s_a, s_b, s_c = conv.sigma_enu2geo(a, b, c, s_a, s_b, s_c)

            point = time_series.Point(a, b, c, t, initype, s_a, s_b, s_c)
            point.anex["sdAB"] = float(fields[10])
            point.anex["sdBC"] = float(fields[11])
            point.anex["sdAC"] = float(fields[12])
            point.anex["Q"] = float(fields[5])

            tsout.add_point(point)

    tsout.meta_set(filein)

    if initype == "ENU":
        tsout.boolENU = True

    try:
        inpfilis = utils.grep(filein, "inp file")
        tsout.anex["rover"] = os.path.basename(inpfilis[0].split()[-1])[0:4].upper()
        tsout.anex["base"] = os.path.basename(inpfilis[1].split()[-1])[0:4].upper()
    except:
        pass

    return tsout

