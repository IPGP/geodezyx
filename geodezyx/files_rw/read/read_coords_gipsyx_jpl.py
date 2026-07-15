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
import gzip

#### Import the logger
import logging
import os
import re

import numpy as np
import pandas as pd

#### geodeZYX modules
from geodezyx import conv
from geodezyx import time_series

log = logging.getLogger("geodezyx")


##########  END IMPORT  ##########
#       _ _____  _          _______ _____ _____   _______     __  ______ _ _
#      | |  __ \| |        / / ____|_   _|  __ \ / ____\ \   / / |  ____(_) |
#      | | |__) | |       / / |  __  | | | |__) | (___  \ \_/ /  | |__   _| | ___  ___
#  _   | |  ___/| |      / /| | |_ | | | |  ___/ \___ \  \   /   |  __| | | |/ _ \/ __|
# | |__| | |    | |____ / / | |__| |_| |_| |     ____) |  | |    | |    | | |  __/\__ \
#  \____/|_|    |______/_/   \_____|_____|_|    |_____/   |_|    |_|    |_|_|\___||___/


def read_gipsy_tdp(filein):
    """
    Read legacy Gipsy TDP (Time Dependent Parameter) file.

    Parameters
    ----------
    filein : str
        Path to the input TDP file.

    Returns
    -------
    TimeSeriePoint
        Time series object containing the extracted coordinates.

    """

    X, Y, Z = np.nan, np.nan, np.nan
    Tx, Ty, Tz, T = np.nan, np.nan, np.nan, np.nan
    sX, sY, sZ = np.nan, np.nan, np.nan

    tsout = time_series.TimeSeriePoint()

    for line in open(filein):

        fields = line.split()

        if fields[4] == "STA" and fields[5] == "Z":
            Tz = conv.tgipsy2dt(fields[0])
            Z = float(fields[2]) * 1000
            sZ = float(fields[3]) * 1000

        if fields[4] == "STA" and fields[5] == "Y":
            Ty = conv.tgipsy2dt(fields[0])
            Y = float(fields[2]) * 1000
            sY = float(fields[3]) * 1000

        if fields[4] == "STA" and fields[5] == "X":
            Tx = conv.tgipsy2dt(fields[0])
            X = float(fields[2]) * 1000
            sX = float(fields[3]) * 1000
            STAT = fields[6]

        if Tx == Ty == Tz:
            T = Tx
            point = time_series.Point(X, Y, Z, T, "XYZ", sX, sY, sZ)
            tsout.add_point(point)

            Tx = 111
            Ty = 222
            Tz = 333

    tsout.meta_set(filein, stat=STAT)

    return tsout


def read_gipsy_tdp_list(filelistin):
    """
    Read multiple legacy Gipsy TDP (Time Dependent Parameter) files.

    Parameters
    ----------
    filelistin : list
        List of input file paths.

    Returns
    -------
    TimeSeriePoint
        Merged time series object containing data from all input files.

    """

    tslist = []
    for fil in filelistin:
        ts = read_gipsy_tdp(fil)
        tslist.append(ts)

    tsout = time_series.merge_ts(tslist)

    stat = list(set([ts.stat for ts in tslist]))[0]

    tsout.meta_set("", stat)

    return tsout


def read_gipsyx_tdp(filein):
    """
    Read GipsyX TDP (Time Dependent Parameter) file.

    Parameters
    ----------
    filein : str
        Path to the input GipsyX TDP file.

    Returns
    -------
    TimeSeriePoint
        Time series object containing the extracted coordinates.

    Raises
    ------
    Exception
        If the file cannot be opened.

    """

    X, Y, Z = np.nan, np.nan, np.nan
    tx, ty, tz, T = np.nan, np.nan, np.nan, np.nan
    s_x, s_y, s_z = np.nan, np.nan, np.nan

    tsout = time_series.TimeSeriePoint()

    try:
        fil = open(filein)
    except Exception as e:
        log.error("unable to open %s", filein)
        raise e

    for line in fil:

        fields = line.split()

        if len(fields) == 0:
            continue

        attribs = fields[-1].split(".")

        if attribs[1] == "Station" and attribs[5] == "Z":
            tz = conv.tgipsy2dt(fields[0])
            Z = float(fields[2])
            s_z = float(fields[3])

        if attribs[1] == "Station" and attribs[5] == "Y":
            ty = conv.tgipsy2dt(fields[0])
            Y = float(fields[2])
            s_y = float(fields[3])

        if attribs[1] == "Station" and attribs[5] == "X":
            tx = conv.tgipsy2dt(fields[0])
            X = float(fields[2])
            s_x = float(fields[3])
            stat = attribs[2]

        if tx == ty == tz:
            T = tx
            point = time_series.Point(X, Y, Z, T, "XYZ", s_x, s_y, s_z)
            tsout.add_point(point)

            tx = np.nan
            ty = np.nan
            tz = np.nan

    tsout.meta_set(filein, stat=stat)
    tsout.sort()

    return tsout


def read_gipsyx_tdp_list(filelistin):
    """
    Read multiple GipsyX TDP (Time Dependent Parameter) files.

    Parameters
    ----------
    filelistin : list
        List of input file paths.

    Returns
    -------
    TimeSeriePoint
        Merged time series object containing data from all input files.

    """

    tslist = []
    for fil in filelistin:
        ts = read_gipsyx_tdp(fil)
        tslist.append(ts)

    tsout = time_series.merge_ts(tslist)

    stat = list(set([ts.stat for ts in tslist]))[0]

    tsout.meta_set("", stat)

    return tsout


def read_gipsy_gdcov(filein):
    """
    Read GIPSY GDCOV (geodetic covariance) file.

    Parameters
    ----------
    filein : str
        Path to the input GDCOV file.

    Returns
    -------
    TimeSeriePoint
        Time series object containing the extracted coordinates with covariance information.

    """
    X, Y, Z = np.nan, np.nan, np.nan
    Tx, Ty, Tz, T = np.nan, np.nan, np.nan, np.nan
    sX, sY, sZ = np.nan, np.nan, np.nan

    tsout = time_series.TimeSeriePoint()

    F = open(filein)
    L = F.readlines()

    ### parameters search
    param = int(L[0].split()[0])

    Lparam = L[1 : param + 1]
    Lcovar = L[param + 2 :]

    for line in Lparam:
        fields = line.split()
        attribs = fields[1].split(".")

        if attribs[1] == "STA" and attribs[-1] == "Z":
            Tz = conv.tgipsy2dt(fields[2])
            Z = float(fields[3])
            sZ = float(fields[4])

        if attribs[1] == "STA" and attribs[-1] == "Y":
            Ty = conv.tgipsy2dt(fields[2])
            Y = float(fields[3])
            sY = float(fields[4])

        if attribs[1] == "STA" and attribs[-1] == "X":
            Tx = conv.tgipsy2dt(fields[2])
            X = float(fields[3])
            sX = float(fields[4])
            STAT = attribs[0]

        if Tx == Ty == Tz:
            T = Tx
            point = time_series.Point(X, Y, Z, T, "XYZ", sX, sY, sZ)
            tsout.add_point(point)

            Tx = np.nan
            Ty = np.nan
            Tz = np.nan

    tsout.meta_set(filein, stat=STAT)

    return tsout


def read_gipsy_gdcov_list(filelistin):
    """
    Read multiple GIPSY GDCOV (geodetic covariance) files.

    Parameters
    ----------
    filelistin : list
        List of input file paths.

    Returns
    -------
    TimeSeriePoint
        Merged time series object containing data from all input files.

    """
    tslist = []
    for fil in filelistin:
        ts = read_gipsy_gdcov(fil)
        tslist.append(ts)

    tsout = time_series.merge_ts(tslist)

    stat = list(set([ts.stat for ts in tslist]))[0]

    tsout.meta_set("", stat)

    return tsout


def read_gipsyx_xfile(filein):
    """
    Read GipsyX X file containing transformation parameters and residuals.

    Parameters
    ----------
    filein : str
        Input file path. Can handle gzip-compressed files (.gz).

    Returns
    -------
    df_trans_out : DataFrame
        Helmert transformation parameters and their standard deviations.
        Columns include transformation parameters (e.g., TX, TY, TZ) and
        their corresponding standard deviations (sTX, sTY, sTZ).
    df_resid_out : DataFrame
        Coordinates residuals (not implemented yet).

    """

    fname = os.path.basename(filein)

    date = conv.date_string_2_dt(fname[:11])

    if filein[-2:] in ("gz", "GZ"):
        F = gzip.open(filein, "r+")
        lines = [e.decode("utf-8") for e in F]
    else:
        F = open(filein, "r+")
        lines = F.readlines()

    df_trans_out = pd.DataFrame()
    df_resid_out = pd.DataFrame()

    df_trans_out.loc[0, "epoch"] = date

    for l in lines:
        #### get transform parameters
        if re.search(" = ", l):
            l2 = l.split()
            label = l2[0]
            val = float(l2[2])
            df_trans_out.loc[0, label] = val

            if len(l2) > 3:
                val_sigma = float(l2[4])
                df_trans_out.loc[0, "s" + label] = val_sigma

        #### get residual values
        # l_resid = []
        # if re.search('^ ( POS| RES)', l):
        #     l_resid.append(l)

        # df_resid_out = pd.read_csv(StringIO("\n".join(l_resid[:-1])))

    return df_trans_out, df_resid_out


def read_gipsyx_xfile_list(filelistin):
    """
    Read multiple GipsyX X files containing transformation parameters and residuals.

    Parameters
    ----------
    filelistin : list
        List of input file paths. Can handle gzip-compressed files (.gz).

    Returns
    -------
    df_trans_out : DataFrame
        Combined Helmert transformation parameters and their standard deviations
        from all input files. Columns include transformation parameters and their
        corresponding standard deviations.
    df_resid_out : DataFrame
        Coordinates residuals (not implemented yet).

    """
    dflist = []
    for fil in filelistin:
        df_trans_mono, df_resid_mono = read_gipsyx_xfile(fil)
        dflist.append(df_trans_mono)

    df_trans_out = pd.concat(dflist)
    df_trans_out.reset_index(drop=True, inplace=True)
    df_resid_out = pd.DataFrame()

    return df_trans_out, df_resid_out


def read_gipsy_bosser(filein):
    """
    Read a GIPSY file in the format used by P. Bosser at ENSTA Brest.

    Parameters
    ----------
    filein : str
        Path to the input file.

    Returns
    -------
    TimeSeriePoint
        Time series object containing the extracted coordinates in FLH format.

    Notes
    -----
    File format contains year, day of year, latitude, longitude, height,
    and RMS values. Latitude and longitude are converted from radians to degrees.

    """

    F, L, H = 0, 0, 0
    T = 0
    sF, sL, sH = 0, 0, 0

    tsout = time_series.TimeSeriePoint()

    for line in open(filein):
        f = line.split()

        y = int(f[0])
        doydec = float(f[1])
        doy = int(doydec)

        T = conv.doy2dt(y, doy) + dt.timedelta(days=(doydec - doy))

        F = np.rad2deg(float(f[2]))
        L = np.rad2deg(float(f[3]))
        H = float(f[4])
        RMS = float(f[8])

        point = time_series.Point(F, L, H, T, "FLH")
        point.anex["RMS"] = RMS
        tsout.add_point(point)

    tsout.meta_set(filein, stat="STAT")

    return tsout


def read_gipsy_apps(filein):
    """
    Read GIPSY APPS (Online processing tool) output file.

    Parameters
    ----------
    filein : str
        Path to the input file.

    Returns
    -------
    TimeSeriePoint
        Time series object containing the extracted coordinates in XYZ format
        with standard deviations.

    Notes
    -----
    The file format contains timestamp, X, Y, Z coordinates and their corresponding
    standard deviations. Comment lines starting with '#' and lines containing
    'Kinematic Processing' are skipped.

    """

    tsout = time_series.TimeSeriePoint()

    for l in open(filein):
        if l[0] == "#" or ("Kinematic Processing" in l):
            continue
        f = l.split()

        date_lis = [int(float(e)) for e in f[1].split(":")]

        T = dt.datetime(*date_lis)

        X = float(f[2])
        sX = float(f[3])
        Y = float(f[4])
        sY = float(f[5])
        Z = float(f[6])
        sZ = float(f[7])

        point = time_series.Point(X, Y, Z, T, "XYZ", sX, sY, sZ)
        tsout.add_point(point)

    tsout.meta_set(filein)
    return tsout


def read_jpl_timeseries_solo(latlonrad_files_list):
    """
    Read JPL time series from separate latitude, longitude, and radial component files.

    Parameters
    ----------
    latlonrad_files_list : list
        List of three file paths containing latitude (.lat), longitude (.lon),
        and radial (.rad) components of the time series.

    Returns
    -------
    TimeSeriePoint
        Time series object containing the extracted coordinates in ENU (East-North-Up)
        format with standard deviations. Station name is extracted from the file names.

    Raises
    ------
    Exception
        If the files have different epoch times or if station names do not match
        across the three input files.

    Notes
    -----
    The three input files are expected to contain the same epochs in the same order.
    Station name is extracted from the filename prefix (before the first dot).
    Coordinates are assumed to be in centimeter units and converted to meter units.

    """
    tsout = time_series.TimeSeriePoint()

    latpath = [f for f in latlonrad_files_list if ".lat" in f][0]
    lonpath = [f for f in latlonrad_files_list if ".lon" in f][0]
    radpath = [f for f in latlonrad_files_list if ".rad" in f][0]

    latfile = open(latpath)
    lonfile = open(lonpath)
    radfile = open(radpath)

    for llat, llon, lrad in zip(latfile, lonfile, radfile):
        flat = [float(e) for e in llat.split()[0:3]]
        flon = [float(e) for e in llon.split()[0:3]]
        frad = [float(e) for e in lrad.split()[0:3]]

        if not (flat[0] == flon[0] == frad[0]):
            log.error("%s %s %s", flat[0], flon[0], frad[0])
            log.error("Time dont corresponds !!!")
            raise Exception

        statlat = os.path.basename(latpath).split(".")[0]
        statlon = os.path.basename(lonpath).split(".")[0]
        statrad = os.path.basename(radpath).split(".")[0]

        if not (statlat == statlon == statrad):
            log.info("%s %s %s", statlat, statlon, statrad)
            log.error("Station name do not corresponds !!!")
            raise Exception

        T = conv.year_decimal2dt(flat[0])

        N = flat[1] * 10**-2
        E = flon[1] * 10**-2
        U = frad[1] * 10**-2

        sN = flat[2] * 10**-2
        sE = flon[2] * 10**-2
        sU = frad[2] * 10**-2

        point = time_series.Point(E, N, U, T, "ENU", sE, sN, sU)

        tsout.boolENU = True
        tsout.add_point(point)

    tsout.stat = statlat

    return tsout

