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

#   _____ ______ ______    ________ _____   ____   _____   ______ _ _
#  / ____|  ____|___  /   / /  ____|  __ \ / __ \ / ____| |  ____(_) |
# | |  __| |__     / /   / /| |__  | |__) | |  | | (___   | |__   _| | ___  ___
# | | |_ |  __|   / /   / / |  __| |  ___/| |  | |\___ \  |  __| | | |/ _ \/ __|
# | |__| | |     / /__ / /  | |____| |    | |__| |____) | | |    | | |  __/\__ \
#  \_____|_|    /_____/_/   |______|_|     \____/|_____/  |_|    |_|_|\___||___/


def read_epos_sta_kinematics(filein):
    """
    read an EPOS kinematic solutions
    """

    F = open(filein)
    Lines_4_DF_stk = []
    for l in F:
        fields = l.split()
        if l[0] != "K" and l[0] != "U" and l[0] != "X":
            continue
        if l[0] == "K" or l[0] == "U" or l[0] == "X":
            namstat = fields[2]
            numstat = int(fields[1])
            MJD_epo = float(fields[3])
            numobs = int(fields[4])

            X = float(fields[6])
            Y = float(fields[7])
            Z = float(fields[8])
            sX = float(fields[10])
            sY = float(fields[11])
            sZ = float(fields[12])

            N = float(fields[14])
            E = float(fields[15])
            U = float(fields[16])
            sN = float(fields[18])
            sE = float(fields[19])
            sU = float(fields[20])

            tup_4_df = (
                namstat,
                numstat,
                MJD_epo,
                numobs,
                X,
                Y,
                Z,
                sX,
                sY,
                sZ,
                N,
                E,
                U,
                sN,
                sE,
                sU,
            )
            Lines_4_DF_stk.append(tup_4_df)

    columns = (
        "site",
        "site_num",
        "MJD_epo",
        "numobs",
        "x",
        "y",
        "z",
        "sx",
        "sy",
        "sz",
        "N",
        "E",
        "U",
        "sN",
        "sE",
        "s_u",
    )

    DFout = pd.DataFrame(Lines_4_DF_stk, columns=columns)
    return DFout


def read_epos_sta_coords_mono(filein, return_df=True):
    """
    Read an GFZ EPOS's coordinate file.
    To read several files at the same time see read_epos_sta_coords_multi

    Parameters
    ----------
    filein : str
        path of the input coordinate file.
    return_df : bool, optional
        if True, returns the coordinates as a Pandas DataFrame.
        if False, return a list of GeodeZYX's point objects (advanced)
        The default is True.

    Returns
    -------
    DataFrame or List of Points
        Ouput coordinates.

    """
    # """
    # read an EPOS's YYYY_DDD_XX_sta_coordinates coordinates files
    # and return a list of Points objects

    # ... TBC ...
    # """

    from unlzw import unlzw

    if type(filein) is str:  ### case 1 : path compressed
        if filein[-2:] in (".Z"):
            with open(filein, "rb") as fh:
                compressed_data = fh.read()
                F = unlzw(compressed_data)

        if filein[-2:] in ("gz", "GZ"):
            F = gzip.open(filein, "r+")
            F = [e.decode("utf-8") for e in F]
        else:  ### case 2 : path uncompressed
            try:
                F = open(filein, "r", encoding="ISO-8859-1")
            except:
                F = open(filein, "r")
    else:  ### case 3 : already a list of lines
        F = open(filein)

    points_list_stk = []
    lines_4_df_stk = []

    for l in F:
        fields = l.split()
        if l[0] != " ":
            continue
        if "SITE" in fields[0]:
            namestat = fields[8]
            numstat = int(fields[2])
            tecto_plate = fields[4]
            MJD_ref = int(fields[5])
            MJD_strt = int(fields[6])
            MJD_end = int(fields[7])
            MJD_mid = np.mean([MJD_strt, MJD_end])
            T = conv.numpy_dt2dt(conv.mjd2dt(MJD_mid))

        if "POS_VEL:XYZ" in fields[0]:
            X = float(fields[4])
            Y = float(fields[5])
            Z = float(fields[6])
            Vx = float(fields[7])
            Vy = float(fields[8])
            Vz = float(fields[9])

        if "SIG_PV_XYZ" in fields[0]:
            sX = float(fields[4].replace("D", "E"))
            sY = float(fields[5].replace("D", "E"))
            sZ = float(fields[6].replace("D", "E"))
            sVx = float(fields[7])
            sVy = float(fields[8])
            sVz = float(fields[9])

            #### Last useful line for the point, store it
            if not return_df:
                point = time_series.Point(X, Y, Z, T, "XYZ", sX, sY, sZ, name=namestat)
                point.anex["Vx"] = sVx
                point.anex["Vy"] = sVy
                point.anex["Vz"] = sVz
                points_list_stk.append(point)

            #### And store for the DataFrame
            else:
                tup_4_DF = (
                    namestat,
                    numstat,
                    tecto_plate,
                    conv.mjd2dt(MJD_strt),
                    MJD_ref,
                    MJD_strt,
                    MJD_end,
                    X,
                    Y,
                    Z,
                    sX,
                    sY,
                    sZ,
                    Vx,
                    Vy,
                    Vz,
                    sVx,
                    sVy,
                    sVz,
                )

                lines_4_df_stk.append(tup_4_DF)

    if return_df:
        columns = (
            "site",
            "site_num",
            "tecto_plate",
            "epoch",
            "MJD_ref",
            "MJD_start",
            "MJD_end",
            "x",
            "y",
            "z",
            "sx",
            "sy",
            "sz",
            "Vx",
            "Vy",
            "Vz",
            "sVx",
            "sVy",
            "sVz",
        )

        DFout = pd.DataFrame(lines_4_df_stk, columns=columns)

        return DFout
    else:
        return points_list_stk


def read_epos_sta_coords_multi(filein_list, output_type="DataFrame"):
    """
    Read several GFZ EPOS's coordinate files.

    Parameters
    ----------
    filein_list : list
        list of input coordinate files inputs.
    output_type : str, optional
        "DataFrame": returns a Pandas DataFrame containing the coordinates
        "TSobjects": returns a dictionary of GeodeZYX's TimeSeries objects (advanced)
        The default is "DataFrame".

    Returns
    -------
    OUT : DataFrame or dict
        See "output_type" input parameter.

    """

    if output_type == "TSobjects":
        OUT = read_epos_sta_coords_multi_legacy(filein_list, return_dict=True)
    elif output_type == "DataFrame":
        DFfil_stk = []
        for fil in filein_list:
            DFfil = read_epos_sta_coords_mono(fil, return_df=True)
            DFfil_stk.append(DFfil)
        DFall = pd.concat(DFfil_stk)
        DFall.reset_index(inplace=True, drop=True)
        OUT = DFall

    return OUT


def read_epos_sta_coords_multi_legacy(filein_list, return_dict=True):
    """
    Read several GFZ EPOS's coordinate files.
    Legacy version

    Parameters
    ----------
    filein_list : list
        list of input coordinate files inputs.
    return_dict : bool, optional
        True: returns a dictionary of GeodeZYX's TimeSeries objects
        "TSobjects": returns a list of GeodeZYX's TimeSeries objects
        The default is True.

    Returns
    -------
    OUT : dict or list
        See "return_dict" input parameter.
    """

    filein_list = sorted(filein_list)
    Points_list = []
    statname_stk = []

    for fil in filein_list:
        Points_daily_list = read_epos_sta_coords_mono(fil, return_df=False)
        Points_list = Points_list + Points_daily_list
        statname_stk = statname_stk + [e.name for e in Points_daily_list]

    statname_uniq = sorted(list(set(statname_stk)))

    ts_dict = dict()

    for point in Points_list:
        if not point.name in ts_dict.keys():
            ts_dict[point.name] = time_series.TimeSeriePoint(stat=point.name)
        ts_dict[point.name].add_point(point)

    if return_dict:
        return ts_dict
    else:
        ts_list = []
        for k, val in ts_dict.items():
            ts_list.append(val)
        return ts_list


def read_epos_slv_times(p, convert_to_time=False):
    """
    convert_to_time : divide by the speed of light to get time-homogene values.
    Values in meter instead
    If convert_to_time : time in sec
    """

    L = utils.extract_text_between_elements_2(
        p, r"\+sum_times/estimates", r"\-sum_times/estimates"
    )

    Lgood_stat = []
    Lgood_sat = []

    for l in L[1:-1]:
        if "EPOCHE" in l:
            cur_epoc_line = l
            cur_epoc_f = cur_epoc_line.split()
            cur_epoc = conv.mjd2dt(int(cur_epoc_f[1])) + dt.timedelta(
                seconds=int(86400 * float(cur_epoc_f[2]))
            )

        if re.match("^   [0-9]{4}.*", l):
            Lgood_stat.append([cur_epoc] + [float(e) for e in l.split()])

        if re.match("^ [A-Z][0-9]{2}.*", l):
            e = l.split()
            Lgood_sat.append([cur_epoc] + [e[0], float(e[1]), float(e[2])])

    ### stations
    DF_stat = pd.DataFrame(
        Lgood_stat, columns=["epoch", "stat", "offset", "offset_sig"]
    )
    DF_stat["stat"] = DF_stat["stat"].astype("int")
    if convert_to_time:
        DF_stat[["offset", "offset_sig"]] = (
            DF_stat[["offset", "offset_sig"]] / 299792458.0
        )

    ### satellites
    DF_sat = pd.DataFrame(Lgood_sat, columns=["epoch", "sat", "offset", "offset_sig"])
    if convert_to_time:
        DF_sat[["offset", "offset_sig"]] = (
            DF_sat[["offset", "offset_sig"]] / 299792458.0
        )

    return DF_stat, DF_sat


def read_epos_tim(tim_file_in, convert_to_sec=False):
    """
    results in microsec
    """
    F = open(tim_file_in)

    head_stop = False

    if convert_to_sec:
        koef = 10**-6
    else:
        koef = 1.0

    Val_stk = []
    for l in F:
        if re.match(r"^\*  [0-9]{4} *([0-9]{1,2} *){4}", l):
            head_stop = True
            epoc = conv.datetime_improved(*l[3:30].split())
        if head_stop and re.match("[A-Z][0-9]{2}.* [0-9]*", l):
            val = l.split()
            val[1] = float(val[1]) * koef
            val.insert(0, epoc)

            Val_stk.append(val)

    DF = pd.DataFrame(Val_stk, columns=["epoch", "sat", "offset"])

    return DF

def stations_in_epos_sta_coords_file_mono(coords_file_path):
    """
    Gives stations in a EPOS coords. file (YYYY_DDD_sta_coordinates)

    Parameters
    ----------
    coords_file_path : str
        path of the EPOS coords. file.

    Returns
    -------
    epoch : datetime
        the main mean epoch in the EPOS coords. file.
    stats_list : list
        list of 4 char station list.
    """

    site_line_list = utils.grep(coords_file_path , " SITE            m")

    stats_list = []
    mean_mjd_list = []
    for l in site_line_list:
        stat = l.split()[8].lower()
        stats_list.append(stat)
        mean_mjd = np.mean([float(l.split()[6]) , float(l.split()[7])])
        mean_mjd_list.append(mean_mjd)

    mjd_final = utils.most_common(mean_mjd_list)
    epoch = conv.mjd2dt(mjd_final)
    return epoch , stats_list

