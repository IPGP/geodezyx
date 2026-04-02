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


#   _____ _   _ ______  _____    _______ _____ _   _  _____   ______ _ _
#  / ____| \ | |  ____|/ ____|  / / ____|_   _| \ | |/ ____| |  ____(_) |
# | |    |  \| | |__  | (___   / / |  __  | | |  \| | (___   | |__   _| | ___  ___
# | |    | . ` |  __|  \___ \ / /| | |_ | | | | . ` |\___ \  |  __| | | |/ _ \/ __|
# | |____| |\  | |____ ____) / / | |__| |_| |_| |\  |____) | | |    | | |  __/\__ \
#  \_____|_| \_|______|_____/_/   \_____|_____|_| \_|_____/  |_|    |_|_|\___||___/


def read_gins_solution(filein, mode="cinematic"):
    """
    Read a GINS solution file.

    Parameters
    ----------
    filein : str
        Path of the input file.
    mode : str, optional
        Processing mode. Either 'cinematic' to return a TimeSeriePoint object
        or 'static' to return a Point object. Default is 'cinematic'.

    Returns
    -------
    TimeSeriePoint or Point
        If mode='cinematic', returns a TimeSeriePoint object containing
        the time series. If mode='static', returns the first Point object.
        Returns None if no points are found in the file.

    """

    F = open(filein)

    pts_list_tmp = []
    namestat = "XXXX"

    datexere = re.search(r"[0-9]{6}_[0-9]{6}", os.path.basename(filein))
    if not datexere:
        datexe = "991231_235959"
    else:
        datexe = datexere[0]

    ginsvers = "unknown"

    for l in F:
        f = l.split()

        if "STATION_NAME" in l:
            # get per default 4char name
            namestat = f[1] if len(f) > 1 else f[-1]

            # try to catch the 9char name in the filename
            filnam = os.path.basename(filein)
            reout = re.search(namestat + "[0-9]{2}[A-Z]{3}", filnam)
            if reout:
                namestat = reout.group(0)

        if "GINS_VERSION" in l:
            ginsvers = f[1]
        elif l[0] == "#":
            continue

        # Traw = float(f[2])

        if "XYZ_SOL" in l:
            coordstype = "XYZ"
            X = float(f[4])
            Y = float(f[6])
            Z = float(f[8])

            sX = np.sqrt(float(f[9]))
            sY = np.sqrt(float(f[10]))
            sZ = np.sqrt(float(f[11]))

            if (
                "T24:00:00.000" in f[1]
            ):  # Manage the special case if we are at the border bw 2 days
                Txyz = (
                    conv.string_date2dt(f[1][:10])
                    + dt.timedelta(days=+1)
                    + dt.timedelta(seconds=-19)
                )
            else:
                Txyz = conv.string_date2dt(f[1]) + dt.timedelta(seconds=-19)

            point = time_series.Point(
                X, Y, Z, Txyz, coordstype, sX, sY, sZ, name=namestat
            )

            point.anex["sdXY"] = float(f[12])
            point.anex["sdXZ"] = float(f[13])
            point.anex["sdYZ"] = float(f[14])
            point.anex["sol_path"] = filein
            point.anex["dateofexe"] = datexe
            point.anex["gins_version"] = ginsvers

        elif "FLH_SOL" in l:
            coordstype = "FLH"
            F = float(f[4])
            L = float(f[6])
            H = float(f[8])

            sF = np.rad2deg(np.sqrt(float(f[9])))
            sL = np.rad2deg(np.sqrt(float(f[10])))
            sH = np.sqrt(float(f[11]))

            point.F, point.L, point.H = F, L, H
            point.sF, point.sL, point.sH = sF, sL, sH

            # point.X , point.Y , point.Z = X,Y,Z

            point.anex["sdFL"] = float(f[12])
            point.anex["sdFH"] = float(f[13])
            point.anex["sdLH"] = float(f[14])

            pts_list_tmp.append(point)

    #### End of reading, export
    if not pts_list_tmp:
        log.warning("no point found in: %s. Returns None", filein)
        return None

    if mode == "cinematic":
        tsout = time_series.TimeSeriePoint()
        for point in pts_list_tmp:
            tsout.add_point(point)
        tsout.meta_set(filein, namestat)
        return tsout

    elif mode == "static":
        pt_out = pts_list_tmp[0]
        return pt_out


def read_gins_solution_multi(filein_list, return_dict=True):
    """
    Read multiple GINS solution files and return a dictionary or list of time series.

    Parameters
    ----------
    filein_list : list
        List of input file paths.
    return_dict : bool, optional
        If True, return a dictionary with station names as keys and TimeSeriePoint
        objects as values. If False, return a list of TimeSeriePoint objects.
        Default is True.

    Returns
    -------
    dict or list
        If return_dict=True, a dictionary with station names as keys and
        TimeSeriePoint objects as values. If return_dict=False, a list of
        TimeSeriePoint objects.

    """

    filein_list = sorted(filein_list)
    points_list = []
    statname_stk = []

    for fil in filein_list:
        pt_day = read_gins_solution(fil, mode="static")
        if not pt_day:
            continue
        points_list.append(pt_day)
        statname_stk.append(pt_day.name)

    statname_uniq = sorted(list(set(statname_stk)))
    ts_dict = dict()

    for point in points_list:
        if not point.name in ts_dict.keys():
            ts_dict[point.name] = time_series.TimeSeriePoint(stat=point.name)
            ts_dict[point.name].meta_set(stat=point.name, name=point.name)
        ts_dict[point.name].add_point(point)

    ## return a dictionnary of TS or a list of TS
    if return_dict:
        return ts_dict
    else:
        ts_list = []
        for k, val in ts_dict.items():
            ts_list.append(val)
        return ts_list


def read_gins(
    filein,
    kineorstatic="kine",
    flh_in_rad=True,
    force_get_convergence=False,
    kf_result=False,
):
    """
    Read a GINS listing file and extract coordinate time series.

    Parameters
    ----------
    filein : str
        Path to the input GINS listing file.
    kineorstatic : str, optional
        Type of processing. Either 'kine' for kinematic processing (returns
        TimeSeriePoint) or 'static' for static processing (returns Point).
        Default is 'kine'.
    flh_in_rad : bool, optional
        If True, FLH coordinates are in radians and will be converted to degrees.
        Default is True.
    force_get_convergence : bool, optional
        If True, forces retrieval of 'convergence' part even if
        'COORDONNEES DES STATIONS AJUSTEES EN HAUTE FREQUENCE' field is present.
        Default is False.
    kf_result : bool, optional
        If True, processes Kalman Filter results. Default is False.

    Returns
    -------
    TimeSeriePoint or Point
        If kineorstatic='kine', returns a TimeSeriePoint object containing
        the time series. If kineorstatic='static', returns the first Point object.

    Notes
    -----
    The function handles special cases such as:
    - Double convergence fields (keeps the last one)
    - Final adjustment coordinates ('COORDONNEES DES STATIONS AJUSTEES EN HAUTE FREQUENCE')
    - Midnight boundary (hour=24) handling

    """

    if ".prepars" in filein:
        log.warning(
            "%s seems to be a prepars file, are you sure of what you are doing ?",
            filein,
        )

    # Pour le static, il y a des blancs en fin de ligne ...

    if kineorstatic == "kine":
        # regex = '\[S[PLHXYZ] .*\]$'
        regex = r"\[S[PLHXYZ] .*\]$"
        tsout = time_series.TimeSeriePoint()
        if kf_result:
            regex = r"\[S[PLHXYZ][E ].*\]     $"
    elif kineorstatic == "static":
        # regex = '\[S[PLHXYZ] .*\]     $'
        regex = r"\[S[PLHXYZ][E ].*\]     $"
        tsout = time_series.TimeSeriePoint()
    else:
        log.error("ERR")

    A, B, C = 0, 0, 0
    Ta, Tb, Tc, T = 111, 222, 333, 0
    sA, sB, sC = 0, 0, 0

    if kf_result:
        regex = r"\[S[PLHXYZ][E ].*\]     $"

    # Specific si 2ble convergence
    grep_conv = utils.grep(filein, "c o n v e r g e n c e")
    if len(grep_conv) == 2:
        IPPmode = True
        converg_compt = 0
        log.info("%s have 2  c o n v e r g e n c e  fields", os.path.basename(filein))
        log.info("keeping the last one")
    else:
        IPPmode = False

    # Specific si Ajustement Final
    greped_adj = utils.grep(
        filein, "COORDONNEES DES STATIONS AJUSTEES EN HAUTE FREQUENCE"
    )
    if force_get_convergence:
        FinalAdj_mode = False
    elif len(greped_adj) != 0:
        FinalAdj_mode = True
        FinalAdj_found = False
        log.info("%s have a COORD STAT AJ EN HTE FREQ  field", os.path.basename(filein))
        log.info("     keeping this one (and not the convergence)")
    else:
        FinalAdj_mode = False

    fileopened = open(filein, encoding="ISO-8859-1")

    regex_valid_line_count = 0

    for line in fileopened:

        if re.compile("__Nom__").search(line):
            nextline = next(fileopened).split()
            namestat = nextline[3]
            Xref = float(nextline[4])
            Yref = float(nextline[5])
            Zref = float(nextline[6])

            Fref, Lref, Href = conv.xyz2geo(Xref, Yref, Zref)

        #        if 'angles en deg' in line:
        #            flh_in_rad = False

        # Specific search
        if IPPmode and re.compile("c o n v e r g e n c e").search(line):
            converg_compt = converg_compt + 1
        if FinalAdj_mode and re.compile(
            "COORDONNEES DES STATIONS AJUSTEES EN HAUTE FREQUENCE"
        ).search(line):
            FinalAdj_found = True

        # Specific skip
        if FinalAdj_mode and not FinalAdj_found:
            continue
        if IPPmode and converg_compt != 2:
            continue

        if "real" in line:
            rawexectime = line.split()[-1]

        if re.compile(regex).search(line):
            regex_valid_line_count += 1
            fields = line.split()

            if re.compile("[XYZ]").search(line):
                initype = "XYZ"
                Aref, Bref, Cref = Xref, Yref, Zref
            elif re.compile("[PLH]").search(line):
                initype = "FLH"
                Aref, Bref, Cref = Fref, Lref, Href
            else:
                log.error("wrong initype")

            if float(fields[2]) == 0:
                continue

            if fields[0] == "stations":
                pass

            # securité pour les lignes du type
            #  ------------------------------------------------------------------------------------------------------
            # stations      175  -0.574353783560234E+07  +/-   0.000000000000000E+00   1  [SX  1212001892701M005071]

            jour = int(line[108:110])
            h = int(line[110:112])
            m = int(line[112:114])
            s = int(line[114:116])
            yy = int(line[125:127])

            # gestion des années
            if 80 < yy <= 99:
                yy = yy + 1900
            else:
                yy = yy + 2000

            # pour le mois, si > sept (9), alors lettre ...
            mm = line[127]

            if mm == "O":
                mm = 10
            elif mm == "N":
                mm = 11
            elif mm == "D":
                mm = 12
            else:
                mm = int(mm)

            try:
                if h == 24:
                    # cas exceptionnel ou on doit gerer minuit
                    # on retranche l'heure dans l'int en input et on l'ajoute dans le dt
                    Ttemp = (
                        dt.datetime(yy, mm, jour, h - 1, m, s)
                        + dt.timedelta(seconds=-19)
                        + dt.timedelta(hours=1)
                    )
                else:
                    Ttemp = dt.datetime(yy, mm, jour, h, m, s) + dt.timedelta(
                        seconds=-19
                    )

                if line[105] == "X" or line[105] == "p":
                    Ta = Ttemp
                    A = float(fields[3])
                    sA = float(fields[4])

                if line[105] == "Y" or line[105] == "L":
                    Tb = Ttemp
                    B = float(fields[3])
                    sB = float(fields[4])

                if line[105] == "Z" or line[105] == "H":
                    Tc = Ttemp
                    C = float(fields[3])
                    sC = float(fields[4])
            except ValueError as err:
                log.error("yy,mm,jour,h,m,s, %s %s %s %s %S ", yy, mm, jour, h, m, s)
                raise err

            if Ta == Tb == Tc:
                T = Ta
                if initype == "FLH" and not FinalAdj_mode:
                    if flh_in_rad:
                        A = np.rad2deg(A)
                        B = np.rad2deg(B)
                        sA = np.rad2deg(sA)
                        sB = np.rad2deg(sB)
                if kf_result and 0:
                    A = A + Aref
                    B = B + Bref
                    C = C + Cref
                point = time_series.Point(
                    A, B, C, T, initype, sA, sB, sC, name=namestat
                )

                Ta = 111
                Tb = 222
                Tc = 333

                if kineorstatic == "static":
                    return point
                elif kineorstatic == "kine":
                    tsout.add_point(point)
                else:
                    log.error("ERROR")

    if regex_valid_line_count == 0:
        log.warning("no valid line (with regex check) was found !!!")

    tsout.anex["exec_time"] = rawexectime
    tsout.meta_set(filein, namestat)
    return tsout


def gins_read_time(line):
    """
    Extract and parse time information from a GINS listing line.

    Parameters
    ----------
    line : str
        A line from a GINS listing file containing time information at
        fixed character positions.

    Returns
    -------
    datetime.datetime
        The parsed datetime object. Returns 1970-01-01 if parsing fails.

    Notes
    -----
    Extracts time components from fixed character positions in the line:
    - Characters 108-110: day (jour)
    - Characters 110-112: hour (h)
    - Characters 112-114: minute (m)
    - Characters 114-116: second (s)
    - Characters 125-127: year (yy)

    Month is parsed from character 127 where letters O, N, D represent
    October, November, December respectively.

    """
    jour = int(line[108:110])
    h = int(line[110:112])
    m = int(line[112:114])
    s = int(line[114:116])
    yy = int(line[125:127])

    # gestion des années
    if 80 < yy <= 99:
        yy = yy + 1900
    else:
        yy = yy + 2000

    # pour le mois, si > sept (9), alors lettre ...
    mm = line[127]

    if mm == "O":
        mm = 10
    elif mm == "N":
        mm = 11
    elif mm == "D":
        mm = 12
    else:
        mm = int(mm)

    try:
        if h == 24:
            # cas exceptionnel ou on doit gerer minuit
            # on retranche l'heure dans l'int en input et on l'ajoute dans le dt
            Ttemp = (
                dt.datetime(yy, mm, jour, h - 1, m, s)
                + dt.timedelta(seconds=-19)
                + dt.timedelta(hours=1)
            )
        else:
            Ttemp = dt.datetime(yy, mm, jour, h, m, s) + dt.timedelta(seconds=-19)
    except:
        Ttemp = dt.datetime(1970, 1, 1)

    T = Ttemp

    return T


def gins_read_MZB(filein, return_df=False):
    """
    Read Mean Zonal Bias (MZB) from a GINS listing file.

    Parameters
    ----------
    filein : str
        Path to the input GINS listing file.
    return_df : bool, optional
        If True, return results as a pandas DataFrame. If False, return as
        separate lists. Default is False.

    Returns
    -------
    tuple or DataFrame
        If return_df=False: tuple of (Tstk, MZBstk, sMZBstk, NameStat)
            - Tstk : list of datetime objects
            - MZBstk : list of MZB values
            - sMZBstk : list of MZB standard deviations
            - NameStat : list of station names

        If return_df=True: pandas DataFrame with columns:
            - epoch : datetime
            - mzb : MZB value (float)
            - mzb_std : MZB standard deviation (float)
            - site : station name

    """

    F = open(filein)

    regex = r"\[MZB.*\]     $"

    Tstk = []
    MZBstk = []
    sMZBstk = []
    NameStat = []

    for line in F:
        if re.compile("__Nom__").search(line):
            nextline = next(F).split()
            namestat = nextline[3]
            Xref = float(nextline[4])
            Yref = float(nextline[5])
            Zref = float(nextline[6])

        if re.compile(regex).search(line):
            fields = line.split()

            # [MZB  801980015051301 GPS]

            Traw = fields[6][-8:]
            yy = int(Traw[0:2])
            mm = int(Traw[2:4])
            dd = int(Traw[4:6])
            hh = int(Traw[6:8])

            if 80 < yy <= 99:
                yy = yy + 1900
            else:
                yy = yy + 2000

            T = dt.datetime(yy, mm, dd, hh)

            MZB = float(fields[3])
            sMZB = float(fields[4])

            Tstk.append(T)
            MZBstk.append(MZB)
            sMZBstk.append(sMZB)
            NameStat.append(namestat)

    if not return_df:
        return Tstk, MZBstk, sMZBstk, NameStat
    else:
        DF = pd.DataFrame((Tstk, MZBstk, sMZBstk, NameStat))
        DF = DF.T
        DF.columns = ("epoch", "mzb", "mzb_std", "site")
        DF.mzb = DF.mzb.astype(float)
        DF.mzb_std = DF.mzb_std.astype(float)
        return DF


def gins_readTROPOZ(filein):
    """
    Read TROPOZ (zenith tropospheric delay) from a GINS listing file.

    Parameters
    ----------
    filein : str
        Path to the input GINS listing file.

    Returns
    -------
    DataFrame
        pandas DataFrame with columns:
        - jjul_cnes : CNES Julian day number
        - tropoz_std : TROPOZ standard deviation
        - tropoz : TROPOZ value
        - epoch : datetime of the measurement (rounded to 1 second)

    """

    L = utils.grep(filein, "TROPOZ COR_ZEN_ESTIM")
    DF = pd.DataFrame([e.split()[2:] for e in L]).astype(float)
    DF.columns = ["jjul_cnes", "tropoz_std", "tropoz"]
    DF["epoch"] = conv.jjul_cnes2dt(DF["jjul_cnes"]).dt.round("1s") - dt.timedelta(
        seconds=19
    )

    return DF


def write_ATM_GAMIT(Tstk, MZBstk, sMZBstk, namestat, file_out):
    """
    Write atmospheric data in GAMIT ATM_ZEN format.

    Parameters
    ----------
    Tstk : list
        List of datetime objects.
    MZBstk : list
        List of Mean Zonal Bias (MZB) values.
    sMZBstk : list
        List of MZB standard deviations.
    namestat : str
        Station name.
    file_out : str
        Path to the output file.

    Returns
    -------
    str
        Path to the output file.

    """
    Fout = open(file_out, "w+")
    for T, mzb, smzb in zip(Tstk, MZBstk, sMZBstk):
        yy = T.year
        mm = T.month
        dd = T.day
        hh = T.hour
        Line = "ATM_ZEN X {}  1 {:4} {:2} {:2} {:2}  0  {:6.4f} +-   {:6.4f}    {:6.4f}\n".format(
            namestat.upper(), yy, mm, dd, hh, mzb, smzb, mzb
        )
        Fout.write(Line)
    Fout.close()
    return file_out


def MZB_GINS_2_ATM_GAMIT(listing_in, path_out):
    """
    Convert Mean Zonal Bias (MZB) from GINS listing to GAMIT ATM format.

    Parameters
    ----------
    listing_in : str
        Path to the input GINS listing file.
    path_out : str
        Path to the output directory.

    Returns
    -------
    str
        Path to the generated output file.

    """
    Tstk, MZBstk, sMZBstk, namestat = gins_read_MZB(listing_in)
    doy, yy = conv.dt2doy_year(Tstk[0], str)
    file_out = os.path.join(
        path_out, "_".join(("MZB_GINS_2_ATM_GAMIT", namestat, doy, yy, ".txt"))
    )
    write_ATM_GAMIT(Tstk, MZBstk, sMZBstk, namestat, file_out)

    return file_out


def read_gins_wrapper(input_list_or_path, flh_in_rad=True):
    """
    Read multiple GINS listing files and return a list of time series.

    Parameters
    ----------
    input_list_or_path : str or list
        Either a list of file paths or a path glob pattern to match GINS files.
    flh_in_rad : bool, optional
        If True, FLH coordinates are in radians. Default is True.

    Returns
    -------
    list
        List of TimeSeriePoint objects, one per GINS listing file.

    Notes
    -----
    Files without the .gins extension are automatically excluded from processing.

    """

    if type(input_list_or_path) is str:
        gins_listings_list = glob.glob(input_list_or_path)
    else:
        gins_listings_list = input_list_or_path

    tslis = []
    for f in gins_listings_list:
        if not ".gins" in f:
            log.warning("WARN : no .gins ext, skipping")
            continue

        ts = read_gins(f, flh_in_rad=flh_in_rad)
        tslis.append(ts)

    return tslis


def convert_sp3_clk_2_GINS_clk(
    sp3_path_in, clk_gins_out, interpo_30sec=True, return_as_DF=True
):
    """
    Convert SP3 satellite clock data to GINS clock format.

    Parameters
    ----------
    sp3_path_in : str
        Path to the input SP3 file.
    clk_gins_out : str
        Path to the output GINS clock file.
    interpo_30sec : bool, optional
        If True, interpolate clock data to 30-second intervals. Default is True.
    return_as_DF : bool, optional
        If True, return results as a pandas DataFrame. If False, return the
        path to the output file. Default is True.

    Returns
    -------
    str or DataFrame
        If return_as_DF=True, returns a pandas DataFrame with columns:
            - epoch : datetime
            - sv : satellite vehicle number
            - clk : clock value

        If return_as_DF=False, returns the path to the output file.

    Notes
    -----
    This function is in beta status and currently only supports GPS clocks.

    """
    DF = files_rw.read_sp3(sp3_path_in)

    Fout = open(clk_gins_out, "w+")

    def write_GINS_signaletic_elt_clk(
        dt_in, sv, signaletic_name="MNG", add_19sec_to_dt_in=True
    ):
        """
        very beta, only for GPS clk
        """

        if add_19sec_to_dt_in:
            dt_work = dt_in + dt.timedelta(seconds=19)
        else:
            dt_work = dt_in

        jjul = conv.dt2jjul_cnes(dt_work)
        sec_in_day = (dt_work - conv.jjul_cnes2dt(jjul)).seconds

        # MNG0000000jjjjjcccccnnnn

        outstr = (
            "[MNG0000000"
            + str(jjul)
            + str(sec_in_day).zfill(5)
            + "GP"
            + str(sv).zfill(2)
            + "]"
        )

        return outstr

    c = 299792458

    if interpo_30sec:
        Epoc_work = []
        Sv_work = []
        Clk_work = []

        for sv in sorted(DF["sv"].unique()):
            DFsv = DF[DF["sv"] == sv]

            Epoc_inp = np.array(conv.dt2posix(DFsv["epoch"]))
            Clk_inp = np.array(DFsv["clk"])

            Epoc_interp = np.arange(np.min(Epoc_inp), np.max(Epoc_inp), 30)
            Epoc_interp_dt = conv.posix2dt(Epoc_interp)

            I = scipy.interpolate.interp1d(Epoc_inp, Clk_inp)

            Clk_interp = I(Epoc_interp)

            Sv_work = Sv_work + [sv] * len(Epoc_interp)
            Epoc_work = Epoc_work + list(Epoc_interp_dt)
            Clk_work = Clk_work + list(Clk_interp)

    else:
        Epoc_work = DF["epoch"]
        Sv_work = DF["sv"]
        Clk_work = DF["clk"]

    DF_work = pd.DataFrame(
        list(zip(Epoc_work, Sv_work, Clk_work)), columns=("epoch", "sv", "clk")
    )

    DF_work.sort_values(["epoch", "sv"], inplace=True)

    for epoc, sv, clk in zip(DF_work["epoch"], DF_work["sv"], DF_work["clk"]):
        signaletik = write_GINS_signaletic_elt_clk(epoc, sv, clk)

        str_final = " 0 0 {:}  {:+17.15e} {:+17.15e}\n".format(
            signaletik, clk * 10**-6 * c, 0
        )

        Fout.write(str_final)

    if not return_as_DF:
        return clk_gins_out
    else:
        return DF_work


def read_gins_multi_raw_listings(filelistin, kineorstatic="static", flh_in_rad=True):
    """
    Read multiple raw GINS listing files and return a single time series.

    Parameters
    ----------
    filelistin : list
        List of input GINS listing file paths.
    kineorstatic : str, optional
        Type of processing. Either 'static' or 'kine' for kinematic.
        Default is 'static'.
    flh_in_rad : bool, optional
        If True, FLH coordinates are in radians. Default is True.

    Returns
    -------
    TimeSeriePoint
        A single TimeSeriePoint object containing points from all input files.

    Notes
    -----
    Files must contain the "c o n v e r g e n c e" field to be processed.
    All points should have the same station name.

    """

    tsout = time_series.TimeSeriePoint()
    refname = "RIEN"
    if kineorstatic == "static":
        for filein in filelistin:
            log.info(filein)
            if not utils.check_regex(filein, "c o n v e r g e n c e"):
                continue
            pt = read_gins(filein, kineorstatic="static", flh_in_rad=flh_in_rad)
            if refname == "RIEN":
                refname = pt.name
            if refname != pt.name:
                log.warning("nom de stat. != reference")
            tsout.add_point(pt)
        tsout.meta_set(stat=refname)

    elif kineorstatic == "kine":
        tsoutlis = []
        for filein in filelistin:
            log.info(filein)
            if not utils.check_regex(filein, "c o n v e r g e n c e"):
                continue
            ts = read_gins(filein, kineorstatic="kine", flh_in_rad=flh_in_rad)
            tsoutlis.append(ts)
        tsout = time_series.merge_ts(tsoutlis)

    else:
        log.error("check kineorstatic keyword")

    tsout.sort()
    return tsout


def read_gins_multi_extracted(filelistin, flh_in_rad=True):
    """
    Read extracted GINS listing files and return a time series.

    Parameters
    ----------
    filelistin : list
        List of input file paths. Must contain exactly 3 files, one for each
        coordinate component (X/p, Y/L, Z/H).
    flh_in_rad : bool, optional
        If True, FLH coordinates are in radians and will be converted to degrees.
        Default is True.

    Returns
    -------
    TimeSeriePoint or None
        A TimeSeriePoint object containing the time series if successful.
        Returns None if the list does not contain exactly 3 files.

    Notes
    -----
    The input files are expected to have the format: S<HLP>__ddhhiissxxxxxxxxxyym.HOUE
    where one file contains X/p, one contains Y/L, and one contains Z/H coordinates.

    """
    tsout = time_series.TimeSeriePoint()
    if len(filelistin) != 3:
        log.error("listfilein != 3 elts")
        return None
    statnameset = list(set([f.split(".")[-1] for f in filelistin]))
    statname = statnameset[0]
    if len(statnameset) != 1:
        log.error("len(statnameset) != 1")
    fileopenedlist = [open(f, "r+") for f in filelistin]
    coortypelist = [os.path.basename(f)[1] for f in filelistin]
    log.info("%s %s", coortypelist, filelistin)
    if "X" in coortypelist:
        initype = "XYZ"
        ia = coortypelist.index("X")
        ib = coortypelist.index("Y")
        ic = coortypelist.index("Z")
    elif "p" in coortypelist:
        initype = "FLH"
        ia = coortypelist.index("p")
        ib = coortypelist.index("L")
        ic = coortypelist.index("H")

    for lf in zip(*fileopenedlist):
        # OLD MODE : MAVAIS GESTION DE LA DATE
        #        t1 = lf[0].split()[0]
        #        t2 = lf[1].split()[0]
        #        t3 = lf[2].split()[0]
        #        if not (t1 == t2 == t3):
        #            print "WARN : read_gins_multi_extracted : not (t1 == t2 == t3)"
        #
        fields = lf[0].split()
        line = lf[0]
        blok = fields[-1]
        jour = int(blok[0:2])
        h = int(blok[2:4])
        m = int(blok[4:6])
        s = int(blok[6:8])
        yy = int(blok[-4:-2])
        # gestion des années
        if 80 < yy <= 99:
            yy = yy + 1900
        else:
            yy = yy + 2000

        # pour le mois, si > sept (9), alors lettre ...
        mm = blok[-2]

        if mm == "O":
            mm = 10
        elif mm == "N":
            mm = 11
        elif mm == "D":
            mm = 12
        else:
            mm = int(mm)

        t1 = dt.datetime(yy, mm, jour, h, m, s) + dt.timedelta(seconds=-19)

        lfiasplit = lf[ia].split()
        lfibsplit = lf[ib].split()
        lficsplit = lf[ic].split()

        A = float(lfiasplit[4])
        B = float(lfibsplit[4])
        C = float(lficsplit[4])
        sA = float(lfiasplit[5])
        sB = float(lfibsplit[5])
        sC = float(lficsplit[5])
        T = t1

        if initype == "FLH" and flh_in_rad:
            A = np.rad2deg(A)
            B = np.rad2deg(B)
            sA = np.rad2deg(sA)
            sB = np.rad2deg(sB)

        pt = time_series.Point(A, B, C, T, initype, sA, sB, sC)

        tsout.add_point(pt)
    tsout.meta_set(stat=statname)
    tsout.sort()

    return tsout


def read_gins_double_diff(filein):
    """
    Extract point objects from a GINS double difference listing file.

    Parameters
    ----------
    filein : str
        Path to the input GINS double difference listing file.

    Returns
    -------
    list or None
        List of Point objects extracted from the file.
        Returns None if the file has no convergence section.

    """

    if utils.grep(filein, "c o n v e r g e n c e") == "":
        log.error("%s have no convergence, return None", filein)
        return None

    fileopened = open(filein)

    statlist = []
    rawdatalist = []
    timelist = []

    regstat = re.compile(r"[0-9]{5}[A-Z][0-9]{3}  [0-9]{7} .*$")
    regresult = re.compile(r"\[S[PLHXYZ][E ].*\]     $")

    for line in fileopened:
        # finding
        if regstat.search(line):
            stat = line.split()[3]
            statlist.append(stat)

        if regresult.search(line):
            fields = line.split()
            rawdatalist.append([float(e) for e in fields[:-2]])
            timelist.append(gins_read_time(line))

    rawdatatab = np.vstack(rawdatalist)

    Pstk = []

    for istat, stat in enumerate(statlist):
        P = time_series.Point(
            rawdatatab[3 * istat, 3],
            rawdatatab[3 * istat + 1, 3],
            rawdatatab[3 * istat + 2, 3],
            timelist[3 * istat],
            "XYZ",
            rawdatatab[3 * istat, 4],
            rawdatatab[3 * istat, 4],
            rawdatatab[3 * istat, 4],
            name=stat,
        )

        Pstk.append(P)

    return Pstk


def read_gins_double_diff_multi(filelistin):
    """
    Read multiple GINS double difference listing files and return a dictionary of time series.

    Parameters
    ----------
    filelistin : list
        List of input GINS double difference listing file paths.

    Returns
    -------
    dict
        Dictionary with station names as keys and TimeSeriePoint objects as values.

    """

    Ptsstk = []
    for filein in filelistin:
        Pstk = read_gins_double_diff(filein)
        if Pstk != None:
            Ptsstk = Ptsstk + Pstk

    statlist = set([e.name for e in Ptsstk])

    tsdico = dict()

    for stat in statlist:
        tsdico[stat] = time_series.TimeSeriePoint(stat=stat)

    for pt in Ptsstk:
        tsdico[pt.name].add_point(pt)

    for ts in tsdico.values():
        ts.sort()

    return tsdico


def read_spotgins_quick(p):
    """
    Read SpotGINS quick format file.

    Parameters
    ----------
    p : str
        Path to the input SpotGINS file.

    Returns
    -------
    DataFrame
        pandas DataFrame with columns:
        - E : East coordinate
        - N : North coordinate
        - U : Up coordinate
        Indexed by date (datetime).

    """
    df = pd.read_csv(p, comment="#", header=None, sep=r"\s+")

    with open(p) as f:
        header_line = [line for line in f if line.startswith("#")][-1]
    df.columns = header_line.strip().split()
    df = df[df.columns[:4]]
    df.columns = ["mjd", "E", "N", "U"]
    df["date"] = conv.mjd2dt(df["mjd"])
    df.set_index("date", inplace=True)
    df = df.drop("mjd", axis=1)

    return df


def diff_spotgins_quick(df1, df2):
    """
    Compute the difference between two SpotGINS coordinate dataframes.

    Parameters
    ----------
    df1 : DataFrame
        First SpotGINS dataframe with E, N, U columns indexed by date.
    df2 : DataFrame
        Second SpotGINS dataframe with E, N, U columns indexed by date.

    Returns
    -------
    DataFrame
        DataFrame containing the differences (df1 - df2) with rounded hourly indices,
        with NaN values removed.

    """
    def _rndidx(d):
        duse = d.index.to_series().dt.round("1h").values
        dout = d.set_index(duse)
        return dout

    def _jjul(d):
        d = _rndidx(d)
        duse = d.index.to_series()
        duse = conv.numpy_dt2dt(duse)
        jjul = conv.dt2jjul_cnes(duse, True)

        dout = d.set_index(jjul)  # + jjul[1]/86400)
        return dout

    lbda = _rndidx

    dfout = (lbda(df1) - lbda(df2)).dropna()
    # print(dfout.to_string())
    return dfout
