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


def read_all_points(filein):
    """selectionne automatiquement le type de fichier brut en entrée
    INPUT  : chemin du fichier brut de POINTS
    OUTPUT : Une TimeSeriePoint"""

    firstline = open(filein).readline()

    if re.compile("RTKLIB").search(firstline):
        tsout = read_rtklib(filein)

    elif re.compile("Kinematic Processing").search(firstline):
        tsout = read_gipsy_apps(filein)

    elif re.compile("tdp.llh").search(filein):
        tsout = read_gipsy_bosser(filein)

    elif (
        re.compile("STA").search(firstline)
        or re.compile("tdp").search(filein)
        or re.compile("TRPAZ").search(firstline)
    ):
        tsout = read_gipsy_tdp(filein)

    elif re.compile("YY  MM DD HR MIN").search(firstline):
        tsout = read_track(filein)

    # elif re.compile('Latitude').search(firstline):
    #     tsout = read_sonardyne_posi(filein)

    # elif re.compile('Heading').search(firstline):
    #     tsout = read_sonardyne_attitude(filein)

    elif re.compile(r"\*\*\* warning").search(firstline):
        tsout = read_gins(filein, "kine")

    elif re.compile("#GINS_VERSION").search(firstline):
        tsout = read_gins_solution(filein)

    # elif re.compile('OCTANS_ATTITUDE').search(firstline):
    #     tsout = read_qinsy(filein, 2014, 0o4, 0o4)

    elif re.compile("PBO Station Position Time Series").search(firstline):
        tsout = read_pbo_pos(filein)

    elif re.compile("latitude_degre_decimal").search(firstline):
        tsout = read_nrcan_csv(filein)

    elif re.compile("--------------------------------------------------").search(
        firstline
    ):  # NDLR : Best parser ever
        tsout = read_nrcan_pos(filein)

    else:
        log.error("pas de motif valide pour lect. auto")
        log.info(filein)
        log.info(firstline)
        raise Exception

    return tsout


def read_all_obs(filein):
    """selectionne automatiquement le type de fichier brut en entrée
    INPUT  : chemin du fichier brut de observations génériques
    OUTPUT : Une LISTE de TimeSerieObs"""

    firstline = open(filein).readline()

    # if re.compile('Heading').search(firstline):
    #     tsout = read_sonardyne_attitude(filein)
    # else:
    log.error("pas de motif valide pour lect. auto")

    return None


#  __  __ _____ _______  _______          __  __ _____ _______   ______ _ _
# |  \/  |_   _|__   __|/ / ____|   /\   |  \/  |_   _|__   __| |  ____(_) |
# | \  / | | |    | |  / / |  __   /  \  | \  / | | |    | |    | |__   _| | ___  ___
# | |\/| | | |    | | / /| | |_ | / /\ \ | |\/| | | |    | |    |  __| | | |/ _ \/ __|
# | |  | |_| |_   | |/ / | |__| |/ ____ \| |  | |_| |_   | |    | |    | | |  __/\__ \
# |_|  |_|_____|  |_/_/   \_____/_/    \_\_|  |_|_____|  |_|    |_|    |_|_|\___||___/


def read_track_2(filein, site_name=None):
    """
    Read a kinematic track file

    Parameters
    ----------
    filein : str
        path of the file.

    Returns
    -------
    DF : Pandas DataFrame

    """
    DF = pd.read_csv(filein, delim_whitespace=True, skiprows=[1, 2])
    DF.columns = [
        "year",
        "month",
        "day",
        "hour",
        "minute",
        "second",
        "dX",
        "dX_std",
        "dY",
        "dY_std",
        "dZ",
        "dZ_std",
        "rms",
        "dd",
        "atm",
        "atm_std",
        "fract_doy",
        "n_epoch",
        "BF",
        "not",
        "f",
        "rho_ua",
        "null",
    ]

    DF.drop("null", axis=1, inplace=True)

    Epoch = conv.ymdhms2dt(DF.year, DF.month, DF.day, DF.hour, DF.minute, DF.second)
    DF["epoch"] = Epoch

    if site_name:
        DF["site"] = site_name()

    return DF


def read_track(filein):
    """
    Read GAMIT/TRACK File

    Parameters
    ----------
    filein : str
        input file path.

    Returns
    -------
    tsout : TimeSeries Object
        output TimeSerie.
    """

    tsout = time_series.TimeSeriePoint()

    for line in open(filein):

        if re.compile("dNorth").search(line):
            initype = "ENU"
        elif re.compile("dX").search(line):
            initype = "XYZ"

        if line[0] != " ":
            continue

        fields = line.split()

        T = dt.datetime(
            int(fields[0]),
            int(fields[1]),
            int(fields[2]),
            int(fields[3]),
            int(fields[4]),
            int(fields[5].split(".")[0]),
            int(np.round(float(fields[5].split(".")[1]), 4)),
        )
        A = float(fields[6])
        B = float(fields[8])
        C = float(fields[10])

        sA = float(fields[7])
        sB = float(fields[9])
        sC = float(fields[11])

        # On inverse A et B car NEU => ENU
        if initype == "ENU":
            point = time_series.Point(B, A, C, T, initype, sB, sA, sC)
        elif initype == "XYZ":
            point = time_series.Point(A, B, C, T, initype, sA, sB, sC)
        else:
            "ERR : track_read : bad initype"

        tsout.add_point(point)

    if initype == "ENU":
        tsout.boolENU = True

    tsout.meta_set(filein)

    tsout.anex["rover"] = tsout.name.split(".")[-2].upper()

    # recherche de la base
    try:
        if ".LC" in filein:
            stat_n_base_lc_files = glob.glob(filein[:-8] + "*")
            stat_n_base_lc_files.remove(filein)
            statbase = os.path.basename(stat_n_base_lc_files[0])[-7:-3]
            tsout.anex["base"] = statbase
        elif ".L1+L2" in filein:
            stat_n_base_l1l2_files = glob.glob(filein[:-11] + "*")
            stat_n_base_l1l2_files.remove(filein)
            statbase = os.path.basename(stat_n_base_l1l2_files[0])[-10:-6]
            tsout.anex["base"] = statbase
    except:
        log.warning("unable to find the base for the TRACK experience")
        log.info(filein)
        pass

    return tsout

def read_nevada(filein, input_coords="enu"):
    """
    input_coords="enu" or "xyz"
    """

    tsout = time_series.TimeSeriePoint()

    envfile = open(filein)

    if input_coords == "enu":
        for l in envfile:
            f = l.split()

            if "site YYMMMDD" in l:
                continue
            if len(l) == 0:
                continue

            stat = f[0]

            T = conv.year_decimal2dt(float(f[2]))

            N = float(f[10])
            E = float(f[8])
            U = float(f[12])

            sN = float(f[15])
            sE = float(f[14])
            sU = float(f[16])

            point = time_series.Point(E, N, U, T, "ENU", sE, sN, sU)

            # tsout.refENU = time_series.Point()

            tsout.boolENU = True
            tsout.add_point(point)

    if input_coords == "xyz":
        for l in envfile:
            f = l.split()

            if "site YYMMMDD" in l:
                continue
            if len(l) == 0:
                continue

            stat = f[0]

            T = conv.year_decimal2dt(float(f[2]))

            X = float(f[3])
            Y = float(f[4])
            Z = float(f[5])

            sX = float(f[6])
            sY = float(f[7])
            sZ = float(f[8])

            point = time_series.Point(X, Y, Z, T, "XYZ", sX, sY, sZ)

            point.anex["Rxy"] = float(f[9])
            point.anex["Rxz"] = float(f[10])
            point.anex["Ryz"] = float(f[11])

            tsout.add_point(point)

    tsout.stat = stat

    return tsout


def read_IGS_coords(filein, initype="auto"):
    Tstk, Astk, Bstk, Cstk = [], [], [], []
    tsout = time_series.TimeSeriePoint()
    for l in open(filein):
        f = l.split()

        if initype == "auto":
            if "plh" in os.path.basename(filein):
                initype = "FLH"
            elif "xyz" in os.path.basename(filein):
                initype = "XYZ"

        T = conv.mjd2dt(float(f[3]))
        A = float(f[6])
        B = float(f[7])
        C = float(f[8])
        sA = float(f[9])
        sB = float(f[10])
        sC = float(f[11])

        initype = "FLH"

        pt = time_series.Point(A, B, C, T, initype, sA, sB, sC, f[0])
        tsout.add_point(pt)

    tsout.meta_set(filein, f[0])
    return tsout


def sorting_a_calais_file(openedfile):
    openedfile.seek(0)
    T, A, sA, STAT = [], [], [], []
    for l in openedfile:
        f = l.split()
        T.append(float(f[0]))
        A.append(float(f[1]))
        sA.append(float(f[2]))
        STAT.append(str(f[3]))

    DATA = [T, A, sA, STAT]
    DATA2 = utils.sort_table(DATA, 0)

    return DATA2


def read_calais(filelist):
    """filelistin est une liste de 3 fichier E N & U"""

    if filelist == []:
        log.warning("files list empty , exiting ...")
        return None
    filelist.sort()

    fileopenedlist = [open(f) for f in filelist]

    sorted_data_lis = [sorting_a_calais_file(fil) for fil in fileopenedlist]

    statnameset = list(set([f.split(".")[0] for f in filelist]))
    statname = os.path.basename(statnameset[0])

    # LOADING ALL DATA IN A BIG MATRIX
    bigT = np.array(sorted(set(np.hstack([np.array(sd[0]) for sd in sorted_data_lis]))))
    bigDATA = np.empty((len(bigT), 3))
    bigDATA.fill(np.nan)
    bigDATAsigma = np.empty((len(bigT), 3))
    bigDATAsigma.fill(np.nan)

    for i, bt in enumerate(bigT):
        for j, data in enumerate(sorted_data_lis):
            for t, d, sd, stat in zip(*data):
                if t == bt:
                    bigDATA[i, j] = d
                    bigDATAsigma[i, j] = sd

    DATA = np.hstack((bigDATA / 1000.0, bigDATAsigma / 100.0))

    ptslist = []

    # MAKING POINTS
    for i in range(DATA.shape[0]):
        pt = time_series.Point(
            conv.year_decimal2dt(bigT[i]), "ENU", DATA[i, 3], DATA[i, 4], DATA[i, 5]
        )
        ptslist.append(pt)

    tsout = time_series.TimeSeriePoint()

    for pt in ptslist:
        tsout.add_point(pt)

    # FINDING DISCONT
    # finding discont directly in files

    # Finding the composant with max of data
    lendata = [len(sd[0]) for sd in sorted_data_lis]
    ii = lendata.index(max(lendata))

    discont = []
    T = sorted_data_lis[ii][0]
    STAT = sorted_data_lis[ii][-1]
    for i in range(len(STAT) - 1):
        if STAT[i + 1] != STAT[i]:
            discont.append(T[i + 1])

    tsout.set_discont(discont)
    tsout.meta_set(stat=statname)
    tsout.boolENU = True
    tsout.sort()

    return tsout


def read_renag_synthetic(filein, discont_file_in=None):
    tsout = time_series.TimeSeriePoint()

    fil = open(filein)

    for l in fil:

        if l[0] == "#":
            continue

        f = l.split()

        T = conv.year_decimal2dt(float(f[0]))

        N = float(f[1])
        E = float(f[2])
        U = float(f[3])

        sN = float(f[4])
        sE = float(f[5])
        sU = float(f[6])

        point = time_series.Point(E, N, U, T, "ENU", sE, sN, sU)

        # tsout.refENU = time_series.Point()

        tsout.boolENU = True
        tsout.add_point(point)

    if discont_file_in:
        DiscontInp = open(discont_file_in)

        Discont = []
        for l in DiscontInp:
            if l[0] == "#":
                continue
            f = l.split()
            try:
                Discont.append(conv.doy2dt(int(f[0]), int(f[1])))
            except:
                log.warning("something went wrong during discont. file reading")
                pass

        Discont = sorted(Discont)
        tsout.set_discont(Discont)

    stat_name = os.path.basename(filein).split(".")[0].split(".")[0]
    tsout.meta_set(path=filein, stat=stat_name, name=stat_name)

    return tsout


def read_jump_file(filein, returned_events=("S", "E", "D")):
    """
    From a "Jump" File (p. Sakic internal file)
    Return a dictionnairy with events

    Parameters
    ----------
    filein : str
        path of the Jump File

    returned_events : tuple or list
        contains the inital letter of the event type which will be stored in the dico
        (See below)

    Returns
    -------
    jump_dico : dict of dict of datetime
        outputed events, in the form jump_dico["STAT"]["L"] = datetime
        where L is the inital letter of the event type

    Note
    ----
    A jump file contains infos like this :

    #>>> STAT S 2000 001
    #>>> STAT E 2001 001
    #>>> STAT D 2000 06 01

    it can manage YEAR DOY or YEAR MM DD or DECIMAL YEAR

    A non-blank 1st column is a commented line

    After a #, it is a commentary

    event type letters :
        S : Start

        E : End

        D : Discontinuity

    """
    import pytz

    F = open(filein)
    jump_dico = dict()
    for l in F:
        l = l.split("#")[0]
        f = l.split()

        ## Skip comment lines
        if (not f) or (l[0] != " ") or ("#" in l):
            continue
        else:
            stat = f[0]
            event = f[1]
            ## Create the key for the station
            if not f[0] in jump_dico.keys():
                jump_dico[stat] = dict()
                for rtn_evt in returned_events:
                    jump_dico[stat][rtn_evt] = []

            ## Fill the dico
            if event in returned_events:
                if len(f) == 4:  ### DOY
                    date = conv.doy2dt(int(f[2]), int(f[3]))
                elif len(f) == 5:  ### YYYY MM DD
                    date = dt.datetime(*[int(e) for e in f[2:]])
                elif len(f) == 3:  ### DECIMAL YEAR
                    date = conv.year_decimal2dt(float(f[2]))

                date_tz = date.replace(tzinfo=pytz.UTC)
                jump_dico[stat][event].append(date_tz)
    return jump_dico


def read_nav_step1_geodesea(filein):
    M = np.loadtxt(filein)
    tsout = time_series.TimeSeriePoint()

    for m in M:
        pt = time_series.Point(
            np.rad2deg(m[1]), np.rad2deg(m[2]), m[3], m[0], initype="FLH"
        )
        tsout.add_point(pt)

    return tsout


#  _   _ _____   _____          _   _   ______ _ _
# | \ | |  __ \ / ____|   /\   | \ | | |  ____(_) |
# |  \| | |__) | |       /  \  |  \| | | |__   _| | ___  ___
# | . ` |  _  /| |      / /\ \ | . ` | |  __| | | |/ _ \/ __|
# | |\  | | \ \| |____ / ____ \| |\  | | |    | | |  __/\__ \
# |_| \_|_|  \_\\_____/_/    \_\_| \_| |_|    |_|_|\___||___/


def read_nrcan_csv(filein, associated_ps_file="", statname=""):
    """
    associated_ps_file is highly recommanded
    because of the time managing

    WARN : Must be avoided b/c of the weak decimal precision of the angles !!!
    """

    if statname == "":
        statname = os.path.basename(filein)[0:4]

    pdcsv = pd.read_csv(filein)

    F = np.array(pdcsv["latitude_degre_decimal"])
    L = np.array(pdcsv["longitude_degre_decimal"])
    H = np.array(pdcsv["hauteur_ellipsoidale_m"])
    heure = np.array(pdcsv["heure_decimal"])
    doy = np.array(pdcsv["jour_de_l_annee"])
    year = np.array(pdcsv["annee"])

    T = conv.doy2dt(year, doy, heure)

    if associated_ps_file != "":
        T = []
        for l in open(associated_ps_file):
            if "BWD" in l:
                f = l.split()
                t = conv.date_string_2_dt(f[4] + " " + f[5])
                T.append(t)
        if statname == "":
            statname = f[2]

    tsout = time_series.TimeSeriePoint()
    for f, l, h, t in zip(F, L, H, T):
        point = time_series.Point(f, l, h, t, "FLH", name=statname)
        tsout.add_point(point)

    tsout.meta_set(filein, statname)

    return tsout


def read_nrcan_pos(filein):
    """
    .pos file are more precise than .csv, should be used !
    """
    tsout = time_series.TimeSeriePoint()
    start_read = False

    for l in open(filein):
        if l[0:3] == "DIR":
            start_read = True
            lhead = l.split()
            i_lat_d = lhead.index("LATDD")
            i_lat_m = lhead.index("LATMN")
            i_lat_s = lhead.index("LATSS")

            i_lon_d = lhead.index("LONDD")
            i_lon_m = lhead.index("LONMN")
            i_lon_s = lhead.index("LONSS")

            i_h = lhead.index("HGT(m)")

            i_slat = lhead.index("SDLAT(95%)")
            i_slon = lhead.index("SDLON(95%)")
            i_sh = lhead.index("SDHGT(95%)")

            continue
        elif not start_read:
            continue
        else:
            f = l.split()
            lat = (
                np.abs(float(f[i_lat_d]))
                + 1 / 60.0 * float(f[i_lat_m])
                + 1 / 3600.0 * float(f[i_lat_s])
            ) * np.sign(float(f[i_lat_d]))
            lon = (
                np.abs(float(f[i_lon_d]))
                + 1 / 60.0 * float(f[i_lon_m])
                + 1 / 3600.0 * float(f[i_lon_s])
            ) * np.sign(float(f[i_lon_d]))
            h = float(f[i_h])

            ### old and useless conversion (2021-01)
            # sE = float(f[15])
            # sN = float(f[16])
            # sU = float(f[17])
            # slat , slon , sh = conv.sigma_enu2geo(lat,lon,h,sE,sN,sU)

            slat, slon, sh = float(f[i_slat]), float(f[i_slon]), float(f[i_sh])

            t = conv.date_string_2_dt(f[4] + " " + f[5])

            pt = time_series.Point(lat, lon, h, t, "FLH", slat, slon, sh, name=f[2])
            tsout.add_point(pt)

    tsout.meta_set(filein, f[2])
    return tsout


def read_pbo_pos(filein):
    filobj = open(filein, "r")
    tsout = time_series.TimeSeriePoint()
    header = True
    for line in filobj:
        if not header:
            f = line.split()
            f2 = [float(e) for e in f[:-1]]
            t = dt.datetime(
                int(f[0][0:4]),
                int(f[0][4:6]),
                int(f[0][6:]),
                int(f[1][0:2]),
                int(f[1][2:4]),
                int(f[1][4:]),
            )
            pt = time_series.Point(f2[3], f2[4], f2[5], t, "XYZ", f2[6], f2[7], f2[8])
            # pt.FLHset(f2[12], f2[13], f2[14]) # useless (251018)
            pt.ENUset(
                f2[16], f2[15], f2[17], f2[19], f2[18], f2[20]
            )  # useless? (251018)
            pt.anex["sdXY"] = f2[9]
            pt.anex["sdXZ"] = f2[10]
            pt.anex["sdYZ"] = f2[11]
            pt.anex["sdEN"] = f2[-4]
            pt.anex["sdNU"] = f2[-3]
            pt.anex["sdEU"] = f2[-2]
            tsout.add_point(pt)
        if line[0] == "*":
            header = False
    tsout.boolENU = True
    tsout.meta_set(filein)
    tsout.anex["refXYZ"] = (0, 0, 0)

    tsout.stat = os.path.basename(filein)[0:4]
    return tsout





def read_hector_neu(filein):
    log.warning("XYZ/FLH conversion not implemented")
    M = np.loadtxt(filein)
    stat = utils.grep(filein, "Site :", only_first_occur=True).split()[3]
    tsout = time_series.ts_from_list(
        M[:, 2],
        M[:, 1],
        M[:, 3],
        conv.year_decimal2dt(M[:, 0]),
        "ENU",
        M[:, 4],
        M[:, 5],
        M[:, 6],
        stat=stat,
        name=stat,
    )

    return tsout


#  _______ _    _  _____     _______ _____   ____   ____  _____   _____   ______ _ _
# |__   __| |  | |/ ____|   / / ____|  __ \ / __ \ / __ \|  __ \ / ____| |  ____(_) |
#    | |  | |  | | |  __   / / |  __| |__) | |  | | |  | | |__) | (___   | |__   _| | ___  ___
#    | |  | |  | | | |_ | / /| | |_ |  _  /| |  | | |  | |  ___/ \___ \  |  __| | | |/ _ \/ __|
#    | |  | |__| | |__| |/ / | |__| | | \ \| |__| | |__| | |     ____) | | |    | | |  __/\__ \
#    |_|   \____/ \_____/_/   \_____|_|  \_\\____/ \____/|_|    |_____/  |_|    |_|_|\___||___/


def read_groops_position(Filesin):
    if not utils.is_iterable(Filesin):
        Filesin = [Filesin]

    Statnames = list(set([os.path.basename(f)[-8:-4] for f in Filesin]))

    if len(Statnames) > 1:
        log.warn("several stations for the same TimeSerie!:" + str(Statnames))

    statname = Statnames[0]

    tsout = time_series.TimeSeriePoint()

    for filein in Filesin:
        DF = pd.read_csv(filein, skiprows=6, header=None, sep=r"\s+")
        T = conv.dt2posix(conv.mjd2dt(DF[0].values))
        X, Y, Z = DF[1], DF[2], DF[3]

        for t, x, y, z in zip(T, X, Y, Z):
            point = time_series.Point(x, y, z, t, "XYZ", name=statname)
            tsout.add_point(point)

    tsout.meta_set(stat=statname)
    tsout.sort()

    return tsout


def _pride_pppar_end_header(filein):
    colheader = 0
    stat = "XXXX"

    with open(filein) as F:
        L = F.readlines()

    for i, l in enumerate(L):
        if "STATION" in l:
            stat = l.split()[0]

        if "END OF HEADER" in l:
            colheader = i + 1
            break
    return colheader, stat


def read_pride_pppar_pos_mono(filein):
    colheader, stat_header = _pride_pppar_end_header(filein)

    df = pd.read_csv(
        filein,
        skiprows=colheader + 1,
        # delim_whitespace=True,
        sep=r"\s?\*?\s+",
        engine="python",
        header=None,
    )

    df.columns = [
        "stat",
        "Mjd",
        "X",
        "Y",
        "Z",
        "Sx",
        "Sy",
        "Sz",
        "Rxy",
        "Rxz",
        "Ryz",
        "Sig0",
        "Nobs",
    ]

    df = df.squeeze()

    T = conv.dt2posix(conv.mjd2dt(df["Mjd"]))

    fuv = df["Sig0"] ** 2  # variance of unit weight

    anex = dict()
    anex["sdXY"] = df["Rxy"] * fuv
    anex["sdXZ"] = df["Rxz"] * fuv
    anex["sdYZ"] = df["Ryz"] * fuv

    # not sure for the sigma computation
    pt = time_series.Point(
        df["X"],
        df["Y"],
        df["Z"],
        T,
        "XYZ",
        np.sqrt(df["Sx"]) * fuv,
        np.sqrt(df["Sy"]) * fuv,
        np.sqrt(df["Sz"]) * fuv,
        name=df["stat"],
        anex=anex,
    )

    return pt


def read_pride_pppar_pos(files_list_in):
    tsout = time_series.TimeSeriePoint()

    for file in files_list_in:
        try:
            pt = read_pride_pppar_pos_mono(file)
        except pd.errors.EmptyDataError as e:
            log.error("%s, %s skipped", e, file)
            continue

        tsout.add_point(pt)
    tsout.meta_set(stat=pt.name)
    tsout.sort()

    return tsout


def read_pride_pppar_kin(filein):
    stat = "XXXX"

    colheader, stat = _pride_pppar_end_header(filein)

    df = pd.read_csv(
        filein,
        skiprows=colheader + 1,
        # delim_whitespace=True,
        sep=r"\s?\*?\s+",
        engine="python",
        header=None,
    )

    t_arr = conv.mjd2dt(df[0]) + df[1].apply(lambda x: dt.timedelta(seconds=x))

    tsout = time_series.TimeSeriePoint()

    tsout = time_series.ts_from_list(
        df[2].values, df[3].values, df[4].values, t_arr, "XYZ", stat=stat, name=stat
    )

    return tsout


def read_webobs(filein, typein="txt", coordtreat=False, dropna=False):
    if coordtreat:
        lbda_colname = lambda c: c + "_treat"
    else:
        lbda_colname = lambda c: c + "ern" if c != "Up" else "Up"

    if typein == "txt":
        header = utils.grep(filein, "#")[-1][1:].strip().split()
        DF = pd.read_csv(
            filein, sep=" ", comment="#", names=header, on_bad_lines="warn"
        )
        unit_suffix = "(m)"
        ### Time  conversion
        DFtime = DF[["yyyy", "mm", "dd", "HH", "MM", "SS"]].copy()
        DFtime.columns = ["year", "month", "day", "h", "m", "s"]
        DF["T"] = pd.to_datetime(DFtime)

    elif typein == "csv":
        DF = pd.read_csv(filein, sep=";", on_bad_lines="warn")
        unit_suffix = ""
        ### Time  conversion
        DFtimedelta = pd.to_timedelta(DF.HH * 3600 + DF.MM * 60 + DF.SS, unit="S")
        DF["T"] = pd.to_datetime(DF["yyyy-mm-dd"]) + DFtimedelta

    if dropna:
        DF = DF.dropna()

    tsout = time_series.TimeSeriePoint()

    T = conv.dt2posix(DF["T"].values)
    A = DF[lbda_colname("East") + unit_suffix].values
    B = DF[lbda_colname("North") + unit_suffix].values
    C = DF[lbda_colname("Up") + unit_suffix].values

    tsout.from_list(T, A, B, C, coortype="UTM")

    return tsout

def read_spotgins_masterfile(master_inp):
    return pd.read_fwf(master_inp, infer_nrows=11, skiprows=5, comment="#")

