#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.files_rw deals with file format conversion
between different file standards .

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
import itertools

#### Import the logger
import logging
import os
import re
import shutil
import textwrap
from io import BytesIO, StringIO

import dateutil
import numpy as np
import pandas as pd

#### geodeZYX modules
from geodezyx import conv
from geodezyx import operational
from geodezyx import utils
from geodezyx.files_rw import read_logsheets

log = logging.getLogger("geodezyx")

##########  END IMPORT  ##########


# _____          __  __ _____ _______                                   _
# / ____|   /\   |  \/  |_   _|__   __|                                 | |
# | |  __   /  \  | \  / | | |    | |      __ _  ___ _ __   ___ _ __ __ _| |
# | | |_ | / /\ \ | |\/| | | |    | |     / _` |/ _ \ '_ \ / _ \ '__/ _` | |
# | |__| |/ ____ \| |  | |_| |_   | |    | (_| |  __/ | | |  __/ | | (_| | |
# \_____/_/    \_\_|  |_|_____|  |_|     \__, |\___|_| |_|\___|_|  \__,_|_|
# __/ |
# |___/
#


def list_stat_in_statinfo(statinfoin):
    print(
        "WARN : DISCONTINUED : list_stat_in_statinfo doesnt work well, should be replaced by stat_list_in_station_info !!!!"
    )
    listtemp = []
    for l in open(statinfoin):
        if l[0] != " ":
            continue
        if l[1:5] == "\n":
            continue
    listtemp.append(l[1:5])
    listout = utils.uniqify_list(listtemp)
    return listout


def read_station_info_solo(filein, stat, column_type="ulr"):
    """
    column_type : str
        "ulr" or "sopac"
    """
    dicout = {
        k: []
        for k in [
            "Start",
            "End",
            "Rec",
            "Ant",
            "AntHt",
            "Dome",
            "Ant N",
            "Ant E",
            "Vers",
            "SwVers",
            "Rec SN",
            "Ant SN",
        ]
    }

    if column_type == "sopac":
        rg_start_1 = slice(25, 29)
        rg_start_2 = slice(30, 33)
        rg_start_3 = slice(34, 36)
        rg_start_4 = slice(37, 39)
        rg_start_5 = slice(40, 42)

        rg_end_1 = slice(44, 48)
        rg_end_2 = slice(49, 52)
        rg_end_3 = slice(53, 55)
        rg_end_4 = slice(56, 58)
        rg_end_5 = slice(59, 61)

        rg_rec = slice(96, 114)
        rg_antht = slice(64, 70)
        rg_ant = slice(170, 187)
        rg_dome = slice(187, 191)
        rg_antN = slice(80, 86)
        rg_antE = slice(89, 95)

        rg_vers = slice(119, 141)
        rg_sw_vers = slice(141, 148)
        rg_rec_sn = slice(148, 170)
        rg_ant_sn = slice(170, 187)

    elif column_type == "ulr":
        rg_start_1 = slice(25, 29)
        rg_start_2 = slice(30, 33)
        rg_start_3 = slice(34, 36)
        rg_start_4 = slice(37, 39)
        rg_start_5 = slice(40, 42)

        rg_end_1 = slice(44, 48)
        rg_end_2 = slice(49, 52)
        rg_end_3 = slice(53, 55)
        rg_end_4 = slice(56, 58)
        rg_end_5 = slice(59, 61)

        rg_rec = slice(97, 103)
        rg_antht = slice(64, 70)
        rg_ant = slice(134, 140)
        rg_dome = slice(142, 146)
        rg_antN = slice(80, 86)
        rg_antE = slice(89, 95)

    else:
        print("ERR : check column_type")
        return None

    for l in open(filein):
        if l[1:5] == stat:
            # cas specifique du temps (repris d'une autre fct)
            f = l
            start = conv.doy2dt(
                int(f[rg_start_1]),
                int(f[rg_start_2]),
                int(f[rg_start_3]),
                int(f[rg_start_4]),
                int(f[rg_start_5]),
            )
            dicout["Start"].append(start)
            if int(f[rg_end_1]) == 9999 or int(f[rg_end_2]) == 999:
                end = conv.doy2dt(2099, 1, 0, 0, 0)
            else:
                end = conv.doy2dt(
                    int(f[rg_end_1]),
                    int(f[rg_end_2]),
                    int(f[rg_end_3]),
                    int(f[rg_end_4]),
                    int(f[rg_end_5]),
                )
            dicout["End"].append(end)
            # Autres Donnees importantes
            dicout["Rec"].append(l[rg_rec].strip())
            dicout["AntHt"].append(float(l[rg_antht].strip()))
            dicout["Ant"].append(l[rg_ant].strip())
            dicout["Dome"].append(l[rg_dome].strip())
            dicout["Ant N"].append(l[rg_antN].strip())
            dicout["Ant E"].append(l[rg_antE].strip())

            if column_type == "sopac":
                dicout["Vers"].append(l[rg_vers].strip())
                dicout["SwVers"].append(l[rg_sw_vers].strip())
                dicout["Rec SN"].append(l[rg_rec_sn].strip())
                dicout["Ant SN"].append(l[rg_ant_sn].strip())

    return dicout


def read_station_info_solo_date(filein, stat, date, column_type="ulr"):
    DIC = read_station_info_solo(filein, stat, column_type=column_type)

    Start = DIC["Start"]
    End = DIC["End"]

    i_good = None
    for i, (s, e) in enumerate(zip(Start, End)):
        if (s <= date) and (date <= e):
            if i_good:
                print(
                    "WARN : index of the good period has been already defined, something is weird ..."
                )
            i_good = i

    if i_good is None:
        print(
            "ERR : no corresponding metadata have been found in the station.info for :"
        )
        log.info(stat, date)
        return None

    dicout = dict()

    for k, v in DIC.items():
        dicout[k] = v[i_good]

    return dicout


# dic = read_station_info_solo('/home/pierre/Documents/stationinfo2gins/station.info.ovsg','HOUE')


def read_lfile_solo(filein, stat):
    """
    WEAK : should be improved with discontinuites management
    """
    for l in open(filein):
        X, Y, Z = 0.0, 0.0, 0.0
        vX, vY, vZ = 0.0, 0.0, 0.0

        T = dt.datetime.now()

        if l[1:5] == stat:
            f = l[1:].split()

            X = float(f[1])
            Y = float(f[2])
            Z = float(f[3])
            Ttmp = float(f[7])

            if np.isclose(Ttmp, 0.0):
                T = conv.year_decimal2dt(2000.0)
            else:
                T = conv.year_decimal2dt(Ttmp)

            if len(l) > 225:  # if velocities are given (125 is arbitrary)
                vX = float(f[8])
                vY = float(f[9])
                vZ = float(f[10])
            else:  # if no velocityis given
                vX = 0.0
                vY = 0.0
                vZ = 0.0

            return X, Y, Z, T, vX, vY, vZ

    log.warning("no coords found for : %s ")
    return None


def read_pbo_vel_file_solo(velfilein, stat):
    """
    return a LIST of panda dataframe (corrsponding to different discontinuites),
    where the keywords of the pbo .vel can be used directly

    exemple :
    ftp://data-out.unavco.org/pub/products/velocity/pbo.final_nam08.vel
    """
    fil = open(velfilein)
    skiplinelis = []
    for i, l in enumerate(fil):
        if l[0] not in (" ", "*"):
            skiplinelis.append(i)

    data = pd.read_table(
        open(velfilein), skiprows=skiplinelis, sep=" *", header=max(skiplinelis) + 1
    )
    # data['Ref_epoch']

    i_lis = []
    for i, s in enumerate(list(data["*Dot#"])):
        if s == stat.upper():
            i_lis.append(i)

    data2 = data.transpose()

    data_lis = []
    for i in i_lis:
        data_lis.append(data2[i])
    return data_lis


def read_globk_vel_file(velfile_in):

    p = velfile_in
    i = 0
    for l in open(p):
        i += 1
        if "(deg)" in l:
            break

    D = pd.read_table(p, skiprows=i, header=-1, delim_whitespace=True)

    D.rename(
        columns={
            0: "Long",
            1: "Lat",
            2: "E",
            3: "N",
            4: "E_adj",
            5: "N_adj",
            6: "E_sig",
            7: "N_sig",
            8: "Rho",
            9: "H",
            10: "H_adj",
            11: "H_sig",
            12: "SITE",
        },
        inplace=True,
    )

    return D


def read_sinex_discontinuity_solo(snxfile, stat, PorV="p"):

    flagxyz = False

    soln_lis = []
    start_lis = []
    end_lis = []
    PorV_lis = []
    notes_lis = []

    for line in open(snxfile):

        if re.compile("SOLUTION/DISCONTINUITY").search(line):
            flagxyz = not flagxyz
            continue

        if line[0] != " ":
            continue

        if flagxyz == True:
            fields = line.split()

            if line[0] != " ":
                continue

            if fields[0] == stat.upper():

                if not fields[6] == PorV:
                    continue

                soln_lis.append(int(fields[2]))
                start_lis.append(conv.datestr_sinex_2_dt(fields[4]))
                end = conv.datestr_sinex_2_dt(fields[5])
                if end == dt.datetime(1970, 1, 1):
                    end = dt.datetime(2099, 1, 1)
                end_lis.append(end)
                PorV_lis.append(fields[6])
                notes_lis.append(line.split("-")[-1].strip())

    return start_lis, end_lis, notes_lis, PorV_lis


def stat_list_in_station_info(filein):
    statlist = []
    fstatinfo = open(filein, "r+")
    for line in fstatinfo:
        if len(line) == 0:
            log.warning("empty line in the station.info !!!")
            continue
        elif line[0] != " ":
            continue
        elif " " == line[0]:
            f = line.split()
            statlist.append(f[0])
    out = sorted(list(set(statlist)))
    return out


def convert_statinfo2eqfile(statinfoin, eqfileout):
    outfileobj = open(eqfileout, "w+")
    stat_list = stat_list_in_station_info(statinfoin)
    for stat in stat_list:
        i = 0
        startlis, endlis = read_station_info_time_solo(statinfoin, stat)
        statbisproto = stat + "_GPS"
        for s, e in zip(startlis, endlis):
            if i != 0:
                statbis_list = list(statbisproto)
                if i > 9:
                    statbis_list[5] = str(i)[0]
                    statbis_list[6] = str(i)[1]
                else:
                    statbis_list[5] = str(i)
                statbis = "".join(statbis_list)
            else:
                statbis = statbisproto
            sy = s.year
            smo = s.month
            sd = s.day
            sh = s.hour
            smi = s.minute
            ss = s.second
            ey = e.year
            emo = e.month
            ed = e.day
            eh = e.hour
            emi = e.minute
            es = e.second
            ligne = " rename {}     {} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2}\n".format(
                stat, statbis, sy, smo, sd, sh, smi, ey, emo, ed, eh, emi
            )
            outfileobj.write(ligne)
            i = i + 1
    outfileobj.close()
    return None


def read_station_info_time_solo(filein, stat):
    """BROUILLON
    pour une station, donne les startdate et enddate d'une periode
    SORTIE : 2 x liste de dt"""

    startlis = []
    endlis = []

    for l in open(filein):

        if l[1:5] == stat:

            f = l[25:].split()

            start = conv.doy2dt(int(f[0]), int(f[1]), int(f[2]), int(f[3]), int(f[4]))
            startlis.append(start)

            if int(f[5]) == 9999 or int(f[6]) == 999:
                end = conv.doy2dt(2099, 1, 0, 0, 0)
            else:
                end = conv.doy2dt(int(f[5]), int(f[6]), int(f[7]), int(f[8]), int(f[9]))
            endlis.append(end)

    return startlis, endlis


def read_eqfile_time_solo(filein, stat):
    """BROUILLON
    pour une station, donne les startdate et enddate d'une periode
    SORTIE : 2 x liste de dt"""

    startlis = []
    endlis = []

    for l in open(filein):
        f = l.split()
        if f[1] == stat:

            start = dt.datetime(int(f[3]), int(f[4]), int(f[5]), int(f[6]), int(f[7]))
            startlis.append(start)

            if int(f[7]) == 9999 or int(f[8]) == 999:
                end = dt.datetime(2099, 1, 1, 0, 0)
            else:
                end = dt.datetime(
                    int(f[8]), int(f[9]), int(f[10]), int(f[11]), int(f[12])
                )
            endlis.append(end)

    return startlis, endlis


def read_eqfile_as_dico(filein, delGPSfield=False):
    """
    for an eqfile return a dico
    dic[STAT] = [disc1 , disc2 ...] where disc are
    the discont of the FIRST column

    delGPSfield manage the elimination of the line of type _GPS
    which are not regular
    EDIT : 150817 WHY ???

    """

    startlis = []
    endlis = []

    outdico = dict()

    for l in open(filein):
        f = l.split()
        stat = f[1]
        if stat not in outdico:
            outdico[stat] = []
        if f[2][-3] == "G" and delGPSfield:
            continue
        datelis = [f[3], f[4], f[5], f[6], f[7], f[8]]
        datelis = [int(e) for e in datelis]
        outdico[stat].append(dt.datetime(*datelis))

    return outdico


def write_eqfile_from_dico(dicoin, outdir, outfilename):
    outfilepath = os.path.join(outdir, outfilename)
    fil = open(outfilepath, "w+")

    #    dicoin = OrderedDict(sorted(dicoin.items(), key=lambda dicoin: dicoin[1]))
    #    dicoin = sorted(dicoin, key=dicoin.get)
    for k, v in dicoin.items():
        for i in range(len(v)):
            stat = k
            if i + 2 > 9:
                statbis = stat + "_" + str(i + 2) + "S"
            else:
                statbis = stat + "_" + str(i + 2) + "PS"

            s = v[i]
            if i == len(v) - 1:
                e = dt.datetime(2099, 1, 1)
            else:
                e = v[i + 1]

            sy = s.year
            smo = s.month
            sd = s.day
            sh = s.hour
            smi = s.minute
            ss = s.second
            ey = e.year
            emo = e.month
            ed = e.day
            eh = e.hour
            emi = e.minute
            es = e.second

            ligne = " rename {}     {} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2}\n".format(
                stat, statbis, sy, smo, sd, sh, smi, ey, emo, ed, eh, emi
            )
            fil.write(ligne)
    return None


#      _        _   _             _        __        ___     _____       _______ _____
#     | |      | | (_)           (_)      / _|      |__ \   / ____|   /\|__   __/ ____|
#  ___| |_ __ _| |_ _  ___  _ __  _ _ __ | |_ ___      ) | | |       /  \  | | | (___
# / __| __/ _` | __| |/ _ \| '_ \| | '_ \|  _/ _ \    / /  | |      / /\ \ | |  \___ \
# \__ \ || (_| | |_| | (_) | | | | | | | | || (_) |  / /_  | |____ / ____ \| |  ____) |
# |___/\__\__,_|\__|_|\___/|_| |_|_|_| |_|_| \___/  |____|  \_____/_/    \_\_| |_____/
#


def statname_of_catsfile(catsneu_path):
    catsneu_f = open(catsneu_path, "r")
    for line in catsneu_f:
        if "Site" in line:
            return line.split(" ")[2]
    return "XXXX"


def statinfo_2_cats(statinfo_path, catsneu_path):
    stat = statname_of_catsfile(catsneu_path)
    log.info("stat = ", stat)
    if stat == "ABMF":
        return None
    stat_dic = read_station_info_solo(statinfo_path, stat)
    start = stat_dic["Start"]
    antht = stat_dic["AntHt"]

    outcats_path = (
        os.path.dirname(catsneu_path) + os.path.basename(catsneu_path) + ".offset"
    )
    #    shutil.copyfile(catsneu_path,outcats_path)

    antht = [a for (s, a) in sorted(zip(start, antht))]
    start = sorted(start)
    log.info(stat_dic, stat)
    s0 = start[0]
    h0 = antht[0]
    s_stk = []
    for s, h in zip(start, antht):
        if h != h0:
            s_stk.append(s)
            h0 = h
    log.info(outcats_path)
    catsneuout_f = open(outcats_path, "w")
    catsneu_f = open(catsneu_path, "a+")
    lastline = ""
    for line in catsneu_f:
        if lastline == line:
            continue
        elif "Height" in line:
            catsneuout_f.write(line)
            for s in s_stk:
                catsneuout_f.write(
                    "# offsets : " + str(conv.dt2year_decimal(s)) + " 1\n"
                )
        # cleaning outlier
        elif line[0] != "#":
            if abs(float(line.split()[1])) > 1 or abs(float(line.split()[2])) > 1:
                pass
            else:
                catsneuout_f.write(line)
        else:
            catsneuout_f.write(line)
        lastline = line

    catsneuout_f.close()
    # type 1 for only vertical (cf doc)
    return None


#      _        _   _             _        __        ___          _
#     | |      | | (_)           (_)      / _|      |__ \        (_)
#  ___| |_ __ _| |_ _  ___  _ __  _ _ __ | |_ ___      ) |   __ _ _ _ __  ___
# / __| __/ _` | __| |/ _ \| '_ \| | '_ \|  _/ _ \    / /   / _` | | '_ \/ __|
# \__ \ || (_| | |_| | (_) | | | | | | | | || (_) |  / /_  | (_| | | | | \__ \
# |___/\__\__,_|\__|_|\___/|_| |_|_|_| |_|_| \___/  |____|  \__, |_|_| |_|___/
#                                                            __/ |
#                                                           |___/

# def oceload_associed(statginsfile):
#    import os
#
#    aaa = os.path.abspath(statginsfile)
#    bbb = 'outocload'
#    filearg = open('temparg','wr')
#    filearg.write(aaa + '\n')
#    filearg.write(os.path.basename(aaa) + '/' + bbb)
#
#    os.system('bash /home/psakicki/scripts/exe_loadoce < temparg')
#


def receptor_gins_corrector(inRec):
    outRec = inRec
    outRec = outRec.replace("ASHTECH", "AS")
    outRec = outRec.replace("TRIMBLE", "TR")
    outRec = outRec.replace("LEICA GRX1200+GNS", "LEICA GRX12")
    outRec = outRec.replace("LEICA GRX1200GGPRO", "LEICA GRX12")
    outRec = outRec.replace("ROGUE", "RG")
    return outRec


def station_info_2_gins(
    statinfoin,
    coordfilein,
    outfile,
    coordfile_type="pbovelfile",
    specific_stats_lis=[],
    ellipsoid="GRS80",
    station_info_columns_type="sopac",
):
    """
    Convert a GAMIT station.info to a GINS stations files

    Parameters
    ----------
    statinfoin : str
        path of input station.info.
    coordfilein : str
        path of input coordinates files.
    outfile : str
         path of output GINS station file.
    coordfile_type : str, optional
        GAMIT coordinates files type : 'lfile' or 'pbovelfile'.
        The default is 'pbovelfile'.
    specific_stats_lis : list, optional
         list of specific stations. The default is [].
    ellipsoid : str, optional
        ellipsoid. The default is 'GRS80'.
    station_info_columns_type : TYPE, optional
        it exists different subtypes of station.info ...
        handles "sopac" or "ulr" yet. The default is "sopac".


    Returns
    -------
    outfile : str
         path of GINS outfile.

    """

    # Info de referentiel de la 1ère ligne
    refiel = "i08i08"
    # Info de referentiel des lignes de réoccupations
    refchange = "igs10"
    # code 'GRG3'
    grg = "GRG3"
    # sigma sur les offsets
    sigoffset = 0.000
    # sigma sur les positions
    sigposi = 0.030
    # sigma sur les vitesses
    sigvit = 0.005
    # Vitesses

    # ID DOMES du pays
    # ex : 971 pour la Guadeloupe
    # mettre 999 dans la plupart des cas
    # ATTENTION : string entre quotes
    # idprefix='971'
    # idprefix='990'
    # EDIT 1805, for very big station.info now idprfix is disabled
    # idprefix='990'

    # id initial, les stations auront un numero atribué de maniere décroissante
    # 97198 97197 97196 ...
    # idstatinit=99
    # EDIT 1805, for very big station.info now idstatinit is now at 99999
    idstatinit = 99999

    # ======== FIN DES MODIFIABLES ========

    outEnd = "000000"

    statgins_filobj = open(outfile, "w")

    if specific_stats_lis == []:
        listat = stat_list_in_station_info(statinfoin)
    else:
        listat = specific_stats_lis

    log.info(listat)

    sigvit = "{:.4f}".format(sigvit).lstrip("0")

    header = header_from_ellipsoid(ellipsoid)
    statgins_filobj.write(header)

    for stat in listat:

        log.info(stat)

        # On degresse l'ID quoiqu'il arrive après
        # pour conserver une homogenéité dans les pseudoDOMES
        idstatinit = idstatinit - 1

        if coordfile_type == "lfile":
            output_lfile = None
            output_lfile = read_lfile_solo(coordfilein, stat)

            if not output_lfile:
                log.info("WARN : station_info_2_gins : skip station ", stat)
                continue
            else:
                X, Y, Z, T, vX, vY, vZ = output_lfile

        elif coordfile_type == "pbovelfile":
            data = read_pbo_vel_file_solo(coordfilein, stat)
            if len(data) != 0:
                data = data[0]
                X = data["Ref_X"]
                Y = data["Ref_Y"]
                Z = data["Ref_Z"]
                T = conv.mjd2dt(data["Ref_jday"])
                vX = data["dX/dt"]
                vY = data["dY/dt"]
                vZ = data["dZ/dt"]
            else:
                log.warning("no %s found in %s", stat, coordfilein)
                continue
        else:
            raise Exception("Wrong coordfile_type")

        if X == 0:
            continue

        # EDIT 1805, for very big station.info now idprefix is disabled
        # idstat = str(idprefix) + str(idstatinit)
        idstat = str(idstatinit)

        T = T.strftime("%d%m%y")

        dicout = read_station_info_solo(
            statinfoin, stat, column_type=station_info_columns_type
        )

        n = len(dicout["Start"])
        log.info("N = %s", n)

        void = "             "
        if stat == "ABMF":
            idstat = str(97103)

        fortupsta = (
            idstat,
            "M001",
            stat,
            void,
            X,
            sigposi,
            Y,
            Z,
            vX,
            sigvit,
            vY,
            vZ,
            T,
            refiel,
        )

        outputline = "   {0} {1} {2}{3}{4:>-12.3f} ~{5:5.3f} {6:>-12.3f} ~{5:5.3f} {7:>-12.3f} ~{5:5.3f} {8:>-7.4f} ~{9} {10:>-7.4f} ~{9} {11:>-7.4f} ~{9} {12}        {13}".format(
            *fortupsta
        )

        outputline = outputline + "\n"
        statgins_filobj.write(outputline)

        for i in range(n):
            idligne = idstat + "{:0>2d}".format(i + 1)
            outAnt = dicout["Ant"][i]
            outAntHt = dicout["AntHt"][i]
            Xant, Yant, Zant = conv.enu2xyz_legacy(0, 0, outAntHt, X, Y, Z)
            lastEnd = outEnd
            outStart = dicout["Start"][i].strftime("%d%m%y")
            outEnd = dicout["End"][i].strftime("%d%m%y")

            if outStart == lastEnd:
                outStart = (dicout["Start"][i] + dt.timedelta(days=1)).strftime(
                    "%d%m%y"
                )

            if dicout["End"][i] == dt.datetime(2099, 1, 1, 0, 0, 0):
                outEnd = "      "

            outDome = dicout["Dome"][i]
            outRec = dicout["Rec"][i]

            if "---------------" in outRec:
                continue
                # methode curieuse dans le stat.info sopac : on note
                # l'installation de l'antenne sans le rec : useless

            outRec = receptor_gins_corrector(outRec)[:11]
            outAnt = outAnt[:11]
            void = "                     "

            formtuple = (
                idligne,
                grg,
                stat,
                outRec,
                Xant,
                sigoffset,
                Yant,
                sigoffset,
                Zant,
                sigoffset,
                outAnt,
                outDome,
                void,
                outStart,
                outEnd,
                refchange,
            )

            # outputline ="%s %4s %4s %20s %5.3f %s %5.3f %s %5.3f %s %15s %s %s %s %s" \
            # %(idligne,grg,stat,outRec,Xant,sigoffset,Yant,sigoffset,Zant,sigoffset,outAnt,outDome,outStart,outEnd,refchange)

            outputline = " {0} {1} {2} {3:<11} {4:>-12.3f}  {5:5.3f} {6:>-12.3f}  {7:5.3f} {8:>-12.3f}  {9:5.3f}   {10:15} {11:4} {12} {13} {14} {15}".format(
                *formtuple
            )

            outputline = outputline + "\n"
            statgins_filobj.write(outputline)

        statgins_filobj.write("\n")

    if station_info_columns_type == "ulr":
        log.info(
            "NOTA : input station.info is a 'ulr' type, the antenna/reciever type has to be corrected manually in the returned GINS file"
        )

    return outfile

    # oceload_associed(outfile)


def stat_file_gins_new_fmt(
    file_out_path,
    STAT="STAT",
    xyz=[0.0, 0.0, 0.0],
    ecc_une=[0.0, 0.0, 0.0],
    rec="NONE",
    ant="NONE",
    radom="NONE",
    ant_id="0000",
):

    f = open(file_out_path, "w+")

    ellipsoid = """#VERSION 1.0
#RADIUS    0.6378137000000000E+07
#1/FLAT    0.2982572221010000E+03
# ITRF2014"""

    x, y, z = xyz[0], xyz[1], xyz[2]
    ecc_u, ecc_n, ecc_e = ecc_une[0], ecc_une[1], ecc_une[2]

    l1 = "   99999 M001 {:4}                             mar_xyz {:15.6f} 0.000000 {:15.6f} 0.000000 {:15.6f} 0.000000 vit_xyz  0.00000 0.0000  0.00000 0.0000  0.00000 0.0000  20100101_000000   EURA          i14i14"
    l1_fmt = l1.format(STAT, x, y, z)
    l2 = " 9999999 GRG3 {:4}      {:23}ecc_une {:15.6f} 0.000000 {:15.6f} 0.000000 {:15.6f} 0.000000            {:16}{:4} {:25}20000101_000000                 rinex"
    l2_fmt = l2.format(STAT, rec, ecc_u, ecc_n, ecc_e, ant, radom, ant_id)

    out = "\n".join((ellipsoid, l1_fmt, l2_fmt))

    f.write(out)
    f.close()
    return file_out_path


def write_station_info_from_datalists(
    period_lis_lis, site_lis, location_lis, station_info_out_path
):
    """


    Parameters
    ----------
    period_lis_lis : list of lists
        list of lists of periods.
    site_lis : list
        list of sites.
    location_lis : list
        list of locations.
    station_info_out_path : str
        the full path (directory+name) of the generated station.info.

    Returns
    -------
    None.

    Note
    ----
    The different data lists (i.e. period_lis_lis , stat_lis , loc_lis)
    are produced by:
        * mono_logsheet_read
              do not forget to activate return_lists = True

        * multi_logsheet_read
              do not forget to activate return_dico = False
    """

    si_file = open(station_info_out_path, "w")
    si_file.write(
        "*SITE  Station Name      Session Start      Session Stop       Ant Ht   HtCod  Ant N    Ant E    Receiver Type         Vers                  SwVer  Receiver SN           Antenna Type     Dome   Antenna SN          \n"
    )
    # writing the period in a station info
    for period_lis, sit, loc in zip(period_lis_lis, site_lis, location_lis):
        for period in period_lis:
            d1, d2, a, r = period
            # simplifiy name, removing space and non ascii char. ,
            # and fill with spaces to have 16 chars
            name = loc.City_or_Town[0:16].replace(" ", "_")
            name = "".join([i if ord(i) < 128 else "_" for i in name])
            if len(name) != 16:
                name = name.ljust(16)
            stat = sit.Four_Character_ID
            # managing the 4 chars
            stat = stat[0:4]
            if len(stat) != 4:
                stat = stat.ljust(4)

            strtup = (
                stat,
                name,
                d1.year,
                conv.dt2doy(d1),
                int(d1.hour),
                int(d1.minute),
                int(d1.second),
                d2.year,
                conv.dt2doy(d2),
                int(d2.hour),
                int(d2.minute),
                int(d2.second),
                float(a.Up_Ecc),
                a.ARPSmart(),
                float(a.North_Ecc),
                float(a.East_Ecc),
                r.Receiver_Type,
                r.Firmware_Version,
                r.FirmwareSmart(),
                str(r.Serial_Number)[:20],
                a.AntTypSmart(),
                a.Antenna_Radome_Type,
                str(a.Serial_Number)[:20],
            )

            outline = " {}  {:16}  {} {} {:0>2d} {:0>2d} {:0>2d}  {} {} {:0>2d} {:0>2d} {:0>2d}  {:+.4f}  {}  {:+.4f}  {:+.4f}  {:20}  {:<20}  {:>5.2f}  {:<20}  {:15}  {:5}  {}\n".format(
                *strtup
            )
            si_file.write(outline)
    si_file.close()
    return None


def write_lfile_from_datalists(site_lis, location_lis, lfile_out_path):
    """datalists (period_lis_lis , stat_lis , loc_lis  ) are produced by :
    * mono_logsheet_read
    * multi_logsheet_read"""
    lf_file = open(lfile_out_path, "w")

    proto_str = " {}_GPS {:.4f} {:.4f} {:.4f} {:+6.4f} {:+6.4f} {:+6.4f} {:6.1f} {:5.3f} {:5.3f} {:5.3f} {:6.4f} {:6.4f} {:6.4f}\n"

    for sit, loc in zip(site_lis, location_lis):
        vx, vy, vz = 0, 0, 0
        sx, sy, sz = 0, 0, 0
        svx, svy, svz = 0, 0, 0

        yr = conv.dt2year_decimal(loc.Reference_epoch)

        outline = proto_str.format(
            sit.Four_Character_ID,
            loc.X_coordinate_m,
            loc.Y_coordinate_m,
            loc.Z_coordinate_m,
            vx,
            vy,
            vz,
            yr,
            sx,
            sy,
            sz,
            svx,
            svy,
            svz,
        )

        lf_file.write(outline)
    lf_file.close()
    return None


def header_from_ellipsoid(ellipsoid):
    if ellipsoid in ("GRS80", "WGS84"):
        header = textwrap.dedent(
            """APRES ITRF2008 TOUTES STATIONS       ,  R= 6378137.00,   f= 298.257222101
 station type site                    X (m)   sig.        Y (m)   sig.        Z(m)    sig.    Xp (m/an)      Yp (m/an)      Zp (m/an)    date  plaque source
   -9999 M000 Earth c. of mass        0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                              010100
   -1998 EGB  Bias       (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                              010198 311298
    -199 EGB  Bias       (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                              010199 310199
     -61 EGC1 1/year Cos (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                                            dy ecc
     -51 EGS1 1/year Sin (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                                            dy ecc
     -62 EGC2 2/year Cos (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                                            dy ecc
     -52 EGS2 2/year Sin (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                                            dy ecc

   14001 M004 ZIMM              4331297.063 ~0.001   567555.878 ~0.001  4633133.936 ~0.001 -0.0135 ~.0001  0.0181 ~.0000  0.0121 ~.0001 010105        i08i08
 1400101 GRG3 ZIMM TR 4000SSE         0.000  0.000        0.000  0.000        0.000  0.000   TRM14532.00     NONE 3311A                 010593 050897 igs06
 1400102 GRG3 ZIMM TR 4000SSI         0.000  0.000        0.000  0.000        0.000  0.000   TRM14532.00     NONE 3311A                 060897 051198 igs06

   14001 M004 ZIMM              4331297.062 ~0.001   567555.882 ~0.001  4633133.935 ~0.001 -0.0135 ~.0001  0.0181 ~.0000  0.0121 ~.0001 010105        i08i08
 1400106 GRG3 ZIMM TR 4000SSI         0.000  0.000        0.000  0.000        0.000  0.000   TRM29659.00     NONE 99390                 061198 110803 igs06
 1400107 GRG3 ZIMM TR 4700            0.000  0.000        0.000  0.000        0.000  0.000   TRM29659.00     NONE 99390                 120803 210206 igs06
 1400108 GRG3 ZIMM TR NETRS           0.000  0.000        0.000  0.000        0.000  0.000   TRM29659.00     NONE 99390                 220206        igs06

            """
        )
    elif ellipsoid in ("GRIM", "EIGEN"):
        header = textwrap.dedent(
            """APRES ITRF2008 TOUTES STATIONS       ,  R= 6378136.46,   f= 298.257650
 station type site                    X (m)   sig.        Y (m)   sig.        Z(m)    sig.    Xp (m/an)      Yp (m/an)      Zp (m/an)    date  plaque source
   -9999 M000 Earth c. of mass        0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                              010100
   -1998 EGB  Bias       (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                              010198 311298
    -199 EGB  Bias       (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                              010199 310199
     -61 EGC1 1/year Cos (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                                            dy ecc
     -51 EGS1 1/year Sin (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                                            dy ecc
     -62 EGC2 2/year Cos (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                                            dy ecc
     -52 EGS2 2/year Sin (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                                            dy ecc

   14001 M004 ZIMM              4331297.063 ~0.001   567555.878 ~0.001  4633133.936 ~0.001 -0.0135 ~.0001  0.0181 ~.0000  0.0121 ~.0001 010105        i08i08
 1400101 GRG3 ZIMM TR 4000SSE         0.000  0.000        0.000  0.000        0.000  0.000   TRM14532.00     NONE 3311A                 010593 050897 igs06
 1400102 GRG3 ZIMM TR 4000SSI         0.000  0.000        0.000  0.000        0.000  0.000   TRM14532.00     NONE 3311A                 060897 051198 igs06

   14001 M004 ZIMM              4331297.062 ~0.001   567555.882 ~0.001  4633133.935 ~0.001 -0.0135 ~.0001  0.0181 ~.0000  0.0121 ~.0001 010105        i08i08
 1400106 GRG3 ZIMM TR 4000SSI         0.000  0.000        0.000  0.000        0.000  0.000   TRM29659.00     NONE 99390                 061198 110803 igs06
 1400107 GRG3 ZIMM TR 4700            0.000  0.000        0.000  0.000        0.000  0.000   TRM29659.00     NONE 99390                 120803 210206 igs06
 1400108 GRG3 ZIMM TR NETRS           0.000  0.000        0.000  0.000        0.000  0.000   TRM29659.00     NONE 99390                 220206        igs06

            """
        )
    else:
        log.error("Check ellipsoid Name")
        return None
    return header


def write_station_file_gins_from_datalists(
    period_lis_lis, site_lis, location_lis, station_info_out_path, ellipsoid="GRS80"
):
    """datalists (period_lis_lis , stat_lis , loc_lis  ) are produced by :
    * mono_logsheet_read
    * multi_logsheet_read"""

    # Info de referentiel de la 1ère ligne
    refiel = "i08i08"
    # Info de referentiel des lignes de réoccupations
    refchange = "igs10"
    # code 'GRG3'
    grg = "GRG3"
    # sigma sur les offsets
    sigoffset = 0.000
    # sigma sur les positions
    sigposi = 0.030
    # sigma sur les vitesses
    sigvit = 0.005
    # Vitesses

    # ID DOMES du pays
    # ex : 971 pour la Guadeloupe
    # mettre 999 dans la plupart des cas
    # ATTENTION : string entre quotes
    idprefix = "971"
    idprefix = "990"

    # id initial, les stations auront un numero atribué de maniere décroissante
    # 97198 97197 97196 ...
    idstatinit = 99

    void = "             "
    void2 = "                     "

    statgins_filobj = open(station_info_out_path, "w+")

    header = header_from_ellipsoid(ellipsoid)
    statgins_filobj.write(header)

    i = 0
    for i, (period_lis, site, loca) in enumerate(
        zip(period_lis_lis, site_lis, location_lis)
    ):
        T = loca.Reference_epoch.strftime("%d%m%y")
        idstatA = site.IERS_DOMES_Number[0:5]
        idstatB = site.IERS_DOMES_Number[5:9]

        # case where an auto DOMES must be given
        if int(idstatA) == 99999:
            idstatA = str(int(idstatA) - i)

        stat = site.Four_Character_ID
        # managing the 4 chars
        stat = stat[0:4]
        if len(stat) != 4:
            stat = stat.ljust(4)

        X = float(loca.X_coordinate_m)
        Y = float(loca.Y_coordinate_m)
        Z = float(loca.Z_coordinate_m)

        vX = float(loca.X_velocity)
        vY = float(loca.Y_velocity)
        vZ = float(loca.Z_velocity)

        sigposi = float(loca.X_coordinate_sigma)
        sigvit = float(loca.X_velocity_sigma)

        fortupsta = (
            idstatA,
            idstatB,
            stat,
            void,
            X,
            sigposi,
            Y,
            Z,
            vX,
            sigvit,
            vY,
            vZ,
            T,
            refiel,
        )
        outputline = "   {0} {1} {2}{3}{4:>-12.3f} ~{5:5.3f} {6:>-12.3f} ~{5:5.3f} {7:>-12.3f} ~{5:5.3f} {8:>-7.4f} ~{9} {10:>-7.4f} ~{9} {11:>-7.4f} ~{9} {12}        {13}".format(
            *fortupsta
        )

        outputline = outputline + "\n"
        statgins_filobj.write(outputline)

        for j, peri in enumerate(period_lis):
            if j == 0:
                d_new_peri = 0
            else:
                d_new_peri = 1
            start, end, ant, rec = peri
            idligne = idstatA + str(j).zfill(2)
            outRec = receptor_gins_corrector(rec.Receiver_Type)

            Xant, Yant, Zant = conv.enu2xyz_legacy(
                ant.East_Ecc, ant.North_Ecc, ant.Up_Ecc, X, Y, Z
            )

            sigoffset = 0.0
            outAnt = ant.AntTypSmart()
            outDome = ant.Antenna_Radome_Type[0:4].ljust(4)
            outStart = (start + dt.timedelta(days=d_new_peri)).strftime("%d%m%y")
            if end == dt.datetime(2099, 1, 1, tzinfo=dateutil.tz.tzutc()):
                outEnd = "      "
            else:
                outEnd = end.strftime("%d%m%y")

            formtuple = (
                idligne,
                grg,
                stat,
                outRec,
                Xant,
                sigoffset,
                Yant,
                sigoffset,
                Zant,
                sigoffset,
                outAnt,
                outDome,
                void2,
                outStart,
                outEnd,
                refchange,
            )

            outputline = " {0} {1} {2} {3:<11} {4:>-12.3f}  {5:5.3f} {6:>-12.3f}  {7:5.3f} {8:>-12.3f}  {9:5.3f}   {10:15} {11} {12} {13} {14} {15}".format(
                *formtuple
            )
            outputline = outputline + "\n"
            statgins_filobj.write(outputline)

        statgins_filobj.write("\n")
    return None


def smart_elt_list(list_raw, n_elt, replacement=""):
    try:
        return list_raw[n_elt]
    except:
        return replacement


def read_rinex_2_dataobjts(rinex_path):

    if utils.empty_file_check(rinex_path):
        log.error("the RINEX file is empty ...")
        log.error(rinex_path)

        return None, None, None, None

    ant_raw = utils.grep(rinex_path, "ANT #", True)
    rec_raw = utils.grep(rinex_path, "REC #", True)
    xyz_raw = utils.grep(rinex_path, "APPROX POSITION XYZ", True).split()
    stat_raw = utils.grep(rinex_path, "MARKER NAME", True).split()
    domes_raw = utils.grep(rinex_path, "MARKER NUMBER", True).split()
    d_hen_raw = utils.grep(rinex_path, "ANTENNA: DELTA H/E/N", True).split()
    t_raw = utils.grep(rinex_path, "TIME OF FIRST OBS", True).split()

    Antobj, Recobj, Siteobj, Locobj = (
        read_logsheets.Antenna(),
        read_logsheets.Receiver(),
        read_logsheets.Site(),
        read_logsheets.Location(),
    )

    Locobj.X_coordinate_m = float(smart_elt_list(xyz_raw, 0))
    Locobj.Y_coordinate_m = float(smart_elt_list(xyz_raw, 1))
    Locobj.Z_coordinate_m = float(smart_elt_list(xyz_raw, 2))

    Antobj.Up_Ecc = float(smart_elt_list(d_hen_raw, 0))
    Antobj.East_Ecc = float(smart_elt_list(d_hen_raw, 1))
    Antobj.North_Ecc = float(smart_elt_list(d_hen_raw, 2))

    Siteobj.Four_Character_ID = smart_elt_list(stat_raw, 0)

    if re.search("[0-9]{5}[A-Z][0-9]{3}", smart_elt_list(domes_raw, 0)):
        Siteobj.IERS_DOMES_Number = smart_elt_list(domes_raw, 0)
    else:
        Siteobj.IERS_DOMES_Number = (
            str(np.random.randint(99999)).zfill(5) + "M001"
        )  # smart_elt_list(domes_raw,0) A CHANGER !!!!!!!

    Antobj.Antenna_Radome_Type = ant_raw[36:40].strip()
    Antobj.Antenna_Type = ant_raw[20:36].strip()

    Recobj.Receiver_Type = rec_raw[20:40].strip()

    # ======== OBSOLETE car utilisation de la fct rinex_start_end ========
    ##securité pour les secondes
    # if not (0 <= t_raw[5] < 60):
    # if t_raw[5] > 30:
    # t_raw[5] = 59
    # else:
    # t_raw[5] = 0

    # Date_Installed = dt.datetime(*[int(float(e)) for e in t_raw[0:6]],
    # tzinfo=dateutil.tz.tzutc())
    # Date_Removed = Date_Installed + dt.timedelta(days=1)

    date_installed, date_removed = operational.rinex_start_end(
        rinex_path, add_tzinfo=1, verbose=0
    )

    Antobj.Date_Removed = date_removed
    Antobj.Date_Installed = date_installed

    Recobj.Date_Removed = date_removed
    Recobj.Date_Installed = date_installed

    Locobj.Reference_epoch = date_installed

    return Antobj, Recobj, Siteobj, Locobj


def write_station_file_gins_from_rinex(
    rinex_path, station_file_out, stat_code_filename_prioritary=True
):
    antobj, recobj, siteobj, locobj = read_rinex_2_dataobjts(rinex_path)

    rnx_name = os.path.basename(rinex_path)
    stat_code_filename = os.path.basename(rinex_path)[0:4]

    if siteobj.Four_Character_ID.upper() != stat_code_filename.upper():
        print("WARN : different 4-char. code in RINEX header and filename")
        print("       for", rnx_name, "(", siteobj.Four_Character_ID, ")")
        if stat_code_filename_prioritary:
            print("       keeping the 4-char. code in filename")
            siteobj.Four_Character_ID = stat_code_filename.upper()

    write_station_file_gins_from_datalists(
        [
            [
                (
                    antobj.Date_Installed,
                    antobj.Date_Removed + dt.timedelta(days=1),
                    antobj,
                    recobj,
                )
            ]
        ],
        [siteobj],
        [locobj],
        station_file_out,
    )

    return None


#                     _   _   _ __  __ ______             _____  _____
#                    | | | \ | |  \/  |  ____|   /\      / ____|/ ____|   /\
#  _ __ ___  __ _  __| | |  \| | \  / | |__     /  \    | |  __| |  __   /  \
# | '__/ _ \/ _` |/ _` | | . ` | |\/| |  __|   / /\ \   | | |_ | | |_ | / /\ \
# | | |  __/ (_| | (_| | | |\  | |  | | |____ / ____ \  | |__| | |__| |/ ____ \
# |_|  \___|\__,_|\__,_| |_| \_|_|  |_|______/_/    \_\  \_____|\_____/_/    \_\
#


def read_nmea(
    file_path,
    outtype="FLH",
    df_out=True,
    use_altitude=False,
    startdate=dt.datetime(1980, 1, 1),
    enddate=dt.datetime(2099, 1, 1),
    export_path="",
):
    """if export_path != '', export a Matlab readable file to this path
    WARNING !!! la coord de ref est codée en dur, à coriger !!!!"""
    T = []
    Lat = []
    Long = []
    Haut = []
    Qual = []
    day = 999
    month = 999
    year = 999

    if use_altitude:
        field_h_alt = 9
    else:
        field_h_alt = 11

    ## open the file as binary because can be mixed with bin data
    for l in open(file_path, "rb"):
        # print(l)
        try:
            l = (l).decode("ascii")
            # print(l)

        except:
            continue

        if not l[0] == "$":
            continue
        if "GPZDA" in l or "GNZDA" in l:
            f = l.split(",")
            day = int(f[2])
            month = int(f[3])
            year = int(f[4])
        if "GPGGA" in l or "GNGGA" in l:
            f = l.split(",")
            h = int(f[1][0:2])
            m = int(f[1][2:4])
            s = int(f[1][4:6])
            qual = int(f[6])
            if f[5] == "W":
                EW = -1
            elif f[5] == "E":
                EW = 1
            else:
                print("WARN: problems....")
            if day == 999:
                continue
            if h == 0 and m == 0 and s == 0:
                continue  # c'est sale de zapper minuit mais bon ...
            t = dt.datetime(year, month, day, h, m, s)
            if not (startdate <= t <= enddate):
                continue
            T.append(t)
            Lat.append(float(f[2][0:2]) + float(f[2][2:]) / 60.0)
            Long.append(EW * (float(f[4][0:3]) + float(f[4][3:]) / 60.0))
            Haut.append(float(f[field_h_alt]))
            Qual.append(qual)

    X, Y, Z = conv.geo2xyz(Lat, Long, Haut)
    f, l, h = np.mean(Lat), np.mean(Long), np.mean(Haut)
    x0, y0, z0 = conv.geo2xyz(f, l, h)

    f, l, h = 43.4417981389, 7.83481522597, 6.59449264956
    x0, y0, z0 = 4595047.79934, 632288.017869, 4363273.52335
    E, N, U = conv.xyz2enu(X, Y, Z, x0, y0, z0)

    if export_path != "":
        outf = open(export_path, "w+")
        outf.write(" ".join(("#lat0,long0,h0 :", str(f), str(l), str(h), "\n")))
        outf.write(" ".join(("#x0 ,y0 , z0   :", str(x0), str(y0), str(z0), "\n")))
        for i in range(len(T)):
            datalis = [
                conv.dt2posix(T[i]),
                T[i].year,
                T[i].month,
                T[i].day,
                T[i].hour,
                T[i].minute,
                T[i].second,
                Lat[i],
                Long[i],
                Haut[i],
                X[i],
                Y[i],
                Z[i],
                E[i],
                N[i],
                U[i],
            ]
            datalis = [str(e) for e in datalis] + ["\n"]
            outf.write(",".join(datalis))  # )
        outf.close()

    if outtype == "ENU":
        OUTTUP = T, E, N, U, Qual
        colnames = ("T", "E", "N", "U", "Qual")
    elif outtype == "XYZ":
        OUTTUP = T, X, Y, Z, Qual
        colnames = ("T", "X", "Y", "Z", "Qual")
    else:
        OUTTUP = T, Lat, Long, Haut, Qual
        colnames = ("T", "F", "L", "H", "Qual")

    if df_out:
        DFout = pd.DataFrame(np.column_stack(OUTTUP))
        DFout.columns = colnames
        DFout = DFout.infer_objects()
        return DFout
    else:
        return OUTTUP


# strtd = dt.datetime(2015,06,19,16,00)
# T,E,N,U,Q = read_nmea(file_path='/home/psakicki/geodesea_nav_final.dat',enuout=1,export_path='/home/psakicki/Documents/geodesea_nav_matrix.dat') #,startdate=strtd)


def plot_nmea(T, E, N, U, Q):
    """a quick and dirty fct to plot the NMEA"""

    import matplotlib
    import matplotlib.pyplot as plt

    plt.close("all")

    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=1)

    plt.figure()
    plt.plot(T, U, "+")
    plt.title("U")
    ax = plt.gca()
    ax.yaxis.set_major_formatter(y_formatter)

    plt.figure()
    plt.plot(T, N, "+")
    plt.title("N")
    ax = plt.gca()
    ax.yaxis.set_major_formatter(y_formatter)

    plt.figure()
    plt.title("E")
    col = ["r"] * len(E)
    plt.scatter(T, E, Q, marker="+")
    ax = plt.gca()
    ax.yaxis.set_major_formatter(y_formatter)

    plt.figure()
    plt.title("E-N")
    plt.plot(E, N, "+")
    ax = plt.gca()
    ax.yaxis.set_major_formatter(y_formatter)

    plt.figure()
    plt.plot(T, Q, "+")
    ax = plt.gca()
    ax.yaxis.set_major_formatter(y_formatter)


def write_latlontime_file_4_OTPS_tide(
    outfilepath,
    lat,
    lon,
    strt,
    end=None,
    sec_step=1,
    generate_inp_file=True,
    tide_model_path="./DATA/Model_atlas",
    set_relative_path_in_inp_file=True,
):
    """


    Parameters
    ----------
    outfilepath : TYPE
        DESCRIPTION.
    lat : TYPE
        DESCRIPTION.
    lon : TYPE
        DESCRIPTION.
    strt : datetime
        start as a datetime.
    end : TYPE, optional
        length of the period in sec (integer)
        OR
        if None, strt is a list. The default is None.
    sec_step : TYPE, optional
        DESCRIPTION. The default is 1.
    generate_inp_file : TYPE, optional
        DESCRIPTION. The default is True.
    tide_model_path : TYPE, optional
        DESCRIPTION. The default is './DATA/Model_atlas'.
    set_relative_path_in_inp_file : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    outdir : TYPE
        DESCRIPTION.


    Note
    ----
    lat>0 - degrees North, lon>0 - degrees East
    lat<0 - degrees South, lon<0 - degrees West

    lat/lon can be a list, with the same size as dates_lis

    """

    fil = open(outfilepath, "w+")

    if not end:
        dates_lis = strt

    elif type(end) is dt.datetime:
        dates_lis = [strt]
        while dates_lis[-1] <= end:
            dates_lis.append(dates_lis[-1] + dt.timedelta(seconds=sec_step))

    else:
        dates_lis = [
            strt + dt.timedelta(seconds=x) for x in np.arange(0, end + 1, sec_step)
        ]

    if not utils.is_iterable(lat):
        lat = [float(lat)] * len(dates_lis)
        lon = [float(lon)] * len(dates_lis)
    else:
        if len(lat) != len(dates_lis) or len(lon) != len(dates_lis):
            print(
                "ERR : lat/lon are list, but not with the same size as the dates list, please check"
            )
            raise Exception

    print("INFO : first & last dates")
    print(dates_lis[0], dates_lis[-1])

    for d, lat_i, lon_i in zip(dates_lis, lat, lon):
        y = d.year
        m = d.month
        da = d.day
        h = d.hour
        mi = d.minute
        s = d.second
        lin = "{:12.8}{:12.8}{:7}{:7}{:7}{:7}{:7}{:7}\n".format(
            lat_i, lon_i, y, m, da, h, mi, s
        )
        fil.write(lin)

    if generate_inp_file:
        outdir = os.path.dirname(outfilepath)
        outname = os.path.basename(outfilepath) + ".inp"
        outoutname = os.path.basename(outfilepath) + ".out"

        fil2 = open(os.path.join(outdir, outname), "w+")

        fil2.write(tide_model_path + "\n")
        if not set_relative_path_in_inp_file:
            fil2.write(outfilepath + " \n")
        else:
            fil2.write("./" + os.path.basename(outfilepath) + " \n")

        fil2.write("z" + " \n")
        fil2.write("" + " \n")
        fil2.write("AP" + " \n")
        fil2.write("geo" + " \n")
        fil2.write("1" + " \n")
        if not set_relative_path_in_inp_file:
            fil2.write(os.path.join(outdir, outoutname) + " \n")
        else:
            fil2.write("./" + outoutname + " \n")

        fil2.close()

    return outdir


def read_OTPS_tide_file(pathin, print_line=False, return_array=True):
    """
    Return :
        latlis , lonlis , datlis , hlis
    """
    filout = open(pathin)

    latlis, lonlis, datlis, hlis = [], [], [], []

    for i, l in enumerate(filout):
        f = l.split()
        if np.mod(i, 1000) == 0 and print_line:
            print("line", i)
        try:
            latlis.append(float(f[0]))
            lonlis.append(float(f[1]))
            datlis.append(dateutil.parser.parse(f[2] + " " + f[3]))
            hlis.append(float(f[4]))
        except:
            continue
    if return_array:
        return np.array(latlis), np.array(lonlis), np.array(datlis), np.array(hlis)
    else:
        return latlis, lonlis, datlis, hlis


def unzip_gz_z(
    inp_gzip_file, out_gzip_file="", remove_inp=False, force=False, out_gzip_dir=None
):
    """
    frontend function to unzip files (gzip and legacy Z compression)

    Parameters
    ----------
    inp_gzip_file : str
        the path of the input file.
    out_gzip_file : str, optional
        if out_gzip_file AND out_gzip_dir not precised
        file will be extracted in the same folder
        as inp_gzip_file. The default is ''.
    remove_inp : bool, optional
        remove the input file. The default is False.
    force : bool, optional
        force the decompression. The default is False.
    out_gzip_dir : str, optional
        output directory. will be used only if out_gzip_file == ''
        unzipped file will keep the same


    Returns
    -------
    out_gzip_file :
        path of the uncompressed file.

    Warning
    -------

    .Z decompression is implemented, but is very unstable (avoid .Z, prefer .gz)
    """

    import gzip

    if inp_gzip_file.endswith(".gz"):
        is_gz = True
    elif inp_gzip_file.endswith(".Z"):
        is_gz = False
    else:
        is_gz = True

    if out_gzip_file == "":
        if not out_gzip_dir:
            out_gzip_dir_use = os.path.dirname(inp_gzip_file)
        else:
            out_gzip_dir_use = out_gzip_dir

        if not os.path.isdir(out_gzip_dir_use):
            log.error("%s does not exists, unzipping is cancelled", out_gzip_dir_use)
            raise Exception

        out_gzip_file = os.path.join(
            out_gzip_dir_use, ".".join(os.path.basename(inp_gzip_file).split(".")[:-1])
        )

    if os.path.isfile(out_gzip_file) and not force:
        log.info("%s already exists, skiping (use force option)", out_gzip_file)
        pass
    else:
        if is_gz:
            with gzip.open(inp_gzip_file, "rb") as f_in, open(
                out_gzip_file, "wb"
            ) as f_out:
                shutil.copyfileobj(f_in, f_out)
        else:
            log.warning(
                "zlib decompress is unstable!! and .Z should be definitly avoided..."
            )
            out_gzip_file = utils.uncompress(inp_gzip_file)

        log.info("uncompress %s to  %s", inp_gzip_file, out_gzip_file)

    ### Removing part
    if remove_inp and type(remove_inp) is bool and os.path.getsize(out_gzip_file) > 0:
        log.info("removing %s")
        os.remove(inp_gzip_file)

    return out_gzip_file


#  ______                _   _                _____                                         _
# |  ____|              | | (_)              / ____|                                       | |
# | |__ _   _ _ __   ___| |_ _  ___  _ __   | |  __ _ __ __ ___   _____ _   _  __ _ _ __ __| |
# |  __| | | | '_ \ / __| __| |/ _ \| '_ \  | | |_ | '__/ _` \ \ / / _ \ | | |/ _` | '__/ _` |
# | |  | |_| | | | | (__| |_| | (_) | | | | | |__| | | | (_| |\ v /  __/ |_| | (_| | | | (_| |
# |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|  \_____|_|  \__,_| \_/ \___|\__, |\__,_|_|  \__,_|
#                                                                        __/ |
#                                                                       |___/


########################################################################'

# ## Read the blocs
# DFcoor = files_rw.read_sinex_versatile(sinex_path_in,"SOLUTION/ESTIMATE",
#                                        header_line_idx=None)
# DFepoc = files_rw.read_sinex_versatile(sinex_path_in,"SOLUTION/EPOCHS",
#                                        header_line_idx=None)

# ## Rename the Coord DF
# DFcoor.rename(columns={0: "INDEX",
#                        1: "TYPE",
#                        2: "CODE",
#                        3: "PT",
#                        4: "SOLN",
#                        5: "REF_EPOCH",
#                        6: "UNIT",
#                        7: "N",
#                        8: "ESTIMATED_VALUE",
#                        9: "STD_DEV"},inplace=True)

# ## Rename the Epoch DF
# DFepoc.rename(columns={0: "CODE",
#                        1: "PT",
#                        2: "SOLN",
#                        3: "T",
#                        4: "DATA_START",
#                        5: "DATA_END",
#                        6: "MEAN_EPOCH"},inplace=True)


# try:
#     DFepoc[DFepoc["SOLN"].str.contains("----")] == 1
# except:
#     pass

# DFcoor = DFcoor[DFcoor["TYPE"].str.contains("STA")]

# AAA = DFcoor.groupby(["CODE",
#                       "SOLN"])
# for a in AAA:
#     print(a[1]["ESTIMATED_VALUE"].T)


# Convert sp3 2 gins
# brouillon
#
# satdic = OrderedDict()
#
# sp3in = "/media/psakicki/AF0E-DA43/CHAL_BUGS/TEMP_DATA/jp211236.sp3"
# orbginsout = '/home/psakicki/aaa_FOURBI/outgins'
#
# for l in open(sp3in):
#    if l[0] in ('#','+','%','/'):
#        continue
#
#    f = l.split()
#
#    if l[0] == '*':
#        ll = [int(float(e)) for e in f[1:]]
#        epoch = dt.datetime(*ll)
#        continue
#    if l[0] == 'p':
#        prn = int(f[0][2:])
#        x = float(f[1]) * 10**3
#        y = float(f[2]) * 10**3
#        z = float(f[3]) * 10**3
#
#        print epoch,prn,x,y,z
#
#        if not satdic.has_key(prn):
#            satdic[prn] = []
#
#        jjul = conv.dt2jjul_cnes(epoch)
#        sec = conv.dt2secinday(epoch) + 19
#
#        satdic[prn].append((jjul,sec,x,y,z))
#
# fil = open(orbginsout,'w+')
#
# for k,vv in satdic.iteritems():
#    for v in vv:
#        if k in (11,13,14,18,20,28):
#            finalprn = 66600 + k
#        elif k in (2,15,17,21):
#            finalprn = 88800 + k
#        else:
#            finalprn = 77700 + k
#
#        jjul = v[0]
#        sec  = v[1]
#        x = v[2]
#        y = v[3]
#        z = v[4]
#
#        finalline =  '  {:5} {:6} {:12.6f} tai xyz ine  0 {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:.3f}  {:.3E}\n'.format(finalprn,jjul,float(sec),0,0,0,0,0,0,0,0,0,x,y,z,0,0,0,0.,0)
#        fil.write(finalline)
#
# fil.close()


#                   ___    _____   ____  __  __ ______  _____                                  _
#                  |__ \  |  __ \ / __ \|  \/  |  ____|/ ____|                                | |
#   ___ _____   __    ) | | |  | | |  | | \  / | |__  | (___    _ __ ___  __ _ _   _  ___  ___| |_
#  / __/ __\ \ / /   / /  | |  | | |  | | |\/| |  __|  \___ \  | '__/ _ \/ _` | | | |/ _ \/ __| __|
# | (__\__ \\ v /   / /_  | |__| | |__| | |  | | |____ ____) | | | |  __/ (_| | |_| |  __/\__ \ |_
#  \___|___/ \_/   |____| |_____/ \____/|_|  |_|______|_____/  |_|  \___|\__, |\__,_|\___||___/\__|
#                                                                           | |
#                                                                           |_|

# def csv2DOMES(csvin):

# filein = '/home/psakicki/gin/TP/GWADA/gps_continu_full.csv'
# import pandas as pd
# df = pd.read_csv(open(filein))
# saved_column = list(df['Nom'])
#
#
# for i,row in df.iterrows():
#    name = eval(row['Nom'])
#    stat = row['Alias']
#    outdomesfile = open('/home/psakicki/gin/TP/GWADA/' + stat + '_DOMES_request.txt' ,'w+' )
#    outdomesfile.write(' 1. Request from (full name) : ' + '\n')
#    outdomesfile.write('Agency                   : '+ 'Institut de Pysique du Globe de Paris' +  '\n')
#    outdomesfile.write('E-mail                   : '+ '\n')
#    outdomesfile.write('Date                     : '+ '\n')
#    outdomesfile.write(' '+ '\n')
#    outdomesfile.write('2. Site Name                : ' + name + '\n')
#    outdomesfile.write('3. Country                  : ' + 'Guadeloupe , French West Indies' + '\n'  )
#    outdomesfile.write('4. Point Description        : ' + row['Type'] + '\n' )
#    outdomesfile.write('5. DOMES Number             : ' + 'N/A' + '\n' )
#    outdomesfile.write('6. Local Number             : ' + 'N/A' + '\n')
#    outdomesfile.write('7. 4-Char Code              : ' + row['Alias'] + '\n' )
#    outdomesfile.write(' '+ '\n')
#    outdomesfile.write('8. Approximate Position ' + '\n')
#    outdomesfile.write('   Latitude (deg min)       : '  + str(row['Lat. (WGS84)']) + '\n' )
#    outdomesfile.write('   Longitude (deg min)      : '  + str(row['Lon. (WGS84)']) + '\n' )
#    outdomesfile.write('   Elevation (m)            : '  + str(row['Elev. (m)']) + '\n')
#    outdomesfile.write(' '+ '\n')
#    outdomesfile.write('9. Instrument               : '  + str('GPS on ' + row['Type']) + '\n' )
#    outdomesfile.write('10. Date of Installation     : ' + str(row['Début / Installation']) + '\n')
#    outdomesfile.write(' '+ '\n')
#    outdomesfile.write('11. Operation Contact Name   : ' + '\n')
#    outdomesfile.write('    Agency                   : ' + '\n')
#    outdomesfile.write('    E-mail                   : ' + '\n')
#    outdomesfile.write(' ' + '\n')
#    outdomesfile.write('12. Site Contact Name        : ' + '\n')
#    outdomesfile.write('    Agency                   : ' + '\n')
#    outdomesfile.write('    E-mail                   : ' + '\n')


#  ________   ________ __  __ _____  _      ______  _____
# |  ____\ \ / /  ____|  \/  |  __ \| |    |  ____|/ ____|
# | |__   \ v /| |__  | \  / | |__) | |    | |__  | (___
# |  __|   > < |  __| | |\/| |  ___/| |    |  __|  \___ \
# | |____ / . \| |____| |  | | |    | |____| |____ ____) |
# |______/_/ \_\______|_|  |_|_|    |______|______|_____/
#
#
## CATS
#
# neulist = ['/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wABMF_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wADE0_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wASF0_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wCBE0_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wDHS0_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wDSD0_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wFNA0_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wFSDC_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wHOUE_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wLAM0_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wMGL0_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wPDB0_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wPSA1_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wSOUF_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wTDB0_GINS_PS.neu']
#
# statinfo_path = '/home/psakicki/THESE/CODES/stationinfo2gins/station.info.ovsg'
#
# for neu in neulist:
#    print neu
#    statinfo_2_cats(statinfo_path,neu)
#
# read_station_info_solo(statinfo_path,'ABMF')

# read_station_info_solo('/home/psakicki/THESE/CODES/stationinfo2gins/station.info.ovsg','HOUE')
# statname_of_catsfile('/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wASF0_GINS_PS.neu')


# GINS
# station_info_2_gins('/home/psakicki/gin/TP/GWADA/station.info.ovsg.volobsis','/home/psakicki/gin/TP/GWADA/lfile_mod_mk2','/home/psakicki/gin/TP/GWADA/GWADA_mk2.stat')
# coordfile = '/home/psakicki/gin/TP/GWADA/RINX2/KARIB/pbo.snaps_igs08.vel'
# statinfo = '/home/psakicki/gin/TP/GWADA/RINX2/KARIB/station.info'

# specific = ['CN40' , 'CN19' , 'CN38' , 'SAN0' , 'CN35' , 'ROA0' , 'CN29' , 'CN18' , 'CBMD' , 'LCSB' , 'GCFS' , 'GCEA' , 'CN10' , 'CN12' , 'PRCG' , 'SMRT' , 'BARA']

# station_info_2_gins(statinfo,coordfile,'/home/psakicki/gin/TP/GWADA/KARIB_mk3.stat',specific_stats_lis=specific)
