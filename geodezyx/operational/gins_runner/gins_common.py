# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.operational contains functions to run the 
GNSS processing software GINS. 

it can be imported directly with:
from geodezyx import operational

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""


########## BEGIN IMPORT ##########
#### External modules
import collections
import copy
import datetime as dt
import glob
#### Import the logger
import logging
import multiprocessing as mp
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import yaml

#### geodeZYX modules
from geodezyx import conv
from geodezyx import files_rw
from geodezyx import operational
from geodezyx import time_series
from geodezyx import utils

log = logging.getLogger('geodezyx')

##########  END IMPORT  ##########


def get_gins_path(extended=False):
    try:
        gs_user = os.environ["GS_USER"]
    except:
        print("ERR : env. var. $GS_USER dont exists !!!")
        return None
        # gs_user =  os.environ['HOME']
    if not extended:
        return gs_user
    else:
        return os.path.join(os.environ["GS_USER"], "gin")


def check_if_stat_in_stationfile(stat, stationfile):
    boolout = utils.check_regex(stationfile, stat)
    if not boolout:
        print("WARN :", stat, "not in", stationfile)
        print("       check your RINEX header and particularly")
        print("       the MARKER NAME field (station 4-char. code)")
    return boolout


def check_if_DOMES_in_oceanloadfile(domes, oclofile):
    boolout = utils.check_regex(oclofile, "^   " + str(domes))
    if not boolout:
        print("WARN :", domes, "not in", oclofile)
    return boolout


def find_DOMESstat_in_stationfile(stat, stationfile):
    fil = open(stationfile)
    for l in fil:
        if stat in l:
            f = l.split()
            return f[0], f[1]
    return 00000, "M000"


def check_if_file_in_gin_folder(a_file_path, gins_path=None):
    if not gins_path:
        gins_path = os.path.join(get_gins_path(), "gin")

    real_path_file = os.path.realpath(a_file_path)
    real_path_gins = os.path.realpath(gins_path)

    if real_path_gins in real_path_file:
        boolout = True
    else:
        boolout = False

    if not boolout:
        log.warning("%s not in %s", a_file_path, gins_path)

    return boolout

def copy_file_in_gin_folder(
    a_file_path, repository_folder=""
):
    if repository_folder == "":
        repository_folder = os.path.join(get_gins_path(), "TEMP_DATA")
    if not os.path.exists(repository_folder):
        os.makedirs(repository_folder)

    shutil.copy2(a_file_path, repository_folder)

    out = os.path.join(repository_folder, os.path.basename(a_file_path))

    return out


def bring_to_gin_folder(a_file_path, temp_data_folder, gins_path=None):
    if not gins_path:
        gins_path = os.path.join(get_gins_path(), "gin")

    if not check_if_file_in_gin_folder(a_file_path, gins_path):
        log.info("copy %s > %s",a_file_path, temp_data_folder)
        file_path_out = copy_file_in_gin_folder(a_file_path, temp_data_folder)
    else:
        file_path_out = a_file_path

    return file_path_out


def check_good_exec_of_gins(streamin, director_name):
    if "Exécution terminée du fichier" in streamin.read():
        print("INFO : happy end for " + director_name + " :)")
        return True
    else:
        print("WARN : bad end for " + director_name + " :(")
        return False


def make_path_ginsstyle(pathin):
    """input path must be an absolute path with /gin/ inside
    output will be .temp.gin/<rest of the path>"""

    if pathin.startswith(".temp.gin"):
        log.info("already a director compatible path: %s", pathin)
        return pathin

    if not "/gin/" in pathin:
        log.error("not /gin/ in %s", pathin)
        return None

    rnx_path_lis = pathin.split("/")
    i_gin = rnx_path_lis.index("gin")
    rnx_path_lis2 = rnx_path_lis[i_gin:]
    rnx_path_lis2[0] = ".temp.gin"
    pathout = os.path.join(*rnx_path_lis2)
    return pathout


def write_oceanload_file(station_file, oceanload_out_file, fes_yyyy=2004):
    temp_cmd_fil = os.path.join(os.path.dirname(oceanload_out_file), "oclo.cmd.tmp")
    temp_cmd_filobj = open(temp_cmd_fil, "w")

    temp_cmd_filobj.write(station_file + "\n")
    temp_cmd_filobj.write(oceanload_out_file + "\n")
    if fes_yyyy == 2012:
        exe_loadoce_cmd = "exe_loadoce_fes2012"
    elif fes_yyyy == 2004:
        exe_loadoce_cmd = "exe_loadoce"

    p = subprocess.Popen(
        exe_loadoce_cmd + " < " + temp_cmd_fil,
        shell=True,
        executable="/bin/bash",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    #    os.remove(temp_cmd_fil)
    temp_cmd_filobj.close()

    if True:
        return oceanload_out_file
    else:
        return ""


def get_temp_data_gins_path():

    gins_path = get_gins_path()

    # be sure there is a TEMP DATA folder
    temp_data_folder = os.path.join(gins_path, "gin", "TEMP_DATA")

    if not os.path.exists(temp_data_folder):
        os.makedirs(temp_data_folder)

    return temp_data_folder


def dirs_copy_generik_2_working(
    director_name_prefix, director_generik_path, temp_data_folder=None
):

    if not temp_data_folder:
        temp_data_folder = get_temp_data_gins_path()

    shutil.copy(director_generik_path, temp_data_folder)
    director_generik_path_tmp = os.path.join(
        temp_data_folder, os.path.basename(director_generik_path)
    )
    director_generik_path_out = os.path.join(
        temp_data_folder, director_name_prefix + ".yml"
    )
    os.rename(director_generik_path_tmp, director_generik_path_out)

    return director_generik_path_out


def get_director_list(wildcard_dir):
    """with a wildcard (e.g. 'GWADA_MK2*') and return a list of corresponding
    directors found in gin/data/directeur folder"""
    gins_path = get_gins_path()
    di_run_lis = [
        os.path.basename(e)
        for e in glob.glob(os.path.join(gins_path, "data", "directeur", wildcard_dir))
    ]
    return di_run_lis


def get_rinex_list(
    parent_folder,
    specific_stats=[],
    invert=False,
    compressed=True,
    start=dt.datetime(1980, 1, 1),
    end=dt.datetime(2099, 1, 1),
):
    """return :
        all RINEXs found in a parent folder and his subfolders
        (compressed or not)

    parent_folder :
        can be a the path of the partent folder (the rinex archive path)
        but also a RINEX list already found
        (modification of 170524 to gain speed)

    specific_stats :
        MUST BE a list ['STA1','STA2','STA3']
        So, if only one elt in specific_stats, use a syntax as : ['STA1']

    invert :
        False = keeping the specific stats
        True  = removing the specific stats
        or for all stations, leave a empty
        tuple in specific_stats

    NB : the end date is included
    """
    if type(parent_folder) is str:
        wholefilelist = []
        for root, dirs, files in os.walk(parent_folder, topdown=False):
            for name in files:
                wholefilelist.append(os.path.join(root, name))

        wholefilelist = list(set(wholefilelist))

        if compressed:
            rnxregex = operational.rinex_regex()
        else:
            rnxregex = operational.rinex_regex(False)

        rinexfilelist = [fil for fil in wholefilelist if re.search(rnxregex, fil)]
    elif utils.is_iterable(parent_folder):
        # is parent_folder is a list containing already found RINEXs
        rinexfilelist = parent_folder

    specific_stats = [e.lower() for e in specific_stats]

    if specific_stats != []:
        if not invert:
            goodrnxfilelist = [
                fil
                for fil in rinexfilelist
                if os.path.basename(fil)[0:4] in specific_stats
            ]
        else:
            goodrnxfilelist = [
                fil
                for fil in rinexfilelist
                if os.path.basename(fil)[0:4] not in specific_stats
            ]
    else:
        goodrnxfilelist = rinexfilelist

    goodrnxfilelist = sorted(goodrnxfilelist)

    if goodrnxfilelist == []:
        print(
            "WARN : get_rinex_list : no RINEX found, check the path of the parent folder !"
        )
    else:
        print("INFO : get_rinex_list :", len(goodrnxfilelist), "RINEXs found")

    goodrnxfilelist_date = []
    end = end + dt.timedelta(seconds=86399)
    for rnx in goodrnxfilelist:
        dtrnx = conv.rinexname2dt(os.path.basename(rnx))
        if start <= dtrnx <= end:
            goodrnxfilelist_date.append(rnx)
    if len(goodrnxfilelist_date) != len(goodrnxfilelist):
        dif = len(goodrnxfilelist) - len(goodrnxfilelist_date)
        print(
            "INFO : get_rinex_list :",
            dif,
            "RINEXs removed bc. not in the date interval",
        )

    if goodrnxfilelist_date == []:
        print("WARN : get_rinex_list : all RINEXs removed bc. of date interval")

    return goodrnxfilelist_date

def sort_by_stations(archive_path, wildcard, i):
    """i is the indice of the first character of the station name in
    eg : i = 18 for filename KARIB_MK3_FLH_v2__bara_22282_2011_003

    archive path is ABSOLUTE b.c. it can be outside of the gins folder"""

    path = os.path.join(archive_path, wildcard)
    filis = glob.glob(path)

    firstfil = os.path.basename(filis[0])
    print("sure the stat name is good ? ")
    print("e.g.", firstfil[i : i + 4], "for", firstfil)
    print("(you have 5sec to abort if not)")
    time.sleep(5)

    for f in filis:
        stat = os.path.basename(f)[i : i + 4]
        STAT = stat.upper()
        archiv = os.path.join(archive_path, STAT)
        if not os.path.exists(archiv):
            os.makedirs(archiv)
        shutil.move(f, archiv)

    return None

def merge_yaml(yaml1, yaml2, yaml_out=None):

    dic1 = yaml.load(open(yaml1))
    dic2 = yaml.load(open(yaml2))

    def merge(dic1, dic2):
        if isinstance(dic1, dict) and isinstance(dic2, dict):
            for k, v in dic2.items():
                if k not in dic1:
                    dic1[k] = v
                else:
                    dic1[k] = merge(dic1[k], v)
        return dic1

    dic3 = merge(dic1, dic2)

    if yaml_out:
        with open(yaml_out, "w+") as outfile:
            outfile.write(yaml.dump(dic3, default_flow_style=False))

    return dic3


#    return dic1

# yaml1    = '/home/psakicki/GINS/gin/data/EXE_PPP/DIR_REF_PPP.yml'
# yaml2    = '/home/psakicki/GINS/gin/TP/TP_RELAX/DIR_MC0.yml'
# yaml_out = '/home/psakicki/GINS/gin/TP/TP_RELAX/DIR_MC0_perso2.yml'
#
# merge_yaml(yaml1,yaml2,yaml_out)
