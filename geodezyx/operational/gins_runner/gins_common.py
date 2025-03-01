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
import datetime as dt
import glob


import os
import re
import shutil
import subprocess
import time
import yaml

#### geodeZYX modules
from geodezyx import conv
from geodezyx import operational
from geodezyx import utils

#### Import the logger
import logging
log = logging.getLogger("geodezyx")

##########  END IMPORT  ##########


def get_gin_path(extended=False):
    """
    Get the GIN path from the environment variable GS_USER.

    Parameters
    ----------
    extended : bool, optional
        If True, returns the path to the 'gin' directory inside the GIN path.
        Defaults to False.

    Returns
    -------
    str or None
        The GIN path or None if the environment variable GS_USER does not exist.
    """
    try:
        gs_user = os.environ["GS_USER"]
    except:
        log.error("env. var. $GS_USER dont exists !!!")
        return None
        # gs_user =  os.environ['HOME']
    if not extended:
        return gs_user
    else:
        return os.path.join(os.environ["GS_USER"], "gin")


def get_spotgins_path():
    try:
        return os.environ["SPOTGINS_DIR"]
    except:
        log.error("env. var. SPOTGINS_DIR dont exists !!!")
        return None

def check_stat_in_statfile(stat, stationfile):
    """
    Check if a station code is present in a station file.

    Parameters
    ----------
    stat : str
        The station code to search for.
    stationfile : str
        The path to the station file.

    Returns
    -------
    bool
        True if the station code is found in the station file, False otherwise.
    """
    boolout = utils.check_regex(stationfile, stat)
    if not boolout:
        log.warning("%s not in %s", stat, stationfile)
        log.warning("check your RINEX header and its the MARKER NAME field (station 4-char. code)")
    return boolout

def check_domes_in_oclo(domes, oclofile):
    """
    Check if a DOMES number is present in an OCLO file.

    Parameters
    ----------
    domes : str
        The DOMES number to search for.
    oclofile : str
        The path to the OCLO file.

    Returns
    -------
    bool
        True if the DOMES number is found in the OCLO file, False otherwise.
    """
    boolout = utils.check_regex(oclofile, "^   " + str(domes))
    if not boolout:
        log.warning("%s not in %s", domes, oclofile)

    return boolout



def find_domes_in_statfile(stat_code, stationfile):
    fil = open(stationfile)
    for l in fil:
        if stat_code in l:
            f = l.split()
            return f[0], f[1]
    return "00000", "M000"


def check_if_in_gin(file_path_inp, gins_path=None):
    """
    Check if a file is within the GIN path.

    Parameters
    ----------
    file_path_inp : str
        The path of the file to check.
    gins_path : str, optional
        The GIN path to check against. If not provided, it defaults to the 'gin' directory inside the GIN path.

    Returns
    -------
    bool
        True if the file is within the GIN path, False otherwise.
    """
    if not gins_path:
        gins_path = os.path.join(get_gin_path(), "gin")

    real_path_file = os.path.realpath(file_path_inp)
    real_path_gins = os.path.realpath(gins_path)

    if real_path_gins in real_path_file:
        boolout = True
    else:
        boolout = False

    if not boolout:
        log.warning("%s not in %s", file_path_inp, gins_path)

    return boolout


def copy_in_gin(file_path_inp, temp_data_folder=None):
    """
    Copy a file to the GIN temporary data folder.

    Parameters
    ----------
    file_path_inp : str
        The path of the file to be copied.
    temp_data_folder : str, optional
        The path to the temporary data folder. If not provided, it defaults to 'TEMP_DATA' inside the GIN path.

    Returns
    -------
    str
        The path of the copied file in the temporary data folder.
    """
    if not temp_data_folder:
        temp_data_folder = os.path.join(get_gin_path(), "TEMP_DATA")
    if not os.path.exists(temp_data_folder):
        os.makedirs(temp_data_folder)

    log.info("copy %s > %s", file_path_inp, temp_data_folder)
    shutil.copy2(file_path_inp, temp_data_folder)

    file_path_out = os.path.join(temp_data_folder, os.path.basename(file_path_inp))

    return file_path_out

def bring_to_gin(file_path_inp, temp_data_folder=None, gins_path=None):
    """
    Ensure a file is within the GIN path, copying it to the temporary data folder if necessary.

    Combo of check_if_in_gin and copy_in_gin functions.

    Parameters
    ----------
    file_path_inp : str
        The path of the file to be checked and possibly copied.
    temp_data_folder : str, optional
        The path to the temporary data folder. If not provided, it defaults to 'TEMP_DATA' inside the GIN path.
    gins_path : str, optional
        The GIN path to check against. If not provided, it defaults to the 'gin' directory inside the GIN path.

    Returns
    -------
    str
        The path of the file within the GIN path.
    """
    if not check_if_in_gin(file_path_inp, gins_path=gins_path):
        file_path_out = copy_in_gin(file_path_inp, temp_data_folder=temp_data_folder)
    else:
        file_path_out = file_path_inp

    return file_path_out


def check_gins_exe(streamin, director_name):
    """
    Check the execution status of a GINS director.

    Parameters
    ----------
    streamin : file-like object
        The input stream to read the execution output from.
    director_name : str
        The name of the director being checked.

    Returns
    -------
    bool
        True if the execution was successful, False otherwise.
    """
    if "Exécution terminée du fichier" in streamin.read():
        print("INFO : happy end for " + director_name + " :)")
        return True
    else:
        print("WARN : bad end for " + director_name + " :(")
        return False


def make_path_ginsstyle(pathin):
    """
    Convert an absolute path containing '/gin/' to a GINS director-compatible path.

    Parameters
    ----------
    pathin : str
        The input path, which must be an absolute path
        containing '/gin/'.

    Returns
    -------
    str or None
        The converted path with '.temp.gin' as the root directory,
         or None if the input path does not contain '/gin/'.
    """
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


def write_oclo_file(station_file, oceanload_out_file, fes_yyyy=2004):
    temp_cmd_fil = os.path.join(os.path.dirname(oceanload_out_file), "oclo.cmd.tmp")
    temp_cmd_filobj = open(temp_cmd_fil, "w")

    temp_cmd_filobj.write(station_file + "\n")
    temp_cmd_filobj.write(oceanload_out_file + "\n")
    if fes_yyyy == 2012:
        exe_loadoce_cmd = "exe_loadoce_fes2012"
    elif fes_yyyy == 2004:
        exe_loadoce_cmd = "exe_loadoce"
    else:
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

    gins_path = get_gin_path()

    # be sure there is a TEMP DATA folder
    temp_data_folder = os.path.join(gins_path, "gin", "TMP_GYNS")

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
    gins_path = get_gin_path()
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


def check_solution(dir_name_inp, gin_path_inp=None):

    dir_name_inp = os.path.basename(dir_name_inp)

    if not gin_path_inp:
        gin_path = get_gin_path()
    else:
        gin_path = gin_path_inp
    # check if the solution already exists
    sol_folder = os.path.join(gin_path, "gin", "batch", "solution")
    sols_matching = glob.glob(sol_folder + '/' + dir_name_inp + '*')
    if len(sols_matching) > 0:
        log.info("Solution %s exists for %s", sols_matching, dir_name_inp)
    return sols_matching

#    return dic1

# yaml1    = '/home/psakicki/GINS/gin/data/EXE_PPP/DIR_REF_PPP.yml'
# yaml2    = '/home/psakicki/GINS/gin/TP/TP_RELAX/DIR_MC0.yml'
# yaml_out = '/home/psakicki/GINS/gin/TP/TP_RELAX/DIR_MC0_perso2.yml'
#
# merge_yaml(yaml1,yaml2,yaml_out)
