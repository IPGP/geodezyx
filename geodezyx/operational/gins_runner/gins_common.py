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


def check_site_stfl(stat, stationfile):
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
        log.warning("check RINEX header and its MARKER NAME field (4-char. code)")
    return boolout


def chek_domes_oclo(domes, oclofile):
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
        log.error("%s not in %s, %s", domes, oclofile, "be sure of what you are doing!")

    return boolout


def find_domes_in_stfl(stat_code, stationfile):
    """
    Find the DOMES number for a station code in a station file.
    """
    fil = open(stationfile)
    for l in fil:
        if stat_code in l:
            f = l.split()
            return f[0], f[1]
    return "00000", "M000"


def check_if_in_gin(file_path_inp, gins_path=None, verbose=True):
    """
    Check if a file is within the GIN path.

    Parameters
    ----------
    file_path_inp : str
        The path of the file to check.
    gins_path : str, optional
        The GIN path to check against. If not provided, it defaults to the 'gin' directory inside the GIN path.
    verbose : bool
        verbose messages. Default is True.

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

    if not boolout and verbose:
        log.warning("%s not in %s", file_path_inp, gins_path)

    return boolout


def copy_in_gin(file_path_inp, temp_data_folder=None, verbose=True):
    """
    Copy a file to the GIN temporary data folder.

    Parameters
    ----------
    file_path_inp : str
        The path of the file to be copied.
    temp_data_folder : str, optional
        The path to the temporary data folder.
        If not provided, it defaults to 'TEMP_DATA' inside the GIN path.
    verbose : bool
        verbose messages. Default is True.

    Returns
    -------
    str
        The path of the copied file in the temporary data folder.
    """
    if not temp_data_folder:
        temp_data_folder = os.path.join(get_gin_path(), "TEMP_DATA")
    if not os.path.exists(temp_data_folder):
        os.makedirs(temp_data_folder)
    if verbose:
        log.info("copy %s > %s", file_path_inp, temp_data_folder)
    shutil.copy2(file_path_inp, temp_data_folder)

    file_path_out = os.path.join(temp_data_folder, os.path.basename(file_path_inp))

    return file_path_out


def bring_to_gin(file_path_inp, temp_data_folder=None, gins_path=None, verbose=True):
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
    verbose : bool
        verbose messages. Default is True.

    Returns
    -------
    str
        The path of the file within the GIN path.
    """
    if not check_if_in_gin(file_path_inp, gins_path=gins_path, verbose=verbose):
        file_path_out = copy_in_gin(file_path_inp, temp_data_folder=temp_data_folder, verbose=verbose)
    else:
        file_path_out = file_path_inp

    return file_path_out

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

def rm_old_tmp_shel(directory='/tmp', age_sec=600):
    """
    Removes files in the specified directory older than the given age in minutes.

    Parameters:
    directory (str): The directory to search in.
    age_sec (int): Age of files in minutes to delete. Default is 10 minutes.
    """
    try:
        # Call the `find` command to delete files older than the specified age
        subprocess.run(
            ['find', directory, '-type', 'f', '-mmin', f'+{age_sec}', '-delete'],
            check=True
        )
    except subprocess.CalledProcessError as e:
        log.error(f"Error: {e.stderr}")


def get_temp_data_gins_path():
    """
    Get the path to the temporary data folder within the GIN directory.
    """

    gins_path = get_gin_path()

    # be sure there is a TEMP DATA folder
    temp_data_folder = os.path.join(gins_path, "gin", "TMP_GYNS")

    if not os.path.exists(temp_data_folder):
        os.makedirs(temp_data_folder)

    return temp_data_folder

def check_solution(dir_name_inp, gin_path_inp=None, verbose=True):
    """
    Check if a solution already exists for the given directory name.
    """

    dir_name_inp = os.path.basename(dir_name_inp)

    if not gin_path_inp:
        gin_path = get_gin_path()
    else:
        gin_path = gin_path_inp
    # check if the solution already exists
    sol_folder = os.path.join(gin_path, "gin", "batch", "solution")
    sols_matching = glob.glob(sol_folder + "/" + dir_name_inp + "*")
    if len(sols_matching) > 0 and verbose:
        log.info("Solution %s exists for %s", sols_matching, dir_name_inp)
    return sols_matching


