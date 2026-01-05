#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.operational contains functions to manipulate
RINEX files (version 2, quite obsolate).

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

#### Import the logger
import logging
import os
import re
import shutil
import string
import subprocess

import dateutil
import hatanaka
import numpy as np
import pandas as pd

#### geodeZYX modules
from geodezyx import conv
from geodezyx import utils

log = logging.getLogger("geodezyx")


##########  END IMPORT  ##########


def read_rnx_epoch_line(line, rnx2=True):
    """
    read the epoch line of a RINEX (2 or 3)
    and extract it as a datetime

    Parameters
    ----------
    line : string
        input line, assuming it's a RINEX epoch line.
    rnx2 : TYPE, optional
        True if RINEX2, False if RINEX3/4.
        The default is True.

    Returns
    -------
    epochdt : datetime
        The epoch of the line.
    """

    if rnx2:
        strtidx = 0
    else:
        strtidx = 1
    fraw = line[strtidx:].split()  # fraw est la pour gerer les fracion de sec ...
    fraw = fraw[0:6]
    fraw = [float(e) for e in fraw]

    f = [int(e) for e in fraw]
    msec = fraw[5] - np.floor(fraw[5])
    msec = np.round(msec, 4)
    msec = int(msec * 10**6)
    f.append(msec)  # ajour des fractions de sec

    if rnx2:
        if f[0] < 50:
            f[0] = f[0] + 2000
        else:
            f[0] = f[0] + 1900

    if f[5] == 60:  # cas particulier rencontré dans des rinex avec T = 60sec
        rinex_60sec = True
        f[5] = 59
    else:
        rinex_60sec = False

    epochdt = dt.datetime(*f)

    if rinex_60sec:
        epochdt = epochdt + dt.timedelta(seconds=1)
    return epochdt


#  _____  _____ _   _ ________   __   _____       _ _ _
# |  __ \|_   _| \ | |  ____\ \ / /  / ____|     | (_) |
# | |__) | | | |  \| | |__   \ v /  | (___  _ __ | |_| |_ ___ _ __
# |  _  /  | | | . ` |  __|   > <    \___ \| '_ \| | | __/ _ \ '__|
# | | \ \ _| |_| |\  | |____ / . \   ____) | |_) | | | ||  __/ |
# |_|  \_\_____|_| \_|______/_/ \_\ |_____/| .__/|_|_|\__\___|_|
#                                          | |
#                                          |_|


def check_if_compressed_rinex(rinex_path):
    """
    Check if a RINEX file is compressed (CRZ or O.Z / D.Z)

    Parameters
    ----------
    rinex_path : string
        path of the RINEX file.

    Returns
    -------
    bool
        True if compressed RINEX, False else.
    """
    boolout = bool(re.search(r".*((d|o)\.(Z)|(gz))$", rinex_path))
    return boolout


def crz2rnx(rinex_path, outdir="", force=True, path_of_crz2rnx="CRZ2RNX", verbose=True):
    """assuming that CRZ2RNX is in the system PATH per default"""

    if not os.path.isfile(rinex_path):
        raise Exception(rinex_path + " dont exists !")

    if force:
        forcearg = "-f"
    else:
        forcearg = ""

    curdir = os.getcwd()
    command = path_of_crz2rnx + " -c " + forcearg + " " + rinex_path
    if outdir == "":
        outdir = os.path.dirname(rinex_path)

    out_rinex_name_splited = os.path.basename(rinex_path).split(".")

    if conv.rinex_regex_search_tester(rinex_path, short_name=False, long_name=True):
        out_rinex_name_splited = os.path.basename(rinex_path).split(".")
        out_rinex_name = out_rinex_name_splited[0]
        out_rinex_name = out_rinex_name + ".rnx"
        out_rinex_path = os.path.join(outdir, out_rinex_name)
    else:
        out_rinex_name_splited = os.path.basename(rinex_path).split(".")
        out_rinex_name = ".".join(out_rinex_name_splited[:-1])
        out_rinex_name = out_rinex_name[:-1] + "o"
        out_rinex_path = os.path.join(outdir, out_rinex_name)

    if os.path.isfile(out_rinex_path) and not force:
        log.warning(out_rinex_path + "already exists, skiping ...")
    else:
        os.chdir(outdir)
        # stream = os.popen(command)

        proc = subprocess.Popen(
            command, shell=True, stdout=subprocess.PIPE, executable="/bin/bash"
        )
        status = proc.wait()

        if status != 0 or not os.path.isfile(out_rinex_path):
            log.error(out_rinex_path + " not created!")
        elif verbose:
            log.info(out_rinex_path + " created")
        else:
            pass
        #    print stream.read()
        os.chdir(curdir)
    return out_rinex_path


def rnx2crz(rinex_path, outdir="", force=True, path_of_rnx2crz="RNX2CRZ"):
    """assuming that RNX2CRZ is in the system PATH per default"""

    if not os.path.isfile(rinex_path):
        raise Exception(rinex_path + " dont exists !")

    if force:
        forcearg = "-f"
    else:
        forcearg = ""

    curdir = os.curdir
    command = path_of_rnx2crz + " " + forcearg + " " + rinex_path
    if outdir == "":
        outdir = os.path.dirname(rinex_path)

    out_rinex_name_o = os.path.basename(rinex_path)
    out_rinex_name_d_z = out_rinex_name_o[:-1] + "d.Z"
    out_rinex_path = os.path.join(outdir, out_rinex_name_d_z)

    if os.path.isfile(out_rinex_path) and not force:
        log.info(out_rinex_path + "already exists, skiping ...")
    else:
        os.chdir(outdir)
        stream = os.popen(command)
        log.info(command + " output in " + outdir)
        #    print stream.read()
    os.chdir(curdir)
    return out_rinex_path


def rinex_read_epoch(
    input_rinex_path_or_string,
    interval_out=False,
    add_tzinfo=False,
    out_array=True,
    out_index=False,
):
    """
    Read the epochs contained in a RINEX File. Can handle RINEX 2 and 3

    Parameters
    ----------
    input_rinex_path_or_string : see below
        input RINEX.
        can be the path of a RINEX file as string or as Path object,
        or directly the RINEX content as a string, bytes, StringIO object or a
        list of lines

    interval_out : bool, optional
        output also the intervals. The default is False.

    add_tzinfo : bool, optional
        add timezone information in the datetime's Epoches.
        The default is False.

    out_array : bool, optional
        output results as array.
        The default is True.

    out_index : bool, optional
        output also the index of the epoch line.
        The default is False.

    Returns
    -------
    array or list
        the epochs in the RINEX.
    """

    ##161019 : dirty copier coller de rinex start end
    epochs_list = []
    rinex_60sec = False

    try:
        input_rinex_path_or_string = hatanaka.decompress(input_rinex_path_or_string)
    except:
        pass

    rnx_lines = utils.open_readlines_smart(input_rinex_path_or_string)

    index_list = []
    for iline, line in enumerate(rnx_lines):
        epoch_rnx2 = re.search("^ {1,2}([0-9]{1,2} * ){5}", line)
        epoch_rnx3 = re.search("^>", line)

        if epoch_rnx2:
            try:
                epochdt = read_rnx_epoch_line(line, rnx2=True)
                epochs_list.append(epochdt)
                index_list.append(iline)
            except:
                log.error(line)

        elif epoch_rnx3:
            try:
                epochdt = read_rnx_epoch_line(line, rnx2=False)
                epochs_list.append(epochdt)
                index_list.append(iline)
            except:
                log.error(line)

    if add_tzinfo:
        epochs_list = [e.replace(tzinfo=dateutil.tz.tzutc()) for e in epochs_list]

    if out_index:
        output = [epochs_list, index_list]
    else:
        output = epochs_list

    if out_array:
        output = np.array(output).T

    return output


def same_day_rinex_check(rinex1, rinex2):
    if os.path.basename(rinex1)[4:7] == os.path.basename(rinex2)[4:7]:
        return True
    else:
        return False


def rinex_start_end(
    input_rinex_path,
    interval_out=False,
    add_tzinfo=False,
    verbose=True,
    safety_mode=True,
):
    """
    Return the first and the last epoch of a RINEX file
    (based on the actual content of the file, not the header)

    Can handle RINEX 2 and 3

    Parameters
    ----------
    input_rinex_path : TYPE
        path of the rinex file.
        can be the path of a RINEX or directly
        the RINEX content as a string

    interval_out : bool, optional
        output also the intervals. The default is False.

    add_tzinfo : bool, optional
        add timezone information in the datetime's Epoches.
        The default is False.

    verbose : bool, optional
        verbose mode. The default is True.

    safety_mode : TYPE, optional
        if the epoch reading fails (e.g. in case of a compressed RINEX)
        activate a reading of the header and the file name as backup.
        The default is True.

    Returns
    -------
    first_epoch , last_epoch , [interval]
        First, last epochs and interval if asked.

    """

    # une liste d'epochs en début et fin de fichier
    # => en trouver le min et le max
    # NB : FAIRE UN FONCTION READ EPOCH A L'OCCAZ
    # NBsuite : c'est fait au 161018 mais par contre c'est un dirty copier coller

    epochs_list = []
    head = utils.head(input_rinex_path, 1500)
    epochs_list_head = rinex_read_epoch(
        head, interval_out=interval_out, add_tzinfo=add_tzinfo, out_array=False
    )

    tail = utils.tail(input_rinex_path, 1500)
    epochs_list_tail = rinex_read_epoch(
        tail, interval_out=interval_out, add_tzinfo=add_tzinfo, out_array=False
    )

    epochs_list = epochs_list_head + epochs_list_tail

    if len(epochs_list) > 0:
        first_epoch = np.min(epochs_list)
        last_epoch = np.max(epochs_list)
    else:
        ## here is the save mode where everything can goes wrong ;)
        first_epoch = conv.rinexname2dt(input_rinex_path)
        if first_epoch:
            last_epoch = first_epoch + dt.timedelta(seconds=85399)
        else:
            ### here the return is None, your file is corrupted
            last_epoch = None
            return first_epoch, last_epoch

        ##### this legacy block if for short name only
        # alphabet = list(string.ascii_lowercase)
        # if os.path.basename(input_rinex_path)[7] in alphabet:
        #    last_epoch = first_epoch + dt.timedelta(hours=1)
        # else:
        #    last_epoch = first_epoch + dt.timedelta(hours=24,seconds=-1)

    if add_tzinfo:
        first_epoch = first_epoch.replace(tzinfo=dateutil.tz.tzutc())
        last_epoch = last_epoch.replace(tzinfo=dateutil.tz.tzutc())

    if verbose:
        log.debug("first & last epochs: %s %s", first_epoch, last_epoch)

    if not interval_out:
        return first_epoch, last_epoch
    else:
        interv_lis = np.diff(epochs_list)
        interv_lis = [e.seconds + e.microseconds * 10**-6 for e in interv_lis]
        interval = utils.most_common(interv_lis)
        log.debug("interval, last epoch: %s %s", interval, last_epoch)

        # return interv_lis , epochs_list
        return first_epoch, last_epoch, interval


def rinex_session_id(first_epoch, last_epoch, full_mode=False):
    """
    full_mode:
        gives the letter of the starting session  & the length in hour
        of the session
    """

    alphabet = list(string.ascii_lowercase)
    rinex_length = last_epoch - first_epoch
    if rinex_length.seconds > 86000:
        rnx_interval_ext = "0"
    else:
        if not full_mode:
            rnx_interval_ext = alphabet[first_epoch.hour]
        else:
            rnx_interval_ext = alphabet[first_epoch.hour] + str(
                int(rinex_length.seconds / 3600.0)
            )

    return rnx_interval_ext


def rinex_spliter_gfzrnx(
    input_rinex_path,
    output_directory,
    stat_out_name="",
    interval_size=86400,
    shift=0,
    inclusive=False,
    gfzrnx_cmd="GFZRNX",
    output_name="::RX3::",
    custom_cmds="",
):
    interval_size_hour = interval_size / 3600.0

    if stat_out_name == "":
        stat_out_name = os.path.basename(input_rinex_path)[0:4]

    # check if the RINEX is compressed ...
    bool_comp_rnx = check_if_compressed_rinex(input_rinex_path)

    # ... if not crz2rnx !
    if bool_comp_rnx:
        input_rinex_path = crz2rnx(input_rinex_path)

    inp_rinex_obj = open(input_rinex_path, "r+")
    out_dir = output_directory  # nom redondant mais j'ai la flemme d'aller corriger le nom de la variable

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    os.chdir(out_dir)

    first_epoch, last_epoch = rinex_start_end(input_rinex_path)

    # In this function, Date are truncated epochs
    # (only the day if interval == 24 , + the hour else)
    first_date = dt.datetime(
        first_epoch.year, first_epoch.month, first_epoch.day, first_epoch.hour
    )
    last_date = dt.datetime(
        last_epoch.year, last_epoch.month, last_epoch.day, last_epoch.hour
    ) + dt.timedelta(hours=1)

    log.info("first & last dates (truncated): %s %s", first_date, last_date)

    time_interval = utils.get_interval(
        first_date, last_date, dt.timedelta(hours=interval_size_hour)
    )

    alphabet = list(string.ascii_lowercase)

    rinex_out_name_lis = []

    for i, curr_date in enumerate(time_interval):

        if shift != 0:
            log.warning("shifted mode on, be careful !!!")

        if bool(shift) and i != 0:
            curr_date = curr_date + dt.timedelta(seconds=shift)
            interval_size_ope = interval_size
        # elif bool(shift) and i == 0 :
        #    interval_size_ope = interval_size + float(shift) / 60.
        else:
            interval_size_ope = interval_size

        if not inclusive:
            interval_size_ope = interval_size_ope - 1.0

        interval_size_ope = int(np.floor(interval_size_ope))

        if not bool(shift):
            if interval_size == 24:
                rnx_interval_ext = "0."
            else:
                rnx_interval_ext = alphabet[curr_date.hour] + "."
        else:
            if interval_size == 24:
                rnx_interval_ext = "0."
            else:
                rnx_interval_ext = alphabet[curr_date.hour] + "."

        p = subprocess.Popen(
            "",
            executable="/bin/bash",
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        # - 1/3600. # one_sec2h
        command = (
            gfzrnx_cmd
            + " -finp "
            + input_rinex_path
            + " -fout "
            + os.path.join(output_directory, output_name)
            + " -site "
            + stat_out_name.upper()
            + " -epo_beg "
            + curr_date.strftime("%Y%m%d_%H%M%S")
            + " --duration "
            + str(interval_size_ope)
            + " "
            + custom_cmds
        )
        log.info(command)

        # rinex_out_name = stat_out_name + curr_date.strftime('%j') + rnx_interval_ext + curr_date.strftime('%y') + 'o'
        std_log_name = curr_date.strftime("%j") + rnx_interval_ext + "err.log"
        err_log_name = curr_date.strftime("%j") + rnx_interval_ext + "err.log"

        stdout, stderr = p.communicate(command.encode())

        std_file = open(std_log_name, "w")
        std_file.write(stdout.decode("utf-8"))
        std_file.close()

        if stderr:
            log.info(err_log_name + " output:")
            log.info(stderr.decode("utf-8"))
            err_file = open(err_log_name, "w")
            err_file.write(stderr.decode("utf-8"))
            err_file.close()

        ### get the latest files
        list_of_files = glob.glob(
            output_directory + "/*"
        )  # * means all if need specific format then *.csv
        latest_file = sorted(list_of_files, key=os.path.getctime)[
            -3:
        ]  ### -2 and -1 are the logs

        latest_file = [e for e in latest_file if ".rnx" in e]

        if len(latest_file) == 0:
            latest_file = ""
            rnx_splt_path = ""
        else:
            latest_file = latest_file[0]
            rnx_splt_path = os.path.join(out_dir, latest_file)

        rinex_out_name_lis.append(rnx_splt_path)


# def rinex_renamer_gfz_odc_2_igs(rinex_in,output_dir=None):
