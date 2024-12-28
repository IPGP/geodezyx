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

log = logging.getLogger(__name__)

##########  END IMPORT  ##########


def rinexs_table_from_list(rnxs_inp, site9_col=False, round_date=False, path_col=True):
    """
    From a simple RINEX list, summarize the data in an ad-hoc DataFrame

    Parameters
    ----------
    rnxs_inp : iterable or str
        if iterable, a list of RINEX files
        if str, path of an input RINEX list as text file.

    Returns
    -------
    DF : Pandas DataFrame
        A DataFrame with the RINEX info in it.
    """

    if utils.is_iterable(rnxs_inp):
        DF = pd.DataFrame(rnxs_inp, columns=["path"])
    else:
        DF = pd.read_csv(rnxs_inp, names=["path"])

    DF["name"] = DF["path"].apply(os.path.basename)
    DF["site4"] = DF["name"].str[:4].str.lower()

    if site9_col:
        DF["site9"] = DF["name"].str[:9].str.lower()
        #### set generic 9-char name if 4 char name file
        BoolNewName = DF["name"].str.match(conv.rinex_regex_long_name())
        BoolOldName = np.logical_not(BoolNewName)
        Site9_generic = DF["site4"].str.upper() + "00XXX"
        DF["site9"].loc[BoolOldName] = Site9_generic[BoolOldName]

    DF["date"] = DF["name"].apply(conv.rinexname2dt)
    if round_date:
        DF["date"] = conv.round_dt(DF["date"].values, "1D", True, "floor")
    DoyYear = DF["date"].apply(conv.dt2doy_year)
    DF["doy"] = DoyYear.apply(lambda x: x[0])
    DF["year"] = DoyYear.apply(lambda x: x[1])

    cols = DF.columns.tolist()

    cols = cols[1:] + [cols[0]]

    DF = DF[cols]

    if not path_col:
        DF.drop("path", axis=1, inplace=True)

    pd.options.display.max_info_columns = 150

    return DF


def rinex_sats_checker(p_rnx):
    """
    Check the consistency of a RINEX's sat list for each epoch

    Designed for checking the consistency of IPGP-OVPF corrupted RINEXs
    RINEX-2 only

    Parameters
    ----------
    p_rnx : str
        input RINEX path.

    Returns
    -------
    bad RINEX path or None.

    """

    import hatanaka

    ### This bloc must  be replaced by open_readlines_smart
    if p_rnx.endswith(".Z") or p_rnx.endswith(".gz"):
        O = hatanaka.decompress(p_rnx)
        O = str(O)
        L = O.split("\\n")
    else:
        F = open(p_rnx, "r+")
        L = F.readlines()

    iline_bloc = 0
    nlines_bloc = -1

    Blocs_stack = []

    ###### This loop read all lines
    for l in L:
        re_epoch = "^ {1,2}([0-9]{1,2} * ){5}"
        # re_sat="[A-Z][0-9][0-9]"
        bool_epoch = re.search(re_epoch, l)
        # bool_sat=re.search(re_sat,l)

        ### we found an epoch line
        if bool_epoch:
            in_epoch = True
            nsat = int(l[30:32])
            iline_bloc = 0
            nlines_bloc = int(np.ceil(nsat / 12))
            LineBloc = []
            lineconcat = ""
            date = read_rnx_epoch_line(l, rnx2=True)

        ### we read the sat lines based on the number of sat
        if iline_bloc <= nlines_bloc:
            LineBloc.append(l[32:].strip())
            lineconcat = lineconcat + l[32:].strip()
            iline_bloc += 1

        ### we stack everything when the sat block is over
        if iline_bloc == nlines_bloc:
            in_epoch = False
            bloc_tuple = (date, nsat, lineconcat)
            Blocs_stack.append(bloc_tuple)

    #### we do a DF with all the found line
    DF = pd.DataFrame(Blocs_stack)

    DF[3] = None

    for irow, row in DF.iterrows():
        nsat = str(row[1])
        satstr = row[2]
        re_sat = "([A-Z][0-9][0-9]){" + nsat + "}"
        bool_sat = re.search(re_sat, satstr)

        #### the the right number of sat is found True = Good
        if bool_sat:
            DF.iloc[irow, 3] = True
        else:
            DF.iloc[irow, 3] = False

    if DF[3].sum() != len(DF):
        return rinex_sats_checker
    else:
        return None


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
# | |__) | | | |  \| | |__   \ V /  | (___  _ __ | |_| |_ ___ _ __
# |  _  /  | | | . ` |  __|   > <    \___ \| '_ \| | | __/ _ \ '__|
# | | \ \ _| |_| |\  | |____ / . \   ____) | |_) | | | ||  __/ |
# |_|  \_\_____|_| \_|______/_/ \_\ |_____/| .__/|_|_|\__\___|_|
#                                          | |
#                                          |_|


def check_if_compressed_rinex(rinex_path):
    boolout = bool(re.search(r".*((d|o)\.(Z)|(gz))$", rinex_path))
    return boolout


def crz2rnx(rinex_path, outdir="", force=True, path_of_crz2rnx="CRZ2RNX"):
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
    out_rinex_name = ".".join(out_rinex_name_splited[:-1])
    out_rinex_name = out_rinex_name[:-1] + "o"
    out_rinex_path = os.path.join(outdir, out_rinex_name)

    if os.path.isfile(out_rinex_path) and not force:
        log.info(out_rinex_path + "already exists, skiping ...")
    else:
        os.chdir(outdir)
        # stream = os.popen(command)

        proc = subprocess.Popen(
            command, shell=True, stdout=subprocess.PIPE, executable="/bin/bash"
        )
        status = proc.wait()

        if status != 0 or not os.path.isfile(out_rinex_path):
            log.info(out_rinex_path + " not created ! :(")
        else:
            log.info(out_rinex_path + " created :)")
        #    print stream.read()
        os.chdir(curdir)
    return out_rinex_path


def crz2rnx_bad(crinex_in_path, outdir="", force=True, path_of_crz2rnx="CRZ2RNX"):
    """assuming that CRZ2RNX is in the system PATH per default"""

    if not os.path.isfile(crinex_in_path):
        raise Exception(crinex_in_path + " dont exists !")

    if force:
        forcearg = "-f"
    else:
        forcearg = ""

    command = path_of_crz2rnx + " " + forcearg + " " + crinex_in_path
    convert_dir = os.path.dirname(crinex_in_path)

    if outdir == "":
        outdir = convert_dir

    out_rinex_name_splited = os.path.basename(crinex_in_path).split(".")
    out_rinex_name = ".".join(out_rinex_name_splited[:-1])
    out_rinex_name = out_rinex_name[:-1] + "o"

    out_rinex_path = os.path.join(outdir, out_rinex_name)
    temp_rinex_path = os.path.join(convert_dir, out_rinex_name)

    if os.path.isfile(out_rinex_path) and not force:
        log.info(out_rinex_path + "already exists, skiping ...")

    else:
        stream = os.popen(command)
        log.info("output in " + outdir)
        try:
            shutil.move(temp_rinex_path, outdir)
        except:
            pass

        #    print stream.read()
    #    os.chdir(curdir)
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
    out_rinex_name_d_Z = out_rinex_name_o[:-1] + "d.Z"
    out_rinex_path = os.path.join(outdir, out_rinex_name_d_Z)

    if os.path.isfile(out_rinex_path) and not force:
        log.info(out_rinex_path + "already exists, skiping ...")
    else:
        os.chdir(outdir)
        stream = os.popen(command)
        log.info(command + " output in " + outdir)
        #    print stream.read()
    os.chdir(curdir)
    return out_rinex_path


# def rinex_regex(compressed=True,compiled=False):
#     if compressed:
#         regexstr = "....[0-9]{3}.\.[0-9]{2}((d\.(Z)|(gz))|(o)|(d))"
#     else:
#         regexstr = "....[0-9]{3}.\.[0-9]{2}o"

#     if compiled:
#         return re.compile(regexstr)
#     else:
#         return regexstr


# def rinex_regex_new_name(compressed=True,compiled=False):
#     if compressed:
#         regexstr = ".{4}[0-9]{2}.{3}_(R|S|U)_[0-9]{11}_[0-9]{2}\w_[0-9]{2}\w_\w{2}\.\w{3}\.gz"
#     else:
#         regexstr = ".{4}[0-9]{2}.{3}_(R|S|U)_[0-9]{11}_[0-9]{2}\w_[0-9]{2}\w_\w{2}\.\w{3}"

#     if compiled:
#         return re.compile(regexstr)
#     else:
#         return regexstr


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

    RnxLines = utils.open_readlines_smart(input_rinex_path_or_string)

    Index_list = []
    for iline, line in enumerate(RnxLines):
        epoch_rnx2 = re.search("^ {1,2}([0-9]{1,2} * ){5}", line)
        epoch_rnx3 = re.search("^>", line)

        if epoch_rnx2:
            try:
                epochdt = read_rnx_epoch_line(line, rnx2=True)
                epochs_list.append(epochdt)
                Index_list.append(iline)
            except:
                log.error(line)

        elif epoch_rnx3:
            try:
                epochdt = read_rnx_epoch_line(line, rnx2=False)
                epochs_list.append(epochdt)
                Index_list.append(iline)
            except:
                log.error(line)

    if add_tzinfo:
        epochs_list = [e.replace(tzinfo=dateutil.tz.tzutc()) for e in epochs_list]

    if out_index:
        OUTPUT = [epochs_list, Index_list]
    else:
        OUTPUT = epochs_list

    if out_array:
        OUTPUT = np.array(OUTPUT).T

    return OUTPUT


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
        First, las epoches and interval if asked.

    """

    # une liste d'epochs en début et fin de fichier
    # => en trouver le min et le max
    # NB : FAIRE UN FONCTION READ EPOCH A L'OCCAZ
    # NBsuite : c'est fait au 161018 mais par contre c'est un dirty copier coller

    epochs_list = []
    Head = utils.head(input_rinex_path, 1500)
    epochs_list_head = rinex_read_epoch(
        Head, interval_out=interval_out, add_tzinfo=add_tzinfo, out_array=False
    )

    Tail = utils.tail(input_rinex_path, 1500)
    epochs_list_tail = rinex_read_epoch(
        Tail, interval_out=interval_out, add_tzinfo=add_tzinfo, out_array=False
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


# def rinex_spliter(input_rinex_path,output_directory,stat_out_name='',
#                  interval_size=24,compress=False,shift = 0):
#    """
#    if shift != 0:
#    the start/end of a session is shifted of shift minutes
#    """
#
#    if stat_out_name == '':
#        stat_out_name = os.path.basename(input_rinex_path)[0:4]
#
#    # check if the RINEX is compressed ...
#    bool_comp_rnx = check_if_compressed_rinex(input_rinex_path)
#    # ... if not crz2rnx !
#    if bool_comp_rnx:
#        input_rinex_path = crz2rnx(input_rinex_path)
#
#    inp_rinex_obj=open(input_rinex_path,'r+')
#    out_dir = output_directory # nom redondant mais j'ai la flemme d'aller corriger le nom de la variable
#
#    if not os.path.exists(out_dir):
#        os.makedirs(out_dir)
#    os.chdir(out_dir)
#
#    first_epoch , last_epoch = rinex_start_end(input_rinex_path)
#
#    # In this function, Date are truncated epochs
#    # (only the day if interval == 24 , + the hour else)
#    first_date = datetime.datetime(first_epoch.year,first_epoch.month,
#                                   first_epoch.day,first_epoch.hour)
#    last_date = datetime.datetime(last_epoch.year,last_epoch.month,
#                                  last_epoch.day,last_epoch.hour)+datetime.timedelta(hours=1)
#
#    print "first & last dates (truncated) : " , first_date , last_date
#
#    time_interval = genefun.get_interval(first_date,last_date,
#                                         datetime.timedelta(hours=interval_size))
#
#    alphabet = list(string.ascii_lowercase)
#
#    rinex_out_name_lis = []
#
#    for i,curr_date in enumerate(time_interval):
#        if interval_size == 24:
#            rnx_interval_ext = '0.'
#        else:
#            rnx_interval_ext = alphabet[curr_date.hour] + '.'
#
#        p = subprocess.Popen('',executable='/bin/bash', stdin=subprocess.PIPE , stdout=subprocess.PIPE , stderr=subprocess.PIPE)
#        command = 'teqc ' + '-O.mo ' + stat_out_name.upper() +' -st ' +  curr_date.strftime('%y%m%d%H%M%S') + ' +dh ' + str(interval_size) + ' ' + input_rinex_path
#        print command
#        rinex_out_name = stat_out_name + curr_date.strftime('%j') + rnx_interval_ext + curr_date.strftime('%y') + 'o'
#        err_log_name = curr_date.strftime('%j') + rnx_interval_ext + "err.log"
#        print rinex_out_name
#        stdout,stderr = p.communicate( command )
#        std_file = open(rinex_out_name, "w")
#        std_file.write(stdout)
#        std_file.close()
#        if stderr != '':
#            print err_log_name + " is not empty, must be checked !"
#            print stderr
#            err_file = open(err_log_name, "w")
#            err_file.write(stderr)
#            err_file.close()
#
#        rnx_splt_path = os.path.join(out_dir,rinex_out_name)
#        if compress:
#            rinex_out_final = rnx2crz(rnx_splt_path)
#            os.remove(rnx_splt_path)
#        else:
#            rinex_out_final = rnx_splt_path
#
#        rinex_out_name_lis.append(rinex_out_final)
#
#    return rinex_out_name_lis


def rinex_spliter(
    input_rinex_path,
    output_directory,
    stat_out_name="",
    interval_size=24,
    compress=False,
    shift=0,
    inclusive=False,
    teqc_cmd="teqc",
):
    """
    if shift != 0:
    the start/end of a session is shifted of shift minutes

    inclusive:
    delta of exaclty interval_size => add the 1st epoch of the next sess
    not inclusive:
    delta of interval_size - 1s
    """
    """
    Split RINEX file using teqc

    Parameters
    ----------
    input_rinex_path : str
        Path of the input RINEX

    output_directory : str
        Description param2
        
    stat_out_name : str
        Station name for the output RINEX
        
    interval_size : int
        Size of the splitted interval
    
    compress : bool
        compress the outputed RINEX    

    interval_size : int
        Size of the splitted interval 
        
        
    FINIR CE HEADER !!!!
                
    Returns
    -------
    out1 : float or int or str or dict or n-tuple or bool or list or numpy.array
        Description out1
    
    out2 : float or int or str or dict or n-tuple or bool or list or numpy.array
        Description out2
        
    Note
    ----
    Misc. Notes

    Source
    ------
    www.forum-source.com
    
    Examples
    --------
    >> answer
    42    
    """

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
        first_date, last_date, dt.timedelta(hours=interval_size)
    )

    alphabet = list(string.ascii_lowercase)

    rinex_out_name_lis = []

    for i, curr_date in enumerate(time_interval):

        if shift != 0:
            log.warning("shifted mode on, be careful !!!")

        if bool(shift) and i != 0:
            curr_date = curr_date + dt.timedelta(minutes=shift)
            interval_size_ope = interval_size
        elif bool(shift) and i == 0:
            interval_size_ope = interval_size + float(shift) / 60.0
        else:
            interval_size_ope = interval_size

        if not inclusive:
            interval_size_ope = interval_size_ope - 1 / 3600.0

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
            teqc_cmd
            + " -O.mo "
            + stat_out_name.upper()
            + " -st "
            + curr_date.strftime("%y%m%d%H%M%S")
            + " +dh "
            + str(interval_size_ope)
            + " "
            + input_rinex_path
        )
        log.info(command)
        rinex_out_name = (
            stat_out_name
            + curr_date.strftime("%j")
            + rnx_interval_ext
            + curr_date.strftime("%y")
            + "o"
        )
        err_log_name = curr_date.strftime("%j") + rnx_interval_ext + "err.log"
        log.info(rinex_out_name)
        stdout, stderr = p.communicate(command.encode())
        std_file = open(rinex_out_name, "w")
        std_file.write(stdout.decode("utf-8"))
        std_file.close()
        if stderr:
            log.warning(err_log_name + " is not empty, must be checked !")
            log.info(stderr.decode("utf-8"))
            err_file = open(err_log_name, "w")
            err_file.write(stderr.decode("utf-8"))
            err_file.close()

        rnx_splt_path = os.path.join(out_dir, rinex_out_name)
        if compress:
            rinex_out_final = rnx2crz(rnx_splt_path)
            os.remove(rnx_splt_path)
        else:
            rinex_out_final = rnx_splt_path

        rinex_out_name_lis.append(rinex_out_final)

    return rinex_out_name_lis


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


def teqc_qc(rinex_path, quick_mode=False, optional_args=""):
    """
    quick mode : reduced qc and no summary file written
    """
    if quick_mode:
        qcmode = "+qcq"
    else:
        qcmode = "+qc"

    # check if the RINEX is compressed ...
    bool_comp_rnx = check_if_compressed_rinex(rinex_path)

    # ... if not crz2rnx !
    if bool_comp_rnx:
        rinex_path_work = crz2rnx(rinex_path)
    else:
        rinex_path_work = rinex_path

    kommand = " ".join(("teqc", qcmode, optional_args, rinex_path))

    log.info(kommand)

    proc = subprocess.Popen(
        kommand, shell=True, stdout=subprocess.PIPE, executable="/bin/bash"
    )
    status = proc.wait()

    if bool_comp_rnx:
        os.remove(rinex_path_work)

    if status:
        log.error("crash for " + rinex_path + ", code " + str(proc.poll()))
        return None

    output = proc.stdout.read()

    return output.decode("ASCII")


def rinex_renamer(input_rinex_path, output_directory, stat_out_name="", remove=False):

    if stat_out_name == "":
        stat_out_name = os.path.basename(input_rinex_path)[0:4]

    stat_out_name = stat_out_name.lower()

    inp_rinex_obj = open(input_rinex_path, "r+")
    out_dir = output_directory  # nom redondant mais j'ai la flemme d'aller corriger le nom de la variable

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    os.chdir(out_dir)

    first_epoch, last_epoch = rinex_start_end(input_rinex_path)
    rnx_interval_ext = rinex_session_id(first_epoch, last_epoch) + "."
    rinex_out_name = (
        stat_out_name
        + first_epoch.strftime("%j")
        + rnx_interval_ext
        + first_epoch.strftime("%y")
        + "o"
    )

    log.info(rinex_out_name)

    output_rinex_path = os.path.join(out_dir, rinex_out_name)

    if input_rinex_path != output_rinex_path:
        log.info("copy of " + input_rinex_path + " to " + output_rinex_path)
        shutil.copy(input_rinex_path, output_rinex_path)

        if remove and os.isfile(output_rinex_path):
            log.info("removing ", input_rinex_path)
            os.remove(input_rinex_path)
    else:
        log.info(input_rinex_path)
        log.info("and " + output_rinex_path + " are the same file")
        log.info("nothing's done ...")

    return output_rinex_path


# def rinex_renamer_gfz_odc_2_igs(rinex_in,output_dir=None):
