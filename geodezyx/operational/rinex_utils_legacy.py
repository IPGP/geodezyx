#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 23/09/2025 12:17:54

@author: psakic

This module contains some utilities to handle RINEX files
in legacy mode i.e. RINEX2 using teqc, crz2rnx, rnx2crz

"""

import os
import re
import pandas as pd
import datetime as dt
import subprocess
import numpy as np
import geodezyx.utils as utils
import geodezyx.conv as conv
import string
import logging
import shutil

log = logging.getLogger("geodezyx")


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

    blocs_stack = []

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
            line_bloc = []
            lineconcat = ""
            date = read_rnx_epoch_line(l, rnx2=True)

        ### we read the sat lines based on the number of sat
        if iline_bloc <= nlines_bloc:
            line_bloc.append(l[32:].strip())
            lineconcat = lineconcat + l[32:].strip()
            iline_bloc += 1

        ### we stack everything when the sat block is over
        if iline_bloc == nlines_bloc:
            in_epoch = False
            bloc_tuple = (date, nsat, lineconcat)
            blocs_stack.append(bloc_tuple)

    #### we do a DF with all the found line
    df = pd.DataFrame(blocs_stack)

    df[3] = None

    for irow, row in df.iterrows():
        nsat = str(row[1])
        satstr = row[2]
        re_sat = "([A-Z][0-9][0-9]){" + nsat + "}"
        bool_sat = re.search(re_sat, satstr)

        #### the the right number of sat is found True = Good
        if bool_sat:
            df.iloc[irow, 3] = True
        else:
            df.iloc[irow, 3] = False

    if df[3].sum() != len(df):
        return rinex_sats_checker
    else:
        return None


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