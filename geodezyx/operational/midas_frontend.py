#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.operational contains functions to run the 
time series velocities estimation software MIDAS. 

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
import glob
#### Import the logger
import logging
import os
import shutil
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#### geodeZYX modules
from geodezyx import conv
from geodezyx import stats
from geodezyx import utils

log = logging.getLogger(__name__)

##########  END IMPORT  ##########


#
#  __  __ _____ _____           _____
# |  \/  |_   _|  __ \   /\    / ____|
# | \  / | | | | |  | | /  \  | (___
# | |\/| | | | | |  | |/ /\ \  \___ \
# | |  | |_| |_| |__| / ____ \ ____) |
# |_|  |_|_____|_____/_/    \_\_____/
#
#


def midas_run(
    tenu_file_path,
    work_dir="",
    path_midas_soft="",
    step_file_path="",
    with_plot=True,
    keep_plot_open=True,
):
    ### Prepare paths of facultative arguments
    if not work_dir:
        work_dir = os.path.dirname(tenu_file_path)
    if not path_midas_soft:
        path_midas_soft = "midas.e"

    ### SECURITY : rm all previous MIDAS files
    MIDAS_files_rm = [
        "MIDAS.STEPIN",
        "MIDAS.STEPOUT",
        "MIDAS.ERR",
        "MIDAS.TENV",
        "MIDAS.RENV",
        "MIDAS.TENU",
        "MIDAS.RENU",
        "MIDAS.VEL",
        "fort.91",
    ]

    for fil in MIDAS_files_rm:
        fil_full_path = os.path.join(work_dir, fil)
        if os.path.isfile(fil_full_path):
            log.info("old %s removed", fil)
            os.remove(fil_full_path)

    ### Find the Station Name
    DF = pd.read_table(tenu_file_path, comment="*", header=-1, delim_whitespace=True)
    stat = list(set(DF[0]))[0]

    ### Prepare the tmp input
    ### The file is called MIDAS.TENU in any case
    os.chdir(work_dir)
    work_file_path_tenu = work_dir + "/MIDAS.TENU"
    shutil.copy(tenu_file_path, work_file_path_tenu)

    ### Prepare the tmp STEP FILE input
    step_OK = False
    if step_file_path and not utils.empty_file_check(step_file_path):
        work_file_path_step = work_dir + "/MIDAS.STEP"
        shutil.copy(step_file_path, work_file_path_step)
        step_OK = True
    else:
        work_file_path_step = ""

    ### Prepare the name of the outputed file
    work_file_path_vel = os.path.join(work_dir, "MIDAS.VEL")

    ### LETS RUN !
    command = path_midas_soft
    log.info("LAUNCHING : %s for stat. %s", command, stat)
    p = subprocess.Popen(
        "",
        executable="/bin/bash",
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    stdout, stderr = p.communicate(command.encode())
    # logs files
    vel_file = os.path.join(work_dir, "MIDAS.VEL")
    vel_file_obj = open(vel_file, "w")
    vel_file_obj.write(stdout.decode("utf-8"))
    vel_file_obj.close()

    exec_OK = False

    if not utils.empty_file_check(work_file_path_vel):
        exec_OK = True
        log.info("MIDAS.VEL exists for %s, should be OK :)", stat)

    if not utils.empty_file_check(os.path.join(work_dir, "MIDAS.ERR")):
        log.error("MIDAS.ERR is not empty for, must be checked :(", stat)

    if stderr:
        log.warning("err.log is not empty, must be checked !")
        log.warning(stderr.decode("utf-8"))
        log.warning("")
        err_file = os.path.join(work_dir, stat + ".err.log")
        err_file_obj = open(err_file, "w")
        err_file_obj.write(stderr.decode("utf-8"))
        err_file_obj.close()

    #### PLOT
    if exec_OK and with_plot:
        fig = midas_plot(work_file_path_tenu, work_file_path_vel, work_file_path_step)
        work_file_path_plot = os.path.join(work_dir, "plot." + stat.upper())
        for ext in (".pdf", ".png"):
            plt.savefig(work_file_path_plot + ext)
        if not keep_plot_open:
            plt.close(fig)

    #### RENAME ALL MIDAS FILES WITH A STATION PREFIX
    for fname in glob.glob("*MIDAS*"):
        os.rename(fname, fname.lower() + "_" + stat.upper())

    return None


def midas_vel_files_2_pandas_DF(vel_files_in):
    """
    Convert MIDAS Velocity files to a Pandas DataFrame

    Parameters
    ----------
    vel_files_in : str or list of str
        if list of str, will consider directly the files inside the list
        if str, can be the path of a single file, or a generic (wilcard) path
        to consider several files
    Returns
    -------
    DF : Pandas DataFrame
    """

    if utils.is_iterable(vel_files_in):
        L_vel_files = vel_files_in
    else:
        L_vel_files = glob.glob(vel_files_in)

    if len(L_vel_files) > 1:
        work_dir = os.path.dirname(L_vel_files[0])
        vel_file_opera = utils.cat(work_dir + "/MERGED_VEL", *L_vel_files)
    else:
        vel_file_opera = L_vel_files[0]

    DF = pd.read_table(vel_file_opera, comment="*", header=-1, delim_whitespace=True)

    DF.columns = [
        "Station",
        "soft",
        "epoch_first",
        "epoch_last",
        "duration",
        "nb_epoch_all",
        "nb_epoch_good",
        "nb_pairs",
        "V_East",
        "V_North",
        "V_Up",
        "sV_East",
        "sV_North",
        "sV_Up",
        "offset_e_1st_epoch",
        "offset_n_1st_epoch",
        "offset_u_1st_epoch",
        "outlier_ratio_e",
        "outlier_ratio_n",
        "outlier_ratio_u",
        "std_velo_pair_e",
        "std_velo_pair_n",
        "std_velo_pair_u",
        "nb_steps",
    ]

    return DF


def midas_plot(path_tenu, path_vel="", path_step=""):
    """
    based on a plot for TimeSeriePoint
    """
    import geoclass as gcls

    DF = pd.read_table(path_tenu, header=-1, delim_whitespace=True)
    stat = list(set(DF[0]))[0]

    if path_vel:
        DFvel = midas_vel_files_2_pandas_DF(path_vel)
        plot_vel = True
    else:
        plot_vel = False

    if path_step:
        DFstep = pd.read_table(path_step, header=-1, delim_whitespace=True)
        Discont = conv.year_decimal2dt(list(DFstep[1]))

    T = DF[1]
    Tdt = conv.year_decimal2dt(T)
    E = DF[2]
    N = DF[3]
    U = DF[4]

    TS = gcls.TimeSeriePoint()

    TS.from_list(Tdt, E, N, U, "ENU")
    TS.meta_set(path_tenu, stat, stat)

    if path_step:
        TS.set_discont(Discont)
        log.info("Discontinuity", Discont)

    #### PLOT ########
    TSplot_fig = TS.plot(fig=plt.figure())
    if path_step:
        TS.plot_discont(fig=TSplot_fig)
    ##################
    TSplot_axe = TSplot_fig.axes

    a_keys_list = ["V_East", "V_North", "V_Up"]

    b_keys_list = ["offset_e_1st_epoch", "offset_n_1st_epoch", "offset_u_1st_epoch"]

    sigma_key_list = ["sV_East", "sV_North", "sV_Up"]

    if plot_vel:
        for i, (a_key, b_key, sigma_key) in enumerate(
            zip(a_keys_list, b_keys_list, sigma_key_list)
        ):
            ax = TSplot_axe[i]

            a = DFvel[a_key][0]
            b = DFvel[b_key][0]
            sigma = DFvel[sigma_key][0]

            _, Vlin = stats.linear_reg_getvalue(T - T[0], a, b)

            ax.plot(Tdt, Vlin, c="xkcd:dark orange")
            text_vel = (
                str(np.round(a * 1000, 4))
                + " +/- "
                + str(np.round(sigma * 1000, 4))
                + " mm/yr"
            )

            ax.text(
                0.8,
                0.1,
                text_vel,
                horizontalalignment="center",
                verticalalignment="center",
                transform=ax.transAxes,
            )

    return TSplot_fig
