#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.operational contains functions to send runs 
to the GFZ's cluster in a batch mode. 

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
import getpass
#### Import the logger
import logging
import os
import re
import subprocess
import time

import numpy as np

#### geodeZYX modules
from geodezyx import utils

#### Import star style
# from geodezyx import *                   # Import the GeodeZYX modules
# from geodezyx.externlib import *         # Import the external modules
# from geodezyx.megalib.megalib import *   # Import the legacy modules names
log = logging.getLogger('geodezyx')

##########  END IMPORT  ##########
def cluster_GFZ_run(
    commands_list,
    bunch_on_off=True,
    bunch_job_nbr=10,
    bunch_wait_time=600,
    bj_check_on_off=True,
    bj_check_mini_nbr=2,
    bj_check_wait_time=120,
    bj_check_user="auto",
    add_cjob_cmd_prefix=True,
    wait_sleeping_before_launch=5,
):
    """
    Parameters
    ----------
    commands_list : list of str
        List of commands you want to run (one command per job).
    bunch_on_off : bool, optional
        If False, send all the jobs to the cluster at once.
        If True, just send <bunch_job_nbr> of them.
        See the options below.
        The default is True.
    bunch_job_nbr : int, optional
        number of jobs which will be sent to the cluser at once.
        The default is 10.
    bunch_wait_time : int, optional
        Minimal time between two bunch runs in sec.
        The default is 600.
    bj_check_on_off : bool, optional
        Do a check before a new bunch run, if previous jobs are still running.
        The default is True.
    bj_check_mini_nbr : int, optional
        The fuctions will wait <bj_check_wait_time> sec
        more if <bj_check_mini_nbr> or more jobs are running.
        The default is 2.
    bj_check_wait_time : int, optional
        wait time between two checks in sec. The default is 120.
    bj_check_user : str, optional
        username you want to check in the batchjobs.
        The default is "auto" i.e. your own username.
    add_cjob_cmd_prefix : bool, optional
        If, the input commands in command_list do not contain
        the cjob prefix command, it will be added automatically
        The default is True.
    wait_sleeping_before_launch : int, optional
        Waiting time between two successive runs.
        (To cancel the run if necessary)
        The default is 5.

    Returns
    -------
    None.

    """

    username = getpass.getuser()

    history_file_path = None

    i_bunch = 0

    if bj_check_user == "auto":
        bj_check_user = utils.get_username()

    log_path = "/home/" + bj_check_user + "/test_tmp.log"
    LOGobj = open(log_path, "w+")

    if not add_cjob_cmd_prefix:
        commands_list_opera = commands_list
    else:
        commands_list_opera = []
        for cmd in commands_list:
            cmd_ope = "cjob -c '" + cmd + "'"
            commands_list_opera.append(cmd_ope)

    log.info("****** JOBS THAT WILL BE LAUNCHED ******")
    log.info("Number of jobs : " + str(len(commands_list_opera)))
    log.info("****************************************")

    for ikommand, kommand in enumerate(commands_list_opera):

        ########## LOG/PRINT command
        log.info(kommand)
        LOGobj.write(kommand + "\n")

        ########## LOG/PRINT sleep
        info_sleep = (
            "script is sleeping for "
            + str(wait_sleeping_before_launch)
            + "sec (so you can cancel it) "
        )
        log.info(info_sleep)
        LOGobj.write(info_sleep + "\n")
        time.sleep(wait_sleeping_before_launch)

        ########## LOG/PRINT start
        info_start = "script starts @ " + str(dt.datetime.now())
        log.info(info_start)
        LOGobj.write(info_start + "\n")

        ########## RUN command here !!
        if True:
            p = subprocess.Popen(
                "",
                executable="/bin/csh",
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            stdout, stderr = p.communicate(kommand.encode())
        else:
            os.system(kommand)

        if history_file_path:
            with open(history_file_path, "a") as myfile:
                myfile.write(kommand + "\n")

        i_bunch += 1

        ########## Bunch On/Off : check if a bunch of job has been launched
        if bunch_on_off and np.mod(i_bunch, bunch_job_nbr) == 0:
            info_bunch = (
                "sleeping @ "
                + str(dt.datetime.now())
                + " for "
                + str(bunch_wait_time)
                + "s b.c. a bunch of "
                + str(bunch_job_nbr)
                + " jobs has been launched"
            )
            log.info(info_bunch)
            time.sleep(bunch_wait_time)
            LOGobj.write(info_bunch + "\n")

            ########## Bunch On/Off Check : check if the bunch is finished (experimental but on is better)
            ###### THIS PART MUST BE MERGE WITH THE SMALLER FCT BELLOW
            if log.info:
                log.info(
                    "All jobs should be finished now, let's see if there is some latecomers"
                )
                bj_check_tigger = False
            while not bj_check_tigger:
                bj_check_tigger = sleep_job_user(
                    minjob=bj_check_mini_nbr, bj_check_wait_time=bj_check_wait_time
                )

    return None


def number_job_user(bj_check_user=None, verbose=True):
    """
    Internal function for sleep_job_user
    """
    username = getpass.getuser()

    if not bj_check_user:
        bj_check_user = utils.get_username()

    bj_command = (
        "perl /dsk/igs2/soft_wrk/"
        + username
        + "/SOFT_EPOS8_BIN_TOOLS/SCRIPTS/e8_bjobs_local.pl"
    )

    bj_list = subprocess.check_output(bj_command, shell="/bin/csh")
    bj_list = bj_list.decode("utf-8")
    bj_list = bj_list.split("\n")

    bj_pattern_checked = bj_check_user + " *" + bj_check_user

    bj_list_checked = [bool(re.search(bj_pattern_checked, l)) for l in bj_list]
    bj_list_checked_sum = np.sum(bj_list_checked)

    if verbose:
        log.info(bj_list_checked_sum, "running jobs found for", bj_check_user)

    return bj_list_checked_sum, bj_list_checked, bj_pattern_checked


def sleep_job_user(bj_check_user=None, minjob=20, bj_check_wait_time=20):
    """
    Internal function for cluster_GFZ_run
    """
    if not bj_check_user:
        bj_check_user = utils.get_username()

    n_job, bj_list_checked, bj_pattern_checked = number_job_user(bj_check_user)

    if n_job >= minjob:
        check_tigger = False
        log.info(
            "sleeping @ "
            + str(dt.datetime.now())
            + " for "
            + str(bj_check_wait_time)
            + "s b.c."
            + str(n_job)
            + " jobs are runing"
        )
        time.sleep(bj_check_wait_time)
    else:
        check_tigger = True
        log.info("let's continue, no job matchs the pattern " + bj_pattern_checked)

    return check_tigger
