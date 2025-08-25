#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 26/02/2025 11:07:24

@author: psakic
"""

#### Import the logger
import os
import time
import datetime as dt
import subprocess
import numpy as np

#### Import geodezyx GINS submodules
import geodezyx.operational.gins_runner.gins_common as gynscmn

#### geodeZYX modules
from geodezyx import utils

#### Import the logger
import logging

log = logging.getLogger("geodezyx")


def run_directors(
    dir_paths_inp,
    opts_gins_pc="",
    opts_gins_90="",
    version="OPERA",
    cmd_mode="exe_gins_dir",
    force=False,
    verbose=True,
    sleep_time_max=10,
):
    """
    Executes GINS commands for a single or multiple directories, handling different modes and retry logic.

    Parameters
    ----------
    dir_paths_inp : str or list
        Path(s) to the directory/directories to process. Can be a single directory (str) or a list of directories (list).
    opts_gins_pc : str, optional
        Options for the GINS PC command. Default is an empty string.
    opts_gins_90 : str, optional
        Options for the GINS 90 command. Default is an empty string.
    version : str, optional
        Version of the GINS software to use. Default is "OPERA".
    cmd_mode : str, optional
        Mode of execution. Can be "ginspc", "exe_gins_dir", or "exe_gins_fic". Default is "exe_gins_dir".
    force : bool, optional
        If True, forces execution even if solutions already exist. Default is False.
    verbose : bool, optional
        If True, enables verbose logging. Default is True.
    sleep_time_max : int, optional
        Maximum sleep time (in seconds) to avoid conflicts during parallel execution. Default is 10.

    Returns
    -------
    None
        The function does not return any value.

    Notes
    -----
    - Handles both single and multi-directory modes.
    - Automatically detects `.fic` files and adjusts the command mode accordingly.
    - Implements retry logic for failed executions.
    - Logs execution details, including start and end times, and handles parallel execution conflicts.
    """
    # Multi or Single Mode ?
    if type(dir_paths_inp) is list:
        multimode = True
        director_path_lis = dir_paths_inp
        log.info("******** DIRECTORS RUNS *******")
    elif type(dir_paths_inp) is str:
        multimode = False
        director_path_lis = [dir_paths_inp]
    else:
        log.error("check the rinex_paths_in !!!")
        return None

    director_path_lis = list(sorted(director_path_lis))

    ndir = len(director_path_lis)

    for i, dir_path in enumerate(director_path_lis):
        dir_nam = os.path.basename(dir_path)

        if multimode:
            log.info(f"======== {i + 1} / {ndir} ======== ")
        log.info("start %s at %s", dir_nam, dt.datetime.now())
        if dir_path[-4:] == ".fic":
            cmd_mode = "exe_gins_fic"
            log.debug("input file ends with .fic, fic_mode is activated")
        start = time.time()

        if dir_path[-1] == "~":
            log.warning("geany ~ temp file, skiping this dir. ")
            continue

        sols_exist = gynscmn.check_solution(os.path.basename(dir_path), verbose=verbose)
        if len(sols_exist) > 0 and not force:
            log.info("solution %s already exists, skipping ...", sols_exist)
            continue

        ogpc = "-F" + opts_gins_pc
        og90 = opts_gins_90
        vers = version

        # log.info("options ginsPC: %s / gins90: %s", ogpc, og90)

        if "IPPP" in og90 and cmd_mode != "ginspc":
            _check_dir_keys(dir_path)

        ### cmd must be a string here... (and not a list of small str...)
        if cmd_mode == "ginspc":
            cmd = " ".join(("ginspc.bash", ogpc, dir_nam, og90, "-v", vers, "-f"))
        elif cmd_mode == "exe_gins_fic":
            cmd = " ".join(("exe_gins", "-fic", dir_nam, "-v", vers, og90))
        elif cmd_mode == "exe_gins_dir":
            cmd = " ".join(("exe_gins", "-dir", dir_nam, "-v", vers, og90))
        else:
            log.error("mode not recognized !!!")
            return None

        gins_path = gynscmn.get_gin_path()
        log_path = utils.create_dir(os.path.join(gins_path, "python_logs"))
        log_path = os.path.join(gins_path, "python_logs", dir_nam + ".log")

        ## artifical sleep to avoid conflict when parallelizing
        if sleep_time_max > 0:
            time.sleep(np.random.randint(1, sleep_time_max * 1000) * 10**-3)

        ## parameters for retry
        itry = 0
        itry_max = 2
        retry_str = ""
        timeout = 300

        while itry <= itry_max:
            itry += 1

            if verbose:
                log.info("submit. command : %s", cmd)

            process = subprocess.run(
                [cmd],
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                executable="/bin/bash",
                timeout=timeout
            )

            # gynscmn.check_gins_exe(open(log_path), director_name)
            #
            # with open(log_path + ".exec", "w") as f:
            #     f.write("exec time : " + str(time.time() - start))
            #
            ok_sol = gynscmn.check_solution(dir_nam, verbose=verbose)
            if ok_sol:
                log.info("solution ok for: %s :) %s", dir_nam, retry_str)
            else:
                log.error("no solution for: %s :( %s", dir_nam, retry_str)
                #log.error("STDOUT: %s", process.stdout)
                #log.error("STDERR: %s", process.stderr)

            if not ok_sol and check_for_retry(process.stderr) and itry < itry_max:
                log.warning(f"retryable directeur {dir_nam} (try {itry} / {itry_max})")
                itry = itry + 1
                retry_str = f"(try {itry} / {itry_max})"
                time.sleep(3)
            else:
                itry = itry_max + 1

        log.info(
            "run %s end at %s (exec: %8.3f s)",
            dir_nam,
            dt.datetime.now(),
            time.time() - start,
        )

    return None


def check_for_retry(stderr_inp):
    if "exe_gins_recup_fic.sh: 67: [: =: unexpected operator" in stderr_inp:
        return True
    else:
        return False

def _check_dir_keys(director_path_inp):
    for grepstr in (
        "userext_gps__qualiteorb",
        "userext_gps__haute_freq",
        "userext_gps__hor_interp",
        "GPS__QUALITEORB",
        "GPS__HAUTE_FREQ",
        "GPS__HOR_INTERP",
    ):
        grep_out = utils.grep(director_path_inp, grepstr)
        if grep_out == "":
            log.warning("IPPP mode on, but no %sin the dir !!!", grepstr)


#######################################################################################################

#  ______                _   _                                                             _
# |  ____|              | | (_)                                                           | |
# | |__ _   _ _ __   ___| |_ _  ___  _ __     __ _ _ __ __ ___   _____ _   _  __ _ _ __ __| |
# |  __| | | | '_ \ / __| __| |/ _ \| '_ \   / _` | '__/ _` \ \ / / _ \ | | |/ _` | '__/ _` |
# | |  | |_| | | | | (__| |_| | (_) | | | | | (_| | | | (_| |\ v /  __/ |_| | (_| | | | (_| |
# |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|  \__, |_|  \__,_| \_/ \___|\__, |\__,_|_|  \__,_|
#                                             __/ |                     __/ |
#                                            |___/                     |___/



########## OLD FUNCTIONS ##########


