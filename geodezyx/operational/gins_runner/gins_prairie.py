#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 26/02/2025 11:09:30

@author: psakic
"""

########## BEGIN IMPORT ##########
#### External modules
#### Import the logger
import logging
import os
import subprocess

import geodezyx.operational.gins_runner.gins_common as ginscmn

#### geodeZYX modules
from geodezyx import operational

log = logging.getLogger('geodezyx')

def prairie_manual(
    rinex_paths_in,
    temp_data_folder="",
    force=True,
    with_historik=True,
    with_wsb=True,
    argsdict=dict(),
):
    """
    Process RINEX files using the GINS prairie manual method.

    This function processes one or more RINEX files through the GINS prairie
    executable, handling file management, decompression, and command execution
    automatically. It supports various prairie processing options and ensures
    proper file placement within the GINS environment.

    Parameters
    ----------
    rinex_paths_in : str or list of str
        Path(s) to the RINEX file(s) to process. Can be a single string path
        or a list of string paths for batch processing.
    temp_data_folder : str, optional
        Path to the temporary data folder where files will be processed.
        If empty string (default), uses the default GINS TEMP_DATA folder.
    force : bool, optional
        If True, forces processing even if the output file already exists.
        Default is True.
    with_historik : bool, optional
        If True, includes GLONASS historik processing using the default
        historik file. Default is True.
    with_wsb : bool, optional
        If True, includes WSB (satellite clock bias) processing using the
        default WSB reference file. Default is True.
    argsdict : dict, optional
        Additional arguments for the prairie command. Keys should be argument
        names (e.g., '-options'), and values should be their corresponding
        file paths or values. Arguments cc2noncc, wsb, options, and out are
        automatically managed if not provided.

    Returns
    -------
    str or list of str
        Path(s) to the output prairie file(s). Returns a single string if
        one file is processed, or a list of strings if multiple files are
        processed.

    Notes
    -----
    - Automatically copies RINEX files to the temporary data folder if needed
    - Decompresses compressed RINEX files (.crx, .gz) automatically
    - Manages default prairie command arguments for standard processing
    - Executes the prairie command using subprocess in the temp data folder
    - Validates output file creation and size before returning paths

    Examples
    --------
    >>> # Process a single RINEX file
    >>> output = prairie_manual('/path/to/rinex.obs')

    >>> # Process multiple files with custom options
    >>> files = ['/path/to/file1.obs', '/path/to/file2.obs']
    >>> args = {'-options': '/custom/options.dat'}
    >>> outputs = prairie_manual(files, argsdict=args)

    argsdict :
        a dictionnary so as argsdict[argument] = val
        e.g.
        argsdict['-options'] = "/home/psakicki/THESE/SOFTWARES/GINSv2B/gin/data/prairie/options_GPS.dat"
        argsdict['-wsb']     = "/home/psakicki/THESE/SOFTWARES/GINSv2B/gin/data/prairie/WSB_NEW_1996.dat"

    argument cc2noncc, wsb, options, and out are automatically managed
    """

    if type(rinex_paths_in) is list:
        multimode = True
        rinex_path_lis = rinex_paths_in
    elif type(rinex_paths_in) is str:
        multimode = False
        rinex_path_lis = [rinex_paths_in]
    else:
        log.error("check the rinex_paths_in !!!")
        return None

    # be sure there is a TEMP DATA folder
    if temp_data_folder == "":
        temp_data_folder = os.path.join(ginscmn.get_gin_path(1), "TEMP_DATA")
    if not os.path.exists(temp_data_folder):
        os.makedirs(temp_data_folder)

    out_pra_path_lis = []
    for rinex_path in rinex_path_lis:
        #        # be sure the RINEX is in gins folder ...
        #        bool_rinex_in_gin = check_if_in_gin(rinex_path)
        #        # ... and copy it otherwise
        #        if not bool_rinex_in_gin:
        expected_rinex_path = os.path.join(
            temp_data_folder, os.path.basename(rinex_path)
        )
        if not os.path.isfile(expected_rinex_path):
            print("INFO : will be copied in ", temp_data_folder)
            rinex_path = ginscmn.copy_in_gin(rinex_path, temp_data_folder)
        else:
            print("INFO : ", expected_rinex_path, "already exists")
            rinex_path = expected_rinex_path

        # check if the RINEX is compressed ...
        bool_comp_rnx = operational.check_if_compressed_rinex(rinex_path)
        # ... if not crz2rnx !
        if bool_comp_rnx:
            crinex_path = rinex_path
            rinex_path = operational.crz2rnx(crinex_path, temp_data_folder)
            if not os.path.isfile(rinex_path):
                print("ERR : something went wrong during CRZ2RNX, skiping the file")
                print(crinex_path)

        rinex_name = os.path.basename(rinex_path)
        out_pra_name = rinex_name + ".pra"
        print("INFO : input & output name :", rinex_name, out_pra_name)
        out_pra_path = os.path.join(temp_data_folder, out_pra_name)

        if os.path.isfile(out_pra_path) and not force:
            print(
                "INFO : prairie_manual : ", out_pra_path, "already exists, skiping ..."
            )

        command = "exe_prairie "

        if force:
            forcearg = " -f "
        else:
            forcearg = ""

        command = command + " " + rinex_name + forcearg

        argsdict["-out"] = out_pra_name

        if not "-cc2noncc" in list(argsdict.keys()):
            argsdict["-cc2noncc"] = os.path.join(
                ginscmn.get_gin_path(1), "archives/p1c1bias.2000p"
            )
        if with_historik and not "-historik" in list(argsdict.keys()):
            argsdict["-historik"] = os.path.join(
                ginscmn.get_gin_path(1), "data/constell/historik_glonass"
            )
        if with_wsb and not "-wsb" in list(argsdict.keys()):
            argsdict["-wsb"] = os.path.join(ginscmn.get_gin_path(1), "archives/WSBREF.res.dat")
        if not "-options" in list(argsdict.keys()):
            argsdict["-options"] = os.path.join(
                ginscmn.get_gin_path(1), "data/prairie/options_GPS.dat"
            )

        for k, v in argsdict.items():
            if "-" == k[0]:
                moins = ""
            else:
                moins = "-"

            command = command + " " + moins + k + " " + v

        #        stream = os.popen(command)
        curdir = os.getcwd()
        os.chdir(temp_data_folder)
        p = subprocess.Popen(
            [command], shell=True
        )  # , stdout=subprocess.PIPE , stderr=subprocess.PIPE)
        p.wait()
        os.chdir(curdir)

        print("")
        print(command)
        print("")
        #
        if os.path.isfile(out_pra_path) and os.stat(out_pra_path).st_size > 0:
            print("INFO : prairie_manual : OK for manual prairie :)")
            print("output in ", out_pra_path)
            out_pra_path_lis.append(out_pra_path)
        else:
            print("ERR : prairie_manual : ", out_pra_path, "dont exist or empty :(")

    if len(out_pra_path_lis) == 1:
        return out_pra_path
    else:
        return out_pra_path_lis
