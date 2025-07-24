#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.operational contains functions to run the 
GNSS processing software RTKLIB. 

it can be imported directly with:
from geodezyx import operational

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""


import collections
########## BEGIN IMPORT ##########
#### External modules
import datetime as dt
#### Import the logger
import logging
import os
import subprocess
import numpy as np

#### geodeZYX modules
from geodezyx import files_rw
from geodezyx import operational
from geodezyx import utils

log = logging.getLogger('geodezyx')

##########  END IMPORT  ##########


def read_conf_file(filein):
    outdic = collections.OrderedDict()
    for l in open(filein):
        if l[0] == "#":
            continue
        if l.strip() == "":
            continue

        f = l.split("=")
        key = f[0]
        if len(f) == 1:
            val = ""
        else:
            try:
                val = f[1].split("#")[0]
            except IndexError:
                val = f[1]
        outdic[key.strip()] = val.strip()
    return outdic


def rtklib_run_from_rinex(
    rnx_rover,
    rnx_base,
    generik_conf,
    working_dir,
    experience_prefix="",
    rover_auto_conf=False,
    base_auto_conf=True,
    XYZbase=[0, 0, 0],
    outtype="auto",
    calc_center="IGS0OPSFIN",
    exe_path = '/home/psakicki/SOFTWARE/RTKLIB_explorer/RTKLIB/app/consapp/rnx2rtkp/gcc/rnx2rtkp'
):
    """
    auto_conf :
        read the header of the rinex and write the conf. file according to it
        if the mode is disabled, the antenna/rec part of the
        conf. file will be the same as the generic one

        NB : RTKLIB "core" have it's own reading header option.
        Thus, my advice is to disable the auto mode for the rover and leave
        ``ant1-postype=rinexhead`` in the generic conf file
        and enable it for the base with a XYZbase vector initialized or
        the good XYZ in header on the rinex (the XYZbase vector is prioritary
        over the base RINEX header)

        (but in anycase, prepars good rinex's headers ;) )

    outtype :
        'auto' (as defined in the generic config file) or
        'dms' 'deg' 'xyz' 'enu'
        can manage the upper case XYZ or FLH
    """

    import shutil

    # paths & files
    working_dir = utils.create_dir(working_dir)
    out_dir = utils.create_dir(os.path.join(working_dir, "OUTPUT"))

    # paths & files
    temp_dir = os.path.join(working_dir,'TEMP')
    clean_temp_dir = False
    if clean_temp_dir:
        shutil.rmtree(temp_dir)
        temp_dir = os.path.join(working_dir,'TEMP')
    out_dir  = os.path.join(working_dir,'OUTPUT')

    # uncompressing rinex if compressed
    if operational.check_if_compressed_rinex(rnx_rover):
        rnx_rover = operational.crz2rnx(rnx_rover, temp_dir)
    if operational.check_if_compressed_rinex(rnx_base):
        rnx_base = operational.crz2rnx(rnx_base, temp_dir)

    # RINEX START & END
    rov_srt, rov_end, rov_itv = operational.rinex_start_end(rnx_rover, 1)
    bas_srt, bas_end, bas_itv = operational.rinex_start_end(rnx_base, 1)

    # RINEX NAMES
    rov_name = os.path.basename(rnx_rover)[0:4]
    bas_name = os.path.basename(rnx_base)[0:4]

    srt_str = rov_srt.strftime("%Y_%j")
    exp_full_name = "_".join((experience_prefix, rov_name, bas_name, srt_str))

    out_conf_fil = os.path.join(out_dir, exp_full_name + ".conf")
    out_result_fil = os.path.join(out_dir, exp_full_name + ".out")

    dicoconf = read_conf_file(generik_conf)

    if rover_auto_conf:
        antobj_rov, recobj_rov, siteobj_rov, locobj_rov = (
            files_rw.read_rinex_2_dataobjts(rnx_rover)
        )
        dicoconf["ant1-postype"] = "xyz"
        dicoconf["ant1-anttype"] = antobj_rov.Antenna_Type
        dicoconf["ant1-pos1"] = locobj_rov.X_coordinate_m
        dicoconf["ant1-pos2"] = locobj_rov.Y_coordinate_m
        dicoconf["ant1-pos3"] = locobj_rov.Z_coordinate_m
        dicoconf["ant1-antdelu"] = antobj_rov.Up_Ecc
        dicoconf["ant1-antdeln"] = antobj_rov.North_Ecc
        dicoconf["ant1-antdele"] = antobj_rov.East_Ecc

    if not outtype.lower() == "auto":
        dicoconf["out-solformat"] = outtype.lower()
        log.info("out-solformat", dicoconf["out-solformat"])

    if base_auto_conf:
        antobj_bas, recobj_bas, siteobj_bas, locobj_bas = (
            files_rw.read_rinex_2_dataobjts(rnx_base)
        )
        dicoconf["ant2-postype"] = "xyz"
        dicoconf["ant2-anttype"] = antobj_bas.Antenna_Type
        if XYZbase[0] != 0:
            dicoconf["ant2-pos1"] = XYZbase[0]
            dicoconf["ant2-pos2"] = XYZbase[1]
            dicoconf["ant2-pos3"] = XYZbase[2]
        else:
            dicoconf["ant2-pos1"] = locobj_bas.X_coordinate_m
            dicoconf["ant2-pos2"] = locobj_bas.Y_coordinate_m
            dicoconf["ant2-pos3"] = locobj_bas.Z_coordinate_m
        dicoconf["ant2-antdelu"] = antobj_bas.Up_Ecc
        dicoconf["ant2-antdeln"] = antobj_bas.North_Ecc
        dicoconf["ant2-antdele"] = antobj_bas.East_Ecc

    if not (bas_srt <= rov_srt <= rov_end <= bas_end):
        log.warning("not bas_srt <= rov_srt <= rov_end <= bas_end !!!")

    outconffilobj = open(out_conf_fil, "w+")
    for k, v in dicoconf.items():
        lin = k.ljust(20) + "=" + str(v) + "\n"
        outconffilobj.write(lin)
    outconffilobj.close()

    ###### ORBITS
    ### SP3
    orblis = operational.download_gnss_products(
        archive_dir=temp_dir,
        startdate=bas_srt,
        enddate=bas_end,
        archtype="/",
        AC_names=(calc_center , )
    )

    unzip = lambda o: files_rw.unzip_gz_z(o) if o.endswith(".gz") or o.endswith(".Z") else o

    sp3_z = orblis[0]
    sp3 = unzip(sp3_z)

    ### BRDC
    statdic = dict()
    statdic["nav"] = ["BRDC"]
    nav_srt = dt.datetime(bas_srt.year, bas_srt.month, bas_srt.day)
    brdclis = operational.download_gnss_rinex(
        statdic, temp_dir, nav_srt, bas_end, archtype="/"
    )

    if len(brdclis) > 0 and brdclis[0][1]:
        nav_z = brdclis[0][0]
        nav = unzip(nav_z)

    else:
        log.error("No BRDC nav file found in the archive")
        raise FileNotFoundError("No BRDC nav file found in the archive")

    # Command
    com_config = "-k " + out_conf_fil
    com_interval = "-ti " + str(np.round(rov_itv,6))
    com_mode = ""
    # com_mode="-p 4"
    com_resultfile = "-o " + out_result_fil
    # com_combinsol="-c"

    # exe_path = "rnx2rtkp"
    #    exe_path = "/home/pierre/install_softs/RTKLIB/rnx2rtkp"
    # exe_path = "/home/psakicki/SOFTWARE/RTKLIB/RTKLIB/app/rnx2rtkp/gcc/rnx2rtkp"
    # exe_path = ""

    bigcomand = " ".join(
        (
            exe_path,
            com_config,
            com_interval,
            com_mode,
            com_resultfile,
            rnx_rover,
            rnx_base,
            nav,
            sp3,
        )
    )

    log.info(bigcomand)

    subprocess.call([bigcomand], executable="/bin/bash", shell=True)
    log.info("RTKLIB RUN FINISHED")

    return None
