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
from geodezyx import conv

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
    xyz_base=[0, 0, 0],
    outtype="auto",
    calc_center="IGS0OPSFIN",
    force=False,
    exe_path = '/home/psakicki/SOFTWARE/RTKLIB_explorer/RTKLIB/app/consapp/rnx2rtkp/gcc/rnx2rtkp'
):
    r"""
    Run RTKLIB \`rnx2rtkp\` from rover/base RINEX observations using a generic RTKLIB
    configuration file, optionally overriding antenna/receiver metadata from RINEX
    headers and downloading required GNSS products (SP3 precise orbits and BRDC nav).

    The function:
    \- creates a working directory structure (notably \`OUTPUT\` and \`TEMP\`)
    \- uncompresses compressed RINEX inputs if needed
    \- reads start/end times and sampling interval from the rover/base RINEX
    \- builds a run\-specific RTKLIB config file from \`generik_conf\`
    \- optionally fills rover/base antenna positions and eccentricities from RINEX headers
    \- downloads required products into \`TEMP\`:
      \- SP3 precise orbit from \`calc_center\`
      \- BRDC navigation RINEX
    \- calls the external executable \`rnx2rtkp\` with assembled arguments
    \- writes outputs to the \`OUTPUT\` directory

    Parameters
    ----------
    rnx_rover : str | os.PathLike
        Path to rover observation RINEX. May be compressed (e.g. \`.gz\`, \`.Z\` and
        other formats supported by \`geodezyx.operational.check_if_compressed_rinex\`).
    rnx_base : str | os.PathLike
        Path to base observation RINEX. May be compressed.
    generik_conf : str | os.PathLike
        Path to a generic RTKLIB configuration file (key=value format). This file is
        parsed by \`read_conf_file\` and then overridden according to function options.
    working_dir : str | os.PathLike
        Directory used for temporary files and outputs. Subdirectories created:
        \`TEMP\` for downloads/uncompressed files and \`OUTPUT\` for results.
    experience_prefix : str, default=""
        Prefix added to the output file stem used to name the generated \`.conf\`
        and \`.out\` files.
    rover_auto_conf : bool, default=False
        If True, reads rover RINEX header and injects rover antenna type, XYZ position,
        and antenna eccentricities into the produced config (keys \`ant1\-\*\`).
        If False, rover settings remain those defined in \`generik_conf\`.
    base_auto_conf : bool, default=True
        If True, reads base RINEX header and injects base antenna type, XYZ position,
        and antenna eccentricities into the produced config (keys \`ant2\-\*\`).
    xyz_base : list[float], default=[0, 0, 0]
        Base station ECEF XYZ coordinates in meters \([X, Y, Z]\). If \`xyz_base[0] != 0\`,
        these coordinates override the base RINEX header coordinates; otherwise the base
        header coordinates are used.
        Note: this is a simple sentinel check; a real X coordinate of 0 would be treated
        as "not provided".
    outtype : str, default="auto"
        Output solution format. If not \`"auto"\`, sets RTKLIB \`out\-solformat\` to one of:
        \`"dms"\`, \`"deg"\`, \`"xyz"\`, \`"enu"\` (case\-insensitive).
        If \`"auto"\`, uses whatever is defined in \`generik_conf\`.
    calc_center : str, default="IGS0OPSFIN"
        Analysis center/product identifier used when downloading GNSS precise products
        via \`geodezyx.operational.download_gnss_products\`. Must match what that
        downloader expects.
    exe_path : str, default=...
        Filesystem path to the RTKLIB \`rnx2rtkp\` executable.

    Returns
    -------
    out_result_fil : str
        Path to the RTKLIB output solution file (\`.out\`).

    Raises
    ------
    FileNotFoundError
        If no SP3 orbit file or no BRDC navigation file can be found locally or downloaded.

    Notes
    -----
    \- The generated config file is written to \`working_dir/OUTPUT\` and named using:
      \`<experience_prefix>\_<rover4>\_<base4>\_<YYYY\_DOY>.conf\`.
    \- The RTKLIB command is executed via \`subprocess.call(..., shell=True)\` using
      \`/bin/bash\`, with a single concatenated command string.
    \- A warning is emitted if the base RINEX time span does not fully cover the rover span.
    """
    import shutil
    # \(... implementation continues ...\)
    import shutil

    # paths & files
    working_dir = utils.create_dir(working_dir)
    out_dir = utils.create_dir(os.path.join(working_dir, "OUTPUT"))

    # paths & files
    temp_dir = os.path.join(working_dir,'TEMP')
    utils.create_dir(temp_dir)
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

    srt_str = rov_srt.strftime("%Y_%j_%H%M")
    exp_full_name = "_".join((experience_prefix, rov_name, bas_name, srt_str))

    out_conf_fil = os.path.join(out_dir, exp_full_name + ".conf")
    out_result_fil = os.path.join(out_dir, exp_full_name + ".out")

    if os.path.isfile(out_result_fil) and not force:
        log.info(f"RTKLIB output file {out_result_fil} already exists. Skipping...")
        return out_result_fil

    dicoconf = read_conf_file(generik_conf)

    if not outtype.lower() == "auto":
        dicoconf["out-solformat"] = outtype.lower()
        log.info(f"out-solformat {dicoconf['out-solformat']}")

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

    if base_auto_conf:
        antobj_bas, recobj_bas, siteobj_bas, locobj_bas = (
            files_rw.read_rinex_2_dataobjts(rnx_base)
        )
        dicoconf["ant2-postype"] = "xyz"
        dicoconf["ant2-anttype"] = antobj_bas.Antenna_Type
        if xyz_base[0] != 0:
            dicoconf["ant2-pos1"] = xyz_base[0]
            dicoconf["ant2-pos2"] = xyz_base[1]
            dicoconf["ant2-pos3"] = xyz_base[2]
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

    ##### ORBITS
    ## SP3
    if "FIN" in calc_center:
        orb_srt = conv.round_dt(bas_srt,"1D", mode="floor")
        orb_end = conv.round_dt(bas_end,"1D", mode="ceil")
    elif "RAP" in calc_center:
        orb_srt = conv.round_dt(bas_srt,"1D", mode="floor")
        orb_end = conv.round_dt(bas_end,"1D", mode="floor")
    elif "ULT" in calc_center:
        orb_srt = conv.round_dt(bas_srt,"6H", mode="floor")
        orb_end = conv.round_dt(bas_end,"6H", mode="floor")
    else:
        orb_srt = bas_srt
        orb_end = bas_end

    orblis = operational.download_gnss_products(
        archive_dir=temp_dir,
        startdate=orb_srt,
        enddate=orb_end,
        archtype="/",
        AC_names=(calc_center , )
    )

    unzip = lambda o: files_rw.unzip_gz_z(o) if o.endswith(".gz") or o.endswith(".Z") else o

    if not orblis or not orblis[0]:
        log.error("No SP3 orbit file found remotely nor locally")
        log.error(f"Is analysis center/latency correct?: {calc_center}")
        raise FileNotFoundError("No SP3 orbit file found remotely nor locally")
    else:
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
        log.error("No BRDC nav file found remotely nor locally")
        raise FileNotFoundError("No BRDC nav file found remotely nor locally")

    # Command
    com_config = "-k " + out_conf_fil
    com_interval = "-ti " + str(np.round(rov_itv,6))
    com_mode = ""
    # com_mode="-p 4"
    com_resultfile = "-o " + out_result_fil
    # com_combinsol="-c"


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

    return out_result_fil
