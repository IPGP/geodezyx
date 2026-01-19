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
import shutil as shutils

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

def write_conf_file(dicoconf, fpath_out):
    f = open(fpath_out, "w+")
    for k, v in dicoconf.items():
        lin = k.ljust(20) + "=" + str(v) + "\n"
        f.write(lin)
    f.close()


def prods2tmp(fpath_inp, dir_out):
    """
    Decompresses a file into the temporary directory if compressed,
    or copies it to the temporary directory if not already there.
    """
    if fpath_inp.endswith((".gz", ".Z")):
        # Decompress into dir_out
        fpath_out = files_rw.unzip_gz_z(fpath_inp, out_gzip_dir=dir_out)
    else:
        # Copy to dir_out if not already there
        fpath_out = os.path.join(dir_out, os.path.basename(fpath_inp))
        if not os.path.exists(fpath_out):
            shutils.copy2(fpath_inp, fpath_out)
    return fpath_out

def rtklib_run_from_rinex(
    rnx_rover,
    rnx_base,
    generik_conf,
    out_dir,
    experience_prefix="",
    rover_auto_conf=False,
    base_auto_conf=True,
    xyz_rover = [0, 0, 0],
    xyz_base=[0, 0, 0],
    outtype="auto",
    calc_center="IGS0OPSFIN",
    force=False,
    tmp_dir=None,
    prod_dir=None,
    clean_tmp=False,
    exe_path = '/home/psakicki/SOFTWARE/RTKLIB_explorer/RTKLIB/app/consapp/rnx2rtkp/gcc/rnx2rtkp'
):
    r"""
    Run RTKLIB \`rnx2rtkp\` from rover/base RINEX observations using a generic RTKLIB
    configuration file, optionally overriding antenna/receiver metadata from RINEX
    headers and downloading required GNSS products (SP3 precise orbits and BRDC nav).

    The function:
    \- creates output and temporary directory structure
    \- uncompresses compressed RINEX inputs if needed
    \- reads start/end times and sampling interval from the rover/base RINEX
    \- builds a run\-specific RTKLIB config file from \`generik_conf\`
    \- optionally fills rover/base antenna positions and eccentricities from RINEX headers
    \- downloads required products into \`prod_dir\`:
      \- SP3 precise orbit from \`calc_center\`
      \- BRDC navigation RINEX
    \- calls the external executable \`rnx2rtkp\` with assembled arguments
    \- writes outputs to the \`out_dir\` directory

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
    out_dir : str | os.PathLike
        Directory where results are saved. This parameter is mandatory.
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
    xyz_rover : list[float], default=[0, 0, 0]
        Rover station ECEF XYZ coordinates in meters \([X, Y, Z]\). If \`xyz_rover[0] != 0\`,
        these coordinates override the rover RINEX header coordinates; otherwise the rover
        header coordinates are used.
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
    force : bool, default=False
        If True, forces reprocessing even if output file already exists.
    tmp_dir : str | os.PathLike | None, default=None
        Temporary directory for intermediate files. Optional. If not provided,
        defaults to \`out_dir/TMP\`.
    prod_dir : str | os.PathLike | None, default=None
        Directory where to search for orbits, clocks, and BRDC files. Optional.
        If not provided, defaults to \`tmp_dir\`.
    clean_tmp : bool, default=False
        If True, deletes all contents of \`tmp_dir\` at the beginning of the
        function execution.
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
    \- The generated config file is written to \`out_dir\` and named using:
      \`<experience_prefix>\_<rover4>\_<base4>\_<YYYY\_DOY>.conf\`.
    \- The RTKLIB command is executed via \`subprocess.call(..., shell=True)\` using
      \`/bin/bash\`, with a single concatenated command string.
    \- A warning is emitted if the base RINEX time span does not fully cover the rover span.
    """
    # Directory structure setup
    # out_dir: mandatory, where results are saved
    out_dir = utils.create_dir(out_dir)

    # tmp_dir: optional, defaults to out_dir/TMP if not provided
    if tmp_dir is None:
        tmp_dir = os.path.join(out_dir, 'TMP')
    tmp_dir = utils.create_dir(tmp_dir)

    if clean_tmp:
        shutils.rmtree(tmp_dir + "/*")

    # prod_dir: optional, defaults to tmp_dir if not provided
    # This is where we search for orbits/clocks/BRDC files
    if prod_dir is None:
        prod_dir = tmp_dir
    else:
        prod_dir = utils.create_dir(prod_dir)

    # uncompressing rinex if compressed
    if operational.check_if_compressed_rinex(rnx_rover):
        rnx_rover = operational.crz2rnx(rnx_rover, tmp_dir)
    if operational.check_if_compressed_rinex(rnx_base):
        rnx_base = operational.crz2rnx(rnx_base, tmp_dir)

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

    def _edit_dicoconf(rnx_inp, xyz_inp, ant_n=1):
        antobj, recobj, siteobj, locobj = (
            files_rw.read_rinex_2_dataobjts(rnx_inp)
        )
        n = str(ant_n)
        dicoconf[f"ant{n}-postype"] = "xyz"
        dicoconf[f"ant{n}-anttype"] = antobj.Antenna_Type
        if xyz_inp[0] != 0:
            dicoconf[f"ant{n}-pos1"] = xyz_inp[0]
            dicoconf[f"ant{n}-pos2"] = xyz_inp[1]
            dicoconf[f"ant{n}-pos3"] = xyz_inp[2]
        else:
            dicoconf[f"ant{n}-pos1"] = locobj.X_coordinate_m
            dicoconf[f"ant{n}-pos2"] = locobj.Y_coordinate_m
            dicoconf[f"ant{n}-pos3"] = locobj.Z_coordinate_m
        dicoconf[f"ant{n}-antdelu"] = antobj.Up_Ecc
        dicoconf[f"ant{n}-antdeln"] = antobj.North_Ecc
        dicoconf[f"ant{n}-antdele"] = antobj.East_Ecc

    # Edit conf file dic
    _edit_dicoconf(rnx_rover, xyz_rover, ant_n=1)
    _edit_dicoconf(rnx_base, xyz_base, ant_n=2)

    if not (bas_srt <= rov_srt <= rov_end <= bas_end):
        log.warning("rover/base epoch inconsistency: not bas_srt <= rov_srt <= rov_end <= bas_end !!!")

    # write conf file
    write_conf_file(dicoconf, out_conf_fil)

    ##### ORBITS
    # Function to decompress files into tmp_dir
    ## SP3
    if "FIN" in calc_center or "RAP" in calc_center:
        orb_srt = conv.round_dt(bas_srt,"1D", mode="floor")
        orb_end = conv.round_dt(bas_end,"1D", mode="floor")
    elif "ULT" in calc_center:
        orb_srt = conv.round_dt(bas_srt,"6H", mode="floor")
        orb_end = conv.round_dt(bas_end,"6H", mode="floor")
    else:
        orb_srt = bas_srt
        orb_end = bas_end

    orblis = operational.download_gnss_products(
        archive_dir=prod_dir,
        startdate=orb_srt,
        enddate=orb_end,
        archtype="/",
        AC_names=(calc_center , ),
        archive_center="ign"
    )

    if not orblis or not orblis[0]:
        log.error("No SP3 orbit file found remotely nor locally")
        log.error(f"Is analysis center/latency correct?: {calc_center}")
        raise FileNotFoundError("No SP3 orbit file found remotely nor locally")
    else:
        sp3lis = [prods2tmp(orb,tmp_dir) for orb in orblis]


    ### BRDC
    statdic = dict()
    statdic["nav"] = ["BRDC"]
    nav_srt = conv.round_dt(bas_srt, "1D", mode="floor")
    nav_end = conv.round_dt(bas_end, "1D", mode="floor")
    brdclis = operational.download_gnss_rinex(
        statdic, prod_dir, nav_srt, nav_end, archtype="year/doy"
    )

    brdc_path_lis , brdc_bool_lis = zip(*brdclis)
    if len(brdc_path_lis) == 0 or sum(brdc_bool_lis) == 0:
        log.error("No BRDC nav file found remotely nor locally")
        raise FileNotFoundError("No BRDC nav file found remotely nor locally")
    else:
        navlis = [prods2tmp(n,tmp_dir) for n,b in brdclis if b]

    # Command
    arg_config = "-k " + out_conf_fil
    arg_interval = "-ti " + str(np.round(rov_itv,6))
    arg_mode = ""
    # arg_mode="-p 4"
    arg_resultfile = "-o " + out_result_fil
    # com_combinsol="-c"

    bigcomand = " ".join(
        (
            exe_path,
            arg_config,
            arg_interval,
            arg_mode,
            arg_resultfile,
            rnx_rover,
            rnx_base,
            " ".join(navlis),
            " ".join(sp3lis),
        )
    )

    log.info(bigcomand)

    subprocess.call([bigcomand], executable="/bin/bash", shell=True)
    log.info("RTKLIB RUN FINISHED")

    return out_result_fil
