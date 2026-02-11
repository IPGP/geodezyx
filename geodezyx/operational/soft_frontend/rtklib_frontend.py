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
import concurrent.futures
import datetime as dt
#### Import the logger
import logging
import os
import subprocess
import numpy as np
import shutil as shutils
from threading import Lock

#### geodeZYX modules
from geodezyx import files_rw
from geodezyx import operational
from geodezyx import utils
from geodezyx import conv

log = logging.getLogger('geodezyx')

##########  END IMPORT  ##########

def read_conf_file(filein):
    outdic = collections.OrderedDict()
    with open(filein) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            key, *val = line.split("=")
            outdic[key.strip()] = val[0].split("#")[0].strip() if val else ""
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
    tmp_dir=None,
    prod_dir=None,
    experience_prefix="",
    rover_auto_conf=False,
    base_auto_conf=True,
    xyz_rover = [0, 0, 0],
    xyz_base=[0, 0, 0],
    outtype="auto",
    calc_center="IGS0OPSFIN",
    force=False,
    clean_tmp=False,
    exe_path = '/home/psakicki/SOFTWARE/RTKLIB_explorer/RTKLIB/app/consapp/rnx2rtkp/gcc/rnx2rtkp'
):
    r"""
    Run RTKLIB `rnx2rtkp` from rover/base RINEX observations using a generic RTKLIB
    configuration file, optionally overriding antenna/receiver metadata from RINEX
    headers and downloading required GNSS products (SP3 precise orbits and BRDC nav).

    The function:
      - creates output and temporary directory structure
      - uncompresses compressed RINEX inputs if needed
      - reads start/end times and sampling interval from the rover/base RINEX
      - builds a run-specific RTKLIB config file from `generik_conf`
      - optionally fills rover/base antenna positions and eccentricities from RINEX headers
      - downloads required products into `prod_dir`:
        - SP3 precise orbit from `calc_center`
        - BRDC navigation RINEX
      - calls the external executable `rnx2rtkp` with assembled arguments
      - writes outputs to the `out_dir` directory

    Parameters
    ----------
    rnx_rover : str | os.PathLike
        Path to rover observation RINEX. May be compressed (e.g. `.gz`, `.Z` and
        other formats supported by `geodezyx.operational.check_if_compressed_rinex`).
    rnx_base : str | os.PathLike
        Path to base observation RINEX. May be compressed.
    generik_conf : str | os.PathLike
        Path to a generic RTKLIB configuration file (key=value format). This file is
        parsed by `read_conf_file` and then overridden according to function options.
    out_dir : str | os.PathLike
        Directory where results are saved. This parameter is mandatory.
    tmp_dir : str | os.PathLike | None, default=None
        Temporary directory for intermediate files. Optional. If not provided,
        defaults to `out_dir/TMP`.
    prod_dir : str | os.PathLike | None, default=None
        Directory where to search for orbits, clocks, and BRDC files. Optional.
        If not provided, defaults to `tmp_dir`.
    experience_prefix : str, default=""
        Prefix added to the output file stem used to name the generated `.conf`
        and `.out` files.
    rover_auto_conf : bool, default=False
        If True, reads rover RINEX header and injects rover antenna type, XYZ position,
        and antenna eccentricities into the produced config (keys `ant1-\*`).
        If False, rover settings remain those defined in `generik_conf`.
    base_auto_conf : bool, default=True
        If True, reads base RINEX header and injects base antenna type, XYZ position,
        and antenna eccentricities into the produced config (keys `ant2-\*`).
    xyz_rover : list[float], default=[0, 0, 0]
        Rover station ECEF XYZ coordinates in meters \([X, Y, Z]\). If `xyz_rover[0] != 0`,
        these coordinates override the rover RINEX header coordinates; otherwise the rover
        header coordinates are used.
    xyz_base : list[float], default=[0, 0, 0]
        Base station ECEF XYZ coordinates in meters \([X, Y, Z]\). If `xyz_base[0] != 0`,
        these coordinates override the base RINEX header coordinates; otherwise the base
        header coordinates are used.
        Note: this is a simple sentinel check; a real X coordinate of 0 would be treated
        as "not provided".
    outtype : str, default="auto"
        Output solution format. If not `"auto"`, sets RTKLIB `out-solformat` to one of:
        `"dms"`, `"deg"`, `"xyz"`, `"enu"` (case-insensitive).
        If `"auto"`, uses whatever is defined in `generik_conf`.
    calc_center : str, default="IGS0OPSFIN"
        Analysis center/product identifier used when downloading GNSS precise products
        via `geodezyx.operational.download_gnss_products`. Must match what that
        downloader expects.
    force : bool, default=False
        If True, forces reprocessing even if output file already exists.
    clean_tmp : bool, default=False
        If True, deletes all contents of `tmp_dir` at the beginning of the
        function execution.
    exe_path : str, default=...
        Filesystem path to the RTKLIB `rnx2rtkp` executable.

    Returns
    -------
    out_result_fil : str
        Path to the RTKLIB output solution file (`.out`).

    Raises
    ------
    FileNotFoundError
        If no SP3 orbit file or no BRDC navigation file can be found locally or downloaded.

    Notes
    -----
    - The generated config file is written to `out_dir` and named using:
      `<experience_prefix>\_<rover4>\_<base4>\_<YYYY\_DOY>.conf`.
    - The RTKLIB command is executed via `subprocess.call(..., shell=True)` using
      `/bin/bash`, with a single concatenated command string.
    - A warning is emitted if the base RINEX time span does not fully cover the rover span.
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

    # RINEX START & END - FAST
    rov_srt_fast = conv.rinexname2dt(rnx_rover)
    bas_srt_fast = conv.rinexname2dt(rnx_base)

    # RINEX NAMES
    rov_name = os.path.basename(rnx_rover)[0:4]
    bas_name = os.path.basename(rnx_base)[0:4]

    # EXPERIENCE / OUTPUT NAMES
    srt_str = rov_srt_fast.strftime("%Y_%j_%H%M")
    exp_full_name = "_".join((experience_prefix, rov_name, bas_name, srt_str))

    out_conf_fil = os.path.join(out_dir, exp_full_name + ".conf")
    out_result_fil = os.path.join(out_dir, exp_full_name + ".out")

    if not force and os.path.isfile(out_result_fil) and os.path.getsize(out_result_fil) > 2000:
        log.info(f"RTKLIB output file {out_result_fil} already exists. Skipping...")
        return out_result_fil

    #### START HEAVY PROCESSING ####
    log.info("Starting RTKLIB RUN for {}".format(exp_full_name))
    # uncompressing rinex if compressed
    if operational.check_if_compressed_rinex(rnx_rover):
        rnx_rover = operational.crz2rnx(rnx_rover, tmp_dir)
    if operational.check_if_compressed_rinex(rnx_base):
        rnx_base = operational.crz2rnx(rnx_base, tmp_dir)

    # RINEX START & END - PRECISE
    rov_srt, rov_end, rov_itv = operational.rinex_start_end(rnx_rover, True)
    bas_srt, bas_end, bas_itv = operational.rinex_start_end(rnx_base, True)

    log.info("Rover: {} | Base: {} | Start: {} | End: {} | Interval: {} s".format(
        rov_name, bas_name, rov_srt, rov_end, rov_itv)
    )

    # READ GENERIC CONF FILE
    dicoconf = read_conf_file(generik_conf)

    if not outtype.lower() == "auto":
        dicoconf["out-solformat"] = outtype.lower()
        log.info(f"out-solformat was 'auto', set to: {dicoconf['out-solformat']}")

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
    if rover_auto_conf:
        _edit_dicoconf(rnx_rover, xyz_rover, ant_n=1)
    if base_auto_conf:
        _edit_dicoconf(rnx_base, xyz_base, ant_n=2)

    if not (bas_srt <= rov_srt <= rov_end <= bas_end):
        log.warning("rover/base epoch inconsistency: not bas_srt <= rov_srt <= rov_end <= bas_end !!!")

    # write conf file
    write_conf_file(dicoconf, out_conf_fil)

    ##### ORBITS
    # Function to decompress files into tmp_dir
    ## SP3
    # if "FIN" in calc_center or "RAP" in calc_center:
    #     orb_srt = conv.round_dt(bas_srt,"1D", mode="floor")
    #     orb_end = conv.round_dt(bas_end,"1D", mode="floor")
    # elif "ULT" in calc_center:
    #     orb_srt = conv.round_dt(bas_srt,"6H", mode="floor")
    #     orb_end = conv.round_dt(bas_end,"6H", mode="floor")
    # else:
    #     orb_srt = bas_srt
    #     orb_end = bas_end
    #
    # prodlis = operational.download_gnss_products(
    #     archive_dir=prod_dir,
    #     startdate=orb_srt,
    #     enddate=orb_end,
    #     archtype="year/doy",
    #     AC_names=(calc_center , ),
    #     archive_center="ign"
    # )

    prodlis = operational.dl_prods(prod_dir,
                         (bas_srt, bas_end),
                         calc_center,
                         prod_types=("sp3","clk"),
                         data_centers=("ign",))

    if not prodlis or not prodlis[0]:
        log.error("No SP3/CLK product files found remotely nor locally")
        log.error(f"Is analysis center/latency correct?: {calc_center}")
        log.warning(f"We continue without precise orbits/only with broadcast ones")
        prodlis_ok = []
    else:
        prodlis_ok = [prods2tmp(orb,tmp_dir) for orb in prodlis]


    ### BRDC
    #statdic = dict()
    #statdic["nav"] = ["BRDC"]
    #nav_srt = conv.round_dt(bas_srt, "1D", mode="floor")
    #nav_end = conv.round_dt(bas_end, "1D", mode="floor")
    brdclis = operational.dl_brdc(prod_dir,(bas_srt, bas_end))

    #brdclis = operational.download_gnss_rinex(
    #    statdic, prod_dir, nav_srt, nav_end, archtype="year/doy"
    #)


    #brdc_path_lis , brdc_bool_lis = zip(*brdclis)
    if len(brdclis) == 0:#  or sum(brdc_bool_lis) == 0:
        log.error("No BRDC nav file found remotely nor locally")
        raise FileNotFoundError("No BRDC nav file found remotely nor locally")
    else:
        #brdclis_ok = [prods2tmp(n,tmp_dir) for n,b in brdclis if b]
        brdclis_ok = [prods2tmp(n,tmp_dir) for n in brdclis]

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
            " ".join(brdclis_ok),
            " ".join(prodlis_ok),
        )
    )
    log.info("Rover: {} | Base: {} | Start: {} | End: {} | Interval: {} s".format(
        rov_name, bas_name, rov_srt, rov_end, rov_itv)
    )
    log.info(bigcomand)
    subprocess.call([bigcomand], executable="/bin/bash", shell=True)
    log.info("RTKLIB RUN FINISHED for {}".format(exp_full_name))

    return out_result_fil


def _prepare_shared_products(prod_dir, time_span_list, calc_center, download_lock=None):
    """
    Download shared products (SP3/CLK and BRDC) once for all parallel runs.

    Parameters
    ----------
    prod_dir : str
        Directory where products will be downloaded.
    time_span_list : list of tuples
        List of (start_time, end_time) tuples for all runs.
    calc_center : str
        Analysis center identifier.
    download_lock : threading.Lock, optional
        Lock to protect concurrent downloads (if needed).

    Returns
    -------
    tuple
        (prodlis_ok, brdclis_ok) - lists of product and broadcast file paths
    """
    # Find the overall time span covering all runs
    all_starts = [span[0] for span in time_span_list]
    all_ends = [span[1] for span in time_span_list]
    global_start = min(all_starts)
    global_end = max(all_ends)

    log.info(f"Downloading shared products for period {global_start} to {global_end}")

    # Download SP3/CLK products
    if download_lock:
        download_lock.acquire()
    try:
        prodlis = operational.dl_prods(prod_dir,
                             (global_start, global_end),
                             calc_center,
                             prod_types=("sp3","clk"),
                             data_centers=("ign",))
    finally:
        if download_lock:
            download_lock.release()

    if not prodlis or not prodlis[0]:
        log.warning("No SP3/CLK product files found remotely nor locally")
        prodlis_paths = []
    else:
        prodlis_paths = prodlis

    # Download BRDC files
    if download_lock:
        download_lock.acquire()
    try:
        brdclis = operational.dl_brdc(prod_dir, (global_start, global_end))
    finally:
        if download_lock:
            download_lock.release()

    if len(brdclis) == 0:
        log.error("No BRDC nav file found remotely nor locally")
        raise FileNotFoundError("No BRDC nav file found remotely nor locally")

    brdclis_paths = brdclis

    log.info(f"Shared products downloaded: {len(prodlis_paths)} SP3/CLK, {len(brdclis_paths)} BRDC")

    return prodlis_paths, brdclis_paths


def _rtklib_run_worker(
    rnx_rover, rnx_base, generik_conf, out_dir, tmp_dir,
    prodlis_shared, brdclis_shared,
    experience_prefix="", rover_auto_conf=False, base_auto_conf=True,
    xyz_rover=[0, 0, 0], xyz_base=[0, 0, 0], outtype="auto",
    force=False, exe_path='/home/psakicki/SOFTWARE/RTKLIB_explorer/RTKLIB/app/consapp/rnx2rtkp/gcc/rnx2rtkp'
):
    """
    Worker function for parallel RTKLIB processing.
    Uses pre-downloaded shared products to avoid download concurrency issues.

    This is the function that gets parallelized.
    """
    # RINEX START & END - FAST
    rov_srt_fast = conv.rinexname2dt(rnx_rover)
    bas_srt_fast = conv.rinexname2dt(rnx_base)

    # RINEX NAMES
    rov_name = os.path.basename(rnx_rover)[0:4]
    bas_name = os.path.basename(rnx_base)[0:4]

    # EXPERIENCE / OUTPUT NAMES
    srt_str = rov_srt_fast.strftime("%Y_%j_%H%M")
    exp_full_name = "_".join((experience_prefix, rov_name, bas_name, srt_str))

    out_conf_fil = os.path.join(out_dir, exp_full_name + ".conf")
    out_result_fil = os.path.join(out_dir, exp_full_name + ".out")

    if not force and os.path.isfile(out_result_fil) and os.path.getsize(out_result_fil) > 2000:
        log.info(f"RTKLIB output file {out_result_fil} already exists. Skipping...")
        return out_result_fil

    #### START HEAVY PROCESSING ####
    log.info("Starting RTKLIB RUN for {}".format(exp_full_name))

    # Create worker-specific tmp directory to avoid conflicts
    worker_tmp_dir = os.path.join(tmp_dir, exp_full_name)
    utils.create_dir(worker_tmp_dir)

    # uncompressing rinex if compressed
    if operational.check_if_compressed_rinex(rnx_rover):
        rnx_rover = operational.crz2rnx(rnx_rover, worker_tmp_dir)
    if operational.check_if_compressed_rinex(rnx_base):
        rnx_base = operational.crz2rnx(rnx_base, worker_tmp_dir)

    # RINEX START & END - PRECISE
    rov_srt, rov_end, rov_itv = operational.rinex_start_end(rnx_rover, True)
    bas_srt, bas_end, bas_itv = operational.rinex_start_end(rnx_base, True)

    log.info("Rover: {} | Base: {} | Start: {} | End: {} | Interval: {} s".format(
        rov_name, bas_name, rov_srt, rov_end, rov_itv)
    )

    # READ GENERIC CONF FILE
    dicoconf = read_conf_file(generik_conf)

    if not outtype.lower() == "auto":
        dicoconf["out-solformat"] = outtype.lower()

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
    if rover_auto_conf:
        _edit_dicoconf(rnx_rover, xyz_rover, ant_n=1)
    if base_auto_conf:
        _edit_dicoconf(rnx_base, xyz_base, ant_n=2)

    if not (bas_srt <= rov_srt <= rov_end <= bas_end):
        log.warning("rover/base epoch inconsistency: not bas_srt <= rov_srt <= rov_end <= bas_end !!!")

    # write conf file
    write_conf_file(dicoconf, out_conf_fil)

    # Copy shared products to worker tmp directory
    prodlis_ok = [prods2tmp(orb, worker_tmp_dir) for orb in prodlis_shared] if prodlis_shared else []
    brdclis_ok = [prods2tmp(n, worker_tmp_dir) for n in brdclis_shared]

    # Command
    arg_config = "-k " + out_conf_fil
    arg_interval = "-ti " + str(np.round(rov_itv, 6))
    arg_mode = ""
    arg_resultfile = "-o " + out_result_fil

    bigcomand = " ".join(
        (
            exe_path,
            arg_config,
            arg_interval,
            arg_mode,
            arg_resultfile,
            rnx_rover,
            rnx_base,
            " ".join(brdclis_ok),
            " ".join(prodlis_ok),
        )
    )

    log.info("Rover: {} | Base: {} | Start: {} | End: {} | Interval: {} s".format(
        rov_name, bas_name, rov_srt, rov_end, rov_itv)
    )
    log.info(bigcomand)

    # THIS IS THE PART THAT GETS PARALLELIZED
    subprocess.call([bigcomand], executable="/bin/bash", shell=True)

    log.info("RTKLIB RUN FINISHED for {}".format(exp_full_name))

    return out_result_fil


def rtklib_run_parallel(
    rinex_pairs,
    generik_conf,
    out_dir,
    tmp_dir=None,
    prod_dir=None,
    calc_center="IGS0OPSFIN",
    max_workers=4,
    experience_prefix="",
    rover_auto_conf=False,
    base_auto_conf=True,
    xyz_rovers=None,
    xyz_bases=None,
    outtype="auto",
    force=False,
    clean_tmp=False,
    exe_path='/home/psakicki/SOFTWARE/RTKLIB_explorer/RTKLIB/app/consapp/rnx2rtkp/gcc/rnx2rtkp'
):
    """
    Run RTKLIB processing in parallel for multiple rover/base RINEX pairs.

    This function:
      1. Downloads all required products (SP3/CLK and BRDC) once in a preliminary step
      2. Runs RTKLIB processing in parallel using ThreadPoolExecutor

    This avoids download concurrency issues by separating the download phase from
    the processing phase.

    Parameters
    ----------
    rinex_pairs : list of tuples
        List of (rnx_rover, rnx_base) file path tuples to process.
    generik_conf : str | os.PathLike
        Path to generic RTKLIB configuration file.
    out_dir : str | os.PathLike
        Directory where results are saved.
    tmp_dir : str | os.PathLike | None, default=None
        Temporary directory for intermediate files.
    prod_dir : str | os.PathLike | None, default=None
        Directory where to search for orbits, clocks, and BRDC files.
    calc_center : str, default="IGS0OPSFIN"
        Analysis center/product identifier.
    max_workers : int, default=4
        Maximum number of parallel workers.
    experience_prefix : str, default=""
        Prefix added to output file stems.
    rover_auto_conf : bool, default=False
        If True, reads rover RINEX header for antenna configuration.
    base_auto_conf : bool, default=True
        If True, reads base RINEX header for antenna configuration.
    xyz_rovers : list of lists, optional
        List of [X, Y, Z] coordinates for each rover station.
        If None, uses [0, 0, 0] for all.
    xyz_bases : list of lists, optional
        List of [X, Y, Z] coordinates for each base station.
        If None, uses [0, 0, 0] for all.
    outtype : str, default="auto"
        Output solution format.
    force : bool, default=False
        If True, forces reprocessing even if output exists.
    clean_tmp : bool, default=False
        If True, deletes tmp_dir contents at start.
    exe_path : str
        Filesystem path to the RTKLIB rnx2rtkp executable.

    Returns
    -------
    list
        List of output file paths for all processed pairs.

    Examples
    --------
    >>> pairs = [
    ...     ('/path/to/rover1.rnx', '/path/to/base1.rnx'),
    ...     ('/path/to/rover2.rnx', '/path/to/base2.rnx'),
    ... ]
    >>> results = rtklib_run_parallel(
    ...     pairs,
    ...     '/path/to/config.conf',
    ...     '/path/to/output',
    ...     max_workers=4
    ... )
    """
    # Directory structure setup
    out_dir = utils.create_dir(out_dir)

    if tmp_dir is None:
        tmp_dir = os.path.join(out_dir, 'TMP')
    tmp_dir = utils.create_dir(tmp_dir)

    if clean_tmp:
        shutils.rmtree(tmp_dir + "/*", ignore_errors=True)

    if prod_dir is None:
        prod_dir = tmp_dir
    else:
        prod_dir = utils.create_dir(prod_dir)

    # Prepare coordinate lists
    if xyz_rovers is None:
        xyz_rovers = [[0, 0, 0]] * len(rinex_pairs)
    if xyz_bases is None:
        xyz_bases = [[0, 0, 0]] * len(rinex_pairs)

    log.info(f"Starting parallel RTKLIB processing for {len(rinex_pairs)} pairs with {max_workers} workers")

    # STEP 1: Preliminary - determine time spans and download products once
    log.info("STEP 1: Determining time spans for all RINEX pairs...")
    time_spans = []
    for rnx_rover, rnx_base in rinex_pairs:
        # Use fast method to get approximate time range
        rov_srt = conv.rinexname2dt(rnx_rover)
        bas_srt = conv.rinexname2dt(rnx_base)
        # Assume 24h sessions - adjust if needed
        rov_end = rov_srt + dt.timedelta(days=1)
        bas_end = bas_srt + dt.timedelta(days=1)
        time_spans.append((min(rov_srt, bas_srt), max(rov_end, bas_end)))

    log.info("STEP 2: Downloading shared products (SP3/CLK and BRDC)...")
    try:
        prodlis_shared, brdclis_shared = _prepare_shared_products(
            prod_dir, time_spans, calc_center
        )
    except Exception as e:
        log.error(f"Failed to download shared products: {e}")
        raise

    # STEP 3: Parallel processing
    log.info(f"STEP 3: Running {len(rinex_pairs)} RTKLIB processes in parallel...")
    results = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_pair = {}
        for i, (rnx_rover, rnx_base) in enumerate(rinex_pairs):
            future = executor.submit(
                _rtklib_run_worker,
                rnx_rover=rnx_rover,
                rnx_base=rnx_base,
                generik_conf=generik_conf,
                out_dir=out_dir,
                tmp_dir=tmp_dir,
                prodlis_shared=prodlis_shared,
                brdclis_shared=brdclis_shared,
                experience_prefix=experience_prefix,
                rover_auto_conf=rover_auto_conf,
                base_auto_conf=base_auto_conf,
                xyz_rover=xyz_rovers[i],
                xyz_base=xyz_bases[i],
                outtype=outtype,
                force=force,
                exe_path=exe_path
            )
            future_to_pair[future] = (rnx_rover, rnx_base)

        # Collect results as they complete
        for future in concurrent.futures.as_completed(future_to_pair):
            pair = future_to_pair[future]
            try:
                result = future.result()
                results.append(result)
                log.info(f"Completed: {pair[0]} + {pair[1]} -> {result}")
            except Exception as exc:
                log.error(f"Failed processing {pair}: {exc}")
                results.append(None)

    log.info(f"Parallel RTKLIB processing completed. {sum(1 for r in results if r)} / {len(rinex_pairs)} successful.")

    return results

