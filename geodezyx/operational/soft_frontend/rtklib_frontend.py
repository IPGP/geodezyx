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

########## BEGIN IMPORT ##########
#### External modules
import concurrent.futures

# import datetime as dt
import collections

#### Import the logger
import logging
import os
import subprocess
import numpy as np
import shutil as shutils

# from threading import Lock

#### geodeZYX modules
from geodezyx import files_rw
from geodezyx import operational
from geodezyx import utils
from geodezyx import conv
from geodezyx.operational.soft_frontend.common_frontend import (
    dl_prods,
    get_best_prods,
)

log = logging.getLogger("geodezyx")


##########  END IMPORT  ##########


def _read_cfg(filein):
    """
    Reads a RTKLIB configuration file into an ordered dictionary.
    """
    outdic = collections.OrderedDict()
    with open(filein) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            key, *val = line.split("=")
            outdic[key.strip()] = val[0].split("#")[0].strip() if val else ""
    return outdic


def _write_cfg(dicoconf, fpath_out):
    """
    Writes a RTKLIB configuration file from an ordered dictionary.
    """
    f = open(fpath_out, "w+")
    for k, v in dicoconf.items():
        lin = k.ljust(20) + "=" + str(v) + "\n"
        f.write(lin)
    f.close()


def _prods2tmp(fpath_inp, dir_out):
    """
    Decompresses a file into the temporary directory if compressed,
    or copies it to the temporary directory if not already there.
    """
    if fpath_inp.endswith(".gz"):
        # Decompress into dir_out
        fpath_out = files_rw.unzip_gz(fpath_inp, out_gzip_dir=dir_out, verbose=False)
    elif fpath_inp.endswith(".Z"):
        fpath_out = files_rw.unzip_gz_z(fpath_inp, out_gzip_dir=dir_out)
    else:
        # Copy to dir_out if not already there
        fpath_out = os.path.join(dir_out, os.path.basename(fpath_inp))
        if not os.path.exists(fpath_out):
            shutils.copy2(fpath_inp, fpath_out)
    return fpath_out


def rtklib_run_mono(
    rnx_rover,
    rnx_base,
    cfgfile_generik,
    out_dir,
    tmp_dir,
    prod_dir=None,
    igs_prods="GRG0OPSFIN",
    download_prods=True,
    orbclklis_inp=[],
    brdclis_inp=[],
    exp_prefix="",
    rover_auto_conf=True,
    base_auto_conf=True,
    xyz_rover=[0, 0, 0],
    xyz_base=[0, 0, 0],
    posmode=None,
    solformat=None,
    sateph=None,
    force=False,
    keep_tmp=False,
    exe_path="rnx2rtkp",
):
    """
    Worker function for parallel RTKLIB processing.
    Uses pre-downloaded shared products to avoid download concurrency issues.

    This is the function that gets parallelized.
    """

    if not download_prods and (not orbclklis_inp or not brdclis_inp):
        log.error("download_prods is False but orbclklis_inp or brdclis_inp is empty!")
        raise ValueError(
            "Must provide orbclklis_inp and brdclis_inp if not downloading products."
        )

    # RINEX START & END - FAST
    rov_srt_fast = conv.rinexname2dt(rnx_rover)
    # bas_srt_fast = conv.rinexname2dt(rnx_base)

    # RINEX NAMES
    rov_name = os.path.basename(rnx_rover)[0:4]
    bas_name = os.path.basename(rnx_base)[0:4]

    # EXPERIENCE / OUTPUT NAMES
    srt_str = rov_srt_fast.strftime("%Y_%j_%H%M")
    exp_full_name = "_".join((exp_prefix, rov_name, bas_name, srt_str))

    year_doy_str = rov_srt_fast.strftime("%Y/%j")
    out_dir_year_doy = str(os.path.join(out_dir, year_doy_str))
    utils.create_dir(out_dir_year_doy)

    out_conf_fil = os.path.join(out_dir_year_doy, exp_full_name + ".conf")
    out_res_fil = os.path.join(out_dir_year_doy, exp_full_name + ".out")

    if (
        not force
        and os.path.isfile(out_res_fil)
        and os.path.getsize(out_res_fil) > 2000
    ):
        log.info(f"Output file already exists, skipping: {out_res_fil}")
        return out_res_fil

    #### START HEAVY PROCESSING ####
    log.info("Starting RTKLIB RUN for {}".format(exp_full_name))

    # Create worker-specific tmp directory to avoid conflicts
    tmp_dir_wrk = os.path.join(tmp_dir, exp_full_name)
    utils.create_dir(tmp_dir_wrk)

    # uncompressing rinex if compressed
    if operational.check_if_compressed_rinex(rnx_rover):
        rnx_rover = operational.crz2rnx(rnx_rover, tmp_dir_wrk)
    if operational.check_if_compressed_rinex(rnx_base):
        rnx_base = operational.crz2rnx(rnx_base, tmp_dir_wrk)

    # RINEX START & END - PRECISE
    rov_srt, rov_end, rov_itv = operational.rinex_start_end(rnx_rover, True)
    bas_srt, bas_end, bas_itv = operational.rinex_start_end(rnx_base, True)

    log.info(
        "Rover: {} | Base: {} | Start: {} | End: {} | Interval: {} s".format(
            rov_name, bas_name, rov_srt, rov_end, rov_itv
        )
    )

    # READ GENERIC CONF FILE & MODIFY IT
    cfgdic = _read_cfg(cfgfile_generik)

    if solformat:
        cfgdic["out-solformat"] = solformat.lower()

    if sateph:
        cfgdic["pos1-sateph"] = sateph.lower()

    if posmode:
        cfgdic["pos1-posmode"] = str(posmode).lower()

    def _edit_cfgdic_posi(rnx_inp, xyz_inp, ant_n=1):
        antobj, recobj, siteobj, locobj = files_rw.read_rinex_2_dataobjts(rnx_inp)
        n = str(ant_n)
        cfgdic[f"ant{n}-postype"] = "xyz"
        cfgdic[f"ant{n}-anttype"] = antobj.Antenna_Type
        if xyz_inp[0] != 0:
            cfgdic[f"ant{n}-pos1"] = xyz_inp[0]
            cfgdic[f"ant{n}-pos2"] = xyz_inp[1]
            cfgdic[f"ant{n}-pos3"] = xyz_inp[2]
        else:
            cfgdic[f"ant{n}-pos1"] = locobj.X_coordinate_m
            cfgdic[f"ant{n}-pos2"] = locobj.Y_coordinate_m
            cfgdic[f"ant{n}-pos3"] = locobj.Z_coordinate_m
        cfgdic[f"ant{n}-antdelu"] = antobj.Up_Ecc
        cfgdic[f"ant{n}-antdeln"] = antobj.North_Ecc
        cfgdic[f"ant{n}-antdele"] = antobj.East_Ecc

    # Edit position in config
    if rover_auto_conf:
        _edit_cfgdic_posi(rnx_rover, xyz_rover, ant_n=1)
    if base_auto_conf:
        _edit_cfgdic_posi(rnx_base, xyz_base, ant_n=2)

    if not (bas_srt <= rov_srt <= rov_end <= bas_end):
        log.warning(
            "rover/base epoch inconsistency: not bas_srt <= rov_srt <= rov_end <= bas_end !!!"
        )

    # write conf file
    _write_cfg(cfgdic, out_conf_fil)

    #### PRODUCTS - DOWNLOAD OR SELECT BEST FROM INPUT LISTS (MULTIPROCESSING) ####
    if download_prods:  ## direct download: best method when no multiprocessing
        orbclklis_use, brdclis_use = dl_prods(prod_dir, bas_srt, igs_prods)
    else:  ### use input lists and find best products: best method when
        # multiprocessing to avoid download concurrency issues
        orbclklis_use = get_best_prods(orbclklis_inp, bas_srt, igs_prods)
        brdclis_use = get_best_prods(brdclis_inp, bas_srt, "BRDC", brdc_mode=True)
        # brdclis_use = brdclis_inp

    # Copy best products to tmp directory
    orbclklis_ok = [_prods2tmp(orb, tmp_dir_wrk) for orb in orbclklis_use]
    brdclis_ok = [_prods2tmp(n, tmp_dir_wrk) for n in brdclis_use]

    # Command
    ### INTERNAL CHOICE: customs args are passed through moddified conf file
    ### her are mandatory args only
    arg_config = "-k " + out_conf_fil
    arg_interval = "-ti " + str(np.round(rov_itv, 6))
    arg_mode = ""
    arg_resultfile = "-o " + out_res_fil

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
            " ".join(orbclklis_ok),
        )
    )

    log.info(
        "Rover: {} | Base: {} | Start: {} | End: {} | Interval: {} s".format(
            rov_name, bas_name, rov_srt, rov_end, rov_itv
        )
    )
    log.info(bigcomand)

    # THIS IS THE PART THAT GETS PARALLELIZED
    subprocess.call([bigcomand], executable="/bin/bash", shell=True)
    if not os.path.isfile(out_res_fil) or os.path.getsize(out_res_fil) < 2000:
        log.error(f"RTKLIB failed for {exp_full_name} :(")
    else:
        utils.gzip_compress(out_res_fil + ".stat", rm_inp=True)
        out_prq_fil = out_res_fil.replace(".out", ".parquet")
        if not os.path.isfile(out_prq_fil):
            df_out2prq = files_rw.read_rtklib(out_res_fil, return_df=True)
            df_out2prq.to_parquet(out_prq_fil, engine='auto')
        log.info("RTKLIB RUN OK for {} :)".format(exp_full_name))

    if not keep_tmp:
        shutils.rmtree(tmp_dir_wrk, ignore_errors=True)
        os.remove(out_res_fil.replace(".out", "") + "_events.pos")

    return out_res_fil


def rtklib_run_pair(
    rinex_pairs,
    cfgfile_generik,
    out_dir,
    tmp_dir=None,
    prod_dir=None,
    orbclklis_inp=None,
    brdclis_inp=None,
    igs_prods="GRG0OPSFIN",
    exp_prefix="",
    xyz_dic=None,
    posmode=None,
    solformat=None,
    sateph=None,
    force=False,
    keep_tmp=False,
    procs=4,
    exe_path="rnx2rtkp",
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
        List of (rnx_rov, rnx_bas) file path tuples to process.
    cfgfile_generik : str | os.PathLike
        Path to generic RTKLIB configuration file.
    out_dir : str | os.PathLike
        Directory where results are saved.
    tmp_dir : str | os.PathLike | None, default=None
        Temporary directory for intermediate files.
    prod_dir : str | os.PathLike | None, default=None
        Directory where to search for orbits, clocks, and BRDC files.
    orbclklis_inp : list of str, optional
        List of orbit/clock files to select from if not downloading.
    brdclis_inp : list of str, optional
        List of BRDC files to select from if not downloading.
    igs_prods : str, default="GRG0OPSFIN"
        Analysis center/product identifier.
    procs : int, default=4
        Maximum number of parallel workers.
    exp_prefix : str, default=""
        Prefix added to output file stems.
    rover_auto_conf : bool, default=False
        If True, reads rover RINEX header for antenna configuration.
    base_auto_conf : bool, default=True
        If True, reads base RINEX header for antenna configuration.
    xyz_dic : dict, optional
        [X, Y, Z] coordinates (dict value) for each station (dict key).
        If None, uses [0, 0, 0] for all.
    solformat : str, default=None
        Output solution format.
    force : bool, default=False
        If True, forces reprocessing even if output exists.
    keep_tmp : bool, default=False
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
    >>> results = rtklib_run_pair(
    ...     pairs,
    ...     '/path/to/config.conf',
    ...     '/path/to/output',
    ...     procs=4
    ... )
    """
    # Directory structure setup
    out_dir = utils.create_dir(out_dir)

    if tmp_dir is None:
        tmp_dir = os.path.join(out_dir, "TMP")
    tmp_dir = utils.create_dir(tmp_dir)

    if prod_dir is None:
        prod_dir = tmp_dir
    else:
        prod_dir = utils.create_dir(prod_dir)

    log.info(
        f"Starting parallel RTKLIB processing for {len(rinex_pairs)} pairs with {procs} workers"
    )

    if len(rinex_pairs) == 0:
        log.error("No RINEX pairs provided. Exiting.")
        return []

    xyz_rovers = []
    xyz_bases = []

    #  Preliminary - determine xyz dictionnaries
    log.info("STEP 1: Determining time spans for all RINEX pairs...")
    for rnx_rov, rnx_bas in rinex_pairs:
        # Prepare coordinate lists
        site_rov = os.path.basename(rnx_rov)[0:9]
        site_bas = os.path.basename(rnx_bas)[0:9]
        if not xyz_dic:
            xyz_dic = {}  # create dummy dict if None provided
        xyz_rov = xyz_dic[site_rov] if site_rov in xyz_dic.keys() else [0, 0, 0]
        xyz_bas = xyz_dic[site_bas] if site_bas in xyz_dic.keys() else [0, 0, 0]
        xyz_rovers.append(xyz_rov)
        xyz_bases.append(xyz_bas)

    #  Parallel processing
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=procs) as executor:
        # Submit all tasks
        future_to_pair = {}
        for i, (rnx_rov, rnx_bas) in enumerate(rinex_pairs):
            future = executor.submit(
                rtklib_run_mono,
                rnx_rover=rnx_rov,
                rnx_base=rnx_bas,
                cfgfile_generik=cfgfile_generik,
                out_dir=out_dir,
                tmp_dir=tmp_dir,
                download_prods=False,  # Products already downloaded in STEP 2
                orbclklis_inp=orbclklis_inp,
                brdclis_inp=brdclis_inp,
                igs_prods=igs_prods,
                exp_prefix=exp_prefix,
                rover_auto_conf=True,
                base_auto_conf=True,
                xyz_rover=xyz_rovers[i],
                xyz_base=xyz_bases[i],
                posmode=posmode,
                solformat=solformat,
                sateph=sateph,
                force=force,
                keep_tmp=keep_tmp,
                exe_path=exe_path,
            )
            future_to_pair[future] = (rnx_rov, rnx_bas)

        # Collect results as they complete
        for future in concurrent.futures.as_completed(future_to_pair):
            pair = future_to_pair[future]
            try:
                result = future.result()
                results.append(result)
                log.info(
                    f"Pair done: {os.path.basename(pair[0])} + {os.path.basename(pair[1])} -> {os.path.basename(result)}"
                )
            except Exception as exc:
                log.error(f"Failed {pair}: {exc}")
                results.append(None)

    log.info(
        f"Parallel RTKLIB processing completed. {sum(1 for r in results if r)} / {len(rinex_pairs)} successful."
    )

    if keep_tmp:
        shutils.rmtree(tmp_dir + "/*", ignore_errors=True)

    return results


def make_pairs(
    rnx_dir,
    sites_rovers,
    sites_bases,
    date_srt=None,
    date_end=None,
):
    """
    Find RINEX files and match rover/base pairs based on time coverage.

    Parameters
    ----------
    rnx_dir : str | os.PathLike
        Directory containing RINEX files.
    sites_rovers : str or list of str
        List of rover site codes.
    sites_bases : str or list of str
        Base site code.
    date_srt : datetime.datetime, optional
        Start date for RINEX file search.
    date_end : datetime.datetime, optional
        End date for RINEX file search.

    Returns
    -------
    tuple
        (rnxs_pairs, df_all) where:
        - rnxs_pairs : list of tuples (rover_path, base_path)
        - df_all : pandas DataFrame with all RINEX file information
    """

    sites_rovers = utils.listify(sites_rovers)
    sites_bases = utils.listify(sites_bases)

    rnxs_pairs = []

    sites_all = list(set(sites_rovers + sites_bases))

    rnxs_all = operational.rinex_finder(
        rnx_dir,
        specific_sites=sites_all,
        start_epoch=date_srt,
        end_epoch=date_end,
    )

    df_all = operational.rinex_table_from_list(rnxs_all, site9_col=True)
    df_all["date_end"] = df_all["date"] + df_all["per"]
    df_rovers = df_all[df_all["site9"].isin(sites_rovers)]
    df_bases = df_all[df_all["site9"].isin(sites_bases)]

    if len(df_bases) == 0:
        log.error(f"No base RINEX files found for site: {sites_bases}")
        log.error(f"All sites found: {str(list(df_all["site9"].unique()))}")
        return [], df_all

    bas_prev, bas_next = "", ""

    for _, row_rov in df_rovers.iterrows():
        rov = row_rov["site9"]
        for bas, df_bas in df_bases.groupby("site9"):
            if bas == rov:
                # this test is to silent multiple warning messages
                if bas != bas_prev or bas != bas_next:
                    log.warning(f"Rover '{rov}' is the same as base '{bas}', skipping pair")
                    base_prev, bas_next = bas, rov
                continue
            rov_srt = row_rov["date"]
            rov_end = row_rov["date_end"]

            sel = (df_bas["date"] <= rov_srt) & (rov_end <= df_bas["date_end"])
            df_bas_sel = df_bas[sel]
            if len(df_bas_sel) == 0:
                log.warning(f"no base for rover {rov} at {rov_srt}")
                continue
            elif len(df_bas_sel) > 1:
                log.warning(f"multi. bases for {rov} at {rov_srt}, get 1st:")
                log.warning(f"\n{df_bas_sel.to_string()}")

            row_bas = df_bas_sel.iloc[0]

            rnxs_pair = (row_rov["path"], row_bas["path"])
            rnxs_pairs.append(rnxs_pair)

    return rnxs_pairs, df_all


def rtklib_run(
    rnx_dir,
    cfgfile_generik,
    sites_rovers,
    sites_bases,
    out_dir,
    tmp_dir=None,
    prod_dir=None,
    igs_prods="GRG0OPSFIN",
    exp_prefix="",
    date_srt=None,
    date_end=None,
    xyz_dic=None,
    posmode=None,
    solformat=None,
    sateph=None,
    force=False,
    keep_tmp=False,
    procs=8,
    exe_path="rnx2rtkp",
):

    log.info("STEP 1: Finding RINEX files and matching rover/base pairs")
    rnxs_pairs, df_all = make_pairs(
        rnx_dir,
        sites_rovers,
        sites_bases,
        date_srt=date_srt,
        date_end=date_end,
    )

    log.info("STEP 2: Downloading shared products (SP3/CLK and BRDC)")
    try:
        orbclklis_shr, brdclis_shr = dl_prods(
            prod_dir, list(df_all["date"].unique()), igs_prods
        )
    except Exception as e:
        log.error(f"Failed to download shared products: {e}")
        raise

    log.info(f"STEP 3: Running {len(rnxs_pairs)} RTKLIB processes in parallel")
    return operational.rtklib_run_pair(
        rnxs_pairs,
        cfgfile_generik,
        out_dir=out_dir,
        tmp_dir=tmp_dir,
        prod_dir=prod_dir,
        orbclklis_inp=orbclklis_shr,
        brdclis_inp=brdclis_shr,
        igs_prods=igs_prods,
        exp_prefix=exp_prefix,
        xyz_dic=xyz_dic,
        posmode=posmode,
        solformat=solformat,
        sateph=sateph,
        force=force,
        keep_tmp=keep_tmp,
        procs=procs,
        exe_path=exe_path,
    )
