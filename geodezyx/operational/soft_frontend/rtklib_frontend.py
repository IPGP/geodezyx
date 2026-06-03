#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.operational contains functions to run the
GNSS processing software RTKLIB.

it can be imported directly with:
from geodezyx import operational

The geodezyx toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/IPGP/geodezyx
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
from os import PathLike

import numpy as np
import shutil as shutils
import pandas as pd
import pyarrow.parquet as pq
from pandas.core.interchange.dataframe_protocol import DataFrame
from tqdm import tqdm

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
        and os.path.getsize(out_res_fil) > 25000  # empiric failed file limit
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
        orbclklis_use = get_best_prods(
            orbclklis_inp, bas_srt, prod_ac_name=igs_prods, brdc_mode=False
        )
        brdclis_use = get_best_prods(
            brdclis_inp, bas_srt, prod_ac_name="BRDC", brdc_mode=True
        )

    # Copy best products to tmp directory
    orbclklis_ok = [_prods2tmp(orb, tmp_dir_wrk) for orb in orbclklis_use]
    brdclis_ok = [_prods2tmp(nav, tmp_dir_wrk) for nav in brdclis_use]

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
        df_out2prq = files_rw.read_rtklib(out_res_fil, return_df=True)
        df_out2prq.to_parquet(out_prq_fil, engine="auto")

        resmpl = "01min"
        out_csv_fil = out_res_fil.replace(".out", "_" + resmpl + ".csv")
        _resample_df(df_out2prq, resmpl).to_csv(out_csv_fil, index=False)

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
    cfgfile_generik : str or os.PathLike
        Path to generic RTKLIB configuration file.
    out_dir : str or os.PathLike
        Directory where results are saved.
    tmp_dir : str or os.PathLike or None, default=None
        Temporary directory for intermediate files.
    prod_dir : str or os.PathLike or None, default=None
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
    rnx_dir : str or os.PathLike
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
    rnxs_pairs : list of tuples
        rover_path, base_path
    df_all : pandas DataFrame
        all RINEX file information
    """

    if sites_bases is not None:
        sit_rov_use = utils.listify(sites_rovers)
        sit_bas_use = utils.listify(sites_bases)
        pairs_use = [(r, b) for r in sit_rov_use for b in sit_bas_use]
    else:
        sit_rov_use, sit_bas_use = zip(*sites_rovers)
        pairs_use = sites_rovers

    rnxs_pairs = []

    sit_all = list(set(sit_rov_use + sit_bas_use))

    rnxs_all = operational.rinex_finder(
        rnx_dir,
        specific_sites=sit_all,
        start_epoch=date_srt,
        end_epoch=date_end,
    )

    df_all = operational.rinex_table_from_list(rnxs_all, site9_col=True)
    df_all["date_end"] = df_all["date"] + df_all["per"]
    df_rovers = df_all[df_all["site9"].isin(sit_rov_use)]
    df_bases = df_all[df_all["site9"].isin(sit_bas_use)]

    if len(df_bases) == 0:
        log.error(f"No base RINEX files found for site: {sit_bas_use}")
        sites_fnd = str(list(df_all["site9"].unique()))
        log.error(f"All sites found: {sites_fnd}")
        return [], df_all

    bas_prev, rov_prev = "", ""

    def _find_bas4rov(row_rov, df_bas_inp):
        """Find matching base for a rover row using time coverage."""
        rov_srt = row_rov["date"]
        rov_end = row_rov["date_end"]
        rov_sit = row_rov["site9"]
        sel = (df_bas_inp["date"] <= rov_srt) & (rov_end <= df_bas_inp["date_end"])
        df_bas_sel = df_bas_inp[sel]

        if len(df_bas_sel) == 0:
            log.warning(f"no base for rover {rov_sit} at {rov_srt}")
            return None
        elif len(df_bas_sel) > 1:
            log.warning(f"multi. bases for {rov_sit} at {rov_srt}, get 1st")
            log.warning(f"\n{df_bas_sel.to_string()}")

        row_bas = df_bas_sel.iloc[0]
        return row_rov["path"], row_bas["path"]

    for rov, bas in pairs_use:
        if rov == bas:
            # this test is to silent multiple warning messages
            if bas != bas_prev or rov != rov_prev:
                wm = f"Rover '{rov}' is the same as base '{bas}', skipping pair"
                log.warning(wm)
                bas_prev, rov_prev = bas, rov
            continue

        df_rov = df_rovers[df_rovers["site9"] == rov]
        df_bas = df_bases[df_bases["site9"] == bas]

        # Use apply to process each rover row
        pairs_for_rov = df_rov.apply(lambda r: _find_bas4rov(r, df_bas), axis=1)
        # Add non-None pairs to the list
        rnxs_pairs.extend([pair for pair in pairs_for_rov if pair is not None])

    return rnxs_pairs, df_all


def rtklib_merge_parquet(
    parquet_inp,
    exp_prefix="",
    fast_merge=False,
    rtklib_out_files=None,
):
    """
    Merge individual RTKLIB parquet files into a single consolidated parquet file.

    Parameters
    ----------
    parquet_inp : str or os.PathLike or list of str
        Either a directory path (all ``*.parquet`` files inside are collected
        recursively) **or** an explicit list of parquet file paths.
        The merged output file is written to the directory (or, for a list,
        to the directory of the first file in the list).
    exp_prefix : str, default=""
        Prefix used to name the merged output file (<exp_prefix>_all.parquet).
    fast_merge : bool, default=False
        If True, only merges the parquet files corresponding to
        ``rtklib_out_files`` (or those in the explicit list) and appends them
        to an already-existing ``<exp_prefix>_all.parquet`` file.
        If False, scans the whole directory recursively for parquet files.
    rtklib_out_files : list of str, optional
        List of ``.out`` file paths produced by a previous RTKLIB run.
        Only used when ``fast_merge=True`` and ``parquet_inp`` is a directory,
        to avoid a full recursive scan.

    Returns
    -------
    str
        Path to the merged parquet file.
    """
    # --- resolve source files and output directory ---
    if isinstance(parquet_inp, (str, os.PathLike)) and os.path.isdir(parquet_inp):
        prq_out_dir = str(parquet_inp)
        if fast_merge and rtklib_out_files:
            l_prq = [f.replace(".out", ".parquet") for f in rtklib_out_files]
            l_prq = [f for f in l_prq if os.path.exists(f)]
        else:
            l_prq = utils.find_recursive(prq_out_dir, "*parquet")
    else:
        # parquet_inp is an explicit list of parquet files
        l_prq = list(parquet_inp)
        prq_out_dir = os.path.dirname(os.path.abspath(l_prq[0])) if l_prq else "."

    prq_path_out = os.path.join(prq_out_dir, exp_prefix + "_all.parquet")
    prq_path_tmp = prq_path_out + ".tmp"

    # Exclude the output file itself and stray temp files from the source list
    l_prq = [
        f for f in l_prq if not f.endswith("_all.parquet") and not f.endswith(".tmp")
    ]

    # When fast-merging, prepend the existing merged file so it is streamed
    # first; write to a temp path to avoid reading and writing the same file.
    if fast_merge and os.path.exists(prq_path_out):
        l_prq_merge = [prq_path_out] + l_prq
        prq_path_wrk = prq_path_tmp
    else:
        l_prq_merge = l_prq
        prq_path_wrk = prq_path_out

    def _drop_pandas_meta(tbl):
        """Drop the 'pandas' metadata key so all tables share the same schema."""
        meta = {k: v for k, v in tbl.schema.metadata.items() if k != b"pandas"}
        return tbl.replace_schema_metadata(meta)

    # Stream each source table directly through a ParquetWriter —
    # no pandas conversion, no in-memory concat.
    writer = None
    try:
        pbar = tqdm(l_prq_merge, desc="Merging parquet", unit="file")
        for f in pbar:
            pbar.set_postfix_str(os.path.basename(f), refresh=False)
            tbl = _drop_pandas_meta(pq.read_table(f))
            if tbl.num_columns == 0:
                log.warning(f"Skipping empty/corrupt parquet file: {f}")
                continue
            if writer is None:
                writer = pq.ParquetWriter(prq_path_wrk, tbl.schema)
            writer.write_table(tbl)
    finally:
        if writer:
            writer.close()

    # Atomically replace the previous merged file when using a temp path
    if prq_path_wrk == prq_path_tmp and os.path.exists(prq_path_tmp):
        os.replace(prq_path_tmp, prq_path_out)

    log.info(f"Merged parquet saved to {prq_path_out}")
    return prq_path_out


def _resample_df(df_inp: DataFrame, sample: str = "15min"):
    """
    Resample a DataFrame with GNSS position data to a specified time interval.

    This helper function resamples DataFrame containing GNSS solution positions (x, y, z)
    to a coarser time resolution using median aggregation. It removes any duplicate
    entries resulting from the resampling operation.

    Parameters
    ----------
    df_inp : pandas.DataFrame
        Input DataFrame with an 'epoch' column (datetime or datetime-like) and
        position columns 'x', 'y', 'z' (numeric).
    sample : str, default="15min"
        Resampling interval as a pandas time offset string.
        Examples: "1min", "15min", "1H" (1 hour), "1D" (1 day).

    Returns
    -------
    pandas.DataFrame
        Resampled DataFrame with:
        - 'epoch' column (reset from index)
        - 'x', 'y', 'z' columns containing median values over each resampling interval
        - No duplicate rows

    Notes
    -----
    - Uses median aggregation to provide robust resampling (resistant to outliers)
    - Expects the input DataFrame to have an 'epoch' column with datetime values
    - The original epoch index is reset in the output
    """
    # ...existing code...
    df_epo = df_inp.set_index("epoch")
    df_med = df_epo[["x", "y", "z"]].resample(sample).median()
    df_out = df_med.reset_index(inplace=False)
    df_out = df_out.drop_duplicates(inplace=False)
    return df_out


def parquet2csv(
    prq_inp: str | PathLike, out_dir: str | PathLike, sample: str = "15min"
):
    """
    Convert merged RTKLIB parquet file to CSV format, processed by rover/base pairs.

    Reads a merged parquet file containing GNSS solutions from multiple rover/base
    station pairs, resamples each pair's data independently, and exports to CSV files.
    This function uses PyArrow filters to minimize memory usage by reading only the
    data relevant to each rover/base pair.

    Parameters
    ----------
    prq_inp : str or os.PathLike
        Path to the merged parquet file containing rover/base pair GNSS solutions.
        The file must contain columns: 'epoch', 'rover', 'base', 'x', 'y', 'z'.
    out_dir : str or os.PathLike
        Output directory where CSV files will be saved. Created if it doesn't exist.
    sample : str, default="15min"
        Resampling interval for position data. Passed to _resample_df().
        Examples: "1min", "15min", "1H", "1D"

    Returns
    -------
    None

    Output Files
    ------------
    CSV files in `out_dir` with naming pattern:
        {rover}_{base}_{sample}.csv

    Each CSV contains columns:
        - epoch: datetime of the resampled position
        - x, y, z: median position coordinates for the resampling interval

    Notes
    -----
    - Uses PyArrow filter expressions to read only necessary data from the parquet file
    - Efficiently handles large parquet files by filtering at read time (minimal RAM usage)
    - Processes each rover/base pair sequentially
    - Prints rover/base pair names to console as they are processed
    """
    # Read the merged parquet file to get unique rover/base combinations
    df_rovbas = pd.read_parquet(
        prq_inp, engine="auto", columns=["rover", "base"]
    ).drop_duplicates()

    # Process each unique rover/base pair
    for irow, (rov, bas) in df_rovbas.iterrows():
        log.info("loading rover/base: %s/%s", rov, bas)

        # Define PyArrow filters to read only this rover/base pair
        # Filters minimize memory usage by selecting data at read time
        filters = [
            ("rover", "==", rov),
            ("base", "==", bas),
        ]

        # Read only the filtered data for this pair
        df_grp = pd.read_parquet(prq_inp, engine="auto", filters=filters)

        # Resample to the specified time interval
        df_out = _resample_df(df_grp, sample)

        # Export to CSV file with naming convention: rover_base_sample.csv
        out_path = f"{out_dir}/{rov}_{bas}_{sample}.csv"
        log.info("saving %s resampled data to CSV: %s", sample, out_path)
        df_out.to_csv(out_path, index=False)

    return None


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
    fast_parquet_merge=False,
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
    out_run_pairs = operational.rtklib_run_pair(
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

    # merge all parquet files into one
    log.info("STEP 4: Merging individual parquet files into one")
    rtklib_merge_parquet(
        out_dir,
        exp_prefix=exp_prefix,
        fast_merge=fast_parquet_merge,
        rtklib_out_files=out_run_pairs,
    )

    return out_run_pairs
