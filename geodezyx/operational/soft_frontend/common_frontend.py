#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 18:47:44 2026

@author: psakic
"""

import datetime as dt
from geodezyx import utils

#### Import the logger
import logging
import subprocess
import os

import hatanaka

from geodezyx import conv
from geodezyx import files_rw
from geodezyx import operational

log = logging.getLogger("geodezyx")

import re


def run_command(command):
    """
    Runs a shell command and captures both stdout and stderr.

    Parameters
    ----------
    command : str
        The shell command to be executed.

    Notes
    -----
    This function uses subprocess.Popen to run the command in a new process.
    It continuously reads and prints stdout and stderr until the process finishes.
    The function prints the return code of the process once it completes.
    """
    # Run the command and capture both stdout and stderr
    process = subprocess.Popen(
        [command],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        executable="/bin/bash",
        shell=True,
    )

    # Continuously read and print stdout and stderr
    while True:
        # Read a line from stdout
        stdout_line = process.stdout.read().decode("utf-8")
        if stdout_line:
            log.info("STDOUT: %s", stdout_line.strip())
        # Read a line from stderr
        # stderr_line = process.stderr.read1().decode("utf-8")
        # if stderr_line:
        #    print(f"STDERR: {stderr_line.strip()}", end='', flush=True)

        # Check if the process has finished
        return_code = process.poll()
        if return_code is not None:
            log.info("RETURN CODE: %s", return_code)
            break


def _get_prod_date(date_inp, prod_ac_name=""):

    now = dt.datetime.now()

    ult_delta = dt.timedelta(hours=24)

    #### the latency is because the data center does not provide the products before a certain time after the epoch of the products
    latency_dt = dt.timedelta(hours=24)
    is_during_latency = now - date_inp < latency_dt

    if "ULT" in prod_ac_name:
        if is_during_latency:
            before_delta = dt.timedelta(hours=12)
        else:
            before_delta = dt.timedelta(hours=0)
        date_out_srt = conv.round_dt(date_inp - ult_delta - before_delta, "6h", mode="floor")
    elif "NRT" in prod_ac_name:
        if is_during_latency:
            before_delta = dt.timedelta(hours=1)
        else:
            before_delta = dt.timedelta(hours=0)
        date_out_srt = conv.round_dt(date_inp - ult_delta - before_delta, "1h", mode="floor")
    elif "FIN" in prod_ac_name or "RAP" in prod_ac_name:
        date_out_srt = conv.round_dt(date_inp, "1d", mode="floor")
    elif "BRDC" in prod_ac_name:
        date_out_srt = conv.round_dt(date_inp, "1d", mode="floor")
    else:
        date_out_srt = date_inp

    return date_out_srt


def _get_dates_fmt(dates_inp, prod_date=True, prod_ac_name=""):
    if type(dates_inp) is tuple:
        date_srt, date_end = dates_inp
        dates_lis = conv.dt_range(date_srt, date_end, dt.timedelta(days=1))
    elif type(dates_inp) is list:
        dates_lis = dates_inp
    elif not utils.is_iterable(dates_inp):
        dates_lis = [dates_inp]
    else:
        log.error(
            "dates_lis should be a list, a tuple, or a single datetime, but got %s",
            type(dates_inp),
        )
        raise Exception

    if prod_date:
        dates_lis = [_get_prod_date(date, prod_ac_name) for date in dates_lis]
        dates_lis = list(sorted(list(set(dates_lis))))  # remove duplicates

    return dates_lis


def dl_brdc(prod_parent_dir, dates_inp, redwld_delta=4):
    """
    Downloads BRDC (Broadcast Ephemeris) files for PRIDE PPPAR from a given directory and date list.

    Parameters
    ----------
    prod_parent_dir : str
        The parent directory where the products are stored.
    dates_inp : iterable of datetime
        If list, the list of dates for which the BRDC files are to be downloaded.
        If tuple, the start and end dates for which the BRDC files are to be downloaded.

    Returns
    -------
    brdc_path_lis : list
        A list of downloaded BRDC files.

    Notes
    -----
    This function first rounds the dates in the date list to the nearest day and removes duplicates.
    It then downloads the BRDC files for each unique date using the operational.download_gnss_rinex function.
    The downloaded files are appended to a list which is returned at the end.
    """
    brdc_lis_out = []

    def _force_redwld(date_inp):
        date_floor_re = conv.round_dt(date_inp, "1d", mode="floor")
        now = dt.datetime.now()

        ## We go for GOP non-real time BRDC
        if now - date_floor_re >= dt.timedelta(hours=25):
            return False, "nav"

        ## We go for BKG real time BRDC
        fct = conv.statname_dt2rinexname_long
        brdc_fname = fct(
            "BRDC",
            date_floor_re,
            country="WRD",
            data_source="S",
            data_type="MN",
            format_compression="rnx.gz",
        )

        prod_dir_doy = os.path.join(
            prod_parent_dir, *conv.dt2doy_year(date_floor_re, reverse_order=True)
        )
        brdc_fnd = utils.find_recursive(prod_dir_doy, brdc_fname)

        if len(brdc_fnd) > 0 and brdc_fnd[0] < now.timestamp() - redwld_delta * 3600:
            log.info(
                "BRDC %s found but older than %sh, re-downloading...",
                brdc_fnd[0],
                redwld_delta,
            )
            return True, "nav_rt"
        else:
            return False, "nav_rt"

    ######### BROADCAST
    dates_inp = _get_dates_fmt(dates_inp, prod_date=True, prod_ac_name="BRDC")
    dates_inp = list(set(dates_inp))  # remove duplicates
    brdc_tmp_stk = []
    for date in dates_inp:
        force, nav_type = _force_redwld(date)
        date_floor = conv.round_dt(date, "1d", mode="floor")
        brdc_tmp = operational.download_gnss_rinex(
            # nav_type steers the download to either "nav" or "nav_rt" (real-time),
            # depending on the age of the existing file
            {nav_type: ["BRDC"]},
            prod_parent_dir,
            date_floor,
            date_floor,
            archtype="year/doy",
            parallel_download=1,
            force=force,
        )
        brdc_tmp_stk.extend(brdc_tmp)
    # manage the weird output of download_gnss_rinex,
    # which is a list of tuples (path, bool)
    brdc_path_lis, brdc_bool_lis = zip(*brdc_tmp_stk)

    return brdc_path_lis


def dl_orbclk(
    prod_parent_dir,
    dates_inp,
    prod_ac_name,
    prod_types=("sp3", "clk", "bia", "obx", "erp"),
    data_centers=("ign", "whu"),
):
    """
    Downloads GNSS products for PRIDE PPPAR from a given directory and date list.

    Parameters
    ----------
    prod_parent_dir : str
        The parent directory where the products are stored.
    dates_inp : iterable of datetime
        If list, the list of dates for which the products are to be downloaded.
        If tuple, the start and end dates for which the products are to be downloaded.
    prod_ac_name : str
        The name of the analysis center providing the products.
    prod_types : tuple of str, optional
        The types of GNSS products to be downloaded
        (e.g., "sp3", "clk", "bia", "obx", "erp").
        Defaults to ("sp3", "clk", "bia", "obx", "erp").
    data_centers : tuple of str, optional
        The data centers from which to download the products
        (e.g., "ign", "whu").

    Returns
    -------
    list
        A list of downloaded GNSS products.

    Notes
    -----
    This function downloads various GNSS products such as orbits, clocks, biases, etc.
    It iterates over specified data centers and attempts to download the products.
    If at least 5 products are found, the function stops further downloads.
    """

    dates_lis = _get_dates_fmt(dates_inp, prod_ac_name=prod_ac_name)
    dates_lis = list(set(dates_lis))  # remove duplicates

    dl_prods_fct = operational.download_gnss_products

    if type(data_centers) is not tuple:
        data_centers = (data_centers,)

    prods = []
    ######### ORBITS CLOCKS ETC...
    for data_centers in data_centers:
        mgex = True if "MGX" in prod_ac_name else False
        prods = dl_prods_fct(
            prod_parent_dir,
            dates_lis,
            None,
            AC_names=(prod_ac_name,),
            prod_types=prod_types,
            remove_patterns=("ULA",),
            archtype="year/doy",
            new_name_conv=True,
            parallel_download=1,
            archive_center=data_centers,
            mgex=mgex,
            repro=0,
            sorted_mode=False,
            return_also_uncompressed_files=True,
            ftp_download=False,
            dow_manu=False,
        )

        if len(prods) >= 5:
            log.info("enougth products found: %s", len(prods))
            break

    return prods


def dl_prods(prod_dir, dates_inp, calc_center, download_lock=None):
    """
    Download shared products (SP3/CLK and BRDC) once for all parallel runs.

    Parameters
    ----------
    prod_dir : str
        Directory where products will be downloaded.
    dates_inp : iterable of datetime
        If list, the list of dates for which the products are to be downloaded.
        If tuple, the start and end dates for which the products are to be downloaded.
    calc_center : str
        Analysis center identifier.
    download_lock : threading.Lock, optional
        Lock to protect concurrent downloads (if needed).

    Returns
    -------
    tuple
        A tuple containing two lists: (orbclklis_out, brdclis_out)
        - orbclklis_out: List of downloaded SP3/CLK product files.
        - brdclis_out: List of downloaded BRDC files.
    """

    ######### Download SP3/CLK products
    if download_lock:
        download_lock.acquire()
    try:
        prodlis = operational.dl_orbclk(
            prod_dir,
            dates_inp,
            calc_center,
            prod_types=("sp3", "clk"),
            data_centers=("ign",),
        )
    finally:
        if download_lock:
            download_lock.release()

    if not prodlis or not prodlis[0]:
        log.warning("No SP3/CLK product files found remotely nor locally")
        orbclklis_out = []
    else:
        orbclklis_out = prodlis

    ######### Download BRDC files
    if download_lock:
        download_lock.acquire()
    try:
        brdclis = operational.dl_brdc(prod_dir, dates_inp)
    finally:
        if download_lock:
            download_lock.release()

    if len(brdclis) == 0:
        log.error("No BRDC nav file found remotely nor locally")
        raise FileNotFoundError("No BRDC nav file found remotely nor locally")

    brdclis_out = brdclis

    log.info(
        f"Products downloaded: {len(orbclklis_out)} SP3/CLK, {len(brdclis_out)} BRDC"
    )

    return orbclklis_out, brdclis_out

##### bellow are the function for PRIDE PPPAR only


def get_right_brdc(brdc_lis_inp, tmp_dir_inp):
    """
    Selects the appropriate BRDC (Broadcast Ephemeris) file from a list and unzips it.

    Parameters
    ----------
    brdc_lis_inp : list
        A list of BRDC file paths.
    tmp_dir_inp : str
        The directory where the unzipped BRDC file will be stored.

    Returns
    -------
    tuple
        A tuple containing the original BRDC file path and the unzipped BRDC file path.
    """
    if len(brdc_lis_inp) == 1:  ## normal case
        brdc_ori = brdc_lis_inp[0]
        brdc_unzip = files_rw.unzip_gz_z(brdc_ori, out_gzip_dir=tmp_dir_inp)
    elif len(brdc_lis_inp) == 0:
        brdc_ori = None
        brdc_unzip = None
        log.warning("no brdc. found")
    elif len(brdc_lis_inp) > 1:
        log.warning("several brdc found, keep the last one")
        log.warning(brdc_lis_inp)
        brdc_ori = brdc_lis_inp[-1]
        brdc_unzip = files_rw.unzip_gz_z(brdc_ori, out_gzip_dir=tmp_dir_inp)
    else:
        #### this should never happend
        brdc_ori = None
        brdc_unzip = None
        pass

    return brdc_ori, brdc_unzip


def get_best_latency(prod_lis_inp):
    """
    Selects the best latency from a list of product files.

    Internal function for get_right_prod.

    Parameters
    ----------
    prod_lis_inp : list
        A list of product file paths.

    Returns
    -------
    list
        The best latency found in the list.
        But can be several if the latency is the same for several products.

    """
    latency_priority_lis = ["FIN", "RAP", "ULT", "NRT"]

    out_prods_lis = []

    for lat in latency_priority_lis:
        for prod in prod_lis_inp:
            if lat in prod:
                out_prods_lis.append(prod)
                break

    if len(out_prods_lis) == 0:
        log.warning("no prod. found with a known latency")
        out_prods_lis = prod_lis_inp

    return out_prods_lis


def get_right_prod(prod_lis_inp, tmp_dir_inp, prod_name, default_fallback):
    """
    Selects the appropriate product file from a list and unzips it if necessary.

    Parameters
    ----------
    prod_lis_inp : list
        A list of product file paths.
    tmp_dir_inp : str
        The directory where the unzipped product file will be stored.
    prod_name : str
        The name of the product.
    default_fallback : bool
        If True, uses default values if products are not found.

    Returns
    -------
    tuple
        A tuple containing the unzipped product file path and the original product file path.
    """
    if len(prod_lis_inp) == 1:  ## normal case
        prod_ori = prod_lis_inp[0]
        prod_out = files_rw.unzip_gz_z(prod_ori, out_gzip_dir=tmp_dir_inp)
    elif len(prod_lis_inp) == 0 and default_fallback:
        log.warning("no prod. %s found, fallback to 'Default' in cfg file", prod_name)
        prod_out = "Default"
        prod_ori = "Default"
    elif len(prod_lis_inp) == 0 and not default_fallback:
        log.warning(
            "no prod. %s found, no fallback to Default set (default_fallback option), aborting...",
            prod_name,
        )
        prod_out = None
        prod_ori = None
    elif len(prod_lis_inp) > 1:
        log.warning("several prod found, search for the best latency")
        log.warning(prod_lis_inp)
        prod_lis_bst = get_best_latency(prod_lis_inp)

        if len(prod_lis_bst) == 1:
            prod_ori = prod_lis_bst[0]
        else:
            log.warning("several prod. found with the same latency, keep the last one")
            log.warning(prod_lis_bst)
            prod_ori = prod_lis_bst[-1]

        prod_out = files_rw.unzip_gz_z(prod_ori, out_gzip_dir=tmp_dir_inp)
    else:
        #### this should never happend
        prod_out = None
        prod_ori = None
        pass

    return prod_out, prod_ori
