#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 18:47:44 2026

@author: psakic
"""

import datetime as dt
import shutil

from geodezyx import utils

#### Import the logger
import logging
import subprocess
import os

from geodezyx import conv
from geodezyx import operational

log = logging.getLogger("geodezyx")


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


def get_prod_date(date_inp, prod_ac_name="", period_stepback=0, latency=3):
    """
    Adjusts product dates based on data center latency and product type.

    Accounts for the delay before data centers provide products after the epoch time.
    Automatically increments latency stepback if the current time is within the latency
    window, then rounds the date according to the product type's cadence.

    Parameters
    ----------
    date_inp : datetime.datetime
        The input date to be adjusted.
    prod_ac_name : str, optional
        The product analysis center name. Determines rounding interval:
        - "ULT": Ultra-rapid (6-hour intervals)
        - "NRT": Near Real-Time (1-hour intervals)
        - "FIN" or "RAP": Final or Rapid (1-day intervals)
        - "BRDC": Broadcast ephemeris (1-day intervals)
        Default is empty string.
    period_stepback : int, optional
        Number of periods to step back from the current time.
        Incremented internally if within the latency window.
        Defaults to 0.
    latency : int, optional
        The latency in hours to consider for stepping back.
        the latency is because the data center does not provide
        the products before a certain time after the epoch of
        the products
        IGS Website: latency is offically 3 hours
        GRG on CDDIS: empircally, ~2 hours and 5 minutes
        If we are within the latency window, we step back one
        period to be able to get a substitute product
        Defaults to 3 hours.

    Returns
    -------
    datetime.datetime
        The adjusted and floored date according to the product cadence.

    Notes
    -----
    - The latency window is disabled for BRDC.
    - Dates are floored (rounded down) to the nearest interval boundary.
    - If `prod_ac_name` is unrecognized, the original date is returned with a warning.
    - ULT and NRT products are designed to be 2 days long,
      so an extra day is subtracted to get the correct product date.

    Examples
    --------
    >>> get_prod_date(dt.datetime(2026, 2, 15, 10, 30), "ULT")
    datetime(2026, 2, 14, 6, 0)
    >>> get_prod_date(dt.datetime(2026, 2, 15, 10, 30), "FIN")
    datetime(2026, 2, 15, 0, 0)
    """
    now = dt.datetime.now(dt.timezone.utc)
    date_inp_utc = date_inp.replace(tzinfo=dt.timezone.utc)

    # the latency is because the data center does not provide the products before a
    # certain time after the epoch of the products
    # IGS Website: latency is 3 hours offically
    # GRG on CDDIS: empircally, ~2 hours and 5 minutes
    # If we are within the latency window, we need to step back one
    # period to be able to get a substitute product
    # Only apply latency adjustment on the initial call (period_stepback == 0)
    # to avoid double-counting latency in subsequent iterations
    if not "BRDC" in prod_ac_name and period_stepback == 0:
        laten_td = dt.timedelta(hours=latency)
        is_during_latency = now - date_inp_utc < laten_td
        if is_during_latency:
            period_stepback += 1

    # rnd_def tuple : (period value, period unit, extra_delta)
    if "ULT" in prod_ac_name:
        rnd_def = (6, "h", 1)
    elif "NRT" in prod_ac_name:
        rnd_def = (1, "h", 1)
    elif "FIN" in prod_ac_name or "RAP" in prod_ac_name:
        rnd_def = (1, "d", 0)
    elif "BRDC" in prod_ac_name:
        rnd_def = (1, "d", 0)
    else:
        log.warning("Unknown prods. '%s', unable to find the right date", prod_ac_name)
        return date_inp

    period_val = rnd_def[0]
    period_unit = rnd_def[1]
    # ULT and NRT are by design 2 days long (1 day is calculated, the other is extrapolated)
    # so we need to substract an extra day to the date to be able to get the right product
    extra_delta = dt.timedelta(days=rnd_def[2])

    rnd_use = str(period_val * (1 + period_stepback)) + str(period_unit)
    date_out_srt = conv.round_dt(date_inp, rnd_use, mode="floor") - extra_delta
    return date_out_srt


get_prod_date(dt.datetime(2026,4,8,12,30,0), "BRDC", period_stepback=1)

def get_dates_fmt(dates_inp, prod_date=True, prod_ac_name=""):
    """
    Normalizes and formats input dates for GNSS product downloading.

    Converts various date input formats to a standardized sorted list. Optionally applies
    product-specific date adjustments to account for data center latency.

    Parameters
    ----------
    dates_inp : datetime, list of datetime, or tuple of datetime
        Input dates in one of three formats:
        - Single datetime: converted to a single-element list
        - List of datetime: used as-is
        - Tuple of (start_date, end_date): generates daily range from start to end
    prod_date : bool, optional
        If True, applies `get_prod_date` to adjust dates for product latency.
        Defaults to True.
    prod_ac_name : str, optional
        Analysis center name passed to `get_prod_date` for latency adjustments.
        Defaults to empty string.

    Returns
    -------
    list of datetime
        Sorted list of unique datetime objects with duplicates removed.

    Raises
    ------
    Exception
        If `dates_inp` is neither a datetime, list, tuple, nor iterable.

    Examples
    --------
    >>> get_dates_fmt(dt.datetime(2026, 2, 15))
    [datetime(2026, 2, 15, 0, 0)]
    >>> get_dates_fmt([dt.datetime(2026, 2, 15), dt.datetime(2026, 2, 15)])
    [datetime(2026, 2, 15, 0, 0)]
    >>> get_dates_fmt((dt.datetime(2026, 2, 15), dt.datetime(2026, 2, 17)))
    [datetime(2026, 2, 15, 0, 0), datetime(2026, 2, 16, 0, 0), datetime(2026, 2, 17, 0, 0)]
    """
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
        dates_lis = [get_prod_date(date, prod_ac_name) for date in dates_lis]
        dates_lis = list(sorted(list(set(dates_lis))))  # remove duplicates

    return dates_lis


def get_best_prods(
    prod_list_inp, date_inp, prod_ac_name="", brdc_mode=False, period_stepback_max=None
):
    """
    Finds the best matching product file(s) for a given date with progressive latency tolerance.

    Searches through available product files to match the requested date, progressively
    increasing latency tolerance if no exact match is found. Falls back to returning all
    products if no match is found after maximum latency attempts.

    Parameters
    ----------
    prod_list_inp : list of str
        List of available product file names.
    date_inp : datetime.datetime
        The target date to match.
    prod_ac_name : str, optional
        Analysis center name (e.g., "ULT", "NRT", "FIN"). Determines maximum latency steps
        if not explicitly provided. Defaults to empty string.
    brdc_mode : bool, optional
        If True, treats files as BRDC (Broadcast), and uses `rinexname2dt` for parsing.
        If False, treats files as SP3/CLK, and uses `sp3name_v3_2dt`.
        Defaults to False.
    period_stepback_max : int, optional
        Maximum number of latency steps to attempt. If None, defaults are:
        - "ULT": 3 steps
        - "NRT": 23 steps
        - Others: 0 steps
        Defaults to None.

    Returns
    -------
    list of str
        List of matching product file names. If no matches found, returns the entire input list.
        For BRDC mode with 2+ matches, attempts to filter for "GOP" products.

    Notes
    -----
    - Performs iterative search: starts with exact date match, then incrementally steps back
      through latency periods to find products.
    - File naming convention interpretation varies by `brdc_mode` setting.
    - Contains potential bug in BRDC filtering logic: `[e for e in best_prod_out if "GOP"][0]`
      may raise IndexError if no products contain "GOP" substring.

    """
    if brdc_mode:
        pan = "BRDC"
    else:
        pan = prod_ac_name

    if period_stepback_max is None:
        if "ULT" in prod_ac_name:
            period_stepback_max = 3
        elif "NRT" in prod_ac_name:
            period_stepback_max = 23
        else:
            period_stepback_max = 0

    best_prod_out = []
    psb = 0
    debug_print = True
    while len(best_prod_out) == 0 and psb <= period_stepback_max:
        date_best = get_prod_date(date_inp, prod_ac_name=pan, period_stepback=psb)
        for prod in prod_list_inp:
            if brdc_mode:
                date_prod = conv.rinexname2dt(prod)
            else:
                date_prod = conv.sp3name_v3_2dt(prod)

            if debug_print:
                log.info(f"current prod: {prod}")
                log.info(f"input date: {date_inp}")
                log.info(f"wished best date: {date_best}")
                log.info(f"current prod date: {date_prod}")
                log.info(f"period stepback/max: {psb}/{period_stepback_max}")

            if date_prod == date_best:
                best_prod_out.append(prod)
        psb += 1

        if debug_print:
            log.info("LOOP END, selected prods: %s", str(best_prod_out))

    if len(best_prod_out) == 0:
        log.warning(
            "No optimal prod. found for %s, latency stepback up to %s, returning all prods.",
            date_inp,
            period_stepback_max,
        )
        best_prod_out = prod_list_inp

    elif brdc_mode and len(best_prod_out) >= 2:
        best_prod_out = [e for e in best_prod_out if "GOP"][0]
    else:
        pass

    return best_prod_out


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
    redwld_delta : int, optional
        The age threshold in hours for re-downloading BRDC files. If a found BRDC
        file is older than this threshold, it will be re-downloaded. Defaults to 4 hours.

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
        now = dt.datetime.now(dt.timezone.utc)
        date_floor_re_utc = date_floor_re.replace(tzinfo=dt.timezone.utc)

        ## We go for GOP non-real time BRDC
        if now - date_floor_re_utc >= dt.timedelta(hours=25):
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

        redl_lim = redwld_delta * 3600
        brdc_mtime = os.path.getmtime(brdc_fnd[0])
        now = now.timestamp()

        if len(brdc_fnd) > 0 and now - brdc_mtime > redl_lim:
            log.info(
                "BRDC %s found but older than %sh, re-downloading...",
                brdc_fnd[0],
                redwld_delta,
            )
            return True, "nav_rt"
        else:
            return False, "nav_rt"

    ######### BROADCAST
    dates_inp = get_dates_fmt(dates_inp, prod_date=True, prod_ac_name="BRDC")
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
    data_centers=("cddis", "esa", "ign", "whu"),
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
        Defaults to ("cddis", "esa", "ign", "whu").

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

    dates_lis = get_dates_fmt(dates_inp, prod_ac_name=prod_ac_name)
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
            ac_names=(prod_ac_name,),
            prod_types=prod_types,
            remove_patterns=("ULA",),
            archtype="year/doy",
            new_name_conv=True,
            parallel_download=1,
            archive_center=data_centers,
            mgex=mgex,
            repro=0,
            return_also_uncompressed_files=True,
            dow_manu=False,
        )

        # if len(prods) >= 5:
        #    log.info("enougth products found: %s", len(prods))
        #    break

    return prods


def dl_orbclk_tite(
    prod_parent_dir, dates_inp, prod_ac_name, login="sakic", password=""
):
    """
    Downloads GRG SP3 and CLK products from CNES's TITE server.
    """

    import geodezyx.operational.gins_runner.gins_bd_update as gins_bd_update

    dates_lis = get_dates_fmt(dates_inp, prod_ac_name=prod_ac_name)
    dates_lis = list(set(dates_lis))  # remove duplicates

    files_list = []
    remote_prod_subdir = ""

    for date in dates_lis:
        year = str(date.year)
        doy = str(conv.dt2doy(date)).zfill(3)
        hour = str(date.hour).zfill(2)

        if "ULT" in prod_ac_name:
            remote_prod_subdir = "GRU"
            file_orb = f"{prod_ac_name}_{year}{doy}{hour}00_02D_05M_ORB.SP3.gz"
            file_clk = f"{prod_ac_name}_{year}{doy}{hour}00_02D_05M_CLK.CLK.gz"
        elif "RAP" in prod_ac_name:
            remote_prod_subdir = "GRR"
            file_orb = f"{prod_ac_name}_{year}{doy}0000_01D_01M_ORB.SP3.gz"
            file_clk = f"{prod_ac_name}_{year}{doy}0000_01D_01M_CLK.CLK.gz"
        elif "FIN" in prod_ac_name:
            remote_prod_subdir = "GRG"
            file_orb = f"{prod_ac_name}_{year}{doy}0000_01D_01M_ORB.SP3.gz"
            file_clk = f"{prod_ac_name}_{year}{doy}0000_01D_01M_CLK.CLK.gz"
        else:
            log.warning("Unknown AC name '%s'", prod_ac_name)
            continue

        files_list.append(file_orb)
        files_list.append(file_clk)

    remote_prod_fulldir = os.path.join(
        "/home/gins/MIROIR_STAF/mesures/gps/orbites/", remote_prod_subdir
    )
    gins_bd_update.download_rsync(
        files_list,
        remote_prod_fulldir,
        prod_parent_dir,
        login,
        "tite.get.obs-mip.fr",
        password,
    )

    loc_files_out = []
    for file in files_list:
        tmp_loc_file = str(os.path.join(prod_parent_dir, file))
        if os.path.isfile(tmp_loc_file):
            date = conv.sp3name_v3_2dt(file)
            year, doy = conv.dt2doy_year(date, reverse_order=True)
            year_doy_dir = str(os.path.join(prod_parent_dir, year, doy))
            utils.create_dir(year_doy_dir)
            final_loc_file = os.path.join(year_doy_dir, file)
            shutil.move(tmp_loc_file, final_loc_file)
            loc_files_out.append(final_loc_file)
        else:
            log.warning("File '%s' not found after download attempt.", file)

    return loc_files_out


def dl_prods(prod_dir, dates_inp, prod_ac_name, download_lock=None, use_tite=False):
    """
    Download shared products (SP3/CLK and BRDC) once for all parallel runs.

    Parameters
    ----------
    prod_dir : str
        Directory where products will be downloaded.
    dates_inp : iterable of datetime
        If list, the list of dates for which the products are to be downloaded.
        If tuple, the start and end dates for which the products are to be downloaded.
    prod_ac_name : str
        Analysis center identifier.
    download_lock : threading.Lock, optional
        Lock to protect concurrent downloads (if needed).
    use_tite : bool, optional
        If True, attempts to download SP3/CLK products from CNES's TITE server
        instead of the default method. Defaults to False.

    Returns
    -------
    tuple
        A tuple containing two lists: (orbclklis_out, brdclis_out)
        - orbclklis_out: List of downloaded SP3/CLK product files.
        - brdclis_out: List of downloaded BRDC files.
    """

    orbclklis_out = []
    brdclis_out = []
    ######### Download SP3/CLK products
    if not use_tite:
        if download_lock:
            download_lock.acquire()
        try:
            prodlis = operational.dl_orbclk(
                prod_dir,
                dates_inp,
                prod_ac_name,
                prod_types=("sp3", "clk"),
                data_centers=("esa",),
            )
        finally:
            if download_lock:
                download_lock.release()

        if not prodlis or not prodlis[0]:
            log.warning("No SP3/CLK product files found remotely nor locally")
            orbclklis_out = []
        else:
            orbclklis_out = prodlis

    ######### Download SP3/CLK products from TITE
    if use_tite:
        log.info("Trying to download SP3/CLK products from CNES's TITE server...")
        try:
            orbclklis_out = dl_orbclk_tite(
                prod_dir, dates_inp, prod_ac_name=prod_ac_name
            )
            if orbclklis_out:
                log.info("SP3/CLK products successfully downloaded from TITE.")
            else:
                log.warning(
                    "No SP3/CLK product files found on TITE after download attempt."
                )
        except Exception as e:
            log.error(f"Error while downloading SP3/CLK products from TITE: {e}")
            orbclklis_out = []

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
        f"Products available: {len(orbclklis_out)} SP3/CLK, {len(brdclis_out)} BRDC"
    )

    return orbclklis_out, brdclis_out
