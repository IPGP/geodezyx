#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 05/01/2026 19:09:20

@author: psakic
"""

import os
import geodezyx.operational.download_gnss.download_utils as dlutils
import geodezyx.operational.download_gnss.download_rinex as dlrinex
import geodezyx.conv as conv
import re
from datetime import datetime

#### Import the logger
import logging

log = logging.getLogger("geodezyx")


def list_rnx_eost(year_start=None, year_end=None, excl_rnx2_post2022=True):
    """
    Lists the paths of RINEX files available on the EOST Strasbourg server within a specified year range.

    This is a wrapper around `list_rnx_http` specifically for the EOST server.

    Parameters
    ----------
    year_start : int, optional
        The starting year for filtering RINEX files. If None, no lower bound is applied.
        Default is None.
    year_end : int, optional
        The ending year for filtering RINEX files. If None, no upper bound is applied.
        Default is None.
    excl_rnx2_post2022 : bool, optional
        If True, excludes RINEX version 2 files for years 2022 and later
        (as RINEX version 3 is available). Default is True.

    Returns
    -------
    list of str
        A list of URLs pointing to the RINEX files that match the specified year criteria.

    Notes
    -----
    - Uses the generic `list_rnx_http` function internally.
    - The function fetches directory listings recursively from the server.
    - Logs the number of RINEX files found for each day.

    Examples
    --------
    >>> rnx_files = list_rnx_eost(year_start=2020, year_end=2021)
    >>> len(rnx_files) > 0
    True
    """
    url_eost = "http://loading.u-strasbg.fr/SPOTGINS/TEST/rinex/"

    # Convert year_start and year_end to dates
    date_start = None
    date_end = None
    if year_start is not None:
        date_start = datetime(year_start, 1, 1).date()
    if year_end is not None:
        date_end = datetime(year_end, 12, 31).date()

    # Get the list of files from the generic function
    rnx_path_lst = dlrinex.list_rnx_http(url_eost, date_start, date_end)

    # Apply EOST-specific filtering for RINEX v2 post-2022
    if excl_rnx2_post2022:
        filtered_lst = []
        for rnx_path in rnx_path_lst:
            # Extract filename from URL path
            filename = rnx_path.split("/")[-1]
            # Extract year from URL path (parent directory)
            year = int(rnx_path.split("/")[-3])

            isrnx2 = conv.rinex_regex_search_tester(
                filename, short_name=True, long_name=False
            )
            if year >= 2022 and isrnx2:
                log.info(f"Exclude RNX2 {filename} (>= 2022, RNX3 is available).")
                continue
            filtered_lst.append(rnx_path)
        return filtered_lst

    return rnx_path_lst


def get_rnx_eost(outdir, year_start=None, year_end=None):
    """
    Downloads RINEX files from the EOS Strasbourg server to a local directory.

    Parameters
    ----------
    outdir : str
        The output directory where RINEX files will be saved. Subdirectories
        will be created following the year/day structure from the server.
    year_start : int, optional
        The starting year for filtering RINEX files. If None, no lower bound is applied.
        Default is None.
    year_end : int, optional
        The ending year for filtering RINEX files. If None, no upper bound is applied.
        Default is None.

    Returns
    -------
    None

    Notes
    -----
    - Creates subdirectories in the output directory matching the year/day structure.
    - Uses `list_rnx_eost` to retrieve the list of RINEX files.
    - Downloads files using HTTP protocol via `dlutils.download_http`.

    Examples
    --------
    >>> get_rnx_eost("/path/to/output", year_start=2020, year_end=2021)
    """
    lst_rnx = list_rnx_eost(year_start, year_end)
    for url in lst_rnx:
        day_suffix = "/".join(url.split("/")[-3:-1])
        outdir_day = os.path.join(outdir, day_suffix)
        if not os.path.isdir(outdir_day):
            os.makedirs(outdir_day)
        _ = dlutils.download_http(url, outdir_day)
    return None


def list_rnx_igs_brdc(date_start=None, date_end=None):
    """
    Lists the paths of broadcast ephemeris (BRDC) files available on the IGS BKG server.

    This is a wrapper around `list_rnx_http` specifically for the IGS BRDC repository.

    Parameters
    ----------
    date_start : str or datetime.date, optional
        The starting date for filtering files in format 'YYYY-MM-DD' or as a datetime.date object.
        If None, no lower bound is applied. Default is None.
    date_end : str or datetime.date, optional
        The ending date for filtering files in format 'YYYY-MM-DD' or as a datetime.date object.
        If None, no upper bound is applied. Default is None.

    Returns
    -------
    list of str
        A list of URLs pointing to the broadcast ephemeris files.

    Notes
    -----
    - Uses the IGS BRDC (Broadcast Ephemeris) repository at BKG.
    - Works with both RINEX 2 and RINEX 3 formats.
    - Files are organized by year/day-of-year structure.

    Examples
    --------
    >>> brdc_files = list_rnx_igs_brdc(date_start="2023-01-01", date_end="2023-01-31")
    >>> len(brdc_files) > 0
    True
    """
    url_igs_brdc = "https://igs.bkg.bund.de/root_ftp/IGS/BRDC/"
    return dlrinex.list_rnx_http(url_igs_brdc, date_start, date_end)

