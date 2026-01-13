#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 05/01/2026 19:09:20

@author: psakic
"""

import os
import htmllistparse
import geodezyx.operational.download_utils as dlutils
import re

#### Import the logger
import logging
log = logging.getLogger("geodezyx")

def list_rnx_eost(year_start=None, year_end=None):
    """
    Lists the paths of RINEX files available on the EOS Strasbourg server within a specified year range.

    Parameters
    ----------
    year_start : int, optional
        The starting year for filtering RINEX files. If None, no lower bound is applied.
        Default is None.
    year_end : int, optional
        The ending year for filtering RINEX files. If None, no upper bound is applied.
        Default is None.

    Returns
    -------
    list of str
        A list of URLs pointing to the RINEX files that match the specified year criteria.

    Notes
    -----
    - The function fetches directory listings recursively from the server.
    - Only directories matching the pattern of a 4-digit year are considered.
    - Logs the number of RINEX files found for each day.

    Examples
    --------
    >>> rnx_files = list_rnx_eost(year_start=2020, year_end=2021)
    >>> len(rnx_files) > 0
    True
    """
    url_main = "http://loading.u-strasbg.fr/SPOTGINS/TEST/rinex/"
    cwd, listing = htmllistparse.fetch_listing(url_main, timeout=30)
    rnx_path_lst = []
    for item in listing:
        url_year = os.path.join(url_main,item.name)
        year = item.name.strip('/')
        if not re.match(r'^\d{4}/$', item.name):
            continue
        else:
            year = int(year)
        if year_start is not None:
            if year < year_start:
                continue
        if year_end is not None:
            if year > year_end:
                continue
        cwd, listing = htmllistparse.fetch_listing(url_year, timeout=30)
        for subitem in listing:
            url_day = os.path.join(url_year,subitem.name)
            day = subitem.name.strip('/')
            cwd, listing = htmllistparse.fetch_listing(url_day, timeout=30)
            for file in listing:
                rnx_path = os.path.join(url_day,file.name)
                rnx_path_lst.append(rnx_path)
            log.info(f"{len(listing)} RINEXs found for day {year}-{day}")

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

