#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 23/04/2024 21:31:45

@author: psakic
"""
import datetime as dt
import glob
import itertools
import logging
import os
import pathlib
import re

import pandas as pd

import geodezyx.operational.download_utils as dlutils
from geodezyx import conv

log = logging.getLogger("geodezyx")


def _rnx_obs_rgx(stat, date):
    """
    Generate RINEX observation file regex patterns for both RINEX2 and RINEX3 formats.

    This internal function creates regex patterns to match RINEX observation files
    for a given station and date, supporting both the legacy 8.3 filename convention
    (RINEX2) and the long filename convention (RINEX3).

    Parameters
    ----------
    stat : str
        4-character GNSS station name (e.g., 'ZIMM', 'TLSE').
        Will be converted to lowercase for RINEX2 pattern.
    date : datetime.datetime
        Date for which to generate the regex patterns.

    Returns
    -------
    tuple of (str, str)
        A tuple containing:
        - rnx2rgx : str
            Regex pattern for RINEX2 observation files (8.3 format)
        - rnx3rgx : str
            Regex pattern for RINEX3 observation files (long format)

    Notes
    -----
    The RINEX2 pattern uses the station name in lowercase with a wildcard
    for the file type extension. The RINEX3 pattern follows the standard
    long filename convention with:
    - Country code: wildcard ("...")
    - Data source: wildcard (".")
    - File period: 01D (daily)
    - Data frequency: wildcard ("...")
    - Data type: .O (observation)
    - Format/compression: wildcard (".*")

    Examples
    --------
    ```
    >>> import datetime as dt
    >>> rnx2_pattern, rnx3_pattern = _rnx_obs_rgx('ZIMM', dt.datetime(2020, 1, 1))
    >>> # rnx2_pattern might be: 'zimm001a.20o.*'
    >>> # rnx3_pattern might be: 'ZIMM...._R_20200010000_01D_....O.*'
    ```
    """
    rnx2rgx = conv.statname_dt2rinexname(stat.lower(), date, rnxtype=".*")
    rnx3rgx = conv.statname_dt2rinexname_long(
        stat,
        date,
        country="...",
        data_source=".",
        file_period="01D",
        data_freq="...",
        data_type=".O",
        format_compression=".*",
    )
    return rnx2rgx, rnx3rgx


def _rnx_nav_rgx(stat, date, sys=".", data_source="."):
    """
    Generate RINEX navigation file regex patterns for both RINEX2 and RINEX3 formats.

    This internal function creates regex patterns to match RINEX navigation files
    for a given station and date, supporting both the legacy 8.3 filename convention
    (RINEX2) and the long filename convention (RINEX3). Navigation files contain
    satellite ephemeris and clock correction data.

    Parameters
    ----------
    stat : str
        4-character GNSS station name (e.g., 'BRDC', 'ZIMM').
        For broadcast navigation files, typically use 'BRDC'.
    date : datetime.datetime
        Date for which to generate the regex patterns.
    sys : str, optional
        GNSS system identifier for RINEX3 navigation files. Default is "." (wildcard).
        Common values:
        - "G" : GPS navigation data
        - "R" : GLONASS navigation data
        - "E" : Galileo navigation data
        - "C" : BeiDou navigation data
        - "M" : Mixed/Multi-GNSS navigation data
        - "." : Wildcard for any system
    data_source : str, optional
        Data source identifier for RINEX3 long filenames. Default is "." (wildcard).
        Common values:
        - "R" : Real-time data
        - "S" : Survey/static data
        - "U" : Ultra-rapid data
        - "." : Wildcard for any source

    Returns
    -------
    tuple of (str, str)
        A tuple containing:
        - rnx2rgx : str
            Regex pattern for RINEX2 navigation files (8.3 format)
        - rnx3rgx : str
            Regex pattern for RINEX3 navigation files (long format)

    Notes
    -----
    The RINEX2 pattern uses the station name in lowercase with a wildcard
    for the file type extension. The RINEX3 pattern follows the standard
    long filename convention for navigation files with:
    - Country code: wildcard ("...")
    - Data source: configurable (default wildcard)
    - File period: 01D (daily)
    - Data frequency: empty ("") for navigation files
    - Data type: sys + "N" (Navigation)
    - Format/compression: wildcard (".*")

    Examples
    --------
    >>> import datetime as dt
    >>> rnx2_pattern, rnx3_pattern = _rnx_nav_rgx('BRDC', dt.datetime(2020, 1, 1))
    >>> # rnx2_pattern might be: 'brdc001a.20n.*'
    >>> # rnx3_pattern might be: 'BRDC...._._20200010000_01D_.N.*'

    >>> # For GPS-specific navigation files
    >>> rnx2_gps, rnx3_gps = _rnx_nav_rgx('BRDC', dt.datetime(2020, 1, 1), sys="G", data_source="R")
    >>> # rnx3_gps might be: 'BRDC...._R_20200010000_01D_GN.*'
    """
    rnx2rgx = conv.statname_dt2rinexname(stat.lower(), date, rnxtype=".*")
    rnx3rgx = conv.statname_dt2rinexname_long(
        stat,
        date,
        country="...",
        data_source=data_source,
        file_period="01D",
        data_freq="",
        data_type=sys + "N",
        format_compression=".*",
    )

    return rnx2rgx, rnx3rgx


def _generic_server(stat, date, urlserver, urlsuffix=None):
    """
    Generate RINEX file URLs for a generic FTP server structure.

    This internal function creates download URLs for both RINEX2 and RINEX3 observation
    files following a standard FTP server directory structure: server/year/doy/filename.
    This is the most common layout used by many GNSS data centers.

    Parameters
    ----------
    stat : str
        4-character GNSS station name (e.g., 'ZIMM', 'TLSE').
        Station name case will be handled appropriately for each RINEX format.
    date : datetime.datetime
        Date for which to generate the URLs.
        Used to construct the year/doy directory path and filename patterns.
    urlserver : str
        Base FTP server URL (e.g., 'ftp://example.com/data/').
        Should include the protocol and base path to the RINEX data directory.

    Returns
    -------
    dict
        Dictionary mapping RINEX version numbers to their respective URLs:
        - urldic[2] : str - URL for RINEX2 observation file
        - urldic[3] : str - URL for RINEX3 observation file

    Notes
    -----
    The function follows the standard GNSS data center directory structure:
    - urlserver/YYYY/DOY/filename
    Where:
    - YYYY is the 4-digit year
    - DOY is the 3-digit day of year (001-366)
    - filename follows RINEX2 or RINEX3 naming conventions

    This generic structure is used by servers like IGS SOPAC, IGN, and others.
    Some servers (like CDDIS) have custom structures and require specialized functions.

    Examples
    --------
    >>> import datetime as dt
    >>> urls = _generic_server('ZIMM', dt.datetime(2020, 1, 15), 'ftp://example.com/data/')
    >>> # urls[2] might be: 'ftp://example.com/data/2020/015/zimm015a.20o.*'
    >>> # urls[3] might be: 'ftp://example.com/data/2020/015/ZIMM...._R_20200150000_01D_....O.*'
    """
    rnx2rgx, rnx3rgx = _rnx_obs_rgx(stat, date)

    if not urlsuffix:
        urlsuffix = ""

    ### generate urls
    urldir = str(os.path.join(urlserver, str(date.year), conv.dt2doy(date), urlsuffix))
    rnx2url = os.path.join(urldir, rnx2rgx)
    rnx3url = os.path.join(urldir, rnx3rgx)

    ### generate output urldic, key 2 and 3 are for rinex version
    urldic = dict()
    urldic[2] = rnx2url
    urldic[3] = rnx3url

    return urldic

def igs_sopac_server(stat, date):
    # plante si trop de requete
    urlserver = "ftp://garner.ucsd.edu/pub/rinex/"
    urldic = _generic_server(stat, date, urlserver)
    return urldic


def igs_cddis_server(stat, date):
    # plante si trop de requete
    urlserver = "ftp://gdc.cddis.eosdis.nasa.gov/gps/data/daily/"

    # we can not use _generic_server here because of the specific server path structure

    ### generate regex
    rnx2rgx, rnx3rgx = _rnx_obs_rgx(stat, date)

    ### generate urls
    urldir = os.path.join(
        urlserver, str(date.year), conv.dt2doy(date), date.strftime("%y") + "d"
    )
    rnx2url = os.path.join(urldir, rnx2rgx)
    rnx3url = os.path.join(urldir, rnx3rgx)

    ### generate output urldic, key 2 and 3 are for rinex version
    urldic = dict()
    urldic[2] = rnx2url
    urldic[3] = rnx3url

    return urldic


def igs_ign_server(stat, date):
    # plante si trop de requete
    urlserver = "ftp://igs.ign.fr/pub/igs/data/"
    urldic = _generic_server(stat, date, urlserver)
    return urldic


def igs_ign_ensg_server(stat, date):
    # plante si trop de requete
    urlserver = "ftp://igs.ensg.eu/pub/igs/data/"
    urldic = _generic_server(stat, date, urlserver)
    return urldic


def igs_bkg_server(stat, date):
    urlserver = "ftp://igs-ftp.bkg.bund.de/IGS/obs/"
    # ftp://igs-ftp.bkg.bund.de/IGS/obs/2024/082/IGS00WRD_R_20240820000_01D_MN.rnx.gz

    urldic = _generic_server(stat, date, urlserver)
    return urldic


def nav_rob_server(stat, date):
    urlserver = "ftp://epncb.oma.be/pub/obs/BRDC/"

    # can not use _generic_server here because of the specific server path structure / file name

    ### generate regex
    rnx2rgx, rnx3rgx = _rnx_nav_rgx(stat, date)  ### NAV RNX HERE !!!

    ### generate urls
    urldir = os.path.join(urlserver, str(date.year))  ## NO DOY FOR THIS ONE !!!
    rnx2url = os.path.join(urldir, rnx2rgx)
    rnx3url = os.path.join(urldir, rnx3rgx)

    ### generate output urldic, key 2 and 3 are for rinex version
    urldic = {}
    urldic[2] = rnx2url
    urldic[3] = rnx3url

    return urldic


def sonel_server(stat, date):
    urlserver = "ftp://ftp.sonel.org/gps/data/"
    urldic = _generic_server(stat, date, urlserver)
    return urldic


def rgp_server(stat, date):
    urlserver = "ftp://rgpdata.ign.fr/pub/data/"
    urldic = _generic_server(stat, date, urlserver, "data_30")
    return urldic


def rgp_ensg_server(stat, date):
    urlserver = "ftp://rgpdata.ensg.eu/pub/data/"
    urldic = _generic_server(stat, date, urlserver, "data_30")
    return urldic

def spotgins_eost_server(stat, date):
    urlserver = "http://loading.u-strasbg.fr/SPOTGINS/TEST/rinex/"
    urldic = _generic_server(stat, date, urlserver)
    return urldic

def euref_server(stat, date):
    urlserver = "ftp://epncb.oma.be/pub/obs/"

    # can not use _generic_server here because of upper case RINEX 2 names
    rnx2rgx, rnx3rgx = _rnx_obs_rgx(stat, date)
    rnx2rgx = rnx2rgx.upper()

    ### generate urls
    urldir = str(os.path.join(urlserver, str(date.year), conv.dt2doy(date)))
    rnx2url = os.path.join(urldir, rnx2rgx)
    rnx3url = os.path.join(urldir, rnx3rgx)

    ### generate output urldic, key 2 and 3 are for rinex version
    urldic = dict()
    urldic[2] = rnx2url
    urldic[3] = rnx3url

    return urldic


def nav_bkg_server(stat, date):
    urlserver = "ftp://igs-ftp.bkg.bund.de/IGS/BRDC/"
    # ftp://igs-ftp.bkg.bund.de/IGS/BRDC/2024/082/BRDC00WRD_S_20240820000_01D_MN.rnx.gz

    ### generate regex
    rnx2rgx, rnx3rgx = _rnx_nav_rgx(
        stat, date, sys="M", data_source="S"
    )  ### NAV RNX HERE !!!

    ### generate urls
    urldir = os.path.join(urlserver, str(date.year), conv.dt2doy(date))
    rnx3url = os.path.join(urldir, rnx3rgx)

    urldic = dict()
    urldic[3] = rnx3url

    return urldic


############ not adapted yet after april 2024 mods
def igs_cddis_nav_server_legacy(stat, date):
    # table_proto privilegier
    urlserver = "ftp://cddis.gsfc.nasa.gov/gps/data/daily/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date, "n.Z")
    url = os.path.join(
        urlserver, str(date.year), conv.dt2doy(date), date.strftime("%y") + "n", rnxname
    )
    return url


def rgp_ign_smn_1_hz_server_legacy(stat, date):
    urlserver = "ftp://rgpdata.ign.fr/pub/data/"

    urls = []

    for h in range(24):
        date_session = date
        date_session = date_session.replace(hour=h)

        log.info("%s session %s", date_session, h)
        rnxname = conv.statname_dt2rinexname(
            stat.lower(), date_session, session_a_instead_of_daily_session=True
        )
        url = os.path.join(
            urlserver, str(date.year), conv.dt2doy(date), "data_1", rnxname
        )

        urls.append(url)

    return urls


def unavco_server_legacy(stat, date):
    urlserver = "ftp://data-out.unavco.org/pub/rinex"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, "obs", str(date.year), conv.dt2doy(date), rnxname)
    return url


def renag_server_legacy(stat, date):
    urlserver = "ftp://renag.unice.fr/data/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), rnxname)
    return url


def uwiseismic_server_legacy(stat, date, user="", passwd=""):
    urlserver = "ftp://www2.uwiseismic.com/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, "rinex", str(date.year), conv.dt2doy(date), rnxname)
    return url, user, passwd


def orpheon_server_legacy(stat, date, user="", passwd=""):
    urlserver = "ftp://renag.unice.fr/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), rnxname)
    return url, user, passwd


def ovsg_server_legacy(stat, date, user="", passwd=""):
    if dt.datetime(2009, 1, 1) <= date <= dt.datetime(2014, 2, 10):
        urlserver = (
            "http://webobs.ovsg.univ-ag.fr/rawdata/GPS-GPSDATA.backtemp_20140210/"
        )
    else:
        urlserver = "http://webobs.ovsg.univ-ag.fr/rawdata/GPS/GPSDATA/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), "rinex", rnxname)
    return url, user, passwd


def geoaus_server_legacy(stat, date):
    """Geosciences Australia
    ex : ftp://ftp.ga.gov.au/geodesy-outgoing/gnss/data/daily/2010/10063/"""
    urlserver = "ftp://ftp.ga.gov.au/geodesy-outgoing/gnss/data/daily/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(
        urlserver, str(date.year), date.strftime("%y") + conv.dt2doy(date), rnxname
    )
    return url


def ens_fr_legacy(stat, date):
    urlserver = "ftp://gnss.ens.fr/pub/public/crl/GPS/rinex/"
    rnxname = conv.statname_dt2rinexname(stat.lower(), date)
    url = os.path.join(urlserver, str(date.year), conv.dt2doy(date), rnxname)
    return url


def _server_select(datacenter, site, curdate):
    mode1hz = False
    secure_ftp = False
    urldic = dict()
    if datacenter in ("igs_cddis", "igs"):
        urldic = igs_cddis_server(site, curdate)
        secure_ftp = True
    elif datacenter == "igs_sopac":
        urldic = igs_sopac_server(site, curdate)
    elif datacenter == "igs_ign":
        urldic = igs_ign_server(site, curdate)
    elif datacenter == "igs_ign_ensg":
        urldic = igs_ign_ensg_server(site, curdate)
    elif datacenter == "sonel":
        urldic = sonel_server(site, curdate)
    elif datacenter == "euref":
        urldic = euref_server(site, curdate)
    elif datacenter == "igs_bkg":
        urldic = igs_bkg_server(site, curdate)
    elif datacenter in ("nav", "brdc"):
        urldic = nav_rob_server(site, curdate)
    elif datacenter in ("nav_rt", "brdc_rt"):
        urldic = nav_bkg_server(site, curdate)
    elif datacenter == "rgp":
        urldic = rgp_server(site, curdate)
    elif datacenter == "rgp_ensg":
        urldic = rgp_ensg_server(site, curdate)
    elif datacenter == "spotgins_eost":
        urldic = spotgins_eost_server(site, curdate)
    # elif datacenter == 'rgp':
    #     urldic = rgp_ign_smn_server_legacy(site, curdate)
    # elif datacenter == 'rgp_mlv':
    #     urldic = rgp_ign_mlv_server(site, curdate)
    # elif datacenter == 'rgp_1Hz':
    #     urls = rgp_ign_smn_1_hz_server_legacy(site, curdate)
    #     mode1hz = True
    # elif datacenter == 'renag':
    #     urldic = renag_server(site, curdate)
    # elif datacenter == 'orpheon':
    #     urldic = orpheon_server_legacy(site, curdate)
    # elif datacenter == 'uwiseismic':
    #     urldic = uwiseismic_server(site, curdate)
    # elif datacenter == 'ovsg':
    #     urldic = ovsg_server_legacy(site, curdate)
    # elif datacenter == 'unavco':
    #     urldic = unavco_server_legacy(site, curdate)
    # elif datacenter == 'geoaus':
    #     urldic = geoaus_server_legacy(site, curdate)
    # elif datacenter == 'ens_fr':
    #     urldic = ens_fr_legacy(site, curdate)
    else:
        log.warning("unkwn server dic in the dico, skip ...")
        return None, None, None

    return urldic, secure_ftp, mode1hz

def effective_save_dir(parent_archive_dir, stat, date, archtype="stat"):
    """
    INTERNAL_FUNCTION

    archtype =
        stat
        stat/year
        stat/year/doy
        year/doy
        year/stat
        week/dow
        OR only '/' for a dirty saving in the parent folder
        ... etc ..."""
    if archtype == "/":
        return parent_archive_dir

    if len(archtype) > 0 and archtype.startswith("/"):
        log.warning("The archive type starts with /, remove it to avoid error")

    out_save_dir = parent_archive_dir
    fff = archtype.split("/")
    year = str(date.year)
    doy = conv.dt2doy(date)
    _, _ = year, doy  ## simply to remove the unused linter warning...
    week, dow = conv.dt2gpstime(date)
    for f in fff:
        out_save_dir = os.path.join(out_save_dir, eval(f))
    return out_save_dir


def rnx_regex_indir(rnx_regex, dir_files_list):
    """
    Match files in a directory against a given regex pattern.

    Parameters
    ----------
    rnx_regex : str
        Regex pattern to match filenames.
    dir_files_list : list of str
        List of filenames in the directory.

    Returns
    -------
    str or None
        The first matching filename, or None if no match is found.
    """
    matches = [file for file in dir_files_list if re.search(rnx_regex, file)]
    return matches[0] if matches else None


def crawl_ftp_files(
    table,
    sftp="auto",
    user=None,
    passwd=None,
    path_ftp_crawled_files_save=None,
    path_all_ftp_files_save=None,
    force=False,
):
    """
    Crawl FTP servers to find available RINEX files and update download table.

    This function performs an optimized FTP crawl by reusing connections and
    minimizing directory changes. It checks for existing local files, connects
    to FTP servers, lists remote files, and updates the table with availability
    status and actual file URLs.

    Parameters
    ----------
    table : pd.DataFrame
        Input table containing RINEX download metadata with columns:
        - 'host': FTP server hostname
        - 'dir': Remote directory path
        - 'outdir': Local output directory
        - 'rnxrgx': RINEX filename regex pattern
        - 'sftp': Boolean indicating if SFTP should be used
        - 'crawled': Boolean indicating if already crawled
    sftp : str or bool, optional
        SFTP mode setting. Default is 'auto'.
        - 'auto': Use the 'sftp' column value from each table row
        - True/False: Force SFTP on/off for all connections
    user : str, optional
        FTP username. Default is None (anonymous).
    passwd : str, optional
        FTP password. Default is None (anonymous).
    path_ftp_crawled_files_save : str, optional
        Path to save the crawled files table as CSV.
        If None, no file is saved.
    path_all_ftp_files_save : str, optional
        Path to save all discovered FTP files as CSV.
        If None, no file is saved.
    force : bool, optional
        Force re-download even if files exist locally. Default is False.

    Returns
    -------
    tuple of (pd.DataFrame, pd.Series, pd.Series)
        - table_use : pd.DataFrame
            Updated table with crawl results, including new columns:
            - 'ok_dwl': Boolean indicating file is available for download
            - 'ok_loc': Boolean indicating file exists locally
            - 'rnxnam': Actual filename found on server
            - 'url_true': Complete FTP URL for download
        - all_ftp_files : pd.Series
            All files discovered on FTP servers with full URLs
        - all_loc_files : pd.Series
            Local file paths for files that already exist

    Notes
    -----
    The function implements several optimizations:
    - Reuses FTP connections when possible (same host)
    - Reconnects every 50 operations to avoid timeouts
    - Caches local and remote directory listings
    - Only changes directories when necessary
    - Saves intermediate results for recovery

    The crawling process:
    1. Checks for existing local files first
    2. Connects to FTP server when host changes
    3. Lists remote directory contents when directory changes
    4. Matches files using regex patterns
    5. Updates table with availability status
    6. Generates download URLs for available files

    Examples
    --------
    >>> import pandas as pd
    >>> table = pd.DataFrame({
    ...     'host': ['ftp.example.com'],
    ...     'dir': ['/data/2020/001'],
    ...     'outdir': ['/local/data'],
    ...     'rnxrgx': ['station001a.20o.*'],
    ...     'sftp': [False],
    ...     'crawled': [False]
    ... })
    >>> crawled_table, all_files, local_files = crawl_ftp_files(table)
    """

    def _save_crawled_files(table_inp):
        """Save crawled files table to CSV if path is provided."""
        if path_ftp_crawled_files_save:
            table_inp.to_csv(path_ftp_crawled_files_save)

    def _get_and_save_all_ftp_files(all_ftp_files_stk_inp):
        """Concatenate and save all discovered FTP files."""
        if all_ftp_files_stk_inp:
            all_ftp_files_out = pd.concat(all_ftp_files_stk_inp)
            all_ftp_files_out.reset_index(drop=True, inplace=True)
        else:
            all_ftp_files_out = pd.Series([], dtype=str)

        if path_all_ftp_files_save:
            all_ftp_files_out.to_csv(path_all_ftp_files_save)

        return all_ftp_files_out

    table_use = table.copy()

    # Initialize loop variables
    prev_host = ""  # Track previous host to reuse connections
    prev_dir = ""  # Track previous directory to avoid unnecessary changes
    ftp_files_list = []  # Cache of current directory file listing
    all_ftp_files_stk = []  # Stack to collect all discovered files
    local_files_lis = []  # Cache of current local directory listing
    count_loop = 0  # Counter for connection refresh
    count_nmax = 50  # Maximum operations before reconnecting
    ftpobj = None  # Current FTP connection object

    for irow, row in table_use.iterrows():
        # Skip rows already crawled
        if row["crawled"]:
            continue

        # Check local files when directory changes
        if prev_dir != row["outdir"] or irow == 0:
            local_files_lis = glob.glob(row["outdir"] + "/*")

        # Check if file already exists locally
        rnxlocal = rnx_regex_indir(row["rnxrgx"], local_files_lis)

        # Skip if file exists locally and not forcing download
        if rnxlocal and os.path.getsize(rnxlocal) > 0:
            rnxloc_bn = os.path.basename(rnxlocal)
            if not force:
                log.info("%s already exists locally ;)", rnxloc_bn)
                table_use.loc[irow, "ok_loc"] = True
                table_use.loc[irow, "ok_dwl"] = False
                table_use.loc[irow, "rnxnam"] = os.path.basename(rnxlocal)
                continue
            else:  # force mode
                log.info(
                    "%s already exists locally, but re-download forced",
                    rnxloc_bn,
                )
                table_use.loc[irow, "ok_loc"] = False
                table_use.loc[irow, "ok_dwl"] = True
                table_use.loc[irow, "rnxnam"] = os.path.basename(rnxlocal)

        count_loop += 1

        # Create new FTP connection if needed
        if row["host"] != prev_host or count_loop > count_nmax or count_loop == 1:
            if ftpobj:
                ftpobj.close()

            # Determine SFTP mode: use row value if 'auto', otherwise use parameter
            sftp_use = bool(row["sftp"]) if sftp == "auto" else sftp
            ftpobj, _ = dlutils.ftp_objt_create(
                secure_ftp_inp=sftp_use,
                host=row["host"],
                user=user,
                passwd=passwd,
            )
            prev_host = row["host"]

            # Save intermediate results and reset counter on reconnection
            if count_loop > count_nmax:
                count_loop = 0
                _save_crawled_files(table_use)
                _get_and_save_all_ftp_files(all_ftp_files_stk)

        # Get file list when directory changes
        if prev_dir != row["dir"] or irow == 0:
            log.info("chdir " + row["dir"])
            ftpobj.cwd("/")  # Reset to root directory first

            try:
                ftpobj.cwd(row["dir"])
                ftp_files_list = dlutils.ftp_dir_list_files(ftpobj)

                # Save all FTP files for reporting with full URLs
                ftp_files_series = pd.Series(ftp_files_list)
                ftp_files_series = (
                    os.path.join("ftp://", row["host"], row["dir"])
                    + "/"
                    + ftp_files_series
                )
                all_ftp_files_stk.append(ftp_files_series)

            except Exception as e:
                log.warning("unable to chdir to %s, exception %s", row["dir"], e)
                ftp_files_list = []

            prev_dir = row["dir"]

        # Check if file exists on server using regex pattern
        rnx_return = rnx_regex_indir(row["rnxrgx"], ftp_files_list)

        # Update table based on file availability
        if rnx_return:
            table_use.loc[irow, "rnxnam"] = rnx_return
            table_use.loc[irow, "ok_dwl"] = True
            log.info(rnx_return + " found on server :)")
        else:
            table_use.loc[irow, "rnxnam"] = ""
            table_use.loc[irow, "ok_dwl"] = False
            log.warning(row["rnxrgx"] + " not found on server :(")

        table_use.loc[irow, "crawled"] = True

    # Generate URLs for downloadable files (available remotely but not locally)
    rnx_ok_dwl = table_use["ok_dwl"] & ~table_use["ok_loc"]
    table_use.loc[rnx_ok_dwl, "url_true"] = table_use.loc[rnx_ok_dwl].apply(
        lambda x: os.path.join("ftp://", x["host"], x["dir"], x["rnxnam"]), axis=1
    )

    # Save final results
    _save_crawled_files(table_use)
    all_ftp_files = _get_and_save_all_ftp_files(all_ftp_files_stk)

    # Clean up FTP connection
    if ftpobj:
        ftpobj.close()

    # Generate local file paths for files that exist locally
    rnx_ok_loc = table_use["ok_loc"]
    if rnx_ok_loc.sum() > 0:
        all_loc_files = pd.Series(
            table_use.loc[rnx_ok_loc].apply(
                lambda e: os.path.join(e["outdir"], e["rnxnam"]), axis=1
            )
        )
    else:
        all_loc_files = pd.Series([])

    return table_use, all_ftp_files, all_loc_files


def download_gnss_rinex(
    statdico,
    output_dir,
    startdate,
    enddate,
    archtype='year/doy',
    parallel_download=4,
    user="anonymous",
    passwd="anonymous@isp.com",
    path_ftp_crawled_files_save=None,
    path_ftp_crawled_files_load=None,
    skip_crawl=False,
    path_all_ftp_files_save=None,
    quiet_mode=False,
    final_archive_for_sup_check=None,
    force=False,
    no_rnx2=False,
    no_rnx3=False,
):
    """
    Parameters
    ----------
    statdico : dict
        a statdico is a dictionary associating Archives Centers to list of stations

        Exemple:
            >>> statdico['archive center 1'] = ['STA1','STA200XXX','STA3', ...]
            >>> statdico['archive center 2'] = ['STA200XXX','STA1','STA4', ...]

        the supported archive center are (april 2024):
            * igs_cddis or igs (CDDIS data center)
            * igs_sopac (for the SOPAC/UCSD/SIO data center, but not very reliable)
            * igs_ign (IGN's data center, main server at St Mandé)
            * igs_ign_ensg (IGN's data center, secondary server at ENSG, Marne-la-Vallée)
            * igs_bkg (BKG's IGS data center, for the IGS stations)
            * sonel
            * euref (EPN data center hosted at ROB)
            * nav or brdc as archive center allows to download nav files (using 'BRDC' as station name)
            * from the ROB server, using GOP files
            * nav_rt or brdc_rt as archive center allows to download *real time* nav files from the BKG server
            * rgp (IGN's RGP main server at St Mandé)
            * rgp_ensg (IGN's RGP, secondary server at ENSG, Marne-la-Vallée)
        _ not reimplemented yet _:
            * rgp_1Hz (IGN's RGP, all the 24 hourly rinex for the day will be downloaded)
            * renag
            * ovsg
            * unavco
            * geoaus (Geosciences Australia)
            * ens_fr

    output_dir : str
        the root directory on your local drive were to store the RINEXs

    archtype : str
        string describing how the archive directory is structured, e.g:
            * stat
            * stat/year
            * stat/year/doy
            * year/doy
            * year/stat
            * week/dow/stat
            * etc ...

    user & passwd : str
        user & password for a secure server

    path_ftp_crawled_files_save : str
        will save at the given path (directory+filename) in a CSV file containing
        the list of the existing RINEXs found on the server by the FTP crawler.
        It allows to use this list directly if one face a timeout during
        the download part.

    path_ftp_crawled_files_load : str
        load and use the list of the existing RINEXs found on the FTP server,
        generated by a previous run of the FTP crawler (called by
        download_gnss_rinex or directly by crawl_ftp_files).
        a new call of crawl_ftp_files can be bypassed with
        `skip_crawl`.

    skip_crawl : bool
        when a list of the existing RINEXs found on the FTP server
        is provided with path_ftp_crawled_files_load, skip
        a new call of crawl_ftp_files.
        Then, we assume the provided list as a complete one,
        ready for the download step

    path_all_ftp_files_save : str
        will save at the given path (directory+filename) in a CSV file
        ALL the remote files found on the FTP server.

    quiet_mode : bool
        List the available RINEXs without downloading them.
        Useful only if path_ftp_crawled_files_save is given

    final_archive_for_sup_check : str
        The final archive path or a file containing the archived RINEXs in
        their final destination.
        useful if the final archive is different from archive_dir
        *** not re-implemented yet ***

    force : bool
        Force the download even if the file already exists locally

    no_rnx2 & no_rnx3 : bool
        limit the search/download to RINEX2 (short names) and/or
        RINEX3 (long names) depending on the boolean given

    Returns
    -------
    out_tup_lis : List of tuples
        Returns a list of tuples containing
        the local path of the downloaded file and
        a boolean indicating whether the download was successful.
        e.g. [(local_path1, True), (local_path2, False), ...]

    Minimal exemple
    ---------------
        >>> statdic = dict()
        >>> statdic['igs_cddis'] = ['ZIMM','tlse']
        >>> archive_dir = '/home/USER/test_dl_rnx'
        >>> startdate = dt.datetime(2020,1,1)
        >>> enddate = dt.datetime(2020,1,31)
        >>> geodezyx.operational.download_gnss_rinex(statdic, output_dir, startdate, enddate)
    """

    date_range = conv.dt_range(startdate, enddate)

    log.info("dates: %s to %s", startdate, enddate)

    if path_ftp_crawled_files_load:
        table = pd.read_csv(path_ftp_crawled_files_load)
    else:
        table = gen_crawl_table(
            statdico, date_range, output_dir, archtype, no_rnx2, no_rnx3
        )

    if len(table) == 0:
        log.error("No RINEX files found for the given criteria.")
        return None

    if skip_crawl:
        table_crawl = table
        files_all = pd.Series([], dtype=str)
        files_loc = pd.Series([], dtype=str)
    else:
        table_crawl, files_all, files_loc = crawl_ftp_files(
            table,
            sftp="auto",
            user=user,
            passwd=passwd,
            path_ftp_crawled_files_save=path_ftp_crawled_files_save,
            path_all_ftp_files_save=path_all_ftp_files_save,
            force=force,
        )

    #### get only the valid (true) url
    table_dl = table_crawl.loc[table_crawl["url_true"].dropna().index]

    # Initialize output list
    out_tup_lis = []

    # Only download if files are available and not in quiet mode
    if len(table_dl) == 0:
        log.error(
            "no valid RINEX URL found/selected on the FTP server, check your inputs"
        )
    elif quiet_mode:
        log.warning("quiet mode, no download was performed")
    else:
        out_tup_lis = dlutils.ftp_downld_front(
            table_dl["url_true"].values,
            table_dl["outdir"].values,
            parallel_download=parallel_download,
            secure_ftp=table_dl["sftp"].values,
            user=user,
            passwd=passwd,
            force=force,
        )

    ### add the local paths to the output tuples
    if len(files_loc) > 0:
        loc_tup_lis = [(f, True) for f in files_loc]
        out_tup_lis_fin = loc_tup_lis + out_tup_lis
    else:
        out_tup_lis_fin = out_tup_lis

    log.info(
        "RINEX files fetched: total: %d, downloaded: %d, already here: %d",
        len(out_tup_lis_fin),
        len(out_tup_lis),
        len(files_loc),
    )

    return out_tup_lis


def gen_crawl_table(statdico, date_range, output_dir, archtype, no_rnx2, no_rnx3):
    """
    Generate a crawl table for RINEX file downloads.

    Parameters
    ----------
    statdico : dict
        Dictionary mapping data centers to station lists.
    date_range : list
        List of datetime objects for the date range.
    output_dir : str
        Root output directory.
    archtype : str
        Archive directory structure type.
    no_rnx2 : bool
        Skip RINEX2 files if True.
    no_rnx3 : bool
        Skip RINEX3 files if True.

    Returns
    -------
    pd.DataFrame
        Table with download metadata for each RINEX file.
    """
    table_proto = []

    for datacenter, site_lis in statdico.items():
        log.info("datacenter/stations: %s/%s", datacenter, " ".join(site_lis))

        for date, site in itertools.product(date_range, site_lis):
            urldic, sftp, _ = _server_select(datacenter, site, date)
            if not urldic:
                continue

            outdir = effective_save_dir(output_dir, site, date, archtype)

            for rnxver, rnxurl in urldic.items():
                if (rnxver == 2 and no_rnx2) or (rnxver == 3 and no_rnx3):
                    continue
                table_proto.append((date, site, outdir, rnxver, rnxurl, sftp))

    # Create DataFrame with all collected data
    table = pd.DataFrame(
        table_proto, columns=["date", "site", "outdir", "ver", "url_theo", "sftp"]
    )

    # Add status columns
    table["crawled"] = False
    table["ok_dwl"] = False
    table["ok_loc"] = False
    table["url_true"] = None
    table["rnxnam"] = ""

    # Parse URL components
    urlpaths = table["url_theo"].apply(pathlib.Path)
    table["rnxrgx"] = urlpaths.apply(lambda p: p.name)
    table["host"] = urlpaths.apply(lambda p: p.parts[1])
    table["dir"] = urlpaths.apply(lambda p: os.path.join(*p.parts[2:-1]))

    return table
