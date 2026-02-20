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


def _rnx_obs_rgx(date, site=None):
    """
    Generate RINEX observation file regex patterns for both RINEX2 and RINEX3 formats.

    This internal function creates regex patterns to match RINEX observation files
    for a given station and date, supporting both the legacy 8.3 filename convention
    (RINEX2) and the long filename convention (RINEX3).

    Parameters
    ----------
    date : datetime.datetime
        Date for which to generate the regex patterns.
    site : str, optional
        4-character GNSS station name (e.g., 'ZIMM', 'TLSE').
        Will be converted to lowercase for RINEX2 pattern.
        If None or empty string, will match all stations.

    Returns
    -------
    rnx2rgx : str
        Regex pattern for RINEX2 observation files (8.3 format)
    rnx3rgx : str
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
    """
    # If site is None or empty, use wildcards to match all stations

    if date.hour > 0:
        per = "01H"  # daily file
    else:
        per = "01D"

    if not site:
        site = "...."  # 4-character wildcard for station

    rnx2rgx = conv.statname_dt2rinexname(site.lower(), date, rnxtype=".*")
    rnx3rgx = conv.statname_dt2rinexname_long(
        site,
        date,
        country="...",
        data_source=".",
        file_period=per,
        data_freq="...",
        data_type=".O",
        format_compression=".*",
    )
    return str(rnx2rgx), str(rnx3rgx)


def _rnx_nav_rgx(date, site=None, sys=".", data_source="."):
    """
    Generate RINEX navigation file regex patterns for both RINEX2 and RINEX3 formats.

    This internal function creates regex patterns to match RINEX navigation files
    for a given station and date, supporting both the legacy 8.3 filename convention
    (RINEX2) and the long filename convention (RINEX3). Navigation files contain
    satellite ephemeris and clock correction data.

    Parameters
    ----------
    date : datetime.datetime
        Date for which to generate the regex patterns.
    site : str, optional
        4-character GNSS station name (e.g., 'BRDC', 'ZIMM').
        For broadcast navigation files, typically use 'BRDC'.
        If None or empty string, will match all stations.
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
    rnx2rgx : str
        Regex pattern for RINEX2 navigation files (8.3 format)
    rnx3rgx : str
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
    >>> rnx2_pattern, rnx3_pattern = _rnx_nav_rgx(dt.datetime(2020, 1, 1), 'BRDC')
    >>> # rnx2_pattern might be: 'brdc001a.20n.*'
    >>> # rnx3_pattern might be: 'BRDC...._._20200010000_01D_.N.*'

    >>> # For GPS-specific navigation files
    >>> rnx2_gps, rnx3_gps = _rnx_nav_rgx(dt.datetime(2020, 1, 1), 'BRDC', sys="G", data_source="R")
    >>> # rnx3_gps might be: 'BRDC...._R_20200010000_01D_GN.*'
    """
    # If site is None or empty, use wildcards to match all stations
    if not site:
        site = "...."  # 4-character wildcard for station

    rnx2rgx = conv.statname_dt2rinexname(site.lower(), date, rnxtype=".*")
    rnx3rgx = conv.statname_dt2rinexname_long(
        site,
        date,
        country="...",
        data_source=data_source,
        file_period="01D",
        data_freq="",
        data_type=sys + "N",
        format_compression=".*",
    )

    return str(rnx2rgx), str(rnx3rgx)


def _generic_server(date, site=None, urlserver=None, urlsuffix=None):
    """
    Generate RINEX file URLs for a generic FTP server structure.

    This internal function creates download URLs for both RINEX2 and RINEX3 observation
    files following a standard FTP server directory structure: server/year/doy/filename.
    This is the most common layout used by many GNSS data centers.

    Parameters
    ----------
    date : datetime.datetime
        Date for which to generate the URLs.
        Used to construct the year/doy directory path and filename patterns.
    site : str, optional
        4-character GNSS station name (e.g., 'ZIMM', 'TLSE').
        Station name case will be handled appropriately for each RINEX format.
        If None or empty string, will match all stations in the directory.
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
    >>> urls = _generic_server(dt.datetime(2020, 1, 15), 'ZIMM', 'ftp://example.com/data/')
    >>> # urls[2] might be: 'ftp://example.com/data/2020/015/zimm015a.20o.*'
    >>> # urls[3] might be: 'ftp://example.com/data/2020/015/ZIMM...._R_20200150000_01D_....O.*'
    """
    rnx2rgx, rnx3rgx = _rnx_obs_rgx(date, site)

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

def igs_sopac_server(date, site=None):
    # plante si trop de requete
    urlserver = "ftp://garner.ucsd.edu/pub/rinex/"
    urldic = _generic_server(date, site, urlserver)
    return urldic


def igs_cddis_server(date, site=None):
    # plante si trop de requete
    urlserver = "ftp://gdc.cddis.eosdis.nasa.gov/gps/data/daily/"

    # we can not use _generic_server here because of the specific server path structure

    ### generate regex
    rnx2rgx, rnx3rgx = _rnx_obs_rgx(date, site)

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


def igs_ign_server(date, site=None):
    # plante si trop de requete
    urlserver = "ftp://igs.ign.fr/pub/igs/data/"
    urldic = _generic_server(date, site, urlserver)
    return urldic


def igs_ign_ensg_server(date, site=None):
    # plante si trop de requete
    urlserver = "ftp://igs.ensg.eu/pub/igs/data/"
    urldic = _generic_server(date, site, urlserver)
    return urldic


def igs_bkg_server(date, site=None):
    urlserver = "ftp://igs-ftp.bkg.bund.de/IGS/obs/"
    # ftp://igs-ftp.bkg.bund.de/IGS/obs/2024/082/IGS00WRD_R_20240820000_01D_MN.rnx.gz

    urldic = _generic_server(date, site, urlserver)
    return urldic


def nav_rob_server(date, site=None):
    urlserver = "ftp://epncb.oma.be/pub/obs/BRDC/"

    # can not use _generic_server here because of the specific server path structure / file name

    ### generate regex
    rnx2rgx, rnx3rgx = _rnx_nav_rgx(date, site, sys="M", data_source="R")  ### NAV RNX HERE !!!

    ### generate urls
    urldir = os.path.join(urlserver, str(date.year))  ## NO DOY FOR THIS ONE !!!
    rnx2url = os.path.join(urldir, rnx2rgx)
    rnx3url = os.path.join(urldir, rnx3rgx)

    ### generate output urldic, key 2 and 3 are for rinex version
    urldic = {2: rnx2url, 3: rnx3url}

    return urldic


def sonel_server(date, site=None):
    urlserver = "ftp://ftp.sonel.org/gps/data/"
    urldic = _generic_server(date, site, urlserver)
    return urldic


def rgp_server(date, site=None):
    urlserver = "ftp://rgpdata.ign.fr/pub/data/"
    urldic = _generic_server(date, site, urlserver, "data_30")
    return urldic


def rgp_ensg_server(date, site=None):
    urlserver = "ftp://rgpdata.ensg.eu/pub/data/"
    urldic = _generic_server(date, site, urlserver, "data_30")
    return urldic

def spotgins_eost_server(date, site=None):
    urlserver = "http://loading.u-strasbg.fr/SPOTGINS/TEST/rinex/"
    urldic = _generic_server(date, site, urlserver)
    return urldic

def euref_server(date, site=None):
    urlserver = "ftp://epncb.oma.be/pub/obs/"

    # can not use _generic_server here because of upper case RINEX 2 names
    rnx2rgx, rnx3rgx = _rnx_obs_rgx(date, site)
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


def nav_bkg_server(date, site=None):
    urlserver = "ftp://igs-ftp.bkg.bund.de/IGS/BRDC/"
    # ftp://igs-ftp.bkg.bund.de/IGS/BRDC/2024/082/BRDC00WRD_S_20240820000_01D_MN.rnx.gz

    ### generate regex
    rnx2rgx, rnx3rgx = _rnx_nav_rgx(
        date, site, sys="M", data_source="S"
    )  ### NAV RNX HERE !!!

    ### generate urls
    urldir = os.path.join(urlserver, str(date.year), conv.dt2doy(date))
    rnx3url = os.path.join(urldir, rnx3rgx)

    urldic = dict()
    urldic[3] = rnx3url

    return urldic

def renag_server_crtk(date, site=None, smp_01s=False):
    if smp_01s:
        urlserver = "ftp://renag.unice.fr/centipede_1s/"
    else:
        urlserver = "ftp://renag.unice.fr/centipede_30s/"
    urldic = _generic_server(date, site, urlserver)
    if len(site) != 4:
        urldic.pop(2, None)

    return urldic


def unavco_server(date, site=None):
    """
    UNAVCO server for RINEX observation files.

    Adapted from legacy function to standard format.
    """
    urlserver = "ftp://data-out.unavco.org/pub/rinex/obs/"
    urldic = _generic_server(date, site, urlserver)
    return urldic


def renag_server(date, site=None):
    """
    RENAG server for RINEX observation files.

    Adapted from legacy function to standard format.
    """
    urlserver = "ftp://renag.unice.fr/data/"
    urldic = _generic_server(date, site, urlserver)
    return urldic


def uwiseismic_server(date, site=None):
    """
    UWISeismic server for RINEX observation files.

    Adapted from legacy function to standard format.
    Note: Authentication may be required, use user/passwd parameters in download_gnss_rinex.
    """
    urlserver = "ftp://www2.uwiseismic.com/rinex/"
    urldic = _generic_server(date, site, urlserver)
    return urldic


def orpheon_server(date, site=None):
    """
    Orpheon server for RINEX observation files (RENAG hosted).

    Adapted from legacy function to standard format.
    Note: Authentication may be required, use user/passwd parameters in download_gnss_rinex.
    """
    urlserver = "ftp://renag.unice.fr/"
    urldic = _generic_server(date, site, urlserver)
    return urldic


def ovsg_server(date, site=None):
    """
    OVSG (Observatoire Volcanologique et Sismologique de Guadeloupe) server.

    Adapted from legacy function to standard format.
    Note: Uses HTTP, authentication may be required.
    """
    if dt.datetime(2009, 1, 1) <= date <= dt.datetime(2014, 2, 10):
        urlserver = "http://webobs.ovsg.univ-ag.fr/rawdata/GPS-GPSDATA.backtemp_20140210/"
    else:
        urlserver = "http://webobs.ovsg.univ-ag.fr/rawdata/GPS/GPSDATA/"

    ### generate regex
    rnx2rgx, rnx3rgx = _rnx_obs_rgx(date, site)

    ### generate urls (with rinex subdirectory specific to OVSG)
    urldir = os.path.join(urlserver, str(date.year), conv.dt2doy(date), "rinex")
    rnx2url = os.path.join(urldir, rnx2rgx)
    rnx3url = os.path.join(urldir, rnx3rgx)

    ### generate output urldic
    urldic = dict()
    urldic[2] = rnx2url
    urldic[3] = rnx3url

    return urldic


def geoaus_server(date, site=None):
    """
    Geosciences Australia GNSS server.

    Adapted from legacy function to standard format.
    Example: ftp://ftp.ga.gov.au/geodesy-outgoing/gnss/data/daily/2010/10063/
    """
    urlserver = "ftp://ftp.ga.gov.au/geodesy-outgoing/gnss/data/daily/"

    ### generate regex
    rnx2rgx, rnx3rgx = _rnx_obs_rgx(date, site)

    ### generate urls (with specific year+doy format)
    urldir = os.path.join(urlserver, str(date.year), date.strftime("%y") + conv.dt2doy(date))
    rnx2url = os.path.join(urldir, rnx2rgx)
    rnx3url = os.path.join(urldir, rnx3rgx)

    ### generate output urldic
    urldic = dict()
    urldic[2] = rnx2url
    urldic[3] = rnx3url

    return urldic


def ens_fr_server(date, site=None):
    """
    ENS (École Normale Supérieure) France GNSS server.

    Adapted from legacy function to standard format.
    """
    urlserver = "ftp://gnss.ens.fr/pub/public/crl/GPS/rinex/"
    urldic = _generic_server(date, site, urlserver)
    return urldic


def rgp_ign_smn_01s_server(date, site=None):
    """
    RGP IGN 01S/1Hz data server (hourly RINEX files).

    Adapted from legacy function to standard format.
    Returns URLs for all 24 hourly sessions.

    Note
    ----
    This function returns a list of urldic instead of a single urldic,
    one for each hourly session of the day.
    """
    urlserver = "ftp://rgpdata.ign.fr/pub/data/"

    urls_list = []

    for h in range(24):
        date_session = date.replace(hour=h)

        log.info("%s session %s", date_session, h)

        ### generate regex
        rnx2rgx, rnx3rgx = _rnx_obs_rgx(date_session, site)

        ### For 1Hz data, use session_a format for RINEX2
        if site:
            rnx2name = conv.statname_dt2rinexname(
                site.lower(), date_session, session_a_instead_of_daily_session=True
            )
            rnx2rgx = rnx2name  # Override with specific session name
        # else keep the wildcard pattern

        ### generate urls
        urldir = os.path.join(urlserver, str(date.year), conv.dt2doy(date), "data_1")
        rnx2url = os.path.join(urldir, rnx2rgx)
        rnx3url = os.path.join(urldir, rnx3rgx)

        ### generate output urldic
        urldic = dict()
        urldic[2] = rnx2url
        urldic[3] = rnx3url

        urls_list.append(urldic)

    return urls_list


def igs_cddis_nav_server(date, site=None):
    """
    IGS CDDIS navigation file server.

    Adapted from legacy function to standard format.
    """
    urlserver = "ftp://gdc.cddis.eosdis.nasa.gov/gps/data/daily/"

    ### generate regex
    rnx2rgx, rnx3rgx = _rnx_nav_rgx(date, site)

    ### generate urls (with specific subdirectory structure)
    urldir = os.path.join(
        urlserver, str(date.year), conv.dt2doy(date), date.strftime("%y") + "n"
    )
    rnx2url = os.path.join(urldir, rnx2rgx)
    rnx3url = os.path.join(urldir, rnx3rgx)

    ### generate output urldic
    urldic = dict()
    urldic[2] = rnx2url
    urldic[3] = rnx3url

    return urldic


# Remove legacy function definitions
# def igs_cddis_nav_server_legacy(stat, date):
# def rgp_ign_smn_1_hz_server_legacy(stat, date):
# def unavco_server_legacy(stat, date):
# def renag_server_legacy(stat, date):
# def uwiseismic_server_legacy(stat, date, user="", passwd=""):
# def orpheon_server_legacy(stat, date, user="", passwd=""):
# def ovsg_server_legacy(stat, date, user="", passwd=""):
# def geoaus_server_legacy(stat, date):
# def ens_fr_legacy(stat, date):

def _server_select(datacenter, date, site=None):
    mode1hz = False
    protocol = "ftp"
    urldic = dict()
    if datacenter in ("igs_cddis", "igs"):
        urldic = igs_cddis_server(date, site)
        protocol = "sftp"
    elif datacenter == "igs_sopac":
        urldic = igs_sopac_server(date, site)
    elif datacenter == "igs_ign":
        urldic = igs_ign_server(date, site)
    elif datacenter == "igs_ign_ensg":
        urldic = igs_ign_ensg_server(date, site)
    elif datacenter == "sonel":
        urldic = sonel_server(date, site)
    elif datacenter == "euref":
        urldic = euref_server(date, site)
    elif datacenter == "igs_bkg":
        urldic = igs_bkg_server(date, site)
    elif datacenter in ("nav", "brdc"):
        urldic = nav_rob_server(date, site)
    elif datacenter in ("nav_rt", "brdc_rt"):
        urldic = nav_bkg_server(date, site)
    elif datacenter == "nav_cddis":
        urldic = igs_cddis_nav_server(date, site)
    elif datacenter == "rgp":
        urldic = rgp_server(date, site)
    elif datacenter == "rgp_ensg":
        urldic = rgp_ensg_server(date, site)
    elif datacenter == "spotgins_eost":
        urldic = spotgins_eost_server(date, site)
    elif datacenter == "renag_crtk":
        urldic = renag_server_crtk(date, site)
    elif datacenter == "renag_crtk_01s":
        urldic = renag_server_crtk(date, site, smp_01s=True)
        mode1hz = True
    elif datacenter == "rgp_01s":
        urldic = rgp_ign_smn_01s_server(date, site)
        mode1hz = True
    elif datacenter == "renag":
        urldic = renag_server(date, site)
    elif datacenter == "orpheon":
        urldic = orpheon_server(date, site)
    elif datacenter == "uwiseismic":
        urldic = uwiseismic_server(date, site)
    elif datacenter == "ovsg":
        urldic = ovsg_server(date, site)
    elif datacenter == "unavco":
        urldic = unavco_server(date, site)
    elif datacenter == "geoaus":
        urldic = geoaus_server(date, site)
    elif datacenter == "ens_fr":
        urldic = ens_fr_server(date, site)
    else:
        log.warning("unkwn server dic in the dico, skip ...")
        return None, None, None

    return urldic, protocol, mode1hz

def effective_save_dir(parent_archive_dir, site, date, archtype="stat"):
    """
    INTERNAL_FUNCTION

    archtype =
        site
        site/year
        site/year/doy
        year/doy
        year/site
        week/dow
        OR only '/' for a dirty saving in the parent folder
        ... etc ...

    If site is None or empty, it will be replaced with "ALL_STATIONS" in the path.
    """
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

    # If site is None or empty, replace with a placeholder
    if not site:
        site = "ALL_STATIONS"

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


def rnx_regex_indir_all(rnx_regex, dir_files_list):
    """
    Match all files in a directory against a given regex pattern.

    Parameters
    ----------
    rnx_regex : str
        Regex pattern to match filenames.
    dir_files_list : list of str
        List of filenames in the directory.

    Returns
    -------
    list of str
        List of all matching filenames, or empty list if no match is found.
    """
    matches = [file for file in dir_files_list if re.search(rnx_regex, file)]
    return matches


def _check_local_file_exists(rnxrgx, outdir, local_files_cache, force=False, all_files_mode=False):
    """
    Check if a RINEX file exists locally.

    Parameters
    ----------
    rnxrgx : str
        Regex pattern to match RINEX filenames.
    outdir : str
        Local output directory path.
    local_files_cache : list of str
        Cached list of local files in the directory.
    force : bool, optional
        Force re-download even if file exists. Default is False.
    all_files_mode : bool, optional
        If True, check for all matching files. Default is False.

    Returns
    -------
    tuple of (bool, bool, str)
        - ok_loc : bool
            True if file exists locally and should not be re-downloaded
        - ok_dwl : bool
            True if file should be downloaded (not local or force=True)
        - rnxnam : str
            Filename if found locally (semicolon-separated if all_files_mode=True),
            empty string otherwise
    """
    if all_files_mode:
        # Find all matching local files
        rnxlocal_list = rnx_regex_indir_all(rnxrgx, local_files_cache)

        if rnxlocal_list:
            # Filter out empty or zero-size files
            valid_files = [f for f in rnxlocal_list if os.path.getsize(f) > 0]

            if valid_files:
                rnxloc_bns = [os.path.basename(f) for f in valid_files]
                rnxnam = ";".join(rnxloc_bns)

                if not force:
                    log.info("%d file(s) already exist locally ;)", len(valid_files))
                    return True, False, rnxnam
                else:
                    log.info(
                        "%d file(s) already exist locally, but re-download forced",
                        len(valid_files)
                    )
                    return False, True, rnxnam
    else:
        # Single file mode - original behavior
        rnxlocal = rnx_regex_indir(rnxrgx, local_files_cache)

        if rnxlocal and os.path.getsize(rnxlocal) > 0:
            rnxloc_bn = os.path.basename(rnxlocal)
            if not force:
                log.info("%s already exists locally ;)", rnxloc_bn)
                return True, False, rnxloc_bn
            else:
                log.info("%s already exists locally, but re-download forced", rnxloc_bn)
                return False, True, rnxloc_bn

    return False, False, ""


def _get_ftp_connection(ftpobj, host, protocol, sftp, user, passwd, prev_host, count_loop, count_nmax):
    """
    Get or create an FTP connection.

    Creates a new connection if the host has changed, counter exceeds max,
    or this is the first connection. Otherwise returns the existing connection.

    Parameters
    ----------
    ftpobj : FTP or None
        Current FTP connection object.
    host : str
        Target FTP server hostname.
    protocol : str
        Protocol string ("ftp" or "sftp").
    sftp : str or bool
        SFTP mode setting ('auto', True, or False).
    user : str or None
        FTP username.
    passwd : str or None
        FTP password.
    prev_host : str
        Previous host name.
    count_loop : int
        Current operation counter.
    count_nmax : int
        Maximum operations before reconnection.

    Returns
    -------
    tuple of (FTP, str)
        - ftpobj : FTP connection object
        - new_host : str
            The host name for connection tracking
    """
    if host != prev_host or count_loop > count_nmax or count_loop == 1:
        if ftpobj:
            ftpobj.close()

        # Determine SFTP mode: use protocol value if 'auto', otherwise use parameter
        if sftp == "auto":
            sftp_use = (protocol == "sftp")
        else:
            sftp_use = bool(sftp)

        ftpobj, _ = dlutils.ftp_objt_create(
            secure_ftp_inp=sftp_use,
            host=host,
            user=user,
            passwd=passwd,
        )
        return ftpobj, host

    return ftpobj, prev_host


def _get_ftp_directory_listing(ftpobj, directory, host, prev_dir):
    """
    Get FTP directory file listing.

    Changes to the specified directory and retrieves the list of files.
    Returns empty list if directory change fails.

    Parameters
    ----------
    ftpobj : FTP
        FTP connection object.
    directory : str
        Remote directory path.
    host : str
        FTP server hostname (for logging/URL generation).
    prev_dir : str
        Previous directory path.

    Returns
    -------
    tuple of (list, list)
        - ftp_files_list : list of str
            List of filenames in the directory
        - ftp_files_urls : list of str
            List of complete FTP URLs for all files
    """
    if prev_dir != directory:
        log.info("chdir " + directory)
        ftpobj.cwd("/")  # Reset to root directory first

        try:
            ftpobj.cwd(directory)
            ftp_files_list = dlutils.ftp_dir_list_files(ftpobj)

            # Generate full URLs for all files using string formatting
            ftp_files_urls = [
                f"ftp://{host}{directory}/{f}" for f in ftp_files_list
            ]
            return ftp_files_list, ftp_files_urls

        except Exception as e:
            log.warning("unable to chdir to %s, exception %s", directory, e)
            return [], []

    return None, None  # Directory unchanged, use cached values


def _match_files_in_directory(rnxrgx, ftp_files_list, all_files_mode=False):
    """
    Match files in FTP directory using regex pattern.

    Parameters
    ----------
    rnxrgx : str
        Regex pattern to match filenames.
    ftp_files_list : list of str
        List of filenames in the FTP directory.
    all_files_mode : bool, optional
        If True, return all matching files. If False, return only first match.
        Default is False.

    Returns
    -------
    str or list of str or None
        - If all_files_mode is False: first matching filename or None
        - If all_files_mode is True: list of all matching filenames (may be empty)
    """
    if all_files_mode:
        return rnx_regex_indir_all(rnxrgx, ftp_files_list)
    else:
        return rnx_regex_indir(rnxrgx, ftp_files_list)


def _update_table_row_with_match(table, irow, rnx_match, all_files_mode=False):
    """
    Update table row based on file match results.

    Parameters
    ----------
    table : pd.DataFrame
        Table to update.
    irow : int
        Row index to update.
    rnx_match : str or list of str or None
        Matched filename(s) from FTP directory.
    all_files_mode : bool, optional
        If True, rnx_match is a list of files. Default is False.

    Returns
    -------
    None
        Modifies table in place.
    """
    if rnx_match:
        table.loc[irow, "ok_dwl"] = True
        if all_files_mode:
            # In all_files_mode, rnx_match is a list
            # Join all filenames with semicolon separator
            table.loc[irow, "rnxnam"] = ";".join(rnx_match)
            log.info(f"{len(rnx_match)} file(s) found on server :)")
        else:
            # In single file mode, rnx_match is a string or None
            table.loc[irow, "rnxnam"] = rnx_match
            log.info(rnx_match + " found on server :)")
    else:
        table.loc[irow, "ok_dwl"] = False
        table.loc[irow, "rnxnam"] = ""
        log.warning(f"{table.loc[irow, "rnxrgx"]} not found on server :(")

    table.loc[irow, "crawled"] = True


def _generate_download_urls(table, all_files_mode=False):
    """
    Generate complete FTP URLs for files to be downloaded.

    Creates URL strings for files that are available remotely but not locally.

    Parameters
    ----------
    table : pd.DataFrame
        Table with crawl results.
    all_files_mode : bool, optional
        If True, handle semicolon-separated filenames. Default is False.

    Returns
    -------
    None
        Modifies table in place by adding 'url_true' column.
    """
    rnx_ok_dwl = table["ok_dwl"] & ~table["ok_loc"]

    if all_files_mode:
        # Handle multiple files separated by semicolons
        def make_urls(row):
            filenames = row["rnxnam"].split(";")
            # Use string formatting instead of os.path.join for FTP URLs
            urls = [f"ftp://{row['host']}/{row['dir']}/{fname}"
                    for fname in filenames]
            return ";".join(urls)

        table.loc[rnx_ok_dwl, "url_true"] = table.loc[rnx_ok_dwl].apply(
            make_urls, axis=1
        )
    else:
        # Single file mode - use string formatting for FTP URLs
        table.loc[rnx_ok_dwl, "url_true"] = table.loc[rnx_ok_dwl].apply(
            lambda x: f"ftp://{x['host']}/{x['dir']}/{x['rnxnam']}", axis=1
        )


def _collect_local_files(table):
    """
    Collect local file paths for files that exist locally.

    Parameters
    ----------
    table : pd.DataFrame
        Table with crawl results containing 'ok_loc', 'outdir', and 'rnxnam' columns.

    Returns
    -------
    pd.Series
        Series of local file paths, or empty Series if no local files.
    """
    rnx_ok_loc = table["ok_loc"]
    if rnx_ok_loc.sum() > 0:
        return pd.Series(
            table.loc[rnx_ok_loc].apply(
                lambda e: os.path.join(e["outdir"], e["rnxnam"]), axis=1
            )
        )
    else:
        return pd.Series([])


def crawl_ftp_files(
    table,
    sftp="auto",
    user=None,
    passwd=None,
    path_ftp_crawled_files_save=None,
    path_all_ftp_files_save=None,
    force=False,
    all_files_mode=False,
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
        - 'protocol': String indicating protocol ("ftp" or "sftp")
        - 'crawled': Boolean indicating if already crawled
    sftp : str or bool, optional
        SFTP mode setting. Default is 'auto'.
        - 'auto': Use the 'protocol' column value from each table row
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
    all_files_mode : bool, optional
        If True, download ALL files matching the regex pattern in each directory.
        If False (default), download only the first matching file (legacy behavior).
        When True, the 'rnxnam' column will contain semicolon-separated filenames.
        Default is False.

    Returns
    -------
    tuple of (pd.DataFrame, pd.Series, pd.Series)
        - table_use : pd.DataFrame
            Updated table with crawl results, including new columns:
            - 'ok_dwl': Boolean indicating file is available for download
            - 'ok_loc': Boolean indicating file exists locally
            - 'rnxnam': Actual filename(s) found on server
                       (semicolon-separated if all_files_mode=True)
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
    >>> # Single file mode (default)
    >>> table = pd.DataFrame({
    ...     'host': ['ftp.example.com'],
    ...     'dir': ['/data/2020/001'],
    ...     'outdir': ['/local/data'],
    ...     'rnxrgx': ['station001a.20o.*'],
    ...     'protocol': ['ftp'],
    ...     'crawled': [False]
    ... })
    >>> crawled_table, all_files, local_files = crawl_ftp_files(table)
    >>>
    >>> # All files mode - download all matching files in directory
    >>> crawled_table, all_files, local_files = crawl_ftp_files(
    ...     table, all_files_mode=True
    ... )
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
    prev_outdir = ""  # Track previous local output directory
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

        # Check local files when output directory changes
        if prev_outdir != row["outdir"]:
            local_files_lis = glob.glob(row["outdir"] + "/*")
            prev_outdir = row["outdir"]

        # Check if file already exists locally
        ok_loc, ok_dwl, rnxnam = _check_local_file_exists(
            row["rnxrgx"], row["outdir"], local_files_lis, force, all_files_mode
        )

        if ok_loc:
            table_use.loc[irow, "ok_loc"] = True
            table_use.loc[irow, "ok_dwl"] = False
            table_use.loc[irow, "rnxnam"] = rnxnam
            continue
        elif ok_dwl and force:
            table_use.loc[irow, "ok_loc"] = False
            table_use.loc[irow, "ok_dwl"] = True
            table_use.loc[irow, "rnxnam"] = rnxnam

        count_loop += 1

        # Get or create FTP connection
        ftpobj, prev_host = _get_ftp_connection(
            ftpobj, row["host"], row["protocol"], sftp, user, passwd,
            prev_host, count_loop, count_nmax
        )

        # Save intermediate results and reset counter on reconnection
        if count_loop > count_nmax:
            count_loop = 0
            _save_crawled_files(table_use)
            _get_and_save_all_ftp_files(all_ftp_files_stk)

        # Get file list when directory changes
        ftp_result = _get_ftp_directory_listing(ftpobj, row["dir"], row["host"], prev_dir)
        if ftp_result[0] is not None:  # Directory changed
            ftp_files_list, ftp_files_urls = ftp_result
            prev_dir = row["dir"]

            # Save all FTP files for reporting
            if ftp_files_urls:
                all_ftp_files_stk.append(pd.Series(ftp_files_urls))

        # Match files on server using regex pattern
        rnx_match = _match_files_in_directory(row["rnxrgx"], ftp_files_list, all_files_mode)

        # Update table based on file availability
        _update_table_row_with_match(table_use, irow, rnx_match, all_files_mode)

    # Generate URLs for downloadable files
    _generate_download_urls(table_use, all_files_mode)

    # Save final results
    _save_crawled_files(table_use)
    all_ftp_files = _get_and_save_all_ftp_files(all_ftp_files_stk)

    # Clean up FTP connection
    if ftpobj:
        ftpobj.close()

    # Collect local file paths
    all_loc_files = _collect_local_files(table_use)

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
    intuitive_output = False,
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
            * rgp_01s (IGN's RGP, all the 24 hourly rinex for the day will be downloaded)
            * renag
            * renag_crtk & renag_crtk_01s (Centipede RTK data hosted at RENAG)
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
        NB: files that already exist locally are marked
        as (local_path, True)

    Minimal exemple
    ---------------
    ```
    >>> statdic = dict()
    >>> statdic['igs_cddis'] = ['ZIMM','tlse']
    >>> archive_dir = '/home/USER/test_dl_rnx'
    >>> startdate = dt.datetime(2020,1,1)
    >>> enddate = dt.datetime(2020,1,31)
    >>> geodezyx.operational.download_gnss_rinex(statdic, output_dir, startdate, enddate)
    ```
    """

    if "01s" in list(statdico.keys())[0]:
        day_step = 0
        sec_step = 3600
    else:
        day_step = 1
        sec_step = 0

    date_range = conv.dt_range(startdate, enddate,day_step=day_step, sec_step=sec_step)

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
            secure_ftp=(table_dl["protocol"] == "sftp").values,
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

    if not intuitive_output:
        return out_tup_lis_fin
    else:
        return zip(*out_tup_lis_fin)


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

    all_sites = False
    for datacenter, site_lis in statdico.items():
        if site_lis is None:
            site_lis = [None]
            log.info("datacenter/stations: %s/all", datacenter)
        else:
            log.info("datacenter/stations: %s/%s", datacenter, " ".join(site_lis))
        for date, site in itertools.product(date_range, site_lis):
            if site is None:
                all_sites = True

            urldic, protocol, _ = _server_select(datacenter, date, site)
            if not urldic:
                continue

            outdir = effective_save_dir(output_dir, site, date, archtype)

            for rnxver, rnxurl in urldic.items():
                if (rnxver == 2 and no_rnx2) or (rnxver == 3 and no_rnx3):
                    continue
                table_proto.append((date, site, outdir, rnxver, rnxurl, protocol))

    # Create DataFrame with all collected data
    table = pd.DataFrame(
        table_proto, columns=["date", "site", "outdir", "ver", "url_theo", "protocol"]
    )

    # Add status columns
    table["crawled"] = False
    table["ok_dwl"] = False
    table["ok_loc"] = False
    table["all"] = all_sites
    table["url_true"] = None
    table["rnxnam"] = ""

    # Parse URL components
    urlpaths = table["url_theo"].apply(pathlib.Path)
    table["rnxrgx"] = urlpaths.apply(lambda p: p.name)
    table["host"] = urlpaths.apply(lambda p: p.parts[1])
    table["dir"] = urlpaths.apply(lambda p: os.path.join(*p.parts[2:-1]))

    return table

# date_range = conv.dt_range(dt.datetime(2020,1,1),
#                     dt.datetime(2020,12,31))
# t = gen_crawl_table({"igs_ign":None},
#                     date_range,
#                     "","/",
#                 False, False)
#
# c = crawl_ftp_files(t,all_files_mode=True)
