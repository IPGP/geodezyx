#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Server definitions and dispatcher for GNSS RINEX downloads.

Contains the individual ``*_server`` helper functions, the generic
server builder ``_generic_server``, regex helpers ``_rnx_obs_rgx`` /
``_rnx_nav_rgx``, and the central dispatcher ``_server_select``.
All used exclusively by :mod:`download_rinex`.
"""

import datetime as dt
import logging
import os

from geodezyx import conv

log = logging.getLogger("geodezyx")


def _rnx_obs_rgx(date, site="....", file_period="...", data_freq="..."):
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
    file_period : str, optional
        File period identifier for RINEX3 long filenames.
        Default is "01D" (daily).
    data_freq : str, optional
        Data frequency identifier for RINEX3 long filenames.
        Default is "..." (wildcard).

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

    rnx2rgx = conv.statname_dt2rinexname(site.lower(), date, rnxtype=".*")
    rnx3rgx = conv.statname_dt2rinexname_long(
        site,
        date,
        country="...",
        data_source=".",
        file_period=file_period,
        data_freq=data_freq,
        data_type=".O",
        format_compression=".*",
    )
    return str(rnx2rgx), str(rnx3rgx)


def _rnx_nav_rgx(date, site="....", sys=".", data_source="."):
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
    data_source : str, optional
        Data source identifier for RINEX3 long filenames. Default is "." (wildcard).

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
    - !!!! File period: 01D (daily) for navigation files !!!!
    - !!!! Data frequency: empty ("") for navigation files !!!!
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


def _generic_server(
    date, site=None, urlserver=None, urlsuffix=None, file_period="...", data_freq="..."
):
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
    file_period : str, optional
        File period identifier for RINEX3 long filenames. Default is "01D" (daily).
    data_freq : str, optional
        Data frequency identifier for RINEX3 long filenames. Default is "..." (wildcard).

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
    rnx2rgx, rnx3rgx = _rnx_obs_rgx(
        date, site, file_period=file_period, data_freq=data_freq
    )

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
    rnx2rgx, rnx3rgx = _rnx_nav_rgx(
        date, site, sys="M", data_source="R"
    )  ### NAV RNX HERE !!!

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


def rgp_rnx3_server(date, site=None):
    urlserver = "ftp://rgpdata.ign.fr/pub/data_v3/"
    urldic = _generic_server(date, site, urlserver, "data_30")
    return urldic


def rgp_ensg_server(date, site=None):
    urlserver = "ftp://rgpdata.ensg.eu/pub/data/"
    urldic = _generic_server(date, site, urlserver, "data_30")
    return urldic


def rgp_ensg_rnx3_server(date, site=None):
    urlserver = "ftp://rgpdata.ensg.eu/pub/data_v3/"
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
        freq = "01S"
    else:
        urlserver = "ftp://renag.unice.fr/centipede_30s/"
        freq = "30S"
    urldic = _generic_server(date, site, urlserver, data_freq=freq)
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
        urlserver = (
            "http://webobs.ovsg.univ-ag.fr/rawdata/GPS-GPSDATA.backtemp_20140210/"
        )
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
    urldir = os.path.join(
        urlserver, str(date.year), date.strftime("%y") + conv.dt2doy(date)
    )
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


def _server_select_rnx(datacenter, date, site=None):
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
