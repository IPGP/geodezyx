# -*- coding: utf-8 -*-
"""
Server definitions and dispatcher for GNSS RINEX downloads (refactored).

This refactored version uses configuration from a YAML file instead of
hardcoded server parameters. The configuration is loaded via
`servers_config_loader` module.

Contains:
- Regex helpers: `_rnx_obs_rgx` / `_rnx_nav_rgx`
- Generic server builder: `_generic_server`
- Server-specific builders for non-generic servers
- Central dispatcher: `_server_select`
"""

import datetime as dt
import logging
import os
from typing import Dict, Any, Optional, Tuple, Union, List

from geodezyx import conv
from .servers_config_loader import get_config, ServerConfigError

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
    """
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

    Parameters
    ----------
    date : datetime.datetime
        Date for which to generate the regex patterns.
    site : str, optional
        4-character GNSS station name (e.g., 'BRDC', 'ZIMM').
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

    This function creates download URLs for both RINEX2 and RINEX3 observation
    files following a standard FTP server directory structure: server/year/doy/filename.

    Parameters
    ----------
    date : datetime.datetime
        Date for which to generate the URLs.
    site : str, optional
        4-character GNSS station name (e.g., 'ZIMM', 'TLSE').
    urlserver : str
        Base FTP server URL (e.g., 'ftp://example.com/data/').
    urlsuffix : str, optional
        Additional path suffix between date path and filename
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
    urldic = {2: rnx2url, 3: rnx3url}
    return urldic


def _server_from_config(datacenter: str, date, site=None) -> Tuple[Dict, str, bool]:
    """
    Build server URLs using configuration from YAML file.

    This function handles servers that use the generic URL builder pattern
    and can be fully configured via YAML.

    Parameters
    ----------
    datacenter : str
        Server name/alias
    date : datetime.datetime
        Date for which to generate URLs
    site : str, optional
        4-character station name

    Returns
    -------
    tuple
        (urldic, protocol, mode1hz) where:
        - urldic : dict - URL dictionary (keys 2 and/or 3 for RINEX version)
        - protocol : str - Protocol to use ('ftp', 'https', 'sftp', 'http')
        - mode1hz : bool - True if returns list of URLs (hourly mode)

    Raises
    ------
    ServerConfigError
        If server configuration is not found
    """
    config = get_config()
    
    try:
        server_cfg = config.get_server(datacenter)
    except ServerConfigError as e:
        log.warning("Unknown server '%s': %s", datacenter, e)
        return None, None, None

    protocol = server_cfg.get('protocol', 'ftp')
    urlserver = server_cfg.get('url')
    urlsuffix = server_cfg.get('urlsuffix', '')
    data_freq = server_cfg.get('data_freq', '...')
    mode1hz = server_cfg.get('mode') == 'hourly'

    if not urlserver:
        log.warning("No URL defined for server '%s'", datacenter)
        return None, None, None

    # Handle special cases
    if mode1hz:
        # Hourly mode - return list of urldic
        urls_list = []
        for h in range(24):
            date_session = date.replace(hour=h)
            rnx2rgx, rnx3rgx = _rnx_obs_rgx(date_session, site, data_freq=data_freq)

            if site:
                rnx2name = conv.statname_dt2rinexname(
                    site.lower(), date_session, session_a_instead_of_daily_session=True
                )
                rnx2rgx = rnx2name

            urldir = os.path.join(urlserver, str(date.year), conv.dt2doy(date), urlsuffix)
            rnx2url = os.path.join(urldir, rnx2rgx)
            rnx3url = os.path.join(urldir, rnx3rgx)

            urldic = {2: rnx2url, 3: rnx3url}
            if len(site) != 4:
                urldic.pop(2, None)
            urls_list.append(urldic)

        # Return the list directly - will be handled by dispatcher
        return urls_list, protocol, mode1hz

    # Standard generic server
    urldic = _generic_server(date, site, urlserver, urlsuffix, data_freq=data_freq)
    return urldic, protocol, mode1hz


# ============================================================================
# SPECIAL CASE SERVERS (non-generic or with custom logic)
# ============================================================================

def _igs_cddis_server(date, site=None):
    """IGS CDDIS server with custom directory structure."""
    urlserver = "ftp://gdc.cddis.eosdis.nasa.gov/gps/data/daily/"
    rnx2rgx, rnx3rgx = _rnx_obs_rgx(date, site)
    urldir = os.path.join(
        urlserver, str(date.year), conv.dt2doy(date), date.strftime("%y") + "d"
    )
    rnx2url = os.path.join(urldir, rnx2rgx)
    rnx3url = os.path.join(urldir, rnx3rgx)
    return {2: rnx2url, 3: rnx3url}


def _euref_server(date, site=None):
    """EUREF server with uppercase RINEX2 names."""
    urlserver = "ftp://epncb.oma.be/pub/obs/"
    rnx2rgx, rnx3rgx = _rnx_obs_rgx(date, site)
    rnx2rgx = rnx2rgx.upper()
    urldir = str(os.path.join(urlserver, str(date.year), conv.dt2doy(date)))
    rnx2url = os.path.join(urldir, rnx2rgx)
    rnx3url = os.path.join(urldir, rnx3rgx)
    return {2: rnx2url, 3: rnx3url}


def _nav_rob_server(date, site=None):
    """ROB navigation server (year-only structure)."""
    urlserver = "https://epncb.oma.be/pub/RINEX/BRDC/"
    rnx2rgx, rnx3rgx = _rnx_nav_rgx(date, site, sys="M", data_source="R")
    urldir = os.path.join(urlserver, str(date.year))
    rnx2url = os.path.join(urldir, rnx2rgx)
    rnx3url = os.path.join(urldir, rnx3rgx)
    return {2: rnx2url, 3: rnx3url}


def _nav_bkg_server(date, site=None):
    """BKG navigation server."""
    urlserver = "https://igs.bkg.bund.de/root_ftp/IGS/BRDC/"
    rnx2rgx, rnx3rgx = _rnx_nav_rgx(date, site, sys="M", data_source="S")
    urldir = os.path.join(urlserver, str(date.year), conv.dt2doy(date))
    rnx3url = os.path.join(urldir, rnx3rgx)
    return {3: rnx3url}


def _nav_cddis_server(date, site=None):
    """CDDIS navigation server."""
    urlserver = "ftp://gdc.cddis.eosdis.nasa.gov/gps/data/daily/"
    rnx2rgx, rnx3rgx = _rnx_nav_rgx(date, site)
    urldir = os.path.join(
        urlserver, str(date.year), conv.dt2doy(date), date.strftime("%y") + "n"
    )
    rnx2url = os.path.join(urldir, rnx2rgx)
    rnx3url = os.path.join(urldir, rnx3rgx)
    return {2: rnx2url, 3: rnx3url}


def _ovsg_server(date, site=None):
    """OVSG server with date-dependent URL."""
    if dt.datetime(2009, 1, 1) <= date <= dt.datetime(2014, 2, 10):
        urlserver = "http://webobs.ovsg.univ-ag.fr/rawdata/GPS-GPSDATA.backtemp_20140210/"
    else:
        urlserver = "http://webobs.ovsg.univ-ag.fr/rawdata/GPS/GPSDATA/"

    rnx2rgx, rnx3rgx = _rnx_obs_rgx(date, site)
    urldir = os.path.join(urlserver, str(date.year), conv.dt2doy(date), "rinex")
    rnx2url = os.path.join(urldir, rnx2rgx)
    rnx3url = os.path.join(urldir, rnx3rgx)
    return {2: rnx2url, 3: rnx3url}


def _geoaus_server(date, site=None):
    """Geosciences Australia server with yy+doy format."""
    urlserver = "ftp://ftp.ga.gov.au/geodesy-outgoing/gnss/data/daily/"
    rnx2rgx, rnx3rgx = _rnx_obs_rgx(date, site)
    urldir = os.path.join(
        urlserver, str(date.year), date.strftime("%y") + conv.dt2doy(date)
    )
    rnx2url = os.path.join(urldir, rnx2rgx)
    rnx3url = os.path.join(urldir, rnx3rgx)
    return {2: rnx2url, 3: rnx3url}


# Mapping of special servers to their builder functions
_SPECIAL_SERVERS = {
    'igs_cddis': (_igs_cddis_server, 'sftp'),
    'igs': (_igs_cddis_server, 'sftp'),  # Alias
    'euref': (_euref_server, 'ftp'),
    'nav_rob': (_nav_rob_server, 'https'),
    'nav': (_nav_rob_server, 'https'),  # Alias
    'brdc': (_nav_rob_server, 'https'),  # Alias
    'nav_bkg': (_nav_bkg_server, 'https'),
    'nav_rt': (_nav_bkg_server, 'https'),  # Alias
    'brdc_rt': (_nav_bkg_server, 'https'),  # Alias
    'nav_cddis': (_nav_cddis_server, 'ftp'),
    'ovsg': (_ovsg_server, 'http'),
    'geoaus': (_geoaus_server, 'ftp'),
}


def _server_select_rnx(datacenter, date, site=None):
    """
    Central dispatcher for server selection.

    This function routes server requests either to special-case handlers
    (for servers with custom logic) or to the generic configuration-based
    builder.

    Parameters
    ----------
    datacenter : str
        Server name or alias
    date : datetime.datetime
        Date for data request
    site : str, optional
        4-character station name

    Returns
    -------
    tuple
        (urldic, protocol, mode1hz) where:
        - urldic : dict or list - URL dictionary/dictionaries
        - protocol : str - Protocol to use
        - mode1hz : bool - Hourly mode flag
    """
    # Check if this is a special server
    if datacenter in _SPECIAL_SERVERS:
        server_func, protocol = _SPECIAL_SERVERS[datacenter]
        urldic = server_func(date, site)
        mode1hz = False
        return urldic, protocol, mode1hz

    # Try configuration-based builder
    urldic, protocol, mode1hz = _server_from_config(datacenter, date, site)

    if urldic is None:
        log.warning("Unknown server '%s', skip", datacenter)
        return None, None, None

    return urldic, protocol, mode1hz
