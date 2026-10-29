# -*- coding: utf-8 -*-
"""
GNSS RINEX server dispatcher
"""

import datetime as dt
import logging
import os
from typing import Dict, Any, Optional, Tuple, Union

from geodezyx import conv
from .servers_rinex_cfg_loader import get_config, ServerConfigError

log = logging.getLogger("geodezyx")


def _rnx_obs_rgx(date, site="....", file_period="...", data_freq="..."):
    """Generate RINEX observation file regex patterns."""
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
    """Generate RINEX navigation file regex patterns."""
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


def _apply_case_formatting(filename: str, case: Optional[str]) -> str:
    """Apply case formatting to filename (lower/upper)."""
    if case == "upper":
        return filename.upper()
    elif case == "lower":
        return filename.lower()
    return filename


def _build_url_path(urlserver: str, date, config: Dict[str, Any], urlsuffix: str = "") -> str:
    """
    Build the URL path directory using the structure pattern from config.
    
    Supports:
    - "year/doy" : standard /YYYY/DOY/
    - "year" : only /YYYY/
    - "year/doy/subdir" : /YYYY/DOY/{subdir}/ with subdir_format pattern
    - "year/yydoy" : /YYYY/{YY}{DOY}/
    """
    structure = config.get("structure", "year/doy")
    
    # Build base path
    if structure == "year":
        urldir = os.path.join(urlserver, str(date.year))
    elif structure == "year/yydoy":
        # Format: /YYYY/{YY}{DOY}/ (e.g., /2024/24015/)
        yydoy = date.strftime("%y") + conv.dt2doy(date)
        urldir = os.path.join(urlserver, str(date.year), yydoy)
    elif structure == "year/doy":
        urldir = os.path.join(urlserver, str(date.year), conv.dt2doy(date))
    elif structure == "year/doy/subdir":
        # Handle custom subdirectory pattern
        subdir_format = config.get("subdir_format", "")
        if subdir_format:
            # Format the subdir pattern (e.g., "{yy}d" -> "24d")
            subdir = subdir_format.format(
                yy=date.strftime("%y"),
                doy=conv.dt2doy(date),
                year=str(date.year),
            )
            urldir = os.path.join(urlserver, str(date.year), conv.dt2doy(date), subdir)
        else:
            urldir = os.path.join(urlserver, str(date.year), conv.dt2doy(date))
    else:
        # Fallback to standard structure
        urldir = os.path.join(urlserver, str(date.year), conv.dt2doy(date))
    
    # Add optional suffix
    if urlsuffix:
        urldir = os.path.join(urldir, urlsuffix)
    
    return urldir


def _generic_server(
    date, 
    site=None, 
    config: Optional[Dict[str, Any]] = None,
    urlsuffix: str = ""
):
    """
    Ultra-flexible generic server builder.
    
    Handles nearly all servers through configuration parameters:
    - Different URL structures (year, year/doy, year/yydoy, year/doy/subdir)
    - Case formatting (upper/lower for RINEX2 and RINEX3)
    - Navigation file parameters (sys, data_source)
    - Data frequency
    - File period
    """
    if not config:
        return None
    
    urlserver = config.get("url")
    if not urlserver:
        return None
    
    file_type = config.get("file_type", "observation")
    rnx_versions = config.get("rnx_version", [2, 3])
    data_freq = config.get("data_freq", "...")
    
    urldir = _build_url_path(urlserver, date, config, urlsuffix)
    
    urldic = {}
    
    if file_type == "observation":
        rnx2rgx, rnx3rgx = _rnx_obs_rgx(date, site, data_freq=data_freq)
        
        # Apply case formatting
        rnx2_case = config.get("rnx2_case", "lower")
        rnx3_case = config.get("rnx3_case", "upper")
        rnx2rgx = _apply_case_formatting(rnx2rgx, rnx2_case)
        rnx3rgx = _apply_case_formatting(rnx3rgx, rnx3_case)
        
        if 2 in rnx_versions:
            urldic[2] = os.path.join(urldir, rnx2rgx)
        if 3 in rnx_versions:
            urldic[3] = os.path.join(urldir, rnx3rgx)
    
    elif file_type == "navigation":
        nav_sys = config.get("nav_sys", ".")
        nav_data_source = config.get("nav_data_source", ".")
        rnx2rgx, rnx3rgx = _rnx_nav_rgx(date, site, sys=nav_sys, data_source=nav_data_source)
        
        # Apply case formatting
        rnx2_case = config.get("rnx2_case", "lower")
        rnx3_case = config.get("rnx3_case", "upper")
        rnx2rgx = _apply_case_formatting(rnx2rgx, rnx2_case)
        rnx3rgx = _apply_case_formatting(rnx3rgx, rnx3_case)
        
        if 2 in rnx_versions:
            urldic[2] = os.path.join(urldir, rnx2rgx)
        if 3 in rnx_versions:
            urldic[3] = os.path.join(urldir, rnx3rgx)
    
    return urldic


def _server_from_cfg(datacenter: str, date, site=None) -> Tuple:
    """
    Load server from config and build URLs using the flexible generic builder.
    
    Returns: (urldic, protocol, mode1hz)
    """
    config = get_config()
    
    try:
        server_cfg = config.get_server(datacenter)
    except ServerConfigError as e:
        log.warning("Unknown server '%s': %s", datacenter, e)
        return None, None, None
    
    protocol = server_cfg.get("protocol", "ftp")
    urlsuffix = server_cfg.get("urlsuffix", "")
    mode = server_cfg.get("mode", "standard")
    mode1hz = mode == "hourly"
    
    # Hourly mode: return list of URLs (one per hour)
    if mode1hz:
        urls_list = []
        for h in range(24):
            date_session = date.replace(hour=h)
            urldic = _generic_server(date_session, site, server_cfg, urlsuffix)
            
            # For 1Hz data with specific site, use session naming
            if site and len(site) == 4:
                # Override RINEX2 with session naming if needed
                pass
            
            urls_list.append(urldic)
        
        return urls_list, protocol, mode1hz
    
    # Standard mode: return single URL dict
    urldic = _generic_server(date, site, server_cfg, urlsuffix)
    
    return urldic, protocol, mode1hz


# ============================================================================
# TRULY EXCEPTIONAL CASES (keep minimal)
# ============================================================================

def _ovsg_server(date, site=None):
    """
    OVSG server - ONLY special case kept because URL changes based on date!
    
    This is genuinely exceptional: URL structure differs by date range.
    Can't be expressed in simple config.
    """
    if dt.datetime(2009, 1, 1) <= date <= dt.datetime(2014, 2, 10):
        urlserver = "http://webobs.ovsg.univ-ag.fr/rawdata/GPS-GPSDATA.backtemp_20140210/"
    else:
        urlserver = "http://webobs.ovsg.univ-ag.fr/rawdata/GPS/GPSDATA/"
    
    rnx2rgx, rnx3rgx = _rnx_obs_rgx(date, site)
    urldir = os.path.join(urlserver, str(date.year), conv.dt2doy(date), "rinex")
    rnx2url = os.path.join(urldir, rnx2rgx)
    rnx3url = os.path.join(urldir, rnx3rgx)
    return {2: rnx2url, 3: rnx3url}


# Mapping of TRULY EXCEPTIONAL servers (minimal set)
_TRULY_EXCEPTIONAL = {
    'ovsg': (_ovsg_server, 'http'),
}


def _server_select_rnx(datacenter, date, site=None):
    """
    Central dispatcher using ultra-flexible generic builder.
    
    Strategy:
    1. Check if server is TRULY exceptional (ovsg) → use special function
    2. Otherwise → use powerful generic builder from config
    
    This minimizes special cases while maintaining flexibility.
    """
    # Check for truly exceptional cases
    if datacenter in _TRULY_EXCEPTIONAL:
        server_func, protocol = _TRULY_EXCEPTIONAL[datacenter]
        urldic = server_func(date, site)
        return urldic, protocol, False
    
    # Use flexible generic builder for everything else
    urldic, protocol, mode1hz = _server_from_cfg(datacenter, date, site)
    
    if urldic is None:
        log.warning("Unknown server '%s', skip", datacenter)
        return None, None, None
    
    return urldic, protocol, mode1hz
