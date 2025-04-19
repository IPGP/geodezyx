#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 15:11:12 2020

@author: psakicki

This sub-module of geodezyx.conv deals with rinex name handeling conversion.

it can be imported directly with:
from geodezyx import conv

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""


#### Import the logger
import logging

########## BEGIN IMPORT ##########
#### External modules
import os
import re

### Imported in the corresponding function to avoid cyclic import
### https://stackoverflow.com/questions/1250103/attributeerror-module-object-has-no-attribute
log = logging.getLogger('geodezyx')


def rinex_regex_search_tester(
    str_in,
    short_name=True,
    long_name=True,
    brdc_long_name=False,
    gfz_godc_name=False,
    compressed=None,
):
    """
    Frontend function for a RINEX regex search/match

    Parameters
    ----------
    str_in : str
        input string.
        If it is a path, will extract its basename.
    short_name : bool, optional
        check if the pattern matches a short name RINEX. The default is True.
    long_name : bool, optional
        check if the pattern matches a long name RINEX. The default is True.
    brdc_long_name : bool, optional
        check if the pattern matches a brdc long name RINEX. The default is False.
    gfz_godc_name: bool, optional
        check if the pattern matches a GFZ's GODC (GNSS Operational Data Center)
        internal long name RINEX. The default is False.
    compressed : bool or None
        return the regex for a compressed rinex
        if None, does not matter (return both compressed or not)

    Returns
    -------
    None or re object
        match result.

    """

    str_in = os.path.basename(str_in)

    match = None

    search_out = re.search(rinex_regex(compressed), str_in)
    if short_name and search_out:
        match = search_out

    search_out = re.search(rinex_regex_long_name(compressed), str_in)
    if long_name and search_out:
        match = search_out

    search_out = re.search(rinex_regex_long_name_brdc(compressed), str_in)
    if brdc_long_name and search_out:
        match = search_out

    search_out = re.search(rinex_regex_long_name_gfz_godc(compressed), str_in)
    if gfz_godc_name and search_out:
        match = search_out

    return match


def rinex_regex(compressed=None, compiled=False):
    """
    Return a regex corresponding to a RINEX name (v2 - short convention)

    Parameters
    ----------
    compressed : bool or None
        return the regex for a compressed rinex
        if None, does not matter (return both compressed or not)


    compiled : bool
        return a Python regex object already compliled

    Returns
    -------
    out : string or python's regex
        a regex
    """
    if compressed is None:
        regexstr = r"^....[0-9]{3}.\.[0-9]{2}((d\.(Z|z|gz))|o|d)$"
    elif not compressed:
        regexstr = r"^....[0-9]{3}.\.[0-9]{2}o$"
    else:
        regexstr = r"^....[0-9]{3}.\.[0-9]{2}((d\.(Z|z|gz))|d)$"

    if compiled:
        return re.compile(regexstr)

    return regexstr


def rinex_regex_long_name(compressed=None, compiled=False):
    """
    Return a regex corresponding to a RINEX name (v3/4 - long convention)

    Parameters
    ----------
    compressed : bool or None
        return the regex for a compressed rinex
        if None, does not matter (return both compressed or not)


    compiled : bool
        return a Python regex object already compliled

    Returns
    -------
    out : string or python's regex
        a regex
    """
    if compressed: # compressed only
        regexstr = r"^.{4}[0-9]{2}.{3}_(R|S|U)_[0-9]{11}_[0-9]{2}\w_[0-9]{2}\w_\w{2}\.\w{3}\.gz$"
    elif compressed is None: # compressed or not
        regexstr = r"^.{4}[0-9]{2}.{3}_(R|S|U)_[0-9]{11}_[0-9]{2}\w_[0-9]{2}\w_\w{2}\.\w{3}(\.gz)?$"
    else: # not compressed only
        regexstr = (
            r"^.{4}[0-9]{2}.{3}_(R|S|U)_[0-9]{11}_[0-9]{2}\w_[0-9]{2}\w_\w{2}\.\w{3}$"
        )

    if compiled:
        return re.compile(regexstr)

    return regexstr


def rinex_regex_long_name_brdc(compressed=None, compiled=False):
    """
    Return a regex corresponding to a BROADCAST RINEX name (v3/4 - new convention)

    Parameters
    ----------
    compressed : bool or None
        return the regex for a compressed rinex
        if None, does not matter (return both compressed or not)


    compiled : bool
        return a Python regex object already compliled

    Returns
    -------
    out : string or python's regex
        a regex
    """
    if compressed:
        regexstr = r"^.{4}[0-9]{2}.{3}_(R|S|U)_[0-9]{11}_[0-9]{2}\w_\w{2}\.\w{3}\.gz$"
    elif compressed is None:
        regexstr = (
            r"^.{4}[0-9]{2}.{3}_(R|S|U)_[0-9]{11}_[0-9]{2}\w_\w{2}\.\w{3}(\.gz)?$"
        )
    else:
        regexstr = r"^.{4}[0-9]{2}.{3}_(R|S|U)_[0-9]{11}_[0-9]{2}\w_\w{2}\.\w{3}$"
    if compiled:
        return re.compile(regexstr)

    return regexstr


def rinex_regex_long_name_gfz_godc(compressed=None, compiled=False):
    """
    Return a regex corresponding to a RINEX name (new convention)
    from the GFZ's GODC archive (both OBS and NAV types)

    Parameters
    ----------
    compressed : bool or None
        return the regex for a compressed rinex
        if None, does not matter (return both compressed or not)


    compiled : bool
        return a Python regex object already compliled

    Returns
    -------
    out : string or python's regex
        a regex
    """
    ### Some exemples for tests
    # BRDC00GFZ_00000000_FRO_RX3_MN_20210114_000000_01D_00U_GFZ.rnx
    # AIRA00JPN_00001047_FRO_RX3_MO_20220105_000000_01D_30S_IGS.crx.gz
    # BRDC00GFZ_00000000_FRO_RX3_MN_20220111_000000_01D_00U_GFZ.rnx
    # YEL200CAN_00002760_FRO_RX3_MO_20220105_000000_01D_30S_IGS.crx.gz

    if compressed:
        regexstr = r"^.{4}[0-9]{2}.{3}_[0-9]{8}_.{3}_.{3}_.{2}_[0-9]{8}_[0-9]{6}_[0-9]{2}\w_[0-9]{2}\w_.{3}\.\w{3}\.gz$"
    elif compressed is None:
        regexstr = r"^.{4}[0-9]{2}.{3}_[0-9]{8}_.{3}_.{3}_.{2}_[0-9]{8}_[0-9]{6}_[0-9]{2}\w_[0-9]{2}\w_.{3}\.\w{3}(\.gz)?$"
    else:
        regexstr = r"^.{4}[0-9]{2}.{3}_[0-9]{8}_.{3}_.{3}_.{2}_[0-9]{8}_[0-9]{6}_[0-9]{2}\w_[0-9]{2}\w_.{3}\.\w{3}$"
    # Nouvelle version : peut  normalement digerer tout RINEX du GFZ
    # le suffixe GFZ peut etre nimporte quelle chaine  de upper char !
    # egalement implemene dans RINEXMOD
    # .{4}[0-9]{2}.{3}_[0-9]{8}_.{3}_.{3}_.{2}_[0-9]{8}_[0-9]{6}_[0-9]{2}\w_[0-9]{2}\w_[A-Z]*\.\w{3}(\.gz)?

    if compiled:
        return re.compile(regexstr)

    return regexstr


def rinex_regex_tester():
    """
    Returns two list of string to test the efficiency of the RINEX regexs.

    true_positive_rnxs is a list of strings describing actual RINEX filenames.
    They must be matched

    false_positive_rnxs is a list of strings describing other filenames.
    But which have been matched at one moment because of too weak regexs
    They must NOT be matched
    """

    true_positive_rnxs = [
        "psa13400.21o",
        "gps11000.21d.Z",
        "tar11320.00d.Z",
        "gosi0010.14d.Z",
        "YARR00AUS_R_20190100000_01D_30S_MO.rnx.gz",
        "REYK00ISL_R_20200420000_01D_30S_MO.crx.gz",
        "CFNG00XXX_R_20240241500_01H_01S_MO.crx.gz",
        "BRDC00GOP_R_20190450000_01D_MN.rnx.gz",
        "BRDC00GFZ_00000000_FRO_RX3_MN_20210114_000000_01D_00U_GFZ.rnx",
        "AIRA00JPN_00001047_FRO_RX3_MO_20220105_000000_01D_30S_IGS.crx.gz",
        "BRDC00GFZ_00000000_FRO_RX3_MN_20220111_000000_01D_00U_GFZ.rnx",
        "YEL200CAN_00002760_FRO_RX3_MO_20220105_000000_01D_30S_IGS.crx.gz",
    ]

    false_positive_rnxs = [
        "LROC00FRA_00001225_FRO_RX3_MO_20220613_000000_01D_30S_IGS.rnx.json",
        "SUN600SWE_00002886_FRO_RX3_MO_20220613_000000_01D_30S_EPN.rnx.json",
        "_._YELL00CAN_00001441_FRO_RX3_MO_20220613_000000_01D_30S_IGS.rnx",
    ]

    return true_positive_rnxs, false_positive_rnxs


def interval_from_rinex_name(rnx_name_inp):
    """
    from a RINEX file name (long name) determines the nominal data interval

    Parameters
    ----------
    rnx_name_inp : str
        long RINEX filename.

    Returns
    -------
    interval_ok : int
        nominal data interval in seconds.
    """

    rnx_name = os.path.basename(rnx_name_inp)

    if not rinex_regex_search_tester(rnx_name, False, True):
        log.error(
            "the input RINEX %s does not have a long name, unable to detect interval",
            rnx_name,
        )
        raise Exception

    interval_str = rnx_name[28:31]

    if interval_str[-1] == "D":
        coef = 86400
    elif interval_str[-1] == "H":
        coef = 3600
    elif interval_str[-1] == "M":
        coef = 60
    else:
        coef = 1

    interval_ok = int(interval_str[:2]) * coef

    return interval_ok


def period_from_rinex_name(rnx_name_inp):
    """
    from a RINEX file name (long naming) determines the nominal file period

    Parameters
    ----------
    rnx_name_inp : str
        long RINEX filename.

    Returns
    -------
    interval_ok : int
        nominal data interval in seconds.
    """

    rnx_name = os.path.basename(rnx_name_inp)

    if not rinex_regex_search_tester(rnx_name, False, True):
        log.error(
            "the input RINEX %s does not have a long name, unable to detect period",
            rnx_name,
        )
        raise Exception

    period_str = rnx_name[24:27]

    if period_str[-1] == "D":
        coef = 86400
    elif period_str[-1] == "H":
        coef = 3600
    elif period_str[-1] == "M":
        coef = 60
    else:
        coef = 1

    period_ok = int(period_str[:2]) * coef

    return period_ok
