#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 23/09/2025 12:24:23

@author: psakic
"""

import hatanaka


import geodezyx.operational as opera
import geodezyx.conv as conv
import datetime as dt


def rinex_check(rnx_path_inp):
    """
    Check the validity of a RINEX file.

    Parameters
    ----------
    rnx_path_inp : str
        Path to the RINEX file.

    Returns
    -------
    bool
        True if the RINEX file is valid, False otherwise.
    """

    ####Open RINEX, decompress if needed
    rnx_use = opera.rinex_open(rnx_path_inp)

    #### RINEX's start/end
    first, last, itrvl = opera.rinex_start_end(rnx_use, interval_out=True)
    first_rnd = conv.round_dt(first, "1D")
    last_rnd = conv.round_dt(last, "1D")

    ######## BOOLEAN TESTS
    ## Check if the RINEX file is included in a full day
    bool_srt_end = True if last_rnd == first_rnd + dt.timedelta(days=1) else False

    ######## FINAL BOOLEAN
    bool_final = bool_srt_end

    return bool_final

p = "/home/psakicki/aaa_FOURBI/ZIM3/ZIM300CHE_R_20243660000_01D_30S_MO.crx.gz"
p = "/home/psakicki/aaa_FOURBI/ZIM3/ZIM300CHE_R_20243660000_01D_30S_MO.crx.gz"
p = "/home/psakicki/aaa_FOURBI/ZIM3/ZIM300CHE_R_20243660000_01D_30S_MO.rnx"

rinex_check(p)
