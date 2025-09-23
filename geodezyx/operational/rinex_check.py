#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 23/09/2025 12:24:23

@author: psakic
"""


import geodezyx.operational as opera
import geodezyx.conv as conv

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

    #### start/end test
    first, last, itrvl = opera.rinex_start_end(rnx_path_inp, interval_out=True)
    first_rnd = conv.round_dt(first, "1D")
    last_rnd = conv.round_dt(last, "1D")

    if last_rnd == first_rnd:
        bool_srt_end = True
    else:
        bool_srt_end = False


    bool_final = bool_srt_end

    return bool_final














