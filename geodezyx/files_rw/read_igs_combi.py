#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 11:42:14 2021

@author: psakic
"""

#### Import the logger
import logging
import os
import re

########## BEGIN IMPORT ##########
#### External modules
import numpy as np
import pandas as pd

#### geodeZYX modules
from geodezyx import conv
from geodezyx import utils

log = logging.getLogger('geodezyx')

#  _____ _____  _____    _____                _     _             _   _                _____        __ _       ______ _ _
# |_   _/ ____|/ ____|  / ____|              | |   (_)           | | (_)              / ____|      / _| |     |  ____(_) |
#   | || |  __| (___   | |     ___  _ __ ___ | |__  _ _ __   __ _| |_ _  ___  _ __   | (___   ___ | |_| |_    | |__   _| | ___  ___
#   | || | |_ |\___ \  | |    / _ \| '_ ` _ \| '_ \| | '_ \ / _` | __| |/ _ \| '_ \   \___ \ / _ \|  _| __|   |  __| | | |/ _ \/ __|
#  _| || |__| |____) | | |___| (_) | | | | | | |_) | | | | | (_| | |_| | (_) | | | |  ____) | (_) | | | |_ _  | |    | | |  __/\__ \
# |_____\_____|_____/   \_____\___/|_| |_| |_|_.__/|_|_| |_|\__,_|\__|_|\___/|_| |_| |_____/ \___/|_|  \__(_) |_|    |_|_|\___||___/


def prn_int_2_prn_str(prn_int, full_out=False):
    """
    for read_combi_sum_full

    if full_out : return e.g. "G04","G",4
    """

    const = "X"

    prn_int = int(prn_int)

    prn_int_out = prn_int

    if prn_int >= 400:
        prn_int_out = prn_int - 400
        const = "J"
    elif 300 <= prn_int < 400:
        prn_int_out = prn_int - 300
        const = "C"
    elif 200 <= prn_int < 300:
        prn_int_out = prn_int - 200
        const = "E"
    elif 100 <= prn_int < 200:
        prn_int_out = prn_int - 100
        const = "R"
    else:
        prn_int_out = prn_int
        const = "G"

    prn_str = const + str(prn_int_out).zfill(2)

    if not full_out:
        return prn_str
    else:
        return prn_str, const, prn_int_out


def read_combi_sum_full(sum_full_file, RMS_lines_output=True, set_PRN_as_index=True):
    Vals_stk = []

    for l in open(sum_full_file):

        F = l.split()
        if "|" in F:
            F.remove("|")

        ### Find date line
        if "MJD:" in l:
            date_line = l

        ### Skip useless lines
        if not "|" in l or "------" in l:
            continue

        ### Find AC list
        if "PRN" in l:
            acs_list = F
            acs_list.append("RMS_sat")
            acs_list.append("PRN_str")
            acs_list.append("CONST")

        elif F[0].isnumeric():
            Fout = [float(f) for f in F]
            Fout[0] = int(Fout[0])
            # Add the PRN string and the constellation
            Fout.append(prn_int_2_prn_str(int(Fout[0])))
            Fout.append(Fout[-1][0])

            Vals_stk.append(Fout)

        elif "RMS" in F[0] and RMS_lines_output:
            Fout = [float(f) for f in F[1:]]
            Fout.append(np.nan)
            Fout.insert(0, F[0])
            # Add FAKE the PRN string and the constellation
            Fout.append(F[0])
            Fout.append(None)

            Vals_stk.append(Fout)

    DF = pd.DataFrame(Vals_stk, columns=acs_list)

    ### Date management
    mjd = float(date_line.split("MJD:")[1].split()[0])
    date_dt = conv.mjd2dt(mjd)

    DF.date_mjd = mjd
    DF.date_dt = date_dt

    DF.date_gps = utils.join_improved("", conv.dt2gpstime(date_dt))

    if set_PRN_as_index:
        DF.set_index("PRN_str", inplace=True)

    return DF


def read_combi_sum_exclu(sum_file, return_as_df=True, use_intuitive_bool=True):

    t_dt = conv.sp3name2dt(sum_file)

    with open(sum_file) as f:
        cont = f.readlines()

    excluded_dic = dict()
    useful_ssection = False
    useful_ssection_k = 0
    for l in cont:
        f = l.split()
        if "---|---" in l and useful_ssection_k < 2:
            useful_ssection = not useful_ssection
            useful_ssection_k += 1
            continue

        if not useful_ssection:
            continue

        prn_raw = f[0]

        if "X" in prn_raw:
            exclu = True
        else:
            exclu = False

        if use_intuitive_bool:
            exclu = not exclu

        prn_int = int(f[0].replace("X", "").split()[0])

        prn_good = prn_int_2_prn_str(prn_int)

        excluded_dic[prn_good] = exclu

    if return_as_df:
        return pd.DataFrame(excluded_dic, index=[t_dt])
    else:
        return excluded_dic


def read_combi_clk_rms(
    sum_file,
    return_as_df=True,
    clk_ref_cen_gal="com",
    index_useful_col=-4,
    convert_to_int=True,
):
    """
    based on : read_good_clk_rms_one
    """

    strt = " RESULTS OF FINAL WEIGHTED COMBINATION"
    end = " CLK_REF_CEN_GAL: " + clk_ref_cen_gal

    l = utils.extract_text_between_elements_2(sum_file, strt, end)

    l = l[:-2]

    lres = [e for e in l if re.search(r"^ [a-z]{3} \|", e)]

    lres_splited = [e.split() for e in lres]

    filnam = os.path.basename(sum_file)
    if "log" in filnam:
        week = int(filnam[4:8])
        dow = int(filnam[9])
        tdt = conv.gpstime2dt(week, dow)
    elif "cls" in filnam:
        week = int(filnam[3:7])
        dow = int(filnam[7])
        tdt = conv.gpstime2dt(week, dow)

    rms_dict = dict()

    for e in lres_splited:
        try:
            if convert_to_int:
                rms_dict[e[0]] = int(float(e[index_useful_col]))
            else:
                rms_dict[e[0]] = float(e[index_useful_col])
        except:
            log.warning("WARN : %s not handeled", e[index_useful_col])
            log.warning("replaced with NaN")
            rms_dict[e[0]] = np.nan

    if return_as_df:
        return pd.DataFrame(rms_dict, index=[tdt])
    else:
        return tdt, rms_dict


def read_combi_clk_rms_full_table(path_in, with_stats_rms=False, detailed_df=False):
    """
    recommended for .out file

    detailed_df: the outlier values are more detailled
    X (excuded) => np.inf
    - (not proivided) => np.nan
    >>> (too big for a print, but still kept) => 999999
    """
    strt = r"RMS \(ps\) OF AC CLOCK COMPARED TO COMBINATION"
    end = "---+---"

    if with_stats_rms:
        nth_occur = 2
    else:
        nth_occur = 1

    lines = utils.extract_text_between_elements_2(
        path_in, strt, end, nth_occur_elt_end=nth_occur
    )

    lines_good = []

    for l in lines[1:]:
        if "---+---" in l or "bad" in l:
            continue
        else:
            lines_good.append(l)

    lines_good = [e.replace("|", "") for e in lines_good]
    lines_good = [e.replace("         ", "  SAT    ") for e in lines_good]

    str = "".join(lines_good)

    import io

    ### Simple Mode
    if not detailed_df:
        df = pd.read_table(
            io.StringIO(str),
            na_values=["-", "X", ">>>"],
            delim_whitespace=True,
            error_bad_lines=False,
        )

    #### Mone detailled mode
    else:
        Cols = lines_good[0].split()[1:-2]

        #### We need an ad hoc fct to convert the values
        def conv_detailed(inp_val):
            if "-" in inp_val:
                out_val = np.nan
            elif "X" in inp_val:
                out_val = np.inf
            elif ">>>" in inp_val:
                out_val = 999999
            else:
                out_val = np.int64(inp_val)
            return out_val

        ### and then each column has to have its own convert fct...
        ### (quite stupid but it's the only way...)
        conv_dict = dict()
        for col in Cols:
            conv_dict[col] = conv_detailed

        df = pd.read_table(
            io.StringIO(str),
            delim_whitespace=True,
            error_bad_lines=False,
            converters=conv_dict,
        )

    df = df.set_index("SAT")

    return df


def read_combi_REPORT(Path_list):
    STK = []
    for p in Path_list:
        F = open(p)
        for l in F:
            f = l.split()
            if "epoch" in l:
                epoch = conv.gpstime2dt(int(f[2]), int(f[3]))
            if "orb_flag_x" in l:
                prn_str, const, prn_int = prn_int_2_prn_str((f[2]), True)
                STK.append((epoch, prn_str, const, prn_int, "all"))
            if "orb_excl_sat" in l:
                prn_str, const, prn_int = prn_int_2_prn_str((f[3]), True)
                STK.append((epoch, prn_str, const, prn_int, f[2]))

    DF = pd.DataFrame(STK, columns=("epoch", "PRN_str", "CONST", "PRN", "AC"))

    return DF
