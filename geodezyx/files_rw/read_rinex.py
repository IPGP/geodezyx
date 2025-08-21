#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.files_rw contains functions to 
read RINEX files observation files.

it can be imported directly with:
from geodezyx import files_rw

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""

#### Import the logger
import logging
import os
import pathlib
import re
from io import StringIO

import hatanaka

########## BEGIN IMPORT ##########
#### External modules
import numpy as np
import pandas as pd
from tqdm import tqdm

#### geodeZYX modules
from geodezyx import operational, utils, conv

log = logging.getLogger('geodezyx')


def read_rinex_obs(rnx_in, set_index=None , return_header = False):
    """
    Frontend function to read a RINEX Observation, version 2 or 3

    Parameters
    ----------
    rnx_in : see below
        input RINEX.
        can be the path of a RINEX file as string or as Path object,
        or directly the RINEX content as a string, bytes, StringIO object or a
        list of lines
    set_index : str or list of str, optional
        define the columns for the index.
        If None, the output DataFrame is "flat", with integer index
        ["epoch","prn"] for instance set the epoch and the prn as MultiIndex
        The default is None.
    return_header : bool, optional
        If True, the function returns a tuple (df_rnx_obs, header_lines),
        where header_lines is a list of strings containing the header lines.

    Returns
    -------
    df_rnx_obs : Pandas DataFrame / GeodeZYX's RINEX format
    header_lines : list of str, optional
        If return_header is True, this contains the header lines of the RINEX file.
        Otherwise, it is not returned.
    """

    if conv.rinex_regex_search_tester(rnx_in,short_name=True,long_name=False):
        return read_rinex2_obs(rnx_in, set_index=set_index, return_header=return_header)
    elif conv.rinex_regex_search_tester(rnx_in,short_name=False,long_name=True):
        return read_rinex3_obs(rnx_in, set_index=set_index, return_header=return_header)
    else:
        log.error("RINEX version not recognized: %s", rnx_in)
        log.error("use read_rinex2_obs or read_rinex3_obs functions manually")
        return None


def read_rinex2_obs(rnx_in, set_index=None, return_header = False):
    """
    Read a RINEX Observation, version 2

    Parameters
    ----------
    rnx_in : see below
        input RINEX.
        can be the path of a RINEX file as string or as Path object,
        or directly the RINEX content as a string, bytes, StringIO object or a
        list of lines
    set_index : str or list of str, optional
        define the columns for the index.
        If None, the output DataFrame is "flat", with integer index
        ["epoch","prn"] for instance set the epoch and the prn as MultiIndex
        The default is None.
    return_header : bool, optional
        If True, the function returns a tuple (df_rnx_obs, header_lines),
        where header_lines is a list of strings containing the header lines.

    Returns
    -------
    df_rnx_obs : Pandas DataFrame / GeodeZYX's RINEX format
    header_lines : list of str, optional
        If return_header is True, this contains the header lines of the RINEX file.
        Otherwise, it is not returned.
    """

    #### open rinex
    lines, epochs_dt, filename = open_rinex(rnx_in)

    ### get the header
    l_head = _header_reader(lines)

    # get the rinex type (mono or multi GNSS)
    rnx_typ = l_head[0][40]

    mixed_rnx = rnx_typ == "M"
    # True if M, False if not

    #### get the observables
    l_obs = [l for l in l_head if "# / TYPES OF OBSERV" in l]
    ## clean SYS / # / OBS TYPES
    l_obs = [l[:60] for l in l_obs]

    obs_all_list_raw = " ".join(l_obs).split()
    obs_all_list = obs_all_list_raw[1:]
    obs_all_list = [
        e
        for sublist in [(e, e + "_LLI", e + "_SSI") for e in obs_all_list]
        for e in sublist
    ]

    nobs = int(obs_all_list_raw[0])
    nl_obs = int(np.ceil(nobs / 5))
    ## 5 is the max num of obs in the RINEX specs
    col_width = nobs * [14, 1, 1]

    df_all_stk = []

    #### reading the epochs
    for iepoc in tqdm(range(len(epochs_dt)), desc="Reading " + filename):
        epoch = epochs_dt[iepoc, 0]
        ## define the start/end indices of the epoch block
        il_srt = epochs_dt[iepoc, 1]
        if iepoc == len(epochs_dt) - 1:
            il_end = None
        else:
            il_end = epochs_dt[iepoc + 1, 1]

        ### get the lines of the epoch block
        lines_epoc = lines[il_srt:il_end]

        ### get the satellites for this epoch block
        sats_find_out = _sats_find(lines_epoc)
        flag, epoc, nsat, _, sats_split, iline_sats_end = sats_find_out        

        if flag != 0: # skip the epoch (which is not an observable one)
            log.warning("epoch %s skipped, flag %i",epoc, flag)
            continue

        ### for each sat, merge the breaked lines
        resub80 = lambda x: re.sub(r"[\r\n]+", "", x.ljust(80))
        # .replace("\r", "").replace("\n", "") i.e. re.sub(r"[\r\n]+", "", s) IS NOT .strip() !!
        # .strip() removes the spaces at the beginning and end of the string, and they should be kept!
        l_obs = [ resub80(e) for e in lines_epoc[iline_sats_end + 1: il_end]]
        l_obs_merg = ["".join(l_obs[n * nl_obs:(n + 1) * nl_obs]) for n in range(nsat)]

        ## read the epoch block using pandas' fixed width reader
        b = StringIO("\n".join(l_obs_merg))
        df_epoch = pd.read_fwf(b, header=None, widths=col_width)
        df_epoch.columns = obs_all_list

        # slow
        # df_epoch["prn"] = sats_split
        # df_epoch["epoch"] = epoch
        # faster
        col_epoch = pd.Series([epoch] * len(df_epoch), name="epoch")
        col_prn = pd.Series(sats_split, name="prn")

        if np.any(col_prn.str.len() == 0):
            log.critical("bad epoch %s" , epoc)
            raise Exception()

        df_epoch = pd.concat([col_epoch, col_prn, df_epoch], axis=1)

        df_all_stk.append(df_epoch)

    ## final concat and cosmetic (reorder columns, sort)
    df_rnx_obs = pd.concat(df_all_stk)

    # strip to have just 2 or 3 prn characters
    df_rnx_obs["prn"] = df_rnx_obs["prn"].str.strip()

    # mixed, system is in the epoch line
    if mixed_rnx or df_rnx_obs["prn"].str.len().max() == 3:
        df_rnx_obs["sys"] = df_rnx_obs["prn"].str[0]
    else:  # not mixed, system is the header RINEX type
        df_rnx_obs["sys"] = rnx_typ
        df_rnx_obs["prn"] = rnx_typ + df_rnx_obs["prn"].str.strip()

    df_rnx_obs["prni"] = df_rnx_obs["prn"].str[-2:].astype(int)

    main_cols = ["epoch", "sys", "prn", "prni"]
    df_rnx_obs = df_rnx_obs.reindex(main_cols + list(sorted(obs_all_list)), axis=1)
    df_rnx_obs.sort_values(["epoch", "prn"], inplace=True)
    df_rnx_obs.reset_index(drop=True, inplace=True)

    if set_index:
        df_rnx_obs.set_index(set_index, inplace=True)
        df_rnx_obs.sort_index(inplace=True)

    if return_header:
        return df_rnx_obs, l_head
    else:
        return df_rnx_obs

def read_rinex3_obs(rnx_in, set_index=None, return_header = False):
    """

    Read a RINEX Observation, version 3 or 4

    Parameters
    ----------
    rnx_in : see below
        input RINEX.
        can be the path of a RINEX file as string or as Path object,
        or directly the RINEX content as a string, bytes, StringIO object or a
        list of lines
    set_index : str or list of str, optional
        define the columns for the index.
        If None, the output DataFrame is "flat", with integer index
        ["epoch","prn"] for instance set the epoch and the prn as MultiIndex
        The default is None.
    return_header : bool, optional
        If True, the function returns a tuple (df_rnx_obs, header_lines),
        where header_lines is a list of strings containing the header lines.

    Returns
    -------
    df_rnx_obs : Pandas DataFrame / GeodeZYX's RINEX format
    header_lines : list of str, optional
        If return_header is True, this contains the header lines of the RINEX file.
        Otherwise, it is not returned.
    """

    #### open rinex
    lines, epochs, filename = open_rinex(rnx_in)

    ### get the header
    l_head = _header_reader(lines)

    #### get the systems and observations
    l_sys = [l for l in l_head if "SYS / # / OBS TYPES" in l]
    ## clean SYS / # / OBS TYPES
    l_sys = [l[:60] for l in l_sys]

    ## manage the 2 lines systems
    for il, l in enumerate(l_sys):
        if l[0] == " ":
            l_sys[il - 1] = l_sys[il - 1] + l
            l_sys.remove(l)

    #### store system and observables in a dictionnary
    dict_sys_nobs = dict()
    dict_sys_cols = dict()  # adds the prn, and LLI and SSI indicators
    obs_all_list = []

    for il, l in enumerate(l_sys):
        sysobs = l.split()
        # sysobs[0] = system letter
        # sysobs[1] = observables number
        # sysobs[2:] = observables names

        dict_sys_nobs[sysobs[0]] = int(sysobs[1])
        ## adds the LLI and SSI indicators
        obs_w_lli_ssi = [(e, e + "_LLI", e + "_SSI") for e in sysobs[2:]]
        obs_all_list.extend([e for sublist in obs_w_lli_ssi for e in sublist])

        dict_sys_cols_tmp = [("sys","epoch","prn",)] + obs_w_lli_ssi
        dict_sys_cols[sysobs[0]] = [e for sublist in dict_sys_cols_tmp for e in sublist]

        if len(sysobs[2:]) != int(sysobs[1]):
            log.warning(
                "difference between theorectical and actual obs nbr in the header"
            )

    ## the max number of observable (for the reading)
    nobs_max = max(dict_sys_nobs.values())

    df_all_stk = []
    error_sys_in_obs = False

    #### reading the epochs
    for iepoc in tqdm(range(len(epochs)), desc="Reading " + filename):
        epoch = epochs[iepoc, 0]
        
        ## check the flag
        flag = int(lines[epochs[iepoc, 1]][31])
        if flag != 0: # skip the epoch (which is not an observable one)
            log.warning("epoch %s skipped, flag %i", epoch, flag)
            continue

        ## define the start/end indices of the epoch block
        iline_start = epochs[iepoc, 1] + 1
        if iepoc == len(epochs) - 1:
            iline_end = None
        else:
            iline_end = epochs[iepoc + 1, 1]

        l_epoc = [re.sub(r"[\r\n]+", "", e) for e in lines[iline_start:iline_end]]
        # .replace("\r", "").replace("\n", "") i.e. re.sub(r"[\r\n]+", "", s) IS NOT .strip() !!
        # .strip() removes the spaces at the beginning and end of the string, and they should be kept!

        for line in l_epoc:
            if "SYS / # / OBS TYPES" in line:
                error_sys_in_obs = True
                break

        ## read the epoch block using pandas' fixed width reader
        b = StringIO("\n".join(l_epoc))

        columns_width = [3] + nobs_max * [14, 1, 1]
        df_epoch = pd.read_fwf(b, header=None, widths=columns_width)

        col_epoch = pd.Series([epoch] * len(df_epoch), name="epoch")
        df_epoch_ok = pd.concat([col_epoch, df_epoch], axis=1)

        df_all_stk.append(df_epoch_ok)

    if error_sys_in_obs:
        log.error(
            "SYS / # / OBS TYPES in body, observable's order is modified, read_rinex does not manage this yet")
        return None

    ## final concat and cosmetic (reorder columns, sort)
    df_rnx_obs = pd.concat(df_all_stk)

    ## extract the system
    col_sys = pd.Series(df_rnx_obs[0].str.strip().str[0], name="sys")
    df_rnx_obs = pd.concat([col_sys, df_rnx_obs], axis=1)

    #### assign the correct observable names for each system
    df_rnx_obs_ok_stk = []
    for sys, df_rnx_obs_sys in df_rnx_obs.groupby("sys"):
        df_rnx_obs_sys_clean = df_rnx_obs_sys.iloc[:, : len(dict_sys_cols[sys])]
        df_rnx_obs_sys_clean.columns = dict_sys_cols[sys] # sys, epoch & prn are already in dict_sys_cols
        df_rnx_obs_ok_stk.append(df_rnx_obs_sys_clean)

    df_rnx_obs_ok = pd.concat(df_rnx_obs_ok_stk)

    col_prni = pd.Series(df_rnx_obs_ok["prn"].str.strip().str[-2:], name="prni", dtype=int)
    df_rnx_obs_ok = pd.concat([df_rnx_obs_ok, col_prni], axis=1)

    main_cols = ["epoch", "sys", "prn", "prni"]
    obs_used_list = list(set(df_rnx_obs_ok.columns).difference(set(main_cols)))
    df_rnx_obs_ok = df_rnx_obs_ok.reindex(main_cols + list(sorted(obs_used_list)), axis=1)

    df_rnx_obs_ok.sort_values(["epoch", "prn"], inplace=True)
    df_rnx_obs_ok.reset_index(drop=True, inplace=True)

    if set_index:
        df_rnx_obs_ok.set_index(set_index, inplace=True)
        df_rnx_obs_ok.sort_index(inplace=True)


    if return_header:
        return df_rnx_obs_ok, l_head
    else:
        return df_rnx_obs_ok


############ UTILITY FUNCTIONS


def open_rinex(rnx_inp,verbose=False):
    #### open block
    try:
        rnx_wrk = hatanaka.decompress(rnx_inp)
    except:
        rnx_wrk = rnx_inp
        pass

    lines = utils.open_readlines_smart(rnx_wrk, verbose=verbose)
    epochs_dt = operational.rinex_read_epoch(rnx_wrk, out_index=True)
    if type(rnx_inp) is str or type(rnx_inp) is pathlib.Path:
        filename = os.path.basename(rnx_inp)
    else:
        filename = "unknown filename"

    return lines, epochs_dt, filename


def df_rnx_clean_lli_ssi(df_rnx_in):
    """
    Remove the Loss of Lock Indicator (LLI) and Signal Strength Indicator (SSI)
    columns in a DataFrame RINEX

    Parameters
    ----------
    df_rnx_in : Pandas DataFrame / GeodeZYX's RINEX format
        A RINEX DataFrame with LLI/SSI columns.

    Returns
    -------
    df_rnx_out : Pandas DataFrame / GeodeZYX's RINEX format
        A RINEX DataFrame without LLI/SSI columns.
    """
    cols = df_rnx_in.columns
    cols_clean = [e for e in cols if not "LLI" in e and not "SSI" in e]
    df_rnx_out = df_rnx_in[cols_clean]
    return df_rnx_out


def observables_dict_per_sys(df_rnx_in):
    """
    Gives the GNSS observables for each GNSS system in a dictionnary


    Parameters
    ----------
    df_rnx_in : Pandas DataFrame / GeodeZYX's RINEX format
        A RINEX DataFrame.

    Returns
    -------
    dict_sys_obs : dict
        A dictionnary with GNSS system as key (G,R,E...).
        And the observalbes for each system as values

    Note
    ----

    Use dict_sys_obs_clean_LLI_SSI if your want to remove the LLI & SSI values
    """

    dict_sys_obs = dict()

    for sys in df_rnx_in["sys"].unique():
        df_sys = df_rnx_in[df_rnx_in["sys"] == sys]
        df_sys_mini = df_rnx_clean_lli_ssi(df_sys)

        obs_sys0 = pd.Series(df_sys_mini.isna().sum() == len(df_sys)).apply(np.logical_not)
        obs_sys = obs_sys0.index[obs_sys0][
            3:
        ]  ### strating from 3 to clean epoch sys prn

        # init_tup = ("epoch","sys","prn")
        # init_tup = []
        # obs_sys_full = [init_tup] + [(e,e+"_LLI",e+"_SSI") for e in obs_sys]
        obs_sys_full = [(e, e + "_LLI", e + "_SSI") for e in obs_sys]
        obs_sys_full = [e for sublist in obs_sys_full for e in sublist]

        dict_sys_obs[sys] = obs_sys_full

    return dict_sys_obs


def dict_sys_obs_clean_lli_ssi(dict_sys_obs_in):
    """
    Clean a `dict_sys_obs` (generated by `observables_dict_per_sys`)
    of its LLI and SSI values

    Parameters
    ----------
    dict_sys_obs_in : dict
        A dictionnary with GNSS system as key (G,R,E...).
        And the observalbes for each system as values.

    Returns
    -------
    dict_sys_obs_out : dict
        Same dictionnary cleanned of its LLI and SSI values.

    """
    dict_sys_obs_out = dict()

    for sys, obs in dict_sys_obs_in.items():
        obs_clean = [e for e in obs if not "LLI" in e and not "SSI" in e]
        dict_sys_obs_out[sys] = obs_clean

    return dict_sys_obs_out


############ INTERNAL FUNCTIONS


def _sats_find(lines_inp):
    """
    For RINEX2 only
    search for the satellites for each epoch block by reading the
    EPOCH/SAT record

    Parameters
    ----------
    lines_inp : List of str
        the lines of one epoch block (EPOCH/SAT + OBSERVATION records).

    Returns
    -------
    bloc_tuple : tuple
        a 5-tuple containing:
            epoc : datetime, the epoch of block
            nsat : int, the number of satellites
            lineconcat : str, the satellites as a concatenated string
            sats_split : list of str, the satellites as a list
            il : int, the index of the last line of the EPOCH/SAT record

    """
    iline_bloc = 0
    nlines_bloc = -1

    re_epoch = "^ {1,2}([0-9]{1,2} * ){5}"

    ### dummy initialisation of variables
    line_bloc = []
    lineconcat_stk = []
    epoc = None
    nsat = 0

    for il, l in enumerate(lines_inp):
        #####
        # re_sat="[A-Z][0-9][0-9]"
        # bool_epoch = re.search(re_epoch, l)
        # bool_sat=re.search(re_sat,l)

        ### we found an epoch line
        if re.search(re_epoch, l):
            flag = int(l[28])
            nsat = int(l[30:32])
            iline_bloc = 0
            nlines_bloc = int(np.ceil(nsat / 12))
            line_bloc = []
            #lineconcat = ""
            lineconcat_stk = []
            epoc = operational.read_rnx_epoch_line(l, rnx2=True)

        ### we read the sat lines based on the number of sat
        if iline_bloc <= nlines_bloc:
            line_bloc.append(l[32:].strip())
            # lineconcat = lineconcat + l[32:].strip()
            lineconcat_stk.append(l[32:].strip())
            iline_bloc += 1

        ### we stack everything when the sat block is over
        if iline_bloc == nlines_bloc:
            #in_epoch = False
            lineconcat = "".join(lineconcat_stk)
            sats_split = [lineconcat[3 * n : 3 * n + 3] for n in range(nsat)]
            bloc_tuple = (flag, epoc, nsat, lineconcat, sats_split, il)
            return bloc_tuple


def _header_reader(lines_inp):
    #### Split header and Observation body
    i_end_header = 0
    for il, l in enumerate(lines_inp):
        if "END OF HEADER" in l:
            i_end_header = il
            break
    return lines_inp[: i_end_header + 1]


### FUNCTION GRAVEYARD


def _line_reader(linein, nobs):
    """
    DISCONTINUED

    For RINEX3 Only

    read the content of an observation line. pd.read_fwf does the job
    """
    columns_width = [3] + nobs * [14, 1, 1]
    columns_cumul = [0] + list(np.cumsum(columns_width))

    out_stk = []
    for i in range(len(columns_cumul) - 1):
        slic = slice(columns_cumul[i], columns_cumul[i + 1])
        out = linein[slic]
        if i == 0:
            out_stk.append(out)
        else:
            try:
                out_stk.append(float(out))
            except ValueError:
                out_stk.append(np.nan)

    return out_stk

def read_rinex3_obs_legacy(rnx_in, set_index=None):
    """
    Read a RINEX Observation, version 3 or 4

    This legacy version is much slower but maybe mor relaible than the new one

    Parameters
    ----------
    rnx_in : see below
        input RINEX.
        can be the path of a RINEX file as string or as Path object,
        or directly the RINEX content as a string, bytes, StringIO object or a
        list of lines
    set_index : str or list of str, optional
        define the columns for the index.
        If None, the output DataFrame is "flat", with integer index
        ["epoch","prn"] for instance set the epoch and the prn as MultiIndex
        The default is None.

    Returns
    -------
    df_rnxobs : Pandas DataFrame / GeodeZYX's RINEX format
    """

    #### open rinex
    lines, epochs, filename = open_rinex(rnx_in)

    ### get the header
    lines_header = _header_reader(lines)

    #### get the systems and observations
    lines_sys = [l for l in lines_header if "SYS / # / OBS TYPES" in l]
    ## clean SYS / # / OBS TYPES
    lines_sys = [l[:60] for l in lines_sys]

    ## manage the 2 lines systems
    for il, l in enumerate(lines_sys):
        if l[0] == " ":
            lines_sys[il - 1] = lines_sys[il - 1] + l
            lines_sys.remove(l)

    #### store system and observables in a dictionnary
    dict_sys_nobs = dict()
    dict_sys_cols = dict()  # adds the prn, and LLI and SSI indicators
    obs_all_list = []

    for il, l in enumerate(lines_sys):
        sysobs = l.split()
        # sysobs[0] = system letter
        # sysobs[1] = observables number
        # sysobs[2:] = observables names

        dict_sys_nobs[sysobs[0]] = int(sysobs[1])
        ## adds the LLI and SSI indicators
        obs_w_lli_ssi = [(e, e + "_LLI", e + "_SSI") for e in sysobs[2:]]
        obs_all_list.extend([e for sublist in obs_w_lli_ssi for e in sublist])

        dict_sys_cols_tmp = [("prn",)] + obs_w_lli_ssi
        dict_sys_cols[sysobs[0]] = [e for sublist in dict_sys_cols_tmp for e in sublist]

        if len(sysobs[2:]) != int(sysobs[1]):
            log.warning(
                "difference between theorectical and actual obs nbr in the header"
            )

    ## the max number of observable (for the reading)
    nobs_max = max(dict_sys_nobs.values())

    df_all_stk = []
    #### reading the epochs
    for iepoc in tqdm(range(len(epochs)), desc="Reading " + filename):
        epoch = epochs[iepoc, 0]
        ## define the start/end indices of the epoch block
        iline_start = epochs[iepoc, 1] + 1
        if iepoc == len(epochs) - 1:
            iline_end = None
        else:
            iline_end = epochs[iepoc + 1, 1]

        lines_epoc = [re.sub(r"[\r\n]+", "", e) for e in lines[iline_start:iline_end]]

        # .replace("\r", "").replace("\n", "") i.e. re.sub(r"[\r\n]+", "", s) IS NOT .strip() !!
        # .strip() removes the spaces at the beginning and end of the string, and they should be kept!

        ## read the epoch block using pandas' fixed width reader
        b = StringIO("\n".join(lines_epoc))

        columns_width = [3] + nobs_max * [14, 1, 1]
        df_epoch = pd.read_fwf(b, header=None, widths=columns_width)

        # one add the system to use groupby
        #df_epoch_sysadded = df_epoch.copy()
        #df_epoch_sysadded["systmp"] = df_epoch_sysadded[0].str[0]

        df_epoch_ok_stk = []
        #### assign the correct observable names for each system
        for sys in dict_sys_cols.keys():
        #for sys, df_epoch_sys in df_epoch_sysadded.groupby("systmp"):
            df_epoch_sys = df_epoch[df_epoch[0].str[0] == sys]
            df_epoch_sys_clean = df_epoch_sys.iloc[:, : len(dict_sys_cols[sys])]
            df_epoch_sys_clean.columns = dict_sys_cols[sys]
            df_epoch_ok_stk.append(df_epoch_sys_clean)

        df_epoch_ok = pd.concat(df_epoch_ok_stk)
        # An epoch column is created to fasten the process
        col_epoch = pd.Series([epoch] * len(df_epoch_ok), name="epoch")
        col_sys = pd.Series(df_epoch_ok["prn"].str[0], name="sys")
        col_prni = pd.Series(df_epoch_ok["prn"].str[1:], name="prni", dtype=int)
        df_epoch_ok = pd.concat([col_epoch, col_sys, col_prni, df_epoch_ok], axis=1)

        df_all_stk.append(df_epoch_ok)

    ## final concat and cosmetic (reorder columns, sort)
    df_rnx_obs = pd.concat(df_all_stk)
    main_cols = ["epoch", "sys", "prn", "prni"]
    obs_used_list = list(set(df_rnx_obs.columns).difference(set(main_cols)))
    df_rnx_obs = df_rnx_obs.reindex(main_cols + list(sorted(obs_used_list)), axis=1)

    df_rnx_obs.sort_values(["epoch", "prn"], inplace=True)
    df_rnx_obs.reset_index(drop=True, inplace=True)

    if set_index:
        df_rnx_obs.set_index(set_index, inplace=True)
        df_rnx_obs.sort_index(inplace=True)

    return df_rnx_obs


if __name__ == "__main__" and False:

    # p3 = "/home/psakicki/aaa_FOURBI/RINEXexemple/ENOG00REU_R_20243030000_01D_30S_MO.rnx"
    # df3 = read_rinex3_obs(p3)
    # print(df3)

    time_readrinex = """
from __main__ import read_rinex2_obs
p = "/home/psakicki/aaa_FOURBI/RINEXexemple/smne176z.18o"
p = "/home/psakicki/aaa_FOURBI/RINEXexemple/mlvl176z.18o"
p = "/home/psakicki/Desktop/daej2220.04o"
read_rinex2_obs(p)
"""

    # import timeit
    #
    # t = timeit.timeit(time_readrinex, number=10)
    # print(t)


    #p2 = "/home/psakicki/Desktop/daej2220.04o"
    #df2 = read_rinex2_obs(p2)
    #print(df2)

    p3 = "/home/psakicki/aaa_FOURBI/RINEXexemple/ENOG00REU_R_20243030000_01D_30S_MO.rnx"
    p3 = "/home/psakicki/aaa_FOURBI/RINEXexemple/CFNG00REU_R_20240430000_15M_01S_MO.rnx"
    p3 = "/home/psakicki/aaa_FOURBI/RINEXexemple/REYK00ISL_R_20200420000_01D_30S_MO.rnx"
    df31 = read_rinex3_obs(p3)
    df32 = read_rinex3_obs_legacy(p3)
    print(df31)
    print(df32)
