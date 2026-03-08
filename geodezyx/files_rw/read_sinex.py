#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 27/08/2025 22:26:15

@author: psakic
"""

import re
import os
import itertools
import numpy as np
import pandas as pd
from io import StringIO
from geodezyx import conv
from geodezyx import utils


def read_snx_trop(snxfile, dataframe_output=True):
    """
    Read troposphere solutions from Troposphere SINEX
    """

    stat, epoc = [], []
    tro, stro, tgn, stgn, tge, stge = [], [], [], [], [], []

    flagtrop = False

    for line in open(snxfile, "r", encoding="ISO-8859-1"):

        if re.compile("TROP/SOLUTION").search(line):
            flagtrop = True
            continue

        if re.compile("-TROP/SOLUTION").search(line):
            flagtrop = False
            continue

        if line[0] != " ":
            continue
        else:
            fields = line.split()

        if flagtrop:

            stat.append(fields[0].upper())
            if not ":" in fields[1]:
                epoc.append(conv.year_decimal2dt(fields[1]))
            else:
                date_elts_lis = fields[1].split(":")
                yy = int(date_elts_lis[0]) + 2000
                doy = int(date_elts_lis[1])
                sec = int(date_elts_lis[2])

                epoc.append(conv.doy2dt(yy, doy, seconds=sec))

            tro.append(float(fields[2]))
            stro.append(float(fields[3]))
            tgn.append(float(fields[4]))
            stgn.append(float(fields[5]))
            tge.append(float(fields[6]))
            stge.append(float(fields[7]))

    outtuple = list(zip(*sorted(zip(stat, epoc, tro, stro, tgn, stgn, tge, stge))))

    if dataframe_output:
        return tropsinex_data_frame(outtuple)
    else:
        return outtuple


def tropsinex_data_frame(read_sinex_result):
    df_sinex = pd.DataFrame.from_records(list(read_sinex_result)).transpose()
    colnam = ["STAT", "epoc", "tro", "stro", "tgn", "stgn", "tge", "stge"]
    df_sinex.columns = colnam
    cols_numeric = ["tro", "stro", "tgn", "stgn", "tge", "stge"]
    df_sinex[cols_numeric] = df_sinex[cols_numeric].apply(
        pd.to_numeric, errors="coerce"
    )

    return df_sinex


def read_sinex(sinex_path_in, keep_sites_as_index=False, drop_nan=False):
    """
    Read a *coordinate* SINEX file as DataFrame

    Parameters
    ----------
    sinex_path_in : str
        path of the SINEX.
    keep_sites_as_index : Bool, optional
        use the site names as index. The default is False.
    drop_nan : bool or str, optional
        remove the NaN in the DataFrame.
        False: do not remove them
        'all': remove the row if all values are NaN
        'any':  remove the row if at least one value is NaN (dangerous)
        The default is False.

    Returns
    -------
    DFout : DataFrame
        SINEX DataFrame.

    """
    ## Read the blocs
    df_coor = read_sinex_versatile(
        sinex_path_in, "SOLUTION/ESTIMATE", header_line_idx=None
    )
    df_epoc = read_sinex_versatile(
        sinex_path_in, "SOLUTION/EPOCHS", header_line_idx=None
    )

    ## Rename the Epoch DF
    df_epoc.rename(
        columns={
            0: "STAT",
            1: "pt",
            2: "soln",
            3: "?",
            4: "start",
            5: "end",
            6: "mean",
        },
        inplace=True,
    )
    df_epoc.set_index(["STAT", "pt", "soln"], inplace=True)

    ### Rearange the coords DF
    index = df_coor[[2, 3, 4]].drop_duplicates()

    stk_df_coor2_line = []

    for stat, pt, soln in zip(index[2], index[3], index[4]):
        df2 = df_coor[(df_coor[2] == stat) & (df_coor[3] == pt) & (df_coor[4] == soln)]

        dicttmp = dict()

        epoc = df2[5].unique()[0]
        dicttmp["epoc"] = epoc
        dicttmp["STAT"] = stat
        dicttmp["pt"] = pt
        dicttmp["soln"] = soln

        for typ, comp in itertools.product(("STA", "VEL"), ("X", "Y", "Z")):

            typcomp = typ + comp
            try:
                val = df2[df2[1] == typcomp][8].values[0]
                sig = df2[df2[1] == typcomp][9].values[0]
            except:
                val, sig = np.nan, np.nan

            if typ == "VEL":
                typdic = "v"
            else:
                typdic = ""

            dicttmp[typdic + comp.lower()] = val
            dicttmp["s" + typdic + comp.lower()] = sig

        df_coor2_line = pd.DataFrame(dicttmp, index=[0])
        stk_df_coor2_line.append(df_coor2_line)

    df_coor2 = pd.concat(stk_df_coor2_line)
    df_coor2.set_index(["STAT", "pt", "soln"], inplace=True)

    ### Concatenate both
    DFout = pd.concat((df_epoc, df_coor2), axis=1)

    ### remove NaN (might be more epochs than coords...)
    if drop_nan:
        DFout.dropna(inplace=True, how=drop_nan)

    if not keep_sites_as_index:
        DFout.reset_index(inplace=True)

    return DFout


def read_sinex_legacy(snxfile, dataframe_output=True):
    """
    This function is depreciated !!!!
    """

    stat, soln, epoc, ac = [], [], [], []
    x, y, z, sx, sy, sz = [], [], [], [], [], []
    vx, vy, vz, svx, svy, svz = [], [], [], [], [], []

    start, end = [], []

    flagxyz = False
    flagepochs = False

    for line in open(snxfile, "r", encoding="ISO-8859-1"):

        if re.compile("SOLUTION/ESTIMATE").search(line):
            flagxyz = not flagxyz
            continue

        if re.compile("SOLUTION/EPOCHS").search(line):
            flagepochs = not flagepochs
            continue

        if line[0] != " ":
            continue
        else:
            fields = line.split()

        if flagxyz:

            if fields[1] == "STAX":
                split_sp3_name = os.path.basename(snxfile)
                ac.append(split_sp3_name[0:3])
                stat.append(fields[2].upper())
                soln.append(fields[4])
                if not ":" in fields[5]:
                    epoc.append(conv.year_decimal2dt(fields[5]))
                else:
                    date_elts_lis = fields[5].split(":")
                    yy = int(date_elts_lis[0]) + 2000
                    doy = int(date_elts_lis[1])
                    sec = int(date_elts_lis[2])

                    epoc.append(conv.doy2dt(yy, doy, seconds=sec))

                x.append(float(fields[8]))
                sx.append(float(fields[9]))

            if fields[1] == "STAY":
                y.append(float(fields[8]))
                sy.append(float(fields[9]))

            if fields[1] == "STAZ":
                z.append(float(fields[8]))
                sz.append(float(fields[9]))

            if fields[1] == "VELX":
                vx.append(float(fields[8]))
                svx.append(float(fields[9]))

            if fields[1] == "VELY":
                vy.append(float(fields[8]))
                svy.append(float(fields[9]))

            if fields[1] == "VELZ":
                vz.append(float(fields[8]))
                svz.append(float(fields[9]))

        if flagepochs:
            start_val = conv.datestr_sinex_2_dt(fields[4])
            end_val = conv.datestr_sinex_2_dt(fields[5])
            start.append(start_val)
            end.append(end_val)

    if not vx:
        nanlist = [np.nan] * len(x)
        vx, vy, vz, svx, svy, svz = (
            list(nanlist),
            list(nanlist),
            list(nanlist),
            list(nanlist),
            list(nanlist),
            list(nanlist),
        )

    outtuple = list(
        zip(
            *sorted(
                zip(
                    ac,
                    stat,
                    soln,
                    epoc,
                    x,
                    y,
                    z,
                    sx,
                    sy,
                    sz,
                    vx,
                    vy,
                    vz,
                    svx,
                    svy,
                    svz,
                    start,
                    end,
                )
            )
        )
    )

    if dataframe_output:
        return sinex_data_frame(outtuple)
    else:
        stat, soln, epoc, x, y, z, sx, sy, sz, vx, vy, vz, svx, svy, svz, start, end = (
            outtuple
        )
        print(
            "TIPS : this output can be converted directly to a Pandas DataFrame using sinex_data_frame function"
        )
        return (
            stat,
            soln,
            epoc,
            x,
            y,
            z,
            sx,
            sy,
            sz,
            vx,
            vy,
            vz,
            svx,
            svy,
            svz,
            start,
            end,
        )


def sinex_data_frame(read_sinex_result):
    DF_Sinex = pd.DataFrame.from_records(list(read_sinex_result)).transpose()
    colnam = [
        "AC",
        "STAT",
        "soln",
        "epoc",
        "x",
        "y",
        "z",
        "sx",
        "sy",
        "sz",
        "vx",
        "vy",
        "vz",
        "svx",
        "svy",
        "svz",
        "start",
        "end",
    ]
    DF_Sinex.columns = colnam

    cols_numeric = [
        "x",
        "y",
        "z",
        "sx",
        "sy",
        "sz",
        "vx",
        "vy",
        "vz",
        "svx",
        "svy",
        "svz",
    ]

    DF_Sinex[cols_numeric] = DF_Sinex[cols_numeric].apply(
        pd.to_numeric, errors="coerce"
    )

    return DF_Sinex


def read_sinex_versatile(
    sinex_path_in,
    id_block,
    convert_date_2_dt=True,
    header_line_idx=-1,
    improved_header_detection=True,
    verbose=True,
):
    """
    Read a block from a SINEX and return the data as a DataFrame

    Parameters
    ----------
    sinex_path_in : str
        Description param1

    id_block : str
        Name of the block (without "+" or "-")

    convert_date_2_dt : bool
        Try to convert a SINEX formated date as a python datetime

    header_line_idx : int or None
        If the block header contains several lines, use this line index
        Per default, the last (-1)
        For the first line, use 0
        If no header is properly defined, use None

    improved_header_detection : bool
        Improved header detection.
        Works for most cases but sometime the simple version works better.
        (advenced usage)
        Default is True

    verbose : bool
        print the header and its field size
        Default is True

    Returns
    -------
    df : Pandas DataFrame
        Returned DataFrame
    """

    ### remove the + or - if any
    if id_block in ("+", "-"):
        id_block = id_block[1:]

    id_block_strt = r"\+" + id_block
    id_block_end = r"\-" + id_block

    lines_list = utils.extract_text_between_elements_2(
        sinex_path_in, id_block_strt, id_block_end
    )
    lines_list = lines_list[1:-1]

    if not lines_list:
        print("ERR : read_sinex_versatile : no block found, ", id_block)

    #### Remove commented lines
    lines_list_header = []
    lines_list_ok = []
    header_lines = True
    for i_l, l in enumerate(lines_list):
        if not l[0] in (" ", "\n") and header_lines:
            ## here we store the 1st commented lines i.e. the header
            lines_list_header.append(l)
        elif l[0] in (" ", "\n"):
            ## here we store the data lines (not commented)
            header_lines = False
            lines_list_ok.append(l)
        else:
            ## here we skip the final commented lines (can happend)
            continue

    if len(lines_list_header) > 0 and header_line_idx:
        ### define the header
        header_line = lines_list_header[header_line_idx]

        header_split = header_line.split()
        if not improved_header_detection:
            ### Simple case when the columns are splitted with only a single
            fields_size = [len(e) + 1 for e in header_split]
        else:
            ### Smarter case : we search for the n spaces after the column name
            fields_size = []

            for fld_head_split in header_split:
                if fld_head_split[0] == "*":
                    ## test to manage a * as 1st character
                    ## which can be wrongly interpreted in the regex
                    fld_head_regex = re.compile(r"\*" + fld_head_split[1:] + " *")
                else:
                    fld_head_regex_str = fld_head_split + " *"
                    ### update 202110: now the brackets must be escaped (normal btw...)
                    fld_head_regex_str = fld_head_regex_str.replace("[", r"\[")
                    fld_head_regex_str = fld_head_regex_str.replace("]", r"\]")
                    fld_head_regex_str = fld_head_regex_str.replace("(", r"\(")
                    fld_head_regex_str = fld_head_regex_str.replace(")", r"\)")

                    fld_head_regex = re.compile(fld_head_regex_str)

                fld_head_space = fld_head_regex.search(header_line)

                fields_size.append(len(fld_head_space.group()))
                # print(fld_head_space.group())

                # # weak method (210216) archived for legacy
                # fld_head_regex = re.compile(fld_head_split[1:] + " *") #trick:
                # #1st char is removed, because it can be a *
                # #and then screw the regex. This char is re-added at the end
                # #when the len is stored (the "+1" below)
                # fld_head_space = fld_head_regex.search(header_line)
                # fields_size.append(len(fld_head_space.group()) + 1)
                # ### !!!!! something is weird here !!!!!
                # print(fld_head_space.group())
                # ### and you will see a bug !!!
                # ### PS 210216

        if verbose:
            print("INFO : read_sinex_versatile : Auto detected column names/sizes")
            print("**** Raw header line in the file:")
            print(header_line)
            print("**** Splited header for the DataFrame:")
            print(header_split)
            print("**** Size of the fields")
            print(fields_size)

        ### Add the header in the big string
        lines_str_w_head = header_line + "".join(lines_list_ok)

        ### Read the file
        try:
            df = pd.read_fwf(StringIO(lines_str_w_head), widths=fields_size)
        except pd.errors.EmptyDataError as ee:
            print("ERR: something goes wrong in the header index position")
            print("     try to give its right position manually with header_line_idx")
            raise (ee)

        df.set_axis(header_split, axis=1, inplace=True)

        ### Rename the 1st column (remove the comment marker)
        df.rename(columns={df.columns[0]: df.columns[0][1:]}, inplace=True)

    else:  # no header in the SINEX
        lines_str = "".join(lines_list_ok)
        df = pd.read_csv(StringIO(lines_str), header=None, delim_whitespace=True)

    regex_time = "(([0-9]{2}|[0-9]{4}):[0-9]{3}|[0-9]{7}):[0-9]{5}"
    for col in df.columns:
        if convert_date_2_dt and re.match(regex_time, str(df[col].iloc[0])):
            try:
                df[col] = df[col].apply(lambda x: conv.datestr_sinex_2_dt(x))
            except Exception as e:
                print(
                    "WARN : read_sinex_versatile : convert date string to datetime failed"
                )
                print(e)
                pass

    return df


def read_sinex_bench_antenna(sinex_in):
    f = open(sinex_in, "r")

    line_good_stk = []
    for l in f:
        if l[0] == " ":
            line_good_stk.append(l)

    T = StringIO("".join(line_good_stk))

    header_cols = [
        "TYP",
        "STA_",
        "Occ",
        "Code",
        "PT",
        "T",
        "___x/up___",
        "___y/n____",
        "___z/e____",
        "_Data_Start",
        "_Data_End__",
        "__Antenna_type______",
        "Radome",
        "__S/N__",
    ]
    df_antenna = pd.read_table(
        T, delim_whitespace=True, error_bad_lines=False, names=header_cols
    )

    return df_antenna