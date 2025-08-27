#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.reffram contains functions for operations
related to GNSS-products


it can be imported directly with:
from geodezyx import reffram

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""

########## BEGIN IMPORT ##########
#### External modules
import datetime as dt
import itertools

#### Import the logger
import logging
import os
import re

import matplotlib
import matplotlib.pyplot as plt
import natsort
import numpy as np
import pandas as pd

#### geodeZYX modules
from geodezyx import conv
from geodezyx import files_rw
from geodezyx import stats
from geodezyx import utils
from geodezyx.reffram import kepler_gzyx

### disabled and imported directly in the needed fct
## import geodezyx.reffram.sofa18 as sofa
log = logging.getLogger("geodezyx")

##########  END IMPORT  ##########


def compar_orbit(
    data_inp_1,
    data_inp_2,
    step_data=900,
    sats_used_list=["G"],
    name1="",
    name2="",
    use_name_1_2_for_table_name=False,
    rtn_output=True,
    convert_ecef_eci=True,
    clean_null_values=True,
    conv_coef=10**3,
    return_sat_null=False,
):
    """
    Compares 2 GNSS orbits files (SP3), and gives a summary plot and a
    statistics table

    Parameters
    ----------
    data_inp_1, data_inp_2 : str or Pandas DataFrame
        contains the orbits or path (string) to the sp3

    step_data : int
        per default data sampling

    sats_used_list : list of str
        used constellation or satellite : G E r c ... E01 , G02 ...
        Individuals satellites are prioritary on whole constellations
        e.g. ['G',"E04"]

    rtn_output : bool
        select output, Radial Transverse Normal or XYZ

    convert_ecef_eci : bool
        convert sp3 ECEF => ECI (Terrestrial => Celestrial)
        must be True in operational to avoid artifacts.

    name1, name2 : str (optionals)
        optional custom names for the 2 orbits

    use_name_1_2_for_table_name : bool
        False : use name 1 and 2 for table name, use datafile instead

    clean_null_values : bool or str
        if True or "all" remove sat position in all X,Y,Z values
        are null (0.000000)
        if "any", remove sat position if X or Y or Z is null
        if False, keep everything

    conv_coef : int
        conversion coefficient, km to m 10**3, km to mm 10**6

    Returns
    -------
    diff_sat_all : Pandas DataFrame
    contains differences b/w data_inp_1 & data_inp_2
    in Radial Transverse Normal OR XYZ frame

        Attributes of diff_sat_all :
            diff_sat_all.name : title of the table

    Note
    ----
    clean_null_values if useful (and necessary) only if
    convert_ECEF_ECI = False
    if convert_ECEF_ECI = True, the cleaning will be done by
    a side effect trick : the convertion ECEF => ECI will generate NaN
    for a zero-valued position
    But, nevertheless, activating  clean_null_values = True is better
    This Note is in fact usefull if you want to see bad positions on a plot
    => Then convert_ECEF_ECI = False and clean_null_values = False

    References
    ----------
    "Coordinate Systems", ASEN 3200 1/24/06 George h. Born

    """

    # selection of both used Constellations AND satellites
    sys_used_list = []
    prn_used_list = []
    for e in sats_used_list:
        if len(e) == 1:  ## it is a constellation
            sys_used_list.append(e)
        elif len(e) == 3:  ## it is a satellite
            prn_used_list.append(e)
            if not e[0] in sys_used_list:
                sys_used_list.append(e[0])

    # Read the files or DataFrames
    # metadata attributes are not copied
    # Thus, manual copy ...
    # (Dirty way, should be impoved without so many lines ...)
    if type(data_inp_1) is str:
        d1orig = files_rw.read_sp3(data_inp_1, epoch_as_pd_index=True)
    else:
        d1orig = data_inp_1.copy(True)
        try:
            d1orig.name = data_inp_1.name
        except:
            d1orig.name = "no_name"
        try:
            d1orig.path = data_inp_1.path
        except:
            d1orig.path = "no_path"
        try:
            d1orig.filename = data_inp_1.filename
        except:
            d1orig.filename = "no_filename"

    if type(data_inp_2) is str:
        d2orig = files_rw.read_sp3(data_inp_2, epoch_as_pd_index=True)
    else:
        d2orig = data_inp_2.copy(True)
        try:
            d2orig.name = data_inp_2.name
        except:
            d2orig.name = "no_name"
        try:
            d2orig.path = data_inp_2.path
        except:
            d2orig.path = "no_path"
        try:
            d2orig.filename = data_inp_2.filename
        except:
            d2orig.filename = "no_filename"

    #### NB : It has been decided with GM that the index of a SP3 dataframe
    ####      will be integers, not epoch datetime anymore
    ####      BUT here, for legacy reasons, the index has to be datetime

    if isinstance(d1orig.index[0], (int, np.integer)):
        d1orig.set_index("epoch", inplace=True)

    if isinstance(d2orig.index[0], (int, np.integer)):
        d2orig.set_index("epoch", inplace=True)

    diff_sat_stk = []

    # This block is for removing null values
    if clean_null_values:
        if clean_null_values == "all":
            all_or_any = np.all
        elif clean_null_values == "any":
            all_or_any = np.any
        else:
            all_or_any = np.all

        xyz_lst = ["x", "y", "z"]

        d1_null_bool = all_or_any(np.isclose(d1orig[xyz_lst], 0.0), axis=1)
        d2_null_bool = all_or_any(np.isclose(d2orig[xyz_lst], 0.0), axis=1)

        d1 = d1orig[np.logical_not(d1_null_bool)]
        d2 = d2orig[np.logical_not(d2_null_bool)]

        if np.any(d1_null_bool) or np.any(d2_null_bool):
            sat_nul = utils.join_improved(" ", *list(set(d1orig[d1_null_bool]["prn"])))
            log.warning("Null values contained in SP3 files : ")
            log.warning(
                "f1: %s %s",
                np.sum(d1_null_bool),
                utils.join_improved(" ", *list(set(d1orig[d1_null_bool]["prn"]))),
            )
            log.warning(
                "f2: %s %s",
                np.sum(d2_null_bool),
                utils.join_improved(" ", *list(set(d2orig[d2_null_bool]["prn"]))),
            )
        else:
            sat_nul = []

    else:
        d1 = d1orig.copy()
        d2 = d2orig.copy()

    d1sys_grp = d1.groupby("sys")
    d2sys_grp = d2.groupby("sys")

    for sysuse in sys_used_list:
        d1sys = d1sys_grp.get_group(sysuse)  # d1sys = d1[d1['sys'] == sysuse]
        d2sys = d2sys_grp.get_group(sysuse)  # d2sys = d2[d2['sys'] == sysuse]

        # checking if the data correspond to the step
        bool_step1 = np.mod((d1sys.index - np.min(d1.index)).seconds, step_data) == 0
        bool_step2 = np.mod((d2sys.index - np.min(d2.index)).seconds, step_data) == 0

        d1win = d1sys[bool_step1]
        d2win = d2sys[bool_step2]

        # find common sats and common epochs
        prni_set = sorted(list(set(d1win["prni"]).intersection(set(d2win["prni"]))))
        epoc_set = sorted(list(set(d1win.index).intersection(set(d2win.index))))

        # if special selection of sats, then apply it
        # (it is late and this selection is incredibely complicated ...)
        if np.any([True if sysuse in e else False for e in prn_used_list]):
            # first find the selected sats for the good constellation
            prni_used_select_list = [int(e[1:]) for e in prn_used_list if sysuse in e]
            # and apply it
            prni_set = sorted(
                list(set(prni_set).intersection(set(prni_used_select_list)))
            )

        d1win_prni_grp = d1win.groupby("prni")
        d2win_prni_grp = d2win.groupby("prni")

        for prni in prni_set:
            # First research : find corresponding epoch for the SV
            # this one is sufficent if there is no gaps (e.g. with 0.00000) i.e.
            # same nb of obs in the 2 files
            # NB : .reindex() is smart, it fills the DataFrame
            # with NaN
            try:
                d1prni_orig = d1win_prni_grp.get_group(prni).reindex(epoc_set)
                # d1win[d1win['prni'] == prni].reindex(epoc_set)
                d2prni_orig = d2win_prni_grp.get_group(prni).reindex(epoc_set)
                # d2win[d2win['prni'] == prni].reindex(epoc_set)
            except Exception as exce:
                log.info("ERR : Unable to re-index with an unique epoch")
                log.info(
                    "      are you sure there is no multiple-defined epochs for the same sat ?"
                )
                log.info(
                    "      it happens e.g. when multiple ACs are in the same DataFrame "
                )
                log.info(
                    "TIP : Filter the input Dataframe before calling this fct with"
                )
                log.info("      DF = DF[DF['AC'] == 'gbm']")

                dtmp1 = d1orig[d1orig["prni"] == prni]
                dtmp2 = d2orig[d2orig["prni"] == prni]

                dupli1 = np.sum(dtmp1.duplicated(["epoch", "prn"]))
                dupli2 = np.sum(dtmp2.duplicated(["epoch", "prn"]))

                log.info(
                    "FWIW: duplicated epoch/sat in DF1 & DF2: %s %s", dupli1, dupli2
                )

                raise exce

            # Second research, it is a security in case of gap
            # This step is useless, because .reindex() will fill the DataFrame
            # with NaN
            if len(d1prni_orig) != len(d2prni_orig):
                log.info(
                    "different epochs nbr for SV %s %s %s",
                    prni,
                    len(d1prni_orig),
                    len(d2prni_orig),
                )
                epoc_prni_set = sorted(
                    list(set(d1prni_orig.index).intersection(set(d2prni_orig.index)))
                )
                d1prni = d1prni_orig.loc[epoc_prni_set]
                d2prni = d2prni_orig.loc[epoc_prni_set]
            else:
                d1prni = d1prni_orig
                d2prni = d2prni_orig

            p1 = d1prni[["x", "y", "z"]]
            p2 = d2prni[["x", "y", "z"]]

            # Start ECEF => ECI
            if convert_ecef_eci:
                # Backup because the columns xyz will be reaffected
                # D1sv_bkp = d1prni.copy()
                # D2sv_bkp = d2prni.copy()

                p1b = conv.ecef2eci(
                    np.array(p1),
                    conv.dt_gpstime2dt_utc(p1.index.to_pydatetime(), out_array=True),
                )
                p2b = conv.ecef2eci(
                    np.array(p2),
                    conv.dt_gpstime2dt_utc(p2.index.to_pydatetime(), out_array=True),
                )

                d1prni[["x", "y", "z"]] = p1b
                d2prni[["x", "y", "z"]] = p2b

                p1 = d1prni[["x", "y", "z"]]
                p2 = d2prni[["x", "y", "z"]]
            # End ECEF => ECI

            if not rtn_output:
                # Compatible with the documentation +
                # empirically tested with OV software
                # it is  p1 - p2 (and not p2 - p1)
                delta_p = p1 - p2

                diff_sat = delta_p.copy()
                diff_sat.columns = ["dx", "dy", "dz"]

            else:
                rnorm = np.linalg.norm(p1, axis=1)

                from geodezyx.utils_xtra import pandas_utils

                vx = pandas_utils.diff_pandas(d1prni, "x", use_np_diff=True)
                vy = pandas_utils.diff_pandas(d1prni, "y", use_np_diff=True)
                vz = pandas_utils.diff_pandas(d1prni, "z", use_np_diff=True)

                v = pd.concat((vx, vy, vz), axis=1)
                v.columns = ["vx", "vy", "vz"]

                r = p1.divide(rnorm, axis=0)
                r.columns = ["xnorm", "ynorm", "znorm"]

                h = pd.DataFrame(np.cross(r, v), columns=["hx", "hy", "hz"])
                hnorm = np.linalg.norm(h, axis=1)

                c = h.divide(hnorm, axis=0)
                c.columns = ["hxnorm", "hynorm", "hznorm"]

                i = pd.DataFrame(np.cross(c, r), columns=["ix", "iy", "iz"])

                r_ar = np.array(r)
                i_ar = np.array(i)
                c_ar = np.array(c)

                # r_ar[1]
                beta = np.stack((r_ar, i_ar, c_ar), axis=1)

                # Compatible with the documentation +
                # empirically tested with OV software
                # it is  p1 - p2 (and not p2 - p1)
                delta_p = p1 - p2

                # Final determination
                astk = []

                for i in range(len(delta_p)):
                    a = np.dot(beta[i, :, :], np.array(delta_p)[i])
                    astk.append(a)

                diff_sat = pd.DataFrame(
                    np.vstack(astk), index=p1.index, columns=["dr", "dt", "dn"]
                )

            diff_sat = diff_sat * conv_coef  # metrer conversion

            diff_sat["sys"] = [sysuse] * len(diff_sat.index)
            diff_sat["prni"] = [prni] * len(diff_sat.index)
            diff_sat["prn"] = [sysuse + str(prni).zfill(2)] * len(diff_sat.index)

            diff_sat_stk.append(diff_sat)

    diff_sat_all = pd.concat(diff_sat_stk)
    date = diff_sat.index[0]

    # Attribute definition
    if rtn_output:
        diff_sat_all.frame_type = "RTN"

        # Pandas donesn't manage well iterable as attribute
        # So, it is separated
        diff_sat_all.frame_col_name1 = "dr"
        diff_sat_all.frame_col_name2 = "dt"
        diff_sat_all.frame_col_name3 = "dn"

    else:
        # Pandas donesn't manage well iterable as attribute
        # So, it is separated
        diff_sat_all.frame_col_name1 = "dx"
        diff_sat_all.frame_col_name2 = "dy"
        diff_sat_all.frame_col_name3 = "dz"

        if convert_ecef_eci:
            diff_sat_all.frame_type = "ECI"
        else:
            diff_sat_all.frame_type = "ECEF"

    # Name definitions
    if name1:
        diff_sat_all.name1 = name1
    else:
        diff_sat_all.name1 = d1orig.name

    if name2:
        diff_sat_all.name2 = name2
    else:
        diff_sat_all.name2 = d2orig.name

    diff_sat_all.filename1 = d1orig.filename
    diff_sat_all.filename2 = d2orig.filename

    diff_sat_all.path1 = d1orig.path
    diff_sat_all.path2 = d2orig.path

    diff_sat_all.name = " ".join(
        (
            "Orbits comparison (" + diff_sat_all.frame_type + ") b/w",
            diff_sat_all.name1,
            "(ref.) and",
            diff_sat_all.name2,
            ",",
            date.strftime("%Y-%m-%d"),
            ", doy",
            str(conv.dt2doy(date)),
        )
    )

    if return_sat_null:
        return diff_sat_all, sat_nul
    else:
        return diff_sat_all


def compar_orbit_plot(
    diff_sat_all_df_in,
    save_plot=False,
    save_plot_dir="",
    save_plot_name="auto",
    save_plot_name_suffix=None,
    save_plot_ext=(".pdf", ".png", ".svg"),
    yaxis_limit=None,
    yaxis_label_unit="m",
):
    """
    General description

    Parameters
    ----------
    diff_sat_all_df_in : DataFrame
        a DataFrame produced by compar_orbit

    yaxis_limit : 3-tuple iterable or 2-element tuple
        force the y axis limits. must look like
        [(ymin_r,ymax_r),(ymin_t,ymax_t),(ymin_n,ymax_n)]
        to control all the axis independely
        OR
        (ymin,ymax)
        to set all th axis at the same limits

    Returns
    -------
    the Figure and the 3 Axes if no save is asked
    export path (str) if save is asked
    but plot a plot anyway
    """

    fig, [axr, axt, axn] = plt.subplots(3, 1, sharex="all")

    satdispo = natsort.natsorted(list(set(diff_sat_all_df_in["prn"])))

    symb_stk = []

    cm = plt.get_cmap("viridis")
    num_colors = len(satdispo)
    colors = [cm(1.0 * i / num_colors) for i in range(num_colors)]

    # Pandas donesn't manage well iterable as attribute
    # So, it is separated
    try:
        col_name0 = diff_sat_all_df_in.frame_col_name1
        col_name1 = diff_sat_all_df_in.frame_col_name2
        col_name2 = diff_sat_all_df_in.frame_col_name3
    except:
        col_name0 = diff_sat_all_df_in.columns[0]
        col_name1 = diff_sat_all_df_in.columns[1]
        col_name2 = diff_sat_all_df_in.columns[2]

    for satuse, color in zip(satdispo, colors):
        diffuse = diff_sat_all_df_in[diff_sat_all_df_in["prn"] == satuse]

        time = diffuse.index
        r = diffuse[col_name0]
        t = diffuse[col_name1]
        n = diffuse[col_name2]

        # fig.fmt_xdata = mdates.DateFormatter('%Y-%m-%d')

        symb = axr.plot(time, r, label=satuse, c=color)
        axt.plot(time, t, label=satuse, c=color)
        axn.plot(time, n, label=satuse, c=color)

        symb_stk.append(symb[0])

        fig.autofmt_xdate()

    ylabuni = " (" + yaxis_label_unit + ")"

    if diff_sat_all_df_in.frame_type == "RTN":
        axr.set_ylabel("Radial diff." + ylabuni)
        axt.set_ylabel("Transverse diff." + ylabuni)
        axn.set_ylabel("Normal diff." + ylabuni)

    else:
        axr.set_ylabel(diff_sat_all_df_in.frame_type + " X diff." + ylabuni)
        axt.set_ylabel(diff_sat_all_df_in.frame_type + " Y diff." + ylabuni)
        axn.set_ylabel(diff_sat_all_df_in.frame_type + " Z diff." + ylabuni)

    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    axr.yaxis.set_major_formatter(y_formatter)
    axt.yaxis.set_major_formatter(y_formatter)
    axn.yaxis.set_major_formatter(y_formatter)

    if yaxis_limit and len(yaxis_limit) == 3:  ### indep. axis limit
        axr.set_ylim(yaxis_limit[0])
        axt.set_ylim(yaxis_limit[1])
        axn.set_ylim(yaxis_limit[2])
    elif yaxis_limit and len(yaxis_limit) == 2:
        axr.set_ylim(yaxis_limit)
        axt.set_ylim(yaxis_limit)
        axn.set_ylim(yaxis_limit)
    else:
        pass

    import matplotlib.dates as mdates

    fig.fmt_xdata = mdates.DateFormatter("%Y-%m-%d")

    lgd = fig.legend(
        tuple(symb_stk), satdispo, loc="lower center", ncol=8, columnspacing=1
    )

    fig.set_size_inches(8.27, 11.69)
    plt.suptitle(diff_sat_all_df_in.name)
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)
    plt.subplots_adjust(bottom=0.15)

    if save_plot:
        if save_plot_name == "auto":
            save_plot_name = "_".join(
                (
                    diff_sat_all_df_in.name1,
                    diff_sat_all_df_in.name2,
                    diff_sat_all_df_in.index.min().strftime("%Y-%m-%d"),
                )
            )

        if save_plot_name_suffix:
            save_plot_name = save_plot_name + "_" + save_plot_name_suffix

        for ext in save_plot_ext:
            save_plot_path = os.path.join(save_plot_dir, save_plot_name)
            plt.savefig(save_plot_path + ext)
            return_val = save_plot_path

    else:
        return_val = fig, (axr, axt, axn)

    return return_val


def compar_orbit_table(diff_sat_all_df_in, rms_style="natural", light_tab=False):
    """
    Generate a table with statistical indicators for an orbit comparison
    (RMS mean, standard dev, ...)

    Parameters
    ----------
    diff_sat_all_df_in : Pandas DataFrame
        a DataFrame produced by compar_orbit

    rms_style : str
        'natural': use the natural definition of the RMS
        'GRGS': RMS calc based on the GRGS definition of the RMS (OV help)
        is actually the standard deviation
        'kouba': RMS as defined in Kouba et al. 1994, p75
        using the degree of freedom (3*Nobs - 7)

    light_tab : bool
        produce a table with only RMS, with min/max/arithmetic instead

    Returns
    -------
    Compar_tab_out : DataFrame
        Statistical results of the comparison

    Note
    ----
    you can pretty print the output DataFrame using tabular module
    here is a template:

    >>> from tabulate import tabulate
    >>> print(tabulate(compar_table,headers="keys",floatfmt=".4f"))
    """

    sat_list = utils.uniq_and_sort(diff_sat_all_df_in["prn"])

    # Pandas donesn't manage well iterable as attribute
    # So, it is separated
    try:
        col_name0 = diff_sat_all_df_in.frame_col_name1
        col_name1 = diff_sat_all_df_in.frame_col_name2
        col_name2 = diff_sat_all_df_in.frame_col_name3
    except:
        col_name0 = diff_sat_all_df_in.columns[0]
        col_name1 = diff_sat_all_df_in.columns[1]
        col_name2 = diff_sat_all_df_in.columns[2]

    rms_stk = []

    for sat in sat_list:
        diffwork = utils.df_sel_val_in_col(diff_sat_all_df_in, "prn", sat)

        if rms_style == "natural":
            rms_a = stats.rms_mean(diffwork[col_name0])
            rms_b = stats.rms_mean(diffwork[col_name1])
            rms_c = stats.rms_mean(diffwork[col_name2])
        elif rms_style == "GRGS":
            rms_a = stats.rms_mean(diffwork[col_name0] - diffwork[col_name0].mean())
            rms_b = stats.rms_mean(diffwork[col_name1] - diffwork[col_name1].mean())
            rms_c = stats.rms_mean(diffwork[col_name2] - diffwork[col_name2].mean())
        elif rms_style == "kouba":
            rms_a = stats.rms_mean_kouba(diffwork[col_name0])
            rms_b = stats.rms_mean_kouba(diffwork[col_name1])
            rms_c = stats.rms_mean_kouba(diffwork[col_name2])

        rms_3d = np.sqrt(rms_a**2 + rms_b**2 + rms_c**2)

        min_a = diffwork[col_name0].min()
        min_b = diffwork[col_name1].min()
        min_c = diffwork[col_name2].min()

        max_a = diffwork[col_name0].max()
        max_b = diffwork[col_name1].max()
        max_c = diffwork[col_name2].max()

        mean_a = diffwork[col_name0].mean()
        mean_b = diffwork[col_name1].mean()
        mean_c = diffwork[col_name2].mean()

        if light_tab:
            rms_stk.append([rms_a, rms_b, rms_c, rms_3d])
        else:
            rms_stk.append(
                [
                    rms_a,
                    rms_b,
                    rms_c,
                    rms_3d,
                    min_a,
                    max_a,
                    mean_a,
                    min_b,
                    max_b,
                    mean_b,
                    min_c,
                    max_c,
                    mean_c,
                ]
            )

    #################################
    # ALL SATS
    if rms_style == "natural":
        rms_a = stats.rms_mean(diff_sat_all_df_in[col_name0])
        rms_b = stats.rms_mean(diff_sat_all_df_in[col_name1])
        rms_c = stats.rms_mean(diff_sat_all_df_in[col_name2])
        rms_3d = np.sqrt(rms_a**2 + rms_b**2 + rms_c**2)
    elif rms_style == "GRGS":
        rms_a = stats.rms_mean(
            diff_sat_all_df_in[col_name0] - diff_sat_all_df_in[col_name0].mean()
        )
        rms_b = stats.rms_mean(
            diff_sat_all_df_in[col_name1] - diff_sat_all_df_in[col_name1].mean()
        )
        rms_c = stats.rms_mean(
            diff_sat_all_df_in[col_name2] - diff_sat_all_df_in[col_name2].mean()
        )
        rms_3d = np.sqrt(rms_a**2 + rms_b**2 + rms_c**2)
    elif rms_style == "kouba":
        rms_a = stats.rms_mean_kouba(diff_sat_all_df_in[col_name0])
        rms_b = stats.rms_mean_kouba(diff_sat_all_df_in[col_name1])
        rms_c = stats.rms_mean_kouba(diff_sat_all_df_in[col_name2])
        rms_3d = np.sqrt(rms_a**2 + rms_b**2 + rms_c**2)

    min_a = diff_sat_all_df_in[col_name0].min()
    min_b = diff_sat_all_df_in[col_name1].min()
    min_c = diff_sat_all_df_in[col_name2].min()

    max_a = diff_sat_all_df_in[col_name0].max()
    max_b = diff_sat_all_df_in[col_name1].max()
    max_c = diff_sat_all_df_in[col_name2].max()

    mean_a = diff_sat_all_df_in[col_name0].mean()
    mean_b = diff_sat_all_df_in[col_name1].mean()
    mean_c = diff_sat_all_df_in[col_name2].mean()

    if light_tab:
        rms_stk.append([rms_a, rms_b, rms_c, rms_3d])
    else:
        rms_stk.append(
            [
                rms_a,
                rms_b,
                rms_c,
                rms_3d,
                min_a,
                max_a,
                mean_a,
                min_b,
                max_b,
                mean_b,
                min_c,
                max_c,
                mean_c,
            ]
        )

        # ALL SATS
    #################################

    if diff_sat_all_df_in.frame_type == "RTN":
        if light_tab:
            cols_nam = ["rmsR", "rmsT", "rmsN", "rms3D"]
        else:
            cols_nam = [
                "rmsR",
                "rmsT",
                "rmsN",
                "rms3D",
                "minR",
                "maxR",
                "meanR",
                "minT",
                "maxT",
                "meanT",
                "minN",
                "maxN",
                "meanN",
            ]

    else:
        if light_tab:
            cols_nam = ["rmsX", "rmsY", "rmsZ", "rms3D"]
        else:
            cols_nam = [
                "rmsX",
                "rmsY",
                "rmsZ",
                "rms3D",
                "minX",
                "maxX",
                "meanX",
                "minY",
                "maxY",
                "meanY",
                "minZ",
                "maxZ",
                "meanZ",
            ]

    compar_tab_out = pd.DataFrame(rms_stk, index=sat_list + ["ALL"], columns=cols_nam)

    return compar_tab_out


def compar_orbit_frontend(data_df1, data_df2, ac1, ac2, sats_used_list=["G"]):
    K = compar_orbit(
        data_df1[data_df1["ac"] == ac1],
        data_df2[data_df2["ac"] == ac2],
        sats_used_list=sats_used_list,
    )
    compar_orbit_plot(K)
    return K


def compar_clock(df_clk_inp_1, df_clk_inp_2, col_name="name", bias_col_name="bias"):
    """
    Compares 2 GNSS clock bias DataFrames (from .clk), to a
    statistics table (with compar_clock_table)


    Parameters
    ----------
    df_clk_inp_1 & df_clk_inp_2 : DataFrame
        Clock DataFrame provided by files_rw.read_clk()

    Returns
    -------
    df_clk_diff : DataFrame
        Clock bias difference DataFrame
    """
    df1idx = df_clk_inp_1.set_index([col_name, "epoch"])
    df1idx.sort_index(inplace=True)

    df2idx = df_clk_inp_2.set_index([col_name, "epoch"])
    df2idx.sort_index(inplace=True)

    i1 = df1idx.index
    i2 = df2idx.index

    iinter = i1.intersection(i2)
    iinter = iinter.sort_values()

    df_diff_bias = df1idx.loc[iinter][bias_col_name] - df2idx.loc[iinter][bias_col_name]

    df_clk_diff = df1idx.loc[iinter].copy()
    df_clk_diff[bias_col_name] = df_diff_bias
    if "ac" in df_clk_diff.columns:
        df_clk_diff.drop("ac", axis=1, inplace=True)
    else:
        df_clk_diff.drop("AC", axis=1, inplace=True)
    df_clk_diff.rename({bias_col_name: bias_col_name + "_diff"}, inplace=True, axis=1)

    # Name definitions
    if "ac" in df_clk_inp_1.columns:
        df_clk_diff["name1"] = df_clk_inp_1.ac.values[0]
    if "AC" in df_clk_inp_1.columns:
        df_clk_diff["name1"] = df_clk_inp_1.AC.values[0]

    if "ac" in df_clk_inp_1.columns:
        df_clk_diff["name2"] = df_clk_inp_2.ac.values[0]
    if "AC" in df_clk_inp_1.columns:
        df_clk_diff["name2"] = df_clk_inp_2.AC.values[0]

    return df_clk_diff


def compar_clock_table(df_clk_diff_in, col_name="name", bias_Col_name="bias_diff"):
    """
    Generate a table with statistical indicators for a clock comparison
    (RMS mean, standard dev, ...)

    Parameters
    ----------
    df_clk_diff_in : DataFrame
        Clock bias difference DataFrame (from compar_clock)

    Returns
    -------
    DFcompar_out : DataFrame
        Statistical results of the comparison.

    """

    df_diff_grp = df_clk_diff_in.groupby(col_name)[bias_Col_name]

    smin = df_diff_grp.min().rename("min", inplace=True)
    smax = df_diff_grp.max().rename("max", inplace=True)
    smean = df_diff_grp.mean().rename("mean", inplace=True)
    sstd = df_diff_grp.std().rename("std", inplace=True)
    srms = df_diff_grp.apply(stats.rms_mean).rename("rms", inplace=True)

    df_compar_out = pd.concat([smin, smax, smean, sstd, srms], axis=1)
    df_compar_out.reset_index()

    return df_compar_out


def compar_clk_plot(
    diff_sat_all_df_in,
    save_plot=False,
    save_plot_dir="",
    save_plot_name="auto",
    save_plot_name_suffix=None,
    save_plot_ext=(".pdf", ".png", ".svg"),
    yaxis_limit=None,
    yaxis_label_unit="psec",
    col_name="name",
    bias_Col_name="bias",
):
    """
    General description

    Parameters
    ----------
    diff_sat_all_df_in: DataFrame
        a DataFrame produced by compar_clk

    yaxis_limit: 3-tuple iterable or 2-element tuple
        force the y axis limits. must look like
        (ymin,ymax) to set all th axis at the same limits

    col_name: Normally the name of the column with the sat names
    bias_Col_name: The column with the clk values
    Default: 'bias'

    Returns
    -------
    the Figure and the 3 Axes if no save is asked
    export path (str) if save is asked
    but plot a plot anyway
    """

    fig, axr = plt.subplots(1, 1, sharex="all")
    diff_sat_all_df_in = diff_sat_all_df_in.reset_index()
    satdispo = natsort.natsorted(list(set(diff_sat_all_df_in[col_name])))
    # satdispo = natsort.natsorted(list(set(diff_sat_all_df_in['sat'])))

    symb_stk = []

    cm = plt.get_cmap("viridis")
    num_colors = len(satdispo)
    colors = [cm(1.0 * i / num_colors) for i in range(num_colors)]

    date = conv.numpy_dt2dt(diff_sat_all_df_in.epoch.values[0])
    diff_sat_all_df_in.name = " ".join(
        (
            "Clock comparison  b/w",
            diff_sat_all_df_in.name1.values[0],
            "(ref.) and",
            diff_sat_all_df_in.name2.values[0],
            ",",
            date.strftime("%Y-%m-%d"),
        )
    )
    # Pandas donesn't manage well iterable as attribute
    # So, it is separated
    try:
        col_name0 = diff_sat_all_df_in.frame_col_name1
        col_name1 = diff_sat_all_df_in.frame_col_name2
        col_name2 = diff_sat_all_df_in.frame_col_name3
    except:
        col_name0 = diff_sat_all_df_in.columns[0]
        col_name1 = diff_sat_all_df_in.columns[1]
        col_name2 = diff_sat_all_df_in.columns[2]

    for satuse, color in zip(satdispo, colors):
        diffuse = diff_sat_all_df_in[diff_sat_all_df_in[col_name] == satuse]

        time = diffuse.epoch
        r = diffuse[bias_Col_name + "_diff"] * 10**12

        # fig.fmt_xdata = mdates.DateFormatter('%Y-%m-%d')

        Symb = axr.plot(time, r, label=satuse, c=color)

        symb_stk.append(Symb[0])

        fig.autofmt_xdate()

    ylabuni = " (" + yaxis_label_unit + ")"

    axr.set_ylabel("Bias Diff." + ylabuni)

    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    axr.yaxis.set_major_formatter(y_formatter)

    if yaxis_limit and len(yaxis_limit) == 3:  ### indep. axis limit
        axr.set_ylim(yaxis_limit[0])

    elif yaxis_limit and len(yaxis_limit) == 2:
        axr.set_ylim(yaxis_limit)

    else:
        pass

    import matplotlib.dates as mdates

    fig.fmt_xdata = mdates.DateFormatter("%Y-%m-%d")

    lgd = fig.legend(
        tuple(symb_stk), satdispo, loc="lower center", ncol=8, columnspacing=1
    )

    fig.set_size_inches(8.27, 11.69)
    plt.suptitle(diff_sat_all_df_in.name)
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)
    plt.subplots_adjust(bottom=0.15)

    if save_plot:
        if save_plot_name == "auto":
            save_plot_name = (
                diff_sat_all_df_in.name1.values[0]
                + "_"
                + diff_sat_all_df_in.name2.values[0]
                + "_"
                + date.strftime("%Y-%m-%d")
            )

        if save_plot_name_suffix:
            save_plot_name = save_plot_name + "_" + save_plot_name_suffix

        for ext in save_plot_ext:
            save_plot_path = os.path.join(save_plot_dir, save_plot_name)
            plt.savefig(save_plot_path + ext)
            return_val = save_plot_path

    else:
        return_val = fig, (axr)

    return return_val


def compar_sinex(
    snx1,
    snx2,
    stat_select=None,
    invert_select=False,
    out_means_summary=True,
    out_meta=True,
    out_dataframe=True,
    manu_wwwwd=None,
):

    if type(snx1) is str:
        week1 = utils.split_improved(os.path.basename(snx1), "_", ".")[:]
        week2 = utils.split_improved(os.path.basename(snx2), "_", ".")[:]
        if week1 != week2:
            log.warning(
                "Dates of 2 input files are differents !!! It might be very bad !!! %s %s",
                week1,
                week2,
            )
        else:
            wwwwd = week1
        d1 = files_rw.read_sinex(snx1, True)
        d2 = files_rw.read_sinex(snx2, True)
    else:
        log.warning(
            "WARN : you are giving the SINEX input as a DataFrame, wwwwd has to be given manually using manu_wwwwd"
        )
        d1 = snx1
        d2 = snx2

    if manu_wwwwd:
        wwwwd = manu_wwwwd

    stat_common = set(d1["STAT"]).intersection(set(d2["STAT"]))

    if stat_select:

        stat_common_init = list(stat_common)

        if invert_select:
            select_fct = lambda x: not x
        else:
            select_fct = lambda x: x

        if type(stat_select) is str:
            stat_common = [
                sta
                for sta in stat_common_init
                if select_fct(re.search(stat_select, sta))
            ]
        elif utils.is_iterable(stat_select):
            stat_common = [
                sta for sta in stat_common_init if select_fct(sta in stat_select)
            ]
        else:
            log.warning("WARN : check type of stat_select")

    d1_common = (
        d1[d1["STAT"].isin(stat_common)].sort_values("STAT").reset_index(drop=True)
    )
    d2_common = (
        d2[d2["STAT"].isin(stat_common)].sort_values("STAT").reset_index(drop=True)
    )

    ddiff = pd.DataFrame()
    ddiff = ddiff.assign(STAT=d1_common["STAT"])

    #### XYZ Part
    for xyz in ("x", "y", "z"):

        dif = pd.to_numeric((d2_common[xyz] - d1_common[xyz]))

        ddiff = ddiff.assign(xyz=dif)
        ddiff = ddiff.rename(columns={"xyz": xyz})

    D3D = np.sqrt(
        (ddiff["x"] ** 2 + ddiff["y"] ** 2 + ddiff["z"] ** 2).astype("float64")
    )

    ddiff = ddiff.assign(d3D_xyz=D3D)

    ### ENU Part
    E, N, U = [], [], []
    enu_stk = []

    for (_, l1), (_, l2) in zip(d1_common.iterrows(), d2_common.iterrows()):
        enu = conv.xyz2enu(l1["x"], l1["y"], l1["z"], l2["x"], l2["y"], l2["z"])
        enu_stk.append(np.array(enu))

    if len(enu_stk) == 0:
        E, N, U = np.array([]), np.array([]), np.array([])
    else:
        ENU = np.hstack(enu_stk)
        E, N, U = ENU[0, :], ENU[1, :], ENU[2, :]

    D2D = np.sqrt((E**2 + N**2).astype("float64"))
    D3D = np.sqrt((E**2 + N**2 + U**2).astype("float64"))

    ddiff = ddiff.assign(e=E)
    ddiff = ddiff.assign(n=N)
    ddiff = ddiff.assign(u=U)
    ddiff = ddiff.assign(d2D_enu=D2D)
    ddiff = ddiff.assign(d3D_enu=D3D)

    #    E,N,U    = conv.xyz2enu((X,Y,Z,x0,y0,z0))
    #    E,N,U    = conv.xyz2enu((X,Y,Z,x0,y0,z0))

    if out_dataframe:
        out_meta = True

    if not out_means_summary:
        log.info("INFO : this is not used operationally and it can be improved")
        return ddiff
    else:
        output = []

        col_names = ("x", "y", "z", "d3D_xyz", "e", "n", "u", "d2D_enu", "d3D_enu")

        for xyz in col_names:
            output.append(stats.rms_mean(ddiff[xyz]))
        for xyz in col_names:
            output.append(np.nanmean(ddiff[xyz]))
        for xyz in col_names:
            output.append(np.nanstd(ddiff[xyz]))

        if out_meta:
            nstat = len(stat_common)
            week = int(wwwwd[:4])
            day = int(wwwwd[4:])
            output = [week, day, nstat] + output

        if not out_dataframe:
            return tuple(output)
        else:

            output_DF = pd.DataFrame(output).transpose()

            output_DF.columns = [
                "week",
                "dow",
                "nbstat",
                "x_rms",
                "y_rms",
                "z_rms",
                "d3D_xyz_rms",
                "e_rms",
                "n_rms",
                "u_rms",
                "d2D_enu_rms",
                "d3D_enu_rms",
                "x_ari",
                "y_ari",
                "z_ari",
                "d3D_xyz_ari",
                "e_ari",
                "n_ari",
                "u_ari",
                "d2D_enu_ari",
                "d3D_enu_ari",
                "x_ari",
                "y_std",
                "z_std",
                "d3D_xyz_std",
                "e_ari",
                "n_std",
                "u_std",
                "d2D_enu_std",
                "d3D_enu_std",
            ]

            return output_DF


#   ____       _     _ _     _____        _        ______
#  / __ \     | |   (_) |   |  __ \      | |      |  ____|
# | |  | |_ __| |__  _| |_  | |  | | __ _| |_ __ _| |__ _ __ __ _ _ __ ___   ___  ___
# | |  | | '__| '_ \| | __| | |  | |/ _` | __/ _` |  __| '__/ _` | '_ ` _ \ / _ \/ __|
# | |__| | |  | |_) | | |_  | |__| | (_| | || (_| | |  | | | (_| | | | | | |  __/\__ \
#  \____/|_|  |_.__/|_|\__| |_____/ \__,_|\__\__,_|_|  |_|  \__,_|_| |_| |_|\___||___/
#

### Orbit DataFrames


def orb_df_velocity_calc(df_orb_in, drop_nan=False):
    """
    Compute the velocity of satellites from a DFOrb dataframe
    (differentiate the position)

    Parameters
    ----------
    df_orb_in : Pandas DataFrame
        an Orbit DataFrame.
    drop_nan : bool, optional
        Remove the nan values
        (the first ones, since it is a numerical differentiation).
        The default is False.
        It is recommended to keep the NaN values, thus the output will keep
        same size and index as input

    Returns
    -------
    df_vel : Pandas DataFrame
        an Orbit DataFrame with velocities (vx, vy ,vz columns).

    """

    dfgrp = df_orb_in.groupby("prn")

    dfprn_stk = []
    for prn, dfprn in dfgrp:
        for coord in ["x", "y", "z"]:
            dcoord = dfprn[coord].diff()
            dtime = np.float64(dfprn["epoch"].diff()) * 10**-9
            # timedelta in ns per defalut
            dfprn["v" + coord] = dcoord / dtime
        dfprn_stk.append(dfprn)

    df_vel = pd.concat(dfprn_stk)

    df_vel.sort_index(inplace=True)

    if drop_nan:
        df_vel.dropna(inplace=True)

    return df_vel


def beta_sun_ra_dec(sun_dec, sun_ra, sat_i, sat_o_lan):
    """
    Compute beta angle based on Sun's right ascension and declination
    Angles are in radians

    Parameters
    ----------
    sun_dec : float
        Sun's declination.
    sun_ra : float
        Sun's right ascension.
    sat_i : float
        Satellite's inclination.
    sat_o_lan : float
        Satellite's Longitude of the ascending node.

    Returns
    -------
    beta : float
        beta angle (in radians).

    Source
    ------
    `Satellites: de Kepler au GPS`, Michel Capderou (2012)

    simpler formula here:
        https://www.fxsolver.com/browse/formulas/Beta+Angle


    Note
    ----
    you can use ``pyorbital.astronomy.sun_ra_dec()`` to compute ``sun_dec``
    and ``sun_ra``

    If so, Sun's right ascension and declination computation polynoms are
    based on: `Astronomical algorithms`, Jean Meeus (1st edition, 1991)
    """
    beta = np.arcsin(
        np.cos(sun_dec) * np.sin(sat_i) * np.sin(sat_o_lan - sun_ra)
        + np.sin(sun_dec) * np.cos(sat_i)
    )
    return beta


def beta_sun_eclip_long(sun_ecl_long, sat_o_lan, sat_i, earth_i):
    """
    Compute beta angle based on Sun's Ecliptic longitude
    Angles are in radians

    Parameters
    ----------
    sun_ecl_long : float
        Sun's Ecliptic longitude.
    sat_i : float
        Satellite's inclination.
    sat_o_lan : float
        Satellite's Longitude of the ascending node.
    earth_i : float
        Earth's inclination.

    Returns
    -------
    beta : float
        beta angle (in radians).

    Source
    ------
    `Calculation of the Eclipse Factor for Elliptical Satellite Orbits` NASA (1962)
    https://ntrs.nasa.gov/citations/19630000622
    `Computation of Eclipse Time for Low-Earth Orbiting Small Satellites` Sumanth R M (2019)
    https://commons.erau.edu/ijaaa/vol6/iss5/15/

    Simpler formula:
        https://en.wikipedia.org/wiki/Beta_angle

    Note
    ----
    you can use ``pyorbital.astronomy.sun_ecliptic_longitude()``
    to compute ``sun_ecl_long``

    If so, Sun's right ascension and declination computation polynoms are
    based on: `Astronomical algorithms`, Jean Meeus (1st edition, 1991)

    """
    p1 = np.cos(sun_ecl_long) * np.sin(sat_o_lan) * np.sin(sat_i)
    p2 = np.sin(sun_ecl_long) * np.cos(earth_i) * np.cos(sat_o_lan) * np.sin(sat_i)
    p3 = np.sin(sun_ecl_long) * np.sin(earth_i) * np.cos(sat_i)
    beta = np.arcsin(p1 - p2 + p3)
    return beta


def beta_angle_calc(
    DFOrb_in,
    calc_beta_sun_ra_dec=True,
    calc_beta_sun_eclip_long=True,
    beta_rad2deg=True,
):
    """
    Compute beta angle for GNSS satellite's orbits stored in an orbit
    DataFrame


    Parameters
    ----------
    DFOrb_in : Pandas DataFrame
        an Orbit DataFrame (ECEF frame).
    calc_beta_sun_ra_dec : bool, optional
        compute beta angle with Sun's right ascension and declination.
        The default is True.
    calc_beta_sun_eclip_long : bool, optional
        compute beta angle with Sun's Ecliptic longitude.
        The default is True.
    beta_rad2deg : bool, optional
        convert beta angle to degrees. The default is True.

    Returns
    -------
    df_out : Pandas DataFrame
        df_orb_in with a new 'beta' column.
    df_wrk : Pandas DataFrame
        intermediate values dataframe for debug.
        here, coordinates are in ECI frame
    """
    import pyorbital.astronomy

    #### convert in ECI
    df_eci = DFOrb_in.copy()
    df_eci[["x", "y", "z"]] = conv.ecef2eci(
        DFOrb_in[["x", "y", "z"]].values, DFOrb_in["epoch"].values
    )

    #### compute velocity
    df_wrk = orb_df_velocity_calc(df_eci, drop_nan=False)
    ##keep drop nan False, thus the DF will keep same size and index as input

    #### compute Kepler's parameters
    p = df_wrk[["x", "y", "z"]].values * 1000
    v = df_wrk[["vx", "vy", "vz"]].values * 1000

    kep_col = ["a", "ecc", "i", "o_peri", "o_lan", "m"]
    kep_params = kepler_gzyx.eci_2_kepler_elts(p, v, rad2deg=False)
    df_wrk[kep_col] = np.column_stack(kep_params)
    ######### COMPUTE BETA

    #### cosmetic changes
    if calc_beta_sun_ra_dec and calc_beta_sun_eclip_long:
        b1 = "1"
        b2 = "2"
    else:
        b1 = ""
        b2 = ""

    if beta_rad2deg:
        r2dfct = np.rad2deg
    else:
        r2dfct = lambda x: x

    ##### Beta computed based on sun declination / right ascension
    if calc_beta_sun_ra_dec:
        ##### sun_ra_dec output in RADIANS
        df_wrk[["sun_ra", "sun_dec"]] = np.column_stack(
            pyorbital.astronomy.sun_ra_dec(df_wrk["epoch"])
        )
        df_wrk["beta" + b1] = beta_sun_ra_dec(
            df_wrk["sun_dec"], df_wrk["sun_ra"], df_wrk["i"], df_wrk["o_lan"]
        ).apply(r2dfct)

    ##### Beta computed based on Ecliptic longitude of the sun
    if calc_beta_sun_eclip_long:
        ### sun_ecliptic_longitude output in RADIANS but NO MODULO !!!
        df_wrk["sun_ecl_long"] = np.mod(
            pyorbital.astronomy.sun_ecliptic_longitude(df_wrk["epoch"]), np.pi * 2
        )
        df_wrk["beta" + b2] = beta_sun_eclip_long(
            df_wrk["sun_ecl_long"], df_wrk["o_lan"], df_wrk["i"], np.deg2rad(23.45)
        ).apply(r2dfct)
    df_out = DFOrb_in.copy()
    df_out["beta"] = df_wrk["beta" + b1]

    return df_out, df_wrk


def orb_df_lagrange_interpolate(
    orb_df_inp, titrp, n=10, append_to_input_df=False, plot=False
):
    """
    High level function to interpolate an orbit DataFrame

    Parameters
    ----------
    orb_df_inp : DataFrame
        an Orbit DataFrame.
    titrp : iterable of datetime
        Epochs of the wished points.
    n : int, optional
        degree of the polynom. Better if even. The default is 10.
    append_to_input_df : bool, optional
        append the interpolated DF to the input DF. The default is False.
    plot : bool, optional
        Plot the values. For debug only. The default is False.

    Returns
    -------
    DForb_out : DataFrame
        Interpolated orbits.

    Tips
    ----
    Use conv.dt_range to generate the wished epochs range

    """
    df_orb_stk = []

    for sat, ac in itertools.product(
        orb_df_inp["prn"].unique(), orb_df_inp["ac"].unique()
    ):

        log.info("process %s %s", ac, sat)

        DForb_use = orb_df_inp[
            (orb_df_inp["prn"] == sat) & (orb_df_inp["ac"] == ac)
        ].copy()

        ### faster but anoying Future Waring
        # Tdata = np.array(DForb_use.epoch.dt.to_pydatetime())
        Tdata = conv.numpy_dt2dt(DForb_use.epoch.values)

        Xitrp = stats.lagrange_interpolate(Tdata, DForb_use["x"], titrp, n=n)
        Yitrp = stats.lagrange_interpolate(Tdata, DForb_use["y"], titrp, n=n)
        Zitrp = stats.lagrange_interpolate(Tdata, DForb_use["z"], titrp, n=n)

        Clkitrp = np.interp(
            conv.dt2posix(np.array(titrp)),
            conv.dt2posix(np.array(Tdata)),
            DForb_use["clk"].values,
        )

        # ClkDummy = np.array([999999.999999] * len(titrp))

        d = {"epoch": titrp, "x": Xitrp, "y": Yitrp, "z": Zitrp, "clk": Clkitrp}
        DForb_tmp = pd.DataFrame(d)

        ### sometihng else must be tested o give the annex val directly in the col of DForb_tmp
        DFannex_vals = DForb_use.drop(
            ["epoch", "x", "y", "z", "clk"], axis=1
        ).drop_duplicates()
        DFannex_vals = pd.concat(
            [DFannex_vals] * (len(titrp)), ignore_index=True, axis=0
        )

        DForb_tmp = pd.concat((DForb_tmp, DFannex_vals), axis=1)
        df_orb_stk.append(DForb_tmp)

        if plot:
            # plt.plot(Tdata,DForb_use.x,'o')
            # plt.plot(titrp,Xitrp,'.')
            ## GUS mod 220322
            fig, axr = plt.subplots(1, 1, sharex="all")
            Symb = axr.plot(Tdata, DForb_use.x, "o")
            Symb = axr.plot(titrp, Xitrp, ".")

    DForb_out = pd.concat(df_orb_stk)

    if append_to_input_df:
        DForb_out = pd.concat((orb_df_inp, DForb_out))

    DForb_out.reset_index(drop=True)
    DForb_out[["x", "y", "z", "clk"]] = DForb_out[["x", "y", "z", "clk"]].astype(float)
    return DForb_out


def orb_df_crf2trf(df_orb_inp, df_eop_inp, time_scale_inp="gps", inv_trf2crf=False):
    """
    Convert an Orbit DataFrame from Celetrial Reference Frame to
    Terrestrial Reference Frame.

    Requires EOP to work. Cf. note below.

    Parameters
    ----------
    df_orb_inp : DataFrame
        Input Orbit DataFrame in Celetrial Reference Frame.
    df_eop_inp : DataFrame
        EOP DataFrame  (C04 format).
    time_scale_inp : str, optional
        The time scale used in. manage 'utc', 'tai' and 'gps'.
        The default is "gps".
    inv_trf2crf : bool, optional
        Provide the inverse transformation TRF => CRF.
        The default is False.

    Returns
    -------
    DForb_out : DataFrame
        Output Orbit DataFrame in Terrestrial Reference Frame.
        (or Celestrial if inv_trf2crf is True)

    Note
    ----
    The EOP can be obtained from the IERS C04 products.
    e.g.
    https://datacenter.iers.org/data/latestVersion/224_EOP_C04_14.62-NOW.IAU2000A224.txt
    To get them as a Compatible DataFrame, use the function
    files_rw.read_eop_C04()
    """

    df_orb = df_orb_inp.copy()

    import geodezyx.reffram.sofa18 as sofa

    ### bring everything to UTC
    if time_scale_inp.lower() == "gps":
        df_orb["epoch_utc"] = conv.dt_gpstime2dt_utc(df_orb["epoch"])
    elif time_scale_inp.lower() == "tai":
        df_orb["epoch_utc"] = conv.dt_tai2dt_utc(df_orb["epoch"])
    elif time_scale_inp.lower() == "utc":
        df_orb["epoch_utc"] = df_orb["epoch"]
    ### TT and UT1 are not implemented (quite unlikely to have them as input)

    ### do the time scale's conversion
    df_orb["epoch_tai"] = conv.dt_utc2dt_tai(df_orb["epoch_utc"])
    df_orb["epoch_tt"] = conv.dt_tai2dt_tt(df_orb["epoch_tai"])
    df_orb["epoch_ut1"] = conv.dt_utc2dt_ut1_smart(df_orb["epoch_utc"], df_eop_inp)

    ### Do the EOP interpolation
    df_eop_intrp = eop_interpotate(df_eop_inp, df_orb["epoch_utc"])
    ### bring the EOP to radians
    Xeop = np.deg2rad(conv.arcsec2deg(df_eop_intrp["x"]))
    Yeop = np.deg2rad(conv.arcsec2deg(df_eop_intrp["y"]))

    TRFstk = []

    for tt, ut1, xeop, yeop, x, y, z in zip(
        df_orb["epoch_tt"],
        df_orb["epoch_ut1"],
        Xeop,
        Yeop,
        df_orb["x"],
        df_orb["y"],
        df_orb["z"],
    ):

        MatCRF22TRF = sofa.iau_c2t06a(
            2400000.5, conv.dt2mjd(tt), 2400000.5, conv.dt2mjd(ut1), xeop, yeop
        )
        if inv_trf2crf:
            MatCRF22TRF = np.linalg.inv(MatCRF22TRF)

        CRF = np.array([x, y, z])
        TRF = np.dot(MatCRF22TRF, CRF)

        TRFstk.append(TRF)

    ### Final stack and replacement
    TRFall = np.vstack(TRFstk)
    DForb_out = df_orb_inp.copy()
    DForb_out[["x", "y", "z"]] = TRFall

    return DForb_out


#### FCT DEF
def orb_df_reg_2_multidx(OrbDFin, index_order=["prn", "epoch"]):
    """
    From an regular Orbit DF generated by read_sp3(), set some columns
    (typically ["prn","epoch"]) as indexes
    The outputed DF is then a Multi-index DF
    """
    OrbDFwrk = OrbDFin.reset_index()
    OrbDFwrk = OrbDFwrk.sort_values(index_order)
    OrbDFwrk = OrbDFwrk.set_index(index_order, inplace=False)
    return OrbDFwrk


def orb_df_multidx_2_reg(OrbDFin, index_order=["prn", "epoch"]):
    """
    Convert a Multi-index formatted OrbitDF to his original form
    """
    OrbDFwrk = OrbDFin.reset_index()
    OrbDFwrk = OrbDFwrk.sort_values(index_order)
    OrbDFwrk["sys"] = OrbDFwrk["prn"].apply(lambda x: x[0])
    OrbDFwrk["prni"] = OrbDFwrk["prn"].apply(lambda x: int(x[1:]))
    return OrbDFwrk


def orb_df_common_epoch_finder(
    OrbDFa_in,
    orb_df_b_in,
    return_index=False,
    supplementary_sort=False,
    order=["prn", "epoch"],
    skip_reg2multidx_OrbDFa=False,
    skip_reg2multidx_OrbDFb=False,
):
    """
    This function finds common satellites and epochs in two Orbit DataFrames and outputs the corresponding Orbit DataFrames.

    Parameters
    ----------
    OrbDFa_in : DataFrame
        The first input Orbit DataFrame.
    orb_df_b_in : DataFrame
        The second input Orbit DataFrame.
    return_index : bool, optional
        If True, the function also returns the common index. Default is False.
    supplementary_sort : bool, optional
        If True, an additional sort is performed. This is useful for multi GNSS where the output DataFrame may not be well sorted. Default is False.
    order : list of str, optional
        The order of the index for the multi-index DataFrame. Default is ["prn","epoch"].
    skip_reg2multidx_OrbDFa : bool, optional
        If True, skips the conversion of the first input DataFrame to a multi-index DataFrame. Default is False.
        The inputs are assumed to be already in multi-index format to optimize execution speed.
        (For advanced use only)
    skip_reg2multidx_OrbDFb : bool, optional
        If True, skips the conversion of the second input DataFrame to a multi-index DataFrame. Default is False.
        The inputs are assumed to be already in multi-index format to optimize execution speed.
        (For advanced use only)

    Returns
    -------
    orb_df_a_out : DataFrame
        The first output Orbit DataFrame with common satellites and epochs.
    orb_df_b_out : DataFrame
        The second output Orbit DataFrame with common satellites and epochs.
    Iinter : Index, optional
        The common index. Only returned if return_index is True.

    Note
    ----
    designed for orbits/sp3 first with sat and epoch as order parmeter,
    but can be used also for instance for snx files with
    STAT and epoch as order parmeter
    """
    if not skip_reg2multidx_OrbDFa:
        orb_df_a = orb_df_reg_2_multidx(OrbDFa_in, index_order=order)
    else:
        orb_df_a = OrbDFa_in
    if not skip_reg2multidx_OrbDFb:
        orb_df_b = orb_df_reg_2_multidx(orb_df_b_in, index_order=order)
    else:
        orb_df_b = orb_df_b_in

    i1 = orb_df_a.index
    i2 = orb_df_b.index

    iinter = i1.intersection(i2)
    iinter = iinter.sort_values()

    orb_df_a_out = orb_df_a.loc[iinter]
    orb_df_b_out = orb_df_b.loc[iinter]

    if supplementary_sort:
        orb_df_a_out = orb_df_a_out.sort_values(order)
        orb_df_b_out = orb_df_b_out.sort_values(order)

    if len(orb_df_a_out) != len(orb_df_b_out):
        log.warning("len(Orb/ClkDFa_out) != len(Orb/ClkDFb_out)")
        log.warning("TIPS : ClkDFa_in and/or ClkDFb_in might contain duplicates")

    if return_index:
        return orb_df_a_out, orb_df_b_out, iinter
    else:
        return orb_df_a_out, orb_df_b_out


def orb_df_const_sv_columns_maker(orb_df_in, inplace=True):
    """
    (re)generate the const and sv columns from the sat one
    """
    if inplace:
        orb_df_in["sys"] = orb_df_in["prn"].str[0]
        orb_df_in["prni"] = orb_df_in["prn"].apply(lambda x: int(x[1:]))
        return None
    else:
        orb_df_out = orb_df_in.copy()
        orb_df_out["sys"] = orb_df_out["prn"].str[0]
        orb_df_out["prni"] = orb_df_out["prn"].apply(lambda x: int(x[1:]))
        return orb_df_out


#   _____ _            _      _____        _        ______
#  / ____| |          | |    |  __ \      | |      |  ____|
# | |    | | ___   ___| | __ | |  | | __ _| |_ __ _| |__ _ __ __ _ _ __ ___   ___  ___
# | |    | |/ _ \ / __| |/ / | |  | |/ _` | __/ _` |  __| '__/ _` | '_ ` _ \ / _ \/ __|
# | |____| | (_) | (__|   <  | |__| | (_| | || (_| | |  | | | (_| | | | | | |  __/\__ \
#  \_____|_|\___/ \___|_|\_\ |_____/ \__,_|\__\__,_|_|  |_|  \__,_|_| |_| |_|\___||___/

### Clock DataFrames


def clk_df_filter(
    clk_df_in,
    typ=("AS", "AR"),
    name=None,
    ac=None,
    epoch_strt=dt.datetime(1980, 1, 1),
    epoch_end=dt.datetime(2099, 1, 1),
    name_regex=False,
):
    """
    Filter the content of a Clock DataFrame

    Parameters
    ----------
    clk_df_in : DataFrame
        Input Clock DataFrame
        (a concatenation of DF generated by files_rw.read_clk.
    typ : iterable of str, optional
        List of the types of clocks: AS (satellite) or AR (receiver).
        The default is ("AS","AR").
    name : iterable of str, optional
        List of wished satellites/stations.
        Can be a regex (see also name_regex)
        The default is None.
    ac : iterable of str, optional
        List of wished ACs. The default is None.
    epoch_strt : datetime, optional
        Start epoch. The default is dt.datetime(1980,1,1).
    epoch_end : datetime, optional
        End epoch (not included). The default is dt.datetime(2099,1,1).
    name_regex : bool, optional
        the given names as 'name' arguments are regular expressions
        Some useful regex are given bellow
        The default is False

    Returns
    -------
    Clock DataFrame
        Output Clock DataFrame.

    Notes
    -----
    '^E[0-9]{2}': Galileo Satellites
    '^G[0-9]{2}': GPS Satellites

    """

    if type(clk_df_in) is str:
        clk_df_wrk = utils.pickle_loader(clk_df_in)
    else:
        clk_df_wrk = clk_df_in

    bool = np.ones(len(clk_df_wrk)).astype(bool)

    if typ:
        bool_tmp = clk_df_wrk.type.isin(typ)
        bool = bool & np.array(bool_tmp)

    if name:
        if not name_regex:  ### full name mode
            bool_tmp = clk_df_wrk.name.isin(name)
            bool = bool & np.array(bool_tmp)
        else:  ### REGEX mode
            bool_tmp = np.zeros(len(clk_df_wrk.name)).astype(bool)
            for rgx in name:
                nam_serie = clk_df_wrk.name
                bool_tmp = bool_tmp | np.array(nam_serie.str.contains(rgx))

            bool = bool & np.array(bool_tmp)

    if ac:
        bool_tmp = clk_df_wrk.ac.isin(ac)
        bool = bool & np.array(bool_tmp)

    ##epoch
    bool_tmp = (epoch_strt <= clk_df_wrk.epoch) & (clk_df_wrk.epoch < epoch_end)
    bool = bool & np.array(bool_tmp)

    return clk_df_wrk[bool]


def clk_df_filter2(
    clk_df_in,
    typ=("AS", "AR"),
    name=None,
    ac=None,
    epoch_strt=dt.datetime(1980, 1, 1),
    epoch_end=dt.datetime(2099, 1, 1),
    name_regex=False,
):
    """
    attempt for a faster version of clk_df_filter, but the original is faster
    """

    if type(clk_df_in) is str:
        clk_df_wrk = utils.pickle_loader(clk_df_in)
    else:
        clk_df_wrk = clk_df_in

    clkdf_stk = []

    for (ityp, iname, iac), clkdf_grp in clk_df_wrk.groupby(["type", "name", "ac"]):
        if typ:
            bool_typ = True if ityp in typ else False
        else:
            bool_typ = True

        if name:
            if not name_regex:
                bool_name = True if iname in name else False
            else:
                bool_name = any([re.search(n, iname) for n in name])
        else:
            bool_name = True

        if ac:
            bool_ac = True if iac in ac else False
        else:
            bool_ac = True

        if not (bool_typ and bool_name and bool_ac):
            continue
        else:
            if epoch_strt > dt.datetime(1980, 1, 1) or epoch_end < dt.datetime(
                2099, 1, 1
            ):
                bool_epoc = (epoch_strt <= clkdf_grp["epoch"]) & (
                    clkdf_grp["epoch"] < epoch_end
                )
                clkdf_stk.append(clkdf_grp[bool_epoc])
            else:
                clkdf_stk.append(clkdf_grp)

    clkdf_out = pd.concat(clkdf_stk)

    return clkdf_out


def clk_df_reg_2_multidx(clk_df_inp, index_order=["name", "epoch"]):
    """
    From an regular Clock DF generated by read_clk(), set some columns
    (typically ["name","epoch"]) as indexes
    The outputed DF is then a Multi-index DF

    It an adapted version of orb_df_reg_2_multidx
    """

    return orb_df_reg_2_multidx(clk_df_inp, index_order)


def clk_df_common_epoch_finder(
    ClkDFa_in,
    ClkDFb_in,
    return_index=False,
    supplementary_sort=False,
    order=["name", "epoch"],
):
    """
    Find common sats/station and epochs in to Clock DF, and output the
    corresponding Clock DFs

    Is an adapted version of orb_df_common_epoch_finder
    """

    return orb_df_common_epoch_finder(
        ClkDFa_in,
        ClkDFb_in,
        return_index=return_index,
        supplementary_sort=supplementary_sort,
        order=order,
    )


def clk_df_common_epoch_finder_multi(
    clk_df_list_in,
    return_index=False,
    supplementary_sort=False,
    order=["name", "epoch"],
):
    """
    Find common sats/station and epochs in to Clock DF, and output the
    corresponding Clock DFs

    Is is the multi version of clk_df_common_epoch_finder
    """

    clk_df_ref = clk_df_list_in[0]

    #### First loop: we find the common epochs
    for ClkDF in clk_df_list_in[1:]:

        OUTTUP = orb_df_common_epoch_finder(
            clk_df_ref,
            ClkDF,
            return_index=True,
            supplementary_sort=supplementary_sort,
            order=order,
        )

        clk_df_ref, _, Iinter = OUTTUP

    #### second loop: we use the common epochs found for the outputed ClkDF
    clk_df_list_out = []
    for ClkDF in clk_df_list_in:
        clk_df_out = ClkDF.set_index(order).loc[Iinter]
        clk_df_list_out.append(clk_df_out)

    if not return_index:
        return clk_df_list_out
    else:
        return clk_df_list_out, Iinter


#   _____ _      _____   __      __   _ _     _       _   _
#  / ____| |    |  __ \  \ \    / /  | (_)   | |     | | (_)
# | (___ | |    | |__) |  \ \  / /_ _| |_  __| | __ _| |_ _  ___  _ __
#  \___ \| |    |  _  /    \ \/ / _` | | |/ _` |/ _` | __| |/ _ \| '_ \
#  ____) | |____| | \ \     \  / (_| | | | (_| | (_| | |_| | (_) | | | |
# |_____/|______|_|  \_\     \/ \__,_|_|_|\__,_|\__,_|\__|_|\___/|_| |_|


def svn_prn_equiv_df(path_meta_snx):
    """
    generate a SVN <> PRN equivalent DataFrame

    Parameters
    ----------
    path_meta_snx : str
        path of the MGEX metadata sinex.
        last version avaiable here
        http://mgex.igs.org/IGS_MGEX_Metadata.php

    Returns
    -------
    DFfin : Pandas DataFrame
        SVN <> PRN equivalent DataFrame.

    """

    df_svn = files_rw.read_sinex_versatile(
        path_meta_snx, "SATELLITE/IDENTIFIER", header_line_idx=-2
    )

    df_prn = files_rw.read_sinex_versatile(
        path_meta_snx, "SATELLITE/PRN", header_line_idx=-2
    )

    df_svn.drop(columns="Comment__________________________________", inplace=True)
    df_prn.drop(columns="Comment_________________________________", inplace=True)

    ## the next lines 1 and 3 seems like they have became useless
    df_svn["SVN_"] = df_svn["SVN_"].apply(lambda x: x[0] + x[1:])
    df_prn.replace(dt.datetime(1970, 1, 1), dt.datetime(2099, 1, 1), inplace=True)
    df_prn["SVN_"] = df_prn["SVN_"].apply(lambda x: x[0] + x[1:])

    df_stk = []

    for isat, sat in df_prn.iterrows():
        svn = sat["SVN_"]

        sat["Block"] = df_svn[df_svn["SVN_"] == svn]["Block__________"].values[0]
        df_stk.append(sat)

    df_fin = pd.concat(df_stk, axis=1).transpose()

    df_fin.rename(
        columns={"SVN_": "SVN", "Valid_From____": "start", "Valid_To______": "end"},
        inplace=True,
    )

    df_fin["const"] = df_fin["SVN"].apply(lambda x: x[0])
    df_fin["SVN_int"] = df_fin["SVN"].apply(lambda x: int(x[1:]))
    df_fin["PRN_int"] = df_fin["PRN"].apply(lambda x: int(x[1:]))

    return df_fin


def svn_prn_equiv(sat_in, date_in, svn_prn_equiv_df, mode="svn2prn", full_output=False):
    """
    Get the equivalence SVN <> PRN for a given epoch

    Parameters
    ----------
    sat_in : str
        Satellite "ID", SVN or PRN.
    date_in : datetime
        wished epoch.
    svn_prn_equiv_df : DataFrame
        Equivalence table generated by svn_prn_equiv_df.
    mode : str, optional
        prn2svn: PRN > SVN
        svn2prn: SVN > PRN.
        The default is "svn2prn".
    full_output : bool, optional
        get the complete Equivalence table row. The default is False.

    Returns
    -------
    str or DataFrame
    """

    svnorprn1 = mode[:3].upper()
    svnorprn2 = mode[-3:].upper()

    df_sat = svn_prn_equiv_df[svn_prn_equiv_df[svnorprn1] == sat_in]
    bool_date = np.logical_and((df_sat.start <= date_in), (date_in < df_sat.end))
    df_out = df_sat[bool_date]

    if len(df_out) != 1:
        log.warning("several or no %s entries !!! %s %s", mode, sat_in, date_in)

    if full_output:
        return df_out
    else:
        return df_out[svnorprn2].values[0]


def get_block_svn(sat_in, svn_prn_equiv_df):
    """
    Get the equivalence SVN block type

    Parameters
    ----------
    sat_in : str
        Satellite SVN.
    svn_prn_equiv_df : DataFrame
        Equivalence table generated by svn_prn_equiv_df.
    Returns
    -------
    str with the block name
    """
    df_sat = svn_prn_equiv_df[svn_prn_equiv_df["SVN"] == sat_in]
    # if df_sat.empty:
    #     print('SVN NOT FOUND')
    # else:
    block = df_sat.Block.values[0]

    return block


def stats_slr(df_in, grpby_keys=["sat"], threshold=0.5):
    """
    computes statistics for SLR Residuals

    Parameters
    ----------
    df_in : Pandas DataFrame
        Input residual Dataframe from read_pdm_res_slr.
    grpby_keys : list of str, optional
        The default is ['sat'].
        per day, per solution, per satellite: ['day','sol','sat']
        per day, per solution, per station: ['day','sol','sta']
        per day, per solution, per satellite, per station: ['day','sol','sta','sat']
    threshold : float
        apply a Threshold

    Returns
    -------
    dd : Output statistics DataFrame
        return the mean, the rms and the std.
    """

    dd = df_in[np.abs(df_in["res"]) < threshold]

    dd_grp = dd.groupby(grpby_keys)
    dd_mean = dd_grp["res"].agg(np.mean).rename("mean") * 1000
    dd_rms = dd_grp["res"].agg(stats.rms_mean).rename("rms") * 1000
    dd_std = dd_grp["res"].agg(np.std).rename("std") * 1000
    dd = pd.concat([dd_mean, dd_std, dd_rms], axis=1)
    dd.reset_index(inplace=True)

    return dd


#  ______           _   _        ____       _            _        _   _               _____                               _
# |  ____|         | | | |      / __ \     (_)          | |      | | (_)             |  __ \                             | |
# | |__   __ _ _ __| |_| |__   | |  | |_ __ _  ___ _ __ | |_ __ _| |_ _  ___  _ __   | |__) |_ _ _ __ __ _ _ __ ___   ___| |_ ___ _ __ ___
# |  __| / _` | '__| __| '_ \  | |  | | '__| |/ _ \ '_ \| __/ _` | __| |/ _ \| '_ \  |  ___/ _` | '__/ _` | '_ ` _ \ / _ \ __/ _ \ '__/ __|
# | |___| (_| | |  | |_| | | | | |__| | |  | |  __/ | | | || (_| | |_| | (_) | | | | | |  | (_| | | | (_| | | | | | |  __/ ||  __/ |  \__ \
# |______\__,_|_|   \__|_| |_|  \____/|_|  |_|\___|_| |_|\__\__,_|\__|_|\___/|_| |_| |_|   \__,_|_|  \__,_|_| |_| |_|\___|\__\___|_|  |___/


### EOP / Earth Oreintation Parameters


def eop_interpotate(df_eop, epochs_intrp, eop_params=["x", "y"]):
    """
    Interopolate the EOP provided in a C04-like DataFrame

    Parameters
    ----------
    df_eop : DataFrame
        Input EOP DataFrame (C04 format).
        Can be generated by files_rw.read_eop_C04
    epochs_intrp : datetime of list of datetimes
        Wished epochs for the interpolation.
    eop_params : list of str, optional
        Wished EOP parameter to be interpolated.
        The default is ["x","y"].

    Returns
    -------
    out : DataFrame or Series
        Interpolated parameters.
        Series if onely one epoch is provided, df_eop elsewere
    """

    from geodezyx import interp

    if not utils.is_iterable(epochs_intrp):
        singleton = True
    else:
        singleton = False

    i_eop = dict()
    out_eop = dict()
    out_eop["epoch"] = epochs_intrp

    for eoppar in eop_params:
        intrp = interp.Interp1dTime(df_eop.epoch, df_eop[eoppar])
        i_eop[eoppar] = intrp
        try:
            out_eop[eoppar] = intrp(epochs_intrp)
        except ValueError as err:
            log.error("in EOP interpolation")
            log.error("param.: %s, epoch: %s", eoppar, epochs_intrp)
            raise err

    if not singleton:
        out = pd.DataFrame(out_eop)
    else:
        out = pd.Series(out_eop)

    return out
