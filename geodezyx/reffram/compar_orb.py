#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

The geodezyx toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/IPGP/geodezyx
"""

########## BEGIN IMPORT ##########
#### External modules
import datetime as dt
import itertools

#### Import the logger
import logging
import os

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
                    conv.dt_gpstime2dt_utc(p1.index.to_pydatetime()),
                )
                p2b = conv.ecef2eci(
                    np.array(p2),
                    conv.dt_gpstime2dt_utc(p2.index.to_pydatetime()),
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
    diff_sat_all_df_inp,
    save_plot=False,
    save_plot_dir="",
    save_plot_name="auto",
    save_plot_name_suffix=None,
    save_plot_ext=(".pdf", ".png", ".svg"),
    yaxis_limit=None,
    yaxis_label_unit="m",
):
    """
    Generates a plot for orbit comparison data.

    Parameters
    ----------
    diff_sat_all_df_inp : DataFrame
        A DataFrame produced by the `compar_orbit` function, containing orbit comparison data.
    save_plot : bool, optional
        If True, saves the plot to a file instead of returning the Figure and Axes. Default is False.
    save_plot_dir : str, optional
        Directory where the plot will be saved. Default is an empty string.
    save_plot_name : str, optional
        Name of the saved plot file. If "auto", the name is generated automatically. Default is "auto".
    save_plot_name_suffix : str, optional
        Suffix to append to the saved plot name. Default is None.
    save_plot_ext : tuple of str, optional
        File extensions for the saved plot. Default is (".pdf", ".png", ".svg").
    yaxis_limit : tuple or list, optional
        Limits for the y-axis. Can be a 3-tuple for independent axis limits or a 2-tuple for uniform limits. Default is None.
    yaxis_label_unit : str, optional
        Unit label for the y-axis. Default is "m".

    Returns
    -------
    tuple or str
        If `save_plot` is False, returns the Figure and Axes as a tuple. If `save_plot` is True, returns the path of the saved plot file.

    Notes
    -----
    - The function creates three subplots for Radial, Transverse, and Normal differences (RTN) or X, Y, Z differences.
    - The y-axis limits can be customized for each subplot or applied uniformly.
    - If `save_plot` is True, the plot is saved in the specified directory with the specified name and extensions.
    - The function uses the `viridis` colormap for distinguishing satellites.

    Example
    -------
    >>> fig, axes = compar_orbit_plot(diff_sat_all_df_inp)
    >>> plt.show()
    """

    fig, [axr, axt, axn] = plt.subplots(3, 1, sharex="all")

    satdispo = natsort.natsorted(list(set(diff_sat_all_df_inp["prn"])))

    symb_stk = []

    cm = plt.get_cmap("viridis")
    num_colors = len(satdispo)
    colors = [cm(1.0 * i / num_colors) for i in range(num_colors)]

    # Pandas donesn't manage well iterable as attribute
    # So, it is separated
    try:
        col_name0 = diff_sat_all_df_inp.frame_col_name1
        col_name1 = diff_sat_all_df_inp.frame_col_name2
        col_name2 = diff_sat_all_df_inp.frame_col_name3
    except:
        col_name0 = diff_sat_all_df_inp.columns[0]
        col_name1 = diff_sat_all_df_inp.columns[1]
        col_name2 = diff_sat_all_df_inp.columns[2]

    for satuse, color in zip(satdispo, colors):
        diffuse = diff_sat_all_df_inp[diff_sat_all_df_inp["prn"] == satuse]

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

    if diff_sat_all_df_inp.frame_type == "RTN":
        axr.set_ylabel("Radial diff." + ylabuni)
        axt.set_ylabel("Transverse diff." + ylabuni)
        axn.set_ylabel("Normal diff." + ylabuni)

    else:
        axr.set_ylabel(diff_sat_all_df_inp.frame_type + " X diff." + ylabuni)
        axt.set_ylabel(diff_sat_all_df_inp.frame_type + " Y diff." + ylabuni)
        axn.set_ylabel(diff_sat_all_df_inp.frame_type + " Z diff." + ylabuni)

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
    plt.suptitle(diff_sat_all_df_inp.name)
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)
    plt.subplots_adjust(bottom=0.15)

    if save_plot:
        if save_plot_name == "auto":
            save_plot_name = "_".join(
                (
                    diff_sat_all_df_inp.name1,
                    diff_sat_all_df_inp.name2,
                    diff_sat_all_df_inp.index.min().strftime("%Y-%m-%d"),
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


def compar_orbit_table(diff_sat_all_df_inp, rms_style="natural", light_tab=False):
    """
    Generate a table with statistical indicators for an orbit comparison
    (RMS mean, standard dev, ...)

    Parameters
    ----------
    diff_sat_all_df_inp : Pandas DataFrame
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
    >>> print(tabulate(compar_table, headers="keys", floatfmt=".4f"))
    """

    sat_list = utils.uniq_and_sort(diff_sat_all_df_inp["prn"])

    # Pandas donesn't manage well iterable as attribute
    # So, it is separated
    try:
        col_name0 = diff_sat_all_df_inp.frame_col_name1
        col_name1 = diff_sat_all_df_inp.frame_col_name2
        col_name2 = diff_sat_all_df_inp.frame_col_name3
    except:
        col_name0 = diff_sat_all_df_inp.columns[0]
        col_name1 = diff_sat_all_df_inp.columns[1]
        col_name2 = diff_sat_all_df_inp.columns[2]

    rms_stk = []

    for sat in sat_list:
        diffwork = utils.df_sel_val_in_col(diff_sat_all_df_inp, "prn", sat)

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
        else:
            raise ValueError("Unknown rms_style : " + str(rms_style))

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
        rms_a = stats.rms_mean(diff_sat_all_df_inp[col_name0])
        rms_b = stats.rms_mean(diff_sat_all_df_inp[col_name1])
        rms_c = stats.rms_mean(diff_sat_all_df_inp[col_name2])
        rms_3d = np.sqrt(rms_a**2 + rms_b**2 + rms_c**2)
    elif rms_style == "GRGS":
        rms_a = stats.rms_mean(
            diff_sat_all_df_inp[col_name0] - diff_sat_all_df_inp[col_name0].mean()
        )
        rms_b = stats.rms_mean(
            diff_sat_all_df_inp[col_name1] - diff_sat_all_df_inp[col_name1].mean()
        )
        rms_c = stats.rms_mean(
            diff_sat_all_df_inp[col_name2] - diff_sat_all_df_inp[col_name2].mean()
        )
        rms_3d = np.sqrt(rms_a**2 + rms_b**2 + rms_c**2)
    elif rms_style == "kouba":
        rms_a = stats.rms_mean_kouba(diff_sat_all_df_inp[col_name0])
        rms_b = stats.rms_mean_kouba(diff_sat_all_df_inp[col_name1])
        rms_c = stats.rms_mean_kouba(diff_sat_all_df_inp[col_name2])
        rms_3d = np.sqrt(rms_a**2 + rms_b**2 + rms_c**2)
    else:
        raise ValueError("Unknown rms_style : " + str(rms_style))

    min_a = diff_sat_all_df_inp[col_name0].min()
    min_b = diff_sat_all_df_inp[col_name1].min()
    min_c = diff_sat_all_df_inp[col_name2].min()

    max_a = diff_sat_all_df_inp[col_name0].max()
    max_b = diff_sat_all_df_inp[col_name1].max()
    max_c = diff_sat_all_df_inp[col_name2].max()

    mean_a = diff_sat_all_df_inp[col_name0].mean()
    mean_b = diff_sat_all_df_inp[col_name1].mean()
    mean_c = diff_sat_all_df_inp[col_name2].mean()

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

    if diff_sat_all_df_inp.frame_type == "RTN":
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
