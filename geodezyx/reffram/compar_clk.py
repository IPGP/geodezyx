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
import re

import matplotlib
import matplotlib.pyplot as plt
import natsort
import numpy as np
import pandas as pd

#### geodeZYX modules
from geodezyx import conv
from geodezyx import stats

### disabled and imported directly in the needed fct
## import geodezyx.reffram.sofa18 as sofa
log = logging.getLogger("geodezyx")

##########  END IMPORT  ##########


def compar_clock(clk_df_inp_1, clk_df_inp_2, col_name="name", bias_col_name="bias"):
    """
    Compares 2 GNSS clock bias DataFrames (from .clk), to a
    statistics table (with compar_clock_table)


    Parameters
    ----------
    clk_df_inp_1 & clk_df_inp_2 : DataFrame
        Clock DataFrame provided by files_rw.read_clk()

    Returns
    -------
    df_clk_diff : DataFrame
        Clock bias difference DataFrame
    """
    df1idx = clk_df_inp_1.set_index([col_name, "epoch"])
    df1idx.sort_index(inplace=True)

    df2idx = clk_df_inp_2.set_index([col_name, "epoch"])
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
    if "ac" in clk_df_inp_1.columns:
        df_clk_diff["name1"] = clk_df_inp_1.ac.values[0]
    if "AC" in clk_df_inp_1.columns:
        df_clk_diff["name1"] = clk_df_inp_1.AC.values[0]

    if "ac" in clk_df_inp_1.columns:
        df_clk_diff["name2"] = clk_df_inp_2.ac.values[0]
    if "AC" in clk_df_inp_1.columns:
        df_clk_diff["name2"] = clk_df_inp_2.AC.values[0]

    return df_clk_diff


def compar_clock_table(df_clk_diff_inp, col_name="name", bias_col_name="bias_diff"):
    """
    Generate a table with statistical indicators for a clock comparison
    (RMS mean, standard dev, ...)

    Parameters
    ----------
    df_clk_diff_inp : DataFrame
        Clock bias difference DataFrame (from compar_clock)

    col_name : str
        the column with the sat names

    bias_col_name : str
        The column with the clk difference values


    Returns
    -------
    compar_df_out : DataFrame
        Statistical results of the comparison.

    """

    df_diff_grp = df_clk_diff_inp.groupby(col_name)[bias_col_name]

    smin = df_diff_grp.min().rename("min", inplace=True)
    smax = df_diff_grp.max().rename("max", inplace=True)
    smean = df_diff_grp.mean().rename("mean", inplace=True)
    sstd = df_diff_grp.std().rename("std", inplace=True)
    srms = df_diff_grp.apply(stats.rms_mean).rename("rms", inplace=True)

    df_compar_out = pd.concat([smin, smax, smean, sstd, srms], axis=1)
    df_compar_out.reset_index()

    return df_compar_out


def compar_clk_plot(
    diff_sat_all_df_inp,
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
    diff_sat_all_df_inp: DataFrame
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
    diff_sat_all_df_inp = diff_sat_all_df_inp.reset_index()
    satdispo = natsort.natsorted(list(set(diff_sat_all_df_inp[col_name])))
    # satdispo = natsort.natsorted(list(set(diff_sat_all_df_inp['sat'])))

    symb_stk = []

    cm = plt.get_cmap("viridis")
    num_colors = len(satdispo)
    colors = [cm(1.0 * i / num_colors) for i in range(num_colors)]

    date = conv.numpy_dt2dt(diff_sat_all_df_inp.epoch.values[0])
    diff_sat_all_df_inp.name = " ".join(
        (
            "Clock comparison  b/w",
            diff_sat_all_df_inp.name1.values[0],
            "(ref.) and",
            diff_sat_all_df_inp.name2.values[0],
            ",",
            date.strftime("%Y-%m-%d"),
        )
    )
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
        diffuse = diff_sat_all_df_inp[diff_sat_all_df_inp[col_name] == satuse]

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
    plt.suptitle(diff_sat_all_df_inp.name)
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)
    plt.subplots_adjust(bottom=0.15)

    if save_plot:
        if save_plot_name == "auto":
            save_plot_name = (
                diff_sat_all_df_inp.name1.values[0]
                + "_"
                + diff_sat_all_df_inp.name2.values[0]
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
