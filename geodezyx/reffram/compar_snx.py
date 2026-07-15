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
    """
    Compare 2 SINEX files and provide statistics on station position differences

    Parameters
    ----------
    snx1, snx2 : str or DataFrame
        Paths to the SINEX files to compare or DataFrames containing SINEX data
    stat_select : str or list, optional
        Regular expression or list of station codes to filter the comparison
    invert_select : bool, optional
        If True, invert the selection criteria for stations
    out_means_summary : bool, optional
        If True, return a summary of means; otherwise, return the full DataFrame
    out_meta : bool, optional
        If True, include metadata (week, day, number of stations) in the output
    out_dataframe : bool, optional
        If True, return the output as a DataFrame; otherwise, return as a tuple
    manu_wwwwd : str, optional
        Manual specification of week/day if not extracted from filename
    """

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

    d3_d = np.sqrt(
        (ddiff["x"] ** 2 + ddiff["y"] ** 2 + ddiff["z"] ** 2).astype("float64")
    )

    ddiff = ddiff.assign(d3D_xyz=d3_d)

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

    d2_d = np.sqrt((E**2 + N**2).astype("float64"))
    d3_d = np.sqrt((E**2 + N**2 + U**2).astype("float64"))

    ddiff = ddiff.assign(e=E)
    ddiff = ddiff.assign(n=N)
    ddiff = ddiff.assign(u=U)
    ddiff = ddiff.assign(d2D_enu=d2_d)
    ddiff = ddiff.assign(d3D_enu=d3_d)

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

            output_df = pd.DataFrame(output).transpose()

            output_df.columns = [
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

            return output_df
