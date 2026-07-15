#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.reffram contains functions for operations
related to GNSS-products


it can be imported directly with:
from geodezyx import reffram

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
from geodezyx import files_rw
from geodezyx import stats
from geodezyx import utils
from geodezyx.reffram import kepler_gzyx
from geodezyx.utils_xtra.pandas_utils import df_common_finder, df_reg_2_multidx

### disabled and imported directly in the needed fct
## import geodezyx.reffram.sofa18 as sofa
log = logging.getLogger("geodezyx")

##########  END IMPORT  ##########







#   ____       _     _ _     _____        _        ______
#  / __ \     | |   (_) |   |  __ \      | |      |  ____|
# | |  | |_ __| |__  _| |_  | |  | | __ _| |_ __ _| |__ _ __ __ _ _ __ ___   ___  ___
# | |  | | '__| '_ \| | __| | |  | |/ _` | __/ _` |  __| '__/ _` | '_ ` _ \ / _ \/ __|
# | |__| | |  | |_) | | |_  | |__| | (_| | || (_| | |  | | | (_| | | | | | |  __/\__ \
#  \____/|_|  |_.__/|_|\__| |_____/ \__,_|\__\__,_|_|  |_|  \__,_|_| |_| |_|\___||___/
#

### Orbit DataFrames


def orb_df_velocity_calc(orb_df_inp, drop_nan=False):
    """
    Compute the velocity of satellites from an orbit dataframe
    (differentiate the position)

    Parameters
    ----------
    orb_df_inp : Pandas DataFrame
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

    dfgrp = orb_df_inp.groupby("prn")

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
    sun_dec : float or pandas Series
        Sun's declination.
    sun_ra : float or pandas Series
        Sun's right ascension.
    sat_i : float or pandas Series
        Satellite's inclination.
    sat_o_lan : float
        Satellite's Longitude of the ascending node.

    Returns
    -------
    beta : float or pandas Series
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
    orb_df_inp,
    calc_beta_sun_ra_dec=True,
    calc_beta_sun_eclip_long=True,
    beta_rad2deg=True,
):
    """
    Compute beta angle for GNSS satellite's orbits stored in an orbit
    DataFrame


    Parameters
    ----------
    orb_df_inp : Pandas DataFrame
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
        df_orb_inp with a new 'beta' column.
    df_wrk : Pandas DataFrame
        intermediate values dataframe for debug.
        here, coordinates are in ECI frame
    """
    import pyorbital.astronomy

    #### convert in ECI
    df_eci = orb_df_inp.copy()
    df_eci[["x", "y", "z"]] = conv.ecef2eci(
        orb_df_inp[["x", "y", "z"]].values, orb_df_inp["epoch"].values
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
    orb_df_out = orb_df_inp.copy()
    orb_df_out["beta"] = df_wrk["beta" + b1]

    return orb_df_out, df_wrk


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
    orb_df_out : DataFrame
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

        df_orb_use = orb_df_inp[
            (orb_df_inp["prn"] == sat) & (orb_df_inp["ac"] == ac)
        ].copy()

        ### faster but anoying Future Waring
        # Tdata = np.array(df_orb_use.epoch.dt.to_pydatetime())

        #t_type = "datetime"
        #t_data = conv.numpy_dt2dt(df_orb_use['epoch'].values)

        t_type = "pandas_timestamp"
        t_data = df_orb_use["epoch"]
        titrp_use = pd.Series(titrp)

        xitrp = stats.lagrange_interpolate(t_data, df_orb_use["x"], titrp_use, n=n, t_type=t_type)
        yitrp = stats.lagrange_interpolate(t_data, df_orb_use["y"], titrp_use, n=n, t_type=t_type)
        zitrp = stats.lagrange_interpolate(t_data, df_orb_use["z"], titrp_use, n=n, t_type=t_type)
        clk_itrp = np.interp(
            conv.pandas_timestamp2posix(titrp_use),
            conv.pandas_timestamp2posix(t_data),
            df_orb_use["clk"].values,
        )

        # ClkDummy = np.array([999999.999999] * len(titrp))

        d = {"epoch": titrp, "x": xitrp, "y": yitrp, "z": zitrp, "clk": clk_itrp}
        orb_df_tmp = pd.DataFrame(d)

        ### sometihng else must be tested o give the annex val directly in the col of orb_df_tmp
        annex_df_vals = df_orb_use.drop(
            ["epoch", "x", "y", "z", "clk"], axis=1
        ).drop_duplicates()
        annex_df_vals = pd.concat(
            [annex_df_vals] * (len(titrp)), ignore_index=True, axis=0
        )

        orb_df_tmp = pd.concat((orb_df_tmp, annex_df_vals), axis=1)
        df_orb_stk.append(orb_df_tmp)

        if plot:
            # plt.plot(Tdata,df_orb_use.x,'o')
            # plt.plot(titrp,Xitrp,'.')
            ## GUS mod 220322
            fig, axr = plt.subplots(1, 1, sharex="all")
            Symb = axr.plot(t_data, df_orb_use.x, "o")
            Symb = axr.plot(titrp, xitrp, ".")

    orb_df_out = pd.concat(df_orb_stk)

    if append_to_input_df:
        orb_df_out = pd.concat((orb_df_inp, orb_df_out))

    orb_df_out.reset_index(drop=True)
    orb_df_out[["x", "y", "z", "clk"]] = orb_df_out[["x", "y", "z", "clk"]].astype(float)
    return orb_df_out


def orb_df_crf2trf(orb_df_inp, eop_df_inp, time_scale_inp="gps", inv_trf2crf=False):
    """
    Convert an Orbit DataFrame from Celetrial Reference Frame to
    Terrestrial Reference Frame.

    Requires EOP to work. Cf. note below.

    Parameters
    ----------
    orb_df_inp : DataFrame
        Input Orbit DataFrame in Celetrial Reference Frame.
    eop_df_inp : DataFrame
        EOP DataFrame  (C04 format).
    time_scale_inp : str, optional
        The time scale used in. manage 'utc', 'tai' and 'gps'.
        The default is "gps".
    inv_trf2crf : bool, optional
        Provide the inverse transformation TRF => CRF.
        The default is False.

    Returns
    -------
    orb_df_out : DataFrame
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

    orb_df_use = orb_df_inp.copy()

    import geodezyx.reffram.sofa as sofa

    ### bring everything to UTC
    if time_scale_inp.lower() == "gps":
        orb_df_use["epoch_utc"] = conv.dt_gpstime2dt_utc(orb_df_use["epoch"])
    elif time_scale_inp.lower() == "tai":
        orb_df_use["epoch_utc"] = conv.dt_tai2dt_utc(orb_df_use["epoch"])
    elif time_scale_inp.lower() == "utc":
        orb_df_use["epoch_utc"] = orb_df_use["epoch"]
    ### TT and UT1 are not implemented (quite unlikely to have them as input)

    ### do the time scale's conversion
    orb_df_use["epoch_tai"] = conv.dt_utc2dt_tai(orb_df_use["epoch_utc"])
    orb_df_use["epoch_tt"] = conv.dt_tai2dt_tt(orb_df_use["epoch_tai"])
    orb_df_use["epoch_ut1"] = conv.dt_utc2dt_ut1_smart(orb_df_use["epoch_utc"], eop_df_inp)

    ### Do the EOP interpolation
    eop_df_intrp = eop_interpotate(eop_df_inp, orb_df_use["epoch_utc"])
    ### bring the EOP to radians
    xeop = np.deg2rad(conv.arcsec2deg(eop_df_intrp["x"]))
    yeop = np.deg2rad(conv.arcsec2deg(eop_df_intrp["y"]))

    trf_stk = []

    for tt, ut1, xeop, yeop, x, y, z in zip(
        orb_df_use["epoch_tt"],
        orb_df_use["epoch_ut1"],
        xeop,
        yeop,
        orb_df_use["x"],
        orb_df_use["y"],
        orb_df_use["z"],
    ):

        mat_crf2trf = sofa.iau_c2t06a(
            2400000.5, conv.dt2mjd(tt), 2400000.5, conv.dt2mjd(ut1), xeop, yeop
        )
        if inv_trf2crf:
            mat_use = np.linalg.inv(mat_crf2trf)
        else: # regular case
            mat_use = mat_crf2trf

        crf = np.array([x, y, z])
        trf = np.dot(mat_use, crf)

        trf_stk.append(trf)

    ### Final stack and replacement
    trf_all = np.vstack(trf_stk)
    orb_df_out = orb_df_inp.copy()
    orb_df_out[["x", "y", "z"]] = trf_all

    return orb_df_out


#### FCT DEF

def orb_df_multidx_2_reg(orb_df_inp, index_order=["prn", "epoch"]):
    """
    Convert a Multi-index formatted OrbitDF to his original form
    """
    orb_df_wrk = orb_df_inp.reset_index()
    orb_df_wrk = orb_df_wrk.sort_values(index_order)
    orb_df_wrk["sys"] = orb_df_wrk["prn"].apply(lambda x: x[0])
    orb_df_wrk["prni"] = orb_df_wrk["prn"].apply(lambda x: int(x[1:]))
    return orb_df_wrk


def orb_df_common_epoch_finder(
    orb_df_a_inp,
    orb_df_b_inp,
    return_index=False,
    supplementary_sort=False,
    order=["prn", "epoch"],
    skip_reg2multidx_orb_df_a=False,
    skip_reg2multidx_orb_df_b=False,
):
    """
    This function finds common satellites and epochs in two Orbit DataFrames and outputs the corresponding Orbit DataFrames.

    Parameters
    ----------
    orb_df_a_inp : DataFrame
        The first input Orbit DataFrame.
    orb_df_b_inp : DataFrame
        The second input Orbit DataFrame.
    return_index : bool, optional
        If True, the function also returns the common index. Default is False.
    supplementary_sort : bool, optional
        If True, an additional sort is performed. This is useful for multi GNSS where the output DataFrame may not be well sorted. Default is False.
    order : list of str, optional
        The order of the index for the multi-index DataFrame. Default is ["prn","epoch"].
    skip_reg2multidx_orb_df_a : bool, optional
        If True, skips the conversion of the first input DataFrame to a multi-index DataFrame. Default is False.
        The inputs are assumed to be already in multi-index format to optimize execution speed.
        (For advanced use only)
    skip_reg2multidx_orb_df_b : bool, optional
        If True, skips the conversion of the second input DataFrame to a multi-index DataFrame. Default is False.
        The inputs are assumed to be already in multi-index format to optimize execution speed.
        (For advanced use only)

    Returns
    -------
    orb_df_a_out : DataFrame
        The first output Orbit DataFrame with common satellites and epochs.
    orb_df_b_out : DataFrame
        The second output Orbit DataFrame with common satellites and epochs.
    iinter : Index, optional
        The common index. Only returned if return_index is True.

    Note
    ----
    designed for orbits/sp3 first with sat and epoch as order parmeter,
    but can be used also for instance for snx files with
    STAT and epoch as order parmeter

    This function is a wrapper for df_common_finder() specialized for orbit DataFrames.
    """
    return df_common_finder(
        orb_df_a_inp,
        orb_df_b_inp,
        return_index=return_index,
        supplementary_sort=supplementary_sort,
        order=order,
        skip_reg2multidx_df_a=skip_reg2multidx_orb_df_a,
        skip_reg2multidx_df_b=skip_reg2multidx_orb_df_b,
    )



def orb_df_const_sv_columns_maker(orb_df_inp, inplace=True):
    """
    (re)generate the const and sv columns from the sat one
    """
    if inplace:
        orb_df_inp["sys"] = orb_df_inp["prn"].str[0]
        orb_df_inp["prni"] = orb_df_inp["prn"].apply(lambda x: int(x[1:]))
        return None
    else:
        orb_df_out = orb_df_inp.copy()
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
    clk_df_inp,
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
    clk_df_inp : DataFrame
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

    if type(clk_df_inp) is str:
        clk_df_wrk = utils.pickle_loader(clk_df_inp)
    else:
        clk_df_wrk = clk_df_inp

    bool_fin = np.ones(len(clk_df_wrk)).astype(bool)

    if typ:
        bool_tmp = clk_df_wrk.type.isin(typ)
        bool_fin = bool_fin & np.array(bool_tmp)

    if name:
        if not name_regex:  ### full name mode
            bool_tmp = clk_df_wrk.name.isin(name)
            bool_fin = bool_fin & np.array(bool_tmp)
        else:  ### REGEX mode
            bool_tmp = np.zeros(len(clk_df_wrk.name)).astype(bool_fin)
            for rgx in name:
                nam_serie = clk_df_wrk.name
                bool_tmp = bool_tmp | np.array(nam_serie.str.contains(rgx))

            bool_fin = bool_fin & np.array(bool_tmp)

    if ac:
        bool_tmp = clk_df_wrk.ac.isin(ac)
        bool_fin = bool_fin & np.array(bool_tmp)

    ##epoch
    bool_tmp = (epoch_strt <= clk_df_wrk.epoch) & (clk_df_wrk.epoch < epoch_end)
    bool_fin = bool_fin & np.array(bool_tmp)

    return clk_df_wrk[bool_fin]


def clk_df_filter2(
    clk_df_inp,
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

    if type(clk_df_inp) is str:
        clk_df_wrk = utils.pickle_loader(clk_df_inp)
    else:
        clk_df_wrk = clk_df_inp

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



def clk_df_common_epoch_finder(
    clk_df_a_inp,
    clk_df_b_inp,
    return_index=False,
    supplementary_sort=False,
    order=["name", "epoch"],
):
    """
    Find common sats/station and epochs in two Clock DataFrames and output the
    corresponding Clock DataFrames.

    Parameters
    ----------
    clk_df_a_inp : DataFrame
        The first input Clock DataFrame.
    clk_df_b_inp : DataFrame
        The second input Clock DataFrame.
    return_index : bool, optional
        If True, the function also returns the common index. Default is False.
    supplementary_sort : bool, optional
        If True, an additional sort is performed. Default is False.
    order : list of str, optional
        The order of the index for the multi-index DataFrame. Default is ["name","epoch"].

    Returns
    -------
    clk_df_a_out : DataFrame
        The first output Clock DataFrame with common stations/satellites and epochs.
    clk_df_b_out : DataFrame
        The second output Clock DataFrame with common stations/satellites and epochs.
    iinter : Index, optional
        The common index. Only returned if return_index is True.

    Note
    ----
    This function is a wrapper for df_common_finder() specialized for clock DataFrames.
    """
    return df_common_finder(
        clk_df_a_inp,
        clk_df_b_inp,
        return_index=return_index,
        supplementary_sort=supplementary_sort,
        order=order,
    )


def clk_df_common_epoch_finder_multi(
    clk_df_list_inp,
    return_index=False,
    supplementary_sort=False,
    order=["name", "epoch"],
):
    """
    Find common sats/station and epochs in to Clock DF, and output the
    corresponding Clock DFs

    Is is the multi version of clk_df_common_epoch_finder
    """

    clk_df_ref = clk_df_list_inp[0]

    #### First loop: we find the common epochs
    for clk_df in clk_df_list_inp[1:]:

        outtup = clk_df_common_epoch_finder(
            clk_df_ref,
            clk_df,
            return_index=True,
            supplementary_sort=supplementary_sort,
            order=order,
        )

        clk_df_ref, _, iinter = outtup

    #### second loop: we use the common epochs found for the outputed ClkDF
    clk_df_list_out = []
    for clk_df in clk_df_list_inp:
        clk_df_out = clk_df.set_index(order).loc[iinter]
        clk_df_list_out.append(clk_df_out)

    if not return_index:
        return clk_df_list_out
    else:
        return clk_df_list_out, iinter


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


def svn_prn_equiv(sat_in, date_in, svn_prn_equiv_df_inp, mode="svn2prn", full_output=False):
    """
    Get the equivalence SVN <> PRN for a given epoch

    Parameters
    ----------
    sat_in : str
        Satellite "ID", SVN or PRN.
    date_in : datetime
        wished epoch.
    svn_prn_equiv_df_inp : DataFrame
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

    df_sat = svn_prn_equiv_df_inp[svn_prn_equiv_df_inp[svnorprn1] == sat_in]
    bool_date = np.logical_and((df_sat.start <= date_in), (date_in < df_sat.end))
    df_out = df_sat[bool_date]

    if len(df_out) != 1:
        log.warning("several or no %s entries !!! %s %s", mode, sat_in, date_in)

    if full_output:
        return df_out
    else:
        return df_out[svnorprn2].values[0]


def get_block_svn(sat_inp, svn_prn_equiv_df_inp):
    """
    Get the equivalence SVN block type

    Parameters
    ----------
    sat_inp : str
        Satellite SVN.
    svn_prn_equiv_df_inp : DataFrame
        Equivalence table generated by svn_prn_equiv_df.
    Returns
    -------
    str with the block name
    """
    df_sat = svn_prn_equiv_df_inp[svn_prn_equiv_df_inp["SVN"] == sat_inp]
    # if df_sat.empty:
    #     print('SVN NOT FOUND')
    # else:
    block = df_sat.Block.values[0]

    return block


def stats_slr(df_inp, grpby_keys=["sat"], threshold=0.5):
    """
    computes statistics for SLR Residuals

    Parameters
    ----------
    df_inp : Pandas DataFrame
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

    dd = df_inp[np.abs(df_inp["res"]) < threshold]

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
