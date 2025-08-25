#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.reffram contains functions related to
 Kepler's element operations.

it can be imported directly with:
from geodezyx import reffram

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""

import datetime as dt

#### Import the logger
import logging

########## BEGIN IMPORT ##########
#### External modules
import numpy as np
import pandas as pd

# from pytwobodyorbit import TwoBodyOrbit

#### geodeZYX modules
from geodezyx import conv

log = logging.getLogger("geodezyx")

##########  END IMPORT  ##########


def eci_2_kepler_elts(pos, vel, rad2deg=True, mu=3.9860044188e14):
    """
    Convert ECI coordinates > Kepler's elements

    Parameters
    ----------
    pos : array
        Position in m.
        Can be a 3-element (x,y,z) array (simple point)
        or a (N,3)-shaped numpy array
    vel : array
        Velocity in m/s.
        Can be a 3-element (vx,vy,vz) array (simple point)
        or a (N,3)-shaped numpy array
    rad2deg : bool, optional
        convert Keplerian's angles (anomalies) to deg. The default is True.
    mu : float, optional
        Standard gravitational parameter. The default is 3.9860044188e14.

    Returns
    -------
    a,ecc,i,o_peri,o_lan,m: floats
        Kepler's elements:

        * a : semi-major axis
        * ecc : orbit eccentricity
        * i : orbit inclination
        * o_peri : argument of periapsis ω
        * o_lan : longitude of the ascending node Ω
        * m : mean anomaly m

    Note
    ----
    Be sure your input is in ECI (SP3 are in ECEF for instance).

    If neeeded, use ``conv.ecef2eci()`` (gross results)
    or ``reffram.orb_df_crf2trf(inv_trf2crf=True)`` (precise results)

    You can get the velocity with ``reffram.orb_df_velocity_calc()``

    Source
    ------
    https://downloads.rene-schwarz.com/download/M002-Cartesian_State_Vectors_to_Keplerian_Orbit_Elements.pdf

    """
    p = np.array(pos)
    v = np.array(vel)

    #### case for one point only
    if len(p.shape) < 2:
        p = np.expand_dims(p, axis=-1).T
    if len(v.shape) < 2:
        v = np.expand_dims(v, axis=-1).T

    pnorm = np.linalg.norm(p, axis=1)
    vnorm = np.linalg.norm(v, axis=1)

    ### orbital momentum vector h
    h = np.cross(p, v)
    hnorm = np.linalg.norm(h, axis=1)

    ### eccentricity vector eccvec & orbit eccentricity ecc
    eccvec = (np.cross(v, h) / mu) - (p / np.expand_dims(pnorm, axis=-1))
    ecc = np.linalg.norm(eccvec, axis=1)

    ### true anomaly ν (nu)
    n = np.cross(np.array([0, 0, 1]), h)
    nnorm = np.linalg.norm(n, axis=1)

    dot_pv_sup0 = np.einsum("ij,ij->i", p, v) >= 0
    nu_sup0 = np.arccos(np.einsum("ij,ij->i", eccvec, p) / (ecc * pnorm))
    nu_inf0 = 2 * np.pi - nu_sup0

    nu = np.zeros(len(dot_pv_sup0))
    nu = nu + nu_sup0 * dot_pv_sup0
    nu = nu + nu_inf0 * np.logical_not(dot_pv_sup0)

    ### orbit inclination i
    i = np.arccos(h[:, 2] / hnorm)

    ### eccentric anomaly E (called e here)
    # e1 = 2*np.arctan(np.tan(nu/2)/np.sqrt((1+ecc)/(1-ecc)))
    e2 = 2 * np.arctan2(np.tan(nu / 2), np.sqrt((1 + ecc) / (1 - ecc)))
    e = e2

    ### longitude of the ascending node Ω (o_lan)
    n1_sup0 = n[:, 1] >= 0

    omega_lan_sup0 = np.arccos(n[:, 0] / nnorm)
    omega_lan_inf0 = 2 * np.pi - omega_lan_sup0

    o_lan = np.zeros(len(n1_sup0))
    o_lan = o_lan + omega_lan_sup0 * n1_sup0
    o_lan = o_lan + omega_lan_inf0 * np.logical_not(n1_sup0)

    ### argument of periapsis ω (o_peri)
    eccvec2_sup0 = eccvec[:, 2] >= 0
    o_peri_sup0 = np.arccos((np.einsum("ij,ij->i", n, eccvec)) / (ecc * nnorm))
    o_peri_inf0 = 2 * np.pi - o_peri_sup0

    o_peri = np.zeros(len(eccvec2_sup0))
    o_peri = o_peri + o_peri_sup0 * eccvec2_sup0
    o_peri = o_peri + o_peri_inf0 * np.logical_not(eccvec2_sup0)

    ### mean anomaly m, (Kepler’s Equation)
    m = e - ecc * np.sin(e)

    ### semi-major axis a
    a = ((2 / pnorm) - (vnorm**2 / mu)) ** -1

    #### case for one point only (bring it back to a scalar)
    sqzflt = lambda x: np.float64(x.squeeze())

    a = sqzflt(a)
    ecc = sqzflt(ecc)
    i = sqzflt(i)
    o_lan = sqzflt(o_lan)
    o_peri = sqzflt(o_peri)
    m = sqzflt(m)

    if rad2deg:
        i = np.rad2deg(i)
        o_peri = np.rad2deg(o_peri)
        o_lan = np.rad2deg(o_lan)
        m = np.rad2deg(m)

    #### o_peri and m seems to be inverted, must be investigated!!!
    return a, ecc, i, o_peri, o_lan, m


def eci_2_kepler_elts_mono(p, v, rad2deg=True, mu=3.9860044188e14):
    """
    **Kept for debug purposes, use** ``eci_2_kepler_elts(pos, vel)`` **instead**

    Convert ECI coordinates > Kepler's elements

    Parameters
    ----------
    p : array
        Position.
    v : array
        Velocity.
    rad2deg : bool, optional
        canvert Keplerian's angles (anomalies) to deg. The default is True.
    mu : float, optional
        Standard gravitational parameter. The default is 3.9860044188e14.

    Returns
    -------
    a,ecc,i,omega_periarg,omega_lan,m : floats
        Kepler's elements.
    """
    log.warning(
        "use this function for debug purposes only, use eci_2_kepler_elts instead"
    )
    pnorm = np.linalg.norm(p)
    vnorm = np.linalg.norm(v)

    ### orbital momentum vector h
    h = np.cross(p, v)
    hnorm = np.linalg.norm(h)

    ### eccentricity vector ecc
    ecc = (np.cross(v, h) / mu) - (p / pnorm)
    ecc_norm = np.linalg.norm(ecc)

    ### true anomaly ν (nu)
    n = np.cross(np.array([0, 0, 1]), h)
    nnorm = np.linalg.norm(n)
    if np.dot(p, v) >= 0:
        nu = np.arccos(np.dot(ecc, p) / (ecc_norm * pnorm))
    else:
        nu = 2 * np.pi - np.arccos(np.dot(ecc, p) / (ecc_norm * pnorm))

    ### orbit inclination i
    i = np.arccos(h[2] / hnorm)

    ### orbit eccentricity ecc
    ecc = np.linalg.norm(ecc)

    ### eccentric anomaly E (called e here)
    # e1 = 2*np.arctan(np.tan(nu/2)/np.sqrt((1+ecc)/(1-ecc)))
    e2 = 2 * np.arctan2(np.tan(nu / 2), np.sqrt((1 + ecc) / (1 - ecc)))
    e = e2

    ### longitude of the ascending node Ω (omega_lan)
    if n[1] >= 0:
        omega_lan = np.arccos(n[0] / nnorm)
    else:
        omega_lan = 2 * np.pi - np.arccos(n[0] / nnorm)

    ### argument of periapsis ω (omega_periarg)
    if ecc[2] >= 0:
        omega_periarg = np.arccos((np.dot(n, ecc)) / (nnorm * ecc_norm))
    else:
        omega_periarg = 2 * np.pi - np.arccos((np.dot(n, ecc)) / (nnorm * ecc_norm))

    ### mean anomaly m, (Kepler’s Equation)
    m = e - ecc * np.sin(e)

    ### semi-major axis a
    a = ((2 / pnorm) - (vnorm**2 / mu)) ** -1

    if rad2deg:
        i = np.rad2deg(i)
        omega_periarg = np.rad2deg(omega_periarg)
        omega_lan = np.rad2deg(omega_lan)
        m = np.rad2deg(m)

    return a, ecc, i, omega_periarg, omega_lan, m


def extrapolate_orbit_kepler(p, v, t, t0, mu=3.9860044188e14):
    from pytwobodyorbit import TwoBodyOrbit

    orbit = TwoBodyOrbit("", mu=mu)  # create an instance
    orbit.setOrbCart(t0, p, v)  # define the orbit
    pout, vout = orbit.posvelatt(t)  # get position and velocity at t1
    kepl = orbit.elmKepl()  # get classical orbital elements

    return pout, vout, kepl


def extrapolate_sp3_data_frame(
    df_sp3,
    step=900,
    n_step=9,
    backward=True,
    forward=True,
    until_backward=None,
    until_forward=None,
    return_all=True,
):
    """
    Extrapolate the positions in a SP3 based on the first/last
    epochs' position

    Parameters
    ----------
    df_sp3 : DataFrame
        Input Orbit DataFrame (i.e. generated by files_rw.read_sp3).
    step : int, optional
        step between two epochs. The default is 900.
    n_step : int, optional
        number of epochs to interpolate. The default is 9.
    backward and forward : bool, optional
        extrapolate for the day before/after. The default is True.
    return_all : bool, optional
        if True, returns the input DF enhanced with the extrapolated values
        if False, returns only the extrapolated values.
        The default is True.
    until_backward & until_backward : datetime, optional
        epoch until then the extrapolation has to be done.
        Override n_step

    Returns
    -------
    df_orb_out : DataFrame
        Orbit DataFrame with extrapolated values
        (see return_all).

    """
    new_epoch_stk = []

    for sat in df_sp3["sat"].unique():
        log.info("extrapolate: %s", sat)
        df_sat = df_sp3[(df_sp3["sat"] == sat) & (df_sp3["type"] == "p")].copy()
        df_sat.sort_values("epoch", inplace=True)

        df_line_dummy = df_sat.iloc[0].copy()
        df_line_dummy["clk"] = 999999.999999

        xyz = df_sat[["x", "y", "z"]].values * 1000
        tgps = conv.numpy_dt2dt(df_sat["epoch"].values)
        tutc = conv.dt_gpstime2dt_utc(tgps)
        tsec = [e.total_seconds() for e in (tutc - tutc[0])]

        p = conv.ecef2eci(xyz, tutc)
        v = np.gradient(p, tsec, axis=0)

        def extrapo_intern_fct(i_ref, coef, until):
            ##
            ## backward :
            ## i_ref , coef = 0,-1
            ##
            ## forward :
            ## i_ref , coef = -1,1
            ##
            ## CORRDS ARE GIVEN IN ECI HERE !!!
            ##
            from pytwobodyorbit import TwoBodyOrbit

            orbit_back = TwoBodyOrbit("orbit", mu=3.9860044188e14)  # create an instance
            orbit_back.setOrbCart(tsec[i_ref], p[i_ref], v[i_ref])  # define the orbit

            if until:
                t_rang_strt, t_rang_end = list(sorted([tutc[i_ref], until]))
                range = conv.dt_range(t_rang_strt, t_rang_end, 0, step)
                n_step_intern = len(range) - 1
            else:
                n_step_intern = n_step

            for t in np.arange(step, n_step_intern * step + 1, step):
                pout, vout = orbit_back.posvelatt(coef * t)
                epoc = tgps[i_ref] + coef * dt.timedelta(seconds=int(t))
                df_line_new = df_line_dummy.copy()
                df_line_new[["x", "y", "z"]] = pout
                df_line_new["epoch"] = pd.Timestamp(epoc)

                new_epoch_stk.append(df_line_new)

            return pout, vout

        ### backward
        if backward:
            extrapo_intern_fct(0, -1, until_backward)

        ### forward
        if forward:
            extrapo_intern_fct(-1, 1, until_forward)

    df_new_epoch = pd.DataFrame(new_epoch_stk)

    tutc_new_epoch = conv.dt_gpstime2dt_utc(
        conv.numpy_dt2dt(df_new_epoch["epoch"].values)
    )
    df_new_epoch[["x", "y", "z"]] = conv.eci2ecef(
        df_new_epoch[["x", "y", "z"]].values, tutc_new_epoch
    )
    df_new_epoch[["x", "y", "z"]] = df_new_epoch[["x", "y", "z"]] * 10**-3

    if return_all:
        df_orb_out = pd.concat((df_sp3, df_new_epoch))
    else:
        df_orb_out = df_new_epoch

    df_orb_out.reset_index(drop=True, inplace=True)
    df_orb_out.sort_values(["sat", "epoch"], inplace=True)
    df_orb_out.reset_index(drop=True, inplace=True)

    return df_orb_out
