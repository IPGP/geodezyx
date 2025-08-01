#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 09:23:54 2019

@author: Chaiyaporn Kitpracha
"""

# from geodezyx import *                   # Import the GeodeZYX modules
# from geodezyx.externlib import *         # Import the external modules
# from geodezyx.megalib.megalib import *   # Import the legacy modules names

import numpy as np
from geodezyx import conv, utils


def PWV_conversion(zwd, Tm):
    """
    This function convert from Zenith Wet delay to Precipitate Water Vapor (PWV)

    Parameters
    ----------
    zwd :
        Zenith wet delay in meters
    Tm:
        Mean temperature of troposphere

    Returns
    ----------
    PWV :
        Precipitate Water Vapor in mm

    Sources
    ----------
    Solution and Constant k2' k3 from Atmospheric effects in Space Geodesy Chapter 3.
    """
    if utils.is_iterable(zwd):
        PWV = []
        for zwd_m, Tm_m in zip(zwd, Tm):
            pwv = PWV_conversion(zwd_m, Tm_m)
            PWV.append(pwv)
        return np.array(PWV)

    else:
        k1 = 77.689
        k2 = 71.2952
        Md = 28.965 / 100  # g/mol -> kg/mol
        Mw = 18.016 / 100  # g/mol -> kg/mol
        k2d = k2 - (k1 * (Mw / Md))
        k3 = 375463
        R = 8.3143
        Rhow = 999.97

        CF = 10e6 * Mw / (Rhow * R * (k2d + (k3 / Tm)))
        PWV = CF * zwd * 1000
        return np.round(PWV, 2)


def Tm_bevis(Ts):
    """
    This function determines mean temperature based on surface temperature using Bevis equation

    Parameters
    ----------
    Ts:
        Surface temperature in Kelvin
    Returns
    ----------
    Tm:
        Mean temperature in Kevlin
    Sources
    ----------
    Atmospheric effects in Space Geodesy Chapter 3.
    """

    Tm = 70.2 + 0.72 * Ts
    return np.round(Tm, 2)


def trop_saast(p, dlat, hell, t=0, e=0, mode="dry"):
    """
    This subroutine determines the zenith total delay based on the equation by Saastamoinen (1972) as refined by Davis et al. (1985)

    Parameters
    ----------
    p :
        pressure in hPa

    dlat :
        ellipsoidal latitude in radians

    t :
         temperature in Celcius

    e :
         water vapor pressure in hPa

    hell :
         ellipsoidal height in meters

    mode :
         dry, wet or total

    Returns
    ----------
    res :
        zenith total delay in m (depend on mode)

    Source
    ----------
     c Reference:
     Saastamoinen, J., Atmospheric correction for the troposphere and
     stratosphere in radio ranging of satellites. The use of artificial
     satellites for geodesy, Geophys. Monogr. Ser. 15, Amer. Geophys. Union,
     pp. 274-251, 1972.
     Davis, J.L, T.A. Herring, I.I. Shapiro, A.E.E. Rogers, and G. Elgered,
     Geodesy by Radio Interferometry: Effects of Atmospheric Modeling Errors
     on Estimates of Baseline Length, Radio Science, Vol. 20, No. 6,
     pp. 1593-1607, 1985.

    """
    if utils.is_iterable(dlat):
        ztd = []
        for p_m, dlat_m, hell_m, t_m, e_m in zip(p, dlat, hell, t, e):
            ztd_m = trop_saast(p_m, dlat_m, hell_m, t_m, e_m, mode)
            ztd.append(ztd_m)
        return np.array(ztd)
    else:
        f = (
            1 - 0.00266 * np.cos(2 * dlat) - 0.00000028 * hell
        )  # calculate denominator f
        t = t + 273.15  # convert celcius to kelvin
        if mode == "dry":
            res = 0.0022768 * p / f
        elif mode == "wet":
            res = (0.22768e-2) * ((0.1255e4) + (t + 0.5e-1)) * (e / t)
        elif mode == "total":
            zhd = 0.0022768 * p / f
            zwd = (0.22768e-2) * ((0.1255e4) + (t + 0.5e-1)) * (e / t)
            res = zhd + zwd

        return np.round(res, 4)


def read_grid_gpt(grid_name, cols=64):
    grid = np.loadtxt(grid_name, usecols=range(cols), comments="%", dtype=float)
    return grid


def calc_stand_ties_gpt3(
    epoc, lat_ref, lon_ref, h_ref, lat_rov, lon_rov, h_rov, grid, unit="mm"
):
    """
    Determine standard atmospheric ties from meteological information from GPT3 with analytical equation from Teke et al. (2011)

    Parameters:
    ----------
        epoc :
            time in Python datetime
        lat_ref :
            Latitude of Ref. station
        lat_rov :
            Latitude of Rov. station
        h_ref :
            Height of Ref. station
        h_rov :
            Height of Rov. station
        grid_file :
            meteological grid file
        unit :
            in meters (m) or milimeters (mm)

    Return:
    ----------
        ties :
            Standard ties of total delay in milimeters or meters

    """
    if utils.is_iterable(lat_ref):
        ties = []
        for epoc_m, lat_ref_m, lon_ref_m, h_ref_m, lat_rov_m, lon_rov_m, h_rov_m in zip(
            epoc, lat_ref, lon_ref, h_ref, lat_rov, lon_rov, h_rov
        ):
            ties_m = calc_stand_ties_gpt3(
                epoc_m,
                lat_ref_m,
                lon_ref_m,
                h_ref_m,
                lat_rov_m,
                lon_rov_m,
                h_rov_m,
                grid,
            )
            ties.append(ties_m)
        return np.array(ties)

    else:
        met_ref = gpt3(epoc, lat_ref, lon_ref, h_ref, grid)

        p0, t0, e0 = met_ref[0], met_ref[1], met_ref[4]
        t0 = t0 + 273.15  # convert to Kelvin

        p = p0 * (
            (1 - (0.0065 * (h_rov - h_ref) / t0)) ** (9.80665 / (0.0065 * 287.058))
        )
        dZHD = (0.0022768 * (p0 - p)) / (
            1 - 0.00266 * np.cos(2 * lat_ref) - 0.28e-6 * h_ref
        )
        dZWD = 2.789 * e0 / (t0**2) * ((5383 / t0) - 0.7803) * 0.0065 * (h_rov - h_ref)

        if unit == "m":
            return np.round(dZHD + dZWD, 4)

        elif unit == "mm":
            return np.round((dZHD + dZWD) * 1000, 4)


def calc_stand_ties(epoc, lat_ref, h_ref, h_rov, p0, t0, e0, unit="mm"):
    """
    Determine standard atmospheric ties with analytical equation from Teke et al. (2011)

    Parameters:
    ----------
        epoc :
            time in Python datetime
        lat_ref :
            Latitude of Ref. station
        h_ref :
            Height of Ref. station
        h_rov :
            Height of Rov. station
        p0:
            Pressure of Ref. station in hPa
        t0:
            Temperature in C of Ref. station
        e0:
            Water vapor pressure of Ref. station in hPa
        unit :
            in meters (m) or milimeters (mm)

    Return:
    ----------
        ties :
            Standard ties of total delay in milimeters or meters

    """
    if utils.is_iterable(lat_ref):
        ties = []
        for epoc_m, lat_ref_m, h_ref_m, h_rov_m, p0_m, t0_m, e0_m in zip(
            epoc, lat_ref, h_ref, h_rov, p0, t0, e0
        ):
            ties_m = calc_stand_ties(
                epoc_m, lat_ref_m, h_ref_m, h_rov_m, p0_m, t0_m, e0_m
            )
            ties.append(ties_m)
        return np.array(ties)

    else:
        t0 = t0 + 273.15  # convert to Kelvin

        p = p0 * (
            (1 - (0.0065 * (h_rov - h_ref) / t0)) ** (9.80665 / (0.0065 * 287.058))
        )
        dZHD = (0.0022768 * (p0 - p)) / (
            1 - 0.00266 * np.cos(2 * lat_ref) - 0.28e-6 * h_ref
        )
        dZWD = 2.789 * e0 / (t0**2) * ((5383 / t0) - 0.7803) * 0.0065 * (h_rov - h_ref)

        if unit == "m":
            return np.round(dZHD + dZWD, 4)

        elif unit == "mm":
            return np.round((dZHD + dZWD) * 1000, 4)


def gpt3(dtin, lat, lon, h_ell, C, it=0):
    """
    This subroutine determines pressure, temperature, temperature lapse rate,
    mean temperature of the water vapor, water vapour pressure, hydrostatic
    and wet mapping function coefficients ah and aw, water vapour decrease
    factor, geoid undulation and empirical tropospheric gradients for
    specific sites near the earth's surface.
    It is based on a 5 x 5 degree external grid file ('gpt3_5.grd') with mean
    values as well as sine and cosine amplitudes for the annual and
    semiannual variation of the coefficients.

    Parameters:
    ----------
    dtin :
        datatime in Python datetime object

    lat:
        ellipsoidal latitude in radians [-pi/2:+pi/2]

    lon:
        longitude in radians [-pi:pi] or [0:2pi]

    h_ell:
        ellipsoidal height in m

    it:
        case 1 no time variation but static quantities, case 0 with time variation (annual and semiannual terms)

    Returns:
    ----------
    p:
        pressure in hPa

    T:
        temperature in degrees Celsius

    dT:
        temperature lapse rate in degrees per km

    Tm:
        mean temperature weighted with the water vapor in degrees Kelvin

    e:
        water vapour pressure in hPa

    ah:
        hydrostatic mapping function coefficient at zero height (VMF3)

    aw:
        wet mapping function coefficient (VMF3)

    la:
        water vapour decrease factor

    undu:
        geoid undulation in m

    Gn_h:
        hydrostatic north gradient in m

    Ge_h:
        hydrostatic east gradient in m

    Gn_w:
        wet north gradient in m

    Ge_w:
        wet east gradient in m

    Notes
    ----------
        Modified for Python by Chaiyaporn Kitpracha

    Source
    ----------
        (c) Department of Geodesy and Geoinformation, Vienna University of
        Technology, 2017

        The copyright in this document is vested in the Department of Geodesy and
        Geoinformation (GEO), Vienna University of Technology, Austria. This document
        may only be reproduced in whole or in part, or stored in a retrieval
        system, or transmitted in any form, or by any means electronic,
        mechanical, photocopying or otherwise, either with the prior permission
        of GEO or in accordance with the terms of ESTEC Contract No.
        4000107329/12/NL/LvH.

        D. Landskron, J. Böhm (2018), VMF3/GPT3: Refined Discrete and Empirical Troposphere Mapping Functions,
        J Geod (2018) 92: 349., doi: 10.1007/s00190-017-1066-2.
        Download at: https://link.springer.com/content/pdf/10.1007%2Fs00190-017-1066-2.pdf
    """
    lat = np.array(lat)
    lon = np.array(lon)
    h_ell = np.array(h_ell)
    # Extract data from grid
    p_grid = C[:, 2:7]  # pressure in Pascal
    T_grid = C[:, 7:12]  # temperature in Kelvin
    Q_grid = C[:, 12:17] / 1000  # specific humidity in kg/kg
    dT_grid = C[:, 17:22] / 1000  # temperature lapse rate in Kelvin/m
    u_grid = C[:, 22]  # geoid undulation in m
    Hs_grid = C[:, 23]  # orthometric grid height in m
    ah_grid = (
        C[:, 24:29] / 1000
    )  # hydrostatic mapping function coefficient, dimensionless
    aw_grid = C[:, 29:34] / 1000  # wet mapping function coefficient, dimensionless
    la_grid = C[:, 34:39]  # water vapor decrease factor, dimensionless
    Tm_grid = C[:, 39:44]  # mean temperature in Kelvin
    Gn_h_grid = C[:, 44:49] / 100000  # hydrostatic north gradient in m
    Ge_h_grid = C[:, 49:54] / 100000  # hydrostatic east gradient in m
    Gn_w_grid = C[:, 54:59] / 100000  # wet north gradient in m
    Ge_w_grid = C[:, 59:64] / 100000  # wet east gradient in m

    # Convert from datetime to doy
    doy = float(conv.dt2doy(dtin)) + conv.dt2fracday(dtin)

    # determine the GPT3 coefficients

    # mean gravity in m/s**2
    gm = 9.80665
    # molar mass of dry air in kg/mol
    dMtr = 28.965e-3
    # universal gas constant in J/K/mol
    Rg = 8.3143

    # factors for amplitudes
    if it == 1:  # then  constant parameters
        cosfy = 0
        coshy = 0
        sinfy = 0
        sinhy = 0
    else:
        cosfy = np.cos(doy / 365.25 * 2 * np.pi)  # coefficient for A1
        coshy = np.cos(doy / 365.25 * 4 * np.pi)  # coefficient for B1
        sinfy = np.sin(doy / 365.25 * 2 * np.pi)  # coefficient for A2
        sinhy = np.sin(doy / 365.25 * 4 * np.pi)  # coefficient for B2

    nstat = lat.size

    # initialization
    p = np.zeros([nstat, 1])
    T = np.zeros([nstat, 1])
    dT = np.zeros([nstat, 1])
    Tm = np.zeros([nstat, 1])
    e = np.zeros([nstat, 1])
    ah = np.zeros([nstat, 1])
    aw = np.zeros([nstat, 1])
    la = np.zeros([nstat, 1])
    undu = np.zeros([nstat, 1])
    Gn_h = np.zeros([nstat, 1])
    Ge_h = np.zeros([nstat, 1])
    Gn_w = np.zeros([nstat, 1])
    Ge_w = np.zeros([nstat, 1])

    if lon < 0:
        plon = (lon + 2 * np.pi) * 180 / np.pi
    else:
        plon = lon * 180 / np.pi

    ppod = (-lat + np.pi / 2) * 180 / np.pi

    ipod = np.floor(ppod + 1)
    ilon = np.floor(plon + 1)

    # changed for the 1 degree grid
    diffpod = ppod - (ipod - 0.5)
    difflon = plon - (ilon - 0.5)

    if ipod == 181:
        ipod = 180

    if ilon == 361:
        ilon = 1

    if ilon == 0:
        ilon = 360

    indx = np.zeros(4)
    indx[0] = (ipod - 1) * 360 + ilon

    # near the poles: nearest neighbour interpolation, otherwise: bilinear
    # with the 1 degree grid the limits are lower and upper
    bilinear = 0
    if ppod > 0.5 and ppod < 179.5:
        bilinear = 1

    if bilinear == 0:
        ix = int(indx[0]) - 1

        # transforming ellipsoidal height to orthometric height
        undu = u_grid[ix]
        hgt = h_ell - undu

        # pressure, temperature at the height of the grid
        T0 = (
            T_grid[ix, 0]
            + T_grid[ix, 1] * cosfy
            + T_grid[ix, 2] * sinfy
            + T_grid[ix, 3] * coshy
            + T_grid[ix, 4] * sinhy
        )
        p0 = (
            p_grid[ix, 0]
            + p_grid[ix, 1] * cosfy
            + p_grid[ix, 2] * sinfy
            + p_grid[ix, 3] * coshy
            + p_grid[ix, 4] * sinhy
        )

        # specific humidity
        Q = (
            Q_grid[ix, 0]
            + Q_grid[ix, 1] * cosfy
            + Q_grid[ix, 2] * sinfy
            + Q_grid[ix, 3] * coshy
            + Q_grid[ix, 4] * sinhy
        )

        # lapse rate of the temperature
        dT = (
            dT_grid[ix, 0]
            + dT_grid[ix, 1] * cosfy
            + dT_grid[ix, 2] * sinfy
            + dT_grid[ix, 3] * coshy
            + dT_grid[ix, 4] * sinhy
        )

        # station height - grid height
        redh = hgt - Hs_grid[ix]

        # temperature at station height in Celsius
        T = T0 + dT * redh - 273.15

        # temperature lapse rate in degrees / km
        dT = dT * 1000

        # virtual temperature in Kelvin
        Tv = T0 * [1 + 0.6077 * Q]

        c = gm * dMtr / [Rg * Tv]

        # pressure in hPa
        p = [p0 * np.exp[-c * redh]] / 100

        # hydrostatic and wet coefficients ah and aw
        ah = (
            ah_grid[ix, 0]
            + ah_grid[ix, 2] * cosfy
            + ah_grid[ix, 3] * sinfy
            + ah_grid[ix, 4] * coshy
            + ah_grid[ix, 5] * sinhy
        )
        aw = (
            aw_grid[ix, 0]
            + aw_grid[ix, 2] * cosfy
            + aw_grid[ix, 3] * sinfy
            + aw_grid[ix, 4] * coshy
            + aw_grid[ix, 5] * sinhy
        )

        # water vapour decrease factor la
        la = (
            la_grid[ix, 0]
            + la_grid[ix, 1] * cosfy
            + la_grid[ix, 2] * sinfy
            + la_grid[ix, 3] * coshy
            + la_grid[ix, 4] * sinhy
        )

        # mean temperature Tm
        Tm = (
            Tm_grid[ix, 0]
            + Tm_grid[ix, 1] * cosfy
            + Tm_grid[ix, 2] * sinfy
            + Tm_grid[ix, 3] * coshy
            + Tm_grid[ix, 4] * sinhy
        )

        # north and east gradients [total, hydrostatic and wet]
        Gn_h = (
            Gn_h_grid[ix, 0]
            + Gn_h_grid[ix, 1] * cosfy
            + Gn_h_grid[ix, 2] * sinfy
            + Gn_h_grid[ix, 3] * coshy
            + Gn_h_grid[ix, 4] * sinhy
        )
        Ge_h = (
            Ge_h_grid[ix, 0]
            + Ge_h_grid[ix, 1] * cosfy
            + Ge_h_grid[ix, 2] * sinfy
            + Ge_h_grid[ix, 3] * coshy
            + Ge_h_grid[ix, 4] * sinhy
        )
        Gn_w = (
            Gn_w_grid[ix, 0]
            + Gn_w_grid[ix, 1] * cosfy
            + Gn_w_grid[ix, 2] * sinfy
            + Gn_w_grid[ix, 3] * coshy
            + Gn_w_grid[ix, 4] * sinhy
        )
        Ge_w = (
            Ge_w_grid[ix, 0]
            + Ge_w_grid[ix, 1] * cosfy
            + Ge_w_grid[ix, 2] * sinfy
            + Ge_w_grid[ix, 3] * coshy
            + Ge_w_grid[ix, 4] * sinhy
        )

        # water vapor pressure in hPa
        e0 = Q * p0 / [0.622 + 0.378 * Q] / 100  # on the grid
        e = e0 * (100 * p / p0) ** (
            la + 1
        )  # on the station height - (14] Askne and Nordius, 1987

    else:
        ipod1 = ipod + 1 * np.sign(diffpod)
        ilon1 = ilon + 1 * np.sign(difflon)

        # changed for the 1 degree grid
        if ilon1 == 361:
            ilon1 = 1

        if ilon1 == 0:
            ilon1 = 360

        # get the number of the line
        # changed for the 1 degree grid
        indx[1] = (ipod1 - 1) * 360 + ilon  # along same longitude
        indx[2] = (ipod - 1) * 360 + ilon1  # along same polar distance
        indx[3] = (ipod1 - 1) * 360 + ilon1  # diagonal
        indx = indx.astype(int)
        indx = indx - 1

        # transforming ellipsoidal height to orthometric height: Hortho = -N + Hell
        undul = u_grid[indx]
        hgt = h_ell - undul

        # pressure, temperature at the height of the grid
        T0 = (
            T_grid[indx, 0]
            + T_grid[indx, 1] * cosfy
            + T_grid[indx, 2] * sinfy
            + T_grid[indx, 3] * coshy
            + T_grid[indx, 4] * sinhy
        )
        p0 = (
            p_grid[indx, 0]
            + p_grid[indx, 1] * cosfy
            + p_grid[indx, 2] * sinfy
            + p_grid[indx, 3] * coshy
            + p_grid[indx, 4] * sinhy
        )

        # humidity
        Ql = (
            Q_grid[indx, 0]
            + Q_grid[indx, 1] * cosfy
            + Q_grid[indx, 2] * sinfy
            + Q_grid[indx, 3] * coshy
            + Q_grid[indx, 4] * sinhy
        )

        # reduction = stationheight - gridheight
        Hs1 = Hs_grid[indx]
        redh = hgt - Hs1

        # lapse rate of the temperature in degree / m
        dTl = (
            dT_grid[indx, 0]
            + dT_grid[indx, 1] * cosfy
            + dT_grid[indx, 2] * sinfy
            + dT_grid[indx, 3] * coshy
            + dT_grid[indx, 4] * sinhy
        )

        # temperature reduction to station height
        Tl = T0 + dTl * redh - 273.15

        # virtual temperature
        Tv = T0 * (1 + 0.6077 * Ql)
        c = gm * dMtr / (Rg * Tv)

        # pressure in hPa
        pl = (p0 * np.exp(-c * redh)) / 100

        # hydrostatic and wet coefficients ah and aw
        ahl = (
            ah_grid[indx, 0]
            + ah_grid[indx, 1] * cosfy
            + ah_grid[indx, 2] * sinfy
            + ah_grid[indx, 3] * coshy
            + ah_grid[indx, 4] * sinhy
        )
        awl = (
            aw_grid[indx, 0]
            + aw_grid[indx, 1] * cosfy
            + aw_grid[indx, 2] * sinfy
            + aw_grid[indx, 3] * coshy
            + aw_grid[indx, 4] * sinhy
        )

        # water vapour decrease factor la
        lal = (
            la_grid[indx, 0]
            + la_grid[indx, 1] * cosfy
            + la_grid[indx, 2] * sinfy
            + la_grid[indx, 3] * coshy
            + la_grid[indx, 4] * sinhy
        )

        # mean temperature of the water vapor Tm
        Tml = (
            Tm_grid[indx, 0]
            + Tm_grid[indx, 1] * cosfy
            + Tm_grid[indx, 2] * sinfy
            + Tm_grid[indx, 3] * coshy
            + Tm_grid[indx, 4] * sinhy
        )

        # north and east gradients [total, hydrostatic and wet]
        Gn_hl = (
            Gn_h_grid[indx, 0]
            + Gn_h_grid[indx, 1] * cosfy
            + Gn_h_grid[indx, 2] * sinfy
            + Gn_h_grid[indx, 3] * coshy
            + Gn_h_grid[indx, 4] * sinhy
        )
        Ge_hl = (
            Ge_h_grid[indx, 0]
            + Ge_h_grid[indx, 1] * cosfy
            + Ge_h_grid[indx, 2] * sinfy
            + Ge_h_grid[indx, 3] * coshy
            + Ge_h_grid[indx, 4] * sinhy
        )
        Gn_wl = (
            Gn_w_grid[indx, 0]
            + Gn_w_grid[indx, 1] * cosfy
            + Gn_w_grid[indx, 2] * sinfy
            + Gn_w_grid[indx, 3] * coshy
            + Gn_w_grid[indx, 4] * sinhy
        )
        Ge_wl = (
            Ge_w_grid[indx, 0]
            + Ge_w_grid[indx, 1] * cosfy
            + Ge_w_grid[indx, 2] * sinfy
            + Ge_w_grid[indx, 3] * coshy
            + Ge_w_grid[indx, 4] * sinhy
        )

        # water vapor pressure in hPa
        e0 = Ql * p0 / (0.622 + 0.378 * Ql) / 100  # on the grid
        el = e0 * (100 * pl / p0) ** (
            lal + 1
        )  # on the station height - [14] Askne and Nordius, 1987

        dnpod1 = abs(diffpod)  # distance nearer point
        dnpod2 = 1 - dnpod1  # distance to distant point
        dnlon1 = abs(difflon)
        dnlon2 = 1 - dnlon1

        # pressure
        R1 = dnpod2 * pl[0] + dnpod1 * pl[1]
        R2 = dnpod2 * pl[2] + dnpod1 * pl[3]
        p = dnlon2 * R1 + dnlon1 * R2

        # temperature
        R1 = dnpod2 * Tl[0] + dnpod1 * Tl[1]
        R2 = dnpod2 * Tl[2] + dnpod1 * Tl[3]
        T = dnlon2 * R1 + dnlon1 * R2

        # temperature in degree per km
        R1 = dnpod2 * dTl[0] + dnpod1 * dTl[1]
        R2 = dnpod2 * dTl[2] + dnpod1 * dTl[3]
        dT = (dnlon2 * R1 + dnlon1 * R2) * 1000

        # water vapor pressure in hPa
        R1 = dnpod2 * el[0] + dnpod1 * el[1]
        R2 = dnpod2 * el[2] + dnpod1 * el[3]
        e = dnlon2 * R1 + dnlon1 * R2

        # ah and aw
        R1 = dnpod2 * ahl[0] + dnpod1 * ahl[1]
        R2 = dnpod2 * ahl[2] + dnpod1 * ahl[3]
        ah = dnlon2 * R1 + dnlon1 * R2
        R1 = dnpod2 * awl[0] + dnpod1 * awl[1]
        R2 = dnpod2 * awl[2] + dnpod1 * awl[3]
        aw = dnlon2 * R1 + dnlon1 * R2

        # undulation
        R1 = dnpod2 * undul[0] + dnpod1 * undul[1]
        R2 = dnpod2 * undul[2] + dnpod1 * undul[3]
        undu = dnlon2 * R1 + dnlon1 * R2

        # water vapor decrease factor la
        R1 = dnpod2 * lal[0] + dnpod1 * lal[1]
        R2 = dnpod2 * lal[2] + dnpod1 * lal[3]
        la = dnlon2 * R1 + dnlon1 * R2

        # gradients
        R1 = dnpod2 * Gn_hl[0] + dnpod1 * Gn_hl[1]
        R2 = dnpod2 * Gn_hl[2] + dnpod1 * Gn_hl[3]
        Gn_h = dnlon2 * R1 + dnlon1 * R2
        R1 = dnpod2 * Ge_hl[0] + dnpod1 * Ge_hl[1]
        R2 = dnpod2 * Ge_hl[2] + dnpod1 * Ge_hl[3]
        Ge_h = dnlon2 * R1 + dnlon1 * R2
        R1 = dnpod2 * Gn_wl[0] + dnpod1 * Gn_wl[1]
        R2 = dnpod2 * Gn_wl[2] + dnpod1 * Gn_wl[3]
        Gn_w = dnlon2 * R1 + dnlon1 * R2
        R1 = dnpod2 * Ge_wl[0] + dnpod1 * Ge_wl[1]
        R2 = dnpod2 * Ge_wl[2] + dnpod1 * Ge_wl[3]
        Ge_w = dnlon2 * R1 + dnlon1 * R2

        # mean temperature of the water vapor Tm
        R1 = dnpod2 * Tml[0] + dnpod1 * Tml[1]
        R2 = dnpod2 * Tml[2] + dnpod1 * Tml[3]
        Tm = dnlon2 * R1 + dnlon1 * R2

    soln = [
        np.round(p, 3),
        np.round(T, 3),
        np.round(dT, 3),
        np.round(Tm, 3),
        np.round(e, 3),
        np.round(ah, 3),
        np.round(aw, 3),
        np.round(la, 3),
        np.round(undu, 3),
        np.round(Gn_h, 3),
        np.round(Ge_h, 3),
        np.round(Gn_w, 3),
        np.round(Ge_w, 3),
    ]
    return soln


def gpt2_5(mjd, lat, lon, HELL, IT, VEC):
    """
     (c) Department of Geodesy and Geoinformation, Vienna University of
     Technology, 2013

     The copyright in this document is vested in the Department of Geodesy and
     Geoinformation (GEO), Vienna University of Technology, Austria. This document
     may only be reproduced in whole or in part, or stored in a retrieval
     system, or transmitted in any form, or by any means electronic,
     mechanical, photocopying or otherwise, either with the prior permission
     of GEO or in accordance with the terms of ESTEC Contract No.
     4000107329/12/NL/LvH.
     ---

     This subroutine determines pressure, temperature, temperature lapse rate,
     mean temperature of the water vapor, water vapour pressure, hydrostatic
     and wet mapping function coefficients ah and aw, water vapour decrease
     factor and geoid undulation for specific sites near the Earth surface.
     It is based on a 5 x 5 degree external grid file ('gpt2_5.grd') with mean
     values as well as sine and cosine amplitudes for the annual and
     semiannual variation of the coefficients.

     The hydrostatic mapping function coefficients have to be used with the
     height dependent Vienna Mapping Function 1 (vmf_ht.f) because the
     coefficients refer to zero height.

     Example 1 (Vienna, 2 August 2012, with time variation):

     dmjd = 56141.d0
     dlat(1) = 48.20d0*pi/180.d0
     dlon(1) = 16.37d0*pi/180.d0
     hell(1) = 156.d0
     nstat = 1
     it = 0

     output:
     p = 1002.56 hPa
     T = 22.12 deg Celsius
     dT = -6.53 deg / km
     Tm = 281.11 K
     e = 16.72 hPa
     ah = 0.0012647
     aw = 0.0005726
     la = 2.6964
     undu = 44.06 m

     Example 2 (Vienna, 2 August 2012, without time variation, i.e. constant values):

     dmjd = 56141.d0
     dlat(1) = 48.20d0*pi/180.d0
     dlon(1) = 16.37d0*pi/180.d0
     hell(1) = 156.d0
     nstat = 1
     it = 1

     output:
     p = 1003.49 hPa
     T = 11.95 deg Celsius
     dT = -5.47 deg / km
     Tm = 273.00 K
     e = 10.23 hPa
     ah = 0.0012395
     aw = 0.0005560
     la = 2.6649
     undu = 44.06 m

     Klemens Lagler, 2 August 2012
     Johannes Boehm, 6 August 2012, revision
     Klemens Lagler, 21 August 2012, epoch change to January 1 2000
     Johannes Boehm, 23 August 2012, adding possibility to determine constant field
     Johannes Boehm, 27 December 2012, reference added
     Johannes Boehm, 10 January 2013, correction for dlat = -90 degrees
      (problem found by Changyong He)
     Johannes Boehm, 21 May 2013, bug with dmjd removed (input parameter dmjd was replaced
     unintentionally; problem found by Dennis Ferguson)
     Gregory Pain,   17 June 2013, adding water vapour decrease factor la
     Gregory Pain,   01 July 2013, adding mean temperature Tm
     Gregory Pain,   30 July 2013, changing the method to calculate the water vapor partial pressure (e)
     Gregory Pain,   31 July 2013, correction for (dlat = -90 degrees, dlon = 360 degrees)
     Johannes Boehm, 27 December 2013, copyright notice added
     Johannes Boehm, 25 August 2014, reference changed to Boehm et al. in GPS
     Solutions

    Source
    ----------
         J. Böhm, G. Möller, M. Schindelegger, G. Pain, R. Weber, Development of an
         improved blind model for slant delays in the troposphere (GPT2w),
         GPS Solutions, 2014, doi:10.1007/s10291-014-0403-7
    Notes:
         Modified for Python by: Chaiyaporn Kitpracha

    Parameters:
    ----------
    mjd:
         modified Julian date (scalar, only one epoch per call is possible)

    lat:
         ellipsoidal latitude in degrees

    lon:
         longitude in degrees

    HELL:
         ellipsoidal height in m

    IT:
         case 1: no time variation but static quantities
    case 0: with time variation (annual and semiannual terms)

    VEC:
        GPT2 grid data in numpy array size 5 x 5 degree ('gpt2_5.grd')

    Returns:
    ----------
     p:
         pressure in hPa

     T:
         temperature in degrees Celsius

     dT:
         temperature lapse rate in degrees per km

     e:
         water vapour pressure in hPa

     undu:
         geoid undulation in m

    """
    GM = 9.80665
    DMTR = 28.965e-3
    RG = 8.3143

    DLAT = lat * np.pi / 180
    DLON = lon * np.pi / 180

    DMJD1 = mjd - 51544.5
    if IT == 1:
        COSFY = 0
        COSHY = 0
        SINFY = 0
        SINHY = 0
    else:
        COSFY = np.cos(DMJD1 / 365.25 * 2 * np.pi)
        COSHY = np.cos(DMJD1 / 365.25 * 4 * np.pi)
        SINFY = np.sin(DMJD1 / 365.25 * 2 * np.pi)
        SINHY = np.sin(DMJD1 / 365.25 * 4 * np.pi)

    PGRID = VEC[0:2592, 2:7]
    TGRID = VEC[0:2592, 7:12]
    QGRID = VEC[0:2592, 12:17] / 1000
    DTGRID = VEC[0:2592, 17:22] / 1000
    U = VEC[0:2592, 22]
    HS = VEC[0:2592, 23]

    if DLON < 0:
        PLON = (DLON + 2 * np.pi) * 180 / np.pi
    else:
        PLON = DLON * 180 / np.pi

    PPOD = (-DLAT + np.pi / 2) * 180 / np.pi

    IPOD = np.floor((PPOD + 5) / 5)
    ILON = np.floor((PLON + 5) / 5)

    DIFFPOD = (PPOD - (IPOD * 5 - 2.5)) / 5
    DIFFLON = (PLON - (ILON * 5 - 2.5)) / 5

    if IPOD == 37:
        IPOD = 36

    INDX = np.zeros(4)
    INDX[0] = (IPOD - 1) * 72 + ILON
    IBILINEAR = 0

    if PPOD > 2.5 and PPOD < 177.5:
        IBILINEAR = 1

    if IBILINEAR == 0:
        UNDU = 0.0
        IX = INDX[0] - 1
        UNDU = U[IX]
        HGT = HELL - UNDU

        T0 = (
            TGRID[IX, 0]
            + TGRID[IX, 1] * COSFY
            + TGRID[IX, 2] * SINFY
            + TGRID[IX, 3] * COSHY
            + TGRID[IX, 4] * SINHY
        )
        P0 = (
            PGRID[IX, 0]
            + PGRID[IX, 1] * COSFY
            + PGRID[IX, 2] * SINFY
            + PGRID[IX, 3] * COSHY
            + PGRID[IX, 4] * SINHY
        )
        Q = (
            QGRID[IX, 0]
            + QGRID[IX, 1] * COSFY
            + QGRID[IX, 2] * SINFY
            + QGRID[IX, 3] * COSHY
            + QGRID[IX, 4] * SINHY
        )
        DT = (
            DTGRID[IX, 0]
            + DTGRID[IX, 1] * COSFY
            + DTGRID[IX, 2] * SINFY
            + DTGRID[IX, 3] * COSHY
            + DTGRID[IX, 4] * SINHY
        )

        REDH = HGT - HS[IX]
        T = T0 + DT * REDH - 273.15
        DT = DT * 1000
        TV = T0 * (1 + 0.6077 * Q)
        C = GM * DMTR / (RG * TV)
        P = (P0 * np.exp(-C * REDH)) / 100
        E = (Q * P) / (0.622 + 0.378 * Q)
    else:
        IPOD1 = IPOD + 1 * np.sign(DIFFPOD)
        ILON1 = ILON + 1 * np.sign(DIFFLON)
        if ILON1 == 73:
            ILON1 = 1

        if ILON1 == 0:
            ILON1 = 72

        INDX[1] = (IPOD1 - 1) * 72 + ILON
        INDX[2] = (IPOD - 1) * 72 + ILON1
        INDX[3] = (IPOD1 - 1) * 72 + ILON1

        INDX = INDX.astype(int)
        UNDUL = np.zeros(4)
        QL = np.zeros(4)
        DTL = np.zeros(4)
        TL = np.zeros(4)
        PL = np.zeros(4)

        INDX = INDX - 1

        for L in range(0, 4):
            UNDUL[L] = U[INDX[L]]
            HGT = HELL - UNDUL[L]
            T0 = (
                TGRID[INDX[L], 0]
                + TGRID[INDX[L], 1] * COSFY
                + TGRID[INDX[L], 2] * SINFY
                + TGRID[INDX[L], 3] * COSHY
                + TGRID[INDX[L], 4] * SINHY
            )
            P0 = (
                PGRID[INDX[L], 0]
                + PGRID[INDX[L], 1] * COSFY
                + PGRID[INDX[L], 2] * SINFY
                + PGRID[INDX[L], 3] * COSHY
                + PGRID[INDX[L], 4] * SINHY
            )
            QL[L] = (
                QGRID[INDX[L], 0]
                + QGRID[INDX[L], 1] * COSFY
                + QGRID[INDX[L], 2] * SINFY
                + QGRID[INDX[L], 3] * COSHY
                + QGRID[INDX[L], 4] * SINHY
            )
            HS1 = HS[INDX[L]]
            REDH = HGT - HS1

            DTL[L] = (
                DTGRID[INDX[L], 0]
                + DTGRID[INDX[L], 1] * COSFY
                + DTGRID[INDX[L], 2] * SINFY
                + DTGRID[INDX[L], 3] * COSHY
                + DTGRID[INDX[L], 4] * SINHY
            )

            TL[L] = T0 + DTL[L] * REDH - 273.15
            TV = T0 * (1 + 0.6077 * QL[L])
            C = GM * DMTR / (RG * TV)
            PL[L] = (P0 * np.exp(-C * REDH)) / 100

        DNPOD1 = np.absolute(DIFFPOD)
        DNPOD2 = 1 - DNPOD1
        DNLON1 = np.absolute(DIFFLON)
        DNLON2 = 1 - DNLON1

        R1 = DNPOD2 * PL[0] + DNPOD1 * PL[1]
        R2 = DNPOD2 * PL[2] + DNPOD1 * PL[3]
        P = DNLON2 * R1 + DNLON1 * R2

        R1 = DNPOD2 * TL[0] + DNPOD1 * TL[1]
        R2 = DNPOD2 * TL[2] + DNPOD1 * TL[3]
        T = DNLON2 * R1 + DNLON1 * R2

        R1 = DNPOD2 * DTL[0] + DNPOD1 * DTL[1]
        R2 = DNPOD2 * DTL[2] + DNPOD1 * DTL[3]
        DT = (DNLON2 * R1 + DNLON1 * R2) * 1000

        R1 = DNPOD2 * QL[0] + DNPOD1 * QL[1]
        R2 = DNPOD2 * QL[2] + DNPOD1 * QL[3]
        Q = DNLON2 * R1 + DNLON1 * R2
        E = (Q * P) / (0.622 + 0.378 * Q)

        R1 = DNPOD2 * UNDUL[0] + DNPOD1 * UNDUL[1]
        R2 = DNPOD2 * UNDUL[2] + DNPOD1 * UNDUL[3]
        UNDU = DNLON2 * R1 + DNLON1 * R2
    return round(P, 2), round(T, 2), round(DT, 2), round(E, 2), round(UNDU, 2)


def vmf1(ah, aw, dt, dlat, zd):
    """
    This subroutine determines the VMF1 (Vienna Mapping Functions 1) for specific sites.

    Parameters
    ----------
    ah:
        hydrostatic coefficient a
    aw:
        wet coefficient a
    dt:
        datetime in python datetime
    dlat:
        ellipsoidal latitude in radians
    zd:
        zenith distance in radians

    Return
    ----------
    vmf1h:
        hydrostatic mapping function
    vmf1w:
        wet mapping function

    Reference
    ----------
    Boehm, J., B. Werl, H. Schuh (2006), Troposphere mapping functions for GPS and very long baseline interferometry
    from European Centre for Medium-Range Weather Forecasts operational analysis data,
    J. Geoph. Res., Vol. 111, B02406, doi:10.1029/2005JB003629.

    Notes
    ----------
    Written by Johannes Boehm, 2005 October 2

    Translated to python by Chaiyaporn Kitpracha
    """

    pi = 3.14159265359
    dmjd = conv.dt2mjd(dt)
    doy = dmjd - 44239.0 + 1 - 28

    bh = 0.0029
    c0h = 0.062

    if dlat < 0:
        phh = pi
        c11h = 0.007
        c10h = 0.002
    else:
        phh = 0
        c11h = 0.005
        c10h = 0.001

    ch = c0h + ((np.cos(doy / 365.25 * 2 * pi + phh) + 1) * c11h / 2 + c10h) * (
        1 - np.cos(dlat)
    )

    sine = np.sin(pi / 2 - zd)
    beta = bh / (sine + ch)
    gamma = ah / (sine + beta)
    topcon = 1 + ah / (1 + bh / (1 + ch))
    vmf1h = topcon / (sine + gamma)

    bw = 0.00146
    cw = 0.04391
    beta = bw / (sine + cw)
    gamma = aw / (sine + beta)
    topcon = 1 + aw / (1 + bw / (1 + cw))
    vmf1w = topcon / (sine + gamma)

    return vmf1h, vmf1w


####################### Compare troposphere delay  ####################################################
