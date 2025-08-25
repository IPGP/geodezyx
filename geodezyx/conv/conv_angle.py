#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 15:11:12 2020

@author: psakic

This sub-module of geodezyx.conv deals with angle conversion.

it can be imported directly with:
from geodezyx import conv

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""

#### Import the logger
import logging
import re

########## BEGIN IMPORT ##########
#### External modules
import numpy as np

#### geodeZYX modules
from geodezyx import utils

log = logging.getLogger("geodezyx")

##########  END IMPORT  ##########

#                       _         _____                              _
#     /\               | |       / ____|                            (_)
#    /  \   _ __   __ _| | ___  | |     ___  _ ____   _____ _ __ ___ _  ___  _ __
#   / /\ \ | '_ \ / _` | |/ _ \ | |    / _ \| '_ \ \ / / _ \ '__/ __| |/ _ \| '_ \
#  / ____ \| | | | (_| | |  __/ | |___| (_) | | | \ v /  __/ |  \__ \ | (_) | | | |
# /_/    \_\_| |_|\__, |_|\___|  \_____\___/|_| |_|\_/ \___|_|  |___/_|\___/|_| |_|
#                  __/ |
#                 |___/


### Angle conversion
def degdec2dms(deg_in, only_dm=False):
    """
    Angle representation conversion

    Convert :
    Decimal Degree => Degree Minute Seconds

    See also dms_num2str for Numeric => String representation conversion

    Parameters
    ----------
    deg_in : float or iterable of floats
        Decimal degrees

    only_dm : bool
        True to get the output string as a Degree Minute Angle only

    Returns
    -------
    deg , minu , sec : np.array of float
        3 arrays for Degrees, Minutes, Seconds

    """
    # bring everything to positive
    sign = np.sign(deg_in)
    deg_pos = np.abs(deg_in)

    deg = np.floor(deg_pos)
    decimal_part = deg_pos - deg
    decimal_part_sec = decimal_part * 3600
    minu = np.floor_divide(decimal_part_sec, 60)
    sec = decimal_part_sec - minu * 60
    sec = np.round(sec, 10)  # ROUNDING is it necessary?

    # back to the right hemisphere
    deg = sign * deg

    if only_dm:
        return deg, minu + sec * (1.0 / 60.0)
    else:
        return deg, minu, sec


def dms_num2str(dms_inp, decimal_fmt="6.3f"):
    """
    Angle representation conversion

    Convert :
    numeric Degree-Minute-Second (DMS) angle => a string representation.

    Parameters
    ----------
    dms_inp : list or tuple of floats
        A list or tuple containing the DMS values. The first element is degrees,
        the second element (optional) is minutes, and the third element (optional) is seconds.
    decimal_fmt : str, optional
        The format for the seconds part of the DMS string. Default is '7.4f'.

    Returns
    -------
    out_str : str
        A string representation of the DMS angle in the format "DD°MM'SS.SSSS\"".
        If only degrees and minutes are provided, the format will be "DD°MM'".
    """

    out_str = "{:02d}".format(int(dms_inp[0])) + "°"

    decimal_fmt_use = "{:0" + decimal_fmt + "}"

    if len(dms_inp) == 3:
        out_str = (
            out_str
            + "{:02d}".format(int(dms_inp[1]))
            + "'"
            + decimal_fmt_use.format(dms_inp[2])
            + '"'
        )
    elif len(dms_inp) == 2:
        out_str = out_str + decimal_fmt_use.format(dms_inp[1]) + "'"

    return out_str


def dms2degdec_num(deg, minn=0, sec=0):
    """
    Angle representation conversion

    Convert :
    Degree Minute Second `float` Angle => decimal Degree `float` Angle

    Parameters
    ----------
    deg & minn & sec : float
        degres, minutes and seconds of the input Angle

    Returns
    -------
    dd_float : float
        Decimal degree Angle

    """
    sig = np.sign(deg)

    return deg + sig * minn * (1.0 / 60.0) + sig * sec * (1.0 / 3600.0)


def dms2degdec_str(dms_str, only_dm=False):
    """
    Angle representation conversion

    Convert :
    Degree Minute Second `string` Angle => decimal Degree `float` Angle

    can manage only DM angle in input

    NB : Too Complicated .... must be simplified

    Parameters
    ----------
    dms_str : str
        string of DMS angle
        e.g.
        "2°20'35.09"E"

    only_dm : bool
        True if the string is only a Degree Minute Angle
        e.g.
        "40°52.0931'N"
        "28°31.4136'E"

    Returns
    -------
    dd_float : float
        Decimal degree Angle

    """

    if not only_dm:
        log.warning("DMS mode not well implemented yet !!! ")

    dms_str = dms_str.strip()

    if re.match(".*[swoSWO].*", dms_str):
        sign = -1
    else:
        sign = 1

    # former regex
    lis = re.findall(r"\d+|\D+", dms_str)

    ipt = lis.index(".")
    decimal = lis[ipt - 1] + lis[ipt] + lis[ipt + 1]
    lis[ipt - 1] = decimal
    lis = lis[:ipt]

    for i, e in enumerate(lis):
        try:
            lis[i] = float(e)
        except:
            lis.remove(e)

    lis = [float(e) for e in lis]

    deg = lis[0]
    minu = lis[1]
    if only_dm:
        sec = 0.0
    else:
        try:
            sec = lis[2]
        except IndexError:
            log.error("did you forgot to activate the DM only mode ? ;) ?")
            return None

    dd_float = sign * (deg + minu / 60.0 + sec / 3600.0)
    return dd_float


###### arcsec


def arcsec2deg(arcsec_in):
    """
    Angle conversion

    Convert :
    Arcsecond => Degrees

    NB : Not really useful, should be avoided

    """

    return arcsec_in / 3600.0


def deg2arcsec(deg_in):
    """
    Angle conversion

    Convert :
    Degrees => Arcsecond

    NB : Not really useful, should be avoided

    """
    return deg_in * 3600.0


def angle2equivalent_earth_radius(angle_in, angtype="deg", earth_radius=6371008.8):
    """
    Quick and simple function which gives the equivalent distance on a
    spherical Earth great circle of an angle

    Useful to determine metric varations in latitude

    angle can be : "deg", "rad", "mas"

    """
    earth_circum = earth_radius * 2 * np.pi
    if angtype == "deg":
        equiv_out = (angle_in * earth_circum) / 360.0
    elif angtype == "rad":
        equiv_out = (angle_in * earth_circum) / (np.pi * 2)
    elif angtype == "mas":
        equiv_out = (angle_in * 10**-3 * earth_circum) / 86400.0
    else:
        log.error("bad angtype %s", angtype)
        return None

    return equiv_out


def angle2equivalent_earth_parallel(
    angle_in, latitude_in, angtype="deg", earth_radius=6371008.8
):
    """
    Quick and simple function which gives the equivalent distance on an
    Earth parallel circle of an angle

    Useful to determine metric varaitions in longitude
    (but the latitude is also necessary to determine
     the radius of the parallel)

    angle can be : "deg", "rad", "mas"

    """

    if angtype == "deg":
        latitude_rad = np.deg2rad(latitude_in)
    elif angtype == "rad":
        latitude_rad = latitude_in
    elif angtype == "mas":
        log.warning("double check the mas conversion for security")
        latitude_rad = (
            2 * np.pi * (86400.0 * 10**3)
        ) / latitude_in  ### a verifier ....
    else:
        log.error("bad angtype %s", angtype)
        return None

    parallel_radius = np.cos(latitude_rad) * earth_radius

    return angle2equivalent_earth_radius(angle_in, angtype, parallel_radius)


def anglesfromvects(xa, ya, xb, yb, angtype="deg"):
    """
    Determines the angle between the points A(xa,ya) and B(xb,yb)

    angle can be : "deg", "rad"

    """
    a = np.array([xa, ya])
    b = np.array([xb, yb])

    ps = np.inner(a, b)
    a = np.arccos(ps / (np.linalg.norm(a) * np.linalg.norm(b)))

    if angtype == "deg":
        return np.rad2deg(a)
    elif angtype == "rad":
        return a
    else:
        log.error("bad angtype")
        return None


def angle_interpolation_quick(a, b, w):
    """
    Determine the interpolation between angle A & B
    by conversion to the cartesian space
    the parameter w € [0,1] define the interpoled angle C(w)
    where c(w=0) = a  &  c(w=1) = b

    References
    ----------
    https://stackoverflow.com/questions/2708476/rotation-interpolation

    """

    cs = (1 - w) * np.cos(a) + w * np.cos(b)
    sn = (1 - w) * np.sin(a) + w * np.sin(b)
    c = np.arctan2(sn, cs)

    return c


def angle_from_3_pts(p1, p2, p3):
    """
    Gives angle between 3 points
    p3 is the vertex (sommet)

    References
    ----------
    http://www.les-mathematiques.net/phorum/read.php?8,596072,596231

    """

    p1, p2, p3 = [np.array(e) for e in (p1, p2, p3)]
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3

    kos = ((x1 - x3) * (x2 - x3) + (y1 - y3) * (y2 - y3)) / (
        np.linalg.norm(p1 - p3) * np.linalg.norm(p2 - p3)
    )

    return np.arccos(kos)


def cartesian2polar(x, y):
    """
    Coordinates conversion

    cartesian => polar conversion

    Parameters
    ----------
    x , y : numpy.array of float
        cartesian coordinates

    Returns
    -------
    r , theta : float
        polar coordinates (radius / angle in radians)
    """
    if utils.is_iterable(x):
        x = np.array(x)
    if utils.is_iterable(y):
        y = np.array(y)

    theta = np.arctan2(y, x)
    r = np.sqrt(x ** +(y**2))
    return r, theta


def polar2cartesian(r, theta, ang="deg"):
    """
    Coordinates conversion

    polar => cartesian conversion

    Parameters
    ----------
    r , theta : float or iterable of floats
        polar coordinates

    ang : string
        'deg' (degrees) or 'rad' (radian)

    Returns
    -------
    x , y : numpy.array of float
        cartesian coordinates
    """
    if utils.is_iterable(r):
        r = np.array(r)
    if utils.is_iterable(theta):
        theta = np.array(theta)

    if ang == "deg":
        theta = np.deg2rad(theta)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y
