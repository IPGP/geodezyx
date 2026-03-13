#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 15:11:12 2020

@author: psakicki

This sub-module of geodezyx.conv deals with low-level geometric operation.

it can be imported directly with:
from geodezyx import conv

The geodezyx toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/IPGP/geodezyx
"""

########## BEGIN IMPORT ##########
#### External modules
import numpy as np

##########  END IMPORT  ##########


#  _                     _                _
# | |                   | |              | |
# | |     _____      __ | | _____   _____| |
# | |    / _ \ \ /\ / / | |/ _ \ \ / / _ \ |
# | |___| (_) \ v  v /  | |  __/\ v /  __/ |
# |______\___/ \_/\_/   |_|\___| \_/ \___|_|_         __                  _   _
#                                 | |      (_)       / _|                | | (_)
#   __ _  ___  ___  _ __ ___   ___| |_ _ __ _  ___  | |_ _   _ _ __   ___| |_ _  ___  _ __
#  / _` |/ _ \/ _ \| '_ ` _ \ / _ \ __| '__| |/ __| |  _| | | | '_ \ / __| __| |/ _ \| '_ \
# | (_| |  __/ (_) | | | | | |  __/ |_| |  | | (__  | | | |_| | | | | (__| |_| | (_) | | | |
#  \__, |\___|\___/|_| |_| |_|\___|\__|_|  |_|\___| |_|  \__,_|_| |_|\___|\__|_|\___/|_| |_|
#   __/ |
#  |___/

### Low level geometric function


def dist(a, b):
    """
    Cartesian distance between 2 points A & B

    Parameters
    ----------
    a,b : numpy.array of float
        3D-points

    Returns
    -------
    d : float
        distance between A & B
    """
    a = np.array(a)
    b = np.array(b)
    return np.linalg.norm(a - b)


def dist_diff(a, b):
    """
    First derivative of cartesian distance between 2 points A & B

    Parameters
    ----------
    a, b : numpy.array of float
        "points", can be 2D or 3D vectors (list, np.array ...)

    Returns
    -------
    diff_a : numpy.array of float
        [ dD/dxa , dD/dya , dD/dza ]

    diff_b : numpy.array of float
        [ dD/dxb , dD/dyb , dD/dzb ] = -diff_a
    """

    d_ab = a - b
    d = np.linalg.norm(d_ab)

    diff_a = d_ab / d
    diff_b = -d_ab / d

    return diff_a, diff_b


def relative_orientation(x1, y1, x2, y2, trigo_orient=True):
    """
    Return the angle between a point 1 and a point 2
    (reference vector : to the right)
    """
    if trigo_orient:
        ang1 = np.mod(360 - np.rad2deg(np.arctan2((x2 - x1), (y2 - y1))), 360)
    else:
        ang1 = np.mod(np.rad2deg(np.arctan2((x2 - x1), (y2 - y1))), 360)
    return ang1


def barycenter(points_list_in):
    """
    Determines the barycenter of a list of points
    """
    points_arr = np.array(points_list_in)
    return np.mean(points_arr[:, -3:], axis=0)


def pythagore(a, b, c=0):
    """
    Computes Pythagore's formula
    """
    return np.sqrt(a**2 + b**2 + c**2)


def equilateral_triangle(side):
    """
    Gives coordinates of an equilateral triangle of a given side
    """
    hauteur_len = np.sqrt(side**2 - (side * 0.5) ** 2)
    a = np.array([0, (2.0 / 3.0) * hauteur_len])
    b = np.array([side / 2.0, -(1.0 / 3.0) * hauteur_len])
    c = np.array([-side / 2.0, -(1.0 / 3.0) * hauteur_len])

    return a, b, c


def vincenty_full(point1, point2, miles=False, full=True, azimuth_in_deg=True):
    """

    Vincenty's formula (inverse method) to calculate the distance (in
    kilometers or miles) between two points on the surface of a spheroid.
    Gives also the Azimuth between 2 points

    Parameters
    ----------
    point1, point2 : iterable of float
        coordinates of the points

    miles : bool
        kilometers if True

    full : bool
        Description param3

    azimuth_in_deg : bool
        azimut in Rad if False

    Returns
    -------
    s : float
        Distance betwwen the 2 points

    fwd_az,rev_az : float
        Forward and Reverse Azimuth between the 2 points

    References
    ----------
    https://github.com/maurycyp/vincenty/blob/master/vincenty/

    Examples
    --------
    >>> vincenty_full((0.0, 0.0), (0.0, 0.0))  # coincident points
    0.0
    >>> vincenty_full((0.0, 0.0), (0.0, 1.0))
    111.319491
    >>> vincenty_full((0.0, 0.0), (1.0, 0.0))
    110.574389
    >>> vincenty_full((0.0, 0.0), (0.5, 179.5))  # slow convergence
    19936.288579
    >>> vincenty_full((0.0, 0.0), (0.5, 179.7))  # failure to converge
    >>> boston = (42.3541165, -71.0693514)
    >>> newyork = (40.7791472, -73.9680804)
    >>> vincenty_full(boston, newyork)
    298.396057
    >>> vincenty_full(boston, newyork, miles=True)
    185.414657
    """

    # WGS 84
    a = 6378137  # meters
    f = 1 / 298.257223563
    b = 6356752.314245  # meters; b = (1 - f)a

    miles_per_kilometer = 0.621371

    max_iterations = 200
    convergence_threshold = 1e-12  # .000,000,000,001

    # short-circuit coincident points
    if point1[0] == point2[0] and point1[1] == point2[1]:
        return 0.0

    import math

    u1 = math.atan((1 - f) * math.tan(math.radians(point1[0])))
    u2 = math.atan((1 - f) * math.tan(math.radians(point2[0])))
    llambda = math.radians(point2[1] - point1[1])

    sin_u1 = math.sin(u1)
    cos_u1 = math.cos(u1)
    sin_u2 = math.sin(u2)
    cos_u2 = math.cos(u2)
    sin_l = math.sin(llambda)
    cos_l = math.cos(llambda)

    for iteration in range(max_iterations):
        sin_lambda = math.sin(llambda)
        cos_lambda = math.cos(llambda)
        sin_sigma = math.sqrt(
            (cos_u2 * sin_lambda) ** 2
            + (cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda) ** 2
        )
        if sin_sigma == 0:
            return 0.0  # coincident points
        cos_sigma = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_lambda
        sigma = math.atan2(sin_sigma, cos_sigma)
        sin_alpha = cos_u1 * cos_u2 * sin_lambda / sin_sigma
        cos_sq_alpha = 1 - sin_alpha**2
        try:
            cos2_sigma_m = cos_sigma - 2 * sin_u1 * sin_u2 / cos_sq_alpha
        except ZeroDivisionError:
            cos2_sigma_m = 0
        c = f / 16 * cos_sq_alpha * (4 + f * (4 - 3 * cos_sq_alpha))
        lambda_prev = llambda
        llambda = llambda + (1 - c) * f * sin_alpha * (
            sigma
            + c
            * sin_sigma
            * (cos2_sigma_m + c * cos_sigma * (-1 + 2 * cos2_sigma_m**2))
        )
        if abs(llambda - lambda_prev) < convergence_threshold:
            break  # successful convergence
    else:
        return None  # failure to converge

    u_sq = cos_sq_alpha * (a**2 - b**2) / (b**2)
    A = 1 + u_sq / 16384 * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)))
    B = u_sq / 1024 * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)))
    delta_sigma = (
        B
        * sin_sigma
        * (
            cos2_sigma_m
            + B
            / 4
            * (
                cos_sigma * (-1 + 2 * cos2_sigma_m**2)
                - B
                / 6
                * cos2_sigma_m
                * (-3 + 4 * sin_sigma**2)
                * (-3 + 4 * cos2_sigma_m**2)
            )
        )
    )
    s = b * A * (sigma - delta_sigma)

    s /= 1000  # meters to kilometers
    if miles:
        s *= miles_per_kilometer  # kilometers to miles

    if not full:
        return s
    else:
        fwd_az = np.arctan2(cos_u2 * sin_l, cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_l)
        rev_az = np.arctan2(cos_u1 * sin_l, -sin_u1 * cos_u2 + cos_u1 * sin_u2 * cos_l)

        if azimuth_in_deg:
            fwd_az = np.rad2deg(fwd_az)
            rev_az = np.rad2deg(rev_az)

        return s, fwd_az, rev_az


def orthogonal_projection(Xa, Xb, Xv):
    """
    Project a point A on a line defined with a vector v and a point B

    Parameters
    ----------
    Xa : list/numpy.array of float
        Coordinates of A point, we want to project

    Xb : list/numpy.array of float
        Coordinates of B point, the origin point of the vector

    Xv : list/numpy.array of float
        Coordinates of the line direction vector

    Returns
    -------
    Xh : numpy.array of float
        Coordinates of H point, projection of A

    D : float
        Distance between A and H

    Notes
    -----
    Misc. Notes

    Source
    ------
    https://fr.wikipedia.org/wiki/Projection_orthogonale
    """

    xa, ya = Xa
    xb, yb = Xb
    xv, yv = Xv

    d = np.sqrt(xv**2 + yv**2)
    bh = ((xa - xb) * xv + (ya - yb) * yv) / d

    xh = xb + (bh / d) * xv
    yh = yb + (bh / d) * yv

    Xh = np.array([xh, yh])
    return Xh, dist(Xa, Xh)


def line_maker(x1, y1, x2, y2, nbpts=10000):
    """
    Determine points of a line easily

    Parameters
    ----------
    x1,y1 : float
        Coordinates of the start point

    x2,y2 : float
        Coordinates of the end point

    nbpts : int
        Number of points to determine on the line

    Returns
    -------
    xout, yout : numpy.array of float
        points of the line
    """

    xout = np.linspace(x1, x2, nbpts)
    yout = np.linspace(y1, y2, nbpts)
    return xout, yout
