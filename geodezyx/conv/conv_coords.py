#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 15:11:12 2020

@author: psakic

This sub-module of geodezyx.conv deals with coordinate conversion.

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

########## BEGIN IMPORT ##########
#### External modules
import numpy as np

#### geodeZYX modules
from geodezyx import utils
from geodezyx.conv import conv_rotation_matrices as rotmat

# import re
log = logging.getLogger("geodezyx")

##########  END IMPORT  ##########
#### Coordinates conversion

#   _____                    _ _             _               _____                              _
#  / ____|                  | (_)           | |             / ____|                            (_)
# | |     ___   ___  _ __ __| |_ _ __   __ _| |_ ___  ___  | |     ___  _ ____   _____ _ __ ___ _  ___  _ __
# | |    / _ \ / _ \| '__/ _` | | '_ \ / _` | __/ _ \/ __| | |    / _ \| '_ \ \ / / _ \ '__/ __| |/ _ \| '_ \
# | |___| (_) | (_) | | | (_| | | | | | (_| | ||  __/\__ \ | |___| (_) | | | \ v /  __/ |  \__ \ | (_) | | | |
#  \_____\___/ \___/|_|  \__,_|_|_| |_|\__,_|\__\___||___/  \_____\___/|_| |_|\_/ \___|_|  |___/_|\___/|_| |_|
#
#


def geo2xyz(lat, lon, h, angle="deg", e2=0.00669438003):
    """
    Coordinates conversion

    FLH Geographic => XYZ ECEF Cartesian

    Parameters
    ----------

    lat,lon,h : numpy.array of floats
        latitude (deg/rad), longitude (deg/rad), height (m)

    angle : string
        describe the angle input type : 'deg' or 'rad'

    e2 : floats
        ellipsoid parameters (WGS84 per default)


    Returns
    -------
    X,Y,Z : numpy.array of floats
        cartesian coordinates X,Y,Z in meters

    References
    ----------
    Based on PYACS of J.-M. Nocquet
    """

    if angle == "deg":
        lon = np.deg2rad(lon)
        lat = np.deg2rad(lat)

    # f=1.0 - sqrt(1-e2)

    wn = wnorm(lat)

    x = (wn + h) * np.cos(lat) * np.cos(lon)
    y = (wn + h) * np.cos(lat) * np.sin(lon)
    z = (wn * (1.0 - e2) + h) * np.sin(lat)

    return x, y, z


def xyz2geo(x, y, z, outdeg=True, a=6378137.0, e2=0.00669438003):
    """
    Coordinates conversion

    XYZ ECEF Cartesian => FLH Geographic

    Parameters
    ----------
    x,y,z : numpy.array of floats
        cartesian coordinates X,Y,Z in meters

    outdeg : bool
        if True, degrees output for the angle, else radian

    a, e2 : floats
        ellipsoid parameters (WGS84 per default)


    Returns
    -------
    latitude, longitude, height : numpy.array of floats
        latitude (deg/rad), longitude (deg/rad), height (m)

    References
    ----------
    Based on PYACS of J.-M. Nocquet
    """

    f = 1.0 - np.sqrt(1 - e2)

    tp = np.sqrt(x**2 + y**2)
    r = np.sqrt(tp**2 + z**2)

    tmu = np.arctan2((z / tp) * ((1.0 - f) + (e2 * a) / r), 1)
    rlbda = np.arctan2(y, x)

    s3 = (np.sin(tmu)) ** 3
    c3 = (np.cos(tmu)) ** 3
    t1 = z * (1 - f) + e2 * a * s3
    t2 = (1 - f) * (tp - e2 * a * c3)

    rphi = np.arctan2(t1, t2)

    rhe = tp * (np.cos(rphi)) + z * (np.sin(rphi))
    rhe = rhe - a * (np.sqrt(1 - e2 * (np.sin(rphi) ** 2)))

    if outdeg:
        rphi = np.rad2deg(rphi)
        rlbda = np.rad2deg(rlbda)

    return rphi, rlbda, rhe

def xyz2enu(x, y, z, x0, y0, z0):
    """
    Coordinates conversion

    XYZ ECEF Geocentric => ENU Topocentic

    Parameters
    ----------
    x,y,z : numpy.array of floats
        cartesian coordinates X,Y,Z in meters

    x0,y0,z0 : floats or numpy.array of floats
        coordinate of the wished topocentric origin point in the geocentric frame
        if they are iterable arrays (of the same size as X,Y,Z)
        a different x0,y0,z0 will be applied for each X,Y,Z element

    Returns
    -------
    e,n,u : numpy.array of floats
        East North Up Component (m) w.r.t. x0,y0,z0

    References
    ----------
    https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
    """

    f0, l0, h0 = xyz2geo(x0, y0, z0, outdeg=True)
    d_x = np.array(x) - x0
    d_y = np.array(y) - y0
    d_z = np.array(z) - z0
    e, n, u = xyz2enu_core(d_x, d_y, d_z, f0, l0)
    return e, n, u


def xyz2enu_around_fix_pos(x, y, z):
    """
    Computes the mean of the x,y,z and
    return the topocentric ENU position around its mean position

    Parameters
    ----------
    x,y,z : numpy.array of floats
        cartesian coordinates x,y,z in meters

    Returns
    -------
    E,N,U : numpy.array of floats
        East North Up Component (m) w.r.t. the mean position

    x0 , y0 , z0 : floats
        coordinates of the mean position in the geocentic frame
    """
    x0 = np.nanmean(x)
    y0 = np.nanmean(y)
    z0 = np.nanmean(z)
    return xyz2enu(x, y, z, x0, y0, z0), np.array([x0, y0, z0])


def enu2xyz(e, n, u, x0, y0, z0, velocity_mode=False):
    """
    Coordinates conversion

    ENU Topocentic => XYZ ECEF Geocentric

    Parameters
    ----------
    e,n,u : numpy.array of floats
        cartesian coordinates E,N,U in meters

    x0,y0,z0 : floats or numpy.array of floats
        coordinate of the topocentric origin point in the geocentric frame
        if they are iterable arrays (of the same size as E,N,U)
        a different x0,y0,z0 will be applied for each E,N,U element
    velocity_mode : bool
        For the Velocity mode, coordinates of the topocentric origin point
        are NOT added at the end of the conversion
        Default is False

    Returns
    -------
    x,y,z : numpy.array of floats
        ECEF X,Y,Z Component (m)

    References
    ----------
    https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
    """

    if utils.is_iterable(e):
        if not utils.is_iterable(x0):  ### The input x0,y0,z0 is a single point
            n = len(e)
            x0, y0, z0 = [x0] * n, [y0] * n, [z0] * n
        else:  ### The input x0,y0,z0 is a vector = different x0,y0,z0 for each point
            x0, y0, z0 = x0, y0, z0

        xlist, ylist, zlist = [], [], []

        for e, n, u, x0use, y0use, z0use in zip(e, n, u, x0, y0, z0):
            x, y, z = enu2xyz(e, n, u, x0use, y0use, z0use, velocity_mode=velocity_mode)
            xlist.append(x)
            ylist.append(y)
            zlist.append(z)
        return np.array(xlist), np.array(ylist), np.array(zlist)

    else:
        fr, lr, hr = xyz2geo(x0, y0, z0)
        f0 = np.deg2rad(fr)
        l0 = np.deg2rad(lr)

        r = np.array(
            [
                [-np.sin(l0), np.cos(l0), 0],
                [-np.sin(f0) * np.cos(l0), -np.sin(f0) * np.sin(l0), np.cos(f0)],
                [np.cos(f0) * np.cos(l0), np.cos(f0) * np.sin(l0), np.sin(f0)],
            ]
        )

        #r3 = r.T
        r3 = np.linalg.inv(r)

        enu = np.vstack((e, n, u))

        xyz = np.dot(r3, enu)

        if velocity_mode:
            x = float(xyz[0])
            y = float(xyz[1])
            z = float(xyz[2])
        else:
            x = float(xyz[0]) + x0
            y = float(xyz[1]) + y0
            z = float(xyz[2]) + z0

        return x, y, z



def xyz2azi_ele(x, y, z, x0, y0, z0, outdeg=False):
    """
    Convert ECEF Cartesian coordinates to azimuth and elevation angles.

    Parameters
    ----------
    x, y, z : numpy.array of floats
        ECEF X, Y, Z coordinates in meters.
    x0, y0, z0 : floats or numpy.array of floats
        Reference point coordinates in ECEF frame.

    Returns
    -------
    azi : numpy.array of floats
        Azimuth angle in radians.
    ele : numpy.array of floats
        Elevation angle in radians.
    d : numpy.array of floats
        Distance from the reference point to the point (x, y, z) in meters.
    outdeg : bool, optional
        If True, output angles in degrees. If False, in radians. Default is False.
    """
    e, n, u = xyz2enu(x, y, z, x0, y0, z0)
    d = np.sqrt(e**2 + n**2 + u**2)
    azi = np.arctan2(e, n)
    ele = np.arcsin(u/d)
    #ele2 = np.arctan2(u, np.sqrt(e**2 + n**2))
    if outdeg:
        azi = np.rad2deg(azi)
        ele = np.rad2deg(ele)
    return azi, ele, d






 # __      __       _              __          __
 # \ \    / /      | |             \ \        / /
 #  \ \  / /__  ___| |_ ___  _ __   \ \  /\  / / __ __ _ _ __  _ __   ___ _ __ ___
 #   \ \/ / _ \/ __| __/ _ \| '__|   \ \/  \/ / '__/ _` | '_ \| '_ \ / _ \ '__/ __|
 #    \  /  __/ (__| || (_) | |       \  /\  /| | | (_| | |_) | |_) |  __/ |  \__ \
 #     \/ \___|\___|\__\___/|_|        \/  \/ |_|  \__,_| .__/| .__/ \___|_|  |___/
 #                                                      | |   | |
 #                                                      |_|   |_|


def geo2xyz_vector(llh, angle="deg", a=6378137.0, e2=0.00669438003):
    """
    Convert an array of geographic coordinates (latitude, longitude, height) to ECEF Cartesian coordinates.

    Parameters
    ----------
    llh : array-like
        Array of shape (N, 3) or (3, N) containing latitude, longitude, and height.
    angle : str, optional
        Specifies if input angles are in degrees ('deg') or radians ('rad'). Default is 'deg'.
    a : float, optional
        Semi-major axis of the ellipsoid in meters. Default is 6378137.0 (WGS84).
    e2 : float, optional
        Square of the first eccentricity of the ellipsoid. Default is 0.00669438003 (WGS84).

    Returns
    -------
    xyz : numpy.ndarray
        Array of shape (N, 3) with ECEF X, Y, Z coordinates.
    """
    #### if Nx3 array => 3xN array
    llh = utils.transpose_vector_array(llh)

    x, y, z = geo2xyz(llh[0], llh[1], llh[2], angle=angle, a=a, e2=e2)
    xyz = np.column_stack((x, y, z))

    return xyz


def xyz2enu_vector(xyz, xyz0):
    """
    Convert an array of ECEF Cartesian coordinates to ENU topocentric coordinates.

    Parameters
    ----------
    xyz : array-like
        Array of shape (N, 3) or (3, N) with ECEF X, Y, Z coordinates.
    xyz0 : array-like
        Reference point as array-like of length 3 (X0, Y0, Z0).

    Returns
    -------
    enu : numpy.ndarray
        Array of shape (N, 3) with East, North, Up coordinates.
    """
    xyz = utils.transpose_vector_array(xyz)

    e, n, u = xyz2enu(xyz[0], xyz[1], xyz[2], xyz0[0], xyz0[1], xyz0[2])
    enu = np.column_stack((e, n, u))

    return enu


def xyz2geo_vector(xyz, outdeg=True, a=6378137.0, e2=0.00669438003):
    """
    Convert an array of ECEF Cartesian coordinates to geographic coordinates (latitude, longitude, height).

    Parameters
    ----------
    xyz : array-like
        Array of shape (N, 3) or (3, N) with ECEF X, Y, Z coordinates.
    outdeg : bool, optional
        If True, output angles are in degrees. If False, in radians. Default is True.
    a : float, optional
        Semi-major axis of the ellipsoid in meters. Default is 6378137.0 (WGS84).
    e2 : float, optional
        Square of the first eccentricity of the ellipsoid. Default is 0.00669438003 (WGS84).

    Returns
    -------
    flh : numpy.ndarray
        Array of shape (N, 3) with latitude, longitude, and height.
    """
    xyz = utils.transpose_vector_array(xyz)

    f, l, h = xyz2geo(xyz[0], xyz[1], xyz[2], outdeg=outdeg, a=a, e2=e2)

    flh = np.column_stack((f, l, h))

    return flh


def enu2xyz_vector(enu, xyz_ref):
    """
    Convert an array of ENU topocentric coordinates to ECEF Cartesian coordinates.

    Parameters
    ----------
    enu : array-like
        Array of shape (N, 3) or (3, N) with East, North, Up coordinates.
    xyz_ref : array-like
        Reference point as array-like of length 3 (X0, Y0, Z0).

    Returns
    -------
    xyz : numpy.ndarray
        Array of shape (N, 3) with ECEF X, Y, Z coordinates.
    """
    enu = utils.transpose_vector_array(enu)

    x, y, z = enu2xyz(enu[0], enu[1], enu[2], xyz_ref[0], xyz_ref[1], xyz_ref[2])

    xyz = np.column_stack((x, y, z))

    return xyz

#      _                             __      _      _      _                                                _
#     (_)                           / /     | |    | |    | |                                              (_)
#  ___ _  __ _ _ __ ___   __ _     / /   ___| |_ __| |  __| | _____   __   ___ ___  _ ____   _____ _ __ ___ _  ___  _ __
# / __| |/ _` | '_ ` _ \ / _` |   / /   / __| __/ _` | / _` |/ _ \ \ / /  / __/ _ \| '_ \ \ / / _ \ '__/ __| |/ _ \| '_ \
# \__ \ | (_| | | | | | | (_| |  / /    \__ \ || (_| || (_| |  __/\ v /  | (_| (_) | | | \ v /  __/ |  \__ \ | (_) | | | |
# |___/_|\__, |_| |_| |_|\__,_| /_/     |___/\__\__,_(_)__,_|\___| \_(_)  \___\___/|_| |_|\_/ \___|_|  |___/_|\___/|_| |_|
#         __/ |
#        |___/

def sigma_xyz2enu(x, y, z, s_x, s_y, s_z, s_xy=0, s_yz=0, s_xz=0):
    """
    Convert standard deviation
    Cartesian ECEF XYZ => Cartesian Topocentric ENU

    Note
    ----
    Inputs values are now assumed correlated (241105)

    References
    ----------
    Linear Algebra, Geodesy, and GPS p332
    https://stackoverflow.com/questions/51162460/converting-ecef-xyz-covariance-matrix-to-enu-covariance-matrix
    https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
    """

    sigma_xyz = np.array(
        [[s_x**2, s_xy, s_xz], [s_xy, s_y**2, s_yz], [s_xz, s_yz, s_z**2]]
    )

    f, l, h = xyz2geo(x, y, z)

    # old and bad rotation matrix (bofore 20241104)
    # C = rotmat.C_ecef2enu(F,L,angtype='deg')
    # new and good rotation matrix (after 20241104)
    c = rotmat.c_ecef2enu_sigma(f, l, angtype="deg")

    sigma_enu = np.dot(np.dot(c, sigma_xyz), c.T)
    sigma_enu2 = np.dot(c, np.dot(sigma_xyz, c.T))

    # lon = np.deg2rad(L)
    # lat = np.deg2rad(F)

    s_e = np.sqrt(sigma_enu[0, 0])
    s_n = np.sqrt(sigma_enu[1, 1])
    s_u = np.sqrt(sigma_enu[2, 2])

    return s_e, s_n, s_u


def sigma_geo2xyz(lat, lon, h, s_f, s_l, s_h, ang="deg"):
    """
    Convert standard deviation
    Geographic FLH => Cartesian ECEF XYZ

    Warning
    -------
    Inputs values are assumed as uncorrelated, which is not accurate
    Must be improved

    References
    ----------
    Linear Algebra, Geodesy, and GPS p332
    """

    log.warning("Inputs values are assumed as uncorrelated, which is not accurate")
    log.warning("Prefer sigma_xyz2enu")

    if ang == "deg":
        lat = np.deg2rad(lat)
        lon = np.deg2rad(lon)
        s_f = np.deg2rad(s_f)
        s_l = np.deg2rad(s_l)
    x, y, z = geo2xyz(lat, lon, h)
    x2, y2, z2 = geo2xyz(lat + s_f, lon + s_l, h + s_h)

    return np.abs(x - x2), np.abs(y - y2), np.abs(z - z2)


def sigma_geo2enu(lat, lon, h, s_f, s_l, s_h, ang="deg"):
    """
    Convert standard deviation
    Geographic FLH => Cartesian Topocentric ENU

    Warning
    -------
    Inputs values are assumed as uncorrelated, which is not accurate
    Must be improved

    References
    ----------
    Linear Algebra, Geodesy, and GPS p332
    """
    # conversion batarde du sigma FLH => sigma ENU
    # Par conversion des angles en distance

    log.warning("Inputs values are assumed as uncorrelated, which is not accurate")
    log.warning("Prefer sigma_xyz2enu")

    if ang == "deg":
        lat = np.deg2rad(lat)
        lon = np.deg2rad(lon)
        s_f = np.deg2rad(s_f)
        s_l = np.deg2rad(s_l)

    e2 = 0.00669438003
    # f=1.0 - np.sqrt(1-e2)
    a = 6378137.0

    # Je présume que Rpsi est le rayon du cercle tangent à l'ellipsoide à une
    # lattitude donnée (comm du 150710)
    psi = np.arctan((1 - e2) * np.tan(lat))
    r_psi = a * np.sqrt(1 - e2) / np.sqrt(1 - e2 * np.cos(psi) ** 2)

    # r est le rayon d'un petit cercle (un parallèle)
    r = np.cos(psi) * (r_psi + h)

    # On pourrait simplifier par 2pi mais autant avoir toute la démarche
    s_e = (np.pi * 2 * r * s_l) / (2 * np.pi)
    s_n = (np.pi * 2 * r_psi * s_f) / (2 * np.pi)
    s_u = s_h
    return s_e, s_n, s_u


def sigma_enu2geo(lat, lon, h, s_e, s_n, s_u, ang="deg", a=6378137.0, e2=0.00669438003):
    """
    Convert standard deviation
    Cartesian Topocentric ENU => Geographic FLH

    Warning
    -------
    Inputs values are assumed as uncorrelated, which is not accurate
    Must be improved

    References
    ----------
    Linear Algebra, Geodesy, and GPS p332
    """

    log.warning("Inputs values are assumed as uncorrelated, which is not accurate")
    log.warning("Prefer sigma_xyz2enu")

    # conversion batarde du sigma ENU => sigma FLH
    # Par conversion des angles en distance

    # f=1.0 - np.sqrt(1-e2)

    # Je présume que Rpsi est le rayon du cercle tangent à l'ellipsoide à une
    # lattitude donnée (comm du 150710)
    psi = np.arctan((1 - e2) * np.tan(lat))
    r_psi = a * np.sqrt(1 - e2) / np.sqrt(1 - e2 * np.cos(psi) ** 2)

    # r est le rayon d'un petit cercle (un parallèle)
    r = np.cos(psi) * (r_psi + h)

    # On pourait simplifier par 2pi mais autant avoir toute la démarche
    s_l = s_e * (2 * np.pi) / (np.pi * 2 * r)
    s_f = s_n * (2 * np.pi) / (np.pi * 2 * r_psi)
    s_h = s_u

    if ang == "deg":
        s_l = np.rad2deg(s_l)
        s_f = np.rad2deg(s_f)
    return s_f, s_l, s_h


def eci2rtn_or_rpy(p, v, c, out_rpy=False, rpy_theo_mode=False):
    """
    convert ECI coordinates in RTN (RIC) or RPY (Roll Pitch Yaw)

    Parameters
    ----------
    p : numpy.array
        3D vector, position of the ref object in ECI frame
    v : numpy.array
        3D vector, velocity of the ref object in ECI frame
    c : numpy.array
        3D Vector, coordinates in ECF frame that will be transformed
    out_rpy : bool
        if True output in RPY frame, RTN instead
    rpy_theo_mode : bool
        use the theoretical matrix composition, but wrong
        ans. empirically, only for debug !!

    Returns
    -------
    c_out : conversion of p in RTN or RPY ref. frame

    References
    ----------
    "Coordinate Systems", ASEN 3200 1/24/06 George H. Born
    """

    c_eci2rtn_mat = rotmat.c_eci2rtn(p, v)

    if not out_rpy:
        trans_mat = c_eci2rtn_mat
    else:

        if not rpy_theo_mode:
            # Pour de très obscures raisons la composition est inversée
            # par rapport à l'ordre standard ... (241017)
            trans_mat = np.dot(rotmat.c_rtn2rpy().T, c_eci2rtn_mat.T)
        else:
            log.warning(
                "using the theoretical mode for RPY conversion, UNSTABLE & WRONG !"
            )
            trans_mat = np.dot(rotmat.c_rtn2rpy(), c_eci2rtn_mat)

        # Mais reste compatible avec Wikipedia
        # https://en.wikipedia.org/wiki/Permutation_matrix#Permutation_of_rows_and_columns

    c_out = trans_mat.dot(c)

    # EXEMPLE POUR DEBUG
    #     AAAAAAAAAAAAA yawerr,aux,nomi
    #       -1.606542300481316
    #       -1.606542300481316
    # CCCCCCCCCCCC sunors orientation du soleil
    #   2.379858191647978E-002   6.654859888132475E-001  -7.460308480163662E-001
    # DDDDDDDDDDDDD Inputs : xsat , vsat , xsun
    #   -15779.437215560150435   -23508.964967410494864     8644.453520237560952
    #        1.296220384065493       -1.917314252639830       -2.847293074880089
    #  -1.342404641197934E+008  -5.934080206107079E+007  -2.572342068927875E+007
    # EEEEEEEEEEEE Rcrs2ors
    #   3.532662370257567E-001  -7.688096793507927E-001   5.330195492501879E-001
    #  -5.225364448455557E-001   3.104619621186746E-001   7.941181766553120E-001
    #  -7.759888076421138E-001  -5.590572803191509E-001  -2.920042493231377E-001
    #    Rcrs2ors = np.array([[ 0.35330134, -0.76880968,  0.53301955],
    #                         [-0.52248415,  0.31046196,  0.79411818],
    #                         [-0.77600804, -0.55905728, -0.29200425]])

    #    Psat = np.array([   -15779.437215560150435 ,  -23508.964967410494864  ,   8644.453520237560952])
    #    Vsat = np.array([        1.296220384065493  ,     -1.917314252639830  ,     -2.847293074880089])
    #    Psun = np.array([  -1.342404641197934E+008 , -5.934080206107079E+007 , -2.572342068927875E+007])

    return c_out


def eci2rtn(p,v,c):
    """
    legacy wrapper of eci2rtn_or_rpy
    """
    return eci2rtn_or_rpy(p,v,c, out_rpy=False)


def ecef2eci(xyz, utc_times):
    """
    Convert ECEF (Earth Centered Earth Fixed) positions to ECI (Earth Centered Inertial)
    positions

    Parameters
    ----------
    xyz : numpy.array of floats
        XYZ are cartesian positions in ECEF. Should have shape (N,3)

    utc_times : numpy.array of floats
        UTC_times are UTC timestamps, as datetime objects. Sould have shape (N)

    Returns
    -------
    eci : numpy.array of floats
        Earth Centered Inertial coordinates. will have shape (N,3)

    Note
    ----
    Requires pyorbital module

     [X]    [C -S 0][X]
     [Y]  = [S  C 0][Y]
     [Z]eci [0  0 1][Z]ecf

     C and S are cos() and sin() of gmst (Greenwich Meridian Sideral Time)

    References
    ----------
    http://ccar.colorado.edu/ASEN5070/handouts/coordsys.doc
    Inspired from satellite-js (https://github.com/shashwatak/satellite-js)
    """

    # XYZ and utc_time must have the same shape
    # if not xyz.shape[:-1] == utc_times.shape:
    #    raise ValueError("shape mismatch for XYZ and utc_times (got {} and {})".format(xyz.shape[:-1],utc_times.shape))

    #    gmst = -1 * astronomy.gmst(utc_times) # EDIT 180430 : Why this -1 ??? removed because wrong ! ...
    from pyorbital import astronomy

    gmst = 1 * astronomy.gmst(utc_times)

    eci = xyz.copy()
    eci[:, 0] = xyz[:, 0] * np.cos(gmst) - xyz[:, 1] * np.sin(gmst)
    eci[:, 1] = xyz[:, 0] * np.sin(gmst) + xyz[:, 1] * np.cos(gmst)
    return eci


def eci2ecef(xyz, utc_times):
    """
    Convert ECI (Earth Centered Inertial) positions to ECEF (Earth Centered Earth Fixed)
    positions

    Parameters
    ----------
    xyz : numpy.array of floats
        XYZ are cartesian positions in Earth Centered Inertial. Should have shape (N,3)

    utc_times : numpy.array of floats
        UTC_times are UTC timestamps, as datetime objects. Sould have shape (N)

    Returns
    -------
    ecef : numpy.array of floats
        Earth Centered Earth Fixed coordinates. will have shape (N,3)

    Note
    ----
    Requires pyorbital module


     [X]          ([C -S 0])[X]
     [Y]     = inv([S  C 0])[Y]
     [Z]ecef      ([0  0 1])[Z]eci


    Empirically:
     [X]       [ C  S  0][X]
     [Y]     = [-S  C  0][Y]
     [Z]ecef   [ 0  0  1][Z]eci



     C and S are cos() and sin() of gmst (Greenwich Meridian Sideral Time)

    References
    ----------
    http://ccar.colorado.edu/ASEN5070/handouts/coordsys.doc
    Inspired from satellite-js (https://github.com/shashwatak/satellite-js)


    Note
    ----
    Quick mode of the reverse fct, can be improved
    """
    # XYZ and utc_time must have the same shape
    # if not xyz.shape[:-1] == utc_times.shape:
    #    raise ValueError("shape mismatch for XYZ and utc_times (got {} and {})".format(xyz.shape[:-1],utc_times.shape))

    #    gmst = -1 * astronomy.gmst(utc_times)
    # EDIT 180430 : Why this -1 ??? removed because wrong ! ...
    from pyorbital import astronomy

    gmst = 1 * astronomy.gmst(utc_times)

    ecef = xyz.copy()
    ecef[:, 0] = +xyz[:, 0] * np.cos(gmst) + xyz[:, 1] * np.sin(gmst)
    ecef[:, 1] = -xyz[:, 0] * np.sin(gmst) + xyz[:, 1] * np.cos(gmst)

    return ecef


#### Core & Legacy


 #   _____                          _
 #  / ____|                 ___    | |
 # | |     ___  _ __ ___   ( _ )   | |     ___  __ _  __ _  ___ _   _
 # | |    / _ \| '__/ _ \  / _ \/\ | |    / _ \/ _` |/ _` |/ __| | | |
 # | |___| (_) | | |  __/ | (_>  < | |___|  __/ (_| | (_| | (__| |_| |
 #  \_____\___/|_|  \___|  \___/\/ |______\___|\__, |\__,_|\___|\__, |
 #                                              __/ |            __/ |
 #                                             |___/            |___/

def vector_separator(abc):
    """
    Split a Nx3 Array/DataFrame in three separated 1-component vectors
    To simplify the usage of the conversion functions
    (which take single component vectors as input)

    Parameters
    ----------
    abc : Array or DataFrame
        Nx3 XYZ, ENU... array.

    Returns
    -------
    a : Array
        1st component.
    b : Array
        2nd component.
    c : Array
        3rd component.

    """
    abc = np.array(abc)
    return abc[:, 0], abc[:, 1], abc[:, 2]


def wnorm(phi, a=6378137.0, e2=0.00669438003):
    """
    Compute the Ellipsoid "Grande Normale"

    References
    ----------
    ALG0021 in NTG_71 (IGN Lambert Projection Tech document)
    Based on PYACS of J.-M. Nocquet
    """
    from numpy import sqrt, sin

    # e=sqrt(e2)
    wd = sqrt(1 - e2 * sin(phi) ** 2)
    result = a / wd
    return result


def normal_vector(phi, llambda, angle="deg", normalized=True):
    """
    Compute the Ellipsoid "Normale"

    References
    ----------
    p. Bosser (2012), Géométrie de l'Ellipsoïde, p27
    """

    if angle == "deg":
        llambda = np.deg2rad(llambda)
        phi = np.deg2rad(phi)

    a = np.cos(llambda) * np.cos(phi)
    b = np.sin(llambda) * np.cos(phi)
    c = np.sin(phi)

    n = np.array([a, b, c])

    if normalized:
        n = n / np.linalg.norm(n)

    return n

def xyz2enu_core(d_x, d_y, d_z, lat0, lon0):
    """
    Coordinates conversion

    XYZ ECEF Geocentric => ENU Topocentic

    **use xyz2enu in priority**

    (xyz2enu_core is a core function for xyz2enu)

    dXYZ = XYZrover - XYZref

    Parameters
    ----------
    d_x,d_y,d_z : floats or numpy.array of floats
        cartesian coordinates difference between the considered point(s)
        and the reference point

    lat0,lon0 : floats or numpy.array of floats
        if they are iterable arrays (of the same size as d_x,dY,dZ)
        a different x0,y0,z0 will be applied for each d_x,dY,dZ element

    Returns
    -------
    E,N,U : numpy.array of floats
        East North Up Component (m) w.r.t. x0,y0,z0

    References
    ----------
    https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates

    Note
    ----
    This recursive fuction should be improved
    """

    ## Case one ref point per dXYZ
    if utils.is_iterable(lat0):
        e, n, u = [], [], []
        for dx_m, dy_m, dz_m, lat0_m, lon0_m in zip(d_x, d_y, d_z, lat0, lon0):
            e_m, n_m, u_m = xyz2enu_core(dx_m, dy_m, dz_m, lat0_m, lon0_m)
            e.append(e_m)
            n.append(n_m)
            u.append(u_m)

        return np.squeeze(np.array(e)), np.squeeze(np.array(n)), np.squeeze(np.array(u))

    # case onle ref point for all dXYZ
    else:
        f0 = np.deg2rad(lat0)
        l0 = np.deg2rad(lon0)

        r = np.array(
            [
                [-np.sin(l0), np.cos(l0), 0],
                [-np.sin(f0) * np.cos(l0), -np.sin(f0) * np.sin(l0), np.cos(f0)],
                [np.cos(f0) * np.cos(l0), np.cos(f0) * np.sin(l0), np.sin(f0)],
            ]
        )

        enu = np.dot(r, np.vstack((d_x, d_y, d_z)))

        e = enu[0, :]
        n = enu[1, :]
        u = enu[2, :]

        return e, n, u

def enu2xyz_legacy(e, n, u, x0, y0, z0):
    """
    KEPT FOR LEGACY REASONS, use enu2xyz

    diffère de enu2xyz pour d'obscure raisons, à investiguer !!!
    est laissé pour des scripts de conversion de GINS (170119)

    this fct compute the dXYZ and not the final XYZ

    """

    fr, lr, hr = xyz2geo(x0, y0, z0)
    f0 = np.deg2rad(fr)
    l0 = np.deg2rad(lr)

    r = np.array(
        [
            [-np.sin(l0), np.cos(l0), 0],
            [-np.sin(f0) * np.cos(l0), -np.sin(f0) * np.sin(l0), np.cos(f0)],
            [np.cos(f0) * np.cos(l0), np.cos(f0) * np.sin(l0), np.sin(f0)],
        ]
    )

    r3 = r.T

    enu = np.vstack((e, n, u))

    xyz = np.dot(r3, enu)  # + np.vstack(( x0, y0, z0))

    d_x = float(xyz[0])
    d_y = float(xyz[1])
    d_z = float(xyz[2])

    return d_x, d_y, d_z

############ deprecated aliases ############
import warnings

def deprec_warn(old, new):
    """
    Helper function to issue deprecation warnings
    """
    warnings.warn(
        f"{old} is deprecated and will be removed in future versions. Use {new} instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return None


def ECI2ECEF(*args, **kwargs):
    """
    Deprecated alias for eci2ecef
    """
    deprec_warn("ECI2ECEF", "eci2ecef")
    return eci2ecef(*args, **kwargs)


def ECEF2ECI(*args, **kwargs):
    """
    Deprecated alias for ecef2eci
    """
    deprec_warn("ECEF2ECI", "ecef2eci")
    return ecef2eci(*args, **kwargs)


def XYZ2ENU_2(*args, **kwargs):
    """
    Deprecated alias for xyz2enu
    """
    deprec_warn("XYZ2ENU_2", "xyz2enu")
    return xyz2enu(*args, **kwargs)


def ENU2XYZ(*args, **kwargs):
    """
    Deprecated alias for enu2xyz
    """
    deprec_warn("ENU2XYZ", "enu2xyz")
    return enu2xyz(*args, **kwargs)


def GEO2XYZ(*args, **kwargs):
    """
    Deprecated alias for geo2xyz
    """
    deprec_warn("GEO2XYZ", "geo2xyz")
    return geo2xyz(*args, **kwargs)


def XYZ2GEO(*args, **kwargs):
    """
    Deprecated alias for xyz2geo
    """
    deprec_warn("XYZ2GEO", "xyz2geo")
    return xyz2geo(*args, **kwargs)


def sXYZ2sENU(*args, **kwargs):
    deprec_warn("sXYZ2sENU", "sigma_xyz2enu")
    """
    Deprecated alias for sigma_xyz2enu
    """
    return sigma_xyz2enu(*args, **kwargs)


def sENU2sFLH(*args, **kwargs):
    """
    Deprecated alias for sigma_enu2geo
    """
    deprec_warn("sENU2sFLH", "sigma_enu2geo")
    return sigma_enu2geo(*args, **kwargs)

def sFLH2sENU(*args, **kwargs):
    """
    Deprecated alias for sigma_geo2enu
    """
    deprec_warn("sFLH2sENU", "sigma_geo2enu")
    return sigma_geo2enu(*args, **kwargs)

def sFLH2sXYZ(*args, **kwargs):
    """
    Deprecated alias for sigma_geo2xyz
    """
    deprec_warn("sFLH2sXYZ", "sigma_geo2xyz")
    return sigma_geo2xyz(*args, **kwargs)