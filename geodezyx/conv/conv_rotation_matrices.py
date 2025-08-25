#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 15:11:12 2020

@author: psakicki

This sub-module of geodezyx.conv deals with rotation, implemented with matrices.

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

log = logging.getLogger('geodezyx')


##########  END IMPORT  ##########

#  _____       _        _   _               __  __       _        _               
# |  __ \     | |      | | (_)             |  \/  |     | |      (_)              
# | |__) |___ | |_ __ _| |_ _  ___  _ __   | \  / | __ _| |_ _ __ _  ___ ___  ___ 
# |  _  // _ \| __/ _` | __| |/ _ \| '_ \  | |\/| |/ _` | __| '__| |/ __/ _ \/ __|
# | | \ \ (_) | || (_| | |_| | (_) | | | | | |  | | (_| | |_| |  | | (_|  __/\__ \
# |_|  \_\___/ \__\__,_|\__|_|\___/|_| |_| |_|  |_|\__,_|\__|_|  |_|\___\___||___/
#                                                                                 

### Rotation matrix

### The rotation matrices conventions are described in details in 
### Sakic (2016, PhD manuscript) p126

### It is based on the book of Grewal, M. S., Weill, L. R., & Andrews, A. p. (2007).
### Global Positioning Systems, Inertial Navigation, and Integration. 
### Hoboken, NJ, USA : John Wiley and Sons, Inc.

def c_2d(theta, angtype='deg'):
    if angtype == 'deg':
        theta = np.deg2rad(theta)

    return np.array([[np.cos(theta), -np.sin(theta)],
                     [np.sin(theta), np.cos(theta)]])


def c_ned2enu():
    return np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])


def c_enu2ned():
    return np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])


def c_ned2ecef(phi, llambda, angtype='deg'):
    if angtype == 'deg':
        phi = np.deg2rad(phi)
        llambda = np.deg2rad(llambda)

    in_ned = np.array([1, 0, 0])
    ie_ned = np.array([0, 1, 0])
    id_ned = np.array([0, 0, 1])

    ix_ned = np.array([-np.cos(llambda) * np.sin(phi), -np.sin(llambda), -np.cos(llambda) * np.cos(phi)])
    iy_ned = np.array([-np.sin(llambda) * np.sin(phi), np.cos(llambda), -np.sin(llambda) * np.cos(phi)])
    iz_ned = np.array([np.cos(phi), 0, -np.sin(phi)])

    c_ned2ecef_out = np.array([[np.dot(ix_ned, in_ned), np.dot(ix_ned, ie_ned), np.dot(ix_ned, id_ned)],
                               [np.dot(iy_ned, in_ned), np.dot(iy_ned, ie_ned), np.dot(iy_ned, id_ned)],
                               [np.dot(iz_ned, in_ned), np.dot(iz_ned, ie_ned), np.dot(iz_ned, id_ned)]])

    return c_ned2ecef_out


def c_ecef2ned(phi, llambda, angtype='deg'):
    return np.linalg.inv(c_ned2ecef(phi, llambda, angtype))


def c_enu2ecef(phi, llambda, angtype='deg'):
    if angtype == 'deg':
        phi = np.deg2rad(phi)
        llambda = np.deg2rad(llambda)

    ie_enu = np.array([1, 0, 0])
    in_enu = np.array([0, 1, 0])
    iu_enu = np.array([0, 0, 1])

    ix_enu = np.array([-np.sin(llambda), -np.cos(llambda) * np.sin(phi), np.cos(llambda) * np.cos(phi)])
    iy_enu = np.array([np.cos(llambda), -np.sin(llambda) * np.sin(phi), np.sin(llambda) * np.cos(phi)])
    iz_enu = np.array([0, np.cos(phi), np.sin(phi)])

    c_enu2ecef_out = np.array([[np.dot(ix_enu, ie_enu), np.dot(ix_enu, in_enu), np.dot(ix_enu, iu_enu)],
                               [np.dot(iy_enu, ie_enu), np.dot(iy_enu, in_enu), np.dot(iy_enu, iu_enu)],
                               [np.dot(iz_enu, ie_enu), np.dot(iz_enu, in_enu), np.dot(iz_enu, iu_enu)]])

    return c_enu2ecef_out


def c_ecef2enu(phi, llambda, angtype='deg'):
    return c_enu2ecef(phi, llambda, angtype).T


def c_ecef2enu_sigma(phi, llambda, angtype='deg'):
    """
    Gives the transformation matrix between ECEF and ENU (East North Up) 
    ref. frame for the sigmas (variance propagation)


    Source
    ------
    https://stackoverflow.com/questions/51162460/converting-ecef-xyz-covariance-matrix-to-enu-covariance-matrix
    """
    if angtype == 'deg':
        phi = np.deg2rad(phi)
        llambda = np.deg2rad(llambda)


    c_ecef2enu_sigma_out = np.array([[- np.sin(llambda), np.cos(llambda), 0],
                                     [-np.sin(phi) * np.cos(llambda), -np.sin(phi) * np.sin(llambda), np.cos(phi)],
                                     [np.cos(phi) * np.cos(llambda), np.cos(phi) * np.sin(llambda), np.sin(phi)]])

    return c_ecef2enu_sigma_out


def c_rpy2enu(roll, pitch, yaw, angtype='deg'):
    if angtype == 'deg':
        roll = np.deg2rad(roll)
        pitch = np.deg2rad(pitch)
        yaw = np.deg2rad(yaw)

    r = roll
    p = pitch
    y = yaw

    ir = np.array([np.sin(y) * np.cos(p), np.cos(y) * np.cos(p), np.sin(p)])
    ip = np.array([np.cos(r) * np.cos(y) + np.sin(r) * np.sin(y) * np.sin(p),
                   -np.cos(r) * np.sin(y) + np.sin(r) * np.cos(y) * np.sin(p),
                   -np.sin(r) * np.cos(p)])
    iy = np.array([-np.sin(r) * np.cos(y) + np.cos(r) * np.sin(y) * np.sin(p),
                   np.sin(r) * np.sin(y) + np.cos(r) * np.cos(y) * np.sin(p),
                   - np.cos(r) * np.cos(p)])

    #    Ie = np.array([ np.sin(Y) * np.cos(p) ,
    #    np.cos(R) * np.cos(Y) + np.sin(R) * np.sin(Y) * np.sin(p) ,
    #    -np.sin(R) * np.cos(Y) + np.cos(R) * np.sin(Y) * np.sin(p) ])
    #
    #    In = np.array([ np.cos(Y) * np.cos(p),
    #     -np.cos(R) * np.sin(Y) + np.sin(R) * np.cos(Y) * np.sin(p),
    #     np.sin(R) * np.sin(Y) + np.cos(R) * np.cos(Y) * np.sin(p) ])
    #
    #    Iu = np.array([ np.sin(p),
    #     -np.sin(R) * np.cos(p),
    #     -np.cos(R) * np.cos(p) ])

    c_rpy2enu_out = np.hstack((ir[np.newaxis].T, ip[np.newaxis].T, iy[np.newaxis].T))

    return c_rpy2enu_out


def c_rpy2enu2(roll, pitch, yaw, angtype='deg'):
    if angtype == 'deg':
        roll = np.deg2rad(roll)
        pitch = np.deg2rad(pitch)
        yaw = np.deg2rad(yaw)

    r = roll
    p = pitch
    y = yaw

    sr = np.sin(r)
    cr = np.cos(r)
    sp = np.sin(p)
    cp = np.cos(p)
    sy = np.sin(y)
    cy = np.cos(y)

    c_rpy2enu_out = np.array([[sy * cp, cr * cy + sr * sy * sp, - sr * cy + cr * sy * sp],
                              [cy * cp, -cr * sy + sr * cy * sp, sr * sy + cr * cy * sp],
                              [sp, -sr * cp, -cr * cp]])

    return c_rpy2enu_out


def c_eci2rtn(p, v):
    """
    gives the transformation matrix between ECI ref. frame and RTN (aka RIC or RSW)

    Parameters
    ----------
    p : numpy.array
        3D vector, position of the object in ECI frame

    v : numpy.array
        3D vector, velocity of the object in ECI frame
        
    Returns
    -------
    C_eci2rtn : numpy.array
        transformation matrix
            
    References
    ----------
    "Coordinate Systems", ASEN 3200 1/24/06 George H. Born
    https://www.colorado.edu/ASEN/asen3200/handouts/Coordinate%20System.pdf
    """

    h = np.cross(p, v)

    r = p / np.linalg.norm(p)
    c = h / np.linalg.norm(h)
    i = np.cross(c, r)

    c_eci2rtn_out = np.column_stack((r, i, c))

    return c_eci2rtn_out


def c_rtn2rpy():
    c_rtn2rpy_out = np.array([[0., 0, -1.],
                              [1., 0., 0.],
                              [0, -1., 0.]])
    return c_rtn2rpy_out


def c_cep2itrs(xpole, ypole):
    """
    xpole and ypole are given in mas, miliarcsecond 
    
    References
    ----------
        Hofmann-Wellenhof,et al. 2008 - GNSS - p21
        Xu - 2007 - GPS - p17
        Xu - Orbits - p15
    """

    xpole2 = np.deg2rad(xpole * (1. / 3600.) * (10 ** -3))
    ypole2 = np.deg2rad(ypole * (1. / 3600.) * (10 ** -3))

    c_cep2itrs_out = np.array([[1., 0., xpole2],
                               [0., 1., -ypole2],
                               [-xpole2, ypole2, 1.]])

    return c_cep2itrs_out


def c_euler(phi, theta, psi):
    """
    Gives the matrix of an Euler rotation
        
    References
    ----------
        https://fr.wikipedia.org/wiki/Angles_d%27Euler
    """
    cphi = np.cos(phi)
    ctheta = np.cos(theta)
    cpsi = np.cos(psi)
    sphi = np.sin(phi)
    stheta = np.sin(theta)
    spsi = np.sin(psi)

    c_euler_out = np.array([[cpsi * cphi - spsi * ctheta * sphi, -cpsi * sphi - spsi * ctheta * cphi, spsi * stheta],
                            [spsi * cphi + cpsi * ctheta * sphi, -spsi * sphi + cpsi * ctheta * cphi, -cpsi * stheta],
                            [stheta * sphi, stheta * cphi, ctheta]])

    return c_euler_out


def c_x(theta):
    """
    Gives the rotation matrix along the X-axis
    
    Manage iterable (list, array) as input
    
    [1,0, 0]
    [0,C,-S]
    [0,S, C]
    
    References
    ----------
        https://fr.wikipedia.org/wiki/Matrice_de_rotation#En_dimension_trois
        
    """

    if not utils.is_iterable(theta):
        theta = np.array([theta])

    c = np.cos(theta)
    s = np.sin(theta)
    z = np.zeros(len(theta))
    i = np.ones(len(theta))

    c_x_out = np.stack([[i, z, z],
                        [z, c, -s],
                        [z, s, c]])

    c_x_out = np.squeeze(c_x_out)

    return c_x_out


def c_y(theta):
    """
    Gives the rotation matrix around the Y-axis
    
    Manage iterable (list, array) as input
    
    [ C,0,S]
    [ 0,1,0]
    [-S,0,C]
    
    References
    ----------
        https://fr.wikipedia.org/wiki/Matrice_de_rotation#En_dimension_trois
    """

    if not utils.is_iterable(theta):
        theta = np.array([theta])

    c = np.cos(theta)
    s = np.sin(theta)
    z = np.zeros(len(theta))
    i = np.ones(len(theta))

    c_y_out = np.stack([[c, z, s],
                        [z, i, z],
                        [-s, z, c]])

    c_y_out = np.squeeze(c_y_out)

    return c_y_out


def c_z(theta):
    """
    Gives the rotation matrix around the Z-axis

    [C,-S,0]
    [S, C,0]
    [0, 0,1]
    
    References
    ----------
        https://fr.wikipedia.org/wiki/Matrice_de_rotation#En_dimension_trois
    """

    if not utils.is_iterable(theta):
        theta = np.array([theta])

    c = np.cos(theta)
    s = np.sin(theta)
    z = np.zeros(len(theta))
    i = np.ones(len(theta))

    c_z_out = np.stack([[c, -s, z],
                        [s, c, z],
                        [z, z, i]])

    c_z_out = np.squeeze(c_z_out)

    return c_z_out


def rot_quelconq(theta, vx, vy, vz, angtype='deg'):
    n = np.linalg.norm([vx, vy, vz])

    vx = vx / n
    vy = vy / n
    vz = vz / n

    if angtype == 'deg':
        theta = np.deg2rad(theta)

    c = np.cos(theta)
    s = np.sin(theta)

    r = np.array([[vx ** 2 + (1 - vx ** 2) * c, vx * vy * (1 - c) - vz * s, vx * vz * (1 - c) + vy * s],
                  [vx * vy * (1 - c) + vz * s, vy ** 2 + (1 - vy ** 2) * c, vy * vz * (1 - c) - vx * s],
                  [vx * vz * (1 - c) - vy * s, vy * vz * (1 - c) + vx * s, vz ** 2 + (1 - vz ** 2) * c]])

    return r


def vector_rpy(V, out_deg=True, ad_hoc_mode=False):
    """
    Donne le "Roll Pitch Yaw" pour un Point/Vecteur
    (les coordonnées spheriques plus exactement, le roll sera tjrs nul)

    le ad_hoc_mode est la pour corriger un bug que je n'explique pas
    pour l'utilisation de ads_offset
    (cette methode a été trouvé par les scripts ENURPY2 et ENURPY3)
    """

    x = V[0]
    y = V[1]
    z = V[2]

    n = np.linalg.norm(V)

    roll = 0
    pitch = np.arcsin(z / n)
    yaw = np.arctan2(y, x)

    if out_deg:
        roll = np.rad2deg(roll)
        pitch = np.rad2deg(pitch)
        yaw = np.rad2deg(yaw)

    if ad_hoc_mode:
        roll = -roll
        pitch = -pitch
        yaw = -yaw

        roll, pitch, yaw = pitch, roll, yaw

    return roll, pitch, yaw


def add_offset(direction, delta, point=None,
               out_delta_enu=False):
    """
    Un mode adhoc a du être implémenté dans la fonction élémentaire
    vector_RPY pour que ça puisse marcher correctement, reste à
    comprendre pourquoi ...

    out_delta_enu est prioritaire sur  les autres modes de sortie
    et renvoie le  delta

    par défaut on a le delta + le vecteur Direction
    et si Point est précisé : Point + Direction
    """

    rpy = vector_rpy(direction, ad_hoc_mode=True)

    crpy2enu = c_rpy2enu(*rpy)
    # Cenu2rpy = np.linalg.inv(C_rpy2enu(*rpy))

    delta_enu = crpy2enu.dot(c_enu2ned().dot(delta))

    if out_delta_enu:
        log.info(out_delta_enu)

    if not np.isclose(np.linalg.norm(delta), np.linalg.norm(delta_enu)):
        log.warning("np.linalg.norm(Delta) != np.linalg.norm(Delta_enu)")
        log.warning(np.linalg.norm(delta), np.linalg.norm(delta_enu), \
                    np.linalg.norm(delta) - np.linalg.norm(delta_enu))

    if not type(point) is None:
        outpoint = point + delta_enu
    else:
        outpoint = direction + delta_enu

    return outpoint
