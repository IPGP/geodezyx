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

# import scipy
# from pyorbital import astronomy
# import re
log = logging.getLogger(__name__)


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

### It is based on the book of Grewal, M. S., Weill, L. R., & Andrews, A. P. (2007). 
### Global Positioning Systems, Inertial Navigation, and Integration. 
### Hoboken, NJ, USA : John Wiley and Sons, Inc.

def C_2D(theta, angtype='deg'):
    if angtype == 'deg':
        theta = np.deg2rad(theta)

    return np.array([[np.cos(theta), -np.sin(theta)],
                     [np.sin(theta), np.cos(theta)]])


def C_ned2enu():
    return np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])


def C_enu2ned():
    return np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])


def C_ned2ecef(phi, llambda, angtype='deg'):
    if angtype == 'deg':
        phi = np.deg2rad(phi)
        llambda = np.deg2rad(llambda)

    In_ned = np.array([1, 0, 0])
    Ie_ned = np.array([0, 1, 0])
    Id_ned = np.array([0, 0, 1])

    Ix_ned = np.array([-np.cos(llambda) * np.sin(phi), -np.sin(llambda), -np.cos(llambda) * np.cos(phi)])
    Iy_ned = np.array([-np.sin(llambda) * np.sin(phi), np.cos(llambda), -np.sin(llambda) * np.cos(phi)])
    Iz_ned = np.array([np.cos(phi), 0, -np.sin(phi)])

    C_ned2ecef_out = np.array([[np.dot(Ix_ned, In_ned), np.dot(Ix_ned, Ie_ned), np.dot(Ix_ned, Id_ned)],
                               [np.dot(Iy_ned, In_ned), np.dot(Iy_ned, Ie_ned), np.dot(Iy_ned, Id_ned)],
                               [np.dot(Iz_ned, In_ned), np.dot(Iz_ned, Ie_ned), np.dot(Iz_ned, Id_ned)]])

    return C_ned2ecef_out


def C_ecef2ned(phi, llambda, angtype='deg'):
    return np.linalg.inv(C_ned2ecef(phi, llambda, angtype))


def C_enu2ecef(phi, llambda, angtype='deg'):
    if angtype == 'deg':
        phi = np.deg2rad(phi)
        llambda = np.deg2rad(llambda)

    Ie_enu = np.array([1, 0, 0])
    In_enu = np.array([0, 1, 0])
    Iu_enu = np.array([0, 0, 1])

    Ix_enu = np.array([-np.sin(llambda), -np.cos(llambda) * np.sin(phi), np.cos(llambda) * np.cos(phi)])
    Iy_enu = np.array([np.cos(llambda), -np.sin(llambda) * np.sin(phi), np.sin(llambda) * np.cos(phi)])
    Iz_enu = np.array([0, np.cos(phi), np.sin(phi)])

    C_enu2ecef_out = np.array([[np.dot(Ix_enu, Ie_enu), np.dot(Ix_enu, In_enu), np.dot(Ix_enu, Iu_enu)],
                               [np.dot(Iy_enu, Ie_enu), np.dot(Iy_enu, In_enu), np.dot(Iy_enu, Iu_enu)],
                               [np.dot(Iz_enu, Ie_enu), np.dot(Iz_enu, In_enu), np.dot(Iz_enu, Iu_enu)]])

    return C_enu2ecef_out


def C_ecef2enu(phi, llambda, angtype='deg'):
    return C_enu2ecef(phi, llambda, angtype).T


def C_ecef2enu_sigma(phi, llambda, angtype='deg'):
    """
    Gives the transformation matrix between ECEF and ENU (East North Up) ref. frame
    for the sigmas (variance propagation)


    Source
    ------
    https://stackoverflow.com/questions/51162460/converting-ecef-xyz-covariance-matrix-to-enu-covariance-matrix
    """
    if angtype == 'deg':
        phi = np.deg2rad(phi)
        llambda = np.deg2rad(llambda)


    C_ecef2enu_sigma_out = np.array([[- np.sin(llambda), np.cos(llambda), 0],
                                     [-np.sin(phi) * np.cos(llambda), -np.sin(phi) * np.sin(llambda), np.cos(phi)],
                                     [np.cos(phi) * np.cos(llambda), np.cos(phi) * np.sin(llambda), np.sin(llambda)]])

    return C_ecef2enu_sigma_out


def C_rpy2enu(roll, pitch, yaw, angtype='deg'):
    if angtype == 'deg':
        roll = np.deg2rad(roll)
        pitch = np.deg2rad(pitch)
        yaw = np.deg2rad(yaw)

    R = roll
    P = pitch
    Y = yaw

    Ir = np.array([np.sin(Y) * np.cos(P), np.cos(Y) * np.cos(P), np.sin(P)])
    Ip = np.array([np.cos(R) * np.cos(Y) + np.sin(R) * np.sin(Y) * np.sin(P),
                   -np.cos(R) * np.sin(Y) + np.sin(R) * np.cos(Y) * np.sin(P),
                   -np.sin(R) * np.cos(P)])
    Iy = np.array([-np.sin(R) * np.cos(Y) + np.cos(R) * np.sin(Y) * np.sin(P),
                   np.sin(R) * np.sin(Y) + np.cos(R) * np.cos(Y) * np.sin(P),
                   - np.cos(R) * np.cos(P)])

    #    Ie = np.array([ np.sin(Y) * np.cos(P) ,
    #    np.cos(R) * np.cos(Y) + np.sin(R) * np.sin(Y) * np.sin(P) ,
    #    -np.sin(R) * np.cos(Y) + np.cos(R) * np.sin(Y) * np.sin(P) ])
    #
    #    In = np.array([ np.cos(Y) * np.cos(P),
    #     -np.cos(R) * np.sin(Y) + np.sin(R) * np.cos(Y) * np.sin(P),
    #     np.sin(R) * np.sin(Y) + np.cos(R) * np.cos(Y) * np.sin(P) ])
    #
    #    Iu = np.array([ np.sin(P),
    #     -np.sin(R) * np.cos(P),
    #     -np.cos(R) * np.cos(P) ])

    C_rpy2enu_out = np.hstack((Ir[np.newaxis].T, Ip[np.newaxis].T, Iy[np.newaxis].T))

    return C_rpy2enu_out


def C_rpy2enu2(roll, pitch, yaw, angtype='deg'):
    if angtype == 'deg':
        roll = np.deg2rad(roll)
        pitch = np.deg2rad(pitch)
        yaw = np.deg2rad(yaw)

    R = roll
    P = pitch
    Y = yaw

    Sr = np.sin(R)
    Cr = np.cos(R)
    Sp = np.sin(P)
    Cp = np.cos(P)
    Sy = np.sin(Y)
    Cy = np.cos(Y)

    C_rpy2enu_out = np.array([[Sy * Cp, Cr * Cy + Sr * Sy * Sp, - Sr * Cy + Cr * Sy * Sp],
                              [Cy * Cp, -Cr * Sy + Sr * Cy * Sp, Sr * Sy + Cr * Cy * Sp],
                              [Sp, -Sr * Cp, -Cr * Cp]])

    return C_rpy2enu_out


def C_eci2rtn(P, V):
    """
    gives the transformation matrix between ECI ref. frame and RTN (aka RIC or RSW)

    Parameters
    ----------
    P : numpy.array
        3D vector, position of the object in ECI frame

    V : numpy.array
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

    H = np.cross(P, V)

    R = P / np.linalg.norm(P)
    C = H / np.linalg.norm(H)
    I = np.cross(C, R)

    C_eci2rtn_out = np.column_stack((R, I, C))

    return C_eci2rtn_out


def C_rtn2rpy():
    C_rtn2rpy_out = np.array([[0., 0, -1.],
                              [1., 0., 0.],
                              [0, -1., 0.]])
    return C_rtn2rpy_out


def C_cep2itrs(xpole, ypole):
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

    C_cep2itrs_out = np.array([[1., 0., xpole2],
                               [0., 1., -ypole2],
                               [-xpole2, ypole2, 1.]])

    return C_cep2itrs_out


def C_euler(phi, theta, psi):
    """
    Gives the matrix of an Euler rotation
        
    References
    ----------
        https://fr.wikipedia.org/wiki/Angles_d%27Euler
    """
    Cphi = np.cos(phi)
    Ctheta = np.cos(theta)
    Cpsi = np.cos(psi)
    Sphi = np.sin(phi)
    Stheta = np.sin(theta)
    Spsi = np.sin(psi)

    C_euler_out = np.array([[Cpsi * Cphi - Spsi * Ctheta * Sphi, -Cpsi * Sphi - Spsi * Ctheta * Cphi, Spsi * Stheta],
                            [Spsi * Cphi + Cpsi * Ctheta * Sphi, -Spsi * Sphi + Cpsi * Ctheta * Cphi, -Cpsi * Stheta],
                            [Stheta * Sphi, Stheta * Cphi, Ctheta]])

    return C_euler_out


def C_x(theta):
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

    C = np.cos(theta)
    S = np.sin(theta)
    Z = np.zeros(len(theta))
    I = np.ones(len(theta))

    C_x_out = np.stack([[I, Z, Z],
                        [Z, C, -S],
                        [Z, S, C]])

    C_x_out = np.squeeze(C_x_out)

    return C_x_out


def C_y(theta):
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

    C = np.cos(theta)
    S = np.sin(theta)
    Z = np.zeros(len(theta))
    I = np.ones(len(theta))

    C_y_out = np.stack([[C, Z, S],
                        [Z, I, Z],
                        [-S, Z, C]])

    C_y_out = np.squeeze(C_y_out)

    return C_y_out


def C_z(theta):
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

    C = np.cos(theta)
    S = np.sin(theta)
    Z = np.zeros(len(theta))
    I = np.ones(len(theta))

    C_z_out = np.stack([[C, -S, Z],
                        [S, C, Z],
                        [Z, Z, I]])

    C_z_out = np.squeeze(C_z_out)

    return C_z_out


def rot_quelconq(theta, Vx, Vy, Vz, angtype='deg'):
    N = np.linalg.norm([Vx, Vy, Vz])

    Vx = Vx / N
    Vy = Vy / N
    Vz = Vz / N

    if angtype == 'deg':
        theta = np.deg2rad(theta)

    c = np.cos(theta)
    s = np.sin(theta)

    R = np.array([[Vx ** 2 + (1 - Vx ** 2) * c, Vx * Vy * (1 - c) - Vz * s, Vx * Vz * (1 - c) + Vy * s],
                  [Vx * Vy * (1 - c) + Vz * s, Vy ** 2 + (1 - Vy ** 2) * c, Vy * Vz * (1 - c) - Vx * s],
                  [Vx * Vz * (1 - c) - Vy * s, Vy * Vz * (1 - c) + Vx * s, Vz ** 2 + (1 - Vz ** 2) * c]])

    return R


def vector_RPY(V, out_deg=True, ad_hoc_mode=False):
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


def add_offset(Direction, Delta, Point=None,
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

    rpy = vector_RPY(Direction, ad_hoc_mode=True)

    Crpy2enu = C_rpy2enu(*rpy)
    # Cenu2rpy = np.linalg.inv(C_rpy2enu(*rpy))

    Delta_enu = Crpy2enu.dot(C_enu2ned().dot(Delta))

    if out_delta_enu:
        log.info(out_delta_enu)

    if not np.isclose(np.linalg.norm(Delta), np.linalg.norm(Delta_enu)):
        log.warning("np.linalg.norm(Delta) != np.linalg.norm(Delta_enu)")
        log.warning(np.linalg.norm(Delta), np.linalg.norm(Delta_enu), \
                    np.linalg.norm(Delta) - np.linalg.norm(Delta_enu))

    if not type(Point) is None:
        outpoint = Point + Delta_enu
    else:
        outpoint = Direction + Delta_enu

    return outpoint
