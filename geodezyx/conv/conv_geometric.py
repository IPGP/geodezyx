#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 15:11:12 2020

@author: psakicki

This sub-module of geodezyx.conv deals with low-level geometric operation.

it can be imported directly with:
from geodezyx import conv

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
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

def dist(A,B):
    """
    Cartesian distance between 2 points A & B

    Parameters
    ----------
    A,B : numpy.array of float
        3D-points
        
    Returns
    -------
    D : float
        distance between A & B
    """
    A = np.array(A)
    B = np.array(B)
    return np.linalg.norm(A - B)

def dist_diff(A,B):
    """
    First derivative of cartesian distance between 2 points A & B

    Parameters
    ----------
    A, B : numpy.array of float
        "points", can be 2D or 3D vectors (list, np.array ...)

    Returns
    -------
    diffA : numpy.array of float 
        [ dD/dxa , dD/dya , dD/dza ]

    diffB : numpy.array of float
        [ dD/dxb , dD/dyb , dD/dzb ] = -diffA
    """

    dAB   = A-B
    dist  = np.linalg.norm(dAB)

    diffA =   dAB / dist
    diffB = - dAB / dist

    return diffA, diffB

def relative_orientation(x1,y1,x2,y2,trigo_orient = True):
    """
    Return the angle between a point 1 and a point 2
    (reference vector : to the right)
    """
    if trigo_orient:
        ang1 = np.mod(360 - np.rad2deg(np.arctan2((x2 - x1),(y2 - y1))),360)
    else:
        ang1 = np.mod(np.rad2deg(np.arctan2((x2 - x1),(y2 - y1))),360)
    return ang1

def barycenter(points_list_in):
    """
    Determines the barycenter of a list of points
    """
    points_arr = np.array(points_list_in)
    return np.mean(points_arr[:,-3:], axis=0)

def pythagore(a,b,c=0):
    """
    Computes Pythagore's formula
    """
    return np.sqrt(a**2 + b**2 + c**2)

def equilateral_triangle(side):
    """
    Gives coordinates of an equilateral triangle of a given side
    """
    hauteur_len = np.sqrt(side**2 - (side * .5) **2)
    A = np.array([0,(2./3.) * hauteur_len])
    B = np.array([side / 2. , - (1./3.) * hauteur_len])
    C = np.array([- side / 2. , - (1./3.) * hauteur_len])

    return A,B,C


def vincenty_full(point1, point2, miles=False,full=True,azimuth_in_deg=True):
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
    
    fwdAz,revAz : float
        Forward and Reverse Azimuth between the 2 points 
        
    References
    ----------
    https://github.com/maurycyp/vincenty/blob/master/vincenty/
    
    Examples
    --------
    >>> vincenty((0.0, 0.0), (0.0, 0.0))  # coincident points
    0.0
    >>> vincenty((0.0, 0.0), (0.0, 1.0))
    111.319491
    >>> vincenty((0.0, 0.0), (1.0, 0.0))
    110.574389
    >>> vincenty((0.0, 0.0), (0.5, 179.5))  # slow convergence
    19936.288579
    >>> vincenty((0.0, 0.0), (0.5, 179.7))  # failure to converge
    >>> boston = (42.3541165, -71.0693514)
    >>> newyork = (40.7791472, -73.9680804)
    >>> vincenty(boston, newyork)
    298.396057
    >>> vincenty(boston, newyork, miles=True)
    185.414657
    """
    
    # WGS 84
    a = 6378137  # meters
    f = 1 / 298.257223563
    b = 6356752.314245  # meters; b = (1 - f)a
    
    MILES_PER_KILOMETER = 0.621371
    
    MAX_ITERATIONS = 200
    CONVERGENCE_THRESHOLD = 1e-12  # .000,000,000,001

    # short-circuit coincident points
    if point1[0] == point2[0] and point1[1] == point2[1]:
        return 0.0

    import math 
    
    U1 = math.atan((1 - f) * math.tan(math.radians(point1[0])))
    U2 = math.atan((1 - f) * math.tan(math.radians(point2[0])))
    L = math.radians(point2[1] - point1[1])
    Lambda = L

    sinU1 = math.sin(U1)
    cosU1 = math.cos(U1)
    sinU2 = math.sin(U2)
    cosU2 = math.cos(U2)
    sinL  = math.sin(L)
    cosL  = math.cos(L)

    for iteration in range(MAX_ITERATIONS):
        sinLambda = math.sin(Lambda)
        cosLambda = math.cos(Lambda)
        sinSigma = math.sqrt((cosU2 * sinLambda) ** 2 +
                             (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) ** 2)
        if sinSigma == 0:
            return 0.0  # coincident points
        cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda
        sigma = math.atan2(sinSigma, cosSigma)
        sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma
        cosSqAlpha = 1 - sinAlpha ** 2
        try:
            cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha
        except ZeroDivisionError:
            cos2SigmaM = 0
        C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha))
        LambdaPrev = Lambda
        Lambda = L + (1 - C) * f * sinAlpha * (sigma + C * sinSigma *
                                               (cos2SigmaM + C * cosSigma *
                                                (-1 + 2 * cos2SigmaM ** 2)))
        if abs(Lambda - LambdaPrev) < CONVERGENCE_THRESHOLD:
            break  # successful convergence
    else:
        return None  # failure to converge

    uSq = cosSqAlpha * (a ** 2 - b ** 2) / (b ** 2)
    A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
    B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))
    deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 * (cosSigma *
                 (-1 + 2 * cos2SigmaM ** 2) - B / 6 * cos2SigmaM *
                 (-3 + 4 * sinSigma ** 2) * (-3 + 4 * cos2SigmaM ** 2)))
    s = b * A * (sigma - deltaSigma)

    s /= 1000  # meters to kilometers
    if miles:
        s *= MILES_PER_KILOMETER  # kilometers to miles


    if not full:
        return s
    else:        
        fwdAz = np.arctan2(cosU2*sinL,  cosU1*sinU2-sinU1*cosU2*cosL)
        revAz = np.arctan2(cosU1*sinL, -sinU1*cosU2+cosU1*sinU2*cosL)
       
        if azimuth_in_deg:
            fwdAz = np.rad2deg(fwdAz)
            revAz = np.rad2deg(revAz)
       
        return s,fwdAz,revAz
    
def orthogonal_projection(Xa,Xb,Xv):
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
        
    Note
    ----
    Misc. Notes

    Source
    ------
    https://fr.wikipedia.org/wiki/Projection_orthogonale
    """
    
    xa , ya = Xa
    xb , yb = Xb    
    xv , yv = Xv
    
    D =np.sqrt(xv**2 + yv **2) 
    BH = ((xa - xb)*xv + (ya - yb)*yv) / D
    
    xh = xb + (BH / D) * xv
    yh = yb + (BH / D) * yv
    
    Xh = np.array([xh,yh])
    return Xh , dist(Xa,Xh)


def line_maker(x1,y1,x2,y2,nbpts=10000):
    """
    Determine points of a line easily
    
    Parameters
    ----------
    x1,y1 : float
        Coordinates of the start point

    x2,y2 : float
        Coordinates of the end point
                
    Returns
    -------
    X,Y : numpy.array of float
        points of the line
    """

    X = np.linspace(x1,x2,nbpts)
    Y = np.linspace(y1,y2,nbpts)
    return X,Y

