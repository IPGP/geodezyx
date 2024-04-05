#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 20:37:05 2023

@author: psakic

This sub-module of geodezyx.conv deals with UTM projections.

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

#### geodeZYX modules
#from geodezyx import utils,conv

#### Import the logger
import logging
log = logging.getLogger(__name__)

##########  END IMPORT  ##########

    
def utm_geo2xy(lat,lon,ellips="wgs84",zone=None): 
    """
    Latitude/Longitude to Universal Transverse Mercator(UTM)
    coordinates precise (mm-level) conversion.
    converts coordinates
    Latitude/Longitude (in degrees) to UTM X and Y (in meters). 
    Default ellipsoid is WGS84.
    
    Parameters
    ----------
    lat : float or numpy array
        Latitude (in degrees).
    lon : float or numpy array
        Longitude (in degrees).
    ellips : str or 2-tuple, optional
        uses specific ellipsoid for conversion. ellipsoid can be one
    	of the following strings
		'wgs84': World Geodetic System 1984 (default)
		'nad27': North American ellips 1927
		'clk66': Clarke 1866
		'nad83': North American ellips 1983
		'grs80': Geodetic Reference System 1980
		'int24': International 1924 / Hayford 1909
        OR
        can be a 2-element tuple [A,F] where A is semimajor axis (in meters)
        F is flattening of the user-defined ellipsoid.
        The default is "wgs84".
    zone : int, optional
        forces the UTM ZONE (scalar integer or same size as
        LAT and LON) instead of automatic set. 
        The default is None.

    Returns
    -------
    x,y : float or numpy array
        UTM X and Y coordinates (in meters).
    f : float or numpy array
        returns also the computed UTM ZONE (negative
        value for southern hemisphere points).
        
    Notes
    -----
    It does not perform cross-ellips conversion.
    Precision is near a millimeter.
    
    This function implements the IGN algorithm.
    This one is *more accurate* than `pyproj` or `utm` python lib 
    because it implements more coefficients.
    
    Source
    ------
    Adapted from WebObs source code, François Beauducel et al. 
    https://github.com/IPGP/webobs/blob/master/CODE/matlab/ll2utm.m
    
	I.G.N., Projection cartographique Mercator Transverse: Algorithmes,
    Notes Techniques NT/G 76, janvier 1995.    
    https://geodesie.ign.fr/contenu/fichiers/documentation/algorithmes/notice/NTG_76.pdf
    
    """   
    
    # Available ellipsoids
    ellips_dic = dict()
    ellips_dic['wgs84'] = np.array([6378137.0,298.257223563])
    ellips_dic['nad83'] = np.array([6378137.0,298.257222101])
    ellips_dic['grs80'] = np.array([6378137.0,298.257222101])
    ellips_dic['nad27'] = np.array([6378206.4,294.978698214])
    ellips_dic['int24'] = np.array([6378388.0,297.0])
    ellips_dic['clk66'] = np.array([6378206.4,294.978698214])

    # constants
    D0 = 180 / np.pi
    K0 = 0.9996
    X0 = 500000
    
    # defaults
    zone = []
    
    if type(ellips) is str:
        # LL2UTM(...,ellips) with ellips as char
        if not ellips in ellips_dic.keys():
            raise Exception('Unkown ellips name "%s"',ellips)
        A1 = ellips_dic[ellips][0]
        F1 = ellips_dic[ellips][1]
    else:
        # LL2UTM(...,ellips) with ellips as [A,F] user-defined
        A1 = ellips[1]
        F1 = ellips[2]
    
    p1 = lat / D0
    l1 = lon / D0
    
    # UTM zone automatic setting
    if not zone:
        F0 = np.round((l1 * D0 + 183) / 6)
    else:
        F0 = zone
    
    B1 = A1 * (1 - 1 / F1)
    E1 = np.sqrt((A1 * A1 - B1 * B1) / (A1 * A1))
    P0 = 0 / D0
    L0 = (6 * F0 - 183) / D0
    
    Y0 = 10000000.0 * (p1 < 0)
    
    N = K0 * A1
    C = utm_coef(E1,0)
    B = C[0] * P0 + C[1] * np.sin(2 * P0) + C[2] * np.sin(4 * P0) + C[3] * np.sin(6 * P0) + C[4] * np.sin(8 * P0)
    YS = Y0 - N * B
    C = utm_coef(E1,2)
    L = np.log(np.multiply(np.tan(np.pi / 4 + p1 / 2),(((1 - E1 * np.sin(p1)) / (1 + E1 * np.sin(p1))) ** (E1 / 2))))
    z = complex(np.arctan(np.sinh(L) / np.cos(l1 - L0)),np.log(np.tan(np.pi / 4 + np.arcsin(np.sin(l1 - L0) / np.cosh(L)) / 2)))
    Z = np.multiply(np.multiply(N,C[0]),z) + np.multiply(N,(C[1] * np.sin(2 * z) + C[2] * np.sin(4 * z) + C[3] * np.sin(6 * z) + C[4] * np.sin(8 * z)))
    xs = np.imag(Z) + X0
    ys = np.real(Z) + YS
    
    f = np.multiply(F0,np.sign(lat))
    fu = np.unique(f)
    if np.isscalar(fu):
        f = fu
    
    x = xs
    y = ys
        
    return x,y,f
    
    ###########################################################################
    
def utm_coef(e = None,m = 0): 
    """
    Projection coefficients
	returns a vector of 5 coefficients based on the input parameters

    Parameters
    ----------
    e : TYPE, optional
        first ellipsoid excentricity. The default is None.
    m : int, optional
        Describes the expected coefficients.
		M = 0 for transverse mercator
		M = 1 for transverse mercator reverse coefficients
		M = 2 for merdian arc
        The default is 0.

    Returns
    -------
    c : numpy array
        Coefficient Matrix.

    """

    if 0 == m:
        c0 = np.array([[- 175 / 16384,0,- 5 / 256,0,- 3 / 64,0,- 1 / 4,0,1],
                       [- 105 / 4096,0,- 45 / 1024,0,- 3 / 32,0,- 3 / 8,0,0],
                       [525 / 16384,0,45 / 1024,0,15 / 256,0,0,0,0],
                       [- 175 / 12288,0,- 35 / 3072,0,0,0,0,0,0],
                       [315 / 131072,0,0,0,0,0,0,0,0]])
    elif 1 == m:
        c0 = np.array([[- 175 / 16384,0,- 5 / 256,0,- 3 / 64,0,- 1 / 4,0,1],
                       [1 / 61440,0,7 / 2048,0,1 / 48,0,1 / 8,0,0],
                       [559 / 368640,0,3 / 1280,0,1 / 768,0,0,0,0],
                       [283 / 430080,0,17 / 30720,0,0,0,0,0,0],
                       [4397 / 41287680,0,0,0,0,0,0,0,0]])
    elif 2 == m:
        c0 = np.array([[- 175 / 16384,0,- 5 / 256,0,- 3 / 64,0,- 1 / 4,0,1],
                       [- 901 / 184320,0,- 9 / 1024,0,- 1 / 96,0,1 / 8,0,0],
                       [- 311 / 737280,0,17 / 5120,0,13 / 768,0,0,0,0],
                       [899 / 430080,0,61 / 15360,0,0,0,0,0,0],
                       [49561 / 41287680,0,0,0,0,0,0,0,0]])


    c = np.polyval(c0.T,e)
    
    return c


