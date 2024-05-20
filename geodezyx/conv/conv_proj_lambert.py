#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 15:11:12 2020

@author: psakic

This sub-module of geodezyx.conv deals with Lambert/Conform Conic projections.

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
from geodezyx import conv

log = logging.getLogger(__name__)

##########  END IMPORT  ##########


#  _____           _           _   _                _____           _        
# |  __ \         (_)         | | (_)              / ____|         | |       
# | |__) | __ ___  _  ___  ___| |_ _  ___  _ __   | |     __ _ _ __| |_ ___  
# |  ___/ '__/ _ \| |/ _ \/ __| __| |/ _ \| '_ \  | |    / _` | '__| __/ _ \ 
# | |   | | | (_) | |  __/ (__| |_| | (_) | | | | | |___| (_| | |  | || (_) |
# |_|   |_|  \___/| |\___|\___|\__|_|\___/|_| |_|  \_____\__,_|_|   \__\___/ 
#                _/ |                                                        
#               |__/        

def latitude_isometric(phi,e):
    """
    ALG0001 in NTG_71 (IGN Lambert Projection Tech document)
    """
    A = np.tan( (np.pi/4) + (phi/2) )
    B1 = 1 - e * np.sin(phi)
    B2 = 1 + e * np.sin(phi)
    
    L = np.log(A * ((B1/B2)**(.5 * e)))
    
    return L

def lambert_projection(long,lat,n,c,e,lc,Xs,Ys):
    """
    ALG0003 in NTG_71 (IGN Lambert Projection Tech document)
    """
    L = latitude_isometric(lat,e)
    X = Xs + c*np.exp(-n * L) * np.sin(n * (long - lc))
    Y = Ys - c*np.exp(-n * L) * np.cos(n * (long - lc))
    return X,Y

def lambert_secante_parameter(a , e , l0 , phi0 , phi1 , phi2 , X0 , Y0):
    """
    ALG0054 in NTG_71 (IGN Lambert Projection Tech document)
    """
    lc = l0
    A = np.log((conv.wnorm(phi2,a,e**2) * np.cos(phi2)) / (conv.wnorm(phi1,a,e**2) * np.cos(phi1)))
    B = latitude_isometric(phi1,e) - latitude_isometric(phi2,e)
    
    n = A/B
    
    C = ((conv.wnorm(phi1,a,e**2) * np.cos(phi1))/n) * np.exp(n * latitude_isometric(phi1,e))
        
    if np.isclose(phi0,np.pi/2):
        Xs = X0
        Ys = Y0  
    else:
        Xs = X0
        Ys = Y0 + C * np.exp(n * latitude_isometric(phi0,e))
            
    return e , n , C , lc , Xs , Ys

def lambert_projection_CC_frontend(long,lat,NZ = 93):
    """
    Project WGS84/ITRF Latitude/Longitude to Lambert Conique Conforme

    Parameters
    ----------
    long,lat : float
        WGS84/ITRF Longitude/Latitude
        
    NZ : int
        Lambert93, NZ = 0 or = 93
                
    Returns
    -------
    X,Y : float
        Projected coordinates
    
    Source
    ------
    https://fr.wikipedia.org/wiki/Projection_conique_conforme_de_Lambert   
    https://geodesie.ign.fr/contenu/fichiers/documentation/algorithmes/notice/NTG_71.pdf
    https://geodesie.ign.fr/contenu/fichiers/documentation/rgf93/Lambert-93.pdf
 
    """
    if NZ in (0,93):
        phi0 = np.deg2rad(46.5)
        phi1 = np.deg2rad(44)
        phi2 = np.deg2rad(49)
        l0 = np.deg2rad(3)
        X0 = 700000.
        Y0 = 6600000.
        a = 6378137.
        e = 0.0818191910435
    else:        
        phi0 = np.deg2rad(41 + NZ)
        phi1 = np.deg2rad(phi0 - 0.75)
        phi2 = np.deg2rad(phi0 + 0.75)
        l0 = np.deg2rad(3)
        X0 = 1700000.
        Y0 = (NZ * 1000000.) + 200000.
        a = 6378137.
        e = 0.0818191910435
    
    e , n , c , lc , Xs , Ys = lambert_secante_parameter(a, e, l0, phi0, phi1, phi2, X0, Y0)
    X,Y = lambert_projection(np.deg2rad(long),np.deg2rad(lat),n,c,e,lc,Xs,Ys)
    
    return X,Y
