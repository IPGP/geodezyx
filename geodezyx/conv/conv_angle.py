#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 15:11:12 2020

@author: psakic

This sub-module of geodezyx.conv deals with angle conversion.

it can be imported directly with:
from geodezyx import conv

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU GPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (GFZ, pierre.sakic@gfz-postdam.de)
GitHub repository :
https://github.com/GeodeZYX/GeodeZYX-Toolbox_v4
"""

########## BEGIN IMPORT ##########
#### External modules
import numpy as np
import scipy
from pyorbital import astronomy
import re

#### geodeZYX modules
from geodezyx import utils

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names

##########  END IMPORT  ##########

#                       _         _____                              _             
#     /\               | |       / ____|                            (_)            
#    /  \   _ __   __ _| | ___  | |     ___  _ ____   _____ _ __ ___ _  ___  _ __  
#   / /\ \ | '_ \ / _` | |/ _ \ | |    / _ \| '_ \ \ / / _ \ '__/ __| |/ _ \| '_ \ 
#  / ____ \| | | | (_| | |  __/ | |___| (_) | | | \ V /  __/ |  \__ \ | (_) | | | |
# /_/    \_\_| |_|\__, |_|\___|  \_____\___/|_| |_|\_/ \___|_|  |___/_|\___/|_| |_|
#                  __/ |                                                           
#                 |___/                                                            

### Angle conversion
    
def dms2dec_num(deg,minn=0,sec=0):
    """
    Angle conversion
    
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
    return deg + minn * (1./60.) +  sec * (1./3600.)


def deg_dec2dms(deg_in,only_dm=False):
    """
    Angle conversion
    
    Convert :
    Decimal Degree => Degree Minute Seconds
        
    Parameters
    ----------
    deg_in : float or iterable of floats
        Decimal degrees

    Returns
    -------
    deg , minu , sec : numpy.array of float
        3 arrays for Degrees, Minutes, Seconds
    """
    deg              = np.floor(deg_in)
    decimal_part     = deg_in - deg 
    decimal_part_sec = decimal_part * 3600
    minu             = np.floor_divide(decimal_part_sec,60)
    sec              = decimal_part_sec - minu * 60
    sec              = np.round(sec,8)
    if not only_dm:
        return deg , minu , sec
    else:
        return deg , minu + sec * (1./60.)
        



def dms2dec(dms_str , onlyDM=False):
    """   
    Angle conversion

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

    onlyDM : bool
        True if the string is only a Degree Minute Angle 
        e.g.
        "40°52.0931'N"
        "28°31.4136'E"
                
    Returns
    -------
    dd_float : float
        Decimal degree Angle
    """
    
    if not onlyDM:
        print("WARN : DMS mode not well implemented yet !!! ")
    
    
    dms_str = dms_str.strip()

    if re.match('.*[swoSWO].*', dms_str):
        sign = -1
    else:
        sign = 1

    # former regex
    lis = re.findall(r'\d+|\D+', dms_str)

    ipt = lis.index('.')
    decimal = lis[ipt-1] +  lis[ipt] + lis[ipt+1]
    lis[ipt-1] = decimal
    lis = lis[:ipt]

    for i,e in enumerate(lis):
        try:
            lis[i] = float(e)
        except:
            lis.remove(e)

    lis = [float(e) for e in lis]
    print(lis)

    deg = lis[0]
    minu = lis[1]
    if onlyDM:
        sec = 0.
    else:
        try:
            sec  = lis[2]
        except IndexError:
            print("ERR : did you forgot to activate the DM only mode ? ;) ?")
            return None

    dd_float = sign * (deg + minu / 60. + sec / 3600.)
    return dd_float

def arcsec2deg(arcsec_in):
    """
    Angle conversion
    
    Convert :
    Arcsecond => Degrees
    
    NB : Not really useful, should be avoided
    """
        
    return arcsec_in / 3600.

def deg2arcsec(deg_in):
    """
    Angle conversion
    
    Convert :
    Degrees => Arcsecond
    
    NB : Not really useful, should be avoided
    """
    return deg_in * 3600.


def angle2equivalent_earth_radius(angle_in,angtype='deg',
                                  earth_radius = 6371008.8):
    """
    Quick and simple function which gives the equivalent distance on a 
    Earth great circle of an angle
    
    angle can be : "deg", "rad", "mas"
    """
    earth_circum = earth_radius * 2 * np.pi
    if angtype == "deg":
        equiv_out = (angle_in * earth_circum) / 360.
    elif angtype == "rad":
        equiv_out = (angle_in * earth_circum) / (np.pi *2)
    elif angtype == "mas":
        equiv_out = (angle_in * 10**-3 * earth_circum) / (86400.)  
        
    return equiv_out
        
def anglesfromvects(xa,ya,xb,yb,angtype='deg'):
    """
    Determines the angle between the points A(xa,ya) and B(xb,yb)
    
    angle can be : "deg", "rad"
    """
    A = np.array([xa,ya])
    B = np.array([xb,yb])

    ps = np.inner(A,B)
    a = np.arccos( ps / (np.linalg.norm(A) * np.linalg.norm(B)) )

    if angtype == 'deg':
        return np.rad2deg(a)
    elif angtype == 'rad':
        return a
    else:
        print('ERR : angfromvects : mauvais angtype')

def angle_interpolation_quick(A,B,w):
    """
    Determine the interpolation between angle A & B
    by conversion to the cartesian space
    the parameter w € [0,1] define the interpoled angle C(w)
    where C(w=0) = A  &  C(w=1) = B

    Source
    ------
    https://stackoverflow.com/questions/2708476/rotation-interpolation
    """
    CS = (1-w)*np.cos(A) + w*np.cos(B)
    SN = (1-w)*np.sin(A) + w*np.sin(B)
    C = np.atan2(SN,CS)

    return C

def angle_from_3_pts(p1,p2,p3):
    """
    Gives angle between 3 points
    p3 is the vertex (sommet)
    
    Source
    ------
    http://www.les-mathematiques.net/phorum/read.php?8,596072,596231
    """

    p1,p2,p3 = [np.array(e) for e in (p1,p2,p3)]
    x1 , y1 = p1
    x2 , y2 = p2
    x3 , y3 = p3

    kos = ((x1 - x3)*(x2 - x3) + (y1 - y3)*(y2 - y3)) / (np.linalg.norm(p1 - p3) * np.linalg.norm(p2 - p3))

    return np.arccos(kos)


def cartesian2polar(x,y):
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

    theta = np.arctan2(y,x)
    r     = np.sqrt(x** + y**2)
    return r , theta

def polar2cartesian(r,theta,ang='deg'):
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

    if ang == 'deg':
        theta = np.deg2rad(theta)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x , y