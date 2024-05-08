#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.reffram contains functions for quaternion 
operations. 

it can be imported directly with:
from geodezyx import reffram

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""


########## BEGIN IMPORT ##########
#### External modules
import numpy as np
from scipy.spatial.transform import Rotation

#### geodeZYX modules
from geodezyx import utils


##########  END IMPORT  ##########


#   ____              _                  _
#  / __ \            | |                (_)
# | |  | |_   _  __ _| |_ ___ _ __ _ __  _  ___  _ __  ___
# | |  | | | | |/ _` | __/ _ \ '__| '_ \| |/ _ \| '_ \/ __|
# | |__| | |_| | (_| | ||  __/ |  | | | | | (_) | | | \__ \
#  \___\_\\__,_|\__,_|\__\___|_|  |_| |_|_|\___/|_| |_|___/
#



def quaternion(r,p,h,angtype='rad',euler_sequence="zyx"):
    """
    Args:
        roll (xaxis) pitch (yaxis) heading (zaxis) angle

        angtype : 'rad' or 'deg'

    Returns:
        the associate quaternion

    Note:
        the input args are in XYZ but
        the rotation is in the order ZYX

        Need of the cgkit quaternion module
        http://cgkit.sourceforge.net/doc2/quat.html

    """
    if angtype == 'deg':
        r,p,h = np.deg2rad(r),np.deg2rad(p),np.deg2rad(h)
        
        
    #### Old style cgkit DISCONTINUED
    #M = cgt.mat3().fromEulerZYX(r,p,h)
    #q = cgt.quat().fromMat(M)

    R = Rotation.from_euler(euler_sequence,[r,p,h],degrees=False)    
    q = R.as_quat()
        
    return q

def quatern_multi(rlist,plist,hlist,angtype='rad'):
    """
    wrapper of quaternion, with lists/arrays of angles
    """

    qlisout = []
    for r,p,h in zip(rlist,plist,hlist):
        q = quaternion(r,p,h,angtype=angtype)
        qlisout.append(q)
    qarrout = np.array(qlisout)
    return qarrout

def rotate_quat2(ptin,quatin):
    """
    DISCONTINUED
    plus long que rotate_quat :
    benchmarking 5.53084397316 vs 0.703444957733 pour 100000
    """
    return np.array(quatin.rotateVec(ptin))

def rotate_quat3(ptin,quatin):
    """
    DISCONTINUED
    use the old quaterinons toolbox
    """
    pt1 = np.concatenate(([0],np.array(ptin)))
    pt2 = cgt.quat(*pt1)
    
    qtmp = (quatin) * pt2 * (quatin.inverse())
    
    return np.array([qtmp.x, qtmp.y , qtmp.z])


def rotate_quat(ptin,quatin):
    Robj = 999999999999999999
    return np.concatenate(([0],np.array(ptin)))


def rotate_quat_multi(pts_list_in, quats_list_in):
    """
    Multi version of rotate_quat

    Args:
        pts_list_in : a list of Points object
        quats_list_in : a list of quaternions

    Returns:
        a list of list of rotated points
    """

    pts_lis_out = []
    for pt in pts_list_in:
        temp_pt_lis = []
        for q in quats_list_in:
            ptmp = rotate_quat(pt,q)
            temp_pt_lis.append(ptmp)
        pts_lis_out.append(temp_pt_lis)
    pts_arr_out = np.array(pts_lis_out)
    return pts_arr_out

def interp_quat(Tlis , Quatlis , t):
    """
    Interpolate quaterinon

    Args:
        Tlis (float) : list of time, in POSIX time

        Quatlis : list of quaternion (produced by quaterion fct)

        t : the instant for the interpolation

    Returns:
        an interpoled Quaternion


    """
    i1 , i2 = utils.find_interval_bound(Tlis,t)
    q1 , q2 = Quatlis[i1] , Quatlis[i2]
    t1 , t2 = Tlis[i1] , Tlis[i2]
    tt = (t - t1) / (t2 - t1)
    qq = cgt.slerp(tt,q1,q2)
    return qq

def interp_quat_multi(Tlis , Quatlis , tlis):
    outlis = []
    for t in tlis:
        q = interp_quat(Tlis , Quatlis , t)
        outlis.append(q)
    return np.array(outlis)
