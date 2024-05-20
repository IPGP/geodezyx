#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.marine contains functions to interpolate 
GEBCO bathymetry.

it can be imported directly with:
from geodezyx import marine

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

#### Import star style
# from geodezyx import *                   # Import the GeodeZYX modules
# from geodezyx.externlib import *         # Import the external modules
# from geodezyx.megalib.megalib import *   # Import the legacy modules names
log = logging.getLogger(__name__)

##########  END IMPORT  ##########

def gebco_bathy_grid_extractor(dataset,latmin,latmax,lonmin,lonmax):
    """
    for safety reasons, lat and lon input MUST BE in the dataset,
    replaced by the closest elsewhere

    return latnew , lonnew , Znew
    """
    
    lon = dataset['lon'][:]
    lat = dataset['lat'][:]
    Z   = dataset['elevation'][:]

    latmin_dec = latmin
    latmax_dec = latmax
    lonmin_dec = lonmin
    lonmax_dec = lonmax

    #latmin_dec = conv.dms2dec_num(*latmin)
    #latmax_dec = conv.dms2dec_num(*latmax)
    #lonmin_dec = conv.dms2dec_num(*lonmin)
    #lonmax_dec = conv.dms2dec_num(*lonmax)


    if not np.any(latmin_dec == lat):
        log.warning("... replacing to the nearest")
        latmin_dec = utils.find_nearest(lat,latmin_dec)[0]

    if not np.any(latmax_dec == lat):
        log.warning("... replacing to the nearest")
        latmax_dec = utils.find_nearest(lat,latmax_dec)[0]

    if not np.any(lonmin_dec == lon):
        log.warning("... replacing to the nearest")
        lonmin_dec = utils.find_nearest(lon,lonmin_dec)[0]

    if not np.any(lonmax_dec == lon):
        log.warning("... replacing to the nearest")
        lonmax_dec = utils.find_nearest(lon,lonmax_dec)[0]

    boollat = np.logical_and(latmin_dec <= lat, lat <= latmax_dec)
    boollon = np.logical_and(lonmin_dec <= lon, lon <= lonmax_dec)


    gridlat = np.tile( boollat[:,np.newaxis]  , (1,len(lon)))
    gridlon = np.tile( boollon , (len(lat),1))

    Znew   = Z[ gridlon * gridlat ]
    latnew = lat[boollat]
    lonnew = lon[boollon]
    Znew   = Znew.reshape(len(latnew),len(lonnew))

    return latnew , lonnew , Znew
