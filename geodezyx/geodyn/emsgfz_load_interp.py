#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 16:43:48 2021

@author: psakicki
"""

import pandas as pd
import numpy as np
import netCDF4 as nc
import itertools
import scipy

from geodezyx import utils



def EMSGFZ_extrapolator(path_or_netcdf_object_in,
                        time_xtrp,
                        lat_xtrp,
                        lon_xtrp,
                        wished_values=("duV","duNS","duEW","dg"),
                        output_type = "DataFrame",
                        debug=False):
    """
    Extrapolate loading values from the EMSGFZ models
    esmdata.gfz-potsdam.de:8080/

    Parameters
    ----------
    path_or_netcdf_object_in : string or NetCDF object
        Input 
        can be a file path (string) or direcly the NetCDF object (faster).
    time_xtrp : float or float iterable
        time for the wished extrapolated values
        for daily files: hours of day [0..23].
        for yearly files: day of years [0..364].
    lat_xtrp : float or float iterable
        latitude component for the wished extrapolated values
        ranging from [-90..90]
    lon_xtrp : float or float iterable
        longitude component for the wished extrapolated values.
        ranging from [-180..180]
    wished_values : tuple of string, optional
        the components of the extrapolated values. 
        The default is ("duV","duNS","duEW","dg").
    output_type : str, optional
        Choose the output type.
        "DataFrame","dict","array","tuple","list"
        The default is "DataFrame".
    debug : bool, optional
        returns the NetCDF object for debug purposes

    Returns
    -------
    Points_out : see output_type
        The extrapolated values.

    """
    
    if not utils.is_iterable(time_xtrp):
        time_xtrp = [time_xtrp]
    if not utils.is_iterable(lat_xtrp):
        lat_xtrp = [lat_xtrp]
    if not utils.is_iterable(lon_xtrp):
        lon_xtrp = [lon_xtrp]
    
    Points_xtrp = (time_xtrp,lat_xtrp,lon_xtrp)

    if type(path_or_netcdf_object_in) is str:
        NC =  nc.Dataset(path_or_netcdf_object_in)
    else:
        NC = path_or_netcdf_object_in
        
    if debug:
        return NC
    
    time = np.array(NC['time'])
    lat  = np.flip(np.array(NC['lat']))
    lon  = np.array(NC['lon'])
    
    Points = (time,lat,lon)    
    
    WishVals_Stk = []
    WishVals_dic = dict()
    
    #### prepare dedicated time, lat, lon columns
    Points_xtrp_array = np.array(list(itertools.product(*Points_xtrp)))
    WishVals_dic['time'] = Points_xtrp_array[:,0]
    WishVals_dic['lat']  = Points_xtrp_array[:,1]
    WishVals_dic['lon']  = Points_xtrp_array[:,2]
    
    
    ### do the interpolation for the wished value
    for wishval in wished_values:
        Val = np.array(NC[wishval])
        Val = np.flip(Val,1)
        Val_xtrp = scipy.interpolate.interpn(Points,Val,
                                             Points_xtrp)
        WishVals_Stk.append(Val_xtrp)
        WishVals_dic[wishval] = Val_xtrp
    
    #### choose the output
    if output_type == "DataFrame":
        Points_out = pd.DataFrame(WishVals_dic)
    elif output_type == "dict":
        Points_out = WishVals_dic
    elif output_type == "array":
        Points_out = np.column_stack(WishVals_Stk)
    elif output_type == "tuple":
        Points_out = tuple(WishVals_Stk)
    elif output_type == "list":
        Points_out = list(WishVals_Stk)
    else:
        Points_out = WishVals_Stk
        
    return Points_out
    