#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 16:43:48 2021
@author: psakic

This sub-module of geodezyx.geodyn contains functions to load 
and use the GFZ's ESM loading model.

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""

########## BEGIN IMPORT ##########
#### External modules
import datetime as dt
import itertools
#### Import the logger
import logging
import os
import urllib.request
from urllib.parse import quote

import numpy as np
import pandas as pd
import scipy

#### geodeZYX modules
from geodezyx import conv
from geodezyx import utils

log = logging.getLogger('geodezyx')

##########  END IMPORT  ##########

def ESMGFZ_downloader(latitude,longitude,output_dir,
                      components=["NTAL","NTOL","HYDL","SLEL"],
                      CM_CF="CF",
                      outputformat = "csv",
                      formatvariables = "duNS,duEW,duV",
                      startdate = dt.datetime(2019,1,1),
                      enddate   = dt.datetime(2020,1,1),
                      outfile_prefix=""):
    """
    Download loading contribution values for a specific point from ESM's website
    
    http://rz-vm115.gfz-potsdam.de:8080/repository/entry/show?entryid=2827909c-6c9d-46fd-ba2c-f806bf215170&output=wiki.view

    Parameters
    ----------
    latitude : float
        latitude in degrees.
    longitude : float
        longitude in degrees.
    output_dir : str
        ouput directory for the downloaded files.
    components : list of str, optional
        list of the wished loading contributions. 
        The default is ["NTAL","NTOL","HYDL","SLEL"].
    CM_CF : str, optional
        Center of Figure (CF) or Center of Mass (CM). The default is "CF".
    outputformat : str, optional
        choose the output format. CSV is strongly recommended.
        'netcdf' and 'timeseries' is also supported
        The default is "csv".
    formatvariables : str, optional
        outputed variables. The default is "duNS,duEW,duV".
    startdate : datetime, optional
        start date. The default is dt.datetime(2010,1,1).
    enddate : datetime, optional
        start end. The default is dt.datetime(2020,1,1).
    outfile_prefix : str, optional
        A custo. The default is "".

    Returns
    -------
    output_path : TYPE
        DESCRIPTION.

    """
    
    url_base = "http://esmdata.gfz-potsdam.de:8080/"+"repository/entry/show/Home/Elastic+Surface+Loading/"
    
    ### base entry
    if startdate >= dt.datetime(2010,1,1) and enddate >= dt.datetime(2010,1,1):
        url_period = '/2010-now'
    elif startdate < dt.datetime(2010,1,1) and enddate < dt.datetime(2010,1,1):
        url_period = '/2000-2009'
    else:
        url_period = ''
        
        
    for icomp in components:
        url_entry = icomp +"/" + CM_CF + "/2000-now" + url_period + "?submit=Get%20Point&output=data.gridaspoint"
        
        ### point
        url_lat = "&location.latitude=" + str(latitude)
        url_lon = "&location.longitude=" + str(longitude)
        url_format = "&format=" + outputformat
        
        ### period 
        startdate_str = quote(conv.dt2str(startdate) + " UTC")
        enddate_str   = quote(conv.dt2str(enddate)   + " UTC")
        url_cal = "&calendar=proleptic_gregorian" + "&fromdate=" + startdate_str + "&todate=" + enddate_str
        
        ### variables
        url_var = "&variable=" + formatvariables
        
        ##concatenation
        Url_tuple = (url_base,url_entry,url_lat,url_lon,url_format,url_cal,url_var)
        url_out = "".join(Url_tuple)
        
        log.info("INFO: downloaded URL")
        log.info(url_out)
        
        output_file = "_".join((outfile_prefix,
                                icomp,
                                CM_CF,
                                str(latitude),
                                str(longitude),
                                conv.dt2str(startdate,"%Y%m%d_%H%M%S"),
                                conv.dt2str(enddate,"%Y%m%d_%H%M%S")))
        
        output_path = os.path.join(output_dir,output_file)
        
        log.info("INFO: output file")
        log.info(output_path)
            
        DownObjt = urllib.request.urlretrieve(url_out, output_path)
        
    return output_path
        
    
def ESMGFZ_DataFrame_reader(path_in,contribution='auto'):
    """
    Reader for CSV files downloaded from EMS's website
    
    Parameters
    ----------
    path_in : str
        path of the input CSV file.
    contribution : str, optional or None
        add a column for the loading contribution name.
        'auto' will determine the contribution based on the filename
        None will not generate such column
        '<contribution name>': manual name
        The default is 'auto'.

    Returns
    -------
    DF : DataFrame
        output DataFrame.

    """
    
    DF = pd.read_csv(path_in,sep=",") #,delim_whitespace=False)
    DF.time = pd.to_datetime(DF.time)
    DF.columns = [e.split("[")[0] for e in  DF.columns]
    
    if contribution:
        if contribution == 'auto':
            name = os.path.basename(path_in)
            DF["contrib"] = 'UNKW'
            for icontrib in ["NTAL","NTOL","HYDL","SLEL"]:
                if icontrib in name:
                    DF["contrib"] = icontrib
                    
        else:
            DF["contrib"] = contribution
            
    return DF


def ESMGFZ_extrapolator(path_or_netcdf_object_in,
                        time_xtrp,
                        lat_xtrp,
                        lon_xtrp,
                        wished_values=("duV","duNS","duEW"),
                        output_type = "DataFrame",
                        debug=False,verbose=True,
                        time_smart=True,
                        interp_method="splinef2d"):
    """
    Extrapolate loading values from the EMSGFZ models
    esmdata.gfz-potsdam.de:8080/

    Parameters
    ----------
    path_or_netcdf_object_in : string, list of strings or NetCDF object
        Input 
        can be a file path (string), a list of string (will be concatenated)
        or direcly the NetCDF object (faster).
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
        The default is ("duV","duNS","duEW").
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
    
    #internal import because this module is complex to install
    import netCDF4 as nc

    
    if not utils.is_iterable(time_xtrp):
        time_xtrp = [time_xtrp]
    if not utils.is_iterable(lat_xtrp):
        lat_xtrp = [lat_xtrp]
    if not utils.is_iterable(lon_xtrp):
        lon_xtrp = [lon_xtrp]
    
    Points_xtrp = (time_xtrp,lat_xtrp,lon_xtrp)

    if type(path_or_netcdf_object_in) is str:
        NC =  nc.Dataset(path_or_netcdf_object_in)
    elif type(path_or_netcdf_object_in) in (nc.Dataset,nc.MFDataset):
        NC = path_or_netcdf_object_in
    else:
        NC = nc.MFDataset(sorted(path_or_netcdf_object_in))
        
    if debug:
        return NC
    
    time_orig = np.array(NC['time'][:])
    
    if time_smart:
        # we work in MJD
        start_date = conv.dt2mjd(conv.str_date2dt(NC['time'].units[11:]))
        if len(time_orig) <= 366:
            time = start_date - .0 + time_orig 
        else:
            time = start_date - .0 + np.array(range(len(time_orig)))
        
    lat  = np.flip(np.array(NC['lat'][:]))  ### we flip the lat because not ascending !
    lon  = np.array(NC['lon'][:])
    
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
        if verbose:
            log.info("INFO: %s start interpolation at %s",wishval,dt.datetime.now())
        
        #### Val = np.array(NC[wishval]) ### Slow
        Val = NC[wishval][:]
        
        if verbose:
            log.info("INFO: %s grid loaded at %s",wishval,dt.datetime.now())
            
        Val = np.flip(Val,1) ### we flip the lat because not ascending !
        Val_xtrp = scipy.interpolate.interpn(Points,Val,
                                             Points_xtrp,
                                             method=interp_method)
        WishVals_Stk.append(Val_xtrp)
        WishVals_dic[wishval] = Val_xtrp
    
    #### choose the output
    if output_type == "DataFrame":
        Points_out = pd.DataFrame(WishVals_dic)
        if time_smart:
            Points_out['time_dt'] = conv.mjd2dt(Points_out['time'])
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
    