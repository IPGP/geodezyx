#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 15:11:12 2020

@author: psakicki

This sub-module of geodezyx.conv deals with time conversion.

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
import datetime as dt
import math
import numpy as np
import os 
import pandas as pd
import scipy
from scipy.spatial.transform import Rotation
import string
import struct
import subprocess
import re
import time
import warnings

#### geodeZYX modules
from geodezyx import utils,stats

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names

##########  END IMPORT  ##########



#import utils.utils as utils

#  _______ _                   _____                              _
# |__   __(_)                 / ____|                            (_)
#    | |   _ _ __ ___   ___  | |     ___  _ ____   _____ _ __ ___ _  ___  _ __
#    | |  | | '_ ` _ \ / _ \ | |    / _ \| '_ \ \ / / _ \ '__/ __| |/ _ \| '_ \
#    | |  | | | | | | |  __/ | |___| (_) | | | \ V /  __/ |  \__ \ | (_) | | | |
#    |_|  |_|_| |_| |_|\___|  \_____\___/|_| |_|\_/ \___|_|  |___/_|\___/|_| |_|
#
### Time conversion


def rinex_regex(compressed=True,compiled=False):
    """
    Return a regex corresponding to a RINEX name (old convention)

    Parameters
    ----------
    compressed : bool
        return a the regex for a compressed rinex

    compiled : bool
        return a Python regex object already compliled
        
    Returns
    -------
    out : string or python's regex
        a regex
    """
    if compressed:
        regexstr = "....[0-9]{3}.\.[0-9]{2}((d\.(Z)|(gz))|(o)|(d))"
    else:
        regexstr = "....[0-9]{3}.\.[0-9]{2}o"

    if compiled:
        return re.compile(regexstr)
    else:
        return regexstr


def rinex_regex_new_name(compressed=None,compiled=False):
    """
    Return a regex corresponding to a RINEX name (new convention)

    Parameters
    ----------
    compressed : bool or None
        return a the regex for a compressed rinex
        if None, does not matter (return both compressed or not)
        

    compiled : bool
        return a Python regex object already compliled
        
    Returns
    -------
    out : string or python's regex
        a regex
    """
    if compressed:
        regexstr = ".{4}[0-9]{2}.{3}_(R|S|U)_[0-9]{11}_[0-9]{2}\w_[0-9]{2}\w_\w{2}\.\w{3}\.gz"
    elif compressed is None:
        regexstr = ".{4}[0-9]{2}.{3}_(R|S|U)_[0-9]{11}_[0-9]{2}\w_[0-9]{2}\w_\w{2}\.\w{3}(\.gz)?"        
    else:
        regexstr = ".{4}[0-9]{2}.{3}_(R|S|U)_[0-9]{11}_[0-9]{2}\w_[0-9]{2}\w_\w{2}\.\w{3}"

    if compiled:
        return re.compile(regexstr)
    else:
        return regexstr



def rinex_regex_new_name_brdc(compressed=None,compiled=False):
    """
    Return a regex corresponding to a BROADCAST RINEX name (new convention)

    Parameters
    ----------
    compressed : bool or None
        return a the regex for a compressed rinex
        if None, does not matter (return both compressed or not)
        

    compiled : bool
        return a Python regex object already compliled
        
    Returns
    -------
    out : string or python's regex
        a regex
    """
    if compressed:
        regexstr = ".{4}[0-9]{2}.{3}_(R|S|U)_[0-9]{11}_[0-9]{2}\w_\w{2}\.\w{3}\.gz"
    elif compressed is None:
        regexstr = ".{4}[0-9]{2}.{3}_(R|S|U)_[0-9]{11}_[0-9]{2}\w_\w{2}\.\w{3}(\.gz)?"        
    else:
        regexstr = ".{4}[0-9]{2}.{3}_(R|S|U)_[0-9]{11}_[0-9]{2}\w_\w{2}\.\w{3}"

    if compiled:
        return re.compile(regexstr)
    else:
        return regexstr    
    
    
def rinex_regex_new_name_brdc_gfz_godc(compressed=None,compiled=False):
    """
    Return a regex corresponding to a BROADCAST RINEX name (new convention)
    from the GFZ's GODC archive

    Parameters
    ----------
    compressed : bool or None
        return a the regex for a compressed rinex
        if None, does not matter (return both compressed or not)
        

    compiled : bool
        return a Python regex object already compliled
        
    Returns
    -------
    out : string or python's regex
        a regex
    """
    ### BRDC00GFZ_00000000_FRO_RX3_MN_20210114_000000_01D_00U_GFZ.rnx
    if compressed:
        regexstr = ".{4}[0-9]{2}.{3}_[0-9]{8}_.{3}_.{3}_.{2}_[0-9]{8}_[0-9]{6}_[0-9]{2}\w_[0-9]{2}\w_.{3}.\w{3}\.gz"
    elif compressed is None:
        regexstr = ".{4}[0-9]{2}.{3}_[0-9]{8}_.{3}_.{3}_.{2}_[0-9]{8}_[0-9]{6}_[0-9]{2}\w_[0-9]{2}\w_.{3}.\w{3}(\.gz)?"        
    else:
        regexstr = ".{4}[0-9]{2}.{3}_[0-9]{8}_.{3}_.{3}_.{2}_[0-9]{8}_[0-9]{6}_[0-9]{2}\w_[0-9]{2}\w_.{3}.\w{3}"

    if compiled:
        return re.compile(regexstr)
    else:
        return regexstr    
    
    
    
    

def tgipsy2dt(tin):
    """
    Time representation conversion
    
    GIPSY Time => Datetime

    Parameters
    ----------
    tin : float or list/numpy.array of float
        GIPSY time(s). Can handle several time float in a list.
                
    Returns
    -------
    dtout : datetime or list/numpy.array of datetime
        Converted Datetime(s)
        
    Note
    ----
    GIPSY time is not the 'real' J2000
    but only counted starting from 1st January 2000 at Noon
    """
    if utils.is_iterable(tin):
        typ=utils.get_type_smart(tin)
        return typ([tgipsy2dt(e) for e in tin])
    else:
        j2000 = dt.datetime(2000, 1, 1, 12, 0)
        tout = j2000 + dt.timedelta(seconds=float(tin))
        return tout

    return tout


def matlab_time2dt(matlab_datenum):
    """
    Time representation conversion
    
    MATLAB Time => Datetime

    Parameters
    ----------
    matlab_datenum : float or list/numpy.array of float
        MATLAB time(s).  Can handle several time floats in a list.
                
    Returns
    -------
    python_datetime : datetime or list/numpy.array of datetime
        Converted Datetime(s)
    """
    if utils.is_iterable(matlab_datenum):
        typ=utils.get_type_smart(matlab_datenum)
        return typ([matlab_time2dt(e) for e in matlab_datenum])
    else:
        python_datetime = dt.datetime.fromordinal(int(matlab_datenum)) + \
        dt.timedelta(days=matlab_datenum%1) - dt.timedelta(days = 366)
    return python_datetime


def round_dt(dtin,round_to,python_dt_out=True):
    """
    Round a datetime object to any time laps in seconds
    
    Parameters
    ----------
    dtin : datetime or list/numpy.array of datetime 
        Datetime you want to round
        Can handle several datetimes in an iterable.
            
    round_to : str
        The way to round the datetime. 
        
        It follows the Pandas' Series conventions, e.g.:
            
        * one-day rounding: '1D'
        * one-minute rounding: '1min'
        * one-second rounding: '1s'
        
        Full list is in the Note's link
        
    python_dt_out : bool, optional
        If True, it returns the date as a legacy Python's DateTime.
        If False, it returns a Pandas Timestamp.
        The default is True.
            
    Returns
    -------  
    dtout : datetime or list/numpy.array of datetime
        Rounded Datetime
        
    Note
    ----
    https://pandas.pydata.org/docs/user_guide/timeseries.html#timeseries-offset-aliases
    """

    ### Here we adopt a new scheme for the recursive approach
    ### because we switch per defaut to a Pandas Series (PSakic 2021-02-22)
    if not utils.is_iterable(dtin):
        singleton = True
        dtin_use = pd.Series([dtin])
    else:
        singleton = False
        typ=utils.get_type_smart(dtin)
        dtin_use = pd.Series(dtin)

    dtin_out = dtin_use.dt.round(round_to)
    
    if python_dt_out:
        dtin_out.dt.to_pydatetime()

    if singleton:
        return dtin_out.iloc[0]
    else:
        return typ(dtin_out)
        
##### Nota Bene
##### dt_round & roundTime moved to the function graveyard (PSakic 2021-02-22)


def dt_ceil(dtin):
    """
    Round a datetime object to the begining of the day
    
    Parameters
    ----------
    dtin : datetime or list/numpy.array of datetime 
        Datetime you want to round, default now.
        Can handle several datetimes in an iterable.
            
    Returns
    -------  
    dtout : datetime or list/numpy.array of datetime
        Rounded Datetime
        
    """
    if utils.is_iterable(dtin):
        typ=utils.get_type_smart(dtin)
        return typ([dt_ceil(e) for e in dtin])
    else:
        return date2dt(dtin.date())
   

def dt_in_local_timezone2posix(dtin):
    if not utils.is_iterable(dtin):
        return time.mktime(dtin.timetuple()) + dtin.microsecond * 0.000001
    else:
        typ=utils.get_type_smart(dtin)
        return typ([ dt_in_local_timezone2posix(e) for e in dtin ])


def posix2dt_in_local_timezone(posixin):
    if not utils.is_iterable(posixin):
        if np.isnan(posixin):
            return dt.datetime(1970,1,1)
        else:
            return dt.datetime.fromtimestamp(posixin)
    else:
        typ=utils.get_type_smart(posixin)
        return typ([ posix2dt_in_local_timezone(e) for e in posixin ])

def dt_range(start_dt,end_dt,
             day_step=1,sec_step=0):
    """
    Range of datetime between a start and end (included)
    
    Parameters
    ----------
    start_dt,end_dt : datetime
        Datetimes
        
    Returns
    -------
    Out_range : list of datetime
        range of dates  
    """
    Out_range = [start_dt]
    while Out_range[-1] < end_dt:
        Out_range.append(Out_range[-1] + dt.timedelta(days=day_step) + dt.timedelta(seconds=sec_step))
    return Out_range

def dt2posix(dtin,out_array=False):       
    """
    Time representation conversion
    
    Python's Datetime => POSIX Time

    Parameters
    ----------
    dtin : datetime or list/numpy.array of datetime.
        Datetime(s).  Can handle several datetimes in an iterable.
    out_array : bool
        if iterable as input, force the output as a Numpy array

    Returns
    -------
    L : float or list/numpy.array of floats
        POSIX Time(s)  
    """
    if not utils.is_iterable(dtin):
        D = dtin - dt.datetime(1970,1,1)
        return D.days * 86400 + D.seconds +  D.microseconds * 10**-6
    else:
        typ=utils.get_type_smart(dtin)
        L = typ([ dt2posix(e) for e in dtin  ])
        if out_array:
            return np.array(L)
        else:
            return L

def posix2dt(posixin,out_array=False):       
    """
    Time representation conversion
    
    POSIX Time => Python's Datetime

    Parameters
    ----------
    posixin : float or list/numpy.array of floats.
        POSIX Time.  Can handle several time in a list.
    out_array : bool
        if iterable as input, force the output as a Numpy array

    Returns
    -------
    L : datetime or list/numpy.array of datetime.
        Converted Datetime(s)
    """

    if not utils.is_iterable(posixin):
        if np.isnan(posixin):
            return dt.datetime(1970,1,1)
        else:
            return dt.datetime(1970,1,1) + dt.timedelta(seconds=posixin)
    else:
        typ=utils.get_type_smart(posixin)
        L = typ([ posix2dt(e) for e in posixin ])
        if out_array:
            return np.array(L)
        else:
            return L


def ymdhms2dt(y=0,mo=0,d=0,h=0,mi=0,s=0,ms=0):
    """
    Improved time representation conversion to Datetime

    can manage list, array, or floats as input
    can manage NaN too (return POSIX 0 epoch if NaN)

    Parameters
    ----------
    y,mo,d,h,mi,s,ms : float or list/numpy.array of floats.
        year , month... . Can handle several times in lists.

    Returns
    -------
    Datetime
        Converted Datetime(s)
    """
    if not utils.is_iterable(y):
        y = int(y)
        mo = int(mo)
        d = int(d)
        h = int(h)
        mi = int(mi)
        s = float(s)
        ms = float(ms)
        
        try:
            ms_from_s  = (s - np.floor(s)) * 10**6
            if ms_from_s != 0:
                if ms != 0:
                    print('WARN : input sec contains microsecs, given microsec are ignored')
                ms = ms_from_s
            return dt.datetime(int(y),int(mo),int(d),int(h),int(mi),int(s),int(ms))
        except:
            return dt.datetime(1970,1,1) # si ca deconne, si on donne un NaN par ex
    else:
        typ=utils.get_type_smart(y)
        if ms == 0:
            ms = [0] * len(y)
        return typ([ ymdhms2dt(*e) for e in zip(y,mo,d,h,mi,s,ms) ])
        
                
    
    

def datetime_improved(*args):
    return ymdhms2dt(*args)
    
dt_improved = datetime_improved
    

def dt2ymdhms(dtin,with_microsec = True):
    """
    Time representation conversion

    Python's Datetime => Lists of Years, Months, Days, Hours, Minutes, Seconds

    Parameters
    ----------
        dtin : datetime or list/numpy.array of datetime
            Python DateTime

        with_microsec : bool 
            if False : rounding the microsecs to the nearest sec

    Returns
    -------
        tuple with year , month , day , hour , minute , second , microsecond
    """
    if not utils.is_iterable(dtin):
        if with_microsec:
            return (dtin.year,dtin.month,dtin.day,dtin.hour,dtin.minute,dtin.second,dtin.microsecond)
        else:
            dt2 = roundTime(dtin,1)
            return (dt2.year,dt2.month,dt2.day,dt2.hour,dt2.minute,dt2.second)
    else:
        typ=utils.get_type_smart(dtin)
        return typ([ dt2ymdhms(e,with_microsec) for e in dtin ])

def ymdhms_vectors2dt(yrlis,mlis,dlis,hlis,minlis,slis):
    """
    Time representation conversion
    
    Lists of Years, Months, Days, Hours, Minutes, Seconds => Python's Datetime
    
    Parameters
    ----------
    yrlis,mlis,dlis,hlis,minlis,slis : float or list/numpy.array of floats.
        Lists of Years, Months, Days, Hours, Minutes, Seconds
        
    Returns
    -------
    L : datetime or list/numpy.array of datetime.
        Datetime(s)
    """
    dtlis = []
    for yr,m,d,minn,h,s in zip(yrlis,mlis,dlis,hlis,minlis,slis):
        dtlis.append(dt.datetime(int(yr),int(m),int(d),int(h),int(minn),int(s)))
    return np.array(dtlis)

def doy2dt(year,days,hours=0,minutes=0,seconds=0):
    """
    Time representation conversion
    
    Day of Year Time => Python's datetime

    Parameters
    ----------
    year, days : float or list/numpy.array of floats.
        year, days of year
    hours, minutes, seconds : float or list/numpy.array of floats, optional
        hours, minutes, seconds

    Returns
    -------
    L : datetime or list/numpy.array of datetime.
        Datetime(s)
    """
    if not utils.is_iterable(year):
        # All this because Python cant handle int with a starting with 0 (like 08)
        # => SyntaxError: invalid token
        year    = int(str(year))
        days    = int(str(days))
        hours   = float(str(hours))
        minutes = float(str(minutes))
        seconds = float(str(seconds))

        tempsecs = seconds + 60 * minutes + 3600 * hours
        #finalsecs     = np.floor(tempsecs)
        finalmicrosec = np.round(tempsecs * 10**6)

        return dt.datetime(year, 1, 1) + dt.timedelta(days - 1) + \
        dt.timedelta(microseconds=finalmicrosec)

    else:
        if not utils.is_iterable(hours):
            hours   = [0] * len(year)
        if not utils.is_iterable(minutes):
            minutes = [0] * len(year)
        if  not utils.is_iterable(seconds):
            seconds = [0] * len(year)

        outlis = []
        typ=utils.get_type_smart(year)
        for y,d,h,m,s in zip(year,days,hours,minutes,seconds):
            outlis.append(doy2dt(y,d,h,m,s))
        return typ(outlis)

def dt2doy(dtin,outputtype=str):  
    return outputtype(dtin.strftime('%j'))

def dt2doy_year(dtin,outputtype=str):
    """
    Time representation conversion
    
    Python's datetime => Day of Year, Year

    Parameters
    ----------
    dtin : datetime or list/numpy.array of datetime
        Datetime(s). Can handle several datetimes in an iterable.
        
    Returns
    -------
    doy , year : tuple of int
        Day of Year and Year
        Is a list of tuple if the input is an iterable
    """
    
    if utils.is_iterable(dtin):
        typ=utils.get_type_smart(dtin)
        return typ([dt2doy_year(e,outputtype=outputtype) for e in dtin])
        
    else:
        return outputtype(dtin.strftime('%j')),outputtype(dtin.strftime('%Y'))
    
def dt2fracday(dtin):
    """
    Time representation conversion

    Python's datetime => Fraction of the day

    Parameters
    ----------
    dtin : datetime or list/numpy.array of datetime
        Datetime(s). Can handle several datetimes in an iterable.
        
    Returns
    -------
    fraction_day : float
        Fractional of the day
        Is a list of int if the input is an iterable
    """
    if utils.is_iterable(dtin):
        typ=utils.get_type_smart(dtin)
        return typ([dt2fracday(e) for e in dtin])
    else:
        return dtin.hour / 24 + dtin.minute / (60*24) + dtin.second / (60*60*24)
    
def dt2secinday(dtin):
    """
    Time representation conversion
    
    Python's datetime => Seconds in days

    Parameters
    ----------
    dtin : datetime or list/numpy.array of datetime
        Datetime(s). Can handle several datetimes in an iterable.
        
    Returns
    -------
    secinday : int
        Seconds in the day
        Is a list of int if the input is an iterable
    """
    if utils.is_iterable(dtin):
        typ=utils.get_type_smart(dtin)
        return typ([dt2secinday(e) for e in dtin])
    else:
        return dtin.hour * 3600 + dtin.minute * 60 + dtin.second

def dt2tuple(dtin):
    return tuple(dtin.timetuple())[:-3]



def tup_or_lis2dt(lisin):
    """
    Time representation conversion
    
    Date-looking strings => Python's datetime

    Parameters
    ----------
    lisin : iterable (tuple/list/numpy.array) of string.
        list of Date-looking strings
        like : ["2018","12","31","12","30","00"]
    Returns
    -------
    L : list of datetime
        converted DateTimes
    """
    
    try:
        return dt.datetime(*[int(float(e)) for e in lisin])
    except:
        return posix2dt(0)


def dt_utc2dt_tai(dtin):
    """
    Time scale conversion
    
    Python's datetime in UTC => Python's datetime in TAI
    
    (Wrapper to correct the leap second)

    Parameters
    ----------
    dtin : datetime or list/numpy.array of datetime
        Datetime(s). Can handle several datetimes in an iterable.
        
    Returns
    -------
    L : list of datetime
        converted DateTimes
        
    Note
    ----
    TAI is (currently) ahead of UTC

    """
    
    if utils.is_iterable(dtin):
        typ=utils.get_type_smart(dtin)
        return typ([dt_utc2dt_tai(e) for e in dtin])
    else:
        return dtin + dt.timedelta(seconds=find_leapsecond(dtin)) 

    
def dt_tai2dt_utc(dtin):
    """
    Time scale conversion
    
    Python's datetime in TAI => Python's datetime in UTC
    
    (Wrapper to correct the leap second)

    Parameters
    ----------
    dtin : datetime or list/numpy.array of datetime
        Datetime(s). Can handle several datetimes in an iterable.
        
    Returns
    -------
    L : list of datetime
        converted DateTimes
        
    Note
    ----
    TAI is (currently) ahead of UTC
    
    """
    
    if utils.is_iterable(dtin):
        typ=utils.get_type_smart(dtin)
        return typ([dt_tai2dt_utc(e) for e in dtin])
    else:
        return dtin + dt.timedelta(seconds=-find_leapsecond(dtin))
    
     
    
def dt_tai2dt_tt(dtin):
    """
    Time scale conversion
    
    Python's datetime in TAI => Python's datetime in Terrestrial Time (TT)
    
    Parameters
    ----------
    dtin : datetime or list/numpy.array of datetime
        Datetime(s). Can handle several datetimes in an iterable.
        
    Returns
    -------
    L : list of datetime
        converted DateTimes
        
    Source
    ------
    https://en.wikipedia.org/wiki/Terrestrial_Time
    
    Note
    ----
    Realization of the TT by the BIPM 
    (restimation and correction by dozen of µsec)
    https://www.bipm.org/en/bipm-services/timescales/time-ftp/ttbipm.html
        
    """
    
    if utils.is_iterable(dtin):
        typ=utils.get_type_smart(dtin)
        return typ([dt_tai2dt_tt(e) for e in dtin])
    else:
        return dtin + dt.timedelta(seconds=32.184)


def dt_utc2dt_ut1(dtin,dUT1):
    """
    Time scale conversion
    
    Python's datetime in UTC => Python's datetime in UT1
    
    Parameters
    ----------
    dtin : datetime or list/numpy.array of datetime
        Datetime(s). Can handle several datetimes in an iterable.
    dUT1 : float
        UT1-UTC in seconds.

    Returns
    -------
    L : list of datetime
        converted DateTimes

    """
    if utils.is_iterable(dtin):
        typ=utils.get_type_smart(dtin)
        return typ([dt_utc2dt_ut1(e) for e in dtin])
    else:
        return dtin + dt.timedelta(seconds=dUT1)
    
    

def dt_utc2dt_ut1_smart(dtin,DF_EOP_in,
                        use_interp1d_obj=True,
                        EOP_interpolator=None):
    """
    Time scale conversion
    
    Python's datetime in UTC => Python's datetime in UT1
    using an EOP DataFrame provided by files_rw.read_eop_C04
    
    Parameters
    ----------
    dtin : datetime or list/numpy.array of datetime
        Datetime(s). Can handle several datetimes in an iterable.
    DF_EOP_in : DataFrame
        EOP DataFrame for the UT1-UTC
        provided by files_rw.read_eop_C04
    use_interp1d_obj : TYPE, optional
        Use an interp1d_time Interpolator object for the EOP determination 
        at the right epoch.
        Faster in recursive mode when dtin is a list/array
        The default is True.
    EOP_interpolator : interp1d_time object, optional
        The interp1d_time Interpolator object for the EOP determination
        Will be determined automatically inside the function
        The default is None.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """    
    #### internally we work with "epoch" as index
    if DF_EOP_in.index.name != "epoch":
        DF_EOP = DF_EOP_in.set_index("epoch")
    else:
        DF_EOP = DF_EOP_in

    if utils.is_iterable(dtin): #### ITERABLE CASE
        typ=utils.get_type_smart(dtin)
        
        ### if iterable, we optimize DF_EOP to speed up the fct        
        dtmin = np.min(dtin)
        dtmax = np.max(dtin)
        
        BOOL = ((DF_EOP.index > dtmin - dt.timedelta(days=2)) & 
                (DF_EOP.index < dtmax + dt.timedelta(days=2)))
    
        DF_EOP = DF_EOP[BOOL]
        
        ### We also use the interpolator class
        if use_interp1d_obj:
            IEOP = interp1d_time(DF_EOP.index.values,DF_EOP['UT1-UTC'])
        else:
            IEOP = None
        
        return typ([dt_utc2dt_ut1_smart(e,DF_EOP,
                                        use_interp1d_obj=use_interp1d_obj,
                                        EOP_interpolator=IEOP) for e in dtin])
    
    else: #### SINGLE ELEMENT CASE
        if (dtin in DF_EOP.index): ## the EOP value is directly in the EOP DF
            dUT1 = DF_EOP.iloc[DF_EOP.index.get_loc(dtin)]['UT1-UTC']
        elif EOP_interpolator and use_interp1d_obj: ### use the interp class, faster in recursive mode
            dUT1 = EOP_interpolator(dtin)
        else: ## the EOP value is interpolated with a "manual" linear interpo.
            dUT1bef = DF_EOP.iloc[DF_EOP.index.get_loc(dtin,'ffill')]
            dUT1aft = DF_EOP.iloc[DF_EOP.index.get_loc(dtin,'bfill')]

    
            a,b,b2 = stats.linear_coef_a_b(dUT1bef['MJD'],
                                           dUT1bef['UT1-UTC'],
                                           dUT1aft['MJD'],
                                           dUT1aft['UT1-UTC'])
            
            dUT1 = stats.linear_reg_getvalue(dt2MJD(dtin), a, b,
                                             full=False)
            
        return dt_utc2dt_ut1(dtin,dUT1)
        

def dt2gpstime(dtin,dayinweek=True,inp_ref="utc",outputtype=int):
    """
    Time scale conversion
    
    Python's datetime (UTC,TAI or GPS time scale) => GPS time

    Parameters
    ----------
    dtin : datetime or list/numpy.array of datetime
        Datetime(s). Can handle several datetimes in an iterable.
        
    inp_ref : str
        "utc" : apply the 19 sec & leap second correction at the epoch 
        "gps" : no correction applied 
        "tai" : apply -19 sec correction
        
    dayinweek : bool
        if True : returns  GPS week, day in GPS week
        
        if False : returns  GPS week, sec in GPS week
    Returns
    -------
    GPS_week, GPS_day/GPS_sec : tuple of int
        GPS week or GPS day/GPS sec
        
        Is a list of tuple if the input is an iterable
    """
    
    if utils.is_iterable(dtin):
        return [dt2gpstime(e) for e in dtin]
        
    else:
        week_raw , secs_raw = utc2gpstime(dtin.year,dtin.month,dtin.day,dtin.hour,
                                          dtin.minute,dtin.second)
        

        if inp_ref == "utc":
            ### utc : utc2gpstime did the job
            week , secs = week_raw , secs_raw

        elif inp_ref == "tai":
            ### tai : utc2gpstime did the job, but needs leap sec correction again
            utc_offset = find_leapsecond(dtin)
            week , secs = week_raw , secs_raw - utc_offset

        elif inp_ref == "gps":
            ### tai : utc2gpstime did the job, but needs leap sec & 19sec correction again
            utc_offset = find_leapsecond(dtin)
            week , secs = week_raw , secs_raw + 19 - utc_offset
            
        if dayinweek:
            day = np.floor(np.divide(secs,86400))
            output1 = outputtype(int(week)) 
            output2 = outputtype(int(day))
        else:
            output1 = outputtype(int(week))             
            output2 = outputtype(int(secs))


        if outputtype == str: ### add an extra 0 if week < 1000
            output1 = output1.zfill(4)
            
        return output1 , output2



def dt2gpsweek_decimal(dtin,return_middle_of_day=True):
    """
    Time representation conversion
    
    Python's datetime => GPS week
    
    Return the GPS week in a decimal mode
    
    (Easier for plots)

    Parameters
    ----------
    dtin : datetime or list/numpy.array of datetime
        Datetime(s). Can handle several datetimes in an iterable.
        
    return_middle_of_day : bool
        if False, return the decimal part of the week at the begining of the day
    
    Returns
    -------
    GPS_week_decimal : float
        decimal GPS week
        
        Is a list of float if the input is an iterable
    """
    if utils.is_iterable(dtin):
        typ=utils.get_type_smart(dtin)
        return typ([dt2gpsweek_decimal(e) for e in dtin])
    else:
        if return_middle_of_day:
            mod = 0.5
        else:
            mod = 0
        week , day = dt2gpstime(dtin)
        return float(week) + (float(day + mod) / 7 )



def dt2list(dtin,return_useful_values=True):
    """
    Time representation conversion
    
    Python's Datetime => Years, Months, Days, Hours, Minutes, Seconds

    Parameters
    ----------
    dtin : datetime
        Datetime
    return_useful_values : bool
        if True, returns only Years, Months, Days

    Returns
    -------
    L : list of int
        Years, Months, Days, [Hours, Minutes, Seconds]
    """
    if return_useful_values:
        return list(dtin.timetuple())[:-3]
    else:
        return list(dtin.timetuple())



def gpstime2dt(gpsweek,gpsdow_or_seconds,dow_input = True,
               output_time_scale="utc"):
    """
    Time scale & representation conversion
    
    GPS Time => Python's datetime

    Parameters
    ----------
    gpsweek : int
        year, days of year
    gpsdow_or_seconds : int
        Day of Week OR Seconds in Weeks
    dow_input : bool
        select if Day of Week (True) OR Seconds in Weeks (False)
    output_time_scale : str
        gives the wished time scale : "utc", "tai", "gps"
        
    Returns
    -------
    L : datetime or list/numpy.array of datetime.
        DateTime  

    Note
    ----
    Only reliable at the day level for the moment

    the leapsecond is found only after a first calc without leapsecond
    it can be some side effects when the input time is close to
    a leap second jump
    
    Source
    ------
    https://gist.github.com/jeremiahajohnson/eca97484db88bcf6b124
    """
    
    if utils.is_iterable(gpsweek):
        typ=utils.get_type_smart(gpsweek)
        return typ([gpstime2dt(w,ds,dow_input,output_time_scale) for w,ds in zip(gpsweek,gpsdow_or_seconds)])
    else:
        if dow_input:
            gpsseconds = gpsdow_or_seconds * 86400 + 86400*.5 # so we are around noon
        else:
            gpsseconds = gpsdow_or_seconds
    
        ## First gross run
        epoch   = dt.datetime(1980,1,6)
        elapsed = dt.timedelta(days=(int(gpsweek)*7),
                               seconds=(int(gpsseconds)+0))
        
        prelim_time = epoch + elapsed
    
        if output_time_scale == "gps":
            final_time = prelim_time
        else:
            leapsec = find_leapsecond(prelim_time)
            
            if output_time_scale == "utc":
                deltasec = leapsec - 19
            elif output_time_scale == "tai":
                deltasec = 19
                
            ## Second run with leap second
            epoch   = dt.datetime(1980,1,6)
            elapsed = dt.timedelta(days=(int(gpsweek)*7),
                                   seconds=(int(gpsseconds)+deltasec))
        
            final_time = epoch + elapsed
    
        if dow_input:
            final_time = final_time.date()
            final_time = dt.datetime(final_time.year,
                                     final_time.month,
                                     final_time.day)
    
        return final_time


def gpsweek_decimal2dt(gpsweekdec_in):
    """
    Time representation conversion
    
    Decimal GPS Week => Python's datetime

    Parameters
    ----------
    gpsweekdec_in : float or list/numpy.array of float
        GPS week in decimal form.

    Returns
    -------
    L : list of datetime
        converted DateTimes

    """
    if utils.is_iterable(gpsweekdec_in):
        typ=utils.get_type_smart(gpsweekdec_in)
        return typ([gpsweek_decimal2dt(e) for e in gpsweekdec_in])
    else:
        week_floor    = np.floor(gpsweekdec_in)
        week_dec_part = gpsweekdec_in - week_floor 
        
        dt_out = gpstime2dt(week_floor,week_dec_part * 7 *86400,
                            dow_input = False)
        return dt_out
    



def dt_gpstime2dt_utc(dtgpsin,out_array=False):
    """
    Time scale conversion
    
    Datetime in GPS Time Scale => Datetime in UTC Time Scale
    
    Correct both the Leap second and the 19sec difference between GPS Time and UTC

    Parameters
    ----------
    dtin : datetime or list/numpy.array of datetime.
        Datetime(s) in GPS Time Scale.  Can handle several datetimes in an iterable.
    out_array : bool
        if iterable as input, force the output as a Numpy array

    Returns
    -------
    L : datetime or list/numpy.array of datetime.
        Datetime(s) in UTC Time Scale
    """
    # on converti le dt gps en dt utc
    if utils.is_iterable(dtgpsin):
        typ=utils.get_type_smart(dtgpsin)
        Out = typ([dt_gpstime2dt_utc(e) for e in dtgpsin])
        if not out_array:
            return Out
        else:
            return np.array(Out)
    else:
        leapsec = find_leapsecond(dtgpsin)
        dtutc = dtgpsin + dt.timedelta(seconds=19) - dt.timedelta(seconds=leapsec)
        return dtutc
    
    
def dt_gpstime2dt_tai(dtgpsin,out_array=False):
    """
    Time scale conversion
    
    Datetime in GPS Time Scale => Datetime in TAI Time Scale
    
    Correct d the 19sec difference between GPS Time and TAI

    Parameters
    ----------
    dtin : datetime or list/numpy.array of datetime.
        Datetime(s) in GPS Time Scale.  Can handle several datetimes in an iterable.
    out_array : bool
        if iterable as input, force the output as a Numpy array

    Returns
    -------
    L : datetime or list/numpy.array of datetime.
        Datetime(s) in TAI Time Scale
    """
    # on converti le dt gps en dt utc
    if utils.is_iterable(dtgpsin):
        typ=utils.get_type_smart(dtgpsin)
        Out = typ([dt_gpstime2dt_tai(e) for e in dtgpsin])
        if not out_array:
            return Out
        else:
            return np.array(Out)
    else:
        dtutc = dtgpsin + dt.timedelta(seconds=19)
        return dtutc

def year_decimal2dt(yearin):
    """
    Time representation conversion
    
    Decimal Year => Python's Datetime 

    Parameters
    ----------
    yearin : float or list/numpy.array of floats.
        Decimal year(s).  Can handle several floats in an iterable.

    Returns
    -------
    L : datetime or list of datetime.
        Datetime(s)  
    """
    if utils.is_iterable(yearin):
        typ=utils.get_type_smart(yearin)
        return typ([year_decimal2dt(e) for e in yearin])
    else:
        return convert_partial_year(yearin)

def dt2year_decimal(dtin):
    """
    Time representation conversion
    
    Python's Datetime => Decimal Year

    Parameters
    ----------
    dtin : datetime or list/numpy.array of datetime.
        Datetime(s).  Can handle several datetimes in an iterable.

    Returns
    -------
    L : float or list/numpy.array of floats
        Decimal Year(s)  
    """
    if utils.is_iterable(dtin):
        typ=utils.get_type_smart(dtin)
        return typ([dt2year_decimal(e) for e in dtin])
    else:
        return toYearFraction(dtin)

def date_string_2_dt(strin):
    """
    Time representation conversion
    
    String => Python's Datetime
    
    Wrapper of dateutil.parser.parse
    
    The input string should looks like a timestamp
    
    Parameters
    ----------
    strin : string or list/numpy.array of strings.
        string(s). Can handle several datetimes in an iterable.

    Returns
    -------
    L : datetime or list of datetime.
        Datetime(s) 
    """
    
    import dateutil.parser
    if utils.is_iterable(strin):
        typ=utils.get_type_smart(strin)
        return typ([date_string_2_dt(e) for e in strin])
    else:
        return dateutil.parser.parse(strin)


### aliases
string_date2dt = date_string_2_dt
str_date2dt    = date_string_2_dt
strdate2dt     = date_string_2_dt




def jjulCNES2dt(jjulin):
    """
    Time representation & scale conversion
    
    Julian Day CNES => Python's Datetime

    Parameters
    ----------
    jjulin : float/int or list/numpy.array of float/int.
        Julian Day CNES.  Can handle several float/int in an iterable.

    Returns
    -------
    L : datetime or list of datetime.
        Datetime(s)

    Note
    ----
    Julian Day CNES starts at 0h Jan 1, 1950 (JD − 2433282.5)
    
    https://www.aviso.altimetry.fr/fr/donnees/outils/jours-calendaires-ou-jours-juliens.html
    """
    
    if utils.is_iterable(jjulin):
        typ=utils.get_type_smart(jjulin)
        return typ([jjulCNES2dt(e) for e in jjulin])
    else:
        return dt.datetime(1950,0o1,0o1,00,00,00) + dt.timedelta(float(jjulin))

def dt2jjulCNES(dtin,onlydays=True):
    """
    Time representation & scale conversion
    
    Python's Datetime => Julian Day CNES 

    Parameters
    ----------
    dtin : datetime or list/numpy.array of datetime.
        Datetime(s).  Can handle several datetimes in an iterable.
        
    only_days : bool
        if False, return also the seconds in the day

    Returns
    -------
    L : int or list of int.
        Julian Day CNES, [seconds in the Day]

    Note
    ----
    Julian Day CNES starts at 0h Jan 1, 1950 (JD − 2433282.5)
    
    https://www.aviso.altimetry.fr/fr/donnees/outils/jours-calendaires-ou-jours-juliens.html
    """
    epok  = (dtin - dt.datetime(1950,0o1,0o1,00,00,00))
    if onlydays:
        return epok.days
    else:
        epok = epok + dt.timedelta(seconds=19)
        return epok.days , epok.seconds

def MJD2dt(mjd_in,seconds=None,round_to='1s'):
    """
    Time representation conversion
    
    Modified Julian Day  => Python's Datetime

    Parameters
    ----------
    mjd_in : float/int or list/numpy.array of float/int.
        Modified Julian Day.  Can handle several float/int in an iterable.
    seconds : float/int or list/numpy.array of float/int.
        Handle a number
    round_to : str or None
        round the output datetime
        see round_dt function for details

    Returns
    -------
    L : datetime or list of datetime.
        Datetime(s)

    Note
    ----
    Modified Julian Day starts at 0h Nov 17, 1858  (JD − 2400000.5)
    
    https://en.wikipedia.org/wiki/Julian_day
    """
    
    if round_to: #### define a lambda fct if we want to round the output dt
        rnd = lambda x: round_dt(x, round_to)
    else: ### indentity lambda if we don't wanna round
        rnd = lambda x: x
        
    # ITERABLE CASE (warn: this fct is not recursive)
    if utils.is_iterable(mjd_in): 
        typ=utils.get_type_smart(mjd_in)
        if seconds:
            if len(seconds) != len(mjd_in):
                print("ERR : MJD2dt : len(seconds) != len(mjd_in) !!")
                raise Exception
        else:
            seconds = np.zeros(len(mjd_in))

        #return typ(rnd([dt.datetime(1858,11,17)+dt.timedelta(days=m,seconds=sec) for m,sec in zip(mjd_in,seconds)]))

        Outdt_Stk = []
        for m,sec in zip(mjd_in,seconds):
            dtout = dt.datetime(1858,11,17) + dt.timedelta(days=m,seconds=sec)
            Outdt_Stk.append(dtout)

        return typ(rnd(Outdt_Stk))
    
    # NON ITERABLE / FLOAT CASE
    else:
        if not seconds:
            seconds = 0    
        return rnd(dt.datetime(1858,11,17) + dt.timedelta(days=mjd_in,seconds=seconds))

def dt2MJD(dtin):
    """
    Time representation conversion
    
    Python's Datetime => Modified Julian Day

    Parameters
    ----------
    dtin : datetime or list/numpy.array of datetime.
        Datetime(s).  Can handle several datetimes in an iterable.

    Returns
    -------
    L : float or list of floats
        Modified Julian Day(s)  

    Note
    ----
    Modified Julian Day starts at 0h Nov 17, 1858  (JD − 2400000.5)
    
    https://en.wikipedia.org/wiki/Julian_day
    """
    # cf http://en.wikipedia.org/wiki/Julian_day
    try:
        typ=utils.get_type_smart(dtin)
        return typ([dt2MJD(t) for t in dtin])
    except:
        delta = (dtin.replace(tzinfo=None) - dt.datetime(1858,11,17))
        return delta.days + (delta.seconds / 86400.) + (delta.microseconds / 864e8)

def dt2str(dtin , str_format="%Y-%m-%d %H:%M:%S"):
    """
    Time representation conversion
    
    Python's Datetime => String

    Just a wrapper of strftime
    
    Parameters
    ----------
    dtin : datetime or list/numpy.array of datetime.
        Datetime(s).  Can handle several datetimes in an iterable.
        
    Returns
    -------
    L : string or list of strings
        Time as string(s)  
        
    Source
    ------
    https://stackoverflow.com/questions/7999935/python-datetime-to-string-without-microsecond-component
    
    http://www.jacksay.com/tutoriaux/bash-shell/bashshell-utilisation-commande-date.html
    """
    
    if utils.is_iterable(dtin):
        typ=utils.get_type_smart(dtin)
        return typ([dt2str(e) for e in dtin])
    else:
        return dtin.strftime(str_format)
    
def rinexname2dt(rinexpath):
    """
    Time representation conversion
    
    RINEX Name (old or new naming convention)  => Python's Datetime
    
    Extract the date in a RINEX name

    Parameters
    ----------
    rinexpath : string
        RINEX path. The RINEX basename will be extracted automatically.

    Returns
    -------
    dt : datetime
        Datetime
    """

    rinexname = os.path.basename(rinexpath)
    #rinexname = rinexpath
        
    
    ##### NEW rinex name
    if re.search(rinex_regex_new_name(),rinexname) or re.search(rinex_regex_new_name_brdc(),rinexname):
        date_str = rinexname.split("_")[2]
        yyyy = int(date_str[:4])
        doy  = int(date_str[4:7])        
        hh   = int(date_str[7:9])        
        mm   = int(date_str[9:11])       
        dt_out = doy2dt(yyyy,doy) + dt.timedelta(seconds=hh*3600 + mm*60)
        return dt_out 

    ##### NEW rinex name -- GFZ's GODC internal name
    ###### OBS RINEX for GFZ GODC not implemented yet !!!!
    if re.search(rinex_regex_new_name_brdc_gfz_godc(),rinexname): 
        date_str = rinexname.split("_")[5]
        time_str = rinexname.split("_")[6]
        yyyy = int(date_str[:4])
        mo   = int(date_str[4:6])     
        dd   = int(date_str[6:8])     
        
        hh   = int(time_str[0:2])        
        mm   = int(time_str[2:4])       
        ss   = int(time_str[4:6])             

        dt_out = dt.datetime(yyyy,mo,dd,hh,mm,ss)
        
        return dt_out 
    
    ##### OLD rinex name    
    elif re.search(rinex_regex(),rinexname):
        alphabet = list(string.ascii_lowercase)
    

        doy  = int(rinexname[4:7])
        yy   = int(rinexname[9:11])
    
        if yy > 80:
            year = yy + 1900
        else:
            year = yy + 2000
    
        if rinexname[7] in alphabet:
            h = alphabet.index(rinexname[7])
        else:
            h = 0
            
        return dt.datetime(year,1,1) + dt.timedelta(days = doy - 1 , seconds = h * 3600)

    else:
        print('ERR: RINEX name is not well formated:',rinexname)
        raise Exception
    

def sp3name2dt(sp3path):
    """
    Time representation conversion
    
    Orbit SP3 Name (legacy/new naming convention)  => Python's Datetime
    
    Extract the date in a Orbit SP3 name

    Parameters
    ----------
    sp3path : string
        Orbit SP3 path. The basename will be extracted automatically.

    Returns
    -------
    dt : datetime
        Datetime
    """

    sp3name = os.path.basename(sp3path)

    if len(sp3name) < 17: ###### A REGEX WOULD BE MUCH BETTER !!!! (PSakic 2021-06)
        return sp3name_leg_2dt(sp3path)
    else:
        return sp3name_v3_2dt(sp3path)
             

def sp3name_leg_2dt(sp3path):
    """
    Time representation conversion
    
    Orbit SP3 Name (old/legacy naming convention)  => Python's Datetime
    
    Extract the date in a Orbit SP3 name

    Parameters
    ----------
    sp3path : string
        Orbit SP3 path. The basename will be extracted automatically.

    Returns
    -------
    dt : datetime
        Datetime
    """

    sp3name = os.path.basename(sp3path)

    week = int(sp3name[3:7])
    dow  = int(sp3name[7])
    return gpstime2dt(week,dow)

def sp3name_v3_2dt(sp3path):
    """
    Time representation conversion
    
    Orbit SP3 Name (new naming convention)  => Python's Datetime
    
    Extract the date in a Orbit SP3 name

    Parameters
    ----------
    sp3path : string
        Orbit SP3 path. The basename will be extracted automatically.

    Returns
    -------
    dt : datetime
        Datetime
    """
    sp3name = os.path.basename(sp3path)
    datestr = sp3name.split("_")
    datestr = datestr[1]
    
    yyyy = int(datestr[0:4])
    doy  = int(datestr[4:7])
    hh   = int(datestr[7:9])
    mm   = int(datestr[9:11])
    
    dt_out = doy2dt(yyyy,doy) + dt.timedelta(hh * 3600 + mm * 60)
    
    return dt_out



def statname_dt2rinexname(statname,datein,rnxtype='d.Z',
                          session_a_instead_of_daily_session = False):
    
    """
    Time representation conversion
    
    Python's Datetime (and station name) => RINEX name
    
    Create a RINEX name from a station name and date

    Parameters
    ----------
    statname : string
        name of the station
    
    datein : datetime
        date of the wished RINEX name
        
    rnxtype : string
        extension of the RINEX ("d.Z",".o") ...
        
    session_a_instead_of_daily_session : bool
        if True, gives an hourly session name (a,b,c ...)
        instead of a daily session (0)

    Returns
    -------
    rinexname : string
        The name of the RINEX
    """

    alphabet = list(string.ascii_lowercase)

    if datein.hour != 0:
        sess = alphabet[datein.hour]
    else:
        if session_a_instead_of_daily_session:
            sess = alphabet[datein.hour]
        else:
            sess = '0'

    return statname + dt2doy(datein) + sess + '.' + datein.strftime('%y') + rnxtype


def datestr_sinex_2_dt(datestrin):
    """
    Time representation conversion
    
    SINEX time format => Python's Datetime 
    
    SINEX time format : 'YY:DDD:SSSSS' or 'YYYY:DDD:SSSSS'
    
    Parameters
    ----------
    datestrin : string or list/numpy.array of string.
        SINEX time format string(s).  Can handle several string in an iterable.
    
    Returns
    -------
    dtout : datetime or list of datetime.
        Datetime(s)
    """
    
    if utils.is_iterable(datestrin):
        return [datestr_sinex_2_dt(e) for e in datestrin]
    else:
        datestr = datestrin
        #### CASE WHERE THE DATE LOOKS LIKE YYDDD:SSSSS
        if re.search("[0-9]{7}:[0-9]{5}",datestr):
            datestr_list = list(datestr)
            datestr_list.insert(4,":")
            datestr = "".join(datestr_list)    
        #### CASE WHERE THE DATE LOOKS LIKE YYYYDDD:SSSSS
        elif re.search("[0-9]{5}:[0-9]{5}",datestr):
            datestr_list = list(datestr)
            datestr_list.insert(2,":")
            datestr = "".join(datestr_list)
            
        ### this test must be independent  
        if '00:000:00000' in datestr:
            return dt.datetime(1970,1,1)
    
        dateint = [int(e) for e in datestr.split(':')]
        yr  = dateint[0]
        doy = dateint[1]
        sec = dateint[2]
        
        ## case for year with only 2 digits
        if re.match("[0-9]{2}:[0-9]{3}:[0-9]{5}",datestr):            
            if yr > 50:
                year = 1900 + yr
            else:
                year = 2000 + yr
        else:
            year = yr
        
        return doy2dt(year,doy,seconds=sec)



def dt_2_sinex_datestr(dtin,short_yy=True,year_sep=":"):
    """
    Time representation conversion
    
    Python's Datetime => SINEX time format 
        
    Parameters
    ----------
    dtin : datetime
    
    Returns
    -------
    dtout : str
        Date in SINEX format
        
    year_sep : str
        year separator, usually :, but can also be empty
    """
    
    if utils.is_iterable(dtin):
        return [dt_2_sinex_datestr(e) for e in dtin]
    else:
        if not short_yy:
            year = str(int(dtin.year))
            year = year.zfill(4)
        else:
            year = str(int(str(dtin.year)[2:]))
            year = year.zfill(2)
        
        doy = str(dt2doy(dtin))
        sec = str(dt2secinday(dtin)) 
        
        strfmt = "{:}" + year_sep + "{:3}:{:5}"
        strout = strfmt.format(year,doy.zfill(3),sec.zfill(5))

        return strout

def dt_2_sp3_datestr(dtin):
    """
    Time representation conversion
    
    Python's Datetime => SP3 time format 
        
    Parameters
    ----------
    dtin : datetime
    
    Returns
    -------
    dtout : str
        Date in SP3 format
    """
    
    if utils.is_iterable(dtin):
        return [dt_2_sp3_datestr(e) for e in dtin]
    else:
        return utils.join_improved("",*dt2gpstime(dtin))


def dt2sp3_timestamp(dtin,start_with_star=True):
    """
    Time representation conversion
    
    Python's Datetime => SP3 Timestamp
    e.g. 
    
    *  2000  5 28  0  0  0.00000000
    """
    yyyy = dtin.year
    mm   = dtin.month
    dd   = dtin.day
    hh   = dtin.hour
    mi   = dtin.minute
    sec  = dtin.second
    
    if start_with_star:
        strt = "*  "
    else:
        strt = ""
        
    strout = strt + "{:4} {:2d} {:2d} {:2d} {:2d} {:11.8f}".format(yyyy,mm,dd,hh,mi,sec)
    return strout



def datestr_gins_filename_2_dt(datestrin):
    """
    Time representation conversion
    
    GINS filename time format => Python's Datetime 
    
    GINS filename time format : 'YYMMDD_HHMMSS'
    
    It is the time format used in the listing names
    
    Parameters
    ----------
    datestrin : string or list/numpy.array of string.
        GINS filename time format string(s).  
        Can handle several string in an iterable.
    
    Returns
    -------
    dtout : datetime or list of datetime.
        Datetime(s)
    """
    
    if utils.is_iterable(datestrin):
        typ=utils.get_type_smart(datestrin)
        return typ([datestr_gins_filename_2_dt(e) for e in datestrin])
    
    else:
        #### CASE WHERE THE DATE LOOKS LIKE 00000:00000

        yr = int(datestrin[0:2])
        mm = int(datestrin[2:4])
        dd = int(datestrin[4:6])
        
        hh   = int(datestrin[7:9])
        mmin = int(datestrin[9:11])
        ss   = int(datestrin[11:13])
    
        if yr > 50:
            year = 1900 + yr
        else:
            year = 2000 + yr
    
        return dt.datetime(year,mm,dd,hh,mmin,ss)

        
        hh   = int(datestrin[7:9])
        mmin = int(datestrin[9:11])
        ss   = int(datestrin[11:13])
    
        if yr > 50:
            year = 1900 + yr
        else:
            year = 2000 + yr
    
        return dt.datetime(year,mm,dd,hh,mmin,ss)






#### LEAP SECONDS MANAGEMENT
    
def leap_seconds(f):
    """
    Return a list of tuples of this format: (timestamp, number_of_seconds)
        timestamp: a 32-bit timestamp, seconds since the UNIX epoch
        number_of_seconds: how many leap-seconds occur at timestamp

    INTERNAL_FUNCTION

    Source
    ------
    http://stackoverflow.com/questions/19332902/extract-historic-leap-seconds-from-tzdata

    """
    TZFILE_MAGIC = 'TZif'.encode('US-ASCII')

    fmt = ">4s c 15x 6l"
    size = struct.calcsize(fmt)
    (tzfile_magic, tzfile_format, ttisgmtcnt, ttisstdcnt, leapcnt, timecnt,
        typecnt, charcnt) =  struct.unpack(fmt, f.read(size))
    #print("DEBUG: tzfile_magic: {} tzfile_format: {} ttisgmtcnt: {} ttisstdcnt: {} leapcnt: {} timecnt: {} typecnt: {} charcnt: {}".format(tzfile_magic, tzfile_format, ttisgmtcnt, ttisstdcnt, leapcnt, timecnt, typecnt, charcnt))

    # Make sure it is a tzfile(5) file
    assert tzfile_magic == TZFILE_MAGIC, (
            "Not a tzfile; file magic was: '{}'".format(tzfile_magic))

    # comments below show struct codes such as "l" for 32-bit long integer
    offset = (timecnt*4  # transition times, each "l"
        + timecnt*1  # indices tying transition time to ttinfo values, each "B"
        + typecnt*6  # ttinfo structs, each stored as "lBB"
        + charcnt*1)  # timezone abbreviation chars, each "c"

    f.seek(offset, 1) # seek offset bytes from current position

    fmt = '>{}l'.format(leapcnt*2)
    #print("DEBUG: leapcnt: {}  fmt: '{}'".format(leapcnt, fmt))
    size = struct.calcsize(fmt)
    data = struct.unpack(fmt, f.read(size))

    lst = [(data[i], data[i+1]) for i in range(0, len(data), 2)]
    assert all(lst[i][0] < lst[i+1][0] for i in range(len(lst)-1))
    assert all(lst[i][1] == lst[i+1][1]-1 for i in range(len(lst)-1))
    return lst

def print_leaps(leap_lst):
    """
    INTERNAL_FUNCTION
    """
    # leap_lst is tuples: (timestamp, num_leap_seconds)
    outlist = []
    for ts, num_secs in leap_lst:
        dtime = (dt.datetime.utcfromtimestamp(ts - num_secs+1))
        outlist.append((dtime,num_secs))
    return outlist

def get_leapsecond_frontend():
    """
    INTERNAL_FUNCTION

    Nota :
        Le temps universel (UT1) et le Temps atomique international (TAI) ont
        été définis comme égaux en 1958. Lors de la mise en place d’UTC en
        1972, UT1 s’était décalé d’environ 10 secondes par rapport au TAI.
        On choisit donc un décalage initial de 10 secondes entre UTC et TAI .

        the initial 10sec are added in find_leapsecond
    """
    zoneinfo_fname = '/usr/share/zoneinfo/right/UTC'

    try:
        ### Linux case : AUTO Mode
        with open(zoneinfo_fname, 'rb') as f:
            leap_lst = leap_seconds(f)
            final_leap_lis = print_leaps(leap_lst)
    except:
        ### Windows case : MANUAL Mode
        final_leap_lis = [(dt.datetime(1972, 7, 1, 0, 0), 1),
         (dt.datetime(1973, 1, 1, 0, 0), 2),
         (dt.datetime(1974, 1, 1, 0, 0), 3),
         (dt.datetime(1975, 1, 1, 0, 0), 4),
         (dt.datetime(1976, 1, 1, 0, 0), 5),
         (dt.datetime(1977, 1, 1, 0, 0), 6),
         (dt.datetime(1978, 1, 1, 0, 0), 7),
         (dt.datetime(1979, 1, 1, 0, 0), 8),
         (dt.datetime(1980, 1, 1, 0, 0), 9),
         (dt.datetime(1981, 7, 1, 0, 0), 10),
         (dt.datetime(1982, 7, 1, 0, 0), 11),
         (dt.datetime(1983, 7, 1, 0, 0), 12),
         (dt.datetime(1985, 7, 1, 0, 0), 13),
         (dt.datetime(1988, 1, 1, 0, 0), 14),
         (dt.datetime(1990, 1, 1, 0, 0), 15),
         (dt.datetime(1991, 1, 1, 0, 0), 16),
         (dt.datetime(1992, 7, 1, 0, 0), 17),
         (dt.datetime(1993, 7, 1, 0, 0), 18),
         (dt.datetime(1994, 7, 1, 0, 0), 19),
         (dt.datetime(1996, 1, 1, 0, 0), 20),
         (dt.datetime(1997, 7, 1, 0, 0), 21),
         (dt.datetime(1999, 1, 1, 0, 0), 22),
         (dt.datetime(2006, 1, 1, 0, 0), 23),
         (dt.datetime(2009, 1, 1, 0, 0), 24),
         (dt.datetime(2012, 7, 1, 0, 0), 25),
         (dt.datetime(2015, 7, 1, 0, 0), 26),
         (dt.datetime(2017, 1, 1, 0, 0), 27)]
    return final_leap_lis

def find_leapsecond(dtin,get_leapsec_lis=[],
                    apply_initial_delta=True):
    """
    Find the TAI-UTC leap second for a given datetime

    Parameters
    ----------
    dtin : datetime
        Epoch for which the leap second is researched

    get_leapsec_lis : list, optional
        A list of leap second, provided by get_leapsecond_frontend() 
        automatically determined if given list is empty
        
    apply_initial_delta : bool, optional
        See note below
                 
    Returns
    -------
    leapsec_out : int
        The leap second for the given epoch
        
    Note
    ----
    Le temps universel (UT1) et le Temps atomique international (TAI) ont
    été définis comme égaux en 1958. Lors de la mise en place d’UTC en
    1972, UT1 s’était décalé d’environ 10 secondes par rapport au TAI.
    On choisit donc un décalage initial de 10 secondes entre UTC et TAI .
    
    
    In 1972, the leap-second system was introduced so that the broadcast UTC 
    seconds could be made exactly equal to the standard SI second, while still
    maintaining the UTC time of day and changes of UTC date synchronized with 
    those of UT1 (the solar time standard that superseded GMT).[11] By then, 
    the UTC clock was already 10 seconds behind TAI, which had been 
    synchronized with UT1 in 1958, but had been counting true SI seconds
    since then. After 1972, both clocks have been ticking in SI seconds, 
    so the difference between their readouts at any time is 10 seconds plus 
    the total number of leap seconds that have been applied to UTC 
    (37 seconds as of January 2017).
    """

    if get_leapsec_lis == []:
        get_leapsec_lis = get_leapsecond_frontend()

    get_leapsec_lis2 = [ e[0] for e in get_leapsec_lis ]

    get_leapsec_lis2.append(dtin)
    i_dtin = sorted(get_leapsec_lis2).index(dtin)

    if apply_initial_delta:
        leapsec_out = i_dtin + 10 
    else:
        leapsec_out = i_dtin
        
    return leapsec_out


 #  _____       _   _                 _       _____       _                        _   _____                                     _        _   _                 
 # |  __ \     | | | |               ( )     |_   _|     | |                      | | |  __ \                                   | |      | | (_)                
 # | |__) |   _| |_| |__   ___  _ __ |/ ___    | |  _ __ | |_ ___ _ __ _ __   __ _| | | |__) |___ _ __  _ __ ___  ___  ___ _ __ | |_ __ _| |_ _  ___  _ __  ___ 
 # |  ___/ | | | __| '_ \ / _ \| '_ \  / __|   | | | '_ \| __/ _ \ '__| '_ \ / _` | | |  _  // _ \ '_ \| '__/ _ \/ __|/ _ \ '_ \| __/ _` | __| |/ _ \| '_ \/ __|
 # | |   | |_| | |_| | | | (_) | | | | \__ \  _| |_| | | | ||  __/ |  | | | | (_| | | | | \ \  __/ |_) | | |  __/\__ \  __/ | | | || (_| | |_| | (_) | | | \__ \
 # |_|    \__, |\__|_| |_|\___/|_| |_| |___/ |_____|_| |_|\__\___|_|  |_| |_|\__,_|_| |_|  \_\___| .__/|_|  \___||___/\___|_| |_|\__\__,_|\__|_|\___/|_| |_|___/
 #         __/ |                                                                                 | |                                                            
 #        |___/                                                                                  |_|                                                                           |_| 

### Python's Internal Representations  
    
def date2dt(date_in):
    """
    Time Python type conversion
    
    Python's Date => Python's Datetime
    
    Parameters
    ----------
    date_in : date or list/numpy.array of date.
        Date(s).  Can handle several dates in an iterable.
        
    Returns
    -------
    L : Datetime or list of Datetime
        Time as Datetime(s)  
    """
    if utils.is_iterable(date_in):
        typ=utils.get_type_smart(date_in)
        return typ([date2dt(e) for e in date_in])
    else:
        return dt.datetime(*tuple(date_in.timetuple())[:3])


def dt2date(dt_in):
    """
    Time Python type conversion
    
    Python's Datetime => Python's Date
    
    Parameters
    ----------
    dt_in : datetime or list/numpy.array of datetime.
        Datetime(s).  Can handle several datetimes in an iterable.
        
    Returns
    -------
    L : Datetime or list of Datetime
        Time as Datetime(s)  
    """
    if utils.is_iterable(dt_in):
        typ=utils.get_type_smart(dt_in)
        return typ([dt2date(e) for e in dt_in])
    else:
        return dt_in.date()


def pandas_timestamp2dt(timestamp_in):
    """
    Time Python type conversion
    
    Pandas's Timestamp => Python's Datetime
    
    Parameters
    ----------
    timestamp_in : Timestamp or list/numpy.array of Timestamp.
        Pandas's Timestamp(s).  Can handle several datetimes in an iterable.
        
    Returns
    -------
    L : Datetime or list of Datetime
        Time as Datetime(s)  
    """
    if utils.is_iterable(timestamp_in):
        typ=utils.get_type_smart(timestamp_in)
        return typ([pandas_timestamp2dt(e) for e in timestamp_in])
    else:
        if isinstance(timestamp_in, pd._libs.tslibs.timedeltas.Timedelta):
            return timestamp_in.to_pytimedelta()
        else:
            return timestamp_in.to_pydatetime()
        
def numpy_dt2dt(numpy_dt_in):
    """
    Time representation conversion

    Numpy datetime64 object => Python's Datetime
    
    Parameters
    ----------
    numpy_dt_in : numpy datetime64 object
        numpy datetime64 object

    Returns
    -------
    dt : datetime or list/numpy.array of datetime
        Converted Datetime(s)
        If the input is a Pandas Series, the output is forced as an array
              
    source
    ------
    https://gist.github.com/blaylockbk/1677b446bc741ee2db3e943ab7e4cabd
    """
    
    if utils.is_iterable(numpy_dt_in):        
        typ=utils.get_type_smart(numpy_dt_in)
        
        ### If the type is a Series, it has to be forced as an array
        ### Otherwise, datetime will be converted as numpy_dt again
        if typ == pd.core.series.Series:
            typ = np.array
            
        return typ([numpy_dt2dt(e) for e in numpy_dt_in])
    
    timestamp = ((numpy_dt_in - np.datetime64('1970-01-01T00:00:00'))
                 / np.timedelta64(1, 's'))
    return dt.datetime.utcfromtimestamp(timestamp)



##### Nota Bene
##### numpy_datetime2dt & datetime64_numpy2dt have been moved
##### to the funtion graveyard (PSakic 2021-02-22)


def date_to_jd(year,month,day):

    if month == 1 or month == 2:
        yearp = year - 1
        monthp = month + 12
    else:
        yearp = year
        monthp = month
    
    # this checks where we are in relation to October 15, 1582, the beginning
    # of the Gregorian calendar.
    if ((year < 1582) or
        (year == 1582 and month < 10) or
        (year == 1582 and month == 10 and day < 15)):
        # before start of Gregorian calendar
        B = 0
    else:
        # after start of Gregorian calendar
        A = math.trunc(yearp / 100.)
        B = 2 - A + math.trunc(A / 4.)
        
    if yearp < 0:
        C = math.trunc((365.25 * yearp) - 0.75)
    else:
        C = math.trunc(365.25 * yearp)
        
    D = math.trunc(30.6001 * (monthp + 1))
    
    jd = B + C + D + day + 1720994.5
    
    return jd
    

def jd_to_mjd(jd):
    """
    Convert Julian Day to Modified Julian Day
    
    Parameters
    ----------
    jd : float
        Julian Day
        
    Returns
    -------
    mjd : float
        Modified Julian Day
    
    """
    return jd - 2400000.5    

def mjd_to_jd(mjd):
    """
    Convert Modified Julian Day to Julian Day.
        
    Parameters
    ----------
    mjd : float
        Modified Julian Day
        
    Returns
    -------
    jd : float
        Julian Day
    
        
    """
    return mjd + 2400000.5

def jd_to_date(jd):
    """
    Convert Julian Day to date.
    
    
    Parameters
    ----------
    jd : float
        Julian Day
        
    Returns
    -------
    year : int
        Year as integer. Years preceding 1 A.D. should be 0 or negative.
        The year before 1 A.D. is 0, 10 B.C. is year -9.
        
    month : int
        Month as integer, Jan = 1, Feb. = 2, etc.
    
    day : float
        Day, may contain fractional part.
        
    Examples
    --------
    Convert Julian Day 2446113.75 to year, month, and day.
    
    >>> jd_to_date(2446113.75)
    (1985, 2, 17.25)
    
    """
    jd = jd + 0.5
    
    F, I = math.modf(jd)
    I = int(I)
    
    A = math.trunc((I - 1867216.25)/36524.25)
    
    if I > 2299160:
        B = I + 1 + A - math.trunc(A / 4.)
    else:
        B = I
        
    C = B + 1524
    
    D = math.trunc((C - 122.1) / 365.25)
    
    E = math.trunc(365.25 * D)
    
    G = math.trunc((C - E) / 30.6001)
    
    day = C - E + F - math.trunc(30.6001 * G)
    
    if G < 13.5:
        month = G - 1
    else:
        month = G - 13
        
    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715
        
    return year, month, day    


def hr_to_Day(hr,minu,sec):
    """
    Convert Julian Day to date.
    
    
    Parameters
    ----------
    jd : float
        Julian Day
        
    Returns
    -------
    year : int
        Year as integer. Years preceding 1 A.D. should be 0 or negative.
        The year before 1 A.D. is 0, 10 B.C. is year -9.
        
    month : int
        Month as integer, Jan = 1, Feb. = 2, etc.
    
    day : float
        Day, may contain fractional part.
        
    Examples
    --------
    Convert Julian Day 2446113.75 to year, month, and day.
    
    >>> jd_to_date(2446113.75)
    (1985, 2, 17.25)
    
    """
    hr = hr
    minu_hr = minu/60
    sec_hr = sec/3600
    hr_fim = hr+minu_hr+sec_hr
    dia_fim = hr_fim/24
        
    return dia_fim


def epo_epos_converter(inp,inp_type="mjd",out_type="yyyy",verbose=False):
    """
    Frontend for the GFZ EPOS epo converter
    
    
    Parameters
    ----------
    inp : string or int
        the input day, like 58773 2019290 20754 ...
        
    inp_type : str
        An output format managed by epo command :
        wwwwd,yyddd,yyyyddd,yyyymmdd
    
    out_type : str
        An output format managed by epo command :
        mjd,yyyy,yy,ddd,mon,dmon,hour,min,sec,wwww,wd


    Returns
    -------
    year : int
        Year as integer. Years preceding 1 A.D. should be 0 or negative.
        The year before 1 A.D. is 0, 10 B.C. is year -9.
    """

    epo_cmd = "-epo "  + str(inp)
    inp_cmd = "-type " + str(inp_type)
    out_cmd = "-o "    + str(out_type)

    cmd = " ".join(("perl $EPOS8_BIN_TOOLS/SCRIPTS/get_epoch.pl",epo_cmd,inp_cmd,out_cmd))

    if verbose:
        print(cmd)
    
    result = subprocess.run(cmd,shell=True,stdout=subprocess.PIPE,executable='/bin/bash')
            
    out = int(result.stdout.decode('utf-8'))
    
    if verbose:
        print(out)
    
    return out


 #  _______ _                   _____           _             _____       _                        _       _   _             
 # |__   __(_)                 / ____|         (_)           |_   _|     | |                      | |     | | (_)            
 #    | |   _ _ __ ___   ___  | (___   ___ _ __ _  ___  ___    | |  _ __ | |_ ___ _ __ _ __   ___ | | __ _| |_ _  ___  _ __  
 #    | |  | | '_ ` _ \ / _ \  \___ \ / _ \ '__| |/ _ \/ __|   | | | '_ \| __/ _ \ '__| '_ \ / _ \| |/ _` | __| |/ _ \| '_ \ 
 #    | |  | | | | | | |  __/  ____) |  __/ |  | |  __/\__ \  _| |_| | | | ||  __/ |  | |_) | (_) | | (_| | |_| | (_) | | | |
 #    |_|  |_|_| |_| |_|\___| |_____/ \___|_|  |_|\___||___/ |_____|_| |_|\__\___|_|  | .__/ \___/|_|\__,_|\__|_|\___/|_| |_|
 #                                                                                    | |                                    
 #                                                                                    |_|                                    



class interp1d_time(scipy.interpolate.interp1d):
    """
    Interpolation with datetime as inputs

    
    This class inherites from scipy.interpolate.interp1d
    and can take as input datetime as X
    
    Pierre Sakic 2020-01
    """
    def __init__(self, x, y, kind='linear', axis=-1,
                 copy=True, bounds_error=None, fill_value=np.nan,
                 assume_sorted=False):
        
        ### x (time) is converted as an array
        x = np.array(x)
        
        ### the datetime is converted to posix
        if isinstance(x[0],dt.datetime):
            xposix = dt2posix(x)
        elif isinstance(x[0],np.datetime64):
            xposix = dt2posix(numpy_dt2dt(x))
        else:
            xposix = x
        
        super().__init__(xposix, y, kind=kind, axis=axis,
                 copy=copy, bounds_error=bounds_error, fill_value=fill_value,
                 assume_sorted=assume_sorted)
        
    def __call__(self,x):
        ### x (time) is converted as an array if it is a mono-elt
        if not utils.is_iterable(x):
            singleton = True
            x = np.array([x])
        else:
            singleton = False
            x = np.array(x)

        ### the datetime is converted to posix
        if isinstance(x[0],dt.datetime):
            xposix = dt2posix(x)
        elif isinstance(x[0],np.datetime64):
            xposix = dt2posix(numpy_dt2dt(x))
        else:
            xposix = x
        
        #### Interpolation via the heritagle
        out = super().__call__(xposix)
        
        if singleton:
            return out[0]
        else:
            return out
    

   
class Slerp_time(scipy.spatial.transform.Slerp):
    """
    Slerp interpolation (for quaterinons) with datetime as inputs
    
    This class inherites from scipy.spatial.transform.Slerp
    and can take as input datetime as X
    
    P. Sakic 2020-01
    """
    
    def __init__(self, times, rotations,extrapolate=True):   
        
        ### time is converted as an array
        times = np.array(times)
        
        ### the datetime is converted to posix
        if isinstance(times[0],dt.datetime):
            times_posix = dt2posix(times)
        elif isinstance(times[0],np.datetime64):
            times_posix = dt2posix(numpy_dt2dt(times))
        else:
            times_posix = times
        
        #### For the extrapolation 
        # first value => begining of posix era
        # last value => end of posix era 
        if extrapolate:
            times_posix = np.array(times_posix)
            times_posix = np.insert(times_posix,0,0.)
            times_posix = np.append(times_posix,2147483646.)
            
            rotations_list = list(rotations)
            rotations_list.insert(0,rotations_list[0])
            rotations_list.append(rotations_list[-1])
            rotations = Rotation.from_quat([r.as_quat() for r in rotations_list])
            
        super().__init__(times_posix, rotations)
        
    def __call__(self,times):
        ### time is converted as an array
        times = np.array(times)
        
        if isinstance(times[0],dt.datetime):
            times_posix = dt2posix(times)
        elif isinstance(times[0],np.datetime64):
            times_posix = dt2posix(numpy_dt2dt(times))
        else:
            times_posix = times
            
        return super().__call__(times_posix)
    
    
 #  ______                _   _                _____                                         _ 
 # |  ____|              | | (_)              / ____|                                       | |
 # | |__ _   _ _ __   ___| |_ _  ___  _ __   | |  __ _ __ __ ___   _____ _   _  __ _ _ __ __| |
 # |  __| | | | '_ \ / __| __| |/ _ \| '_ \  | | |_ | '__/ _` \ \ / / _ \ | | |/ _` | '__/ _` |
 # | |  | |_| | | | | (__| |_| | (_) | | | | | |__| | | | (_| |\ V /  __/ |_| | (_| | | | (_| |
 # |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|  \_____|_|  \__,_| \_/ \___|\__, |\__,_|_|  \__,_|
 #                                                                        __/ |                
 #                                                                       |___/
  
# function graveyard  

def datetime64_numpy2dt(npdt64_in):
    """
    Time Python type conversion
    
    Numpy's datetime64 => Python's Datetime
    
    **This function is depreciated !!!!!**
    
    **use numpy_dt2dt instead !!!!!**
    
    Parameters
    ----------
    npdt64_in : datetime64 or list/numpy.array of datetime64.
        Numpy's datetime64(s). Can handle several datetimes in an iterable.
        
    Returns
    -------
    L : Datetime or list of Datetime
        Time as Datetime(s)  
        
    Source
    ------
    
    https://stackoverflow.com/questions/13703720/converting-between-datetime-timestamp-and-datetime64
    """
    print("WARN: datetime64_numpy2dt is depreciated, use numpy_dt2dt instead")
    warnings.warn("WARN: datetime64_numpy2dt is depreciated, use numpy_dt2dt instead",DeprecationWarning)


    if utils.is_iterable(npdt64_in):
        typ=utils.get_type_smart(npdt64_in)
        return typ([datetime64_numpy2dt(e) for e in npdt64_in])
    else:
        return pd.Timestamp(npdt64_in).to_pydatetime()



def numpy_datetime2dt(npdtin):
    """
    Time representation conversion
    
    Numpy Datetime => Datetime
    
    **This function is depreciated !!!**
    **Use numpy_dt2dt instead      !!!**

    Parameters
    ----------
    npdtin : np.datetime64 or list/numpy.array of np.datetime64 
        Numpy Datetime.  Can handle several time in a list.
                
    Returns
    -------
    python_datetime : datetime or list/numpy.array of datetime
        Converted Datetime(s)
        
    Source
    ------
    https://stackoverflow.com/questions/29753060/how-to-convert-numpy-datetime64-into-datetime/29755657
    """
    
    print("WARN: numpy_datetime2dt is depreciated, use numpy_dt2dt instead")
    warnings.warn("WARN: numpy_datetime2dt is depreciated, use numpy_dt2dt instead",DeprecationWarning)
    if utils.is_iterable(npdtin):        
        typ=utils.get_type_smart(npdtin)
        return typ([numpy_datetime2dt(e) for e in npdtin])
    else:
        python_datetime = npdtin.astype('M8[us]').astype('O') 
    return python_datetime


def dt_round(dtin=None, roundTo=60):
    """
    Round a datetime object to any time laps in seconds
    
    **This function is depreciated !!!**
    **Use round_dt instead         !!!**

    
    Parameters
    ----------
    dtin : datetime or list/numpy.array of datetime 
        Datetime you want to round, default now.
        Can handle several datetimes in an iterable.
            
    roundTo : int
        Closest number of seconds to round to, default 1 minute.
            
    Returns
    -------  
    dtout : datetime or list/numpy.array of datetime
        Rounded Datetime
        
    Note
    ---- 
    Based on :
    http://stackoverflow.com/questions/3463930/how-to-round-the-minute-of-a-datetime-object-python
    """
    import datetime as dtmod

    if utils.is_iterable(dtin):
        typ=utils.get_type_smart(dtin)
        return typ([dt_round(e,roundTo=roundTo) for e in dtin])
    else:
        if dtin == None :
            dtin = dtmod.datetime.now()
        seconds = (dtin - dtin.min).seconds
        # // is a floor division, not a comment on following line:
        rounding = (seconds+roundTo/2) // roundTo * roundTo
        return dtin + dtmod.timedelta(0,rounding-seconds,-dtin.microsecond)

def roundTime(*args):
    """
    Wrapper of dt_round for legacy reasons
    """
    return dt_round(*args)

    
def utc2gpstime(year,month,day,hour,min,sec):
    """
    looks useless and discontinued
    (210219)
    """
    
    date_utc = dt.datetime(year,month,day,hour,min,sec)
    
    utc_offset = find_leapsecond(date_utc)
    
    start_gps_time   = dt.datetime(1980,1,6)
    date_diff  = date_utc - start_gps_time + dt.timedelta(seconds=utc_offset) - dt.timedelta(seconds=19)
    
    gpsweek_decimal = date_diff.days / 7.
    gpsweek = np.floor(gpsweek_decimal)
    gpsweek_decimal_part = np.modf(gpsweek_decimal)[0]
    gpssecs = np.round(86400 * 7 * gpsweek_decimal_part) + date_diff.seconds
    
    return int(gpsweek),int(gpssecs)


def utc2gpstime_bad(year,month,day,hour,min,sec):
    from numpy import mod,array,floor
    """
    [gpsweek,gpssecs] = utc2gpstime(year,month,day,hour,min,sec)
    Converts UTC time to GPS week and seconds, Python version, Nov 2008
    
    UNSTABLE, DISCONTINUED !
    """

    # Number of days into the year at the start of each month (ignoring leap years).
    days_in_year = array([0,31,59,90,120,151,181,212,243,273,304,334,365])

    # year referenced from 1980
    ye = year - 1980;

    # Compute num of leap days
    leap_days = ye/4 + 1
    if mod(ye,4) == 0 and month <= 2:
        leap_days = leap_days - 1;

    # days elapsed since midnight jan 5, 1980.
    de = ye*365 + days_in_year[month-1] + day + leap_days - 6;

    # Leap Seconds, good after 1998
    # GPS time is ahead of UTC by this many Leap Seconds

    utc_offset = find_leapsecond(dt.datetime(year,month,day))

    # Desactivated because poor

    #    utc_offset = 13 # 1999 through 2005
    #    if year >= 2009:
    #        utc_offset = 15  # starting in 2009 to ???
    #    elif year >= 2006:
    #        utc_offset = 14 # 2006 through 2008
    #
    #    if year < 1999:
    #        print "Leap seconds invalid for this year!"
    #        return 0,0

    # Convert time to GPS weeks and seconds.
    gpsweek = floor(de/7.0)
    gpssecs = mod(de,7)*86400.0 + hour*3600.0 + min*60.0 + sec + utc_offset

    # Adjust GPS weeks/seconds to guarantee that seconds is in the range 0-604800.0. */
    if gpssecs < 0.0:
        gpsweek -= 1
        gpssecs += 604800.0

    if gpssecs >= 604800.0:
        gpsweek += 1
        gpssecs -= 604800.0

    return gpsweek,gpssecs


def gpstime2utc_bad(gpsweek,gpssecs,utc_offset):

    """
    DISCONTINUED FUNCTION, TOO UNSTABLE (171017)
    use gpstime2dt instead

    [year,month,day,hour,minute,sec] = gpstime2utc(gpsweek,gpssecs,utc_offset)
    Converts GPS week and seconds into UTC time, Python version, March 2009
    """

    from numpy import array,floor,fmod

    days_in_year = array([0,31,59,90,120,151,181,212,243,273,304,334,365])
    days_in_leap_year = array([0,31,60,91,121,152,182,213,244,274,305,335,366])
    # leap dayys since 1980
    leap_day_table = array([1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,
                            0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1]);

    # adjust for utc offset
    secs = gpssecs - utc_offset;

    # time as integer and fractional seconds since midnight Jan 5 1980.
    int_secs = floor(secs)
    fract_secs = secs - int_secs
    secs_1980 = 604800*gpsweek + int_secs

    # Express current GPS time as elapsed UTC days and seconds.
    day_number = floor(secs_1980/86400)
    sec_of_day = secs_1980 - day_number*86400

    leap_days1 = floor((day_number+1406) / 1461)

    # calculate UTC year
    total_years = floor((day_number + 5 - leap_days1) / 365)
    year = 1980 + total_years

    # day of utc year
    leap_days2 = sum(leap_day_table[0:total_years])

    day_of_utc_year = (day_number + 5 - 365*total_years - leap_days2)

    # determine month and day
    month = -1
    if fmod(year,4) != 0:
        # Not a leap year.
        for i in range(12):
            if day_of_utc_year > days_in_year[i]:
                month = i+1
                day = day_of_utc_year - days_in_year[i] + 1
    else:
        # A leap year.
        for i in range(12):
            if day_of_utc_year > days_in_leap_year[i]:
                month = i+1
                day = day_of_utc_year - days_in_leap_year[i] + 1

    # calculate utc hour
    hour = floor(sec_of_day/3600)

    if hour > 23:
        hour = 23

    # calculate utc minute
    minute = floor((sec_of_day - hour*3600) / 60)

    if minute > 59:
        minute = 59

    # calculate utc seconds
    sec = (sec_of_day - hour*3600 - minute*60) + fract_secs

    utc_time = [year,month,day,hour,minute,sec]

    return utc_time


def toYearFraction(date):
    """
    DISCONTINUED
    use dt2year_decimal instead (which is the same)
    """
    #
    # Give the decimal year
    # source :
    # http://stackoverflow.com/questions/6451655/python-how-to-convert-datetime-dates-to-decimal-years

    from datetime import datetime as dt
    import time

    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

def convert_partial_year(number):
    """
    DISCONTINUED
    use year_decimal2dt instead (which is the same)
    """

    # Source
    # http://stackoverflow.com/questions/19305991/convert-fractional-years-to-a-real-date-in-python

    import calendar
    from datetime import timedelta, datetime

    year = int(number)
    d = timedelta(days=(number - year)*(365 + calendar.isleap(year)))
    day_one = datetime(year,1,1)
    date = d + day_one
    return date
