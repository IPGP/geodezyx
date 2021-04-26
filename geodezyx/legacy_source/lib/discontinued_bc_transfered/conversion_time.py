# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 14:12:44 2019

@author: psakicki

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics

Copyright (C) 2019 Pierre Sakic (GFZ, pierre.sakic@gfz-postdam.de)
GitHub repository :
https://github.com/PierreS1/GeodeZYX-Toolbox-Lite

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
"""
import datetime as dt
import time
import numpy as np
import os 
import string
import re
import struct

import conversion_general as cnv_gen

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


def rinex_regex_new_name(compressed=True,compiled=False):
    """
    Return a regex corresponding to a RINEX name (new convention)

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
        regexstr = ".{4}[0-9]{2}.{3}_(R|S|U)_[0-9]{11}_[0-9]{2}\w_[0-9]{2}\w_\w{2}\.\w{3}\.gz"
    else:
        regexstr = ".{4}[0-9]{2}.{3}_(R|S|U)_[0-9]{11}_[0-9]{2}\w_[0-9]{2}\w_\w{2}\.\w{3}"

    if compiled:
        return re.compile(regexstr)
    else:
        return regexstr
    
    

def tgipsy2dt(tin):
    """
    Time conversion
    
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
    if cnv_gen.is_iterable(tin):
        return [tgipsy2dt(e) for e in tin]
    else:
        j2000 = dt.datetime(2000, 1, 1, 12, 0)
        tout = j2000 + dt.timedelta(seconds=float(tin))
        return tout

    return tout

def numpy_datetime2dt(npdtin):
    """
    Time conversion
    
    Numpy Datetime => Datetime

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
    if cnv_gen.is_iterable(npdtin):
        return [numpy_datetime2dt(e) for e in npdtin]
    else:
        python_datetime = npdtin.astype('M8[ms]').astype('O') 
    return python_datetime



def matlab_time2dt(matlab_datenum):
    """
    Time conversion
    
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
    if cnv_gen.is_iterable(matlab_datenum):
        return [matlab_time2dt(e) for e in matlab_datenum]
    else:
        python_datetime = dt.datetime.fromordinal(int(matlab_datenum)) + \
        dt.timedelta(days=matlab_datenum%1) - dt.timedelta(days = 366)
    return python_datetime

def dt_round(dtin=None, roundTo=60):
    """
    Round a datetime object to any time laps in seconds
    
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

    if cnv_gen.is_iterable(dtin):
        return [dt_round(e,roundTo=roundTo) for e in dtin]
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
    

def dt_in_local_timezone2posix(dtin):
    if not cnv_gen.is_iterable(dtin):
        return time.mktime(dtin.timetuple()) + dtin.microsecond * 0.000001
    else:
        return [ dt_in_local_timezone2posix(e) for e in dtin ]


def posix2dt_in_local_timezone(posixin):
    if not cnv_gen.is_iterable(posixin):
        if np.isnan(posixin):
            return dt.datetime(1970,1,1)
        else:
            return dt.datetime.fromtimestamp(posixin)
    else:
        return [ posix2dt_in_local_timezone(e) for e in posixin ]


def dt2posix(dtin,out_array=False):       
    """
    Time conversion
    
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
    if not cnv_gen.is_iterable(dtin):
        D = dtin - dt.datetime(1970,1,1)
        return D.days * 86400 + D.seconds +  D.microseconds * 10**-6
    else:
        L = [ dt2posix(e) for e in dtin  ]
        if out_array:
            return np.array(L)
        else:
            return L

def posix2dt(posixin,out_array=False):       
    """
    Time conversion
    
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

    if not cnv_gen.is_iterable(posixin):
        if np.isnan(posixin):
            return dt.datetime(1970,1,1)
        else:
            return dt.datetime(1970,1,1) + dt.timedelta(seconds=posixin)
    else:
        L = [ posix2dt(e) for e in posixin ]
        if out_array:
            return np.array(L)
        else:
            return L


def datetime_improved(y=0,mo=0,d=0,h=0,mi=0,s=0,ms=0):
    """
    Improved time conversion of Datetime

    can manage when you give list, array, or floats as input
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
    try:
        ms_from_s  = (s - np.floor(s)) * 10**6
        if ms_from_s != 0:
            if ms != 0:
                print('WARN : input sec contains microsecs, given microsec are ignored')
            ms = ms_from_s
        return dt.datetime(int(y),int(mo),int(d),int(h),int(mi),int(s),int(ms))
    except:
        return dt.datetime(1970,1,1) # si ca deconne, si on donne un NaN par ex
            
def dt2ymdhms(dtin,with_microsec = True):
    """
    Time conversion

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
    if not cnv_gen.is_iterable(dtin):
        if with_microsec:
            return (dtin.year,dtin.month,dtin.day,dtin.hour,dtin.minute,dtin.second,dtin.microsecond)
        else:
            dt2 = roundTime(dtin,1)
            return (dt2.year,dt2.month,dt2.day,dt2.hour,dt2.minute,dt2.second)
    else:
        return [ dt2ymdhms(e,with_microsec) for e in dtin ]

def ymdhms_vectors2dt(yrlis,mlis,dlis,hlis,minlis,slis):
    """
    Time conversion
    
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
    Time conversion
    
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
    if not cnv_gen.is_iterable(year):
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
        if not cnv_gen.is_iterable(hours):
            hours   = [0] * len(year)
        if not cnv_gen.is_iterable(minutes):
            minutes = [0] * len(year)
        if  not cnv_gen.is_iterable(seconds):
            seconds = [0] * len(year)

        outlis = []
        for y,d,h,m,s in zip(year,days,hours,minutes,seconds):
            outlis.append(doy2dt(y,d,h,m,s))
        return outlis

def dt2doy(dtin,outputtype=str):  
    return outputtype(dtin.strftime('%j'))

def dt2doy_year(dtin,outputtype=str):
    """
    Time conversion
    
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
    
    if cnv_gen.is_iterable(dtin):
        return [dt2doy_year(e,outputtype=outputtype) for e in dtin] 
        
    else:
        return outputtype(dtin.strftime('%j')),outputtype(dtin.strftime('%Y'))
    
def dt2fracday(dtin):
    """
    Python's datetime => Seconds in days

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
    if cnv_gen.is_iterable(dtin):
        return [dt2fracday(e) for e in dtin]
    else:
        return dtin.hour / 24 + dtin.minute / (60*24) + dtin.second / (60*60*24)
    
def dt2secinday(dtin):
    """
    Time conversion
    
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
    if cnv_gen.is_iterable(dtin):
        return [dt2secinday(e) for e in dtin]
    else:
        return dtin.hour * 3600 + dtin.minute * 60 + dtin.second

def dt2tuple(dtin):
    return tuple(dtin.timetuple())[:-3]

def tup_or_lis2dt(lisin):
    """
    Time conversion
    
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

def gpstime2utc(gpsweek,gpssecs,utc_offset):

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
    leap_day_table = array([1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1]);

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



def utc2gpstime(year,month,day,hour,min,sec):
    
    date_utc = dt.datetime(year,month,day,hour,min,sec)
    
    utc_offset = find_leapsecond(date_utc)
    
    start_gps_time   = dt.datetime(1980,1,6)
    date_diff  = date_utc - start_gps_time + dt.timedelta(seconds=utc_offset) - dt.timedelta(seconds=19)
    
    gpsweek_decimal = date_diff.days / 7.
    gpsweek = np.floor(gpsweek_decimal)
    gpsweek_decimal_part = np.modf(gpsweek_decimal)[0]
    gpssecs = np.round(86400 * 7 * gpsweek_decimal_part) + date_diff.seconds
    
    return int(gpsweek),int(gpssecs)


def dt2gpstime(dtin,dayinweek=True,inp_ref="utc"):
    
    """
    Time conversion
    
    Python's datetime => GPS time

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
    
    if cnv_gen.is_iterable(dtin):
        return [dt2gpstime(e) for e in dtin]
        
    else:
        week_raw , secs_raw = utc2gpstime(dtin.year,dtin.month,dtin.day,dtin.hour,
                                          dtin.minute,dtin.second)
        
        utc_offset = find_leapsecond(dtin)

        if inp_ref == "utc":
            ### utc : utc2gpstime did the job
            week , secs = week_raw , secs_raw

        elif inp_ref == "tai":
            ### tai : utc2gpstime did the job, but needs leap sec correction again
            week , secs = week_raw , secs_raw - utc_offset

        elif inp_ref == "gps":
            ### tai : utc2gpstime did the job, but needs leap sec & 19sec correction again
            week , secs = week_raw , secs_raw + 19 - utc_offset
            
        if dayinweek:
            day = np.floor(np.divide(secs,86400))
            return int(week) , int(day)
        else:
            return int(week) , int(secs)


def dt2gpsweek_decimal(dtin,return_middle_of_day=True):
    """
    Time conversion
    
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
    if cnv_gen.is_iterable(dtin):
        return np.array([dt2gpsweek_decimal(e) for e in dtin])
    else:
        if return_middle_of_day:
            mod = 0.5
        else:
            mod = 0
        week , day = dt2gpstime(dtin)
        return float(week) + (float(day + mod) / 7 )



def dt2list(dtin,return_useful_values=True):
    """
    Time conversion
    
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



def gpstime2dt(gpsweek,gpsdow_or_seconds,dow_input = True):
    """
    Time conversion
    
    Day of Year Time => Python's datetime

    Parameters
    ----------
    gpsweek : int
        year, days of year
    gpsdow_or_seconds : int
        Day of Week OR Seconds in Weeks
    dow_input : bool
        select if Day of Week (True) OR Seconds in Weeks (False)
        
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
    
    if dow_input:
        gpsseconds = gpsdow_or_seconds * 86400 + 86400*.5 # so we are around noon
    else:
        gpsseconds = gpsdow_or_seconds

    ## First gross run
    epoch   = dt.datetime(1980,0o1,0o6)
    elapsed = dt.timedelta(days=(gpsweek*7),seconds=(gpsseconds+0))

    prelim_time = epoch + elapsed


    leapsec = find_leapsecond(prelim_time)

    ## Second run with leap second
    epoch   = dt.datetime(1980,0o1,0o6)
    elapsed = dt.timedelta(days=(gpsweek*7),seconds=(gpsseconds+leapsec))

    final_time = epoch + elapsed

    if dow_input:
        final_time = final_time.date()
        final_time = dt.datetime(final_time.year, final_time.month, final_time.day)

    return final_time


def gpsweek_decimal2dt(gpsweekdec_in):
    if cnv_gen.is_iterable(gpsweekdec_in):
        return [gpsweek_decimal2dt(e) for e in gpsweekdec_in]
    else:
        week_floor    = np.floor(gpsweekdec_in)
        week_dec_part = gpsweekdec_in - week_floor 
        
        dt_out = gpstime2dt(week_floor,week_dec_part * 7)
        return dt_out
    



def dt_gpstime2dt_utc(dtgpsin,out_array=False):
    """
    Time conversion
    
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
    if cnv_gen.is_iterable(dtgpsin):
        Out = [dt_gpstime2dt_utc(e) for e in dtgpsin]
        if not out_array:
            return Out
        else:
            return np.array(Out)
    else:
        leapsec = find_leapsecond(dtgpsin)
        dtutc = dtgpsin + dt.timedelta(seconds=19) - dt.timedelta(seconds=leapsec)
        return dtutc

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


def year_decimal2dt(yearin):
    """
    Time conversion
    
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
    if cnv_gen.is_iterable(yearin):
        return [year_decimal2dt(e) for e in yearin]
    else:
        return convert_partial_year(yearin)

def dt2year_decimal(dtin):
    """
    Time conversion
    
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
    if cnv_gen.is_iterable(dtin):
        return np.array([dt2year_decimal(e) for e in dtin])
    else:
        return toYearFraction(dtin)

def date_string_2_dt(strin):
    """
    Time conversion
    
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
    if cnv_gen.is_iterable(strin):
        return np.array([date_string_2_dt(e) for e in strin])
    else:
        return dateutil.parser.parse(strin)

string_date2dt = date_string_2_dt



def jjulCNES2dt(jjulin):
    """
    Time conversion
    
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
    
    if cnv_gen.is_iterable(jjulin):
        return [jjulCNES2dt(e) for e in jjulin]
    else:
        return dt.datetime(1950,0o1,0o1,00,00,00) + dt.timedelta(float(jjulin))

def dt2jjulCNES(dtin,onlydays=True):
    """
    Time conversion
    
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

def MJD2dt(mjd_in,seconds=None):
    # cf http://en.wikipedia.org/wiki/Julian_day
    """
    Time conversion
    
    Modified Julian Day  => Python's Datetime

    Parameters
    ----------
    mjd_in : float/int or list/numpy.array of float/int.
        Modified Julian Day.  Can handle several float/int in an iterable.

    Returns
    -------
    L : datetime or list of datetime.
        Datetime(s)

    Note
    ----
    Modified Julian Day starts at 0h Nov 17, 1858  (JD − 2400000.5)
    
    https://en.wikipedia.org/wiki/Julian_day
    """

    # ITERABLE CASE
    if cnv_gen.is_iterable(mjd_in): 
        if seconds:
            if len(seconds) != len(mjd_in):
                print("ERR : MJD2dt : len(seconds) != len(mjd_in) !!")
                raise Exception
        else:
            seconds = np.zeros(len(mjd_in))

        return [dt.datetime(1858,11,17) + dt.timedelta(days=m,seconds=sec) for m,sec in zip(mjd_in,seconds)]
    
    # NON ITERABLE / FLOAT CASE
    else:
        if not seconds:
            seconds = 0    
        return dt.datetime(1858,11,17) + dt.timedelta(days=mjd_in,seconds=seconds)


def dt2MJD(dtin):
    """
    Time conversion
    
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
        return [dt2MJD(t) for t in dtin]
    except:
        delta = (dtin - dt.datetime(1858,11,17))
        return delta.days + delta.seconds / 86400.

def dt2str(dtin , str_format="%Y-%m-%d %H:%M:%S"):
    """
    Time conversion
    
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
    
    if cnv_gen.is_iterable(dtin):
        return [dt2str(e) for e in dtin]
    else:
        return dtin.strftime(str_format)

def rinexname2dt(rinexpath):
    """
    Time conversion
    
    RINEX Name (old naming convention)  => Python's Datetime
    
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
        
    if re.search(rinex_regex_new_name(),rinexname):
        date_str = rinexname.split("_")[2]
        yyyy = int(date_str[:4])
        doy  = int(date_str[4:7])        
        hh   = int(date_str[7:9])        
        mm   = int(date_str[9:11])       
        dt_out = doy2dt(yyyy,doy) + dt.timedelta(seconds=hh*3600 + mm*60)
        return dt_out 
    else:
        alphabet = list(string.ascii_lowercase)
    
        if not re.search(rinex_regex(),rinexname):
            raise Exception('RINEX name is not well formated !!!')
        else:
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

        

def sp3name2dt(sp3path):
    """
    Time conversion
    
    Orbit SP3 Name (old naming convention)  => Python's Datetime
    
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
    Time conversion
    
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
    Time conversion
    
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
    Time conversion
    
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
    
    if cnv_gen.is_iterable(datestrin):
        return [datestr_sinex_2_dt(e) for e in datestrin]
    else:
        #### CASE WHERE THE DATE LOOKS LIKE 00000:00000
        if re.search("[0-9]{5}:[0-9]{5}",datestrin):
            datestr_list = list(datestrin)
            datestr_list.insert(2,":")
            datestrin = "".join(datestr_list)
        elif '00:000:00000' in datestrin:
            return dt.datetime(1970,1,1)
    
        dateint = [int(e) for e in datestrin.split(':')]
        yr = dateint[0]
        doy = dateint[1]
        sec = dateint[2]
        
        ## case for year with only 2 digits
        if re.match("[0-9]{2}:[0-9]{3}:[0-9]{5}",datestrin):            
            if yr > 50:
                year = 1900 + yr
            else:
                year = 2000 + yr
        else:
            year = yr
        
        return doy2dt(year,doy,seconds=sec)


def datestr_gins_filename_2_dt(datestrin):
    """
    Time conversion
    
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
    
    if cnv_gen.is_iterable(datestrin):
        return [datestr_gins_filename_2_dt(e) for e in datestrin]
    
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


def dt2sp3_timestamp(dtin,start_with_star=True):
    """
    Time conversion
    
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

def numpy_dt2dt(numpy_dt_in):
    """
    Time conversion

    numpy datetime64 object => Python's Datetime
    
    Parameters
    ----------
    numpy_dt_in : numpy datetime64 object
        numpy datetime64 object

    Returns
    -------
    dt : datetime
        Datetime
              
    source
    ------
    
    https://gist.github.com/blaylockbk/1677b446bc741ee2db3e943ab7e4cabd
    """
    timestamp = ((numpy_dt_in - np.datetime64('1970-01-01T00:00:00'))
                 / np.timedelta64(1, 's'))
    return dt.datetime.utcfromtimestamp(timestamp)




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
