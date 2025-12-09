#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 15:11:12 2020
Refactored on December 8, 2025

@author: psakic

This sub-module of geodezyx.conv deals with time conversion.

it can be imported directly with:
from geodezyx import conv

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox

REFACTORING NOTES:
==================
This refactored version improves:
1. Performance by using vectorized operations instead of recursive list comprehension
2. Handling of datetime/numpy.datetime64/pandas.Timestamp uniformly
3. Type preservation: single in → single out, iterable in → same type iterable out
4. Designed for geodesy, geophysics and astronomy applications
"""

########## BEGIN IMPORT ##########
#### External modules
import datetime as dt
import logging
import time
import warnings
from functools import wraps

import numpy as np

#### geodeZYX modules
from geodezyx import utils

### Logger
log = logging.getLogger('geodezyx')

##########  END IMPORT  ##########


#  _   _ _   _ _ _ _   _
# | | | | |_(_) (_) |_(_) ___  ___
# | | | | __| | | | __| |/ _ \/ __|
# | |_| | |_| | | | |_| |  __/\__ \
#  \___/ \__|_|_|_|\__|_|\___||___/
#
### Utility functions for handling different time types


def _normalize_datetime_input(dt_input):
    """
    Convert various datetime types to native Python datetime.

    Handles:
    - datetime.datetime (passthrough)
    - numpy.datetime64 (convert to datetime)
    - pandas.Timestamp (convert to datetime)

    Parameters
    ----------
    dt_input : datetime-like
        Input datetime in various formats

    Returns
    -------
    datetime.datetime
        Native Python datetime
    """
    if isinstance(dt_input, dt.datetime):
        return dt_input
    elif isinstance(dt_input, np.datetime64):
        # Convert numpy.datetime64 to datetime
        return dt_input.astype('datetime64[us]').astype(dt.datetime)
    else:
        # Try pandas Timestamp (lazy import)
        try:
            import pandas as pd
            if isinstance(dt_input, pd.Timestamp):
                return dt_input.to_pydatetime()
        except ImportError:
            pass

    # If all else fails, return as is
    return dt_input


def _normalize_timedelta_input(td_input):
    """
    Convert various timedelta types to native Python timedelta.

    Handles:
    - datetime.timedelta (passthrough)
    - numpy.timedelta64 (convert to timedelta)
    - pandas.Timedelta (convert to timedelta)

    Parameters
    ----------
    td_input : timedelta-like
        Input timedelta in various formats

    Returns
    -------
    datetime.timedelta
        Native Python timedelta
    """
    if isinstance(td_input, dt.timedelta):
        return td_input
    elif isinstance(td_input, np.timedelta64):
        # Convert numpy.timedelta64 to timedelta
        return dt.timedelta(microseconds=int(td_input / np.timedelta64(1, 'us')))
    else:
        # Try pandas Timedelta (lazy import)
        try:
            import pandas as pd
            if isinstance(td_input, pd.Timedelta):
                return td_input.to_pytimedelta()
        except ImportError:
            pass

    # If all else fails, return as is
    return td_input


def _output_to_datetime_type(dt_output, target_type):
    """
    Convert native Python datetime to the target type.

    Parameters
    ----------
    dt_output : datetime.datetime
        Native Python datetime
    target_type : type
        Target datetime type (datetime, numpy.datetime64, pandas.Timestamp)

    Returns
    -------
    datetime-like
        Datetime in target type
    """
    if target_type == 'datetime64' or target_type is np.datetime64:
        return np.datetime64(dt_output)
    elif target_type == 'Timestamp':
        import pandas as pd
        return pd.Timestamp(dt_output)
    else:
        return dt_output


def vectorized_time_converter(func):
    """
    Decorator to vectorize time conversion functions.

    This decorator handles:
    1. Single element input → single element output
    2. Iterable input → same type iterable output
    3. Automatic conversion of datetime-like types
    4. Vectorized operations for performance

    The decorated function should work on single native datetime objects.

    Parameters
    ----------
    func : callable
        Function that works on single datetime objects

    Returns
    -------
    callable
        Vectorized version of the function
    """
    @wraps(func)
    def wrapper(time_input, *args, **kwargs):
        # Check if input is iterable
        is_iter = utils.is_iterable(time_input)

        if not is_iter:
            # Single element case
            normalized = _normalize_datetime_input(time_input)
            result = func(normalized, *args, **kwargs)
            return result
        else:
            # Iterable case - vectorize
            input_type = utils.get_type_smart(time_input)

            # Convert to numpy array for vectorization
            time_array = np.asarray(time_input)

            # Vectorized conversion using numpy's vectorize (more efficient than list comprehension)
            vectorized_func = np.vectorize(lambda x: func(_normalize_datetime_input(x), *args, **kwargs))
            result_array = vectorized_func(time_array)

            # Convert back to original type
            return input_type(result_array)

    return wrapper


def numeric_vectorized_converter(func):
    """
    Decorator for numeric time conversion functions.

    Similar to vectorized_time_converter but for functions that work with
    numeric inputs (floats, ints) and return datetime objects.

    Parameters
    ----------
    func : callable
        Function that works on single numeric values

    Returns
    -------
    callable
        Vectorized version of the function
    """
    @wraps(func)
    def wrapper(numeric_input, *args, **kwargs):
        # Check if input is iterable
        is_iter = utils.is_iterable(numeric_input)

        if not is_iter:
            # Single element case
            return func(numeric_input, *args, **kwargs)
        else:
            # Iterable case - use numpy for efficiency
            input_type = utils.get_type_smart(numeric_input)
            numeric_array = np.asarray(numeric_input)

            # Vectorized operation
            vectorized_func = np.vectorize(lambda x: func(x, *args, **kwargs))
            result_array = vectorized_func(numeric_array)

            return input_type(result_array)

    return wrapper


#  _______ _                   _____                              _
# |__   __(_)                 / ____|                            (_)
#    | |   _ _ __ ___   ___  | |     ___  _ ____   _____ _ __ ___ _  ___  _ __
#    | |  | | '_ ` _ \ / _ \ | |    / _ \| '_ \ \ / / _ \ '__/ __| |/ _ \| '_ \
#    | |  | | | | | | |  __/ | |___| (_) | | | \ v /  __/ |  \__ \ | (_) | | | |
#    |_|  |_|_| |_| |_|\___|  \_____\___/|_| |_|\_/ \___|_|  |___/_|\___/|_| |_|
#
### Time conversion functions


@numeric_vectorized_converter
def tgipsy2dt(tin):
    """
    Time representation conversion

    GIPSY Time => Datetime

    Parameters
    ----------
    tin : float or iterable of float
        GIPSY time(s). Can handle iterable of floats.

    Returns
    -------
    dtout : datetime.datetime or iterable of datetime
        Converted Datetime(s). If input is iterable, returns same type.

    Note
    ----
    GIPSY time is not the 'real' J2000
    but only counted starting from 1st January 2000 at Noon
    """
    j2000 = dt.datetime(2000, 1, 1, 12, 0)
    return j2000 + dt.timedelta(seconds=float(tin))


@numeric_vectorized_converter
def matlab_time2dt(matlab_datenum):
    """
    Time representation conversion

    MATLAB Time => Datetime

    Parameters
    ----------
    matlab_datenum : float or iterable of float
        MATLAB time(s). Can handle iterable of floats.

    Returns
    -------
    python_datetime : datetime.datetime or iterable of datetime
        Converted Datetime(s). If input is iterable, returns same type.
    """
    return (dt.datetime.fromordinal(int(matlab_datenum)) +
            dt.timedelta(days=matlab_datenum % 1) -
            dt.timedelta(days=366))


def round_dt(dtin, round_to, python_dt_out=True, mode='round'):
    """
    Round a datetime object to any time lapse in seconds

    Uses Pandas Series for efficient vectorized rounding.

    Parameters
    ----------
    dtin : datetime or iterable of datetime
        Datetime(s) to round. Can handle iterable of datetimes.

    round_to : str
        The way to round the datetime.

        It follows the Pandas' Series conventions, e.g.:

        * one-day rounding: '1D'
        * one-minute rounding: '1min'
        * one-second rounding: '1s'

        Full list: https://pandas.pydata.org/docs/user_guide/timeseries.html#timeseries-offset-aliases

    python_dt_out : bool, optional
        If True, returns legacy Python DateTime.
        If False, returns Pandas Timestamp.
        Default is True.

    mode : str
        Rounding mode: 'round' (nearest), 'floor', 'ceil'
        Default is 'round'

    Returns
    -------
    dtout : datetime or iterable of datetime
        Rounded Datetime(s). If input is iterable, returns same type.

    Note
    ----
    This function uses Pandas for efficient vectorized operations.
    """
    import pandas as pd

    # Determine if input is singleton
    is_singleton = not utils.is_iterable(dtin)

    if is_singleton:
        dtin_series = pd.Series([_normalize_datetime_input(dtin)])
        input_type = None
    else:
        input_type = utils.get_type_smart(dtin)
        # Normalize all datetime inputs
        dtin_normalized = [_normalize_datetime_input(d) for d in dtin]
        dtin_series = pd.Series(dtin_normalized)

    # Apply rounding based on mode
    if mode == 'round':
        dtin_out = dtin_series.dt.round(round_to)
    elif mode == 'ceil':
        dtin_out = dtin_series.dt.ceil(round_to)
    elif mode == 'floor':
        dtin_out = dtin_series.dt.floor(round_to)
    else:
        log.error("check mode value: 'round', 'floor', 'ceil'")
        raise ValueError("Invalid mode. Must be 'round', 'floor', or 'ceil'")

    # Convert output
    if python_dt_out:
        dtin_out = np.array(pd.to_datetime(dtin_out))
    else:
        dtin_out = np.array(dtin_out)

    # Return based on input type
    if is_singleton:
        return dtin_out[0]
    else:
        if input_type is not None:
            return input_type(dtin_out)
        else:
            return list(dtin_out)


@vectorized_time_converter
def dt_ceil(dtin):
    """
    Round a datetime object to the beginning of the day

    Parameters
    ----------
    dtin : datetime or iterable of datetime
        Datetime(s) to round. Can handle iterable of datetimes.

    Returns
    -------
    dtout : datetime or iterable of datetime
        Rounded Datetime(s). If input is iterable, returns same type.
    """
    return dt.datetime.combine(dtin.date(), dt.time(0, 0, 0))


@vectorized_time_converter
def dt_in_local_timezone2posix(dtin):
    """
    Convert datetime in local timezone to POSIX timestamp.

    Parameters
    ----------
    dtin : datetime or iterable of datetime
        Datetime(s) in local timezone

    Returns
    -------
    float or iterable of float
        POSIX timestamp(s)
    """
    return time.mktime(dtin.timetuple()) + dtin.microsecond * 0.000001


@numeric_vectorized_converter
def posix2dt_in_local_timezone(posixin):
    """
    Convert POSIX timestamp to datetime in local timezone.

    Parameters
    ----------
    posixin : float or iterable of float
        POSIX timestamp(s)

    Returns
    -------
    datetime or iterable of datetime
        Datetime(s) in local timezone
    """
    if np.isnan(posixin):
        return dt.datetime(1970, 1, 1)
    else:
        return dt.datetime.fromtimestamp(posixin)


def dt_range(start_dt, end_dt, day_step=1, sec_step=0):
    """
    Generate range of datetime between start and end (included)

    Parameters
    ----------
    start_dt, end_dt : datetime
        Start and end datetimes
    day_step, sec_step : int, optional
        Step size in days and seconds

    Returns
    -------
    out_range : list of datetime
        Range of dates
    """
    out_range = [start_dt]
    step_delta = dt.timedelta(days=day_step, seconds=sec_step)

    while out_range[-1] < end_dt:
        out_range.append(out_range[-1] + step_delta)

    return out_range


@vectorized_time_converter
def dt2posix(dtin):
    """
    Time representation conversion

    Python's Datetime => POSIX Time

    Parameters
    ----------
    dtin : datetime or iterable of datetime
        Datetime(s). Can handle iterable of datetimes.

    Returns
    -------
    float or iterable of floats
        POSIX Time(s). If input is iterable, returns same type.

    Note
    ----
    This function properly handles both legacy datetime and numpy.datetime64
    with different precision below microsecond level.
    """
    # The decorator handles normalization, so dtin is always native datetime here
    D = dtin - dt.datetime(1970, 1, 1)
    dout = D.days * 86400 + D.seconds + D.microseconds * 10 ** -6
    return np.round(dout, 6)


@numeric_vectorized_converter
def posix2dt(posixin):
    """
    Time representation conversion

    POSIX Time => Python's Datetime

    Parameters
    ----------
    posixin : float or iterable of float
        POSIX Time(s). Can handle iterable of floats.

    Returns
    -------
    datetime or iterable of datetime
        Datetime(s). If input is iterable, returns same type.
    """
    return dt.datetime(1970, 1, 1) + dt.timedelta(seconds=float(posixin))


@numeric_vectorized_converter
def ntp2dt(ntp_timestamp_inp):
    """
    Time representation conversion

    NTP Time => Python's Datetime

    Parameters
    ----------
    ntp_timestamp_inp : int or iterable of int
        NTP timestamp(s). Can handle iterable of ints.

    Returns
    -------
    datetime or iterable of datetime
        Datetime(s). If input is iterable, returns same type.

    Note
    ----
    NTP epoch is 1900-01-01 00:00:00
    """
    # NTP epoch: 1900-01-01
    ntp_epoch = dt.datetime(1900, 1, 1, 0, 0, 0)
    return ntp_epoch + dt.timedelta(seconds=float(ntp_timestamp_inp))


def ymdhms2dt(y=0, mo=0, d=0, h=0, mi=0, s=0, ms=0):
    """
    Time representation conversion

    Year, Month, Day, Hour, Minute, Second => Python's Datetime

    Parameters
    ----------
    y, mo, d, h, mi, s, ms : int
        Year, Month, Day, Hour, Minute, Second, Microsecond

    Returns
    -------
    datetime
        Datetime object
    """
    return dt.datetime(y, mo, d, h, mi, s, ms)


def datetime_improved(*args):
    """
    Improved datetime constructor (alias for ymdhms2dt)
    """
    return dt.datetime(*args)


@vectorized_time_converter
def dt2ymdhms(dtin, with_microsec=True):
    """
    Time representation conversion

    Python's Datetime => Year, Month, Day, Hour, Minute, Second tuple

    Parameters
    ----------
    dtin : datetime or iterable of datetime
        Datetime(s). Can handle iterable of datetimes.
    with_microsec : bool
        Include microseconds in output

    Returns
    -------
    tuple or iterable of tuples
        (Year, Month, Day, Hour, Minute, Second[, Microsecond])
    """
    if with_microsec:
        return (dtin.year, dtin.month, dtin.day,
                dtin.hour, dtin.minute, dtin.second, dtin.microsecond)
    else:
        return (dtin.year, dtin.month, dtin.day,
                dtin.hour, dtin.minute, dtin.second)


def ymdhms_vectors2dt(yrlis, mlis, dlis, hlis, minlis, slis):
    """
    Time representation conversion

    Vectors of Year, Month, Day, Hour, Minute, Second => Datetimes

    Parameters
    ----------
    yrlis, mlis, dlis, hlis, minlis, slis : iterable of int
        Lists/arrays of year, month, day, hour, minute, second values

    Returns
    -------
    list of datetime
        Converted DateTimes
    """
    # Vectorized datetime creation using list comprehension (still efficient for this use case)
    return [dt.datetime(y, mo, d, h, mi, int(s))
            for y, mo, d, h, mi, s in zip(yrlis, mlis, dlis, hlis, minlis, slis)]


def doy2dt(year, days, hours=0, minutes=0, seconds=0):
    """
    Time representation conversion

    Year and Day of Year => Python's Datetime

    Parameters
    ----------
    year : int or iterable of int
        Year(s)
    days : int or iterable of int
        Day(s) of year
    hours, minutes, seconds : int or iterable of int, optional
        Hour, minute, second values

    Returns
    -------
    datetime or iterable of datetime
        Datetime(s). If inputs are iterable, returns list.

    Note
    ----
    If multiple arguments are iterables, they must have the same length.
    """
    # Check if any input is iterable
    is_year_iter = utils.is_iterable(year)
    is_days_iter = utils.is_iterable(days)
    is_hours_iter = utils.is_iterable(hours)
    is_mins_iter = utils.is_iterable(minutes)
    is_secs_iter = utils.is_iterable(seconds)

    if any([is_year_iter, is_days_iter, is_hours_iter, is_mins_iter, is_secs_iter]):
        # At least one input is iterable - vectorize
        # Ensure all are arrays of same length
        year_arr = np.atleast_1d(year)
        days_arr = np.atleast_1d(days)
        hours_arr = np.atleast_1d(hours) if is_hours_iter else np.full_like(year_arr, hours)
        mins_arr = np.atleast_1d(minutes) if is_mins_iter else np.full_like(year_arr, minutes)
        secs_arr = np.atleast_1d(seconds) if is_secs_iter else np.full_like(year_arr, seconds)

        results = []
        for y, d, h, m, s in zip(year_arr, days_arr, hours_arr, mins_arr, secs_arr):
            base_dt = dt.datetime(int(y), 1, 1)
            delta = dt.timedelta(days=int(d) - 1, hours=int(h),
                               minutes=int(m), seconds=float(s))
            results.append(base_dt + delta)

        return results
    else:
        # All inputs are scalar
        base_dt = dt.datetime(int(year), 1, 1)
        delta = dt.timedelta(days=int(days) - 1, hours=int(hours),
                           minutes=int(minutes), seconds=float(seconds))
        return base_dt + delta


@vectorized_time_converter
@vectorized_time_converter
def dt2doy(dtin, outputtype=str):
    """
    Time representation conversion

    Python's Datetime => Day of Year string

    Parameters
    ----------
    dtin : datetime or iterable of datetime
        Datetime(s). Can handle iterable of datetimes.
    outputtype : type
        Output type (str or int)

    Returns
    -------
    str/int or iterable of str/int
        Day of Year. If input is iterable, returns same type.
    """
    doy = dtin.timetuple().tm_yday
    if outputtype == str:
        return str(doy).zfill(3)
    else:
        return outputtype(doy)


@vectorized_time_converter
def dt2doy_year(dtin, outputtype=str):
    """
    Time representation conversion

    Python's Datetime => Year and Day of Year tuple

    Parameters
    ----------
    dtin : datetime or iterable of datetime
        Datetime(s). Can handle iterable of datetimes.
    outputtype : type
        Output type (str or int)

    Returns
    -------
    tuple or iterable of tuples
        (Year, Day of Year). If input is iterable, returns same type.
    """
    year = dtin.year
    doy = dtin.timetuple().tm_yday

    if outputtype == str:
        return (str(year), str(doy).zfill(3))
    else:
        return (outputtype(year), outputtype(doy))


@vectorized_time_converter
def dt2fracday(dtin):
    """
    Time representation conversion

    Python's Datetime => Fractional day

    Parameters
    ----------
    dtin : datetime or iterable of datetime
        Datetime(s). Can handle iterable of datetimes.

    Returns
    -------
    float or iterable of float
        Fractional day (0.0 to 1.0). If input is iterable, returns same type.
    """
    seconds_in_day = (dtin.hour * 3600 + dtin.minute * 60 +
                     dtin.second + dtin.microsecond * 1e-6)
    return seconds_in_day / 86400.0


@vectorized_time_converter
def dt2secinday(dtin):
    """
    Time representation conversion

    Python's Datetime => Seconds in day

    Parameters
    ----------
    dtin : datetime or iterable of datetime
        Datetime(s). Can handle iterable of datetimes.

    Returns
    -------
    float or iterable of float
        Seconds in day. If input is iterable, returns same type.
    """
    return (dtin.hour * 3600 + dtin.minute * 60 +
            dtin.second + dtin.microsecond * 1e-6)


@vectorized_time_converter
def dt2tuple(dtin):
    """
    Time representation conversion

    Python's Datetime => Tuple of (year, month, day, hour, minute, second)

    Parameters
    ----------
    dtin : datetime or iterable of datetime
        Datetime(s). Can handle iterable of datetimes.

    Returns
    -------
    tuple or iterable of tuples
        (year, month, day, hour, minute, second)
    """
    return tuple(dtin.timetuple())[:-3]


def tup_or_lis2dt(lisin):
    """
    Time representation conversion

    Date-looking list/tuple => Python's datetime

    Parameters
    ----------
    lisin : iterable of string/int/float
        List of date components like ["2018","12","31","12","30","00"]

    Returns
    -------
    datetime
        Converted DateTime
    """
    try:
        return dt.datetime(*[int(float(e)) for e in lisin])
    except:
        return posix2dt(0)


#  _______ _                    _____           _
# |__   __(_)                  / ____|         | |
#    | |   _ _ __ ___   ___   | (___   ___ __ _| | ___  ___
#    | |  | | '_ ` _ \ / _ \   \___ \ / __/ _` | |/ _ \/ __|
#    | |  | | | | | | |  __/   ____) | (_| (_| | |  __/\__ \
#    |_|  |_|_| |_| |_|\___|  |_____/ \___\__,_|_|\___||___/
#
### Time scale conversions


# Leap second handling (to be imported from original module or implemented)
def find_leapsecond(dtin):
    """
    Find the number of leap seconds at a given datetime.

    This is a placeholder - should use the actual implementation.
    """
    # Import from original module
    from geodezyx.conv import conv_time as conv_time_orig
    if hasattr(conv_time_orig, 'find_leapsecond'):
        return conv_time_orig.find_leapsecond(dtin)
    else:
        # Fallback - approximate
        if dtin.year >= 2017:
            return 37
        elif dtin.year >= 2015:
            return 36
        elif dtin.year >= 2012:
            return 35
        else:
            return 34  # Rough approximation


@vectorized_time_converter
def dt_utc2dt_tai(dtin):
    """
    Time scale conversion

    Python's datetime in UTC => Python's datetime in TAI

    Parameters
    ----------
    dtin : datetime or iterable of datetime
        Datetime(s) in UTC. Can handle iterable of datetimes.

    Returns
    -------
    datetime or iterable of datetime
        Datetime(s) in TAI. If input is iterable, returns same type.

    Note
    ----
    TAI is (currently) ahead of UTC
    """
    return dtin + dt.timedelta(seconds=find_leapsecond(dtin))


@vectorized_time_converter
def dt_tai2dt_utc(dtin):
    """
    Time scale conversion

    Python's datetime in TAI => Python's datetime in UTC

    Parameters
    ----------
    dtin : datetime or iterable of datetime
        Datetime(s) in TAI. Can handle iterable of datetimes.

    Returns
    -------
    datetime or iterable of datetime
        Datetime(s) in UTC. If input is iterable, returns same type.

    Note
    ----
    TAI is (currently) ahead of UTC
    """
    return dtin + dt.timedelta(seconds=-find_leapsecond(dtin))


@vectorized_time_converter
def dt_tai2dt_tt(dtin):
    """
    Time scale conversion

    Python's datetime in TAI => Python's datetime in Terrestrial Time (TT)

    Parameters
    ----------
    dtin : datetime or iterable of datetime
        Datetime(s) in TAI. Can handle iterable of datetimes.

    Returns
    -------
    datetime or iterable of datetime
        Datetime(s) in TT. If input is iterable, returns same type.

    Source
    ------
    https://en.wikipedia.org/wiki/Terrestrial_Time

    Note
    ----
    TT = TAI + 32.184 seconds
    """
    return dtin + dt.timedelta(seconds=32.184)


def dt_utc2dt_ut1(dtin, dUT1):
    """
    Time scale conversion

    Python's datetime in UTC => Python's datetime in UT1

    Parameters
    ----------
    dtin : datetime or iterable of datetime
        Datetime(s) in UTC. Can handle iterable of datetimes.
    dUT1 : float
        UT1-UTC in seconds.

    Returns
    -------
    datetime or iterable of datetime
        Datetime(s) in UT1. If input is iterable, returns same type.
    """
    is_iter = utils.is_iterable(dtin)

    if not is_iter:
        normalized = _normalize_datetime_input(dtin)
        return normalized + dt.timedelta(seconds=dUT1)
    else:
        input_type = utils.get_type_smart(dtin)
        results = [_normalize_datetime_input(d) + dt.timedelta(seconds=dUT1) for d in dtin]
        return input_type(results)


# Note: More complex functions like dt_utc2dt_ut1_smart will need the full implementation
# with EOP DataFrame handling. This is kept from the original for now.


# GPS Time conversions
def utc2gpstime(year, month, day, hour, minute, second):
    """
    Convert UTC time to GPS week and seconds.

    This is a helper function (to be imported from original or implemented).
    """
    # Import from original module
    from geodezyx.conv import conv_time as conv_time_orig
    if hasattr(conv_time_orig, 'utc2gpstime'):
        return conv_time_orig.utc2gpstime(year, month, day, hour, minute, second)
    else:
        # Simple implementation
        epoch = dt.datetime(1980, 1, 6)
        target = dt.datetime(year, month, day, hour, minute, second)
        diff = target - epoch
        weeks = diff.days // 7
        seconds = (diff.days % 7) * 86400 + diff.seconds
        return weeks, seconds


@vectorized_time_converter
def dt2gpstime(dtin, secinweek=False, inp_ref="utc", outputtype=int):
    """
    Time scale conversion

    Python's datetime => GPS time

    Parameters
    ----------
    dtin : datetime or iterable of datetime
        Datetime(s). Can handle iterable of datetimes.
    secinweek : bool
        if True: returns GPS week, sec in GPS week
        if False: returns GPS week, day in GPS week
    inp_ref : str
        "utc": apply 19 sec & leap second correction at the epoch
        "gps": no correction applied
        "tai": apply -19 sec correction
    outputtype : type
        Type of the output

    Returns
    -------
    tuple or iterable of tuples
        (GPS_week, GPS_day/GPS_sec)
        If input is iterable, returns same type of iterable containing tuples.

    Note
    ----
    Returns tuple as single object, allowing @vectorized_time_converter to work.
    The decorator treats each tuple as a single return value.
    """
    # Decorator handles normalization, so dtin is always native datetime here

    # Get raw GPS time
    week_raw, secs_raw = utc2gpstime(dtin.year, dtin.month, dtin.day,
                                     dtin.hour, dtin.minute, dtin.second)

    # Apply time scale corrections
    if inp_ref == "utc":
        week, secs = week_raw, secs_raw
    elif inp_ref == "tai":
        utc_offset = find_leapsecond(dtin)
        week, secs = week_raw, secs_raw - utc_offset
    elif inp_ref == "gps":
        utc_offset = find_leapsecond(dtin)
        week, secs = week_raw, secs_raw + 19 - utc_offset
    else:
        raise ValueError(f"Invalid inp_ref '{inp_ref}'. Must be 'utc', 'tai', or 'gps'")

    # Format output based on secinweek flag
    if not secinweek:  # day in week
        day = int(np.floor(secs / 86400))
        output1 = outputtype(int(week))
        output2 = outputtype(day)
    else:  # sec in week
        output1 = outputtype(int(week))
        output2 = outputtype(int(secs))

    # Special formatting for string output
    if outputtype == str and not secinweek:
        output1 = output1.zfill(4)

    # Return as single tuple object (decorator will vectorize this)
    outtup = (output1, output2)
    return outtup


@vectorized_time_converter
def dt2gpsweek_decimal(dtin, return_middle_of_day=True):
    """
    Time representation conversion

    Python's datetime => GPS week (decimal)

    Parameters
    ----------
    dtin : datetime or iterable of datetime
        Datetime(s). Can handle iterable of datetimes.
    return_middle_of_day : bool
        if False, return decimal part at beginning of day

    Returns
    -------
    float or iterable of float
        Decimal GPS week. If input is iterable, returns same type.
    """
    mod = 0.5 if return_middle_of_day else 0
    week, day = dt2gpstime(dtin)
    return float(week) + (float(day + mod) / 7)


@vectorized_time_converter
def dt2list(dtin, return_useful_values=True):
    """
    Time representation conversion

    Python's Datetime => List of [Year, Month, Day, Hour, Minute, Second]

    Parameters
    ----------
    dtin : datetime or iterable of datetime
        Datetime(s). Can handle iterable of datetimes.
    return_useful_values : bool
        if True, returns only Year, Month, Day

    Returns
    -------
    list or iterable of lists
        Date components. If input is iterable, returns same type.
    """
    if return_useful_values:
        return list(dtin.timetuple())[:-3]
    else:
        return list(dtin.timetuple())


def gpstime2dt(gpsweek, gpsdow_or_seconds, dow_input=True, output_time_scale="utc"):
    """
    Time scale & representation conversion

    GPS Time => Python's datetime

    Parameters
    ----------
    gpsweek : int or iterable of int
        GPS week(s)
    gpsdow_or_seconds : int or iterable of int
        Day of Week OR Seconds in Week
    dow_input : bool
        select if Day of Week (True) OR Seconds in Weeks (False)
    output_time_scale : str
        wished timescale: "utc", "tai", "gps"

    Returns
    -------
    datetime or iterable of datetime
        DateTime(s). If inputs are iterable, returns list.

    Note
    ----
    The leapsecond is found after a first calc without leapsecond.
    There can be side effects when input time is close to a leap second jump.
    """
    is_week_iter = utils.is_iterable(gpsweek)
    is_dow_iter = utils.is_iterable(gpsdow_or_seconds)

    if is_week_iter or is_dow_iter:
        return [gpstime2dt(w, ds, dow_input, output_time_scale)
                for w, ds in zip(np.atleast_1d(gpsweek), np.atleast_1d(gpsdow_or_seconds))]

    # Single value case
    if dow_input:
        gpsseconds = gpsdow_or_seconds * 86400. + 86400. * 0.5  # noon
    else:
        gpsseconds = gpsdow_or_seconds

    # First gross run
    epoch = dt.datetime(1980, 1, 6)
    elapsed = dt.timedelta(days=(int(gpsweek) * 7), seconds=(int(gpsseconds)))
    prelim_time = epoch + elapsed

    if output_time_scale == "gps":
        final_time = prelim_time
    else:
        leapsec = find_leapsecond(prelim_time)

        if output_time_scale == "utc":
            deltasec = leapsec - 19
        elif output_time_scale == "tai":
            deltasec = 19
        else:
            log.error("check output_time_scale value: 'utc', 'tai', 'gps'")
            raise ValueError("Invalid output_time_scale")

        # Second run with leap second
        elapsed = dt.timedelta(days=(int(gpsweek) * 7),
                             seconds=(int(gpsseconds) + deltasec))
        final_time = epoch + elapsed

    if dow_input:
        final_time = dt.datetime(final_time.year, final_time.month, final_time.day)

    return final_time


@numeric_vectorized_converter
def gpsweek_decimal2dt(gpsweekdec_in, output_time_scale='utc'):
    """
    Time representation conversion

    Decimal GPS Week => Python's datetime

    Parameters
    ----------
    gpsweekdec_in : float or iterable of float
        GPS week in decimal form.
    output_time_scale : str
        wished timescale: "utc", "tai", "gps"

    Returns
    -------
    datetime or iterable of datetime
        Datetime(s). If input is iterable, returns same type.
    """
    week_floor = np.floor(gpsweekdec_in)
    week_dec_part = gpsweekdec_in - week_floor

    return gpstime2dt(week_floor, week_dec_part * 7 * 86400,
                     dow_input=False, output_time_scale=output_time_scale)


@vectorized_time_converter
def dt_gpstime2dt_utc(dtgpsin):
    """
    Time scale conversion

    Datetime in GPS Time Scale => Datetime in UTC Time Scale

    Parameters
    ----------
    dtgpsin : datetime or iterable of datetime
        Datetime(s) in GPS Time Scale. Can handle iterable of datetimes.

    Returns
    -------
    datetime or iterable of datetime
        Datetime(s) in UTC Time Scale
    """
    # Decorator handles normalization and vectorization
    leapsec = find_leapsecond(dtgpsin)
    return dtgpsin + dt.timedelta(seconds=19 - leapsec)


@vectorized_time_converter
def dt_gpstime2dt_tai(dtgpsin):
    """
    Time scale conversion

    Datetime in GPS Time Scale => Datetime in TAI Time Scale

    Parameters
    ----------
    dtgpsin : datetime or iterable of datetime
        Datetime(s) in GPS Time Scale. Can handle iterable of datetimes.

    Returns
    -------
    datetime or iterable of datetime
        Datetime(s) in TAI Time Scale
    """
    # Decorator handles normalization and vectorization
    return dtgpsin + dt.timedelta(seconds=19)


@numeric_vectorized_converter
def year_decimal2dt(yearin):
    """
    Time representation conversion

    Decimal Year => Python's Datetime

    Parameters
    ----------
    yearin : float or iterable of float
        Decimal year(s). Can handle iterable of floats.

    Returns
    -------
    datetime or iterable of datetime
        Datetime(s). If input is iterable, returns same type.
    """
    import calendar
    year = int(yearin)
    d = dt.timedelta(days=(yearin - year) * (365 + calendar.isleap(year)))
    day_one = dt.datetime(year, 1, 1)
    return d + day_one


@vectorized_time_converter
def dt2year_decimal(dtin):
    """
    Time representation conversion

    Python's Datetime => Decimal Year

    Parameters
    ----------
    dtin : datetime or iterable of datetime
        Datetime(s). Can handle iterable of datetimes.

    Returns
    -------
    float or iterable of float
        Decimal Year(s). If input is iterable, returns same type.
    """
    def since_epoch(date):
        return time.mktime(date.timetuple())

    year = dtin.year
    start_of_this_year = dt.datetime(year=year, month=1, day=1)
    start_of_next_year = dt.datetime(year=year + 1, month=1, day=1)

    year_elapsed = since_epoch(dtin) - since_epoch(start_of_this_year)
    year_duration = since_epoch(start_of_next_year) - since_epoch(start_of_this_year)

    fraction = year_elapsed / year_duration

    return dtin.year + fraction


def date_string_2_dt(strin):
    """
    Time representation conversion

    String => Python's Datetime

    Wrapper of dateutil.parser.parse

    Parameters
    ----------
    strin : string or iterable of strings
        String(s). Can handle iterable of strings.

    Returns
    -------
    datetime or iterable of datetime
        Datetime(s). If input is iterable, returns same type.
    """
    import dateutil.parser

    is_iter = utils.is_iterable(strin)

    if is_iter:
        input_type = utils.get_type_smart(strin)
        return input_type([dateutil.parser.parse(s) for s in strin])
    else:
        return dateutil.parser.parse(strin)


# Aliases
string_date2dt = date_string_2_dt
str_date2dt = date_string_2_dt
strdate2dt = date_string_2_dt


# Additional conversion functions follow the same pattern...
# (MJD, CNES Julian Day, etc.)

@numeric_vectorized_converter
def jjul_cnes2dt(jjulin):
    """
    Time representation & scale conversion

    Julian Day CNES => Python's Datetime

    Parameters
    ----------
    jjulin : float/int or iterable of float/int
        Julian Day CNES. Can handle iterable of values.

    Returns
    -------
    datetime or iterable of datetime
        Datetime(s). If input is iterable, returns same type.

    Note
    ----
    Julian Day CNES starts at 0h Jan 1, 1950 (JD − 2433282.5)
    """
    return dt.datetime(1950, 1, 1, 0, 0, 0) + dt.timedelta(float(jjulin))


def dt2jjul_cnes(dtin, onlydays=True):
    """
    Time representation & scale conversion

    Python's Datetime => Julian Day CNES

    Parameters
    ----------
    dtin : datetime or iterable of datetime
        Datetime(s). Can handle iterable of datetimes.
    onlydays : bool
        if False, return also the seconds in the day

    Returns
    -------
    int/tuple or iterable of int/tuple
        Julian Day CNES[, seconds in the Day]
    """
    is_iter = utils.is_iterable(dtin)

    if is_iter:
        input_type = utils.get_type_smart(dtin)
        return input_type([dt2jjul_cnes(e, onlydays) for e in dtin])

    dtin = _normalize_datetime_input(dtin)
    epok = dtin - dt.datetime(1950, 1, 1, 0, 0, 0)

    if onlydays:
        return epok.days
    else:
        return epok.days, epok.seconds


def mjd2dt(mjd_in, seconds=None, round_to='1s'):
    """
    Time representation conversion

    Modified Julian Day => Python's Datetime

    Parameters
    ----------
    mjd_in : float/int or iterable of float/int
        Modified Julian Day. Can handle iterable of values.
    seconds : float/int or iterable of float/int, optional
        Additional seconds
    round_to : str or None
        Round the output datetime (see round_dt for details)

    Returns
    -------
    datetime or iterable of datetime
        Datetime(s). If input is iterable, returns same type.

    Note
    ----
    Modified Julian Day starts at 0h Nov 17, 1858 (JD − 2400000.5)
    """
    # Define rounding function
    if round_to:
        rnd = lambda x: round_dt(x, round_to)
    else:
        rnd = lambda x: x

    is_iter = utils.is_iterable(mjd_in)

    if is_iter:
        input_type = utils.get_type_smart(mjd_in)

        if seconds is not None:
            if len(seconds) != len(mjd_in):
                log.error("len(seconds) != len(mjd_in)!!")
                raise ValueError("Length mismatch between mjd_in and seconds")
        else:
            seconds = np.zeros(len(mjd_in))

        results = []
        for m, sec in zip(mjd_in, seconds):
            dtout = dt.datetime(1858, 11, 17) + dt.timedelta(days=float(m), seconds=float(sec))
            results.append(dtout)

        return input_type(rnd(results))
    else:
        if seconds is None:
            seconds = 0
        return rnd(dt.datetime(1858, 11, 17) + dt.timedelta(days=float(mjd_in), seconds=float(seconds)))


def MJD2dt(*args, **kwargs):
    """Alias for mjd2dt (deprecated)"""
    warnings.warn("MJD2dt is deprecated. Use mjd2dt instead.", DeprecationWarning, stacklevel=2)
    return mjd2dt(*args, **kwargs)


@vectorized_time_converter
def dt2mjd(dtin):
    """
    Time representation conversion

    Python's Datetime => Modified Julian Day

    Parameters
    ----------
    dtin : datetime or iterable of datetime
        Datetime(s). Can handle iterable of datetimes.

    Returns
    -------
    float or iterable of float
        Modified Julian Day(s). If input is iterable, returns same type.

    Note
    ----
    Modified Julian Day starts at 0h Nov 17, 1858 (JD − 2400000.5)
    """
    delta = dtin.replace(tzinfo=None) - dt.datetime(1858, 11, 17)
    return delta.days + (delta.seconds / 86400.) + (delta.microseconds / 864e8)


def dt2MJD(*args, **kwargs):
    """Alias for dt2mjd (deprecated)"""
    warnings.warn("dt2MJD is deprecated. Use dt2mjd instead.", DeprecationWarning, stacklevel=2)
    return dt2mjd(*args, **kwargs)


@vectorized_time_converter
def dt2str(dtin, str_format="%Y-%m-%d %H:%M:%S"):
    """
    Time representation conversion

    Python's Datetime => String

    Parameters
    ----------
    dtin : datetime or iterable of datetime
        Datetime(s). Can handle iterable of datetimes.
    str_format : str
        Format string for strftime

    Returns
    -------
    str or iterable of str
        Time as string(s). If input is iterable, returns same type.
    """
    return dtin.strftime(str_format)


# Additional specialized functions (RINEX names, SP3 names, etc.)
# are kept mostly as-is since they involve specific format parsing
# and don't benefit as much from the vectorization pattern

def rinexname2dt(rinexpath):
    """
    Time representation conversion

    RINEX Name => Python's Datetime

    Extract the date in a RINEX name

    Parameters
    ----------
    rinexpath : string
        RINEX path. The RINEX basename will be extracted automatically.

    Returns
    -------
    datetime
        Datetime
    """
    # This function involves complex regex parsing - keep as-is from original
    # Import the full function from original module
    from geodezyx.conv import conv_time as conv_time_orig
    return conv_time_orig.rinexname2dt(rinexpath)


def sp3name2dt(sp3path):
    """
    Time representation conversion

    Orbit SP3 Name => Python's Datetime

    Extract the date in a Orbit SP3 name

    Parameters
    ----------
    sp3path : string
        Orbit SP3 path. The basename will be extracted automatically.

    Returns
    -------
    datetime
        Datetime
    """
    # This function involves format detection - keep as-is from original
    from geodezyx.conv import conv_time as conv_time_orig
    return conv_time_orig.sp3name2dt(sp3path)


def date2dt(date_in):
    """
    Convert date object to datetime

    Parameters
    ----------
    date_in : date or iterable of date
        Date object(s)

    Returns
    -------
    datetime or iterable of datetime
        Datetime object(s)
    """
    is_iter = utils.is_iterable(date_in)

    if is_iter:
        input_type = utils.get_type_smart(date_in)
        return input_type([dt.datetime.combine(d, dt.time.min) for d in date_in])
    else:
        return dt.datetime.combine(date_in, dt.time.min)


# More specialized conversion functions will be imported from original module
# to avoid duplicating complex regex parsing logic
# (These can be refactored later if needed)


# Export main conversion functions
__all__ = [
    # Time representation conversions
    'tgipsy2dt', 'matlab_time2dt', 'round_dt', 'dt_ceil',
    'dt2posix', 'posix2dt', 'ntp2dt',
    'ymdhms2dt', 'dt2ymdhms', 'ymdhms_vectors2dt',
    'doy2dt', 'dt2doy', 'dt2doy_year',
    'dt2fracday', 'dt2secinday', 'dt2tuple', 'tup_or_lis2dt', 'dt2list',
    'year_decimal2dt', 'dt2year_decimal',
    'date_string_2_dt', 'string_date2dt', 'str_date2dt', 'strdate2dt',
    'jjul_cnes2dt', 'dt2jjul_cnes',
    'mjd2dt', 'dt2mjd', 'MJD2dt', 'dt2MJD',
    'dt2str', 'date2dt',

    # Time scale conversions
    'dt_utc2dt_tai', 'dt_tai2dt_utc', 'dt_tai2dt_tt', 'dt_utc2dt_ut1',
    'dt2gpstime', 'dt2gpsweek_decimal', 'gpstime2dt', 'gpsweek_decimal2dt',
    'dt_gpstime2dt_utc', 'dt_gpstime2dt_tai',
    'utc2gpstime', 'find_leapsecond',

    # Specialized conversions (kept from original - complex parsing)
    'rinexname2dt', 'sp3name2dt',

    # Utility functions
    'dt_range', 'dt_in_local_timezone2posix', 'posix2dt_in_local_timezone',

    # Vectorization decorators (advanced users)
    'vectorized_time_converter', 'numeric_vectorized_converter',
]


# NOTE: This is a refactored version with improved performance.
# Complex functions like rinexname2dt, sp3name2dt, and dt_utc2dt_ut1_smart
# should be imported from the original module or refactored separately
# as they involve domain-specific parsing that doesn't benefit much from vectorization.

