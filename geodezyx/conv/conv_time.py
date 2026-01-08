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
import re
import os
import string

import numpy as np

#### geodeZYX modules
from geodezyx import utils
from geodezyx.conv import conv_time_leap_sec as leapsec
from geodezyx.conv import conv_rinex

### Logger
log = logging.getLogger("geodezyx")

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
        return dt_input.astype("datetime64[us]").astype(dt.datetime)
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
        return dt.timedelta(microseconds=int(td_input / np.timedelta64(1, "us")))
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
    if target_type == "datetime64" or target_type is np.datetime64:
        return np.datetime64(dt_output)
    elif target_type == "Timestamp":
        import pandas as pd

        return pd.Timestamp(dt_output)
    else:
        return dt_output


def _vectorize_1arg_core(func, input_val, *args, preprocessor=None, **kwargs):
    """
    Helper function to vectorize a function with a single vectorizable argument.

    Parameters
    ----------
    func : callable
        Function to vectorize
    input_val : scalar or iterable
        Input value(s) to process
    *args : tuple
        Additional positional arguments to pass to func
    preprocessor : callable, optional
        Function to preprocess each element before passing to func
    **kwargs : dict
        Additional keyword arguments to pass to func

    Returns
    -------
    scalar or iterable
        Result(s) in the same type as input
    """
    is_iter = utils.is_iterable(input_val)

    if not is_iter:
        # Single element case
        if preprocessor:
            input_val = preprocessor(input_val)
        return func(input_val, *args, **kwargs)
    else:
        # Iterable case - vectorize
        input_type = utils.get_type_smart(input_val)
        input_array = np.asarray(input_val)

        # Create vectorized function
        if preprocessor:
            vectorized_func = np.vectorize(
                lambda x: func(preprocessor(x), *args, **kwargs)
            )
        else:
            vectorized_func = np.vectorize(lambda x: func(x, *args, **kwargs))

        result_array = vectorized_func(input_array)
        return input_type(result_array)


def vector_datetime_conv(func):
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
        return _vectorize_1arg_core(
            func, time_input, *args, preprocessor=_normalize_datetime_input, **kwargs
        )

    return wrapper


def vector_numeric_conv(func):
    """
    Decorator for numeric time conversion functions.

    Similar to vector_datetime_conv but for functions that work with
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
        return _vectorize_1arg_core(func, numeric_input, *args, **kwargs)

    return wrapper

def vector_string_conv(func):
    """
    Decorator for string time conversion functions.

    Similar to vector_numeric_conv but for functions that work with
    string inputs and return datetime objects.

    Note: utils.is_iterable() treats strings as non-iterable by default,
    so we can use the common _vectorize_1arg_core helper.

    Parameters
    ----------
    func : callable
        Function that works on single string values

    Returns
    -------
    callable
        Vectorized version of the function
    """

    @wraps(func)
    def wrapper(string_input, *args, **kwargs):
        return _vectorize_1arg_core(func, string_input, *args, **kwargs)

    return wrapper

def _vectorize_2args_core(func, arg1, arg2, preprocessor_arg1=None, preprocessor_arg2=None, *args, **kwargs):
    """
    Core logic for vectorizing functions with 2 arguments.

    Parameters
    ----------
    func : callable
        Function to vectorize
    arg1 : scalar or iterable
        First argument
    arg2 : scalar or iterable or None
        Second argument (can be None)
    preprocessor_arg1 : callable, optional
        Preprocessor for first argument elements
    preprocessor_arg2 : callable, optional
        Preprocessor for second argument elements
    *args, **kwargs
        Additional arguments to pass to func

    Returns
    -------
    scalar or iterable
        Result(s) in the same type as input
    """
    input_type = lambda x: x  # Default passthrough
    # Check if arguments are iterable
    is_arg1_iter = utils.is_iterable(arg1)
    is_arg2_iter = arg2 is not None and utils.is_iterable(arg2)

    if not is_arg1_iter and not is_arg2_iter:
        # Single element case for both args
        if preprocessor_arg1:
            arg1 = preprocessor_arg1(arg1)
        if preprocessor_arg2 and arg2 is not None:
            arg2 = preprocessor_arg2(arg2)
        return func(arg1, arg2, *args, **kwargs)

    # At least one argument is iterable
    # Determine output type from the first iterable argument
    if is_arg1_iter:
        input_type = utils.get_type_smart(arg1)
        arg1_array = np.atleast_1d(arg1)
    else:
        arg1_array = None

    if is_arg2_iter:
        if not is_arg1_iter:
            # Only arg2 is iterable, use its type for output
            input_type = utils.get_type_smart(arg2)
        arg2_array = np.atleast_1d(arg2)
    else:
        # arg2 is not iterable (single value or None)
        arg2_array = None

    # Determine the iteration pairs
    if is_arg1_iter and is_arg2_iter:
        # Both are iterable - check length compatibility
        if len(arg1_array) != len(arg2_array):
            raise ValueError(
                f"Length mismatch: arg1 has {len(arg1_array)} elements, "
                f"arg2 has {len(arg2_array)} elements"
            )
        pairs = zip(arg1_array, arg2_array)
    elif is_arg1_iter:
        # Only arg1 is iterable, broadcast arg2
        pairs = [(a1, arg2) for a1 in arg1_array]
    else:
        # Only arg2 is iterable, broadcast arg1
        pairs = [(arg1, a2) for a2 in arg2_array]

    # Vectorized operation
    results = []
    for a1, a2 in pairs:
        # Apply preprocessors
        if preprocessor_arg1:
            a1 = preprocessor_arg1(a1)
        elif hasattr(a1, 'item'):
            # Convert numpy types to Python types
            a1 = a1.item()

        if preprocessor_arg2 and a2 is not None:
            a2 = preprocessor_arg2(a2)
        elif a2 is not None and hasattr(a2, 'item'):
            # Convert numpy types to Python types
            a2 = a2.item()

        result = func(a1, a2, *args, **kwargs)
        results.append(result)

    return input_type(results)


def vector_2args_numeric_conv(func):
    """
    Decorator for numeric time conversion functions with 2 vectorizable arguments.

    This decorator allows both the first and second arguments to be vectorized.
    If either or both arguments are iterable, the function is applied element-wise.
    The second argument can be None (optional).

    Parameters
    ----------
    func : callable
        Function that works on single numeric values for first 2 args

    Returns
    -------
    callable
        Vectorized version of the function

    Examples
    --------
    Typical use cases:
    - gpstime2dt(gpsweek, gpsdow_or_seconds, ...)
    - mjd2dt(mjd_in, seconds=None, ...)
    """

    @wraps(func)
    def wrapper(arg1, arg2=None, *args, **kwargs):
        return _vectorize_2args_core(func, arg1, arg2, *args, **kwargs)

    return wrapper


def vector_2args_datetime_numeric_conv(func):
    """
    Decorator for time conversion functions with datetime first arg and numeric second arg.

    This decorator allows both the first (datetime) and second (numeric) arguments
    to be vectorized. If either or both arguments are iterable, the function is
    applied element-wise.

    Parameters
    ----------
    func : callable
        Function that works on single datetime and numeric values

    Returns
    -------
    callable
        Vectorized version of the function

    Examples
    --------
    Typical use cases:
    - dt_utc2dt_ut1(dtin, d_ut1)
    - Functions with datetime + offset pattern
    """

    @wraps(func)
    def wrapper(arg1, arg2=None, *args, **kwargs):
        return _vectorize_2args_core(
            func, arg1, arg2,
            preprocessor_arg1=_normalize_datetime_input,
            *args, **kwargs
        )

    return wrapper


#  _______ _                   _____                              _
# |__   __(_)                 / ____|                            (_)
#    | |   _ _ __ ___   ___  | |     ___  _ ____   _____ _ __ ___ _  ___  _ __
#    | |  | | '_ ` _ \ / _ \ | |    / _ \| '_ \ \ / / _ \ '__/ __| |/ _ \| '_ \
#    | |  | | | | | | |  __/ | |___| (_) | | | \ v /  __/ |  \__ \ | (_) | | | |
#    |_|  |_|_| |_| |_|\___|  \_____\___/|_| |_|\_/ \___|_|  |___/_|\___/|_| |_|
#
### Time conversion functions


@vector_numeric_conv
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


@vector_numeric_conv
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
    return (
        dt.datetime.fromordinal(int(matlab_datenum))
        + dt.timedelta(days=matlab_datenum % 1)
        - dt.timedelta(days=366)
    )


def round_dt(dtin, round_to, python_dt_out=True, mode="round"):
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

    input_type = lambda x: x
    if is_singleton:
        dtin_series = pd.Series([_normalize_datetime_input(dtin)])
    else:
        input_type = utils.get_type_smart(dtin)
        # Normalize all datetime inputs
        dtin_norma = [_normalize_datetime_input(d) for d in dtin]
        dtin_series = pd.Series(dtin_norma)

    # Apply rounding based on mode
    if mode == "round":
        dt_use = dtin_series.dt.round(round_to)
    elif mode == "ceil":
        dt_use = dtin_series.dt.ceil(round_to)
    elif mode == "floor":
        dt_use = dtin_series.dt.floor(round_to)
    else:
        log.error("check mode value: 'round', 'floor', 'ceil'")
        raise ValueError("Invalid mode. Must be 'round', 'floor', or 'ceil'")

    # Convert output
    if python_dt_out:
        # Convert to Python datetime objects
        # Use list comprehension to ensure native datetime, not Timestamp
        dt_out = [x.to_pydatetime() for x in dt_use]
    else:
        # Keep as Pandas Timestamps
        dt_out = list(dt_use)

    # Return based on input type
    if is_singleton:
        return dt_out[0]
    else:
        if input_type is not None:
            return input_type(dt_out)
        else:
            return dt_out


@vector_datetime_conv
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


@vector_datetime_conv
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


@vector_numeric_conv
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


@vector_datetime_conv
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
    dout = D.days * 86400 + D.seconds + D.microseconds * 10**-6
    return np.round(dout, 6)


@vector_numeric_conv
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


@vector_numeric_conv
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


@vector_datetime_conv
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
        return (
            dtin.year,
            dtin.month,
            dtin.day,
            dtin.hour,
            dtin.minute,
            dtin.second,
            dtin.microsecond,
        )
    else:
        return dtin.year, dtin.month, dtin.day, dtin.hour, dtin.minute, dtin.second


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
    return [
        dt.datetime(y, mo, d, h, mi, int(s))
        for y, mo, d, h, mi, s in zip(yrlis, mlis, dlis, hlis, minlis, slis)
    ]


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
        hours_arr = (
            np.atleast_1d(hours) if is_hours_iter else np.full_like(year_arr, hours)
        )
        mins_arr = (
            np.atleast_1d(minutes) if is_mins_iter else np.full_like(year_arr, minutes)
        )
        secs_arr = (
            np.atleast_1d(seconds) if is_secs_iter else np.full_like(year_arr, seconds)
        )

        results = []
        for y, d, h, m, s in zip(year_arr, days_arr, hours_arr, mins_arr, secs_arr):
            base_dt = dt.datetime(int(y), 1, 1)
            delta = dt.timedelta(
                days=int(d) - 1, hours=int(h), minutes=int(m), seconds=float(s)
            )
            results.append(base_dt + delta)

        return results
    else:
        # All inputs are scalar
        base_dt = dt.datetime(int(year), 1, 1)
        delta = dt.timedelta(
            days=int(days) - 1,
            hours=int(hours),
            minutes=int(minutes),
            seconds=float(seconds),
        )
        return base_dt + delta


@vector_datetime_conv
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


@vector_datetime_conv
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
        return str(year), str(doy).zfill(3)
    else:
        return outputtype(year), outputtype(doy)


@vector_datetime_conv
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
    seconds_in_day = (
        dtin.hour * 3600 + dtin.minute * 60 + dtin.second + dtin.microsecond * 1e-6
    )
    return seconds_in_day / 86400.0


@vector_datetime_conv
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
    return dtin.hour * 3600 + dtin.minute * 60 + dtin.second + dtin.microsecond * 1e-6


@vector_datetime_conv
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


@vector_datetime_conv
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
    return dtin + dt.timedelta(seconds=leapsec.find_leapsecond(dtin))


@vector_datetime_conv
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
    return dtin + dt.timedelta(seconds=-leapsec.find_leapsecond(dtin))


@vector_datetime_conv
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


@vector_2args_datetime_numeric_conv
def dt_utc2dt_ut1(dtin, d_ut1):
    """
    Time scale conversion

    Python's datetime in UTC => Python's datetime in UT1

    Parameters
    ----------
    dtin : datetime or iterable of datetime
        Datetime(s) in UTC. Can handle iterable of datetimes.
    d_ut1 : float or iterable of float
        UT1-UTC in seconds. Can be iterable to match dtin.

    Returns
    -------
    datetime or iterable of datetime
        Datetime(s) in UT1. If input is iterable, returns same type.
    """
    return dtin + dt.timedelta(seconds=d_ut1)


def dt_utc2dt_ut1_smart(
    dtin, df_eop_in=None, use_interp1d_obj=True, eop_interpolator=None
):
    """
    Time scale conversion

    Python's datetime in UTC => Python's datetime in UT1
    using an EOP DataFrame provided by files_rw.read_eop_C04

    Parameters
    ----------
    dtin : datetime.datetime or list/numpy.array of datetime.datetime
        Datetime(s). Can handle several datetimes in an iterable.
    df_eop_in : DataFrame
        EOP DataFrame for the UT1-UTC
        provided by files_rw.read_eop_C04
    use_interp1d_obj : TYPE, optional
        Use an Interp1dTime Interpolator object for the EOP determination
        at the right epoch.
        Faster in recursive mode when dtin is a list/array
        The default is True.
    eop_interpolator : Interp1dTime object, optional
        The Interp1dTime Interpolator object for the EOP determination
        Will be determined automatically inside the function
        The default is None.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """

    if not df_eop_in:
        log.error("You must provide an EOP DataFrame")
        raise ValueError("EOP DataFrame is required")

    #### internally we work with "epoch" as index
    if df_eop_in.index.name != "epoch":
        df_eop = df_eop_in.set_index("epoch")
    else:
        df_eop = df_eop_in

    if utils.is_iterable(dtin):  #### ITERABLE CASE
        typ = utils.get_type_smart(dtin)

        ### if iterable, we optimize df_eop to speed up the fct
        dtmin = np.min(dtin)
        dtmax = np.max(dtin)

        boolt = (df_eop.index > dtmin - dt.timedelta(days=2)) & (
            df_eop.index < dtmax + dt.timedelta(days=2)
        )

        df_eop = df_eop[boolt]

        ### We also use the interpolator class
        if use_interp1d_obj:
            from geodezyx import interp

            ieop = interp.Interp1dTime(df_eop.index.values, df_eop["UT1-UTC"])
        else:
            ieop = None

        return typ(
            [
                dt_utc2dt_ut1_smart(
                    e, df_eop, use_interp1d_obj=use_interp1d_obj, eop_interpolator=ieop
                )
                for e in dtin
            ]
        )

    else:  #### SINGLE ELEMENT CASE
        from geodezyx import stats

        if dtin in df_eop.index:  ## the EOP value is directly in the EOP DF
            d_ut1 = df_eop.iloc[df_eop.index.get_loc(dtin)]["UT1-UTC"]
        elif (
            eop_interpolator and use_interp1d_obj
        ):  ### use the interp class, faster in recursive mode
            d_ut1 = eop_interpolator(dtin)
        else:  ## the EOP value is interpolated with a "manual" linear interpo.
            d_ut1bef = df_eop.iloc[df_eop.index.get_loc(dtin, "ffill")]
            d_ut1aft = df_eop.iloc[df_eop.index.get_loc(dtin, "bfill")]

            a, b, b2 = stats.linear_coef_a_b(
                d_ut1bef["MJD"],
                d_ut1bef["UT1-UTC"],
                d_ut1aft["MJD"],
                d_ut1aft["UT1-UTC"],
            )

            d_ut1 = stats.linear_reg_getvalue(dt2mjd(dtin), a, b, full=False)

        return dt_utc2dt_ut1(dtin, d_ut1)


# GPS Time conversions
@vector_datetime_conv
def utc2gpstime(dtin_utc):
    """
    Convert UTC Time to GPS Time

    Parameters
    ----------
    dtin_utc : datetime
        input UTC time.

    Returns
    -------
    gpsweek,gpssecs
        Converted epoch in GPS time, i.e. GPS Week and GPS seconds in week.

    """
    utc_offset = leapsec.find_leapsecond(dtin_utc)

    start_gps_time = dt.datetime(1980, 1, 6)
    date_diff = (
        dtin_utc
        - start_gps_time
        + dt.timedelta(seconds=utc_offset)
        - dt.timedelta(seconds=19)
    )

    gpsweek_decimal = date_diff.days / 7.0
    gpsweek = np.floor(gpsweek_decimal)
    gpsweek_decimal_part = np.modf(gpsweek_decimal)[0]
    gpssecs = np.round(86400 * 7 * gpsweek_decimal_part) + date_diff.seconds

    return int(gpsweek), int(gpssecs)


@vector_datetime_conv
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
    Returns tuple as single object, allowing @vector_datetime_conv to work.
    The decorator treats each tuple as a single return value.
    """
    # Decorator handles normalization, so dtin is always native datetime here

    # Get raw GPS time
    week_raw, secs_raw = utc2gpstime(dtin)

    # Apply time scale corrections
    if inp_ref == "utc":
        week, secs = week_raw, secs_raw
    elif inp_ref == "tai":
        utc_offset = leapsec.find_leapsecond(dtin)
        week, secs = week_raw, secs_raw - utc_offset
    elif inp_ref == "gps":
        utc_offset = leapsec.find_leapsecond(dtin)
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


@vector_datetime_conv
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


@vector_datetime_conv
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

@vector_2args_numeric_conv
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
    # Single value case (decorator handles vectorization)
    if dow_input:
        gpsseconds = gpsdow_or_seconds * 86400.0 + 86400.0 * 0.5  # noon
    else:
        gpsseconds = gpsdow_or_seconds

    # First gross run
    epoch = dt.datetime(1980, 1, 6)
    elapsed = dt.timedelta(days=(int(gpsweek) * 7), seconds=(int(gpsseconds)))
    prelim_time = epoch + elapsed

    if output_time_scale == "gps":
        final_time = prelim_time
    else:
        utc_leapsec = leapsec.find_leapsecond(prelim_time)

        if output_time_scale == "utc":
            deltasec = utc_leapsec - 19
        elif output_time_scale == "tai":
            deltasec = 19
        else:
            log.error("check output_time_scale value: 'utc', 'tai', 'gps'")
            raise ValueError("Invalid output_time_scale")

        # Second run with leap second
        elapsed = dt.timedelta(
            days=(int(gpsweek) * 7), seconds=(int(gpsseconds) + deltasec)
        )
        final_time = epoch + elapsed

    if dow_input:
        final_time = dt.datetime(final_time.year, final_time.month, final_time.day)

    return final_time


@vector_numeric_conv
def gpsweek_decimal2dt(gpsweekdec_in, output_time_scale="utc"):
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

    return gpstime2dt(
        week_floor,
        week_dec_part * 7 * 86400,
        dow_input=False,
        output_time_scale=output_time_scale,
    )


@vector_datetime_conv
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
    utc_leapsec = leapsec.find_leapsecond(dtgpsin)
    return dtgpsin + dt.timedelta(seconds=19 - utc_leapsec)


@vector_datetime_conv
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


@vector_numeric_conv
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


@vector_datetime_conv
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


@vector_string_conv
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

    return dateutil.parser.parse(strin)


# Aliases
string_date2dt = date_string_2_dt
str_date2dt = date_string_2_dt
strdate2dt = date_string_2_dt


# Additional conversion functions follow the same pattern...
# (MJD, CNES Julian Day, etc.)


@vector_numeric_conv
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


@vector_datetime_conv
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
    # At this point dtin is a native Python datetime (decorator normalizes/input-vectorizes)
    epok = dtin - dt.datetime(1950, 1, 1, 0, 0, 0)

    if onlydays:
        return epok.days
    else:
        return epok.days, epok.seconds


@vector_2args_numeric_conv
def mjd2dt(mjd_in, seconds=None, round_to="1s"):
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
    # Single value case (decorator handles vectorization)
    if seconds is None:
        seconds = 0

    dtout = dt.datetime(1858, 11, 17) + dt.timedelta(
        days=float(mjd_in), seconds=float(seconds)
    )

    # Apply rounding if requested
    if round_to:
        return round_dt(dtout, round_to)
    else:
        return dtout


def MJD2dt(*args, **kwargs):
    """Alias for mjd2dt (deprecated)"""
    warnings.warn(
        "MJD2dt is deprecated. Use mjd2dt instead.", DeprecationWarning, stacklevel=2
    )
    return mjd2dt(*args, **kwargs)


@vector_datetime_conv
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
    return delta.days + (delta.seconds / 86400.0) + (delta.microseconds / 864e8)


def dt2MJD(*args, **kwargs):
    """Alias for dt2mjd (deprecated)"""
    warnings.warn(
        "dt2MJD is deprecated. Use dt2mjd instead.", DeprecationWarning, stacklevel=2
    )
    return dt2mjd(*args, **kwargs)


@vector_datetime_conv
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

@vector_string_conv
def rinexname2dt(rinexpath):
    """
    Time representation conversion

    RINEX Name (short or long naming convention)  => Python's Datetime

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
    # rinexname = rinexpath

    ##### LONG rinex name
    if re.search(conv_rinex.rinex_regex_long_name(), rinexname) or re.search(
        conv_rinex.rinex_regex_long_name_brdc(), rinexname
    ):
        date_str = rinexname.split("_")[2]
        yyyy = int(date_str[:4])
        doy = int(date_str[4:7])
        hh = int(date_str[7:9])
        mm = int(date_str[9:11])
        dt_out = doy2dt(yyyy, doy) + dt.timedelta(seconds=hh * 3600 + mm * 60)
        return dt_out

        ##### LONG rinex name -- GFZ's GODC internal name
    elif re.search(conv_rinex.rinex_regex_long_name_gfz_godc(), rinexname):
        date_str = rinexname.split("_")[5]
        time_str = rinexname.split("_")[6]
        yyyy = int(date_str[:4])
        mo = int(date_str[4:6])
        dd = int(date_str[6:8])

        hh = int(time_str[0:2])
        mm = int(time_str[2:4])
        ss = int(time_str[4:6])

        dt_out = dt.datetime(yyyy, mo, dd, hh, mm, ss)

        return dt_out

        ##### SHORT rinex name
    elif re.search(
        conv_rinex.rinex_regex(), rinexname.lower()
    ):  ##EUREF are upper case...
        alphabet = list(string.ascii_lowercase)

        doy = int(rinexname[4:7])
        yy = int(rinexname[9:11])

        if yy > 80:
            year = yy + 1900
        else:
            year = yy + 2000

        if rinexname[7] in alphabet:
            h = alphabet.index(rinexname[7])
        else:
            h = 0

        return dt.datetime(year, 1, 1) + dt.timedelta(days=doy - 1, seconds=h * 3600)

    else:
        log.error("RINEX name is not well formated: %s", rinexname)
        return None

@vector_string_conv
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

    if len(sp3name) < 17:  ###### A REGEX WOULD BE MUCH BETTER !!!! (PSakic 2021-06)
        return sp3name_leg_2dt(sp3path)
    else:
        return sp3name_v3_2dt(sp3path)

@vector_string_conv
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
    dow = int(sp3name[7])
    return gpstime2dt(week, dow)

@vector_string_conv
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
    doy = int(datestr[4:7])
    hh = int(datestr[7:9])
    mm = int(datestr[9:11])

    dt_out = doy2dt(yyyy, doy) + dt.timedelta(hh * 3600 + mm * 60)

    return dt_out


def statname_dt2rinexname(
    statname, datein, rnxtype="d.Z", session_a_instead_of_daily_session=False
):
    """
    Time representation conversion

    Python's Datetime (and station name) => RINEX name

    Create a RINEX name from a station name and date

    Parameters
    ----------
    statname : string
        name of the station
        if is more than 4 characters long (typically for 9-chars. names),
        only the 4 first characters will be used

    datein : datetime.datetime
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

    if session_a_instead_of_daily_session:
        sess = alphabet[datein.hour]
    else:
        sess = "0"

    out_rnx_name = (
        statname[:4] + dt2doy(datein) + sess + "." + datein.strftime("%y") + rnxtype
    )
    out_rnx_name = out_rnx_name.lower()

    return out_rnx_name


def statname_dt2rinexname_long(
    statname,
    datein,
    country="XXX",
    data_source="R",
    file_period="00U",
    data_freq="00U",
    data_type="MO",
    format_compression="crx.gz",
    preset_type=None,
):
    """

    Time representation conversion

    Python's Datetime (and station name) => RINEX long name convention

    Create a RINEX name from a station name, date and ancillary parameters

    Parameters
    ----------
    statname : string
        name of the station
        can be a 4 or 9 char.

    datein : datetime.datetime
        date of the wished RINEX name

    country : string
        optional. the 3char. country code
        if given, will replace the one in an input 9-char statname
        default value is "XXX"

    data_source : str, optional
        Data Source.
        R – From Receiver, data using vendor or other software
        S – From data Stream (RTCM or other)
        U – Unknown
        The default is "R".

    file_period : str, optional
        File Period
        15M–15 Minutes
        01H–1 Hour
        01D–1 Day
        01Y–1 Year
        00U-Unspecified
        The default is "00U".

    data_freq : str, optional
        data frequency.
        None is allowed (for Navigation RINEX)

        XXC – 100 Hertz
        XXZ – Hertz,
        XXS – Seconds,
        XXM – Minutes,
        XXH – Hours,
        XXD – Days
        XXU – Unspecified
        The default is "00U".

    data_type : str, optional
        Two characters represent the data type:
        GO - GPS Obs.
        RO - GLONASS Obs.
        EO - Galileo Obs.
        JO - QZSS Obs.
        CO - BDS Obs.
        IO – NavIC/IRNSS Obs.
        SO - SBAS Obs.
        MO - Mixed Obs.
        GN - Nav. GPS
        RN - GLONASS Nav.
        EN - Galileo Nav.
        JN - QZSS Nav.
        CN - BDS Nav.
        IN – NavIC/IRNSS Nav.
        SN - SBAS Nav.
        MN – Mixed Nav. (All GNSS
        Constellations)
        MM-Meteorological Observation

        The default is "MO".

    format_compression : string
        extension of the RINEX ("rnx","crz","crz.gz") ...
        The default is 'crz.gz'.

    preset_type : str, optional
        takes "daily" or "hourly" values.
        set the most common data_freq and file_period
        for an hourly or daily session.
        The default is None.

    Returns
    -------
    out_rnx_name : str
        The name of the RINEX.

    """

    if not len(statname) in (4, 9):
        log.error("statname not 4 or 9 chars")
        raise Exception

    if len(statname) == 4:
        statname_ok = statname + "00" + country
    elif len(statname) == 9 and country != "XXX":
        statname_ok = statname[:4] + "00" + country
    else:
        statname_ok = statname

    if preset_type == "daily":
        file_period = "01D"
        data_freq = "30S"

    elif preset_type == "hourly":
        file_period = "01H"
        data_freq = "01S"

    date_ok = datein.strftime("%Y") + dt2doy(datein) + datein.strftime("%H%M")

    data_source_ok = "_" + data_source + "_"

    # for nav RINEX, data_freq can be None, thus we filter it
    elts_period_freq = [e for e in (file_period, data_freq, data_type) if e]
    period_freq_ok = "_" + "_".join(elts_period_freq)

    out_rnx_name = statname_ok + data_source_ok + date_ok + period_freq_ok
    out_rnx_name = out_rnx_name.upper() + "." + format_compression

    return out_rnx_name


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
        if re.search("[0-9]{7}:[0-9]{5}", datestr):
            datestr_list = list(datestr)
            datestr_list.insert(4, ":")
            datestr = "".join(datestr_list)
            #### CASE WHERE THE DATE LOOKS LIKE YYYYDDD:SSSSS
        elif re.search("[0-9]{5}:[0-9]{5}", datestr):
            datestr_list = list(datestr)
            datestr_list.insert(2, ":")
            datestr = "".join(datestr_list)

        ### this test must be independent
        if "00:000:00000" in datestr:
            return dt.datetime(1970, 1, 1)

        dateint = [int(e) for e in datestr.split(":")]
        yr = dateint[0]
        doy = dateint[1]
        sec = dateint[2]

        ## case for year with only 2 digits
        if re.match("[0-9]{2}:[0-9]{3}:[0-9]{5}", datestr):
            if yr > 50:
                year = 1900 + yr
            else:
                year = 2000 + yr
        else:
            year = yr

        return doy2dt(year, doy, seconds=sec)


def dt_2_sinex_datestr(dtin, short_yy=True, year_sep=":"):
    """
    Time representation conversion

    Python's Datetime => SINEX time format

    Parameters
    ----------
    dtin : datetime.datetime or list/numpy.array of datetime.datetime

    short_yy : bool
        if True, use a 2-digit year (e.g. 00:123:45678)
        if False, use a 4-digit year (e.g. 2000:123:45678)
        default is True.

    year_sep : str
        separator between the year and the day of year.
        default is ":".

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
        strout = strfmt.format(year, doy.zfill(3), sec.zfill(5))

        return strout


def dt_2_sp3_datestr(dtin):
    """
    Time representation conversion

    Python's Datetime => SP3 time format

    Parameters
    ----------
    dtin : datetime.datetime or list/numpy.array of datetime.datetime

    Returns
    -------
    dtout : str
        Date in SP3 format
    """

    if utils.is_iterable(dtin):
        return [dt_2_sp3_datestr(e) for e in dtin]
    else:
        return utils.join_improved("", *dt2gpstime(dtin))


def dt2sp3_timestamp(dtin, start_with_star=True):
    """
    Time representation conversion

    Python's Datetime => SP3 Timestamp
    e.g.

    *  2000  5 28  0  0  0.00000000
    """
    yyyy = dtin.year
    mm = dtin.month
    dd = dtin.day
    hh = dtin.hour
    mi = dtin.minute
    sec = dtin.second

    if start_with_star:
        strt = "*  "
    else:
        strt = ""

    strout = strt + "{:4} {:2d} {:2d} {:2d} {:2d} {:11.8f}".format(
        yyyy, mm, dd, hh, mi, sec
    )
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
    dtout : datetime.datetime or list/numpy.array of datetime.datetime.
        Datetime(s)
    """

    if utils.is_iterable(datestrin):
        typ = utils.get_type_smart(datestrin)
        return typ([datestr_gins_filename_2_dt(e) for e in datestrin])

    else:
        #### CASE WHERE THE DATE LOOKS LIKE 00000:00000
        yr = int(datestrin[0:2])
        mm = int(datestrin[2:4])
        dd = int(datestrin[4:6])

        hh = int(datestrin[7:9])
        mmin = int(datestrin[9:11])
        ss = int(datestrin[11:13])

        if yr > 50:
            year = 1900 + yr
        else:
            year = 2000 + yr

        return dt.datetime(year, mm, dd, hh, mmin, ss)


def trimble_file2dt(trmfile_in):
    """
    Time representation conversion

    Trimble raw data file => Python's Datetime

    Trimble file exemple : CASG202101070000A.T02

    Parameters
    ----------
    trmfile_in : string or list/numpy.array of string.
        Trimble raw data file as string.
        Can handle several string in an iterable.
        Will extract the basename of a path.

    Returns
    -------
    dtout : datetime or list of datetime.
        Datetime(s)
    """
    if utils.is_iterable(trmfile_in):
        typ = utils.get_type_smart(trmfile_in)
        return typ([trimble_file2dt(e) for e in trmfile_in])
    else:
        trmfile = os.path.basename(trmfile_in)
        year = int(trmfile[4:8])
        month = int(trmfile[8:10])
        day = int(trmfile[10:12])
        hh = int(trmfile[12:14])
        mm = int(trmfile[14:16])

        return dt.datetime(year, month, day, hh, mm, 0)


def dt2epoch_rnx3(dt_in, epoch_flag=0, nsats=0, rec_clk_offset=0):
    """
    Time representation conversion

    Python's Datetime => RINEX3/4 epoch representation
    e.g.
    > 2022 11 29 17 15 25.0000000  0  0       0.000000000000

    Parameters
    ----------
    dt_in : datetime.datetime or list/numpy.array of datetime.datetime.
        Datetime(s).
    epoch_flag : int, optional
        Epoch flag:
            0 : OK
            1 : power failure between previous and current epoch
            >1 : Special event (see RINEX documentation).
        The default is 0.
    nsats : int, optional
        Number of satellites observed in current epoch. The default is 0.
    rec_clk_offset : float, optional
        Receiver clock offset correction (seconds). The default is 0.

    Returns
    -------
    epoch_out : str or list of str
        output RINEX3/4 epoch representation.
    """

    if utils.is_iterable(dt_in):
        typ = utils.get_type_smart(dt_in)
        return typ([dt2epoch_rnx3(e, epoch_flag, nsats, rec_clk_offset) for e in dt_in])

    else:
        fmt = "> {:4} {:2n} {:2n} {:2n} {:2n}{:11.7f}  {:1n}{:3n}      {:15.12f}"
        epoch_out = fmt.format(
            dt_in.year,
            dt_in.month,
            dt_in.day,
            dt_in.hour,
            dt_in.minute,
            dt_in.second,
            epoch_flag,
            nsats,
            rec_clk_offset,
        )

        return epoch_out




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
        typ = utils.get_type_smart(date_in)
        return typ([date2dt(e) for e in date_in])
    else:
        return dt.datetime(*tuple(date_in.timetuple())[:3])


def dt2date(dt_in):
    """
    Time Python type conversion

    Python's Datetime => Python's Date

    Parameters
    ----------
    dt_in : datetime.datetime or list/numpy.array of datetime.datetime.
        Datetime(s).  Can handle several datetimes in an iterable.

    Returns
    -------
    L : Datetime or list of Datetime
        Time as Datetime(s)
    """
    if utils.is_iterable(dt_in):
        typ = utils.get_type_smart(dt_in)
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

    import pandas as pd

    if utils.is_iterable(timestamp_in):
        typ = utils.get_type_smart(timestamp_in)
        return typ([pandas_timestamp2dt(e) for e in timestamp_in])
    else:
        if isinstance(timestamp_in, pd.Timedelta):
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
    dt : datetime.datetime or list/numpy.array of datetime.datetime
        Converted Datetime(s)
        If the input is a Pandas Series, the output is forced as an array

    Source
    ------
    https://gist.github.com/blaylockbk/1677b446bc741ee2db3e943ab7e4cabd
    """

    if utils.is_iterable(numpy_dt_in):
        typ = utils.get_type_smart(numpy_dt_in)

        ### If the type is a pd.core.series.Series or a
        ### pd.core.indexes.datetimes.DatetimeIndex,
        ### it has to be forced as an array
        ### Otherwise, datetime will be converted as numpy_dt again
        if not typ in (list, tuple, np.array):
            typ = np.array
        return typ([numpy_dt2dt(e) for e in numpy_dt_in])

    #timestamp = ((numpy_dt_in - np.datetime64('1970-01-01T00:00:00'))
    #             / np.timedelta64(1, 's'))
    # return dt.datetime.fromtimestamp(timestamp)

    ### better implementation because of the timezone bug (PS 241124)

    else:
        import pandas as pd
        return pd.Timestamp(numpy_dt_in).to_pydatetime()

def time_obj_tester(delta=False, out_iterable=list, start=None):
    """
    Generate test time objects in different formats
    (datetime, numpy.datetime64, pandas.Timestamp).

    Parameters
    ----------
    delta : bool, optional
        If True, return time deltas relative to the `start` time.
        Default is False.
    out_iterable : callable, optional
        A callable that determines the type of iterable to return (e.g., list, tuple).
        Default is list.
    start : datetime.datetime, optional
        The starting datetime for the generated time objects.
        If None, defaults to January 1, 2020, 00:00:00.

    Returns
    -------
    tuple
        A tuple containing three iterables:
        - times_datetime: Iterable of datetime.datetime objects.
        - times_datetime_numpy: Iterable of numpy.datetime64 objects.
        - times_pandas_timestamp: Iterable of pandas.Timestamp objects.

    Notes
    -----
    - The function generates 10 time objects, each separated by 60 seconds.
    - If `delta` is True, the returned values are time deltas relative to the `start` time.
    - The `out_iterable` parameter allows customization of the output container type (e.g., list, tuple).

    Examples
    --------
    >>> time_obj_tester()
    ([datetime.datetime(2020, 1, 1, 0, 0),
      datetime.datetime(2020, 1, 1, 0, 1),
      ...],
     [numpy.datetime64('2020-01-01T00:00:00'),
      numpy.datetime64('2020-01-01T00:01:00'),
      ...],
     [Timestamp('2020-01-01 00:00:00'),
      Timestamp('2020-01-01 00:01:00'),
      ...])

    >>> time_obj_tester(delta=True, out_iterable=tuple)
    ((datetime.timedelta(0),
      datetime.timedelta(seconds=60),
      ...),
     (numpy.timedelta64(0,'s'),
      numpy.timedelta64(60,'s'),
      ...),
     (Timedelta('0 days 00:00:00'),
      Timedelta('0 days 00:01:00'),
      ...))
    """

    import pandas as pd

    if start is None:
        start = dt.datetime(2020, 1, 1, 0, 0, 0)

    # Generate a list of datetime objects, each separated by 60 seconds
    times_datetime = []
    for i in range(10):
        times_datetime.append(start + dt.timedelta(seconds=i * 60))

    # Convert the list of datetime objects to the specified iterable type
    times_datetime = out_iterable(times_datetime)

    # Convert the datetime objects to numpy.datetime64 and pandas.Timestamp formats
    times_datetime_numpy = out_iterable([np.datetime64(t) for t in times_datetime])
    times_pandas_timestamp = out_iterable([pd.Timestamp(t) for t in times_datetime])

    # If delta is True, calculate time deltas relative to the start time
    if delta:
        times_datetime = out_iterable([t - start for t in times_datetime])
        times_datetime_numpy = out_iterable([t - np.datetime64(start) for t in times_datetime_numpy])
        times_pandas_timestamp = out_iterable([t - pd.Timestamp(start) for t in times_pandas_timestamp])

    # Return the generated time objects in the specified formats
    return times_datetime, times_datetime_numpy, times_pandas_timestamp


##### Nota Bene
##### numpy_datetime2dt & datetime64_numpy2dt have been moved
##### to the funtion graveyard (PSakic 2021-02-22)

def epo_epos_converter(inp, inp_type="mjd", out_type="yyyy", verbose=False):
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

    verbose : bool
        verbose mode

    Returns
    -------
    year : int
        Year as integer. Years preceding 1 A.D. should be 0 or negative.
        The year before 1 A.D. is 0, 10 B.C. is year -9.
    """

    epo_cmd = "-epo " + str(inp)
    inp_cmd = "-type " + str(inp_type)
    out_cmd = "-o " + str(out_type)

    cmd = " ".join(("perl $EPOS8_BIN_TOOLS/SCRIPTS/get_epoch.pl", epo_cmd, inp_cmd, out_cmd))

    if verbose:
        log.debug(cmd)

    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, executable='/bin/bash')

    out = int(result.stdout.decode('utf-8'))

    if verbose:
        log.debug(out)

    return out