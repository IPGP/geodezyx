# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.stats contains functions for low-level
statistics.

it can be imported directly with:
from geodezyx import stats

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""

########## BEGIN IMPORT ##########
#### External modules
import datetime as dt

#### Import the logger
import logging

import numpy as np
import pandas as pd
import scipy

#### geodeZYX modules
from geodezyx import conv
from geodezyx import utils

log = logging.getLogger("geodezyx")

##########  END IMPORT  ##########


def linear_regression(x, y, fulloutput=False, simple_lsq=False, alpha=0.95):
    """
    Performs linear regression on two vectors, X and Y, and returns the coefficients a (slope) and b (intercept).

    Parameters
    ----------
    x : list or numpy.array
        The X values.
    y : list or numpy.array
        The Y values.
    simple_lsq : bool, optional
        If True, performs a basic, low-level least square inversion (faster, but less outputs).
        If False, calls scipy's linregress. Default is False.
    fulloutput : bool, optional
        If True, returns additional outputs (confidence interval for the slope and standard deviation).
        Default is False.
    alpha : float, optional
        The alpha value for the confidence interval. Default is .95.

    Returns
    -------
    slope : float
        The slope (a) of the linear regression.
    intercept : float
        The intercept (b) of the linear regression.
    confid_interval_slope : tuple of float, optional
        The confidence interval for the slope. Only returned if fulloutput is True.
    std_err : float, optional
        The standard deviation. Only returned if fulloutput is True.

    Notes
    -----
    This function performs a similar job to scipy.stats.linregress.

    Regarding computation speed: low-level least square inversion is faster for small datasets.
    For larger datasets, scipy's linregress is faster (n points > ~13000).

    """
    # Ensure x and y are numpy arrays
    if not isinstance(x, np.ndarray):
        x = np.array(x)
    if not isinstance(y, np.ndarray):
        y = np.array(y)

    # Check if lengths of x and y are equal
    if len(x) != len(y):
        log.error("len(x) (%i) != len(y) (%i)", len(x), len(y))
        return 0, 0

    # Perform least squares regression if simple_lsq is True
    if simple_lsq:
        a_arr = np.array(
            [x, np.ones(len(x))]
        )  # x2 faster than np.column_stack([x, np.ones(len(x))])
        slope, intercept = np.linalg.lstsq(a_arr.T, y, rcond=None)[
            0
        ]  # obtaining the parameters
        std_err = np.nan
    else:
        # Perform scipy's linregress if simple_lsq is False
        outtup = scipy.stats.linregress(x, y)
        slope, intercept, r_value, p_value, std_err = outtup

    # Return only the slope and intercept if fulloutput is False
    if simple_lsq or not fulloutput:
        return slope, intercept
    else:
        # Return the slope, intercept, confidence interval, and standard deviation if fulloutput is True
        return slope, intercept, confid_interval_slope(x, y, alpha), std_err


def linear_reg_getvalue(x, a, b, full=True):
    """
    Compute Y = a*X + b from a vector X and coefficients a and b.

    Parameters
    ----------
    x : list or numpy.ndarray
        Input values.
    a : float
        Linear regression slope coefficient.
    b : float
        Linear regression intercept coefficient.
    full : bool, optional
        If True, return both X and Y = aX + b.
        If False, return only Y = aX + b.
        Default is True.

    Returns
    -------
    y : numpy.ndarray
        Computed values Y = aX + b.
        Only returned if `full` is False.
    x : numpy.ndarray
        Input values (unchanged).
        Only returned if `full` is True.
    y : numpy.ndarray
        Computed values Y = aX + b.
        Only returned if `full` is True.

    Notes
    -----
    This function may be unstable when working with POSIX Time as X-data due to 
    large values. Decimal Years are recommended for better numerical stability.

    """

    if full:
        return np.array(x), a * np.array(x) + b
    else:
        return a * np.array(x) + b


def linear_coef_a_b(x1, y1, x2, y2):
    """
    Calculate line coefficients from two points.

    Calculates the slope (a) and intercept (b) coefficients of a line passing 
    through two points (x1, y1) and (x2, y2). Input values can be scalars or iterables.

    Parameters
    ----------
    x1 : float or list or numpy.ndarray
        X-coordinate(s) of the first point.
    y1 : float or list or numpy.ndarray
        Y-coordinate(s) of the first point.
    x2 : float or list or numpy.ndarray
        X-coordinate(s) of the second point.
    y2 : float or list or numpy.ndarray
        Y-coordinate(s) of the second point.

    Returns
    -------
    a : float or numpy.ndarray
        Slope coefficient of the line.
    b1 : float or numpy.ndarray
        Y-intercept using first point (should equal b2).
    b2 : float or numpy.ndarray
        Y-intercept using second point (should equal b1).

    """

    if utils.is_iterable(x1):
        x1 = np.array(x1, dtype=np.float64)
        x2 = np.array(x2, dtype=np.float64)
        y1 = np.array(y1, dtype=np.float64)
        y2 = np.array(y2, dtype=np.float64)
    else:
        x1 = float(x1)
        x2 = float(x2)
        y1 = float(y1)
        y2 = float(y2)

    a = (y2 - y1) / (x2 - x1)
    b1 = y1 - a * x1
    b2 = y2 - a * x2
    return a, b1, b2


def detrend_timeseries(x, y):
    """
    Remove linear trend from a time series.

    Removes the linear trend from Y(X) by subtracting the fitted line and 
    restoring the original starting value.

    Parameters
    ----------
    x : list or numpy.ndarray
        Independent variable (time or similar).
    y : list or numpy.ndarray
        Dependent variable (data values).

    Returns
    -------
    x : numpy.ndarray
        Independent variable (unchanged).
    yout : numpy.ndarray
        Detrended dependent variable.

    """

    x = np.array(x)
    y = np.array(y)
    a, b = linear_regression(x, y)

    ylinear = a * x + b
    yout = y - ylinear + y[0]

    return x, yout

def confid_interval_slope(x, y, alpha=0.95):
    """
    Calculate a confidence interval on the slope of a linear trend.

    Parameters
    ----------
    x : array_like
        Independent variable.
    y : array_like
        Dependent variable.
    alpha : float, optional
        Confidence level (default is 0.95 for 95% confidence).

    Returns
    -------
    mi : float
        Lower bound of the confidence interval for the slope.
    ma : float
        Upper bound of the confidence interval for the slope.

    Source
    -------
    Based on methods from:
    http://www.i4.auc.dk/borre/matlab
    http://kom.aau.dk/~borre/matlab/
    """

    sux = np.sum(x)
    xb = np.mean(x)
    suy = np.sum(y)
    yb = np.mean(y)
    n = len(x)
    s1 = np.sum(x * y)
    s2 = sux * suy / n
    sxy = s1 - s2
    s4 = np.sum(x**2)
    s5 = (sux**2) / n
    sxx = s4 - s5
    s7 = np.sum(y**2)
    s8 = (suy**2) / n
    syy = s7 - s8
    b1 = sxy / sxx
    b0 = yb - b1 * xb
    s14 = (sxy**2) / sxx
    s2y = (syy - s14) / (n - 2)
    sy = np.sqrt(s2y)
    s2b1 = s2y / sxx
    s2b0 = s2y * (1 / n + (xb**2) / sxx)
    # t=tq(1-alpha/2,n-2)
    t = scipy.stats.t.ppf(1 - alpha / 2, n - 2)
    mi = b1 - t * np.sqrt(s2b1)
    ma = b1 + t * np.sqrt(s2b1)
    return mi, ma


def running_mean(data_in, window):
    """
    Compute running mean (moving average) of data.

    Parameters
    ----------
    data_in : list or numpy.ndarray
        Input data values.
    window : float or int
        Size of the window for the running mean.

    Returns
    -------
    data_run : numpy.ndarray
        Running mean of the input data.

    """
    return pd.DataFrame(data_in).mean(window=window).values.flatten()


def running_mean_convolut(data_in, window, convolve_mode="same"):
    """
    Compute running mean using convolution mode.

    Computes the moving average of input data using convolution with a 
    specified mode. The result is mean-centered to avoid bias.

    Parameters
    ----------
    data_in : list or numpy.ndarray
        Input data values.
    window : float or int
        Size of the window for the running mean.
    convolve_mode : str, optional
        Mode for underlying convolution operation. Default is 'same'.

    Returns
    -------
    data_run : numpy.ndarray
        Running mean of data_in with same size as input (not shifted).

    Notes
    -----
    After stress testing, this implementation provides output with the same 
    size as input without shifting, though it is slower than some alternatives.
    The subtraction of the mean is an empirical trick to center the result.

    This is a wrapper based on running_mean_4. For more details on convolution 
    modes, see:
    - https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
    - https://stackoverflow.com/questions/11352047/finding-moving-average-from-data-points-in-python

    """

    data_mean = np.nanmean(data_in)
    data_zero_centered = data_in - data_mean

    data_run = running_mean_4(data_zero_centered, window, convolve_mode)

    return data_run + data_mean

# Running Means functions with INTERNAL ID 1 3 5
# return an Y output shorter than the input : not very convenient
# to align it on a X vector ...
#
# INTERNAL ID 2 gives Y with same size as X
# but result is shifted
# Y[0] should be aligned with the middle of the 1st window
#
# INTERNAL ID 4 is selected. It's declared as slow on StackOverflow,
# but at least the job is done.
# https://stackoverflow.com/questions/11352047/finding-moving-average-from-data-points-in-python
#
# About the speed, the ID4 is not the slowest in fact,
# it is the ID2 which is totally slow
# so this weakness has to be relativised
# (the answer is pretty old, so the convolve fct might have been improved)
#
# About the convolution mode, it is detailed here :
# https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
# "valid" mode is advised BUT do the same jobs as the others fcts (smaller output)
# The "same" mode is the best one for our applications.
# And "full" mode is not working either.
#
# BUT all the fcts are maintened here, because they can be usefull for some
# other cases !
#
# ##### EXEMPLE STRESS TEST CODE
# X = np.arange(1,1000,1) / (np.pi * 6)
#
# plt.clf()
# Ytrue = np.sin(X) * 100
# Y = Ytrue + np.random.randn(len(X)) * 100
#
# plt.plot(X,Y)
#
# N = 50
#
# Y1 = stats.movingaverage(Y,N)
# Y2 = stats.runningMean(Y,N)
# Y3 = stats.running_mean_core(Y,N)
# Y4a = stats.movingaverage_bis(Y,N,"same")
# Y4b = stats.movingaverage_bis(Y,N,"full")
# Y5 = stats.movingaverage_ter(Y,N)
#
# plt.clf()
# plt.plot(Y)
# plt.plot(Ytrue)
# plt.plot(Y1,"r.")
# plt.plot(Y2,"b.")
# plt.plot(Y3,"r.")
# plt.plot(Y4a,"y.")
# plt.plot(Y5,"g.")
#
#
# plt.clf()
# plt.plot(X,Y)
# plt.plot(X,Ytrue)
# #plt.plot(X,Y1,"r.")
# plt.plot(X,Y2,"b.")
# #plt.plot(X,Y3,"r.")
# plt.plot(X,Y4a,"y.")
# #plt.plot(X,Y4b,"r.")
# #plt.plot(X,Y5,"g.")


def running_mean_1(values, window):
    """
    Compute running mean using convolution with 'valid' mode.

    This method includes the 'valid' mode which requires enough datapoints. 
    For example, without 'valid', it would start at the first point with no prior 
    points, resulting in (1+0+0)/3 = 0.3333.

    Parameters
    ----------
    values : array-like
        Input data values.
    window : int
        Window size for the running mean.

    Returns
    -------
    smas : numpy.ndarray
        Running mean values using 'valid' convolution mode.

    Notes
    -----
    Internal ID: 1
    See: https://sentdex.com/sentiment-analysisbig-data-and-python-tutorials-algorithmic-trading/

    """
    weigths = np.repeat(1.0, window) / window
    smas = np.convolve(values, weigths, "valid")
    return smas  # as a numpy array


def running_mean_2(x, n):
    """
    Compute running mean by explicit windowing.

    Parameters
    ----------
    x : array-like
        Input data values.
    n : int
        Window size for the running mean.

    Returns
    -------
    y : numpy.ndarray
        Running mean values.

    Notes
    -----
    Internal ID: 2
    See: https://stackoverflow.com/questions/13728392/moving-average-or-running-mean

    """
    y = np.zeros((len(x),))
    for ctr in range(len(x)):
        y[ctr] = np.sum(x[ctr : (ctr + n)])
    return y / n


def running_mean_3(x, n):
    """
    Compute running mean using cumulative sum (moyenne glissante).

    Parameters
    ----------
    x : array-like
        Input data values.
    n : int
        Window size for the running mean.

    Returns
    -------
    xout : numpy.ndarray
        Running mean values.

    Notes
    -----
    Internal ID: 3
    See: https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
    (Alleo's answer)

    """
    cumsum = np.cumsum(np.insert(x, 0, 0))
    xout = (cumsum[n:] - cumsum[:-n]) / n

    return xout


def running_mean_4(interval, window_size, convolve_mode="same"):
    """
    Compute running mean using convolution with configurable mode.

    This method is slower than some alternatives but provides output with the same 
    size as input without shifting.

    Parameters
    ----------
    interval : array-like
        Input data values.
    window_size : int
        Size of the convolution window.
    convolve_mode : str, optional
        Convolution mode ('same', 'valid', 'full'). Default is 'same'.

    Returns
    -------
    result : numpy.ndarray
        Running mean values.

    Notes
    -----
    Internal ID: 4
    See: https://stackoverflow.com/questions/11352047/finding-moving-average-from-data-points-in-python

    """
    window = np.ones(int(window_size)) / float(window_size)
    return np.convolve(interval, window, convolve_mode)


def running_mean_5(data, window_width):
    """
    Compute running mean using cumulative sum method.

    Parameters
    ----------
    data : array-like
        Input data values.
    window_width : int
        Width of the running mean window.

    Returns
    -------
    ma_vec : numpy.ndarray
        Moving average vector.

    Notes
    -----
    Internal ID: 5
    See: https://stackoverflow.com/questions/11352047/finding-moving-average-from-data-points-in-python
    (Roman Kh's answer)

    """
    cumsum_vec = np.cumsum(np.insert(data, 0, 0))
    ma_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width

    return ma_vec


def sinusoide(t, a, omega, phi=0, f=None):
    """
    Generate a sinusoidal waveform.

    Produces a sinusoidal waveform of the form: A * sin(ω*T + φ)

    Parameters
    ----------
    t : float or array-like
        Time variable.
    a : float
        Amplitude, the peak deviation of the function from zero.
    omega : float
        Angular frequency (ω = 2πf), the rate of change of the function 
        argument in units of radians per second.
    phi : float, optional
        Phase offset in radians, specifying where in its cycle the 
        oscillation is at t = 0. Default is 0.
    f : float, optional
        Ordinary frequency (cycles per second). If provided, it overrides 
        the `omega` parameter. To use this parameter, set omega=0. 
        Default is None.

    Returns
    -------
    result : float or numpy.ndarray
        The sinusoidal waveform value(s).

    Notes
    -----
    See: https://en.wikipedia.org/wiki/Sine_wave

    """

    if f:
        omega_use = 2 * np.pi * f
    else:
        omega_use = omega

    return a * np.sin(omega_use * t + phi)


def butter_lowpass(cutoff, fs, order=5):
    """
    Design a Butterworth lowpass digital filter.

    Parameters
    ----------
    cutoff : float
        Critical frequency (Hertz). Frequencies above this will be attenuated.
    fs : float
        Sampling frequency (Hertz).
    order : int, optional
        Order of the filter. Default is 5.

    Returns
    -------
    b : numpy.ndarray
        Numerator (zeros) of the IIR filter.
    a : numpy.ndarray
        Denominator (poles) of the IIR filter.

    See Also
    --------
    scipy.signal.butter : Design IIR filters
    butter_lowpass_filter : Apply the designed filter to data

    """
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = scipy.signal.butter(order, normal_cutoff, btype="low", analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order=5):
    """
    Apply Butterworth lowpass digital filter to data.

    Parameters
    ----------
    data : array-like
        Input signal to be filtered.
    cutoff : float
        Critical frequency (Hertz). Frequencies above this will be attenuated.
    fs : float
        Sampling frequency (Hertz).
    order : int, optional
        Order of the filter. Default is 5.

    Returns
    -------
    y : numpy.ndarray
        Filtered signal.

    See Also
    --------
    butter_lowpass : Design the lowpass filter
    scipy.signal.butter : Design IIR filters
    scipy.signal.lfilter : Apply IIR filters

    """
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = scipy.butter.lfilter(b, a, data)
    return y

def gaussian_filter_gfz(tim_ref, dat_ref, width=7):
    """
    Apply Gaussian filter to smooth data.

    Gaussian filter based on GFZ's GMT_plus.pm/gaussian_kernel. 
    Smooths data by weighted averaging using Gaussian weights.

    Parameters
    ----------
    tim_ref : array-like (list or numpy.ndarray)
        X/T component of the time series (in decimal days).
    dat_ref : array-like (list or numpy.ndarray)
        Y component (the data values).
    width : int, optional
        Size of the smoothing window. An odd number is recommended.
        Default is 7.

    Returns
    -------
    dat_smt2 : numpy.ndarray
        Smoothed Y values with Gaussian weighting.

    Notes
    -----
    For additional smoothing ideas and references, see:
    - https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html
    - https://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way
    - https://stackoverflow.com/questions/32900854/how-to-smooth-a-line-using-gaussian-kde-kernel-in-python-setting-a-bandwidth

    """

    tim_raw = tim_ref
    dat_raw = dat_ref

    dat_smt2 = []

    num_raw = len(tim_raw)

    for ismt in range(num_raw):

        tim_raw_work = np.delete(tim_raw, ismt)
        dat_raw_work = np.delete(dat_raw, ismt)

        x_lag = tim_raw[ismt] - tim_raw_work
        x_fac = np.exp(-((x_lag / width) ** 2) / 2)

        # x_fac[x_fac < 0.01] = 0.

        clean_bool = x_fac > 0.01
        ## It differs a bit of the official fct, because the next element
        ## following this criteria is included anyway

        dat_raw_clean = dat_raw_work[clean_bool]
        x_fac_clean = x_fac[clean_bool]

        dat_raw_clean = dat_raw_clean[~np.isnan(dat_raw_clean)]
        x_fac_clean = x_fac_clean[~np.isnan(x_fac_clean)]

        x_val = np.sum(np.multiply(dat_raw_clean, x_fac_clean))
        x_wht = np.sum(x_fac_clean)

        dat_smt2.append(x_val / x_wht)

    dat_smt2 = np.array(dat_smt2)

    return dat_smt2

def smooth(x, window_len=11, window="hanning"):
    """
    Smooth data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal. 
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends to minimize transient parts at the 
    beginning and end of the output signal.

    Parameters
    ----------
    x : numpy.ndarray
        The input signal. Must be 1-dimensional.
    window_len : int, optional
        The dimension of the smoothing window. Should be an odd integer.
        Default is 11.
    window : str, optional
        The type of window from {'flat', 'hanning', 'hamming', 'bartlett', 'blackman'}.
        'flat' window will produce a moving average smoothing.
        Default is 'hanning'.

    Returns
    -------
    y : numpy.ndarray
        The smoothed signal.

    Raises
    ------
    ValueError
        If input is not 1-dimensional or smaller than window size.

    Notes
    -----
    Works only for equally spaced data.
    
    The output length may differ from input length. To correct this, 
    return `y[(window_len/2-1):-(window_len/2)]` instead of `y`.

    See Also
    --------
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    References
    ----------
    https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html

    Examples
    --------
    >>> t = np.linspace(-2, 2, 0.1)
    >>> x = np.sin(t) + np.random.randn(len(t)) * 0.1
    >>> y = smooth(x)

    """

    if x.ndim != 1:
        raise (ValueError, "smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise (ValueError, "Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if not window in ["flat", "hanning", "hamming", "bartlett", "blackman"]:
        raise (
            ValueError,
            "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'",
        )

    s = np.r_[x[window_len - 1 : 0 : -1], x, x[-2 : -window_len - 1 : -1]]
    if window == "flat":  # moving average
        w = np.ones(window_len, "d")
    else:
        w = eval("numpy." + window + "(window_len)")

    y = np.convolve(w / w.sum(), s, mode="valid")
    return y


def harmonic_mean(a):
    """
    Compute harmonic mean of a list or array.

    Parameters
    ----------
    a : array-like
        Input values.

    Returns
    -------
    result : float
        The harmonic mean of the input values.

    """

    a = np.array(a)
    return len(a) / np.sum(1.0 / a)


def find_intersection(x1, y1, x2, y2):
    """
    Find intersection points of two line plots.

    Parameters
    ----------
    x1 : array-like
        X-coordinates of the first line.
    y1 : array-like
        Y-coordinates of the first line.
    x2 : array-like
        X-coordinates of the second line.
    y2 : array-like
        Y-coordinates of the second line.

    Returns
    -------
    roots : numpy.ndarray
        X-coordinates of the intersection points.
    y_intersect : numpy.ndarray
        Y-coordinates of the intersection points.

    References
    ----------
    http://stackoverflow.com/questions/8094374/python-matplotlib-find-intersection-of-lineplots

    """

    import scipy.interpolate as interpolate
    import scipy.optimize as optimize

    p1 = interpolate.PiecewisePolynomial(x1, y1[:, np.newaxis])
    p2 = interpolate.PiecewisePolynomial(x2, y2[:, np.newaxis])

    def pdiff(x):
        return p1(x) - p2(x)

    xs = np.r_[x1, x2]
    xs.sort()
    x_min = xs.min()
    x_max = xs.max()
    x_mid = xs[:-1] + np.diff(xs) / 2
    roots = set()
    for val in x_mid:
        root, infodict, ier, mesg = optimize.fsolve(pdiff, val, full_output=True)
        # ier==1 indicates a root has been found
        if ier == 1 and x_min < root < x_max:
            roots.add(root[0])
    roots = np.array(list(roots))
    return roots, p1(roots)


def wrap_to360(lon):
    """
    Wrap longitude values to the range [0, 360).

    Parameters
    ----------
    lon : float or array-like
        Longitude values in degrees.

    Returns
    -------
    result : float or numpy.ndarray
        Longitude values wrapped to [0, 360].

    Notes
    -----
    Based on MATLAB's wrap_to360 function.

    """
    lon = np.mod(lon, 360)
    return lon


def wrap_to180(lon):
    """
    Wrap longitude values to the range (-180, 180].

    Parameters
    ----------
    lon : float or array-like
        Longitude values in degrees.

    Returns
    -------
    result : float or numpy.ndarray
        Longitude values wrapped to (-180, 180].

    Notes
    -----
    Based on MATLAB's wrap_to180 function.

    """
    if not (lon is np.array):
        notaarray = True
        lon = np.array([lon])
    else:
        notaarray = False
    q = (lon < -180) + (180 < lon)
    lon[q] = wrap_to360(lon[q] + 180) - 180
    if notaarray:
        lon = lon[0]
    return lon


# Low level statistic function


def rms_mean(a):
    """
    Compute RMS (Root Mean Square) of a list or array.

    Parameters
    ----------
    a : array-like
        Input values.

    Returns
    -------
    result : float
        The RMS value of the input.

    """
    return np.sqrt(np.nanmean(np.square(np.array(a, np.float64))))

def rms_mean_alternativ(a):
    """
    Compute RMS with mean subtraction (standard deviation equivalent).

    Computes: √< (A - Ā)² > instead of √< (A)² >, where Ā is the arithmetic mean.
    This is essentially the standard deviation.

    Parameters
    ----------
    a : array-like
        Input values.

    Returns
    -------
    result : float
        The RMS value with mean subtraction.

    Notes
    -----
    This is equivalent to the standard deviation of the input.

    """
    return np.sqrt(np.nanmean(np.square(a - np.nanmean(a))))


def rms_mean_kouba(a, multipl_coef=3, deg_of_freedom=7):
    """
    Compute weighted RMS (Root Mean Square) with Kouba's method.

    Parameters
    ----------
    a : array-like
        Input values.
    multipl_coef : int, optional
        Multiplication coefficient. Default is 3.
    deg_of_freedom : int, optional
        Degrees of freedom. Default is 7.

    Returns
    -------
    result : float
        The weighted RMS value.

    """
    return np.sqrt(np.sum(np.square(a)) / (multipl_coef * len(a) - deg_of_freedom))


def mad(data, mode="median"):
    """
    Compute Median Absolute Deviation (MAD).

    Parameters
    ----------
    data : array-like
        Input data values.
    mode : str, optional
        Mode for computing deviation center: 'median' or 'mean'.
        Default is 'median'.

    Returns
    -------
    result : float
        The Median Absolute Deviation.

    """

    if mode == "median":
        mad_out = np.nanmedian(np.abs(data - np.nanmedian(data)))
    elif mode == "mean":
        mad_out = np.nanmean(np.abs(data - np.nanmean(data)))

    return mad_out


def outlier_mad(
    data,
    threshold=3.5,
    verbose=False,
    convert_to_np_array=True,
    mad_mode="median",
    seuil=None,
):
    """
    clean the outlier of Ya dataset using the MAD approach

    Parameters
    ----------
    data : list or numpy.array
        Values

    threshold : float
        MAD threshold

    verbose : bool

    convert_to_np_array : bool
        if True returns output as an array, if False as a regular list

    mad_mode : str
        'median' or 'mean' : MAD can also be based on mean (for experimental purposes)

    seuil : float, optional
        legacy name of 'threshold' argument.
        will override threshold value if given


    Returns
    -------
    dataout : numpy.array
        Values cleaned of outliers

    boolbad : numpy.array
        Y-sized booleans

    Source
    ------
    Utilisation de la MAD pour detecter les outliers
    https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm
    https://web.ipac.caltech.edu/staff/fmasci/home/statistics_refs/BetterThanMAD.pdf
    """

    if seuil:
        log.warning("'seuil' argument for outlier_mad is deprecated !!!")
        threshold = seuil

    if convert_to_np_array:
        data = np.array(data)
    nbinp = float(len(data))
    mad_value = mad(data, mode=mad_mode)
    med = np.nanmedian(data)

    if np.isclose(np.sum(np.abs(np.diff(data))), 0.0):
        if verbose:
            log.info("elimination ratio: 0 , all values are equals")
        return data, np.array([True] * len(data))

    if np.isclose(med, 0.0):
        if verbose:
            log.info("elimination ratio: 0 , null median")
        return data, np.array([True] * len(data))

    if np.isclose(mad_value, 0.0):
        if verbose:
            log.info("elimination ratio: 0 , null MAD")
        return data, np.array([True] * len(data))

    diff = data - med
    mzs = 0.6745 * np.abs(diff) / mad_value

    mzs[np.isnan(mzs)] = threshold * 10  ### makes a threshold virtually higher for NaN
    boolbad = mzs < threshold  ## remendir: False means that the values are bad.

    dataout = data[boolbad]
    nbout = float(sum(boolbad))
    ratio = (nbinp - nbout) / nbinp
    if verbose:
        log.info(
            "MAD outlier elimination ratio: %i / %i, %f p.c.",
            nbinp - nbout,
            nbinp,
            ratio * 100,
        )
    return dataout, boolbad


def outlier_mad_binom(
    y, x, threshold=3.5, verbose=False, detrend_first=False, return_booleans=False
):
    """
    clean the outlier of Y using the MAD approach
    and clean the corresponding values in X
    assuming that we have the function : X => Y(X)
    (be carefull, Y is the first argument)

    Parameters
    ----------
    y : list or numpy.array
        Values

    x : list or numpy.array
        X Values so as X => Y(X)

    threshold : float
        MAD threshold

    verbose : bool

    detrend_first : bool
        detrend linear behavior of Y(X) first

    return_booleans : bool
        return good and bad values of Y and X as booleans

    Returns
    -------
    yclean & xclean : numpy.array

    bb : numpy.array (if return_booleans == True)
        Y-sized booleans
    """
    if detrend_first:
        _, y_work = detrend_timeseries(x, y)
    else:
        _, y_work = np.array(x), np.array(y)

    _, bb = outlier_mad(y_work, threshold, verbose)

    x_clean = np.array(x)[bb]
    y_clean = np.array(y)[bb]

    if not return_booleans:
        return y_clean, x_clean
    else:
        return y_clean, x_clean, bb


def outlier_above_below_simple(x, low_bound, upp_bound, return_booleans=True):
    """
    Gives values of X which are between low_bound & upp_bound

    Parameters
    ----------
    x : list or numpy.array
        Values

    low_bound & upp_bound  : float
        lower and upper bound of X values wished

    return_booleans : bool
        return booleans or not

    Returns
    -------
    xout : numpy.array
        X between low_bound & upp_bound

    bbool : bool
         X-sized array of booleans

    """

    xwork = np.array(x)

    if low_bound >= upp_bound:
        log.warning("lower bound >= upper bound !!!")
        log.warning("low_bond : ", low_bound)
        log.warning("upp_bond : ", upp_bound)

    bbool = (low_bound <= xwork) & (xwork <= upp_bound)

    xout = xwork[bbool]

    if return_booleans:
        return xout, bbool
    else:
        return xout


def outlier_above_below(
    x,
    threshold_values,
    reference=np.nanmean,
    theshold_absolute=True,
    return_booleans=True,
    theshold_relative_value="reference",
    verbose=False,
):
    """
    Gives values of X which are between threshold values

    Parameters
    ----------
    threshold_values : single value (float) or a 2-tuple
        (lower bound theshold , upper bound theshold)

        `WARN` : those value(s) have to be positives.
        Minus sign for lower bound and plus sign for upper
        one will be applied internally

    reference : float or callable
        the central reference value
        can be a absolute fixed value (float) or
        a function (e.g. np.mean of np.median)

    theshold_absolute : bool
        if True threshold_values are absolutes values
            >>> low = reference - threshold_values[0]
            >>> upp = reference + threshold_values[1]
        if False they are fractions of theshold_relative_value
            >>> low = reference - threshold_values[0] * theshold_relative_value
            >>> upp = reference + threshold_values[1] * theshold_relative_value
        (see also below)

    theshold_relative_value : str or function
        if the string "reference" or None is given, then it the reference
        value which is used
        if it is a fuction (e.g. np.std()) then it is this value returned
        by this function which is used
        Only useful when theshold_absolute = False

    return_booleans : bool
        return booleans or not

    verbose : bool

    Returns
    -------
    xout : numpy array
        X between low_bound & upp_bound

    bbool : numpy array
        X-sized array of booleans
    """

    if utils.is_iterable(threshold_values):
        ths_input_low = threshold_values[0]
        ths_input_upp = threshold_values[1]
    else:
        ths_input_low = threshold_values
        ths_input_upp = threshold_values

    if ths_input_low < 0.0 or ths_input_upp < 0.0:
        log.warning("threshold_values have to be positive")
        log.warning("minus sign for lower bound will be applied internally")

    if callable(reference):
        ref_val = reference(x)
    else:
        ref_val = reference

    if theshold_relative_value in ("reference", None):
        relativ_val = reference
    elif callable(theshold_relative_value):
        relativ_val = theshold_relative_value(x)
    else:
        relativ_val = reference

    if theshold_absolute:
        ths_low = ref_val - ths_input_low
        ths_upp = ref_val + ths_input_upp
    else:
        ths_low = ref_val - ths_input_low * relativ_val
        ths_upp = ref_val + ths_input_upp * relativ_val

    if verbose:
        log.info("outlier_above_below theshold values")
        log.info("reference: %s", ref_val)
        log.info("effective lower bound: %s", ths_low)
        log.info("effective upper bound: %s", ths_upp)

    xout, bbool = outlier_above_below_simple(x, ths_low, ths_upp)

    if return_booleans:
        return xout, bbool
    else:
        return xout


def outlier_above_below_binom(
    y,
    x,
    threshold_values,
    reference=np.nanmean,
    theshold_absolute=True,
    theshold_relative_value="reference",
    return_booleans=False,
    detrend_first=True,
    verbose=False,
):
    """
    Gives values of Y which are between threshold values, and correct an
    associated X so as X => Y(X)

    Parameters
    ----------
    threshold_values : single value (float) or a 2-tuple
        (lower bound theshold , upper bound theshold)

        `WARN` : those value(s) have to be positives.
        Minus sign for lower bound and plus sign for upper
        one will be applied internally

    reference : float or callable
        the central reference value
        can be a absolute fixed value (float) or
        a function (e.g. np.mean of np.median)

    theshold_absolute : bool
        if True threshold_values are absolutes values
            >>> low = reference - threshold_values[0]
            >>> upp = reference + threshold_values[1]
        if False they are fractions of theshold_relative_value
            >>> low = reference - threshold_values[0] * theshold_relative_value
            >>> upp = reference + threshold_values[1] * theshold_relative_value
        (see also below)

    theshold_relative_value : str or function
        if the string "reference" or None is given, then it the reference
        value which is used
        if it is a fuction (e.g. np.std()) then it is this value returned
        by this function which is used
        Only useful when theshold_absolute = False

    detrend_first : bool
        detrend linear behavior of Y(X) first
        Recommended

    return_booleans : bool
        return booleans or not

    verbose : bool


    Returns
    -------
    x_out : numpy array
        X between low_bound & upp_bound

    bbool : numpy array
        X-sized array of booleans
    """

    if detrend_first:
        _, y_work = detrend_timeseries(x, y)
    else:
        _, y_work = np.array(x), np.array(y)

    _, bb = outlier_above_below(
        y_work,
        threshold_values,
        reference=reference,
        theshold_absolute=theshold_absolute,
        theshold_relative_value=theshold_relative_value,
        return_booleans=True,
        verbose=verbose,
    )

    x_clean = np.array(x)[bb]
    y_clean = np.array(y)[bb]

    if not return_booleans:
        return y_clean, x_clean
    else:
        return y_clean, x_clean, bb


def outlier_sigma(datasigmain, threshold=3):
    """
    Remove outliers based on sigma threshold (legacy method).

    Removes points where sigma > threshold * median(sigmas). This is an old 
    and discontinued method that is not very efficient.

    Parameters
    ----------
    datasigmain : array-like
        Input sigma (uncertainty) values.
    threshold : int, optional
        Multiplier for the sigma threshold. Default is 3.

    Returns
    -------
    datasigmaout : numpy.ndarray
        Sigma values filtered to exclude outliers.
    boolbad : numpy.ndarray
        Boolean array indicating valid values (True = keep, False = outlier).

    Notes
    -----
    This is a legacy function and is rarely used. More modern methods like 
    MAD-based outlier detection are recommended.

    """
    moy = np.median(datasigmain)
    marge = moy * threshold

    log.info("moy,threshold,margin: %f, %i, %f", moy, threshold, marge)

    boolbad = np.abs(datasigmain) < marge

    datasigmaout = datasigmain[boolbad]

    return datasigmaout, boolbad


def lagrange1(points):
    """
    Determine Lagrangian polynomial from points (low-level function).

    Creates a polynomial interpolation function using Lagrange's method. 
    This replaces scipy.interpolate.lagrange which is highly unstable.

    Parameters
    ----------
    points : list of tuple
        List of (x, y) coordinate tuples defining the polynomial.

    Returns
    -------
    p : callable
        Function representing the Lagrangian polynomial.

    Notes
    -----
    More numerically stable than scipy.interpolate.lagrange.

    References
    ----------
    https://gist.github.com/melpomene/2482930

    """

    def poly(x):
        total = 0
        n = len(points)
        for i in range(n):
            xi, yi = points[i]

            def g(ii, nn):

                tot_mul = 1
                for jj in range(nn):
                    if ii == jj:
                        continue
                    xj, yj = points[jj]
                    tot_mul *= (x - xj) / float(xi - xj)

                return tot_mul

            total += yi * g(i, n)
        return total

    return poly


def lagrange2(x, y):
    """
    Determine Lagrangian polynomial from points (more Pythonic version).

    Creates a polynomial interpolation function using Lagrange's method. 
    This version is more Pythonic but slower than lagrange1. Like lagrange1, 
    it replaces scipy.interpolate.lagrange which is highly unstable.

    Parameters
    ----------
    x : array-like
        X-coordinates of the interpolation points.
    y : array-like
        Y-coordinates of the interpolation points.

    Returns
    -------
    p : callable
        Function representing the Lagrangian polynomial.

    Notes
    -----
    More numerically stable than scipy.interpolate.lagrange.
    More Pythonic but slower than lagrange1.

    References
    ----------
    https://gist.github.com/melpomene/2482930

    """

    def poly(x_itrp):
        total = 0
        n = len(x)
        for i in range(n):

            def g(ii, nn):
                x_but_i = np.concatenate((x[:ii], x[ii + 1:]))
                # mask = np.ones(len(x),dtype=bool)
                # mask[i] = False
                # x_but_i = np.concatenate((x[:i],x[i+1:]))
                # x_but_i = x[mask]
                # return np.prod((x_itrp -x[mask])/(x[i] - x[mask]))

                return np.prod((x_itrp - x_but_i) / (x[ii] - x_but_i))

            total += y[i] * g(i, n)
        return total

    return poly


def lagrange_interpolate(tdata, ydata, titrp, n=10, t_type="datetime"):
    """
    Perform temporal Lagrangian polynomial interpolation.

    Interpolates Y values at specified time epochs using Lagrangian polynomial
    fitting. The X-component represents time.

    Parameters
    ----------
    tdata : array-like of datetime
        X/T component (time) of the known interpolation points.
    ydata : array-like of float
        Y component (data values) of the known interpolation points.
    titrp : array-like of datetime
        Epochs at which to compute interpolated values.
    n : int, optional
        Degree of the Lagrangian polynomial. Better if even. Default is 10.
    t_type : str, optional
        type of the time component, can be "datetime", "posix" or "pandas_timestamp".
        The default is "datetime".
        pandas_timestamp is recommended for a more precise applications
        (nanosecond precision instead of microsecond for datetime)

    Returns
    -------
    y_intrp : numpy.ndarray
        Interpolated Y values at the requested epochs.

    See Also
    --------
    lagrange1 : Low-level Lagrangian polynomial function
    conv.dt_range : Generate a range of datetime epochs

    Notes
    -----
    Use conv.dt_range to generate the wished epochs range.
    """
    ydata = np.array(ydata)

    nn = int(n / 2)

    if t_type == "datetime":
        tdata_px = conv.dt2posix(np.array(tdata))
        titrp_px = conv.dt2posix(np.array(titrp))
    elif t_type == "posix":
        tdata_px = np.array(tdata)
        titrp_px = np.array(titrp)
    elif t_type == "pandas_timestamp":
        ## here we work in int64 and nanosecond precision, for better precision and stability
        tdata_px = conv.pandas_timestamp2posix(tdata)
        titrp_px = conv.pandas_timestamp2posix(titrp)
    else:
        log.error("t_type should be 'datetime', 'posix' or 'pandas_timestamp'")
        raise ValueError("t_type should be 'datetime', 'posix' or 'pandas_timestamp'")

    tref = tdata_px[0]

    ### we substract a ref time to avoid numerical instability
    tdata_px = tdata_px - tref
    titrp_px = titrp_px - tref

    sur_val = (np.nan, np.nan)
    sur_idx = (np.nan, np.nan)

    ### some checks
    if np.any(np.diff(tdata_px) == 0):
        log.warning("some Tdata are equals")

    if np.any(np.diff(ydata) == 0):
        log.warning("some Ydata are equals")

    if np.any(titrp_px < 0):
        log.warning("some wanted values are outside the data interval!!!!")

    yintrp = []

    for tintrp in titrp_px:

        if (sur_val[0] <= tintrp) & (tintrp <= sur_val[1]):
            ### the Polynom is alread determined
            pass
        else:
            sur_val, sur_idx = utils.find_surrounding(tdata_px, tintrp)

            if sur_idx[0] - nn < 0:  # manage side effect for first points
                imin = 0
                imax = n + 1
            elif sur_idx[1] + nn > len(ydata):  # manage side effect for last points
                imin = len(ydata) - n - 1
                imax = len(ydata)
            else:  # regular case
                ### if (sur_idx[0] - nn >= 0) and (sur_idx[1] + nn >= len(Ydata)):
                imin = sur_idx[0] - nn
                imax = sur_idx[1] + nn

            tuse = tdata_px[imin:imax]
            yuse = ydata[imin:imax]

            poly = lagrange1(list(zip(tuse, yuse)))
            # poly = lagrange2(tuse,yuse)

        yintrp = poly(tintrp)
        yintrp.append(yintrp)

    return np.array(yintrp)

def dates_middle(start, end):
    """
    Compute the midpoint between two dates.

    Parameters
    ----------
    start : datetime or numeric
        Start date/time.
    end : datetime or numeric
        End date/time.

    Returns
    -------
    middle : datetime or numeric
        The midpoint between start and end.

    """
    return start + (end - start) / 2


def time_win_basic(
    start,
    end,
    t_lis_inp,
    data_lis_inp,
    outposix=True,
    invert=False,
    out_array=False,
    out_boolis=False,
    only_boolis=False,
):
    """
    Filter data within a time window.

    Selects data points that fall within a specified time window. 
    Internally converts to POSIX time for computation.

    Parameters
    ----------
    start : datetime or float
        Start of the time window.
    end : datetime or float
        End of the time window.
    t_lis_inp : array-like of datetime or float
        Time values of the data points.
    data_lis_inp : array-like
        Data values corresponding to the times.
    outposix : bool, optional
        If True, output times in POSIX format.
        If False, output as datetime. Default is True.
    invert : bool, optional
        If True, invert the boolean selection (select data outside window).
        Default is False.
    out_array : bool, optional
        If True, return outputs as numpy arrays. Default is False.
    out_boolis : bool, optional
        If True, also return the boolean selection array. Default is False.
    only_boolis : bool, optional
        If True, skip filtering and only return boolean array.
        Default is False.

    Returns
    -------
    tlisout : array-like
        Filtered time values. None if only_boolis=True.
    datalisout : array-like
        Filtered data values. None if only_boolis=True.
    boolis : numpy.ndarray, optional
        Boolean array indicating selected points.
        Only returned if out_boolis=True.

    """

    if isinstance(t_lis_inp[0], dt.datetime):
        tlis = conv.dt2posix(t_lis_inp)
    else:
        tlis = t_lis_inp

    if isinstance(start, dt.datetime):
        start = conv.dt2posix(start)
    if isinstance(end, dt.datetime):
        end = conv.dt2posix(end)

    if not isinstance(tlis, np.ndarray) or not isinstance(data_lis_inp, np.ndarray):
        tlis = np.array(tlis)
        datalis = np.array(data_lis_inp)
    else:
        tlis = tlis
        datalis = data_lis_inp

    boolis = (start <= tlis) * (tlis <= end)

    if invert:
        boolis = np.logical_not(boolis)

    if only_boolis:
        datalisout = None
        tlisout = None
    else:
        datalisout = datalis[boolis]
        tlisout = tlis[boolis]

        if not outposix:
            tlisout = conv.posix2dt(tlisout)

        if out_array:
            tlisout, datalisout = np.array(tlisout), np.array(datalisout)

    if out_boolis:
        out_tuple = (tlisout, datalisout, boolis)
    else:
        out_tuple = (tlisout, datalisout)

    return out_tuple


def time_win_multi(
    start, end, t_lis, data_lislis, outposix=True, invert=False, out_array=False
):
    """
    Filter multiple datasets within a time window.

    Applies time window filtering to multiple data arrays simultaneously, 
    using the same time array for all datasets.

    Parameters
    ----------
    start : datetime or float
        Start of the time window.
    end : datetime or float
        End of the time window.
    t_lis : array-like
        Time values for filtering.
    data_lislis : list of array-like
        Multiple data arrays to filter.
    outposix : bool, optional
        If True, output times in POSIX format. Default is True.
    invert : bool, optional
        If True, invert the selection. Default is False.
    out_array : bool, optional
        If True, return outputs as numpy arrays. Default is False.

    Returns
    -------
    Tlisout : array-like
        Filtered time values.
    datalislisout : list of array-like
        Filtered data arrays.

    See Also
    --------
    time_win_basic : Filter single dataset within time window

    """
    datalislisout = []
    for i, datalis in enumerate(data_lislis):

        Tlisout, datalisout = time_win_basic(
            start, end, t_lis, datalis, outposix, invert, out_array=out_array
        )
        datalislisout.append(datalisout)
    return Tlisout, datalislisout


def time_win_multi_start_end(
    start_list_inp,
    end_list_inp,
    t_lis_inp,
    data_lis_inp,
    outposix=True,
    invert=False,
    out_array=False,
    out_boolis=False,
):
    """
    Filter data within multiple time windows simultaneously.

    Selects data points that fall within ALL specified time windows. 
    This is useful for selecting intersections of multiple time periods.
    Internally converts to POSIX time for computation.

    Parameters
    ----------
    start_list_inp : list of datetime or float
        Start times of the time windows.
    end_list_inp : list of datetime or float
        End times of the time windows.
    t_lis_inp : array-like of datetime or float
        Time values of the data points.
    data_lis_inp : array-like
        Data values corresponding to the times.
    outposix : bool, optional
        If True, output times in POSIX format. Default is True.
    invert : bool, optional
        If True, invert the selection. Default is False.
    out_array : bool, optional
        If True, return outputs as numpy arrays. Default is False.
    out_boolis : bool, optional
        If True, also return boolean selection arrays. Default is False.

    Returns
    -------
    tlisout : array-like
        Filtered time values.
    datalisout : array-like
        Filtered data values.
    boolis_opera : numpy.ndarray, optional
        Combined boolean array (intersection of all windows).
        Only returned if out_boolis=True.
    boolis_stk : numpy.ndarray, optional
        Stack of individual boolean arrays for each window.
        Only returned if out_boolis=True.

    Raises
    ------
    ValueError
        If len(Start_list_in) != len(End_list_in).

    """

    if len(start_list_inp) != len(end_list_inp):
        log.error("len(Start_list_in) != len(End_list_in) !!")

    boolis_stk = []
    for start, end in zip(start_list_inp, end_list_inp):
        _, _, boolis = time_win_basic(
            start,
            end,
            t_lis_inp,
            data_lis_inp,
            outposix=outposix,
            invert=invert,
            out_boolis=True,
            only_boolis=True,
        )

        boolis_stk.append(boolis)

    boolis_stk = np.stack(boolis_stk)
    boolis_opera = np.all(boolis_stk, axis=0)

    datalis = np.array(data_lis_inp)
    tlis = np.array(t_lis_inp)

    datalisout = datalis[boolis_opera]
    tlisout = tlis[boolis_opera]

    if not outposix:
        tlisout = conv.posix2dt(tlisout)

    if out_array:
        tlisout, datalisout = np.array(tlisout), np.array(datalisout)

    if out_boolis:
        out_tuple = (tlisout, datalisout, boolis_opera, boolis_stk)
    else:
        out_tuple = (tlisout, datalisout)

    return out_tuple


def get_season(now):
    """
    Determine the season for a given date.

    Parameters
    ----------
    now : datetime.date or datetime.datetime
        The date to determine the season for. If datetime is provided, 
        the date component is used.

    Returns
    -------
    season : str
        The season name: 'winter', 'spring', 'summer', or 'autumn'.
        Returns None if the date is outside the defined ranges.

    Notes
    -----
    Season boundaries are:
    - Winter: Dec 21 - Mar 20
    - Spring: Mar 21 - Jun 20
    - Summer: Jun 21 - Sep 22
    - Autumn: Sep 23 - Dec 20

    """

    seasons = [
        ("winter", (dt.date(1, 1, 1), dt.date(1, 3, 20))),
        ("spring", (dt.date(1, 3, 21), dt.date(1, 6, 20))),
        ("summer", (dt.date(1, 6, 21), dt.date(1, 9, 22))),
        ("autumn", (dt.date(1, 9, 23), dt.date(1, 12, 20))),
        ("winter", (dt.date(1, 12, 21), dt.date(1, 12, 31))),
    ]

    # suppressing the year
    if isinstance(now, dt.datetime):
        now = now.date()
    now = now.replace(year=1)

    for season, (start, end) in seasons:
        if start <= now <= end:
            return season


def color_of_season(datein):
    """
    Get color representation for a season.

    Maps seasons to matplotlib-compatible color codes.

    Parameters
    ----------
    datein : datetime.date or datetime.datetime
        The date to get the season color for.

    Returns
    -------
    color : str
        Color code: 'b' (blue) for winter, 'r' (red) for summer,
        'g' (green) for spring, 'k' (black) for autumn.

    See Also
    --------
    get_season : Determine season from date

    """
    season = get_season(datein)
    if season == "winter":
        outcolor = "b"
    elif season == "summer":
        outcolor = "r"
    elif season == "spring":
        outcolor = "g"
    elif season == "autumn":
        outcolor = "k"
    return outcolor


#  ______                _   _                _____                                         _
# |  ____|              | | (_)              / ____|                                       | |
# | |__ _   _ _ __   ___| |_ _  ___  _ __   | |  __ _ __ __ ___   _____ _   _  __ _ _ __ __| |
# |  __| | | | '_ \ / __| __| |/ _ \| '_ \  | | |_ | '__/ _` \ \ / / _ \ | | |/ _` | '__/ _` |
# | |  | |_| | | | | (__| |_| | (_) | | | | | |__| | | | (_| |\ v /  __/ |_| | (_| | | | (_| |
# |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|  \_____|_|  \__,_| \_/ \___|\__, |\__,_|_|  \__,_|
#                                                                        __/ |
#                                                                       |___/


def outlier_mad_binom_legacy(
    x, y, threshold=3.5, verbose=False, detrend_first=False, return_booleans=False
):
    """
    Remove outliers from paired X,Y data using MAD (legacy version).

    Legacy version with different argument order than the main version.
    May be unstable when detrending.

    Parameters
    ----------
    x : array-like
        X values.
    y : array-like
        Y values (dependent on X).
    threshold : float, optional
        MAD threshold. Default is 3.5.
    verbose : bool, optional
        If True, print elimination information. Default is False.
    detrend_first : bool, optional
        If True, detrend before outlier detection. Default is False.
    return_booleans : bool, optional
        If True, also return boolean selection array. Default is False.

    Returns
    -------
    x_clean : numpy.ndarray
        Cleaned X values.
    y_clean : numpy.ndarray
        Cleaned Y values.
    bb : numpy.ndarray, optional
        Boolean array indicating valid points.
        Only returned if return_booleans=True.

    Notes
    -----
    This is a legacy function. Use outlier_mad_binom for the standard version.

    See Also
    --------
    outlier_mad_binom : Main version with correct argument order

    """
    if detrend_first:
        x_work, _ = detrend_timeseries(x, y)
    else:
        x_work, _ = np.array(x), np.array(y)

    _, bb = outlier_mad(x_work, threshold, verbose)

    x_clean = np.array(x)[bb]
    y_clean = np.array(y)[bb]

    if not return_booleans:
        return x_clean, y_clean
    else:
        return x_clean, y_clean, bb


# def plot_vertical_bar(xlis , color='r',linewidth=1):
#     out_bar_list = []
#     for x in xlis:
#         out_bar = plt.axvline(x,color=color,linewidth=linewidth)
#         out_bar_list.append(out_bar)
#     return out_bar_list
#
# def plot_vertical_bar_ax(xlis,ax_in,color='r',linewidth=1):
#     out_bar_list = []
#     for x in xlis:
#         out_bar = ax_in.axvline(x,color=color,linewidth=linewidth)
#         out_bar_list.append(out_bar)
#     return out_bar_list



def gaussian_filter_gfz_legacy(tim_ref, dat_ref, width=7):
    """
    Apply Gaussian filter to smooth data (legacy, slow version).

    Gaussian filter based on GFZ's GMT_plus.pm/gaussian_kernel. 
    This is a legacy version that is VERY SLOW due to a dirty conversion 
    from Perl code. The pythonic version gaussian_filter_gfz should be 
    used instead.

    Parameters
    ----------
    tim_ref : array-like
        X/T component of the time series (in decimal days).
    dat_ref : array-like
        Y component (the data values).
    width : int, optional
        Size of the window. Odd numbers are recommended. Default is 7.

    Returns
    -------
    dat_smt : list
        Smoothed Y values.

    Warnings
    --------
    This function is VERY SLOW. Use gaussian_filter_gfz instead.

    Notes
    -----
    For additional smoothing ideas and references, see:
    - http://scipy-cookbook.readthedocs.io/items/SignalSmooth.html
    - https://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way
    - https://stackoverflow.com/questions/32900854/how-to-smooth-a-line-using-gaussian-kde-kernel-in-python-setting-a-bandwidth

    See Also
    --------
    gaussian_filter_gfz : Pythonic and faster version

    """

    log.warning("WARN : THIS function gaussian_filter_GFZ_style_smoother")
    log.warning("IS VERY SLOW (DIRTY CONVERSION OF A PERL FCT)")
    log.warning("THE PYTHONIC VERSION gaussian_filter_GFZ_style_smoother_improved")
    log.warning("SHOULD BE USED  !!!")

    tim_raw = tim_ref
    dat_raw = dat_ref

    num_raw = len(tim_raw)
    icomp = 0

    dat_smt = [np.nan] * len(dat_ref)

    for ismt in range(num_raw):  # +1

        x_val = 0.0
        x_wht = 0.0

        for iraw in reversed(range(0, ismt)):
            x_lag = tim_raw[ismt] - tim_raw[iraw]
            x_fac = np.exp(-((x_lag / width) ** 2) / 2)
            x_val += dat_raw[iraw] * x_fac
            x_wht += x_fac
            icomp += 1
            if x_fac < 0.01:
                break

        for iraw in range(ismt + 1, num_raw):
            x_lag = tim_raw[ismt] - tim_raw[iraw]
            x_fac = np.exp(-((x_lag / width) ** 2) / 2)
            x_val += dat_raw[iraw] * x_fac
            x_wht += x_fac
            icomp += 1
            if x_fac < 0.01:
                break

        dat_smt[ismt] = x_val / x_wht

    return dat_smt
