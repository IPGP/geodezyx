# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.stats contains functions for low-level
statistics.

it can be imported directly with:
from geodezyx import stats

The geodezyx toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/IPGP/geodezyx
"""

########## BEGIN IMPORT ##########
#### External modules
import datetime as dt

#### Import the logger
import logging

import numpy as np
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
    std_err = None
    if simple_lsq:
        A = np.array(
            [x, np.ones(len(x))]
        )  # x2 faster than np.column_stack([x, np.ones(len(x))])
        slope, intercept = np.linalg.lstsq(A.T, y, rcond=None)[
            0
        ]  # obtaining the parameters
    else:
        # Perform scipy's linregress if simple_lsq is False
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)

    # Return only the slope and intercept if fulloutput is False
    if simple_lsq or not fulloutput:
        return slope, intercept
    else:
        # Return the slope, intercept, confidence interval, and standard deviation if fulloutput is True
        return slope, intercept, confid_interval_slope(x, y, alpha), std_err


def linear_reg_getvalue(X, a, b, full=True):
    """
    From a vector X and coefficients a & b, get Y = a*X + b

    Parameters
    ----------
    X : list or numpy.array
        Values

    a & b : float
        Linear regression coefficients

    full : bool
        True : return X , Y = aX + b , False : return Y = aX + b

    Returns
    -------
    Y : numpy.array
        if full == False

    OR

    X , Y : numpy.array
        if full == True

    Note
    ----
        Unstable while working with POSIX Time as X-data (too heigh values ? ...)
        Decimal Years are recommended

    """

    if full:
        return np.array(X), a * np.array(X) + b
    else:
        return a * np.array(X) + b


def linear_coef_a_b(x1, y1, x2, y2):
    """
    Gives coefficients of the line between two points (x1,y1) & (x2,y2)
    x1,y1,x2,y2 can be iterables

    Parameters
    ----------
    x1,y1,x2,y2 : float or list or numpy.array
        Coordinates of the 1st and the 2nd point

    Returns
    -------
    a : float
        regression coefficient

    b1 & b2 : float
        regression offsets coefficient (b1 must be equal to b2)

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


def detrend_timeseries(X, Y):
    """
    detrend, i.e. remove linear tendence of a timeserie Y(X)

    Parameters
    ----------
    X & Y: list or numpy.array
        Values

    Returns
    -------
    X & Yout: list or numpy.array
        Detrended Y

    """

    X = np.array(X)
    Y = np.array(Y)
    a, b = linear_regression(X, Y)

    # Yout = Y - a * (X - X[0])
    # Yout = Y - ( a * X + b )

    Ylinear = a * X + b
    Yout = Y - Ylinear + Y[0]

    return X, Yout


def confid_interval_slope(x, y, alpha=0.95):
    """
    Calcule un intervalle de confiance sur une tendance

    Parameters
    ----------
    x : array
        la variable indépendante
    y : array
        la variable dépendante
    alpha : float
        la probabilité d'erreur tolérée

    Returns
    -------
    mi : float
        la borne inférieure de l'intervalle
    ma : float
        la borne supérieure de l'intervalle

    References
    ----------
    http://www.i4.auc.dk/borre/matlab
    http://kom.aau.dk/~borre/matlab/
    """

    sux = np.sum(x)
    xb = np.mean(x)
    suy = np.sum(y)
    yb = np.mean(y)
    n = len(x)
    S1 = np.sum(x * y)
    S2 = sux * suy / n
    Sxy = S1 - S2
    S4 = np.sum(x**2)
    S5 = (sux**2) / n
    Sxx = S4 - S5
    S7 = np.sum(y**2)
    S8 = (suy**2) / n
    Syy = S7 - S8
    b1 = Sxy / Sxx
    b0 = yb - b1 * xb
    S14 = (Sxy**2) / Sxx
    s2y = (Syy - S14) / (n - 2)
    sy = np.sqrt(s2y)
    s2b1 = s2y / Sxx
    s2b0 = s2y * (1 / n + (xb**2) / Sxx)
    # t=tq(1-alpha/2,n-2)
    t = scipy.stats.t.ppf(1 - alpha / 2, n - 2)
    mi = b1 - t * np.sqrt(s2b1)
    ma = b1 + t * np.sqrt(s2b1)
    return mi, ma


def running_mean(data_in, window, convolve_mode="same"):
    """
    Gives running mean / moving average of data

    Parameters
    ----------
    data_in : list or numpy.array
        Values

    window :  float or int
        Size of the window for the running mean

    convolve_mode : str
        (expert) mode for the underlying convolution

    Returns
    -------
    data_run : numpy.array
        running mean of data_in (sane size as data_in)
        should stay "same"

    Note
    ----
    Nota :
        After a stress test, this one is the only one to
        provide an output with same size as input
        AND not shifted
        This fct is slow but at leat, do the job

        See running_mean_help for more details

        convolve_mode should stay fixed as "same"

    Nota 2 (for developpers) :
        Wrapper based on fct movingaverage_bis

        The substraction of the mean is an empirical trick
    """

    data_mean = np.nanmean(data_in)
    data_zero_centered = data_in - data_mean

    data_run = movingaverage_bis(data_zero_centered, window)

    return data_run + data_mean


def running_mean_help():
    help_str = """
    
Running Means functions with INTERNAL ID 1 3 5
return an Y output shorter than the input : not very convenient 
to align it on a X vector ...

INTERNAL ID 2 gives Y with same size as X
but result is shifted
Y[0] should be aligned with the middle of the 1st window

INTERNAL ID 4 is selected. It's declared as slow on StackOverflow,
but at least the job is done.
https://stackoverflow.com/questions/11352047/finding-moving-average-from-data-points-in-python

About the speed, the ID4 is not the slowest in fact,
it is the ID2 which is totally slow
so this weakness has to be relativised
(the answer is pretty old, so the convolve fct might have been improved)

About the convolution mode, it is detailed here : 
https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
"valid" mode is advised BUT do the same jobs as the others fcts (smaller output)
The "same" mode is the best one for our applications.
And "full" mode is not working either.

BUT all the fcts are maintened here, because they can be usefull for some 
other cases !
    
##### EXEMPLE STRESS TEST CODE
X = np.arange(1,1000,1) / (np.pi * 6)

plt.clf()
Ytrue = np.sin(X) * 100
Y = Ytrue + np.random.randn(len(X)) * 100

plt.plot(X,Y)

N = 50

Y1 = stats.movingaverage(Y,N)
Y2 = stats.runningMean(Y,N)
Y3 = stats.running_mean_core(Y,N)
Y4a = stats.movingaverage_bis(Y,N,"same")
Y4b = stats.movingaverage_bis(Y,N,"full")
Y5 = stats.movingaverage_ter(Y,N)

plt.clf()
plt.plot(Y)
plt.plot(Ytrue)
plt.plot(Y1,"r.")
plt.plot(Y2,"b.")
plt.plot(Y3,"r.")
plt.plot(Y4a,"y.")
plt.plot(Y5,"g.")


plt.clf()
plt.plot(X,Y)
plt.plot(X,Ytrue)
#plt.plot(X,Y1,"r.")
plt.plot(X,Y2,"b.")
#plt.plot(X,Y3,"r.")
plt.plot(X,Y4a,"y.")
#plt.plot(X,Y4b,"r.")
#plt.plot(X,Y5,"g.")
    """

    return help_str


def movingaverage(values, window):
    """
    Calculate moving average using convolution (INTERNAL_ID_1).

    Parameters
    ----------
    values : array-like
        Input values.
    window : int
        Window size for the moving average.

    Returns
    -------
    smas : numpy.ndarray
        Moving average calculated with 'valid' mode.
        The output size will be len(values) - window + 1.

    Notes
    -----
    The 'valid' mode requires enough datapoints. For example, if you remove 'valid',
    it will start at point one without any prior points.
    This implementation returns a shorter array than the input.

    References
    ----------
    http://sentdex.com/sentiment-analysisbig-data-and-python-tutorials-algorithmic-trading/how-to-chart-stocks-and-forex-doing-your-own-financial-charting/calculate-simple-moving-average-sma-python/
    """
    weigths = np.repeat(1.0, window) / window
    smas = np.convolve(values, weigths, "valid")
    return smas


def runningMean(x, N):
    """
    Calculate running mean (INTERNAL_ID_2).

    Parameters
    ----------
    x : array-like
        Input values.
    N : int
        Window size.

    Returns
    -------
    y : numpy.ndarray
        Running mean with same size as input.

    Notes
    -----
    This implementation returns an output with the same size as input,
    but the result is shifted. Y[0] should be aligned with the middle of the first window.

    References
    ----------
    http://stackoverflow.com/questions/13728392/moving-average-or-running-mean
    """
    y = np.zeros((len(x),))
    for ctr in range(len(x)):
        y[ctr] = np.sum(x[ctr : (ctr + N)])
    return y / N


def running_mean_core(x, N):
    """
    Calculate running mean using cumulative sum (INTERNAL_ID_3).

    Parameters
    ----------
    x : array-like
        Input values.
    N : int
        Window size.

    Returns
    -------
    xout : numpy.ndarray
        Running mean with reduced size.

    Notes
    -----
    More efficient than runningMean but returns output shorter than input.

    References
    ----------
    https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
    """
    cumsum = np.cumsum(np.insert(x, 0, 0))
    xout = (cumsum[N:] - cumsum[:-N]) / N
    return xout


def movingaverage_bis(interval, window_size, convolve_mode="same"):
    """
    Calculate moving average using convolution with specified mode (INTERNAL_ID_4).

    Parameters
    ----------
    interval : array-like
        Input values.
    window_size : int
        Size of the convolution window.
    convolve_mode : str, optional
        Mode for convolution: 'same', 'valid', or 'full'.
        The default is 'same', which returns output with same size as input.

    Returns
    -------
    numpy.ndarray
        Moving average with specified convolution mode applied.

    Notes
    -----
    Slower than other methods but gives output of the same size as input.
    The 'same' mode is recommended for most applications.

    References
    ----------
    https://stackoverflow.com/questions/11352047/finding-moving-average-from-data-points-in-python
    """
    window = np.ones(int(window_size)) / float(window_size)
    return np.convolve(interval, window, convolve_mode)


def movingaverage_ter(data, window_width):
    """
    Calculate moving average using cumulative sum (INTERNAL_ID_5).

    Parameters
    ----------
    data : array-like
        Input values.
    window_width : int
        Width of the moving window.

    Returns
    -------
    ma_vec : numpy.ndarray
        Moving average vector with reduced size compared to input.

    Notes
    -----
    Returns output shorter than the input.

    References
    ----------
    https://stackoverflow.com/questions/11352047/finding-moving-average-from-data-points-in-python
    """


def sinusoide(T, A, omega, phi=0, f=None):
    """
    produce a sinusoidal waveform

    Parameters
    ----------
    T : float
        time variable.
    A : float
        amplitude, the peak deviation of the function from zero.
    omega : float, optional
        ω = 2πf, angular frequency, the rate of change of the
        function argument in units of radians per second.
    phi : float, optional
        phase, specifies (in radians) where in its cycle
        the oscillation is at t = 0.
        The default is 0.
    f : float
        ordinary frequency, the number of oscillations (cycles)
        that occur each second of time.
        If given, it overrides the angular frequency omega.
        Thus, to use it, declare also omega = 0
        The default is None.

    Returns
    -------
    float
        a sinusoidal waveform.

    Notes
    -----
    https://en.wikipedia.org/wiki/Sine_wave

    """

    if f:
        omega_use = 2 * np.pi * f
    else:
        omega_use = omega

    return A * np.sin(omega_use * T + phi)


def butter_lowpass(cutoff, fs, order=5):
    """
    Design a Butterworth low-pass filter.

    Parameters
    ----------
    cutoff : float
        Cutoff frequency.
    fs : float
        Sampling frequency.
    order : int, optional
        Order of the filter. The default is 5.

    Returns
    -------
    b : numpy.ndarray
        Numerator coefficients of the IIR filter.
    a : numpy.ndarray
        Denominator coefficients of the IIR filter.

    Notes
    -----
    This function uses scipy.signal.butter to design the filter.
    """
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = scipy.signal.butter(order, normal_cutoff, btype="low", analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order=5):
    """
    Apply a Butterworth low-pass filter to data.

    Parameters
    ----------
    data : array-like
        Input data to be filtered.
    cutoff : float
        Cutoff frequency.
    fs : float
        Sampling frequency.
    order : int, optional
        Order of the filter. The default is 5.

    Returns
    -------
    y : numpy.ndarray
        Filtered data.

    See Also
    --------
    butter_lowpass : Design the Butterworth low-pass filter.
    """
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = scipy.signal.lfilter(b, a, data)
    return y


def gaussian_filter_GFZ_style_smoother(tim_ref, dat_ref, width=7):
    """
    Gaussian filter to smooth data based on GFZ's GMT_plus.pm/gaussian_kernel.

    This is a deprecated slow version. Use gaussian_filter_GFZ_style_smoother_improved instead.

    Parameters
    ----------
    tim_ref : array-like
        The X/T component of the time series (in decimal days).
    dat_ref : array-like
        The Y component (the data).
    width : int, optional
        Size of the window. An odd number is best. The default is 7.

    Returns
    -------
    dat_smt : list
        Smoothed Y values.

    Warnings
    --------
    This version is very slow (dirty conversion of a Perl function).
    Use gaussian_filter_GFZ_style_smoother_improved instead.

    Notes
    -----
    Some other nice ideas for smoothing:
    http://scipy-cookbook.readthedocs.io/items/SignalSmooth.html
    https://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way
    https://stackoverflow.com/questions/32900854/how-to-smooth-a-line-using-gaussian-kde-kernel-in-python-setting-a-bandwidth

    See Also
    --------
    gaussian_filter_GFZ_style_smoother_improved : Pythonic version that should be used.
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


def gaussian_filter_GFZ_style_smoother_improved(tim_ref, dat_ref, width=7):
    """
    Gaussian filter to smooth data,
    based on GFZ's GMT_plus.pm/gaussian_kernel

    Parameters
    ----------
    tim_ref : iterable (list or array)
        the X/T component of the time serie (in decimal days!)
    dat_ref : iterable (list or array)
        the Y component (the data).
    width : int, optional
        size of the window
        (odd number is best ?).
        The default is 7.

    Returns
    -------
    dat_smt2 : array
        smoothed Y.

    Note
    ----
    Some other nice ideas here
    http://scipy-cookbook.readthedocs.io/items/SignalSmooth.html
    https://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way
    https://stackoverflow.com/questions/32900854/how-to-smooth-a-line-using-gaussian-kde-kernel-in-python-setting-a-bandwidth
    """

    tim_raw = tim_ref
    dat_raw = dat_ref

    dat_smt2 = []

    num_raw = len(tim_raw)

    for ismt in range(num_raw):

        tim_raw_work = np.delete(tim_raw, ismt)
        dat_raw_work = np.delete(dat_raw, ismt)

        X_lag = tim_raw[ismt] - tim_raw_work
        X_fac = np.exp(-((X_lag / width) ** 2) / 2)

        # X_fac[X_fac < 0.01] = 0.

        clean_bool = X_fac > 0.01
        ## It differs a bit of the official fct, because the next element
        ## following this criteria is included anyway

        dat_raw_clean = dat_raw_work[clean_bool]
        X_fac_clean = X_fac[clean_bool]

        dat_raw_clean = dat_raw_clean[~np.isnan(dat_raw_clean)]
        X_fac_clean = X_fac_clean[~np.isnan(X_fac_clean)]

        X_val = np.sum(np.multiply(dat_raw_clean, X_fac_clean))
        X_wht = np.sum(X_fac_clean)

        dat_smt2.append(X_val / X_wht)

    dat_smt2 = np.array(dat_smt2)

    return dat_smt2


# def interpolate_gaps(values, limit=None):
#    """
#    Fill gaps using linear interpolation, optionally only fill gaps up to a
#    size of `limit`.
#    """
#    values = np.asarray(values)
#    i = np.arange(values.size)
#    valid = np.isfinite(values)
#    filled = np.interp(i, i[valid], values[valid])
#
#    if limit is not None:
#        invalid = ~valid
#        for n in range(1, limit+1):
#            invalid[:-n] &= invalid[n:]
#        filled[invalid] = np.nan
#
#    return filled


def smooth(x, window_len=11, window="hanning"):
    """
    Smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the beginning and end part of the output signal.

    Note: works only for equally spaced data

    Parameters
    ----------
    x : array
        the input signal
    window_len : int
        the dimension of the smoothing window; should be an odd integer
    window : str
        the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
        flat window will produce a moving average smoothing.

    Returns
    -------
    array
        the smoothed signal

    Examples
    --------
    >>> t = np.linspace(-2, 2, 0.1)
    >>> x = np.sin(t) + np.random.randn(len(t))*0.1
    >>> y = smooth(x)

    See Also
    --------
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    Notes
    -----
    The window parameter could be the window itself if an array instead of a string
    length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.

    References
    ----------
    http://scipy-cookbook.readthedocs.io/items/SignalSmooth.html
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


def harmonic_mean(A):
    """
    Calculate the harmonic mean of a list or array.

    Parameters
    ----------
    A : array-like
        Input values.

    Returns
    -------
    float
        The harmonic mean of A.

    Notes
    -----
    The harmonic mean is defined as the reciprocal of the arithmetic mean of the reciprocals.
    """
    A = np.array(A)
    return len(A) / np.sum(1.0 / A)


def find_intersection(x1, y1, x2, y2):
    """
    Find intersection points of two curves defined by points (x1, y1) and (x2, y2).

    Parameters
    ----------
    x1 : array-like
        X coordinates of the first curve.
    y1 : array-like
        Y coordinates of the first curve.
    x2 : array-like
        X coordinates of the second curve.
    y2 : array-like
        Y coordinates of the second curve.

    Returns
    -------
    roots : numpy.ndarray
        X coordinates of intersection points.
    y_vals : numpy.ndarray
        Y coordinates of intersection points.

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


def wrapTo360(lon):
    """
    Wrap longitude values to the range [0, 360).

    Parameters
    ----------
    lon : float or array-like
        Longitude value(s) in degrees.

    Returns
    -------
    float or numpy.ndarray
        Wrapped longitude value(s).

    Notes
    -----
    This function is based on MATLAB's wrapTo360 function.
    """
    lon = np.mod(lon, 360)
    return lon


def wrapTo180(lon):
    """
    Wrap longitude values to the range [-180, 180).

    Parameters
    ----------
    lon : float or array-like
        Longitude value(s) in degrees.

    Returns
    -------
    float or numpy.ndarray
        Wrapped longitude value(s).

    Notes
    -----
    This function is based on MATLAB's wrapTo180 function.
    """
    if not (lon is np.array):
        notaarray = True
        lon = np.array([lon])
    else:
        notaarray = False
    q = (lon < -180) + (180 < lon)
    lon[q] = wrapTo360(lon[q] + 180) - 180
    if notaarray:
        lon = lon[0]
    return lon


# Low level statistic function


def rms_mean(A):
    """
    Calculate the root mean square (RMS) mean of a list or array.

    Parameters
    ----------
    A : array-like
        Input values.

    Returns
    -------
    float
        The RMS mean value.

    Notes
    -----
    The RMS mean is defined as: sqrt(mean(A^2))
    NaN values are ignored in the calculation.
    """
    return np.sqrt(np.nanmean(np.square(np.array(A, np.float64))))


def RMSmean(indata):
    """
    Calculate the root mean square (RMS) mean of a list or array.

    Parameters
    ----------
    indata : array-like
        Input values.

    Returns
    -------
    float
        The RMS mean value.

    Warnings
    --------
    This function is redundant with rms_mean. Use rms_mean instead.

    Notes
    -----
    The RMS mean is defined as: sqrt(mean(indata^2))
    NaN values are ignored in the calculation.
    """
    rms = np.sqrt(np.nanmean(np.square(indata)))
    return rms


def rms_mean_alternativ(A):
    """
    Calculate the "GRGS style" RMS of a list or array.

    This calculates the standard deviation (RMS of deviations from the mean).

    Parameters
    ----------
    A : array-like
        Input values.

    Returns
    -------
    float
        The RMS value of deviations from the mean.

    Notes
    -----
    This function calculates: sqrt(mean((A - mean(A))^2))
    which is essentially the standard deviation.
    As of 1808, this is basically the standard deviation.
    """
    return np.sqrt(np.nanmean(np.square(A - np.nanmean(A))))


def rms_mean_kouba(A, multipl_coef=3, deg_of_freedom=7):
    """
    Calculate RMS mean following Kouba's method.

    Parameters
    ----------
    A : array-like
        Input values.
    multipl_coef : int, optional
        Multiplication coefficient. The default is 3.
    deg_of_freedom : int, optional
        Degrees of freedom. The default is 7.

    Returns
    -------
    float
        The RMS value calculated using Kouba's method.

    Notes
    -----
    Formula: sqrt(sum(A^2) / (multipl_coef * len(A) - deg_of_freedom))
    """
    return np.sqrt(np.sum(np.square(A)) / (multipl_coef * len(A) - deg_of_freedom))


def mad(data, mode="median"):
    """
    Calculate the Median Absolute Deviation (MAD) of a list or array.

    Parameters
    ----------
    data : array-like
        Input values.
    mode : str, optional
        Calculation mode. Can be 'median' or 'mean'. The default is 'median'.

    Returns
    -------
    float
        The Median Absolute Deviation value.

    Notes
    -----
    When mode is 'median': MAD = median(|data - median(data)|)
    When mode is 'mean': MAD = mean(|data - mean(data)|)
    NaN values are ignored in the calculation.

    References
    ----------
    http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm
    """

    if mode == "median":
        MAD = np.nanmedian(np.abs(data - np.nanmedian(data)))
    elif mode == "mean":
        MAD = np.nanmean(np.abs(data - np.nanmean(data)))
    return MAD


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
    http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm
    http://web.ipac.caltech.edu/staff/fmasci/home/statistics_refs/BetterThanMAD.pdf
    """

    if seuil:
        log.warning("'seuil' argument for outlier_mad is deprecated !!!")
        threshold = seuil

    if convert_to_np_array:
        data = np.array(data)
    nbinp = float(len(data))
    MAD = mad(data, mode=mad_mode)
    med = np.nanmedian(data)

    if np.isclose(np.sum(np.abs(np.diff(data))), 0.0):
        if verbose:
            log.info("elimination ratio: 0 , all values are equals")
        return data, np.array([True] * len(data))

    if np.isclose(med, 0.0):
        if verbose:
            log.info("elimination ratio: 0 , null median")
        return data, np.array([True] * len(data))

    if np.isclose(MAD, 0.0):
        if verbose:
            log.info("elimination ratio: 0 , null MAD")
        return data, np.array([True] * len(data))

    diff = data - med
    MZS = 0.6745 * np.abs(diff) / MAD

    MZS[np.isnan(MZS)] = threshold * 10  ### makes a threshold virtually higher for NaN
    boolbad = MZS < threshold  ## remendir: False means that the values are bad.

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
    Y, X, threshold=3.5, verbose=False, detrend_first=False, return_booleans=False
):
    """
    clean the outlier of Y using the MAD approach
    and clean the corresponding values in X
    assuming that we have the function : X => Y(X)
    (be carefull, Y is the first argument)

    Parameters
    ----------
    Y : list or numpy.array
        Values

    X : list or numpy.array
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
    Yclean & Xclean : numpy.array

    bb : numpy.array (if return_booleans == True)
        Y-sized booleans
    """
    if detrend_first:
        _, Ywork = detrend_timeseries(X, Y)
    else:
        _, Ywork = np.array(X), np.array(Y)

    _, bb = outlier_mad(Ywork, threshold, verbose)

    Xclean = np.array(X)[bb]
    Yclean = np.array(Y)[bb]

    if not return_booleans:
        return Yclean, Xclean
    else:
        return Yclean, Xclean, bb


def outlier_above_below_simple(X, low_bound, upp_bound, return_booleans=True):
    """
    Gives values of X which are between low_bound & upp_bound

    Parameters
    ----------
    X : list or numpy.array
        Values

    low_bound & upp_bound  : float
        lower and upper bound of X values wished

    return_booleans : bool
        return booleans or not

    Returns
    -------
    Xout : numpy.array
        X between low_bound & upp_bound

    bbool : bool
         X-sized array of booleans

    """

    Xwork = np.array(X)

    if low_bound >= upp_bound:
        log.warning("lower bound >= upper bound !!!")
        log.warning("low_bond : ", low_bound)
        log.warning("upp_bond : ", upp_bound)

    bbool = (low_bound <= Xwork) & (Xwork <= upp_bound)

    Xout = Xwork[bbool]

    if return_booleans:
        return Xout, bbool
    else:
        return Xout


def outlier_above_below(
    X,
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
    X : array
        Input array
    threshold_values : single value (float) or a 2-tuple
        (lower bound threshold, upper bound threshold)

        WARN: those value(s) have to be positives.
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
        print debug information

    Returns
    -------
    Xout : numpy array
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
        ref_val = reference(X)
    else:
        ref_val = reference

    if theshold_relative_value in ("reference", None):
        relativ_val = reference
    elif callable(theshold_relative_value):
        relativ_val = theshold_relative_value(X)
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

    Xout, bbool = outlier_above_below_simple(X, ths_low, ths_upp)

    if return_booleans:
        return Xout, bbool
    else:
        return Xout


def outlier_above_below_binom(
    Y,
    X,
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
    Y : array
        The dependent variable
    X : array
        The independent variable
    threshold_values : single value (float) or a 2-tuple
        (lower bound threshold, upper bound threshold)
    reference : float or callable
        the central reference value
    theshold_absolute : bool
        if True threshold_values are absolutes values
    theshold_relative_value : str or function
        relative threshold value specification
    return_booleans : bool
        return booleans or not
    detrend_first : bool
        remove trend before outlier detection
    verbose : bool
        print debug information

    Returns
    -------
    Yout : array
        Filtered Y values
    Xout : array
        Corresponding filtered X values
    Xout : array
        Corresponding filtered X values
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
    Xout : numpy array
        X between low_bound & upp_bound

    bbool : numpy array
        X-sized array of booleans
    """

    if detrend_first:
        _, Ywork = detrend_timeseries(X, Y)
    else:
        _, Ywork = np.array(X), np.array(Y)

    _, bb = outlier_above_below(
        Ywork,
        threshold_values,
        reference=reference,
        theshold_absolute=theshold_absolute,
        theshold_relative_value=theshold_relative_value,
        return_booleans=True,
        verbose=verbose,
    )

    Xclean = np.array(X)[bb]
    Yclean = np.array(Y)[bb]

    if not return_booleans:
        return Yclean, Xclean
    else:
        return Yclean, Xclean, bb


def outlier_sigma(datasigmain, threshold=3):
    """
    Remove outliers from sigma values above a threshold.

    This is a very old and discontinued function, not very efficient.

    Parameters
    ----------
    datasigmain : array-like
        Input sigma values.
    threshold : float, optional
        Threshold multiplier of median sigma. The default is 3.

    Returns
    -------
    datasigmaout : numpy.ndarray
        Filtered sigma values.
    boolbad : numpy.ndarray
        Boolean array indicating which values were kept (not outliers).

    Notes
    -----
    This function is deprecated and not recommended for new code.
    A point with sigma > threshold * median(sigmas) is considered an outlier.

    Deprecation
    -----------
    This is a very old and discontinued function.
    """
    moy = np.median(datasigmain)
    marge = moy * threshold

    log.info("moy,threshold,margin: %f, %i, %f", moy, threshold, marge)

    boolbad = np.abs(datasigmain) < marge

    datasigmaout = datasigmain[boolbad]

    return datasigmaout, boolbad


def lagrange1(points):
    """
    Low level function to determine a lagrangian polynom

    Replace scipy.interpolate.lagrange which is HIGHLY instable

    Parameters
    ----------
    points : list of n-interable
        point list.

    Returns
    -------
    p : function
        function representing the polynom.

    Source
    ------
    from : https://gist.github.com/melpomene/2482930

    """

    def P(x):
        total = 0
        n = len(points)
        for i in range(n):
            xi, yi = points[i]

            def g(i, n):

                tot_mul = 1
                for j in range(n):
                    if i == j:
                        continue
                    xj, yj = points[j]
                    tot_mul *= (x - xj) / float(xi - xj)

                return tot_mul

            total += yi * g(i, n)
        return total

    return P


def lagrange2(X, Y):
    """
    Low level function to determine a lagrangian polynom

    Replace scipy.interpolate.lagrange which is HIGHLY instable

    this function is more pythonic, but slower thant lagrange1....

    Parameters
    ----------
    points : list of n-interable
        point list.

    Returns
    -------
    p : function
        function representing the polynom.

    Source
    ------
    from : https://gist.github.com/melpomene/2482930

    """

    def P(x_itrp):
        total = 0
        n = len(X)
        for i in range(n):

            def g(i, n):
                X_but_i = np.concatenate((X[:i], X[i + 1 :]))
                # mask = np.ones(len(X),dtype=bool)
                # mask[i] = False
                # X_but_i = np.concatenate((X[:i],X[i+1:]))
                # X_but_i = X[mask]
                # return np.product((x_itrp -X[mask])/(X[i] - X[mask]))

                return np.product((x_itrp - X_but_i) / (X[i] - X_but_i))

            total += Y[i] * g(i, n)
        return total

    return P


def lagrange_interpolate(tdata, ydata, titrp, n=10, t_type="datetime"):
    """
    Perform a temporal lagrangian interpolation
    the X-component is a time

    Parameters
    ----------
    tdata : iterable of datetime
        X/T component of the known points.
    ydata : iterable of floats
        Y component of the known points..
    titrp : iterable of datetime
        Epochs of the wished points.
    n : int, optional
        degree of the polynom. Better if even. The default is 10.
    t_type : str, optional
        type of the time component, can be "datetime", "posix" or "pandas_timestamp".
        The default is "datetime".
        pandas_timestamp is recommended for a more precise applications
        (nanosecond precision instead of microsecond for datetime)

    Returns
    -------
    Yintrp : float array
        output interpolated data.

    Tips
    ----
    Use conv.dt_range to generate the wished epochs range

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

    Yintrp = []

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

            Tuse = tdata_px[imin:imax]
            Yuse = ydata[imin:imax]

            Poly = lagrange1(list(zip(Tuse, Yuse)))
            # Poly = lagrange2(Tuse,Yuse)

        yintrp = Poly(tintrp)
        Yintrp.append(yintrp)

    return np.array(Yintrp)


def dates_middle(start, end):
    """
    Calculate the middle date between two dates.

    Parameters
    ----------
    start : datetime
        Start date.
    end : datetime
        End date.

    Returns
    -------
    datetime
        The middle date between start and end.
    """
    return start + (end - start) / 2


def time_win_basic(
    start,
    end,
    Tlisin,
    Datalisin,
    outposix=True,
    invert=False,
    out_array=False,
    out_boolis=False,
    only_boolis=False,
):
    """
    Filter time series data within a specified time window.

    Parameters
    ----------
    start : datetime or float
        Start time of the window (POSIX if float).
    end : datetime or float
        End time of the window (POSIX if float).
    Tlisin : array-like of datetime or float
        Time values of the input time series.
    Datalisin : array-like
        Data values of the input time series.
    outposix : bool, optional
        If True, output times are in POSIX format. The default is True.
    invert : bool, optional
        If True, invert the boolean mask (select values outside the window).
        The default is False.
    out_array : bool, optional
        If True, convert outputs to numpy arrays. The default is False.
    out_boolis : bool, optional
        If True, also return the boolean mask. The default is False.
    only_boolis : bool, optional
        If True, only compute and return the boolean mask (faster).
        The default is False.

    Returns
    -------
    Tlisout : array or None
        Time values within the window (or outside if invert=True).
        None if only_boolis=True.
    Datalisout : array or None
        Data values within the window (or outside if invert=True).
        None if only_boolis=True.
    boolis : numpy.ndarray, optional
        Boolean mask array if out_boolis=True or only_boolis=True.

    Notes
    -----
    Internally works with POSIX time format.
    """

    if isinstance(Tlisin[0], dt.datetime):
        Tlis = conv.dt2posix(Tlisin)
    else:
        Tlis = Tlisin

    if isinstance(start, dt.datetime):
        start = conv.dt2posix(start)
    if isinstance(end, dt.datetime):
        end = conv.dt2posix(end)

    if not isinstance(Tlis, np.ndarray) or not isinstance(Datalisin, np.ndarray):
        Tlis = np.array(Tlis)
        Datalis = np.array(Datalisin)
    else:
        Tlis = Tlis
        Datalis = Datalisin

    boolis = (start <= Tlis) * (Tlis <= end)

    if invert:
        boolis = np.logical_not(boolis)

    if only_boolis:
        Datalisout = None
        Tlisout = None
    else:
        Datalisout = Datalis[boolis]
        Tlisout = Tlis[boolis]

        if not outposix:
            Tlisout = conv.posix2dt(Tlisout)

        if out_array:
            Tlisout, Datalisout = np.array(Tlisout), np.array(Datalisout)

    if out_boolis:
        out_tuple = (Tlisout, Datalisout, boolis)
    else:
        out_tuple = (Tlisout, Datalisout)

    return out_tuple


def time_win_multi(
    start, end, Tlist, Datalislis, outposix=True, invert=False, out_array=False
):
    """
    Filter multiple time series data within a specified time window.

    Parameters
    ----------
    start : datetime or float
        Start time of the window.
    end : datetime or float
        End time of the window.
    Tlist : array-like
        Time values (common to all data series).
    Datalislis : list of array-like
        List of data arrays to filter.
    outposix : bool, optional
        If True, output times are in POSIX format. The default is True.
    invert : bool, optional
        If True, invert the selection. The default is False.
    out_array : bool, optional
        If True, convert outputs to numpy arrays. The default is False.

    Returns
    -------
    Tlisout : array
        Time values within the window.
    Datalislisout : list of arrays
        Filtered data arrays.

    See Also
    --------
    time_win_basic : Filter a single time series.
    """
    Datalislisout = []
    for i, datalis in enumerate(Datalislis):

        Tlisout, datalisout = time_win_basic(
            start, end, Tlist, datalis, outposix, invert, out_array=out_array
        )
        Datalislisout.append(datalisout)
    return Tlisout, Datalislisout


def time_win_multi_start_end(
    Start_list_in,
    End_list_in,
    Tlisin,
    Datalisin,
    outposix=True,
    invert=False,
    out_array=False,
    out_boolis=False,
):
    """
    Filter time series data within multiple specified time windows.

    All windows must be satisfied simultaneously (logical AND operation).

    Parameters
    ----------
    Start_list_in : array-like
        List of start times for each window.
    End_list_in : array-like
        List of end times for each window.
    Tlisin : array-like
        Time values of the input time series.
    Datalisin : array-like
        Data values of the input time series.
    outposix : bool, optional
        If True, output times are in POSIX format. The default is True.
    invert : bool, optional
        If True, invert the selection. The default is False.
    out_array : bool, optional
        If True, convert outputs to numpy arrays. The default is False.
    out_boolis : bool, optional
        If True, also return boolean masks. The default is False.

    Returns
    -------
    Tlisout : array
        Time values satisfying all windows.
    Datalisout : array
        Data values satisfying all windows.
    boolis_opera : numpy.ndarray, optional
        Combined boolean mask if out_boolis=True.
    boolis_stk : numpy.ndarray, optional
        Stack of individual window boolean masks if out_boolis=True.

    Notes
    -----
    Internally works in POSIX time format.
    A point is included only if it falls within ALL specified time windows (AND operation).

    See Also
    --------
    time_win_basic : Filter using a single time window.
    """

    if len(Start_list_in) != len(End_list_in):
        log.error("len(Start_list_in) != len(End_list_in) !!")

    boolis_stk = []
    for start, end in zip(Start_list_in, End_list_in):
        _, _, boolis = time_win_basic(
            start,
            end,
            Tlisin,
            Datalisin,
            outposix=outposix,
            invert=invert,
            out_boolis=True,
            only_boolis=True,
        )

        boolis_stk.append(boolis)

    boolis_stk = np.stack(boolis_stk)
    boolis_opera = np.all(boolis_stk, axis=0)

    Datalis = np.array(Datalisin)
    Tlis = np.array(Tlisin)

    Datalisout = Datalis[boolis_opera]
    Tlisout = Tlis[boolis_opera]

    if not outposix:
        Tlisout = conv.posix2dt(Tlisout)

    if out_array:
        Tlisout, Datalisout = np.array(Tlisout), np.array(Datalisout)

    if out_boolis:
        out_tuple = (Tlisout, Datalisout, boolis_opera, boolis_stk)
    else:
        out_tuple = (Tlisout, Datalisout)

    return out_tuple


def get_season(now):
    """
    Determine the season for a given date.

    Parameters
    ----------
    now : datetime or date
        The date to determine the season for.

    Returns
    -------
    str
        The season name: 'winter', 'spring', 'summer', or 'autumn'.

    Notes
    -----
    Season boundaries are based on the Northern Hemisphere calendar:
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
    Get a color representation of the season.

    Parameters
    ----------
    datein : datetime or date
        The date to get the season color for.

    Returns
    -------
    str
        Color code: 'b' (blue) for winter, 'r' (red) for summer,
        'g' (green) for spring, 'k' (black) for autumn.

    Notes
    -----
    Useful for plotting time series with color-coded seasons.

    See Also
    --------
    get_season : Get the season name for a given date.
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


def outlier_mad_binom_legacy(
    X, Y, threshold=3.5, verbose=False, detrend_first=False, return_booleans=False
):
    """
    Remove outliers from both X and Y using MAD approach (legacy version).

    Parameters
    ----------
    X : array-like
        X values (dependent variable in this legacy version).
    Y : array-like
        Y values (independent variable in this legacy version).
    threshold : float, optional
        MAD threshold. The default is 3.5.
    verbose : bool, optional
        Print elimination statistics. The default is False.
    detrend_first : bool, optional
        Remove linear trend before outlier detection. The default is False.
    return_booleans : bool, optional
        Return boolean mask along with cleaned data. The default is False.

    Returns
    -------
    Xclean : numpy.ndarray
        Cleaned X values.
    Yclean : numpy.ndarray
        Cleaned Y values.
    bb : numpy.ndarray, optional
        Boolean mask if return_booleans=True.

    Warnings
    --------
    This is a legacy function where the order of X and Y is different from the main version.
    The detrend operation might be unstable.
    Use outlier_mad_binom instead for new code.

    Deprecation
    -----------
    This is the legacy version with inverted argument order.
    Use the main version outlier_mad_binom for new code.

    See Also
    --------
    outlier_mad_binom : The recommended version with correct argument order.
    """


#  ______                _   _                _____                                         _
# |  ____|              | | (_)              / ____|                                       | |
# | |__ _   _ _ __   ___| |_ _  ___  _ __   | |  __ _ __ __ ___   _____ _   _  __ _ _ __ __| |
# |  __| | | | '_ \ / __| __| |/ _ \| '_ \  | | |_ | '__/ _` \ \ / / _ \ | | |/ _` | '__/ _` |
# | |  | |_| | | | | (__| |_| | (_) | | | | | |__| | | | (_| |\ v /  __/ |_| | (_| | | | (_| |
# |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|  \_____|_|  \__,_| \_/ \___|\__, |\__,_|_|  \__,_|
#                                                                        __/ |
#                                                                       |___/


def outlier_mad_binom_legacy(
    X, Y, threshold=3.5, verbose=False, detrend_first=False, return_booleans=False
):
    """
    clean the outlier of X and clean the corresponding values in Y

    legacy : order of X Y is different than in the main version, and here
    it might be unstable for the detrend
    """
    if detrend_first:
        Xwork, _ = detrend_timeseries(X, Y)
    else:
        Xwork, _ = np.array(X), np.array(Y)

    _, bb = outlier_mad(Xwork, threshold, verbose)

    Xclean = np.array(X)[bb]
    Yclean = np.array(Y)[bb]

    if not return_booleans:
        return Xclean, Yclean
    else:
        return Xclean, Yclean, bb


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
