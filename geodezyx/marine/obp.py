#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 21:08:05 2023

@author: psakic

This module regroups the functions for the exploitation of the A0A pressure
sensors in the context of the REVOSIMA network

It is based on the work of Yann-Terden Tranchant, LIENSs, La Rochelle, France
"""


import numpy as np
import pandas as pd
import xarray as xr
# import gsw
#from scipy.signal import butter, filtfilt
import scipy.signal as signal
import scipy.optimize as optimize
import datetime as dt


# from datetime import datetime, timedelta

#import numpy as np
#from scipy.signal import butter

GRAVITY = 9.81

def butter_highpass(cutoff, fs, order=5):
    """
    Designs a highpass Butterworth filter.

    Parameters
    ----------
    cutoff : float
        The cutoff frequency of the filter.
    fs : float
        The sampling frequency of the signal.
    order : int, optional
        The order of the filter. Default is 5.

    Returns
    -------
    b : ndarray
        The numerator coefficient vector of the filter.
    a : ndarray
        The denominator coefficient vector of the filter.
    """
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    out_butter = signal.butter(order, normal_cutoff, btype='high', analog=False, output='ba')
    b, a = out_butter[0], out_butter[1]
    return b, a

def butter_lowpass(cutoff, fs, order=5):
    """
    Designs a lowpass Butterworth filter.

    Parameters
    ----------
    cutoff : float
        The cutoff frequency of the filter.
    fs : float
        The sampling frequency of the signal.
    order : int, optional
        The order of the filter. Default is 5.

    Returns
    -------
    b : ndarray
        The numerator coefficient vector of the filter.
    a : ndarray
        The denominator coefficient vector of the filter.
    """
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    out_butter = signal.butter(order, normal_cutoff, btype='low', analog=False, output='ba')
    b, a = out_butter[0], out_butter[1]

    return b, a

def butter_highpass_filtfilt(data, cutoff, fs, order=4):
    """
    Applies a highpass Butterworth filter to the data using forward and backward filtering.

    Parameters
    ----------
    data : array_like
        The input data to be filtered.
    cutoff : float
        The cutoff frequency of the filter.
    fs : float
        The sampling frequency of the signal.
    order : int, optional
        The order of the filter. Default is 4.

    Returns
    -------
    y : ndarray
        The filtered data.
    """
    b, a = butter_highpass(cutoff, fs, order=order)
    y = signal.filtfilt(b, a, data)
    return y

def butter_lowpass_filtfilt(data, cutoff, fs, order=4):
    """
    Applies a lowpass Butterworth filter to the data using forward and backward filtering.

    Parameters
    ----------
    data : array_like
        The input data to be filtered.
    cutoff : float
        The cutoff frequency of the filter.
    fs : float
        The sampling frequency of the signal.
    order : int, optional
        The order of the filter. Default is 4.

    Returns
    -------
    y : ndarray
        The filtered data.
    """
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = signal.filtfilt(b, a, data)
    return y

def butter_bandpass(lowcut, highcut, fs, order=3):
    """
    Designs a bandpass Butterworth filter.

    Parameters
    ----------
    lowcut : float
        The lower cutoff frequency of the filter.
    highcut : float
        The upper cutoff frequency of the filter.
    fs : float
        The sampling frequency of the signal.
    order : int, optional
        The order of the filter. Default is 3.

    Returns
    -------
    b : ndarray
        The numerator coefficient vector of the filter.
    a : ndarray
        The denominator coefficient vector of the filter.
    """
    return signal.butter(order, [lowcut, highcut], fs=fs, btype='band')

def butter_bandpass_filtfilt(data, lowcut, highcut, fs, order=4):
    """
    Applies a bandpass Butterworth filter to the data using forward and backward filtering.

    Parameters
    ----------
    data : array_like
        The input data to be filtered.
    lowcut : float
        The lower cutoff frequency of the filter.
    highcut : float
        The upper cutoff frequency of the filter.
    fs : float
        The sampling frequency of the signal.
    order : int, optional
        The order of the filter. Default is 4.

    Returns
    -------
    y : ndarray
        The filtered data.
    """
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = signal.filtfilt(b, a, data)
    return y

def butterworth(df, t0=3*24*3600, t1=10*24*3600, kind='bandpass', order=4):
    """
    Applies a Butterworth filter to the data.

    Parameters
    ----------
    df : pandas.DataFrame
        The input data to be filtered.
    t0 : int, optional
        The lower cutoff period in seconds. Default is 3 days (3*24*3600 seconds).
    t1 : int, optional
        The upper cutoff period in seconds. Default is 10 days (10*24*3600 seconds).
    kind : str, optional
        The type of filter to apply. Options are 'highpass', 'lowpass', and 'bandpass'. Default is 'bandpass'.
    order : int, optional
        The order of the filter. Default is 4.

    Returns
    -------
    pandas.Series
        The filtered data.
    """
    difftime = np.diff(df.index)
    if not df.index[0].tzname():
        pass  # directly nanosec
    else:
        difftime = np.array([e.total_seconds() * 10**9 for e in difftime])  # TimeDelta to be conv in nanosec

    if not np.all(difftime == difftime[0]):
        print('The sampling is not regular, return')
        return
    else:
        fs = 1/(float(difftime[0])/1e9)

    if kind == 'highpass':
        cut = 1/t0
        y = butter_highpass_filtfilt(df.values, cut, fs, order=order)
    elif kind == 'lowpass':
        cut = 1/t0
        y = butter_lowpass_filtfilt(df.values, cut, fs, order=order)
    elif kind == 'bandpass':
        lowcut = 1/t1
        highcut = 1/t0
        y = butter_bandpass_filtfilt(df.values, lowcut, highcut, fs, order=order)
    else:
        print(f'Invalid filter type : {kind}. Must be highpass, lowpass, or bandpass.')
        return

    return pd.DataFrame(y, index=df.index)[0]

def read_hycom(file):
    """
    Reads a HYCOM dataset and renames its variables.

    Parameters
    ----------
    file : str
        The path to the HYCOM dataset file.

    Returns
    -------
    xarray.Dataset
        The renamed dataset.
    """
    ds = xr.open_dataset(file)
    return ds.rename(lon='longitude', lat='latitude', water_temp='theta', salinity='salt', surf_el='ssh')

def read_ecco2(file):
    """
    Reads an ECCO2 dataset and renames its variables.

    Parameters
    ----------
    file : str
        The path to the ECCO2 dataset file.

    Returns
    -------
    xarray.Dataset
        The renamed dataset.
    """
    ds = xr.open_dataset(file)
    return ds.rename(LEV='depth')

def read_glorys(file):
    """
    Reads a GLORYS dataset and renames its variables.

    Parameters
    ----------
    file : str
        The path to the GLORYS dataset file.

    Returns
    -------
    xarray.Dataset
        The renamed dataset.
    """
    ds = xr.open_dataset(file, decode_times=True)
    return ds.rename(thetao='theta', so='salt', zos='ssh')

def read_duacs(file):
    """
    Reads a DUACS dataset and renames its variables.

    Parameters
    ----------
    file : str
        The path to the DUACS dataset file.

    Returns
    -------
    xarray.Dataset
        The renamed dataset.
    """
    ds = xr.open_dataset(file)
    return ds.rename(sla='ssh')

def interp_xy(ds, x=None, y=None, method = 'linear'):
    """
    frontend for xarray spatial interpolation

    Parameters
    ----------
    ds : xarray.Dataset
        input xarray grid.
    x : float, optional
        longitude in deg (but depends on ds units). The default is None.
    y : float, optional
        latitude in deg (but depends on ds units). The default is None.
    method : str or Litteral, optional
        interpolation method. The default is 'linear'.

    Returns
    -------
    ds
        Interpolated values.

    """
    if not x and not y:
        print('Must give x or/and y coord')
        return
    if not x:
        return ds.interp(latitude = y, method = method)
    elif not y:
        return ds.interp(longitude = x, method = method)
    else:
        return ds.interp(longitude = x, latitude = y, method = method)


def interp_time(ds, time, method='linear'):
    """
    Interpolates the dataset along the time dimension.

    Parameters
    ----------
    ds : xarray.Dataset
        The input dataset to be interpolated.
    time : array-like or scalar
        The new time points to interpolate to.
    method : str or Litteral, optional
        The interpolation method to use. Default is 'linear'.

    Returns
    -------
    xarray.Dataset
        The interpolated dataset.
    """
    return ds.interp(time=time, method=method)

def interp_z(ds, depth, method='linear'):
    """
    Interpolates the dataset along the depth dimension.

    Parameters
    ----------
    ds : xarray.Dataset
        The input dataset to be interpolated.
    depth : array-like or scalar
        The new depth points to interpolate to.
    method : str or Litteral, optional
        The interpolation method to use. Default is 'linear'.

    Returns
    -------
    xarray.Dataset
        The interpolated dataset.
    """
    return ds.interp(depth=depth, method=method)

def get_bottom_depth(ds):
    """
    Computes the bottom depth of the dataset.

    Parameters
    ----------
    ds : xarray.Dataset
        The input dataset containing depth and theta variables.

    Returns
    -------
    numpy.ndarray
        The bottom depth values.
    """
    land = np.isnan(ds.theta[0, 0, :, :])
    isnan = np.isnan(ds.theta[0, 0:, :, :])
    ind_nan = np.count_nonzero(~isnan, axis=0)
    bottom_depth = ds.depth.values[ind_nan]
    bottom_depth[land.values] = np.nan
    return bottom_depth

def extract_profile(ds, x, y, method='linear', zbottom=None):
    """
    Extracts a profile from the dataset at the specified coordinates and interpolates to the given depth.

    Parameters
    ----------
    ds : xarray.Dataset
        The input dataset containing the data.
    x : float
        The longitude coordinate.
    y : float
        The latitude coordinate.
    method : str, optional
        The interpolation method to use. Default is 'linear'.
    zbottom : int, float, list, or numpy.ndarray, optional
        The depth(s) to interpolate to. If None, no depth interpolation is performed. Default is None.

    Returns
    -------
    xarray.Dataset
        The extracted and interpolated profile.
    """
    ds = interp_xy(ds, x, y, method)

    if not isinstance(zbottom, type(None)):
        if isinstance(zbottom, (int, float)):
            valid = ~np.isnan(ds.theta[0])
            max_depth = ds.depth[valid][-1]
            if zbottom > max_depth.values:
                print(f'{zbottom} m is deeper than the maximum depth layer of the profile : {max_depth.values}')
                return None
            k = len((zbottom - ds.depth)[(zbottom - ds.depth) > 0]) - 1
            new_depth = np.hstack([ds.depth[:k], zbottom])
        elif isinstance(zbottom, (np.ndarray, list)):
            new_depth = zbottom
        else:
            print(f'Invalid zbottom parameter type : {type(zbottom)}. Must be int or float depth (in m).')
            return None
        ds = interp_z(ds, new_depth, method='linear')

    return ds

def compute_dens_profile(profile):
    """
    Computes the density profile from the given dataset.

    Parameters
    ----------
    profile : xarray.Dataset
        The input dataset containing temperature (theta) and salinity (salt) variables.

    Returns
    -------
    numpy.ndarray
        The computed density profile.
    """
    import seawater
    ds = profile
    if np.ndim(ds.theta) != 2:
        print(f'Bad dataset dimension : {np.ndim(ds.theta)}. The function must be applied on a profile.')
        return None
    lat = ds.latitude.mean()
    pres = seawater.pres(ds.depth, lat=lat)
    dens = seawater.dens(ds.salt, ds.theta, pres)
    return dens

def compute_steric(profile, integration='forward'):
    """
    Computes the steric height anomaly from the given profile.

    Parameters
    ----------
    profile : xarray.Dataset
        The input dataset containing temperature (theta) and salinity (salt) variables.
    integration : str, optional
        The integration method to use. Options are 'forward', 'backward', and 'averaged'. Default is 'forward'.

    Returns
    -------
    pandas.Series
        The steric height anomaly.
    """
    ds = profile
    #lat = ds.latitude.mean()
    dens = compute_dens_profile(ds)
    if integration == 'forward':
        steric = pd.DataFrame(GRAVITY * np.nansum(np.diff(ds.depth) * dens[:, :-1], axis=1) / 10000, index=ds.time)[0]
    elif integration == 'backward':
        steric = pd.DataFrame(GRAVITY * np.nansum(np.diff(ds.depth) * dens[:, 1:], axis=1) / 10000, index=ds.time)[0]
    elif integration == 'averaged':
        dens = np.vstack([np.sum(dens[:, i:i+2], axis=1) / 2 for i in range(len(ds.depth) - 1)]).T
        steric = pd.DataFrame(GRAVITY * np.nansum(np.diff(ds.depth) * dens, axis=1) / 10000, index=ds.time)[0]
    else:
        print(f'Invalid integration method : {integration}. Must be forward, backward, or averaged.')
        return None

    return steric - steric.median()

def compute_phibot(profile, ssh=None, integration='forward', rho=10.35, remove_median=True):
    """
    Computes the ocean hydrostatic bottom pressure anomaly (phibot).

    Parameters
    ----------
    profile : xarray.Dataset
        The input dataset containing temperature (theta) and salinity (salt) variables.
    ssh : xarray.DataArray, optional
        The sea surface height data. If None, it is taken from the profile. Default is None.
    integration : str, optional
        The integration method to use for steric height computation.
        Options are 'forward', 'backward', and 'averaged'. Default is 'forward'.
    rho : float, optional
        The reference density in kg/m^3. Default is 10.35.
    remove_median : bool, optional
        Whether to remove the median from the computed phibot values. Default is True.

    Returns
    -------
    pandas.DataFrame
        The computed phibot values.

    Notes
    -----
    phibot = Ocean hydrostatic bottom pressure anomaly
    https://cmr.earthdata.nasa.gov/search/concepts/V2028471168-POCLOUD.html
    https://cmr.earthdata.nasa.gov/search/concepts/V2146301108-POCLOUD.html


    PHIBOT		Bottom Pressure Pot. Anomaly (p/rhonil, m^2/s^2)
                To convert to m, divide by g (g=9.81 m/s^2)
		PHIBOT is the anomaly relative to Depth * rhonil * g
		The absolute bottom pressure in Pa is:
		Depth * rhonil * g + PHIBOT * rhonil (rhonil=1027.5 kg/m^3)
    http://apdrc.soest.hawaii.edu/doc/Readme_ecco2_cube92
    """
    if isinstance(ssh, type(None)):
        ssh = profile.ssh
    #lat = profile.latitude.mean()
    #dens = compute_dens_profile(profile)
    phibot = pd.DataFrame(index=profile.time)
    phibot['subsurface'] = compute_steric(profile, integration)
    phibot['surface'] = GRAVITY * rho * ssh / 100
    phibot['phibot'] = phibot.subsurface + phibot.surface
    if remove_median:
        phibot -= phibot.median(axis=0)

    return phibot
    
def compute_spectrogram(df, max_period=45, nchunks=3600*6):
    """
    Computes the spectrogram of the given data.

    Parameters
    ----------
    df : pandas.DataFrame
        The input data to compute the spectrogram for.
    max_period : int, optional
        The maximum period to consider in the spectrogram. Default is 45.
    nchunks : int, optional
        The number of chunks to divide the data into. Default is 3600*6.

    Returns
    -------
    p : numpy.ndarray
        The periods corresponding to the spectrogram.
    t : pandas.DatetimeIndex
        The time points corresponding to the spectrogram.
    sxx : numpy.ndarray
        The spectrogram values.
    """
    from scipy.signal import spectrogram

    f, t, sxx = spectrogram(df.values, nperseg=nchunks, scaling='density')
    p = (1/f[1:])
    sxx = sxx[1:]
    i = np.argmin(np.abs(p - max_period))
    p = p[i:]
    sxx = sxx[i:]

    t = pd.DatetimeIndex([df.index[0] + dt.timedelta(hours=t) for t in t/3600])
    return p, t, sxx

def log_linear(x, a1, a2, a3, a4, **kwargs):
    """
    Log-linear model function.

    Parameters
    ----------
    x : array-like
        The input data.
    a1 : float
        Coefficient for the logarithmic term.
    a2 : float
        Coefficient for the logarithmic term.
    a3 : float
        Coefficient for the linear term.
    a4 : float
        Constant term.
    **kwargs : dict
        Additional keyword arguments.

    Returns
    -------
    float
        The computed log-linear value.
    """
    return a1 * np.log(a2 * (1 + x)) + a3 * x + a4

def log(x, a1, a2, a3, a4, **kwargs):
    """
    Logarithmic model function.

    Parameters
    ----------
    x : array-like
        The input data.
    a1 : float
        Coefficient for the logarithmic term.
    a2 : float
        Coefficient for the logarithmic term.
    a3 : float
        Constant term.
    **kwargs : dict
        Additional keyword arguments.

    Returns
    -------
    float
        The computed logarithmic value.
    """
    x0 = min(x)
    return a1 * np.log(a2 * (1 + x - x0)) + a3 * (x - x0) + a4

def exp(x, a, b, c, d, **kwargs):
    """
    Exponential model function.

    Parameters
    ----------
    x : array-like
        The input data.
    a : float
        Coefficient for the exponential term.
    b : float
        Coefficient for the exponential term.
    c : float
        Constant term.
    **kwargs : dict
        Additional keyword arguments.

    Returns
    -------
    float
        The computed exponential value.
    """
    return a * np.exp(-b * x) + c * x + d 

def linear(x, a1, a2, a3, a4, **kwargs):
    """
    Linear model function.

    Parameters
    ----------
    x : array-like
        The input data.
    a1 : float
        Slope of the linear function.
    a2 : float
        Intercept of the linear function.
    a3 : float
        Additional coefficient.
    a4 : float
        Additional coefficient.
    **kwargs : dict
        Additional keyword arguments.

    Returns
    -------
    float
        The computed linear value.
    """
    x0 = min(x)
    return a1 * (x - x0) + a2

def exp_linear(x, a, b, c, d, **kwargs):
    """
    Exponential-linear model function.

    Parameters
    ----------
    x : array-like
        The input data.
    a : float
        Coefficient for the exponential term.
    b : float
        Coefficient for the exponential term.
    c : float
        Linear term coefficient.
    d : float
        Constant term.
    **kwargs : dict
        Additional keyword arguments.

    Returns
    -------
    float
        The computed exponential-linear value.
    """
    x0 = min(x)
    return a * np.exp(-1/b * (x - x0)) + c * (x - x0) + d


#### Added by PS (2024-10)
def exp_polynomial_deg2(x, a, b, c,  **kwargs):
    """
    Computes an exponential function combined with a second-degree polynomial.

    This function calculates the value of an exponential function combined with a
    second-degree polynomial given the input data and coefficients.

    Parameters
    ----------
    x : array-like
        The input data.
    a : float
        Coefficient for the polynomial term.
    b : float
        Coefficient for the exponential term.
    c : float
        Constant term.
    **kwargs : dict
        Additional keyword arguments.

    Returns
    -------
    float
        The computed value from the exponential function combined with a second-degree polynomial.
    """
    x0 = min(x)
    return (x - x0)**2 * a * np.exp(b * (x-x0) + c)

# Exponential with offset
def exp_with_offset(x, a, b, c, **kwargs):
    """
    Exponential model function with an offset.

    This function calculates an exponential decay with an added constant offset.

    Parameters
    ----------
    x : array-like
        The input data.
    a : float
        Coefficient for the exponential term.
    b : float
        Exponent for the exponential term.
    c : float
        Constant offset term.
    **kwargs : dict
        Additional keyword arguments.

    Returns
    -------
    float
        The computed value from the exponential model with offset.
    """
    return a * np.exp(-b * x) + c

# Higher degree polynomial (3rd degree)
def polynomial_3(x, a, b, c, d, **kwargs):
    """
    Computes a 3rd degree polynomial.

    This function calculates the value of a 3rd degree polynomial given the input data and coefficients.

    Parameters
    ----------
    x : array-like
        The input data.
    a : float
        The constant term of the polynomial.
    b : float
        The coefficient of the linear term.
    c : float
        The coefficient of the quadratic term.
    d : float
        The coefficient of the cubic term.
    **kwargs : dict
        Additional keyword arguments.

    Returns
    -------
    float
        The computed value of the 3rd degree polynomial.
    """
    return a + b*x + c*x**2 + d*x**3
# Generalized nonlinear model (exponential + power law)
def exp_power_combined(x, a, b, c, d, **kwargs):
    """
    Generalized nonlinear model function combining exponential decay and power law.

    This function models data using a combination of an exponential decay term and a power law term.

    Parameters
    ----------
    x : array-like
        The input data.
    a : float
        Coefficient for the exponential term.
    b : float
        Exponent for the exponential term.
    c : float
        Coefficient for the power law term.
    d : float
        Exponent for the power law term.
    **kwargs : dict
        Additional keyword arguments.

    Returns
    -------
    float
        The computed value from the combined exponential and power law model.
    """
    return a * np.exp(-b * x) + c * np.power(x + 1, -d)

def fit_model(data, t_fit, model='log_linear', offset_last=True, maxfev=1000, pn=2,
              normalize_t=False):
    """
    Fits a model to the given data.

    Parameters
    ----------
    data : pandas.Series
        The input data to fit the model to.
        Is a Series with time as Index
    t_fit : pandas.DatetimeIndex
        Time points to get fitted values once the fit model has been determeined.
    model : str, optional
        The model to use for fitting. Options are 'log_linear', 'log', 'exp',
        'linear', 'exp_linear', and 'poly'. Default is 'log_linear'.
    offset_last : bool, optional
        Whether to offset the last value. Default is True.
    maxfev : int, optional
        The maximum number of function evaluations. Default is 1000.
    pn : int, optional
        The degree of the polynomial if the model is 'poly'. Default is 2.
    normalize_t : bool, optional
        Whether to normalize the time values. Default is False.

    Returns
    -------
    pandas.Series
        The fitted model values.
    """

    # initial time implemented in days, but better in seconds
    #t_num_data = (data.index.to_julian_date() - data.index.to_julian_date()[0]).values
    #t_num_fit = (t_fit.to_julian_date() - t_fit.to_julian_date()[0]).values

    t_num_data = (data.index - data.index[0]).total_seconds().values
    t_num_fit = (t_fit - t_fit[0]).total_seconds().values

    ### normalize
    if normalize_t:
        ref = t_num_data[-1]
        t_num_data /= ref
        t_num_fit /= ref

    if model == 'poly':
        z = np.polyfit(t_num_data, data.values, pn)
        p = np.poly1d(z)
        fit = p(t_num_fit)
        ind_min = np.argmin(fit)
        fit[ind_min:] = fit[ind_min]
        out_opti = None
    elif model == 'log_linear':
        out_opti = optimize.curve_fit(log_linear, t_num_data, data.values, maxfev=maxfev)
        popt, pcov = out_opti[0], out_opti[1]
        fit = log_linear(t_num_fit, *popt)
    elif model == 'log':
        out_opti = optimize.curve_fit(log, t_num_data, data.values, maxfev=maxfev)
        popt, pcov = out_opti[0], out_opti[1]
        fit = log(t_num_fit, *popt)
    elif model == 'exp':
        out_opti = optimize.curve_fit(exp, t_num_data, data.values, maxfev=maxfev)
        popt, pcov = out_opti[0], out_opti[1]
        fit = exp(t_num_fit, *popt)
    elif model == 'exp_linear':
        out_opti = optimize.curve_fit(exp_linear, t_num_data, data.values, maxfev=maxfev)
        popt, pcov = out_opti[0], out_opti[1]
        fit = exp_linear(t_num_fit, *popt)
    elif model == 'exp_polynomial_deg2':
        out_opti = optimize.curve_fit(exp_polynomial_deg2, t_num_data, data.values, maxfev=maxfev)
        popt, pcov = out_opti[0], out_opti[1]
        fit = exp_polynomial_deg2(t_num_fit, *popt)
    elif model == 'linear':
        out_opti = optimize.curve_fit(linear, t_num_data, data.values, maxfev=maxfev)
        popt, pcov = out_opti[0], out_opti[1]
        fit = linear(t_num_fit, *popt)
    elif model == 'exp_power_combined':
        out_opti = optimize.curve_fit(exp_power_combined, t_num_data, data.values, maxfev=maxfev)
        popt, pcov = out_opti[0], out_opti[1]
        fit = exp_power_combined(t_num_fit, *popt)
    else:
        print(f'Invalid model type : {model}. Must be log_linear, log, exp, linear, exp_linear, or poly.')
        return None, None

    return pd.DataFrame(fit, index=t_fit)[0], out_opti