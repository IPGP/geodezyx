#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 21:08:05 2023

@author: psakic

This module regroups the functions for the exploitation of the A0A pressure
sensors in the context of the REVOSIMA network

It is based on the work of Yann Terden Tranchant

"""


import xarray as xr
import numpy as np
import seawater
import pandas as pd
#from datetime import datetime, timedelta
import scipy.optimize as optimize

from scipy.signal import butter, lfilter, sosfilt, filtfilt

#import numpy as np
import matplotlib.pyplot as plt
#from scipy.signal import butter

def butter_highpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_highpass_filtfilt(data, cutoff, fs, order=4):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y

def butter_lowpass_filtfilt(data, cutoff, fs, order=4):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y

def butter_bandpass(lowcut, highcut, fs, order=3):
    return butter(order, [lowcut, highcut], fs=fs, btype='band')

def butter_bandpass_filtfilt(data, lowcut, highcut, fs, order=4):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = filtfilt(b, a, data)
    return y

def butterworth(df, t0 = 3*24*3600, t1 = 10*24*3600, kind = 'bandpass', order = 4):
    
    dt = np.diff(df.index)
    if not df.index[0].tzname():
        pass # directly nanosec
    else:
        dt = np.array([e.total_seconds() * 10**9 for e in dt])  # TimeDelta to be conv in nanosec
            
    
    
    if not np.all(dt == dt[0]):
        print('The sampling is not regular, return')
        return
    else:
        fs = 1/(float(dt[0])/1e9)
        
    if kind == 'highpass':
        cut = 1/t0
        y = butter_highpass_filtfilt(df.values, cut, fs, order=order)
    elif kind == 'lowpass':
        cut = 1/t0
        y = butter_lowpass_filtfilt(df.values, cut, fs, order=order)
    elif kind == 'bandpass':
        lowcut = 1/t1
        highcut = 1/t0
        y = butter_bandpass_filtfilt(df.values, lowcut,
                                     highcut, fs, order=order)
    return pd.DataFrame(y, index = df.index)[0]
    
def read_hycom(file):
    ds = xr.open_dataset(file)
    return ds.rename(lon = 'longitude', lat = 'latitude',
                     water_temp = 'theta', salinity = 'salt', surf_el = 'ssh')

def read_ecco2(file):
    ds = xr.open_dataset(file)
    return ds.rename(LEV = 'depth')

def read_glorys(file):
    ds = xr.open_dataset(file, decode_times=True)
    return ds.rename(thetao = 'theta', so = 'salt', zos = 'ssh')
    
def read_duacs(file):
    ds = xr.open_dataset(file)
    return ds.rename(sla = 'ssh')

g = 9.81

def interp_xy(ds, x=None, y=None, method = 'linear'):
    """
    frontend for xarray spatial interpolation

    Parameters
    ----------
    ds : xarray
        input xarray grid.
    x : float, optional
        longitude in deg (but depends on ds units). The default is None.
    y : float, optional
        latitude in deg (but depends on ds units). The default is None.
    method : str, optional
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


def interp_time(ds, time, method = 'linear'):
    return ds.interp(time = time, method = method)

def interp_z(ds, depth, method = 'linear'):
    return ds.interp(depth = depth, method = method)

def get_bottom_depth(ds):
    land = np.isnan(ds.theta[0, 0, :, :])
    isnan = np.isnan(ds.theta[0, 0:, :, :])
    ind_nan = np.count_nonzero(~isnan, axis = 0)
    bottom_depth = ds.depth.values[ind_nan]
    bottom_depth[land.values] = np.nan
    return bottom_depth

def extract_profile(ds, x, y, method = 'linear', zbottom = None):

    ds = interp_xy(ds, x, y, method)
    
    if not isinstance(zbottom, type(None)):
        if isinstance(zbottom, (int, float)):
            valid = ~np.isnan(ds.theta[0])
            max_depth = ds.depth[valid][-1]
            if zbottom > max_depth.values :
                print(f'{zbottom} m is deeper than the maximum depth layer of the profile : {max_depth.values}')
                return
            k = len((zbottom - ds.depth)[(zbottom - ds.depth)>0]) - 1
            new_depth = np.hstack([ds.depth[:k],zbottom])
        elif isinstance(zbottom, (np.ndarray, list)):
            new_depth = zbottom
        else:
            print(f'Invalid zbottom parameter type : {type(zbottom)}. Must be int or float depth (in m).')
        ds = interp_z(ds, new_depth, method = 'linear')

    return ds

def compute_dens_profile(profile):
    ds = profile
    if np.ndim(ds.theta) != 2:
        print(f'Bad dataset dimension : {np.ndim(hycom.theta)}. The function must be applied on a profile.')
        return
    lat = ds.latitude.mean()
    pres = seawater.pres(ds.depth, lat = lat)
    dens = seawater.dens(ds.salt, ds.theta, pres)
    return dens

def compute_steric(profile, integration = 'forward'):
    ds = profile
    lat = ds.latitude.mean()
    dens = compute_dens_profile(ds)
    if integration == 'forward':
        steric = pd.DataFrame(g*np.nansum(np.diff(ds.depth) * dens[:, :-1],
                                          axis = 1)/10000, index = ds.time)[0]
    elif integration == 'backward':
        steric = pd.DataFrame(g*np.nansum(np.diff(ds.depth) * dens[:, 1:], 
                                          axis = 1)/10000, index = ds.time)[0]
    elif integration == 'averaged':
        dens = np.vstack([np.sum(dens[:, i:i+2], axis = 1)/2 for i in range(len(ds.depth)-1)]).T
        steric = pd.DataFrame(g*np.nansum(np.diff(ds.depth) * dens, axis = 1)/10000, index = ds.time)[0]
    return steric - steric.median()

def compute_phibot(profile, ssh=None, integration = 'forward',
                   rho = 10.35, remove_median = True):
    """
    phibot = Ocean hydrostatic bottom pressure anomaly
    https://cmr.earthdata.nasa.gov/search/concepts/V2028471168-POCLOUD.html
    
    
    PHIBOT		Bottom Pressure Pot. Anomaly (p/rhonil, m^2/s^2)
                To convert to m, divide by g (g=9.81 m/s^2)
		PHIBOT is the anomaly relative to Depth * rhonil * g
		The absolute bottom pressure in Pa is:
		Depth * rhonil * g + PHIBOT * rhonil (rhonil=1027.5 kg/m^3)
    http://apdrc.soest.hawaii.edu/doc/Readme_ecco2_cube92
        
    """
    if isinstance(ssh, type(None)):
        ssh = profile.ssh
    lat = profile.latitude.mean()
    dens = compute_dens_profile(profile)
    phibot = pd.DataFrame(index = profile.time)
    phibot['subsurface'] = compute_steric(profile, integration)
    phibot['surface'] = g*rho*ssh/100
    phibot['phibot'] = phibot.subsurface + phibot.surface 
    if remove_median:
        phibot -= phibot.median(axis = 0)
        
    return phibot
    
def compute_spectrogram(df, max_period = 45, nchunks = 3600*6):
    from scipy.signal import spectrogram

    f, t, Sxx = spectrogram(df.values, nperseg=nchunks, scaling = 'density')
    T = (1/f[1:])
    Sxx = Sxx[1:]
    i = np.argmin(np.abs(T - max_period))
    T = T[i:]
    Sxx = Sxx[i:]

    t = pd.DatetimeIndex([df.index[0] + timedelta(hours = t) for t in t/3600])
    return T, t, Sxx
    
import scipy.optimize as optimize
def log_linear(x, A1, A2, A3, A4, **kwargs):
    #x0 = min(x)
    return A1 * np.log(A2 * (1 + x))  +A4 +  A3 * (x)

def log(x, A1, A2, A3, **kwargs):
    x0 = min(x)
    return A1 * np.log(A2 * (1 + x -x0))  + A3  # +A3 * (x - x0)

def exp(x, a, b, c, d, **kwargs):
    #x0 = min(x)
    return a * np.exp(-b * (x )) + c #+ d * (x - x0)

def linear(x, A1, A2, A3, A4, **kwargs):
    x0 = min(x)
    return A1 * (x - x0) + A2

def exp_linear(x, a, b, c, d, **kwargs):
    x0 = min(x)
    return a * np.exp(-1/b * (x - x0)) + c + d * (x-x0)

def fit_model(data, t_fit, model = 'log_linear', offset_last = True, maxfev = 1000, pn = 2):
    
    t_num = (data.index.to_julian_date() - data.index.to_julian_date()[0]).values
    t_fit_num = (t_fit.to_julian_date() - t_fit.to_julian_date()[0]).values

    ref = t_num[-1]
    t_num /= ref
    t_fit_num/=ref
    
    if model == 'poly':

        z = np.polyfit(t_num, data.values, pn)
        p = np.poly1d(z)

        fit = p(t_fit_num)
        ind_min = np.argmin(fit)
        fit[ind_min:] = fit[ind_min]
        
    if model == 'log_linear':
        popt, pcov = optimize.curve_fit(log_linear, t_num, data.values, maxfev = maxfev)
        fit = log_linear(t_fit_num, *popt)
    elif model == 'log':
        popt, pcov = optimize.curve_fit(log, t_num, data.values, maxfev = maxfev)
        fit = log(t_fit_num, *popt)
    elif model == 'exp':
        popt, pcov = optimize.curve_fit(exp, t_num, data.values, maxfev = maxfev)
        fit = exp(t_fit_num, *popt)
    elif model == 'exp_linear':
        popt, pcov = optimize.curve_fit(exp_linear, t_num, data.values, maxfev = maxfev)
        fit = exp_linear(t_fit_num, *popt)
    elif model == 'linear':
        popt, pcov = optimize.curve_fit(linear, t_num, data.values, maxfev = maxfev)
        fit = linear(t_fit_num, *popt)
    return pd.DataFrame(fit, index = t_fit)[0]