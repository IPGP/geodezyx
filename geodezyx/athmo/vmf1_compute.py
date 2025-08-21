#!/usr/bin/env python3
"""
Python translation of the VMF1G computation Fortran modules.

This module provides atmospheric delay computation functions for geodetic applications
using Vienna Mapping Function 1 Gridded (VMF1G) data.

Translated from Fortran modules:
- get_vmf1g_values.f90
- compute_vmf1g_values.f90
- get_orography.f90
- sd_vmf1_grid_elevation.f90
- sd_gmf.f90
- sd_gpt.f90
- jour_julien_ff.f90
- calend1.f90

Author: Translated to Python
Original Fortran authors: FF, FR, Johannes Boehm
"""

import os
import sys
import math
import numpy as np
from typing import Dict, Tuple, Union
from geodezyx import conv
import datetime as dt

def get_orography(orography_file: str, lat: float, lon: float) -> np.ndarray:
    """
    Get orographic heights for the 4 grid nodes surrounding a location.

    Parameters:
    -----------
    orography_file : str
        File containing orography data
    lat, lon : float
        Latitude and longitude in degrees

    Returns:
    --------
    np.ndarray : Orographic heights for 4 grid nodes
    """
    # Initialize height array
    h_tmp = np.zeros((91, 145), dtype=int)

    # Read orography file
    fic1 = os.path.join(orography_file, 'orography_ell.dat')

    with open(fic1, 'r') as f:
        # Skip header
        f.readline()

        # Read data
        i = 1
        j = 1
        for line in f:
            if i <= 14:
                values = list(map(int, line.strip().split()))
                start_idx = (i-1) * 10
                end_idx = start_idx + len(values)
                h_tmp[91-j, start_idx:end_idx] = values
                i += 1
            else:
                values = list(map(int, line.strip().split()))
                start_idx = (i-1) * 10
                end_idx = start_idx + len(values)
                h_tmp[91-j, start_idx:end_idx] = values
                i = 1
                j += 1
                if j == 92:
                    break

    # Calculate grid coordinates
    latg4 = lat - (lat % 2.0)
    if lat < 0:
        latg4 = latg4 - 2.0
    latg1 = latg4 + 2.0

    long4 = lon - (lon % 2.5)
    long3 = long4 + 2.5

    if latg1 > 90.0:  # North pole overflow
        latg1 = 90.0

    if long3 > 357.5:  # Greenwich meridian overflow
        long3 = 0.0

    latg2 = latg1
    latg3 = latg4
    long1 = long4
    long2 = long3

    # Calculate indices
    ind_lat1 = int((latg1 - (-90) + 2.0) / 2.0)
    ind_lat2 = int((latg2 - (-90) + 2.0) / 2.0)
    ind_lat3 = int((latg3 - (-90) + 2.0) / 2.0)
    ind_lat4 = int((latg4 - (-90) + 2.0) / 2.0)

    ind_lon1 = int((long1 + 2.5) / 2.5)
    ind_lon2 = int((long2 + 2.5) / 2.5)
    ind_lon3 = int((long3 + 2.5) / 2.5)
    ind_lon4 = int((long4 + 2.5) / 2.5)

    # Extract orographic heights
    oro_zhdg = np.array([
        float(h_tmp[ind_lat1-1, ind_lon1-1]),  # Convert to 0-based indexing
        float(h_tmp[ind_lat2-1, ind_lon2-1]),
        float(h_tmp[ind_lat3-1, ind_lon3-1]),
        float(h_tmp[ind_lat4-1, ind_lon4-1])
    ])

    return oro_zhdg


def sd_vmf1_grid_elevation(ah: float, aw: float, dmjd: float, dlat: float,
                          ht: float, zd: float) -> Tuple[float, float]:
    """
    Determine VMF1 (Vienna Mapping Functions 1) with height correction.

    Parameters:
    -----------
    ah : float
        Hydrostatic coefficient a
    aw : float
        Wet coefficient a
    dmjd : float
        Modified Julian date
    dlat : float
        Latitude in radians
    ht : float
        Ellipsoidal height in meters
    zd : float
        Zenith distance in radians

    Returns:
    --------
    tuple : (vmf1h, vmf1w) - hydrostatic and wet mapping functions
    """
    pi = math.pi

    # Reference day is 28 January (consistent with Niell 1996)
    doy = dmjd + 33282.0 - 44239.0 + 1.0 - 28.0

    # Hydrostatic mapping function
    bh = 0.0029
    c0h = 0.062

    if dlat < 0.0:  # Southern hemisphere
        phh = pi
        c11h = 0.007
        c10h = 0.002
    else:  # Northern hemisphere
        phh = 0
        c11h = 0.005
        c10h = 0.001

    ch = c0h + ((math.cos(doy/365.25*2.0*pi + phh) + 1) * c11h/2.0 + c10h) * (1.0 - math.cos(dlat))

    sine = math.sin(pi/2.0 - zd)
    beta = bh / (sine + ch)
    gamma = ah / (sine + beta)
    topcon = (1.0 + ah/(1 + bh/(1.0 + ch)))
    vmf1h = topcon / (sine + gamma)

    # Height correction (Niell, 1996)
    a_ht = 2.53e-5
    b_ht = 5.49e-3
    c_ht = 1.14e-3
    hs_km = ht / 1000.0

    beta = b_ht / (sine + c_ht)
    gamma = a_ht / (sine + beta)
    topcon = (1.0 + a_ht/(1.0 + b_ht/(1.0 + c_ht)))
    ht_corr_coef = 1.0/sine - topcon/(sine + gamma)
    ht_corr = ht_corr_coef * hs_km
    vmf1h = vmf1h + ht_corr

    # Wet mapping function
    bw = 0.00146
    cw = 0.04391
    beta = bw / (sine + cw)
    gamma = aw / (sine + beta)
    topcon = (1.0 + aw/(1.0 + bw/(1.0 + cw)))
    vmf1w = topcon / (sine + gamma)

    return vmf1h, vmf1w


def sd_gpt_simple(dmjd: float, dlat: float, dlon: float, dhgt: float) -> Tuple[float, float, float]:
    """
    Determine Global Pressure and Temperature based on Spherical Harmonics.

    Parameters:
    -----------
    dmjd : float
        Modified Julian date
    dlat : float
        Latitude in radians
    dlon : float
        Longitude in radians
    dhgt : float
        Ellipsoidal height in m

    Returns:
    --------
    tuple : (pressure in hPa, temperature in Celsius, geoid undulation in m)
    """
    pi = math.pi

    # Reference day is 28 January
    doy = dmjd + 33282 - 44239 + 1 - 28

    # Spherical harmonic coefficients (truncated for brevity - full arrays needed for complete implementation)
    # These would need to be the complete arrays from the Fortran code
    # For now, using simplified GPT model

    # Simplified GPT implementation - pressure at sea level
    pres0 = 1013.25  # Standard atmospheric pressure
    temp0 = 15.0     # Standard temperature

    # Height correction for pressure
    h_corr = dhgt / 1000.0  # Convert to km
    pres = pres0 * math.exp(-h_corr / 8.5)  # Simplified barometric formula

    # Temperature lapse rate
    temp = temp0 - 6.5 * h_corr  # Standard lapse rate

    # Simplified geoid undulation (would need full spherical harmonic implementation)
    undu = 0.0

    return pres, temp, undu


def sd_gmf_simple(dmjd: float, dlat: float, dlon: float, dhgt: float, zd: float) -> Tuple[float, float]:
    """
    Determine Global Mapping Functions (GMF).

    Parameters:
    -----------
    dmjd : float
        Modified Julian date
    dlat : float
        Latitude in radians
    dlon : float
        Longitude in radians
    dhgt : float
        Height in m
    zd : float
        Zenith distance in radians

    Returns:
    --------
    tuple : (gmfh, gmfw) - hydrostatic and wet mapping functions
    """
    pi = math.pi

    # Reference day is 28 January
    doy = dmjd + 33282 - 44239 + 1 - 28

    # Simplified GMF implementation
    # This would need the full spherical harmonic coefficients from the Fortran code

    # Basic mapping function calculation
    sine = math.sin(pi/2.0 - zd)

    # Simplified coefficients (full implementation would use spherical harmonics)
    ah = 0.00127
    bh = 0.0029
    ch = 0.062

    aw = 0.00058
    bw = 0.00146
    cw = 0.04391

    # Hydrostatic mapping function
    beta_h = bh / (sine + ch)
    gamma_h = ah / (sine + beta_h)
    topcon_h = (1.0 + ah/(1 + bh/(1.0 + ch)))
    gmfh = topcon_h / (sine + gamma_h)

    # Wet mapping function
    beta_w = bw / (sine + cw)
    gamma_w = aw / (sine + beta_w)
    topcon_w = (1.0 + aw/(1.0 + bw/(1.0 + cw)))
    gmfw = topcon_w / (sine + gamma_w)

    return gmfh, gmfw


def sd_gmf(dmjd, dlat, dlon, dhgt, zd):
    """
    Python translation of sd_gmf subroutine
    Determines Global Mapping Functions GMF

    Parameters:
    -----------
    dmjd : float
        modified julian date
    dlat : float
        latitude in radians
    dlon : float
        longitude in radians
    dhgt : float
        height in m
    zd : float
        zenith distance in radians

    Returns:
    --------
    gmfh : float
        hydrostatic mapping function
    gmfw : float
        wet mapping function
    """

    pi = 3.14159265359

    # reference day is 28 January
    doy = dmjd + 33282 - 44239.0 + 1 - 28

    # Coefficients for GMF (complete arrays)
    ah_mean = np.array(
        [
            1.2517e02,
            8.503e-01,
            6.936e-02,
            -6.760e00,
            1.771e-01,
            1.130e-02,
            5.963e-01,
            1.808e-02,
            2.801e-03,
            -1.414e-03,
            -1.212e00,
            9.300e-02,
            3.683e-03,
            1.095e-03,
            4.671e-05,
            3.959e-01,
            -3.867e-02,
            5.413e-03,
            -5.289e-04,
            3.229e-04,
            2.067e-05,
            3.000e-01,
            2.031e-02,
            5.900e-03,
            4.573e-04,
            -7.619e-05,
            2.327e-06,
            3.845e-06,
            1.182e-01,
            1.158e-02,
            5.445e-03,
            6.219e-05,
            4.204e-06,
            -2.093e-06,
            1.540e-07,
            -4.280e-08,
            -4.751e-01,
            -3.490e-02,
            1.758e-03,
            4.019e-04,
            -2.799e-06,
            -1.287e-06,
            5.468e-07,
            7.580e-08,
            -6.300e-09,
            -1.160e-01,
            8.301e-03,
            8.771e-04,
            9.955e-05,
            -1.718e-06,
            -2.012e-06,
            1.170e-08,
            1.790e-08,
            -1.300e-09,
            1.000e-10,
        ]
    )

    bh_mean = np.array(
        [
            0.000e00,
            0.000e00,
            3.249e-02,
            0.000e00,
            3.324e-02,
            1.850e-02,
            0.000e00,
            -1.115e-01,
            2.519e-02,
            4.923e-03,
            0.000e00,
            2.737e-02,
            1.595e-02,
            -7.332e-04,
            1.933e-04,
            0.000e00,
            -4.796e-02,
            6.381e-03,
            -1.599e-04,
            -3.685e-04,
            1.815e-05,
            0.000e00,
            7.033e-02,
            2.426e-03,
            -1.111e-03,
            -1.357e-04,
            -7.828e-06,
            2.547e-06,
            0.000e00,
            5.779e-03,
            3.133e-03,
            -5.312e-04,
            -2.028e-05,
            2.323e-07,
            -9.100e-08,
            -1.650e-08,
            0.000e00,
            3.688e-02,
            -8.638e-04,
            -8.514e-05,
            -2.828e-05,
            5.403e-07,
            4.390e-07,
            1.350e-08,
            1.800e-09,
            0.000e00,
            -2.736e-02,
            -2.977e-04,
            8.113e-05,
            2.329e-07,
            8.451e-07,
            4.490e-08,
            -8.100e-09,
            -1.500e-09,
            2.000e-10,
        ]
    )

    ah_amp = np.array(
        [
            -2.738e-01,
            -2.837e00,
            1.298e-02,
            -3.588e-01,
            2.413e-02,
            3.427e-02,
            -7.624e-01,
            7.272e-02,
            2.160e-02,
            -3.385e-03,
            4.424e-01,
            3.722e-02,
            2.195e-02,
            -1.503e-03,
            2.426e-04,
            3.013e-01,
            5.762e-02,
            1.019e-02,
            -4.476e-04,
            6.790e-05,
            3.227e-05,
            3.123e-01,
            -3.535e-02,
            4.840e-03,
            3.025e-06,
            -4.363e-05,
            2.854e-07,
            -1.286e-06,
            -6.725e-01,
            -3.730e-02,
            8.964e-04,
            1.399e-04,
            -3.990e-06,
            7.431e-06,
            -2.796e-07,
            -1.601e-07,
            4.068e-02,
            -1.352e-02,
            7.282e-04,
            9.594e-05,
            2.070e-06,
            -9.620e-08,
            -2.742e-07,
            -6.370e-08,
            -6.300e-09,
            8.625e-02,
            -5.971e-03,
            4.705e-04,
            2.335e-05,
            4.226e-06,
            2.475e-07,
            -8.850e-08,
            -3.600e-08,
            -2.900e-09,
            0.000e00,
        ]
    )

    bh_amp = np.array(
        [
            0.000e00,
            0.000e00,
            -1.136e-01,
            0.000e00,
            -1.868e-01,
            -1.399e-02,
            0.000e00,
            -1.043e-01,
            1.175e-02,
            -2.240e-03,
            0.000e00,
            -3.222e-02,
            1.333e-02,
            -2.647e-03,
            -2.316e-05,
            0.000e00,
            5.339e-02,
            1.107e-02,
            -3.116e-03,
            -1.079e-04,
            -1.299e-05,
            0.000e00,
            4.861e-03,
            8.891e-03,
            -6.448e-04,
            -1.279e-05,
            6.358e-06,
            -1.417e-07,
            0.000e00,
            3.041e-02,
            1.150e-03,
            -8.743e-04,
            -2.781e-05,
            6.367e-07,
            -1.140e-08,
            -4.200e-08,
            0.000e00,
            -2.982e-02,
            -3.000e-03,
            1.394e-05,
            -3.290e-05,
            -1.705e-07,
            7.440e-08,
            2.720e-08,
            -6.600e-09,
            0.000e00,
            1.236e-02,
            -9.981e-04,
            -3.792e-05,
            -1.355e-05,
            1.162e-06,
            -1.789e-07,
            1.470e-08,
            -2.400e-09,
            -4.000e-10,
        ]
    )

    aw_mean = np.array(
        [
            5.640e01,
            1.555e00,
            -1.011e00,
            -3.975e00,
            3.171e-02,
            1.065e-01,
            6.175e-01,
            1.376e-01,
            4.229e-02,
            3.028e-03,
            1.688e00,
            -1.692e-01,
            5.478e-02,
            2.473e-02,
            6.059e-04,
            2.278e00,
            6.614e-03,
            -3.505e-04,
            -6.697e-03,
            8.402e-04,
            7.033e-04,
            -3.236e00,
            2.184e-01,
            -4.611e-02,
            -1.613e-02,
            -1.604e-03,
            5.420e-05,
            7.922e-05,
            -2.711e-01,
            -4.406e-01,
            -3.376e-02,
            -2.801e-03,
            -4.090e-04,
            -2.056e-05,
            6.894e-06,
            2.317e-06,
            1.941e00,
            -2.562e-01,
            1.598e-02,
            5.449e-03,
            3.544e-04,
            1.148e-05,
            7.503e-06,
            -5.667e-07,
            -3.660e-08,
            8.683e-01,
            -5.931e-02,
            -1.864e-03,
            -1.277e-04,
            2.029e-04,
            1.269e-05,
            1.629e-06,
            9.660e-08,
            -1.015e-07,
            -5.000e-10,
        ]
    )

    bw_mean = np.array(
        [
            0.000e00,
            0.000e00,
            2.592e-01,
            0.000e00,
            2.974e-02,
            -5.471e-01,
            0.000e00,
            -5.926e-01,
            -1.030e-01,
            -1.567e-02,
            0.000e00,
            1.710e-01,
            9.025e-02,
            2.689e-02,
            2.243e-03,
            0.000e00,
            3.439e-01,
            2.402e-02,
            5.410e-03,
            1.601e-03,
            9.669e-05,
            0.000e00,
            9.502e-02,
            -3.063e-02,
            -1.055e-03,
            -1.067e-04,
            -1.130e-04,
            2.124e-05,
            0.000e00,
            -3.129e-01,
            8.463e-03,
            2.253e-04,
            7.413e-05,
            -9.376e-05,
            -1.606e-06,
            2.060e-06,
            0.000e00,
            2.739e-01,
            1.167e-03,
            -2.246e-05,
            -1.287e-04,
            -2.438e-05,
            -7.561e-07,
            1.158e-06,
            4.950e-08,
            0.000e00,
            -1.344e-01,
            5.342e-03,
            3.775e-04,
            -6.756e-05,
            -1.686e-06,
            -1.184e-06,
            2.768e-07,
            2.730e-08,
            5.700e-09,
        ]
    )

    aw_amp = np.array(
        [
            1.023e-01,
            -2.695e00,
            3.417e-01,
            -1.405e-01,
            3.175e-01,
            2.116e-01,
            3.536e00,
            -1.505e-01,
            -1.660e-02,
            2.967e-02,
            3.819e-01,
            -1.695e-01,
            -7.444e-02,
            7.409e-03,
            -6.262e-03,
            -1.836e00,
            -1.759e-02,
            -6.256e-02,
            -2.371e-03,
            7.947e-04,
            1.501e-04,
            -8.603e-01,
            -1.360e-01,
            -3.629e-02,
            -3.706e-03,
            -2.976e-04,
            1.857e-05,
            3.021e-05,
            2.248e00,
            -1.178e-01,
            1.255e-02,
            1.134e-03,
            -2.161e-04,
            -5.817e-06,
            8.836e-07,
            -1.769e-07,
            7.313e-01,
            -1.188e-01,
            1.145e-02,
            1.011e-03,
            1.083e-04,
            2.570e-06,
            -2.140e-06,
            -5.710e-08,
            2.000e-08,
            -1.632e00,
            -6.948e-03,
            -3.893e-03,
            8.592e-04,
            7.577e-05,
            4.539e-06,
            -3.852e-07,
            -2.213e-07,
            -1.370e-08,
            5.800e-09,
        ]
    )

    bw_amp = np.array(
        [
            0.000e00,
            0.000e00,
            -8.865e-02,
            0.000e00,
            -4.309e-01,
            6.340e-02,
            0.000e00,
            1.162e-01,
            6.176e-02,
            -4.234e-03,
            0.000e00,
            2.530e-01,
            4.017e-02,
            -6.204e-03,
            4.977e-03,
            0.000e00,
            -1.737e-01,
            -5.638e-03,
            1.488e-04,
            4.857e-04,
            -1.809e-04,
            0.000e00,
            -1.514e-01,
            -1.685e-02,
            5.333e-03,
            -7.611e-05,
            2.394e-05,
            8.195e-06,
            0.000e00,
            9.326e-02,
            -1.275e-02,
            -3.071e-04,
            5.374e-05,
            -3.391e-05,
            -7.436e-06,
            6.747e-07,
            0.000e00,
            -8.637e-02,
            -3.807e-03,
            -6.833e-04,
            -3.861e-05,
            -2.268e-05,
            1.454e-06,
            3.860e-07,
            -1.068e-07,
            0.000e00,
            -2.658e-02,
            -1.947e-03,
            7.131e-04,
            -3.506e-05,
            1.885e-07,
            5.792e-07,
            3.990e-08,
            2.000e-08,
            -5.700e-09,
        ]
    )

    # Degree n and order m
    nmax = 9
    mmax = 9

    # Unit vector
    x = math.cos(dlat) * math.cos(dlon)
    y = math.cos(dlat) * math.sin(dlon)
    z = math.sin(dlat)

    # Legendre polynomials
    V = np.zeros((10, 10))
    W = np.zeros((10, 10))

    V[0, 0] = 1
    W[0, 0] = 0
    V[1, 0] = z * V[0, 0]
    W[1, 0] = 0

    for n in range(2, nmax + 1):
        V[n, 0] = ((2 * n - 1) * z * V[n - 1, 0] - (n - 1) * V[n - 2, 0]) / n
        W[n, 0] = 0.0

    for m in range(1, nmax + 1):
        V[m, m] = (2 * m - 1) * (x * V[m - 1, m - 1] - y * W[m - 1, m - 1])
        W[m, m] = (2 * m - 1) * (x * W[m - 1, m - 1] + y * V[m - 1, m - 1])
        if m < nmax:
            V[m + 1, m] = (2 * m + 1) * z * V[m, m]
            W[m + 1, m] = (2 * m + 1) * z * W[m, m]
        for n in range(m + 2, nmax + 1):
            V[n, m] = ((2 * n - 1) * z * V[n - 1, m] - (n + m - 1) * V[n - 2, m]) / (
                n - m
            )
            W[n, m] = ((2 * n - 1) * z * W[n - 1, m] - (n + m - 1) * W[n - 2, m]) / (
                n - m
            )

    # Hydrostatic
    bh = 0.0029
    c0h = 0.062
    if dlat < 0:  # southern hemisphere
        phh = pi
        c11h = 0.007
        c10h = 0.002
    else:  # northern hemisphere
        phh = 0
        c11h = 0.005
        c10h = 0.001

    ch = c0h + ((math.cos((doy * (2 * pi)) / 365.25 + phh) + 1) * c11h / 2 + c10h) * (
        1 - math.cos(dlat)
    )

    ahm = 0
    aha = 0
    i = 0
    for n in range(nmax + 1):
        for m in range(n + 1):
            ahm = ahm + (ah_mean[i] * V[n, m] + bh_mean[i] * W[n, m])
            aha = aha + (ah_amp[i] * V[n, m] + bh_amp[i] * W[n, m])
            i = i + 1

    ah = (ahm + aha * math.cos(doy / 365.25 * 2 * pi)) * 1e-5

    sine = math.sin(pi / 2.0 - zd)
    beta = bh / (sine + ch)
    gamma = ah / (sine + beta)
    topcon = 1 + ah / (1 + bh / (1 + ch))
    gmfh = topcon / (sine + gamma)

    # Height correction for hydrostatic mapping function from Niell (1996)
    a_ht = 2.53e-5
    b_ht = 5.49e-3
    c_ht = 1.14e-3
    hs_km = dhgt / 1000.0

    beta = b_ht / (sine + c_ht)
    gamma = a_ht / (sine + beta)
    topcon = 1 + a_ht / (1 + b_ht / (1 + c_ht))
    ht_corr_coef = 1.0 / sine - topcon / (sine + gamma)
    ht_corr = ht_corr_coef * hs_km
    gmfh = gmfh + ht_corr

    # Wet
    bw = 0.00146
    cw = 0.04391

    awm = 0
    awa = 0
    i = 0
    for n in range(nmax + 1):
        for m in range(n + 1):
            awm = awm + (aw_mean[i] * V[n, m] + bw_mean[i] * W[n, m])
            awa = awa + (aw_amp[i] * V[n, m] + bw_amp[i] * W[n, m])
            i = i + 1

    aw = (awm + awa * math.cos(doy * (2 * pi) / 365.25)) * 1e-5

    beta = bw / (sine + cw)
    gamma = aw / (sine + beta)
    topcon = 1 + aw / (1 + bw / (1 + cw))
    gmfw = topcon / (sine + gamma)

    return gmfh, gmfw


def sd_gpt(dmjd, dlat, dlon, dhgt):
    """
    Python translation of sd_gpt subroutine
    Determines Global Pressure and Temperature

    Parameters:
    -----------
    dmjd : float
        modified julian date
    dlat : float
        latitude in radians
    dlon : float
        longitude in radians
    dhgt : float
        ellipsoidal height in m

    Returns:
    --------
    pres : float
        pressure in hPa
    temp : float
        temperature in Celsius
    undu : float
        Geoid undulation in m
    """

    pi = 3.14159265359

    # reference day is 28 January
    doy = dmjd + 33282 - 44239 + 1 - 28

    # Geoid coefficients
    a_geoid = np.array(
        [
            -5.6195e-01,
            -6.0794e-02,
            -2.0125e-01,
            -6.4180e-02,
            -3.6997e-02,
            +1.0098e01,
            +1.6436e01,
            +1.4065e01,
            +1.9881e00,
            +6.4414e-01,
            -4.7482e00,
            -3.2290e00,
            +5.0652e-01,
            +3.8279e-01,
            -2.6646e-02,
            +1.7224e00,
            -2.7970e-01,
            +6.8177e-01,
            -9.6658e-02,
            -1.5113e-02,
            +2.9206e-03,
            -3.4621e00,
            -3.8198e-01,
            +3.2306e-02,
            +6.9915e-03,
            -2.3068e-03,
            -1.3548e-03,
            +4.7324e-06,
            +2.3527e00,
            +1.2985e00,
            +2.1232e-01,
            +2.2571e-02,
            -3.7855e-03,
            +2.9449e-05,
            -1.6265e-04,
            +1.1711e-07,
            +1.6732e00,
            +1.9858e-01,
            +2.3975e-02,
            -9.0013e-04,
            -2.2475e-03,
            -3.3095e-05,
            -1.2040e-05,
            +2.2010e-06,
            -1.0083e-06,
            +8.6297e-01,
            +5.8231e-01,
            +2.0545e-02,
            -7.8110e-03,
            -1.4085e-04,
            -8.8459e-06,
            +5.7256e-06,
            -1.5068e-06,
            +4.0095e-07,
            -2.4185e-08,
        ]
    )

    b_geoid = np.array(
        [
            +0.0000e00,
            +0.0000e00,
            -6.5993e-02,
            +0.0000e00,
            +6.5364e-02,
            -5.8320e00,
            +0.0000e00,
            +1.6961e00,
            -1.3557e00,
            +1.2694e00,
            +0.0000e00,
            -2.9310e00,
            +9.4805e-01,
            -7.6243e-02,
            +4.1076e-02,
            +0.0000e00,
            -5.1808e-01,
            -3.4583e-01,
            -4.3632e-02,
            +2.2101e-03,
            -1.0663e-02,
            +0.0000e00,
            +1.0927e-01,
            -2.9463e-01,
            +1.4371e-03,
            -1.1452e-02,
            -2.8156e-03,
            -3.5330e-04,
            +0.0000e00,
            +4.4049e-01,
            +5.5653e-02,
            -2.0396e-02,
            -1.7312e-03,
            +3.5805e-05,
            +7.2682e-05,
            +2.2535e-06,
            +0.0000e00,
            +1.9502e-02,
            +2.7919e-02,
            -8.1812e-03,
            +4.4540e-04,
            +8.8663e-05,
            +5.5596e-05,
            +2.4826e-06,
            +1.0279e-06,
            +0.0000e00,
            +6.0529e-02,
            -3.5824e-02,
            -5.1367e-03,
            +3.0119e-05,
            -2.9911e-05,
            +1.9844e-05,
            -1.2349e-06,
            -7.6756e-09,
            +5.0100e-08,
        ]
    )

    # Pressure coefficients
    ap_mean = np.array(
        [
            +1.0108e03,
            +8.4886e00,
            +1.4799e00,
            -1.3897e01,
            +3.7516e-03,
            -1.4936e-01,
            +1.2232e01,
            -7.6615e-01,
            -6.7699e-02,
            +8.1002e-03,
            -1.5874e01,
            +3.6614e-01,
            -6.7807e-02,
            -3.6309e-03,
            +5.9966e-04,
            +4.8163e00,
            -3.7363e-01,
            -7.2071e-02,
            +1.9998e-03,
            -6.2385e-04,
            -3.7916e-04,
            +4.7609e00,
            -3.9534e-01,
            +8.6667e-03,
            +1.1569e-02,
            +1.1441e-03,
            -1.4193e-04,
            -8.5723e-05,
            +6.5008e-01,
            -5.0889e-01,
            -1.5754e-02,
            -2.8305e-03,
            +5.7458e-04,
            +3.2577e-05,
            -9.6052e-06,
            -2.7974e-06,
            +1.3530e00,
            -2.7271e-01,
            -3.0276e-04,
            +3.6286e-03,
            -2.0398e-04,
            +1.5846e-05,
            -7.7787e-06,
            +1.1210e-06,
            +9.9020e-08,
            +5.5046e-01,
            -2.7312e-01,
            +3.2532e-03,
            -2.4277e-03,
            +1.1596e-04,
            +2.6421e-07,
            -1.3263e-06,
            +2.7322e-07,
            +1.4058e-07,
            +4.9414e-09,
        ]
    )

    bp_mean = np.array(
        [
            +0.0000e00,
            +0.0000e00,
            -1.2878e00,
            +0.0000e00,
            +7.0444e-01,
            +3.3222e-01,
            +0.0000e00,
            -2.9636e-01,
            +7.2248e-03,
            +7.9655e-03,
            +0.0000e00,
            +1.0854e00,
            +1.1145e-02,
            -3.6513e-02,
            +3.1527e-03,
            +0.0000e00,
            -4.8434e-01,
            +5.2023e-02,
            -1.3091e-02,
            +1.8515e-03,
            +1.5422e-04,
            +0.0000e00,
            +6.8298e-01,
            +2.5261e-03,
            -9.9703e-04,
            -1.0829e-03,
            +1.7688e-04,
            -3.1418e-05,
            +0.0000e00,
            -3.7018e-01,
            +4.3234e-02,
            +7.2559e-03,
            +3.1516e-04,
            +2.0024e-05,
            -8.0581e-06,
            -2.3653e-06,
            +0.0000e00,
            +1.0298e-01,
            -1.5086e-02,
            +5.6186e-03,
            +3.2613e-05,
            +4.0567e-05,
            -1.3925e-06,
            -3.6219e-07,
            -2.0176e-08,
            +0.0000e00,
            -1.8364e-01,
            +1.8508e-02,
            +7.5016e-04,
            -9.6139e-05,
            -3.1995e-06,
            +1.3868e-07,
            -1.9486e-07,
            +3.0165e-10,
            -6.4376e-10,
        ]
    )

    ap_amp = np.array(
        [
            -1.0444e-01,
            +1.6618e-01,
            -6.3974e-02,
            +1.0922e00,
            +5.7472e-01,
            -3.0277e-01,
            -3.5087e00,
            +7.1264e-03,
            -1.4030e-01,
            +3.7050e-02,
            +4.0208e-01,
            -3.0431e-01,
            -1.3292e-01,
            +4.6746e-03,
            -1.5902e-04,
            +2.8624e00,
            -3.9315e-01,
            -6.4371e-02,
            +1.6444e-02,
            -2.3403e-03,
            +4.2127e-05,
            +1.9945e00,
            -6.0907e-01,
            -3.5386e-02,
            -1.0910e-03,
            -1.2799e-04,
            +4.0970e-05,
            +2.2131e-05,
            -5.3292e-01,
            -2.9765e-01,
            -3.2877e-02,
            +1.7691e-03,
            +5.9692e-05,
            +3.1725e-05,
            +2.0741e-05,
            -3.7622e-07,
            +2.6372e00,
            -3.1165e-01,
            +1.6439e-02,
            +2.1633e-04,
            +1.7485e-04,
            +2.1587e-05,
            +6.1064e-06,
            -1.3755e-08,
            -7.8748e-08,
            -5.9152e-01,
            -1.7676e-01,
            +8.1807e-03,
            +1.0445e-03,
            +2.3432e-04,
            +9.3421e-06,
            +2.8104e-06,
            -1.5788e-07,
            -3.0648e-08,
            +2.6421e-10,
        ]
    )

    bp_amp = np.array(
        [
            +0.0000e00,
            +0.0000e00,
            +9.3340e-01,
            +0.0000e00,
            +8.2346e-01,
            +2.2082e-01,
            +0.0000e00,
            +9.6177e-01,
            -1.5650e-02,
            +1.2708e-03,
            +0.0000e00,
            -3.9913e-01,
            +2.8020e-02,
            +2.8334e-02,
            +8.5980e-04,
            +0.0000e00,
            +3.0545e-01,
            -2.1691e-02,
            +6.4067e-04,
            -3.6528e-05,
            -1.1166e-04,
            +0.0000e00,
            -7.6974e-02,
            -1.8986e-02,
            +5.6896e-03,
            -2.4159e-04,
            -2.3033e-04,
            -9.6783e-06,
            +0.0000e00,
            -1.0218e-01,
            -1.3916e-02,
            -4.1025e-03,
            -5.1340e-05,
            -7.0114e-05,
            -3.3152e-07,
            +1.6901e-06,
            +0.0000e00,
            -1.2422e-02,
            +2.5072e-03,
            +1.1205e-03,
            -1.3034e-04,
            -2.3971e-05,
            -2.6622e-06,
            +5.7852e-07,
            +4.5847e-08,
            +0.0000e00,
            +4.4777e-02,
            -3.0421e-03,
            +2.6062e-05,
            -7.2421e-05,
            +1.9119e-06,
            +3.9236e-07,
            +2.2390e-07,
            +2.9765e-09,
            -4.6452e-09,
        ]
    )

    # Temperature coefficients
    at_mean = np.array(
        [
            +1.6257e01,
            +2.1224e00,
            +9.2569e-01,
            -2.5974e01,
            +1.4510e00,
            +9.2468e-02,
            -5.3192e-01,
            +2.1094e-01,
            -6.9210e-02,
            -3.4060e-02,
            -4.6569e00,
            +2.6385e-01,
            -3.6093e-02,
            +1.0198e-02,
            -1.8783e-03,
            +7.4983e-01,
            +1.1741e-01,
            +3.9940e-02,
            +5.1348e-03,
            +5.9111e-03,
            +8.6133e-06,
            +6.3057e-01,
            +1.5203e-01,
            +3.9702e-02,
            +4.6334e-03,
            +2.4406e-04,
            +1.5189e-04,
            +1.9581e-07,
            +5.4414e-01,
            +3.5722e-01,
            +5.2763e-02,
            +4.1147e-03,
            -2.7239e-04,
            -5.9957e-05,
            +1.6394e-06,
            -7.3045e-07,
            -2.9394e00,
            +5.5579e-02,
            +1.8852e-02,
            +3.4272e-03,
            -2.3193e-05,
            -2.9349e-05,
            +3.6397e-07,
            +2.0490e-06,
            -6.4719e-08,
            -5.2225e-01,
            +2.0799e-01,
            +1.3477e-03,
            +3.1613e-04,
            -2.2285e-04,
            -1.8137e-05,
            -1.5177e-07,
            +6.1343e-07,
            +7.8566e-08,
            +1.0749e-09,
        ]
    )

    bt_mean = np.array(
        [
            +0.0000e00,
            +0.0000e00,
            +1.0210e00,
            +0.0000e00,
            +6.0194e-01,
            +1.2292e-01,
            +0.0000e00,
            -4.2184e-01,
            +1.8230e-01,
            +4.2329e-02,
            +0.0000e00,
            +9.3312e-02,
            +9.5346e-02,
            -1.9724e-03,
            +5.8776e-03,
            +0.0000e00,
            -2.0940e-01,
            +3.4199e-02,
            -5.7672e-03,
            -2.1590e-03,
            +5.6815e-04,
            +0.0000e00,
            +2.2858e-01,
            +1.2283e-02,
            -9.3679e-03,
            -1.4233e-03,
            -1.5962e-04,
            +4.0160e-05,
            +0.0000e00,
            +3.6353e-02,
            -9.4263e-04,
            -3.6762e-03,
            +5.8608e-05,
            -2.6391e-05,
            +3.2095e-06,
            -1.1605e-06,
            +0.0000e00,
            +1.6306e-01,
            +1.3293e-02,
            -1.1395e-03,
            +5.1097e-05,
            +3.3977e-05,
            +7.6449e-06,
            -1.7602e-07,
            -7.6558e-08,
            +0.0000e00,
            -4.5415e-02,
            -1.8027e-02,
            +3.6561e-04,
            -1.1274e-04,
            +1.3047e-05,
            +2.0001e-06,
            -1.5152e-07,
            -2.7807e-08,
            +7.7491e-09,
        ]
    )

    at_amp = np.array(
        [
            -1.8654e00,
            -9.0041e00,
            -1.2974e-01,
            -3.6053e00,
            +2.0284e-02,
            +2.1872e-01,
            -1.3015e00,
            +4.0355e-01,
            +2.2216e-01,
            -4.0605e-03,
            +1.9623e00,
            +4.2887e-01,
            +2.1437e-01,
            -1.0061e-02,
            -1.1368e-03,
            -6.9235e-02,
            +5.6758e-01,
            +1.1917e-01,
            -7.0765e-03,
            +3.0017e-04,
            +3.0601e-04,
            +1.6559e00,
            +2.0722e-01,
            +6.0013e-02,
            +1.7023e-04,
            -9.2424e-04,
            +1.1269e-05,
            -6.9911e-06,
            -2.0886e00,
            -6.7879e-02,
            -8.5922e-04,
            -1.6087e-03,
            -4.5549e-05,
            +3.3178e-05,
            -6.1715e-06,
            -1.4446e-06,
            -3.7210e-01,
            +1.5775e-01,
            -1.7827e-03,
            -4.4396e-04,
            +2.2844e-04,
            -1.1215e-05,
            -2.1120e-06,
            -9.6421e-07,
            -1.4170e-08,
            +7.8720e-01,
            -4.4238e-02,
            -1.5120e-03,
            -9.4119e-04,
            +4.0645e-06,
            -4.9253e-06,
            -1.8656e-06,
            -4.0736e-07,
            -4.9594e-08,
            +1.6134e-09,
        ]
    )

    bt_amp = np.array(
        [
            +0.0000e00,
            +0.0000e00,
            -8.9895e-01,
            +0.0000e00,
            -1.0790e00,
            -1.2699e-01,
            +0.0000e00,
            -5.9033e-01,
            +3.4865e-02,
            -3.2614e-02,
            +0.0000e00,
            -2.4310e-02,
            +1.5607e-02,
            -2.9833e-02,
            -5.9048e-03,
            +0.0000e00,
            +2.8383e-01,
            +4.0509e-02,
            -1.8834e-02,
            -1.2654e-03,
            -1.3794e-04,
            +0.0000e00,
            +1.3306e-01,
            +3.4960e-02,
            -3.6799e-03,
            -3.5626e-04,
            +1.4814e-04,
            +3.7932e-06,
            +0.0000e00,
            +2.0801e-01,
            +6.5640e-03,
            -3.4893e-03,
            -2.7395e-04,
            +7.4296e-05,
            -7.9927e-06,
            -1.0277e-06,
            +0.0000e00,
            +3.6515e-02,
            -7.4319e-03,
            -6.2873e-04,
            -8.2461e-05,
            +3.1095e-05,
            -5.3860e-07,
            -1.2055e-07,
            -1.1517e-07,
            +0.0000e00,
            +3.1404e-02,
            +1.5580e-02,
            -1.1428e-03,
            +3.3529e-05,
            +1.0387e-05,
            -1.9378e-06,
            -2.7327e-07,
            +7.5833e-09,
            -9.2323e-09,
        ]
    )

    # Parameter t
    t = math.sin(dlat)

    # Degree n and order m
    n = 9
    m = 9

    # Determine n! (factorial) moved by 1
    dfac = np.ones(2 * n + 2)
    for i in range(1, 2 * n + 2):
        dfac[i] = dfac[i - 1] * i

    # Determine Legendre functions (Heiskanen and Moritz, Physical Geodesy, 1967, eq. 1-62)
    P = np.zeros((n + 1, n + 1))
    for i in range(n + 1):
        for j in range(min(i, m) + 1):
            ir = int((i - j) / 2)
            sum_val = 0
            for k in range(ir + 1):
                sum_val = sum_val + (
                    (-1) ** k
                    * dfac[2 * i - 2 * k]
                    / dfac[k]
                    / dfac[i - k]
                    / dfac[i - j - 2 * k]
                    * t ** (i - j - 2 * k)
                )
            # Legendre functions moved by 1
            P[i, j] = 1 / 2**i * math.sqrt((1 - t**2) ** j) * sum_val

    # Spherical harmonics
    aP = np.zeros(55)
    bP = np.zeros(55)
    i = 0
    for n_idx in range(10):
        for m_idx in range(n_idx + 1):
            aP[i] = P[n_idx, m_idx] * math.cos(m_idx * dlon)
            bP[i] = P[n_idx, m_idx] * math.sin(m_idx * dlon)
            i = i + 1

    # Geoidal height
    undu = 0
    for i in range(55):
        undu = undu + (a_geoid[i] * aP[i] + b_geoid[i] * bP[i])

    # Orthometric height
    hort = dhgt - undu

    # Surface pressure on the geoid
    apm = 0
    apa = 0
    for i in range(55):
        apm = apm + (ap_mean[i] * aP[i] + bp_mean[i] * bP[i])
        apa = apa + (ap_amp[i] * aP[i] + bp_amp[i] * bP[i])
    pres0 = apm + apa * math.cos(doy / 365.25 * 2 * pi)

    # Height correction for pressure
    pres = pres0 * (1 - 0.0000226 * hort) ** 5.225

    # Surface temperature on the geoid
    atm = 0
    ata = 0
    for i in range(55):
        atm = atm + (at_mean[i] * aP[i] + bt_mean[i] * bP[i])
        ata = ata + (at_amp[i] * aP[i] + bt_amp[i] * bP[i])
    temp0 = atm + ata * math.cos(doy / 365.25 * 2 * pi)

    # Height correction for temperature
    temp = temp0 - 0.0065 * hort

    return pres, temp, undu


def get_vmf1g_values(data_dir: str, year: int, month: int, day: int, sec: float,
                     lat: float, lon: float) -> Dict[str, Union[float, np.ndarray]]:
    """
    Get VMF1G values for a given location and time.

    Parameters:
    -----------
    data_dir : str
        Directory containing VMF1G data files
    year, month, day : int
        Calendar date
    sec : float
        Seconds since midnight
    lat, lon : float
        Latitude and longitude in degrees

    Returns:
    --------
    dict : Dictionary containing VMF1G grid data and interpolation parameters
    """
    # Initialize output arrays
    ahg = np.zeros((4, 2))
    awg = np.zeros((4, 2))
    zhdg = np.zeros((4, 2))
    zwdg = np.zeros((4, 2))

    fic0_tmp = 'VMFG_'

    # Determine time windows
    dh0 = sec / 3600.0 - 0.0
    dh6 = sec / 3600.0 - 6.0
    dh12 = sec / 3600.0 - 12.0
    dh18 = sec / 3600.0 - 18.0

    if dh0 >= 0.0 and dh0 < 6.0:
        dh = dh0
        h1 = 0.0
        h2 = 6.0
    elif dh6 >= 0.0 and dh6 < 6.0:
        dh = dh6
        h1 = 6.0
        h2 = 12.0
    elif dh12 >= 0.0 and dh12 < 6.0:
        dh = dh12
        h1 = 12.0
        h2 = 18.0
    elif dh18 >= 0.0 and dh18 < 6.0:
        dh = dh18
        h1 = 18.0
        h2 = 24.0
    else:
        # Default case
        dh = dh0
        h1 = 0.0
        h2 = 6.0

    # Format date/time strings
    year_str = f"{int(year):04d}"
    month_str = f"{int(month):02d}"
    day_str = f"{int(day):02d}"
    h1_str = f"{int(h1):02d}"

    fic1 = os.path.join(data_dir, f"{fic0_tmp}{year_str}{month_str}{day_str}.H{h1_str}")

    # Handle day rollover for h2=24
    jy = year // 4
    ix = year - jy * 4

    if h2 == 24.0:
        day_tmp = day + 1
        month_tmp = month

        # Handle month rollover
        if month in [1, 3, 5, 7, 8, 10, 12] and day_tmp == 32:
            day_tmp = 1
            month_tmp = month + 1
        elif month == 2 and ix == 0 and day_tmp == 30:  # Leap year
            day_tmp = 1
            month_tmp = month + 1
        elif month == 2 and ix != 0 and day_tmp == 29:  # Non-leap year
            day_tmp = 1
            month_tmp = month + 1

        day_str = f"{int(day_tmp):02d}"
        month_str = f"{int(month_tmp):02d}"
        h2_str = "00"
    else:
        day_str = f"{int(day):02d}"
        h2_str = f"{int(h2):02d}"

    fic2 = os.path.join(data_dir, f"{fic0_tmp}{year_str}{month_str}{day_str}.H{h2_str}")

    print('----------------------------------------------------------')
    print(f'fic1 : {fic1}')
    print(f'fic2 : {fic2}')

    # Calculate grid coordinates
    latg4 = lat - (lat % 2.0)
    if lat < 0:
        latg4 = latg4 - 2.0
    latg1 = latg4 + 2.0

    long4 = lon - (lon % 2.5)
    long3 = long4 + 2.5

    if latg1 > 90.0:  # North pole overflow
        latg1 = 90.0

    if long3 > 357.5:  # Greenwich meridian overflow
        long3 = 0.0

    latg2 = latg1
    latg3 = latg4
    long1 = long4
    long2 = long3

    def read_vmf_file(filename: str, time_index: int):
        """Read VMF data from file and extract grid point values"""
        if not os.path.exists(filename):
            print(f"Warning: File {filename} not found. Using default values.")
            return

        with open(filename, 'r') as f:
            # Skip header (7 lines)
            for _ in range(7):
                f.readline()

            # Read data lines
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 6:
                    lat_tmp = float(parts[0])
                    lon_tmp = float(parts[1])
                    ah_tmp = float(parts[2])
                    aw_tmp = float(parts[3])
                    zhd_tmp = float(parts[4])
                    zwd_tmp = float(parts[5])

                    # Check which grid point this corresponds to
                    if abs(lat_tmp - latg1) < 0.1 and abs(lon_tmp - long1) < 0.1:  # Node #1
                        ahg[0, time_index] = ah_tmp
                        awg[0, time_index] = aw_tmp
                        zhdg[0, time_index] = zhd_tmp
                        zwdg[0, time_index] = zwd_tmp
                    elif abs(lat_tmp - latg2) < 0.1 and abs(lon_tmp - long2) < 0.1:  # Node #2
                        ahg[1, time_index] = ah_tmp
                        awg[1, time_index] = aw_tmp
                        zhdg[1, time_index] = zhd_tmp
                        zwdg[1, time_index] = zwd_tmp
                    elif abs(lat_tmp - latg4) < 0.1 and abs(lon_tmp - long4) < 0.1:  # Node #4
                        ahg[3, time_index] = ah_tmp
                        awg[3, time_index] = aw_tmp
                        zhdg[3, time_index] = zhd_tmp
                        zwdg[3, time_index] = zwd_tmp
                    elif abs(lat_tmp - latg3) < 0.1 and abs(lon_tmp - long3) < 0.1:  # Node #3
                        ahg[2, time_index] = ah_tmp
                        awg[2, time_index] = aw_tmp
                        zhdg[2, time_index] = zhd_tmp
                        zwdg[2, time_index] = zwd_tmp

        # Handle special case for polar limits (North pole only)
        if ahg[2, time_index] == 0 and ahg[3, time_index] == 0:
            # Copy values from nodes 1,2 to nodes 4,3 respectively
            ahg[2, time_index] = ahg[1, time_index]
            ahg[3, time_index] = ahg[0, time_index]
            awg[2, time_index] = awg[1, time_index]
            awg[3, time_index] = awg[0, time_index]
            zhdg[2, time_index] = zhdg[1, time_index]
            zhdg[3, time_index] = zhdg[0, time_index]
            zwdg[2, time_index] = zwdg[1, time_index]
            zwdg[3, time_index] = zwdg[0, time_index]

    # Read both files
    read_vmf_file(fic1, 0)
    read_vmf_file(fic2, 1)

    return {
        'dh': dh,
        'latg1': latg1, 'latg2': latg2, 'latg3': latg3, 'latg4': latg4,
        'long1': long1, 'long2': long2, 'long3': long3, 'long4': long4,
        'ahg': ahg,
        'awg': awg,
        'zhdg': zhdg,
        'zwdg': zwdg
    }


def compute_vmf1g_values(data_dir: str, year: int, month: int, day: int,
                        hh: int, min_val: int, sec: int, lat: float, lon: float,
                        h: float, el: float) -> Dict[str, float]:
    """
    Main function to compute VMF1G atmospheric delay values.

    Parameters:
    -----------
    data_dir : str
        Directory containing VMF1G data files
    year, month, day : int
        Calendar date
    hh, min_val, sec : int
        Time (hour, minute, second)
    lat, lon : float
        Latitude and longitude in degrees
    h : float
        Ellipsoidal height in meters
    el : float
        Elevation angle in degrees

    Returns:
    --------
    dict : Dictionary containing computed atmospheric delay values
    """
    pi = math.pi

    # Convert to radians
    phi = lat * pi / 180.0
    lambda_rad = lon * pi / 180.0

    # Handle day rollover
    jy = year // 4
    ix = year - jy * 4
    if hh == 24:
        hh = 0
        day = day + 1
        if month in [1, 3, 5, 7, 8, 10, 12] and day == 32:
            day = 1
            month = month + 1
        elif month == 2 and ix == 0 and day == 30:  # Leap year
            day = 1
            month = month + 1
        elif month == 2 and ix != 0 and day == 29:  # Non-leap year
            day = 1
            month = month + 1

    sec_long = float(hh * 3600 + min_val * 60 + sec)

    # Convert calendar date to Julian days (CNES & regular
    jj_sec_cnes = conv.dt2jjul_cnes(dt.datetime(year, month, day, hh, min_val, sec))
    jj_sec = conv.dt2mjd(dt.datetime(year, month, day, hh, min_val, sec))


    # Get VMF1G values
    vmf_data = get_vmf1g_values(data_dir, year, month, day, sec_long, lat, lon)

    # Get orographic heights
    oro_zhdg = get_orography(data_dir, lat, lon)

    # Temporal interpolation
    dh = vmf_data['dh']
    ahg = vmf_data['ahg']
    awg = vmf_data['awg']
    zhdg = vmf_data['zhdg']
    zwdg = vmf_data['zwdg']

    # Interpolate in time
    ahg_date = np.zeros(4)
    awg_date = np.zeros(4)
    zhdg_date = np.zeros(4)
    zwdg_date = np.zeros(4)

    for i in range(4):
        ahg_date[i] = ahg[i, 0] + dh * (ahg[i, 1] - ahg[i, 0]) / 6.0
        awg_date[i] = awg[i, 0] + dh * (awg[i, 1] - awg[i, 0]) / 6.0
        zhdg_date[i] = zhdg[i, 0] + dh * (zhdg[i, 1] - zhdg[i, 0]) / 6.0
        zwdg_date[i] = zwdg[i, 0] + dh * (zwdg[i, 1] - zwdg[i, 0]) / 6.0

    # Geographic interpolation for ah and aw
    latg1, latg2, latg3, latg4 = vmf_data['latg1'], vmf_data['latg2'], vmf_data['latg3'], vmf_data['latg4']
    long1, long2, long3, long4 = vmf_data['long1'], vmf_data['long2'], vmf_data['long3'], vmf_data['long4']

    # Interpolate along parallels and meridian for ah
    ahg_tmp1 = ahg_date[0] + (lon - long1) * (ahg_date[1] - ahg_date[0]) / 2.5
    ahg_tmp2 = ahg_date[3] + (lon - long4) * (ahg_date[2] - ahg_date[3]) / 2.5
    ahg_interp = ahg_tmp2 + (lat - latg4) * (ahg_tmp1 - ahg_tmp2) / 2.0

    # Interpolate along parallels and meridian for aw
    awg_tmp1 = awg_date[0] + (lon - long1) * (awg_date[1] - awg_date[0]) / 2.5
    awg_tmp2 = awg_date[3] + (lon - long4) * (awg_date[2] - awg_date[3]) / 2.5
    awg_interp = awg_tmp2 + (lat - latg4) * (awg_tmp1 - awg_tmp2) / 2.0

    # Compute mapping functions
    vmf1h, vmf1w = sd_vmf1_grid_elevation(ahg_interp, awg_interp, jj_sec, phi, h,
                                          (pi/2 - el*pi/180.0))

    # GMF mapping functions for comparison
    sd_gmf_use = sd_gmf
    gmfh, gmfw = sd_gmf_use(jj_sec, phi, lambda_rad, h, (pi / 2 - el * pi / 180.0))

    # GPT pressure and temperature
    sd_gpt_use = sd_gpt
    pres0, temp0, undu0 = sd_gpt_use(jj_sec, phi, lambda_rad, h)

    # Calculate pressures at grid nodes
    pressures = []
    temps = []
    for i in range(4):
        lat_grid = [latg1, latg2, latg3, latg4][i] * pi / 180.0
        lon_grid = [long1, long2, long3, long4][i] * pi / 180.0
        pres, temp, undu = sd_gpt_use(jj_sec, lat_grid, lon_grid, oro_zhdg[i])
        pressures.append(pres)
        temps.append(temp)

    # Calculate ZHD using Saastamoinen formula
    zhd0 = 2.2768e-3 * pres0 / (1. - 0.00266*math.cos(2*phi) - 0.00028*h*1.e-3)

    zhds = []
    for i in range(4):
        lat_grid = [latg1, latg2, latg3, latg4][i] * pi / 180.0
        zhd = 2.2768e-3 * pressures[i] / (1. - 0.00266*math.cos(2*lat_grid) - 0.00028*oro_zhdg[i]*1.e-3)
        zhds.append(zhd)

    # Correct ZHDG values
    zhdg_date_init = zhdg_date.copy()
    for i in range(4):
        zhdg_date[i] = zhdg_date[i] + (zhd0 - zhds[i])

    # Geographic interpolation for corrected ZHDG and ZWDG
    zhdg_tmp1 = zhdg_date[0] + (lon - long1) * (zhdg_date[1] - zhdg_date[0]) / 2.5
    zhdg_tmp2 = zhdg_date[3] + (lon - long4) * (zhdg_date[2] - zhdg_date[3]) / 2.5
    zhdg_interp = zhdg_tmp2 + (lat - latg4) * (zhdg_tmp1 - zhdg_tmp2) / 2.0

    zwdg_tmp1 = zwdg_date[0] + (lon - long1) * (zwdg_date[1] - zwdg_date[0]) / 2.5
    zwdg_tmp2 = zwdg_date[3] + (lon - long4) * (zwdg_date[2] - zwdg_date[3]) / 2.5
    zwdg_interp = zwdg_tmp2 + (lat - latg4) * (zwdg_tmp1 - zwdg_tmp2) / 2.0

    # Print results
    print('----------------------------------------------------------')
    print(f'JJ CNES   : {jj_sec:14.8f}')
    print(f'date cal  : {day:02d} {month:02d} {year:4d}')
    print(f'sec       : {sec:02d}')
    print(f'el        : {el:7.4f}')
    print(f'dt/grd    : {dh:7.1f}')
    print(f'position  : {lat:7.1f} {lon:7.1f} {h:7.1f}')
    print(f'lat_nodes : {latg1:5.1f} {latg2:5.1f} {latg3:5.1f} {latg4:5.1f}')
    print(f'lon_nodes : {long1:5.1f} {long2:5.1f} {long3:5.1f} {long4:5.1f}')
    print(f'oro.      : {oro_zhdg[0]:7.1f} {oro_zhdg[1]:7.1f} {oro_zhdg[2]:7.1f} {oro_zhdg[3]:7.1f}')
    print(f'ahg_date  : {ahg_date[0]:10.8f} {ahg_date[1]:10.8f} {ahg_date[2]:10.8f} {ahg_date[3]:10.8f}')
    print(f'awg_date  : {awg_date[0]:10.8f} {awg_date[1]:10.8f} {awg_date[2]:10.8f} {awg_date[3]:10.8f}')
    print(f'zhdg_date0: {zhdg_date_init[0]:7.4f} {zhdg_date_init[1]:7.4f} {zhdg_date_init[2]:7.4f} {zhdg_date_init[3]:7.4f}')
    print(f'zhdg_redu : {zhd0-zhds[0]:7.4f} {zhd0-zhds[1]:7.4f} {zhd0-zhds[2]:7.4f} {zhd0-zhds[3]:7.4f}')
    print(f'zhdg_date : {zhdg_date[0]:7.4f} {zhdg_date[1]:7.4f} {zhdg_date[2]:7.4f} {zhdg_date[3]:7.4f}')
    print(f'zwdg_date : {zwdg_date[0]:7.4f} {zwdg_date[1]:7.4f} {zwdg_date[2]:7.4f} {zwdg_date[3]:7.4f}')
    print('----------------------------------------------------------')
    print(f'GMF/mf    : {gmfh:14.8f} {gmfw:14.8f}')
    print(f'GPT/ZD    : {zhd0:14.8f}')
    print(f'GPT*GMF/SD: {zhd0*gmfh:14.8f}')
    print(f'VMF1G/mf  : {vmf1h:14.8f} {vmf1w:14.8f}')
    print(f'VMF1G/ZD  : {zhdg_interp:14.8f} {zwdg_interp:14.8f}')
    print(f'VMF1G/SD  : {zhdg_interp*vmf1h:14.8f} {zwdg_interp*vmf1w:14.8f}')
    print('----------------------------------------------------------')

    return {
        'jj_sec': jj_sec,
        'vmf1h': vmf1h,
        'vmf1w': vmf1w,
        'gmfh': gmfh,
        'gmfw': gmfw,
        'zhd_gpt': zhd0,
        'zhdg_interp': zhdg_interp,
        'zwdg_interp': zwdg_interp,
        'slant_delay_h': zhdg_interp * vmf1h,
        'slant_delay_w': zwdg_interp * vmf1w,
        'slant_delay_total': zhdg_interp * vmf1h + zwdg_interp * vmf1w
    }

def main():
    """Main function to run VMF1G computation from command line or interactive input."""
    import argparse

    parser = argparse.ArgumentParser(description='Compute VMF1G atmospheric delay values')
    parser.add_argument('data_dir', nargs='?', help='Data directory')
    parser.add_argument('year', type=int, nargs='?', help='Year')
    parser.add_argument('month', type=int, nargs='?', help='Month')
    parser.add_argument('day', type=int, nargs='?', help='Day')
    parser.add_argument('hour', type=int, nargs='?', help='Hour')
    parser.add_argument('minute', type=int, nargs='?', help='Minute')
    parser.add_argument('second', type=int, nargs='?', help='Second')
    parser.add_argument('latitude', type=float, nargs='?', help='Latitude (degrees)')
    parser.add_argument('longitude', type=float, nargs='?', help='Longitude (degrees, 0-360)')
    parser.add_argument('height', type=float, nargs='?', help='Height above ellipsoid (m)')
    parser.add_argument('elevation', type=float, nargs='?', help='Elevation angle (degrees)')

    args = parser.parse_args()

    # Get parameters from command line or interactive input
    if args.data_dir is None:
        data_dir = input("Data directory ? : ")
        year = int(input("Year ? : "))
        month = int(input("Month ? : "))
        day = int(input("Day ? : "))
        hour = int(input("Hour(s) ? : "))
        minute = int(input("Minute(s) ? : "))
        second = int(input("Second(s) ? : "))
        lat = float(input("Latitude(deg) ? : "))
        lon = float(input("Longitude(deg;[0-360]) ? : "))
        h = float(input("Height above ell.(m) ? : "))
        el = float(input("Elevation angle(deg) ? : "))
    else:
        data_dir = args.data_dir
        year = args.year
        month = args.month
        day = args.day
        hour = args.hour
        minute = args.minute
        second = args.second
        lat = args.latitude
        lon = args.longitude
        h = args.height
        el = args.elevation

    # Compute VMF1G values
    results = compute_vmf1g_values(data_dir, year, month, day, hour, minute, second,
                                  lat, lon, h, el)

    return results


if __name__ == "__main__":
    main()
