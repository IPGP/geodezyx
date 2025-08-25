# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd


def gmf(dmjd=None, dlat=None, dlon=None, dhgt=None, zd=None):
    """
    Determine the Global Mapping Functions (GMF) for hydrostatic and wet components.

    This function calculates mapping functions used in atmospheric delay modeling
    for satellite geodesy applications. The GMF is based on numerical weather model
    data and provides improved accuracy over earlier mapping functions.

    Parameters
    ----------
    dmjd : float
        Modified Julian Date.
    dlat : float
        Ellipsoidal latitude in radians.
    dlon : float
        Longitude in radians.
    dhgt : float
        Height above sea level in meters.
    zd : float
        Zenith distance in radians.

    Returns
    -------
    gmfh : float
        Hydrostatic mapping function value.
    gmfw : float
        Wet mapping function value.

    References
    ----------
    Boehm, J., A.E. Niell, p. Tregoning, H. Schuh (2006),
    Global Mapping Functions (GMF): A new empirical mapping function based on
    numerical weather model data, Geoph. Res. Letters, Vol. 33, L07304,
    doi:10.1029/2005GL025545.

    Notes
    -----
    Original implementation by Johannes Boehm, 2005 August 30.
    Recursions for Legendre polynomials updated 2006 Aug. 14 (O. Montenbruck).
    Latitude changed to ellipsoidal latitude 2011 Jul. 21 (J. Boehm).
    Converted to Python by Zohreh Adavi (zohreh.adavi@tuwien.ac.at), 2023-02-20.
    """
    ah_mean = np.array(
        [
            +125.17,
            +0.8503,
            +0.06936,
            -6.76,
            +0.1771,
            +0.0113,
            +0.5963,
            +0.01808,
            +0.002801,
            -0.001414,
            -1.212,
            +0.093,
            +0.003683,
            +0.001095,
            +4.671e-05,
            +0.3959,
            -0.03867,
            +0.005413,
            -0.0005289,
            +0.0003229,
            +2.067e-05,
            +0.3,
            +0.02031,
            +0.0059,
            +0.0004573,
            -7.619e-05,
            +2.327e-06,
            +3.845e-06,
            +0.1182,
            +0.01158,
            +0.005445,
            +6.219e-05,
            +4.204e-06,
            -2.093e-06,
            +1.54e-07,
            -4.28e-08,
            -0.4751,
            -0.0349,
            +0.001758,
            +0.0004019,
            -2.799e-06,
            -1.287e-06,
            +5.468e-07,
            +7.58e-08,
            -6.3e-09,
            -0.116,
            +0.008301,
            +0.0008771,
            +9.955e-05,
            -1.718e-06,
            -2.012e-06,
            +1.17e-08,
            +1.79e-08,
            -1.3e-09,
            +1e-10,
        ]
    )
    bh_mean = np.array(
        [
            +0.0,
            +0.0,
            +0.03249,
            +0.0,
            +0.03324,
            +0.0185,
            +0.0,
            -0.1115,
            +0.02519,
            +0.004923,
            +0.0,
            +0.02737,
            +0.01595,
            -0.0007332,
            +0.0001933,
            +0.0,
            -0.04796,
            +0.006381,
            -0.0001599,
            -0.0003685,
            +1.815e-05,
            +0.0,
            +0.07033,
            +0.002426,
            -0.001111,
            -0.0001357,
            -7.828e-06,
            +2.547e-06,
            +0.0,
            +0.005779,
            +0.003133,
            -0.0005312,
            -2.028e-05,
            +2.323e-07,
            -9.1e-08,
            -1.65e-08,
            +0.0,
            +0.03688,
            -0.0008638,
            -8.514e-05,
            -2.828e-05,
            +5.403e-07,
            +4.39e-07,
            +1.35e-08,
            +1.8e-09,
            +0.0,
            -0.02736,
            -0.0002977,
            +8.113e-05,
            +2.329e-07,
            +8.451e-07,
            +4.49e-08,
            -8.1e-09,
            -1.5e-09,
            +2e-10,
        ]
    )
    ah_amp = np.array(
        [
            -0.2738,
            -2.837,
            +0.01298,
            -0.3588,
            +0.02413,
            +0.03427,
            -0.7624,
            +0.07272,
            +0.0216,
            -0.003385,
            +0.4424,
            +0.03722,
            +0.02195,
            -0.001503,
            +0.0002426,
            +0.3013,
            +0.05762,
            +0.01019,
            -0.0004476,
            +6.79e-05,
            +3.227e-05,
            +0.3123,
            -0.03535,
            +0.00484,
            +3.025e-06,
            -4.363e-05,
            +2.854e-07,
            -1.286e-06,
            -0.6725,
            -0.0373,
            +0.0008964,
            +0.0001399,
            -3.99e-06,
            +7.431e-06,
            -2.796e-07,
            -1.601e-07,
            +0.04068,
            -0.01352,
            +0.0007282,
            +9.594e-05,
            +2.07e-06,
            -9.62e-08,
            -2.742e-07,
            -6.37e-08,
            -6.3e-09,
            +0.08625,
            -0.005971,
            +0.0004705,
            +2.335e-05,
            +4.226e-06,
            +2.475e-07,
            -8.85e-08,
            -3.6e-08,
            -2.9e-09,
            +0.0,
        ]
    )
    bh_amp = np.array(
        [
            +0.0,
            +0.0,
            -0.1136,
            +0.0,
            -0.1868,
            -0.01399,
            +0.0,
            -0.1043,
            +0.01175,
            -0.00224,
            +0.0,
            -0.03222,
            +0.01333,
            -0.002647,
            -2.316e-05,
            +0.0,
            +0.05339,
            +0.01107,
            -0.003116,
            -0.0001079,
            -1.299e-05,
            +0.0,
            +0.004861,
            +0.008891,
            -0.0006448,
            -1.279e-05,
            +6.358e-06,
            -1.417e-07,
            +0.0,
            +0.03041,
            +0.00115,
            -0.0008743,
            -2.781e-05,
            +6.367e-07,
            -1.14e-08,
            -4.2e-08,
            +0.0,
            -0.02982,
            -0.003,
            +1.394e-05,
            -3.29e-05,
            -1.705e-07,
            +7.44e-08,
            +2.72e-08,
            -6.6e-09,
            +0.0,
            +0.01236,
            -0.0009981,
            -3.792e-05,
            -1.355e-05,
            +1.162e-06,
            -1.789e-07,
            +1.47e-08,
            -2.4e-09,
            -4e-10,
        ]
    )
    aw_mean = np.array(
        [
            +56.4,
            +1.555,
            -1.011,
            -3.975,
            +0.03171,
            +0.1065,
            +0.6175,
            +0.1376,
            +0.04229,
            +0.003028,
            +1.688,
            -0.1692,
            +0.05478,
            +0.02473,
            +0.0006059,
            +2.278,
            +0.006614,
            -0.0003505,
            -0.006697,
            +0.0008402,
            +0.0007033,
            -3.236,
            +0.2184,
            -0.04611,
            -0.01613,
            -0.001604,
            +5.42e-05,
            +7.922e-05,
            -0.2711,
            -0.4406,
            -0.03376,
            -0.002801,
            -0.000409,
            -2.056e-05,
            +6.894e-06,
            +2.317e-06,
            +1.941,
            -0.2562,
            +0.01598,
            +0.005449,
            +0.0003544,
            +1.148e-05,
            +7.503e-06,
            -5.667e-07,
            -3.66e-08,
            +0.8683,
            -0.05931,
            -0.001864,
            -0.0001277,
            +0.0002029,
            +1.269e-05,
            +1.629e-06,
            +9.66e-08,
            -1.015e-07,
            -5e-10,
        ]
    )
    bw_mean = np.array(
        [
            +0.0,
            +0.0,
            +0.2592,
            +0.0,
            +0.02974,
            -0.5471,
            +0.0,
            -0.5926,
            -0.103,
            -0.01567,
            +0.0,
            +0.171,
            +0.09025,
            +0.02689,
            +0.002243,
            +0.0,
            +0.3439,
            +0.02402,
            +0.00541,
            +0.001601,
            +9.669e-05,
            +0.0,
            +0.09502,
            -0.03063,
            -0.001055,
            -0.0001067,
            -0.000113,
            +2.124e-05,
            +0.0,
            -0.3129,
            +0.008463,
            +0.0002253,
            +7.413e-05,
            -9.376e-05,
            -1.606e-06,
            +2.06e-06,
            +0.0,
            +0.2739,
            +0.001167,
            -2.246e-05,
            -0.0001287,
            -2.438e-05,
            -7.561e-07,
            +1.158e-06,
            +4.95e-08,
            +0.0,
            -0.1344,
            +0.005342,
            +0.0003775,
            -6.756e-05,
            -1.686e-06,
            -1.184e-06,
            +2.768e-07,
            +2.73e-08,
            +5.7e-09,
        ]
    )
    aw_amp = np.array(
        [
            +0.1023,
            -2.695,
            +0.3417,
            -0.1405,
            +0.3175,
            +0.2116,
            +3.536,
            -0.1505,
            -0.0166,
            +0.02967,
            +0.3819,
            -0.1695,
            -0.07444,
            +0.007409,
            -0.006262,
            -1.836,
            -0.01759,
            -0.06256,
            -0.002371,
            +0.0007947,
            +0.0001501,
            -0.8603,
            -0.136,
            -0.03629,
            -0.003706,
            -0.0002976,
            +1.857e-05,
            +3.021e-05,
            +2.248,
            -0.1178,
            +0.01255,
            +0.001134,
            -0.0002161,
            -5.817e-06,
            +8.836e-07,
            -1.769e-07,
            +0.7313,
            -0.1188,
            +0.01145,
            +0.001011,
            +0.0001083,
            +2.57e-06,
            -2.14e-06,
            -5.71e-08,
            +2e-08,
            -1.632,
            -0.006948,
            -0.003893,
            +0.0008592,
            +7.577e-05,
            +4.539e-06,
            -3.852e-07,
            -2.213e-07,
            -1.37e-08,
            +5.8e-09,
        ]
    )
    bw_amp = np.array(
        [
            +0.0,
            +0.0,
            -0.08865,
            +0.0,
            -0.4309,
            +0.0634,
            +0.0,
            +0.1162,
            +0.06176,
            -0.004234,
            +0.0,
            +0.253,
            +0.04017,
            -0.006204,
            +0.004977,
            +0.0,
            -0.1737,
            -0.005638,
            +0.0001488,
            +0.0004857,
            -0.0001809,
            +0.0,
            -0.1514,
            -0.01685,
            +0.005333,
            -7.611e-05,
            +2.394e-05,
            +8.195e-06,
            +0.0,
            +0.09326,
            -0.01275,
            -0.0003071,
            +5.374e-05,
            -3.391e-05,
            -7.436e-06,
            +6.747e-07,
            +0.0,
            -0.08637,
            -0.003807,
            -0.0006833,
            -3.861e-05,
            -2.268e-05,
            +1.454e-06,
            +3.86e-07,
            -1.068e-07,
            +0.0,
            -0.02658,
            -0.001947,
            +0.0007131,
            -3.506e-05,
            +1.885e-07,
            +5.792e-07,
            +3.99e-08,
            +2e-08,
            -5.7e-09,
        ]
    )

    # reference day is 28 January
    # this is taken from Niell (1996) to be consistent
    doy = dmjd - 44239 + 1 - 28

    # degree n and order m
    nmax = 9
    mmax = 9
    # unit vector
    x = np.cos(dlat) * np.cos(dlon)
    y = np.cos(dlat) * np.sin(dlon)
    z = np.sin(dlat)
    # Legendre polynomials
    v = pd.DataFrame(np.zeros((nmax + 1, nmax + 1)))
    w = pd.DataFrame(np.zeros((nmax + 1, nmax + 1)))

    v.loc[0, 0] = 1
    w.loc[0, 0] = 0
    v.loc[1, 0] = z * v.loc[0, 0]
    w.loc[1, 0] = 0
    for n in range(2, nmax + 1):
        v.loc[n, 0] = (
            (2 * n - 1) * z * v.loc[n - 1, 0] - (n - 1) * v.loc[n - 2, 0]
        ) / n
        w.loc[n, 0] = 0

    for m in range(1, nmax + 1):
        v.loc[m, m] = (2 * m - 1) * (x * v.loc[m - 1, m - 1] - y * w.loc[m - 1, m - 1])
        w.loc[m, m] = (2 * m - 1) * (x * w.loc[m - 1, m - 1] + y * v.loc[m - 1, m - 1])
        if m < nmax:
            v.loc[m + 1, m] = (2 * m + 1) * z * v.loc[m, m]
            w.loc[m + 1, m] = (2 * m + 1) * z * w.loc[m, m]
        for n in range(m + 2, nmax + 1):
            v.loc[n, m] = (
                (2 * n - 1) * z * v.loc[n - 1, m] - (n + m - 1) * v.loc[n - 2, m]
            ) / (n - m)
            w.loc[n, m] = (
                (2 * n - 1) * z * w.loc[n - 1, m] - (n + m - 1) * w.loc[n - 2, m]
            ) / (n - m)

    # (1) hydrostatic mf

    bh = 0.0029
    c0h = 0.062
    if dlat < 0:
        phh = np.pi
        c11h = 0.007
        c10h = 0.002
    else:
        phh = 0
        c11h = 0.005
        c10h = 0.001

    ch = c0h + ((np.cos(doy / 365.25 * 2 * np.pi + phh) + 1) * c11h / 2 + c10h) * (
        1 - np.cos(dlat)
    )
    ahm = 0
    aha = 0
    i = 0
    for n in range(0, nmax + 1):
        for m in range(0, n + 1):

            ahm = ahm + (ah_mean[i] * v.loc[n, m] + bh_mean[i] * w.loc[n, m])
            aha = aha + (ah_amp[i] * v.loc[n, m] + bh_amp[i] * w.loc[n, m])
            i = i + 1

    ah = (ahm + aha * np.cos(doy / 365.25 * 2 * np.pi)) * 1e-05
    sine = np.sin(np.pi / 2 - zd)
    cose = np.cos(np.pi / 2 - zd)
    beta = bh / (sine + ch)
    gamma = ah / (sine + beta)
    topcon = 1 + ah / (1 + bh / (1 + ch))
    gmfh = topcon / (sine + gamma)
    # height correction for hydrostatic mapping function from Niell (1996) in order to reduce the coefficients to sea level
    a_ht = 2.53e-05
    b_ht = 0.00549
    c_ht = 0.00114
    hs_km = dhgt / 1000
    beta = b_ht / (sine + c_ht)
    gamma = a_ht / (sine + beta)
    topcon = 1 + a_ht / (1 + b_ht / (1 + c_ht))
    ht_corr_coef = 1 / sine - topcon / (sine + gamma)
    ht_corr = ht_corr_coef * hs_km
    gmfh = gmfh + ht_corr
    # (2) wet mf

    bw = 0.00146
    cw = 0.04391
    awm = 0
    awa = 0
    i = 0
    for n in range(0, nmax + 1):
        for m in range(0, n + 1):

            awm = awm + (aw_mean[i] * v.loc[n, m] + bw_mean[i] * w.loc[n, m])
            awa = awa + (aw_amp[i] * v.loc[n, m] + bw_amp[i] * w.loc[n, m])
            i = i + 1

    aw = (awm + awa * np.cos(doy / 365.25 * 2 * np.pi)) * 1e-05
    beta = bw / (sine + cw)
    gamma = aw / (sine + beta)
    topcon = 1 + aw / (1 + bw / (1 + cw))
    gmfw = topcon / (sine + gamma)
    return gmfh, gmfw
