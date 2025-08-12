#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 07:43:47 2024

@author: snahmani
"""


# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import pandas as pd   





def gmf(dmjd = None,dlat = None,dlon = None,dhgt = None,zd = None):
# =============================================================================  
#     This subroutine determines the Global Mapping Functions GMF
#     Reference: Boehm, J., A.E. Niell, P. Tregoning, H. Schuh (2006),
#     Global Mapping Functions (GMF): A new empirical mapping function based on numerical weather model data,
#     Geoph. Res. Letters, Vol. 33, L07304, doi:10.1029/2005GL025545.
    
#    input data
#    ----------
#    dmjd: modified julian date
#    dlat: ellipsoidal latitude in radians
#    dlon: longitude in radians
#    dhgt: height in m
#    zd:   zenith distance in radians
    
#    output data
#     -----------
#    gmfh: hydrostatic mapping function
#    gmfw: wet mapping function
    
#    Johannes Boehm, 2005 August 30
    
#    ref 2006 Aug. 14: recursions for Legendre polynomials (O. Montenbruck)
#    ref 2011 Jul. 21: latitude -> ellipsoidal latitude (J. Boehm)

#    File Converted to python script by Zohreh Adavi zohreh.adavi@tuwien.ac.at, 2023-02-20
# =============================================================================  
    
    ah_mean = np.array([+ 125.17,+ 0.8503,+ 0.06936,- 6.76,+ 0.1771,+ 0.0113,+ 0.5963,+ 0.01808,+ 0.002801,- 0.001414,- 1.212,+ 0.093,+ 0.003683,+ 0.001095,+ 4.671e-05,+ 0.3959,- 0.03867,+ 0.005413,- 0.0005289,+ 0.0003229,+ 2.067e-05,+ 0.3,+ 0.02031,+ 0.0059,+ 0.0004573,- 7.619e-05,+ 2.327e-06,+ 3.845e-06,+ 0.1182,+ 0.01158,+ 0.005445,+ 6.219e-05,+ 4.204e-06,- 2.093e-06,+ 1.54e-07,- 4.28e-08,- 0.4751,- 0.0349,+ 0.001758,+ 0.0004019,- 2.799e-06,- 1.287e-06,+ 5.468e-07,+ 7.58e-08,- 6.3e-09,- 0.116,+ 0.008301,+ 0.0008771,+ 9.955e-05,- 1.718e-06,- 2.012e-06,+ 1.17e-08,+ 1.79e-08,- 1.3e-09,+ 1e-10])
    bh_mean = np.array([+ 0.0,+ 0.0,+ 0.03249,+ 0.0,+ 0.03324,+ 0.0185,+ 0.0,- 0.1115,+ 0.02519,+ 0.004923,+ 0.0,+ 0.02737,+ 0.01595,- 0.0007332,+ 0.0001933,+ 0.0,- 0.04796,+ 0.006381,- 0.0001599,- 0.0003685,+ 1.815e-05,+ 0.0,+ 0.07033,+ 0.002426,- 0.001111,- 0.0001357,- 7.828e-06,+ 2.547e-06,+ 0.0,+ 0.005779,+ 0.003133,- 0.0005312,- 2.028e-05,+ 2.323e-07,- 9.1e-08,- 1.65e-08,+ 0.0,+ 0.03688,- 0.0008638,- 8.514e-05,- 2.828e-05,+ 5.403e-07,+ 4.39e-07,+ 1.35e-08,+ 1.8e-09,+ 0.0,- 0.02736,- 0.0002977,+ 8.113e-05,+ 2.329e-07,+ 8.451e-07,+ 4.49e-08,- 8.1e-09,- 1.5e-09,+ 2e-10])
    ah_amp = np.array([- 0.2738,- 2.837,+ 0.01298,- 0.3588,+ 0.02413,+ 0.03427,- 0.7624,+ 0.07272,+ 0.0216,- 0.003385,+ 0.4424,+ 0.03722,+ 0.02195,- 0.001503,+ 0.0002426,+ 0.3013,+ 0.05762,+ 0.01019,- 0.0004476,+ 6.79e-05,+ 3.227e-05,+ 0.3123,- 0.03535,+ 0.00484,+ 3.025e-06,- 4.363e-05,+ 2.854e-07,- 1.286e-06,- 0.6725,- 0.0373,+ 0.0008964,+ 0.0001399,- 3.99e-06,+ 7.431e-06,- 2.796e-07,- 1.601e-07,+ 0.04068,- 0.01352,+ 0.0007282,+ 9.594e-05,+ 2.07e-06,- 9.62e-08,- 2.742e-07,- 6.37e-08,- 6.3e-09,+ 0.08625,- 0.005971,+ 0.0004705,+ 2.335e-05,+ 4.226e-06,+ 2.475e-07,- 8.85e-08,- 3.6e-08,- 2.9e-09,+ 0.0])
    bh_amp = np.array([+ 0.0,+ 0.0,- 0.1136,+ 0.0,- 0.1868,- 0.01399,+ 0.0,- 0.1043,+ 0.01175,- 0.00224,+ 0.0,- 0.03222,+ 0.01333,- 0.002647,- 2.316e-05,+ 0.0,+ 0.05339,+ 0.01107,- 0.003116,- 0.0001079,- 1.299e-05,+ 0.0,+ 0.004861,+ 0.008891,- 0.0006448,- 1.279e-05,+ 6.358e-06,- 1.417e-07,+ 0.0,+ 0.03041,+ 0.00115,- 0.0008743,- 2.781e-05,+ 6.367e-07,- 1.14e-08,- 4.2e-08,+ 0.0,- 0.02982,- 0.003,+ 1.394e-05,- 3.29e-05,- 1.705e-07,+ 7.44e-08,+ 2.72e-08,- 6.6e-09,+ 0.0,+ 0.01236,- 0.0009981,- 3.792e-05,- 1.355e-05,+ 1.162e-06,- 1.789e-07,+ 1.47e-08,- 2.4e-09,- 4e-10])
    aw_mean = np.array([+ 56.4,+ 1.555,- 1.011,- 3.975,+ 0.03171,+ 0.1065,+ 0.6175,+ 0.1376,+ 0.04229,+ 0.003028,+ 1.688,- 0.1692,+ 0.05478,+ 0.02473,+ 0.0006059,+ 2.278,+ 0.006614,- 0.0003505,- 0.006697,+ 0.0008402,+ 0.0007033,- 3.236,+ 0.2184,- 0.04611,- 0.01613,- 0.001604,+ 5.42e-05,+ 7.922e-05,- 0.2711,- 0.4406,- 0.03376,- 0.002801,- 0.000409,- 2.056e-05,+ 6.894e-06,+ 2.317e-06,+ 1.941,- 0.2562,+ 0.01598,+ 0.005449,+ 0.0003544,+ 1.148e-05,+ 7.503e-06,- 5.667e-07,- 3.66e-08,+ 0.8683,- 0.05931,- 0.001864,- 0.0001277,+ 0.0002029,+ 1.269e-05,+ 1.629e-06,+ 9.66e-08,- 1.015e-07,- 5e-10])
    bw_mean = np.array([+ 0.0,+ 0.0,+ 0.2592,+ 0.0,+ 0.02974,- 0.5471,+ 0.0,- 0.5926,- 0.103,- 0.01567,+ 0.0,+ 0.171,+ 0.09025,+ 0.02689,+ 0.002243,+ 0.0,+ 0.3439,+ 0.02402,+ 0.00541,+ 0.001601,+ 9.669e-05,+ 0.0,+ 0.09502,- 0.03063,- 0.001055,- 0.0001067,- 0.000113,+ 2.124e-05,+ 0.0,- 0.3129,+ 0.008463,+ 0.0002253,+ 7.413e-05,- 9.376e-05,- 1.606e-06,+ 2.06e-06,+ 0.0,+ 0.2739,+ 0.001167,- 2.246e-05,- 0.0001287,- 2.438e-05,- 7.561e-07,+ 1.158e-06,+ 4.95e-08,+ 0.0,- 0.1344,+ 0.005342,+ 0.0003775,- 6.756e-05,- 1.686e-06,- 1.184e-06,+ 2.768e-07,+ 2.73e-08,+ 5.7e-09])
    aw_amp = np.array([+ 0.1023,- 2.695,+ 0.3417,- 0.1405,+ 0.3175,+ 0.2116,+ 3.536,- 0.1505,- 0.0166,+ 0.02967,+ 0.3819,- 0.1695,- 0.07444,+ 0.007409,- 0.006262,- 1.836,- 0.01759,- 0.06256,- 0.002371,+ 0.0007947,+ 0.0001501,- 0.8603,- 0.136,- 0.03629,- 0.003706,- 0.0002976,+ 1.857e-05,+ 3.021e-05,+ 2.248,- 0.1178,+ 0.01255,+ 0.001134,- 0.0002161,- 5.817e-06,+ 8.836e-07,- 1.769e-07,+ 0.7313,- 0.1188,+ 0.01145,+ 0.001011,+ 0.0001083,+ 2.57e-06,- 2.14e-06,- 5.71e-08,+ 2e-08,- 1.632,- 0.006948,- 0.003893,+ 0.0008592,+ 7.577e-05,+ 4.539e-06,- 3.852e-07,- 2.213e-07,- 1.37e-08,+ 5.8e-09])
    bw_amp = np.array([+ 0.0,+ 0.0,- 0.08865,+ 0.0,- 0.4309,+ 0.0634,+ 0.0,+ 0.1162,+ 0.06176,- 0.004234,+ 0.0,+ 0.253,+ 0.04017,- 0.006204,+ 0.004977,+ 0.0,- 0.1737,- 0.005638,+ 0.0001488,+ 0.0004857,- 0.0001809,+ 0.0,- 0.1514,- 0.01685,+ 0.005333,- 7.611e-05,+ 2.394e-05,+ 8.195e-06,+ 0.0,+ 0.09326,- 0.01275,- 0.0003071,+ 5.374e-05,- 3.391e-05,- 7.436e-06,+ 6.747e-07,+ 0.0,- 0.08637,- 0.003807,- 0.0006833,- 3.861e-05,- 2.268e-05,+ 1.454e-06,+ 3.86e-07,- 1.068e-07,+ 0.0,- 0.02658,- 0.001947,+ 0.0007131,- 3.506e-05,+ 1.885e-07,+ 5.792e-07,+ 3.99e-08,+ 2e-08,- 5.7e-09])
    np.pi = 3.14159265359
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
    V=pd.DataFrame(np.zeros((nmax+1,nmax+1)))
    W=pd.DataFrame(np.zeros((nmax+1,nmax+1)))
    
    
    V.loc[0,0] = 1
    W.loc[0,0] = 0
    V.loc[1,0] = z * V.loc[0,0]
    W.loc[1,0] = 0
    for n in range(2,nmax+1):
        V.loc[n ,0] = ((2 * n - 1) * z * V.loc[n-1,0] - (n - 1) * V.loc[n-2,0]) / n
        W.loc[n ,0] = 0
    
    for m in range(1,nmax+1):
        V.loc[m ,m ] = (2 * m - 1) * (x * V.loc[m-1,m-1]  - y * W.loc[m-1,m-1] )
        W.loc[m ,m ] = (2 * m - 1) * (x * W.loc[m-1,m-1] +  y * V.loc[m-1,m-1])
        if (m < nmax):
            V.loc[m + 1,m ] = (2 * m + 1) * z * V.loc[m,m]
            W.loc[m + 1,m ] = (2 * m + 1) * z * W.loc[m,m]
        for n in range(m+2,nmax+1):
            V.loc[n ,m ] = ((2 * n - 1) * z * V.loc[n-1,m] - (n + m - 1) * V.loc[n-2,m]) / (n - m)
            W.loc[n ,m ] = ((2 * n - 1) * z * W.loc[n-1,m] - (n + m - 1) * W.loc[n-2,m]) / (n - m)
    
    # (1) hydrostatic mf
    
    bh = 0.0029
    c0h = 0.062
    if (dlat < 0):
        phh = np.pi
        c11h = 0.007
        c10h = 0.002
    else:
        phh = 0
        c11h = 0.005
        c10h = 0.001
    
    ch = c0h + ((np.cos(doy / 365.25 * 2 * np.pi + phh) + 1) * c11h / 2 + c10h) * (1 - np.cos(dlat))
    ahm = 0
    aha = 0
    i = 0
    for n in range(0,nmax+1):
        for m in range(0,n+1):
            
            ahm = ahm + (ah_mean[i] * V.loc[n,m] + bh_mean[i]  * W.loc[n,m])
            aha = aha + (ah_amp[i]  * V.loc[n,m] + bh_amp[i]  * W.loc[n,m])
            i = i + 1
    
    ah = (ahm + aha * np.cos(doy / 365.25 * 2 * np.pi)) * 1e-05
    sine = np.sin(np.pi / 2 - zd)
    cose = np.cos(np.pi / 2 - zd)
    beta = bh / (sine + ch)
    gamma = ah / (sine + beta)
    topcon = (1 + ah / (1 + bh / (1 + ch)))
    gmfh = topcon / (sine + gamma)
    # height correction for hydrostatic mapping function from Niell (1996) in order to reduce the coefficients to sea level
    a_ht = 2.53e-05
    b_ht = 0.00549
    c_ht = 0.00114
    hs_km = dhgt / 1000
    beta = b_ht / (sine + c_ht)
    gamma = a_ht / (sine + beta)
    topcon = (1 + a_ht / (1 + b_ht / (1 + c_ht)))
    ht_corr_coef = 1 / sine - topcon / (sine + gamma)
    ht_corr = ht_corr_coef * hs_km
    gmfh = gmfh + ht_corr
    # (2) wet mf
    
    bw = 0.00146
    cw = 0.04391
    awm = 0
    awa = 0
    i = 0
    for n in range(0,nmax+1):
        for m in range(0,n+1):
            
            awm = awm + (aw_mean[i] * V.loc[n,m] + bw_mean[i] * W.loc[n,m])
            awa = awa + (aw_amp [i] * V.loc[n,m] + bw_amp[i] * W.loc[n,m])
            i = i + 1
    
    aw = (awm + awa * np.cos(doy / 365.25 * 2 * np.pi)) * 1e-05
    beta = bw / (sine + cw)
    gamma = aw / (sine + beta)
    topcon = (1 + aw / (1 + bw / (1 + cw)))
    gmfw = topcon / (sine + gamma)
    return gmfh,gmfw
 
def gpt3_5_fast_readGrid(filename ='gpt3_5.grd'):
# =============================================================================  
#     gpt3_5_fast_readGrid.py
#     This routine reads the grid 'gpt3_5.grd' and overgives the respective
#     parameters in a cell. It must be placed ahead of the for loop, which runs
#    through all observations in order to save a huge time amount. The related
#    function "gpt3_5_fast.py" must then replace the default
#   "gpt3_5.m" at the respective location in the text.
    
    
#    File created by Daniel Landskron, 2017-04-12
    
#    File Converted to python script by Zohreh Adavi zohreh.adavi@tuwien.ac.at, 2023-02-20
# =============================================================================  
    
# read gridfile

    
    file=open(filename)
    C = pd.read_csv(file ,skiprows=[i for i in range(0,1)],sep='\s+',skipfooter =0, skip_blank_lines=True,header = None, index_col= False )
    file.close()   


    
    p_grid = C.loc[:,np.arange(2,7)]
    
    T_grid = C.loc[:,np.arange(7,12)]
    
    Q_grid = C.loc[:,np.arange(12,17)] / 1000
    
    dT_grid= C.loc[:,np.arange(17,22)] / 1000
    
    u_grid = C.loc[:,22]
    
    Hs_grid = C.loc[:,23]
    
    ah_grid = C.loc[:,np.arange(24,29)] / 1000 
    
    aw_grid = C.loc[:,np.arange(29,34)] / 1000 
    
    la_grid = C.loc[:,np.arange(34,39)]
    
    Tm_grid = C.loc[:,np.arange(39,44)]
    
    Gn_h_grid = C.loc[:,np.arange(44,49)] / 100000
    
    Ge_h_grid = C.loc[:,np.arange(49,54)] / 100000
    
    Gn_w_grid = C.loc[:,np.arange(54,59)] / 100000 
    
    Ge_w_grid = C.loc[:,np.arange(59,64)] / 100000
    
    # combine all data to on cell grid
    gpt3_grid={}
    gpt3_grid[0]=p_grid
    gpt3_grid[1]=T_grid
    gpt3_grid[2]=Q_grid
    gpt3_grid[3]=dT_grid
    gpt3_grid[4]=u_grid
    gpt3_grid[5]=Hs_grid
    gpt3_grid[6]=ah_grid
    gpt3_grid[7]=aw_grid
    gpt3_grid[8]=la_grid
    gpt3_grid[9]=Tm_grid
    gpt3_grid[10]=Gn_h_grid
    gpt3_grid[11]=Ge_h_grid
    gpt3_grid[12]=Gn_w_grid
    gpt3_grid[13]=Ge_w_grid
            
    return gpt3_grid



    
def gpt3_5_fast(mjd = None,lat = None,lon = None,h_ell = None,it = None,grid = None): 
    
#    (c) Department of Geodesy and Geoinformation, Vienna University of Technology, 2017
    
# =============================================================================  
#     This subroutine determines pressure, temperature, temperature lapse rate,
#     mean temperature of the water vapor, water vapour pressure, hydrostatic
#     and wet mapping function coefficients ah and aw, water vapour decrease
#     factor, geoid undulation and empirical tropospheric gradients for
#     specific sites near the earth's surface. All output values are valid for
#     the specified ellipsoidal height h_ell.
#     GPT3_5 is based on a 5°x5° external grid file ('gpt3_5.grd') with mean
#     values as well as sine and cosine amplitudes for the annual and
#     semiannual variation of the coefficients.
#     This file goes together with 'gpt3_5_fast_readGrid.m', which reads the
#     grid and saves it in cell structures, what, especially for large numbers
#     of stations, significantly reduces runtime compared to gpt3_5.m.
#     gpt3_5_fast_readGrid' must be placed ahead of the for loop, in which
#     gpt3_5_fast is positioned.
    
    
#     Reference:
#     D. Landskron, J. Böhm (2018), VMF3/GPT3: Refined Discrete and Empirical Troposphere Mapping Functions,
#     J Geod (2018) 92: 349., doi: 10.1007/s00190-017-1066-2.
#     Download at: https://link.springer.com/content/pdf/10.1007#2Fs00190-017-1066-2.pdf
    
    
#     Input parameters:
#     -------------------    
#     mjd:   modified Julian date (scalar, only one epoch per call is possible)
#     lat:   ellipsoidal latitude in radians [-pi/2:+pi/2] (vector)
#     lon:   longitude in radians [-pi:pi] or [0:2pi] (vector)
#     h_ell: ellipsoidal height in m (vector)
#     it:    case 1: no time variation but static quantities
#     case 0: with time variation (annual and semiannual terms)
    
#     Output parameters:
#     -------------------
#     p:    pressure in hPa (vector)
#     T:    temperature in degrees Celsius (vector)
#     dT:   temperature lapse rate in degrees per km (vector)
#     Tm:   mean temperature weighted with the water vapor in degrees Kelvin (vector)
#     e:    water vapour pressure in hPa (vector)
#     ah:   hydrostatic mapping function coefficient (VMF3) (vector)
#     aw:   wet mapping function coefficient (VMF3) (vector)
#     la:   water vapour decrease factor (vector)
#     undu: geoid undulation in m (vector)
#     Gn_h: hydrostatic north gradient in m (vector)
#     Ge_h: hydrostatic east gradient in m (vector)
#     Gn_w: wet north gradient in m (vector)
#     Ge_w: wet east gradient in m (vector)
    
#     File created by Daniel Landskron, 2017-04-12
#     File Converted to python script by Zohreh Adavi zohreh.adavi@tuwien.ac.at, 2023-02-20
# =============================================================================  
    
# convert mjd to doy
    
    hour = int(np.floor((mjd - int(np.floor(mjd))) * 24))
    
    minu = int(np.floor((((mjd - int(np.floor(mjd))) * 24) - hour) * 60))
    
    sec = (((((mjd - int(np.floor(mjd))) * 24) - hour) * 60) - minu) * 60
    
    # change secs, min hour whose sec==60
    
    
    if sec==60:
        minu=minu+1
        sec=0
        
    if minu==60:
        hour=hour+1
        minu=0
    # calc jd (yet wrong for hour==24)
    jd = mjd + 2400000.5
    # if hr==24, correct jd and set hour==0
    if hour==24:
        jd=jd+1
        hour=0
    
    # integer julian date
    jd_int = int(np.floor(jd + 0.5))
    aa = jd_int + 32044
    bb = int(np.floor((4 * aa + 3) / 146097))
    cc = aa - int(np.floor((bb * 146097) / 4))
    dd = int(np.floor((4 * cc + 3) / 1461))
    ee = cc - int(np.floor((1461 * dd) / 4))
    mm = int(np.floor((5 * ee + 2) / 153))
    day = ee - int(np.floor((153 * mm + 2) / 5)) + 1
    month = mm + 3 - 12 * int(np.floor(mm / 10))
    year = bb * 100 + dd - 4800 + int(np.floor(mm / 10))
    # first check if the specified year is leap year or not (logical output)
    leapYear = (np.logical_or((np.mod(year,4) == np.logical_and(0,np.mod(year,100)) != 0),np.mod(year,400)) == 0)
    days = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    doy = sum(days[range(0,month-1)]) + day
    if leapYear == 1 and month > 2:
        doy = doy + 1
    
    doy = doy + mjd - int(np.floor(mjd))
    
    # determine the GPT3 coefficients
    
    # mean gravity in m/s**2
    gm = 9.80665
    # molar mass of dry air in kg/mol
    dMtr = 28.965 * 10 ** - 3
    # universal gas constant in J/K/mol
    Rg = 8.3143
    # factors for amplitudes
    if (it == 1):
        cosfy = 0
        coshy = 0
        sinfy = 0
        sinhy = 0
    else:
        cosfy = np.cos(doy / 365.25 * 2 * np.pi)
        coshy = np.cos(doy / 365.25 * 4 * np.pi)
        sinfy = np.sin(doy / 365.25 * 2 * np.pi)
        sinhy = np.sin(doy / 365.25 * 4 * np.pi)
    
    # read the respective data from the grid
    p_grid = grid[0]
    
    T_grid = grid[1]
    
    Q_grid = grid[2]
    
    dT_grid = grid[3]
    
    u_grid = grid[4]
    
    Hs_grid = grid[5]
    
    ah_grid = grid[6]
    
    aw_grid = grid[7]
    
    la_grid = grid[8]
    
    Tm_grid = grid[9]
    
    Gn_h_grid = grid[10]
    
    Ge_h_grid = grid[11]
    
    Gn_w_grid = grid[12]
    
    Ge_w_grid = grid[13]
    
    # determine the number of stations
    nstat = len(lat)
    # initialization
    p    = np.zeros(np.array([nstat,1]))
    T    = np.zeros(np.array([nstat,1]))
    dT   = np.zeros(np.array([nstat,1]))
    Tm   = np.zeros(np.array([nstat,1]))
    e    = np.zeros(np.array([nstat,1]))
    ah   = np.zeros(np.array([nstat,1]))
    aw   = np.zeros(np.array([nstat,1]))
    la   = np.zeros(np.array([nstat,1]))
    undu = np.zeros(np.array([nstat,1]))
    Gn_h = np.zeros(np.array([nstat,1]))
    Ge_h = np.zeros(np.array([nstat,1]))
    Gn_w = np.zeros(np.array([nstat,1]))
    Ge_w = np.zeros(np.array([nstat,1]))
    # loop over stations
    for k in range(0,nstat):
        # only positive longitude in degrees
        if lon[k] < 0:
            plon = (lon[k] + 2 * np.pi) * 180 / np.pi
        else:
            plon = lon[k] * 180 / np.pi
        # transform to polar distance in degrees
        ppod = (- lat[k] + np.pi / 2) * 180 / np.pi
        # find the index (line in the grid file) of the nearest point
        ipod = int(np.floor((ppod + 5) / 5))
        ilon = int(np.floor((plon + 5) / 5))
        # normalized (to one) differences, can be positive or negative
        diffpod = (ppod - (ipod * 5 - 2.5)) / 5
        difflon = (plon - (ilon * 5 - 2.5)) / 5
        if ipod == 36:
            ipod = 35
        if ilon == 72:
            ilon = 0
        if ilon == -1:
            ilon = 71
        indx=[]
        # get the number of the corresponding line
        indx.append((ipod - 1) * 72 + ilon-1)
        # near the poles: nearest neighbour interpolation, otherwise: bilinear
        bilinear = 0
        if ppod > 2.5 and ppod < 177.5:
            bilinear = 1
        # case of nearest neighbourhood
        if bilinear == 0:
            ix = indx[0]
            # transforming ellipsoidal height to orthometric height
            undu[k] = u_grid[ix]
            hgt = h_ell[k] - undu[k]
            # pressure, temperature at the height of the grid
            T0 = T_grid.iloc[ix,0] + T_grid.iloc[ix,1] * cosfy + T_grid.iloc[ix,2] * sinfy + T_grid.iloc[ix,3] * coshy + T_grid.iloc[ix,4] * sinhy
            p0 = p_grid.iloc[ix,0] + p_grid.iloc[ix,1] * cosfy + p_grid.iloc[ix,2] * sinfy + p_grid.iloc[ix,3] * coshy + p_grid.iloc[ix,4] * sinhy
            # specific humidity
            Q = Q_grid.iloc[ix,0] + Q_grid.iloc[ix,1] * cosfy + Q_grid.iloc[ix,2] * sinfy + Q_grid.iloc[ix,3] * coshy + Q_grid.iloc[ix,4] * sinhy
            # lapse rate of the temperature
            dT[k] = dT_grid.iloc[ix,0] + dT_grid.iloc[ix,1] * cosfy + dT_grid.iloc[ix,2] * sinfy + dT_grid.iloc[ix,3] * coshy + dT_grid.iloc[ix,4] * sinhy
            # station height - grid height
            redh = hgt - Hs_grid[ix]
            # temperature at station height in Celsius
            T[k] = T0 + dT[k] * redh - 273.15
            # temperature lapse rate in degrees / km
            dT[k] = dT[k] * 1000
            # virtual temperature in Kelvin
            Tv = T0 * (1 + 0.6077 * Q)
            c = gm * dMtr / (Rg * Tv)
            # pressure in hPa
            p[k] = (p0 * np.exp(- c * redh)) / 100
            # hydrostatic and wet coefficients ah and aw
            ah[k] = ah_grid.iloc[ix,0] + ah_grid.iloc[ix,1] * cosfy + ah_grid.iloc[ix,2] * sinfy + ah_grid.iloc[ix,3] * coshy + ah_grid.iloc[ix,4] * sinhy
            aw[k] = aw_grid.iloc[ix,0] + aw_grid.iloc[ix,1] * cosfy + aw_grid.iloc[ix,2] * sinfy + aw_grid.iloc[ix,3] * coshy + aw_grid.iloc[ix,4] * sinhy
            # water vapour decrease factor la
            la[k] = la_grid.iloc[ix,0] + la_grid.iloc[ix,1] * cosfy + la_grid.iloc[ix,2] * sinfy + la_grid.iloc[ix,3] * coshy + la_grid.iloc[ix,4] * sinhy
            # mean temperature Tm
            Tm[k] = Tm_grid.iloc[ix,0] + Tm_grid.iloc[ix,1] * cosfy + Tm_grid.iloc[ix,2] * sinfy + Tm_grid.iloc[ix,3] * coshy + Tm_grid.iloc[ix,4] * sinhy
            # north and east gradients (total, hydrostatic and wet)
            Gn_h[k] = Gn_h_grid.iloc[ix,0] + Gn_h_grid.iloc[ix,1] * cosfy + Gn_h_grid.iloc[ix,2] * sinfy + Gn_h_grid.iloc[ix,3] * coshy + Gn_h_grid.iloc[ix,4] * sinhy
            Ge_h[k] = Ge_h_grid.iloc[ix,0] + Ge_h_grid.iloc[ix,1] * cosfy + Ge_h_grid.iloc[ix,2] * sinfy + Ge_h_grid.iloc[ix,3] * coshy + Ge_h_grid.iloc[ix,4] * sinhy
            Gn_w[k] = Gn_w_grid.iloc[ix,0] + Gn_w_grid.iloc[ix,1] * cosfy + Gn_w_grid.iloc[ix,2] * sinfy + Gn_w_grid.iloc[ix,3] * coshy + Gn_w_grid.iloc[ix,4] * sinhy
            Ge_w[k] = Ge_w_grid.iloc[ix,0] + Ge_w_grid.iloc[ix,1] * cosfy + Ge_w_grid.iloc[ix,2] * sinfy + Ge_w_grid.iloc[ix,3] * coshy + Ge_w_grid.iloc[ix,4] * sinhy
            # water vapor pressure in hPa
            e0 = Q * p0 / (0.622 + 0.378 * Q) / 100
            e[k] = e0 * (100 * p[k] / p0) ** (la[k] + 1)
        else:
            ipod1 = int(ipod + np.sign(diffpod))
            ilon1 = int(ilon + np.sign(difflon))
            if ilon1 == 72:
                ilon1 = 0
            if ilon1 == -1:
                ilon1 = 71
            # get the number of the line
            indx.append((ipod1 - 1) * 72 + ilon-1)
            indx.append((ipod - 1) * 72 + ilon1-1)
            indx.append((ipod1 - 1) * 72 + ilon1-1)
            # transforming ellipsoidal height to orthometric height: Hortho = -N + Hell
            undul = u_grid[indx]
            hgt = h_ell[k] - undul
            # pressure, temperature at the height of the grid
            T0 = T_grid.iloc[indx,0] + T_grid.iloc[indx,1] * cosfy + T_grid.iloc[indx,2] * sinfy + T_grid.iloc[indx,3] * coshy + T_grid.iloc[indx,4] * sinhy
            p0 = p_grid.iloc[indx,0] + p_grid.iloc[indx,1] * cosfy + p_grid.iloc[indx,2] * sinfy + p_grid.iloc[indx,3] * coshy + p_grid.iloc[indx,4] * sinhy
            # humidity
            Ql = Q_grid.iloc[indx,0] + Q_grid.iloc[indx,1] * cosfy + Q_grid.iloc[indx,2] * sinfy + Q_grid.iloc[indx,3] * coshy + Q_grid.iloc[indx,4] * sinhy
            # reduction = stationheight - gridheight
            Hs1 = Hs_grid[indx]
            redh = hgt - Hs1
            # lapse rate of the temperature in degree / m
            dTl = dT_grid.iloc[indx,0] + dT_grid.iloc[indx,1] * cosfy + dT_grid.iloc[indx,2] * sinfy + dT_grid.iloc[indx,3] * coshy + dT_grid.iloc[indx,4] * sinhy
            # temperature reduction to station height
            Tl = T0 + np.multiply(dTl,redh) - 273.15
            # virtual temperature
            Tv = np.multiply(T0,(1 + 0.6077 * Ql))
            c = gm * dMtr / (Rg * Tv)
            # pressure in hPa
            pl = (np.multiply(p0,np.exp(np.multiply(- c,redh)))) / 100
            # hydrostatic and wet coefficients ah and aw
            ahl = ah_grid.iloc[indx,0] + ah_grid.iloc[indx,1] * cosfy + ah_grid.iloc[indx,2] * sinfy + ah_grid.iloc[indx,3] * coshy + ah_grid.iloc[indx,4] * sinhy
            awl = aw_grid.iloc[indx,0] + aw_grid.iloc[indx,1] * cosfy + aw_grid.iloc[indx,2] * sinfy + aw_grid.iloc[indx,3] * coshy + aw_grid.iloc[indx,4] * sinhy
            # water vapour decrease factor la
            lal = la_grid.iloc[indx,0] + la_grid.iloc[indx,1] * cosfy + la_grid.iloc[indx,2] * sinfy + la_grid.iloc[indx,3] * coshy + la_grid.iloc[indx,4] * sinhy
            # mean temperature of the water vapor Tm
            Tml = Tm_grid.iloc[indx,0] + Tm_grid.iloc[indx,1] * cosfy + Tm_grid.iloc[indx,2] * sinfy + Tm_grid.iloc[indx,3] * coshy + Tm_grid.iloc[indx,4] * sinhy
            # north and east gradients (total, hydrostatic and wet)
            Gn_hl = Gn_h_grid.iloc[indx,0] + Gn_h_grid.iloc[indx,1] * cosfy + Gn_h_grid.iloc[indx,2] * sinfy + Gn_h_grid.iloc[indx,3] * coshy + Gn_h_grid.iloc[indx,4] * sinhy
            Ge_hl = Ge_h_grid.iloc[indx,0] + Ge_h_grid.iloc[indx,1] * cosfy + Ge_h_grid.iloc[indx,2] * sinfy + Ge_h_grid.iloc[indx,3] * coshy + Ge_h_grid.iloc[indx,4] * sinhy
            Gn_wl = Gn_w_grid.iloc[indx,0] + Gn_w_grid.iloc[indx,1] * cosfy + Gn_w_grid.iloc[indx,2] * sinfy + Gn_w_grid.iloc[indx,3] * coshy + Gn_w_grid.iloc[indx,4] * sinhy
            Ge_wl = Ge_w_grid.iloc[indx,0] + Ge_w_grid.iloc[indx,1] * cosfy + Ge_w_grid.iloc[indx,2] * sinfy + Ge_w_grid.iloc[indx,3] * coshy + Ge_w_grid.iloc[indx,4] * sinhy
            # water vapor pressure in hPa
            e0 = np.multiply(Ql,p0) / (0.622 + 0.378 * Ql) / 100
            el = np.multiply(e0,(100.0 * pl / p0) ** (lal + 1))
            dnpod1 = np.abs(diffpod)
            dnpod2 = 1 - dnpod1
            dnlon1 = np.abs(difflon)
            dnlon2 = 1 - dnlon1
            # pressure
            R1 = dnpod2 * pl.iloc[0] + dnpod1 * pl.iloc[1]
            R2 = dnpod2 * pl.iloc[2] + dnpod1 * pl.iloc[3]
            p[k] = dnlon2 * R1 + dnlon1 * R2
            # temperature
            R1 = dnpod2 * Tl.iloc[0] + dnpod1 * Tl.iloc[1]
            R2 = dnpod2 * Tl.iloc[2] + dnpod1 * Tl.iloc[3]
            T[k] = dnlon2 * R1 + dnlon1 * R2
            # temperature in degree per km
            R1 = dnpod2 * dTl.iloc[0] + dnpod1 * dTl.iloc[1]
            R2 = dnpod2 * dTl.iloc[2] + dnpod1 * dTl.iloc[3]
            dT[k] = (dnlon2 * R1 + dnlon1 * R2) * 1000
            # water vapor pressure in hPa
            R1 = dnpod2 * el.iloc[0] + dnpod1 * el.iloc[1]
            R2 = dnpod2 * el.iloc[2] + dnpod1 * el.iloc[3]
            e[k] = dnlon2 * R1 + dnlon1 * R2
            # ah and aw
            R1 = dnpod2 * ahl.iloc[0] + dnpod1 * ahl.iloc[1]
            R2 = dnpod2 * ahl.iloc[2] + dnpod1 * ahl.iloc[3]
            ah[k] = dnlon2 * R1 + dnlon1 * R2
            R1 = dnpod2 * awl.iloc[0] + dnpod1 * awl.iloc[1]
            R2 = dnpod2 * awl.iloc[2] + dnpod1 * awl.iloc[3]
            aw[k] = dnlon2 * R1 + dnlon1 * R2
            # undulation
            R1 = dnpod2 * undul.iloc[0] + dnpod1 * undul.iloc[1]
            R2 = dnpod2 * undul.iloc[2] + dnpod1 * undul.iloc[3]
            undu[k] = dnlon2 * R1 + dnlon1 * R2
            # water vapor decrease factor la
            R1 = dnpod2 * lal.iloc[0] + dnpod1 * lal.iloc[1]
            R2 = dnpod2 * lal.iloc[2] + dnpod1 * lal.iloc[3]
            la[k] = dnlon2 * R1 + dnlon1 * R2
            # gradients
            R1 = dnpod2 * Gn_hl.iloc[0] + dnpod1 * Gn_hl.iloc[1]
            R2 = dnpod2 * Gn_hl.iloc[2] + dnpod1 * Gn_hl.iloc[3]
            Gn_h[k] = (dnlon2 * R1 + dnlon1 * R2)
            R1 = dnpod2 * Ge_hl.iloc[0] + dnpod1 * Ge_hl.iloc[1]
            R2 = dnpod2 * Ge_hl.iloc[2] + dnpod1 * Ge_hl.iloc[3]
            Ge_h[k] = (dnlon2 * R1 + dnlon1 * R2)
            R1 = dnpod2 * Gn_wl.iloc[0] + dnpod1 * Gn_wl.iloc[1]
            R2 = dnpod2 * Gn_wl.iloc[2] + dnpod1 * Gn_wl.iloc[3]
            Gn_w[k] = (dnlon2 * R1 + dnlon1 * R2)
            R1 = dnpod2 * Ge_wl.iloc[0] + dnpod1 * Ge_wl.iloc[1]
            R2 = dnpod2 * Ge_wl.iloc[2] + dnpod1 * Ge_wl.iloc[3]
            Ge_w[k] = (dnlon2 * R1 + dnlon1 * R2)
            # mean temperature of the water vapor Tm
            R1 = dnpod2 * Tml.iloc[0] + dnpod1 * Tml.iloc[1]
            R2 = dnpod2 * Tml.iloc[2] + dnpod1 * Tml.iloc[3]
            Tm[k] = dnlon2 * R1 + dnlon1 * R2
    
    return p,T,dT,Tm,e,ah,aw,la,undu,Gn_h,Ge_h,Gn_w,Ge_w

