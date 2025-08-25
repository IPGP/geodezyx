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

import geodezyx.gnss_edu as gnss_edu


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
    C = pd.read_csv(file ,skiprows=[i for i in range(0,1)],sep=r'\s+',skipfooter =0, skip_blank_lines=True,header = None, index_col= False )
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

