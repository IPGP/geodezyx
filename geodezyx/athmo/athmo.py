#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 09:23:54 2019

@author: Chaiyaporn Kitpracha
"""

from geodezyx import np,conv,re,df


def trop_saast(p,dlat,hell,t=0,e=0,mode="dry"):
    """
    This subroutine determines the zenith total delay based on the equation by Saastamoinen (1972) as refined by Davis et al. (1985)

    Parameters
    ----------
    p :
        pressure in hPa

    dlat :
        ellipsoidal latitude in radians 
        
    t :
         temperature in Celcius
         
    e : 
         water vapor pressure in hPa
                
    hell :
         ellipsoidal height in meters
         
    mode :
         dry, wet or total
         
    Returns
    -------
    res :
        zenith total delay in m (depend on mode) 
    
    Source
    ------
     c Reference:
     Saastamoinen, J., Atmospheric correction for the troposphere and 
     stratosphere in radio ranging of satellites. The use of artificial 
     satellites for geodesy, Geophys. Monogr. Ser. 15, Amer. Geophys. Union, 
     pp. 274-251, 1972.
     Davis, J.L, T.A. Herring, I.I. Shapiro, A.E.E. Rogers, and G. Elgered, 
     Geodesy by Radio Interferometry: Effects of Atmospheric Modeling Errors 
     on Estimates of Baseline Length, Radio Science, Vol. 20, No. 6, 
     pp. 1593-1607, 1985.
    
    """
    f = 1-0.00266*np.cos(2*dlat) - 0.00000028*hell # calculate denominator f
    t = t + 273.15 #convert celcius to kelvin 
    if mode == "dry":
        res = 0.0022768*p/f
    elif mode == "wet":
        res = (0.22768e-2) * ((0.1255e+4)+(t+0.5e-1)) * (e/t)
    elif mode == "total":
        zhd = 0.0022768*p/f
        zwd = (0.22768e-2) * ((0.1255e+4)+(t+0.5e-1)) * (e/t)
        res = zhd + zwd
        
    return np.round(res,4)

def read_grid_gpt(grid_name,cols=64):
    grid = np.loadtxt(grid_name,usecols=range(cols),comments="%",dtype=float)
    return grid

def calc_stand_ties(epoc, lat_ref , lon_ref , h_ref , lat_rov , lon_rov , h_rov ,grid,unit="mm"):
    """
    Determine standard atmospheric ties from meteological information from GPT3 with Saastamonien model
    
    Parameters:
    -------
        epoc : 
            time in Python datetime
        lat_ref : 
            Latitude of Ref. station
        lat_rov : 
            Latitude of Rov. station
        h_ref : 
            Height of Ref. station
        h_rov : 
            Height of Rov. station
        grid_file : 
            meteological grid file
        unit : 
            in meters (m) or milimeters (mm)
        
    Return:
    -------
        ties : 
            Standard ties of total delay in milimeters
        
    """
    
    
    met_ref = gpt3(epoc,lat_ref,lon_ref,h_ref,grid)
    met_rov = gpt3(epoc,lat_rov,lon_rov,h_rov,grid)
    
    ztd_ref = trop_saast(met_ref[0],lat_ref,h_ref,met_ref[1],met_ref[4],"total")
    ztd_rov = trop_saast(met_rov[0],lat_rov,h_ref,met_rov[1],met_rov[4],"total")
    if unit == "m":
        return ztd_ref - ztd_rov
    elif unit == "mm":
        return (ztd_ref - ztd_rov) * 1000

def gpt3(dtin,lat,lon,h_ell,C,it=0):
    
   
    """ 
    This subroutine determines pressure, temperature, temperature lapse rate, 
    mean temperature of the water vapor, water vapour pressure, hydrostatic 
    and wet mapping function coefficients ah and aw, water vapour decrease
    factor, geoid undulation and empirical tropospheric gradients for 
    specific sites near the earth's surface.
    It is based on a 5 x 5 degree external grid file ('gpt3_5.grd') with mean
    values as well as sine and cosine amplitudes for the annual and
    semiannual variation of the coefficients.
    
    Parameters:
    -------  
    dtin : 
        datatime in Python datetime object
    
    lat:   
        ellipsoidal latitude in radians [-pi/2:+pi/2]
    
    lon:   
        longitude in radians [-pi:pi] or [0:2pi]
    
    h_ell: 
        ellipsoidal height in m
    
    it:    
        case 1 no time variation but static quantities, case 0 with time variation (annual and semiannual terms)
           
    Returns:
    -------    
    p:    
        pressure in hPa  
    
    T:    
        temperature in degrees Celsius 
    
    dT:   
        temperature lapse rate in degrees per km
    
    Tm:   
        mean temperature weighted with the water vapor in degrees Kelvin 
    
    e:    
        water vapour pressure in hPa 
    
    ah:   
        hydrostatic mapping function coefficient at zero height (VMF3) 
    
    aw:   
        wet mapping function coefficient (VMF3) 
    
    la:   
        water vapour decrease factor 
    
    undu: 
        geoid undulation in m 
    
    Gn_h: 
        hydrostatic north gradient in m 
    
    Ge_h: 
        hydrostatic east gradient in m 
    
    Gn_w: 
        wet north gradient in m 
    
    Ge_w: 
        wet east gradient in m 
    
    Notes:
        Modified for Python by : Chaiyaporn Kitpracha
        
    Source:
        
        (c) Department of Geodesy and Geoinformation, Vienna University of
        Technology, 2017
      
        The copyright in this document is vested in the Department of Geodesy and
        Geoinformation (GEO), Vienna University of Technology, Austria. This document
        may only be reproduced in whole or in part, or stored in a retrieval
        system, or transmitted in any form, or by any means electronic,
        mechanical, photocopying or otherwise, either with the prior permission
        of GEO or in accordance with the terms of ESTEC Contract No.
        4000107329/12/NL/LvH.
        
        D. Landskron, J. Böhm (2018), VMF3/GPT3: Refined Discrete and Empirical Troposphere Mapping Functions, 
        J Geod (2018) 92: 349., doi: 10.1007/s00190-017-1066-2. 
        Download at: https://link.springer.com/content/pdf/10.1007%2Fs00190-017-1066-2.pdf
    """
    lat = np.array(lat)
    lon = np.array(lon)
    h_ell = np.array(h_ell)
    # Extract data from grid
    p_grid    = C[:,2:7]          # pressure in Pascal
    T_grid    = C[:,7:12]         # temperature in Kelvin
    Q_grid    = C[:,12:17]/1000   # specific humidity in kg/kg
    dT_grid   = C[:,17:22]/1000   # temperature lapse rate in Kelvin/m
    u_grid    = C[:,22]           # geoid undulation in m
    Hs_grid   = C[:,23]           # orthometric grid height in m
    ah_grid   = C[:,24:29]/1000   # hydrostatic mapping function coefficient, dimensionless
    aw_grid   = C[:,29:34]/1000   # wet mapping function coefficient, dimensionless
    la_grid   = C[:,34:39]     	  # water vapor decrease factor, dimensionless
    Tm_grid   = C[:,39:44]        # mean temperature in Kelvin
    Gn_h_grid = C[:,44:49]/100000    # hydrostatic north gradient in m
    Ge_h_grid = C[:,49:54]/100000    # hydrostatic east gradient in m
    Gn_w_grid = C[:,54:59]/100000    # wet north gradient in m
    Ge_w_grid = C[:,59:64]/100000    # wet east gradient in m
    
    # Convert from datetime to doy 
    doy = float(conv.dt2doy(dtin)) + conv.dt2fracday(dtin)
    
    # determine the GPT3 coefficients

    # mean gravity in m/s**2
    gm = 9.80665
    # molar mass of dry air in kg/mol
    dMtr = 28.965e-3
    # universal gas constant in J/K/mol
    Rg = 8.3143
    
    # factors for amplitudes
    if it==1: # then  constant parameters
        cosfy = 0
        coshy = 0
        sinfy = 0
        sinhy = 0
    else:
        cosfy = np.cos(doy/365.25*2*np.pi)   # coefficient for A1
        coshy = np.cos(doy/365.25*4*np.pi)   # coefficient for B1
        sinfy = np.sin(doy/365.25*2*np.pi)   # coefficient for A2
        sinhy = np.sin(doy/365.25*4*np.pi)   # coefficient for B2
    
    nstat = lat.size
    
    # initialization
    p    = np.zeros([nstat , 1])
    T    = np.zeros([nstat , 1])
    dT   = np.zeros([nstat , 1])
    Tm   = np.zeros([nstat , 1])
    e    = np.zeros([nstat , 1])
    ah   = np.zeros([nstat , 1])
    aw   = np.zeros([nstat , 1])
    la   = np.zeros([nstat , 1])
    undu = np.zeros([nstat , 1])
    Gn_h = np.zeros([nstat , 1])
    Ge_h = np.zeros([nstat , 1])
    Gn_w = np.zeros([nstat , 1])
    Ge_w = np.zeros([nstat , 1])
        
    if lon < 0:
        plon = (lon + 2*np.pi)*180/np.pi
    else:
        plon = lon*180/np.pi
        
    ppod = (-lat + np.pi/2)*180/np.pi
    
    ipod = np.floor(ppod+1) 
    ilon = np.floor(plon+1)
    
    # changed for the 1 degree grid
    diffpod = (ppod - (ipod - 0.5))
    difflon = (plon - (ilon - 0.5))
    
    if ipod == 181:
        ipod = 180
        
    if ilon == 361:
        ilon = 1

    if ilon == 0:
        ilon = 360
    
    indx = np.zeros(4)
    indx[0] = (ipod - 1)*360 + ilon
   
    # near the poles: nearest neighbour interpolation, otherwise: bilinear
    # with the 1 degree grid the limits are lower and upper
    bilinear = 0
    if ppod > 0.5 and ppod < 179.5:
           bilinear = 1
    
    
    
    if bilinear == 0:
        ix = int(indx[0]) - 1
        
        # transforming ellipsoidal height to orthometric height
        undu = u_grid[ix]
        hgt = h_ell-undu
            
        # pressure, temperature at the height of the grid
        T0 = T_grid[ix,0] + T_grid[ix,1]*cosfy + T_grid[ix,2]*sinfy + T_grid[ix,3]*coshy + T_grid[ix,4]*sinhy
        p0 = p_grid[ix,0] + p_grid[ix,1]*cosfy + p_grid[ix,2]*sinfy + p_grid[ix,3]*coshy + p_grid[ix,4]*sinhy
         
        # specific humidity
        Q = Q_grid[ix,0] + Q_grid[ix,1]*cosfy + Q_grid[ix,2]*sinfy + Q_grid[ix,3]*coshy + Q_grid[ix,4]*sinhy
            
        # lapse rate of the temperature
        dT = dT_grid[ix,0] + dT_grid[ix,1]*cosfy + dT_grid[ix,2]*sinfy + dT_grid[ix,3]*coshy + dT_grid[ix,4]*sinhy 

        # station height - grid height
        redh = hgt - Hs_grid[ix]

        # temperature at station height in Celsius
        T = T0 + dT*redh - 273.15
        
        # temperature lapse rate in degrees / km
        dT = dT*1000

        # virtual temperature in Kelvin
        Tv = T0*[1+0.6077*Q]
        
        c = gm*dMtr/[Rg*Tv]
        
        # pressure in hPa
        p = [p0*np.exp[-c*redh]]/100
        
            
        # hydrostatic and wet coefficients ah and aw 
        ah = ah_grid[ix,0] + ah_grid[ix,2]*cosfy + ah_grid[ix,3]*sinfy + ah_grid[ix,4]*coshy + ah_grid[ix,5]*sinhy
        aw = aw_grid[ix,0] + aw_grid[ix,2]*cosfy + aw_grid[ix,3]*sinfy + aw_grid[ix,4]*coshy + aw_grid[ix,5]*sinhy
		
		# water vapour decrease factor la
        la = la_grid[ix,0] + \
                la_grid[ix,1]*cosfy + la_grid[ix,2]*sinfy + \
                la_grid[ix,3]*coshy + la_grid[ix,4]*sinhy 
		
		# mean temperature Tm
        Tm = Tm_grid[ix,0] + \
                Tm_grid[ix,1]*cosfy + Tm_grid[ix,2]*sinfy + \
                Tm_grid[ix,3]*coshy + Tm_grid[ix,4]*sinhy
            
        # north and east gradients [total, hydrostatic and wet]
        Gn_h = Gn_h_grid[ix,0] + Gn_h_grid[ix,1]*cosfy + Gn_h_grid[ix,2]*sinfy + Gn_h_grid[ix,3]*coshy + Gn_h_grid[ix,4]*sinhy
        Ge_h = Ge_h_grid[ix,0] + Ge_h_grid[ix,1]*cosfy + Ge_h_grid[ix,2]*sinfy + Ge_h_grid[ix,3]*coshy + Ge_h_grid[ix,4]*sinhy
        Gn_w = Gn_w_grid[ix,0] + Gn_w_grid[ix,1]*cosfy + Gn_w_grid[ix,2]*sinfy + Gn_w_grid[ix,3]*coshy + Gn_w_grid[ix,4]*sinhy
        Ge_w = Ge_w_grid[ix,0] + Ge_w_grid[ix,1]*cosfy + Ge_w_grid[ix,2]*sinfy + Ge_w_grid[ix,3]*coshy + Ge_w_grid[ix,4]*sinhy
		
		# water vapor pressure in hPa
        e0 = Q*p0/[0.622+0.378*Q]/100 # on the grid
        e = e0*(100*p/p0)**(la+1)   # on the station height - (14] Askne and Nordius, 1987
        
    else:
        ipod1 = ipod + 1*np.sign(diffpod)
        ilon1 = ilon + 1*np.sign(difflon)
       
		# changed for the 1 degree grid
        if ilon1 == 361:
            ilon1 = 1
        
        if ilon1 == 0:
            ilon1 = 360
        
        # get the number of the line
		# changed for the 1 degree grid
        indx[1] = (ipod1 - 1)*360 + ilon  # along same longitude
        indx[2] = (ipod  - 1)*360 + ilon1 # along same polar distance
        indx[3] = (ipod1 - 1)*360 + ilon1 # diagonal
        indx = indx.astype(int)
        indx = indx - 1
       
        # transforming ellipsoidal height to orthometric height: Hortho = -N + Hell
        undul = u_grid[indx]
        hgt = h_ell-undul
        
        # pressure, temperature at the height of the grid
        T0 = T_grid[indx,0] + T_grid[indx,1]*cosfy + T_grid[indx,2]*sinfy + T_grid[indx,3]*coshy + T_grid[indx,4]*sinhy
        p0 = p_grid[indx,0] + p_grid[indx,1]*cosfy + p_grid[indx,2]*sinfy + p_grid[indx,3]*coshy + p_grid[indx,4]*sinhy
        
        # humidity
        Ql = Q_grid[indx,0] + Q_grid[indx,1]*cosfy + Q_grid[indx,2]*sinfy + Q_grid[indx,3]*coshy + Q_grid[indx,4]*sinhy
        
        # reduction = stationheight - gridheight
        Hs1 = Hs_grid[indx]
        redh = hgt - Hs1
        
        # lapse rate of the temperature in degree / m
        dTl = dT_grid[indx,0] + dT_grid[indx,1]*cosfy + dT_grid[indx,2]*sinfy + dT_grid[indx,3]*coshy + dT_grid[indx,4]*sinhy
        
        # temperature reduction to station height
        Tl = T0 + dTl * redh - 273.15
        
        # virtual temperature
        Tv = T0 *(1+0.6077*Ql)
        c = gm*dMtr / (Rg*Tv)
        
        # pressure in hPa
        pl = (p0 *np.exp(-c * redh))/100
        
        # hydrostatic and wet coefficients ah and aw
        ahl = ah_grid[indx,0] + ah_grid[indx,1]*cosfy + ah_grid[indx,2]*sinfy + ah_grid[indx,3]*coshy + ah_grid[indx,4]*sinhy
        awl = aw_grid[indx,0] + aw_grid[indx,1]*cosfy + aw_grid[indx,2]*sinfy + aw_grid[indx,3]*coshy + aw_grid[indx,4]*sinhy
        
        # water vapour decrease factor la
        lal = la_grid[indx,0] + la_grid[indx,1]*cosfy + la_grid[indx,2]*sinfy + la_grid[indx,3]*coshy + la_grid[indx,4]*sinhy
        
        # mean temperature of the water vapor Tm
        Tml = Tm_grid[indx,0] + Tm_grid[indx,1]*cosfy + Tm_grid[indx,2]*sinfy + Tm_grid[indx,3]*coshy + Tm_grid[indx,4]*sinhy
        
        # north and east gradients [total, hydrostatic and wet]
        Gn_hl = Gn_h_grid[indx,0] + Gn_h_grid[indx,1]*cosfy + Gn_h_grid[indx,2]*sinfy + Gn_h_grid[indx,3]*coshy + Gn_h_grid[indx,4]*sinhy
        Ge_hl = Ge_h_grid[indx,0] + Ge_h_grid[indx,1]*cosfy + Ge_h_grid[indx,2]*sinfy + Ge_h_grid[indx,3]*coshy + Ge_h_grid[indx,4]*sinhy
        Gn_wl = Gn_w_grid[indx,0] + Gn_w_grid[indx,1]*cosfy + Gn_w_grid[indx,2]*sinfy + Gn_w_grid[indx,3]*coshy + Gn_w_grid[indx,4]*sinhy
        Ge_wl = Ge_w_grid[indx,0] + Ge_w_grid[indx,1]*cosfy + Ge_w_grid[indx,2]*sinfy + Ge_w_grid[indx,3]*coshy + Ge_w_grid[indx,4]*sinhy
        
        # water vapor pressure in hPa
        e0 = Ql *p0 /(0.622+0.378*Ql)/100 # on the grid
        el =  e0 * (100 *pl /p0)**(lal+1)  # on the station height - [14] Askne and Nordius, 1987
            
        dnpod1 = abs(diffpod) # distance nearer point
        dnpod2 = 1 - dnpod1   # distance to distant point
        dnlon1 = abs(difflon)
        dnlon2 = 1 - dnlon1
        
        # pressure
        R1 = dnpod2*pl[0]+dnpod1*pl[1]
        R2 = dnpod2*pl[2]+dnpod1*pl[3]
        p = dnlon2*R1+dnlon1*R2
            
        # temperature
        R1 = dnpod2*Tl[0]+dnpod1*Tl[1]
        R2 = dnpod2*Tl[2]+dnpod1*Tl[3]
        T = dnlon2*R1+dnlon1*R2
        
        # temperature in degree per km
        R1 = dnpod2*dTl[0]+dnpod1*dTl[1]
        R2 = dnpod2*dTl[2]+dnpod1*dTl[3]
        dT = (dnlon2*R1+dnlon1*R2)*1000
            
        # water vapor pressure in hPa
        R1 = dnpod2*el[0]+dnpod1*el[1]
        R2 = dnpod2*el[2]+dnpod1*el[3]
        e = dnlon2*R1+dnlon1*R2
            
        # ah and aw
        R1 = dnpod2*ahl[0]+dnpod1*ahl[1]
        R2 = dnpod2*ahl[2]+dnpod1*ahl[3]
        ah = dnlon2*R1+dnlon1*R2
        R1 = dnpod2*awl[0]+dnpod1*awl[1]
        R2 = dnpod2*awl[2]+dnpod1*awl[3]
        aw = dnlon2*R1+dnlon1*R2
        
        # undulation
        R1 = dnpod2*undul[0]+dnpod1*undul[1]
        R2 = dnpod2*undul[2]+dnpod1*undul[3]
        undu = dnlon2*R1+dnlon1*R2
		
		# water vapor decrease factor la
        R1 = dnpod2*lal[0]+dnpod1*lal[1]
        R2 = dnpod2*lal[2]+dnpod1*lal[3]
        la = dnlon2*R1+dnlon1*R2
        
        # gradients
        R1 = dnpod2*Gn_hl[0]+dnpod1*Gn_hl[1]
        R2 = dnpod2*Gn_hl[2]+dnpod1*Gn_hl[3]
        Gn_h = (dnlon2*R1 + dnlon1*R2)
        R1 = dnpod2*Ge_hl[0]+dnpod1*Ge_hl[1]
        R2 = dnpod2*Ge_hl[2]+dnpod1*Ge_hl[3]
        Ge_h = (dnlon2*R1 + dnlon1*R2)
        R1 = dnpod2*Gn_wl[0]+dnpod1*Gn_wl[1]
        R2 = dnpod2*Gn_wl[2]+dnpod1*Gn_wl[3]
        Gn_w = (dnlon2*R1 + dnlon1*R2)
        R1 = dnpod2*Ge_wl[0]+dnpod1*Ge_wl[1]
        R2 = dnpod2*Ge_wl[2]+dnpod1*Ge_wl[3]
        Ge_w = (dnlon2*R1 + dnlon1*R2)
		
		# mean temperature of the water vapor Tm
        R1 = dnpod2*Tml[0]+dnpod1*Tml[1]
        R2 = dnpod2*Tml[2]+dnpod1*Tml[3]
        Tm = dnlon2*R1+dnlon1*R2
    
    soln = [np.round(p,3),np.round(T,3),np.round(dT,3),np.round(Tm,3),np.round(e,3), \
            np.round(ah,3),np.round(aw,3),np.round(la,3),np.round(undu,3),np.round(Gn_h,3),np.round(Ge_h,3), \
            np.round(Gn_w,3),np.round(Ge_w,3)]   
    return soln

def read_snx_trop(snxfile,dataframe_output=True):
    """
    Read troposphere solutions from Troposphere SINEX
    """
    
    STAT , epoc = [] , []
    tro , stro , tgn , stgn , tge , stge = [] , [] , [] , [] , [] , []
    
    flagtrop = False
    
    for line in open(snxfile,"r",encoding = "ISO-8859-1"):
        
        if re.compile('TROP/SOLUTION').search(line):
            flagtrop = True
            continue
        
        if re.compile('-TROP/SOLUTION').search(line):
            flagtrop = False
            continue
        
        if line[0] != ' ':
            continue
        else:
            fields = line.split()
        
        if flagtrop ==True:
            
            STAT.append(fields[0].upper())
            if not ':' in fields[1]:
                epoc.append(conv.convert_partial_year(fields[1]))
            else:
                date_elts_lis = fields[1].split(':')
                yy =  int(date_elts_lis[0]) + 2000
                doy = int(date_elts_lis[1])
                sec = int(date_elts_lis[2])

                epoc.append(conv.doy2dt(yy,doy,seconds=sec))
            
            tro.append(float(fields[2]))
            stro.append(float(fields[3]))
            tgn.append(float(fields[4]))
            stgn.append(float(fields[5]))
            tge.append(float(fields[6]))
            stge.append(float(fields[7]))
            
    outtuple = \
    list(zip(*sorted(zip(STAT , epoc , tro , stro , tgn , stgn , tge , stge))))
    
    if dataframe_output:
        return Tropsinex_DataFrame(outtuple)
                
def Tropsinex_DataFrame(read_sinex_result):
     """
      General description

        Parameters
        ----------
        read_sinex_result : 
            List values from read_snx_trop function
                    
        Returns
        -------
        DF_Sinex : 
            troposphere information from SINEX in Dataframe
            
     """
     DF_Sinex = df.DataFrame.from_records(list(read_sinex_result)).transpose()
     colnam = ['STAT','epoc','tro','stro','tgn','stgn','tge','stge']
     DF_Sinex.columns = colnam
     cols_numeric = ['tro','stro','tgn','stgn','tge','stge']
     DF_Sinex[cols_numeric] = DF_Sinex[cols_numeric].apply(df.to_numeric, errors='coerce')
     
     return DF_Sinex



####################### Compare troposphere delay  ####################################################
def compare_trop_ties(input_file,STA1,STA2,coord_file="",grid_met="",apply_ties=False,mode="DF",mode_coor="sinex",coord_t="static"):
    """
    Calculate differences of tropospheric delay and gradients between selected stations (Atmospheric ties)
    Args  :
          sinex_file : troposphere solutions from sinex
          STA1 : Reference station
          STA2 : Rover station
          coord_file : Station coordinates file path
          grid_met : Grid file for meteological information from Global Pressture Temperature (GPT)
          apply_ties : Apply height coorections from standard ties (Default No)
          mode : data in DataFrame (DF) or SINEX (SINEX) (Default DataFrame)
          mode_coor : Coordinates file format in SINEX or EPOS or DataFrame format (Argruments : sinex , epos , df)
          coord_t : static (default) or kinematic
    Return :
         trop_diff : difference of tropospheric delay and gradients between selected stations (Atmospheric ties)
                     and uncertainty of atmospheric ties and gradients ties
    """

    if mode == "SINEX":
        trop_pd = gfc.read_snx_trop(input_file)
    elif mode == "DF":
        trop_pd = input_file
    else:
        import sys
        print("No this option for troposphere solution")
        sys.exit()

    trop_ref = trop_pd[trop_pd.STAT == STA1]
    trop_rov = trop_pd[trop_pd.STAT == STA2]

    if trop_ref.empty:
        print("No solution for "+STA1+" Reference station")
        return None,None

    if trop_rov.empty:
        print("No solution for "+STA2+" Rover station")
        return None,None

    diff_pd = pd.merge(trop_ref,trop_rov,how='outer',on='epoc')
    diff_pd = diff_pd.dropna()
    # Tropospheric ties
    diff_pd['Trop_ties'] = diff_pd['tro_x'] - diff_pd['tro_y']
    diff_pd['STrop_ties'] = np.nan # add blank column before input values
    diff_pd['STrop_ties'] = diff_pd.apply(lambda x: np.round(np.sqrt(x.stro_x**2 + x.stro_y**2),2),axis=1)

    # North gradients ties
    diff_pd['Ngra_ties'] = diff_pd['tgn_x'] - diff_pd['tgn_y']
    diff_pd['SNgra_ties'] = np.nan # add blank column before input values
    diff_pd['SNgra_ties'] = diff_pd.apply(lambda x: np.round(np.sqrt(x.stgn_x**2 + x.stgn_y**2),2),axis=1)

    # East gradients ties
    diff_pd['Egra_ties'] = diff_pd['tge_x'] - diff_pd['tge_y']
    diff_pd['SEgra_ties'] = np.nan # add blank column before input values
    diff_pd['SEgra_ties'] = diff_pd.apply(lambda x: np.round(np.sqrt(x.stge_x**2 + x.stge_y**2),2),axis=1)

     # Apply height corrections from Standard ties
    if apply_ties == True:

        if isinstance(grid_met,str):
            import sys
            print("Please read grid file before use this function")
            sys.exit()

        if mode_coor == "epos" and coord_t == "static":
            coord = read_epos_sta_coords_mono(coord_file)
            # Extract coordinates Ref station in Lat. Lon. height
            coord_ref = coord[coord.site==STA1]
            lat_ref , lon_ref , h_ref = geok.XYZ2GEO(coord_ref.x,coord_ref.y,coord_ref.z,False)
            # Extract coordinates Rov station in Lat. Lon. height
            coord_rov = coord[coord.site==STA2]
            lat_rov , lon_rov , h_rov = geok.XYZ2GEO(coord_rov.x,coord_rov.y,coord_rov.z,False)

            #Merge coordinates results
            coord_res = pd.merge(coord_ref,coord_rov,how='outer',on='MJD_ref')
            coord_res['lat_ref'] , coord_res['lon_ref'] , coord_res['h_ref'] = float(lat_ref) , float(lon_ref) , float(h_ref)
            coord_res['lat_rov'] , coord_res['lon_rov'] , coord_res['h_rov'] = float(lat_rov) , float(lon_rov) , float(h_rov)

            #drop unnecessary column
            coord_res = coord_res.drop([ 'site_num_x', 'tecto_plate_x', 'MJD_start_x',
                                        'MJD_end_x', 'Vx_x',
                                        'Vy_x', 'Vz_x', 'sVx_x', 'sVy_x', 'sVz_x', 'site_num_y',
                                        'tecto_plate_y', 'MJD_start_y', 'MJD_end_y', 'Vx_y', 'Vy_y', 'Vz_y', 'sVx_y', 'sVy_y',
                                        'sVz_y'],axis=1)

        elif mode_coor == "sinex" and coord_t == "static":
            coord = gfc.read_sinex(coord_file,True)
             # Extract coordinates Ref station in Lat. Lon. height
            coord_ref = coord[coord.STAT==STA1]
            lat_ref , lon_ref , h_ref = geok.XYZ2GEO(coord_ref.x,coord_ref.y,coord_ref.z,False)
            # Extract coordinates Rov station in Lat. Lon. height
            coord_rov = coord[coord.STAT==STA2]
            lat_rov , lon_rov , h_rov = geok.XYZ2GEO(coord_rov.x,coord_rov.y,coord_rov.z,False)

            #Merge coordinates results
            coord_res = pd.merge(coord_ref,coord_rov,how='outer',on='epoc')
            coord_res['lat_ref'] , coord_res['lon_ref'] , coord_res['h_ref'] = float(lat_ref) , float(lon_ref) , float(h_ref)
            coord_res['lat_rov'] , coord_res['lon_rov'] , coord_res['h_rov'] = float(lat_rov) , float(lon_rov) , float(h_rov)

             #drop unnecessary columns
            coord_res = coord_res.drop(['AC_x', 'soln_x', 'vx_x', 'vy_x', 'vz_x', 'svx_x', 'svy_x', 'svz_x', 'start_x',
                                       'end_x', 'AC_y', 'soln_y', 'vx_y', 'vy_y', 'vz_y', 'svx_y', 'svy_y', 'svz_y',
                                       'start_y', 'end_y'],axis=1)

        elif mode_coor == "epos" and coord_t == "kinematic":
            coord = read_epos_sta_kinematics(coord_file)
            # Extract coordinates Ref station in Lat. Lon. height
            coord_ref = coord[coord.site==STA1.lower()]
            f_x_ref , f_y_ref , f_z_ref = interpolator_with_extrapolated(coord_ref.MJD_epo,coord_ref.x,coord_ref.y,coord_ref.z)
            x_ref_new = f_x_ref(geok.dt2MJD(diff_pd.epoc))
            y_ref_new = f_y_ref(geok.dt2MJD(diff_pd.epoc))
            z_ref_new = f_z_ref(geok.dt2MJD(diff_pd.epoc))
            lat_ref , lon_ref , h_ref = geok.XYZ2GEO(x_ref_new,y_ref_new,z_ref_new,False)

            # Extract coordinates Rov station in Lat. Lon. height
            coord_rov = coord[coord.site==STA2.lower()]
            f_x_rov , f_y_rov , f_z_rov = interpolator_with_extrapolated(coord_rov.MJD_epo,coord_rov.x,coord_rov.y,coord_rov.z)
            x_rov_new = f_x_rov(geok.dt2MJD(diff_pd.epoc))
            y_rov_new = f_y_rov(geok.dt2MJD(diff_pd.epoc))
            z_rov_new = f_z_rov(geok.dt2MJD(diff_pd.epoc))
            lat_rov , lon_rov , h_rov = geok.XYZ2GEO(x_rov_new,y_rov_new,z_rov_new,False)

            diff_pd['lat_ref'] , diff_pd['lon_ref'] , diff_pd['h_ref'] = lat_ref , lon_ref , h_ref
            diff_pd['lat_rov'] , diff_pd['lon_rov'] , diff_pd['h_rov'] = lat_rov , lon_rov , h_rov

            #Merge coordinates results
            coord_res = diff_pd[['STAT_x','STAT_y','epoc','lat_ref','lon_ref','h_ref','lat_rov','lon_rov','h_rov']].copy()
            coord_res = coord_res.rename(index=str,columns={"STAT_x":"STAT_ref","STAT_y":"STAT_rov"})
        else:
            import sys
            print("No this option for coordinates")
            sys.exit()
        #Extract standard ties

        grid = grid_met
        if coord_t == "kinematic":
            diff_pd['stand_ties'] = diff_pd.apply(lambda x:athmo.calc_stand_ties(x['epoc'], x.lat_ref , x.lon_ref , x.h_ref,x.lat_rov , x.lon_rov , x.h_rov,grid),axis=1)
        elif coord_t == "static":
            diff_pd['stand_ties'] = diff_pd.apply(lambda x: gtro.calc_stand_ties(x['epoc'], float(lat_ref) , float(lon_ref) , float(h_ref) , float(lat_rov) , float(lon_rov) , float(h_rov),grid),axis=1)


        diff_pd['Trop_ties_corr'] = diff_pd.apply(lambda x: x['tro_x'] - (x['tro_y']+x['stand_ties']),axis=1)

    #drop unnecessary column
    if coord_t == "kinematic":
        diff_pd = diff_pd.drop(['tro_x', 'stro_x', 'tgn_x', 'stgn_x', 'tge_x', 'stge_x','tro_y', 'stro_y', 'tgn_y', 'stgn_y', 'tge_y','stge_y','lat_ref','lon_ref','h_ref','lat_rov','lon_rov','h_rov'],axis=1)
    else:
        diff_pd = diff_pd.drop(['tro_x', 'stro_x', 'tgn_x', 'stgn_x', 'tge_x', 'stge_x','tro_y', 'stro_y', 'tgn_y', 'stgn_y', 'tge_y','stge_y'],axis=1)

    #Change column name of station
    diff_pd = diff_pd.rename(index=str,columns={"STAT_x":"STAT_ref","STAT_y":"STAT_rov"})

    return diff_pd,coord_res

def stat_summary_trop_ties(df):
    wmean_no_ties = np.round(np.average(df.Trop_ties,weights = 1/df.STrop_ties),3)
    wmean_wt_ties = np.round(np.average(df.Trop_ties_corr,weights = 1/df.STrop_ties),3)
    rms_mean_no_ties = np.round(geok.rms_mean(df.Trop_ties),3)
    rms_mean_wt_ties = np.round(geok.rms_mean(df.Trop_ties_corr),3)

    return [wmean_no_ties,wmean_wt_ties,rms_mean_no_ties,rms_mean_wt_ties]

def plot_trop_ties(df,ref_sta,rov_sta,analy_coor=False,analy_num_obs=False,df_coord="",savePlot=False,filePath="",fileName=""):
    """
    Plot tropospheric ties function
    Input:
        df : DataFrame from "compare_trop_ties" function
        ref_sta : Refernce station
        rov_sta : Rover station
        analy_coor : Add height difference information from "compare_trop_ties"
        analy_num_obs : Add plot number of observation from "compare_trop_ties"
        df_coord : list of height difference
        savePlot : save figure
        filePath : Directory to save
        fileName : Filename of figure
    """
    epo_plt = geok.dt2year_decimal(df.epoc)
    h_diff = df_coord.h_ref - df_coord.h_rov
    if analy_coor:

        fig, ax = plt.subplots(3,1,sharex=True)
        axA = ax[0]
        axB = ax[1]
        axC = ax[2]
        axA.plot(epo_plt,df.Trop_ties , marker="P",linestyle="--",label="Trop.ties")
        axA.plot(np.unique(epo_plt), np.poly1d(np.polyfit(epo_plt, df.Trop_ties, 1))(np.unique(epo_plt)),label="Trop. ties fitline")
        if 'Trop_ties_corr' in df.columns:
            axA.plot(epo_plt,df.Trop_ties_corr,marker="*",linestyle="-.",label="Trop. ties apply height corr.")
            axA.plot(np.unique(epo_plt), np.poly1d(np.polyfit(epo_plt, df.Trop_ties_corr, 1))(np.unique(epo_plt)),label="Trop. ties apply height fitline")
        axA.grid()
        axA.legend()
        axA.set_ylabel("Trop. ties (mm)")
        axA.set_xlabel("Time")
        axA.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))

        axB.plot(epo_plt,h_diff,marker="P",linestyle="--",label="Height difference")
        axB.plot(np.unique(epo_plt), np.poly1d(np.polyfit(epo_plt, h_diff, 1))(np.unique(epo_plt)),label="Height diff. fitline")
        axB.grid()
        axB.legend()
        axB.set_ylabel("Height difference (m)")
        axB.set_xlabel("Time")
        axB.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
        if analy_num_obs:
            if len(epo_plt) != len(df.num_obs_ref):
                print("No plot number of observations")
                fig.delaxes(axC)
            else:
                axC.plot(epo_plt,df.num_obs_ref,marker="P",linestyle="--",label="Num obs. Ref sta")
                axC.plot(epo_plt,df.num_obs_rov,marker="*",linestyle="-.",label="Num obs. Rov sta")
                axC.grid()
                axC.legend()
                axC.set_ylabel("Num Obs.")
                axC.set_xlabel("Time")
                axC.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
    else:
        fig , ax = plt.subplots(2,1,sharex=True)
        axA = ax[0]
        axC = ax[1]
        axA.plot(epo_plt,df.Trop_ties , marker="P",linestyle="--",label="Trop.ties")
        axA.plot(np.unique(epo_plt), np.poly1d(np.polyfit(epo_plt, df.Trop_ties, 1))(np.unique(epo_plt)),label="Trop. ties fitline")
        if 'Trop_ties_corr' in df.columns:
            axA.plot(epo_plt,df.Trop_ties_corr,marker="*",linestyle="-.",label="Trop. ties apply height corr.")
            axA.plot(np.unique(epo_plt), np.poly1d(np.polyfit(epo_plt, df.Trop_ties_corr, 1))(np.unique(epo_plt)),label="Trop. ties apply height fitline")
        axA.grid()
        axA.legend()
        axA.set_ylabel("Trop. ties (mm)")
        axA.set_xlabel("Time")
        axA.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))

        if analy_num_obs:
            if len(epo_plt) != len(df.num_obs_ref):
                print("No plot number of observations")
                fig.delaxes(axC)
            else:
                axC.plot(epo_plt,df.num_obs_ref,marker="P",linestyle="--",label="Num obs. Ref sta")
                axC.plot(epo_plt,df.num_obs_rov,marker="*",linestyle="-.",label="Num obs. Rov sta")
                axC.grid()
                axC.legend()
                axC.set_ylabel("Num Obs.")
                axC.set_xlabel("Time")
                axC.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))

    plt.tight_layout()
    plt.suptitle("Total delay ties of " + ref_sta + "-" + rov_sta)
    if savePlot:
        export_ts_figure_pdf(fig,filePath,fileName,True)

    plt.show()

    return None




