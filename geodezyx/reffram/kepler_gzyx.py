#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.reffram contains functions related to
 Kepler's element operations. 

it can be imported directly with:
from geodezyx import reffram

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""

import datetime as dt
#### Import the logger
import logging

########## BEGIN IMPORT ##########
#### External modules
import numpy as np
import pandas as pd
# from pytwobodyorbit import TwoBodyOrbit
from pytwobodyorbit import TwoBodyOrbit

#### geodeZYX modules
from geodezyx import conv

log = logging.getLogger(__name__)

##########  END IMPORT  ##########


def ECI_2_kepler_elts(pos,vel,rad2deg=True,
                      mu=3.9860044188e14):
    """
    Convert ECI coordinates > Kepler's elements

    Parameters
    ----------
    pos : array
        Position in m. 
        Can be a 3-element (x,y,z) array (simple point) 
        or a (N,3)-shaped numpy array
    vel : array
        Velocity in m/s.
        Can be a 3-element (vx,vy,vz) array (simple point) 
        or a (N,3)-shaped numpy array
    rad2deg : bool, optional
        convert Keplerian's angles (anomalies) to deg. The default is True.
    mu : float, optional
        Standard gravitational parameter. The default is 3.9860044188e14.

    Returns
    -------
    a,ecc,i,o_peri,o_lan,m: floats
        Kepler's elements:  
            
        * a : semi-major axis
        * ecc : orbit eccentricity
        * i : orbit inclination
        * o_peri : argument of periapsis ω
        * o_lan : longitude of the ascending node Ω
        * m : mean anomaly m

    Note
    ----
    Be sure your input is in ECI (SP3 are in ECEF for instance).  

    If neeeded, use ``conv.ECEF2ECI()`` (gross results) 
    or ``reffram.OrbDF_crf2trf(inv_trf2crf=True)`` (precise results)
    
    You can get the velocity with ``reffram.DFOrb_velocity_calc()``
    
    Source
    ------
    https://downloads.rene-schwarz.com/download/M002-Cartesian_State_Vectors_to_Keplerian_Orbit_Elements.pdf

    """
    p = np.array(pos)
    v = np.array(vel)
    
    #### case for one point only
    if len(p.shape) < 2:
        p = np.expand_dims(p,axis=-1).T
    if len(v.shape) < 2:
        v = np.expand_dims(v,axis=-1).T
        
    pnorm = np.linalg.norm(p,axis=1)
    vnorm = np.linalg.norm(v,axis=1)
    
    ### orbital momentum vector h
    h = np.cross(p,v)
    hnorm = np.linalg.norm(h,axis=1)
    
    ### eccentricity vector eccvec & orbit eccentricity ecc
    eccvec = (np.cross(v,h)/mu) - (p/np.expand_dims(pnorm,axis=-1))
    ecc = np.linalg.norm(eccvec,axis=1)
    
    ### true anomaly ν (nu)
    n = np.cross(np.array([0,0,1]),h)
    nnorm = np.linalg.norm(n,axis=1)
    
    dot_pv_sup0 = np.einsum('ij,ij->i',p,v) >= 0
    nu_sup0 = np.arccos(np.einsum('ij,ij->i',eccvec,p)/(ecc*pnorm))
    nu_inf0 = 2*np.pi - nu_sup0
        
    nu = np.zeros(len(dot_pv_sup0))
    nu = nu + nu_sup0 * dot_pv_sup0
    nu = nu + nu_inf0 * np.logical_not(dot_pv_sup0)
    
    ### orbit inclination i
    i = np.arccos(h[:,2]/hnorm)

    ### eccentric anomaly E (called e here)
    #e1 = 2*np.arctan(np.tan(nu/2)/np.sqrt((1+ecc)/(1-ecc)))
    e2 = 2*np.arctan2(np.tan(nu/2),np.sqrt((1+ecc)/(1-ecc)))
    e=e2
    
    ### longitude of the ascending node Ω (o_lan)
    n1_sup0 = n[:,1] >= 0
        
    omega_lan_sup0 = np.arccos(n[:,0]/nnorm)
    omega_lan_inf0 = 2*np.pi - omega_lan_sup0
    
    o_lan = np.zeros(len(n1_sup0))
    o_lan = o_lan + omega_lan_sup0 * n1_sup0
    o_lan = o_lan + omega_lan_inf0 * np.logical_not(n1_sup0)
    
    ### argument of periapsis ω (o_peri)
    eccvec2_sup0 = eccvec[:,2] >= 0
    o_peri_sup0 = np.arccos((np.einsum('ij,ij->i',n, eccvec))/(ecc*nnorm))
    o_peri_inf0 = 2 * np.pi - o_peri_sup0
    
    o_peri = np.zeros(len(eccvec2_sup0))
    o_peri = o_peri + o_peri_sup0 * eccvec2_sup0
    o_peri = o_peri + o_peri_inf0 * np.logical_not(eccvec2_sup0)
    
    ### mean anomaly m, (Kepler’s Equation)
    m = e - ecc*np.sin(e)
    
    ### semi-major axis a
    a = ((2/pnorm) - (vnorm**2 / mu))**-1

    #### case for one point only (bring it back to a scalar)
    sqzflt = lambda x: np.float64(x.squeeze())

    a = sqzflt(a)
    ecc = sqzflt(ecc)
    i = sqzflt(i)
    o_lan = sqzflt(o_lan)
    o_peri = sqzflt(o_peri)
    m = sqzflt(m)

    if rad2deg:
        i=np.rad2deg(i)
        o_peri=np.rad2deg(o_peri)
        o_lan=np.rad2deg(o_lan)
        m=np.rad2deg(m)
        
    #### o_peri and m seems to be inverted, must be investigated!!!
    return a,ecc,i,o_peri,o_lan,m

def ECI_2_kepler_elts_mono(P,V,rad2deg=True,
                           mu=3.9860044188e14):
    """
    **Kept for debug purposes, use** ``ECI_2_kepler_elts(pos, vel)`` **instead**
    
    Convert ECI coordinates > Kepler's elements

    Parameters
    ----------
    P : array
        Position.
    V : array
        Velocity.
    rad2deg : bool, optional
        canvert Keplerian's angles (anomalies) to deg. The default is True.
    mu : float, optional
        Standard gravitational parameter. The default is 3.9860044188e14.

    Returns
    -------
    a,ecc,i,omega_periarg,omega_LAN,m : floats
        Kepler's elements.
    """
    log.warning('use this function for debug purposes only, use ECI_2_kepler_elts instead')
    Pnorm = np.linalg.norm(P)
    Vnorm = np.linalg.norm(V)
    
    ### orbital momentum vector H
    H = np.cross(P,V)
    Hnorm = np.linalg.norm(H)
    
    ### eccentricity vector Ecc
    Ecc = (np.cross(V,H)/mu) - (P/Pnorm)
    Ecc_norm = np.linalg.norm(Ecc)
    
    ### true anomaly ν (nu)
    N = np.cross(np.array([0,0,1]),H)
    Nnorm = np.linalg.norm(N)
    if np.dot(P,V) >= 0:
        nu = np.arccos(np.dot(Ecc,P)/(Ecc_norm*Pnorm))
    else:
        nu = 2*np.pi - np.arccos(np.dot(Ecc,P)/(Ecc_norm*Pnorm))        
    
    ### orbit inclination i
    i = np.arccos(H[2]/Hnorm)
    
    ### orbit eccentricity ecc
    ecc = np.linalg.norm(Ecc)
    
    ### eccentric anomaly E (called e here)
    #e1 = 2*np.arctan(np.tan(nu/2)/np.sqrt((1+ecc)/(1-ecc)))
    e2 = 2*np.arctan2(np.tan(nu/2),np.sqrt((1+ecc)/(1-ecc)))
    e=e2
    
    ### longitude of the ascending node Ω (omega_LAN)
    if N[1] >= 0:
        omega_LAN = np.arccos(N[0]/Nnorm)
    else:
        omega_LAN = 2*np.pi - np.arccos(N[0]/Nnorm)
    
    ### argument of periapsis ω (omega_periarg)
    if Ecc[2] >= 0:
        omega_periarg = np.arccos((np.dot(N,Ecc))/(Nnorm * Ecc_norm))
    else:
        omega_periarg = 2*np.pi - np.arccos((np.dot(N,Ecc))/(Nnorm * Ecc_norm))
    
    ### mean anomaly m, (Kepler’s Equation)
    m = e - ecc*np.sin(e)
    
    ### semi-major axis a
    a = ((2/Pnorm) - (Vnorm**2 / mu))**-1

    if rad2deg:
        i=np.rad2deg(i)
        omega_periarg=np.rad2deg(omega_periarg)
        omega_LAN=np.rad2deg(omega_LAN)
        m=np.rad2deg(m)

    return a,ecc,i,omega_periarg,omega_LAN,m

def extrapolate_orbit_kepler(P,V,t,t0,mu=3.9860044188e14):
    orbit = TwoBodyOrbit("", mu=mu)   # create an instance
    orbit.setOrbCart(t0, P, V)        # define the orbit
    Pout, Vout = orbit.posvelatt(t)   # get position and velocity at t1
    kepl = orbit.elmKepl()            # get classical orbital elements
    
    return Pout,Vout,kepl



def extrapolate_sp3_DataFrame(DFsp3,step=900,n_step=9,
                              backward=True,forward=True,
                              until_backward=None,
                              until_forward=None,
                              return_all=True):
    """
    Extrapolate the positions in a SP3 based on the first/last
    epochs' position

    Parameters
    ----------
    DFsp3 : DataFrame
        Input Orbit DataFrame (i.e. generated by files_rw.read_sp3).
    step : int, optional
        step between two epochs. The default is 900.
    n_step : int, optional
        number of epochs to interpolate. The default is 9.
    backward and forward : bool, optional
        extrapolate for the day before/after. The default is True.
    return_all : bool, optional
        if True, returns the input DF enhanced with the extrapolated values
        if False, returns only the extrapolated values.
        The default is True.
    until_backward & until_backward : datetime, optional
        epoch until then the extrapolation has to be done.
        Override n_step
        
    Returns
    -------
    DForb_out : DataFrame
        Orbit DataFrame with extrapolated values
        (see return_all).

    """
    NewEpoch_stk = []
    
    for sat in DFsp3["sat"].unique():
        log.info("extrapolate: %s",sat)
        DFsat = DFsp3[(DFsp3["sat"] == sat) & (DFsp3["type"] == "P")].copy()
        DFsat.sort_values("epoch",inplace=True)
        
        DFline_dummy = DFsat.iloc[0].copy()
        DFline_dummy["clk"] = 999999.999999
        
        XYZ   = DFsat[["x","y","z"]].values *1000
        Tgps = conv.numpy_dt2dt(DFsat["epoch"].values)
        Tutc = conv.dt_gpstime2dt_utc(Tgps)
        Tsec = [e.total_seconds() for e in (Tutc-Tutc[0])]
        
        P=conv.ECEF2ECI(XYZ,Tutc)
        V=np.gradient(P,Tsec,axis=0)
        
        def extrapo_intern_fct(i_ref,coef,until):
            ##
            ## backward :
            ## i_ref , coef = 0,-1
            ##
            ## forward :
            ## i_ref , coef = -1,1
            ##
            ## CORRDS ARE GIVEN IN ECI HERE !!!
            ##
            orbit_back = TwoBodyOrbit("orbit", mu=3.9860044188e14)  # create an instance    
            orbit_back.setOrbCart(Tsec[i_ref], P[i_ref], V[i_ref])  # define the orbit
            
            if until:
                t_rang_strt,t_rang_end = list(sorted([Tutc[i_ref],until]))
                Range = conv.dt_range(t_rang_strt,t_rang_end,0,step)
                n_step_intern = len(Range) - 1 
            else:
                n_step_intern = n_step
            
            for t in np.arange(step,n_step_intern*step +1,step):
                Pout, Vout = orbit_back.posvelatt(coef * t)
                epoc = Tgps[i_ref] + coef * dt.timedelta(seconds=int(t))
                DFline_new = DFline_dummy.copy()
                DFline_new[["x","y","z"]] = Pout
                DFline_new["epoch"] = pd.Timestamp(epoc)
        
                NewEpoch_stk.append(DFline_new)
                
            return Pout, Vout
    
        ### backward
        if backward:
            extrapo_intern_fct(0,-1,until_backward)
    
        ### forward
        if forward:
            extrapo_intern_fct(-1,1,until_forward)
        
    DFNewEpoch = pd.DataFrame(NewEpoch_stk)
    
    Tutc_NewEpoch = conv.dt_gpstime2dt_utc(conv.numpy_dt2dt(DFNewEpoch["epoch"].values))
    DFNewEpoch[["x","y","z"]] = conv.ECI2ECEF(DFNewEpoch[["x","y","z"]].values,Tutc_NewEpoch)
    DFNewEpoch[["x","y","z"]] = DFNewEpoch[["x","y","z"]] * 10**-3
    
    if return_all:
        DForb_out = pd.concat((DFsp3,DFNewEpoch))
    else:
        DForb_out = DFNewEpoch
    
    DForb_out.reset_index(drop=True,inplace=True)
    DForb_out.sort_values(["sat","epoch"],inplace=True)
    DForb_out.reset_index(drop=True,inplace=True)
    
    return DForb_out


