#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 13:39:48 2019

@author: psakicki
"""

########## BEGIN IMPORT ##########
#### External modules
import numpy as np
import pandas as pd
import datetime as dt
#from pytwobodyorbit import TwoBodyOrbit
#### geodeZYX modules
from geodezyx import conv
from geodezyx import utils
from geodezyx import files_rw
from geodezyx import operational
from geodezyx import time_series

##########  END IMPORT  ##########


def ECI_2_kepler_elts(P,V,rad2deg=True,
                      mu=3.9860044188e14):
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



def extrapolate_sp3_DataFrame(DFsp3,step=900,n_step=9):
    NewEpoch_stk = []
    
    for sat in DFsp3["sat"].unique():
        print(sat)
        DFsat = DFsp3[(DFsp3["sat"] == sat) & (DFsp3["type"] == "P")]
        
        DFline_dummy = DFsat.iloc[0].copy()
        DFline_dummy["clk"] = 999999.999999
        
        XYZ   = DFsat[["x","y","z"]].values *1000
        Tgps = conv.datetime64_numpy2dt(DFsat["epoch"].values)
        Tutc = conv.dt_gpstime2dt_utc(Tgps)
        Tsec = [e.total_seconds() for e in (Tutc-Tutc[0])]
        
        P=conv.ECEF2ECI(XYZ,Tutc)
        V=np.gradient(P,Tsec,axis=0)
        
    
        def extrapo_intern_fct(i_ref,coef):
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
            for t in np.arange(step,n_step*step +1,step):
                Pout, Vout = orbit_back.posvelatt(coef * t)
                epoc = Tgps[i_ref] + coef * dt.timedelta(seconds=int(t))
                DFline_new = DFline_dummy.copy()
                DFline_new[["x","y","z"]] = Pout
                DFline_new["epoch"] = pd.Timestamp(epoc)
        
                NewEpoch_stk.append(DFline_new)
                
            return Pout, Vout
    
        ### backward
        extrapo_intern_fct(0,-1)
    
        ### forward
        extrapo_intern_fct(-1,1)
        
    DFNewEpoch = pd.DataFrame(NewEpoch_stk)
    
    Tutc_NewEpoch = conv.dt_gpstime2dt_utc(conv.datetime64_numpy2dt(DFNewEpoch["epoch"].values))
    DFNewEpoch[["x","y","z"]] = conv.ECI2ECEF(DFNewEpoch[["x","y","z"]].values,Tutc_NewEpoch)
    DFNewEpoch[["x","y","z"]] = DFNewEpoch[["x","y","z"]] * .001
    
    DFout = pd.concat((DFsp3,DFNewEpoch))
    DFout.reset_index(drop=True,inplace=True)
    DFout.sort_values(["sat","epoch"],inplace=True)
    DFout.reset_index(drop=True,inplace=True)
    
    return DFout


