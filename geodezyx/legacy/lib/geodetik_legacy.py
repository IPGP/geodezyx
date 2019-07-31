# -*- coding: utf-8 -*-
"""
Created on Fri May 30 12:50:26 2014

@author: psakicki


legacy version of geodetik, before
split in smaller module (20190406)
"""

import genefun
import softs_runner
import transformations as trans

import numpy as np
import datetime as dt
import time
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
import re
import os
import multiprocessing as mp
import struct
import string
import itertools
import platform
from scipy.signal import butter, lfilter, freqz
import warnings



#if platform.node() in ('calipso' , "diamant"):
    #import cgkitmod.cgtypes as cgt
    ##import cgtypes as cgt

import sys
sys.dont_write_bytecode = True



### Rotation matrix

def C_2D(theta,angtype='deg'):

    if angtype == 'deg':
        theta = np.deg2rad(theta)

    return np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta), np.cos(theta)]])

def C_ned2enu():
    return np.array([[0,1,0],[1,0,0],[0,0,-1]])

def C_enu2ned():
    return np.array([[0,1,0],[1,0,0],[0,0,-1]])

def C_ned2ecef(phi,llambda,angtype='deg'):

    if angtype == 'deg':
        phi = np.deg2rad(phi)
        llambda = np.deg2rad(llambda)

    In_ned = np.array([1,0,0])
    Ie_ned = np.array([0,1,0])
    Id_ned = np.array([0,0,1])

    Ix_ned = np.array([-np.cos(llambda) * np.sin(phi),-np.sin(llambda),-np.cos(llambda)*np.cos(phi)])
    Iy_ned = np.array([-np.sin(llambda) * np.sin(phi),np.cos(llambda),-np.sin(llambda)*np.cos(phi)])
    Iz_ned = np.array([np.cos(phi),0,-np.sin(phi)])

    C_ned2ecef = np.array([[np.dot(Ix_ned,In_ned),np.dot(Ix_ned,Ie_ned),np.dot(Ix_ned,Id_ned)],
                    [np.dot(Iy_ned,In_ned),np.dot(Iy_ned,Ie_ned),np.dot(Iy_ned,Id_ned)],
                    [np.dot(Iz_ned,In_ned),np.dot(Iz_ned,Ie_ned),np.dot(Iz_ned,Id_ned)]])

    return C_ned2ecef

def C_ecef2ned(phi,llambda,angtype='deg'):
    return np.linalg.inv(C_ned2ecef(phi,llambda,angtype))


def C_enu2ecef(phi,llambda,angtype='deg'):
    if angtype == 'deg':
        phi = np.deg2rad(phi)
        llambda = np.deg2rad(llambda)


    Ie_enu = np.array([1,0,0])
    In_enu = np.array([0,1,0])
    Iu_enu = np.array([0,0,1])

    Ix_enu = np.array([-np.sin(llambda),-np.cos(llambda) * np.sin(phi),np.cos(llambda)*np.cos(phi)])
    Iy_enu = np.array([np.cos(llambda),-np.sin(llambda)*np.sin(phi),np.sin(llambda)*np.cos(phi)])
    Iz_enu = np.array([0,np.cos(phi),np.sin(phi)])


    C_enu2ecef = np.array([[np.dot(Ix_enu,Ie_enu),np.dot(Ix_enu,In_enu),np.dot(Ix_enu,Iu_enu)],
                [np.dot(Iy_enu,Ie_enu),np.dot(Iy_enu,In_enu),np.dot(Iy_enu,Iu_enu)],
                [np.dot(Iz_enu,Ie_enu),np.dot(Iz_enu,In_enu),np.dot(Iz_enu,Iu_enu)]])

    return C_enu2ecef

def C_ecef2enu(phi,llambda,angtype='deg'):
    return C_enu2ecef(phi,llambda,angtype).T


def C_rpy2enu(roll,pitch,yaw,angtype='deg'):

    if angtype == 'deg':
        roll = np.deg2rad(roll)
        pitch = np.deg2rad(pitch)
        yaw = np.deg2rad(yaw)

    R = roll
    P = pitch
    Y = yaw

    Ir = np.array([np.sin(Y)*np.cos(P),np.cos(Y)*np.cos(P),np.sin(P)])
    Ip = np.array([np.cos(R) * np.cos(Y) + np.sin(R)*np.sin(Y)*np.sin(P),
                   -np.cos(R)*np.sin(Y)+np.sin(R)*np.cos(Y)*np.sin(P),
                   -np.sin(R)*np.cos(P)])
    Iy = np.array([-np.sin(R) *  np.cos(Y) + np.cos(R) * np.sin(Y) * np.sin(P) ,
                   np.sin(R) * np.sin(Y) + np.cos(R) * np.cos(Y) * np.sin(P) ,
                   - np.cos(R) * np.cos(P)    ])

#    Ie = np.array([ np.sin(Y) * np.cos(P) ,
#    np.cos(R) * np.cos(Y) + np.sin(R) * np.sin(Y) * np.sin(P) ,
#    -np.sin(R) * np.cos(Y) + np.cos(R) * np.sin(Y) * np.sin(P) ])
#
#    In = np.array([ np.cos(Y) * np.cos(P),
#     -np.cos(R) * np.sin(Y) + np.sin(R) * np.cos(Y) * np.sin(P),
#     np.sin(R) * np.sin(Y) + np.cos(R) * np.cos(Y) * np.sin(P) ])
#
#    Iu = np.array([ np.sin(P),
#     -np.sin(R) * np.cos(P),
#     -np.cos(R) * np.cos(P) ])

    C_rpy2enu = np.hstack((Ir[np.newaxis].T,Ip[np.newaxis].T,Iy[np.newaxis].T))

    return C_rpy2enu


def C_rpy2enu2(roll,pitch,yaw,angtype='deg'):

    if angtype == 'deg':
        roll = np.deg2rad(roll)
        pitch = np.deg2rad(pitch)
        yaw = np.deg2rad(yaw)

    R = roll
    P = pitch
    Y = yaw

    Sr = np.sin(R)
    Cr = np.cos(R)
    Sp = np.sin(P)
    Cp = np.cos(P)
    Sy = np.sin(Y)
    Cy = np.cos(Y)

    C_rpy2enu = np.array([[Sy*Cp,  Cr*Cy + Sr*Sy*Sp , - Sr*Cy + Cr*Sy*Sp],
                         [Cy*Cp , -Cr*Sy + Sr*Cy*Sp ,   Sr*Sy + Cr*Cy*Sp],
                         [Sp    , -Sr*Cp , -Cr * Cp]])

    return C_rpy2enu


def C_eci2rtn(P,V):
    """
    gives the transformation matrix between ECI ref. frame and RTN (aka RIC or RSW)

    Parameters
    ----------
    P : numpy.array
        3D vector, position of the object in ECI frame

    V : numpy.array
        3D vector, velocity of the object in ECI frame
        
    Returns
    -------
    C_eci2rtn : numpy.array
        transformation matrix
            
    Source
    ------
    "Coordinate Systems", ASEN 3200 1/24/06 George H. Born
    https://www.colorado.edu/ASEN/asen3200/handouts/Coordinate%20System.pdf
    """


    H = np.cross(P,V)

    R = P / np.linalg.norm(P)
    C = H / np.linalg.norm(H)
    I = np.cross(C,R)

    C_eci2rtn = np.column_stack((R,I,C))

    return C_eci2rtn


def C_rtn2rpy():
    C_rtn2rpy = np.array([[ 0.,  0  ,  -1.],
                          [ 1.,  0. ,   0.],
                          [ 0,  -1. ,   0.]])
    return C_rtn2rpy


def C_cep2itrs(xpole , ypole):
    """
    xpole and ypole are given in mas, miliarcsecond 
    
    Source
    ------
        Hofmann-Wellenhof,et al. 2008 - GNSS - p21
        Xu - 2007 - GPS - p17
        Xu - Orbits - p15
    """
    
    xpole2 = np.deg2rad(xpole * (1./3600.) * (10** -3))
    ypole2 = np.deg2rad(ypole * (1./3600.) * (10** -3))
    
    C_cep2itrs = np.array([[1.      , 0.     ,  xpole2],
                           [0.      , 1.     , -ypole2],
                           [-xpole2 , ypole2 , 1.     ]])
    
    return C_cep2itrs



def rot_quelconq(theta,Vx,Vy,Vz,angtype='deg'):

    N = np.linalg.norm([Vx,Vy,Vz])

    Vx = Vx / N
    Vy = Vy / N
    Vz = Vz / N

    if angtype == 'deg':
        theta = np.deg2rad(theta)

    c = np.cos(theta)
    s = np.sin(theta)

    R = np.array([[Vx ** 2 + (1-Vx**2)*c , Vx * Vy * (1 - c) - Vz * s , Vx * Vz * (1-c) +Vy * s],
                 [Vx * Vy * (1-c) + Vz * s , Vy**2 + (1-Vy**2)*c , Vy*Vz * (1-c) - Vx * s ],
                 [Vx * Vz * (1-c) - Vy * s , Vy * Vz * ( 1-c) +Vx * s , Vz**2 + (1-Vz**2) * c]])

    return R



def vector_RPY(V,out_deg=True,ad_hoc_mode=False):
    """
    Donne le "Roll Pitch Yaw" pour un Point/Vecteur
    (les coordonnées spheriques plus exactement, le roll sera tjrs nul)

    le ad_hoc_mode est la pour corriger un bug que je n'explique pas
    pour l'utilisation de ads_offset
    (cette methode a été trouvé par les scripts ENURPY2 et ENURPY3)
    """

    x = V[0]
    y = V[1]
    z = V[2]

    n = np.linalg.norm(V)

    roll  = 0
    pitch = np.arcsin(z/n)
    yaw   = np.arctan2(y,x)

    if out_deg:
        roll  = np.rad2deg(roll)
        pitch = np.rad2deg(pitch)
        yaw   = np.rad2deg(yaw)

    if ad_hoc_mode:
        roll  = -roll
        pitch = -pitch
        yaw   = -yaw

        roll , pitch , yaw = pitch , roll , yaw

    return roll ,  pitch , yaw

def add_offset(Direction , Delta , Point = None ,
               out_delta_enu = False ):
    """
    Un mode adhoc a du être implémenté dans la fonction élémentaire
    vector_RPY pour que ça puisse marcher correctement, reste à
    comprendre pourquoi ...

    out_delta_enu est prioritaire sur  les autres modes de sortie
    et renvoie le  delta

    par défaut on a le delta + le vecteur Direction
    et si Point est précisé : Point + Direction
    """

    rpy = vector_RPY(Direction,ad_hoc_mode=True)

    Crpy2enu = C_rpy2enu(*rpy)
    Cenu2rpy = np.linalg.inv(C_rpy2enu(*rpy))

    Delta_enu = Crpy2enu.dot(C_enu2ned().dot(Delta))

    if out_delta_enu:
        print(out_delta_enu)

    if not np.isclose(np.linalg.norm(Delta) , np.linalg.norm(Delta_enu)):
        print("WARN : np.linalg.norm(Delta) != np.linalg.norm(Delta_enu)")
        print(np.linalg.norm(Delta) , np.linalg.norm(Delta_enu) , \
        np.linalg.norm(Delta) - np.linalg.norm(Delta_enu))


    if not type(Point) is None:
        outpoint = Point + Delta_enu
    else:
        outpoint = Direction + Delta_enu

    return outpoint

#### Coordinates conversion

def wnorm(phi,a=6378137.,e2=0.00669438003):
    """
    D'après PYACS - JMN
    
    ALG0021 in NTG_71 (IGN Lambert Projection Tech document)
    """
    from numpy import sqrt,sin

    e=sqrt(e2)

    wd=sqrt(1 - e2*sin(phi)**2)
    result=a/wd
    return result

def normal_vector(phi,llambda,angle='deg',normalized=True):
    #[Bosser 2012] p27

    if angle == 'deg':
        llambda = np.deg2rad(llambda)
        phi = np.deg2rad(phi)

    a = np.cos(llambda) * np.cos(phi)
    b = np.sin(llambda) * np.cos(phi)
    c = np.sin(phi)

    n = np.array([a,b,c])

    if normalized:
        n = n / np.linalg.norm(n)

    return n


def GEO2XYZ(phi,llambda,he,angle='deg'):
    # D'après PYACS - JMN
    from numpy import sqrt,sin,cos

    if angle == 'deg':
        llambda = np.deg2rad(llambda)
        phi = np.deg2rad(phi)

    a=6378137.
    e2=0.00669438003
    f=1.0 - sqrt(1-e2)

    wn=wnorm(phi)

    xx=(wn+he)*cos(phi)*cos(llambda)
    yy=(wn+he)*cos(phi)*sin(llambda)
    zz=(wn*(1.0-e2)+he)*sin(phi)

    return xx,yy,zz


def XYZ2GEO(x,y,z,outdeg=True):
    # D'après PYACS - JMN
    """
    Coordinates conversion

    XYZ Cartesian => FLH Geographic

    Args :
        x,y,z : cartesian coordinates

        outdeg : if True, degrees output for the angle, else radian

    Returns :
        latitude, longitude, height
    """

    from numpy import sqrt,arctan2,sin,cos,rad2deg

    A=6378137.
    E2=0.00669438003
    F=1.0 - sqrt(1-E2)

    TP = sqrt (x**2 + y**2)
    R = sqrt (TP**2 + z**2)

    TMU = arctan2( (z/TP)*((1. - F) + (E2*A)/R) ,1)
    RLBDA = arctan2(y,x)

    S3 = (sin(TMU))**3
    C3 = (cos(TMU))**3
    T1 = z*(1 - F) + E2*A*S3
    T2 = (1 - F)*( TP - E2*A*C3 )

    RPHI = arctan2(T1,T2)

    RHE = TP*(cos(RPHI)) + z*(sin(RPHI))
    RHE = RHE - A*( sqrt(1-E2*(sin(RPHI)**2)) )

    if outdeg:
        RPHI = rad2deg(RPHI)
        RLBDA = rad2deg(RLBDA)

    return RPHI,RLBDA,RHE


def XYZ2ENU(dX,dY,dZ,lat0,lon0):
    '''
    core function for XYZ2ENU_2

    dXYZ = XYZrover - XYZref
    '''
    f0 = np.deg2rad(lat0)
    l0 = np.deg2rad(lon0)

    R = np.array([[ -np.sin(l0) , np.cos(l0) , 0 ] ,
    [ -np.sin(f0)*np.cos(l0) , -np.sin(f0)*np.sin(l0) , np.cos(f0) ] ,
    [  np.cos(f0)*np.cos(l0) , np.cos(f0)*np.sin(l0)  , np.sin(f0)]])

    enu = np.dot(R,np.vstack((dX,dY,dZ)))

    E = enu[0,:]
    N = enu[1,:]
    U = enu[2,:]

    return E,N,U

def XYZ2ENU_2(X,Y,Z,x0,y0,z0):
    '''
    Coordinates conversion

    XYZ ECEF Geocentric => ENU Topocentic

    Args :
        X,Y,Z    : geocentric coordinates to be converted (may be numpy array)

        x0,y0,z0 : coordinate of the topocentric origin point
        in the geocentric frame
    Returns :
        E,N,U    :  topocentric coordinates
    '''
    f0,l0,h0 = XYZ2GEO(x0,y0,z0,outdeg=1)
    dX = np.array(X) - x0
    dY = np.array(Y) - y0
    dZ = np.array(Z) - z0
    E,N,U = XYZ2ENU(dX,dY,dZ,f0,l0)
    return E,N,U


def XYZ2ENU_around_fix_pos(X,Y,Z):
    """
    calc the mean of the X,Y,Z and
    return the ENU pos to this mean pos

    Returns
        Earr , Narr , Uarr ,
        [x0 , y0 , z0] :

    """
    x0 = np.nanmean(X)
    y0 = np.nanmean(Y)
    z0 = np.nanmean(Z)
    return XYZ2ENU_2(X,Y,Z,x0,y0,z0) , np.array([x0,y0,z0])


def ENU2XYZ(E,N,U,Xref,Yref,Zref):
    if genefun.is_iterable(E):
        Xlist , Ylist , Zlist = [] , [] , []
        for e,n,u in zip(E,N,U):
            x,y,z = ENU2XYZ(e,n,u,Xref,Yref,Zref)
            Xlist.append(x)
            Ylist.append(y)
            Zlist.append(z)
        return np.array(Xlist) , np.array(Ylist) , np.array(Zlist)
    
    else:
        fr,lr,hr = XYZ2GEO(Xref,Yref,Zref)
        f0 = np.deg2rad(fr)
        l0 = np.deg2rad(lr)
    
        R = np.array([[ -np.sin(l0)            , np.cos(l0)             , 0          ],
                      [ -np.sin(f0)*np.cos(l0) , -np.sin(f0)*np.sin(l0) , np.cos(f0) ],
                      [  np.cos(f0)*np.cos(l0) , np.cos(f0)*np.sin(l0)  , np.sin(f0)]])
    
        R3 = R.T
        R3 = scipy.linalg.inv(R)
    
        ENU = np.vstack((E,N,U))
    
        xyz = np.dot(R3,ENU) #+ np.vstack((Xref,Yref,Zref))
    
        X = float(xyz[0]) + Xref
        Y = float(xyz[1]) + Yref
        Z = float(xyz[2]) + Zref
    
        return X,Y,Z


def ENU2XYZ_legacy(E,N,U,Xref,Yref,Zref):
    """
    KEPT FOR LEGACY REASONS, use ENU2XYZ

    differe de ENU2XYZ pour d'obscure raisons, à investiguer !!!
    est laissé pour des scripts de conversion de GINS (170119)
    """

    fr,lr,hr = XYZ2GEO(Xref,Yref,Zref)
    f0 = np.deg2rad(fr)
    l0 = np.deg2rad(lr)

    R = np.array([[ -np.sin(l0) , np.cos(l0) , 0 ] ,
    [ -np.sin(f0)*np.cos(l0) , -np.sin(f0)*np.sin(l0) , np.cos(f0) ] ,
    [  np.cos(f0)*np.cos(l0) , np.cos(f0)*np.sin(l0) , np.sin(f0)]])

    R3 = R.T

    ENU = np.vstack((E,N,U))

    xyz = np.dot(R3,ENU) #+ np.vstack((Xref,Yref,Zref))

    X = float(xyz[0])
    Y = float(xyz[1])
    Z = float(xyz[2])

    return X,Y,Z


def sFLH2sXYZ(F,L,H,sF,sL,sH,ang='deg'):
    # conversion batarde du sigma FLH => sigma XYZ
    if ang == 'deg':
        F  = np.deg2rad(F)
        L  = np.deg2rad(L)
        sF = np.deg2rad(sF)
        sL = np.deg2rad(sL)
    X,Y,Z  = GEO2XYZ(F,L,H)
    X2,Y2,Z2 = GEO2XYZ(F+sF,L+sL,H+sH)

    return np.abs(X-X2) , np.abs(Y-Y2) , np.abs(Z-Z2)

#def sXYZ2sFLH(X,Y,Z,sX,sY,sZ):
#    F,L,H    = XYZ2GEO(X,Y,Z)
#    F2,L2,H2 = XYZ2GEO(X+sX,Y+sY,Z+sZ)
#
#    return np.abs(F-F2) , np.abs(L-L2) , np.abs(H-H2)
#

def sFLH2sENU(F,L,H,sF,sL,sH,ang='deg'):
    # conversion batarde du sigma FLH => sigma ENU
    # Par conversion des angles en distance
    if ang == 'deg':
        F = np.deg2rad(F)
        L = np.deg2rad(L)
        sF = np.deg2rad(sF)
        sL = np.deg2rad(sL)

    E2=0.00669438003
    f=1.0 - np.sqrt(1-E2)
    A=6378137.

    # Je présume que Rpsi est le rayon du cercle tangent à l'ellipsoide à une
    # lattitude donnée (comm du 150710)
    Psi = np.arctan((1 - E2) * np.tan(F))
    Rpsi = A * np.sqrt(1 - E2) / np.sqrt(1 - E2 * np.cos(Psi)**2)

    # r est le rayon d'un petit cercle (un parallèle)
    r = np.cos(Psi) * (Rpsi + H)

    # On pourait simplifier par 2pi mais autant avoir toute la démarche
    sE = (np.pi * 2 * r    * sL) / (2 * np.pi)
    sN = (np.pi * 2 * Rpsi * sF) / (2 * np.pi)
    sU = sH
    return sE,sN,sU

def sENU2sFLH(F,L,H,sE,sN,sU,ang='deg'):
    # conversion batarde du sigma ENU => sigma FLH
    # Par conversion des angles en distance

    E2=0.00669438003
    f=1.0 - np.sqrt(1-E2)
    A=6378137.

    # Je présume que Rpsi est le rayon du cercle tangent à l'ellipsoide à une
    # lattitude donnée (comm du 150710)
    Psi = np.arctan((1 - E2) * np.tan(F))
    Rpsi = A * np.sqrt(1 - E2) / np.sqrt(1 - E2 * np.cos(Psi)**2)

    # r est le rayon d'un petit cercle (un parallèle)
    r = np.cos(Psi) * (Rpsi + H)

    # On pourait simplifier par 2pi mais autant avoir toute la démarche
    sL = sE *  (2 * np.pi) / (np.pi * 2 * r    )
    sF = sN *  (2 * np.pi) / (np.pi * 2 * Rpsi )
    sH = sU

    if ang=='deg':
        sL = np.rad2deg(sL)
        sF = np.rad2deg(sF)
    return sF,sL,sH


def sXYZ2sENU(X,Y,Z,sX,sY,sZ,sXY=0,sYZ=0,sXZ=0):
    """
    conversion "batarde" en assumant par défaut
    que XYZ ne sont pas corrélé (pas bien ...)
    Linear Algebra, Geodesy, and GPS p332
    """

    SIGMAxyz = np.array([[sX**2,sXY,sXZ],[sXY,sY**2,sYZ],[sXZ,sYZ,sZ**2]])

    F,L,H = XYZ2GEO(X,Y,Z)

    C = C_ecef2enu(F,L,angtype='deg')

    SIGMAenu = np.dot(np.dot(C,SIGMAxyz),C.T)

    sE = np.sqrt(SIGMAenu[0,0])
    sN = np.sqrt(SIGMAenu[1,1])
    sU = np.sqrt(SIGMAenu[2,2])

    return sE,sN,sU


def ECI2RTN_or_RPY(P,V,C,out_rpy=False,rpy_theo_mode=False):
    """
    convert ECI coordinates in RTN (RIC) or RPY (Roll Pitch Yaw)

    Parameters
    ----------
    P : numpy.array
        3D vector, position of the ref object in ECI frame
    V : numpy.array
        3D vector, velocity of the ref object in ECI frame
    C : numpy.array
        3D Vector, coordinates in ECF frame that will be transformed
    out_rpy : bool
        if True output in RPY frame, RTN instead
    rpy_theo_mode : bool
        use the theoretical matrix composition, but wrong ans.
        empirically, only for debug !!

    Returns
    -------
    Cout : conversion of P in RTN or RPY ref. frame

    Source
    ------
    "Coordinate Systems", ASEN 3200 1/24/06 George H. Born
    """

    C_eci2rtn_mat  = C_eci2rtn(P,V)

    if not out_rpy:
        TransMat = C_eci2rtn_mat
    else:

        if not rpy_theo_mode:
            # Pour de très obscures raisons la composition est inversée
            # par rapport à l'ordre standard ... (241017)
            TransMat = np.dot(C_rtn2rpy().T,C_eci2rtn_mat.T)
        else:
            print("WARN : using the theoretical mode for RPY conversion, UNSTABLE & WRONG !")
            TransMat = np.dot(C_rtn2rpy(),C_eci2rtn_mat)

        # Mais reste compatible avec Wikipedia
        # https://en.wikipedia.org/wiki/Permutation_matrix#Permutation_of_rows_and_columns

    Cout = TransMat.dot( C )

    # EXEMPLE POUR DEBUG
    #     AAAAAAAAAAAAA yawerr,aux,nomi
    #       -1.606542300481316
    #       -1.606542300481316
    # CCCCCCCCCCCC sunors orientation du soleil
    #   2.379858191647978E-002   6.654859888132475E-001  -7.460308480163662E-001
    # DDDDDDDDDDDDD Inputs : xsat , vsat , xsun
    #   -15779.437215560150435   -23508.964967410494864     8644.453520237560952
    #        1.296220384065493       -1.917314252639830       -2.847293074880089
    #  -1.342404641197934E+008  -5.934080206107079E+007  -2.572342068927875E+007
    # EEEEEEEEEEEE Rcrs2ors
    #   3.532662370257567E-001  -7.688096793507927E-001   5.330195492501879E-001
    #  -5.225364448455557E-001   3.104619621186746E-001   7.941181766553120E-001
    #  -7.759888076421138E-001  -5.590572803191509E-001  -2.920042493231377E-001
    #    Rcrs2ors = np.array([[ 0.35330134, -0.76880968,  0.53301955],
    #                         [-0.52248415,  0.31046196,  0.79411818],
    #                         [-0.77600804, -0.55905728, -0.29200425]])

    #    Psat = np.array([   -15779.437215560150435 ,  -23508.964967410494864  ,   8644.453520237560952])
    #    Vsat = np.array([        1.296220384065493  ,     -1.917314252639830  ,     -2.847293074880089])
    #    Psun = np.array([  -1.342404641197934E+008 , -5.934080206107079E+007 , -2.572342068927875E+007])

    return Cout


def ECI2RTN(P,V,C):
    """
    legacy wrapper of ECI2RTN_or_RPY
    """
    return ECI2RTN_or_RPY(P,V,C,out_rpy=False)


def ECEF2ECI(xyz,utc_times):
    """ 
    Convert ECEF (Earth Centered Earth Fixed) positions to ECI (Earth Centered Inertial)
    positions:
        XYZ are cartesian positions in ECEF. Should have shape (N,3)
        UTC_times are UTC_times, as datetime objects. Sould have shape (N)


     [X]    [C -S 0][X]
     [Y]  = [S  C 0][Y]
     [Z]eci [0  0 1][Z]ecf

     C and S are cos() and sin() of gmst (Greenwich Meridian Sideral Time)

    Source
    ------
    http://ccar.colorado.edu/ASEN5070/handouts/coordsys.doc
    Inspired from satellite-js (https://github.com/shashwatak/satellite-js)
    """
    from pyorbital import astronomy
    # XYZ and utc_time must have the same shape
    #if not xyz.shape[:-1] == utc_times.shape:
    #    raise ValueError("shape mismatch for XYZ and utc_times (got {} and {})".format(xyz.shape[:-1],utc_times.shape))

    #    gmst = -1 * astronomy.gmst(utc_times) # EDIT 180430 : Why this -1 ??? removed because wrong ! ...
    gmst = 1 * astronomy.gmst(utc_times)

    eci = xyz.copy()
    eci[:,0] = xyz[:,0]*np.cos(gmst) - xyz[:,1]*np.sin(gmst)
    eci[:,1] = xyz[:,0]*np.sin(gmst) + xyz[:,1]*np.cos(gmst)
    return eci

    from pyorbital import astronomy




### Angle conversion
    
def arcsec2deg(arcsec_in):
    return arcsec_in / 3600.

def deg2arcsec(deg_in):
    return deg_in * 3600.


def dms2dec_num(deg,minn=0,sec=0):
    """
    Convert :
    Degree Minute Second `float` Angle => decimal Degree `float` Angle
        
    Parameters
    ----------
    deg & minn & sec : float
        degres, minutes and seconds of the input Angle

    Returns
    -------
    dd_float : float
        Decimal degree Angle
    """
    return deg + minn * (1./60.) +  sec * (1./3600.)


def dms2dec(dms_str , onlyDM=False):
    """   
    Convert :
    Degree Minute Second `string` Angle => decimal Degree `float` Angle
    
    can manage only DM angle in input  
    
    Too Complicated .... must be simplified
    
    Parameters
    ----------
    dms_str : str
        string of DMS angle
        e.g.
        "2°20'35.09"E"

    onlyDM : bool
        True if the string is only a Degree Minute Angle 
        e.g.
        "40°52.0931'N"
        "28°31.4136'E"
                
    Returns
    -------
    dd_float : float
        Decimal degree Angle
    """
    
    if not onlyDM:
        print("WARN : DMS mode not well implemented yet !!! ")
    
    
    dms_str = dms_str.strip()

    if re.match('[swoSWO]', dms_str):
        sign = -1
    else:
        sign = 1

    # former regex
    lis = re.findall(r'\d+|\D+', dms_str)

    ipt = lis.index('.')
    decimal = lis[ipt-1] +  lis[ipt] + lis[ipt+1]
    lis[ipt-1] = decimal
    lis = lis[:ipt]

    for i,e in enumerate(lis):
        try:
            lis[i] = float(e)
        except:
            lis.remove(e)

    lis = [float(e) for e in lis]
    print(lis)

    deg = lis[0]
    minu = lis[1]
    if onlyDM:
        sec = 0.
    else:
        try:
            sec  = lis[2]
        except IndexError:
            print("ERR : did you forgot to activate the DM only mode ? ;) ?")
            return None

    dd_float = sign * (deg + minu / 60. + sec / 3600.)
    return dd_float

def deg_dec2dms(deg_in):
    deg              = np.floor(deg_in)
    decimal_part     = deg_in - deg 
    decimal_part_sec = decimal_part * 3600
    minu             = np.floor_divide(decimal_part_sec,60)
    sec              = decimal_part_sec - minu * 60
    sec              = np.round(sec,8)
    
    return deg , minu , sec

def angle2equivalent_earth_radius(angle_in,angtype='deg'):
    earth_radius = 6371008.8
    earth_circum = earth_radius * 2 * np.pi
    if angtype == "deg":
        equiv_out = (angle_in * earth_circum) / 360.
    elif angtype == "rad":
        equiv_out = (angle_in * earth_circum) / (np.pi *2)
    elif angtype == "mas":
        equiv_out = (angle_in * 10**-3 * earth_circum) / (86400.)  
        
    return equiv_out
        

def anglesfromvects(xa,ya,xb,yb,angtype='deg'):

    A = np.array([xa,ya])
    B = np.array([xb,yb])

    ps = np.inner(A,B)
    a = np.arccos( ps / (np.linalg.norm(A) * np.linalg.norm(B)) )

    if angtype == 'deg':
        return np.rad2deg(a)
    elif angtype == 'rad':
        return a
    else:
        print('ERR : angfromvects : mauvais angtype')

def angle_interpolation_quick(A,B,w):
    """
    Determine the interpolation between angle A & B
    by conversion to the cartesian space
    the parameter w € [0,1] define the interpoled angle C(w)
    where C(w=0) = A  &  C(w=1) = B

    Source
    ------
    https://stackoverflow.com/questions/2708476/rotation-interpolation
    """
    CS = (1-w)*np.cos(A) + w*np.cos(B)
    SN = (1-w)*np.sin(A) + w*np.sin(B)
    C = np.atan2(SN,CS)

    return C

def angle_from_3_pts(p1,p2,p3):
    """
    Gives angle between 3 points
    p3 is the vertex (sommet)
    
    Source
    ------
    http://www.les-mathematiques.net/phorum/read.php?8,596072,596231
    """

    p1,p2,p3 = [np.array(e) for e in (p1,p2,p3)]
    x1 , y1 = p1
    x2 , y2 = p2
    x3 , y3 = p3

    kos = ((x1 - x3)*(x2 - x3) + (y1 - y3)*(y2 - y3)) / (np.linalg.norm(p1 - p3) * np.linalg.norm(p2 - p3))

    return np.arccos(kos)


def cartesian2polar(x,y):
    """
    cartesian => polar conversion
    Args:
        x , y (float or list of floats) : cartesian coordinates
    Returns:
        r , theta (float or list of floats) : polar coordinates
    """
    if genefun.is_iterable(x):
        x = np.array(x)
    if genefun.is_iterable(y):
        y = np.array(y)

    theta = np.arctan2(y,x)
    r     = np.sqrt(x** + y**2)
    return r , theta

def polar2cartesian(r,theta,ang='deg'):
    """
    polar => cartesian conversion
    Args:
        r , theta (float or list of floats) : polar coordinates

        ang (string) : 'deg' (degrees) or 'rad' (radian)

    Returns:
        x , y (float or list of floats) : cartesian coordinates

    """
    if genefun.is_iterable(r):
        r = np.array(r)
    if genefun.is_iterable(theta):
        theta = np.array(theta)

    if ang == 'deg':
        theta = np.deg2rad(theta)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x , y

### Low level geometric function

def dist(A,B):
    """
    Cartesian distance between 2 points A & B

    Args:
        A, B (vectors) : "points", can be 2D or 3D vectors (list, np.array ...)

    Returns :
        distance between A & B
    """
    A = np.array(A)
    B = np.array(B)
    return np.linalg.norm(A - B)

def dist_diff(A,B):
    """
    First derivative of cartesian distance between 2 points A & B

    Args:
        A, B (vectors) : "points", can be 2D or 3D vectors (list, np.array ...)

    Returns :
        diffA (vector) : [ dD/dxa , dD/dya , dD/dza ]

        diffB (vector) : [ dD/dxb , dD/dyb , dD/dzb ] = -diffA
    """

    dAB   = A-B
    dist  = scipy.linalg.norm(dAB)

    diffA =   dAB / dist
    diffB = - dAB / dist

    return diffA, diffB

def relative_orientation(x1,y1,x2,y2,trigo_orient = True):
    if trigo_orient:
        ang1 = np.mod(360 - np.rad2deg(np.arctan2((x2 - x1),(y2 - y1))),360)
    else:
        ang1 = np.mod(np.rad2deg(np.arctan2((x2 - x1),(y2 - y1))),360)
    return ang1

def barycenter(points_list_in):
    points_arr = np.array(points_list_in)
    return np.mean(points_arr[:,-3:], axis=0)

def pythagore(a,b,c=0):
    return np.sqrt(a**2 + b**2 + c**2)

def equilateral_triangle(side):
    """
    Gives coordinates of an equilateral triangle of a given side
    """
    hauteur_len = np.sqrt(side**2 - (side * .5) **2)
    A = np.array([0,(2./3.) * hauteur_len])
    B = np.array([side / 2. , - (1./3.) * hauteur_len])
    C = np.array([- side / 2. , - (1./3.) * hauteur_len])

    return A,B,C


def vincenty_full(point1, point2, miles=False,full=True,azimuth_in_deg=True):
    """
    gives the Azimuth between 2 points
    
    Vincenty's formula (inverse method) to calculate the distance (in
    kilometers or miles) between two points on the surface of a spheroid
    Doctests:
    >>> vincenty((0.0, 0.0), (0.0, 0.0))  # coincident points
    0.0
    >>> vincenty((0.0, 0.0), (0.0, 1.0))
    111.319491
    >>> vincenty((0.0, 0.0), (1.0, 0.0))
    110.574389
    >>> vincenty((0.0, 0.0), (0.5, 179.5))  # slow convergence
    19936.288579
    >>> vincenty((0.0, 0.0), (0.5, 179.7))  # failure to converge
    >>> boston = (42.3541165, -71.0693514)
    >>> newyork = (40.7791472, -73.9680804)
    >>> vincenty(boston, newyork)
    298.396057
    >>> vincenty(boston, newyork, miles=True)
    185.414657
    
    https://github.com/maurycyp/vincenty/blob/master/vincenty/
    """
    
    # WGS 84
    a = 6378137  # meters
    f = 1 / 298.257223563
    b = 6356752.314245  # meters; b = (1 - f)a
    
    MILES_PER_KILOMETER = 0.621371
    
    MAX_ITERATIONS = 200
    CONVERGENCE_THRESHOLD = 1e-12  # .000,000,000,001

    # short-circuit coincident points
    if point1[0] == point2[0] and point1[1] == point2[1]:
        return 0.0

    import math 
    
    U1 = math.atan((1 - f) * math.tan(math.radians(point1[0])))
    U2 = math.atan((1 - f) * math.tan(math.radians(point2[0])))
    L = math.radians(point2[1] - point1[1])
    Lambda = L

    sinU1 = math.sin(U1)
    cosU1 = math.cos(U1)
    sinU2 = math.sin(U2)
    cosU2 = math.cos(U2)
    sinL  = math.sin(L)
    cosL  = math.cos(L)

    for iteration in range(MAX_ITERATIONS):
        sinLambda = math.sin(Lambda)
        cosLambda = math.cos(Lambda)
        sinSigma = math.sqrt((cosU2 * sinLambda) ** 2 +
                             (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) ** 2)
        if sinSigma == 0:
            return 0.0  # coincident points
        cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda
        sigma = math.atan2(sinSigma, cosSigma)
        sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma
        cosSqAlpha = 1 - sinAlpha ** 2
        try:
            cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha
        except ZeroDivisionError:
            cos2SigmaM = 0
        C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha))
        LambdaPrev = Lambda
        Lambda = L + (1 - C) * f * sinAlpha * (sigma + C * sinSigma *
                                               (cos2SigmaM + C * cosSigma *
                                                (-1 + 2 * cos2SigmaM ** 2)))
        if abs(Lambda - LambdaPrev) < CONVERGENCE_THRESHOLD:
            break  # successful convergence
    else:
        return None  # failure to converge

    uSq = cosSqAlpha * (a ** 2 - b ** 2) / (b ** 2)
    A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
    B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))
    deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 * (cosSigma *
                 (-1 + 2 * cos2SigmaM ** 2) - B / 6 * cos2SigmaM *
                 (-3 + 4 * sinSigma ** 2) * (-3 + 4 * cos2SigmaM ** 2)))
    s = b * A * (sigma - deltaSigma)

    s /= 1000  # meters to kilometers
    if miles:
        s *= MILES_PER_KILOMETER  # kilometers to miles


    if not full:
        return round(s, 6)
    else:        
        fwdAz = np.arctan2(cosU2*sinL,  cosU1*sinU2-sinU1*cosU2*cosL)
        revAz = np.arctan2(cosU1*sinL, -sinU1*cosU2+cosU1*sinU2*cosL)
       
        if azimuth_in_deg:
            fwdAz = np.rad2deg(fwdAz)
            revAz = np.rad2deg(revAz)
       
        return round(s, 6),fwdAz,revAz
    
def orthogonal_projection(Xa,Xb,Xv):
    """
    General description

    Parameters
    ----------
    param1 : float or int or str or dict or n-tuple or bool or list or numpy.array
        Description param1

    param2 : float or int or str or dict or n-tuple or bool or list or numpy.array
        Description param2
        
    param3 : float or int or str or dict or n-tuple or bool or list or numpy.array
        Description param3
                
    Returns
    -------
    out1 : float or int or str or dict or n-tuple or bool or list or numpy.array
        Description out1
    
    out2 : float or int or str or dict or n-tuple or bool or list or numpy.array
        Description out2
        
    Note
    ----
    Misc. Notes

    Source
    ------
    https://fr.wikipedia.org/wiki/Projection_orthogonale
    
    Examples
    --------
    >>> answer
    42 
    """
    xa , ya = Xa
    xb , yb = Xb    
    xv , yv = Xv
    
    D =np.sqrt(xv**2 + yv **2) 
    BH = ((xa - xb)*xv + (ya - yb)*yv) / D
    
    xh = xb + (BH / D) * xv
    yh = yb + (BH / D) * yv
    
    Xh = np.array([xh,yh])
    return Xh , dist(Xa,Xh)


def line_maker(X1,Y1,X2,Y2,nbpts=10000):
    X = np.linspace(X1,X2,nbpts)
    Y = np.linspace(Y1,Y2,nbpts)
    return X,Y
    

#  _____           _           _   _                _____           _        
# |  __ \         (_)         | | (_)              / ____|         | |       
# | |__) | __ ___  _  ___  ___| |_ _  ___  _ __   | |     __ _ _ __| |_ ___  
# |  ___/ '__/ _ \| |/ _ \/ __| __| |/ _ \| '_ \  | |    / _` | '__| __/ _ \ 
# | |   | | | (_) | |  __/ (__| |_| | (_) | | | | | |___| (_| | |  | || (_) |
# |_|   |_|  \___/| |\___|\___|\__|_|\___/|_| |_|  \_____\__,_|_|   \__\___/ 
#                _/ |                                                        
#               |__/        

def latitude_isometric(phi,e):
    """
    ALG0001 in NTG_71 (IGN Lambert Projection Tech document)
    """
    A = np.tan( (np.pi/4) + (phi/2) )
    B1 = 1 - e * np.sin(phi)
    B2 = 1 + e * np.sin(phi)
    
    L = np.log(A * ((B1/B2)**(.5 * e)))
    
    return L

def lambert_projection(long,lat,n,c,e,lc,Xs,Ys):
    """
    ALG0003 in NTG_71 (IGN Lambert Projection Tech document)
    """
    L = latitude_isometric(lat,e)
    X = Xs + c*np.exp(-n * L) * np.sin(n * (long - lc))
    Y = Ys - c*np.exp(-n * L) * np.cos(n * (long - lc))
    return X,Y

def lambert_secante_parameter(a , e , l0 , phi0 , phi1 , phi2 , X0 , Y0):
    """
    ALG0054 in NTG_71 (IGN Lambert Projection Tech document)
    """
    lc = l0
    A = np.log((wnorm(phi2,a,e**2) * np.cos(phi2)) / (wnorm(phi1,a,e**2) * np.cos(phi1)))
    B = latitude_isometric(phi1,e) - latitude_isometric(phi2,e)
    
    n = A/B
    
    C = ((wnorm(phi1,a,e**2) * np.cos(phi1))/n) * np.exp(n * latitude_isometric(phi1,e))
        
    if np.isclose(phi0,np.pi/2):
        Xs = X0
        Ys = Y0  
    else:
        Xs = X0
        Ys = Y0 + C * np.exp(n * latitude_isometric(phi0,e))
            
    return e , n , C , lc , Xs , Ys

def lambert_projection_CC_frontend(long,lat,NZ = 93):
    """
    General description

    Parameters
    ----------
    param1 : float or int or str or dict or n-tuple or bool or list or numpy.array
        Description param1

    param2 : float or int or str or dict or n-tuple or bool or list or numpy.array
        Description param2
        
    NZ : int
        Lambert93, NZ = 0 or = 93
                
    Returns
    -------
    out1 : float or int or str or dict or n-tuple or bool or list or numpy.array
        Description out1
    
    out2 : float or int or str or dict or n-tuple or bool or list or numpy.array
        Description out2
        
    Note
    ----
    Misc. Notes

    Source
    ------
    https://fr.wikipedia.org/wiki/Projection_conique_conforme_de_Lambert   
    geodesie.ign.fr/contenu/fichiers/documentation/algorithmes/notice/NTG_71.pdf
    https://geodesie.ign.fr/contenu/fichiers/documentation/rgf93/Lambert-93.pdf
 
    """
    if NZ in (0,93):
        phi0 = np.deg2rad(46.5)
        phi1 = np.deg2rad(44)
        phi2 = np.deg2rad(49)
        l0 = np.deg2rad(3)
        X0 = 700000.
        Y0 = 6600000.
        a = 6378137.
        e = 0.0818191910435
    else:        
        phi0 = np.deg2rad(41 + NZ)
        phi1 = np.deg2rad(phi0 - 0.75)
        phi2 = np.deg2rad(phi0 + 0.75)
        l0 = np.deg2rad(3)
        X0 = 1700000.
        Y0 = (NZ * 1000000.) + 200000.
        a = 6378137.
        e = 0.0818191910435
    
    e , n , c , lc , Xs , Ys = lambert_secante_parameter(a, e, l0, phi0, phi1, phi2, X0, Y0)
    X,Y = lambert_projection(np.deg2rad(long),np.deg2rad(lat),n,c,e,lc,Xs,Ys)
    
    return X,Y
    
    
    
#  _______ _                   _____                              _
# |__   __(_)                 / ____|                            (_)
#    | |   _ _ __ ___   ___  | |     ___  _ ____   _____ _ __ ___ _  ___  _ __
#    | |  | | '_ ` _ \ / _ \ | |    / _ \| '_ \ \ / / _ \ '__/ __| |/ _ \| '_ \
#    | |  | | | | | | |  __/ | |___| (_) | | | \ V /  __/ |  \__ \ | (_) | | | |
#    |_|  |_|_| |_| |_|\___|  \_____\___/|_| |_|\_/ \___|_|  |___/_|\___/|_| |_|
#
### Time conversion

def leap_seconds(f):
    """
    Return a list of tuples of this format: (timestamp, number_of_seconds)
        timestamp: a 32-bit timestamp, seconds since the UNIX epoch
        number_of_seconds: how many leap-seconds occur at timestamp

    INTERNAL_FUNCTION

    Source
    ------
    http://stackoverflow.com/questions/19332902/extract-historic-leap-seconds-from-tzdata

    """
    TZFILE_MAGIC = 'TZif'.encode('US-ASCII')

    fmt = ">4s c 15x 6l"
    size = struct.calcsize(fmt)
    (tzfile_magic, tzfile_format, ttisgmtcnt, ttisstdcnt, leapcnt, timecnt,
        typecnt, charcnt) =  struct.unpack(fmt, f.read(size))
    #print("DEBUG: tzfile_magic: {} tzfile_format: {} ttisgmtcnt: {} ttisstdcnt: {} leapcnt: {} timecnt: {} typecnt: {} charcnt: {}".format(tzfile_magic, tzfile_format, ttisgmtcnt, ttisstdcnt, leapcnt, timecnt, typecnt, charcnt))

    # Make sure it is a tzfile(5) file
    assert tzfile_magic == TZFILE_MAGIC, (
            "Not a tzfile; file magic was: '{}'".format(tzfile_magic))

    # comments below show struct codes such as "l" for 32-bit long integer
    offset = (timecnt*4  # transition times, each "l"
        + timecnt*1  # indices tying transition time to ttinfo values, each "B"
        + typecnt*6  # ttinfo structs, each stored as "lBB"
        + charcnt*1)  # timezone abbreviation chars, each "c"

    f.seek(offset, 1) # seek offset bytes from current position

    fmt = '>{}l'.format(leapcnt*2)
    #print("DEBUG: leapcnt: {}  fmt: '{}'".format(leapcnt, fmt))
    size = struct.calcsize(fmt)
    data = struct.unpack(fmt, f.read(size))

    lst = [(data[i], data[i+1]) for i in range(0, len(data), 2)]
    assert all(lst[i][0] < lst[i+1][0] for i in range(len(lst)-1))
    assert all(lst[i][1] == lst[i+1][1]-1 for i in range(len(lst)-1))
    return lst

def print_leaps(leap_lst):
    """
    INTERNAL_FUNCTION
    """
    # leap_lst is tuples: (timestamp, num_leap_seconds)
    outlist = []
    for ts, num_secs in leap_lst:
        dtime = (dt.datetime.utcfromtimestamp(ts - num_secs+1))
        outlist.append((dtime,num_secs))
    return outlist

def get_leapsecond_frontend():
    """
    INTERNAL_FUNCTION

    Nota :
        Le temps universel (UT1) et le Temps atomique international (TAI) ont
        été définis comme égaux en 1958. Lors de la mise en place d’UTC en
        1972, UT1 s’était décalé d’environ 10 secondes par rapport au TAI.
        On choisit donc un décalage initial de 10 secondes entre UTC et TAI .

        the initial 10sec are added in find_leapsecond
    """
    zoneinfo_fname = '/usr/share/zoneinfo/right/UTC'

    try:
        with open(zoneinfo_fname, 'rb') as f:
            leap_lst = leap_seconds(f)
            final_leap_lis = print_leaps(leap_lst)
    except:
        final_leap_lis = [(datetime.datetime(1972, 7, 1, 0, 0), 1),
         (datetime.datetime(1973, 1, 1, 0, 0), 2),
         (datetime.datetime(1974, 1, 1, 0, 0), 3),
         (datetime.datetime(1975, 1, 1, 0, 0), 4),
         (datetime.datetime(1976, 1, 1, 0, 0), 5),
         (datetime.datetime(1977, 1, 1, 0, 0), 6),
         (datetime.datetime(1978, 1, 1, 0, 0), 7),
         (datetime.datetime(1979, 1, 1, 0, 0), 8),
         (datetime.datetime(1980, 1, 1, 0, 0), 9),
         (datetime.datetime(1981, 7, 1, 0, 0), 10),
         (datetime.datetime(1982, 7, 1, 0, 0), 11),
         (datetime.datetime(1983, 7, 1, 0, 0), 12),
         (datetime.datetime(1985, 7, 1, 0, 0), 13),
         (datetime.datetime(1988, 1, 1, 0, 0), 14),
         (datetime.datetime(1990, 1, 1, 0, 0), 15),
         (datetime.datetime(1991, 1, 1, 0, 0), 16),
         (datetime.datetime(1992, 7, 1, 0, 0), 17),
         (datetime.datetime(1993, 7, 1, 0, 0), 18),
         (datetime.datetime(1994, 7, 1, 0, 0), 19),
         (datetime.datetime(1996, 1, 1, 0, 0), 20),
         (datetime.datetime(1997, 7, 1, 0, 0), 21),
         (datetime.datetime(1999, 1, 1, 0, 0), 22),
         (datetime.datetime(2006, 1, 1, 0, 0), 23),
         (datetime.datetime(2009, 1, 1, 0, 0), 24),
         (datetime.datetime(2012, 7, 1, 0, 0), 25),
         (datetime.datetime(2015, 7, 1, 0, 0), 26),
         (datetime.datetime(2017, 1, 1, 0, 0), 27)]
    return final_leap_lis

def find_leapsecond(dtin,get_leapsec_lis=[],
                    apply_initial_delta=True):
    """
    Find the TAI-UTC leap second for a given datetime

    Parameters
    ----------
    dtin : datetime
        Epoch for which the leap second is researched

    get_leapsec_lis : list, optional
        A list of leap second, provided by get_leapsecond_frontend() 
        automatically determined if given list is empty
        
    apply_initial_delta : bool, optional
        See note below
                 
    Returns
    -------
    leapsec_out : int
        The leap second for the given epoch
        
    Note
    ----
    Le temps universel (UT1) et le Temps atomique international (TAI) ont
    été définis comme égaux en 1958. Lors de la mise en place d’UTC en
    1972, UT1 s’était décalé d’environ 10 secondes par rapport au TAI.
    On choisit donc un décalage initial de 10 secondes entre UTC et TAI .
    
    
    In 1972, the leap-second system was introduced so that the broadcast UTC 
    seconds could be made exactly equal to the standard SI second, while still
    maintaining the UTC time of day and changes of UTC date synchronized with 
    those of UT1 (the solar time standard that superseded GMT).[11] By then, 
    the UTC clock was already 10 seconds behind TAI, which had been 
    synchronized with UT1 in 1958, but had been counting true SI seconds
    since then. After 1972, both clocks have been ticking in SI seconds, 
    so the difference between their readouts at any time is 10 seconds plus 
    the total number of leap seconds that have been applied to UTC 
    (37 seconds as of January 2017).
    """

    if get_leapsec_lis == []:
        get_leapsec_lis = get_leapsecond_frontend()

    get_leapsec_lis2 = [ e[0] for e in get_leapsec_lis ]

    get_leapsec_lis2.append(dtin)
    i_dtin = sorted(get_leapsec_lis2).index(dtin)

    if apply_initial_delta:
        leapsec_out = i_dtin + 10 
    else:
        leapsec_out = i_dtin
        
    return leapsec_out

def tgipsy2dt(tin):
    """
    Time conversion
    GIPSY Time => Datetime

    Parameters
    ----------
    tin : float or list  or numpy.array
        GIPSY time. Can handle several time float in a list.
                
    Returns
    -------
    dtout : datetime or list
        converted Datetime
        
    Note
    ----
    le temps GIPSY n'est pas le 'vrai' J2000 mais
    seulement compté a partir du 1er janvier 2000 midi
    """
    if genefun.is_iterable(tin):
        return [tgipsy2dt(e) for e in tin]
    else:
        j2000 = dt.datetime(2000, 1, 1, 12, 0)
        tout = j2000 + dt.timedelta(seconds=float(tin))
        return tout

    return tout

def matlab_time2dt(matlab_datenum):
    """
    Time conversion
    MATLAB Time => Datetime

    Parameters
    ----------
    matlab_datenum : float or list or numpy.array
        MATLAB time.  Can handle several time float in a list.
                
    Returns
    -------
    python_datetime : datetime of list
        converted Datetime 
    """
    if genefun.is_iterable(matlab_datenum):
        return [matlab_time2dt(e) for e in matlab_datenum]
    else:
        python_datetime = dt.datetime.fromordinal(int(matlab_datenum)) + \
        dt.timedelta(days=matlab_datenum%1) - dt.timedelta(days = 366)
    return python_datetime

def roundTime(dtin=None, roundTo=60):
    """Round a datetime object to any time laps in seconds
    dt : datetime.datetime object, default now.
    roundTo : Closest number of seconds to round to, default 1 minute.
    Author: Thierry Husson 2012 - Use it as you want but don't blame me.
    http://stackoverflow.com/questions/3463930/how-to-round-the-minute-of-a-datetime-object-python
    """
    import datetime as dtmod

    if dtin == None :
        dtin = dtmod.datetime.now()
    seconds = (dtin - dtin.min).seconds
    # // is a floor division, not a comment on following line:
    rounding = (seconds+roundTo/2) // roundTo * roundTo
    return dtin + dtmod.timedelta(0,rounding-seconds,-dtin.microsecond)


def dt_in_local_timezone2posix(dtin):
    if not genefun.is_iterable(dtin):
        return time.mktime(dtin.timetuple()) + dtin.microsecond * 0.000001
    else:
        return [ dt_in_local_timezone2posix(e) for e in dtin ]


def posix2dt_in_local_timezone(posixin):
    if not genefun.is_iterable(posixin):
        if np.isnan(posixin):
            return dt.datetime(1970,1,1)
        else:
            return dt.datetime.fromtimestamp(posixin)
    else:
        return [ posix2dt_in_local_timezone(e) for e in posixin ]

old_cnv_fcts = False
# old functions are bad !!!
# elles convertissent dans la timezone d'arrivée
if old_cnv_fcts:
    def dt2posix(dtin,out_array=False):
        if not genefun.is_iterable(dtin):
            return time.mktime(dtin.timetuple()) + dtin.microsecond * 0.000001
        else:
            L = [ dt2posix(e) for e in dtin ]
            if out_array:
                return np.array(L)
            else:
                return L

    def posix2dt(posixin,out_array=False):
        if not genefun.is_iterable(posixin):
            if np.isnan(posixin):
                return dt.datetime(1970,1,1)
            else:
                return dt.datetime.fromtimestamp(posixin)
        else:
            L = [ posix2dt(e) for e in posixin  ]
            if out_array:
                return np.array(L)
            else:
                return L

else:

    def dt2posix(dtin,out_array=False):       
        """
        Time conversion
        Python's Datetime => POSIX Time
    
        Parameters
        ----------
        dtin : datetime or list/numpy.array of datetime.
            Datetime object(s).  Can handle several datetimes in a list.
        out_array : bool
            if iterable as input, force the output as a Numpy array

        Returns
        -------
        L : float or list/numpy.array of floats
            POSIX Time  
        """
        if not genefun.is_iterable(dtin):
            D = dtin - dt.datetime(1970,1,1)
            return D.days * 86400 + D.seconds +  D.microseconds * 10**-6
        else:
            L = [ dt2posix(e) for e in dtin  ]
            if out_array:
                return np.array(L)
            else:
                return L

    def posix2dt(posixin,out_array=False):       
        """
        Time conversion
        POSIX Time => Python's datetime
    
        Parameters
        ----------
        posixin : float or list/numpy.array of floats.
            POSIX Time.  Can handle several time in a list.
        out_array : bool
            if iterable as input, force the output as a Numpy array

        Returns
        -------
        L : datetime or list/numpy.array of datetime.
            DateTime  
        """

        if not genefun.is_iterable(posixin):
            if np.isnan(posixin):
                return dt.datetime(1970,1,1)
            else:
                return dt.datetime(1970,1,1) + dt.timedelta(seconds=posixin)
        else:
            L = [ posix2dt(e) for e in posixin  ]
            if out_array:
                return np.array(L)
            else:
                return L


def datetime_improved(y=0,mo=0,d=0,h=0,mi=0,s=0,ms=0):
    """
    Improved time conversion of  Datetime

    can manage when you give list, array, or floats in input
    can manage NaN too (return POSIX 0 epoch if)

    Parameters
    ----------
    y,mo,d,h,mi,s,ms : float or list/numpy.array of floats.
        year , month... . Can handle several time in lists.

    Returns
    -------
    L : datetime or list/numpy.array of datetime.
        DateTime  
    """
    try:
        ms_from_s  = (s - np.floor(s)) * 10**6
        if ms_from_s != 0:
            if ms != 0:
                print('WARN : input sec contains microsecs, given microsec are ignored')
            ms = ms_from_s
        return dt.datetime(int(y),int(mo),int(d),int(h),int(mi),int(s),int(ms))
    except:
        return dt.datetime(1970,1,1) # si ca deconne, si on donne un NaN par ex
            
def dt2ymdhms(dtin,with_microsec = True):
    """
    Time conversion

    POSIX Time => Python's datetime

    Args :
        dtin  (datetime object or list of dt) : Python DateTime

        if with_microsec : if False rounding the microsecs to the nearest sec

    Returns :
        tuple with year , month , day , hour , minute , second , microsecond
    """
    if not genefun.is_iterable(dtin):
        if with_microsec:
            return (dtin.year,dtin.month,dtin.day,dtin.hour,dtin.minute,dtin.second,dtin.microsecond)
        else:
#            roundsec = int(np.round(dtin.second + dtin.microsecond * (10**-6)))
#            dt2 = dtin.replace(second=0)
#            dt2 = dt2.replace(microsecond=0)
#            dt2 = dt2 + dt.timedelta(seconds=roundsec)
            dt2 = roundTime(dtin,1)
            return (dt2.year,dt2.month,dt2.day,dt2.hour,dt2.minute,dt2.second)
    else:
        return [ dt2ymdhms(e,with_microsec) for e in dtin ]

def ymdhms_vectors2dt(yrlis,mlis,dlis,minlis,hlis,slis):
    """
    """
    dtlis = []
    for yr,m,d,minn,h,s in zip(yrlis,mlis,dlis,minlis,hlis,slis):
        dtlis.append(dt.datetime(int(yr),int(m),int(d),int(minn),int(h),int(s)))
    return np.array(dtlis)

# OBSOLETE
def dt2posix_list(Lin):
    return [ dt2posix(e) for e in Lin ]

def posix2dt_list(Lin):
    return [ posix2dt(e) for e in Lin ]

#
# Autore: Scott Gleason
# License: BSD, see bsd.txt
# cf google Doc pour source
#

def doy2dt(year,days,hours=0,minutes=0,seconds=0):
    """
    Time conversion
    Day of Year Time => Python's datetime

    Parameters
    ----------
    year, days : float or list/numpy.array of floats.
        year, days of year
    hours, minutes, seconds : float or list/numpy.array of floats, optional
        hours, minutes, seconds

    Returns
    -------
    L : datetime or list/numpy.array of datetime.
        DateTime  
    """
    if not genefun.is_iterable(year):
        # All this because Python cant handle int with a starting with 0 (like 08)
        # => SyntaxError: invalid token
        year    = int(str(year))
        days    = int(str(days))
        hours   = float(str(hours))
        minutes = float(str(minutes))
        seconds = float(str(seconds))

        tempsecs = seconds + 60 * minutes + 3600 * hours
        finalsecs     = np.floor(tempsecs)
        finalmicrosec = np.round(tempsecs * 10**6)

        return dt.datetime(year, 1, 1) + dt.timedelta(days - 1) + \
        dt.timedelta(microseconds=finalmicrosec)

    else:
        if not genefun.is_iterable(hours):
            hours   = [0] * len(year)
        if not genefun.is_iterable(minutes):
            minutes = [0] * len(year)
        if  not genefun.is_iterable(seconds):
            seconds = [0] * len(year)

        outlis = []
        for y,d,h,m,s in zip(year,days,hours,minutes,seconds):
            outlis.append(doy2dt(y,d,h,m,s))
        return outlis

def dt2doy(dtin,outputtype=str):
    return outputtype(dtin.strftime('%j'))

def dt2doy_year(dtin,outputtype=str):
    return outputtype(dtin.strftime('%j')),outputtype(dtin.strftime('%Y'))


def dt2secinday(dtin):
    return dtin.hour * 3600 + dtin.minute * 60 + dtin.second

def dt2tuple(dtin):
    return tuple(dtin.timetuple())[:-3]

def tup_or_lis2dt(lisin):
    """
    Time conversion
    Date-looking strings iterable => Python's datetime

    Parameters
    ----------
    lisin : iterable (tuple/list/numpy.array) of string.
        list of Date-looking strings
        
        like : ["2018","12","31","12","30","00"]
    Returns
    -------
    L : datetime
        converted DateTime  
    """
    
    try:
        return dt.datetime(*[int(float(e)) for e in lisin])
    except:
        return posix2dt(0)

def gpstime2utc(gpsweek,gpssecs,utc_offset):

    """
    DISCONTINUED FUNCTION, TOO UNSTABLE (171017)

    [year,month,day,hour,minute,sec] = gpstime2utc(gpsweek,gpssecs,utc_offset)
    Converts GPS week and seconds into UTC time, Python version, March 2009
    """

    from numpy import array,floor,fmod

    days_in_year = array([0,31,59,90,120,151,181,212,243,273,304,334,365])
    days_in_leap_year = array([0,31,60,91,121,152,182,213,244,274,305,335,366])
    # leap dayys since 1980
    leap_day_table = array([1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1]);

    # adjust for utc offset
    secs = gpssecs - utc_offset;

    # time as integer and fractional seconds since midnight Jan 5 1980.
    int_secs = floor(secs)
    fract_secs = secs - int_secs
    secs_1980 = 604800*gpsweek + int_secs

    # Express current GPS time as elapsed UTC days and seconds.
    day_number = floor(secs_1980/86400)
    sec_of_day = secs_1980 - day_number*86400

    leap_days1 = floor((day_number+1406) / 1461)

    # calculate UTC year
    total_years = floor((day_number + 5 - leap_days1) / 365)
    year = 1980 + total_years

    # day of utc year
    leap_days2 = sum(leap_day_table[0:total_years])

    day_of_utc_year = (day_number + 5 - 365*total_years - leap_days2)

    # determine month and day
    month = -1
    if fmod(year,4) != 0:
        # Not a leap year.
        for i in range(12):
            if day_of_utc_year > days_in_year[i]:
                month = i+1
                day = day_of_utc_year - days_in_year[i] + 1
    else:
        # A leap year.
        for i in range(12):
            if day_of_utc_year > days_in_leap_year[i]:
                month = i+1
                day = day_of_utc_year - days_in_leap_year[i] + 1

    # calculate utc hour
    hour = floor(sec_of_day/3600)

    if hour > 23:
        hour = 23

    # calculate utc minute
    minute = floor((sec_of_day - hour*3600) / 60)

    if minute > 59:
        minute = 59

    # calculate utc seconds
    sec = (sec_of_day - hour*3600 - minute*60) + fract_secs

    utc_time = [year,month,day,hour,minute,sec]

    return utc_time


def utc2gpstime_bad(year,month,day,hour,min,sec):
    from numpy import mod,array,floor
    """
    [gpsweek,gpssecs] = utc2gpstime(year,month,day,hour,min,sec)
    Converts UTC time to GPS week and seconds, Python version, Nov 2008
    
    UNSTABLE, DISCONTINUED !
    """

    # Number of days into the year at the start of each month (ignoring leap years).
    days_in_year = array([0,31,59,90,120,151,181,212,243,273,304,334,365])

    # year referenced from 1980
    ye = year - 1980;

    # Compute num of leap days
    leap_days = ye/4 + 1
    if mod(ye,4) == 0 and month <= 2:
        leap_days = leap_days - 1;

    # days elapsed since midnight jan 5, 1980.
    de = ye*365 + days_in_year[month-1] + day + leap_days - 6;

    # Leap Seconds, good after 1998
    # GPS time is ahead of UTC by this many Leap Seconds

    utc_offset = find_leapsecond(dt.datetime(year,month,day))

    # Desactivated because poor

    #    utc_offset = 13 # 1999 through 2005
    #    if year >= 2009:
    #        utc_offset = 15  # starting in 2009 to ???
    #    elif year >= 2006:
    #        utc_offset = 14 # 2006 through 2008
    #
    #    if year < 1999:
    #        print "Leap seconds invalid for this year!"
    #        return 0,0

    # Convert time to GPS weeks and seconds.
    gpsweek = floor(de/7.0)
    gpssecs = mod(de,7)*86400.0 + hour*3600.0 + min*60.0 + sec + utc_offset

    # Adjust GPS weeks/seconds to guarantee that seconds is in the range 0-604800.0. */
    if gpssecs < 0.0:
        gpsweek -= 1
        gpssecs += 604800.0

    if gpssecs >= 604800.0:
        gpsweek += 1
        gpssecs -= 604800.0

    return gpsweek,gpssecs


def utc2gpstime(year,month,day,hour,min,sec):
    
    date_utc = dt.datetime(year,month,day,hour,min,sec)
    
    utc_offset = find_leapsecond(date_utc)
    
    start_gps_time   = dt.datetime(1980,1,6)
    date_diff  = date_utc - start_gps_time + dt.timedelta(seconds=utc_offset) - dt.timedelta(seconds=19)
    
    gpsweek_decimal = date_diff.days / 7.
    gpsweek = np.floor(gpsweek_decimal)
    gpsweek_decimal_part = np.modf(gpsweek_decimal)[0]
    gpssecs = np.round(86400 * 7 * gpsweek_decimal_part) + date_diff.seconds
    
    return int(gpsweek),int(gpssecs)




def dt2gpstime(dtin,dayinweek=True):
    """
    Returns:
        dayinweek=True
            week GPS , day in week
        dayinweek=False
            week GPS , sec in week
    """
    week , secs = utc2gpstime(dtin.year,dtin.month,dtin.day,
                              dtin.hour,dtin.minute,dtin.second)
    if dayinweek:
        day = np.floor(np.divide(secs,86400))
        return int(week) , int(day)
    else:
        return int(week) , int(secs)


def dt2gpsweek_decimal(dtin,return_middle_of_day=True):
    if genefun.is_iterable(dtin):
        return np.array([dt2gpsweek_decimal(e) for e in dtin])
    else:
        if return_middle_of_day:
            mod = 0.5
        else:
            mod = 0
        week , day = dt2gpstime(dtin)
        return float(week) + (float(day + mod) / 7 )



def dt2list(dtin,return_useful_values=True):
    if return_useful_values:
        return list(dtin.timetuple())[:-3]
    else:
        return list(dtin.timetuple())



def gpstime2dt(gpsweek,gpsdow_or_seconds,dow_input = True):
    """
    Time conversion
    Day of Year Time => Python's datetime

    Parameters
    ----------
    gpsweek : int
        year, days of year
    gpsdow_or_seconds : int
        Day of Week OR Seconds in Weeks
    dow_input : bool
        select if Day of Week (True) OR Seconds in Weeks (False)
        
    Returns
    -------
    L : datetime or list/numpy.array of datetime.
        DateTime  

    Note
    ----
    Only reliable at the day level for the moment

    the leapsecond is found only after a first calc without leapsecond
    it can be some side effects when the input time is close to
    a leap second jump
    
    Source
    ------
    https://gist.github.com/jeremiahajohnson/eca97484db88bcf6b124
    """
    
    if dow_input:
        gpsseconds = gpsdow_or_seconds * 86400 + 86400*.5 # so we are around noon
    else:
        gpsseconds = gpsdow_or_seconds

    ## First gross run
    epoch   = dt.datetime(1980,0o1,0o6)
    elapsed = dt.timedelta(days=(gpsweek*7),seconds=(gpsseconds+0))

    prelim_time = epoch + elapsed


    leapsec = find_leapsecond(prelim_time)

    ## Second run with leap second
    epoch   = dt.datetime(1980,0o1,0o6)
    elapsed = dt.timedelta(days=(gpsweek*7),seconds=(gpsseconds+leapsec))

    final_time = epoch + elapsed

    if dow_input:
        final_time = final_time.date()
        final_time = dt.datetime(final_time.year, final_time.month, final_time.day)

    return final_time

def dt_gpstime2dt_utc(dtgpsin,out_array=False):
    # on converti le dt gps en dt utc
    if genefun.is_iterable(dtgpsin):
        Out = [dt_gpstime2dt_utc(e) for e in dtgpsin]
        if not out_array:
            return Out
        else:
            return np.array(Out)
    else:
        leapsec = find_leapsecond(dtgpsin)
        dtutc = dtgpsin + dt.timedelta(seconds=19) - dt.timedelta(seconds=leapsec)
        return dtutc

def toYearFraction(date):
    # Give the decimal year
    # source :
    # http://stackoverflow.com/questions/6451655/python-how-to-convert-datetime-dates-to-decimal-years

    from datetime import datetime as dt
    import time

    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

def date_string_2_dt(strin):
    import dateutil.parser
    return dateutil.parser.parse(strin)

string_date2dt = date_string_2_dt

def convert_partial_year(number):
    """
    DISCONTINUED
    use year_decimal2dt instead (which is the same)
    """

    # Source
    # http://stackoverflow.com/questions/19305991/convert-fractional-years-to-a-real-date-in-python

    import calendar
    from datetime import timedelta, datetime

    year = int(number)
    d = timedelta(days=(number - year)*(365 + calendar.isleap(year)))
    day_one = datetime(year,1,1)
    date = d + day_one
    return date

def year_decimal2dt(yearin):
    #just a wrapper
    if genefun.is_iterable(yearin):
        return [year_decimal2dt(e) for e in yearin]
    else:
        return convert_partial_year(yearin)

def dt2year_decimal(dtin):
    #just a wrapper
    if genefun.is_iterable(dtin):
        return np.array([dt2year_decimal(e) for e in dtin])
    else:
        return toYearFraction(dtin)

def jjulCNES2dt(jjulin):
    if genefun.is_iterable(jjulin):
        return [jjulCNES2dt(e) for e in jjulin]
    else:
        return dt.datetime(1950,0o1,0o1,00,00,00) + dt.timedelta(float(jjulin))

def dt2jjulCNES(dtin,onlydays=True):
    epok  = (dtin - dt.datetime(1950,0o1,0o1,00,00,00))
    if onlydays:
        return epok.days
    else:
        epok = epok + dt.timedelta(seconds=19)
        return epok.days , epok.seconds

def MJD2dt(mjd_in,seconds=None):
    # cf http://en.wikipedia.org/wiki/Julian_day
    """
    If seconds is precised, it has to be the same type/size as mjd_in
    """

    # ITERABLE CASE
    if genefun.is_iterable(mjd_in): 
        if seconds:
            if len(seconds) != len(mjd_in):
                print("ERR : MJD2dt : len(seconds) != len(mjd_in) !!")
                raise Exception
        else:
            seconds = np.zeros(len(mjd_in))

        return [dt.datetime(1858,11,17) + dt.timedelta(days=m,seconds=sec) for m,sec in zip(mjd_in,seconds)]
    
    # NON ITERABLE / FLOAT CASE
    else:
        if not seconds:
            seconds = 0    
        return dt.datetime(1858,11,17) + dt.timedelta(days=mjd_in,seconds=seconds)


def dt2MJD(dtin):
    # cf http://en.wikipedia.org/wiki/Julian_day
    try:
        return [dt2MJD(t) for t in dtin]
    except:
        delta = (dtin - dt.datetime(1858,11,17))
        return delta.days + delta.seconds / 86400.

def dt2str(dtin , str_format="%Y-%m-%d %H:%M:%S"):
    """
    just a wrapper of strftime

    source :
        https://stackoverflow.com/questions/7999935/python-datetime-to-string-without-microsecond-component
        http://www.jacksay.com/tutoriaux/bash-shell/bashshell-utilisation-commande-date.html
    """

    return dtin.strftime(str_format)

def rinexname2dt(rinexpath):

    rinexname = os.path.basename(rinexpath)
    #rinexname = rinexpath
        
    if re.search(softs_runner.rinex_regex_new_name(),rinexname):
        date_str = rinexname.split("_")[2]
        yyyy = int(date_str[:4])
        doy  = int(date_str[4:7])        
        hh   = int(date_str[7:9])        
        mm   = int(date_str[9:11])       
        dt_out = doy2dt(yyyy,doy) + dt.timedelta(seconds=hh*3600 + mm*60)
        return dt_out 
    else:
        alphabet = list(string.ascii_lowercase)
    
        if not re.search(softs_runner.rinex_regex(),rinexname):
            raise Exception('RINEX name is not well formated !!!')
        else:
            doy  = int(rinexname[4:7])
            yy   = int(rinexname[9:11])
    
        if yy > 80:
            year = yy + 1900
        else:
            year = yy + 2000
    
        if rinexname[7] in alphabet:
            h = alphabet.index(rinexname[7])
        else:
            h = 0
    
        return dt.datetime(year,1,1) + dt.timedelta(days = doy - 1 , seconds = h * 3600)

        

def sp3name2dt(sp3path):
    sp3name = os.path.basename(sp3path)

    week = int(sp3name[3:7])
    dow  = int(sp3name[7])
    return gpstime2dt(week,dow)

def sp3name_v3_2dt(sp3path):
    sp3name = os.path.basename(sp3path)
    datestr = sp3name.split("_")
    datestr = datestr[1]
    
    yyyy = int(datestr[0:4])
    doy  = int(datestr[4:7])
    hh   = int(datestr[7:9])
    mm   = int(datestr[9:11])
    
    dt_out = doy2dt(yyyy,doy) + dt.timedelta(hh * 3600 + mm * 60)
    
    return dt_out



def statname_dt2rinexname(statname,datein,rnxtype='d.Z',
                          session_a_instead_of_daily_session = False):

    alphabet = list(string.ascii_lowercase)

    if datein.hour != 0:
        sess = alphabet[datein.hour]
    else:
        if session_a_instead_of_daily_session:
            sess = alphabet[datein.hour]
        else:
            sess = '0'

    return statname + dt2doy(datein) + sess + '.' + datein.strftime('%y') + rnxtype


def datestr_sinex_2_dt(datestrin):
    
    if genefun.is_iterable(datestrin):
        return [datestr_sinex_2_dt(e) for e in datestrin]
    else:
        #### CASE WHERE THE DATE LOOKS LIKE 00000:00000
        if re.search("[0-9]{5}:[0-9]{5}",datestrin):
            datestr_list = list(datestrin)
            datestr_list.insert(2,":")
            datestrin = "".join(datestr_list)
            
    
        if datestrin == '00:000:00000':
            return dt.datetime(1970,1,1)
    
        dateint = [int(e) for e in datestrin.split(':')]
        yr = dateint[0]
        doy = dateint[1]
        sec = dateint[2]
    
        if yr > 50:
            year = 1900 + yr
        else:
            year = 2000 + yr
    
        return doy2dt(year,doy,seconds=sec)


def datestr_gins_filename_2_dt(datestrin):
    """
    datestr : 181211_124642
    """
    
    if genefun.is_iterable(datestrin):
        return [datestr_gins_filename_2_dt(e) for e in datestrin]
    
    else:
        #### CASE WHERE THE DATE LOOKS LIKE 00000:00000

        yr = int(datestrin[0:2])
        mm = int(datestrin[2:4])
        dd = int(datestrin[4:6])
        
        hh   = int(datestrin[7:9])
        mmin = int(datestrin[9:11])
        ss   = int(datestrin[11:13])
    
        if yr > 50:
            year = 1900 + yr
        else:
            year = 2000 + yr
    
        return dt.datetime(year,mm,dd,hh,mmin,ss)

        
        hh   = int(datestrin[7:9])
        mmin = int(datestrin[9:11])
        ss   = int(datestrin[11:13])
    
        if yr > 50:
            year = 1900 + yr
        else:
            year = 2000 + yr
    
        return dt.datetime(year,mm,dd,hh,mmin,ss)

#date_gins_whole_line2dt = gcls.gins_read_time
#def find_time_gap(Tin,pas=0,marge=1):
#
#    # marge = si dt > pas x marge, alors dt est un gap
#    # pas = intervalle de temps nominale entre 2 mesures
#
#    # On travaille en temps POSIX
#
#    if isinstance(Tin[0],dt.datetime):
#        T = dt2posix_list(Tin)
#    else:
#        T = Tin
#
#    dT = np.round(np.diff(T),1)
#
#    if pas == 0:
#        pas=np.round(np.min(dT),1)
#
#    print "pas nominal : %f s" %(pas)
#
#    # On cherche le pas optimal si pas précisé
#    # mais c'est mieux de l'imposer en manu
#
#    ind_list = []
#    win_list = []
#    bool_list = []
#
#    for idt, vdt in enumerate(dT):
#        ibool = True
#
#        if vdt > pas * marge :
#            print "gap trouvé : %f s (i : %i)" %(vdt,idt)
#            ind_list.append(idt)
#            win_list.append([T[idt],T[idt+1]])
#            ibool = False
#
#        bool_list.append(ibool)
#
#    return ind_list , win_list
#

def dates_middle(start,end):
    return start + (end - start)/2


def time_win_basic(start,end,Tlisin,Datalisin,outposix=True,invert=False,
                   out_array=False , out_boolis = False , only_boolis = False):
    """
    In Intern, we works in POSIX

    only_boolis : To gain speed, no operation on Tlis & Datalisin is be done
                  None is outputed for Tlisout , Datalisout 

    Outputs :
        If out_boolis == True:
            Tlisout , Datalisout , boolis
        If out_boolis == False:
            Tlisout , Datalisout
    """

    if isinstance(Tlisin[0],dt.datetime):
        Tlis  = dt2posix(Tlisin)
    else:
        Tlis = Tlisin


    if isinstance(start,dt.datetime):
        start = dt2posix(start)
    if isinstance(end,dt.datetime):
        end   = dt2posix(end)
        


    if not isinstance(Tlis,np.ndarray) or not isinstance(Datalisin,np.ndarray):    
        Tlis    = np.array(Tlis)
        Datalis = np.array(Datalisin)
    else:
        Tlis    = Tlis
        Datalis = Datalisin

    boolis = (start <= Tlis) * (Tlis <= end)

    if invert:
        boolis = np.logical_not(boolis)

    if only_boolis:
        Datalisout = None
        Tlisout    = None  
    else:
        Datalisout = Datalis[boolis]
        Tlisout    = Tlis[boolis]
    
        if not outposix:
            Tlisout = posix2dt(Tlisout)
            
        if out_array:
            Tlisout , Datalisout = np.array(Tlisout) , np.array(Datalisout)


    if out_boolis:
        out_tuple = (Tlisout , Datalisout , boolis)
    else:
        out_tuple = (Tlisout , Datalisout)
        
    return out_tuple


def time_win_multi(start,end,Tlist,Datalislis,outposix=True,
                   invert=False,out_array=False):
    Datalislisout = []
    for i,datalis in enumerate(Datalislis):
#        print 'INFO : time_win_multi : list no' , i
        Tlisout , datalisout = time_win_basic(start,end,Tlist,datalis,outposix,
                                              invert,out_array=out_array)
        Datalislisout.append(datalisout)
    return Tlisout , Datalislisout


def time_win_multi_start_end(Start_list_in,End_list_in,Tlisin,Datalisin,
                             outposix=True,invert = False,
                             out_array=False , out_boolis = False):
    
    """
    In Intern, we works in POSIX

    Outputs :
        If out_boolis == True:
            Tlisout , Datalisout , boolis_opera , boolis_stk (4 values !!)
        If out_boolis == False:
            Tlisout , Datalisout
    """
    
    if len(Start_list_in) != len(End_list_in):
        print("ERR : time_win_multi_start_end : len(Start_list_in) != len(End_list_in) !!")
    
    
    boolis_stk = []
    for start , end in zip(Start_list_in , End_list_in):
        _ , _ , boolis = time_win_basic(start,end,Tlisin,Datalisin,
                                        outposix=outposix,invert=invert,
                                        out_boolis = True , only_boolis = True)
        
        boolis_stk.append(boolis)
        
    boolis_stk   = np.stack(boolis_stk)     
    boolis_opera = np.all(boolis_stk,axis=0)
    
    Datalis = np.array(Datalisin)
    Tlis    = np.array(Tlisin)

    Datalisout = Datalis[boolis_opera]
    Tlisout    = Tlis[boolis_opera]
    
    if not outposix:
        Tlisout = posix2dt(Tlisout)
        
    if out_array:
        Tlisout , Datalisout = np.array(Tlisout) , np.array(Datalisout)

    if out_boolis:
        out_tuple = (Tlisout , Datalisout , boolis_opera , boolis_stk)
    else:
        out_tuple = (Tlisout , Datalisout)
    
    return out_tuple


def get_season(now):
    from datetime import date, datetime

    seasons = [('winter', (date(1,  1,  1),  date(1,  3, 20))),
               ('spring', (date(1,  3, 21),  date(1,  6, 20))),
               ('summer', (date(1,  6, 21),  date(1,  9, 22))),
               ('autumn', (date(1,  9, 23),  date(1, 12, 20))),
               ('winter', (date(1, 12, 21),  date(1, 12, 31)))]

    # suppressing the year
    if isinstance(now, datetime):
        now = now.date()
    now = now.replace(year=1)

    for season, (start, end) in seasons:
        if start <= now <= end:
            return season


def color_of_season(datein):
    season = get_season(datein)
    if season == 'winter':
        outcolor = 'b'
    elif season == 'summer':
        outcolor = 'r'
    elif season == 'spring':
        outcolor = 'g'
    elif season == 'autumn':
        outcolor = 'k'
    return outcolor

# Low level statistic function

def rms_mean(A):
    """
    returns RMS mean of a list/array
    """
    return np.sqrt(np.nanmean(np.square(A)))

def RMSmean(indata):
    """
    returns RMS mean of a list/array

    useless redundancy with rms_mean
    this function use shall be avoided
    """
    rms = np.sqrt(np.nanmean(np.square(indata)))
    return rms


def rms_mean_alternativ(A):
    """
    returns "GRGS style" RMS of a list/array
    the arithmetic mean of the values is substracted from the values
    NB 1808 : It is basically the standard deviation ...

    i.e.    _
    √< (A - A)^2 > instead of √< (A)^2 >
    """
    return np.sqrt(np.nanmean(np.square(A - np.nanmean(A))))


def mad(data,mode='median'):
    """
    returns Median Absolute Deviation (MAD) a list/array
    """

    if mode == 'median':
        MAD = np.nanmedian(np.abs(data - np.nanmedian(data)))
    elif mode == 'mean':
        MAD = np.nanmean(np.abs(data - np.nanmean(data)))
    return MAD


def outlier_mad(data,seuil=3.5,verbose=False,convert_to_np_array=True,
                mad_mode = 'median' ):
    
    """    
    clean the outlier of Y usind MAD approach 
    and clean the corresponding values in X
    assuming that we have the function : X => Y(X)
    (be carefull, Y is the first argument)
    
    Parameters
    ----------
    data : list or numpy.array
        Values

    seuil : float
        MAD threshold        
        
    verbose : bool
        
    convert_to_np_array : bool
        if True returns output as an array, if False as a regular list
        
    mad_mode : str
        'median' or 'mean' : MAD can also be based on mean (for experimental purposes)
  
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

    if convert_to_np_array:
        data = np.array(data)
    nbinp = float(len(data))
    MAD = mad(data,mode=mad_mode)
    med = np.nanmedian(data)

    if np.isclose(np.sum(np.abs(np.diff(data))) , 0.):
        if verbose:
            print("ratio d'elimination : 0 , données toutes égales")
        return data , np.array([True] * len(data))

    if np.isclose(med , 0.):
        if verbose:
            print("ratio d'elimination : 0 , mediane nulle")
        return data , np.array([True] * len(data))

    if np.isclose(MAD,0.):
        if verbose:
            print("ratio d'elimination : 0 , MAD nulle")
        return data , np.array([True] * len(data))


    diff = data - med
    MZS = 0.6745 * np.abs(diff) / MAD
    MZS[np.isnan(MZS)] = seuil * 10
    boolbad = MZS < seuil
    dataout = data[boolbad]
    nbout = float(sum(boolbad))
    ratio = (nbinp-nbout)/nbinp
    if verbose:
        print("ratio d'elimination : %i / %i, %f p.c." %(nbinp-nbout,nbinp,ratio * 100))
    return dataout , boolbad

def outiler_mad(data,seuil=3.5,verbose=False,convert_to_np_array=True,
                mad_mode = 'median' ):
    """
    wrapper of outlier_mad, maintened for legacy with a typo
    """
    return outlier_mad(data,seuil , verbose , convert_to_np_array , mad_mode )


def outlier_mad_binom(Y,X,seuil=3.5,verbose=False,detrend_first=False,
                      return_booleans = False):   
    """    
    clean the outlier of Y usind MAD approach 
    and clean the corresponding values in X
    assuming that we have the function : X => Y(X)
    (be carefull, Y is the first argument)
    
    Parameters
    ----------
    Y : list or numpy.array
        Values

    X : list or numpy.array
        X Values so as X => Y(X)

    seuil : float
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
        Xwork , Ywork = detrend_timeseries(X,Y)
    else:
        Xwork , Ywork = np.array(X) , np.array(Y)
        
    _ , bb = outiler_mad(Ywork,seuil,verbose)
    
    Xclean = np.array(X)[bb]    
    Yclean = np.array(Y)[bb]
    
    if not return_booleans:
        return Yclean , Xclean
    else:
        return Yclean , Xclean , bb

def outlier_mad_binom_legacy(X,Y,seuil=3.5,verbose=False,detrend_first=False,
                      return_booleans = False):
    """
    clean the outlier of X and clean the corresponding values in Y
    
    legacy : order of X Y is different than in the main version, and here 
    it might be unstable for the detrend
    """
    if detrend_first:
        Xwork , _ = detrend_timeseries(X,Y)
    else:
        Xwork , _ = np.array(X) , np.array(Y)
        
    _ , bb = outiler_mad(Xwork,seuil,verbose)
    
    Xclean = np.array(X)[bb]    
    Yclean = np.array(Y)[bb]
    
    if not return_booleans:
        return Xclean , Yclean
    else:
        return Xclean , Yclean , bb
        

def outlier_above_below_simple(X , low_bound  , upp_bound,
                        return_booleans = True):    
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
        print("WARN : outlier_above_below_simple : lower bound >= upper bound !!!")
        print("      low_bond : " , low_bound)
        print("      upp_bond : " , upp_bound)
    
    bbool = (low_bound <= Xwork) & (Xwork <= upp_bound)
    
    Xout = Xwork[bbool]
    
    if return_booleans:
        return Xout, bbool
    else:
        return Xout
        

def outlier_above_below(X , threshold_values ,
                        reference = np.nanmean  , 
                        theshold_absolute = True,
                        return_booleans   = True,
                        theshold_relative_value = "reference",
                        verbose = False):

    """    
    Gives values of X which are between threshold values

    Parameters
    ----------
    threshold_values : single value (float) or a 2-tuple 
        (lower bound theshold , upper bound theshold)
        
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
    
    if genefun.is_iterable(threshold_values):
        ths_input_low = threshold_values[0]
        ths_input_upp = threshold_values[1]
    else:
        ths_input_low = threshold_values
        ths_input_upp = threshold_values        
        
    if ths_input_low < 0. or ths_input_upp < 0.:
        print("WARN : outlier_above_below : threshold_values have to be positive")
        print("       minus sign for lower bound will be applied internally")
        
    
    if callable(reference):
        ref_val = reference(X)
    else:
        ref_val = reference
        
    if theshold_relative_value in ("reference" , None):
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
        print("INFO : outlier_above_below theshold values")
        print("       reference : " , ref_val )
        print("       effective lower bound : " , ths_low )
        print("       effective upper bound : " , ths_upp )
                
    Xout , bbool = outlier_above_below_simple(X , ths_low , ths_upp)
    
    if return_booleans:
        return Xout , bbool
    else:
        return Xout
    
    
def outlier_above_below_binom(Y , X , 
                              threshold_values ,
                              reference = np.nanmean  , 
                              theshold_absolute = True,
                              theshold_relative_value = "reference",
                              return_booleans   = False,
                              detrend_first     = True,
                              verbose           = False):
    
    
    """    
    Gives values of Y which are between threshold values, and correct an 
    associated X so as X => Y(X)

    Parameters
    ----------
    threshold_values : single value (float) or a 2-tuple 
        (lower bound theshold , upper bound theshold)
        
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
        Xwork , Ywork = detrend_timeseries(X,Y)
    else:
        Xwork , Ywork = np.array(X) , np.array(Y)
        

    _ , bb = outlier_above_below(Ywork , threshold_values ,
                        reference = reference , 
                        theshold_absolute = theshold_absolute,
                        theshold_relative_value=theshold_relative_value,
                        return_booleans   = True,
                        verbose = verbose)
    
    Xclean = np.array(X)[bb]    
    Yclean = np.array(Y)[bb]
    
    if not return_booleans:
        return Yclean , Xclean
    else:
        return Yclean , Xclean , bb


def outlier_sigma(datasigmain,seuil=3):
    """
    si un point a un sigma > seuil * moy(sigmas) on le vire
    
    really old and discontinued, and not really efficient
    """
    moy = np.median(datasigmain)
    marge = moy * seuil

    print("INFO : outlier_sigma : moy,seuil,marge",  moy,seuil,marge)

    boolbad = np.abs(datasigmain) < marge

    datasigmaout = datasigmain[boolbad]

    return datasigmaout,boolbad


def outlier_overmean(Xin,Yin,marge=0.1):
    """
    really old and discontinued, use outlier_above_below instead
    """

    # elimine les points qui sont au dela d'une certaine marge au dessus de la moyenne

    nbinp = float(len(Yin))

    # Xin sont des array pas des listes

    moy = np.abs(np.nanmean(Yin))

    boolbad = np.abs(np.abs(Yin) - np.abs(moy)) < marge

    Yout = Yin[boolbad]
    Xout = Xin[boolbad]


    nbout = float(sum(boolbad))

    ratio = (nbinp-nbout)/nbinp
    print("ratio d'elimination : %i / %i, %f" %(nbinp-nbout,nbinp,ratio))
    print("moyenne : %f" %(moy))

    plt.figure(12)
    plt.clf()
    plt.plot(Xin,[moy + marge] * len(Xin))
    plt.plot(Xin,[moy - marge] * len(Xin))
    plt.plot(Xin,[moy] * len(Xin))
    plt.plot(Xin,Yin,'*')
    plt.plot(Xout,Yout,'+')


    return Xout , Yout, boolbad


def BL_from_points(listpointin):
    ''' 
    A partir d'une liste de points,
    retourne les baselines entre ces points dans une matrice 
    '''
    """    
    From a list of 2-D or 3-dD points, returns the a matrix with distance 
    between each points 
    
    Parameters
    ----------
    listpointin : list or numpy.array
        List of N 2D or 3D points [[x1,y1,z1] ... [xn , yn , zn]]
                
    Returns
    -------
    BL : numpy.array
        matrix with distances between each points
        
    """

    N = len(listpointin)
    BL = np.empty((N,N))

    for i,pt1 in enumerate(listpointin):
        for j,pt2 in enumerate(listpointin):

            if i == j:
                BL[i,j] = 0
            else:
                BL[i,j] = np.linalg.norm(pt1 - pt2)

    return BL



def mat_poids(Sinp,Ninp,fuvinp=1):
    """
    discontinued
    """
    # Sinp : liste des Sigmas sig = sqrt(var)
    # Ninp : liste de la taille de chaque blocs (obs)
    # fuvinp = 1 : facteur unitaire de variance

    if len(Sinp) != len(Ninp):
        raise Exception("S et N de taille differente")

    Ktemp = []

    for i in range(len(Sinp)):
        print(Sinp[i])
        Ktemp.append(np.eye(Ninp[i]) * Sinp[i]**2)

    K = scipy.linalg.block_diag(*Ktemp)
    Q = (1/fuvinp) * K
    P = scipy.linalg.inv(Q)

    return K , Q , P

def rotmat2(theta,angtype='deg'):

    if angtype == 'deg':
        theta = np.deg2rad(theta)

    rotmat = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])

    return rotmat


def rotmat3(alpha,beta,gamma,xyzreftuple = ([1, 0, 0], [0, 1, 0], [0, 0, 1]),angtype='deg'):

    xaxis, yaxis, zaxis = xyzreftuple

    if angtype == 'deg':
        alpha = np.deg2rad(alpha)
        beta  = np.deg2rad(beta)
        gamma = np.deg2rad(gamma)

    Rx = trans.rotation_matrix(alpha, xaxis)
    Ry = trans.rotation_matrix(beta, yaxis)
    Rz = trans.rotation_matrix(gamma, zaxis)
    R = trans.concatenate_matrices(Rz, Ry, Rx)[:3,:3]

    return R

def rotate_points(alphal,betal,gammal,pointlin,Rtype='R1',
                  xyzreftuple = ([1, 0, 0], [0, 1, 0], [0, 0, 1]),
                  angtype='deg',fullout = False):
    '''
    R1  = Rz(g) * Ry(b) * Rx(a)
         si les RPY sont donnés dans le NED
         alors les positions résultantes sont dans le NED

    R2  =  matrice RPY2ENU
        si les RPY sont donnés dans le NED
        alors les  résultantes sont DANS LE ENU
        pas besoin de rotation NED2ENU

        Grewal et al. 2007

    Entrée :
        Angles n = A
        liste de listes de P * [ points ]

    Sortie :
        liste de listes [ [ xA ] [ xA ] ... xP [ xA ] ]  '''

    xaxis, yaxis, zaxis = xyzreftuple

    if not genefun.is_iterable(alphal):
        alphal = np.array([alphal])
        betal = np.array([betal])
        gammal = np.array([gammal])
        boolnotiterable = True
    else:
        boolnotiterable = False

    pointlout = []
    R_out = []


    for pt in pointlin:

        if not genefun.is_iterable(pt) or len(pt) != 3:
            print("ERR : rotate_points : pts != 3 coords")
            return 0

        pointltmp = []

        for a,b,g in zip(alphal,betal,gammal):

            R1 = rotmat3(a,b,g,angtype=angtype,xyzreftuple=xyzreftuple)
            R2 = C_rpy2enu(a,b,g,angtype=angtype)

            if Rtype == 'R1':
                R = R1
            elif Rtype == 'R2':
                R = R2
            R_out.append(R)

            pointltmp.append(np.dot(R,pt))

        pointlout.append(pointltmp)

        if boolnotiterable:
            pointlout = pointltmp

        pointlout = np.array(pointlout)

    if fullout:
        return pointlout , R_out
    else:
        return pointlout


def guess_seq_len(seq):
    #source
    #http://stackoverflow.com/questions/11385718/python-finding-repeating-sequence-in-list-of-integers
    guess = 1

    if len(set(seq)) == 1:
        return 1

    max_len = len(seq) / 2
    for x in range(2, max_len):
        if seq[0:x] == seq[x:2*x] :
            return x

    return guess

def wrapTo2Pi(lon):
    """
     wrapTo2Pi Wrap angle in radians to [0 2*pi]

    lambdaWrapped = wrapTo2Pi(LAMBDA) wraps angles in LAMBDA, in radians,
    to the interval [0 2*pi] such that zero maps to zero and 2*pi maps
    to 2*pi. (In general, positive multiples of 2*pi map to 2*pi and
    negative multiples of 2*pi map to zero.)

    See also wrapToPi, wrapTo180, wrapTo360.

    """
    lon = np.array(lon)
    positiv = lon > 0
    outlon = np.mod(lon , 2*np.pi)
    outlon[np.logical_and(outlon == 0 , positiv)] = 2 * np.pi
    return outlon


def wrapToPi(lon):
    """
    wrapToPi Wrap angle in radians to [-pi pi]

       lambdaWrapped = wrapToPi(LAMBDA) wraps angles in LAMBDA, in radians,
       to the interval [-pi pi] such that pi maps to pi and -pi maps to
       -pi.  (In general, odd, positive multiples of pi map to pi and odd,
       negative multiples of pi map to -pi.)

       See also wrapTo2Pi, wrapTo180, wrapTo360.

    """

    outlon = np.array(lon)
    q =  np.logical_and((outlon < -np.pi) , (np.pi < outlon))
    outlon[q] = wrapTo2Pi(outlon[q] + np.pi) - np.pi
    return outlon


def wrapTo180(lonin):
    """
    wrapTo180 Wrap angle in degrees to [-180 180]

    lonWrapped = wrapTo180(LON) wraps angles in LON, in degrees, to the
    interval [-180 180] such that 180 maps to 180 and -180 maps to -180.
    (In general, odd, positive multiples of 180 map to 180 and odd,
    negative multiples of 180 map to -180.)

    See also wrapTo360, wrapTo2Pi, wrapToPi.
    """
    lon = np.array(lonin)
    q = (lon < -180) and (180 < lon)
    lon[q] = wrapTo360(lon[q] + 180) - 180

    return lon

def wrapTo360(lonin):
    """
    wrapTo360 Wrap angle in degrees to [0 360]

    lonWrapped = wrapTo360(LON) wraps angles in LON, in degrees, to the
    interval [0 360] such that zero maps to zero and 360 maps to 360.
    (In general, positive multiples of 360 map to 360 and negative
    multiples of 360 map to zero.)

    See also wrapTo180, wrapToPi, wrapTo2Pi.
    """

    lon = np.array(lonin)

    positiveInput = (lon > 0)
    lon = np.mod(lon, 360)
    lon[(lon == 0) & positiveInput] = 360
    return lon



# Pas convaincu de son utilité
def unwrap180(anglist,angtype='deg'):

    if angtype == 'deg':
        seuil = 360

    angout = []

    for a in anglist:
        if a > seuil / 2:
            a = a - seuil
        angout.append(a)

    return angout


def wrap360(anglist,angtype='deg'):

    angout = []

    if angtype == 'deg':
        seuil = 360
    elif angtype == 'rad':
        seuil = 2*np.pi

    for a in anglist:
        if a < 0:
            a = a + seuil

        angout.append(a)

    return angout

class interp1d_ang():

    def __init__(self,T,A,angtype='deg',kind='linear',bounds_error=False):

        if angtype == 'deg':
            A = np.deg2rad(A)

        self.A = A
        self.T = T
        self.C = np.cos(A)
        self.S = np.sin(A)

        self.CfT = interpolate.interp1d(T,self.C,kind=kind,bounds_error=bounds_error)
        self.SfT = interpolate.interp1d(T,self.S,kind=kind,bounds_error=bounds_error)


    def __call__(self,T,angtype='deg'):

        I = np.arctan2(self.SfT(T) ,self.CfT(T) )
        I = wrap360(I,angtype='rad')

        if angtype == 'deg':
            return np.rad2deg(I)
        else:
            return I

def group_consecutives(vals, step=1):
    """
    Return list of consecutive lists of numbers from vals (number list).
    """
    run = []
    result = [run]
    expect = None
    for v in vals:
        if (v == expect) or (expect is None):
            run.append(v)
        else:
            run = [v]
            result.append(run)
        expect = v + step

        result2 = []
        for r in result:
            if len(r) > 1:
                result2.append([r[0],r[-1]])
            else:
                result2.append(r)

    return result2

def linear_regression(x,y,fulloutput=False,alpha=.95):
    """    
    From 2 vectors X and Y, returns linear regression coefficients a and b

    Parameters
    ----------
    X & Y : list or numpy.array
        Values

    fulloutput : bool
        full output
        
    alpha : float
        alpha value for the confidence interval
                
    Returns
    -------
    a & b : float
        Linear regression coefficients
    
    If fulloutput == True:
        
    confid_interval_slope : float
        confid_interval_slope
        
    std_err : float
        standard deviation
        
    Note
    ----
    http://glowingpython.blogspot.fr/2012/03/linear-regression-with-numpy.html 
    
    This function is doing more or less the same job as scipy.stats.linregress
    """

    # On bosse avec des arrays
    x = np.array(x)
    y = np.array(y)

    if len(x) != len(y):
        print("ERR : linear_regression : len(x) != len(y)")
        print("      len(x) : " , len(x))
        print("      len(y) : " , len(y))

        return 0,0

    A = np.array([x, np.ones(len(x))])
    # linearly generated sequence
    w = np.linalg.lstsq(A.T,y,rcond=None)[0] # obtaining the parameters

    if not fulloutput:
        return w[0],w[1]
    else:
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x,y)
        return w[0],w[1],confid_interval_slope(x,y,alpha),std_err


def linear_reg_getvalue(X,a,b,full=True):
    """    
    From 2 vector X and coefficients a & b, get Y = a*X + b

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

def linear_coef_a_b(x1,y1,x2,y2):

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

    if genefun.is_iterable(x1):
        x1 = np.array(x1,dtype=np.float64)
        x2 = np.array(x2,dtype=np.float64)
        y1 = np.array(y1,dtype=np.float64)
        y2 = np.array(y2,dtype=np.float64)
    else:
        x1 = float(x1)
        x2 = float(x2)
        y1 = float(y1)
        y2 = float(y2)

    a = (y2 - y1) / (x2 - x1)
    b1 = y1 - a*x1
    b2 = y2 - a*x2
    return a , b1 , b2

def detrend_timeseries(X,Y):
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
    a,b = linear_regression(X,Y)

    #Yout = Y - a * (X - X[0])
    #Yout = Y - ( a * X + b )  
    
    Ylinear =  ( a * X + b )
    Yout    = Y - Ylinear + Y[0] 
    
    return X , Yout


def confid_interval_slope(x,y,alpha=.95):
    """
     Calcule un intervalle de confiance sur une tendance
     En entrée: x     = la variable indépendante
                y     = la variable dépendante
                alpha = la probabilité d'erreur tolérée
     En sortie: mi    = la borne inférieure de l'intervalle
                ma    = la borne supérieure de l'intervalle
                
    Source (???? => En fait non ...)
    http://www.i4.auc.dk/borre/matlab
    http://kom.aau.dk/~borre/matlab/
    """

    sux=np.sum(x)
    xb=np.mean(x)
    suy=np.sum(y)
    yb=np.mean(y)
    n=len(x)
    S1=np.sum(x*y)
    S2=sux*suy/n
    Sxy=S1-S2
    S4=np.sum(x**2)
    S5=(sux**2)/n
    Sxx=S4-S5
    S7=np.sum(y**2)
    S8=(suy**2)/n
    Syy=S7-S8
    b1=Sxy/Sxx
    b0=yb-b1*xb
    S14=(Sxy**2)/Sxx
    s2y=(Syy-S14)/(n-2)
    sy=np.sqrt(s2y)
    s2b1=s2y/Sxx
    s2b0=s2y*(1/n+(xb**2)/Sxx)
    #t=tq(1-alpha/2,n-2)
    t = scipy.stats.t.ppf(1-alpha/2,n-2)
    mi=b1-t*np.sqrt(s2b1)
    ma=b1+t*np.sqrt(s2b1)
    return mi,ma

def plot_vertical_bar(xlis , color='r',linewidth=1):
    out_bar_list = []
    for x in xlis:
        out_bar = plt.axvline(x,color=color,linewidth=linewidth)
        out_bar_list.append(out_bar)
    return out_bar_list

def plot_vertical_bar_ax(xlis,ax_in,color='r',linewidth=1):
    out_bar_list = []
    for x in xlis:
        out_bar = ax_in.axvline(x,color=color,linewidth=linewidth)
        out_bar_list.append(out_bar)
    return out_bar_list


def running_mean(data_in , window , convolve_mode="same"):   
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

Y1 = geok.movingaverage(Y,N)
Y2 = geok.runningMean(Y,N)
Y3 = geok.running_mean_core(Y,N)
Y4a = geok.movingaverage_bis(Y,N,"same")
Y4b = geok.movingaverage_bis(Y,N,"full")
Y5 = geok.movingaverage_ter(Y,N)

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

def movingaverage(values,window):
    """
    including valid will REQUIRE there to be enough datapoints.
    for example, if you take out valid, it will start @ point one,
    not having any prior points, so itll be 1+0+0 = 1 /3 = .3333
    http://sentdex.com/sentiment-analysisbig-data-and-python-tutorials-algorithmic-trading/how-to-chart-stocks-and-forex-doing-your-own-financial-charting/calculate-simple-moving-average-sma-python/
    
    INTERNAL_ID_1
    """
    weigths = np.repeat(1.0, window)/window
    smas = np.convolve(values, weigths, 'valid')
    return smas # as a numpy array

def runningMean(x, N):
    """
    http://stackoverflow.com/questions/13728392/moving-average-or-running-mean
    
    INTERNAL_ID_2
    """
    y = np.zeros((len(x),))
    for ctr in range(len(x)):
         y[ctr] = np.sum(x[ctr:(ctr+N)])
    return y/N

def running_mean_core(x, N):
    """
    moyenne glissante
    https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
    Alleo answer
    
    INTERNAL_ID_3
    """
    cumsum  = np.cumsum(np.insert(x, 0, 0))
    xout = (cumsum[N:] - cumsum[:-N]) / N

    return xout

def movingaverage_bis(interval, window_size , convolve_mode="same"):
    """
    moyenne glissante, plus lente mais donne une sortie de meme taille que l'entrée
    https://stackoverflow.com/questions/11352047/finding-moving-average-from-data-points-in-python
    
    INTERNAL_ID_4
    """
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, convolve_mode)


def movingaverage_ter(data, window_width):
    """
    https://stackoverflow.com/questions/11352047/finding-moving-average-from-data-points-in-python
    Roman Kh ans
    
    INTERNAL_ID_5  
    """
    cumsum_vec = np.cumsum(np.insert(data, 0, 0)) 
    ma_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
    
    return ma_vec


def harmonic_mean(A):
    """
    harmonic mean of a list/array A
    """

    A = np.array(A)
    return len(A) / np.sum(1.0/A)

def find_intersection(x1,y1,x2,y2):
    #http://stackoverflow.com/questions/8094374/python-matplotlib-find-intersection-of-lineplots

    import scipy.interpolate as interpolate
    import scipy.optimize as optimize

    p1=interpolate.PiecewisePolynomial(x1,y1[:,np.newaxis])
    p2=interpolate.PiecewisePolynomial(x2,y2[:,np.newaxis])

    def pdiff(x):
        return p1(x)-p2(x)

    xs=np.r_[x1,x2]
    xs.sort()
    x_min=xs.min()
    x_max=xs.max()
    x_mid=xs[:-1]+np.diff(xs)/2
    roots=set()
    for val in x_mid:
        root,infodict,ier,mesg = optimize.fsolve(pdiff,val,full_output=True)
        # ier==1 indicates a root has been found
        if ier==1 and x_min<root<x_max:
            roots.add(root[0])
    roots=np.array(list(roots))
    return roots,p1(roots)

def wrapTo360(lon):
    # according to the MATLAB fct
    lon = np.mod(lon, 360)
    return lon

def wrapTo180(lon):
    # according to the MATLAB fct
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


def rinex_lister(path,add_new_names=True):
    """
    find all rinex in a folder and his subfolders
    path can be a string or a tuple of string => manage multi files :)

    is very similar with softs_runner.multi_finder_rinex and
    gins_runner.get_rinex_list
    """

    if type(path) is str:
        path = [path]

    paths_walked = []
    for p in path:
        paths_walked = paths_walked + list(os.walk(p))

    wholefilelist = []
    for tup in paths_walked:
        wholefilelist = wholefilelist + tup[-1]

    wholefilelist = list(set(wholefilelist))
    
    rinexfilelist           = [fil for fil in wholefilelist if re.search( softs_runner.rinex_regex() , os.path.basename(fil))]
    if add_new_names:
        rinexfilelist_new_names = [fil for fil in wholefilelist if re.search( softs_runner.rinex_regex_new_name() , os.path.basename(fil))]
        rinexfilelist = rinexfilelist + rinexfilelist_new_names
    
    print('INFO : ' , len(rinexfilelist) , 'RINEXs found')

    return rinexfilelist


def rinex_timeline(inputlist_or_paths,start = dt.datetime(1980,1,1) ,
                   end = dt.datetime(2099,1,1),use_rinex_lister = True,
                   dots_plot=False,jul_date_plot=False,return_figure=False):
    """
    if use_rinex_lister = True :
    inputlist_or_paths is a path OR a list of path
    where RINEX can be found (use the subfunction rinex_lister for that)
    else if use_rinex_lister = False:
    it's a list of RINEXs
    """

#    if use_rinex_lister:
#        filelist = rinex_lister(inputlist_or_paths)
#    else:
#        filelist = inputlist_or_paths
#
#    rinexfilelist = [fil for fil in filelist if re.match( '.*' + softs_runner.rinex_regex() + '$', fil)]
#
#    if not use_rinex_lister:
#        rinexfilelist = [os.path.basename(e) for e in rinexfilelist]
#
#    print('INFO : ', len(rinexfilelist), 'RINEXs will be ploted on the timeline')
#
#    statname_lis = sorted(list(set([rin[0:4] for rin in rinexfilelist])))
#
#    print('INFO : ', len(statname_lis), 'stations will be ploted on the timeline')
#
#    datadico = dict()
#
#    for stat in statname_lis:
#        datadico[stat] = []
#
#    for rnx in rinexfilelist:
#        try:
#            datadico[rnx[0:4]].append((rnx,rinexname2dt(rnx)))
#        except:
#            print('error with : ', rnx)
#

    datadico = rinex_timeline_datadico(inputlist_or_paths,
                                       use_rinex_lister = use_rinex_lister)

    fig  = timeline_plotter(datadico,start = start ,
                     end=end ,dots_plot=dots_plot,
                     jul_date_plot=jul_date_plot)

#    ax.xaxis.grid(True)
#    plotstat_lis = []
#    for i,stat in enumerate(reversed(sorted(datadico.keys()))):
#        print(stat)
#        T = [ e[-1] for e in datadico[stat] ]
#        T = [ t for t in T if start <= t <= end ]
#        T = sorted(T)
#
#        plotstat_lis.append(stat)
#        if dots_plot:
#            #old old with dots
#            ax.plot(T,i*np.ones(len(T)), '.')
#        else:
#            TGrp = genefun.consecutive_groupIt(dt2MJD(T),True)
#            for tgrp in TGrp:
#                if not jul_date_plot:
#                    tgrp = MJD2dt(tgrp)
#                if tgrp[0] == tgrp[1]:
#                    ax.plot(tgrp[0],i, '.')
#                else:
#                    ax.plot(tgrp,[i]*2, '-')
#
##    ax.set_yticks(np.arange(0,len(plotstat_lis)-1),plotstat_lis)
#    plt.yticks(np.arange(0,len(plotstat_lis)),plotstat_lis)
#    fig.autofmt_xdate()
#    fig.set_size_inches(11.69,i * 0.28) #16.53
#    ax.set_ylim([-1 , len(plotstat_lis) + 1])
##    return fig
    if not return_figure:
        return datadico
    else:
        return datadico , fig

def rinex_timeline_datadico(inputlist_or_paths,use_rinex_lister = True,
                            optional_info=''):
    """
    convention for RINEX datadico :
        datadico[stat] = [(rinexname1,optional1,date1) ... (rinexnameN,optionalN,dateN)]
    """
    if use_rinex_lister:
        filelist = rinex_lister(inputlist_or_paths)
    else:
        filelist = inputlist_or_paths

    #rinexfilelist = [fil for fil in filelist if re.search( '.*' + softs_runner.rinex_regex() + '$', fil)]
    rinexfilelist_old = [fil for fil in filelist if re.search( softs_runner.rinex_regex() , fil)]
    rinexfilelist_new = [fil for fil in filelist if re.search( softs_runner.rinex_regex() , fil)]

    rinexfilelist = rinexfilelist_old + rinexfilelist_new

    if not use_rinex_lister:
        rinexfilelist = [os.path.basename(e) for e in rinexfilelist]

    print('INFO : ', len(rinexfilelist), 'RINEXs will be ploted on the timeline')

    statname_lis = sorted(list(set([rin[0:4] for rin in rinexfilelist])))

    print('INFO : ', len(statname_lis), 'stations will be ploted on the timeline')

    datadico = dict()

    for stat in statname_lis:
        datadico[stat] = []

    for rnx in rinexfilelist:
        try:
            datadico[rnx[0:4]].append((rnx,optional_info,rinexname2dt(rnx)))
        except:
            print('error with : ', rnx)

    return datadico

def rinex_timeline_datadico_merge_not_very_smart(datadico_list,priority_list):
    """
    Merge different RINEXs datadico, produced by rinex_timeline_datadico
    coming from different archives
    Args :
        rinex_timeline_datadico : list of RINEX datadico
        priority_list : priority list of 'optional_info' (archive ID)
                        it will erase optional_info of lower priority
    Returns :
        datadico_out : a merged datadico
    """

    datadico_out  = dict()

    datadico_merged = genefun.dicts_of_list_merge(*datadico_list)

    for k , dataval in datadico_merged.items():

        rnxname_list = [e[0]  for e in dataval]
        archive_list = [e[1]  for e in dataval]
        date_list    = [e[-1] for e in dataval]

        out_date_list , out_all_list = [] , []
        for r,a,d in zip(rnxname_list,archive_list,date_list):
            if d not in out_date_list:
                out_date_list.append(d)
                out_all_list.append((r,a,d))
            else:
                ind_existing   = out_date_list.index(d)
                archd_existing = out_all_list[ind_existing][1]
                if priority_list.index(a) < priority_list.index(archd_existing):
                    out_date_list.remove(d)
                    out_all_list.remove(out_all_list[ind_existing])
                    out_date_list.append(d)
                    out_all_list.append((r,a,d))

        datadico_out[k] = out_all_list

    return datadico_out


def rinex_timeline_datadico_merge(datadico_list,priority_list=None):
    """
    Merge different RINEXs datadico, produced by rinex_timeline_datadico
    coming from different archives
    Args :
        rinex_timeline_datadico : list of RINEX datadico
        priority_list : priority list of 'optional_info' (archive ID)
                        it will erase optional_info of lower priority
                        NB : it is not very useful, just sort
                             datadico_list in the right order ...
    Returns :
        datadico_out : a merged datadico
    """

    datadico_out  = dict()

    datadico_merged = genefun.dicts_of_list_merge(*datadico_list)

    for k , dataval in datadico_merged.items():

        rnxname_list = [e[0]  for e in dataval]
        archive_list = [e[1]  for e in dataval]
        date_list    = [e[-1] for e in dataval]

        out_date_list , out_all_list = [] , []
        for r,a,d in zip(rnxname_list,archive_list,date_list):
            if d not in out_date_list:
                out_date_list.append(d)
                out_all_list.append((r,a,d))
            elif priority_list:
                ind_existing   = out_date_list.index(d)
                archd_existing = out_all_list[ind_existing][1]
                if priority_list.index(a) < priority_list.index(archd_existing):
                    out_date_list.remove(d)
                    out_all_list.remove(out_all_list[ind_existing])
                    out_date_list.append(d)
                    out_all_list.append((r,a,d))
            else:
                pass

        datadico_out[k] = out_all_list

    return datadico_out

def timeline_plotter(datadico,start = dt.datetime(1980,1,1) ,
                     end = dt.datetime(2099,1,1),dots_plot=False,
                     jul_date_plot=False,datadico_anex_list = [],
                     use_only_stats_of_main_datadico=False,
                     colordico_for_main_datadico=None):
    """
    A simpler version has been commited to geodezyx toolbox for archive
    on 20180118 15:59A
    """

    fig , ax = plt.subplots()
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)

    if not use_only_stats_of_main_datadico:
        stats_concat_list = list(datadico.keys()) + sum([list(e.keys()) for e in datadico_anex_list], [])
        stats_concat_list = list(reversed(sorted(list(set(stats_concat_list)))))
    else:
        stats_concat_list = list(reversed(sorted(list(datadico.keys()))))

    # the plot has not the same behavior if it is the morning or not 
    # (rinexs timelines wont be ploted if it is the morning)
    if dt.datetime.now().hour < 12:
        morning_shift = dt.timedelta(days=1)
    else:
        morning_shift = dt.timedelta(days=0)

    legend_list = [] # must be here before the loop
    for i,stat in enumerate(stats_concat_list):
        # PART 1 : PLOT MAIN DATADICO
        if not stat in datadico.keys():
            continue

        #T = Time, O = Station name (Observation)
        Torig = [ e[-1] for e in datadico[stat] ]
        Oorig = [ e[1]  for e in datadico[stat] ]

        # Time windowing
        T , O = [] , []
        for t , o in zip(Torig , Oorig):
            if ( start - morning_shift ) <= t <= end: 
                T.append(t)
                O.append(o)

        T,O = genefun.sort_binom_list(T,O)

        TMJD=dt2MJD(T)
        TGrpAll = genefun.consecutive_groupIt(TMJD,True) # Tuples (start,end) of continue period

        for tgrp in TGrpAll:
            color1 = ''
            color2 = ''
            extra_archive = False

            ###*** managing colors
            if colordico_for_main_datadico:
                igrpstart    = TMJD.index(tgrp[0])
                igrpend      = TMJD.index(tgrp[1])
                Ogrp         = O[igrpstart:igrpend+1]
                Tgrp         = TMJD[igrpstart:igrpend+1] # all dates in the current continue period

                opt_set      = list(set(Ogrp))

                if len(opt_set) == 1: # Regular case : only one archive
                    if opt_set[0] in colordico_for_main_datadico:
                        color1 = colordico_for_main_datadico[opt_set[0]]
                else: # several archives ... so the line has to be splited in segments
                    extra_archive = True
                    OSubgrp = genefun.identical_groupIt(Ogrp)
                    TSubgrp = genefun.sublistsIt(Tgrp , [len(e) for e in OSubgrp])

                    Tgrp_plt      = [ (e[0] , e[-1] + 1)    for e in TSubgrp ] # +1 because the end boundary day is not included
                    Ogrp_plt      = [ list(set(e))[0] for e in OSubgrp     ]
                    Color_grp_plt = [ colordico_for_main_datadico[e] if e in colordico_for_main_datadico else color2 for e in Ogrp_plt ]
            ###*** End of managing colors

            tgrp = (tgrp[0] , tgrp[1] + 1) # +1 because the end boundary day is not included
                                           # must stay there, in case of not colordico_for_main_datadico

            if not jul_date_plot:
                tgrp = MJD2dt(tgrp)
                if extra_archive:
                    Tgrp_plt = [ (MJD2dt(e[0]) , MJD2dt(e[1])) for e in Tgrp_plt ]

            #PLOT part
            if tgrp[0] == tgrp[1] + dt.timedelta(days=1): # CASE NO PERIOD, ONLY ONE DAY
                ax.plot(tgrp[0],i , '.' + color1)
            else:                  # CASE PERIOD
                if not extra_archive: # ONE ARCHIVE REGULAR CASE
                    ax.plot(tgrp,[i]*2, '-' + color1)
                else:              # SUBCASE "EXTRA" : SEVERAL ARCHIVE
                    for tt , cc in zip(Tgrp_plt , Color_grp_plt):
                        ax.plot(tt,[i]*2 , '-' + cc)

        # PART 2 : PLOT ANEX DATADICO
        for idatadico_anex,datadico_anex in enumerate(datadico_anex_list):
            if stat in list(datadico_anex.keys()):
                T = datadico_anex[stat]
                T = [ t for t in T if start <= t <= end  ]
                T = sorted(T)

                if jul_date_plot:
                    T = MJD2dt(T)

                pale_blue_dot , = ax.plot(T,i*np.ones(len(T)), 'o', color='skyblue',label="final SNX")
                legend_list     = [pale_blue_dot]

    #### LEGEND
    if colordico_for_main_datadico:
        import matplotlib.lines as mlines
        for arcnam , col in colordico_for_main_datadico.items():
            legend_list.append(mlines.Line2D([], [], color=col,
                                                     label=arcnam))
            plt.legend(handles=legend_list,loc='upper left')


#    ax.set_yticks(np.arange(0,len(plotstat_lis)-1),plotstat_lis)
    plt.yticks(np.arange(0,len(stats_concat_list)),stats_concat_list)
    fig.autofmt_xdate()
    koef = np.sqrt(2) * 1
    fig.set_size_inches(koef*11.69,koef*len(stats_concat_list) * 0.28) #16.53
    ax.set_ylim([-1 , len(stats_concat_list) + 1])

    return fig


# def timeline_plotter_old_bkp(datadico,start = dt.datetime(1980,1,1) ,
#                      end = dt.datetime(2099,1,1),dots_plot=False,
#                      jul_date_plot=False,datadico_anex_list = [],
#                      use_only_stats_of_main_datadico=False,
#                      colordico_for_main_datadico=None):
#     """
#     A simpler version has been commited to geodezyx toolbox for archive
#     on 20180118 15:59A
#     """
#
#     fig , ax = plt.subplots()
#     ax.xaxis.grid(True)
#     ax.yaxis.grid(True)
#
#     if not use_only_stats_of_main_datadico:
#         stats_concat_list = list(datadico.keys()) + sum([list(e.keys()) for e in datadico_anex_list], [])
#         stats_concat_list = list(reversed(sorted(list(set(stats_concat_list)))))
#     else:
#         stats_concat_list = list(reversed(sorted(list(datadico.keys()))))
#
#     legend_list = [] # must be here before the loop
#     for i,stat in enumerate(stats_concat_list):
#         # PART 1 : PLOT MAIN DATADICO
#         if stat in datadico.keys():
#             T = [ e[-1] for e in datadico[stat] ]
#             T = [ t for t in T if start <= t <= end ]
#             T = sorted(T)
#
#             O = [ e[1] for e in datadico[stat] ]
#
#
#             if dots_plot:
#                 #old old with dots
#                 ax.plot(T,i*np.ones(len(T)), '.')
#             else:
#                 TMJD=dt2MJD(T)
#                 TGrpAll = genefun.consecutive_groupIt(TMJD,True)
#
#                 for tgrp in TGrpAll:
#
#                     color1 = ''
#                     color2 = ''
#                     extra_archive = False
#
#                     # managing colors
#                     if colordico_for_main_datadico:
#                         igrpstart    = TMJD.index(tgrp[0])
#                         igrpend      = TMJD.index(tgrp[1])
#                         Ogrp         = O[igrpstart:igrpend+1]
#                         Tgrp         = TMJD[igrpstart:igrpend+1]
#                         opt_set      = list(set(Ogrp))
#
#                         if len(opt_set) == 1: # Regular case : only one archive
#                             if opt_set[0] in colordico_for_main_datadico:
#                                 color1 = colordico_for_main_datadico[opt_set[0]]
#                         else: # several archives ... so the line has to be splited in segments
#                             extra_archive = True
#                             OSubgrp = genefun.identical_groupIt(Ogrp)
#                             TSubgrp = genefun.sublistsIt(Tgrp , [len(e) for e in OSubgrp])
#
#                             Tgrp_plt      = [ (e[0],e[-1] + 1)    for e in TSubgrp ] # +1 because the end boundary day is not included
#                             Ogrp_plt      = [ list(set(e))[0] for e in OSubgrp     ]
#                             Color_grp_plt = [ colordico_for_main_datadico[e] if e in colordico_for_main_datadico else color2 for e in Ogrp_plt ]
#
#                             print(Ogrp)
#
#                     tgrp = (tgrp[0] , tgrp[1] + 1) # +1 because the end boundary day is not included
#
#                     if not jul_date_plot:
#                         tgrp = MJD2dt(tgrp)
#                         if extra_archive:
#                             Tgrp_plt = [ (MJD2dt(e[0]) , MJD2dt(e[1])) for e in Tgrp_plt ]
#
#
#                     #PLOT part
#                     if tgrp[0] == tgrp[1]: # CASE NO PERIOD, ONLY ONE DAY
#                         ax.plot(tgrp[0],i , '.' + color1)
#                     else:                  # CASE PERIOD
#                         if not extra_archive: # ONE ARCHIVE REGULAR CASE
#                             ax.plot(tgrp,[i]*2, '-' + color1)
#                         else:              # SUBCASE "EXTRA" : SEVERAL ARCHIVE
#                             for tt , cc in zip(Tgrp_plt , Color_grp_plt):
#                                 print(i,stat,tt,cc)
#                                 ax.plot(tt,[i]*2 , '-' + cc)
#
#         # PART 2 : PLOT ANEX DATADICO
#         for idatadico_anex,datadico_anex in enumerate(datadico_anex_list):
#             if stat in list(datadico_anex.keys()):
#                 T = datadico_anex[stat]
#                 T = [ t for t in T if start <= t <= end  ]
#                 T = sorted(T)
#
#                 if jul_date_plot:
#                     T = MJD2dt(T)
#
#                 pale_blue_dot , = ax.plot(T,i*np.ones(len(T)), 'o', color='skyblue',label="final SNX")
#                 legend_list = [pale_blue_dot]
#
#     #### LEGEND
#     if colordico_for_main_datadico:
#         import matplotlib.lines as mlines
#         for arcnam , col in colordico_for_main_datadico.items():
#             legend_list.append(mlines.Line2D([], [], color=col,
#                                                      label=arcnam))
#             plt.legend(handles=legend_list,loc='upper left')
#
#
# #    ax.set_yticks(np.arange(0,len(plotstat_lis)-1),plotstat_lis)
#     plt.yticks(np.arange(0,len(stats_concat_list)),stats_concat_list)
#     fig.autofmt_xdate()
#     koef = np.sqrt(2) * 1
#     fig.set_size_inches(koef*11.69,koef*len(stats_concat_list) * 0.28) #16.53
#     ax.set_ylim([-1 , len(stats_concat_list) + 1])
#
#     return fig


# THIS FIRST SOLUTION PRINT DASHED LINES, NOT GOOD
#                        if len(opt_set) == 1: # Regular case : only one archive
#                            if opt_set[0] in colordico_for_main_datadico:
#                                color1 = colordico_for_main_datadico[opt_set[0]]
#                        else:
#                            dashed_line = True
#                            if opt_set[0] in colordico_for_main_datadico:
#                                color1 = colordico_for_main_datadico[opt_set[0]]
#                            if opt_set[1] in colordico_for_main_datadico:
#
# THIS SECOND SOLUTION PRINT X, NOT GOOD
#                            o_major = genefun.most_common(Ogrp)
#                            if o_major in colordico_for_main_datadico:
#                                color1 = colordico_for_main_datadico[o_major]
#
#                            Textra , Color_extra = [],[]
#                            Tgrp      = T[igrpstart:igrpend+1]
#
#                            for oo , tt in zip(Ogrp,Tgrp):
#                                if oo == o_major:
#                                    continue
#                                else:
#                                    if oo in colordico_for_main_datadico:
#                                        color2 = colordico_for_main_datadico[oo]
#                                    Textra.append(tt)
#                                    Color_extra.append(color2)


def listing_gins_timeline(path,stat_strt,stat_end,date_strt,date_end,suffix_regex=''):
    """ find all gins listings in a folder and his subfolders
        and plot timeline of the avaiable listings

        stat_strt,stat_end,date_strt,stat_end : where to find
        in the name the statname and the date"""

    fig , ax = plt.subplots()

    aaa = list(os.walk(path))
    wholefilelist = []
    for tup in aaa:
        wholefilelist = wholefilelist + tup[-1]

    wholefilelist = list(set(wholefilelist))
    lifilelist = [fil for fil in wholefilelist if re.search(suffix_regex +'.*\gins', fil)]

    statname_lis = sorted(list(set([li[stat_strt:stat_end] for li in lifilelist])))

    list(set(statname_lis))

    datadico = dict()

    for stat in statname_lis:
        datadico[stat] = []

    for li in lifilelist:
        try:
            tup = (li , jjulCNES2dt(li[date_strt:date_end] ))
            datadico[li[stat_strt:stat_end]].append(tup)
        except:
            print('error with : ', li)
    plotstat_lis = []
    for i,stat in enumerate(sorted(datadico.keys())):
        print(i , stat)
        T = [e[-1] for e in datadico[stat]]
        plotstat_lis.append(stat)
        ax.plot(T,i*np.ones(len(T)), '.')

    i_list = np.arange(0,len(plotstat_lis))

    print(i_list)
    plt.yticks(i_list,plotstat_lis)

    return fig

def error_ellipse(xm,ym,sigx,sigy,sigxy,nsig=1,ne = 100,scale=1):
    """
    from matlab fct
    http://kom.aau.dk/~borre/matlab/geodesy/errell.m
    It works but don't ask why ...

    (X,Y) orientation convention is inverted  => (Y,X) ...
    so in a practical way you must invert X ,Y
    (it is not important for the axis but it is for the orientation)
    AND
    sigx,sigy,sigxy must be first normalized with the fuv

    sigx, sigy, sigxy :
        so as we can generate a covariance matrix
        cov = np.array([[sigx ** 2,sigxy],[sigxy,sigy ** 2]])

    ne :
        nb of segements for the ellipse

    RETURNS :
        xe,ye,dx2,dy2


    DEBUG:

    si on a
    xe1,ye1,_,_ = geok.error_ellipse(pxp[0],pxp[1], sigxB , sigyB , sigxyB, scale= 10000)
    xe2,ye2,_,_ = geok.error_ellipse(pxp[0],pxp[1], sigyB , sigxB , sigxyB, scale= 10000)
    et PAS les - à D et dxy0 => on a 2 ellipses differentes

    si on a
    xe1,ye1,_,_ = geok.error_ellipse(pxp[0],pxp[1], sigxB , sigyB , sigxyB, scale= 10000)
    xe2,ye2,_,_ = geok.error_ellipse(pxp[0],pxp[1], sigyB , sigxB , sigxyB, scale= 10000)
    et AVEC les - à D et dxy0 => on a 2 ellipses differentes
    au moins une ellipse coincide avec celle de Ghiliani

    A investiguer, en attendant, à éviter
    """

    cov = np.array([[sigx ** 2,sigxy],[sigxy,sigy ** 2]])
    V,D = np.linalg.eig(cov)

    D = -D.T

    std1 = np.sqrt(V[0])
    std2 = np.sqrt(V[1])

    if std1 < std2:
        z1 = np.linspace(-std1,std1,ne)
        z2 = np.sqrt(V[1] * (np.ones(ne) - (z1 / std1) **2))
        dxy = np.dot(D,np.vstack((np.hstack((z1,np.flipud(z1))),np.hstack((z2,-z2)))))
    else:
        z2 = np.linspace(-std2,std2,ne)
        z1 = np.sqrt(V[0] * (np.ones(ne) - (z2 / std2) **2))
        dxy = np.dot(D,np.vstack((np.hstack((z1,-z1)),np.hstack((z2,np.flipud(z2))))))

    dx =  -dxy[0,:]
    dy =  dxy[1,:]

    dx2 =  scale*nsig*dx
    dy2 =  scale*nsig*dy
    xe = xm * np.ones(2*ne) + dx2
    ye = ym * np.ones(2*ne) + dy2

    return xe,ye,dx2,dy2

def error_ellipse_parameters(qxx,qyy,qxy,fuv,out_t=False):
    """
    INPUT :
        qxx,qyy,qxy  : factors as in the varcovar matrix (no normalisation with the fuv or other)
        fuv
    OUTPUT :
        Su/a Sb/b    : semimajor and semiminor axis
        t            : angle that the u/a axis makes with the y axis in clockwise direction
        OR
        phi          : angle that the u/a axis makes with the x axis in trigo direction

    (this one is the natural way, perfect to test those of the Strang & Borre)
    """
    S0   = np.sqrt(fuv)
    t    = 0.5 * np.arctan2( 2*qxy , (qyy - qxx))
    tdeg = np.rad2deg(t)

    quu = qxx * np.sin(t)**2 + 2*qxy * np.cos(t) * np.sin(t) + qyy * np.cos(t)**2
    qvv = qxx * np.cos(t)**2 - 2*qxy * np.cos(t) * np.sin(t) + qyy * np.sin(t)**2
    Su = S0 * np.sqrt(quu)
    Sv = S0 * np.sqrt(qvv)

    a   = Su
    b   = Sv
    phi = 90 - tdeg
    if out_t:
        return a,b,tdeg
    else:
        return a,b,phi


def error_ellipse_parameters_2(sigx,sigy,sigxy,out_deg=True):
    """
    ref : Strang & Borre p 337
    (X,Y) orientation convention is inverted  => (Y,X) ...
    so in a practical way you must invert X ,Y
    (it is not important for the axis but it is for the orientation)
    AND
    sigx,sigy,sigxy must be normalized with the fuv
    """

    cov = np.array([[sigx ** 2,sigxy],[sigxy,sigy ** 2]])

    sigx2 = sigx ** 2
    sigy2 = sigy ** 2
    sigxy2 = sigxy ** 2

    cov = np.array([[sigx ** 2,sigxy],[sigxy,sigy ** 2]])

    V,D   = np.linalg.eig(cov)
    V2,D2 = np.linalg.eig(np.linalg.inv(cov))

    Dsqrt = np.sqrt(V)

    l1 = 0.5 * (sigx2 + sigy2 + np.sqrt((sigx2 + sigy2)**2 - 4*(sigx2 * sigy2 - sigxy2)))
    l2 = 0.5 * (sigx2 + sigy2 - np.sqrt((sigx2 + sigy2)**2 - 4*(sigx2 * sigy2 - sigxy2)))

    a = np.sqrt(l1)
    b = np.sqrt(l2)
    phi = ( 0.5 * np.arctan2(2*sigxy,sigx2-sigy2))

    if out_deg:
        phi = np.rad2deg(phi)
    return a,b,phi


def ellipse_get_coords(a=0.0, b=0.0, x=0.0, y=0.0, angle=0.0, k=2 ,
                       out_separate_X_Y = True , trigo = True):
    """ Draws an ellipse using (360*k + 1) discrete points; based on pseudo code
    given at http://en.wikipedia.org/wiki/Ellipse
    k = 1 means 361 points (degree by degree)
    a = major axis distance,
    b = minor axis distance,
    x = offset along the x-axis
    y = offset along the y-axis
    angle = trigo/clockwise rotation [in degrees] of the ellipse;
        * angle=0  : the ellipse is aligned with the positive x-axis
        * angle=30 : rotated 30 degrees trigo/clockwise from positive x-axis

    trigo sense is the standard convention

    NB : clockwise is the internal convention, but we prefer trigo
         convention for the Ghiliani ellipses made by error_ellipse_parameters

    source : scipy-central.org/item/23/1/plot-an-ellipse
    """

    if trigo:
        angle = - angle



    pts = np.zeros((360*k+1, 2))

    beta = -angle * np.pi/180.0
    sin_beta = np.sin(beta)
    cos_beta = np.cos(beta)
    alpha = np.radians(np.r_[0.:360.:1j*(360*k+1)])

    sin_alpha = np.sin(alpha)
    cos_alpha = np.cos(alpha)

    pts[:, 0] = x + (a * cos_alpha * cos_beta - b * sin_alpha * sin_beta)
    pts[:, 1] = y + (a * cos_alpha * sin_beta + b * sin_alpha * cos_beta)

    if not out_separate_X_Y:
        return pts
    else:
        return pts[:, 0] , pts[:, 1]


def ellipse_center(a):
    """
    http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
    """
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])


def ellipse_angle_of_rotation( a , outdeg = True ):
    """
    core fct for ellipse_fit
    http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
    """
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    if not outdeg:
        return 0.5*np.arctan(2*b/(a-c))
    else:
        return np.rad2deg(0.5*np.arctan(2*b/(a-c)))


def ellipse_axis_length( a ):
    """
    core fct for ellipse_fit
    http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
    """
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(up/down1)
    res2=np.sqrt(up/down2)
    return np.array([res1, res2])

def fitEllipse_core(x,y):
    """
    core fct for ellipse_fit
    http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
    """
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  np.linalg.eig(np.dot(np.linalg.inv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    return a

def ellipse_fit(x,y):
    """
    find the parameters
    a,b,phi,x0,y0 of an ellipse

    from
    http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
    """
    tmp   = fitEllipse_core(x,y)
    a,b   = ellipse_axis_length( tmp )
    phi   = ellipse_angle_of_rotation(tmp)
    x0,y0 = ellipse_center(tmp)
    return a,b,phi,x0,y0


#  _                    _      _____                  _    _ _   _ _
# | |                  | |    / ____|                | |  | | | (_) |
# | |     ___  __ _ ___| |_  | (___   __ _ _ __ ___  | |  | | |_ _| |___
# | |    / _ \/ _` / __| __|  \___ \ / _` | '__/ __| | |  | | __| | / __|
# | |___|  __/ (_| \__ \ |_   ____) | (_| | |  \__ \ | |__| | |_| | \__ \
# |______\___|\__,_|___/\__| |_____/ \__, |_|  |___/  \____/ \__|_|_|___/
#                                       | |
#                                       |_|

def get_accur_coeff(i):
    accur_coeff_mat = np.array([[0,0,0,-1/2.,0,1/2.,0,0,0],
    [0,0,1/12. ,-2/3.,0,2/3.,-1/12.,0,0],
    [0,-1/60.,3/20.,-3/4.,0,3/4.,-3/20.,1/60.,0],
    [1/280.,-4/105.,1/5.,-4/5.,0,4/5.,-1/5.,4/105.,-1/280.]])
    if i > 3:
        return accur_coeff_mat[-1]
    else:
        return accur_coeff_mat[i]

def partial_derive_old(f,var_in,var_out=0,kwargs_f={},args_f=[],h=0):
    ''' var_in :
            detrivation with respect to this variable
            can be a int (starts with 0) or a string descirbing the name of
            the var in f
        var_out :
            the output of f which needs to be considerated as the output
            ** must be a int **
        args_f & kwargs_f :
            tuple/list & dict describing the arguments of f
        h :
            derivation step, if h == 0 give x * sqrt(epsilon)
            (source : http://en.wikipedia.org/wiki/Numerical_differentiation) '''

    # tuple => list pour plus d'aisance
    args_f = list(args_f)

    # operational arguments
    args_f_m = list(args_f)
    args_f_p = list(args_f)
    kwargs_f_m = dict(kwargs_f)
    kwargs_f_p = dict(kwargs_f)

    args_name_list = list(f.__code__.co_varnames)
    # var in is a int
    if type(var_in) is int:
        var_ind = var_in
        var_name = args_name_list[var_ind]
    # var in is a str
    else:
        var_name = var_in
        try:
            var_ind = args_name_list.index(var_name)
        except ValueError:
            print(args_name_list)
            raise Exception('wrong var_in name (not in args name list)')

#    if var_ind < len(args_f):
#        x = args_f[var_ind]
#        if h == 0:
#            h = x * np.sqrt(np.finfo(float).eps)
#            if h == 0:
#                print 'WARN : h == 0 ! setting @ 10**-6 '
#                h = 10**-6
#        args_f_m[var_ind] = x - h
#        args_f_p[var_ind] = x + h
#    else:
#        x = kwargs_f[var_name]
#        if h == 0:
#            h = x * np.sqrt(np.finfo(float).eps)
#            if h == 0:
#                print 'WARN : h == 0 ! setting @ 10**-6 '
#                h = 10**-6
#        kwargs_f_m[var_name] = x - h
#        kwargs_f_p[var_name] = x + h


    if var_ind < len(args_f):
        x = args_f[var_ind]
    else:
        x = kwargs_f[var_name]

    if h == 0:
        h = x * np.sqrt(np.finfo(float).eps)
        if h == 0:
            print('WARN : h == 0 ! setting @ 10**-6 ')
            h = 10**-6

    if var_ind < len(args_f):
        args_f_m[var_ind] = x - h
        args_f_p[var_ind] = x + h
    else:
        kwargs_f_m[var_name] = x - h
        kwargs_f_p[var_name] = x + h

    m = f(*args_f_m,**kwargs_f_m)
    p = f(*args_f_p,**kwargs_f_p)

    if genefun.is_iterable(m):
        m = m[var_out]
        p = p[var_out]
#    print p,m,h,x
#    print h
#    print h == 0
    dout = (p - m) / (2. * float(h))

    return dout


def partial_derive(f,var_in,var_out=0,kwargs_f={},args_f=[],h=0,accur=-1):
    ''' var_in :
            detrivation with respect to this variable
            can be a int (starts with 0) or a string describing the name of
            the var in f
        var_out :
            the output of f which needs to be considerated as the output
            ** must be a int **

        args_f & kwargs_f :
            tuple/list & dict describing the arguments of f
        h :
            derivation step, if h == 0 give x * sqrt(epsilon)
            (source : http://en.wikipedia.org/wiki/Numerical_differentiation) '''

    # tuple => list pour plus d'aisance
    args_f = list(args_f)

    # operational arguments
    args_f_i = list(args_f)
    kwargs_f_i = dict(kwargs_f)

    args_name_list = list(f.__code__.co_varnames)
    # var in is a int
    if type(var_in) is int:
        var_ind = var_in
        var_name = args_name_list[var_ind]
    # var in is a str
    else:
        var_name = var_in
        try:
            var_ind = args_name_list.index(var_name)
        except ValueError:
            print(args_name_list)
            raise Exception('wrong var_in name (not in args name list)')

    if var_ind < len(args_f):
        x = args_f[var_ind]
    else:
        x = kwargs_f[var_name]

    if h == 0:
        h = x * np.sqrt(np.finfo(float).eps)
        if h == 0:
            print('WARN : h == 0 ! setting @ 10**-6 ')
            h = 10**-6

    res_stk = []
    accur_coeff = get_accur_coeff(accur)
    for i,k in enumerate(accur_coeff):
        if k == 0:
            res_stk.append(0.)
        else:
            if var_ind < len(args_f):
                args_f_i[var_ind] = x + h * (i-4)
            else:
                kwargs_f_i[var_name]   = x + h * (i-4)

            res = f(*args_f_i,**kwargs_f_i)

            if genefun.is_iterable(res):
                res = res[var_out]

            res_stk.append(res)

    dout = np.dot(np.array(res_stk) , accur_coeff) / h

    return dout


def jacobian_line(f,var_in_list,var_out=0,kwargs_f={},args_f=[],h=0,aray=True):
    """ same argument as ``partial_derive`` except var_in becomes var_in_list :
    it's an iterable of all variables the derivation must be performed
    """
    line_out =  []
    for var_in in var_in_list:
        line_out.append(partial_derive(f,var_in,var_out,kwargs_f,args_f,h))
    if aray:
        return np.array(line_out)
    else:
        return line_out


#def jacobian_old_bkp(f,var_in_list,var_out,kwargs_f_list=[],args_f_list=[],h=10**-6):
#    if args_f_list == [] and kwargs_f_list == []:
#        raise Exception('jacobian : args_f_list == [] and kwargs_f_list ==  []')
#    elif args_f_list == [] and kwargs_f_list != []:
#        args_f_list =  [ [] * len(kwargs_f_list)]
#    elif args_f_list != [] and kwargs_f_list == []:
#        kwargs_f_list =  [ [] * len(args_f_list)]
#
#    if len(args_f_list) != len(kwargs_f_list):
#        raise Exception("Jacobian : len(args_f_list) != len(kwargs_f_list)")
#    jacob_temp = []
#    for kwargs_f , args_f in zip(kwargs_f_list,args_f_list):
#        line = jacobian_line(f,var_in_list,var_out,kwargs_f,args_f,h)
#        jacob_temp.append(line)
#
#    return np.vstack(*jacob_temp)

def kwargs_for_jacobian(kwdic_generik,kwdic_variables):
    """ Building a list of kwargs for the jacobian function
        kwdic_generik : parameters which not gonna change
        kwdic_variable : parameters which gonna change, so must be associated
        with iterable
        """
    keys_list = list(kwdic_variables.keys())

    for k,v in kwdic_variables.items():
        if not genefun.is_iterable(v):
            print('WARN : key',k,'val',v,'is not iterable !!!')

    values_combined = itertools.product(*list(kwdic_variables.values()))

    kwdic_list_out =[]
    for values in values_combined:
        kwdic_out = dict(kwdic_generik)
        for k,v in zip(keys_list,values):
            kwdic_out[k] = v
        kwdic_list_out.append(kwdic_out)

    return kwdic_list_out

def jacobian(f,var_in_list,var_out,kwargs_f_list=[],h=10**-6,nproc=4):
    """ il n'y a que les kwargs qui sont gérés """
    args_f_list = []
    for i in range(len(kwargs_f_list)):
        args_f_list.append([])


    if len(args_f_list) != len(kwargs_f_list):
        print(len(args_f_list) , len(kwargs_f_list))
        raise Exception("Jacobian : len(args_f_list) != len(kwargs_f_list)")
    jacob_temp = []
    args_list = []
    for kwargs_f , args_f in zip(kwargs_f_list,args_f_list):
        args_list.append((f,var_in_list,var_out,kwargs_f,args_f,h))

    pool = mp.Pool(processes=nproc)
    results = [pool.apply(jacobian_line, args=x) for x in args_list]
    pool.close()
    jacob_temp = results

    return np.vstack(jacob_temp)


def nan_cleaner(Ain,Bin):
    """
    remove A & B of their respective NaN

    Args :
        Ain , Bin : lists/arrays

    Returns:
        Aout , Bout : A & B withour NaN
    """
    Ain = np.array(Ain)
    Bin = np.array(Bin)

    notnanA = np.logical_not(np.isnan(Ain))
    notnanB = np.logical_not(np.isnan(Bin))

    notnan = notnanA * notnanB

    Aout = Ain[notnan]
    Bout = Bin[notnan]

    return Aout , Bout


def clean_nan(A,L):
    """
    DISCONTINUED

    A est un array bi dimentionnel
    L est un array mono dimentionnel
    renvoie un A et un L nettoyé mutuellement
    de leurs NaN respectifs
    return np.sqrt(a + a.T - np.diag(a.diagonal()))
    """
    if A.ndim != 2  or L.ndim != 1:
        raise Exception("ERREUR : revoir les dimensions de A (2) et L (1)")

    if A.shape[0] != L.shape[0]:
        raise Exception("ERREUR : A et L ont des longeurs differentes")

    # 11) Detecter les nan de L
    indLnonan = np.where(~np.isnan(L))
    nbnanL = np.isnan(L).sum()
    nbnanAtot = np.isnan(A).sum()
    nbnanA = np.isnan(A).any(1).sum()

    # 12) Virer les nan de L dans A
    A2 = A[indLnonan[0],:]
    # 13) Virer les nan de L dans L
    L2 = L[indLnonan]
    # 21) Detecter les nan de A
    # (même si cetains ont potentiellement déjà été supprimé)
    boolA2arenan = np.isnan(A2)
    nbnanA2 = boolA2arenan.any(1).sum()
    # 22) Virer les nan de A dans A
    Aout = A2[~boolA2arenan.any(1)]

    # 23) Virer les nan de A dans L
    Lout = L2[~boolA2arenan.any(1)]

    if  np.isnan(Lout).sum() != 0 or  np.isnan(Aout).sum() != 0 :
        print("pour info, isnan(Lout).sum(), isnan(Aout).sum()")
        print(np.isnan(Lout).sum(), np.isnan(Aout).sum())
        raise Exception("ERREUR : A_out et L_out contiennent des NaN")

    if Aout.shape[0] != Lout.shape[0]:
        raise Exception("ERREUR : A_out et L_out ont des longeurs differentes")

    print("%i NaN dans le V" %nbnanL)
    print("%i NaN dans la M vert., APRES suppr. des NaN du V" %nbnanA2)
    print("%i lignes supprimées dans le V et la M (en théorie)" %(nbnanA2 + nbnanL))
    print("")
    print("pour info :")
    print("%i NaN dans la M vert., AVANT suppr. des NaN du V" %nbnanA)
    print("%i NaN dans la M AU TOTAL (vert. et horiz.)" %nbnanAtot)

    return Aout , Lout

#def weight_mat_OLD(Sinp,Ninp,fuvinp=1):
#    """ Sinp : liste des Sigmas sig = sqrt(var)
#        Ninp : liste de la taille de chaque blocs (obs)
#        fuvinp = 1 : facteur unitaire de variance
#        inspiré de mat_poids , fct écrite dans la lib resolution de GPShApy
#
#        retourne :
#        K : matrice de var-covar
#        Q : matrice des cofacteurs
#        P : matrice des poids inv(Q)
#    """
#
#    if len(Sinp) != len(Ninp):
#        raise Exception("S et N de taille differente")
#
#    Ktemp = []
#
#    for i in range(len(Sinp)):
#        try:
#            Ktemp.append(np.eye(Ninp[i]) * Sinp[i]**2)
#        except DeprecationWarning:
#            print "weight_mat : Are you sure you don't invert Sinp <> Ninp ?"
#    K = scipy.linalg.block_diag(*Ktemp)
#    Q = (1/fuvinp) * K
#    # normalement P = scipy.linalg.inv(Q)
#    # mais pour que ce soit + rapide
#    Qdiag = np.diagonal(Q)
#    Pdiag = 1/Qdiag
#    P = scipy.linalg.block_diag(*Pdiag)
#
#    return K , Q , P

def weight_mat(Sinp,Ninp=[],fuvinp=1,sparsediag=False):
    """
    Args :
        Sinp : liste des Sigmas sig = sqrt(var)
        Ninp : liste de la taille de chaque blocs (obs)
        fuvinp = 1 : facteur unitaire de variance
        inspiré de mat_poids , fct écrite dans la lib resolution de GPShApy

    Returns :
        K : matrice de var-covar
        Q : matrice des cofacteurs
        P : matrice des poids inv(Q)
    """

    if Ninp == []:
        Ninp = [1] * len(Sinp)

    Sinp = np.array(Sinp)
    if np.any(Sinp == 0):
        print("WARN : some sigma inputs are == 0 !!!")

    if len(Sinp) != len(Ninp):
        raise Exception("S et N de taille differente")

    Ktemp = []
    for i in range(len(Sinp)):
        try:
            Ktemp = Ktemp + Ninp[i] * [Sinp[i]**2]
        except DeprecationWarning:
            print("weight_mat : Are you sure you don't invert Sinp <> Ninp ?")
    Ktemp = np.array(Ktemp).astype(np.float64)
    Qtemp = (1. /fuvinp) * Ktemp
    Ptemp = 1. / Qtemp

    if sparsediag:
        K = scipy.sparse.diags(Ktemp,0)
        Q = scipy.sparse.diags(Qtemp,0)
        P = scipy.sparse.diags(Ptemp,0)
    else:
        K = np.diag(Ktemp)
        Q = np.diag(Qtemp,0)
        P = np.diag(Ptemp,0)

    #    #K = scipy.linalg.block_diag(*Ktemp)
    #    Q = (1/fuvinp) * K
    #    # normalement P = scipy.linalg.inv(Q)
    #    # mais pour que ce soit + rapide
    #    Qdiag = np.diagonal(Q)
    #    Pdiag = 1/Qdiag
    #    P = scipy.linalg.block_diag(*Pdiag)

    return K , Q , P



def weight_mat_simple(Pinp,Ninp=[],sparsediag=False,
                      return_digaonal_only=False):
    """
    Simple version of weight_mat : takes directly the weights (Pinp)
    and the size for each weigths blocks (Ninp)
    
    Pinp and Ninp have to have the same length
    
    Args :
        Pinp : list of weigths
        Ninp : list of the size of each block (obs number)
        fuvinp = 1 : facteur unitaire de variance

    Returns :
        P : weigth matrix
    """

    if Ninp == []:
        Ninp = [1] * len(Pinp)

    if len(Pinp) != len(Ninp):
        raise Exception("different sizes for Pinp and Ninp !!!")

    Ptemp = []
    for i in range(len(Pinp)):
        try:
            Ptemp = Ptemp + Ninp[i] * [Pinp[i]]
        except DeprecationWarning:
            print("weight_mat_simple : Are you sure you don't invert Pinp <> Ninp ?")

    if return_digaonal_only:
        P = np.array(Ptemp)
    
    else:
        if sparsediag:
            P = scipy.sparse.diags(Ptemp,0)
        else:
            P = np.diag(Ptemp,0)

    return P

def fuv_calc_OLD(V,A):
    fuv = np.dot(V.T,V) / np.abs(A.shape[0] - A.shape[1])
    return fuv

def fuv_calc_OLD2(V,A,P=None):
    if P is None:
        P = np.eye(len(V))

    P = np.matrix(P)
    V = np.matrix(V)
    A = np.matrix(A)

    if V.shape[0] == 1:
        V = V.T

    fuv = (V.T * P * V) / np.abs(A.shape[0] - A.shape[1])
    return fuv[0,0]


def fuv_calc(V,A,P=1,normafuv=1):
    """
    Args :
        V : residuals vector

        A : Jacobian matrix

        P : weight matrix

        Can manage both standard arrays and sparse array

    Returns :
        fuv : Facteur unitaire de variance (unitary variance factor)

    Notes :
        le fuv dépend de la martice de poids
        mais les sigmas non
        ex :
        poids de 10**-6
        fuv    :  439828.260843
        sigma  :  [ 5.21009306  5.09591568  0.04098106]
        poids de 1
        fuv    :  4.39828260843e-07
        sigma  :  [ 5.21009306  5.09591568  0.04098106]
    """
    V = np.squeeze(np.array(V))

    if not genefun.is_iterable(P):
        P = np.ones(len(V)) * P
    elif scipy.sparse.issparse(P):
        P = P.diagonal()
    else:
        print("DEPRECIATION : modification done in fuv_calc, P should be given as Matrix-shaped now")
        P = P.diagonal()

    P = P * (1 / normafuv)

    VPV = np.column_stack((V,P,V))
    numera = np.sum(np.product(VPV,1))

    # nummera is just an adaptation of VT * P * V
    # EDIT 1806 : Je veux bien ... mais c'est une adaptation pourrie !
    # ou alors il faut bien s'assurer que l'on a extrait la diagonale de P
    
    
    if A.ndim == 1:
      A = np.expand_dims(A,0)

    fuv = numera / np.abs(A.shape[0] - A.shape[1])
    return fuv


def sigmas_formal_calc(N,V,A,fuv=None,P=None):
    if not fuv:
        fuv      = fuv_calc(V,A,P)
        
    varcovar = (fuv) * scipy.linalg.inv(N)
    sigmas   = np.sqrt(varcovar.diagonal())
    return sigmas , varcovar
  

def smart_i_giver(subgrp_len_list,i_in_sublis,sublis_id,
                  advanced=False,sublis_id_list=[]):
    """
    eg
    subgrp_len_list = [4201, 4186, 4157, 4041, 4058, 4204, 4204, 4204, 4204]
    i_in_sublis = 2
    sublis_id = 3
    return 4201 + 4186 + 4157 + 2

    advanced = True:
    the sublis_id is not a int but and generic identifier ( str ,int , set ... )
    present in sublis_id_list
    else
    must be an int
    """

    if sublis_id_list == [] and not type(sublis_id) is int:
        print("ERR : smart_i_giver")
        return None

    if advanced:
        sublis_id = sublis_id_list.index(sublis_id)

    if sublis_id == 0:
        decaleur = 0
    else:
        decaleur = -1

    before = np.sum(subgrp_len_list[:sublis_id])
    return int(before + decaleur + i_in_sublis)


def constraint_improve_N(N,C,trans=False,outsparsetype = 'csc'):
    """ give N normal matrix and C constraints matrix
        returns N compined with C
        trans is a (dirty) way to transpose C
        if made in wrong shape

        convention Ghilani 2011 p424 :

        N    C.T
        C     0
    """

    if trans:
        C = C.T

    O = np.zeros((C.shape[0],C.shape[0]))

    if scipy.sparse.issparse(C):
        L1   = scipy.sparse.hstack((N,C.T))
        L2   = scipy.sparse.hstack((C,O))
        Nout = scipy.sparse.vstack((L1,L2))

        if outsparsetype == 'csc':
            Nout = scipy.sparse.csc_matrix(Nout)
        elif outsparsetype == 'csr':
            Nout = scipy.sparse.csr_matrix(Nout)

    else:
        L1 = np.hstack((N,C.T))
        L2 = np.hstack((C,O))
        Nout = np.vstack((L1,L2))

    return Nout

def triangle_arr2vect(triarrin,k=1):
    ind = np.triu_indices_from(triarrin,k)
    return triarrin[ind]


#   ____              _                  _
#  / __ \            | |                (_)
# | |  | |_   _  __ _| |_ ___ _ __ _ __  _  ___  _ __  ___
# | |  | | | | |/ _` | __/ _ \ '__| '_ \| |/ _ \| '_ \/ __|
# | |__| | |_| | (_| | ||  __/ |  | | | | | (_) | | | \__ \
#  \___\_\\__,_|\__,_|\__\___|_|  |_| |_|_|\___/|_| |_|___/
#


def quaternion(r,p,h,angtype='rad'):
    """
    Args:
        roll (xaxis) pitch (yaxis) heading (zaxis) angle

        angtype : 'rad' or 'deg'

    Returns:
        the associate quaternion

    Note:
        the input args are in XYZ but
        the rotation is in the order ZYX

        Need of the cgkit quaternion module
        http://cgkit.sourceforge.net/doc2/quat.html

    """
    if angtype == 'deg':
        r,p,h = np.deg2rad(r),np.deg2rad(p),np.deg2rad(h)
    M = cgt.mat3().fromEulerZYX(r,p,h)
    q = cgt.quat().fromMat(M)
    return q

def quatern_multi(rlist,plist,hlist,angtype='rad'):
    """
    wrapper of quaternion, with lists/arrays of angles
    """

    qlisout = []
    for r,p,h in zip(rlist,plist,hlist):
        q = quaternion(r,p,h,angtype=angtype)
        qlisout.append(q)
    qarrout = np.array(qlisout)
    return qarrout

def rotate_quat2(ptin,quatin):
    """
    DISCONTINUED
    plus long que rotate_quat :
    benchmarking 5.53084397316 vs 0.703444957733 pour 100000
    """
    return np.array(quatin.rotateVec(ptin))

def rotate_quat(ptin,quatin):

    pt1 = np.concatenate(([0],np.array(ptin)))
    pt2 = cgt.quat(*pt1)
    qtmp = (quatin) * pt2 * (quatin.inverse())
    return np.array([qtmp.x, qtmp.y , qtmp.z])

def rotate_quat_multi(pts_list_in, quats_list_in):
    """
    Multi version of rotate_quat

    Args:
        pts_list_in : a list of Points object
        quats_list_in : a list of quaternions

    Returns:
        a list of list of rotated points
    """

    pts_lis_out = []
    for pt in pts_list_in:
        temp_pt_lis = []
        for q in quats_list_in:
            ptmp = rotate_quat(pt,q)
            temp_pt_lis.append(ptmp)
        pts_lis_out.append(temp_pt_lis)
    pts_arr_out = np.array(pts_lis_out)
    return pts_arr_out

def interp_quat(Tlis , Quatlis , t):
    """
    Interpolate quaterinon

    Args:
        Tlis (float) : list of time, in POSIX time

        Quatlis : list of quaternion (produced by quaterion fct)

        t : the instant for the interpolation

    Returns:
        an interpoled Quaternion


    """
    i1 , i2 = genefun.find_interval_bound(Tlis,t)
    q1 , q2 = Quatlis[i1] , Quatlis[i2]
    t1 , t2 = Tlis[i1] , Tlis[i2]
    tt = (t - t1) / (t2 - t1)
    qq = cgt.slerp(tt,q1,q2)
    return qq

def interp_quat_multi(Tlis , Quatlis , tlis):
    outlis = []
    for t in tlis:
        q = interp_quat(Tlis , Quatlis , t)
        outlis.append(q)
    return np.array(outlis)


def randomwalk_normal(N=100, d=2 , moy = 0 , sigma = 1):
    """
    d = dimension
    """
    return np.cumsum(moy + np.random.randn(N,d) * sigma)

def randomwalk_uniform(N=100, d=2 , bound = 0.5):
    """
    d = dimension
    bound = contraint of the random walk
    """
    return np.cumsum(np.random.uniform(-bound,bound,(N,d)))


def circle_draw(xc,yc,R,N):
    theta = np.linspace(0,2 * np.pi,N)
    X = np.cos(theta) * R + xc
    Y = np.sin(theta) * R + yc
    return X,Y


def random_walk_in_a_circle(x0 , y0 , xc , yc ,
                            R , N , step_size ,  param = 1 ,
                            polar = True , uniform_or_normal = 'n',
                            rand_seed = -1):

    """
    random : normal ou uniform
    coords : polaire ou cartésien

    param est un paramètre très versatile pour controler le random :
    plage pour le uniform
    sigma pour le normal

    on recommande plutot le polaire normal : on a un pas constant & une derive réaliste sur le cap

    Returns :
        X,Y, Xcircle , Ycircle

    Exemple :

    for un in ('u','n'):
        for pol in range(2):
            X,Y , Xcircle , Ycircle = random_walk_in_a_circle(10,10,0,0,50,10000,polar = pol,uniform_or_normal=un)

            plt.figure()
            plt.plot(Xcircle,Ycircle)
            plt.plot(X,Y)
            plt.axis('equal')
            plt.suptitle(un + str(pol))
    """

    X = [x0]
    Y = [y0]

    if rand_seed > -1:
        RAND = np.random.RandomState(rand_seed)
    else:
        RAND = np.random.RandomState(np.random.randint(10**6))

    Xcircle,Ycircle = circle_draw(xc,yc,R,500)

    for i in range(N-1):
        D = R+1
        iwhil = 0
        while D > R:
            iwhil += 1
            if iwhil > 500:
                print('WARN : infinite loop in random_walk_in_a_circle ...' , iwhil)
            if polar:
                if uniform_or_normal == 'u':
                    dalpha = RAND.uniform(-param,param) * 2 * np.pi
                else:
                    dalpha = RAND.normal(0,param)       * 2 * np.pi
                drho = step_size
                dx = drho * np.cos(dalpha)
                dy = drho * np.sin(dalpha)
                #print dx,dy
            else:
                if uniform_or_normal == 'u':
                    dx = np.random.uniform(-param,param)
                    dy = np.random.uniform(-param,param)
                else:
                    dx = np.random.normal(0,param)
                    dy = np.random.normal(0,param)
                print(dx , dy)
            xtemp = X[-1] + dx
            ytemp = Y[-1] + dy
            D = np.sqrt((xtemp - xc)**2 + (ytemp - yc)**2)

        X.append(xtemp)
        Y.append(ytemp)

    X = np.array(X)
    Y = np.array(Y)

    return X,Y, Xcircle , Ycircle



def randn_bool(N,true_ratio = 0.5,RandGene = None):
    if RandGene is None:
        RandGene = np.random.RandomState()
    if type(RandGene) is int:
        RandGene = np.random.RandomState(RandGene)
    try:
        randlis = RandGene.uniform(size=N)
    except AttributeError:
        "ERR : AttributeError : RandGene  may be an int32/int64, but it's not an authentic int as required ..."
    boolout_lis = []
    for r in randlis:
        if r < true_ratio:
            boolout_lis.append(True)
        else:
            boolout_lis.append(False)

    return boolout_lis


def points_circle_border(Npts,r,r_sigma,az_type_normal=True,
                         main_dir=3.14159,dir_range=3.14159,seed=None):
    if not seed:
        seed = np.random.randint(10000)

    S = np.random.RandomState(seed)

    if not az_type_normal:
        Az = S.rand(Npts)  * 2 * np.pi
    else:
        Az = S.randn(Npts) * dir_range + main_dir


    R = np.array(Npts * [r]) - np.abs(S.randn(Npts) * r_sigma)

    X,Y = polar2cartesian(R,Az,'rad')

    return X , Y


def estimated_autocorrelation(x):
    """
    http://stackoverflow.com/q/14297012/190597
    http://en.wikipedia.org/wiki/Autocorrelation#Estimation
    """
    n = len(x)
    variance = x.var()
    x = x-x.mean()
    r = np.correlate(x, x, mode = 'full')[-n:]
    assert np.allclose(r, np.array([(x[:n-k]*x[-(n-k):]).sum() for k in range(n)]))
    result = r/(variance*(np.arange(n, 0, -1)))
    return result

def savage_buford_formula(Vs,X,d):
    """
    X : distance à la faille , un iterable pour toutes le profil,
    un nombre pour la longeur max

    d : profondeur de la faille
    retourne X , et Vdeform(X)

    X et d doivent être dans la même unité, Vs pas forcément
    """

    if not genefun.is_iterable(X):
        X = np.arange(-X,X,1)
    return X , ( Vs / np.pi ) * np.arctan2(X,d)




def R2_calc(y_obs,y_fit,with_r2_bis=False):
    #https://en.wikipedia.org/wiki/Coefficient_of_determination
    ybar = np.mean(y_obs)
    SStot = np.sum((y_obs - ybar)**2)
    SSreg = np.sum((y_fit - ybar)**2)
    SSres = np.sum((y_obs - y_fit)**2)

    r2 = 1. - ( SSres / SStot)
    r2bis = ( SSreg / SStot)

    if not with_r2_bis:
        return r2
    else:
        return r2 , r2bis

def R2_from_a_line_regress(Xobs,Yobs,a,b):
    #https://en.wikipedia.org/wiki/Coefficient_of_determination
    Xfit , Yfit = geok.linear_reg_getvalue(Xobs,a,b)
    r2 = R2_calc(Yobs,Yfit)
    return r2


def bins_middle(bin_edges):
    bin_edges2 = []
    for i in range(len(bin_edges)-1):
        bin_edges2.append((bin_edges[i]+bin_edges[i+1])/2.)
    return bin_edges2

def chi2_test_frontend(dist_inp,nbins=10,ddof=2,debug=0,mode2=False,aaa=1):
    """
    mode1 (par def) : on normalise la theorique et pas la observée
    mode2 : on normalise la distribution observée est pas la théorique
    INCOHERENT AVEC MATLAB => A EVITER

    Enfin bon, la manière dont on fabrique les valeurs théoriques est quand même
    un peu vaseuse ... penser à porter le code matlab chi2gof.m l.185
    Et comprendre aussi pourquoi ils ont un ddof de 2 (quon ajoute ici aussi
    par bete copiage) ...

    en debug mode bin_edges,bin_edges2,hist,gauss,chi2

    """
    hist , bin_edges =  np.histogram(dist_inp,nbins,density=mode2)
    bin_edges2 = bins_middle(bin_edges)
    #    gauss = scipy.stats.norm.pdf(bin_edges2,
    #                         np.mean(dist_inp),np.std(dist_inp))
    #    if not mode2:
    #        koefnorm = scipy.integrate.trapz(hist,bin_edges2)
    #    else:
    #        koefnorm = 1.

    cdf = scipy.stats.norm.cdf(bin_edges[1:-1],np.mean(dist_inp),np.std(dist_inp)*aaa)
    cdf = np.concatenate(([0],cdf,[1]))
    dif = np.diff(cdf)

    gauss = len(dist_inp) * dif

    #    gauss1 = gauss
    #    gauss = gauss * koefnorm
    try:
        chi2 = scipy.stats.chisquare(hist , gauss , ddof)
    except:
        chi2 = (np.nan,np.nan)
    if not debug:
        return chi2
    else:
        plt.plot(bin_edges2,gauss,'+')
        plt.plot(bin_edges2,hist,'x')
        return bin_edges,bin_edges2,hist,gauss,chi2


def chi2_test_lsq(V , A ,  P = None , fuvin = None , risk = 0.05,
                  cleaning_std = False , cleaning_normalized = False,
                  koefP = 1):

    """
    P est uniquement la diagonale de la matrice des poids

    les cleaning sont des astuces pour se rapprocher de 1 (en nettoyant les plus mauvaises valeurs)
    cleaning_normalized est à privilégier (et override cleaning_std)

    koefP est coefficient qu'on donne a P pour trouver une solution viable

    si koefP != 1, le nouveau P est donné en avant dernier argument

    """

    if fuvin is None and P is None:
        print("ERR : chi2_test_lsq : fuvin == None and P == None")
        return None

    ddl = (np.max(A.shape) - np.min(A.shape))

    if cleaning_std:
        bb = Vwork < np.std(Vwork) * 3

    if cleaning_normalized:
        Vnorma = V / np.sqrt(1/P)
        bb = Vnorma < 3

    if cleaning_std or cleaning_normalized:
        V = V[bb]
        P = P[bb]

    if koefP != 1:
        P = P * koefP

    if fuvin is not None:
        esfuv = fuvin
    else:
        esfuv = np.sum( V.T * P * V ) / ddl

        pmin = scipy.stats.chi2.ppf(risk/2,ddl) / ddl
        pmax = scipy.stats.chi2.ppf(1 - risk/2,ddl) / ddl

    if pmin < esfuv and esfuv < pmax:
        boolchi2 = True
    else:
        boolchi2 = False

    if koefP == 1:
        return esfuv , pmin , pmax , boolchi2
    else:
        return esfuv , pmin , pmax , P , boolchi2

def greek_alphabet(num=None,maj=False):
    if not num:
        a = ['\u03B1',
        '\u03B2',
        '\u03B3',
        '\u03B4',
        '\u03B5',
        '\u03B6',
        '\u03B7',
        '\u03B8',
        '\u03B9',
        '\u03BA',
        '\u03BB',
        '\u03BC',
        '\u03BD',
        '\u03BE',
        '\u03BF',
        '\u03C0',
        '\u03C1',
        '\u03C3',
        '\u03C4',
        '\u03C5',
        '\u03C6',
        '\u03C7',
        '\u03C8',
        '\u03C9']

        A=['\u0391',
        '\u0392',
        '\u0393',
        '\u0394',
        '\u0395',
        '\u0396',
        '\u0397',
        '\u0398',
        '\u0399',
        '\u039A',
        '\u039B',
        '\u039C',
        '\u039D',
        '\u039E',
        '\u039F',
        '\u03A0',
        '\u03A1',
        '\u03A3',
        '\u03A4',
        '\u03A5',
        '\u03A6',
        '\u03A7',
        '\u03A8',
        '\u03A9']

        if maj:
            return A
        else:
            return a
    else:
        return greek_alphabet()[num]

def sinusoide(T,A,omega,phi=0):
    """
    amplitude de la grandeur, appelée aussi valeur de crête, dans l'unité de la grandeur mesurée
    omega : pulsation de la grandeur en rad⋅s-1
    phi : phase instantanée en rad
    phi : phase à l'origine en rad (souvent fixée par l'expérimentateur)
    """
    return A * np.sin(omega * T + phi)

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y


def gaussian_filter_GFZ_style_smoother(tim_ref, dat_ref, width=7):
    """
    Gaussian filter to smooth data, based on
    GFZ's GMT_plus.pm/gaussian_kernel
      
    Args :
        tim_ref : the X/T component of the time serie (in decimal days !)
        
        dat_ref : the Y component (the data)
        
        width   : size of the window (odd number is best ?)
        
    Returns :
        dat_smt : smoothed Y
    
    NB :
        Some other nice ideas here
        http://scipy-cookbook.readthedocs.io/items/SignalSmooth.html
        https://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way
        https://stackoverflow.com/questions/32900854/how-to-smooth-a-line-using-gaussian-kde-kernel-in-python-setting-a-bandwidth
        
    NB2 : THIS VERSION IS VERY SLOW (DIRTY CONVERSION OF A PERL FCT)
          THE PYTHONIC VERSION gaussian_filter_GFZ_style_smoother_improved
          BELOW SHOULD BE USED  !!!
    """
        
    tim_raw = tim_ref
    dat_raw = dat_ref
      
    num_raw = len(tim_raw)
    icomp=0
    
    dat_smt = [np.nan] * len(dat_ref)
      
    for ismt in range(num_raw ): #+1
        
        x_val  = 0.0;
        x_wht  = 0.0;
          
        for iraw in reversed(range( 0 , ismt )): 
            x_lag = tim_raw[ismt] - tim_raw[iraw]
            x_fac = np.exp( - (x_lag / width ) ** 2 / 2 )
            x_val += dat_raw[iraw] * x_fac
            x_wht += x_fac
            icomp+=1
            if x_fac < 0.01:
                break
              
        for iraw in range(ismt+1 , num_raw):
            x_lag = tim_raw[ismt] - tim_raw[iraw]
            x_fac = np.exp( - ( x_lag / width ) ** 2 / 2 )
            x_val += dat_raw[iraw] * x_fac
            x_wht += x_fac
            icomp+=1;
            if x_fac < 0.01:
                break
          
        dat_smt[ismt] = x_val/x_wht
    
    return dat_smt



def gaussian_filter_GFZ_style_smoother_improved(tim_ref, dat_ref, width=7):
    """
    Gaussian filter to smooth data, based on
    GFZ's GMT_plus.pm/gaussian_kernel
      
    Args :
        tim_ref : the X/T component of the time serie (in decimal days !)
        
        dat_ref : the Y component (the data)
        
        width   : size of the window (odd number is best ?)
        
    Returns :
        dat_smt : smoothed Y
    
    NB :
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
        tim_raw_work = np.delete(tim_raw , ismt)
        dat_raw_work = np.delete(dat_raw , ismt)
        
        X_lag = tim_raw[ismt] - tim_raw_work
        X_fac = np.exp( - (X_lag / width ) ** 2 / 2 )
        
        #X_fac[X_fac < 0.01] = 0.

        clean_bool = X_fac > 0.01
        ## It differs a bit of the official fct, because the next element
        ## following this criteria is included anyway
        
        dat_raw_clean = dat_raw_work[clean_bool]
        X_fac_clean   = X_fac[clean_bool]
        
        dat_raw_clean = dat_raw_clean[~np.isnan(dat_raw_clean)]
        X_fac_clean   = X_fac_clean[~np.isnan(X_fac_clean)]

        X_val = np.sum(np.multiply(dat_raw_clean,X_fac_clean))
        X_wht = np.sum(X_fac_clean)
    
        dat_smt2.append(X_val/X_wht)
    
    dat_smt2 = np.array(dat_smt2)
    
    return dat_smt2
        




def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    NOTA PERSO : works only for equaly spaced data ....
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    SOURCE : http://scipy-cookbook.readthedocs.io/items/SignalSmooth.html
    """

    if x.ndim != 1:
        raise(ValueError, "smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise(ValueError, "Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise(ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=numpy.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='valid')
    return y

def gebco_bathy_grid_extractor(dataset,latmin,latmax,lonmin,lonmax):
    """
    for safety reasons, lat and lon input MUST BE in the dataset,
    replaced by the closest elsewhere

    return latnew , lonnew , Znew
    """

    import genefun as gf

    lon = dataset['lon'][:]
    lat = dataset['lat'][:]
    Z   = dataset['elevation'][:]

    latmin_dec = latmin
    latmax_dec = latmax
    lonmin_dec = lonmin
    lonmax_dec = lonmax

    #latmin_dec = geok.dms2dec_num(*latmin)
    #latmax_dec = geok.dms2dec_num(*latmax)
    #lonmin_dec = geok.dms2dec_num(*lonmin)
    #lonmax_dec = geok.dms2dec_num(*lonmax)


    if not np.any(latmin_dec == lat):
        print("ERR : ... replacing to the nearest")
        latmin_dec = gf.find_nearest(lat,latmin_dec)[0]

    if not np.any(latmax_dec == lat):
        print("ERR : ... replacing to the nearest")
        latmax_dec = gf.find_nearest(lat,latmax_dec)[0]

    if not np.any(lonmin_dec == lon):
        print("ERR : ... replacing to the nearest")
        lonmin_dec = gf.find_nearest(lon,lonmin_dec)[0]

    if not np.any(lonmax_dec == lon):
        print("ERR : ... replacing to the nearest")
        lonmax_dec = gf.find_nearest(lon,lonmax_dec)[0]

    boollat = np.logical_and(latmin_dec <= lat, lat <= latmax_dec)
    boollon = np.logical_and(lonmin_dec <= lon, lon <= lonmax_dec)


    gridlat = np.tile( boollat[:,np.newaxis]  , (1,len(lon)))
    gridlon = np.tile( boollon , (len(lat),1))

    Znew   = Z[ gridlon * gridlat ]
    latnew = lat[boollat]
    lonnew = lon[boollon]
    Znew   = Znew.reshape(len(latnew),len(lonnew))

    return latnew , lonnew , Znew


### High level geodetic function

def itrf_speed_calc(x0,y0,z0,t0,vx,vy,vz,t):
    """
    Args :
        x0,y0,z0 (floats) : coordinates at the reference epoch (m)

        t0 (float) : reference epoch (decimal year)

        vx,vy,vz (floats) : speed of the point (m/yr)

        t (float) : output epoch
    Returns :
        xout,yout,zout : coordinates of the point @ the ref. epoch (m)
    """

    xout = x0 + vx * ( t - t0 )
    yout = y0 + vy * ( t - t0 )
    zout = z0 + vz * ( t - t0 )

    return xout,yout,zout


def itrf_psd_fundamuntal_formula(t,A_l,t_l,tau_l,A_e,t_e,tau_e):
    """
    http://itrf.ensg.ign.fr/ITRF_solutions/2014/doc/ITRF2014-PSD-model-eqs-IGN.pdf
    """
    
    dL = A_l * np.log(1 + (t - t_l)/ tau_l) + A_e * (1 + (t - t_e)/ tau_e)







def calc_pos_speed_itrf(x0,y0,z0,t0,vx,vy,vz,t):
    """
    just a wrapper of itrf_speed_calc
    for legacy reasons
    """
    return itrf_speed_calc(x0,y0,z0,t0,vx,vy,vz,t)


def helmert_trans(Xa,params='itrf2008_2_etrf2000',invert=True,workepoc=2009.):
    """
    NB 1 : http://etrs89.ensg.ign.fr/memo-V8.pdf
    NB 2 :
    Transformation inverse : * -1 pour les paramètres
    cf https://en.wikipedia.org/wiki/Helmert_transformation

    optimisé pour RGF93 => ITRF2008 (d'ou le invert = True & workepoc = 2009.)

    NB3 : Attention losque l'on compare avec
          le convertisseur EUREF avec une vitesse
          parce que elle aussi est modifiée dans la conversion ETRS => ITRS.
          Conclusion : le faire en 2 étapes ETRS = > ITRS dans la même epoc
                                            ITRS epoc 1 => ITRS epoc 2

    if manual
    then params is a tuple
    params = (t1,t2,t3,dab,r1,r2,r3)

    NE MARCHE PAS PARCE BESOIN DU RATE(TAUX) DES PARAMS D'HELMERT !!!!!!
    (160923)
    """
    if invert:
        inver = -1.
    else:
        inver = 1.

    mas2rad = 0.0000000048481368
    mas2rad = 4.8481368111e-6 * 1e-3

    if params == 'itrf2008_2_etrf2000':

        t1rate = .1 * 10**-3
        t2rate = .1 * 10**-3
        t3rate = -1.8 * 10**-3
        dabrate = .08 * 10**-9
        r1rate =  .081 * mas2rad
        r2rate =  .490 * mas2rad
        r3rate = -.792 * mas2rad

        t1rate = .1 * 10**-3
        t2rate = .1 * 10**-3
        t3rate = -1.8 * 10**-3
        dabrate = .08 * 10**-9
        r1rate =  .081 * mas2rad
        r2rate =  .490 * mas2rad
        r3rate = -.792 * mas2rad


        t1  =(52.1  * 10**-3   + t1rate  * ( workepoc - 2000.)) *inver
        t2  =(49.3  * 10**-3   + t2rate  * ( workepoc - 2000.)) *inver
        t3  =(-58.5 * 10**-3   + t3rate  * ( workepoc - 2000.)) *inver
        dab =( 1.34 * 10**-9   + dabrate * ( workepoc - 2000.)) *inver
        r1  =( 0.891 * mas2rad + r1rate  * ( workepoc - 2000.)) *inver
        r2  =( 5.39  * mas2rad + r2rate  * ( workepoc - 2000.)) *inver
        r3  =( -8.712* mas2rad + r3rate  * ( workepoc - 2000.)) *inver


    elif params == 'itrf2000_2_etrf2000':
        t1  =54.0  * 10**-3   *inver
        t2  =51.0  * 10**-3   *inver
        t3  =-48.0 * 10**-3   *inver
        dab = 0.0  * 10**-9   *inver
        r1  = 0.891 * mas2rad *inver
        r2  = 5.390 * mas2rad *inver
        r3  = -8.712* mas2rad *inver




    R = np.matrix([[dab,-r3,r2],
                   [r3,dab,-r1],
                   [-r2,r1,dab]])

    Xb = Xa + np.matrix([t1,t2,t3]) + np.dot(R,Xa)
    Xb = np.squeeze(np.array(Xb))

    return Xb


def helmert_trans_estim_matrixs_maker(X1 , X2):
    x1 , y1 , z1 = X1
    x2 , y2 , z2 = X2
    
    
    block_1 = np.eye(3)
    
    block_2 = np.array([[ 0. , -z1,  y1, x1],
                       [ z1,   0., -x1, y1],
                       [-y1,  x1,   0., z1]])
    
    l = X2 - X1
    A = np.hstack((block_1 , block_2))

    return l,A


def helmert_trans_estim(X1list , X2list):
    """
    estimates 7 parameters of a 3D Helmert transformation between a set of points
    X1 and a set of points X2 (compute transformation X1 => X2)
    
    Inputs :
        X1list & X2list : list of N (x,y,z) points ,
        or an numpy array of shape (3, N)
    
    Outputs :
        Res : 7 Helmert params. : x,y,z translations, x,y,z rotations, scale
        A : Design matrix
        l : X2 - X1
        
    Source :
        https://elib.uni-stuttgart.de/bitstream/11682/9661/1/BscThesis_GaoYueqing.pdf
    """

    l_stk = []
    A_stk = []
    
    
    for X1 , X2 in zip(X1list , X2list):
        lmono , Amono = helmert_trans_estim_matrixs_maker(X1,X2)
            
        l_stk.append(lmono)
        A_stk.append(Amono)
        
    A = np.vstack(A_stk)
    l = np.hstack(l_stk)
    
    Res = scipy.linalg.inv((A.T).dot(A)).dot(A.T).dot(l)
    
    return Res , A , l


#  ______      _             _____      _      
# |  ____|    | |           |  __ \    | |     
# | |__  _   _| | ___ _ __  | |__) |__ | | ___ 
# |  __|| | | | |/ _ \ '__| |  ___/ _ \| |/ _ \
# | |___| |_| | |  __/ |    | |  | (_) | |  __/
# |______\__,_|_|\___|_|    |_|   \___/|_|\___|
#                                              
    
### Euler Pole determination


def euler_pole_calc(lat_ref,long_ref,vn_ref,ve_ref,
                    incvn_ref=None,incve_ref=None,Rt=6.3710088e6):
    """
    Compute the Euler pole of a set of reference points
    
    Written by C. Geisert (ENSTA/LIENSs) - 2017

    Parameters
    ----------
    lat_ref,long_ref : list or numpy.array
        latitude and longitude of the reference points (deg)

    vn_ref,ve_ref : list or numpy.array
        north and east velocities of the reference points (m/yr)
        
    incvn_ref,incve_ref : list or numpy.array
        uncertainties on north and east velocities of the reference points (m/yr)

    Rt : float
        Earth Radius
        IUGG value = 6.3710088e6
        IAU value  = 6.378e6
        
                
    Returns
    -------
    w : numpy.array
        Euler vector (rad/yr)
    
    wratedeg : float
        Rate of rotation (deg/Myr)
    
    wlat,wlong : float
        latitude and longitude of the Euler pole (deg)
        
    wwmat : numpy.array
        weight matrix (for debug)

    desmat : numpy.array
        design matrix (for debug)
        
    nrmatinv : 2-tuple
        output of scipy's lstsq fct (for debug)

    Source
    ------
    based on
    Goudarzi, M. A., Cocard, M., & Santerre, R. (2014). EPC: Matlab software to estimate Euler pole parameters. GPS Solutions, 18(1), 153–162. https://doi.org/10.1007/s10291-013-0354-4
        
    
    Notes
    -----
    w is the common element for all Euler pole functions
    
    Should remains in rad/yr
    """
    if (incve_ref is None) or (incvn_ref is None):
        incve_ref = np.ones(len(ve_ref))
        incvn_ref = np.ones(len(vn_ref))
            
    all_pos_ref =np.column_stack((lat_ref,long_ref))
    coords=(all_pos_ref*np.pi/180)  # [rad] position des stations de la plaque de ref lat,long [°] > lat,long [rad]
    dm=topo2dm(coords)  # [rad]
    desmat=Rt*(dm) #design matrix [m] [rad]
    Tdesmat=np.transpose(desmat) # [m] [rad]
    n  = len(all_pos_ref)
    vel_vector = list(np.zeros([2 * n,1]))     # vn,ve ; vn, ve alternes  vecteur (vn,ve)0 puis (vn,ve)1, puis stations suivantes..
    for  index in range (0,2*n-1,2):    
        vel_vector[index]   =  vn_ref[int(index/2)]
        vel_vector[index+1]   =   ve_ref[int(index/2)]#vn_ref, ve_ref > vel_vector    
    s=[0]*len(incvn_ref)                              
    covars=np.transpose(np.array( [incvn_ref,incve_ref,s])) #[m/yr]     
    wwmat=covarvec2wtmat(covars) #weight matrix of observations  #[m/yr]
    wwmat=np.array(wwmat[0])
    nrmat=Tdesmat.dot(wwmat).dot(desmat) # [m] [rad] 
    lennrmat=np.shape(nrmat) # []
    nrmatinv=np.linalg.lstsq(nrmat,np.eye(lennrmat[0]),rcond=None)  # [m] [│rad]         
    A=(Tdesmat.dot(wwmat)).dot((vel_vector))  # [m, rad] * [m/yr]   
    w=nrmatinv[0].dot(A) # taille 3x3  # [m, rad] * [m/yr]  #   [rad/yr]
#    wx=w[0]  #[rad/ yr]  
#    wy=w[1]  #[rad/ yr] 
#    wz=w[2]  #[rad/ yr] 
#    wrate=np.sqrt(wx**2+wy**2+wz**2)*1e6    # [rad/yr]> [rad/Myr]
#    wratedeg=wrate*180/np.pi  # [deg/Myr]
#    wlat=(np.arctan(wz/(np.sqrt(wx**2+wy**2))))*180/np.pi # [°] 
#    wlong=(np.arctan(wy/wx))*180/np.pi # [°]     

    wlat,wlong,wrate,wratedeg = euler_pole_vector_to_latlongrate(w)  
    
    return w,wratedeg,wlat,wlong,wwmat,desmat,nrmatinv


def euler_pole_vector_to_latlongrate(w):   
    """
    Convert Euler pole vector to latitude longitude and rate  
    
    Written by P. Sakic based on C. Geisert work
    
    Parameters
    ----------
    w : numpy.array
        Pole of rotation vector computed by euler_pole_calc (rad/yr)

    Returns
    -------
    wlat,wlong : float
        latitude and longitude of the Euler pole (deg)
        
    wrate : float
        Rate of rotation (rad/Myr)
        
    wratedeg : float
        Rate of rotation (deg/Myr)
        
    Notes
    -----
    w is the common element for all Euler pole functions
    
    Should remains in rad/yr
    """
    wx=w[0]  #[rad/ yr]  
    wy=w[1]  #[rad/ yr] 
    wz=w[2]  #[rad/ yr] 
    wrate=np.sqrt(wx**2+wy**2+wz**2)*1e6    # [rad/yr]> [rad/Myr]
    wratedeg=wrate*180/np.pi  # [deg/Myr]
    wlat=(np.arctan(wz/(np.sqrt(wx**2+wy**2))))*180/np.pi # [°] 
    wlong=(np.arctan(wy/wx))*180/np.pi # [°]   
    return wlat,wlong,wrate,wratedeg
    

def euler_vels_relative_to_ref(w,lat_ITRF,long_ITRF,vn_ITRF,ve_ITRF,
                               incvn_ITRF=None,incve_ITRF=None,Rt=6.378e6):
    """
    Compute relative velocities of points with respect to an reference plate/Euler pole
    
    Written by C. Geisert (ENSTA/LIENSs) - 2017
    
    Parameters
    ----------
    w : numpy.array
        Pole of rotation vector computed by euler_pole_calc (rad/yr)

    lat_ITRF,long_ITRF : list or numpy.array
        latitude and longitude of the points (deg)
        
    vn_ITRF,ve_ITRF : float or int or str or dict or n-tuple or bool or list or numpy.array
        north and east velocities of the points (m/yr)

    incvn_ITRF,incve_ITRF : float or int or str or dict or n-tuple or bool or list or numpy.array
        uncertainties on north and east velocities of the points (m/yr)
                              
    Returns
    -------
    vel_reltoref : numpy.array
        relative velocities of points with respect to the Euler pole w
    
    Notes
    -----
    w is the common element for all Euler pole functions
    
    Should remains in rad/yr
    """
    if (incvn_ITRF is None) or (incve_ITRF is None):
        incvn_ITRF = np.ones(len(vn_ITRF))
        incve_ITRF = np.ones(len(ve_ITRF))
        
    all_pos = np.column_stack((lat_ITRF,long_ITRF))
    nstation=len(all_pos)
    vel_rot=[0]*nstation
    for i in range(nstation): #
        mat=Rt*np.array([np.array([np.sin(all_pos[i][1]*np.pi/180), -np.cos(all_pos[i][1]*np.pi/180),0]),
                         np.array([-np.sin(all_pos[i][0]*np.pi/180)*np.cos(all_pos[i][1]*np.pi/180), -np.sin(all_pos[i][0]*np.pi/180)*np.sin(all_pos[i][1]*np.pi/180),np.cos(all_pos[i][0]*np.pi/180)])]) # [°] > [rad] pour le calcul
        vel_rot[i]=mat.dot(w)  #   [m] cross [rad/yr]
    #vel_rot=np.array(vel_rot)*1000  # [m/an]  > [mm/an]
    v_ITRF=np.transpose(np.array([vn_ITRF,ve_ITRF])) # les vitesses de toutes les stations [m/yr]
    vel_reltoref= v_ITRF-vel_rot  # [m/an]
    
    return vel_reltoref

def euler_pole_vector_from_latlongrate(wlat,wlong,wrate,
                                return_w_in_deg_per_Myr=False):
    """
    Compute the Euler vector from the pole latitude, longitude and rate
    
    Written by C. Geisert (ENSTA/LIENSs) - 2017

    Parameters
    ----------
    wlat,wlong : float
        latitude and longitude of the Euler pole (deg)

    wrate : float
        rate of the Euler pole (rad/Myr)
        
    Returns
    -------
    w : numpy.array
        Pole of rotation vector (rad/yr)

    Notes
    -----
    w is the common element for all Euler pole functions
    
    Should remains in rad/yr per default
    """
    # OUTPUT : w [rad/year]
     
    ### Conversion deg => rad coordinates
    wlong = np.deg2rad(wlong)    
    wlat  = np.deg2rad(wlat)
    
    wx=(wrate/1e6)*np.cos(wlat)*np.cos(wlong)#*np.pi/180  # [rad/yr]     
    wy=(wrate/1e6)*np.cos(wlat)*np.sin(wlong)#*np.pi/180  # [rad/yr]              
    wz=(wrate/1e6)*np.sin(wlat)#*np.pi/180                # [rad/yr]      
    w=[0]*3
    
    if return_w_in_deg_per_Myr:
        k = (180/np.pi) * 10**6
    else:
        k = 1.
    
    w[0]=wx*k
    w[1]=wy*k 
    w[2]=wz*k 
    
    return w



#%% ____Estimation de l'incertitude lors du calcul des EPP
def euler_pole_quality(w,vn_ref,ve_ref,nrmatinv,desmat,wwmat,
                       pretty_output=True):
    """ 
    Compute the uncertainties of the Euler pole determination

    Parameters
    ----------
    w : numpy.array
        Pole of rotation vector computed by euler_pole_calc 

    vn_ref,ve_ref : list or numpy.array
        north and east velocities of the reference points (m/yr)
        
    nrmatinv : 2-tuple
        output of scipy's lstsq fct from euler_pole_calc 
        
    wwmat : numpy.array
        weight matrix from euler_pole_calc 

    desmat : numpy.array
        design matrix from euler_pole_calc 
        
    pretty_output : bool
        if True, convert sigma_ww_latlon to pertinent units directly,
        returns raw units instead
                
    Returns
    -------
    sigma_ww : numpy.array
        Uncertainty on the Euler vector

    sigma_ww_latlon : numpy.array
        Uncertainty on the Euler pole : [rateSigma, latSigma, longSigma] 
        if pretty_output == True :  [deg/Myr,deg,deg] 
        if pretty_output == False : [rad/yr,rad,rad]
                
    dV_topo3 : numpy.array
        Residual velocities for the references points
        
    wrmse : float
        weigthed RMS on Residual velocities (m)
        
    wrmse_norm : float
        nomalized weigthed RMS on Residual velocities (m)
        
    rmse : float
        unweigthed RMS on Residual velocities (m)
        
    apost_sigma : float
        a-posteriori sigma (m)
        
    Notes
    -----
    w is the common element for all Euler pole functions
    
    Should remains in rad/yr
    """
    v_ref=np.transpose(np.array([vn_ref,ve_ref])) # les vitesses de toutes les stations [mm/yr]
    estm_vel=desmat.dot(w)
    estm_vel=np.transpose(np.array( [estm_vel[0::2],estm_vel[1::2]]))
    estm_vel_diff = v_ref- (estm_vel)  
    dV_topo= np.transpose(estm_vel_diff) 
    dV_topo=np.transpose(dV_topo)
    nd=len(dV_topo) 
    dV_topo=np.append(dV_topo[:,0],dV_topo[:,1])
    dV_topo2=(np.zeros([2*nd])) 
    for  index in range (0,2*nd-1,2):   
            dV_topo2[index]   = dV_topo[int(index/2)]
            dV_topo2[index+1] = dV_topo[int(index/2+nd)]                  
    df_value = len(dV_topo2) - 3;
    s0_2     = (np.transpose(dV_topo2).dot(wwmat).dot(dV_topo2)) / df_value;
    c_ww     = s0_2*(nrmatinv[0])  
    omega_length = np.sqrt(np.transpose(w) .dot( w))  
    #H  
    H00=  w[0]/omega_length
    H01=  w[1]/omega_length
    H02=  w[2]/omega_length
    H10= (-w[0] * w[2]) / (omega_length**2 * np.sqrt(w[0]**2 + w[1]**2))  
    H11= (-w[1] * w[2]) / (omega_length**2 * np.sqrt(w[0]**2 + w[1]**2)) 
    H12= -np.sqrt(w[0]**2 + w[1]**2) / omega_length**2  
    H20= -w[1]/(w[0]**2 + w[1]**2)
    H21=  w[0]/(w[0]**2 + w[1]**2)
    H22= 0
    H=np.array([[H00,H01,H02],
                [H10,H11,H12],
                [H20,H21,H22]])     
    c_ww_latlong = H.dot( c_ww ).dot( np.transpose(H) )  
    
    ### OUTPUT GENERATION
    sigma_ww         = np.sqrt(np.diag(c_ww)) 
    sigma_ww_latlong = np.sqrt(np.diag(c_ww_latlong)) 
    if pretty_output:
        rateSigma=((sigma_ww_latlong[0])*1e6*180/np.pi) # [rad/yr]> [°/Myr]
        latSigma=((sigma_ww_latlong[1])*180/np.pi) # [rad] >[°]
        longSigma=((sigma_ww_latlong[2])*180/np.pi) # [rad]>[°]
        sigma_ww_latlong = np.array([rateSigma,latSigma,longSigma])
    dV_topo3 = np.reshape(dV_topo2,(int(len(dV_topo2)/2),2))
    
    wrmse      = np.sqrt((np.transpose(dV_topo2).dot(wwmat)).dot(dV_topo2) / len(dV_topo2))
    wrmse_norm = np.sqrt((np.transpose(dV_topo2).dot(wwmat/(np.sum(np.diag(wwmat))))) .dot(dV_topo2) / len(dV_topo2))
    rmse       = np.sqrt(np.transpose(dV_topo2).dot(dV_topo2) / len(dV_topo2))
    apost_sigma = np.sqrt(s0_2)

    return sigma_ww,sigma_ww_latlong,dV_topo3,wrmse,wrmse_norm,rmse,apost_sigma

def topo2dm(coords):
    """
    Internal fuction for euler_pole_calc
    """
    n  = len(coords);
    dm = np.zeros([2 * n,3]) 
    for  index in range (0,2*n-1,2):       
        dm[index,0]   =  np.sin(coords[int((index+1)/2)][1]);
        dm[index,1]   = -np.cos(coords[int((index+1)/2)][1]);  
        dm[index+1,0] = -np.sin(coords[int((index+1)/2)][0]) * np.cos(coords[int((index+1)/2)][1]);
        dm[index+1,1] = -np.sin(coords[int((index+1)/2)][0]) * np.sin(coords[int((index+1)/2)][1]);
        dm[index+1,2] =  np.cos(coords[int((index+1)/2)][0]);   
    return dm


def covarvec2wtmat(covars):
    """
    Internal fuction for euler_pole_calc
    """
    wtmat = []
    aa=np.shape(covars)
    cc = np.zeros([aa[0] * 2,aa[0] * 2]);
    a=np.shape(cc)
    for idx in range( 0 ,a[0],2):    
        cc[idx,idx]         = covars[int(idx/2)][0]**2;
        cc[idx + 1,idx + 1] = covars[int(idx/2)][1]**2;

        cc[idx,idx + 1]     = covars[int(idx/2)][2];
        cc[idx + 1,idx]     = covars[int(idx/2)][2];
    
    cc[np.isnan(cc)] = 0;  
    wtmat = np.linalg.lstsq(cc,np.eye(a[0]),rcond=None) 
    return wtmat  








#### ASTRONOMY FUNCTION
def semi_major_axis_from_mean_motion(n):
    """
    source : https://space.stackexchange.com/questions/18289/how-to-get-semi-major-axis-from-tle
    """
    mu = 3.9860044189 * 10**14
    a  = (mu**(1./3.)) / ((2*n*np.pi/86400)**(2./3.))
    return a




# ============================= CIMETIERE DE FONCTIONS =======================
# ============================= CIMETIERE DE FONCTIONS =======================
# ============================= CIMETIERE DE FONCTIONS =======================


    #def exclude_data(Tin,Yin,biggap):
#
#    Tout = []
#    Yout = []
#
#    for i in range(len(Tin)):
#
#        if Tin[i] < biggap[0] or Tin[i] > biggap[1]:
#
#            Tout.append(Tin[i])
#            Yout.append(Yin[i])
#
#        else :
#            print str(Tin[i]) + " in gap"
#
#    return np.asarray(Tout) , np.asarray(Yout)

#def exclude_data(Tin,Yin,gaplist):
#
#    Tout = []
#    Yout = []
#
#    for i in range(len(Tin)):
#        gapbool = False
#
#        print 'aaaa'
#
#        for gap in gaplist:
#            if Tin[i] > gap[0] and Tin[i] < gap[1]:
#                print str(Tin[i]) + " in gap " + str(gap)
#                gapbool = True
#                break
#
#        if gapbool == False:
#            Tout.append(Tin[i])
#            Yout.append(Yin[i])
#        else:
#            print 'bbb'
#
#    return np.asarray(Tout) , np.asarray(Yout)
#
#def find_pattern(listin):
#
#    patlen = 1
#    lenmax = len(listin)
#
#    boolgood = False
#
#    patlenlist = []
#
#    while (( patlen < lenmax )):
#
#        motif1 = listin[:patlen]
#        motif2 = listin[patlen:patlen*2]
#
#        if motif1 == motif2:
#            boolgood = True
#            patlenlist.append(patlen)
#            patlen = patlen + 1
#
#        else:
#            patlen = patlen + 1
#
#    if patlenlist[0] == 1:
#        return patlenlist[1]
#    else:
#        return patlenlist[0]


#def find_timegap0(Tin,pas=0,marge=1.1,advenced=False,margebiggap=100):
#    ### ==== PAS FINI ====
#
#    if isinstance(Tin[0],dt.datetime):
#        T = dt2posix_list(Tin)
#    else:
#        T = Tin
#
#    dT = np.diff(T)
#
#    if pas == 0:
#        pas=np.mean(dT)
#
#    # dtbad_bool = les delta T qui sont trop grands = True p/r au pas nominal
#    dtbad_bool = dT > pas * marge
#
#    dtbad_ind = np.where(dtbad_bool == True)[0]
#    dtbad = dT[np.where(dtbad_bool == True)[0]]
#
#
#    allperiod = T[-1] - T[0]
#
#    gaplist = []
#
#    for ind in dtbad_ind:
#        gaplist.append([T[ind],T[ind+1]])
#
#    biggap = [gaplist[0][0],gaplist[-1][1]]
#
#
#
#    if advenced == True:
#
#        antigaplist = []
#        antidtlist = []
#
#
#        for i in range(len(gaplist)-1):
#            end1stgap = gaplist[i][1]
#            start2ndgap = gaplist[i+1][0]
#
#            antidtlist.append(start2ndgap - end1stgap)
#            antigaplist.append([end1stgap,start2ndgap])
#
#        #for i in range(len(antidtlist)):
#
#    return gaplist , dtbad, biggap

# CHANGEMENT DE REPERE OBSOLETE
#def C_ned2enu():
#    return np.array([[0,1,0],[1,0,0],[0,0,-1]])
#
#
#def C_ned2ecef(phi,llambda,angtype='deg'):
#
#    if angtype == 'deg':
#        phi = np.deg2rad(phi)
#        llambda = np.deg2rad(llambda)
#
#
#    In_ned = np.array([-np.cos(llambda) * np.sin(phi),-np.sin(llambda)*np.sin(phi),np.cos(phi)])
#    Ie_ned = np.array([-np.sin(llambda) , np.cos(llambda) , 0])
#    Id_ned = np.array([-np.cos(llambda) * np.cos(phi),-np.sin(llambda) * np.cos(phi),-np.sin(phi)])
#
#    Ix_ned = np.array([-np.cos(llambda) * np.sin(phi),-np.sin(llambda),-np.cos(llambda)*np.cos(phi)])
#    Iy_ned = np.array([-np.sin(llambda) * np.sin(phi),np.cos(llambda),-np.sin(llambda)*np.cos(phi)])
#    Iz_ned = np.array([np.cos(phi),0,-np.sin(phi)])
#
#
#    C_ned2ecef = np.array([[np.dot(Ix_ned.T,In_ned),np.dot(Ix_ned.T,Ie_ned),np.dot(Ix_ned.T,Id_ned)],
#                    [np.dot(Iy_ned.T,In_ned),np.dot(Iy_ned.T,Ie_ned),np.dot(Iy_ned.T,Id_ned)],
#                    [np.dot(Iz_ned.T,In_ned),np.dot(Iz_ned.T,Ie_ned),np.dot(Iz_ned.T,Id_ned)]])
#
#    return C_ned2ecef
#
#def C_ecef2ned(phi,llambda,angtype='deg'):
#    return np.linalg.inv(C_ned2ecef(phi,llambda,angtype))
#
#
#def C_enu2ecef(phi,llambda,angtype='deg'):
#    if angtype == 'deg':
#        phi = np.deg2rad(phi)
#        llambda = np.deg2rad(llambda)
#
#    Ie_enu = np.array([-np.sin(llambda),np.cos(llambda),0])
#    In_enu = np.array([-np.cos(llambda)*np.sin(phi),-np.sin(llambda)*np.sin(phi),np.cos(phi)])
#    Iu_enu = np.array([np.cos(llambda)*np.cos(phi),np.sin(llambda)*np.cos(phi),np.sin(phi)])
#
#    Ix_enu = np.array([-np.sin(llambda),-np.cos(llambda) * np.sin(phi),np.cos(llambda)*np.cos(phi)])
#    Iy_enu = np.array([np.cos(llambda),-np.sin(llambda)*np.sin(phi),np.sin(llambda)*np.cos(phi)])
#    Iz_enu = np.array([0,np.cos(phi),np.sin(phi)])
#
#    C_enu2ecef = np.array([[np.dot(Ix_enu.T,Ie_enu),np.dot(Ix_enu.T,In_enu),np.dot(Ix_enu.T,Iu_enu)],
#                [np.dot(Iy_enu.T,Ie_enu),np.dot(Iy_enu.T,In_enu),np.dot(Iy_enu.T,Iu_enu)],
#                [np.dot(Iz_enu.T,Ie_enu),np.dot(Iz_enu.T,In_enu),np.dot(Iz_enu.T,Iu_enu)]])
#
#    return C_enu2ecef
#
#def C_ecef2enu(phi,llambda,angtype='deg'):
#    return np.linalg.inv(C_enu2ecef(phi,llambda,angtype))
#
#
#def C_rpy2enu(roll,pitch,yaw,angtype='deg'):
#
#    if angtype == 'deg':
#        roll = np.deg2rad(roll)
#        pitch = np.deg2rad(pitch)
#        yaw = np.deg2rad(yaw)
#
#    R = roll
#    P = pitch
#    Y = yaw
#
#    Ir = np.array([np.sin(Y)*np.cos(P),np.cos(Y)*np.cos(P),np.sin(P)])
#    Ip = np.array([np.cos(R) * np.cos(Y) + np.sin(R)*np.sin(Y)*np.sin(P),
#                   -np.cos(R)*np.sin(Y)+np.sin(R)*np.cos(Y)*np.sin(P),
#                   -np.sin(R)*np.cos(P)])
#    Iy = np.array([-np.sin(R) *  np.cos(Y) + np.cos(R) * np.sin(Y) * np.sin(P) ,
#                   np.sin(R) * np.sin(Y) + np.cos(R) * np.cos(Y) * np.sin(P) ,
#                   - np.cos(R) * np.cos(P)    ])
#
#
#
#    Ie = np.array([ np.sin(Y) * np.cos(P) ,
#    np.cos(R) * np.cos(Y) + np.sin(R) * np.sin(Y) * np.sin(P) ,
#    -np.sin(R) * np.cos(Y) + np.cos(R) * np.sin(Y) * np.sin(P) ])
#
#    In = np.array([ np.cos(Y) * np.cos(P),
#     -np.cos(R) * np.sin(Y) + np.sin(R) * np.cos(Y) * np.sin(P),
#     np.sin(R) * np.sin(Y) + np.cos(R) * np.cos(Y) * np.sin(P) ])
#
#    Iu = np.array([ np.sin(P),
#     -np.sin(R) * np.cos(P),
#     -np.cos(R) * np.cos(P) ])
#
#    C_rpy2enu = np.hstack((Ir[np.newaxis].T,Ip[np.newaxis].T,Iy[np.newaxis].T))
#    #C_rpy2enu2 = np.vstack((Ie,In,Iu))
#
#    return C_rpy2enu
