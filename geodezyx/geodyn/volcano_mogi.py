#!/usr/bin/env python3

"""
Created on Mon Mar 14 12:58:32 2022
@author: psakic

This sub-module of geodezyx.geodyn contains functions to for forward 
volcano-geodesy analytic models.

It is direcly inspired by the work of Scott Henderson
https://github.com/scottyhq/cov9

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU GPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (GFZ, pierre.sakic@gfz-postdam.de)
GitHub repository :
https://github.com/GeodeZYX/GeodeZYX-Toolbox_v4
"""



"""

"""
import numpy as np
import matplotlib.pyplot as plt

# =====================
# Inverse Models
# =====================
def invert(xargs,xcen,ycen,depth,dV):
    """
    Wrapper of mogi.forward to project to LOS and adjust arguments to work
    with scipy.omptimize.curvefit. Assumes UTM input for X and Y
    """
    #NOTE: nu fixed to default 0.25 by leaving out
    X,Y,incidence,heading = xargs
    ux, uy, uz = forward(X,Y,xcen,ycen,depth,dV)
    dataVec = np.dstack([ux, uy, uz])
    cart2los = -get_cart2los(incidence,heading)
    los = np.sum(dataVec * cart2los, axis=2)

    return los.ravel()


# =====================
# Forward Models
# =====================
def forward(x,y,xcen=0,ycen=0,d=3e3,dV=1e6, nu=0.25):
    """
    Calculates surface deformation based on point source

    References: Mogi 1958, Segall 2010 p.203

    Args:
    ------------------
    x: x-coordinate grid (m)
    y: y-coordinate grid (m)

    Kwargs:
    -----------------
    xcen: y-offset of point source epicenter (m)
    ycen: y-offset of point source epicenter (m)
    d: depth to point (m)
    dV: change in volume (m^3)
    nu: poisson's ratio for medium

    Returns:
    -------
    (ux, uy, uz)
    """
    # Center coordinate grid on point source
    x = x - xcen
    y = y - ycen

    # Convert to surface cylindrical coordinates
    th, rho = cart2pol(x,y)
    R = np.hypot(d,rho)

    # Mogi displacement calculation
    C = ((1-nu) / np.pi) * dV
    ur = C * rho / R**3
    uz = C * d / R**3

    ux, uy = pol2cart(th, ur)

    return np.array([ux,uy,uz])


def forward_dp(x,y,xcen=0,ycen=0,d=3e3,a=500,dP=100e6,mu=4e9,nu=0.25):
    """
    dP instead of dV, NOTE: dV = pi * dP * a**3 / mu
    981747.7 ~ 1e6
    """
    dV = np.pi * dP * a**3 / mu
    return forward(x,y,xcen,ycen,d,dV,nu)


# =====================
# Utilities
# =====================
def dP2dV(dP,a,mu=30e9):
    dV = (np.pi * dP * a**3) / mu
    return dV

def dV2dP(dV,a,mu=30e9):
    dP = (dV * mu) / (np.pi * a**3)
    return dP

def cart2pol(x1,x2):
    theta = np.arctan2(x2,x1)
    r = np.hypot(x2,x1)
    return theta, r


def pol2cart(theta,r):
    x1 = r * np.cos(theta)
    x2 = r * np.sin(theta)
    return x1,x2

def get_cart2los(incidence,heading):
    '''
    coefficients for projecting cartesian displacements into LOS vector
    '''
    incidence = np.deg2rad(incidence)
    heading = np.deg2rad(heading)

    EW2los = np.sin(heading) * np.sin(incidence)
    NS2los = np.cos(heading) * np.sin(incidence)
    Z2los = -np.cos(incidence)

    cart2los = np.dstack([EW2los, NS2los, Z2los])

    return cart2los

# =====================
# Benchmark
# =====================
def benchmark(normalize=False,
              params = None,
              mesh_size=15000,
              color="blue",
              figure=None):
    """
    Mogi Source in an elastic halfspace
    (Segall Figure 7.5)
    """
    # Set parameters
    
    if not params:
        params = dict(xcen = 0,
                      ycen = 0,
                      d = 3e3, #m
                      dV = 1e6, #m^3
                      nu = 0.25)
    depth = params['d']

    npts = 100
    
    x = np.linspace(-mesh_size,mesh_size,npts)
    y = np.linspace(-mesh_size,mesh_size,npts)
    X,Y = np.meshgrid(x,y)

    # Run mogi model with delta volume input
    dx,dy,dz = forward(X,Y,**params)
    dr = np.hypot(dx,dy)

    # Normalize results
    
    nptshalf =  int(npts/2)
    
    if normalize:
        z = dz[nptshalf, nptshalf:] / dz.max()
        r = dr[nptshalf, nptshalf:] / dz.max()
        x = x[nptshalf:] / depth
    else:
        z = dz[nptshalf, nptshalf:] 
        r = dr[nptshalf, nptshalf:]
        x = x[nptshalf:] 
        
    # Reproduce the figure
    if not figure:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        fig = figure
        ax = plt.gca()
    
    label = "  " + str(d) 
    ax.plot(x, z,'-', lw=3, label='dz' + label,color=color)
    ax.plot(x, r,'--', lw=3, label='dr' + label,color=color)
    plt.legend()
    plt.grid(True)
    plt.title('Mogi Displacements')
    
    if normalize:
        plt.xlabel('normalized distance (r/d, m)')
        plt.ylabel('normalized displacement (dxi / dz.max, m)')
    else:
        plt.xlabel('distance (m)')
        plt.ylabel('displacement (dxi)')
        
    plt.show()


if __name__ == '__main__':
    
    figure = plt.figure()

    for idd,d in enumerate((500,600,700,800,900,1000)):
        params = dict(xcen = 0,
                      ycen = 0,
                      d = d, #m
                      dV = 1e6, #m^3
                      nu = 0.25)
        
        benchmark(params=params,
                  figure=figure,
                  mesh_size=200,
                  color="C" + str(idd))
