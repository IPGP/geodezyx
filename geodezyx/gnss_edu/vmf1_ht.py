# -*- coding: utf-8 -*-
"""
"""

import numpy as np
    
def vmf1_ht(ah = None,aw = None,dmjd = None,dlat = None,ht = None,zd = None): 
# =============================================================================
#     !!! This is the version with height correction !!!
#     !!! It has to be used with the grid !!!
#
#     This subroutine determines the VMF1 (Vienna Mapping Functions 1)
#     Reference: Boehm, J., B. Werl, H. Schuh (2006),
#     Troposphere mapping functions for GPS and very long baseline interferometry
#     from European Centre for Medium-Range Weather Forecasts operational analysis data,
#     J. Geoph. Res., Vol. 111, B02406, doi:10.1029/2005JB003629.
#
#     Please mind that the coefficients in this paper are wrong. The corrected version of
#     the paper can be found at:
#     http://ggosatm.hg.tuwien.ac.at/DOCS/PAPERS/2006Boehm_etal_VMF1.pdf
#
#     input data
#     ----------
#     ah:   hydrostatic coefficient a (http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/)
#     aw:   wet coefficient a         (http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/)
#     dmjd: modified julian date
#     dlat: ellipsoidal latitude in radians
#     ht:   ellipsoidal height in meter
#     zd:   zenith distance in radians
#
#     output data
#     -----------
#     vmf1h: hydrostatic mapping function
#     vmf1w: wet mapping function
#
#     Johannes Boehm, 2005 October 2
#      Rev 2011 July 21: latitude -> ellipsoidal latitude

#     implicit double precision (a-h,o-z)
    
#     pi = 3.14159265359d0
    
#     reference day is 28 January
#     this is taken from Niell (1996) to be consistent
#     File Converted to python script by Zohreh Adavi zohreh.adavi@tuwien.ac.at, 2023-02-20
# =============================================================================  
    doy = dmjd - 44239.0 + 1 - 28
    bh = 0.0029
    c0h = 0.062
    if (dlat < 0):
        phh = np.pi
        c11h = 0.007
        c10h = 0.002
    else:
        phh = 0.0
        c11h = 0.005
        c10h = 0.001
    
    ch = c0h + ((np.cos(doy / 365.25 * 2.0 * np.pi + phh) + 1.0) * c11h / 2.0 + c10h) * (1.0 - np.cos(dlat))
    sine = np.sin(np.pi / 2.0 - zd)
    beta = bh / (sine + ch)
    gamma = ah / (sine + beta)
    topcon = (1.0 + ah / (1.0 + bh / (1.0 + ch)))
    vmf1h = topcon / (sine + gamma)
    # C  height correction for hydrotatic part [Niell, 1996]
    a_ht = 2.53e-05
    b_ht = 0.00549
    c_ht = 0.00114
    hs_km = ht / 1000.0
    beta = b_ht / (sine + c_ht)
    gamma = a_ht / (sine + beta)
    topcon = (1.0 + a_ht / (1.0 + b_ht / (1.0 + c_ht)))
    ht_corr_coef = 1.0 / sine - topcon / (sine + gamma)
    ht_corr = np.multiply(ht_corr_coef,hs_km)
    vmf1h = vmf1h + ht_corr
    bw = 0.00146
    cw = 0.04391
    beta = bw / (sine + cw)
    gamma = aw / (sine + beta)
    topcon = (1.0 + aw / (1.0 + bw / (1.0 + cw)))
    vmf1w = topcon / (sine + gamma)
    return vmf1h,vmf1w