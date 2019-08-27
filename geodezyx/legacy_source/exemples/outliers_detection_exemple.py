#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 17:02:48 2019

@author: psakicki
"""

import geoclass as gcls
import geok_statistics as gstats


### Let's load a time serie file from M
p = "/home/psakicki/Downloads/RAIL.NA12.tenv3"

TS = gcls.read_nevada(p)

### We extract the interesting ENU component as array
E,N,U,T,sE,sN,sU = TS.to_list("ENU")

## plot the raw time serie
plt.plot(T,E,".",label="Raw")


#clean the outlier of Y usind MAD approach 
#and clean the corresponding values in X
#assuming that we have the function : X => Y(X)
#(be carefull, Y is the first argument)

Emad_clean , Tmad_clean = gstats.outlier_mad_binom(E,T,threshold=3.5,
                                                   verbose=False, 
                                                   detrend_first=True, 
                                                   return_booleans=False)

#Parameters
#----------
#Y : list or numpy.array
#    Values
#
#X : list or numpy.array
#    X Values so as X => Y(X)
#
#threshold : float
#    MAD threshold        
#    
#verbose : bool
#    
#detrend_first : bool
#    detrend linear behavior of Y(X) first
#    
#return_booleans : bool
#    return good and bad values of Y and X as booleans
#  
#Returns
#-------
#Yclean & Xclean : numpy.array
#
#bb : numpy.array (if return_booleans == True)
#    Y-sized booleans

## plot the MAD-cleaned time serie
plt.plot(Tmad_clean,Emad_clean,"+",label="MAD")

#Use the function outlier_above_below to remove residual outliers 
#Here we remove points which are 2mm above or below the median value
#(after a detrend)

Eab_clean , Tab_clean = gstats.outlier_above_below_binom(Emad_clean , Tmad_clean, 
                                                         threshold_values= 0.002 ,
                                                         reference = np.median  , 
                                                         theshold_absolute = True ,
                                                         theshold_relative_value = np.std,
                                                         return_booleans   = False,
                                                         detrend_first     = True,
                                                         verbose           = True)


#    Gives values of Y which are between threshold values, and correct an 
#    associated X so as X => Y(X)
#
#    Parameters
#    ----------
#    threshold_values : single value (float) or a 2-tuple 
#        (lower bound theshold , upper bound theshold)
#        
#        `WARN` : those value(s) have to be positives.
#        Minus sign for lower bound and plus sign for upper 
#        one will be applied internally
#        
#    reference : float or callable
#        the central reference value
#        can be a absolute fixed value (float) or 
#        a function (e.g. np.mean of np.median)
#
#    theshold_absolute : bool
#        if True threshold_values are absolutes values
#            >>> low = reference - threshold_values[0] 
#            >>> upp = reference + threshold_values[1] 
#        if False they are fractions of theshold_relative_value 
#            >>> low = reference - threshold_values[0] * theshold_relative_value 
#            >>> upp = reference + threshold_values[1] * theshold_relative_value
#        (see also below)
#    
#    theshold_relative_value : str or function
#        if the string "reference" or None is given, then it the reference 
#        value which is used
#        if it is a fuction (e.g. np.std()) then it is this value returned
#        by this function which is used
#        Only useful when theshold_absolute = False
#        
#    detrend_first : bool
#        detrend linear behavior of Y(X) first
#        Recommended
#        
#    return_booleans : bool
#        return booleans or not
#
#    verbose : bool
#                
#        
#    Returns
#    -------
#    Xout : numpy array
#        X between low_bound & upp_bound
#        
#    bbool : numpy array
#        X-sized array of booleans


#plot the above_below cleaned points
plt.plot(Tab_clean , Eab_clean , "x", label="above below")

plt.legend()


