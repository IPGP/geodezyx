#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 13:58:51 2021

@author: psakicki
"""

import datetime as dt
#### Import the logger
import logging

import numpy as np
import scipy

from geodezyx import utils, conv

log = logging.getLogger('geodezyx')


#  _______ _                   _____           _             _____       _                        _       _   _             
# |__   __(_)                 / ____|         (_)           |_   _|     | |                      | |     | | (_)            
#    | |   _ _ __ ___   ___  | (___   ___ _ __ _  ___  ___    | |  _ __ | |_ ___ _ __ _ __   ___ | | __ _| |_ _  ___  _ __  
#    | |  | | '_ ` _ \ / _ \  \___ \ / _ \ '__| |/ _ \/ __|   | | | '_ \| __/ _ \ '__| '_ \ / _ \| |/ _` | __| |/ _ \| '_ \ 
#    | |  | | | | | | |  __/  ____) |  __/ |  | |  __/\__ \  _| |_| | | | ||  __/ |  | |_) | (_) | | (_| | |_| | (_) | | | |
#    |_|  |_|_| |_| |_|\___| |_____/ \___|_|  |_|\___||___/ |_____|_| |_|\__\___|_|  | .__/ \___/|_|\__,_|\__|_|\___/|_| |_|
#                                                                                    | |                                    
#                                                                                    |_|                                    

class Interp1dTime(scipy.interpolate.interp1d):
   """
   Interpolation with datetime as inputs

   
   This class inherites from scipy.interpolate.interp1d
   and can take as input datetime as X
   
   Pierre Sakic 2020-01
   """
   def __init__(self, x, y, kind='linear', axis=-1,
                copy=True, bounds_error=None, fill_value=np.nan,
                assume_sorted=False):
       
       ### x (time) is converted as an array
       x = np.array(x)
       
       ### the datetime is converted to posix
       if isinstance(x[0],dt.datetime):
           xposix = conv.dt2posix(x)
       elif isinstance(x[0],np.datetime64):
           ### xposix = conv.dt2posix(conv.numpy_dt2dt(x))
           # dt2posix handles now (202108) the np datetimes 
           # bc of a pb of time resolution
           xposix = conv.dt2posix(x)
       else:
           xposix = x
       
       super().__init__(xposix, y, kind=kind, axis=axis,
                copy=copy, bounds_error=bounds_error, fill_value=fill_value,
                assume_sorted=assume_sorted)
       
   def __call__(self,x):
       ### x (time) is converted as an array if it is a mono-elt
       if not utils.is_iterable(x):
           singleton = True
           x = np.array([x])
       else:
           singleton = False
           x = np.array(x)

       ### the datetime is converted to posix
       if isinstance(x[0],dt.datetime):
           xposix = conv.dt2posix(x)
       elif isinstance(x[0],np.datetime64):
           ### xposix = conv.dt2posix(conv.numpy_dt2dt(x))
           # dt2posix handles now (202108) the np datetimes 
           # bc of a pb of time resolution
           xposix = conv.dt2posix(x)
       else:
           xposix = x
       
       #### Interpolation via the heritagle
       out = super().__call__(xposix)
       
       if singleton:
           return out[0]
       else:
           return out
   

  
class SlerpTime(scipy.spatial.transform.Slerp):
   """
   Slerp interpolation (for quaterinons) with datetime as inputs
   
   This class inherites from scipy.spatial.transform.Slerp
   and can take as input datetime as X
   
   p. Sakic 2020-01
   """
   
   def __init__(self, times, rotations,extrapolate=True):
       from scipy.spatial.transform import Rotation

       ### time is converted as an array
       times = np.array(times)
       
       
       ### the datetime is converted to posix
       if isinstance(times[0],dt.datetime):
           times_posix = conv.dt2posix(times)
       elif isinstance(times[0],np.datetime64):
           ### times_posix = conv.dt2posix(conv.numpy_dt2dt(times))
           # dt2posix handles now (202108) the np datetimes 
           # bc of a pb of time resolution
           times_posix = conv.dt2posix(times)
       else:
           times_posix = times
                      
       #### check if the time_posix are stricly assending
       bool_assending = np.all(np.diff(times_posix) > 0)
       if not bool_assending:
           log.debug(list(times_posix))
           log.debug(list(np.diff(times_posix)))
           IndexBad = np.where(np.diff(times_posix) <= 0)
           log.debug(list(IndexBad))
           log.debug(list(times_posix[IndexBad]))
           log.debug(list(np.diff(times_posix)[IndexBad]))
           
           log.warning("times_posix is not stricly assending, it will crash")
           log.warning("then we show below the times_posix for debug ")
       
       #### For the extrapolation 
       # first value => begining of posix era
       # last value => end of posix era 
       if extrapolate:
           times_posix = np.array(times_posix)
           times_posix = np.insert(times_posix,0,0.)
           times_posix = np.append(times_posix,2147483646.)
           
           rotations_list = list(rotations)
           rotations_list.insert(0,rotations_list[0])
           rotations_list.append(rotations_list[-1])
           rotations = Rotation.from_quat([r.as_quat() for r in rotations_list])
           
       super().__init__(times_posix, rotations)
       
   def __call__(self,times):
       ### time is converted as an array
       times = np.array(times)
       
       if isinstance(times[0],dt.datetime):
           times_posix = conv.dt2posix(times)
       elif isinstance(times[0],np.datetime64):
           ### times_posix = conv.dt2posix(conv.numpy_dt2dt(times))
           # dt2posix handles now (202108) the np datetimes 
           # bc of a pb of time resolution
           times_posix = conv.dt2posix(times)
       else:
           times_posix = times
           
       return super().__call__(times_posix)
   
