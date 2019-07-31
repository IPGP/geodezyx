# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 15:40:53 2015

@author: psakicki
"""

import numpy as np
import scipy
from scipy import interpolate
import scipy.optimize as optimize
import scipy.stats
import scipy.sparse as sprs
import time
import datetime as dt
import dateutil
import os
import inspect
import itertools
import multiprocessing as mp
import glob
import platform
import shutil
import subprocess
import dateutil.parser
import collections
import genefun as gf
import re
import copy
import sympy
import pandas as pd
import tabulate

#if platform.node() in ('calipso','diamant'):
#    import cgkit.cgtypes as cgt


#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#import mayavi.mlab as mlab

#lib persos geodezyx
import genefun
import geodetik as geok
import geoclass as gcls
import geo_files_converter_lib as gfc

# libs persos acoustic
if platform.node() in ('calipso' , "psakicki-MS-16Y1"):
    import acouclass as acls
    import raytrace as rt
    import SSP as ssp
    import MOVEclass as mcls
if platform.node() in ('calipso'):    
    import acoustic_ranging_lib as arl


#if platform.node() in ('calipso' , "diamant"):
#    import acoustic_ranging_lib as arl

np.set_printoptions(suppress=True)



#import matplotlib
## matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
#import matplotlib.pyplot as plt

