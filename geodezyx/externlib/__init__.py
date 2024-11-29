"""
Created on Thu Aug 13 15:11:12 2020

@author: psakic

This module geodezyx.externlib import easily all the useful extrern Python's packages.

it can be imported directly with:
from geodezyx.externlib import *  

!!!! THE EXTERN LIB IS JUST A LAZY TRICK FOR THE DEVELOPPER !!!!
!!!! TO IMPORT ALL THE USEFUL PACKAGES IN ONE LINE !!!!
!!!! IT SHOULD NOT BE USED OPERATIONNALY !!!!
!!!! IT IS NOT A GOOD PRACTICE TO IMPORT ALL THE PACKAGES IN ONE LINE !!!!

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""


from .externlib import *

# #### IMPORT EXTERNAL MODULES
import bisect
from collections import defaultdict
from collections import Counter
import copy
import datetime as dt
import fnmatch
from ftplib import FTP
import glob
import inspect
import itertools
import linecache
import math
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing as mp
from natsort import natsorted, ns
import numpy as np
import numpy as npaa
import operator
import os 
from os import remove, close
import pandas as df
import pandas as pd
import pickle
import scipy
from scipy.signal import butter, lfilter, freqz
import shutil
import string
import struct
import sys
import subprocess
import re
import tempfile
from tempfile import mkstemp
import tabulate
import time
import urllib
import uuid

