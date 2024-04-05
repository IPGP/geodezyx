#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#### IMPORT EXTERNAL MODULES
# import bisect
# from collections import defaultdict
# from collections import Counter
# import copy
# import datetime as dt
# import fnmatch
# from ftplib import FTP
# import glob
# import inspect
# import itertools
# import linecache
import logging
import logging.config
# import math
# import matplotlib
#matplotlib.use("agg") ### avoid tk import error
# import matplotlib.pyplot as plt
# import multiprocessing as mp
## from natsort import natsorted, ns
# import numpy as np
# import numpy as npaa
# import operator
import os 
from os import remove, close, path
# import pandas as df
# import pandas as pd
# import pickle
# import scipy
# from scipy.signal import butter, lfilter, freqz
# import shutil
# import string
# import struct
# import sys
# import subprocess
# import re
# import tempfile
# from tempfile import mkstemp
# import tabulate
# import time
# import urllib
# import uuid


#### IMPORT CONFIG FOR LOGGER
log_file_path = os.path.join(path.dirname(path.abspath(__file__)),'logconfig','loggzyx.py')
from . import logconfig
if os.path.isfile(log_file_path):
    ##print("INFO:",log_file_path,"found")
    ###### old config file format
    ###logging.config.fileConfig(fname=log_file_path, disable_existing_loggers=False)
    ###### new dict
    logging.config.dictConfig(logconfig.loggzyx.log_config_dict)
else:
    print("ERR:logger config file",log_file_path,"is missing")


#### IMPORT GEODEZYX INTERNAL SUBMODULES
from . import externlib
from . import megalib

from . import athmo
from . import conv
from . import files_rw
from . import geodyn
from . import marine
from . import operational
from . import reffram
from . import stats
from . import time_series
from . import utils


__all__ = ['athmo',
           'conv',
           'externlib',
           'files_rw',
           'geodyn',
           'marine',
           'operational',
           'reffram',
           'stats',
           'time_series',
           'utils']




#### Import extern libraires in in the geodezyx namespace
#from geodezyx.extern import *

# NB : if you want extern in the main namespace 
#      don't call with 
#      import geodezyx
#      but with
#      from geodezyx import *




