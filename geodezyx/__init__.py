#### IMPORT EXTERNAL MODULES
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


#### IMPORT GEODEZYX INTERNAL SUBMODULES
from . import externlib
from . import megalib

from . import athmo
from . import conv
from . import files_rw
from . import geodyn
#from . import kepler_core
from . import operational
#from . import legacy
from . import reffram
from . import stats
from . import time_series
from . import utils


__all__ = ['athmo',
           'conv',
           'externlib',
           'files_rw',
           'geodyn',
#           'kepler_core',
           'operational',
#           'legacy',
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




