#### IMPORT EXTERNAL MODULES
import bisect
from collections import defaultdict
from collections import Counter
import copy
import datetime as dt
import fnmatch
import glob
import inspect
import itertools
import math
import matplotlib
import matplotlib.pyplot as plt
from natsort import natsorted, ns
import numpy as np
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
import time
import uuid


#### IMPORT GEODEZYX INTERNAL SUBMODULES
from . import athmo
from . import conv
from . import files_rw
from . import geodyn
#from . import legacy
from . import megalib
from . import reffram
from . import stats
from . import time_series
from . import utils

__all__ = ['athmo',
           'conv',
           'files_rw',
           'geodyn',
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




