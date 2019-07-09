#### Import all subpackages
import numpy as np
import pandas as df
import datetime as dt
import time
import os 
import string
import re
import struct
import math
import sys
import scipy

from . import extern
from . import utils
from . import athmo
from . import conv

__all__ = ['athmo','utils','extern','conv']
#### Import extern libraires in in the geodezyx namespace
#from geodezyx.extern import *

# NB : if you want extern in the main namespace 
#      don't call with 
#      import geodezyx
#      but with
#      from geodezyx import *

