#### Import all subpackages
from . import extern
from . import utils
from . import athmo
import numpy as np
import pandas as pd
__all__ = ['athmo','utils','extern']
#### Import extern libraires in in the geodezyx namespace
#from geodezyx.extern import *

# NB : if you want extern in the main namespace 
#      don't call with 
#      import geodezyx
#      but with
#      from geodezyx import *

