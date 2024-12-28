#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#### IMPORT EXTERNAL MODULES
import logging.config # MUST remain logging.config for < v3.10
# (import logging alone doesn't work for < 3.10)
import os
from os import path

__version__='4.5.3'  ## increase it with bump-my-version !!!

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
# from . import externlib # disable 202406xx
# from . import megalib # disable 202411xx
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

if __name__ == '__main__':
    print("geodezyx version",__version__)
    print("geodezyx.__all__",__all__)

#### Import extern libraires in in the geodezyx namespace
#from geodezyx.extern import *

# NB : if you want extern in the main namespace 
#      don't call with 
#      import geodezyx
#      but with
#      from geodezyx import *

