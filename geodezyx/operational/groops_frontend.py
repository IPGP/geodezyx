#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 10:54:49 2023

@author: psakicki
"""

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names


def groops_basic_runner(xml_config_path="",
                        global_var_dict=dict(),
                        xml_var_dict=dict(),
                        quiet=False,
                        groops_bin_path='/home/psakicki/SOFTWARE/GROOPS/groops/bin/groops'):


    global_args_str = ""
    for key, val in global_var_dict.items():
        global_args_str = global_args_str + " ".join(("--global",key+"="+val))
        
        
    command = " ".join((groops_bin_path,global_args_str,xml_config_path))
    log.info("groops command: %s",command)

    process1 = subprocess.run(command,
                              capture_output=True,
                              text=True,
                              shell=True)

groops_basic_runner()


