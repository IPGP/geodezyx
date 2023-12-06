#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 17:46:14 2023

@author: psakic
"""

import re
from io import StringIO
import pandas as pd 
import numpy as np


def read_ringo_systems(file_in):
    
    F = open(file_in)
    L = F.readlines()
    
    start_pattern = re.compile(r"# (\w+)$")
    end_pattern = re.compile(r"# -{35,}|\Z") 
    
    # Initialize variables to track the current satellite system and its data
    current_system = None
    current_data = []
    
    # Dictionary to store extracted tables
    tables = {}      
            
    # Iterate through the lines of the text
    for line in L:
        # Check if the line matches the start pattern
        start_match = start_pattern.match(line)
        if start_match:
            # Initialize data for the new satellite system
            current_system = start_match.group(1)
        # Check if the line matches the end pattern
        elif end_pattern.match(line):
            # Store the current data for the last satellite system
            if current_system:
                tables[current_system] = current_data
                current_system = None
                current_data = []
        elif current_system:
            # If neither start nor end pattern matches, append the line to the current data
            current_data.append(line)
        else:
            pass
            
    for sys,data in tables.items():
        df = pd.read_csv(StringIO("\n".join(data)),delim_whitespace=True)
        tables[sys] = df.rename(columns={"%":"prn"})
        
    return tables


file_in="/home/psakicki/Downloads/qc_HOUE_20230990000"


def read_ringo_qc(file_in):
    F = open(file_in)
    L = F.readlines()
    
    # Define patterns for identifying the start and end of each table
    end_pattern = re.compile(r"# -{35,}|\Z")  # Updated pattern to also match the end of the input
    
    # Define patterns for the "Summary of each satellite system" and "For each satellite" tables
    summary_pattern = re.compile(r"# Summary of each satellite system")
    for_each_satellite_pattern = re.compile(r"# For each satellite")
    
    current_data = []
    
    # Dictionary to store extracted tables
    tables = {}
    
    # Flag to determine if the current section is within "Summary of each satellite system" or "For each satellite"
    in_summary_section = False
    in_for_each_satellite_section = False
    
    # Iterate through the lines of the text
    for iline,line in enumerate(L):
        # Check if the line matches the start pattern
        
        # Check if the line matches the "Summary of each satellite system" pattern
        if summary_pattern.match(line):
            in_summary_section = True
            current_data = []
    
        # Check if the line matches the "For each satellite" pattern
        elif for_each_satellite_pattern.match(line):
            in_for_each_satellite_section = True
            current_data = []
    
        # Check if the line matches the end pattern
        elif end_pattern.match(line) or iline+1 == len(L):
            # If the current section is within "Summary of each satellite system," save the data in the tables dictionary
            if in_summary_section:
                tables["qc_sys"] = current_data
                in_summary_section = False
            # If the current section is within "For each satellite," save the data in the tables dictionary
            elif in_for_each_satellite_section:
                tables["qc_prn"] = current_data
                in_for_each_satellite_section = False
    
        # If the current line contains data, append it to the current_data list
        else:
            current_data.append(line.strip())
            
        
    ### format the table
    for qc,data in tables.items():
        # interpret it as a DataFrame
        df = pd.read_csv(StringIO("\n".join(data)),sep=",")
        df.replace('----',np.nan,regex=True,inplace=True)
        
        # rename the columns
        # sor stands for slips/observation ratio
        cols = [ qc[3:], 'std_mp1', 'std_mp2', 'std_mp5', 'sor_mp1',
                'sor_mp2', 'sor_mp5', 'sor_gf', 'sor_mw', ' sor_iod']
        df.columns = cols
        
        
        
        for col in cols[4:]:
            # remove the ratio with / and convert it for a list [slips , nobs]
            df[col] = df[col].str.split("/")
            df[col] = df[col].apply(lambda x: [int(x[0]),int(x[1])] if int(x[1]) != 0 else [np.nan]*2)
            so_stk = np.stack(df[col],dtype=float)
            # extract slips  an nobs in dedicated df columns
            df["slips_"+col.split("_")[1]] = so_stk[:,0]
            df["nobs_"+col.split("_")[1]] = so_stk[:,1]
            # finaly repalce the sor columns with the actual ratio
            df[col] = df[col].apply(lambda x: np.divide(*x))
            
        tables[qc] = df
        
    
    return tables


T = read_ringo_qc(file_in)

T['qc_sys']