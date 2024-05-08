#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 17:46:14 2023

@author: psakic
"""

#### Import the logger
import logging
import os
import re
from io import StringIO

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

#                       _     _     
#     /\               | |   (_)    
#    /  \   _ __  _   _| |__  _ ___ 
#   / /\ \ | '_ \| | | | '_ \| / __|
#  / ____ \| | | | |_| | |_) | \__ \
# /_/    \_\_| |_|\__,_|_.__/|_|___/
                        

def read_anubis_xtr_sum(xtr_in):
    
    if os.path.basename(xtr_in).split(".")[-1]  != "xtr":
        log.warn("%s is not a XTR file",xtr_in)
    
    L = open(xtr_in).readlines()
    
    def _date_format(linein,end_epoch=False):
        ltmp = list(linein)
        ltmp[18] = "T"
        if end_epoch:
            ltmp[38] = "T"
        return "".join(ltmp)
    
    sum_pattern = re.compile(r"^(#|=)...SUM ")
    
    totsum_stk = []
    gnssum_stk = []
    
    for line in L:
        sum_match = sum_pattern.match(line)
        if sum_match:
            if "TOTSUM" in line:
                totsum_stk.append(line)
            else:
                gnssum_stk.append(line)
                
        if "Header information" in line:
            break
    
    ########### GNSSUM
    ## remove 1st line of GNSSUM
    gnssum_stk.remove(gnssum_stk[0])
    
    for il,l in enumerate(gnssum_stk):
        gnssum_stk[il] = _date_format(l,False)
    
    df_gnssum = pd.read_csv(StringIO("\n".join(gnssum_stk)),
                            delim_whitespace=True,
                            na_values="-")
    
    df_gnssum.rename({"#GNSSUM":"sys",
                      df_gnssum.columns[1]:"Epoch"},axis=1,inplace=True)
    df_gnssum["sys"] = df_gnssum["sys"].apply(lambda x:x[1:4])
    
    for epocol in ["Epoch"]:
        df_gnssum[epocol] = pd.to_datetime(df_gnssum[epocol])
        
    
    
    ########### TOTSUM
    totsum_stk[1] = _date_format(totsum_stk[1],True)
    
    df_totsum = pd.read_csv(StringIO("\n".join(totsum_stk)),
                            delim_whitespace=True)
    
    for epocol in ["First_Epoch________","Last_Epoch_________"]:
        df_totsum[epocol] = pd.to_datetime(df_totsum[epocol])
        
    df_totsum.drop("#TOTSUM",axis=1,inplace=True)
    
    return df_totsum , df_gnssum


#  _____  _                   
# |  __ \(_)                  
# | |__) |_ _ __   __ _  ___  
# |  _  /| | '_ \ / _` |/ _ \ 
# | | \ \| | | | | (_| | (_) |
# |_|  \_\_|_| |_|\__, |\___/ 
#                  __/ |      
#                 |___/   

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