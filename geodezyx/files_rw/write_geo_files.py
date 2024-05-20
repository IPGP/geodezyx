#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.files_rw contains functions to 
write misc. geodetic data in dedicated files.

it can be imported directly with:
from geodezyx import files_rw

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""

########## BEGIN IMPORT ##########
#### External modules
import datetime as dt
#### Import the logger
import logging
import os
import re

import numpy as np
import pandas as pd

#### geodeZYX modules
from geodezyx import conv
from geodezyx import files_rw
from geodezyx import reffram
from geodezyx import utils

# from geodezyx.megalib.megalib import *   # Import the legacy modules names
log = logging.getLogger(__name__)


#### Import star style
# from geodezyx import *                   # Import the GeodeZYX modules
# from geodezyx.externlib import *         # Import the external modules
# from geodezyx.megalib.megalib import *   # Import the legacy modules names
##########  END IMPORT  ##########
    
def write_sp3(SP3_DF_in,outpath,outname=None,prefix='orb',
              skip_null_epoch=True,force_format_c=False):
    """
    Write a SP3 file from an Orbit DataFrame

    Parameters
    ----------
    SP3_DF_in : DataFrame
        Input Orbit DataFrame.
    outpath : str
        The output path of the file (see also outname).
    outname : None or str, optional
        None = outpath is the full path (directory + filename) of the output.
        A string = a manual name for the file.
        'auto_old_cnv' = automatically generate the filename (old convention)
        'auto_new_cnv' = automatically generate the filename (new convention)
        The default is None.
    prefix : str, optional
        the output 3-char. name of the AC. The default is 'orb'.
    skip_null_epoch : bool, optional
        Do not write an epoch if all sats are null (filtering). 
        The default is True.
    force_format_c : bool, optional
        force SP3's format c. The default is False.

    Returns
    -------
    The string containing the formatted SP3 data.
    """
    
    ################## MAIN DATA
    LinesStk = []

    SP3_DF_wrk = SP3_DF_in.sort_values(["epoch","prn"])

    EpochRawList  = SP3_DF_wrk["epoch"].unique()
    SatList    = sorted(SP3_DF_wrk["prn"].unique())
    SatList    = list(sorted(SP3_DF_wrk["prn"].unique()))
    ## SatList    = list(reversed(SatList)) 
    #### PS 210721
    #### reversed bc the sats are sorted ascending=False, but why???
    #### to have G before E ??
    SatListSet = set(SatList)
    EpochUsedList = []
    
    if not "clk" in SP3_DF_wrk.columns:
        SP3_DF_wrk["clk"] = 999999.999999
    
    for epoc in EpochRawList:
        SP3epoc   = pd.DataFrame(SP3_DF_wrk[SP3_DF_wrk["epoch"] == epoc])
        
        ######## if keep_missing_sat_in_epoch:
        ## manage missing sats for the current epoc
        MissingSats = SatListSet.difference(set(SP3epoc["prn"]))
        
        for miss_sat in MissingSats:
            miss_line = SP3epoc.iloc[0].copy()
            miss_line["prn"]   = miss_sat
            miss_line["sys"] = miss_sat[0]
            ### check the sp3 doc 
            # bad position = 0.000000
            # bad clock    = 999999.9999999
            miss_line["x"]     = 0.000000
            miss_line["y"]     = 0.000000
            miss_line["z"]     = 0.000000
            miss_line["clk"]   = 999999.999999
            
            SP3epoc = SP3epoc.append(miss_line)
        #### end of missing sat bloc

        SP3epoc.sort_values("prn",inplace=True,ascending=True)
        timestamp = conv.dt2sp3_timestamp(conv.numpy_dt2dt(epoc)) + "\n"

        linefmt = "P{:}{:14.6f}{:14.6f}{:14.6f}{:14.6f}\n"

        LinesStkEpoch = []
        sum_val_epoch = 0
        for ilin , lin in SP3epoc.iterrows():
            if not "clk" in lin.index:  # manage case if no clk in columns
                lin["clk"] = 999999.999999
            line_out = linefmt.format(lin["prn"],lin["x"],lin["y"],lin["z"],lin["clk"])
            
            sum_val_epoch += lin["x"]+lin["y"]+lin["z"]

            LinesStkEpoch.append(line_out)


        ### if skip_null_epoch activated, print only if valid epoch 
        if not ( np.isclose(sum_val_epoch,0) and skip_null_epoch):
            LinesStk.append(timestamp)          # stack the timestamp
            LinesStk = LinesStk + LinesStkEpoch # stack the values
            EpochUsedList.append(epoc)          # stack the epoc as dt


    ################## HEADER
    ######### SATELLITE LIST

    Satline_stk   = []
    Sigmaline_stk = []


    if force_format_c:
        nlines = 5
    else:
        div,mod = np.divmod(len(SatList),17)
        
        if div < 5:
            nlines = 5
        else:
            nlines = div

            if mod != 0:
                nlines += 1
        
        
    for i in range(nlines):
        SatLine = SatList[17*i:17*(i+1)]
        SatLineSigma = len(SatLine) * " 01"
        
        if len(SatLine) < 17:
            complem = " 00" * (17 - len(SatLine))
        else:
            complem = ""

        if i == 0:
            nbsat4line = len(SatList)
        else:
            nbsat4line = ''

        satline = "+  {:3}   ".format(nbsat4line) + "".join(SatLine) + complem + "\n"
        sigmaline = "++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n"
        sigmaline = "++       " + SatLineSigma + complem  + "\n"
        
        Satline_stk.append(satline)
        Sigmaline_stk.append(sigmaline)


    ######### 2 First LINES
    start_dt = conv.numpy_dt2dt(np.min(EpochUsedList))
    
    header_line1 = "#cP" + conv.dt2sp3_timestamp(start_dt,False) + "     {:3}".format(len(EpochUsedList)) + "   u+U IGSXX FIT  XXX\n"

    delta_epoch = int(utils.most_common(np.diff(EpochUsedList) * 10**-9))
    MJD  = conv.dt2MJD(start_dt)
    MJD_int = int(np.floor(MJD))
    MJD_dec = MJD - MJD_int
    gps_wwww , gps_sec = conv.dt2gpstime(start_dt,False,"gps")

    header_line2 = "## {:4} {:15.8f} {:14.8f} {:5} {:15.13f}\n".format(gps_wwww,gps_sec,delta_epoch,MJD_int,MJD_dec)


    ######### HEADER BOTTOM
    header_bottom = """%c M  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
%c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
%f  1.2500000  1.025000000  0.00000000000  0.000000000000000
%f  0.0000000  0.000000000  0.00000000000  0.000000000000000
%i    0    0    0    0      0      0      0      0         0
%i    0    0    0    0      0      0      0      0         0
/* PCV:IGSXX_XXXX OL/AL:FESXXXX  NONE     YN CLK:CoN ORB:CoN
/*     GeodeZYX Toolbox Output
/*
/*
"""


    ################## FINAL STACK

    FinalLinesStk = []

    FinalLinesStk.append(header_line1)
    FinalLinesStk.append(header_line2)
    FinalLinesStk = FinalLinesStk + Satline_stk + Sigmaline_stk
    FinalLinesStk.append(header_bottom)
    FinalLinesStk = FinalLinesStk + LinesStk + ["EOF"]

    FinalStr = "".join(FinalLinesStk)


    ### Manage the file path
    prefix_opera = prefix
    
    if not outname:
        outpath_opera = outpath
    elif outname == 'auto_old_cnv':
        week , dow = conv.dt2gpstime(start_dt)
        filename = prefix_opera + str(week) + str(dow) + '.sp3'
        outpath_opera = os.path.join(outpath,filename)
        
    elif outname == 'auto_new_cnv':
        log.error("not implemented yet !!!!!")
        raise Exception
        
    F = open(outpath_opera,"w+")
    F.write(FinalStr)
    
    
def write_clk(DFclk_in,outpath,
              outname=None,prefix='orb',
              header="",output_std_values=False):
    """
    Write a SP3 Clock file from an Clock DataFrame

    Parameters
    ----------
    DFclk_in : DataFrame
        Input Clock DataFrame.
    outpath : str
        The output path of the file (see also outname).
    outname : None or str, optional
        None = outpath is the full path (directory + filename) of the output.
        A string = a manual name for the file.
        'auto_old_cnv' = automatically generate the filename (old convention)
        'auto_new_cnv' = automatically generate the filename (new convention)
        The default is None.
    prefix : str, optional
        the output 3-char. name of the AC. The default is 'orb'.
    header : str, optional
        A string describing the clk file header. The default is "".
    output_std_values : bool, optional
        Add observation sigmas as the last column. The default is False.

    Returns
    -------
    The string containing the formatted clock data.
    """
    
    HEAD = header
    Row_str_stk = []

    if output_std_values:
        row_str_proto = "{:2} {:4} {:4d} {:02d} {:02d} {:02d} {:02d} {:9.6f} {:2d}   {:19.12e} {:19.12e}"
    else:
        row_str_proto = "{:2} {:4} {:4d} {:02d} {:02d} {:02d} {:02d} {:9.6f} {:2d}   {:19.12e}"
        
    for irow, row in DFclk_in.iterrows():

        if output_std_values:
            one_or_two=2
            row_str = row_str_proto.format(row["type"],row["name"],
                                           int(row["year"]),
                                           int(row["month"]),
                                           int(row["day"]),
                                           int(row["hour"]),
                                           int(row["minute"]),
                                           int(row["second"]),
                                           one_or_two,
                                           row["bias"],
                                           row["sigma"])
        else:
            one_or_two=1
            row_str = row_str_proto.format(row["type"],row["name"],
                                           int(row["year"]),
                                           int(row["month"]),
                                           int(row["day"]),
                                           int(row["hour"]),
                                           int(row["minute"]),
                                           int(row["second"]),
                                           one_or_two,
                                           row["bias"])            
        Row_str_stk.append(row_str)
    
    ## Add EOF
    Row_str_stk.append("EOF")
    
    CORPSE = "\n".join(Row_str_stk)
       
    OUT = HEAD + "                                                            END OF HEADER\n" + CORPSE
    
    ### Manage the file path
    prefix_opera = prefix
    
    if not outname:
        outpath_opera = outpath
    elif outname == 'auto_old_cnv':
        start_dt = dt.datetime(int(DFclk_in.iloc[0]["year"]),
                               int(DFclk_in.iloc[0]["month"]),
                               int(DFclk_in.iloc[0]["day"]))
        week , dow = conv.dt2gpstime(start_dt)
        filename = prefix_opera + str(week) + str(dow) + '.clk'
        outpath_opera = os.path.join(outpath,filename)
        
    elif outname == 'auto_new_cnv':
        log.error("not implemented yet !!!!!")
        raise Exception
        
    else:
        outpath_opera = os.path.join(outpath,outname)
        
    OUT = OUT   

    with open(outpath_opera,"w+") as Fout:
        Fout.write(OUT)
        Fout.close()
    
        
    return OUT


def ine_block_mono(sat,dt_in,extra_intrvl_strt=.1,extra_intrvl_end=.4,step=300):
    """
    Write an EPOS INE block
    """
    
    Fields = ['orb____1',
    'orb____2',
    'orb____3',
    'orb____4',
    'orb____5',
    'orb____6',
    'orb___db',
    'orb_s2db',
    'orb_c2db',
    'orb_s4db',
    'orb_c4db',
    'orb___yb',
    'orb___xb',
    'orb_sixb',
    'orb_coxb',
    'orb___cr']
    
    
    mjd = np.floor(conv.dt2MJD(dt_in))
    mjd_strt = mjd - extra_intrvl_strt
    mjd_end  = mjd + extra_intrvl_end + 1
    
    Lines = []
    
    l1 =  " sat_nr  : " + sat + "\n"
    l2 =  " stepsize: {:3}  {:6.2f}\n".format(sat,step)
    
    Lines.append(l1)
    Lines.append(l2)
    
    
    for field in Fields:
        line = " {:}: {:3}  0.000000000000000E+00 {:11.5f} {:11.5f}\n".format(field,sat,mjd_strt,mjd_end)
        Lines.append(line)
        
    Lines.append(" end_sat\n")
        
    str_out = "".join(Lines)
    
    return str_out



def write_ine_dummy_file(Sat_list,dt_in,extra_intrvl_strt=.1,
             extra_intrvl_end=.4,step=300,out_file_path=None):
    """
    Write an EPOS INE dummy (empty values) file
    """

    Lines = []
    
    mjd = np.floor(conv.dt2MJD(dt_in))
    mjd_strt = mjd - extra_intrvl_strt
    mjd_end  = mjd + extra_intrvl_end + 1
    
    datestr = conv.dt2str(dt.datetime.now(),str_format='%Y/%m/%d %H:%M:%S')
    
    mjd_strt_deci = mjd_strt - np.floor(mjd_strt)
    
    
    head_proto="""%=INE 1.00 {:} NEWSE=INE+ORBCOR                                                                                 
+global
 day_info: 
 epoch   :                            {:5}  {:16.14f}
 interval:                            {:11.5f} {:11.5f}
 stepsize:      {:6.2f}
-global
+initial_orbit
"""
    head = head_proto.format(datestr,int(mjd),0,mjd_strt,mjd_end,step)
    
    Lines.append(head)
    
    for sat in Sat_list:
        Lines.append("******************************************************************\n")
        sat_str = ine_block_mono(sat,dt_in,extra_intrvl_strt,extra_intrvl_end,step)
        Lines.append(sat_str)
        Lines.append("******************************************************************\n")
    
    str_end = """-initial_orbit
%ENDINE
"""
    
    Lines.append(str_end)
         
    str_out = "".join(Lines)
    
    if out_file_path:
        with open(out_file_path,"w") as f:
            f.write(str_out)
            f.close()

    return str_out

def sp3_overlap_creator(ac_list,dir_in,dir_out,
                        suffix_out_input = None,
                        overlap_size = 7200,
                        force = False,
                        manage_missing_sats='exclude_missing_epoch',
                        eliminate_null_sat=True,
                        severe=False,
                        separated_systems_export=False,
                        first_date=None,
                        end_date=None,
                        #new_naming = False,
                        exclude_bad_epoch=True,
                        sys = None,
                        new_name = False):
    """
    Generate an SP3 Orbit file with overlap based on the SP3s of the 
    days before and after
    
    Parameters
    ----------
    ac_list : list
        3-character codes of the ACs.
    dir_in : str
        where the input sp3 are.
    dir_out : str
         where the output sp3 will be outputed.
    suffix_out_input : str, optional
        last char of the 3-char. code. if None, then it is the same as input.
    overlap_size : int, optional
        Overlapsize. The default is 7200.
    force : True, optional
        force overwrite. The default is False.
    manage_missing_sats : str, optional
        'exclude_missing_day' : generate a file with only the common sat 
        between the 3 days. Thus, exclude the missing sats for a complete day\n
        'exclude_missing_epoch' : generate a file with only sat with full epochs\n
        'extrapolate' : extrapolate the missing sats based on the first/last epoch\n
        The default is 'exclude_missing_epoch'.
    eliminate_null_sat : bool, optional
        eliminate null sat. The default is True.
    severe : bool, optional
        raise an exception if problem. The default is False.
    separated_systems_export : bool, optional
        export different sp3 for different system. The default is False.
    first_date : datetime, optional
        exclude SP3 before this epoch
    end_date :datetime, optional
        exclude SP3 after this epoch
    exclude_bad_epoch : bool, optional
        remove bad epoch (usually filled with 99999999.9999999 or 0.00000000)
    sys : string, optional
        if just to keep one system (e.g.: G)
    
    Returns
    -------
    None.

    Note
    ----
    start/end date are not implemented
    the force option skips existing files 
    not implemented for new names
    """


    Dict_Lfiles_ac = dict()

    for ac in ac_list:
        Dict_Lfiles_ac[ac] = []
        Lfile = Dict_Lfiles_ac[ac]
        

        if new_name:
            Lfile = softs_runner.find_IGS_products_files(dir_in,["sp3"],
                                                                  ac_list,
                                                                  first_date,
                                                                  date_end=end_date,
                                                                  severe=False,
                                                                  recursive_search=True,
                                                                  regex_old_naming = True,
                                                                  regex_new_naming = True,
                                                                  regex_igs_tfcc_naming = False,
                                                                  compressed="incl") 
            
        else:
            #Extlist = ["sp3","SP3","sp3.gz","SP3.gz"]
            Extlist = ["sp3","sp3.gz",'eph',"SP3","SP3.gz"]
            for ext in Extlist:
                Lfile = Lfile + utils.find_recursive(dir_in,"*" + ac + "*" + ext)
                
        log.info("Nb of SP3 found for %s %s",ac,len(Lfile))
        
        if not suffix_out_input:
            suffix_out = ac
        else:
            suffix_out = ac[:2] + suffix_out_input ## GM  2021-09-06 this is not good for the new namning convention
        
        D     = []
        WWWWD = []
        
        New_name_list = []
        
        regex_suffix = "[0-9]{4}_[0-9]{2}[A-Z]_[0-9]{2}[A-Z]_ORB.SP3"
         
        for sp3 in Lfile:
            #wwwwd_str = os.path.basename(sp3)[3:8]
            #D.append(conv.gpstime2dt(int(wwwwd_str[:4]),int(wwwwd_str[4:])))
            extension = sp3[-3:]
            dat = conv.sp3name2dt(sp3)
            D.append(dat)
        
            if re.search(regex_suffix,sp3):
                New_name_list.append(True)
            else:
                New_name_list.append(False)
                
       
#        ## just a check for the files list. May be excluded in the future. 
#        for item in Lfile:
#            if 'COD0MGXFIN_20200050000_01D_05M_ORB.SP3' in item:
#                print(Lfile.index(item))
       
            
        for dat,newnamebool in zip(D[1:-1],New_name_list[1:-1]): ####if selection manuel, zip > 2lists !!!
#        for dat in D[875:876]: ####if selection manuel, zip > 2lists !!!

            dummy= 1
            try:
                log.info("*********** %s %s",ac,dat)
                
                if first_date:
                    if dat < first_date:
                        log.info("SKIP date",dat)
                        continue
                if end_date:
                    if dat > end_date:
                        log.info("SKIP date",dat)
                        continue
                        
                    
                wwwwd_str = conv.dt_2_sp3_datestr(dat).zfill(5)
            
                dat_bef = dat - dt.timedelta(days=1)
                dat_aft = dat + dt.timedelta(days=1)
                
                wwwwd_str_bef = utils.join_improved("",*conv.dt2gpstime(dat_bef)).zfill(5)
                wwwwd_str_aft = utils.join_improved("",*conv.dt2gpstime(dat_aft)).zfill(5)
                
                ###### check if exists
                dir_out_wk = os.path.join(dir_out,"wk" + str(wwwwd_str)[:4])
                utils.create_dir(dir_out_wk)
                fil_out = dir_out_wk + "/" + suffix_out  + wwwwd_str + ".sp3"
                
                if not force and os.path.isfile(fil_out):
                    log.info("0)) " + fil_out + " exists, skipping...")
                    continue


                ### *************** STEP 1 ***************
                log.info("1)) Search for the days before/after")                
                log.info("1)) %s %s",dat_bef,dat_aft)
    
    
                ## GM 2021-09-06 in case of the new naming conventation this is step needs to be done:
                ## not the best way to handle this, since it is hardcoded and the step needs to be given in a "if"
                ## maybe try to find a better way with regex!
                # if new_naming:
                #     if ac in ['WUM','GRG','SHA']:
                #         step_str = '15'
                #     else:
                #         step_str = '05'
                        
                #     year,day = conv.dt_to_doy(conv.gpstime2dt(int(wwwwd_str[:4]),int(wwwwd_str[-1])))                    
                #     p1    = utils.find_regex_in_list(str(year)+str(day).zfill(3)  + "0000_01D_"+step_str+"M_ORB.SP3",Lfile,True)
                    
                #     year_bef,day_bef = conv.dt_to_doy(conv.gpstime2dt(int(wwwwd_str_bef[:4]),int(wwwwd_str_bef[-1])))
                #     p_bef = utils.find_regex_in_list(str(year_bef)+str(day_bef).zfill(3)  + "0000_01D_"+step_str+"M_ORB.SP3",Lfile,True)
                    
                #     year_aft,day_aft = conv.dt_to_doy(conv.gpstime2dt(int(wwwwd_str_aft[:4]),int(wwwwd_str_aft[-1])))
                #     p_aft = utils.find_regex_in_list(str(year_aft)+str(day_aft).zfill(3)  + "0000_01D_"+step_str+"M_ORB.SP3",Lfile,True)
                    
                #if re.search(regex_suffix,) 
                if newnamebool:
                    day,year = conv.dt2doy_year(dat)                    
                    regex_prefix = str(year)+str(day).zfill(3)
                    p1    = utils.find_regex_in_list(regex_prefix + regex_suffix,Lfile,True)

                    day_bef,year_bef = conv.dt2doy_year(dat_bef)      
                    regex_prefix_bef = str(year_bef)+str(day_bef).zfill(3)
                    p_bef = utils.find_regex_in_list(regex_prefix_bef + regex_suffix ,Lfile,True)
                    
                    day_aft,year_aft = conv.dt2doy_year(dat_aft)      
                    regex_prefix_aft = str(year_aft)+str(day_aft).zfill(3) 
                    p_aft = utils.find_regex_in_list(regex_prefix_aft + regex_suffix ,Lfile,True)
                    
                else: 
                    if extension == "eph":
                        p1    = utils.find_regex_in_list(wwwwd_str     + ".eph",Lfile,True)
                        p_bef = utils.find_regex_in_list(wwwwd_str_bef + ".eph",Lfile,True)
                        p_aft = utils.find_regex_in_list(wwwwd_str_aft + ".eph",Lfile,True)
                    else:
                        p1    = utils.find_regex_in_list(wwwwd_str     + ".sp3",Lfile,True)
                        p_bef = utils.find_regex_in_list(wwwwd_str_bef + ".sp3",Lfile,True)
                        p_aft = utils.find_regex_in_list(wwwwd_str_aft + ".sp3",Lfile,True)

                log.info("1)) Files found for the days before/after")                            
                log.info("0b) %s",p_bef)
                log.info("01) %s",p1)
                log.info("0a) %s",p_aft)
            
                if not p1 or not p_bef or not p_aft:
                    log.error("with day %s",dat)
                    continue
                
                SP3     = files_rw.read_sp3(p1)
                SP3_bef = files_rw.read_sp3(p_bef)
                SP3_aft = files_rw.read_sp3(p_aft)
                
                ### Filtering to keep P only
                SP3 = SP3[SP3.type == "P"]
                SP3_bef = SP3_bef[SP3_bef.type == "P"]
                SP3_aft = SP3_aft[SP3_aft.type == "P"]
                
                if sys:
                    SP3     = SP3[SP3.sys == sys ]
                    SP3_bef = SP3_bef[SP3_bef.sys == sys ]
                    SP3_aft = SP3_aft[SP3_aft.sys == sys ]
                
                
                SP3_bef = SP3_bef[SP3_bef["epoch"] < SP3["epoch"].min()]
                SP3_aft = SP3_aft[SP3_aft["epoch"] > SP3["epoch"].max()]
                
                SP3concat = pd.concat((SP3_bef,SP3,SP3_aft))
                
                dat_filter_bef = dat - dt.timedelta(seconds=overlap_size)
                dat_filter_aft = dat + dt.timedelta(seconds=overlap_size) + dt.timedelta(days=1)

                ### *************** STEP 2 ***************
                log.info("2)) dates of the overlap period before/after")                   
                log.info("2)) %s %s",dat_filter_bef,dat_filter_aft)

                ### *************** STEP 3 *************** 
                log.info("3)) Dates of: SP3 concatenated, before, current, after")                       
                log.info("3)) %s %s",SP3concat["epoch"].min(),SP3concat["epoch"].max())
                log.info("3b) %s %s",SP3_bef["epoch"].min(),SP3_bef["epoch"].max())
                log.info("31) %s %s",SP3["epoch"].min(),SP3["epoch"].max())
                log.info("3a) %s %s",SP3_aft["epoch"].min(),SP3_aft["epoch"].max())
                
                SP3concat = SP3concat[(SP3concat["epoch"] >= dat_filter_bef) & (SP3concat["epoch"] <= dat_filter_aft)]
                
                if exclude_bad_epoch:
                    Good_epochs_bool_999 = SP3concat[['x','y','z']] < 999999.
                    Good_epochs_bool_000 = np.logical_not(np.isclose(SP3concat[['x','y','z']],0.))
                    
                    Good_epochs_bool = np.logical_and(Good_epochs_bool_999,Good_epochs_bool_000)
                    Good_epochs_bool = np.all(Good_epochs_bool,axis=1)
                                                            
                    SP3concat = SP3concat[Good_epochs_bool]
                
                ########## HERE WE MANAGE THE MISSING SATS
                if manage_missing_sats == "exclude_missing_day":     
                    log.info("4))","remove missing sats -- day")                                     
                    common_sats = set(SP3_bef["prn"]).intersection(set(SP3["prn"])).intersection(set(SP3_aft["prn"]))
                    SP3concat = SP3concat[SP3concat["prn"].isin(common_sats)]
                    
                elif manage_missing_sats == "exclude_missing_epoch":
                    log.info("4))","remove missing sats -- epoch")      
                    nepoc = len(SP3concat["epoch"].unique())
                    SP3concat_satgrp = SP3concat.groupby("prn")
                    
                    All_sats = SP3concat["prn"].unique()
                    Good_sats = SP3concat_satgrp.count() == nepoc

                    ###### Good_sats = Good_sats.reset_index()["sat"]
                    ## we get the good sats based one column containing a boolean
                    ## because of the test just before (abitrarily epoch column)
                    ## and after get the corresponding good sats names
                    Good_sats = Good_sats[Good_sats["epoch"]].reset_index()["prn"]
                    
                    Bad_sats = list(set(All_sats) - set(Good_sats))
                    log.info("excluded bad sats: %s", Bad_sats)
                    
                    SP3concat = SP3concat[SP3concat["prn"].isin(Good_sats)]
                    
                
                elif manage_missing_sats == "extrapolate":
                    log.info("4)) extrapolate missing sats ")                                     
                    for iovl,SP3_ovl in enumerate((SP3_bef,SP3_aft)):
                        if iovl == 0:
                            backward = True
                            forward  = False
                            backfor = "backward"
                        elif iovl == 1:
                            backward = False
                            forward  = True
                            backfor = "forward"
                            
                        Sats = set(SP3["prn"])
                        Sats_ovl = set(SP3_ovl["prn"])
                    
                        Sats_miss = Sats.difference(Sats_ovl)
                        if not Sats_miss:
                            continue
                        log.info("4a) extrapolate missing sats %s %s",backfor,Sats_miss)                                     

                        SP3extrapo_in = SP3concat[SP3concat["prn"].isin(Sats_miss)]
                        
                        #step = utils.most_common(SP3concat["epoch"].diff().dropna())
                        #step = step.astype('timedelta64[s]').astype(np.int32)
                        step = 900
                        #print(step)
                        
                        #print("SP3extrapo_in",SP3extrapo_in)
                        
                        SP3extrapo = reffram.extrapolate_sp3_DataFrame(SP3extrapo_in,
                                                                       step=step,
                                                                       n_step=int(overlap_size/step),
                                                                       backward=backward,
                                                                       forward=forward,
                                                                       until_backward=dat_filter_bef,
                                                                       until_forward=dat_filter_aft,
                                                                       return_all=False)
                        
                        SP3concat = pd.concat((SP3concat,SP3extrapo))
                        log.info(SP3extrapo)

                else:
                    log.error("check manage_missing_sats value")
                    raise Exception
                    
                if eliminate_null_sat:
                    GoodSats = []
                    for sat in SP3concat["prn"].unique():
                        XYZvals = SP3concat[SP3concat["prn"] == sat][["x","y","z"]].sum(axis=1)
                        
                        V = np.sum(np.isclose(XYZvals,0)) / len(XYZvals)
                                            
                        if V < 0.50:
                            GoodSats.append(sat)
                        else:
                            log.info("6) eliminate because null position %s",sat)
                        
                    SP3concat = SP3concat[SP3concat["prn"].isin(GoodSats)]

                ### *************** STEP 7 ***************           
                log.info("7)) Start/End Epoch of the concatenated file ")                                     
                log.info("7)) %s %s",SP3concat["epoch"].min(),SP3concat["epoch"].max())

                #### All systems        
                log.info("8)) outputed file")
                log.info(fil_out)
                write_sp3(SP3concat,fil_out)
                
                #### system separated
                if False:
                    for sys in SP3concat["sys"].unique():
                        try:
                            SP3concat_sys = SP3concat[SP3concat["sys"] == sys]
                            fil_out_sys = dir_out_wk + "/" + suffix_out[:2] + sys.lower() + wwwwd_str.zfill(5) + ".sp3"
                            log.info("9)) outputed file")
                            log.info(fil_out_sys)
                            write_sp3(SP3concat_sys,fil_out_sys)
                        except:
                            continue
            
            except KeyboardInterrupt:
                raise KeyboardInterrupt
                
            except Exception as e:
                if severe:
                    log.error(e)
                    raise e
                else:
                    log.warning("Error %s but no severe mode, continue...",e)
                    continue


    """
    sort_wrt="site" or "site_num"
    
    soln_in_DF
    use soln AND pt information in the input DataFrame
    """
    


def write_epos_sta_coords(DF_in,file_out,sort_wrt="site",
                          no_time_limit_for_first_period = True,
                          no_time_limit_for_last_period = True,
                          soln_in_DF=True,
                          TRF_name="xTRFnn"):
    """
    Write an EPOS coordinate file

    Parameters
    ----------
    DF_in : DataFrame
        Input Orbit DataFrame.
    file_out : str
        The output path of the file.
    sort_wrt : bool, optional
        Sort the values with respect to a DF column. 
        The default is "site".
    no_time_limit_for_first_period : bool, optional
        No time limit for the first period. 
        The default is True.
    no_time_limit_for_last_period : bool, optional
        No time limit for the last period. 
        The default is True.
    soln_in_DF : bool, optional
        Soln in DF. 
        The default is True.

    Returns
    -------
    None.

    """
    
    DF_work = DF_in.sort_values([sort_wrt,"MJD_start"])

    Stat_lines_blk_stk = []

    generic_header = """+info
 FLATTENING                  298.2550
 MAJOR_AXIS              6378140.0000
 REFERENCE_FRAME                IGS14
 NUMBER_OF_STATIONS             {:5d}
 REF_MJD                        {:5d}
-info
"""

    generic_header = generic_header.format(len(DF_work["site_num"].unique()),
                                           int(utils.most_common(DF_work["MJD_ref"])))

    Stat_lines_blk_stk.append(generic_header)

    Stat_lines_blk_stk.append("+station_coordinates")

    for site in DF_work[sort_wrt].unique():

        Stat_lines_blk_stk.append("*------------------------- ---- ----- -beg- -end- -**- ------------------------------------------------\n*")

        DF_SiteBlock = DF_work[DF_work[sort_wrt] == site]
        
        DF_SiteBlock.reset_index(inplace=True)

        for i_l ,(_ , l) in enumerate(DF_SiteBlock.iterrows()):

            if soln_in_DF:
                iope = int(l["soln"])
                pt = l["pt"]
            else:
                iope = i_l + 1
                pt = "A"
            
            if no_time_limit_for_first_period and i_l == 0:
                MJD_start = 0
            else:
                MJD_start = l["MJD_start"]
                                            
            if no_time_limit_for_last_period and (i_l+1) == len(DF_SiteBlock):
                MJD_end = 0
            else:
                MJD_end = l["MJD_end"]
                

            line_site_fmt = " SITE            m {:4d}  {:1d} {:} {:5d} {:5d} {:5d} {:}   {:}  {:1d}      LOG_CAR       LOG_CAR"
            line_site_fmt = " SITE            m {:4d}  {:1d} {:} {:5d} {:5d} {:5d} {:}   {:}  {:1d}{:>13}       {:<13}"
            line_valu_fmt = " POS_VEL:XYZ     m {:4d}  {:1d} {:+15.4f} {:+15.4f} {:+15.4f}      {:+6.4f} {:+6.4f} {:+6.4f}"
            line_sigm_fmt = " SIG_PV_XYZ      m {:4d}  {:1d} {:+15.4f} {:+15.4f} {:+15.4f}      {:+6.4f} {:+6.4f} {:+6.4f}"

            line_site = line_site_fmt.format(int(l["site_num"]),
                                             int(iope),
                                             l["tecto_plate"].upper(),
                                             int(l["MJD_ref"]),
                                             int(MJD_start),
                                             int(MJD_end),
                                             l["site"],
                                             pt,
                                             int(iope),
                                             TRF_name,
                                             TRF_name)
            
            line_valu = line_valu_fmt.format(int(l["site_num"]),
                                             int(iope),
                                             l["x"],
                                             l["y"],
                                             l["z"],
                                             l["Vx"],
                                             l["Vy"],
                                             l["Vz"])
            
            line_sigm = line_sigm_fmt.format(int(l["site_num"]),
                                             int(iope),
                                             l["sx"],
                                             l["sy"],
                                             l["sz"],
                                             l["sVx"],
                                             l["sVy"],
                                             l["sVz"])

            Stat_lines_blk_stk.append(line_site)
            Stat_lines_blk_stk.append(line_valu)
            Stat_lines_blk_stk.append(line_sigm)
            Stat_lines_blk_stk.append("*")

    Stat_lines_blk_stk.append("-station_coordinates")

    final_str = "\n".join(Stat_lines_blk_stk)


    with open(file_out,"w+") as f:
        f.write(final_str)

    return final_str



def write_sndy_light_dat(ts_in,outdir,outprefix):
    """ Not properly implemented """
    fil = open(os.path.join(outdir,outprefix),'w+')
    if isinstance(ts_in,TimeSeriePoint):
        if ts_in.initype() == 'FLH':
            for pt in ts_in.pts:
                lin = ' '.join([str(e) for e in [pt.F , pt.L , pt.H , pt.T , pt.sF , pt.sL , pt.sH ]])
                fil.write(lin + '\n')
    elif isinstance(ts_in,TimeSerieObs):
        if ts_in.typeobs == 'RPY':
            for att in ts_in.obs:
                lin = ' '.join([str(e) for e in [att.R , att.P , att.Y , att.T , att.Q.w , att.Q.x , att.Q.y , att.Q.z ]])
                fil.write(lin + '\n')
    fil.close()
