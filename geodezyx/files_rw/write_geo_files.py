#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 18:00:21 2019

@author: psakicki
"""

#### geodeZYX modules
from geodezyx import conv
from geodezyx import operational
from geodezyx import utils

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names
##########  END IMPORT  ##########


def write_sndy_light_dat(ts_in,outdir,outprefix):
    """pas fini"""
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
    
def write_sp3(SP3_DF_in,outpath,outname=None,prefix='orb',
              skip_null_epoch=True,force_format_c=False):
    """
    Write DOCSTRING

    outname : None = outpath is the full path (directory + filename) of the outputfile
             'auto_old_cnv' = automatically generate the filename (old convention)
             'auto_new_cnv' = automatically generate the filename (new convention)
             
    prefix : the output name of the AC
    
    skip_null_epoch: Do not write an epoch if all sats are null (filtering)

    """
    ################## MAIN DATA
    LinesStk = []

    SP3_DF_wrk = SP3_DF_in.sort_values(["epoch","sat"])

    EpochRawList  = SP3_DF_wrk["epoch"].unique()
    SatList    = sorted(SP3_DF_wrk["sat"].unique())
    SatList    = list(reversed(SatList))
    SatListSet = set(SatList)
    EpochUsedList = []
    
    if not "clk" in SP3_DF_wrk.columns:
        SP3_DF_wrk["clk"] = 999999.999999
    
    for epoc in EpochRawList:
        SP3epoc   = pd.DataFrame(SP3_DF_wrk[SP3_DF_wrk["epoch"] == epoc])
        ## Missing Sat
        MissingSats = SatListSet.difference(set(SP3epoc["sat"]))
        
        for miss_sat in MissingSats:
            miss_line = SP3epoc.iloc[0].copy()
            miss_line["sat"]   = miss_sat
            miss_line["const"] = miss_sat[0]
            miss_line["x"]     = 0.000000
            miss_line["y"]     = 0.000000
            miss_line["z"]     = 0.000000
            miss_line["clk"]   = 999999.999999
            
            SP3epoc = SP3epoc.append(miss_line)

        SP3epoc.sort_values("sat",inplace=True,ascending=False)
        timestamp = conv.dt2sp3_timestamp(conv.numpy_datetime2dt(epoc)) + "\n"

        linefmt = "P{:}{:14.6f}{:14.6f}{:14.6f}{:14.6f}\n"

        LinesStkEpoch = []
        sum_val_epoch = 0
        for ilin , lin in SP3epoc.iterrows():
            line_out = linefmt.format(lin["sat"],lin["x"],lin["y"],lin["z"],lin["clk"])
            
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
    start_dt = conv.numpy_datetime2dt(np.min(EpochUsedList))
    
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
        print("ERR: not implemented yet !!!!!")
        raise Exception
        
    F = open(outpath_opera,"w+")
    F.write(FinalStr)
    
    
def write_clk(DFclk_in,clk_file_out,header="",output_std_values=False):
    HEAD = header
    Row_str_stk = []

    if output_std_values:
        row_str_proto = "{:2} {:4} {:4d} {:02d} {:02d} {:02d} {:02d} {:9.6f} {:2d}   {:19.12e} {:19.12e}"
    else:
        row_str_proto = "{:2} {:4} {:4d} {:02d} {:02d} {:02d} {:02d} {:9.6f} {:2d}   {:19.12e}"
        
    for irow, row in DFclk_in.iterrows():

        if output_std_values:
            one_or_two=2
            row_str = row_str_proto.format(row["type"],row["name"],row["year"],
                                           row["month"],row["day"],row["hour"],row["minute"],
                                           row["second"],one_or_two,row["bias"],row["sigma"])
        else:
            one_or_two=1
            row_str = row_str_proto.format(row["type"],row["name"],row["year"],
                                           row["month"],row["day"],row["hour"],row["minute"],
                                           row["second"],one_or_two,row["bias"])
            
        Row_str_stk.append(row_str)
    
    ## Add EOF
    Row_str_stk.append("EOF")
    
    CORPSE = "\n".join(Row_str_stk)
       
    OUT = HEAD + CORPSE

    with open(clk_file_out,"w+") as Fout:
        Fout.write(OUT)
        Fout.close()
        
    return OUT


def ine_block_mono(sat,dt_in,extra_intrvl_strt=.1,extra_intrvl_end=.4,step=300):
    
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
                        common_sats_only=True,
                        eliminate_null_sat=True,
                        severe=False,
                        separated_systems_export=False):
    """
    Parameters
    ----------
    ac_list : list
        3-character codes of the ACs.
    dir_in : str
        where the input sp3 are.
    dir_out : str
         where the output sp3 will be outputed.
    suffix_out_input : str, optional
        last char of the 3-char. code. The default is None, then it is s.
    overlap_size : int, optional
        Overlapsize. The default is 7200.
    force : True, optional
        force overwrite. The default is False.
    common_sats_only : True, optional
        generate a file with only the common sat between the 3 days.
        The default is True.
    eliminate_null_sat : bool, optional
        eliminate null sat. The default is True.
    severe : bool, optional
        raise an exception if problem. The default is False.
    separated_systems_export : bool, optional
        export different sp3 for different system. The default is False.

    Returns
    -------
    None.

    """
    


    for ac in ac_list:
        Lfile = utils.find_recursive(dir_in,"*" + ac + "*sp3")
        
        if not suffix_out_input:
            suffix_out = ac[:2] + "s" 
        else:
            suffix_out = suffix_out_input
        
        D     = []
        WWWWD = []
            
        for sp3 in Lfile:
            wwwwd_str = os.path.basename(sp3)[3:8]
            D.append(conv.gpstime2dt(int(wwwwd_str[:4]),int(wwwwd_str[4:])))
            
            
        for dat in D[1:-1]: ####if selection manuel, zip > 2lists !!!
            try:
                print("******",ac,dat)
                
                if conv.dt2gpstime(dat)[0] < 1800:
                    print("SKIP",dat)
                    continue
                    
                wwwwd_str = conv.dt_2_sp3_datestr(dat)
            
                dat_bef = dat - dt.timedelta(days=1)
                dat_aft = dat + dt.timedelta(days=1)
                
                wwwwd_str_bef = utils.join_improved("",*conv.dt2gpstime(dat_bef))
                wwwwd_str_aft = utils.join_improved("",*conv.dt2gpstime(dat_aft))
                
                
                ###### check if exsists
                dir_out_wk = os.path.join(dir_out,"wk" + str(wwwwd_str)[:4])
                utils.create_dir(dir_out_wk)
                fil_out = dir_out_wk + "/" + suffix_out  + wwwwd_str + ".sp3"
                
                if not force and os.path.isfile(fil_out):
                    print("0))",fil_out,"exsists, skipping...")
                    continue


                ### *************** STEP 1 ***************
                print("1)) Search for the days before/after")                
                print("1))",dat_bef,dat_aft)
                
                p1    = utils.find_regex_in_list(wwwwd_str     + ".sp3",Lfile,True)
                p_bef = utils.find_regex_in_list(wwwwd_str_bef + ".sp3",Lfile,True)
                p_aft = utils.find_regex_in_list(wwwwd_str_aft + ".sp3",Lfile,True)

                print("1)) Files found for the days before/after")                            
                print("0b)",p_bef)
                print("01)",p1)
                print("0a)",p_aft)
            
                if not p1 or not p_bef or not p_aft:
                    print("ERROR with day",dat)
                    continue
                
                SP3     = read_sp3(p1)
                SP3_bef = read_sp3(p_bef)
                SP3_aft = read_sp3(p_aft)
                
                SP3_bef = SP3_bef[SP3_bef["epoch"] < SP3["epoch"].min()]
                SP3_aft = SP3_aft[SP3_aft["epoch"] > SP3["epoch"].max()]
                
                SP3concat = pd.concat((SP3_bef,SP3,SP3_aft))
                
                dat_filter_bef = dat - dt.timedelta(seconds=overlap_size)
                dat_filter_aft = dat + dt.timedelta(seconds=overlap_size) + dt.timedelta(days=1)

                ### *************** STEP 2 ***************
                print("2)) dates of the overlap period before/after")                   
                print("2))",dat_filter_bef,dat_filter_aft)

                ### *************** STEP 3 *************** 
                print("3)) Dates of: SP3 concatenated, before, current, after")                       
                print("3))",SP3concat["epoch"].min(),SP3concat["epoch"].max())
                print("3b)",SP3_bef["epoch"].min(),SP3_bef["epoch"].max())
                print("31)",SP3["epoch"].min(),SP3["epoch"].max())
                print("3a)",SP3_aft["epoch"].min(),SP3_aft["epoch"].max())
                
                SP3concat = SP3concat[(SP3concat["epoch"] >= dat_filter_bef) & (SP3concat["epoch"] <= dat_filter_aft)]
                
                if common_sats_only:           
                    common_sats = set(SP3_bef["sat"]).intersection(set(SP3["sat"])).intersection(set(SP3_aft["sat"]))
                    SP3concat = SP3concat[SP3concat["sat"].isin(common_sats)]
                    
                    
                if eliminate_null_sat:
                    GoodSats = []
                    for sat in SP3concat["sat"].unique():
                        XYZvals = SP3concat[SP3concat["sat"] == sat][["x","y","z"]].sum(axis=1)
                        
                        V = np.sum(np.isclose(XYZvals,0)) / len(XYZvals)
                                            
                        if V < 0.50:
                            GoodSats.append(sat)
                        else:
                            print("6) eliminate because null position",sat)
                        
                    SP3concat = SP3concat[SP3concat["sat"].isin(GoodSats)]


                ### *************** STEP 7 ***************           
                print("7))","Start/End Epoch of the concatenated file ")                                     
                print("7))",SP3concat["epoch"].min(),SP3concat["epoch"].max())

                #### All systems        
                dir_out_wk = os.path.join(dir_out,"wk" + str(wwwwd_str)[:4])
                utils.create_dir(dir_out_wk)
                fil_out = dir_out_wk + "/" + suffix_out  + wwwwd_str + ".sp3"
                print("8)) outputed file")
                print(fil_out)
                write_sp3(SP3concat,fil_out)
                
                #### system separated
                if False:
                    for sys in SP3concat["const"].unique():
                        try:
                            SP3concat_sys = SP3concat[SP3concat["const"] == sys]
                            fil_out_sys = dir_out_wk + "/" + suffix_out[:2] + sys.lower() + wwwwd_str + ".sp3"
                            print("9)) outputed file")
                            print(fil_out_sys)
                            write_sp3(SP3concat_sys,fil_out_sys)
                        except:
                            continue
            
            except KeyboardInterrupt:
                raise KeyboardInterrupt
                
            except Exception as e:
                if severe:
                    raise e
                else:
                    print("ERR:",e)
                    raise e

