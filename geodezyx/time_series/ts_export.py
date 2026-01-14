#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 15:56:33 2023

@author: psakicki
"""

########## BEGIN IMPORT ##########
#### External modules
import copy
import datetime as dt
import numpy as np
#### Import the logger
import logging
import os
import matplotlib.pyplot as plt


#### geodeZYX modules
from geodezyx import conv
from geodezyx import time_series
from geodezyx import utils

log = logging.getLogger('geodezyx')

##########  END IMPORT  ##########

def export_ts_figure_pdf(fig,export_path,filename,close=False):
    """ fig can accept a int (id of a Figure)
         OR the figure Object itself """

    if type(fig) is int:
        f = plt.figure(fig)
    elif type(fig) is plt.Figure:
        f = fig

    f.set_size_inches(16.53,11.69)
    # tralala pour avoir des dates propres
    # parce que dans des cas + simples f.autofmt_xdate() suffit
    for a in f.axes[1:]:
        labels = a.get_xticklabels()
        for l in labels:
            l.set_rotation(40)
            l.set_horizontalalignment('right')
    out_path = os.path.join(export_path,filename+'.pdf')
    f.tight_layout()
    f.subplots_adjust(top=0.94)
    f.savefig(out_path,papertype='a4',format='pdf')

    if close:
        plt.close(f)

    return None

def export_ts_plot(tsin,export_path,coortype='ENU',export_type=("pdf","png"),
                   plot_B = False,close_fig_after_export=True):
    """ Very beta ...
        to be implemented : merge w/ the export_figure_pdf fct """
    # plot A avec les barres de sigma
    plt.clf()
    tsin.plot(coortype,fig=1)
    tsin.plot_discont()
    f = plt.gcf()
    f.set_size_inches(16.53,11.69)
    # tralala pour avoir des dates propres
    # parce que dans des cas + simples f.autofmt_xdate() suffit
    for a in f.axes[1:]:
        labels = a.get_xticklabels()
        for l in labels:
            l.set_rotation(40)
            l.set_horizontalalignment('right')
    f.tight_layout()
    f.subplots_adjust(top=0.92)
    for typ in export_type:
        export_file = os.path.join(export_path,tsin.stat+'a.' + typ)
        f.savefig(export_file,
                  papertype='a4',format=typ)
        log.info("plot exported in %s" , export_file)
    if close_fig_after_export:
        plt.close(f)

    # plot B avec une droite de regression
    if plot_B:
        time_series.linear_regress_ts(tsin)
        f = plt.gcf()
        #    f.set_dpi(300)
        f.set_size_inches(16.53,11.69)
        for a in f.axes:
            labels = a.get_xticklabels()
            for l in labels:
                l.set_rotation(40)
                l.set_horizontalalignment('right')
        f.tight_layout()
        f.subplots_adjust(top=0.92)
        for typ in export_type:
            export_file = os.path.join(export_path,tsin.stat+'b.' + typ)
            f.savefig(export_file,
                      papertype='a4',format=typ)
        log.info("plot exported in %s" , export_file)
        if close_fig_after_export:
            plt.close(f)
    return None


def export_ts(ts, outdir, coordtype='ENU', outprefix='', write_header=False):
    """
    Export a time series to a custom `.ts.dat` format.

    Parameters
    ----------
    ts : TimeSeries
        Input time series object containing geodetic data.
    outdir : str
        Directory where the output file will be saved.
    coordtype : str, optional
        Coordinate type to use for the export (e.g., 'ENU'), by default 'ENU'.
    outprefix : str, optional
        Prefix to add to the output file name, by default an empty string.
    write_header : bool, optional
        If True, writes a header to the output file, by default False.

    Notes
    -----
    - The function writes a `.ts.dat` file containing the time series data.
    - The file includes geodetic coordinates (East, North, Up or other types) and their uncertainties.
    - The header, if enabled, describes the columns in the file.

    Returns
    -------
    None
        The function writes the output file to the specified directory and does not return any value.
        write_header not well implemented !!!
    """

    proto_str = '{:23} ' * 14

    A,B,C,T,sA,sB,sC = ts.to_list(coordtype)
    if outprefix != '':
        outprefix = outprefix + '_'

    outfilenam = outprefix +  ts.stat + "_" +  ts.name + '.' + coordtype + '.ts.dat'
    outpath = os.path.join(outdir,outfilenam)

    filobj = open(outpath,'w+')

    if write_header:
        header = "#" + proto_str.format(*(coordtype[0],"sigma"+coordtype[0],
                                          coordtype[1],"sigma"+coordtype[1],
                                          coordtype[2],"sigma"+coordtype[2],
                                          "year","month","day",
                                          "hour","minute","seconds",
                                          "year_decimal","posix_time"))
        filobj.write(header + "\n")

    for a,b,c,t,sa,sb,sc in zip(A,B,C,T,sA,sB,sC):
        paramstr = [str(e) for e in [a,sa,b,sb,c,sc]]
        #paramstr = [e.ljust(18, '0') for e in paramstr]
        
        tt = conv.posix2dt(t)
        yr_dec_str = str(conv.dt2year_decimal(tt))
        posix_str = str(t)

        paramstr_time = list(tt.strftime("%Y %m %d %H %M %S").split())

        paramstr2 = paramstr + paramstr_time + [yr_dec_str,posix_str]
        
        outlin =  proto_str.format(*paramstr2) + "\n"
        filobj.write(outlin)
    filobj.close()

    log.info('timeserie exported in %s', outpath)
    return None


def export_ts_as_neu(tsin, outdir, outprefix, coordtype='ENU'):
    """
    Export a time series to a HECTOR .neu compatible format.

    Parameters
    ----------
    tsin : TimeSeries
        Input time series object containing geodetic data.
    outdir : str
        Directory where the output .neu file will be saved.
    outprefix : str
        Prefix to add to the output file name.
    coordtype : str, optional
        Coordinate type to use for the export, by default 'ENU'.

    Notes
    -----
    - The function writes a .neu file containing the time series data in HECTOR format.
    - The header includes metadata such as the reference epoch, reference position, and field descriptions.
    - The data section includes geodetic displacements (North, East, Up) relative to the first point.
    - If the time series lacks XYZ attributes, default values are used for the reference position.

    Returns
    -------
    None
        The function writes the output file to the specified directory and does not return any value.
    """
    if not hasattr(tsin[0],'X'):
        log.warning('no XYZ in ts')
        noXYZ = True
    else:
        noXYZ = False

    tswork = copy.deepcopy(tsin)
    #if coordtype == 'XYZ':
    #    mp = tswork.mean_posi()
    #    tswork.ENUcalc(mp)
    outpath = outdir +'/' + outprefix + tswork.stat + '.neu'
    outfile = open(outpath,'w+')
    E,N,U,T,sE,sN,sU = tswork.to_list(coordtype)
    first_pt = tswork.pts[0]
    if noXYZ:
        first_pt.X = 0.
        first_pt.Y = 0.
        first_pt.Z = 0.
        first_pt.L = 0.
        first_pt.H = 0.
        first_pt.F = 0.


    e0,n0,u0,t0 = list(zip(E,N,U,T))[0]
    # write the header
    outfile.write('# Site : {} \n'.format(tswork.stat))
    if 'calc_center' in list(tswork.anex.keys()):
        outfile.write('# Analysis Centre: {} \n'.format(tswork.anex['calc_center']))
    else:
        outfile.write('# Analysis Centre: N/A \n')

    outfile.write('# Solution code: GINS_PS \n')
    outfile.write('# Datum: ITRF2008\n')
    outfile.write('#\n')
    outfile.write('# Reference epoch: {}\n'.format(conv.dt2year_decimal(first_pt.Tdt)))
    outfile.write('# X : {}\n'.format(first_pt.X))
    outfile.write('# Y : {}\n'.format(first_pt.Y))
    outfile.write('# Z : {}\n'.format(first_pt.Z))
    outfile.write('#\n')
    outfile.write('# Longitude : {}\n'.format(first_pt.L))
    outfile.write('# Latitude  : {}\n'.format(first_pt.F))
    outfile.write('# Height    : {}\n'.format(first_pt.H))
    outfile.write('#\n')
    if coordtype == "ENU":
        outfile.write('# Components : ' + "NEU" + "\n")
    elif coordtype == "XYZ":
        outfile.write('# Components : ' + "YXZ" + "\n")
        outfile.write('# Cartesian components are undirect to maintain consistency with NEU\n')                
    outfile.write('#\n')
    if tswork.bool_discont:
        outfile.write('# type_of_offset : from discontinuties got from a station.info\n')
        outfile.write('#\n')
        for disc in sorted(tswork.discont):
            outfile.write('# offset {} 7\n'.format(conv.dt2year_decimal(disc)))
        outfile.write('#\n')
    # write the data
    for e,n,u,t,se,sn,su in zip(E,N,U,T,sE,sN,sU):
        t = conv.dt2year_decimal(conv.posix2dt(t))
        if coordtype == "ENU":        
            outfile.write('{:.5f}   {:+.6f}    {:+.6f}    {:+.6f} {:+.6f} {:+.6f} {:+.6f}\n'.format(t,n-n0,e-e0,u-u0,se,sn,su))
        elif coordtype == "XYZ":
            outfile.write('{:.5f}   {:+.6f}    {:+.6f}    {:+.6f} {:+.6f} {:+.6f} {:+.6f}\n'.format(t,n-n0,e-e0,u-u0,se,sn,su))

    log.info('timeserie exported in %s',outpath)
    return None


def export_ts_as_hector_enu(tsin,outdir,outprefix,coordtype = 'ENU'):
    """
    export to a HECTOR .enu (and not .neu !) compatible format
    This format is simpler : just gives MJD E N U

    This format is necessary to force a sampling period.

    outfile will be writed in
    /outdir/outprefixSTAT.enu
    """

    log.warning("NOT IMPLEMENTED YET !")

    return None

def export_ts_as_midas_tenu(tsin, outdir, outprefix, coordtype='ENU', export_step=True):
    """
    Export a time series to a MIDAS .tneu compatible format.

    Parameters
    ----------
    tsin : TimeSeries
        Input time series object containing geodetic data.
    outdir : str
        Directory where the output .tneu file will be saved.
    outprefix : str
        Prefix to add to the output file name.
    coordtype : str, optional
        Coordinate type to use for the export, by default 'ENU'.
    export_step : bool, optional
        If True, exports a step file containing discontinuities, by default True.

    Notes
    -----
    - The function writes a .tneu file containing the time series data in MIDAS format.
    - If `export_step` is True and the time series contains discontinuities, a separate step file is created.
    - The .tneu file includes geodetic displacements (East, North, Up) relative to the first point.

    Returns
    -------
    None
        The function writes the output file(s) to the specified directory and does not return any value.
    """
    if not hasattr(tsin[0],'X'):
        log.warning('no XYZ in ts')
        noXYZ = True
    else:
        noXYZ = False

    tswork = copy.deepcopy(tsin)
    stat = tswork.stat
    if coordtype == 'XYZ':
        mp = tswork.mean_posi()
        tswork.ENUcalc(mp)
    outpath = outdir +'/' + outprefix + tswork.stat + '.tenu'
    outfile = open(outpath,'w+')
    E,N,U,T,sE,sN,sU = tswork.to_list('ENU')
    first_pt = tswork.pts[0]
    if noXYZ:
        first_pt.X = 0.
        first_pt.Y = 0.
        first_pt.Z = 0.
        first_pt.L = 0.
        first_pt.H = 0.
        first_pt.F = 0.

    e0,n0,u0,t0 = list(zip(E,N,U,T))[0]

    for e,n,u,t in zip(E,N,U,T):
        t = conv.dt2year_decimal(conv.posix2dt(t))
        #outfile.write('{} {:.5f} {:+.6f} {:+.6f} {:+.6f} \n'.format(stat,t,n-n0,e-e0,u-u0))
        outfile.write('{} {:.5f} {:+.6f} {:+.6f} {:+.6f} \n'.format(stat,t,e-e0,n-n0,u-u0))

    log.info('timeserie exported in %s',outpath)

    if export_step and tswork.bool_discont:
        outpath_step = outdir +'/' + outprefix + tswork.stat + '.step'
        outfile_step = open(outpath_step,'w+')
        for d in tswork.discont:
            d = conv.dt2year_decimal(d)
            line = tswork.stat + " " + str(d) + "\n"
            outfile_step.write(line)

            log.info('timeserie discont. (steps) exported in %s',outpath_step)

    return None


def export_ts_as_pbo_pos(tsin, outdir, outprefix='' ,
                         force=None, force_1st_pt_as_ref=True,
                         verbose=False):
    """
    Export a time series to the GAMIT/GLOBK PBO position (.pos) format.

    Parameters
    ----------
    tsin : TimeSeries
        Input time series object containing geodetic data.
    outdir : str
        Directory where the output .pos file will be saved.
    outprefix : str, optional
        Prefix to add to the output file name, by default an empty string.
    force : str, optional
        Determines the data source to use for the export. Can be 'data' or 'data_xyz'.
        If None, the function decides based on the available attributes, by default None.
    force_1st_pt_as_ref : bool, optional
        If True, forces the first point of the time series to be used as the reference position, by default True.
    verbose : bool, optional
        If True, enables verbose logging for debugging purposes, by default False.

    Notes
    -----
    - The function writes a .pos file containing the time series data in the PBO format.
    - The header includes metadata such as the reference epoch, reference position, and field descriptions.
    - The data section includes geodetic coordinates (X, Y, Z) and their uncertainties, along with other metadata.

    Returns
    -------
    outfile_path : str
        Path to the exported PBO file
    """

    # force option
    
    # if force is not None:
    #     if force not in ['data','data_xyz']:
    #         print('!!! ERROR: force keyword must be either data or data_xyz')
    #         return(self)
    
    #     if force == 'data':
    #         self.data_xyz = None

    #     if force == 'data_xyz':
    #         self.data = None

    # # cdata
    
    # if force != 'data_xyz':
    #     if not self.cdata(data=True):
    #         print('! Can not write .pos file. Problem with .data attribute.')
    #         self.cdata(data=True,verbose=True)
    #         return(self)


    ###############################################
    def __print_header_pos_file(f_pos,code,
                                first_epoch,last_epoch,
                                release_date,XYZ_ref,
                                geo_ref,verbose=False):
        
        f_pos.write("PBO Station Position Time Series. Reference Frame : Unknown\n")
        f_pos.write("Format Version: 1.1.0\n")
        f_pos.write("4-character ID: %s\n" % code)
        f_pos.write("Station name  : %s        \n" % code)
        f_pos.write("First Epoch   : %s\n" % first_epoch)
        f_pos.write("Last Epoch    : %s\n" % last_epoch)
        f_pos.write("Release Date  : %s\n" % release_date)
        f_pos.write("XYZ Reference position :   %s\n" % XYZ_ref)
        f_pos.write("NEU Reference position :   %s (Unknown/WGS84)\n" % geo_ref)
        f_pos.write("Start Field Description\n")
        f_pos.write("YYYYMMDD      Year, month, day for the given position epoch\n")
        f_pos.write("HHMMSS        Hour, minute, second for the given position epoch\n")
        f_pos.write("JJJJJ.JJJJJ   Modified Julian day for the given position epoch\n")
        f_pos.write("X             X coordinate, Specified Reference Frame, meters\n")
        f_pos.write("Y             Y coordinate, Specified Reference Frame, meters\n")
        f_pos.write("Z             Z coordinate, Specified Reference Frame, meters\n")
        f_pos.write("Sx            Standard deviation of the X position, meters\n")
        f_pos.write("Sy            Standard deviation of the Y position, meters\n")
        f_pos.write("Sz            Standard deviation of the Z position, meters\n")
        f_pos.write("Rxy           Correlation of the X and Y position\n")
        f_pos.write("Rxz           Correlation of the X and Z position\n")
        f_pos.write("Ryz           Correlation of the Y and Z position\n")
        f_pos.write("Nlat          North latitude, WGS-84 ellipsoid, decimal degrees\n")
        f_pos.write("Elong         East longitude, WGS-84 ellipsoid, decimal degrees\n")
        f_pos.write("Height (Up)   Height relative to WGS-84 ellipsoid, m\n")
        f_pos.write("dN            Difference in North component from NEU reference position, meters\n")
        f_pos.write("dE            Difference in East component from NEU reference position, meters\n")
        f_pos.write("du            Difference in vertical component from NEU reference position, meters\n")
        f_pos.write("Sn            Standard deviation of dN, meters\n")
        f_pos.write("Se            Standard deviation of dE, meters\n")
        f_pos.write("Su            Standard deviation of dU, meters\n")
        f_pos.write("Rne           Correlation of dN and dE\n")
        f_pos.write("Rnu           Correlation of dN and dU\n")
        f_pos.write("Reu           Correlation of dEand dU\n")
        f_pos.write("Soln          corresponding to products  generated with rapid or final orbit products, in supplemental processing, campaign data processing or reprocessing\n")
        f_pos.write("End Field Description\n")
        f_pos.write("*YYYYMMDD HHMMSS JJJJJ.JJJJ         X             Y             Z            Sx        Sy       Sz     Rxy   Rxz    Ryz            NLat         Elong         Height         dN        dE        dU         Sn       Se       Su      Rne    Rnu    Reu  Soln\n")

    ###############################################    

    def __format_data(tsinin):
                
        lformatted_out = []
        
        for pt in tsinin.pts:
            YYYYMMDD = int(conv.dt2str(pt.Tdt,"%Y%m%d"))
            HHMMSS = int(conv.dt2str(pt.Tdt,"%H%M%S"))
            MJD = conv.dt2mjd(pt.Tdt)
            X, Y, Z = pt.X, pt.Y, pt.Z
            Sx, Sy, Sz = pt.sX, pt.sY, pt.sZ  
            Rxy,Rxz,Ryz = pt.anex.get('sdXY', 0.0),pt.anex.get('sdXZ', 0.0),pt.anex.get('sdYZ', 0.0)
            Elong,NLat,Height = pt.F,pt.L,pt.H
            dN,dE,dU = pt.N,pt.E,pt.U
            Sn,Se,Su = pt.sN,pt.sE,pt.sU
            Rne,Rnu,Reu = pt.anex.get('sdEN', 0.0),pt.anex.get('sdNU', 0.0),pt.anex.get('sdEU', 0.0)
            Soln = "XXXXXXXXXXXXXX"
            
            l = [YYYYMMDD,HHMMSS,MJD,X,Y,Z,Sx,Sy,Sz,Rxy,Rxz,Ryz,Elong,NLat,
                 Height,dN,dE,dU,Sn,Se,Su,Rne,Rnu,Reu,Soln]
            
            lformatted_out.append(l)
            
        return lformatted_out
            
    
    ###############################################
    def __print_content_pos_file(f_pos,lformatted):

        for l in lformatted:
                        
            #l = [YYYYMMDD,HHMMSS,MJD,X,Y,Z,Sx,Sy,Sz,Rxy,Rxz,Ryz,Elong,NLat,Height,dN,dE,dU,Sn,Se,Su,Rne,Rnu,Reu,Soln]
            line_proto = " {:8d} {:6d} {:10.4f} {:14.5f} {:14.5f} {:14.5f} {:8.5f} {:8.5f} {:8.5f} {:6.3f} {:6.3f} {:6.3f}     {:14.10f}  {:14.10f} {:10.5f}  {:10.5f}{:10.5f}{:10.5f}   {:8.5f} {:8.5f} {:8.5f} {:6.3f} {:6.3f} {:6.3f} {:s}\n"
            
            f_pos.write(line_proto.format(*l))

    ###############################################


    tswork = copy.deepcopy(tsin)
    
    if not tswork.boolENU:
        tswork.ENUcalc_from_first_posi()
    else:
        first_pt = tswork.pts[0]
        if first_pt.E + first_pt.N + first_pt.U != 0. and force_1st_pt_as_ref:
            tswork.ENUcalc_from_first_posi()


    # output directory
    if not os.path.isdir(outdir):
        try:
            os.mkdir(outdir)
            print("-- Creating directory ",outdir)
        except:
            print("!!! Error : Could not create directory ",outdir)
    

        
    outpath = outdir +'/' +  tswork.stat + outprefix + '.pos'
        
    # open pos file
    f_pos=open(outpath,'w+')
    
    import datetime
    
    current_date=datetime.datetime.now()
    release_date= conv.dt2str(current_date) #("%d%02d%02d %02d%02d%02d"% (current_date.year,current_date.month,current_date.day,current_date.hour,current_date.minute,current_date.second))
    
    refpt = tswork.refENU
    
    XYZ_ref=(" %14.5lf %14.5lf %14.5lf (unknown)" % (refpt.X,refpt.Y,refpt.Z))
    geo_ref=(" %16.10lf %16.10lf %16.10lf (unknown)" % (refpt.F,refpt.L,refpt.H))
    
    __print_header_pos_file(f_pos,
                            tswork.stat,
                            tswork.startdate(),
                            tswork.enddate(),
                            release_date,
                            XYZ_ref,
                            geo_ref)
    
    
    wpos = __format_data(tswork)

    __print_content_pos_file(f_pos,wpos)
    
    f_pos.close()
    
    return outpath


def export_ts_as_spotgins(tsin, outdir, ac, data_src="unknown", version=3):
    """
    Export time series to SPOTGINS format.

    Parameters
    ----------
    tsin : TimeSeries
        Input time series object
    outdir : str
        Output directory path
    ac : str
        Analysis Center acronym (3-4 characters)
    data_src : str, optional
        Data source identifier, by default "unknown"
    version : int, optional
        SPOTGINS format version (2 or 3), by default 2

    Returns
    -------
    outfile_path : str
        Path to the exported SPOTGINS file
    """
    now_date = utils.get_timestamp(utc=True)
    df = tsin.to_dataframe("ENU")
    df["Tdt"] = df["Tdt"] + dt.timedelta(seconds=19)
    name = tsin.stat
    outfile_path = os.path.join(outdir, f"SPOTGINS_{name}.enu")

    datexelis = [p.anex.get('dateofexe', 'NA') for p in tsin.pts]
    ginsverslis = [p.anex.get('gins_version', 'NA') for p in tsin.pts]

    const_lbda = lambda x: "GE" if x >= dt.datetime(2018, 10, 7) else "G"
    const = const_lbda(np.max(df["Tdt"]))

    if version == 2:
        fmtstr2 = " {:14.8f} {:14.6f} {:14.6f} {:14.6f} {:14.6f} {:14.6f} {:14.6f}  {:4s}{:2s}{:2s}{:2s}{:2s}{:2.0f}  {:12.7f}  {:5s}  {:s}   {:s}\n"
        head2 = f"""# SPOTGINS SOLUTION [POSITION] v2
# DATE             : {now_date}
#------------------------------------
# STATION          : {name}
# OPERATOR         : {ac}
# STRATEGY_SUMMARY : https://www.poleterresolide.fr/visualisation-gnss/#/solution/SPOTGINS
# REF_FRAME        : IGS20
# PRODUCTS         : G20/GRG
# CONSTEL          : {const}
# UNITS            : meters
# ELLIPSOID        : GRS80
# DATA_SOURCE      : {data_src}
# ACKNOWLEDGMENTS  :
#------------------------------------
# X_pos            : {tsin.refENU.X:.6f}
# Y_pos            : {tsin.refENU.Y:.6f}
# Z_pos            : {tsin.refENU.Z:.6f}
# Longitude        : {tsin.refENU.L:.6f}
# Latitude         : {tsin.refENU.F:.6f}
# Height           : {tsin.refENU.H:.6f}
#------------------------------------
#jjjjj.jjjjjjjj         _____E         _____N         _____U         ____dE         ____dN         ____dU  yyyymmddHHMMSS  yyyy.yyyyyyy  Const  Dateofexe       GinsVersion
"""

        with open(outfile_path, "w") as outfile:
            outfile.write(head2)

            for ir, r in df.iterrows():
                t = r["Tdt"]
                outstr = fmtstr2.format(
                    np.round(conv.dt2mjd(t), 7),
                    r["E"], r["N"], r["U"],
                    r["sE"], r["sN"], r["sU"],
                    str(t.year),
                    str(t.month).zfill(2),
                    str(t.day).zfill(2),
                    str(t.hour).zfill(2),
                    str(t.minute).zfill(2),
                    t.second,
                    conv.dt2year_decimal(t),
                    const_lbda(t),
                    datexelis[ir],
                    ginsverslis[ir]
                )
                outfile.write(outstr)

    elif version == 3:
        # Extract correlation coefficients from point data
        prairieverslis = [p.anex.get('prairie_version', 'NA') for p in tsin.pts]
        qflaglis = [p.anex.get('quality_flag', 'NA') for p in tsin.pts]

        # Calculate correlation coefficients (placeholder - adjust based on actual data structure)
        corr_en = [p.anex.get('sdEN', 0.0) for p in tsin.pts]
        corr_eu = [p.anex.get('sdEU', 0.0) for p in tsin.pts]
        corr_nu = [p.anex.get('sdNU', 0.0) for p in tsin.pts]

        fmtstr3 = " {:6.1f} {:14.6f} {:14.6f} {:14.6f} {:14.6f} {:14.6f} {:14.6f} {:10.6f} {:10.6f} {:10.6f}  {:4s}-{:2s}-{:2s}T12:00:00  {:11.6f}  {:5s}  {:4s}  {:13s}  {:11s}  {:3s}\n"

        head3 = f"""# SPOTGINS SOLUTION [POSITION] v3
# Creation {now_date}
#----------------------------------------------------------
# STATION          : {name}
# ANALYSIS_CENTRE  : {ac}
# STRATEGY_SUMMARY : https://www.poleterresolide.fr/geodesy-plotter/#/solution/SPOTGINS
# REF_FRAME        : IGS20
# PRODUCTS         : G20/GRG
# CONSTELLATION    : {const}
# UNITS            : meters
# ELLIPSOID        : GRS80
# DATA_SOURCE      : {data_src}
# FLAG             : 0 = Full quality
# FLAG             : 1 = Antenna/radome not calibrated
# FLAG             : 2 = Less than 12 h of observations
# FLAG             : 4 = Less than 80% of fixed ambiguities
# ACKNOWLEDGMENTS  :
#----------------------------------------------------------
# X_pos            : {tsin.refENU.X:15.6f}
# Y_pos            : {tsin.refENU.Y:15.6f}
# Z_pos            : {tsin.refENU.Z:15.6f}
# Longitude        : {tsin.refENU.L:11.6f}
# Latitude         : {tsin.refENU.F:11.6f}
# Height           : {tsin.refENU.H:11.6f}
#----------------------------------------------------------
#MJD           DispEast      DispNorth         DispUp      SigmaEast     SigmaNorth        SigmaUp     CorrEN     CorrEU     CorrNU  yyyy-mm-ddTHH:MM:SS  DecimalYear  Const  Flag  DateOfExe      GinsVersion  PrairieVersion
"""

        with open(outfile_path, "w") as outfile:
            outfile.write(head3)

            for ir, r in df.iterrows():
                t = r["Tdt"]
                mjd = np.round(conv.dt2mjd(t), 7) # Convert to MJD (Modified Julian Date)

                outstr = fmtstr3.format(
                    mjd,
                    r["E"], r["N"], r["U"],
                    r["sE"], r["sN"], r["sU"],
                    corr_en[ir], corr_eu[ir], corr_nu[ir],
                    str(t.year),
                    str(t.month).zfill(2),
                    str(t.day).zfill(2),
                    conv.dt2year_decimal(t),
                    const_lbda(t),
                    str(qflaglis[ir]),
                    datexelis[ir],
                    ginsverslis[ir],
                    str(prairieverslis[ir])
                )
                outfile.write(outstr)

    #log.info('SPOTGINS v%s time series exported to %s', version, outfile_path)
    outfile.close()
    return outfile_path

