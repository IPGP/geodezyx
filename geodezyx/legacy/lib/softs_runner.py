# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 17:59:36 2015

@author: psakicki

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics

Copyright (C) 2019 Pierre Sakic (GFZ, pierre.sakic@gfz-postdam.de)
GitHub repository :
https://github.com/PierreS1/GeodeZYX-Toolbox-Lite

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
"""

import geodetik as geok
import geoclass
import genefun
import geo_files_converter_lib as gfc

import os
import shutil
import subprocess
import glob
import matplotlib
#matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import numpy as np
import datetime
import time
import os
import urllib.request, urllib.parse, urllib.error
import datetime as dt
import urllib.request, urllib.error, urllib.parse
import multiprocessing as mp
import string
import re
import configparser
import collections
import dateutil
import string
import pandas


#  _    _ ______ _____ _______ ____  _____
# | |  | |  ____/ ____|__   __/ __ \|  __ \
# | |__| | |__ | |       | | | |  | | |__) |
# |  __  |  __|| |       | | | |  | |  _  /
# | |  | | |___| |____   | | | |__| | | \ \
# |_|  |_|______\_____|  |_|  \____/|_|  \_\
#

def keeping_specific_stats(listoffiles,specific_stats,invert=False):
    """ if invert = True : NOT keeping BUT removing specific stats"""
    if not (type(specific_stats) is tuple):
        raise Exception('keeping_specific_stats : specific_stats must be a tuple')
    if not invert:
        outlistoffiles = []
    else:
        outlistoffiles = list(listoffiles)

    if specific_stats != ():
        for fil in listoffiles:
            for stat in specific_stats:
                if stat in fil:
                    if not invert:
                        outlistoffiles.append(fil)
                    else:
                        outlistoffiles.remove(fil)
        return outlistoffiles
    else:
        return listoffiles

def neufile_outlier_removing(inp_neufile,generik_conf_file,outdir='',remove_ctl_file=True):
    """ from a NEU file
    => removeoutlier preprocessing
    => 3 MOM files (one per component) """
    neufile_name = os.path.basename(inp_neufile)
    prefix_inp = neufile_name.replace('.neu','_')
    inpdir = os.path.dirname(inp_neufile)
    if outdir == '':
        outdir = inpdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    comp_list = [ 'North' , 'East' , 'Up' ]

    for comp in comp_list:
        work_conf_file = inpdir + '/' + prefix_inp  + comp[0] + '_rmoutlier.ctl'
        prefix_out = prefix_inp + comp[0] + '_pre'
        momfile_name = prefix_out + '.mom'
        momfile_path = outdir + '/' + momfile_name
        shutil.copyfile(generik_conf_file, work_conf_file)
#        work_conf_file_obj = open(work_conf_file,'w+')
        genefun.replace(work_conf_file, 'DataFile',      'DataFile            ' + neufile_name)
        genefun.replace(work_conf_file, 'DataDirectory', 'DataDirectory       ' + inpdir)
        genefun.replace(work_conf_file, 'OutputFile',    'OutputFile          ' + momfile_path)
        genefun.replace(work_conf_file, 'component' ,    'component           ' + comp)

        p = subprocess.Popen('',executable='/bin/bash', stdin=subprocess.PIPE , stdout=subprocess.PIPE , stderr=subprocess.PIPE)
        command = "removeoutliers " +  work_conf_file
        print('LAUNCHING : ',  command)
        stdout,stderr = p.communicate( command.encode() )
        # logs files
        std_file = open(outdir + '/' + prefix_out + ".std.log", "w")
        std_file.write(stdout.decode("utf-8"))
        std_file.close()
        if stderr:
            print("err.log is not empty, must be checked !")
            print(stderr.decode("utf-8"))
            print('')
            err_file = open(outdir + '/' + prefix_out + ".err.log", "w")
            err_file.write(stderr.decode("utf-8"))
            err_file.close()

        # remove the ctl file
        if remove_ctl_file:
            os.remove(work_conf_file)

    return None


def multi_neufile_outlier_removing(inpdir,generik_conf_file,outdir='',
                                   extention='neu',specific_stats=(),
                                   invert_specific=False,remove_ctl_file=True):
    if not os.path.exists(inpdir):
        os.makedirs(inpdir)
    os.chdir(inpdir)

    wildcarded_path = inpdir + '/' + '*.' + extention

    listofenufile = glob.glob(wildcarded_path)

    if len(listofenufile) == 0:
        print("WARN : no files found in specified directory ...")
        print(wildcarded_path)

    listofenufile = keeping_specific_stats(listofenufile,specific_stats,
                                           invert=invert_specific)

    print("removing outliers of", len(listofenufile), ' timeseries')

    for enu in sorted(listofenufile):
        neufile_outlier_removing(enu,generik_conf_file,outdir=outdir,
                                 remove_ctl_file=remove_ctl_file)

    return None


def momfile_trend_processing(inp_momfile,generik_conf_file,outdir='',
                             remove_ctl_file=True):
    momfile_inp_name = os.path.basename(inp_momfile)
    inpdir = os.path.dirname(inp_momfile)
    if outdir == '':
        outdir = inpdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    work_conf_file = inpdir + '/' + momfile_inp_name.replace('.mom','') + '_trend.ctl'
    prefix_out_name = momfile_inp_name.replace('_pre.mom','')
    momfile_out_name = momfile_inp_name.replace('pre.mom','out.mom')
    momfile_out_path = outdir + '/' + momfile_out_name
    shutil.copyfile(generik_conf_file, work_conf_file)
    genefun.replace(work_conf_file, 'DataFile',      'DataFile            ' + momfile_inp_name)
    genefun.replace(work_conf_file, 'DataDirectory', 'DataDirectory       ' + inpdir)
    genefun.replace(work_conf_file, 'OutputFile',    'OutputFile          ' + momfile_out_path)

    p = subprocess.Popen('',executable='/bin/bash', stdin=subprocess.PIPE , stdout=subprocess.PIPE , stderr=subprocess.PIPE)
    command = "estimatetrend " +  work_conf_file
    print('LAUNCHING : ',  command)
    stdout,stderr = p.communicate( command.encode() )
    # logs files
    sumfile = outdir + '/' + prefix_out_name + ".sum"
    std_file = open(sumfile, "w")
    std_file.write(stdout.decode("utf-8") )
    std_file.close()
    if stderr:
        print("err.log is not empty, must be checked !")
        print(stderr.decode("utf-8"))
        print('')
        err_file = open(outdir + '/' + prefix_out_name + ".err.log", "w")
        err_file.write(stderr.decode("utf-8"))
        err_file.close()

    # plot a graph
    # matplotlib.use('Agg')
    fig = plt.figure()
    try:
        ax = fig.gca()
        ax.set_xlabel("Date")
        ax.set_ylabel("Displacement (m)")
        data = np.loadtxt(momfile_out_path)
        ax.plot(MJD2dt(data[:,0]),data[:,1],'b+')
        ax.plot(MJD2dt(data[:,0]),data[:,2],'g.')

        # add offsets
        mom_obj = open(inp_momfile,'r+')
        offset_stk = []
        for line in mom_obj:
            if 'offset' in line:
                offset_stk.append(float(line.split()[-1]))
        print('offset detected for plot :', offset_stk)
        for off in offset_stk:
            plt.axvline(MJD2dt(off),color='r')
        mom_obj.close()

        # add trend & bias
        mom_obj = open(sumfile,'r+')
        for line in mom_obj:
            try:
                if 'bias :' in line[0:7]:
                    biasline = line
                if 'trend:' in line[0:7]:
                    trendline = line
            except:
                continue
        mom_obj.close()

    #    .text(0.1, 0.15, biasline , horizontalalignment='left',verticalalignment='center',transform=ax.transAxes)
        plt.figtext(0.1, 0.05, trendline + biasline , horizontalalignment='left',verticalalignment='center',transform=ax.transAxes)
        plt.suptitle(prefix_out_name)

        #save
        fig.set_size_inches(16.53,11.69*.6667)
        fig.autofmt_xdate()
        fig.savefig(outdir + '/' + prefix_out_name + '.plt.pdf' )
        fig.savefig(outdir + '/' + prefix_out_name + '.plt.png' )

    except:
        plt.close(fig)

    # remove the ctl file
    if remove_ctl_file:
        os.remove(work_conf_file)

    return None

def MJD2dt(mjd_in):
    # cf http://en.wikipedia.org/wiki/Julian_day
    try:
        return [datetime.datetime(1858,11,17) + datetime.timedelta(m )for m in mjd_in]
    except:
        return datetime.datetime(1858,11,17) + datetime.timedelta(mjd_in)

def multi_momfile_trend_processing(inpdir,generik_conf_file,outdir='',extention='pre.mom',remove_ctl_file=True,specific_stats=(),invert_specific=False):
    start = time.time()
    if not os.path.exists(inpdir):
        os.makedirs(inpdir)
    os.chdir(inpdir)
    listofmomfile = glob.glob(inpdir + '/' + '*' + extention)
    listofmomfile = keeping_specific_stats(listofmomfile,specific_stats,invert=invert_specific)

    print("processing", len(listofmomfile), 'files')

    for imom , mom in enumerate(sorted(listofmomfile)):
        print("processing file ", imom+1 , "/" , len(listofmomfile))
        momfile_trend_processing(mom,generik_conf_file,outdir=outdir,remove_ctl_file=remove_ctl_file)

    print('multi_momfile_trend_processing exec time : ' , time.time() - start)
    return None

# ===== FCTS 4 multi_sumfiles_trend_extract =====

def sumfiles_to_statdico(inpdir,specific_stats=(),invert_specific=False):
    ''' this fct search for every sum file in a folder
        a stat dico contains no data
        only the paths to the E,N,U sum files
        statdico[stat] = [path/E.sum,path/N.sum,path/U.sum ]
         for each stat getting the 3 ENU sum files

        Thoses lists will be send in sumfiles_trend_extract
         '''
    if not os.path.exists(inpdir):
        os.makedirs(inpdir)
    os.chdir(inpdir)
    listofsumfile = glob.glob(inpdir + '/' + '*' + '.sum')
    listofsumfile = keeping_specific_stats(listofsumfile,specific_stats,invert=invert_specific)
    listofstat = [os.path.basename(sumfil).split('.')[0].split('_')[-2] for sumfil in listofsumfile]
    statdico = {}
    for stat in listofstat:
        statdico[stat] = []
        for sumfil in listofsumfile:
            if stat in sumfil:
                statdico[stat].append(sumfil)
    return statdico

def sumfiles_trend_extract(listof3ENUsumfile):
    V , sV = {} , {}
    if len(listof3ENUsumfile) != 3:
        raise Exception('listofENUsumfile != 3')
    for fil in listof3ENUsumfile:
        bnfil = os.path.basename(fil)
        coord = bnfil.split('.')[0][-1]
        for l in open(fil):
            if 'trend:' in l:
                f = l.split()
                V[coord]  = float(f[1])
                sV[coord] = float(f[3])

    statname = bnfil.split('.')[0].split('_')[-2]
    return statname , V , sV

def get_FLH_from_NEUfile(neufilepath):
    for l in open(neufilepath):
        f = l.split()
        if 'Longitude' in l:
            lon = float(f[-1])
        if 'Latitude' in l:
            lat = float(f[-1])
        if 'Height' in l:
            hau = float(f[-1])
    return lat , lon , hau

def velfile_from_a_list_of_statVsV_tuple(listoftup,out_dir, out_prefix , raw_neu_dir='',style='epc'):
    if not os.path.exists(out_dir ):
        os.makedirs(out_dir  )
    """
    make a GLOBK style .vel file
    """

    ### HEADER DEFINITION
    if style == 'globk':
        outfile = open(out_dir +'/' + out_prefix + '.vel','w+')
        outfile.write('* Velocity field created with Hector Runner / P. Sakic - La Rochelle Univ. (FRA) \n')
        outfile.write('*  Long         Lat        Evel    Nvel    dEv     dNv    E +-    N +-    Rne      Hvel     dHv    H +-  Site\n')
        outfile.write('*  deg          deg       mm/yr   mm/yr   mm/yr   mm/yr   mm/yr   mm/yr            mm/yr   mm/yr   mm/yr\n')
        k = 1000.
    elif style == 'epc':
        outfile = open(out_dir +'/' + out_prefix + '.vel.neu','w+')
        k = 1
    elif style == "csv_renag":
        outfile = open(out_dir +'/' + out_prefix + '.vel.csv','w+')
        outfile.write('Station,V_North,V_East,V_Up,sV_North,sV_East,sV_Up\n')
        k = 1000.
    elif style == "dataframe":
        column_names = ['Station','Latitude','Longitude','V_North','V_East','V_Up','sV_North','sV_East','sV_Up']
        lines_stk = []
        k = 1


    ### FILLING DEFINITION
    for tup in listoftup:
        stat , V , sV = tup
        try:
            neufile = genefun.regex2filelist(raw_neu_dir,stat)[0]
            lat , lon , hau = get_FLH_from_NEUfile(neufile)
        except:
            lat , lon , hau = 0,0,0
        print(stat)
        if style == 'globk':
            line = '{:11.5f}{:11.5f} {:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f} {:6.3f}  {:8.2f}{:8.2f}{:8.2f} {}\n'.format(lon , lat , V['E']*k , V['N']*k , 0. , 0. , sV['E']*k , sV['N']*k , 0. , V['U']*k ,0. ,  sV['U']*k , stat + '_GPS')
        elif style == 'epc':
            line = '{} {} {} {} {} {} {} {}\n'.format(stat,lat,geok.wrapTo180(lon),V['N']*k,V['E']*k,sV['N']*k,sV['E']*k,0)
        elif style == 'csv_renag':
            line = '{},{:10.5f},{:10.5f},{:10.5f},{:10.5f},{:10.5f},{:10.5f}\n'.format(stat,V['N']*k,V['E']*k,V['U']*k,sV['N']*k,sV['E']*k,sV['U']*k)
        elif style == "dataframe":
            line = [stat,lat,lon,V['N']*k,V['E']*k,V['U']*k,sV['N']*k,sV['E']*k,sV['U']*k]
            lines_stk.append(line)

        if style != "dataframe":
            outfile.write(line)

    #### FINALISATION
    if style == "dataframe":
        DF = pandas.DataFrame(lines_stk,columns=column_names)
        genefun.pickle_saver(DF,out_dir,out_prefix + "_DataFrame")
        return DF
    else:
        outfile.close()
        return None

#def velNEUfile_from_a_list_of_statVsV_tuple(listoftup,out_dir, out_prefix , raw_neu_dir=''):
#    """ make a dirty velocity file compatible with EPC """
#    outfile = open(out_dir +'/' + out_prefix + '.vel.neu','w+')
#    for tup in listoftup:
#        stat , V , sV = tup
#        try:
#            neufile = genefun.regex2filelist(raw_neu_dir,stat)[0]
#            lat , lon , hau = get_FLH_from_NEUfile(neufile)
#        except:
#            lat , lon , hau = 0,0,0
#
#        line = '{} {} {} {} {} {} {} {}\n'.format(stat,lat,geok.wrapTo180(lon),V['N'],V['E'],sV['N'],sV['E'],0)
#        outfile.write(line)
#    return None

def multi_sumfiles_trend_extract(inp_dir , out_dir , out_prefix ,
                                 raw_neu_dir='',specific_stats=(),
                                 invert_specific=False,style='epc'):
    ''' style = globk OR epc
    make a GLOBK style .vel file or
    make a dirty velocity file compatible with EPC '''

    genefun.create_dir(out_dir)

    statdico = sumfiles_to_statdico(inp_dir,specific_stats,invert_specific=invert_specific)
    listoftup = []
    for stat , listof3ENUsumfile in statdico.items():
        tup = sumfiles_trend_extract(listof3ENUsumfile)
        if tup[0] != stat:
            raise Exception('tup[0] != stat')
        listoftup.append(tup)
    listoftup.sort(key=lambda x : x[0])
    velfile_from_a_list_of_statVsV_tuple(listoftup,out_dir, out_prefix , raw_neu_dir,style)



#  __  __ _____ _____           _____
# |  \/  |_   _|  __ \   /\    / ____|
# | \  / | | | | |  | | /  \  | (___
# | |\/| | | | | |  | |/ /\ \  \___ \
# | |  | |_| |_| |__| / ____ \ ____) |
# |_|  |_|_____|_____/_/    \_\_____/
#











# _____  _____ _   _ ________   __  _____                      _                 _
#|  __ \|_   _| \ | |  ____\ \ / / |  __ \                    | |               | |
#| |__) | | | |  \| | |__   \ V /  | |  | | _____      ___ __ | | ___   __ _  __| | ___ _ __
#|  _  /  | | | . ` |  __|   > <   | |  | |/ _ \ \ /\ / / '_ \| |/ _ \ / _` |/ _` |/ _ \ '__|
#| | \ \ _| |_| |\  | |____ / . \  | |__| | (_) \ V  V /| | | | | (_) | (_| | (_| |  __/ |
#|_|  \_\_____|_| \_|______/_/ \_\ |_____/ \___/ \_/\_/ |_| |_|_|\___/ \__,_|\__,_|\___|_|

def igs_garner_server(stat,date):
    # plante si trop de requete
    urlserver = "ftp://garner.ucsd.edu/pub/rinex/"
    rnxname = geok.statname_dt2rinexname(stat.lower(),date)
    url = os.path.join( urlserver , str(date.year) , geok.dt2doy(date) , rnxname )
    return url

def igs_cddis_server(stat,date):
    # a privilegier
    urlserver = "ftp://cddis.gsfc.nasa.gov/gps/data/daily/"
    rnxname = geok.statname_dt2rinexname(stat.lower(),date)
    url = os.path.join(urlserver , str(date.year) , geok.dt2doy(date) , date.strftime('%y') + 'd' , rnxname)
    return url

def igs_cddis_nav_server(stat,date):
    # a privilegier
    urlserver = "ftp://cddis.gsfc.nasa.gov/gps/data/daily/"
    rnxname = geok.statname_dt2rinexname(stat.lower(),date,'n.Z')
    url = os.path.join(urlserver , str(date.year) , geok.dt2doy(date) , date.strftime('%y') + 'n' , rnxname)
    return url

def rgp_ign_smn_server(stat,date):
    urlserver = "ftp://rgpdata.ign.fr/pub/data/"
    rnxname = geok.statname_dt2rinexname(stat.lower(),date)
    url = os.path.join(urlserver , str(date.year) , geok.dt2doy(date) , 'data_30' , rnxname)
    return url

def rgp_ign_smn_1Hz_server(stat,date):
    urlserver = "ftp://rgpdata.ign.fr/pub/data/"

    urls = []

    for h in range(24):
        date_session = date
        date_session = date_session.replace(hour=h)

        print(date_session , 'session' , h)
        rnxname = geok.statname_dt2rinexname(stat.lower(),date_session ,
                                             session_a_instead_of_daily_session = 1)
        url = os.path.join(urlserver , str(date.year) , geok.dt2doy(date) ,
                           'data_1' , rnxname)

        urls.append(url)

    return urls

def unavco_server(stat,date):
    urlserver='ftp://data-out.unavco.org/pub/rinex'
    rnxname = geok.statname_dt2rinexname(stat.lower(),date)
    url = os.path.join(urlserver , 'obs', str(date.year) , geok.dt2doy(date) , rnxname)
    return url

def renag_server(stat,date):
    urlserver = "ftp://renag.unice.fr/data/"
    rnxname = geok.statname_dt2rinexname(stat.lower(),date)
    url = os.path.join(urlserver , str(date.year) , geok.dt2doy(date) , rnxname)
    return url

def orpheon_server(stat,date,user='',passwd=''):
    urlserver = "ftp://renag.unice.fr/"
    rnxname = geok.statname_dt2rinexname(stat.lower(),date)
    url = os.path.join(urlserver , str(date.year) , geok.dt2doy(date) , rnxname)
    return url,user,passwd

def ovsg_server(stat,date,user='',passwd=''):
    if dt.datetime(2009,1,1) <= date <= dt.datetime(2014,2,10):
        urlserver = "http://webobs.ovsg.univ-ag.fr/rawdata/GPS-GPSDATA.backtemp_20140210/"
    else:
        urlserver = "http://webobs.ovsg.univ-ag.fr/rawdata/GPS/GPSDATA/"
    rnxname = geok.statname_dt2rinexname(stat.lower(),date)
    url = os.path.join(urlserver , str(date.year) , geok.dt2doy(date) , 'rinex' , rnxname)
    return url,user,passwd

def geoaus_server(stat,date):
    """ Geosciences Australia
        ex : ftp://ftp.ga.gov.au/geodesy-outgoing/gnss/data/daily/2010/10063/ """
    urlserver = "ftp://ftp.ga.gov.au/geodesy-outgoing/gnss/data/daily/"
    rnxname = geok.statname_dt2rinexname(stat.lower(),date)
    url = os.path.join(urlserver , str(date.year) , date.strftime('%y') + geok.dt2doy(date) , rnxname)
    return url

def sonel_server(stat,date):
    """ex : ftp://ftp.sonel.org/gps/data/2015/001/ """
    urlserver = 'ftp://ftp.sonel.org/gps/data/'
    rnxname = geok.statname_dt2rinexname(stat.lower(),date)
    url = os.path.join(urlserver,str(date.year),geok.dt2doy(date),rnxname)
    return url



def orbclk_cddis_server(date,center='igs', sp3clk = 'sp3', repro=0, mgex=False,
                        longname = False, force_weekly_file=False):
    """
    longname is experimental and only for MGEX yet !! (180426)
    """
    urlserver='ftp://cddis.gsfc.nasa.gov/pub/gps/products/'
    if mgex:
        urlserver = urlserver + 'mgex/'
    if repro == 0:
        rep_fldr = ''
    else:
        rep_fldr = 'repro' + str(repro)
    if repro != 0:
        center     = list(center)
        center[-1] = str(repro)
        center = ''.join(center)
    if center in ("cod","cof","co2","cf2") and sp3clk == "sp3":
        print("INFO : CODE orbit extension changed to eph")
        sp3clk = "eph"

    # date definition
    week, day = geok.dt2gpstime(date)

    if force_weekly_file == False:
        pass

    elif type(force_weekly_file) is str or type(force_weekly_file) is int:
        print("INFO : The weekly file will be downloaded (DoW =",force_weekly_file,")")
        print("       Check force_weekly_file option if you don't want it")
        day = force_weekly_file

    elif sp3clk == "erp" and force_weekly_file == True:
        print("INFO : The weekly file (DoW = 7) will be downloaded for ERP")
        day = 7

    elif sp3clk == "sum" and force_weekly_file == True:
        print("INFO : The weekly file (DoW = 7) will be downloaded for SUM")
        day = 7

    if not longname: # e.g. gbm19903.sp3.Z
        if not 'igu' in center:
            orbname = center + str(week).zfill(4)  + str(day) +'.'+ sp3clk +'.Z'
            url = os.path.join(urlserver, str(week).zfill(4) ,rep_fldr,orbname)
        else:
            igusuffix = center[-2:]
            center    = center[:-2]
            orbname = center + str(week).zfill(4)  + str(day)  + '_' + igusuffix +'.'+ sp3clk + '.Z'
            url = os.path.join(urlserver, str(week).zfill(4) ,rep_fldr,orbname)
    else: # e.g. COD0MGXFIN_20180580000_01D_05M_ORB.SP3.gz
        if len(center) == 3:
            center = center.upper() + "0"
        else:
            center = center.upper()

        datelong = date.strftime("%Y%j")

        if "SHA" in center:
            if sp3clk == "sp3":
                sp3clk_long = "_01D_15M_ORB.SP3"
            elif sp3clk == "clk":
                sp3clk_long = "_01D_05M_CLK.CLK"
            elif sp3clk == "erp":
                sp3clk_long = "_03D_12H_ERP.ERP"
            elif sp3clk == "bia":
                sp3clk_long = "_01D_01D_OSB.BIA"
            elif sp3clk == "snx":
                sp3clk_long = "_01D_000_SOL.SNX"

            orbname = center + "MGXRAP_" + datelong + "0000"  + sp3clk_long + ".gz"
        else:
            if sp3clk == "sp3":
                sp3clk_long = "_01D_05M_ORB.SP3"
            elif sp3clk == "clk":
                sp3clk_long = "_01D_30S_CLK.CLK"
            elif sp3clk == "erp":
                sp3clk_long = "_03D_12H_ERP.ERP"
            elif sp3clk == "bia":
                sp3clk_long = "_01D_01D_OSB.BIA"
            elif sp3clk == "snx":
                sp3clk_long = "_01D_000_SOL.SNX"

            orbname = center + "MGXFIN_" + datelong + "0000"  + sp3clk_long + ".gz"

        url = os.path.join(urlserver, str(week).zfill(4) ,rep_fldr,orbname)

    return url



def orbclk_ign_server(date,center='igs', sp3clk = 'sp3', repro=0, mgex=False,
                        longname = False,force_weekly_file=False):
    """
    longname is experimental and only for MGEX yet !! (180426)
    Must be merged with orbclk_cddis_server to make a big equivalent fct (180523)
    """
    urlserver='ftp://igs.ign.fr/pub/igs/products/'
    if mgex:
        urlserver = urlserver + 'mgex/'
    if repro == 0:
        rep_fldr = ''
    else:
        rep_fldr = 'repro' + str(repro)
    if repro != 0:
        center     = list(center)
        center[-1] = str(repro)
        center = ''.join(center)
    week, day = geok.dt2gpstime(date)


    if force_weekly_file == False:
        pass
    elif type(force_weekly_file) is str or type(force_weekly_file) is int:
        print("INFO : The weekly file will be downloaded (DoW =",force_weekly_file,")")
        day = force_weekly_file

    elif sp3clk == "erp" and force_weekly_file == True:
        print("INFO : The weekly file (DoW = 7) will be downloaded for ERP")
        day = 7

    elif sp3clk == "sum" and force_weekly_file == True:
        print("INFO : The weekly file (DoW = 7) will be downloaded for SUM")
        day = 7

    if not longname: # e.g. gbm19903.sp3.Z
        if not 'igu' in center:
            orbname = center + str(week).zfill(4)  + str(day) +'.'+ sp3clk +'.Z'
            url = os.path.join(urlserver, str(week).zfill(4) ,rep_fldr,orbname)
        else:
            igusuffix = center[-2:]
            center    = center[:-2]
            orbname = center + str(week).zfill(4)  + str(day)  + '_' + igusuffix +'.'+ sp3clk + '.Z'
            url = os.path.join(urlserver, str(week).zfill(4) ,rep_fldr,orbname)
    else: # e.g. COD0MGXFIN_20180580000_01D_05M_ORB.SP3.gz
        if len(center) == 3:
            center = center.upper() + "0"
        else:
            center = center.upper()

        datelong = date.strftime("%Y%j")

        if "SHA" in center:
            if sp3clk == "sp3":
                sp3clk_long = "_01D_15M_ORB.SP3"
            elif sp3clk == "clk":
                sp3clk_long = "_01D_05M_CLK.CLK"
            elif sp3clk == "erp":
                sp3clk_long = "_03D_12H_ERP.ERP"
            elif sp3clk == "bia":
                sp3clk_long = "_01D_01D_OSB.BIA"

            orbname = center + "MGXRAP_" + datelong + "0000"  + sp3clk_long + ".gz"
        else:
            if sp3clk == "sp3":
                sp3clk_long = "_01D_05M_ORB.SP3"
            elif sp3clk == "clk":
                sp3clk_long = "_01D_30S_CLK.CLK"
            elif sp3clk == "erp":
                sp3clk_long = "_03D_12H_ERP.ERP"
            elif sp3clk == "bia":
                sp3clk_long = "_01D_01D_OSB.BIA"

            orbname = center + "MGXFIN_" + datelong + "0000"  + sp3clk_long + ".gz"

        url = os.path.join(urlserver, str(week).zfill(4) ,rep_fldr,orbname)

    return url




def orbclk_igscb_server(date,center='gfz', sp3clk = 'sp3',repro=0):
    urlserver='ftp://igscb.jpl.nasa.gov/pub/product/'
    if repro == 0:
        rep_fldr = ''
    elif repro == 3:
        rep_fldr = 'repro' + str(repro)
    if repro != 0:
        center     = list(center)
        center[-1] = str(repro)
        center = ''.join(center)
    week, day = geok.dt2gpstime(date)
    if not 'igu' in center:
        orbname = center + str(week).zfill(4) + str(day) +'.' + sp3clk + '.Z'
        url = os.path.join(urlserver,str(week).zfill(4) ,rep_fldr,orbname)
    else:
        igusuffix = center[-2:]
        center    = center[:-2]
        orbname = center + str(week).zfill(4) + str(day)  + '_' + igusuffix +'.' + sp3clk + '.Z'
        url = os.path.join(urlserver,str(week).zfill(4) ,rep_fldr,orbname)
    return url

def orbclk_gfz_local_server(date,center='gfz', sp3clk='sp3',repro=0):
    if repro == 0:
        urlserver = '/dsk/igs_archive/IGS/SAVE_PROD_1d/'
    elif repro == 3:
        urlserver = '/dsk/repro3/ARCHIVE/IGS/SAVE_PROD_1d/'
    else:
        print("ERR : check the repro !!!")
        raise Exception

    week, day = geok.dt2gpstime(date)
    
    orbname = center + str(week).zfill(4) + str(day) +'.' + sp3clk + '.Z'
    url = os.path.join(urlserver,str(week).zfill(4),orbname)
    
    return url

def downloader(url,savedir,force = False,
               check_if_file_already_exists_uncompressed=True):
    """
    general function to download a file

    INTERNAL_FUNCTION
    """
#    print url
    if type(url) is tuple:
        need_auth = True
        username = url[1]
        password = url[2]
        url = url[0]
    else:
        need_auth = False

    url_print = str(url)

    rnxname = os.path.basename(url)

    potential_exisiting_files_list = [os.path.join(savedir , rnxname)]

    if check_if_file_already_exists_uncompressed:
        potential_exisiting_files_list.append(os.path.join(savedir , rnxname.replace(".gz","")))
        potential_exisiting_files_list.append(os.path.join(savedir , rnxname.replace(".Z","")))
        potential_exisiting_files_list = list(set(potential_exisiting_files_list))

    for f in potential_exisiting_files_list:
        if os.path.isfile(f) and (not force):
            print("INFO :", os.path.basename(f) , "already exists locally ;)")
            return None
        
    ##### LOCAL FILE (particular case for GFZ)
    if os.path.isfile(url):
        print("INFO : downloader : the is a local file, a simple copy will be used")
        print("       URL :",url)
        shutil.copy(url,savedir)
    
    ##### REMOTE FILE (General case)
    elif ("ftp" in url) or ("http" in url) :
        # managing a authentification
        if need_auth:
            if 'ftp://' in url: # FTP with Auth
                url = url.replace('ftp://','ftp://' + username + ':' + password + '@')
                opener = urllib.request.build_opener()
            else: # HTTP with Auth
                password_mgr = urllib.request.HTTPPasswordMgrWithDefaultRealm()
                top_level_url = url
                password_mgr.add_password(None, top_level_url, username, password)
                handler = urllib.request.HTTPBasicAuthHandler(password_mgr)
                # create "opener" (OpenerDirector instance)
                opener = urllib.request.build_opener(handler)
        else: # FTP or HTTP without Auth
            opener = urllib.request.build_opener()
    
        # use the opener to fetch a URL
        try:
            f = opener.open(url)
        except (urllib.error.HTTPError , urllib.error.URLError):
            print("WARN :",rnxname,"not found on server :(")
            print(url_print)
            return None
        print("INFO :" , rnxname ," found on server :)")
        data = f.read()
        if not os.path.exists(savedir):
            os.makedirs(savedir)
        with open(os.path.join(savedir , rnxname), "wb") as code:
            code.write(data)
    else:
        print("ERR : something goes wrong with the URL")
        print("     ", url)


    # effective downloading (old version)
#    try:
#        f = urllib2.urlopen(url)
#    except urllib2.HTTPError:
#        print rnxname," not found :("
#        print url
#        return None
#    print rnxname," found :)"
#    data = f.read()
#    if not os.path.exists(savedir):
#        os.makedirs(savedir)
#    with open(os.path.join(savedir , rnxname), "wb") as code:
#        code.write(data)
    return None

def start_end_date_easy(start_year,start_doy,end_year,end_doy):
    start = geok.doy2dt(start_year,start_doy)
    end   = geok.doy2dt(end_year,end_doy)
    return start , end

def effective_save_dir(parent_archive_dir,stat,date,archtype ='stat'):
    """
    INTERNAL_FUNCTION

    archtype =
        stat
        stat/year
        stat/year/doy
        year/doy
        year/stat
        week/dow
        OR only '/' for a dirty saving in the parent folder
        ... etc ... """
    if archtype == '/':
        return parent_archive_dir

    out_save_dir = parent_archive_dir
    fff = archtype.split('/')
    year = str(date.year)
    doy = geok.dt2doy(date)
    week, dow = geok.dt2gpstime(date)
    for f in fff:
        out_save_dir = os.path.join(out_save_dir,eval(f))
    return out_save_dir

def effective_save_dir_orbit(parent_archive_dir,calc_center,date,
                             archtype ='year/doy/'):
    """
    INTERNAL_FUNCTION

    archtype =
        stat
        stat/year
        stat/year/doy
        year/doy
        year/stat
        week/dow
        wkwwww : use a GFZ's CF-ORB wk<wwww> naming
        OR only '/' for a dirty saving in the parent folder
        ... etc ...
    """
    if archtype == '/':
        return parent_archive_dir

    out_save_dir = parent_archive_dir
    fff = archtype.split('/')
    year = str(date.year)
    doy = geok.dt2doy(date)
    week, dow = geok.dt2gpstime(date)

    for f in fff:
        if "wkwwww" in f:
            f_evaluated = "wk" + str(week).zfill(4)
        else:
            f_evaluated = eval(f)
        out_save_dir = os.path.join(out_save_dir,f_evaluated)
    return out_save_dir


def downloader_wrap(intup):
    downloader(*intup)
    return None

def multi_downloader_rinex(statdico,archive_dir,startdate,enddate,
                           archtype ='stat',parallel_download=4,
                           sorted_mode=False,user='',passwd=''):
    """
    Parameters
    ----------
    statdico : dict
        a statdico is a dictionary associating Archives Centers to list of stations

        Exemple:
            >>> statdico['archive center 1'] = ['STA1','STA2','STA3', ...]
            >>> statdico['archive center 2'] = ['STA2','STA1','STA4', ...]

        the supported archive center are (july 2015):
            igs (cddis center)

            igs_garner (for the garner center, but not very reliable)

            rgp (St Mandé center)

            rgp_1Hz (all the 24 hourly rinex for the day will be downloaded)

            renag

            ovsg

            unavco

            sonel

            geoaus (Geosciences Australia)

            nav or brdc as archive center allows to download nav files (using 'BRDC' as station name) from the CDDIS server


    archtype : str
        string describing how the archive directory is structured, e.g :

            stat

            stat/year

            stat/year/doy

            year/doy

            year/stat

            week/dow/stat

            ... etc ...


    sorted_mode : bool
        if False:
            using the map multiprocess fct so the download order will
            be scrambled
        if True:
            using the apply multiprocess fct so the download order will be
            in the chrono. order
        The scrambled (False) is better, bc. it doesn't create zombies processes


    user & passwd : str
        user & password for a locked server


    Returns
    -------
    url_list : list of str
        list of URLs

    savedir_list : list of str
        list of downloaded products paths
    """

    curdate = startdate

    pool = mp.Pool(processes=parallel_download)

    urllist = []
    savedirlist = []


    print("generating the list of potential RINEXs ...")
    while curdate <= enddate:
        for netwk , statlis in list(statdico.items()):
            for stat in statlis:
                stat = stat.lower()
                mode1Hz = False

                if netwk == 'igs':
                    url = igs_cddis_server(stat,curdate)
                elif netwk == 'igs_garner':
                    url = igs_garner_server(stat,curdate)
                elif netwk == 'rgp':
                    url = rgp_ign_smn_server(stat,curdate)
                elif netwk == 'rgp_1Hz':
                    urls = rgp_ign_smn_1Hz_server(stat,curdate)
                    mode1Hz = True
                elif netwk == 'renag':
                    url = renag_server(stat,curdate)
                elif netwk == 'orpheon':
                    url = orpheon_server(stat,curdate,user,passwd)
                elif netwk == 'ovsg':
                    url = ovsg_server(stat,curdate,user,passwd)
                elif netwk == 'unavco':
                    url = unavco_server(stat,curdate)
                elif netwk == 'sonel':
                    url = sonel_server(stat,curdate)
                elif netwk == 'geoaus':
                    url = geoaus_server(stat,curdate)
                elif netwk in ('nav' , 'brdc'):
                    url = igs_cddis_nav_server(stat,curdate)
                else:
                    print('WARN : unkwn server dic in the dico, skip ...')
                    continue

                savedir = effective_save_dir(archive_dir,stat,curdate,archtype)

                if not mode1Hz:
                    urllist.append(url)
                    savedirlist.append(savedir)
                else:
                    urllist = urllist + urls
                    savedirlist = savedirlist + [savedir] * len(urls)

        curdate = curdate + dt.timedelta(days=1)

    #savedirlist = [x for (y,x) in sorted(zip(urllist,savedirlist))]
    #urllist     = sorted(urllist)
    
    urllist,savedirlist = genefun.sort_binom_list(urllist,savedirlist)

    print(" ... done")
    print(len(urllist),"potential RINEXs")

    if sorted_mode:
        results = [pool.apply_async(downloader , args=(u,sd)) for u,sd in zip(urllist,savedirlist)]
    else:
        _ = pool.map(downloader_wrap,list(zip(urllist,savedirlist)))


    localfiles_lis = []
    skiped_url = 0
    for url , savedir in zip(urllist,savedirlist):
        try:
            localfile = os.path.join(savedir,os.path.basename(url))
            if os.path.isfile(localfile):
                localfiles_lis.append(localfile)
        except:
            # because of a weird error
            # i = p.rfind('/') + 1
            # AttributeError: 'tuple' object has no attribute 'rfind'
            skiped_url += 1
            continue

    pool.close()
    print('INFO : ' , skiped_url , 'returned url skiped because of a weird error, but it is not important :)')
    return localfiles_lis

#    return zip(urllist,savedirlist)


def orbclk_long2short_name(longname_filepath_in,rm_longname_file=True,
                           center_id_last_letter=None,
                           center_manual_short_name=None):
    """
    Rename a long naming new convention IGS product file to the short old
    convention
    Naming will be done automaticcaly based on the 3 first charaters of the
    long AC id
    e.g. CODE => cod, GRGS => grg, NOAA => noa ...

    Parameters
    ----------
    longname_filepath_in : str
        Full path of the long name product file

    rm_longname_file : bool
        Remove the original long name product file

    center_id_last_letter : str
        replace the last letter of the short AC id by another letter
        (see note below)
        
    center_manual_short_name : str
        replace completely the long name with this one
        overrides center_id_last_letter

    Returns
    -------
    shortname_filepath : str
        Path of the  short old-named product file

    Note
    ----
    if you rename MGEX orbits, we advise to set
    center_id_last_letter="m"
    the AC code name will be changed to keep a MGEX convention
    (but any other caracter can be used too)

    e.g. for Bern's products, the long id is CODE

    if center_id_last_letter=None, it will become cod,
    if center_id_last_letter=m, it will become com

    """

    longname_basename = os.path.basename(longname_filepath_in)
    longname_dirname  = os.path.dirname(longname_filepath_in)

    center = longname_basename[:3]


    if center_manual_short_name:
        center = center_manual_short_name
    elif center_id_last_letter:
        center_as_list = list(center)
        center_as_list[-1] = center_id_last_letter
        center = "".join(center_as_list)

    yyyy   = int(longname_basename.split("_")[1][:4])
    doy    = int(longname_basename.split("_")[1][4:7])

    import geodetik as geok

    day_dt = geok.doy2dt(yyyy,doy)

    wwww , dow = geok.dt2gpstime(day_dt)

    shortname_prefix = center.lower() + str(wwww) + str(dow)

    if   "SP3" in longname_basename:
        shortname = shortname_prefix + ".sp3"
    elif "CLK" in longname_basename:
        shortname = shortname_prefix + ".clk"
    elif "ERP" in longname_basename:
        shortname = shortname_prefix + ".erp"
    elif "BIA" in longname_basename:
        shortname = shortname_prefix + ".bia"

    shortname_filepath = os.path.join(longname_dirname , shortname)

    shutil.copy2(longname_filepath_in , shortname_filepath)

    if rm_longname_file:
        print("INFO : remove " , longname_filepath_in)
        os.remove(longname_filepath_in)

    return shortname_filepath


def multi_downloader_orbs_clks(archive_dir,startdate,enddate,calc_center='igs',
                            sp3clk='sp3',archtype ='year/doy',parallel_download=4,
                            archive_center='cddis',repro=0,sorted_mode=False,
                            force_weekly_file=False, return_also_uncompressed_files=True):

    """
    Download IGS products. Can manage MGEX products too
    (see archive_center argument)

    Parameters
    ----------
    archive_dir : str
        Parent archive directory where files will be stored

    startdate & enddate : datetime
        Start and End of the wished period

    calc_center : str or list of str
        calc_center can be a string or a list, describing the calc center
        e.g. 'igs','grg','cod','jpl' ...

    sp3clk : str
        Product type, can handle :

            'clk'

            'clk_30s'

            'sp3'

            'erp'

            'bia'

    archive_center : str
        server of download, "regular" IGS or MGEX, can handle :

            'cddis'

            'cddis_mgex'

            'cddis_mgex_longname'

            'ign_mgex'

            'ign_mgex_longname'
            
            'gfz_local'

    archtype: str
        string describing how the archive directory is structured, e.g :

            stat

            stat/year

            stat/year/doy

            year/doy

            year/stat

            week/dow/stat

            ... etc ...

    repro : int
        number of the IGS reprocessing (0 = routine processing)

    sorted_mode : bool
        if False:
            using the map multiprocess fct so the download order will
            be scrambled
        if True:
            using the apply multiprocess fct so the download order will be
            in the chrono. order
        The scrambled (False) is better, bc. it doesn't create zombies processes

    Returns
    -------
    localfiles_lis : list of str
        list of downloaded products paths
    """

    if type(calc_center) is str:
        calc_center =  [ ''.join([calc_center]) ] # POURQUOI CETTE LIGNE ?? (150717)
#        calc_center =  [calc_center]

    pool = mp.Pool(processes=parallel_download)
    urllist = []
    savedirlist = []

    if sp3clk == 'clk':
        typ = 'Clocks'
    elif sp3clk == 'clk_30s':
        typ = 'Clocks (30s)'
    elif sp3clk == 'sp3':
        typ = 'Orbits'
    elif sp3clk == 'erp':
        typ = 'ERP'
    elif sp3clk == 'snx':
        typ = 'SINEXs'
    elif sp3clk == 'sum':
        typ = 'SUM files'
    elif sp3clk == 'bia':
        typ = 'ISBs'
    else:
        typ = '????'

    print("generating the list of potential " + typ + " ...")
    
    for cc in calc_center:
            curdate = startdate
            while curdate <= enddate:
                if re.search("igs([0-9][0-9]|yy|YY)P",cc):
                    cc = "igs" + str(curdate.year)[2:] + "P"
                    print("INFO : IGS reference frame snx/ssc, correcting the year : ",cc)

                url = ''
                if archive_center == 'cddis':
                    url = orbclk_cddis_server(curdate,cc,repro=repro,sp3clk=sp3clk,
                                              force_weekly_file=force_weekly_file)
                elif archive_center == 'cddis_mgex':
                    url = orbclk_cddis_server(curdate,cc,repro=repro,sp3clk=sp3clk,
                                              mgex=True,force_weekly_file=force_weekly_file)
                elif archive_center == 'cddis_mgex_longname':
                    url = orbclk_cddis_server(curdate,cc,repro=repro,
                                              sp3clk=sp3clk,mgex=True,longname=True,
                                              force_weekly_file=force_weekly_file)
                elif archive_center == 'ign_mgex':
                    url = orbclk_ign_server(curdate,cc,repro=repro,sp3clk=sp3clk,
                                              mgex=True,force_weekly_file=force_weekly_file)
                elif archive_center == 'ign_mgex_longname':
                    url = orbclk_ign_server(curdate,cc,repro=repro,
                                              sp3clk=sp3clk,mgex=True,longname=True,
                                              force_weekly_file=force_weekly_file)
                elif archive_center == 'gfz_local':
                    url = orbclk_gfz_local_server(curdate,cc,repro=repro,
                                              sp3clk=sp3clk)

                else:
                    print('ERR : Wrong archive_center name !!! :' , archive_center)
                urllist.append(url)
                savedir = effective_save_dir_orbit(archive_dir,cc,curdate,archtype)
                savedirlist.append(savedir)
                curdate = curdate + dt.timedelta(days=1)

    savedirlist = [x for (y,x) in sorted(zip(urllist,savedirlist))]
    urllist = sorted(urllist)

    print(" ... done")
    print(len(urllist),"potential " + typ)

#    if sorted_mode:
#        print 'aaaa'
#        results = [pool.apply_async(downloader , args=(u,sd)) for u,sd in zip(urllist,savedirlist)]
#    else:
#        print 'bb'
#        _ = pool.map(downloader_wrap,zip(urllist,savedirlist))
#
#
    if not sorted_mode:
        _ = pool.map(downloader_wrap,list(zip(urllist,savedirlist)))
    else:
         results = [pool.apply_async(downloader , args=(u,sd)) for u,sd in zip(urllist,savedirlist)]

    localfiles_lis = []
    
    if not return_also_uncompressed_files:
        for url , savedir in zip(urllist,savedirlist):
            localfile = os.path.join(savedir,os.path.basename(url))
            if os.path.isfile(localfile):
                localfiles_lis.append(localfile)
    else:
    
        for url , savedir in zip(urllist,savedirlist):

            localfile = os.path.join(savedir,os.path.basename(url))        
            
            potential_exisiting_files_list = [localfile]
            potential_exisiting_files_list.append(localfile.replace(".gz",""))
            potential_exisiting_files_list.append(localfile.replace(".Z",""))
            potential_exisiting_files_list = list(set(potential_exisiting_files_list))
            
            for potential_exisiting_file in potential_exisiting_files_list:
                if os.path.isfile(potential_exisiting_file):
                    localfiles_lis.append(potential_exisiting_file)
    
    pool.close()
    return localfiles_lis

def multi_finder_rinex(main_dir,rinex_types=('o','d','d.Z','d.z'),
                       specific_stats = [] ):
    """
    from a main_dir, find all the rinexs in this folder and his subfolder

    (corresponding to the rinex_types)

    and return a list of the found rinexs

    is very similar with geodetik.rinex_lister and  gins_runner.get_rinex_list

    But this one is the most elaborated , must be used in priority !!!
    """
    files_raw_lis , _ = genefun.walk_dir(main_dir)

    yylis = [str(e).zfill(2) for e in list(range(80,100)) + list(range(0,dt.datetime.now().year - 2000 + 1))]

    rinex_lis = []


    for f in files_raw_lis:
        for rnxext in rinex_types:
            for yy in yylis:
                if f.endswith(yy + rnxext):
                    rinex_lis.append(f)


    # CASE FOR specific_stats
    if len(specific_stats) > 0:
        rinex_lis2 = []
        for stat in specific_stats:
            for rnx in rinex_lis:
                if stat in os.path.basename(rnx):
                    rinex_lis2.append(rnx)
        rinex_lis = rinex_lis2

    print('INFO : ' , len(rinex_lis) , 'RINEXs found')

    return rinex_lis


def multi_archiver_rinex(rinex_lis,parent_archive_dir,archtype='stat',
                         move=True,force_mv_or_cp=True):
    """
    from rinex_lis, a list of rinex (generated by the function
    multi_finder_rinex)

    move (if move=True) of copy (if move=False) those rinexs in the
    parent_archive_dir according to the archtype,
    string describing how the archive directory is structured, e.g :
            stat

            stat/year

            stat/year/doy

            year/doy

            year/stat

            week/dow/stat

            ... etc ...
    """

    mv_cnt   = 0
    skip_cnt = 0
    print('INFO : RINEXs as input :' , len(rinex_lis))

    if move:
        mv_fct = genefun.move
    else:
        mv_fct = genefun.copy2

    for rnx in rinex_lis:
        date = geok.rinexname2dt(rnx)
        rnxname = os.path.basename(rnx)
        stat = rnxname[0:4]
        savedir = effective_save_dir(parent_archive_dir, stat, date, archtype)
        genefun.create_dir(savedir)

        if not force_mv_or_cp and os.path.isfile(os.path.join(savedir,rnxname)):
            skip_cnt += 1
            continue
        else:
            mv_fct(rnx,savedir)
            mv_cnt += 1

    print('INFO : RINEXs skiped   :' , skip_cnt, '(because already existing)')
    print('INFO : RINEXs moved    :' , mv_cnt)

    return None


def find_IGS_products_files(parent_dir,File_type,ACs,date_start,date_end=None,
                            recursive_search=True,severe=True):
    """
    Find all product files in a parent folder which correspond to file type(s),
    AC(s) and date(s)

    Has to be improved for the new naming convention

    Parameters
    ----------
    parent_dir : str or list of str
        The parent directory (i.e. the archive) where files are stored
        can be a string (path of the archive) or a list of file paths
        (given by the function genefun.find_recursive) in order to gain time


    File_type : str or list of str
        File type(s) researched (sp3, erp, clk ...)
        can be a list of string for several file type paths
        or a string like 'sp3' if only one file type is researched

    ACs : 'all' or str or list of str
        AC(s) researched
        can be a list of string for several ACs
        or a string like 'gfz' if only one AC is researched
        if 'all', search for all the ACs

    date_start : dt.datetime or 2-tuple list of int
        begining of the time period researched
        can be a datetime
        or a 2-tuple (wwww,d) e.g. (1990,0)

    date_end : None or dt.datetime or 2-tuple list of int
        end of the time period researched
        can be a datetime
        or a 2-tuple (wwww,d) e.g. (1990,0)
        if None, then only date_start is researched
    
    severe : bool
        If True, raises an exception if something goes wrong
        
    Naming_conv : str or list of str

    Returns
    -------
    Files_select_cumul_list : list
        List of files found

    """
    
    ###### Prelim Checks   ##############
    if not os.path.exists(parent_dir):
        print("ERR : parent directory doesn't exist")
        print(parent_dir)
        if severe:
            raise Exception

    ###### Time management ##############
    ###### work internally with datetime
    ### date_start
    if type(date_start) is dt.datetime:
        date_start_ok = date_start
    else:
        date_start_ok = geok.gpstime2dt(*date_start)
    ### date_end
    if not date_end:
        date_end_ok = date_start_ok
    elif type(date_end) is dt.datetime:
        date_end_ok = date_end
    else:
        date_end_ok = geok.gpstime2dt(*date_end)
    ### generate time period with a while loop
    Dates_list = [date_start_ok]
    while Dates_list[-1] < date_end_ok:
        Dates_list.append(Dates_list[-1]  + dt.timedelta(days=1))

    Dates_wwwwd_list = [genefun.join_improved("",*geok.dt2gpstime(d)) for d in Dates_list]

    ###### File type / ACs management ##############

    if not genefun.is_iterable(File_type):
        File_type = [File_type]
    if not genefun.is_iterable(ACs):
        ACs = [ACs]

    ###### General file search management ##############
    if genefun.is_iterable(parent_dir):
        FILE_LIST = parent_dir
    elif recursive_search:
        # All the files are listed first
        FILE_LIST = genefun.find_recursive(parent_dir,".*",case_sensitive=False)
    else:
        FILE_LIST = glob.glob(parent_dir + "/*")


    ###### Regex Definition ##############
    
    regex_old_naming = True
    regex_new_naming = True
    regex_igs_tfcc_naming = True
    
    join_regex_and = lambda  L : "(" +  "|".join(L) + ")"
    
    #####
    ##### WORK IN PROGESS HERE FOR THE REGEX DEFINTION FOR THE NEW NAMING CONVENTION
    ##### AND THE FCT ARGUMENTS
    #####
    
    Re_patt_big_stk = []
    
    if regex_old_naming: 
        if ACs[0] == "all":
            re_patt_ac = "\w{3}"
        else:
            re_patt_ac = join_regex_and(ACs)
        re_patt_date   = join_regex_and(Dates_wwwwd_list)
        re_patt_filtyp = join_regex_and(File_type)
        re_patt_big_old_naming = ".*".join((re_patt_ac,re_patt_date,re_patt_filtyp))
        Re_patt_big_stk.append(re_patt_big_old_naming)

    if regex_igs_tfcc_naming:
        Dates_yy_list = list(set([str(geok.gpstime2dt(int(e[0:4]),int(e[4])).year)[2:] for e in Dates_wwwwd_list]))
        Dates_wwww_list = list(set([e[:-1] for e in Dates_wwwwd_list]))
        #Dates_wwww_dot_list = [e + "\." for e in Dates_wwww_list]
        re_patt_year = join_regex_and(Dates_yy_list)
        # 2x re_patt_date : because .sum doesn't the day
        re_patt_date = join_regex_and(Dates_wwwwd_list + Dates_wwww_list)
        re_patt_filtyp = "\." +  join_regex_and(File_type)

        re_patt_big_igs_tfcc_naming = "igs" + re_patt_year + "P" + re_patt_date + ".*" + re_patt_filtyp
        Re_patt_big_stk.append(re_patt_big_igs_tfcc_naming)

    re_patt_big = join_regex_and(Re_patt_big_stk)
        
    print("INFO : REGEX researched :")    
    print(re_patt_big)
    
    ###### Specific file search management ##############
    Files_select_list = []
    for fil in FILE_LIST:
        if re.search(re_patt_big,os.path.basename(fil),re.IGNORECASE):
            Files_select_list.append(fil)
            
    if len(Files_select_list) == 0:
        print("ERR : no products found")
        if severe:
            raise Exception

    print("INFO : number of files found :",len(Files_select_list))    
    print(re_patt_big)
    
    return Files_select_list

#  _____  _____ _   _ ________   __   _____       _ _ _
# |  __ \|_   _| \ | |  ____\ \ / /  / ____|     | (_) |
# | |__) | | | |  \| | |__   \ V /  | (___  _ __ | |_| |_ ___ _ __
# |  _  /  | | | . ` |  __|   > <    \___ \| '_ \| | | __/ _ \ '__|
# | | \ \ _| |_| |\  | |____ / . \   ____) | |_) | | | ||  __/ |
# |_|  \_\_____|_| \_|______/_/ \_\ |_____/| .__/|_|_|\__\___|_|
#                                          | |
#                                          |_|

def check_if_compressed_rinex(rinex_path):
    boolout = bool(re.search('.*((d|o)\.(Z)|(gz))$', rinex_path))
    return boolout

def crz2rnx(rinex_path,outdir='',force=True,path_of_crz2rnx='CRZ2RNX'):
    """ assuming that CRZ2RNX is in the system PATH per default """

    if not os.path.isfile(rinex_path):
        raise Exception( rinex_path + ' dont exists !')

    if force:
        forcearg = '-f'
    else:
        forcearg = ''

    curdir = os.getcwd()
    command = path_of_crz2rnx + ' -c ' + forcearg + ' ' + rinex_path
    if outdir == '':
        outdir = os.path.dirname(rinex_path)

    out_rinex_name_splited = os.path.basename(rinex_path).split('.')
    out_rinex_name = '.'.join(out_rinex_name_splited[:-1])
    out_rinex_name = out_rinex_name[:-1] + 'o'
    out_rinex_path = os.path.join(outdir,out_rinex_name)

    if os.path.isfile(out_rinex_path) and not force:
        print("INFO : crz2rnx : ", out_rinex_path , 'already exists, skiping ...')
    else:
        os.chdir(outdir)
        #stream = os.popen(command)

        proc = subprocess.Popen(command, shell=True,
                            stdout=subprocess.PIPE,executable='/bin/bash')
        status = proc.wait()

        if status != 0 or not os.path.isfile(out_rinex_path):
            print('ERR : crz2rnx : ',out_rinex_path,"not created ! :(")
        else:
            print('INFO : crz2rnx : ',out_rinex_path, 'created :)')
        #    print stream.read()
        os.chdir(curdir)
    return out_rinex_path



def crz2rnx_bad(crinex_in_path,outdir='',force=True,path_of_crz2rnx='CRZ2RNX'):
    """ assuming that CRZ2RNX is in the system PATH per default """

    if not os.path.isfile(crinex_in_path):
        raise Exception( crinex_in_path  + ' dont exists !')

    if force:
        forcearg = '-f'
    else:
        forcearg = ''

    command = path_of_crz2rnx + ' ' + forcearg + ' ' + crinex_in_path
    convert_dir = os.path.dirname(crinex_in_path)

    if outdir == '':
        outdir = convert_dir

    out_rinex_name_splited = os.path.basename(crinex_in_path).split('.')
    out_rinex_name = '.'.join(out_rinex_name_splited[:-1])
    out_rinex_name = out_rinex_name[:-1] + 'o'

    out_rinex_path  = os.path.join(outdir,out_rinex_name)
    temp_rinex_path = os.path.join(convert_dir,out_rinex_name)

    if os.path.isfile(out_rinex_path) and not force:
        print("INFO : crz2rnx : ", out_rinex_path , 'already exists, skiping ...')
    else:
        stream = os.popen(command)
        print('output in',outdir)
        try:
            shutil.move(temp_rinex_path, outdir)
        except:
            pass

        #    print stream.read()
#    os.chdir(curdir)
    return out_rinex_path

def rnx2crz(rinex_path,outdir='',force=True,path_of_rnx2crz='RNX2CRZ'):
    """ assuming that RNX2CRZ is in the system PATH per default """

    if not os.path.isfile(rinex_path):
        raise Exception( rinex_path + ' dont exists !')

    if force:
        forcearg = '-f'
    else:
        forcearg = ''

    curdir = os.curdir
    command = path_of_rnx2crz + ' ' + forcearg + ' ' + rinex_path
    if outdir == '':
        outdir = os.path.dirname(rinex_path)

    out_rinex_name_o = os.path.basename(rinex_path)
    out_rinex_name_d_Z = out_rinex_name_o[:-1] + 'd.Z'
    out_rinex_path = os.path.join(outdir,out_rinex_name_d_Z)

    if os.path.isfile(out_rinex_path) and not force:
        print("INFO : rnx2crz : ", out_rinex_path , 'already exists, skiping ...')
    else:
        os.chdir(outdir)
        stream = os.popen(command)
        print(command,', output in',outdir)
        #    print stream.read()
    os.chdir(curdir)
    return out_rinex_path

def rinex_regex(compressed=True,compiled=False):
    if compressed:
        regexstr = "....[0-9]{3}.\.[0-9]{2}((d\.(Z)|(gz))|(o)|(d))"
    else:
        regexstr = "....[0-9]{3}.\.[0-9]{2}o"

    if compiled:
        return re.compile(regexstr)
    else:
        return regexstr


def rinex_regex_new_name(compressed=True,compiled=False):
    if compressed:
        regexstr = ".{4}[0-9]{2}.{3}_(R|S|U)_[0-9]{11}_[0-9]{2}\w_[0-9]{2}\w_\w{2}\.\w{3}\.gz"
    else:
        regexstr = ".{4}[0-9]{2}.{3}_(R|S|U)_[0-9]{11}_[0-9]{2}\w_[0-9]{2}\w_\w{2}\.\w{3}"

    if compiled:
        return re.compile(regexstr)
    else:
        return regexstr


def rinex_read_epoch(input_rinex_path_or_string,interval_out=False,
                    add_tzinfo=False,out_array=True):
    """
    input_rinex_path_or_string :

        can be the path of a RINEX or directly the RINEX content as a string

    161019 : dirty copier coller de rinex start end
    """
    epochs_list = []
    rinex_60sec = False



    if  genefun.is_iterable(input_rinex_path_or_string):
        Finp = input_rinex_path_or_string
    elif os.path.isfile(input_rinex_path_or_string):
        Finp = open(input_rinex_path_or_string)
    else:
        Finp = input_rinex_path_or_string


    for line in Finp:
        epoch=re.search('^ {1,2}([0-9]{1,2} * ){5}',line)
        if epoch:
            fraw = line.split() #fraw est la pour gerer les fracion de sec ...
            fraw = fraw[0:6]
            fraw = [float(e) for e in fraw]
            f = [int(e) for e in fraw]
            msec = (fraw[5] - np.floor(fraw[5]))
            msec = np.round(msec,4)
            msec = int(msec * 10**6)
            f.append( msec )  # ajour des fractions de sec

            if f[0] < 50:
                f[0] = f[0] + 2000
            else:
                f[0] = f[0] + 1900
            if f[5] == 60: # cas particulier rencontré dans des rinex avec T = 60sec
                rinex_60sec = True
                f[5] = 59
            else:
                rinex_60sec = False
            try:
                epochdt = datetime.datetime(*f)
                if rinex_60sec:
                    epochdt = epochdt + dt.timedelta(seconds=1)
                epochs_list.append(epochdt)

            except:
                print("ERR : rinex_read_epoch : " , f)


    if add_tzinfo:
        epochs_list = [ e.replace(tzinfo=dateutil.tz.tzutc()) for e in epochs_list]

    if out_array:
        return np.array(epochs_list)
    else:
        return epochs_list

def same_day_rinex_check(rinex1,rinex2):
    if os.path.basename(rinex1)[4:7] == os.path.basename(rinex2)[4:7]:
        return True
    else:
        return False


def rinex_start_end(input_rinex_path,interval_out=False,
                    add_tzinfo=False,verbose = True,
                    safety_mode = True):
    """
    safety_mode :

        if the epoch reading fails (e.g. in case of a compressed RINEX)
        activate a reading of the header and the file name as backup.


    une liste d'epochs en début et fin de fichier
    => en trouver le min et le max
    NB : FAIRE UN FONCTION READ EPOCH A L'OCCAZ
    NBsuite : c'est fait au 161018 mais par contre c'est un dirty copier coller
    """
    epochs_list = []
    Head = genefun.head(input_rinex_path,500)
    epochs_list_head = rinex_read_epoch(Head,interval_out=interval_out,
                                        add_tzinfo=add_tzinfo,out_array=False)


    Tail =  genefun.tail(input_rinex_path,500)
    epochs_list_tail = rinex_read_epoch(Tail,interval_out=interval_out,
                                        add_tzinfo=add_tzinfo,out_array=False)

    epochs_list = epochs_list_head + epochs_list_tail

    if len(epochs_list) == 0:
        first_epoch = geok.rinexname2dt(input_rinex_path)
        alphabet = list(string.ascii_lowercase)

        if os.path.basename(input_rinex_path)[7] in alphabet:
            last_epoch = first_epoch + dt.timedelta(hours=1)
        else:
            last_epoch = first_epoch + dt.timedelta(hours=24,seconds=-1)
    else:
        first_epoch = np.min(epochs_list)
        last_epoch = np.max(epochs_list)

    if add_tzinfo:
        first_epoch = first_epoch.replace(tzinfo=dateutil.tz.tzutc())
        last_epoch = last_epoch.replace(tzinfo=dateutil.tz.tzutc())

    if verbose:
        print("first & last epochs : " , first_epoch , last_epoch)
    if not interval_out:
        return first_epoch , last_epoch
    else:
        interv_lis = np.diff(epochs_list)
        interv_lis = [e.seconds + e.microseconds * 10**-6 for e in interv_lis]
        interval   = genefun.most_common(interv_lis)
        print("interval : " , interval , last_epoch)

        #return interv_lis , epochs_list
        return first_epoch , last_epoch , interval

def rinex_session_id(first_epoch,last_epoch,full_mode=False):
    """
    full_mode:
        gives the letter of the starting session  & the length in hour
        of the session
    """

    alphabet = list(string.ascii_lowercase)
    rinex_length = last_epoch - first_epoch
    if rinex_length.seconds > 86000:
        rnx_interval_ext = '0'
    else:
        if not full_mode:
            rnx_interval_ext = alphabet[first_epoch.hour]
        else:
            rnx_interval_ext = alphabet[first_epoch.hour]  \
                               + str(int(rinex_length.seconds / 3600.))

    return rnx_interval_ext

#def rinex_spliter(input_rinex_path,output_directory,stat_out_name='',
#                  interval_size=24,compress=False,shift = 0):
#    """
#    if shift != 0:
#    the start/end of a session is shifted of shift minutes
#    """
#
#    if stat_out_name == '':
#        stat_out_name = os.path.basename(input_rinex_path)[0:4]
#
#    # check if the RINEX is compressed ...
#    bool_comp_rnx = check_if_compressed_rinex(input_rinex_path)
#    # ... if not crz2rnx !
#    if bool_comp_rnx:
#        input_rinex_path = crz2rnx(input_rinex_path)
#
#    inp_rinex_obj=open(input_rinex_path,'r+')
#    out_dir = output_directory # nom redondant mais j'ai la flemme d'aller corriger le nom de la variable
#
#    if not os.path.exists(out_dir):
#        os.makedirs(out_dir)
#    os.chdir(out_dir)
#
#    first_epoch , last_epoch = rinex_start_end(input_rinex_path)
#
#    # In this function, Date are truncated epochs
#    # (only the day if interval == 24 , + the hour else)
#    first_date = datetime.datetime(first_epoch.year,first_epoch.month,
#                                   first_epoch.day,first_epoch.hour)
#    last_date = datetime.datetime(last_epoch.year,last_epoch.month,
#                                  last_epoch.day,last_epoch.hour)+datetime.timedelta(hours=1)
#
#    print "first & last dates (truncated) : " , first_date , last_date
#
#    time_interval = genefun.get_interval(first_date,last_date,
#                                         datetime.timedelta(hours=interval_size))
#
#    alphabet = list(string.ascii_lowercase)
#
#    rinex_out_name_lis = []
#
#    for i,curr_date in enumerate(time_interval):
#        if interval_size == 24:
#            rnx_interval_ext = '0.'
#        else:
#            rnx_interval_ext = alphabet[curr_date.hour] + '.'
#
#        p = subprocess.Popen('',executable='/bin/bash', stdin=subprocess.PIPE , stdout=subprocess.PIPE , stderr=subprocess.PIPE)
#        command = 'teqc ' + '-O.mo ' + stat_out_name.upper() +' -st ' +  curr_date.strftime('%y%m%d%H%M%S') + ' +dh ' + str(interval_size) + ' ' + input_rinex_path
#        print command
#        rinex_out_name = stat_out_name + curr_date.strftime('%j') + rnx_interval_ext + curr_date.strftime('%y') + 'o'
#        err_log_name = curr_date.strftime('%j') + rnx_interval_ext + "err.log"
#        print rinex_out_name
#        stdout,stderr = p.communicate( command )
#        std_file = open(rinex_out_name, "w")
#        std_file.write(stdout)
#        std_file.close()
#        if stderr != '':
#            print err_log_name + " is not empty, must be checked !"
#            print stderr
#            err_file = open(err_log_name, "w")
#            err_file.write(stderr)
#            err_file.close()
#
#        rnx_splt_path = os.path.join(out_dir,rinex_out_name)
#        if compress:
#            rinex_out_final = rnx2crz(rnx_splt_path)
#            os.remove(rnx_splt_path)
#        else:
#            rinex_out_final = rnx_splt_path
#
#        rinex_out_name_lis.append(rinex_out_final)
#
#    return rinex_out_name_lis


def rinex_spliter(input_rinex_path,output_directory,stat_out_name='',
                  interval_size=24,compress=False,shift = 0,
                  inclusive = False,teqc_cmd='teqc'):
    """
    if shift != 0:
    the start/end of a session is shifted of shift minutes

    inclusive:
    delta of exaclty interval_size => add the 1st epoch of the next sess
    not inclusive:
    delta of interval_size - 1s
    """
    """
    Split RINEX file using teqc

    Parameters
    ----------
    input_rinex_path : str
        Path of the input RINEX

    output_directory : str
        Description param2
        
    stat_out_name : str
        Station name for the output RINEX
        
    interval_size : int
        Size of the splitted interval
    
    compress : bool
        compress the outputed RINEX    

    interval_size : int
        Size of the splitted interval 
        
        
    FINIR CE HEADER !!!!
                
    Returns
    -------
    out1 : float or int or str or dict or n-tuple or bool or list or numpy.array
        Description out1
    
    out2 : float or int or str or dict or n-tuple or bool or list or numpy.array
        Description out2
        
    Note
    ----
    Misc. Notes

    Source
    ------
    www.forum-source.com
    
    Examples
    --------
    >> answer
    42    
    """

    if stat_out_name == '':
        stat_out_name = os.path.basename(input_rinex_path)[0:4]

    # check if the RINEX is compressed ...
    bool_comp_rnx = check_if_compressed_rinex(input_rinex_path)

    # ... if not crz2rnx !
    if bool_comp_rnx:
        input_rinex_path = crz2rnx(input_rinex_path)

    inp_rinex_obj=open(input_rinex_path,'r+')
    out_dir = output_directory # nom redondant mais j'ai la flemme d'aller corriger le nom de la variable

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    os.chdir(out_dir)

    first_epoch , last_epoch = rinex_start_end(input_rinex_path)

    # In this function, Date are truncated epochs
    # (only the day if interval == 24 , + the hour else)
    first_date = datetime.datetime(first_epoch.year,first_epoch.month,
                                   first_epoch.day,first_epoch.hour)
    last_date = datetime.datetime(last_epoch.year,last_epoch.month,
                                  last_epoch.day,last_epoch.hour)+datetime.timedelta(hours=1)

    print("first & last dates (truncated) : " , first_date , last_date)

    time_interval = genefun.get_interval(first_date,last_date,
                                         datetime.timedelta(hours=interval_size))

    alphabet = list(string.ascii_lowercase)

    rinex_out_name_lis = []

    for i,curr_date in enumerate(time_interval):

        if shift != 0:
            print('WARN : rinex_spliter : shifted mode on, be careful !!!')


        if bool(shift) and i != 0:
            curr_date = curr_date + dt.timedelta(minutes = shift)
            interval_size_ope = interval_size
        elif bool(shift) and i == 0 :
            interval_size_ope = interval_size + float(shift) / 60.
        else:
            interval_size_ope = interval_size

        if not inclusive:
            interval_size_ope = interval_size_ope - 1/3600.

        if not bool(shift):
            if interval_size == 24:
                rnx_interval_ext = '0.'
            else:
                rnx_interval_ext = alphabet[curr_date.hour] + '.'
        else:
            if interval_size == 24:
                rnx_interval_ext = '0.'
            else:
                rnx_interval_ext = alphabet[curr_date.hour]  +'.'

        p = subprocess.Popen('',executable='/bin/bash', stdin=subprocess.PIPE , stdout=subprocess.PIPE , stderr=subprocess.PIPE)
        # - 1/3600. # one_sec2h
        command = teqc_cmd  + ' -O.mo ' + stat_out_name.upper() +' -st ' +  curr_date.strftime('%y%m%d%H%M%S') + ' +dh ' + str(interval_size_ope) + ' ' + input_rinex_path
        print(command)
        rinex_out_name = stat_out_name + curr_date.strftime('%j') + rnx_interval_ext + curr_date.strftime('%y') + 'o'
        err_log_name = curr_date.strftime('%j') + rnx_interval_ext + "err.log"
        print(rinex_out_name)
        stdout,stderr = p.communicate( command.encode() )
        std_file = open(rinex_out_name, "w")
        std_file.write(stdout.decode("utf-8"))
        std_file.close()
        if stderr:
            print(err_log_name + " is not empty, must be checked !")
            print(stderr.decode("utf-8"))
            err_file = open(err_log_name, "w")
            err_file.write(stderr.decode("utf-8"))
            err_file.close()

        rnx_splt_path = os.path.join(out_dir,rinex_out_name)
        if compress:
            rinex_out_final = rnx2crz(rnx_splt_path)
            os.remove(rnx_splt_path)
        else:
            rinex_out_final = rnx_splt_path

        rinex_out_name_lis.append(rinex_out_final)

    return rinex_out_name_lis



def teqc_qc(rinex_path,quick_mode=False,optional_args=''):
    """
    quick mode : reduced qc and no summary file written
    """
    if quick_mode:
        qcmode = "+qcq"
    else:
        qcmode = "+qc"

    # check if the RINEX is compressed ...
    bool_comp_rnx = check_if_compressed_rinex(rinex_path)

    # ... if not crz2rnx !
    if bool_comp_rnx:
        rinex_path_work = crz2rnx(rinex_path)
    else:
        rinex_path_work = rinex_path



    kommand = " ".join(("teqc" , qcmode , optional_args , rinex_path))

    print(kommand)

    proc = subprocess.Popen(kommand, shell=True,
                            stdout=subprocess.PIPE,executable='/bin/bash')
    status = proc.wait()

    if bool_comp_rnx:
        os.remove(rinex_path_work)

    if status:
        print("ERR : teqc qc : crash for " + rinex_path + ", code " + str(proc.poll()))
        return None

    output = proc.stdout.read()

    return output.decode('ASCII')

def rinex_renamer(input_rinex_path,output_directory,stat_out_name='',remove=False):

    if stat_out_name == '':
        stat_out_name = os.path.basename(input_rinex_path)[0:4]

    stat_out_name = stat_out_name.lower()

    inp_rinex_obj=open(input_rinex_path,'r+')
    out_dir = output_directory # nom redondant mais j'ai la flemme d'aller corriger le nom de la variable

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    os.chdir(out_dir)

    first_epoch , last_epoch = rinex_start_end(input_rinex_path)
    rnx_interval_ext = rinex_session_id(first_epoch,last_epoch) + '.'
    rinex_out_name = stat_out_name + first_epoch.strftime('%j') + rnx_interval_ext + first_epoch.strftime('%y') + 'o'

    print(rinex_out_name)

    output_rinex_path = os.path.join(out_dir,rinex_out_name)

    if input_rinex_path != output_rinex_path:
        print("INFO : copy of ", input_rinex_path , ' to ' , output_rinex_path)
        shutil.copy(input_rinex_path,output_rinex_path)

        if remove and os.isfile(output_rinex_path) :
            print("INFO : removing " , input_rinex_path)
            os.remove(input_rinex_path)
    else:
        print("INFO : " , input_rinex_path)
        print("and", output_rinex_path ,"are the same file")
        print("nothing's done ...")

    return output_rinex_path

#  _____ _______ _  ___      _____ ____
# |  __ \__   __| |/ / |    |_   _|  _ \
# | |__) | | |  | ' /| |      | | | |_) |
# |  _  /  | |  |  < | |      | | |  _ <
# | | \ \  | |  | . \| |____ _| |_| |_) |
# |_|  \_\ |_|  |_|\_\______|_____|____/
#

def read_conf_file(filein):
    outdic = collections.OrderedDict()
    for l  in open(filein):
        if l[0] == "#":
            continue
        if l.strip() == '':
            continue

        f = l.split('=')
        key = f[0]
        if len(f) == 1:
            val = ''
        else:
            try:
                val = f[1].split('#')[0]
            except IndexError:
                val = f[1]
        outdic[key.strip()] = val.strip()
    return outdic

def rtklib_run_from_rinex(rnx_rover,rnx_base,generik_conf,working_dir,
                          experience_prefix="",rover_auto_conf=False,
                          base_auto_conf=True,XYZbase=[0,0,0],outtype = 'auto',
                          calc_center='igs'):

    """
    auto_conf :
        read the header of the rinex and write the conf. file according to it
        if the mode is disabled, the antenna/rec part of the
        conf. file will be the same as the generic one

        NB : RTKLIB "core" have it's own reading header option.
        Thus, my advice is to disable the auto mode for the rover and leave
        ``ant1-postype=rinexhead`` in the generic conf file
        and enable it for the base with a XYZbase vector initialized or
        the good XYZ in header on the rinex (the XYZbase vector is prioritary
        over the base RINEX header)

        (but in anycase, prepars good rinex's headers ;) )

    outtype :
        'auto' (as defined in the generic config file) or
        'dms' 'deg' 'xyz' 'enu'
        can manage the upper case XYZ or FLH
    """

    # paths & files
    working_dir = genefun.create_dir(working_dir)
    temp_dir    = genefun.create_dir(os.path.join(working_dir,'TEMP_' + genefun.get_timestamp()))
    out_dir     = genefun.create_dir(os.path.join(working_dir,'OUTPUT'))

    # uncompressing rinex if compressed
    if check_if_compressed_rinex(rnx_rover):
        rnx_rover = crz2rnx(rnx_rover,temp_dir)
    if check_if_compressed_rinex(rnx_base):
        rnx_base  = crz2rnx(rnx_base,temp_dir)

    # RINEX START & END
    rov_srt, rov_end , rov_itv = rinex_start_end(rnx_rover,1)
    bas_srt, bas_end , bas_itv = rinex_start_end(rnx_base,1)

    # RINEX NAMES
    rov_name = os.path.basename(rnx_rover)[0:4]
    bas_name = os.path.basename(rnx_base)[0:4]

    # paths & files
    #temp_dir = os.path.join(working_dir,'TEMP')
    #if clean_temp_dir:
    #    shutil.rmtree(temp_dir)
    #    temp_dir = os.path.join(working_dir,'TEMP')
    #out_dir  = os.path.join(working_dir,'OUTPUT')

    srt_str = rov_srt.strftime("%Y_%j")
    exp_full_name = '_'.join((experience_prefix,rov_name,bas_name,srt_str))

    out_conf_fil   = os.path.join(out_dir,exp_full_name + '.conf')
    out_result_fil = os.path.join(out_dir,exp_full_name + '.out' )

    dicoconf = read_conf_file(generik_conf)

    if rover_auto_conf:
        Antobj_rov , Recobj_rov , Siteobj_rov , Locobj_rov = \
        gfc.read_rinex_2_dataobjts(rnx_rover)
        dicoconf['ant1-postype'] ='xyz'
        dicoconf['ant1-anttype'] = Antobj_rov.Antenna_Type
        dicoconf['ant1-pos1'] = Locobj_rov.X_coordinate_m
        dicoconf['ant1-pos2'] = Locobj_rov.Y_coordinate_m
        dicoconf['ant1-pos3'] = Locobj_rov.Z_coordinate_m
        dicoconf['ant1-antdelu'] = Antobj_rov.Up_Ecc
        dicoconf['ant1-antdeln'] = Antobj_rov.North_Ecc
        dicoconf['ant1-antdele'] = Antobj_rov.East_Ecc

    if not outtype.lower() == 'auto':
        dicoconf['out-solformat'] = outtype.lower()
        print('out-solformat' , dicoconf['out-solformat'])


    if base_auto_conf:
        Antobj_bas , Recobj_bas , Siteobj_bas , Locobj_bas = \
        gfc.read_rinex_2_dataobjts(rnx_base)
        dicoconf['ant2-postype'] ='xyz'
        dicoconf['ant2-anttype'] = Antobj_bas.Antenna_Type
        if XYZbase[0] != 0:
            dicoconf['ant2-pos1'] = XYZbase[0]
            dicoconf['ant2-pos2'] = XYZbase[1]
            dicoconf['ant2-pos3'] = XYZbase[2]
        else:
            dicoconf['ant2-pos1'] = Locobj_bas.X_coordinate_m
            dicoconf['ant2-pos2'] = Locobj_bas.Y_coordinate_m
            dicoconf['ant2-pos3'] = Locobj_bas.Z_coordinate_m
        dicoconf['ant2-antdelu'] = Antobj_bas.Up_Ecc
        dicoconf['ant2-antdeln'] = Antobj_bas.North_Ecc
        dicoconf['ant2-antdele'] = Antobj_bas.East_Ecc


    if not (bas_srt <= rov_srt <= rov_end <= bas_end):
        print('WARN : not bas_srt <= rov_srt <= rov_end <= bas_end !!!')

    outconffilobj = open(out_conf_fil,'w+')
    for k,v in dicoconf.items():
        lin = k.ljust(20)+'='+str(v)+'\n'
        outconffilobj.write(lin)
    outconffilobj.close()

    # ORBITS
    # SP3
    orblis = multi_downloader_orbs_clks( temp_dir , bas_srt , bas_end , archtype='/',
                                        calc_center = calc_center)
    sp3Z = orblis[0]
    sp3 = genefun.uncompress(sp3Z)

    # BRDC
    statdic = dict()
    statdic['nav'] = ['BRDC']
    nav_srt = dt.datetime(bas_srt.year, bas_srt.month , bas_srt.day )
    orblis = multi_downloader_rinex(statdic,temp_dir , nav_srt , bas_end ,
                                    archtype='/', sorted_mode=0)
    navZ = orblis[0]
    nav = genefun.uncompress(navZ)

    # Command
    com_config  = "-k " + out_conf_fil
    com_interval="-ti " + str(rov_itv)
    com_mode = ""
    #com_mode="-p 4"
    com_resultfile="-o " + out_result_fil
    #com_combinsol="-c"

    exe_path = "rnx2rtkp"
#    exe_path = "/home/pierre/install_softs/RTKLIB/rnx2rtkp"

    bigcomand = ' '.join((exe_path,com_config,com_interval,com_mode,
                          com_resultfile,rnx_rover,rnx_base,nav,sp3))

    print(bigcomand)

    subprocess.call([bigcomand], executable='/bin/bash', shell=True)
    print("RTKLIB RUN FINISHED")

    return None

def track_runner(rnx_rover,rnx_base,working_dir,experience_prefix,
                 XYZbase  = [], XYZrover = [] , outtype = 'XYZ',mode = 'short',
                 interval=None,antmodfile = "~/gg/tables/antmod.dat",
                 calc_center='igs' , forced_sp3_path = ''):

    # paths & files
    working_dir = genefun.create_dir(working_dir)
    temp_dir    = genefun.create_dir(os.path.join(working_dir,'TEMP'))
    out_dir     = genefun.create_dir(os.path.join(working_dir,'OUTPUT'))

    if check_if_compressed_rinex(rnx_rover):
        rnx_rover = crz2rnx(rnx_rover,temp_dir)
    else:
        shutil.copy(rnx_rover,temp_dir)

    if check_if_compressed_rinex(rnx_base):
        rnx_base  = crz2rnx(rnx_base,temp_dir)
    else:
        shutil.copy(rnx_base,temp_dir)

    # RINEX START & END
    rov_srt, rov_end , rov_itv = rinex_start_end(rnx_rover,1)
    bas_srt, bas_end , bas_itv = rinex_start_end(rnx_base,1)

    # RINEX NAMES
    rov_name = os.path.basename(rnx_rover)[0:4]
    bas_name = os.path.basename(rnx_base)[0:4]

    rov_name_uper = rov_name.upper()
    bas_name_uper = bas_name.upper()


    srt_str = rov_srt.strftime("%Y_%j")
    exp_full_name = '_'.join((experience_prefix,rov_name,bas_name,srt_str))

    out_conf_fil   = os.path.join(out_dir,exp_full_name + '.cmd')
    out_result_fil = os.path.join(out_dir,exp_full_name + '.out' )

    print(out_conf_fil)

    confobj = open(out_conf_fil,'w+')


    # Obs Files
    confobj.write(' obs_file' + '\n')
    confobj.write(' '.join((' ',bas_name_uper,os.path.basename(rnx_base) ,'F'))+ '\n')
    confobj.write(' '.join((' ',rov_name_uper,os.path.basename(rnx_rover),'K'))+ '\n')
    confobj.write('\n')

    date = geok.rinexname2dt(os.path.basename(rnx_rover))

    # Nav File
    if forced_sp3_path == '':
        strt_rnd = dt.datetime(*bas_srt.timetuple()[:3])
        end_rnd  = dt.datetime(*bas_end.timetuple()[:3])
        orblis = multi_downloader_orbs_clks( temp_dir , strt_rnd , end_rnd ,
                                            archtype='/',
                                            calc_center = calc_center)
        

        #sp3Z = orblis[0]
        sp3 = [genefun.uncompress(sp3Z) for sp3Z in orblis]
        sp3 = [e  if ".sp3" in e[-5:] else e + ".sp3" for e in sp3]
    else:
        if genefun.is_iterable(forced_sp3_path):
            sp3 = forced_sp3_path
        else:
            sp3 = [forced_sp3_path]
    for sp3_mono in sp3: 
        confobj.write(' '.join((' ','nav_file',sp3_mono ,' sp3'))+ '\n')
    confobj.write('\n')

    # Mode
    confobj.write(' mode ' +  mode + '\n')
    confobj.write('\n')

    # Output
    confobj.write(' pos_root ' + exp_full_name +'.pos' + '\n' )
    confobj.write(' res_root ' + exp_full_name +'.res' + '\n' )
    confobj.write(' sum_file ' + exp_full_name +'.sum' + '\n' )
    confobj.write('\n')

    # Outtype
    confobj.write(' out_type ' + outtype + '\n')
    confobj.write('\n')

    # Interval
    if not interval:
        confobj.write(' interval ' + str(rov_itv) + '\n')
    else:
        confobj.write(' interval ' + str(interval) + '\n')

    confobj.write('\n')

    # Coords
    bool_site_pos = False
    if XYZbase != []:
        if not bool_site_pos:
            confobj.write(' site_pos \n')
            bool_site_pos = True
        XYZbase = [str(e) for e in XYZbase]
        confobj.write(' '.join([' ', bas_name_uper] + XYZbase + ['\n']))

    if XYZrover != []:
        if not bool_site_pos:
            confobj.write(' site_pos \n')
            bool_site_pos = True
        XYZrover = [str(e) for e in XYZrover]
        confobj.write(' '.join([' ', rov_name_uper] + XYZrover + ['\n']))

    if bool_site_pos:
        confobj.write('\n')

    # Offsets
    confobj.write(' ante_off \n')

    Antobj_rov , Recobj_rov , Siteobj_rov , Locobj_rov = \
    gfc.read_rinex_2_dataobjts(rnx_rover)

    confobj.write(' '.join([' ', rov_name_uper ,
                            str(Antobj_rov.North_Ecc) ,
                            str(Antobj_rov.East_Ecc) ,
                            str(Antobj_rov.Up_Ecc) ,
                            Antobj_rov.Antenna_Type , '\n']))

    Antobj_bas , Recobj_bas , Siteobj_bas , Locobj_bas = \
    gfc.read_rinex_2_dataobjts(rnx_base)

    confobj.write(' '.join([' ', bas_name_uper ,
                            str(Antobj_bas.North_Ecc) ,
                            str(Antobj_bas.East_Ecc) ,
                            str(Antobj_bas.Up_Ecc) ,
                            Antobj_bas.Antenna_Type , '\n']))
    confobj.write('\n')

    # Site_stats
    confobj.write(' site_stats \n')
    confobj.write(' ' + bas_name_uper  + " 0.1 0.1 0.1 0 0 0" + '\n')
    confobj.write(' ' + rov_name_uper  + " 20 20 20 0.5 0.5 0.5" + '\n')
    confobj.write('\n')

    # Misc
    confobj.write(" USE_GPTGMF" + '\n')
    confobj.write(" ANTMOD_FILE " + antmodfile + '\n')


    confobj.write(" atm_stats" + '\n')
    confobj.write('  all 0.1 0.00030.00023' + '\n')


    confobj.close()
    #END OF FILE WRITING

    dowstring = ''.join([str(e) for e in geok.dt2gpstime(date)])
    bigcomand = ' '.join(("track -f" ,  out_conf_fil , '-d' , geok.dt2doy(date) ,'-w', dowstring))

    print('INFO : command launched :')
    print(bigcomand)

    # START OF PROCESSING
    os.chdir(temp_dir)
    subprocess.call([bigcomand], executable='/bin/bash', shell=True)

    outfiles = []
    outfiles = outfiles + glob.glob(os.path.join(temp_dir,exp_full_name + '*sum*'))
    outfiles = outfiles + glob.glob(os.path.join(temp_dir,exp_full_name + '*pos*'))
    outfiles = outfiles + glob.glob(os.path.join(temp_dir,exp_full_name + '*cmd*'))

    Antobj_rov , Recobj_rov , Siteobj_rov , Locobj_rov = \
    gfc.read_rinex_2_dataobjts(rnx_rover)

    [shutil.copy(e,out_dir) for e in outfiles]
    [os.remove(e) for e in outfiles]

    print("TRACK RUN FINISHED")
    print('results available in ' , out_dir)
    return None


def run_track(temp_dir,exp_full_name,out_conf_fil,date,rnx_rover):
    dowstring = ''.join([str(e) for e in geok.dt2gpstime(date)])
    bigcomand = ' '.join(("track -f" ,  out_conf_fil , '-d' , geok.dt2doy(date) ,'-w', dowstring))

    print('INFO : command launched :')
    print(bigcomand)

    # START OF PROCESSING
    os.chdir(temp_dir)
    subprocess.call([bigcomand], executable='/bin/bash', shell=True)

    outfiles = []
    outfiles = outfiles + glob.glob(os.path.join(temp_dir,exp_full_name + '*sum*'))
    outfiles = outfiles + glob.glob(os.path.join(temp_dir,exp_full_name + '*pos*'))
    outfiles = outfiles + glob.glob(os.path.join(temp_dir,exp_full_name + '*cmd*'))

    Antobj_rov , Recobj_rov , Siteobj_rov , Locobj_rov = \
    gfc.read_rinex_2_dataobjts(rnx_rover)

    [shutil.copy(e,out_dir) for e in outfiles]
    [os.remove(e) for e in outfiles]

    print("TRACK RUN FINISHED")
    print('results available in ' , out_dir)

    return None



def cluster_GFZ_run(commands_list,
                    bunch_on_off = True,
                    bunch_job_nbr = 10,
                    bunch_wait_time = 600,
                    bj_check_on_off = True,
                    bj_check_mini_nbr = 2,
                    bj_check_wait_time = 120,
                    bj_check_user="auto"):


    history_file_path = None
    wait_sleeping_before_launch=5

    i_bunch = 0

    if bj_check_user == "auto":
        bj_check_user=genefun.get_username()

    log_path = "/home/" + bj_check_user + "/test_tmp.log"
    LOGobj = open(log_path , 'w+')

    print ("****** JOBS THAT WILL BE LAUNCHED ******")
    print ('Number of jobs : ' + str(len(commands_list)))
    print ("****************************************")

    for kommand in commands_list:

        ########## LOG/PRINT command
        print(kommand)
        LOGobj.write(kommand + '\n')

        ########## LOG/PRINT sleep
        info_sleep = "INFO : script is sleeping for " +str(wait_sleeping_before_launch) +"sec (so you can cancel it) "
        print(info_sleep)
        LOGobj.write(info_sleep + '\n')
        time.sleep(wait_sleeping_before_launch)

        ########## LOG/PRINT start
        info_start = "INFO : script starts @ " + str(dt.datetime.now())
        print(info_start)
        LOGobj.write(info_start + '\n')

        ########## RUN command here !!
        p = subprocess.Popen('',executable='/bin/csh', stdin=subprocess.PIPE , stdout=subprocess.PIPE , stderr=subprocess.PIPE)
        stdout,stderr = p.communicate( kommand.encode() )

        if history_file_path:
            with open(history_file_path , "a") as myfile:
                myfile.write(kommand + '\n')

        i_bunch += 1

        ########## Bunch On/Off : check if a bunch of job has been launched
        if bunch_on_off and np.mod(i_bunch , bunch_job_nbr) == 0:
            info_bunch = "INFO : sleeping @ " + str(dt.datetime.now()) + " for " + str(bunch_wait_time) + "s b.c. a bunch of " + str(bunch_job_nbr) + " jobs has been launched"
            print(info_bunch)
            time.sleep(bunch_wait_time)
            LOGobj.write(info_bunch + '\n')

            ########## Bunch On/Off Check : check if the bunch is finished (experimental but on is better)
            if bj_check_on_off:
                print("INFO : BJ Check : All jobs should be finished now, let's see if there is some latecomers")
                bj_check_tigger = False
                bj_command = "perl /dsk/igs2/soft_wrk/psakicki/SOFT_EPOS8_BIN_TOOLS/SCRIPTS/e8_bjobs_local.pl"
            while not bj_check_tigger:
                bj_list = subprocess.check_output(bj_command,shell='/bin/csh')
                bj_pattern_checked =   bj_check_user + ' *' +  bj_check_user

                bj_list_checked = [re.search(bj_pattern_checked,l) for l in bj_list]
                bj_list_checked_sum = np.sum(bj_list_checked)

                if bj_list_checked_sum > bj_check_mini_nbr:
                    print("INFO : sleeping @ " + str(dt.datetime.now()) + " for " + str(bj_check_wait_time) + "s b.c." + str(bj_list_checked_sum) + "job(s) match pattern " + bj_pattern_checked)
                    print(bj_list)
                    time.sleep(bj_check_wait_time)

                else:
                    bj_check_tigger = True
                    print("INFO : let's continue, no job matchs the pattern " + bj_pattern_checked)


#
#  __  __ _____ _____           _____
# |  \/  |_   _|  __ \   /\    / ____|
# | \  / | | | | |  | | /  \  | (___
# | |\/| | | | | |  | |/ /\ \  \___ \
# | |  | |_| |_| |__| / ____ \ ____) |
# |_|  |_|_____|_____/_/    \_\_____/
#
#

def midas_run(tenu_file_path,work_dir="",
              path_midas_soft="",
              step_file_path="",
              with_plot=True,
              keep_plot_open = True):
    ### Prepare paths of facultative arguments
    if not work_dir:
        work_dir=os.path.dirname(tenu_file_path)
    if not path_midas_soft:
        path_midas_soft = "midas.e"

    ### SECURITY : rm all previous MIDAS files
    MIDAS_files_rm = ["MIDAS.STEPIN", "MIDAS.STEPOUT", "MIDAS.ERR",
                      "MIDAS.TENV", "MIDAS.RENV", "MIDAS.TENU",
                      "MIDAS.RENU", "MIDAS.VEL","fort.91"]

    for fil in MIDAS_files_rm:
        fil_full_path = os.path.join(work_dir,fil)
        if os.path.isfile(fil_full_path):
            print("INFO : old",fil,"removed")
            os.remove(fil_full_path)

    ### Find the Station Name
    DF = pandas.read_table(tenu_file_path,comment='*',header=-1,delim_whitespace = True)
    stat = list(set(DF[0]))[0]

    ### Prepare the tmp input
    ### The file is called MIDAS.TENU in any case
    os.chdir(work_dir)
    work_file_path_tenu = work_dir + "/MIDAS.TENU"
    shutil.copy(tenu_file_path , work_file_path_tenu)

    ### Prepare the tmp STEP FILE input
    step_OK = False
    if step_file_path and not genefun.empty_file_check(step_file_path):
        work_file_path_step = work_dir + "/MIDAS.STEP"
        shutil.copy(step_file_path , work_file_path_step)
        step_OK = True
    else:
        work_file_path_step = ""

    ### Prepare the name of the outputed file
    work_file_path_vel  = os.path.join(work_dir,'MIDAS.VEL')

    ### LETS RUN !
    command = path_midas_soft
    print('LAUNCHING : ',  command, "for stat." , stat)
    p = subprocess.Popen('',executable='/bin/bash', stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE , stderr=subprocess.PIPE)
    stdout,stderr = p.communicate( command.encode() )
    # logs files
    vel_file = os.path.join(work_dir ,  "MIDAS.VEL")
    vel_file_obj = open(vel_file, "w")
    vel_file_obj.write(stdout.decode("utf-8") )
    vel_file_obj.close()

    exec_OK = False

    if not genefun.empty_file_check(work_file_path_vel):
        exec_OK = True
        print("INFO : MIDAS.VEL exists for",stat,", should be OK :)")

    if not genefun.empty_file_check(os.path.join(work_dir,'MIDAS.ERR')):
        print("ERR : MIDAS.ERR is not empty for",stat,", must be checked :(")

    if stderr:
        print("err.log is not empty, must be checked !")
        print(stderr.decode("utf-8"))
        print('')
        err_file = os.path.join(work_dir , stat + ".err.log")
        err_file_obj = open(err_file, "w")
        err_file_obj.write(stderr.decode("utf-8"))
        err_file_obj.close()

    #### PLOT
    if exec_OK and with_plot:
        fig = midas_plot(work_file_path_tenu,
                         work_file_path_vel,
                         work_file_path_step)
        work_file_path_plot = os.path.join(work_dir, "plot." + stat.upper())
        for ext in (".pdf",".png"):
            plt.savefig(work_file_path_plot + ext)
        if not keep_plot_open:
            plt.close(fig)


    #### RENAME ALL MIDAS FILES WITH A STATION PREFIX
    for fname in glob.glob('*MIDAS*'):
        os.rename(fname, fname.lower() + "_" + stat.upper())

    return None



def midas_vel_files_2_pandas_DF(vel_files_in):
    """
    Convert MIDAS Velocity files to a Pandas DataFrame

    Parameters
    ----------
    vel_files_in : str or list of str
        if list of str, will consider directly the files inside the list
        if str, can be the path of a single file, or a generic (wilcard) path
        to consider several files
    Returns
    -------
    DF : Pandas DataFrame
    """

    if genefun.is_iterable(vel_files_in):
        L_vel_files = vel_files_in
    else:
        L_vel_files = glob.glob(vel_files_in)


    if len(L_vel_files) > 1:
        work_dir = os.path.dirname(L_vel_files[0])
        vel_file_opera = genefun.cat(work_dir + "/MERGED_VEL" , *L_vel_files)
    else:
        vel_file_opera = L_vel_files[0]

    DF = pandas.read_table(vel_file_opera,comment='*',header=-1,delim_whitespace = True)

    DF.columns = ["Station",
              "soft",
              "epoch_first",
              "epoch_last",
              "duration",
              "nb_epoch_all",
              "nb_epoch_good",
              "nb_pairs",
              "V_East",
              "V_North",
              "V_Up",
              "sV_East",
              "sV_North",
              "sV_Up",
              "offset_e_1st_epoch",
              "offset_n_1st_epoch",
              "offset_u_1st_epoch",
              "outlier_ratio_e",
              "outlier_ratio_n",
              "outlier_ratio_u",
              "std_velo_pair_e",
              "std_velo_pair_n",
              "std_velo_pair_u",
              "nb_steps"]

    return DF


def midas_plot(path_tenu,path_vel="",path_step=""):
    """
    based on a plot for TimeSeriePoint
    """
    import geoclass as gcls

    DF = pandas.read_table(path_tenu,header=-1,delim_whitespace = True)
    stat = list(set(DF[0]))[0]

    if path_vel:
        DFvel = midas_vel_files_2_pandas_DF(path_vel)
        plot_vel = True
    else:
        plot_vel = False


    if path_step:
        DFstep  = pandas.read_table(path_step,header=-1,delim_whitespace = True)
        Discont = geok.year_decimal2dt(list(DFstep[1]))


    T = DF[1]
    Tdt = geok.year_decimal2dt(T)
    E = DF[2]
    N = DF[3]
    U = DF[4]

    TS = gcls.TimeSeriePoint()

    TS.from_list(Tdt,E,N,U,"ENU")
    TS.meta_set(path_tenu,stat,stat)

    if path_step:
        TS.set_discont(Discont)
        print("AAAAAAAA",Discont)

    #### PLOT ########
    TSplot_fig = TS.plot(fig=plt.figure())
    if path_step:
        TS.plot_discont(fig=TSplot_fig)
    ##################
    TSplot_axe = TSplot_fig.axes

    a_keys_list = ["V_East",
             "V_North",
             "V_Up"]

    b_keys_list = ["offset_e_1st_epoch",
             "offset_n_1st_epoch",
             "offset_u_1st_epoch"]

    sigma_key_list = ["sV_East",
             "sV_North",
             "sV_Up"]

    if plot_vel:
        for i , (a_key , b_key , sigma_key) in enumerate(zip(a_keys_list ,
                b_keys_list,
                sigma_key_list)):
            ax = TSplot_axe[i]


            a = DFvel[a_key][0]
            b = DFvel[b_key][0]
            sigma = DFvel[sigma_key][0]

            _ , Vlin = geok.linear_reg_getvalue(T-T[0],a,b)

            ax.plot(Tdt,Vlin,c="xkcd:dark orange")
            text_vel = str(np.round(a*1000,4)) + " +/- " + str(np.round(sigma*1000,4)) + " mm/yr"

            ax.text(0.8, 0.1, text_vel, horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)

    return TSplot_fig
