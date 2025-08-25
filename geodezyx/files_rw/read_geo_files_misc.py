#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.files_rw contains reading functions to 
import files containing geodetic observations/products.

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
import linecache
#### Import the logger
import logging
import os

import numpy as np
import pandas as pd

#### geodeZYX modules
from geodezyx import utils

log = logging.getLogger('geodezyx')


##########  END IMPORT  ##########


 #  ______                _   _                _____                                         _ 
 # |  ____|              | | (_)              / ____|                                       | |
 # | |__ _   _ _ __   ___| |_ _  ___  _ __   | |  __ _ __ __ ___   _____ _   _  __ _ _ __ __| |
 # |  __| | | | '_ \ / __| __| |/ _ \| '_ \  | | |_ | '__/ _` \ \ / / _ \ | | |/ _` | '__/ _` |
 # | |  | |_| | | | | (__| |_| | (_) | | | | | |__| | | | (_| |\ v /  __/ |_| | (_| | | | (_| |
 # |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|  \_____|_|  \__,_| \_/ \___|\__, |\__,_|_|  \__,_|
 #                                                                        __/ |                
 #                                                                       |___/       

def read_erp_bad(path,return_array=False):
    """
    This function is discontinued, use read_erp1 instead
    """
    F = open(path)
    L = []
    for l in F:
        if "MJD" in l:
            keys = l.split()
        try:
            float(l[0])
        except:
            continue
        l=utils.str_2_float_line(l.strip())
        L.append(l)

    M = np.vstack(L)

    if return_array:
        return M
    else:
        return pd.DataFrame(M)


############################# READ CLK 30
#def read_clk1(file_path_in, returns_pandas = True, interval=None):
#    """
#    General description
#
#    Parameters
#    ----------
#    file_path_in :  str
#        Path of the file in the local machine.
#
#    returns_pandas :  bool
#        Define if pandas function will be used to put the data on tables
#
#    interval :  int
#        Defines which interval should be used in the data tables. The interval is always in minutes unit.
#
#    Returns
#    -------
#    out1 : float or int oandas table
#        Returns a panda table format with the data extracted from the file.
#
#    out2 :  list
#       In case of the pandas function are not in use, the return will be a list with the data extract from the file.
#
#
#    """
#
#    file = open(file_path_in,'r')
#    fil = file.readlines()
#    file.close()
#
#
#    types_m, name_m, year_m, month_m, day_m, hr_m, minute_m, seconds_m, epoch_m, offset_m, rms_m = [],[],[],[],[],[],[],[],[],[],[]
#
#    Clk_read  = []
#
#        #####IGNORES HEADER
#
#    le = 0
#    i = 0
#    count = 0
#    for le in range(len(fil)):
#     linhaatual = linecache.getline(file_path_in, le)
#     header_end = (linhaatual[60:73])
#     count +=1
#     if header_end =="END OF HEADER":
#        i = count
##############################################################
#    while i <= len(fil):
#          linhaatual = linecache.getline(file_path_in, i)
#          types = (linhaatual[0:2]) # Column refers to AS or AR
#          name = (linhaatual[3:7])# Name of the station or satellite
#          year = int((linhaatual[8:12]))
#          month = int((linhaatual[13:15]))
#          day = int((linhaatual[16:18]))
#          hr= int((linhaatual[19:22]))# hour
#          minute = int((linhaatual[22:25]))
#          seconds = float((linhaatual[25:33]))
#          offset = float(((linhaatual[40:59])))# Clock offset
#          ##### check if there is a value for thr rms
#          if (linhaatual[65:70]) == '     ' or len(linhaatual[65:70])==0:
#              rms = 0
#          else:
#              rms = float(linhaatual[61:79])
#          ############
#          epoch = dt.datetime(year,month,day,hr,minute,int(seconds),int(((int(seconds)-seconds)*10**-6))) # epoch as date.time function
#
#        ############ Data in a defined interval
#          if interval:
#              if (float(minute)%interval) == 0 and float(seconds) == 0:
#
#                  if returns_pandas:
#                            Clk_dados = [types,name,year,month,day,hr,minute,seconds,epoch,offset,rms]
#                            Clk_read.append(Clk_dados)
#                  else:
#                            types_m.append(types)
#                            name_m.append(name)
#                            year_m.append(year)
#                            month_m.append(month)
#                            day_m.append(day)
#                            hr_m.append(hr)
#                            minute_m.append(minute)
#                            seconds_m.append(seconds)
#                            epoch_m.append(epoch)
#                            offset_m.append(offset)
#                            rms_m.append(rms)
#         #########################################################
#        ########################################################### standart file epochs
#          else:
#              if returns_pandas:
#                            Clk_dados = [types,name,year,month,day,hr,minute,seconds,epoch,offset,rms]
#                            Clk_read.append(Clk_dados)
#              else:
#                            types_m.append(types)
#                            name_m.append(name)
#                            year_m.append(year)
#                            month_m.append(month)
#                            day_m.append(day)
#                            hr_m.append(hr)
#                            minute_m.append(minute)
#                            seconds_m.append(seconds)
#                            epoch_m.append(epoch)
#                            offset_m.append(offset)
#                            rms_m.append(rms)
#
#
#          i +=1
#
#  ################################################################################################
#   ############################################ put on pandas table format
#    if returns_pandas:
#     Clk_readed = pd.DataFrame(Clk_read, columns=['type','name','year','month','day','hr','minutes','seconds','epoch','offset','rms'])
#     Clk_readed.path = file_path_in
#
#     return Clk_readed
#
#          ###############################
#    else:
#           print(" ...")
#           return  types_m, name_m, year_m, month_m, day_m, hr_m, minute_m, seconds_m, offset_m, rms_m
#



############################# END FUNCTION READ CLK 30




################################################READ ERP GUS
# def read_erp1(caminho_arq,ac):

    # le = open(caminho_arq, 'r')
    # letudo = le.readlines()
    # le.close()
    # tamanho = len(letudo) #usado para saber quantas linhas tem o arquivo



    # para = tamanho #onde o arquivo para de ser lido

    # numeros = ['0','1','2','3','4','5','6','7','8','9']
    # le = 0
    # numlin = 0 #numero de linhas da matriz de epocas
    # numcol = 16 #numero de colunas que a matriz final deve ter


    # while le <= para:
     # linhaatual = linecache.getline(caminho_arq, le)
     # if linhaatual[0:1] in numeros:
         # numlin +=1
     # le +=1

    # dt = np.dtype(object)
    # ERP = np.empty((numlin,numcol), dtype=dt)
    # n = 0
    # j = 0
    # l = 0
    # g = 0
# #####################################
    # if ac == 'wum':
        # while l <=para:
         # linhaatual = linecache.getline(caminho_arq, l)
         # if linhaatual[0:1] in numeros:
             # ERP[j,0] = float(linhaatual[0:8])
             # ERP[j,1] = round((float(linhaatual[11:17])*(10**-6)),6)
             # ERP[j,2] = float(linhaatual[18:25])*(10**-6)
             # ERP[j,4] = float(linhaatual[26:34])*(10**-7)
             # ERP[j,3] = float(linhaatual[35:41])*(10**-7)
             # ERP[j,5] = float(linhaatual[44:47])*(10**-6)
             # ERP[j,6] = float(linhaatual[50:53])*(10**-6)
             # ERP[j,7] = float(linhaatual[57:61])*(10**-7)
             # ERP[j,8] = float(linhaatual[63:69])*(10**-7)
             # ERP[j,9] = float(linhaatual[69:73])
             # ERP[j,10] = float(linhaatual[74:76])
             # ERP[j,11] = float(linhaatual[76:80])
             # ERP[j,12] = float(linhaatual[81:86])*(10**-6)
             # ERP[j,13] = float(linhaatual[88:93])*(10**-6)
             # ERP[j,14] = float(linhaatual[96:101])*(10**-6)
             # ERP[j,15] = float(linhaatual[102:107])*(10**-6)

             # j +=1
         # l +=1

        # Erp_end = pd.DataFrame(ERP, columns=['MJD','X-p', 'Y-p', 'UT1UTC','LOD','S-X','S-Y','S-UT','S-LD','NR', 'NF', 'NT',
                                             # 'X-RT','Y-RT','S-XR','S-YR'])
    # #    Erp_end.set_index('MJD',inplace=True)
        # return Erp_end
# #######################################
    # if ac == 'COD':
        # while n <=para:
         # linhaatual = linecache.getline(caminho_arq, n)
         # if linhaatual[0:1] in numeros:
             # ERP[j,0] = float(linhaatual[0:8])
             # ERP[j,1] = round((float(linhaatual[12:17])*(10**-6)),6)
             # ERP[j,2] = float(linhaatual[20:27])*(10**-6)
             # ERP[j,4] = float(linhaatual[38:41])*(10**-7)
             # ERP[j,3] = float(linhaatual[29:35])*(10**-7)
             # ERP[j,5] = float(linhaatual[47:49])*(10**-6)
             # ERP[j,6] = float(linhaatual[53:55])*(10**-6)
             # ERP[j,7] = float(linhaatual[59:61])*(10**-7)
             # ERP[j,8] = float(linhaatual[65:66])*(10**-7)
             # ERP[j,9] = float(linhaatual[67:70])
             # ERP[j,10] = float(linhaatual[71:74])
             # ERP[j,11] = float(linhaatual[75:77])
             # ERP[j,12] = float(linhaatual[81:85])*(10**-6)
             # ERP[j,13] = float(linhaatual[88:92])*(10**-6)
             # ERP[j,14] = float(linhaatual[94:98])*(10**-6)
             # ERP[j,15] = float(linhaatual[100:104])*(10**-6)

             # j +=1
         # n +=1

        # Erp_end = pd.DataFrame(ERP, columns=['MJD','X-p', 'Y-p', 'UT1UTC','LOD','S-X','S-Y','S-UT','S-LD','NR', 'NF', 'NT',
                                             # 'X-RT','Y-RT','S-XR','S-YR'])
    # #    Erp_end.set_index('MJD',inplace=True)
        # return Erp_end
# ################################################
    # if ac == 'gbm':
        # while g <=para:
         # linhaatual = linecache.getline(caminho_arq, g)
         # if linhaatual[0:1] in numeros:
             # ERP[j,0] = float(linhaatual[0:9])
             # ERP[j,1] = round((float(linhaatual[12:18])*(10**-6)),6)
             # ERP[j,2] = float(linhaatual[19:25])*(10**-6)
             # ERP[j,4] = float(linhaatual[26:37])*(10**-7)  ###UT1 - TAI!!!!!!!!!!!!
             # ERP[j,3] = float(linhaatual[39:44])*(10**-7)
             # ERP[j,5] = float(linhaatual[46:49])*(10**-6)
             # ERP[j,6] = float(linhaatual[50:54])*(10**-6)
             # ERP[j,7] = float(linhaatual[56:59])*(10**-7)
             # ERP[j,8] = float(linhaatual[60:64])*(10**-7)
             # ERP[j,9] = float(linhaatual[65:70])
             # ERP[j,10] = float(linhaatual[71:72])
             # ERP[j,11] = float(linhaatual[73:76])
             # ERP[j,12] = float(linhaatual[77:82])*(10**-6)
             # ERP[j,13] = float(linhaatual[83:88])*(10**-6)
             # ERP[j,14] = float(linhaatual[92:96])*(10**-6)
             # ERP[j,15] = float(linhaatual[100:104])*(10**-6)

             # j +=1
         # g +=1

        # Erp_end = pd.DataFrame(ERP, columns=['MJD','X-p', 'Y-p', 'UT1UTC','LOD','S-X','S-Y','S-UT','S-LD','NR', 'NF', 'NT',
                                             # 'X-RT','Y-RT','S-XR','S-YR'])
    # #    Erp_end.set_index('MJD',inplace=True)
        # return Erp_end
        # #out = '/home/mansur/Documents/MGEX/Saida/outERP_'+caminho_arq[28:-4]+'.txt'
        # #np.savetxt(out, ERP, fmt="%02s %04s %05s %05s %05s %05s %05s %10s %10s %10s %10s %10s %10s %10s %10s %10s")
################################################### END READ ERP GUS
##########################################################################################################################

def read_clk_old(file_path_in, returns_pandas = True, interval=None,old_naming=True):
    """
    Read an IGS clk file
    Slow and complex, use read_clk instead
    
    Parameters
    ----------
    file_path_in :  str
        Path of the file in the local machine.

    returns_pandas :  bool
        Define if pandas function will be used to put the data on tables

    interval :  int
        Defines which interval should be used in the data tables. The interval is always in minutes unit.

    Returns
    -------
    out1 : float or int oandas table
        Returns a panda table format with the data extracted from the file.

    out2 :  list
       In case of the pandas function are not in use, the return will be a list with the data extract from the file.


    """

    file = open(file_path_in,'r')
    fil = file.readlines()
    file.close()


    types_m, name_m, year_m, month_m, day_m, hr_m, minute_m, seconds_m, epoch_m, offset_m, rms_m = [],[],[],[],[],[],[],[],[],[],[]

    Clk_read  = []

        #####IGNORES HEADER

    le = 0
    i = 0
    count = 0
    for le in range(len(fil)):
     linhaatual = linecache.getline(file_path_in, le)
     header_end = (linhaatual[60:73])
     count +=1
     if header_end =="END OF HEADER":
        i = count
#############################################################
    while i <= len(fil):
          linhaatual = linecache.getline(file_path_in, i)
          types = (linhaatual[0:2]) # Column refers to AS or AR
          name = (linhaatual[3:7])# Name of the station or satellite
          year = int((linhaatual[8:12]))
          month = int((linhaatual[13:15]))
          day = int((linhaatual[16:18]))
          hr= int((linhaatual[19:22]))# hour
          minute = int((linhaatual[22:25]))
          seconds = float((linhaatual[25:33]))
          offset = float(((linhaatual[40:59])))# Clock offset
          ##### check if there is a value for thr rms
          if (linhaatual[65:70]) == '     ' or len(linhaatual[65:70])==0:
              rms = 0
          else:
              rms = float(linhaatual[61:79])
          ############
          epoch = dt.datetime(year,month,day,hr,minute,int(seconds),int(((int(seconds)-seconds)*10**-6))) # epoch as date.time function

        ############ Data in a defined interval
          if interval:
              if (float(minute)%interval) == 0 and float(seconds) == 0:

                  if returns_pandas:
                            Clk_dados = [types,name,year,month,day,hr,minute,seconds,epoch,offset,rms]
                            Clk_read.append(Clk_dados)
                  else:
                            types_m.append(types)
                            name_m.append(name)
                            year_m.append(year)
                            month_m.append(month)
                            day_m.append(day)
                            hr_m.append(hr)
                            minute_m.append(minute)
                            seconds_m.append(seconds)
                            epoch_m.append(epoch)
                            offset_m.append(offset)
                            rms_m.append(rms)
         #########################################################
        ########################################################### standart file epochs
          else:
              if returns_pandas:
                            Clk_dados = [types,name,year,month,day,hr,minute,seconds,epoch,offset,rms]
                            Clk_read.append(Clk_dados)
              else:
                            types_m.append(types)
                            name_m.append(name)
                            year_m.append(year)
                            month_m.append(month)
                            day_m.append(day)
                            hr_m.append(hr)
                            minute_m.append(minute)
                            seconds_m.append(seconds)
                            epoch_m.append(epoch)
                            offset_m.append(offset)
                            rms_m.append(rms)


          i +=1

  ################################################################################################
   ############################################ put on pandas table format
    if returns_pandas:
        if old_naming:
            name_list = ['offset','rms']
        else:
            name_list = ['bias','sigma']
            
        Clk_readed = pd.DataFrame(Clk_read, columns= ['type','name','year','month','day','h','minutes','seconds','epoch'] + name_list)
        Clk_readed['epoch'] =  pd.to_datetime(Clk_readed[[ 'year' ,'month' ,'day','h','minutes','seconds']])
        Clk_readed.path = file_path_in
    
        return Clk_readed

          ###############################
    else:
           return  types_m, name_m, year_m, month_m, day_m, hr_m, minute_m, seconds_m, offset_m, rms_m




############################# END FUNCTION READ CLK 30





# def read_erp(caminho_arq,ac):
    # """
    # General description
    # Units: ('MJD','X-p (arcsec)', 'Y-p (arcsec)', 'UT1UTC (E-7S)','LOD (E-7S/D)','S-X (E-6" arcsec)','S-Y (E-6" arcsec)',
    # 'S-UT (E-7S)','S-LD (E-7S/D)','NR (E-6" arcsec)', 'NF (E-6" arcsec)', 'NT (E-6" arcsec)',
    # 'X-RT (arcsec/D)','Y-RT (arcsec/D)','S-XR (E-6" arcsec/D)','S-YR (E-6" arcsec/D)', 'C-XY', 'C-XT',
    # 'C-YT', 'DPSI', 'DEPS','S-DP','S-DE')

    # Parameters
    # ----------
    # file_path_in :  str
        # Path of the file in the local machine.

    # which AC :  str
        # The analisys center that will be used


    # Returns
    # -------
    # out1 :  pandas table
        # Returns a panda table format with the data extracted from the file.

    # Obs: this function need to be improved to read also the values of UTC and etc. So far reads only the Pole rates. 

    # """

    # #### FIND DELIVERY DATE
    # name = os.path.basename(caminho_arq)
    
    # if len(name) == 12:
        # dt_delivery = conv.sp3name2dt(caminho_arq)
    # elif len(name) == 38:
        # dt_delivery = conv.sp3name_v3_2dt(caminho_arq)
    # else:
        # dt_delivery = conv.posix2dt(0)


    # le = open(caminho_arq, 'r')
    # letudo = le.readlines()
    # le.close()
    # tamanho = len(letudo) 

    

    # numeros = ['0','1','2','3','4','5','6','7','8','9']


    # ERP=[]


    # if caminho_arq[-3:] in ('snx','ssc'):
        # file = open(caminho_arq)
        # Lines =  file.readlines()
        # XPO_stk  = []
        # XPO_std_stk = []
        # YPO_stk  = []
        # YPO_std_stk = []
        # LOD_stk  = []
        # LOD_std_stk = []
        # MJD_stk = []
        # marker = False

        # for i in range(len(Lines)):

            # if len(Lines[i].strip()) == 0:
                # continue
            # else:

                # if Lines[i].split()[0] == '+SOLUTION/ESTIMATE':
                    # marker = True

                # if Lines[i].split()[0] == '-SOLUTION/ESTIMATE':
                    # marker = False

                # if utils.contains_word(Lines[i],'XPO') and marker:
                    # Doy = (Lines[i][30:33])
                    # Year = (Lines[i][27:29])
                    # Pref_year = '20'
                    # Year = int(Pref_year+Year)
                    # Date = conv.doy2dt(Year,Doy)
                    # XPO = float(Lines[i][47:68])*(10**-3)
                    # XPO_std = float(Lines[i][69:80])*(10**-3)
                    # XPO_stk.append(XPO)
                    # XPO_std_stk.append(XPO_std)
                    # MJD_stk.append(conv.jd_to_mjd(conv.date_to_jd(Date.year,Date.month,Date.day)))
                    
                # if utils.contains_word(Lines[i],'YPO') and marker:
                    # Doy = (Lines[i][30:33])
                    # Year = str(Lines[i][27:29])
                    # Pref_year = '20'
                    # Year = int(Pref_year+Year)
                    # Date = conv.doy2dt(Year,Doy)
                    # YPO = float(Lines[i][47:68])*(10**-3)
                    # YPO_std = float(Lines[i][69:80])*(10**-3)
                    # YPO_stk.append(YPO)
                    # YPO_std_stk.append(YPO_std)
                    # MJD_stk.append(conv.jd_to_mjd(conv.date_to_jd(Date.year,Date.month,Date.day)))
                    
                # if utils.contains_word(Lines[i],'LOD') and marker:
                    # Doy = (Lines[i][30:33])
                    # Year = str(Lines[i][27:29])
                    # Pref_year = '20'
                    # Year = int(Pref_year+Year)
                    # Date = conv.doy2dt(Year,Doy)
                    # LOD = float(Lines[i][47:68])*(10**+4)
                    # LOD_std = float(Lines[i][69:80])*(10**+4)
                    # LOD_stk.append(LOD)
                    # LOD_std_stk.append(LOD_std)
                    # MJD_stk.append(conv.jd_to_mjd(conv.date_to_jd(Date.year,Date.month,Date.day)))
                    
        # MJD = list(set(MJD_stk))
        # if len(LOD_stk) == 0:
                # LOD_stk = ['0']*len(MJD)
                # LOD_std_stk = ['0']*len(MJD)

        # for i in range(len(MJD)):

            # ERP_data = [ac, MJD[i], XPO_stk[i], YPO_stk[i], 0, LOD_stk[i], XPO_std_stk[i], YPO_std_stk[i],
                     # 0, LOD_std_stk[i], 0, 0, 0, 0, 0, 0, 0, dt_delivery]
            # print(ERP_stk)
            # ERP.append(ERP_data)



    # if ac in ('COD','cod','com', 'cof', 'grg', 'mit', 'sio'):
        # for i in range(tamanho+1):
            # linhaatual = linecache.getline(caminho_arq, i)
            # if linhaatual[0:1] in numeros:
                # ERP_data = linhaatual.split()
                # for j in range(len(ERP_data)):
                    # ERP_data[j] = float(ERP_data[j])
                # ERP_data.insert(0,ac)
                # ERP_data[2] =  ERP_data[2]*(10**-6)
                # ERP_data[3] =  ERP_data[3]*(10**-6)
                # ERP_data[13] = ERP_data[13]*(10**-6)
                # ERP_data[14] = ERP_data[14]*(10**-6)
                # del ERP_data[17:]
                # ERP_data.append(dt_delivery)

                # ERP.append(ERP_data)



    # if ac in ('wum','grg','esa', 'mit', 'ngs', 'sio'):
        # for i in range(tamanho+1):
            # linhaatual = linecache.getline(caminho_arq, i)
            # if linhaatual[0:1] in numeros:
                # ERP_data = linhaatual.split()
                # for j in range(len(ERP_data)):
                    # ERP_data[j] = float(ERP_data[j])
                # ERP_data.insert(0,ac)
                # ERP_data[2] =  ERP_data[2]*(10**-6)
                # ERP_data[3] =  ERP_data[3]*(10**-6)
                # ERP_data[13] = ERP_data[13]*(10**-6)
                # ERP_data[14] = ERP_data[14]*(10**-6)
                # ERP_data.append(dt_delivery)

                # ERP.append(ERP_data)

# #
    # if ac in ('gbm', 'gfz'):
        # for i in range(tamanho+1):
            # linhaatual = linecache.getline(caminho_arq, i)
            # if linhaatual[0:1] in numeros:
                # ERP_data = linhaatual.split()
                # for j in range(len(ERP_data)):
                    # ERP_data[j] = float(ERP_data[j])
                # ERP_data.insert(0,ac)
                # ERP_data[2] =  ERP_data[2]*(10**-6)
                # ERP_data[3] =  ERP_data[3]*(10**-6)
                # ERP_data[13] = ERP_data[13]*(10**-6)
                # ERP_data[14] = ERP_data[14]*(10**-6)
                # ERP_data.append(dt_delivery)

                # ERP.append(ERP_data)


    # header = []
    # if ac in ('emr'):
        # for i in range(tamanho+1):
            # linhaatual = linecache.getline(caminho_arq, i)
            # if linhaatual == 'EOP  SOLUTION':
                # del ERP_data[:]
                # header = ['EOP  SOLUTION']
            # if linhaatual[0:1] in numeros and 'EOP  SOLUTION' in header:
                # ERP_data = linhaatual.split()
                # for j in range(len(ERP_data)):
                    # ERP_data[j] = float(ERP_data[j])
                # ERP_data.insert(0,ac)
                # ERP_data[2] =  ERP_data[2]*(10**-6)
                # ERP_data[3] =  ERP_data[3]*(10**-6)
                # ERP_data[13] = ERP_data[13]*(10**-6)
                # ERP_data[14] = ERP_data[14]*(10**-6)
                # del ERP_data[17:]
                # ERP_data.append(dt_delivery)

                # ERP.append(ERP_data)

    # Erp_end = pd.DataFrame(ERP, columns=['AC','MJD','X-p', 'Y-p', 'UT1UTC(UT1 -TAI)','LOD','S-X','S-Y','S-UT','S-LD',
                                         # 'NR', 'NF', 'NT',
                                         # 'X-RT','Y-RT','S-XR','S-YR',
                                         # 'Delivered_date'])
    # return Erp_end


def list_files(dire,file = None):
    """
    written for MGEX combi by GM
    should be discontinued
    """
    List_arq = os.listdir(dire)

    path = []

    for l in range(len(List_arq)):
        atual = List_arq[l]
        if file == 'sp3':
            if atual[-3:] == 'SP3' or atual[-3:] == 'sp3':
                path.append(os.path.join(dire,List_arq[l]))

        if file == 'clk':
            if atual[-3:] == 'CLK' or atual[-3:] == 'clk':
                path.append(os.path.join(dire,List_arq[l]))

        if file == 'erp':
            if atual[-3:] == 'ERP' or atual[-3:] == 'erp':
                path.append(os.path.join(dire,List_arq[l]))

        if file == 'eph':
            if atual[-3:] == 'EPH' or atual[-3:] == 'eph':
                path.append(os.path.join(dire,List_arq[l]))


        if file == 'snx':
            if atual[-3:] == 'SNX' or atual[-3:] == 'snx':
                path.append(os.path.join(dire,List_arq[l]))


    return path

 
def AC_equiv_vals(AC1,AC2):
    """
    Find common values between 2 Orbit DataFrame (from SP3)

    Note
    ----
    Redundant with reffram.orb_df_common_epoch_finder
    this function should not be used anymore
    
    Parameters
    ----------
    AC1 : DataFrame
        Orbit DataFrame 1.
    AC2 : DataFrame
        Orbit DataFrame 2.

    Returns
    -------
    AC1_ok : DataFrame
        Cleaned Orbit DataFrame 1.
    AC2_ok : DataFrame
        Cleaned Orbit DataFrame 2.
    ACmerged : DataFrame
        Merged orbit DataFrame.
    """
    
    log.warning("This function is redundant with reffram.orb_df_common_epoch_finder, use the latter one")
    
    ### 1) Merge the 2 DF to find common lines
    ACmerged = pd.merge(AC1 , AC2 , how='inner', on=['epoch', 'sat'])

    ### 2) Extract merged epoch & sv
    common_epoch = ACmerged["epoch"]
    common_sat   = ACmerged["sat"]
    ### 3) Create a boolean line based on common epoch / sv
    common_sat_epoch_AC1 = (AC1["epoch"].isin(common_epoch)) & (AC1["sat"].isin(common_sat))
    common_sat_epoch_AC2 = (AC2["epoch"].isin(common_epoch)) & (AC2["sat"].isin(common_sat))
    ### 4) Get epoch and sv in the combined sol which correspond to the SP3
    AC1new = AC1[common_sat_epoch_AC1].copy()
    AC2new = AC2[common_sat_epoch_AC2].copy()
    ### 5) A sort to compare the same things
    AC1new.sort_values(by=['sat','sv'],inplace=True)
    AC2new.sort_values(by=['sat','sv'],inplace=True)

    ### Check for > 99999 vals
    AC1_bad_bool_9  = (AC1new["x"] > 9999) & (AC1new["y"] > 9999) & (AC1new["z"] > 9999)
    AC1_bad_bool_9  = np.logical_not(np.array(AC1_bad_bool_9))

    AC2_bad_bool_9  = (AC2new["x"] > 9999) & (AC2new["y"] > 9999) & (AC2new["z"] > 9999)
    AC2_bad_bool_9  = np.logical_not(np.array(AC2_bad_bool_9))

    AC12_bad_bool_9 = np.array(np.logical_and(AC1_bad_bool_9 , AC2_bad_bool_9))

    ### Check for NaN vals
    AC1_bad_bool_nan  = np.isnan(AC1new["x"]) & np.isnan(AC1new["y"]) & np.isnan(AC1new["z"])
    AC1_bad_bool_nan  = np.logical_not(np.array(AC1_bad_bool_nan))

    AC2_bad_bool_nan  = np.isnan(AC2new["x"]) & np.isnan(AC2new["y"]) & np.isnan(AC2new["z"])
    AC2_bad_bool_nan  = np.logical_not(np.array(AC2_bad_bool_nan))

    AC12_bad_bool_nan = np.array(np.logical_and(AC1_bad_bool_nan , AC2_bad_bool_nan))

    AC12_bad_bool = np.logical_and(AC12_bad_bool_9, AC12_bad_bool_nan)

    AC1_ok = AC1new[AC12_bad_bool]
    AC2_ok = AC2new[AC12_bad_bool]

    return AC1_ok , AC2_ok , ACmerged
