#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 10:24:42 2021

@author: psakic
"""

########## BEGIN IMPORT ##########
#### External modules
import linecache
#### Import the logger
import logging
import os

import numpy as np
import pandas as pd

#### geodeZYX modules
from geodezyx import conv
from geodezyx import utils

log = logging.getLogger(__name__)


def read_erp(file_path_in,ac=None):
    """
    
    Read IGS Analysis Center ERP files

    Parameters
    ----------
    file_path_in :  str
        Path of the file in the local machine.

    ac :  str
        The analysis center that will be used. 
        If not precised, will be the first 3 letters of the input name


    Returns
    -------
    out1 :  pandas table
        Returns a panda table format with the data extracted from the file.
        

    Note
    ----
    Columns name
    
    ('MJD','X-P (arcsec)', 'Y-P (arcsec)', 'UT1UTC (E-7S)','LOD (E-7S/D)','S-X (E-6" arcsec)','S-Y (E-6" arcsec)',
    'S-UT (E-7S)','S-LD (E-7S/D)','NR (E-6" arcsec)', 'NF (E-6" arcsec)', 'NT (E-6" arcsec)',
    'X-RT (arcsec/D)','Y-RT (arcsec/D)','S-XR (E-6" arcsec/D)','S-YR (E-6" arcsec/D)', 'C-XY', 'C-XT',
    'C-YT', 'DPSI', 'DEPS','S-DP','S-DE')

    """
    
    caminho_arq = file_path_in

    #### FIND DELIVERY DATE
    name = os.path.basename(caminho_arq)

    if not ac:
        ac = name[:3]

    if len(name) == 12:
        dt_delivery = conv.sp3name2dt(caminho_arq)
    elif len(name) == 38:
        dt_delivery = conv.sp3name_v3_2dt(caminho_arq)
    else:
        dt_delivery = conv.posix2dt(0)


    le = open(caminho_arq, 'r')
    letudo = le.readlines()
    le.close()
    tamanho = len(letudo) #usado para saber quantas linhas tem o arquivo

    #para = tamanho #onde o arquivo para de ser lido

    numeros = ['0','1','2','3','4','5','6','7','8','9']
    #le = 0
    #numlin = 0 #numero de linhas da matriz de epocas
    #numcol = 16 #numero de colunas que a matriz final deve ter

    ERP=[]


    if caminho_arq[-3:] in ('snx','ssc'):
        file = open(caminho_arq)
        Lines =  file.readlines()
        XPO_stk  = []
        XPO_std_stk = []
        YPO_stk  = []
        YPO_std_stk = []
        LOD_stk  = []
        LOD_std_stk = []
        MJD_stk = []
        marker = False

        for i in range(len(Lines)):

            if len(Lines[i].strip()) == 0:
                continue
            else:

                if Lines[i].split()[0] == '+SOLUTION/ESTIMATE':
                    marker = True

                if Lines[i].split()[0] == '-SOLUTION/ESTIMATE':
                    marker = False

                if utils.contains_word(Lines[i],'XPO') and marker:
                    # Doy = (Lines[i][30:33])
                    # Year = (Lines[i][27:29])
                    # Pref_year = '20'
                    # Year = int(Pref_year+Year)
                    # Date = conv.doy2dt(Year,Doy)
                    Date = conv.datestr_sinex_2_dt(Lines[i].split()[5])
                    XPO = float(Lines[i][47:68])*(10**-3)
                    XPO_std = float(Lines[i][69:80])*(10**-3)
                    XPO_stk.append(XPO)
                    XPO_std_stk.append(XPO_std)
                    MJD_stk.append(conv.dt2MJD(Date))
                    #MJD_stk.append(cmg.jd_to_mjd(cmg.date_to_jd(Date.year,Date.month,Date.day)))
                    
                if utils.contains_word(Lines[i],'YPO') and marker:
                    # Doy = (Lines[i][30:33])
                    # Year = (Lines[i][27:29])
                    # Pref_year = '20'
                    # Year = int(Pref_year+Year)
                    # Date = conv.doy2dt(Year,Doy)
                    Date = conv.datestr_sinex_2_dt(Lines[i].split()[5])
                    YPO = float(Lines[i][47:68])*(10**-3)
                    YPO_std = float(Lines[i][69:80])*(10**-3)
                    YPO_stk.append(YPO)
                    YPO_std_stk.append(YPO_std)
                    MJD_stk.append(conv.dt2MJD(Date))
                    #MJD_stk.append(cmg.jd_to_mjd(cmg.date_to_jd(Date.year,Date.month,Date.day)))

                    
                if utils.contains_word(Lines[i],'LOD') and marker:
                    # Doy = (Lines[i][30:33])
                    # Year = (Lines[i][27:29])
                    # Pref_year = '20'
                    # Year = int(Pref_year+Year)
                    # Date = conv.doy2dt(Year,Doy)
                    Date = conv.datestr_sinex_2_dt(Lines[i].split()[5])
                    LOD = float(Lines[i][47:68])*(10**+4)
                    LOD_std = float(Lines[i][69:80])*(10**+4)
                    LOD_stk.append(LOD)
                    LOD_std_stk.append(LOD_std)
                    #MJD_stk.append(cmg.jd_to_mjd(cmg.date_to_jd(Date.year,Date.month,Date.day)))
                    MJD_stk.append(conv.dt2MJD(Date))

        MJD = list(sorted(set(MJD_stk)))
        if len(LOD_stk) == 0:
                LOD_stk = ['0']*len(MJD)
                LOD_std_stk = ['0']*len(MJD)


        for i in range(len(MJD)):

            ERP_data = [ac, MJD[i], XPO_stk[i], YPO_stk[i], 0, LOD_stk[i], XPO_std_stk[i], YPO_std_stk[i],
                     0, LOD_std_stk[i], 0, 0, 0, 0, 0, 0, 0, dt_delivery]

            ERP.append(ERP_data)



    if ac in ('COD','cod','com', 'cof', 'grg', 'mit', 'sio','igs','igr'):
        for i in range(tamanho+1):
            linhaatual = linecache.getline(caminho_arq, i)
            if linhaatual[0:1] in numeros:
                ERP_data = linhaatual.split()
                for j in range(len(ERP_data)):
                    ERP_data[j] = float(ERP_data[j])
                ERP_data.insert(0,ac)
                ERP_data[2] =  ERP_data[2]*(10**-6)
                ERP_data[3] =  ERP_data[3]*(10**-6)
                ERP_data[13] = ERP_data[13]*(10**-6)
                ERP_data[14] = ERP_data[14]*(10**-6)
                del ERP_data[17:]
                ERP_data.append(dt_delivery)

                ERP.append(ERP_data)
                
        linecache.clearcache()
        

#        Erp_end = pd.DataFrame(ERP, columns=['AC','MJD','X-P', 'Y-P', 'UT1UTC(UT1 -TAI)','LOD','S-X','S-Y','S-UT','S-LD','NR', 'NF', 'NT',
#                                                 'X-RT','Y-RT','S-XR','S-YR',"Delivery_date"])
#        return Erp_end
#

    if ac in ('wum','grg','esa', 'mit', 'ngs', 'sio'):
        for i in range(tamanho+1):
            linhaatual = linecache.getline(caminho_arq, i)
            if linhaatual[0:1] in numeros:
                ERP_data = linhaatual.split()
                for j in range(len(ERP_data)):
                    ERP_data[j] = float(ERP_data[j])
                ERP_data.insert(0,ac)
                ERP_data[2] =  ERP_data[2]*(10**-6)
                ERP_data[3] =  ERP_data[3]*(10**-6)
                ERP_data[13] = ERP_data[13]*(10**-6)
                ERP_data[14] = ERP_data[14]*(10**-6)
                ERP_data.append(dt_delivery)

                ERP.append(ERP_data)
                
        linecache.clearcache()

#        Erp_end = pd.DataFrame(ERP, columns=['AC','MJD','X-P', 'Y-P', 'UT1UTC(UT1 -TAI)','LOD','S-X','S-Y','S-UT','S-LD','NR', 'NF', 'NT',
#                                                 'X-RT','Y-RT','S-XR','S-YR'])
#        return Erp_end
#
    if ac in ('gbm', 'gfz','gfr',"p1_","p1r"):
        for i in range(tamanho+1):
            linhaatual = linecache.getline(caminho_arq, i)
            if linhaatual[0:1] in numeros:
                ERP_data = linhaatual.split()
                for j in range(len(ERP_data)):
                    ERP_data[j] = float(ERP_data[j])
                ERP_data.insert(0,ac)
                ERP_data[2] =  ERP_data[2]*(10**-6)
                ERP_data[3] =  ERP_data[3]*(10**-6)
                ERP_data[13] = ERP_data[13]*(10**-6)
                ERP_data[14] = ERP_data[14]*(10**-6)
                ERP_data.append(dt_delivery)

                ERP.append(ERP_data)
        linecache.clearcache()

#        Erp_end = pd.DataFrame(ERP, columns=['AC','MJD','X-P', 'Y-P', 'UT1UTC(UT1 -TAI)','LOD','S-X','S-Y','S-UT','S-LD','NR', 'NF', 'NT',
#                                                 'X-RT','Y-RT','S-XR','S-YR'])  ##EH TBM O RATE XY POR DIA??????
#        return Erp_end

    header = []
    if ac in ('emr'):
        for i in range(tamanho+1):
            linhaatual = linecache.getline(caminho_arq, i)
            if linhaatual == 'EOP  SOLUTION':
                del ERP_data[:]
                header = ['EOP  SOLUTION']
            if linhaatual[0:1] in numeros and 'EOP  SOLUTION' in header:
                ERP_data = linhaatual.split()
                for j in range(len(ERP_data)):
                    ERP_data[j] = float(ERP_data[j])
                ERP_data.insert(0,ac)
                ERP_data[2] =  ERP_data[2]*(10**-6)
                ERP_data[3] =  ERP_data[3]*(10**-6)
                ERP_data[13] = ERP_data[13]*(10**-6)
                ERP_data[14] = ERP_data[14]*(10**-6)
                del ERP_data[17:]
                ERP_data.append(dt_delivery)

                ERP.append(ERP_data)
        linecache.clearcache()



    Erp_end = pd.DataFrame(ERP, columns=['AC','MJD','X-P', 'Y-P', 
                                         'UT1UTC(UT1 -TAI)','LOD',
                                         'S-X','S-Y','S-UT','S-LD',
                                         'NR', 'NF', 'NT',
                                         'X-RT','Y-RT','S-XR','S-YR',
                                         'Delivered_date'])

        
    return Erp_end


### read_erp2 = read_erp
    
def read_erp_multi(path_list, 
                   return_array=False,
                   smart_mode=True,
                   ac=None):
    """
    

    Parameters
    ----------
    path_list : list
        a list of ERP files.
    return_array : ool, optional
        retur results in Array or DataFrame. The default is False.
    smart_mode : bool, optional
        keep only the latest value (True is recommended). 
        The default is True.
    ac : str, optional
        select a specific AC. The default is None.

    Returns
    -------
     Array or DataFrame
        output ERPs.

    """
    path_list = sorted(path_list)
    Lstk = []
    for path in path_list:
        L = read_erp(path,ac)
        Lstk.append(L)

    M = np.vstack(Lstk)

    if smart_mode:
        Msmart = np.array([])
        for ilm , lm in enumerate(np.flipud(M)):
            if ilm == 0:
                Msmart = np.array([lm])
            elif lm[0] in Msmart[:,0]:
                continue
            else:
                Msmart = np.vstack((Msmart,lm))

        M = np.flipud(Msmart)

    if return_array:
        return M
    else:
        return pd.DataFrame(M)



def read_bull_B(file_path_in):
    """
    Read an Bulletin B file

    Parameters
    ----------
    file_path_in :  str
        Path of the file in the local machine.
    Returns
    -------
    DFout : pandas DataFrame
        Returns a panda table format with the data extracted from the file.
    """

    if not utils.is_iterable(file_path_in):
        file_path_in = [file_path_in]

    file_path_in = sorted(file_path_in)

    DFstk = []

    for path_solo in file_path_in:
        S = utils.extract_text_between_elements(path_solo,
                                                "1 - DAILY FINAL VALUES" ,
                                                "2 - DAILY FINAL VALUES" )

        L = S.replace('\t','\n').split("\n")

        L2 = []
        for e in L:
            if len(e) > 0:
                if e[0] !=  " ":
                    L2.append(e)
        L3 = []
        for e in L2:
            L4 = []
            for ee in e.split():
                L4.append(float(ee))
            L3.append(L4)

        DF = pd.DataFrame(np.vstack(L3))
        DFstk.append(DF)

    DFout = pd.concat(DFstk)
    DFout.columns = ["year","month","day","MJD","x","y","UT1-UTC","dX","dY",
                  "x err","y err","UT1 err","X err","Y err"]

    return DFout

    
def read_erp_snx(snx_in):
    """
    Read ERP in SINEX Format

    Parameters
    ----------
    snx_in : str
        path of the SINEX file.

    Returns
    -------
    None.

    """
    
    log.warning("read_erp_snx not properly IMPLEMENTED !!!!!")
    
    utils.extract_text_between_elements_2(snx_in,
                                          r'\+SOLUTION/ESTIMATE',
                                          '-SOLUTION/ESTIMATE',
                                          return_string = True)
    
    

def read_eop_C04(file_path_in):
    """
    read EOP C04 file

    Parameters
    ----------
    file_path_in : str
        path of the EOP C04 file.

    Returns
    -------
    DF : DataFrame
        out DataFrame.

    """
    DF = pd.read_csv(file_path_in,
                     delim_whitespace=True,
                     skip_blank_lines=True,
                     skiprows=13,header=None)
    
    cols = ['yyyy','mm','dd','MJD','x','y','UT1-UTC','LOD','dX','dY',
            'x_Err','y_Err','UT1-UTC_Err','LOD_Err','dX_Err','dY_Err']
    
    DF.columns = cols
    
    DF["epoch"] = conv.MJD2dt(DF.MJD)
    
    return DF

def read_eop_finals(file_path_in,precnut_model = 2000,
                    simplified_EOP_DF="mixed"):
    """
    read EOP finals file


    Parameters
    ----------
    file_path_in : str
        path of the EOP finals file.
    precnut_model : int, optional
        IAU Precession-Nutation Model. Valid values are 1980 and 2000.
        The default is 2000.

    Returns
    -------
    DF : DataFrame
        out DataFrame - Complete content of the finals file.
    DF2 : DataFrame
        out DataFrame - Simplified content of the finals file.
        Based on C04 EOP DataFrame structure.
        If no Bulletin-B values provided, use Bulletin-A.
        
    Note
    ----
    Format description
    https://github.com/marando/phpIERS/blob/master/src/Marando/IERS/default/readme.finals
    """
    
    if precnut_model == 2000:
        x_or_psi = "X"
        y_or_eps = "Y"
    elif precnut_model == 1980:
        x_or_psi = "PSI"
        y_or_eps = "EPS"   
    else:
        raise Exception
    
    cols = ((0,2    ),
    (2,4    ),
    (4,6    ),
    (6,7    ),
    (7,15   ),
    (15,16  ),
    (16,17  ),
    (17,18  ),
    (18,27  ),
    (27,36  ),
    (36,37  ),
    (37,46  ),
    (46,55  ),
    (55,57  ),
    (57,58  ),
    (58,68  ),
    (68,78  ),
    (78,79  ),
    (79,86  ),
    (86,93  ),
    (93,95  ),
    (95,96  ),
    (96,97  ),
    (97,106 ),
    (106,115),
    (115,116),
    (116,125),
    (125,134),
    (134,144),
    (144,154),
    (154,165),
    (165,175),
    (175,185))
        
    col_names = ['yyyy',
    'mm',
    'dd',
    '[blank]',
    'MJD',
    '[blank]',
    'flg_xy',
    '[blank]',
    'x_A',
    'x_A_Err',
    '[blank]',
    'y_A',
    'y_A_Err',
    '[blanks]',
    'flg_UT1-UTC',
    'UT1-UTC_A',
    'UT1-UTC_A_Err',
    '[blank]',
    'LOD_A',
    'LOD_A_Err',
    '[blanks]',
    'flg_nut',
    '[blank]',
    'd' + x_or_psi + '_A',
    'd' + x_or_psi + '_A_Err',
    '[blank]',
    'd' + y_or_eps + '_A',
    'd' + y_or_eps + '_A_Err',
    'x_B',
    'y_B',
    'UT1-UTC_B',
    'd' + x_or_psi + '_B',
    'd' + y_or_eps + '_B']
    
    DF = pd.read_fwf(file_path_in,colspecs=cols,header=None)
    DF.columns = col_names
    ### remove everything with NaN
    DF.dropna(axis=1,how="all",inplace=True)
    ### if not x_A provided, then its a blank line
    DF = DF[np.logical_not(np.isnan(DF['x_A']))]
    
    ### set mjd as int
    DF["MJD"]   = DF["MJD"].astype(np.int64)
    
    #### DF 2 is a simplified version, base on the C04 DF
    DF2 = pd.DataFrame()
    DF2["epoch"] = conv.MJD2dt(DF["MJD"])
    
    if simplified_EOP_DF == 'mixed':
        ### Create epoch and use B-values per default
        DF2[['MJD','x','y','UT1-UTC','LOD','dX','dY']] = DF[['MJD','x_B','y_B','UT1-UTC_B','LOD_A','dX_B','dY_B']]
        
        DF2["Bul"] = "B"
        ### If no B-values, use A-values
        DF2.loc[np.isnan(DF.x_B),'Bul'] = "A"
        DF2.loc[np.isnan(DF.x_B),'x'] = DF.loc[np.isnan(DF.x_B),'x_A']
        DF2.loc[np.isnan(DF.y_B),'y'] = DF.loc[np.isnan(DF.y_B),'y_A']
        DF2.loc[np.isnan(DF['UT1-UTC_B']),'UT1-UTC'] = DF.loc[np.isnan(DF['UT1-UTC_B']),'UT1-UTC_A']
        DF2.loc[np.isnan(DF.dX_B),'dX'] = DF.loc[np.isnan(DF.dX_B),'dX_A']
        DF2.loc[np.isnan(DF.dX_B),'dY'] = DF.loc[np.isnan(DF.dY_B),'dY_A']
                
    elif simplified_EOP_DF == "A":
        DF2[['MJD','x','y','UT1-UTC','LOD','dX','dY']] = DF[['MJD','x_A','y_A','UT1-UTC_A','LOD_A','dX_A','dY_A']]
        DF2["Bul"] = "A"
        DF2.dropna(axis=1,how="all",inplace=True)       

    elif simplified_EOP_DF == "B":
        DF2[['MJD','x','y','UT1-UTC','LOD','dX','dY']] = DF[['MJD','x_B','y_B','UT1-UTC_B','LOD_A','dX_B','dY_B']]
        DF2["Bul"] = "B"
        DF2.dropna(axis=1,how="all",inplace=True)

    else:
        log.error("check the simplified_EOP_DF parameter !!!")
        raise Exception
        
    return DF,DF2

########################################################################################################################################
############################################### READ IERS GUS

def read_IERS(file_path_in):
#    file_path_in =  '/dsk/mansur/MGEX_C/ERP_IERS'

    file = open(file_path_in,'r')
    fil = file.readlines()
    file.close()
    return fil


def read_IERS_info(fil, mjd):
    #    mjd = int(mjd)
    X, Y, ut1utc, dx, dy, xerr, yerr, ut1utcerr, dxerr, dyerr, LOD, sig_lod, OMEGA, sig_ome = [],[],[],[],[],[],[],[],[],[],[],[],[],[]
    PoleXY = []

#    mjd = 58230

    for i in range(len(fil)):
    #    linhaatual = linecache.getline(file_path_in, i)
        linhaatual = fil[i]
        if (linhaatual[15:20]) == str(mjd) and (linhaatual[104:105]) != '':
          X = float(linhaatual[23:30])*(10**-3)
          Y = float(linhaatual[31:39])*(10**-3)
          ut1utc = float(linhaatual[40:49])*(10**-3)
          dx = float(linhaatual[51:58])*(10**-3)
          dy = float(linhaatual[59:65])*(10**-3)
          xerr = float(linhaatual[68:74])*(10**-3)
          yerr = float(linhaatual[77:82])*(10**-3)
          ut1utcerr = float(linhaatual[86:93])*(10**-3)
          dxerr = float(linhaatual[94:100])*(10**-3)
          dyerr = float(linhaatual[101:106])*(10**-3)
        elif (linhaatual[15:20]) == str(mjd) and (linhaatual[66:68]) != '':
          LOD = float(linhaatual[22:29])*(10**-3)
          sig_lod = float(linhaatual[30:37])*(10**-3)
          OMEGA = float(linhaatual[38:53])*(10**-3)
          sig_ome = float(linhaatual[56:70])*(10**-3)

    line_data = [mjd,X,Y,ut1utc,dx,dy,xerr,yerr,ut1utcerr,dxerr,dyerr,LOD,sig_lod,OMEGA,sig_ome]
    PoleXY.append(line_data)

    PoleXY_data = pd.DataFrame(PoleXY, columns=['mjd','X','Y','ut1-utc','dx','dy','xerr','yerr','ut1utcerr','dxerr','dyerr','LOD','sig_lod','OMEGA','sig_ome'])
    PoleXY_data.set_index('mjd',inplace=True)
    return PoleXY_data

################################################## END READ IERS GUS