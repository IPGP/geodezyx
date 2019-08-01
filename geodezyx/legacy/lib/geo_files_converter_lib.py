# -*- coding: utf-8 -*-
"""
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
import genefun
import softs_runner


import datetime as dt
import numpy as np
import textwrap
import shutil
import os
import pandas
import re

import dateutil.parser
import geodetik as geok
import copy
import glob

import matplotlib.pyplot as plt

#for read rinex nav
from io import BytesIO,StringIO
from pandas import DataFrame,Series
from pandas.io.pytables import read_hdf

#
   #_____          __  __ _____ _______                                   _
  #/ ____|   /\   |  \/  |_   _|__   __|                                 | |
 #| |  __   /  \  | \  / | | |    | |      __ _  ___ _ __   ___ _ __ __ _| |
 #| | |_ | / /\ \ | |\/| | | |    | |     / _` |/ _ \ '_ \ / _ \ '__/ _` | |
 #| |__| |/ ____ \| |  | |_| |_   | |    | (_| |  __/ | | |  __/ | | (_| | |
  #\_____/_/    \_\_|  |_|_____|  |_|     \__, |\___|_| |_|\___|_|  \__,_|_|
                                          #__/ |
                                         #|___/
#

def list_stat_in_statinfo(statinfoin):
    print("WARN : DISCONTINUED : list_stat_in_statinfo doesnt work well, should be replaced by stat_list_in_station_info !!!!")
    listtemp = []
    for l in open(statinfoin):
        if l[0] != ' ':
            continue
        if l[1:5] == '\n':
            continue
    listtemp.append(l[1:5])
    listout = genefun.uniqify_list(listtemp)
    return listout


def read_station_info_solo(filein,stat,column_type="ulr"):
    """
    column_type : str
        "ulr" or "sopac"
    """
    dicout = {k: [] for k in ['Start','End','Rec','Ant',
              'AntHt','Dome','Ant N','Ant E','Vers',
              'SwVers','Rec SN','Ant SN']}

    if column_type == "sopac":
        rg_start_1 = slice(25,29)
        rg_start_2 = slice(30,33)
        rg_start_3 = slice(34,36)
        rg_start_4 = slice(37,39)
        rg_start_5 = slice(40,42)

        rg_end_1   = slice(44,48)
        rg_end_2   = slice(49,52)
        rg_end_3   = slice(53,55)
        rg_end_4   = slice(56,58)
        rg_end_5   = slice(59,61)

        rg_rec   = slice(96,114)
        rg_antht = slice(64,70)
        rg_ant   = slice(170,187)
        rg_dome  = slice(187,191)
        rg_antN  = slice(80,86)
        rg_antE  = slice(89,95)

        rg_vers     = slice(119,141)
        rg_sw_vers  = slice(141,148)
        rg_rec_sn   = slice(148,170)
        rg_ant_sn   = slice(170,187)

    elif column_type == "ulr":
        rg_start_1 = slice(25,29)
        rg_start_2 = slice(30,33)
        rg_start_3 = slice(34,36)
        rg_start_4 = slice(37,39)
        rg_start_5 = slice(40,42)

        rg_end_1   = slice(44,48)
        rg_end_2   = slice(49,52)
        rg_end_3   = slice(53,55)
        rg_end_4   = slice(56,58)
        rg_end_5   = slice(59,61)

        rg_rec   = slice(97,103)
        rg_antht = slice(64,70)
        rg_ant   = slice(134,140)
        rg_dome  = slice(142,146)
        rg_antN  = slice(80,86)
        rg_antE  = slice(89,95)

    else:
        print("ERR : check column_type")
        return None

    for l in open(filein):
        if l[1:5] == stat:
            # cas specifique du temps (repris d'une autre fct)
            f = l
            start = geok.doy2dt(int(f[rg_start_1]),
                                int(f[rg_start_2]),
                                int(f[rg_start_3]),
                                int(f[rg_start_4]),
                                int(f[rg_start_5]))
            dicout['Start'].append(start)
            if int(f[rg_end_1]) == 9999 or int(f[rg_end_2]) == 999:
                end = geok.doy2dt(2099,1,0,0,0)
            else:
                end = geok.doy2dt(int(f[rg_end_1]),
                                  int(f[rg_end_2]),
                                  int(f[rg_end_3]),
                                  int(f[rg_end_4]),
                                  int(f[rg_end_5]))
            dicout['End'].append(end)
            # Autres Donnees importantes
            dicout['Rec'].append(l[rg_rec].strip())
            dicout['AntHt'].append(float(l[rg_antht].strip()))
            dicout['Ant'].append(l[rg_ant].strip())
            dicout['Dome'].append(l[rg_dome].strip())
            dicout['Ant N'].append(l[rg_antN].strip())
            dicout['Ant E'].append(l[rg_antE].strip())

            if column_type == "sopac":
                dicout['Vers'].append(l[rg_vers].strip())
                dicout['SwVers'].append(l[rg_sw_vers].strip())
                dicout['Rec SN'].append(l[rg_rec_sn].strip())
                dicout['Ant SN'].append(l[rg_ant_sn].strip())

    return dicout

def read_station_info_solo_date(filein,stat,date,column_type="ulr"):
    DIC = read_station_info_solo(filein,stat,column_type=column_type)

    Start = DIC["Start"]
    End   = DIC["End"]

    i_good = None
    for i,(s,e) in enumerate(zip(Start,End)):
        if (s <= date) and (date <= e):
            if i_good:
                print("WARN : index of the good period has been already defined, something is weird ...")
            i_good = i

    if i_good is None:
        print('ERR : no corresponding metadata have been found in the station.info for :' )
        print(stat,date)
        return None

    dicout = dict()

    for k , v in DIC.items():
        dicout[k] = v[i_good]

    return dicout



#dic = read_station_info_solo('/home/pierre/Documents/stationinfo2gins/station.info.ovsg','HOUE')

def read_lfile_solo(filein,stat):
    """
    WEAK : should be improved with discontinuites management
    """
    for l in open(filein):
        X, Y, Z  = 0.,0.,0.
        vX,vY,vZ = 0.,0.,0.

        T = dt.datetime.now()

        if l[1:5] == stat:
            f = l[1:].split()

            X  = float(f[1])
            Y  = float(f[2])
            Z  = float(f[3])
            Ttmp = float(f[7])

            if np.isclose(Ttmp , 0.):
                T  = geok.convert_partial_year(2000.)
            else:
                T  = geok.convert_partial_year(Ttmp)

            if len(l) > 225: # if velocities are given (125 is arbitrary)
                vX = float(f[8])
                vY = float(f[9])
                vZ = float(f[10])
            else: # if no velocityis given
                vX = 0.
                vY = 0.
                vZ = 0.


            return X,Y,Z,T,vX,vY,vZ

    print("WARN : read_lfile_solo : no coords found for : " , stat)
    return None



def read_pbo_vel_file_solo(velfilein,stat):
    """
    return a LIST of panda dataframe (corrsponding to different discontinuites),
    where the keywords of the pbo .vel can be used directly

    exemple :
    ftp://data-out.unavco.org/pub/products/velocity/pbo.final_nam08.vel
    """
    fil = open(velfilein)
    skiplinelis = []
    for i,l in enumerate(fil):
        if l[0] not in (' ','*'):
            skiplinelis.append(i)

    data = pandas.read_table(open(velfilein),skiprows=skiplinelis,sep= ' *',header=max(skiplinelis)+1)
    #data['Ref_epoch']

    i_lis = []
    for i,s in enumerate(list(data['*Dot#'])):
        if s == stat.upper():
            i_lis.append(i)

    data2 = data.transpose()

    data_lis = []
    for i in i_lis:
        data_lis.append(data2[i])
    return data_lis


def read_globk_vel_file(velfile_in):

    p = velfile_in
    i = 0
    for l in open(p):
        i += 1
        if '(deg)' in l:
            break

    D = pandas.read_table(p,skiprows=i,header=-1,delim_whitespace = True )

    D.rename(columns={0: 'Long',
                      1: 'Lat',
                      2: 'E',
                      3: 'N',
                      4: 'E_adj',
                      5: 'N_adj',
                      6: 'E_sig',
                      7: 'N_sig',
                      8: 'Rho',
                      9: 'H',
                      10: 'H_adj',
                      11: 'H_sig',
                      12: 'SITE'}, inplace=True)

    return D


def read_sinex_discontinuity_solo(snxfile,stat,PorV = 'P'):

    flagxyz = False

    soln_lis  = []
    start_lis = []
    end_lis   = []
    PorV_lis  = []
    notes_lis = []

    for line in open(snxfile):

        if re.compile('SOLUTION/DISCONTINUITY').search(line):
            flagxyz = not flagxyz
            continue

        if line[0] != ' ':
            continue

        if flagxyz == True:
            fields = line.split()

            if line[0] != ' ':
                continue

            if fields[0] == stat.upper():

                if not fields[6] == PorV:
                    continue

                soln_lis.append(int(fields[2]))
                start_lis.append(geok.datestr_sinex_2_dt(fields[4]))
                end = geok.datestr_sinex_2_dt(fields[5])
                if end == dt.datetime(1970,1,1):
                    end = dt.datetime(2099,1,1)
                end_lis.append(end)
                PorV_lis.append(fields[6])
                notes_lis.append(line.split('-')[-1].strip())

    return start_lis , end_lis , notes_lis , PorV_lis

def stat_list_in_station_info(filein):
    statlist = []
    fstatinfo = open(filein,"r+")
    for line in fstatinfo:
        if  len(line) == 0:
            print('WARN : stat_list_in_station_info : empty line in the station.info !!!')
            continue
        elif  line[0] != ' ':
            continue
        elif ' ' == line[0]:
            f = line.split()
            statlist.append(f[0])
    out = sorted(list(set(statlist)))
    return out

def convert_statinfo2eqfile(statinfoin,eqfileout):
    outfileobj = open(eqfileout,'w+')
    stat_list = stat_list_in_station_info(statinfoin)
    for stat in stat_list:
        i = 0
        startlis , endlis = read_station_info_time_solo(statinfoin,stat)
        statbisproto = stat + '_GPS'
        for s,e in zip(startlis,endlis):
            if i !=0:
                statbis_list = list(statbisproto)
                if i > 9:
                    statbis_list[5] = str(i)[0]
                    statbis_list[6] = str(i)[1]
                else:
                    statbis_list[5] = str(i)
                statbis = ''.join(statbis_list)
            else:
                statbis = statbisproto
            sy = s.year
            smo = s.month
            sd = s.day
            sh = s.hour
            smi = s.minute
            ss = s.second
            ey = e.year
            emo = e.month
            ed = e.day
            eh = e.hour
            emi = e.minute
            es = e.second
            ligne = ' rename {}     {} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2}\n'.format(stat,statbis,sy,smo,sd,sh,smi,ey,emo,ed,eh,emi)
            outfileobj.write(ligne)
            i = i+1
    outfileobj.close()
    return None


def read_station_info_time_solo(filein,stat):
    ''' BROUILLON
        pour une station, donne les startdate et enddate d'une periode
        SORTIE : 2 x liste de dt '''

    startlis = []
    endlis = []


    for l in open(filein):

        if l[1:5] == stat:

            f = l[25:].split()

            start = geok.doy2dt(int(f[0]),int(f[1]),int(f[2]),int(f[3]),int(f[4]))
            startlis.append(start)

            if int(f[5]) == 9999 or int(f[6]) == 999:
                end = geok.doy2dt(2099,1,0,0,0)
            else:
                end = geok.doy2dt(int(f[5]),int(f[6]),int(f[7]),int(f[8]),int(f[9]))
            endlis.append(end)

    return startlis,endlis



def read_eqfile_time_solo(filein,stat):
    ''' BROUILLON
        pour une station, donne les startdate et enddate d'une periode
        SORTIE : 2 x liste de dt '''

    startlis = []
    endlis = []

    for l in open(filein):
        f = l.split()
        if f[1] == stat:

            start = dt.datetime(int(f[3]),int(f[4]),int(f[5]),int(f[6]),int(f[7]))
            startlis.append(start)

            if int(f[7]) == 9999 or int(f[8]) == 999:
                end = dt.datetime(2099,1,1,0,0)
            else:
                end = dt.datetime(int(f[8]),int(f[9]),int(f[10]),int(f[11]),int(f[12]))
            endlis.append(end)

    return startlis,endlis


def read_eqfile_as_dico(filein, delGPSfield = False):
    '''
    for an eqfile return a dico
    dic[STAT] = [disc1 , disc2 ...] where disc are
    the discont of the FIRST column

    delGPSfield manage the elimination of the line of type _GPS
    which are not regular
    EDIT : 150817 WHY ???

    '''

    startlis = []
    endlis = []

    outdico = dict()

    for l in open(filein):
        f = l.split()
        stat = f[1]
        if stat not in outdico:
            outdico[stat] = []
        if f[2][-3] == 'G' and delGPSfield:
            continue
        datelis = [f[3],f[4],f[5],f[6],f[7],f[8]]
        datelis = [int(e) for e in datelis]
        outdico[stat].append(dt.datetime(*datelis))

    return outdico

def write_eqfile_from_dico(dicoin,outdir,outfilename):
    outfilepath = os.path.join(outdir,outfilename)
    fil = open(outfilepath,'w+')

    from collections import OrderedDict
#    dicoin = OrderedDict(sorted(dicoin.items(), key=lambda dicoin: dicoin[1]))
#    dicoin = sorted(dicoin, key=dicoin.get)
    for k,v in dicoin.items():
        for i in range(len(v)):
            stat = k
            if i+2 > 9:
                statbis = stat + '_' + str(i+2) + 'S'
            else:
                statbis = stat + '_' + str(i+2) + 'PS'

            s = v[i]
            if i == len(v)-1:
                e = dt.datetime(2099,1,1)
            else:
                e = v[i+1]

            sy = s.year
            smo = s.month
            sd = s.day
            sh = s.hour
            smi = s.minute
            ss = s.second
            ey = e.year
            emo = e.month
            ed = e.day
            eh = e.hour
            emi = e.minute
            es = e.second

            ligne = ' rename {}     {} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2} {:2}\n'.format(stat,statbis,sy,smo,sd,sh,smi,ey,emo,ed,eh,emi)
            fil.write(ligne)
    return None





#      _        _   _             _        __        ___     _____       _______ _____
#     | |      | | (_)           (_)      / _|      |__ \   / ____|   /\|__   __/ ____|
#  ___| |_ __ _| |_ _  ___  _ __  _ _ __ | |_ ___      ) | | |       /  \  | | | (___
# / __| __/ _` | __| |/ _ \| '_ \| | '_ \|  _/ _ \    / /  | |      / /\ \ | |  \___ \
# \__ \ || (_| | |_| | (_) | | | | | | | | || (_) |  / /_  | |____ / ____ \| |  ____) |
# |___/\__\__,_|\__|_|\___/|_| |_|_|_| |_|_| \___/  |____|  \_____/_/    \_\_| |_____/
#


def statname_of_catsfile(catsneu_path):
    catsneu_f = open(catsneu_path,'r')
    for line in catsneu_f:
        if 'Site' in line:
            return line.split(' ')[2]
    return 'XXXX'

def statinfo_2_cats(statinfo_path,catsneu_path):
    stat = statname_of_catsfile(catsneu_path)
    print('stat = ',stat)
    if stat == 'ABMF':
        return None
    stat_dic = read_station_info_solo(statinfo_path,stat)
    start = stat_dic['Start']
    antht = stat_dic['AntHt']

    outcats_path = os.path.dirname(catsneu_path) + os.path.basename(catsneu_path) + '.offset'
#    shutil.copyfile(catsneu_path,outcats_path)

    antht = [a for (s,a) in sorted(zip(start,antht))]
    start = sorted(start)
    print(stat_dic,stat)
    s0 = start[0]
    h0 = antht[0]
    s_stk = []
    for s, h in zip(start,antht):
        if h != h0:
            s_stk.append(s)
            h0 = h
    print(outcats_path)
    catsneuout_f = open(outcats_path,'w')
    catsneu_f = open(catsneu_path,'a+')
    lastline=''
    for line in catsneu_f:
        if lastline == line:
            continue
        elif 'Height' in line:
            catsneuout_f.write(line)
            for s in s_stk:
                catsneuout_f.write('# offsets : ' + str(geok.toYearFraction(s)) + ' 1\n')
        # cleaning outlier
        elif line[0] != '#':
            if (abs(float(line.split()[1])) > 1 or abs(float(line.split()[2])) > 1):
                print("aaaa")
            else:
                catsneuout_f.write(line)
        else:
            catsneuout_f.write(line)
        lastline = line


    catsneuout_f.close()
    # type 1 for only vertical (cf doc)
    return None

#      _        _   _             _        __        ___          _
#     | |      | | (_)           (_)      / _|      |__ \        (_)
#  ___| |_ __ _| |_ _  ___  _ __  _ _ __ | |_ ___      ) |   __ _ _ _ __  ___
# / __| __/ _` | __| |/ _ \| '_ \| | '_ \|  _/ _ \    / /   / _` | | '_ \/ __|
# \__ \ || (_| | |_| | (_) | | | | | | | | || (_) |  / /_  | (_| | | | | \__ \
# |___/\__\__,_|\__|_|\___/|_| |_|_|_| |_|_| \___/  |____|  \__, |_|_| |_|___/
#                                                            __/ |
#                                                           |___/

#def oceload_associed(statginsfile):
#    import os
#
#    aaa = os.path.abspath(statginsfile)
#    bbb = 'outocload'
#    filearg = open('temparg','wr')
#    filearg.write(aaa + '\n')
#    filearg.write(os.path.basename(aaa) + '/' + bbb)
#
#    os.system('bash /home/psakicki/scripts/exe_loadoce < temparg')
#

def receptor_gins_corrector(inRec):
    outRec = inRec
    outRec = outRec.replace('ASHTECH','AS')
    outRec = outRec.replace('TRIMBLE','TR')
    outRec = outRec.replace('LEICA GRX1200+GNS','LEICA GRX12')
    outRec = outRec.replace('LEICA GRX1200GGPRO','LEICA GRX12')
    outRec = outRec.replace('ROGUE','RG')
    return outRec

def station_info_2_gins(statinfoin,coordfilein,outfile,
                        coordfile_type ='pbovelfile',specific_stats_lis=[],
                        ellipsoid='GRS80',station_info_columns_type="sopac"):
    """
    Convert a GAMIT station.info to a GINS stations files

    Inputs :
        statinfoin : path of input station.info
        coodfilein : path of input coordinates files
        outfile : path of output GINS station file
        coordfile_type : GAMIT coordinates files type : 'lfile' or 'pbovelfile'
        specific_stats_lis : list of specific stats
        ellipsoid : 'GRS80'
        station_info_columns_type : it exists different subtypes of station.info ...
                                    handles "sopac" or "ulr" yet

    Output :
          path of GINS outfile
    """

    # Info de referentiel de la 1ère ligne
    refiel = 'i08i08'
    # Info de referentiel des lignes de réoccupations
    refchange = 'igs10'
    # code 'GRG3'
    grg = 'GRG3'
    # sigma sur les offsets
    sigoffset = 0.000
    # sigma sur les positions
    sigposi = 0.030
    # sigma sur les vitesses
    sigvit=0.005
    # Vitesses

    # ID DOMES du pays
    # ex : 971 pour la Guadeloupe
    # mettre 999 dans la plupart des cas
    # ATTENTION : string entre quotes
    #idprefix='971'
    #idprefix='990'
    #EDIT 1805, for very big station.info now idprfix is disabled
    #idprefix='990'

    # id initial, les stations auront un numero atribué de maniere décroissante
    # 97198 97197 97196 ...
    # idstatinit=99
    #EDIT 1805, for very big station.info now idstatinit is now at 99999
    idstatinit=99999

    # ======== FIN DES MODIFIABLES ========

    outEnd='000000'

    statgins_filobj=open(outfile,'w')

    if specific_stats_lis == []:
        listat = stat_list_in_station_info(statinfoin)
    else:
        listat = specific_stats_lis

    print(listat)

    sigvit = '{:.4f}'.format(sigvit).lstrip('0')

    header  = header_from_ellipsoid(ellipsoid)
    statgins_filobj.write(header)

    for stat in listat:

        print(stat)

        # On degresse l'ID quoiqu'il arrive après
        # pour conserver une homogenéité dans les pseudoDOMES
        idstatinit = idstatinit - 1

        if coordfile_type == 'lfile':
            output_lfile = None
            output_lfile = read_lfile_solo(coordfilein,stat)

            if not output_lfile:
                print("WARN : station_info_2_gins : skip station " , stat )
                continue
            else:
                X,Y,Z,T,vX,vY,vZ = output_lfile

        elif coordfile_type == 'pbovelfile':
            data = read_pbo_vel_file_solo(coordfilein,stat)
            if len(data) != 0:
                data = data[0]
                X = data['Ref_X']
                Y = data['Ref_Y']
                Z = data['Ref_Z']
                T = geok.MJD2dt(data['Ref_jday'])
                vX = data['dX/dt']
                vY = data['dY/dt']
                vZ = data['dZ/dt']
            else:
                print('WARN : no', stat,'found in',coordfilein)
                continue
        else:
            raise Exception('Wrong coordfile_type')

        if X == 0:
            continue

        #EDIT 1805, for very big station.info now idprefix is disabled
        #idstat = str(idprefix) + str(idstatinit)
        idstat = str(idstatinit)

        T = T.strftime('%d%m%y')

        dicout = read_station_info_solo(statinfoin,stat,
                                        column_type=station_info_columns_type)

        n = len(dicout['Start'])
        print("N = ", n)

        void='             '
        if stat == 'ABMF':
            idstat=str(97103)

        fortupsta=(idstat,'M001',stat,void,X,sigposi,Y,Z,vX,sigvit,vY,vZ,T,refiel)

        outputline="   {0} {1} {2}{3}{4:>-12.3f} ~{5:5.3f} {6:>-12.3f} ~{5:5.3f} {7:>-12.3f} ~{5:5.3f} {8:>-7.4f} ~{9} {10:>-7.4f} ~{9} {11:>-7.4f} ~{9} {12}        {13}".format(*fortupsta)

        outputline=outputline+'\n'
        statgins_filobj.write(outputline)


        for i in range(n):
            idligne  = idstat + '{:0>2d}'.format(i+1)
            outAnt   = dicout['Ant'][i]
            outAntHt = dicout['AntHt'][i]
            Xant,Yant,Zant = geok.ENU2XYZ_legacy(0,0,outAntHt,X,Y,Z)
            lastEnd  = outEnd
            outStart = dicout['Start'][i].strftime('%d%m%y')
            outEnd   = dicout['End'][i].strftime('%d%m%y')

            if outStart == lastEnd:
                outStart = (dicout['Start'][i] + dt.timedelta(days=1)).strftime('%d%m%y')

            if dicout['End'][i] == dt.datetime(2099,1,1,0,0,0):
                outEnd = '      '

            outDome = dicout['Dome'][i]
            outRec = dicout['Rec'][i]

            if '---------------' in outRec :
                continue
                # methode curieuse dans le stat.info sopac : on note
                # l'installation de l'antenne sans le rec : useless

            outRec = receptor_gins_corrector(outRec)[:11]
            outAnt = outAnt[:11]
            void='                     '

            formtuple=(idligne,grg,stat,outRec,Xant,sigoffset,Yant,sigoffset,Zant,sigoffset,outAnt,outDome,void,outStart,outEnd,refchange)

            #outputline ="%s %4s %4s %20s %5.3f %s %5.3f %s %5.3f %s %15s %s %s %s %s" \
            #%(idligne,grg,stat,outRec,Xant,sigoffset,Yant,sigoffset,Zant,sigoffset,outAnt,outDome,outStart,outEnd,refchange)

            outputline=" {0} {1} {2} {3:<11} {4:>-12.3f}  {5:5.3f} {6:>-12.3f}  {7:5.3f} {8:>-12.3f}  {9:5.3f}   {10:15} {11:4} {12} {13} {14} {15}".format(*formtuple)

            outputline=outputline+'\n'
            statgins_filobj.write(outputline)


        statgins_filobj.write('\n')

    if station_info_columns_type == "ulr":
        print("NOTA : input station.info is a 'ulr' type, the antenna/reciever type has to be corrected manually in the returned GINS file")


    return outfile

    #oceload_associed(outfile)


def stat_file_GINS_new_fmt(file_out_path,
                           STAT="STAT",
                           xyz=[0.,0.,0.],
                           ecc_une=[0.,0.,0.],
                           rec = "NONE",
                           ant = "NONE",
                           radom = "NONE",
                           ant_id = "0000"):

    F = open(file_out_path,"w+")

    Ellipsoid = """#VERSION 1.0
#RADIUS    0.6378137000000000E+07
#1/FLAT    0.2982572221010000E+03
# ITRF2014"""

    x , y , z = xyz[0] , xyz[1] , xyz[2]
    ecc_u , ecc_n , ecc_e = ecc_une[0] , ecc_une[1] , ecc_une[2]

    l1     = "   99999 M001 {:4}                             mar_xyz {:15.6f} 0.000000 {:15.6f} 0.000000 {:15.6f} 0.000000 vit_xyz  0.00000 0.0000  0.00000 0.0000  0.00000 0.0000  20100101_000000   EURA          i14i14"
    l1_fmt = l1.format(STAT,x,y,z)
    l2 = " 9999999 GRG3 {:4}      {:23}ecc_une {:15.6f} 0.000000 {:15.6f} 0.000000 {:15.6f} 0.000000            {:16}{:4} {:25}20000101_000000                 rinex"
    l2_fmt = l2.format(STAT,rec,ecc_u,ecc_n,ecc_e,ant,radom,ant_id)

    OUT = "\n".join((Ellipsoid , l1_fmt , l2_fmt))

    F.write(OUT)

    F.close()

    return file_out_path


#  _                     _               _                            _
# | |                   | |             | |                          | |
# | |     ___   __ _ ___| |__   ___  ___| |_ ___   _ __ ___  __ _  __| | ___ _ __
# | |    / _ \ / _` / __| '_ \ / _ \/ _ \ __/ __| | '__/ _ \/ _` |/ _` |/ _ \ '__|
# | |___| (_) | (_| \__ \ | | |  __/  __/ |_\__ \ | | |  __/ (_| | (_| |  __/ |
# |______\___/ \__, |___/_| |_|\___|\___|\__|___/ |_|  \___|\__,_|\__,_|\___|_|
#               __/ |
#              |___/
#


# ==== OBJECTS ====


class Event(object):
    def __init__(self):
        self.__Date_Installed            =  dt.datetime(1980,1,1)
        self.__Date_Removed              =  dt.datetime(2099,1,1)

    @property
    def Date_Installed(self):
        return self.__Date_Installed

    @Date_Installed.setter
    def Date_Installed(self,indate):
        if isinstance(indate,dt.datetime):
            self.__Date_Installed = indate
        elif 'CCYY-MM-DD' in indate:
            self.__Date_Installed = dt.datetime(1980,1,1,
                                                tzinfo=dateutil.tz.tzutc())
        elif 'XXXX' in indate:
            self.__Date_Installed = dt.datetime(1980,1,1,
                                                tzinfo=dateutil.tz.tzutc())
        else:
            self.__Date_Installed = dateutil.parser.parse(indate).replace(hour=0,
            minute=0, second=0, microsecond=0,tzinfo=dateutil.tz.tzutc())
#        print self.Date_Installed

    @property
    def Date_Removed(self):
        return self.__Date_Removed

    @Date_Removed.setter
    def Date_Removed(self,indate):
        if isinstance(indate,dt.datetime):
            self.__Date_Removed = indate
        elif 'CCYY-MM-DD' in indate:
            self.__Date_Removed = dt.datetime(2099,1,1,
                                              tzinfo=dateutil.tz.tzutc())
        elif 'XXXX' in indate:
            self.__Date_Removed = dt.datetime(2099,1,1,
                                              tzinfo=dateutil.tz.tzutc())
        else:
            self.__Date_Removed = dateutil.parser.parse(indate).replace(hour=0,
            minute=0, second=0, microsecond=0,tzinfo=dateutil.tz.tzutc())
#        print self.__Date_Removed

class Reciever(Event):
    def __init__(self):
        Event.__init__(self)
        self.Receiver_Type             =  'XXXX'
        self.Satellite_System          =  'XXXX'
        self.Serial_Number             =  'XXXX'
        self.Firmware_Version          =  'XXXX'
        self.Elevation_Cutoff_Setting  =  'XXXX'
        self.Temperature_Stabiliz      =  'XXXX'
        self.Additional_Information    =  'XXXX'
#        self.Date_Installed            = 'XXXX'
#        self.Date_Removed              = 'XXXX'

    def __repr__(self):
        return "{},{},{}".format(self.Receiver_Type ,self.Date_Installed,self.Date_Removed)
        return "{},{},{}".format(self.Antenna_Type,self.Date_Installed,self.Date_Removed)

    def FirmwareSmart(self):
        aaa = str(self.Firmware_Version).split()
        for a in aaa:
            try:
                return float(a)
            except:
                continue
        return 0.


class Antenna(Event):
    def __init__(self):
        Event.__init__(self)
        self.Antenna_Type             =  'XXXX'
        self.Serial_Number            =  'XXXX'
        self.Antenna_Reference_Point  =  'XXXX'
        self.Up_Ecc                   =  0.
        self.North_Ecc                =  0.
        self.East_Ecc                 =  0.
        self.Alignment_from_True_N    =  0.
        self.__Antenna_Radome_Type    =  'XXXX'
        self.Radome_Serial_Number     =  'XXXX'
        self.Antenna_Cable_Type       =  'XXXX'
        self.Antenna_Cable_Length     =  'XXXX'
        self.Additional_Information   =  'XXXX'

    @property
    def Antenna_Radome_Type(self):
        return self.__Antenna_Radome_Type

    @Antenna_Radome_Type.setter
    def Antenna_Radome_Type(self,inradome):
        if inradome == '':
            self.__Antenna_Radome_Type = 'NONE'
        else:
            self.__Antenna_Radome_Type = inradome[0:4].upper()

    def __repr__(self):
        return "{},{},{}".format(self.Antenna_Type,self.Date_Installed,self.Date_Removed)

    def AntTypSmart(self):
        #Elimination of the Radome type
        if len(self.Antenna_Type.split()) > 1:
            return self.Antenna_Type.split()[0]
        else:
            return self.Antenna_Type

    def ARPSmart(self):
        #return 'DH'+self.Antenna_Reference_Point
        # see http://rses.anu.edu.au/geodynamics/gps/papers/gamit/apdx2_ps.pdf
        #return 'DHPAB'
        return 'DHARP'


class Site(object):
    def __init__(self):
        self.Site_Name               =  'XXXX'
        self.Four_Character_ID       =  'XXXX'
        self.Monument_Inscription    =  'XXXX'
        self.__IERS_DOMES_Number     =  'XXXX'
        self.CDP_Number              =  'XXXX'
        self.Monument_Description    =  'XXXX'
        self.Height_of_the_Monument  =  'XXXX'
        self.Monument_Foundation     =  'XXXX'
        self.Foundation_Depth        =  'XXXX'
        self.Marker_Description      =  'XXXX'
        self.Date_Installed          =  'XXXX'
        self.Geologic_Characteristic =  'XXXX'
        self.Bedrock_Type            =  'XXXX'
        self.Bedrock_Condition       =  'XXXX'
        self.Fracture_Spacing        =  'XXXX'
        self.Fault_zones_nearby      =  'XXXX'
        self.Distance_activity       =  'XXXX'
        self.Additional_Information  =  'XXXX'

    @property
    def IERS_DOMES_Number(self):
        return  self.__IERS_DOMES_Number
    @IERS_DOMES_Number.setter
    def IERS_DOMES_Number(self,indomes):
        #if indomes == '':
        #    self.__IERS_DOMES_Number = '99999M001'
        #else:
        #    self.__IERS_DOMES_Number = indomes

        if re.search('[0-9]{4}M[0-9]{3}', indomes.strip()):
            self.__IERS_DOMES_Number = indomes.strip()
        else:
            self.__IERS_DOMES_Number = '99999M001'

    def __repr__(self):
        return "{},{},{}".format(self.Site_Name,self.Four_Character_ID,self.Date_Installed)


class Location(object):
    def __init__(self):
        self.City_or_Town            = 'XXXX'
        self.State_or_Province       = 'XXXX'
        self.Country                 = 'XXXX'
        self.Tectonic_Plate          = 'XXXX'
        self.X_coordinate_m          = 0.
        self.Y_coordinate_m          = 0.
        self.Z_coordinate_m          = 0.
        self.Latitude                = 0.
        self.Longitude               = 0.
        self.Elevation_m_ellips      = 0.
        self.Additional_Information  = 'XXXX'

        self.X_velocity              = 0.
        self.Y_velocity              = 0.
        self.Z_velocity              = 0.

        self.X_coordinate_sigma      = 0.008
        self.Y_coordinate_sigma      = 0.008
        self.Z_coordinate_sigma      = 0.008

        self.X_velocity_sigma        = 0.008
        self.Y_velocity_sigma        = 0.008
        self.Z_velocity_sigma        = 0.008

        self.Reference_epoch         = dt.datetime(2005,1,1)

    def __repr__(self):
        return "{},{},{},{}".format(self.X_coordinate_m,self.Y_coordinate_m,
                                 self.Z_coordinate_m,self.Reference_epoch)


    def export_as_string():
        Str_list = []
        Str_list.append('City or Town             : {:}'.format(self.City_or_Town))
        Str_list.append('State or Province        : {:}'.format(self.State_or_Province))
        Str_list.append('Country                  : {:}'.format(self.Country))
        Str_list.append('Tectonic Plate           : {:}'.format(self.Tectonic_Plate))
        Str_list.append('Approximate Position (ITRF)')

        X = self.X_coordinate_m
        Y = self.Y_coordinate_m
        Z = self.Z_coordinate_m

        lat , lon , h = geok.XYZ2GEO(X,Y,Z)
        lat_deg , lat_min , lat_sec = geok.deg2deg_dec2dms(lat)
        lon_deg , lon_min , lon_sec = geok.deg2deg_dec2dms(lon)

        Str_list.append('X coordinate (m)       : {:}'.format(X))
        Str_list.append('Y coordinate (m)       : {:}'.format(Y))
        Str_list.append('Z coordinate (m)       : {:}'.format(Z))

        Str_list.append('Latitude (N is +)      : {:+3d}{:2d}{:5.2d}'.format((lat_deg,lat_min,lat_sec)))
        Str_list.append('Longitude (E is +)     : {:+3d}{:2d}{:5.2d}'.format((lat_deg,lat_min,lat_sec)))
        Str_list.append('Elevation (m,ellips.)  : {:7.1f}'.format(h))
        Str_list.append('Additional Information   : N/A')

        out_str = "\n".join(Str_list)

        return out_str



# ===== functions =====

def read_blocks_logsheet(input_file, block_id):
    proto_objects = [None,Site(),Location(),Reciever(),Antenna()]

    blkstr = str(block_id) + '.'
    blkstrnxt = str(block_id+1) + '.'
    blkstrX = str(block_id) + '.x'

    objlis = []
    inblock = False
    

    for line in open(input_file,'rb'):
        
        ### manage unknown caracters
        line = str(line.decode('utf-8',errors='ignore'))


        if line == '\n':
            continue

        if line.startswith(blkstrnxt) or line.startswith(blkstrX):
            inblock = False
            break

        if line.startswith(blkstr) and (not 'Information' in line or block_id == 2):
            Obj = copy.copy(proto_objects[block_id])
            objlis.append(Obj)
            inblock = True

        if inblock:
            if 'Approximate Position' in line:
                continue
            prop = line[4:30].strip().replace(' ', '_').replace('/','_').replace('(','').replace(')','').replace('.','').replace(',','_')
            if 'Marker->ARP Up Ecc' in prop:
                prop = 'Up_Ecc'
            elif 'Marker->ARP North Ecc' in prop:
                prop = 'North_Ecc'
            elif 'Marker->ARP East Ecc' in prop:
                prop = 'East_Ecc'
            elif 'Latitude (N is +)' in prop:
                prop = 'Latitude'
            elif 'Longitude (E is +)' in prop:
                prop = 'Longitude'
            elif 'Elevation (m,ellips.)' in prop:
                prop = 'Elevation'
            data = line[31:].strip()
            #print prop , data
            try:
                setattr(Obj,prop,float(data))
            except:
                setattr(Obj,prop,data)

    return objlis

def mono_logsheet_read(logsheet_path,return_lists = True):
    '''
    from a logsheet, returns :

    * a "period list" i.e. a list of tuple
    (start date , end date , antenna object, reciever object)

    * a site object

    * a location object (station have new materials, but dont move ;) )

    if return_lists = True

    each object is returned in a list, so they can be managed immediatly by
    write_station_info_from_datalists
    '''
    sit_lis = read_blocks_logsheet(logsheet_path,1)
    loc_lis = read_blocks_logsheet(logsheet_path,2)
    rec_lis = read_blocks_logsheet(logsheet_path,3)
    ant_lis = read_blocks_logsheet(logsheet_path,4)

    ant_install = [e.Date_Installed for e in ant_lis]
    rec_install = [e.Date_Installed for e in rec_lis]

    # merging all install date
    all_install_date = sorted(list(set(ant_install + rec_install )))

    # for each install date, find the Rec/Ant couple
    date_ant_rec_couple_lis = []
    for d in all_install_date:
        # exclude extremal date
        if d > dt.datetime(2098,1,1,tzinfo=dateutil.tz.tzutc()):
            continue

        potential_ant = []
        potential_rec = []

        for a in ant_lis:
            if a.Date_Installed <= d < a.Date_Removed:
                potential_ant.append(a)
        for r in rec_lis:
            if r.Date_Installed <= d < r.Date_Removed:
                potential_rec.append(r)

        if len(potential_ant) == 0:
            print('WARN : missing Antenna info for',sit_lis[0].Four_Character_ID,d,'skip ...')
            continue
        if len(potential_rec) == 0:
            print('WARN : missing Receiver info for',sit_lis[0].Four_Character_ID,d,'skip ...')
            continue
        if len(potential_ant) != 1:
            print('WARN : several Antennas found for',sit_lis[0].Four_Character_ID,d)
        if len(potential_rec) != 1:
            print('WARN : several Receivers found for',sit_lis[0].Four_Character_ID,d)

        date_ant_rec_couple_lis.append((d , potential_ant[0] , potential_rec[0]))

    if len(date_ant_rec_couple_lis) != len(set(date_ant_rec_couple_lis)):
        print('bug 2')

    # construction of a period
    period_lis = []
    for i in range(len(date_ant_rec_couple_lis)):
        d1 = date_ant_rec_couple_lis[i][0]
        if i+1 == len(date_ant_rec_couple_lis):
            d2 = dt.datetime(2099,1,1,tzinfo=dateutil.tz.tzutc())
        else:
            d2 = date_ant_rec_couple_lis[i+1][0]
        a = date_ant_rec_couple_lis[i][1]
        r = date_ant_rec_couple_lis[i][2]

        period_lis.append((d1 , d2 , a , r))

    sit = sit_lis[0]
    loc = loc_lis[0]


    if return_lists:
        #each object is returned in a list, so they can be managed immediatly by
        #write_station_info_from_datalists
        return [period_lis] , [sit] , [loc]
    else:
        return period_lis   , sit , loc



def multi_logsheet_read(pathin,wildcardin='*log',return_dico=False):
    """
    if return_dico = False :

    return period_lis_lis , stat_lis , loc_lis

    this mode is useful for station.info generation

    if return_dico = True  :

    stations_dico['STAT'] =

    More human readable


    """
    fullpath = os.path.join(pathin,wildcardin)
    logsheet_list = sorted(glob.glob(fullpath))
    
    if not logsheet_list:
        print("ERR : no logsheets found, exiting ...")
        print("    ",fullpath)
        return None
    print("logsheet_list", logsheet_list)

    period_lis_lis = []
    stat_lis       = []
    loc_lis        = []

    stations_dico = dict()

    for ls in logsheet_list:
        try:
            p,s,l = mono_logsheet_read(ls)
        except:
            print("WARN :",ls,"skipped for unknown reason ...")
            print("       logsheet must be checked")
            continue

        period_lis_lis = period_lis_lis + p
        stat_lis       = stat_lis       + s
        loc_lis        = loc_lis        + l

        stations_dico[s[0].Four_Character_ID] = (p,s,l)


    if not return_dico:
        return period_lis_lis , stat_lis , loc_lis
    else:
        return stations_dico

def write_station_info_from_datalists(period_lis_lis,site_lis,location_lis,station_info_out_path):
    ''' datalists (period_lis_lis , stat_lis , loc_lis  ) are produced by :
        * mono_logsheet_read
        * multi_logsheet_read  '''

    si_file = open(station_info_out_path,'w')
    si_file.write('*SITE  Station Name      Session Start      Session Stop       Ant Ht   HtCod  Ant N    Ant E    Receiver Type         Vers                  SwVer  Receiver SN           Antenna Type     Dome   Antenna SN          \n')
    # writing the period in a station info
    for period_lis,sit,loc in zip(period_lis_lis,site_lis,location_lis):
        for period in period_lis:
            d1 , d2 , a , r = period
            # simplifiy name, removing space and non ascii char. ,
            # and fill with spaces to have 16 chars
            name = loc.City_or_Town[0:16].replace(' ','_')
            name = ''.join([i if ord(i) < 128 else '_' for i in name])
            if len(name) != 16:
                name = name.ljust(16)
            stat = sit.Four_Character_ID
            #managing the 4 chars
            stat = stat[0:4]
            if len(stat) != 4:
                stat = stat.ljust(4)

            strtup = stat , name , \
            d1.year,geok.dt2doy(d1),int(d1.hour),int(d1.minute),int(d1.second), \
            d2.year,geok.dt2doy(d2),int(d2.hour),int(d2.minute),int(d2.second), \
            float(a.Up_Ecc) , a.ARPSmart() , float(a.North_Ecc) , float(a.East_Ecc) , \
            r.Receiver_Type , r.Firmware_Version , r.FirmwareSmart() , \
            str(r.Serial_Number)[:20] , a.AntTypSmart() , \
            a.Antenna_Radome_Type , str(a.Serial_Number)[:20]

            outline =  " {}  {:16}  {} {} {:0>2d} {:0>2d} {:0>2d}  {} {} {:0>2d} {:0>2d} {:0>2d}  {:+.4f}  {}  {:+.4f}  {:+.4f}  {:20}  {:<20}  {:>5.2f}  {:<20}  {:15}  {:5}  {}\n".format(*strtup)
            si_file.write(outline)
    si_file.close()
    return None


def write_lfile_from_datalists(site_lis,location_lis,lfile_out_path):
    ''' datalists (period_lis_lis , stat_lis , loc_lis  ) are produced by :
        * mono_logsheet_read
        * multi_logsheet_read  '''
    lf_file = open(lfile_out_path,'w')

    
    proto_str = " {}_GPS {:.4f} {:.4f} {:.4f} {:+6.4f} {:+6.4f} {:+6.4f} {:6.1f} {:5.3f} {:5.3f} {:5.3f} {:6.4f} {:6.4f} {:6.4f}\n"
    
    for sit,loc in zip(site_lis,location_lis):
        vx  , vy  , vz  = 0,0,0
        sx  , sy  , sz  = 0,0,0
        svx , svy , svz = 0,0,0
        
        yr = geok.dt2year_decimal(loc.Reference_epoch)

        outline = proto_str.format(sit.Four_Character_ID ,
                                   loc.X_coordinate_m ,
                                   loc.Y_coordinate_m ,
                                   loc.Z_coordinate_m,
                                   vx  , vy  , vz,
                                   yr  ,
                                   sx  , sy  , sz,
                                   svx , svy , svz)
        
        lf_file.write(outline)
    lf_file.close()
    return None

def header_from_ellipsoid(ellipsoid):
    if ellipsoid in ( 'GRS80' , 'WGS84'):
        header=textwrap.dedent("""APRES ITRF2008 TOUTES STATIONS       ,  R= 6378137.00,   f= 298.257222101
 station type site                    X (m)   sig.        Y (m)   sig.        Z(m)    sig.    Xp (m/an)      Yp (m/an)      Zp (m/an)    date  plaque source
   -9999 M000 Earth c. of mass        0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                              010100
   -1998 EGB  Bias       (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                              010198 311298
    -199 EGB  Bias       (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                              010199 310199
     -61 EGC1 1/year Cos (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                                            dy ecc
     -51 EGS1 1/year Sin (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                                            dy ecc
     -62 EGC2 2/year Cos (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                                            dy ecc
     -52 EGS2 2/year Sin (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                                            dy ecc

   14001 M004 ZIMM              4331297.063 ~0.001   567555.878 ~0.001  4633133.936 ~0.001 -0.0135 ~.0001  0.0181 ~.0000  0.0121 ~.0001 010105        i08i08
 1400101 GRG3 ZIMM TR 4000SSE         0.000  0.000        0.000  0.000        0.000  0.000   TRM14532.00     NONE 3311A                 010593 050897 igs06
 1400102 GRG3 ZIMM TR 4000SSI         0.000  0.000        0.000  0.000        0.000  0.000   TRM14532.00     NONE 3311A                 060897 051198 igs06

   14001 M004 ZIMM              4331297.062 ~0.001   567555.882 ~0.001  4633133.935 ~0.001 -0.0135 ~.0001  0.0181 ~.0000  0.0121 ~.0001 010105        i08i08
 1400106 GRG3 ZIMM TR 4000SSI         0.000  0.000        0.000  0.000        0.000  0.000   TRM29659.00     NONE 99390                 061198 110803 igs06
 1400107 GRG3 ZIMM TR 4700            0.000  0.000        0.000  0.000        0.000  0.000   TRM29659.00     NONE 99390                 120803 210206 igs06
 1400108 GRG3 ZIMM TR NETRS           0.000  0.000        0.000  0.000        0.000  0.000   TRM29659.00     NONE 99390                 220206        igs06

            """)
    elif ellipsoid in ('GRIM' , 'EIGEN'):
        header=textwrap.dedent("""APRES ITRF2008 TOUTES STATIONS       ,  R= 6378136.46,   f= 298.257650
 station type site                    X (m)   sig.        Y (m)   sig.        Z(m)    sig.    Xp (m/an)      Yp (m/an)      Zp (m/an)    date  plaque source
   -9999 M000 Earth c. of mass        0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                              010100
   -1998 EGB  Bias       (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                              010198 311298
    -199 EGB  Bias       (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                              010199 310199
     -61 EGC1 1/year Cos (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                                            dy ecc
     -51 EGS1 1/year Sin (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                                            dy ecc
     -62 EGC2 2/year Cos (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                                            dy ecc
     -52 EGS2 2/year Sin (mm)         0.000 ~0.000        0.000 ~0.000        0.000 ~0.000                                                            dy ecc

   14001 M004 ZIMM              4331297.063 ~0.001   567555.878 ~0.001  4633133.936 ~0.001 -0.0135 ~.0001  0.0181 ~.0000  0.0121 ~.0001 010105        i08i08
 1400101 GRG3 ZIMM TR 4000SSE         0.000  0.000        0.000  0.000        0.000  0.000   TRM14532.00     NONE 3311A                 010593 050897 igs06
 1400102 GRG3 ZIMM TR 4000SSI         0.000  0.000        0.000  0.000        0.000  0.000   TRM14532.00     NONE 3311A                 060897 051198 igs06

   14001 M004 ZIMM              4331297.062 ~0.001   567555.882 ~0.001  4633133.935 ~0.001 -0.0135 ~.0001  0.0181 ~.0000  0.0121 ~.0001 010105        i08i08
 1400106 GRG3 ZIMM TR 4000SSI         0.000  0.000        0.000  0.000        0.000  0.000   TRM29659.00     NONE 99390                 061198 110803 igs06
 1400107 GRG3 ZIMM TR 4700            0.000  0.000        0.000  0.000        0.000  0.000   TRM29659.00     NONE 99390                 120803 210206 igs06
 1400108 GRG3 ZIMM TR NETRS           0.000  0.000        0.000  0.000        0.000  0.000   TRM29659.00     NONE 99390                 220206        igs06

            """)
    else:
        print('ERR : Check ellipsoid Name')
        return None
    return header




def write_station_file_gins_from_datalists(period_lis_lis,site_lis,location_lis,
                                           station_info_out_path,ellipsoid="GRS80"):
    ''' datalists (period_lis_lis , stat_lis , loc_lis  ) are produced by :
        * mono_logsheet_read
        * multi_logsheet_read  '''

    # Info de referentiel de la 1ère ligne
    refiel = 'i08i08'
    # Info de referentiel des lignes de réoccupations
    refchange = 'igs10'
    # code 'GRG3'
    grg = 'GRG3'
    # sigma sur les offsets
    sigoffset = 0.000
    # sigma sur les positions
    sigposi = 0.030
    # sigma sur les vitesses
    sigvit=0.005
    # Vitesses

    # ID DOMES du pays
    # ex : 971 pour la Guadeloupe
    # mettre 999 dans la plupart des cas
    # ATTENTION : string entre quotes
    idprefix='971'
    idprefix='990'

    # id initial, les stations auront un numero atribué de maniere décroissante
    # 97198 97197 97196 ...
    idstatinit=99

    void='             '
    void2='                     '

    statgins_filobj= open(station_info_out_path,'w+')

    header  = header_from_ellipsoid(ellipsoid)
    statgins_filobj.write(header)

    i = 0
    for i , ( period_lis , site , loca) in enumerate(zip(period_lis_lis,site_lis,location_lis)):
        T = loca.Reference_epoch.strftime('%d%m%y')
        idstatA = site.IERS_DOMES_Number[0:5]
        idstatB = site.IERS_DOMES_Number[5:9]

        # case where an auto DOMES must be given
        if int(idstatA) == 99999:
            idstatA = str(int(idstatA) - i)

        stat = site.Four_Character_ID
        # managing the 4 chars
        stat = stat[0:4]
        if len(stat) != 4:
            stat = stat.ljust(4)

        X = float(loca.X_coordinate_m)
        Y = float(loca.Y_coordinate_m)
        Z = float(loca.Z_coordinate_m)

        vX = float(loca.X_velocity)
        vY = float(loca.Y_velocity)
        vZ = float(loca.Z_velocity)

        sigposi = float(loca.X_coordinate_sigma)
        sigvit  = float(loca.X_velocity_sigma)

        fortupsta=(idstatA,idstatB,stat,void,X,sigposi,Y,Z,vX,sigvit,vY,vZ,T,refiel)
        outputline="   {0} {1} {2}{3}{4:>-12.3f} ~{5:5.3f} {6:>-12.3f} ~{5:5.3f} {7:>-12.3f} ~{5:5.3f} {8:>-7.4f} ~{9} {10:>-7.4f} ~{9} {11:>-7.4f} ~{9} {12}        {13}".format(*fortupsta)

        outputline=outputline+'\n'
        statgins_filobj.write(outputline)

        for j,peri in enumerate(period_lis):
            if j == 0:
                d_new_peri = 0
            else:
                d_new_peri = 1
            start , end , ant ,rec = peri
            idligne = idstatA + str(j).zfill(2)
            outRec = receptor_gins_corrector( rec.Receiver_Type )

            Xant,Yant,Zant = geok.ENU2XYZ_legacy(ant.East_Ecc,ant.North_Ecc,ant.Up_Ecc,X,Y,Z)

            sigoffset = 0.
            outAnt = ant.AntTypSmart()
            outDome = ant.Antenna_Radome_Type[0:4].ljust(4)
            outStart = (start + dt.timedelta(days=d_new_peri)).strftime('%d%m%y')
            if end == dt.datetime(2099,1,1,tzinfo=dateutil.tz.tzutc()):
                outEnd = '      '
            else:
                outEnd = end.strftime('%d%m%y')

            formtuple=(idligne,grg,stat,outRec,Xant,sigoffset,Yant,sigoffset,Zant,sigoffset,outAnt,outDome,void2,outStart,outEnd,refchange)

            outputline=" {0} {1} {2} {3:<11} {4:>-12.3f}  {5:5.3f} {6:>-12.3f}  {7:5.3f} {8:>-12.3f}  {9:5.3f}   {10:15} {11} {12} {13} {14} {15}".format(*formtuple)
            outputline=outputline+'\n'
            statgins_filobj.write(outputline)

        statgins_filobj.write('\n')
    return None

def smart_elt_list(list_raw,n_elt,replacement=''):
    try:
        return list_raw[n_elt]
    except:
        return replacement

def read_rinex_2_dataobjts(rinex_path):

    if genefun.empty_file_check(rinex_path):
        print('ERR : the RINEX file is empty ...')
        print('      ' , rinex_path)

        return None , None , None , None

    ant_raw   = genefun.grep(rinex_path,'ANT #',True)
    rec_raw   = genefun.grep(rinex_path,'REC #',True)
    xyz_raw   = genefun.grep(rinex_path,'APPROX POSITION XYZ',True).split()
    stat_raw  = genefun.grep(rinex_path,'MARKER NAME',True).split()
    domes_raw = genefun.grep(rinex_path,'MARKER NUMBER',True).split()
    d_hen_raw = genefun.grep(rinex_path,'ANTENNA: DELTA H/E/N',True).split()
    t_raw     = genefun.grep(rinex_path,'TIME OF FIRST OBS',True).split()

    Antobj, Recobj , Siteobj , Locobj = Antenna(),Reciever(),Site(),Location()

    Locobj.X_coordinate_m = float(smart_elt_list(xyz_raw,0))
    Locobj.Y_coordinate_m = float(smart_elt_list(xyz_raw,1))
    Locobj.Z_coordinate_m = float(smart_elt_list(xyz_raw,2))

    Antobj.Up_Ecc    = float(smart_elt_list(d_hen_raw,0))
    Antobj.East_Ecc  = float(smart_elt_list(d_hen_raw,1))
    Antobj.North_Ecc = float(smart_elt_list(d_hen_raw,2))

    Siteobj.Four_Character_ID = smart_elt_list(stat_raw,0)

    if re.search('[0-9]{5}[A-Z][0-9]{3}',smart_elt_list(domes_raw,0)):
        Siteobj.IERS_DOMES_Number =  smart_elt_list(domes_raw,0)
    else:
        Siteobj.IERS_DOMES_Number = str(np.random.randint(99999)).zfill(5) + 'M001' # smart_elt_list(domes_raw,0) A CHANGER !!!!!!!

    Antobj.Antenna_Radome_Type = ant_raw[36:40].strip()
    Antobj.Antenna_Type        = ant_raw[20:36].strip()

    Recobj.Receiver_Type = rec_raw[20:40].strip()

    # ======== OBSOLETE car utilisation de la fct rinex_start_end ========
    ##securité pour les secondes
    #if not (0 <= t_raw[5] < 60):
        #if t_raw[5] > 30:
            #t_raw[5] = 59
        #else:
            #t_raw[5] = 0

    #Date_Installed = dt.datetime(*[int(float(e)) for e in t_raw[0:6]],
                                   #tzinfo=dateutil.tz.tzutc())
    #Date_Removed = Date_Installed + dt.timedelta(days=1)


    Date_Installed , Date_Removed = softs_runner.rinex_start_end(rinex_path,
                                                                 add_tzinfo=1,
                                                                 verbose=0)


    Antobj.Date_Removed = Date_Removed
    Antobj.Date_Installed = Date_Installed

    Recobj.Date_Removed = Date_Removed
    Recobj.Date_Installed = Date_Installed

    Locobj.Reference_epoch = Date_Installed

    return Antobj , Recobj , Siteobj , Locobj


def write_station_file_gins_from_rinex(rinex_path,station_file_out,
                                       stat_code_filename_prioritary = True):
    Antobj , Recobj , Siteobj , Locobj = read_rinex_2_dataobjts(rinex_path)

    rnx_name           = os.path.basename(rinex_path)
    stat_code_filename = os.path.basename(rinex_path)[0:4]

    if Siteobj.Four_Character_ID.upper() != stat_code_filename.upper():
        print('WARN : different 4-char. code in RINEX header and filename')
        print('       for' , rnx_name , '(' , Siteobj.Four_Character_ID , ')')
        if stat_code_filename_prioritary:
            print('       keeping the 4-char. code in filename')
            Siteobj.Four_Character_ID = stat_code_filename.upper()


    write_station_file_gins_from_datalists([[(Antobj.Date_Installed, \
    Antobj.Date_Removed + dt.timedelta(days=1),Antobj,Recobj)]],[Siteobj],[Locobj],station_file_out)

    return None







# Convert sp3 2 gins
# brouillon
#
#satdic = OrderedDict()
#
#sp3in = "/media/psakicki/AF0E-DA43/CHAL_BUGS/TEMP_DATA/jp211236.sp3"
#orbginsout = '/home/psakicki/aaa_FOURBI/outgins'
#
#for l in open(sp3in):
#    if l[0] in ('#','+','%','/'):
#        continue
#
#    f = l.split()
#
#    if l[0] == '*':
#        ll = [int(float(e)) for e in f[1:]]
#        epoch = dt.datetime(*ll)
#        continue
#    if l[0] == 'P':
#        prn = int(f[0][2:])
#        x = float(f[1]) * 10**3
#        y = float(f[2]) * 10**3
#        z = float(f[3]) * 10**3
#
#        print epoch,prn,x,y,z
#
#        if not satdic.has_key(prn):
#            satdic[prn] = []
#
#        jjul = geok.dt2jjulCNES(epoch)
#        sec = geok.dt2secinday(epoch) + 19
#
#        satdic[prn].append((jjul,sec,x,y,z))
#
#fil = open(orbginsout,'w+')
#
#for k,vv in satdic.iteritems():
#    for v in vv:
#        if k in (11,13,14,18,20,28):
#            finalprn = 66600 + k
#        elif k in (2,15,17,21):
#            finalprn = 88800 + k
#        else:
#            finalprn = 77700 + k
#
#        jjul = v[0]
#        sec  = v[1]
#        x = v[2]
#        y = v[3]
#        z = v[4]
#
#        finalline =  '  {:5} {:6} {:12.6f} tai xyz ine  0 {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:+.14E} {:.3f}  {:.3E}\n'.format(finalprn,jjul,float(sec),0,0,0,0,0,0,0,0,0,x,y,z,0,0,0,0.,0)
#        fil.write(finalline)
#
#fil.close()

#                     _   _   _ __  __ ______             _____  _____
#                    | | | \ | |  \/  |  ____|   /\      / ____|/ ____|   /\
#  _ __ ___  __ _  __| | |  \| | \  / | |__     /  \    | |  __| |  __   /  \
# | '__/ _ \/ _` |/ _` | | . ` | |\/| |  __|   / /\ \   | | |_ | | |_ | / /\ \
# | | |  __/ (_| | (_| | | |\  | |  | | |____ / ____ \  | |__| | |__| |/ ____ \
# |_|  \___|\__,_|\__,_| |_| \_|_|  |_|______/_/    \_\  \_____|\_____/_/    \_\
#

def read_nmea(file_path , enuout = True ,
            startdate = dt.datetime(1980,1,1) ,
            enddate = dt.datetime(2099,1,1) , export_path = ''):
    """ if export_path != '', export a Matlab readable file to this path
        WARNING !!! la coord de ref est codée en dur, à coriger !!!!"""
    T    = []
    Lat  = []
    Long = []
    Haut = []
    Qual = []
    day = 999
    month = 999
    year = 999

    for l in open(file_path):
        if not l[0] == '$':
            continue
        if 'GPZDA' in l:
            f = l.split(',')
            day = int(f[2])
            month = int(f[3])
            year = int(f[4])
        if 'GPGGA' in l:
            f = l.split(',')
            h = int(f[1][0:2])
            m = int(f[1][2:4])
            s = int(f[1][4:6])
            qual = int(f[6])
            if  f[5] == 'W':
                EW = -1
            elif f[5] == 'E':
                EW = 1
            else:
                print('ProBleme')
            if day == 999:
                continue
            if h == 0 and m == 0 and s == 0:
                continue # c'est sale de zapper minuit mais bon ...
            t = dt.datetime(year,month,day,h,m,s)
            if not ( startdate <= t <= enddate):
                continue
            T.append(t)
            Lat.append(float(f[2][0:2]) + float(f[2][2:]) / 60.)
            Long.append( EW * (float(f[4][0:3]) + float(f[4][3:]) / 60.))
            Haut.append(float(f[9]))
            Qual.append(qual)

    X,Y,Z = geok.GEO2XYZ(Lat,Long,Haut)
    f,l,h = np.mean(Lat),np.mean(Long),np.mean(Haut)
    x0,y0,z0 = geok.GEO2XYZ(f,l,h)

    f,l,h =   43.4417981389 , 7.83481522597 , 6.59449264956
    x0,y0,z0 = 4595047.79934 , 632288.017869 , 4363273.52335
    E,N,U = geok.XYZ2ENU_2(X,Y,Z,x0,y0,z0)

    if export_path != '':
        outf = open(export_path,'w+')
        outf.write(' '.join(('#lat0,long0,h0 :',str(f),str(l),str(h),'\n')))
        outf.write(' '.join(('#x0 ,y0 , z0   :',str(x0),str(y0),str(z0),'\n')))
        for i in range(len(T)):
            datalis = [geok.dt2posix(T[i]),T[i].year,T[i].month,T[i].day,
                       T[i].hour,T[i].minute,T[i].second,Lat[i],Long[i],
                        Haut[i],X[i],Y[i],Z[i],E[i],N[i],U[i] ]
            datalis = [str(e) for e in datalis] + ['\n']
            outf.write(','.join(datalis)) #)
        outf.close()
    if enuout:
        return T,E,N,U,Qual
    else:
        return T,Lat,Long,Haut,Qual

#strtd = dt.datetime(2015,06,19,16,00)
#T,E,N,U,Q = read_nmea(file_path='/home/psakicki/geodesea_nav_final.dat',enuout=1,export_path='/home/psakicki/Documents/geodesea_nav_matrix.dat') #,startdate=strtd)

def plot_nmea(T,E,N,U,Q):
    """ a quick and dirty fct to plot the NMEA"""
    plt.close('all')

    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=1)

    plt.figure()
    plt.plot(T,U,'+')
    plt.title('U')
    ax = plt.gca()
    ax.yaxis.set_major_formatter(y_formatter)

    plt.figure()
    plt.plot(T,N,'+')
    plt.title('N')
    ax = plt.gca()
    ax.yaxis.set_major_formatter(y_formatter)

    plt.figure()
    plt.title('E')
    col = ['r'] * len(E)
    plt.scatter(T,E,Q,marker='+')
    ax = plt.gca()
    ax.yaxis.set_major_formatter(y_formatter)

    plt.figure()
    plt.title('E-N')
    plt.plot(E,N,'+')
    ax = plt.gca()
    ax.yaxis.set_major_formatter(y_formatter)

    plt.figure()
    plt.plot(T,Q,'+')
    ax = plt.gca()
    ax.yaxis.set_major_formatter(y_formatter)

def write_latlontime_file_4_OTPS_tide(outfilepath , lat , lon ,
                                      strt , end = None, sec_step=1 ,
                                      generate_inp_file = True,
                                      tide_model_path = './DATA/Model_atlas',
                                      set_relative_path_in_inp_file=True):
    """
    strt    : start as a datetime
    end     : it can be :
        length of the period in sec (integer)
        OR
        if None, strt is a list

    lat>0 - degrees North, lon>0 - degrees East
    lat<0 - degrees South, lon<0 - degrees West
    
    lat/lon can be a list, with the same size as dates_lis
    """
    fil = open(outfilepath,'w+')


    if not end:
        dates_lis = strt
        
    elif type(end) is dt.datetime:
        dates_lis  = [strt]
        while dates_lis[-1] <= end:
            dates_lis.append( dates_lis[-1] + dt.timedelta(seconds=sec_step))

    else:
        dates_lis = [strt + dt.timedelta(seconds=x) for x in np.arange(0,end+1,sec_step)]


    if not genefun.is_iterable(lat):
        lat = [float(lat)] * len(dates_lis)
        lon = [float(lon)] * len(dates_lis)
    else:
        if len(lat) != len(dates_lis) or len(lon) != len(dates_lis):
            print("ERR : lat/lon are list, but not with the same size as the dates list, please check")
            raise Exception
            

    print('INFO : first & last dates')
    print(dates_lis[0] , dates_lis[-1])

    for d,lat_i,lon_i in zip(dates_lis,lat,lon):
        y  = d.year
        m  = d.month
        da = d.day
        h  = d.hour
        mi = d.minute
        s  = d.second
        lin = '{:12.8}{:12.8}{:7}{:7}{:7}{:7}{:7}{:7}\n'.format(lat_i,lon_i,y,m,da,h,mi,s)
        fil.write(lin)

    if generate_inp_file:
        outdir   = os.path.dirname(outfilepath)
        outname  = os.path.basename(outfilepath) + '.inp'
        outoutname = os.path.basename(outfilepath) + '.out'

        fil2 = open(os.path.join(outdir,outname),'w+')

        fil2.write(tide_model_path + '\n' )
        if not set_relative_path_in_inp_file:
            fil2.write(outfilepath + ' \n')
        else:
            fil2.write("./" + os.path.basename(outfilepath) + ' \n')
            
        fil2.write('z'         + ' \n')
        fil2.write(''          + ' \n')
        fil2.write('AP'        + ' \n')
        fil2.write('geo'       + ' \n')
        fil2.write('1'         + ' \n')
        if not set_relative_path_in_inp_file:
            fil2.write( os.path.join(outdir,outoutname)  + ' \n')
        else:
            fil2.write("./" + outoutname + ' \n')

        fil2.close()

    return outdir

def read_OTPS_tide_file(pathin,print_line=False,return_array=True):
    """
    Return :
        latlis , lonlis , datlis , hlis
    """
    filout = open(pathin)

    latlis , lonlis , datlis , hlis = [],[],[],[]

    for i,l in enumerate(filout):
        f = l.split()
        if np.mod(i , 1000) == 0 and print_line:
            print('line', i)
        try:
            latlis.append(float(f[0]))
            lonlis.append(float(f[1]))
            datlis.append(dateutil.parser.parse(f[2] + ' ' + f[3]))
            hlis.append(float(f[4]))
        except:
            continue
    if return_array:
        return np.array(latlis) , np.array(lonlis) ,  np.array(datlis) , np.array(hlis)
    else:
        return latlis , lonlis , datlis , hlis

def read_snx_trop(snxfile,dataframe_output=True):
    """
    Read troposphere solutions from Troposphere SINEX
    """
    
    STAT , epoc = [] , []
    tro , stro , tgn , stgn , tge , stge = [] , [] , [] , [] , [] , []
    
    flagtrop = False
    
    for line in open(snxfile,"r",encoding = "ISO-8859-1"):
        
        if re.compile('TROP/SOLUTION').search(line):
            flagtrop = True
            continue
        
        if re.compile('-TROP/SOLUTION').search(line):
            flagtrop = False
            continue
        
        if line[0] != ' ':
            continue
        else:
            fields = line.split()
        
        if flagtrop ==True:
            
            STAT.append(fields[0].upper())
            if not ':' in fields[1]:
                epoc.append(geok.convert_partial_year(fields[1]))
            else:
                date_elts_lis = fields[1].split(':')
                yy =  int(date_elts_lis[0]) + 2000
                doy = int(date_elts_lis[1])
                sec = int(date_elts_lis[2])

                epoc.append(geok.doy2dt(yy,doy,seconds=sec))
            
            tro.append(float(fields[2]))
            stro.append(float(fields[3]))
            tgn.append(float(fields[4]))
            stgn.append(float(fields[5]))
            tge.append(float(fields[6]))
            stge.append(float(fields[7]))
            
    outtuple = \
    list(zip(*sorted(zip(STAT , epoc , tro , stro , tgn , stgn , tge , stge))))
    
    if dataframe_output:
        return Tropsinex_DataFrame(outtuple)
                
def Tropsinex_DataFrame(read_sinex_result):
     DF_Sinex = pandas.DataFrame.from_records(list(read_sinex_result)).transpose()
     colnam = ['STAT','epoc','tro','stro','tgn','stgn','tge','stge']
     DF_Sinex.columns = colnam
     cols_numeric = ['tro','stro','tgn','stgn','tge','stge']
     DF_Sinex[cols_numeric] = DF_Sinex[cols_numeric].apply(pandas.to_numeric, errors='coerce')
     
     return DF_Sinex
     
def read_sinex(snxfile,dataframe_output=False):

    STAT , soln , epoc, AC = [] , [] , [], []
    x , y , z , sx , sy , sz = [] , [] , [] , [] , [] , []
    vx , vy , vz , svx , svy , svz = [] , [] , [] , [] , [] , []
    
    start , end = [] , []

    flagxyz    = False
    flagepochs = False

    for line in open(snxfile,"r",encoding = "ISO-8859-1"):

        if re.compile('SOLUTION/ESTIMATE').search(line):
            flagxyz = not flagxyz
            continue
        
        if re.compile('SOLUTION/EPOCHS').search(line):
            flagepochs = not flagepochs
            continue        

        if line[0] != ' ':
            continue
        else:
            fields = line.split()


        if flagxyz == True:

            if fields[1] == 'STAX':
                split_sp3_name = os.path.basename(snxfile)
                AC.append(split_sp3_name[0:3])
                STAT.append(fields[2].upper())
                soln.append(fields[4])
                if not ':' in fields[5]:
                    epoc.append(geok.convert_partial_year(fields[5]))
                else:
                    date_elts_lis = fields[5].split(':')
                    yy =  int(date_elts_lis[0]) + 2000
                    doy = int(date_elts_lis[1])
                    sec = int(date_elts_lis[2])

                    epoc.append(geok.doy2dt(yy,doy,seconds=sec))

                x.append(float(fields[8]))
                sx.append(float(fields[9]))

            if fields[1] == 'STAY':
                y.append(float(fields[8]))
                sy.append(float(fields[9]))

            if fields[1] == 'STAZ':
                z.append(float(fields[8]))
                sz.append(float(fields[9]))


            if fields[1] == 'VELX':
                vx.append(float(fields[8]))
                svx.append(float(fields[9]))

            if fields[1] == 'VELY':
                vy.append(float(fields[8]))
                svy.append(float(fields[9]))

            if fields[1] == 'VELZ':
                vz.append(float(fields[8]))
                svz.append(float(fields[9]))
                
        if flagepochs:
            start_val = geok.datestr_sinex_2_dt(fields[4])
            end_val   = geok.datestr_sinex_2_dt(fields[5])
            start.append(start_val)
            end.append(end_val)
    
    if not vx:
        nanlist = [np.nan] * len(x)
        vx , vy , vz , svx , svy , svz = list(nanlist) , list(nanlist) , list(nanlist) , \
                                         list(nanlist) , list(nanlist) ,list(nanlist)
                                         

    outtuple = \
    list(zip(*sorted(zip(AC, STAT , soln , epoc , x , y , z , sx , sy , sz ,
                         vx , vy , vz , svx , svy , svz , start , end))))

    if dataframe_output:
        return sinex_DataFrame(outtuple)
    else:
        STAT , soln , epoc , x , y , z , sx , sy , sz , vx , vy , vz , svx , svy , svz , start , end= outtuple
        print('TIPS : this output can be converted directly to a Pandas DataFrame using sinex_DataFrame function')
        return STAT , soln , epoc , x , y , z , sx , sy , sz , vx , vy , vz , svx , svy , svz , start ,  end

def sinex_DataFrame(read_sinex_result):
    DF_Sinex = pandas.DataFrame.from_records(list(read_sinex_result)).transpose()
    colnam = ['AC', 'STAT' , 'soln' , 'epoc' , 'x' , 'y' , 'z' , 'sx' , 'sy' , 'sz' , 'vx' , 'vy' , 'vz' , 'svx' , 'svy' , 'svz' , 'start' , 'end']
    DF_Sinex.columns = colnam
    
    cols_numeric = ["x","y","z",
                "sx","sy","sz",
                "vx","vy","vz",
                "svx","svy","svz"]
        
    DF_Sinex[cols_numeric] = DF_Sinex[cols_numeric].apply(pandas.to_numeric, errors='coerce')
        

    return DF_Sinex


def read_sinex_versatile(sinex_path_in , id_block,
                         convert_date_2_dt = True,
                         header_line_idx = -1):
    """
    Read a block from a SINEX and return the data as a DataFrame

    Parameters
    ----------
    sinex_path_in : str
        Description param1

    id_block : str
        Name of the block (without "+" or "-")
        
    convert_date_2_dt : bool
        Try to convert a SINEX formated date as a python datetime
        
    header_line_idx : int
        If the block header contains several lines, use this line index
        Per default, the last (-1)
        For the first line, use 0
        
                
    Returns
    -------
    DF : Pandas DataFrame
        Returned DataFrame
    """

    id_block_strt = "\+" + id_block
    id_block_end  = "\-" + id_block
    
    Lines_list = genefun.extract_text_between_elements_2(sinex_path_in,
                                                         id_block_strt,
                                                         id_block_end)
    Lines_list = Lines_list[1:-1]
    
    if not Lines_list:
        print("ERR : read_sinex_versatile : no block found, ",id_block)
    
    #### Remove commented lines
    Lines_list_header = []
    Lines_list_OK     = []
    header_lines = True
    for i_l , l in enumerate(Lines_list):
        if not l[0] in (" ","\n") and header_lines:
            Lines_list_header.append(l)
        if l[0] == " ":
            header_lines = False
            Lines_list_OK.append(l)

    Lines_str  = "".join(Lines_list_OK)
                            
    if len(Lines_list_header) > 0:
        ### define the header
        header_line = Lines_list_header[header_line_idx]
        

        Header_split = header_line.split()
        Fields_size = [len(e)+1 for e in Header_split]

        print(header_line , Fields_size)



    
        ### Read the file
        DF = pandas.read_fwf(StringIO(Lines_str),widths=Fields_size)
        DF.set_axis(Header_split, axis=1, inplace=True)
        
        ### Rename the 1st column (remove the comment marker)
        DF.rename(columns={DF.columns[0]:DF.columns[0][1:]}, inplace=True)

    else: # no header in the SINEX
        DF = pandas.read_csv(StringIO(Lines_str),header=-1 ,
                             delim_whitespace=True)

    for col in DF.columns:
        if convert_date_2_dt and re.match("([0-9]{2}|[0-9]{4}):[0-9]{3}:[0-9]{5}",
                                          str(DF[col][0])):
            try:
                DF[col] = DF[col].apply(lambda x : geok.datestr_sinex_2_dt(x))
            except:
                print("WARN : read_sinex_versatile : convert date string to datetime failed")
                pass
        
    return DF




def read_sinex_bench_antenna(sinex_in):
    F = open(sinex_in,"r")

    line_good_stk = []
    for l in F:
        if l[0] == " ":
            line_good_stk.append(l)

    T = StringIO("".join(line_good_stk))

    header_cols = ["TYP" ,"STA_" ,"Occ", "Code","PT", "T",
                   "___x/up___" ,"___y/n____", "___z/e____" ,
                   "_Data_Start", "_Data_End__",
                   "__Antenna_type______" ,"Radome","__S/N__"]
    DFantenna = pandas.read_table(T,delim_whitespace = True,error_bad_lines=False,names=header_cols)

    return DFantenna


def read_rinex_nav(fn,writeh5=None,version=2):
    """
    Based on Michael Hirsch readRinexNav function
    Manage RINEX3
    http://gage14.upc.es/gLAB/HTML/GPS_Navigation_Rinex_v2.11.html
    """

    with open(fn , 'r') as f:
        #find end of header, which has non-constant length
        while True:

            line = f.readline()

            if 'RINEX VERSION' in line:
                if float(line.split()[0]) < 3:
                    version = 2
                else:
                    version = 3


            if 'END OF HEADER' in line: break

            if version == 2:
                startcol = 3 #column where numerical data starts
                nfloat=19 #number of text elements per float data number
                nline=7 #number of lines per record
            elif version == 3:
                startcol = 5 #column where numerical data starts
                nfloat=19 #number of text elements per float data number
                nline=7 #number of lines per record


        #handle frame by frame
        sv = []; epoch=[]; raws=''


        while True:
            headln = f.readline()
            if not headln: break
            #handle the header
            if version == 2:
                sv.append(headln[:2])
                year = int(headln[2:5])
                if 80<= year <=99:
                    year+=1900
                elif year<80: #good till year 2180
                    year+=2000

                epoch.append(dt.datetime(year =year,
                                      month   =int(headln[5:8]),
                                      day     =int(headln[8:11]),
                                      hour    =int(headln[11:14]),
                                      minute  =int(headln[14:17]),
                                      second  =int(headln[17:20]),
                                      microsecond=int(headln[21])*100000))


                raw = (headln[22:].rstrip() +
                       ''.join(f.readline()[startcol:].rstrip() for _ in range(nline-1))
                       +f.readline()[startcol:79].rstrip())




            if version == 3:
                fields = headln[:23].split(' ')
                sv.append(fields[0])
                year = int(fields[1])

                epoch.append(dt.datetime(year = year,
                                      month   =int(fields[2]),
                                      day     =int(fields[3]),
                                      hour    =int(fields[4]),
                                      minute  =int(fields[5]),
                                      second  =int(fields[6])))

                """
                now get the data.
                Use rstrip() to chomp newlines consistently on Windows and Python 2.7/3.4
                Specifically [:-1] doesn't work consistently as .rstrip() does here.
                """
                raw = (headln[23:].rstrip() + ' ' +
                       ''.join(f.readline()[startcol:].rstrip() + ' ' for _ in range(nline-1))
                       +f.readline()[startcol:80].rstrip())

            raws += raw + '\n'

    raws = raws.replace('D','E')

    strio = BytesIO(raws.encode())
    darr = np.genfromtxt(strio,delimiter=nfloat)

    nav= DataFrame(darr, epoch,
               ['SVclockBias','SVclockDrift','SVclockDriftRate','IODE',
                'Crs','DeltaN','M0','Cuc','Eccentricity','Cus','sqrtA','TimeEph',
                'Cic','OMEGA','CIS','Io','Crc','omega','OMEGA DOT','IDOT',
                'CodesL2','GPSWeek','L2Pflag','SVacc','SVhealth','TGD','IODC',
                'TransTime','FitIntvl',"spare1","spare2"])

    if version == 3:
        Const = [e[0]  for e in sv]
        SV    = [int(e[1:]) for e in sv]
        nav['Const'] = Series(np.array(Const), index=nav.index)
        nav['sv'] = Series(np.array(SV), index=nav.index)
    elif version == 2:
        rinexnav_type = os.path.basename(fn)[-1]
        if rinexnav_type == 'n':
            nav['Const'] = Series(['G']*len(nav.index), index=nav.index)
        else:
            nav['Const'] = Series([rinexnav_type.upper()]*len(nav.index), index=nav.index)
        nav['sv'] = Series(np.array(sv), index=nav.index)



    if False:
        h5fn = fn.with_suffix('.h5')
        print(('saving NAV data to {}'.format(h5fn)))
        nav.to_hdf(h5fn,key='NAV',mode='a',complevel=6,append=False)

    return nav


def unzip_gz_Z(inp_gzip_file,out_gzip_file='',remove_inp=False, force = False):
    """
    if out_gzip_file if not precised, file will be extracted in the same folder
    as inp_gzip_file

    .Z decompression is implemented, but is very unstable (avoid .Z, prefer .gz)
    """

    import gzip,zlib

    if inp_gzip_file.endswith('.gz'):
        is_gz = True
    elif inp_gzip_file.endswith('.Z'):
        is_gz = False
    else:
        is_gz = True

    if out_gzip_file == '':
        out_gzip_file = os.path.join(os.path.dirname(inp_gzip_file) ,
                                     '.'.join(os.path.basename(inp_gzip_file).split('.')[:-1]))


    if os.path.isfile(out_gzip_file) and not force:
        print('INFO : ' , out_gzip_file , 'already exists, skiping (use force option)')
        pass
    else:
        if is_gz:
            with gzip.open(inp_gzip_file, 'rb') as f_in, open(out_gzip_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        else:
            print("WARN : zlib decompress is unstable !! and .Z should be definitly avoided ... ")
            #str_object1 = open(inp_gzip_file, 'rb').read()
            #str_object2 = zlib.decompress(str_object1)
            #f = open(out_gzip_file, 'wb')
            #f.write(str_object2)
            #f.close()

            out_gzip_file = genefun.uncompress(inp_gzip_file)

        print('INFO : uncompressing ' + inp_gzip_file + " to " + out_gzip_file )

    if remove_inp and os.path.getsize(out_gzip_file) > 0:
        print("INFO : removing " + inp_gzip_file)
        os.remove(inp_gzip_file)

    return out_gzip_file




#                   ___    _____   ____  __  __ ______  _____                                  _
#                  |__ \  |  __ \ / __ \|  \/  |  ____|/ ____|                                | |
#   ___ _____   __    ) | | |  | | |  | | \  / | |__  | (___    _ __ ___  __ _ _   _  ___  ___| |_
#  / __/ __\ \ / /   / /  | |  | | |  | | |\/| |  __|  \___ \  | '__/ _ \/ _` | | | |/ _ \/ __| __|
# | (__\__ \\ V /   / /_  | |__| | |__| | |  | | |____ ____) | | | |  __/ (_| | |_| |  __/\__ \ |_
#  \___|___/ \_/   |____| |_____/ \____/|_|  |_|______|_____/  |_|  \___|\__, |\__,_|\___||___/\__|
#                                                                           | |
#                                                                           |_|

#def csv2DOMES(csvin):

#filein = '/home/psakicki/gin/TP/GWADA/gps_continu_full.csv'
#import pandas as pd
#df = pd.read_csv(open(filein))
#saved_column = list(df['Nom'])
#
#
#for i,row in df.iterrows():
#    name = eval(row['Nom'])
#    stat = row['Alias']
#    outdomesfile = open('/home/psakicki/gin/TP/GWADA/' + stat + '_DOMES_request.txt' ,'w+' )
#    outdomesfile.write(' 1. Request from (full name) : ' + '\n')
#    outdomesfile.write('Agency                   : '+ 'Institut de Pysique du Globe de Paris' +  '\n')
#    outdomesfile.write('E-mail                   : '+ '\n')
#    outdomesfile.write('Date                     : '+ '\n')
#    outdomesfile.write(' '+ '\n')
#    outdomesfile.write('2. Site Name                : ' + name + '\n')
#    outdomesfile.write('3. Country                  : ' + 'Guadeloupe , French West Indies' + '\n'  )
#    outdomesfile.write('4. Point Description        : ' + row['Type'] + '\n' )
#    outdomesfile.write('5. DOMES Number             : ' + 'N/A' + '\n' )
#    outdomesfile.write('6. Local Number             : ' + 'N/A' + '\n')
#    outdomesfile.write('7. 4-Char Code              : ' + row['Alias'] + '\n' )
#    outdomesfile.write(' '+ '\n')
#    outdomesfile.write('8. Approximate Position ' + '\n')
#    outdomesfile.write('   Latitude (deg min)       : '  + str(row['Lat. (WGS84)']) + '\n' )
#    outdomesfile.write('   Longitude (deg min)      : '  + str(row['Lon. (WGS84)']) + '\n' )
#    outdomesfile.write('   Elevation (m)            : '  + str(row['Elev. (m)']) + '\n')
#    outdomesfile.write(' '+ '\n')
#    outdomesfile.write('9. Instrument               : '  + str('GPS on ' + row['Type']) + '\n' )
#    outdomesfile.write('10. Date of Installation     : ' + str(row['Début / Installation']) + '\n')
#    outdomesfile.write(' '+ '\n')
#    outdomesfile.write('11. Operation Contact Name   : ' + '\n')
#    outdomesfile.write('    Agency                   : ' + '\n')
#    outdomesfile.write('    E-mail                   : ' + '\n')
#    outdomesfile.write(' ' + '\n')
#    outdomesfile.write('12. Site Contact Name        : ' + '\n')
#    outdomesfile.write('    Agency                   : ' + '\n')
#    outdomesfile.write('    E-mail                   : ' + '\n')


#  ________   ________ __  __ _____  _      ______  _____
# |  ____\ \ / /  ____|  \/  |  __ \| |    |  ____|/ ____|
# | |__   \ V /| |__  | \  / | |__) | |    | |__  | (___
# |  __|   > < |  __| | |\/| |  ___/| |    |  __|  \___ \
# | |____ / . \| |____| |  | | |    | |____| |____ ____) |
# |______/_/ \_\______|_|  |_|_|    |______|______|_____/
#
#
## CATS
#
#neulist = ['/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wABMF_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wADE0_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wASF0_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wCBE0_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wDHS0_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wDSD0_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wFNA0_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wFSDC_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wHOUE_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wLAM0_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wMGL0_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wPDB0_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wPSA1_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wSOUF_GINS_PS.neu',
#'/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wTDB0_GINS_PS.neu']
#
#statinfo_path = '/home/psakicki/THESE/CODES/stationinfo2gins/station.info.ovsg'
#
#for neu in neulist:
#    print neu
#    statinfo_2_cats(statinfo_path,neu)
#
#read_station_info_solo(statinfo_path,'ABMF')

#read_station_info_solo('/home/psakicki/THESE/CODES/stationinfo2gins/station.info.ovsg','HOUE')
#statname_of_catsfile('/home/psakicki/gin/TP/GWADA/DATA_CATS/NEU/wASF0_GINS_PS.neu')



# GINS
#station_info_2_gins('/home/psakicki/gin/TP/GWADA/station.info.ovsg.volobsis','/home/psakicki/gin/TP/GWADA/lfile_mod_mk2','/home/psakicki/gin/TP/GWADA/GWADA_mk2.stat')
#coordfile = '/home/psakicki/gin/TP/GWADA/RINX2/KARIB/pbo.snaps_igs08.vel'
#statinfo = '/home/psakicki/gin/TP/GWADA/RINX2/KARIB/station.info'

#specific = ['CN40' , 'CN19' , 'CN38' , 'SAN0' , 'CN35' , 'ROA0' , 'CN29' , 'CN18' , 'CBMD' , 'LCSB' , 'GCFS' , 'GCEA' , 'CN10' , 'CN12' , 'PRCG' , 'SMRT' , 'BARA']

#station_info_2_gins(statinfo,coordfile,'/home/psakicki/gin/TP/GWADA/KARIB_mk3.stat',specific_stats_lis=specific)
