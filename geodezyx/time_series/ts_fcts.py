# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 17:38:41 2019

@author: psakicki
"""

########## BEGIN IMPORT ##########
#### External modules
import copy
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import os 
import scipy

#### geodeZYX modules
from geodezyx import conv
from geodezyx import stats
from geodezyx import utils
from geodezyx import time_series
from geodezyx import reffram

#### Import star style
from geodezyx import *                   # Import the GeodeZYX modules
from geodezyx.externlib import *         # Import the external modules
from geodezyx.megalib.megalib import *   # Import the legacy modules names
##########  END IMPORT  ##########



def print4compar(dA,dB,dC,dD,coortype):

    if coortype == 'ENU':
        Astr,Bstr,Cstr = 'E','N','U'
    elif coortype == 'XYZ':
        Astr,Bstr,Cstr = 'X','Y','Z'
    else:
        Astr,Bstr,Cstr = 'A','B','C'

    print("moyenne aritm & RMS et std composante " + Astr)
    print(str(np.nanmean(dA)))
    print(stats.RMSmean(dA))
    print(str(np.nanstd(dA)))
    print('')

    print("moyenne aritm & RMS et std composante " + Bstr)
    print(str(np.nanmean(dB)))
    print(stats.RMSmean(dB))
    print(str(np.nanstd(dB)))
    print('')

    print("moyenne aritm & RMS et std composante " + Cstr)
    print(str(np.nanmean(dC)))
    print(stats.RMSmean(dC))
    print(str(np.nanstd(dC)))
    print('')

    print("moyenne aritm & RMS et std composante D")
    print(str(np.nanmean(dD)))
    print(stats.RMSmean(dD))
    print(str(np.nanstd(dD)))
    print('')

    print("RMS3D : sqrt((RMS_{}**2 + RMS_{}**2 + RMS_{}**2)/3 ) ".format(Astr,Bstr,Cstr))
    print(stats.RMSmean([stats.RMSmean(dA),stats.RMSmean(dB),stats.RMSmean(dC)]))
    print("RMS2D : uniquement sur les 2 composantes plani")
    print(stats.RMSmean([stats.RMSmean(dA),stats.RMSmean(dB)]))
    print('')

def print4compar_tabular(dicolist,split=0,print_2D3D_if_any=True):

    dAlist,dBlist,dClist,dDlist,dDblist,statlist = [],[],[],[],[],[]
    for dic in dicolist:
        if "dD2D" in list(dic.keys()):
            D2Dn3D = True
        else:
            D2Dn3D = False

        if not print_2D3D_if_any:
            D2Dn3D = False

        coortype = dic['coortype']
        Dtype = dic['Dtype']

        dAlist.append(dic['dA'])
        dBlist.append(dic['dB'])
        dClist.append(dic['dC'])

        if not D2Dn3D:
            dDlist.append(dic['dD'])
        else:
            dDlist.append(dic['dD2D'])
            dDblist.append(dic['dD3D'])

        statlist.append(dic['name'])

    if coortype == 'ENU':
        Astr,Bstr,Cstr = 'E','N','U'
    elif coortype == 'XYZ':
        Astr,Bstr,Cstr = 'X','Y','Z'
    else:
        Astr,Bstr,Cstr = 'A','B','C'

    datacol = ("Arithm. mean", "RMS mean", "std. dev.")
    datacol = ("Arit.", "RMS", utils.greek_alphabet()[17])

    headerline = ['Experiment']

    if not D2Dn3D:
        headervartup = (Astr,Bstr,Cstr,'D ' + Dtype)
    else:
        headervartup = (Astr,Bstr,Cstr,'D 2D', 'D 3D')

    for abc in headervartup:
        for dc in datacol:
            headerline.append(' '.join((dc,abc)))

#    headerline.append("RMS3D : sqrt((RMS_{}**2 + RMS_{}**2 + RMS_{}**2)/3 )".format(Astr,Bstr,Cstr))
#    headerline.append("RMS2D : uniquement sur les 2 composantes plani")

    headerline.append("RMS3D")
    headerline.append("RMS2D")

    LINES_STK = [headerline]

    if not D2Dn3D:
        dDblist = [None] * len(dDlist)

    for (dA,dB,dC,dD,dDb,stat) in zip(dAlist,dBlist,dClist,dDlist,dDblist,statlist):

        if np.isclose(np.sum(dA) + np.sum(dB) + np.sum(dC),0):
            continue

        line = []
        line.append(stat)

        line.append(np.nanmean(dA))
        line.append(stats.RMSmean(dA))
        line.append(np.nanstd(dA))

        line.append(np.nanmean(dB))
        line.append(stats.RMSmean(dB))
        line.append(np.nanstd(dB))

        line.append(np.nanmean(dC))
        line.append(stats.RMSmean(dC))
        line.append(np.nanstd(dC))

        line.append(np.nanmean(dD))
        line.append(stats.RMSmean(dD))
        line.append(np.nanstd(dD))

        if D2Dn3D:
            line.append(np.nanmean(dDb))
            line.append(stats.RMSmean(dDb))
            line.append(np.nanstd(dDb))

        line.append(stats.RMSmean([stats.RMSmean(dA),stats.RMSmean(dB),stats.RMSmean(dC)]))
        line.append(stats.RMSmean([stats.RMSmean(dA),stats.RMSmean(dB)]))

        LINES_STK.append(line)

    if split != 0:
        LINES_arr = np.vstack(LINES_STK)
        LINES1 = LINES_arr[:,:split]
        LINES2 = np.column_stack((LINES_arr[:,0],LINES_arr[:,split:]))
        return LINES1,LINES2
    else:
        return LINES_STK

def compar_plot(dico_list_in, namest = 0, namend    = 10  ,
                alpha = 0.8 , diapt = 1.5 , new_style = True,
                colormap = 'gnuplot'):


    if not new_style:
        fig , ax_arr = plt.subplots(2,2)
        fig.set_size_inches(16.53,11.69)
    else:
        fig , ax_arr = plt.subplots(4,1)
        fig.set_size_inches(11.69,16.53)

    cm = plt.get_cmap(colormap)
    NUM_COLORS = len(dico_list_in)
    color_list = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]


    for color , D in zip(color_list,dico_list_in):

        nametmp  = D['name']
        nametmp  = nametmp.ljust(namend - namest)
        coortype = D['coortype']
        Dtype    = D['Dtype']

        if coortype == 'ENU':
            Ati = '\u0394' + 'East component'
            Bti = '\u0394' + 'North component'
            Cti = '\u0394' + 'Up component'
            Dti = '\u0394' + 'Distance (' + Dtype + ') component'
            ylbl = 'meters'
        elif coortype == 'XYZ':
            Ati = 'delta X'
            Bti = 'delta Y'
            Cti = 'delta Z'
            Dti = 'delta D ' + Dtype
            ylbl = 'meters'
        elif coortype == 'FLH':
            Ati = 'delta Latitude'
            Bti = 'delta Longitude'
            Cti = 'delta Height'
            Dti = 'delta D ' + Dtype + ' (no sens)'
            ylbl = 'degrees'
        else:
            Ati = 'delta A'
            Bti = 'delta B'
            Cti = 'delta C'
            Dti = 'delta D ' + Dtype
            ylbl = '???'

        if not new_style:
            axA = ax_arr[0][0]
            axB = ax_arr[0][1]
            axC = ax_arr[1][0]
            axD = ax_arr[1][1]
        else:
            axA = ax_arr[0]
            axB = ax_arr[1]
            axC = ax_arr[2]
            axD = ax_arr[3]

        axA.plot(D['TA'],D['dA'], '.-' ,label=nametmp[namest:namend] ,
                 markersize=diapt, alpha=alpha,color=color)
        axA.set_ylabel(ylbl)
        axA.legend()
        axA.set_title(Ati)

        axB.plot(D['TB'] ,D['dB'] ,'.-',label=nametmp[namest:namend] ,
                 markersize=diapt, alpha=alpha,color=color)
        axB.set_ylabel(ylbl)
        axB.legend()
        axB.set_title(Bti)

        axC.plot(D['TC'] ,D['dC'],'.-',label=nametmp[namest:namend] ,
                 markersize=diapt, alpha=alpha,color=color)
        axC.set_ylabel(ylbl)
        axC.legend()
        axC.set_title(Cti)

        axD.plot(D['TD'] ,D['dD'],'.-',label=nametmp[namest:namend] ,
                 markersize=diapt, alpha=alpha,color=color)
        axD.set_ylabel(ylbl)
        axD.legend()
        axD.set_title(Dti)

    return fig

def compar(tstup , coortype='ENU' , seuil=3. , win=[] , mode='keep' ,
           Dtype='2D3D', namest=0,namend=10,alpha = 5 , diapt = 5 , verbose=True,
           print_report=True,plot=True):
    """
        si seuil == 0, pas de nettoyage

        On preconise compar en mode keep :
        Ainsi le Tref est strictement cantonné aux bonnes valeurs des Ti
        En mode del, les Ti sont extrapolé aux valeur de Tref => nan => pbs pour la mad


        if 2D3D /3D2D : selectionne 2D ou 3D (selon) en prioritaire pour dD mais
        détermine 2D et 3D et les introduit dans des clés adhoc dD2D dD3D
        '''
    """

    if len(Dtype) == 4:
        D2n3 = True
        Dtype = Dtype[0:2]
    else:
        D2n3 = False

    diapt=10
    alpha=0.8

    if len(tstup) <= 1:
        print("ERR : compar : len(tstup) <= 1 , on ne compar pas un élément seul ;) !!!!")

    if verbose:
        print("")
        print("BILAN DE LA COMPARAISON")
        print("reference : " +  tstup[0].name)
        print("coordonnées : " +  coortype)
        print("")


    dicolist = []

    for ivar,tsvar in enumerate(tstup):


        #dico de sortie
        dicovar = dict()

        if verbose:
            print("===========================================================")
            print(tsvar.name)
            print("------------------------------")


#        tsout = TimeSeriePoint()
#        tsout = copy.copy(tsvar)
#        tsout.del_data()

        if tsvar.bool_interp_uptodate == False:
            tsvar.interp_set()

        A, B, C, T, sA , sB , sC = tsvar.to_list(coortype=coortype)

        tsref = tstup[0].interp_get(T,coortype=coortype)
        Aref, Bref, Cref , Tref , sA , sB , sC = tsref.to_list(coortype=coortype)


        # WTF IS THOSE LINES
        # autowin = time_gap(tstup[0],marge=1.1,mode='keep')
        # win = autowin
        #
        # outbool = []


        # (winlis == [] or (ivar in winlis)) and
        if win != []:
            if verbose:
                print("Application du fenetrage")
                print("------------------------------")
            bl = time_win_T(T,win,mode=mode)
            if verbose:
                print("Nb de pts : " , len(bl))
                print("Pts validés :" , np.sum(bl))
                print('')

        else:
            bl = np.asarray([True] * len(Tref))

        dA = A[bl] - Aref[bl]
        dB = B[bl] - Bref[bl]
        dC = C[bl] - Cref[bl]
        Tref = Tref[bl]

        #print "N of diff close to 0 => Bug ..."
        #print 'dA' , np.sum(np.isclose(dA,0))
        #print 'dB' , np.sum(np.isclose(dB,0))
        #print 'dC' , np.sum(np.isclose(dC,0))


        dD2D = np.sqrt(dA ** 2 + dB ** 2)
        dD3D = np.sqrt(dA ** 2 + dB ** 2 + dC ** 2)

        if Dtype == '2D':
            dD = np.sqrt(dA ** 2 + dB ** 2)
        elif Dtype == '3D':
            dD = np.sqrt(dA ** 2 + dB ** 2 + dC ** 2)

        else:
            print("ERR : compar : wrong Dtype")


#        print "Stats sans Nettoyage"
#        print "------------------------------"
#        print4compar(dA,dB,dC,dD,coortype)

        # nettoyage par la MAD
        if seuil != 0:
            if verbose:
                print("Nettoyage par la MAD")
                print("------------------------------")
            dAin = dA
            dA,bb = stats.outlier_mad(dA,seuil=seuil)
            dAout = dA
            TA = conv.posix2dt(np.array(Tref[bb]))

            dBin = dB
            dB,bb = stats.outlier_mad(dB,seuil=seuil)
            dBout = dB
            TB = conv.posix2dt(np.array(Tref[bb]))

            dCin = dC
            dC,bb = stats.outlier_mad(dC,seuil=seuil)
            dCout = dC
            TC = conv.posix2dt(np.array(Tref[bb]))

            dDin = dD
            dD,bb = stats.outlier_mad(dD,seuil=seuil)
            TD = conv.posix2dt(np.array(Tref[bb]))
            dDout = dD
            if D2n3:
                dD2Din = dD2D
                dD2D,bb = stats.outlier_mad(dD2D,seuil=seuil)
                TD2D = conv.posix2dt(np.array(Tref[bb]))
                dD2Dout = dD2D
                dD3Din = dD3D
                dD3D,bb = stats.outlier_mad(dD3D,seuil=seuil)
                TD3D = conv.posix2dt(np.array(Tref[bb]))
                dD3Dout = dD3D

            if verbose:
                print('')
                print("Stats aprÃ¨s Nettoyage")
                print("------------------------------")
            if print_report:
                print4compar(dA,dB,dC,dD,coortype)

        else:
            Tdt = conv.posix2dt(Tref)
            TA,TB,TC,TD = Tdt,Tdt,Tdt,Tdt
            dAin,dBin,dCin,dDin = dA,dB,dC,dD

        # pour la legende : les indices de debut et fin du nom

        # fabrication du Dico de sortie
        dicovar['name']     = tsvar.name
        dicovar['coortype'] = coortype
        dicovar['Dtype']    = Dtype

        dicovar['TA'] = np.array(TA)
        dicovar['dA'] = np.array(dA)
        dicovar['dAbrut'] = np.array(dAin)
        dicovar['TB'] = np.array(TB)
        dicovar['dB'] = np.array(dB)
        dicovar['dBbrut'] = np.array(dBin)
        dicovar['TC'] = np.array(TC)
        dicovar['dC'] = np.array(dC)
        dicovar['dCbrut'] = np.array(dCin)
        dicovar['TD'] = np.array(TD)
        dicovar['dD'] = np.array(dD)
        dicovar['dDbrut'] = np.array(dDin)
        dicovar['dDtype'] = np.array(Dtype)
        if D2n3:
            dicovar['TD2D'] = np.array(TD2D)
            dicovar['dD2D'] = np.array(dD2D)
            dicovar['dD2Dbrut'] = np.array(dD2Din)
            dicovar['dD2Dtype'] = np.array('2D')

            dicovar['TD3D'] = np.array(TD3D)
            dicovar['dD3D'] = np.array(dD2D)
            dicovar['dD3Dbrut'] = np.array(dD3Din)
            dicovar['dD3Dtype'] = np.array('3D')

        dicolist.append(dicovar)

    if plot:
        compar_plot(dicolist,namest,namend,alpha,diapt)
        try:
            suptit = ' '.join((tstup[-1].stat , str(tstup[-1].startdate()) ,
                               str(tstup[-1].enddate())))
            plt.suptitle(suptit)
        except:
            pass

    return dicolist

def compar2(tstup,coortype='ENU',seuil=3.,win=[],mode='keep',
           Dtype='3D',namest=0,namend=10,alpha = 5 , diapt = 5 ,verbose=True,
           print_report=True,plot=True):
    """
    160903 cette fonction semble discontinuée
    """

    tsref = tstup[0]
    dicolist = []

    for tsvar in tstup:
        Aref,Bref,Cref,Tref,_,_,_  = tsref.to_list()
        Avar,Bvar,Cvar,Tvar,_,_,_  = tsvar.to_list()

        Tref = np.round(Tref,1)
        Tvar = np.round(Tvar,1)

        Tcommon    = np.intersect1d(Tref,Tvar)
        Tcommon_dt = conv.posix2dt(Tcommon)

        ind_ref = np.nonzero(np.in1d(Tref , Tcommon))[0]
        ind_var = np.nonzero(np.in1d(Tvar , Tcommon))[0]

        dA = Avar[ind_var] - Aref[ind_ref]
        dB = Bvar[ind_var] - Bref[ind_ref]
        dC = Cvar[ind_var] - Cref[ind_ref]

        if Dtype == '2D':
            dD = np.sqrt(dA ** 2 + dB ** 2)
        elif Dtype == '3D':
            dD = np.sqrt(dA ** 2 + dB ** 2 + dC ** 2)

        dicovar = dict()

        dicovar['name']     = tsvar.name
        dicovar['coortype'] = coortype
        dicovar['Dtype']    = Dtype

        dicovar['TA'] = Tcommon_dt
        dicovar['dA'] = np.array(dA)
        dicovar['dAbrut'] = np.array(dA)
        dicovar['TB'] = Tcommon_dt
        dicovar['dB'] = np.array(dB)
        dicovar['dBbrut'] = np.array(dB)
        dicovar['TC'] = Tcommon_dt
        dicovar['dC'] = np.array(dC)
        dicovar['dCbrut'] = np.array(dC)
        dicovar['TD'] = Tcommon_dt
        dicovar['dD'] = np.array(dD)
        dicovar['dDbrut'] = np.array(dD)
        dicovar['dDtype'] = np.array(Dtype)

        dicolist.append(dicovar)

    if plot:
        compar_plot(dicolist,namest,namend,alpha,diapt)
        try:
            suptit = ' '.join((tstup[-1].stat , str(tstup[-1].startdate()) ,
                               str(tstup[-1].enddate())))
            plt.suptitle(suptit)
        except:
            pass

    return dicolist


def mad_cleaner(tsin,seuil=3.5,method='dist',coortype='ABC',
                detrend_first=False,output_detrended=False, verbose=False):
    '''
    method : methode d'élimination :
            dist : on élimine les point qu sont trop loin en distance de la posi de ref
            indep : on traite les point independaments

            dist est a privilégier
    output_detrended ne marche que si detrend_first est activé
    '''

    if coortype == 'ABC':
        coortype = tsin.initype()
    print(coortype)

    if detrend_first:
        tswork = detrend_ts(tsin,coortype)
    else:
        tswork = tsin

    A,B,C,T,sA,sB,sC = tswork.to_list(coortype)

    D = np.sqrt(A ** 2 + B ** 2 + C ** 2)

    if method == 'dist':
        Dout , bb = stats.outlier_mad(D,seuil)
    if method == 'indep':
        Aout , bbA = stats.outlier_mad(A,seuil)
        Bout , bbB = stats.outlier_mad(B,seuil)
        Cout , bbC = stats.outlier_mad(C,seuil)

        bb = bbA * bbB * bbC

    if output_detrended:
        tsout = bool_cleaner(tswork,bb)
    else:
        tsout = bool_cleaner(tsin,bb)

    if verbose:
        killratio = np.round(float(tsout.nbpts) / float(tsin.nbpts),4)

        print("INFO : mad_cleaner :" , tsin.stat , "clean ratio" , tsin.nbpts , tsout.nbpts , killratio)

    return tsout

def sigma_cleaner(tsin,seuil=3,coortype='ABC',cleantype='any', verbose=False):

    if coortype == 'ABC':
        coortype = tsin.initype()

    A,B,C,T,sA,sB,sC = tsin.to_list(coortype)

    sAout , bbsA = stats.outiler_sigma(sA,seuil)
    sBout , bbsB = stats.outiler_sigma(sB,seuil)
    sCout , bbsC = stats.outiler_sigma(sC,seuil)

    if cleantype == 'any':
        boolbad = bbsA * bbsB * bbsC
    elif cleantype == 'all':
        boolbad = bbsA + bbsB + bbsC
    tsout = bool_cleaner(tsin,boolbad)

    if verbose:
        killratio = np.round(float(tsout.nbpts) / float(tsin.nbpts),4)

        print("INFO : sigma_cleaner :" , tsin.stat , "clean ratio" , tsin.nbpts , tsout.nbpts , killratio)

    return tsout

def std_dev_cleaner(tsin,stddev_threshold,coortype="ABC",cleantype="any",
                    verbose=False):
    """
    A rebooted (1807) version of sigma_cleaner
    just remove values in a timeserie with a high sigma/std deviation
    """

    if coortype == 'ABC':
        coortype = tsin.initype()

    A,B,C,T,sA,sB,sC = tsin.to_list(coortype)

    bbsA = sA <= stddev_threshold
    bbsB = sB <= stddev_threshold
    bbsC = sC <= stddev_threshold

    if cleantype == 'any':
        boolbad = bbsA * bbsB * bbsC
    elif cleantype == 'all':
        boolbad = bbsA + bbsB + bbsC
    tsout = bool_cleaner(tsin,boolbad)

    if verbose:
        killratio = np.round(float(tsout.nbpts) / float(tsin.nbpts),4)

        print("INFO : std_dev_cleaner :" , tsin.stat , "clean ratio" , tsin.nbpts , tsout.nbpts , killratio)

    return tsout

def bool_cleaner(tsin,boollist, verbose=False):
    ''' A partir d'une liste de bool de meme longeur que le nbre de points
        on ne conserve que les points True '''
    tsout = copy.copy(tsin)
    ptslistin = tsin.pts
    if len(ptslistin) != len(boollist):
        print("ERR : bool_cleaner : len bool != len pts")
        return 0
    ptslistout = []
    for pt,b in zip(ptslistin,boollist):
         if b:
             ptslistout.append(pt)
    tsout.pts = ptslistout

    if verbose:
        killratio = np.round(float(tsout.nbpts) / float(tsin.nbpts),4)

        print("INFO : bool_cleaner :" , tsin.stat , "clean ratio" , tsin.nbpts , tsout.nbpts , killratio)

    return tsout

def linear_regress_find_coeff(tsin,coortype='ENU'):
    A, B, C , T , sA , sB , sC = tsin.to_list(coortype)
    outW = []
    for i,composante in enumerate((A,B,C)):
        M = np.array([ T , np.ones(len(T)) ])
        w  = scipy.linalg.lstsq(M.T,composante)[0]
        outW.append((w[0],w[1]))
    return outW

def detrend_ts(tsin,coortype='ENU'):
    A, B, C , T , sA , sB , sC = tsin.to_list(coortype)
    try:
        Xref , Yref , Zref = tsin.mean_posi()
    except:
        Xref , Yref , Zref = 0,0,0

    W = linear_regress_find_coeff(tsin,coortype)
    UnTr = [] # UnTr
    for i,composante in enumerate((A,B,C)):
        w = W[i]
        comp_line = w[0]*T+w[1]
        comp_untrend = composante - comp_line
        UnTr.append(comp_untrend)

    tsout = ts_from_list(UnTr[0],UnTr[1],UnTr[2],T,coortype,sA,sB,sC,tsin.stat,tsin.name)




    tsout.anex['Xref'] = Xref
    tsout.anex['Yref'] = Yref
    tsout.anex['Zref'] = Zref

    return tsout

def linear_regress_ts(tsin,coortype='ENU',titledetails = ''):
    """ doit être cablé ASAP linear_regress_find_coeff"""
    try:
        A, B, C , T , sA , sB , sC = tsin.to_list(coortype)
        fig , axes = plt.subplots(2,2)
        fig.suptitle(tsin.stat + ' ' + titledetails)
        axes = np.reshape(axes,-1)

        W = linear_regress_find_coeff(tsin,coortype)

        for i,composante in enumerate((A,B,C)):
            M = np.array([ T , np.ones(len(T)) ])
#            w  = scipy.linalg.lstsq(M.T,composante)[0]
            w = W[i]
            xi = np.arange(min(T),max(T),10**5)
            line = w[0]*xi+w[1]
            axes[i].plot(conv.posix2dt(xi),line,'r-')
            axes[i].plot(conv.posix2dt(T),composante,'+')
            axes[i].set_title(coortype[i])
            axes[i].set_xlabel('Date')
            axes[i].set_ylabel('Displacement (m)')
            v = w[0] * 31556926 / 10**-3
            v2 = round(v,2)
            print('composante ' + coortype[i] + ' = ' + str(v) + 'mm/year')
            axes[i].text(0.9, 0.1, str(v2) + 'mm/year' , horizontalalignment='center',verticalalignment='center',transform=axes[i].transAxes)
    except:
        pass

def linear_regress_ts_discont(tsin,coortype = 'ENU') :
    discont_improved = tsin.discont + [dt.datetime(2099,1,1)]
    for j in range(len(discont_improved)-1):
        print('période ',j,':',discont_improved[j],discont_improved[j+1])
        tsbis = time_win(tsin,[[discont_improved[j],discont_improved[j+1]]])
        datetitle = ' : ' + str(discont_improved[j]) + "\u2192" + str(discont_improved[j+1])
        linear_regress_ts(tsbis,coortype,titledetails=datetitle)

#def export_ts_as_pbo_pos(tsin):
#    # IN PROGRESS
#    tswork = copy.deepcopy(tsin)
#    outtfile.write('{:.5f}   {:+.6f}    {:+.6f}    {:+.6f} {:+.6f} {:+.6f} {:+.6f}\n'.format(t,n-n0,e-e0,u-u0,se,sn,su))
#

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
        print("INFO : plot exported in " , export_file)
    if close_fig_after_export:
        plt.close(f)

    # plot B avec une droite de regression
    if plot_B:
        linear_regress_ts(tsin)
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
        print("INFO : plot exported in " , export_file)
        if close_fig_after_export:
            plt.close(f)
    return None


def export_ts(ts,outdir,coordtype = 'ENU',outprefix='',write_header=False):
    """
    export the timeserie

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

    print('INFO : timeserie exported in ' + outpath)
    return None


def export_ts_as_neu(tsin,outdir,outprefix,coordtype = 'ENU'):
    """
    export to a HECTOR .neu compatible format

    outfile will be writed in
    /outdir/outprefixSTAT.neu
    
    NB: The XYZ mode is quite dirty (191001)
    """
    if not hasattr(tsin[0],'X'):
        print('WARN : export_ts_as_neu : no XYZ in ts')
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
    outfile.write('# Reference epoch: {}\n'.format(conv.toYearFraction(first_pt.Tdt)))
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
            outfile.write('# offset {} 7\n'.format(conv.toYearFraction(disc)))
        outfile.write('#\n')
    # write the data
    for e,n,u,t,se,sn,su in zip(E,N,U,T,sE,sN,sU):
        t = conv.toYearFraction(conv.posix2dt(t))
        if coordtype == "ENU":        
            outfile.write('{:.5f}   {:+.6f}    {:+.6f}    {:+.6f} {:+.6f} {:+.6f} {:+.6f}\n'.format(t,n-n0,e-e0,u-u0,se,sn,su))
        elif coordtype == "XYZ":
            outfile.write('{:.5f}   {:+.6f}    {:+.6f}    {:+.6f} {:+.6f} {:+.6f} {:+.6f}\n'.format(t,n-n0,e-e0,u-u0,se,sn,su))

    print('INFO : timeserie exported in ' + outpath)
    return None


def export_ts_as_hector_enu(tsin,outdir,outprefix,coordtype = 'ENU'):
    """
    export to a HECTOR .enu (and not .neu !) compatible format
    This format is simpler : just gives MJD E N U

    This format is necessary to force a sampling period.

    outfile will be writed in
    /outdir/outprefixSTAT.enu
    """

    print("NOT IMPLEMENTED YET !")

    return None



def export_ts_as_midas_tenu(tsin,outdir,outprefix,coordtype = 'ENU',
                            export_step=True):
    """
    export to a MIDAS .tneu compatible format

    outfile will be writed in
    /outdir/outprefixSTAT.tneu

    if export_step == True:
    export a step file as
    /outdir/outprefixSTAT.step
    """
    if not hasattr(tsin[0],'X'):
        print('WARN : export_ts_as_midas_tneu : no XYZ in ts')
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
        t = conv.toYearFraction(conv.posix2dt(t))
        #outfile.write('{} {:.5f} {:+.6f} {:+.6f} {:+.6f} \n'.format(stat,t,n-n0,e-e0,u-u0))
        outfile.write('{} {:.5f} {:+.6f} {:+.6f} {:+.6f} \n'.format(stat,t,e-e0,n-n0,u-u0))

    print('INFO : timeserie exported in ' + outpath)

    if export_step and tswork.bool_discont:
        outpath_step = outdir +'/' + outprefix + tswork.stat + '.step'
        outfile_step = open(outpath_step,'w+')
        for d in tswork.discont:
            d = conv.toYearFraction(d)
            line = tswork.stat + " " + str(d) + "\n"
            outfile_step.write(line)

            print('INFO : timeserie discont. (steps) exported in ' + outpath_step)

    return None




def decimate_cleaner(tsin,minval,in_place = False):
    """
    in_place DOES'NT WORK !!!
    """
    if not in_place:
        tsout = copy.deepcopy(tsin)
    else:
        tsout = tsin
    if minval == 0:
        return tsout

    tsout.del_data()
    for pt in tsin.pts:
        #print pt.T ,minval , np.mod(pt.T,minval)
        if np.mod(np.round(pt.T,2),minval) < 10**-3:
            tsout.add_point(pt)
    tsout.interval_nominal()
    return tsout

def decimate_cleaner_2(tsin,N,in_place = False):
    """ keep a value every N vals """
    if in_place:
        tsout = copy.deepcopy(tsin)
    else:
        tsout = tsin
    tsout.del_data()
    tsout.pts = tsin.pts[:-1:N]
    tsout.interval_nominal()
    return tsout


def mean_list_of_pts(ptslisin):
    """ useful for merge fct
    ONLY IMPLEMENTED FOR ENU coords for the moment """
    E = np.nanmean([p.E for p in ptslisin])
    N = np.nanmean([p.N for p in ptslisin])
    U = np.nanmean([p.U for p in ptslisin])
    sE = np.nanmean([p.sE for p in ptslisin])
    sN = np.nanmean([p.sN for p in ptslisin])
    sU = np.nanmean([p.sU for p in ptslisin])

    Tlis = [p.T for p in ptslisin]
    T = ((np.max(Tlis) - np.min(Tlis)) / 2.) + np.min(Tlis)
    pt = time_series.Point(E,N,U,T,'ENU',sE,sN,sU,ptslisin[0].name)

    return pt


def merge(tsin,N):
    """ merge N points in one """
    tsin.sort()
    tsout = copy.deepcopy(tsin)
    tsout.del_data()
    slic = utils.sliceIt(tsin.pts,N)
    for sl in slic:
        pt = mean_list_of_pts(sl)
        tsout.add_point(pt)
    tsout.interval_nominal()
    return tsout

def merge_ts(ts_list_in):
    pts_list_merged = []
    for ts in ts_list_in:
        pts_list_merged = pts_list_merged + ts.pts

    ts_out = time_series.TimeSeriePoint()
    ts_out.pts = pts_list_merged
    ts_out.anex = ts_list_in[0].anex
    ts_out.interval_nominal()
    ts_out.sort()
    return ts_out


def time_win(tsin,windows,mode='keep',outbool=False):
    if type(windows[0][0]) == dt.datetime:
        Ttype = 'Tdt'
    else:
        Ttype = 'T'

    tsout = copy.copy(tsin)
    tsout.del_data()

    if isinstance(tsout,time_series.TimeSeriePoint):
        data = 'pts'
    elif isinstance(tsout,time_series.TimeSerieObs):
        data = 'obs'
    else:
        print('BUG : time_win : ca chie')

    outboollist = []

    print(tsin)
    for i,e in enumerate(windows):
        print('windows', i)
        print(e[0] , e[1])


    for pt in getattr(tsin,data):
        boollist = []
        for win in windows:
            if win[1] < win[0]:
                print("WARN : time_win : le début de periode est + grand que la fin")
                print(win[0] , win[1])
                print(mode)

            if win[0] < getattr(pt,Ttype) and win[1] > getattr(pt,Ttype):
                boollist.append(True)
            else:
                boollist.append(False)

        if mode == 'keep' :
            finalbool = np.any(boollist)
        elif mode == 'del':
            finalbool = not np.any(boollist)

        if finalbool:
            tsout.add_point(pt)

        outboollist.append(finalbool)

    tsout.bool_interp_uptodate = False

    #tsout.interp_set()
    if outbool:
        outboollist = np.array(outboollist)
        return tsout,outboollist
    else:
        return tsout

def time_win_T(Tin,win,mode='del'):

    boolfinalist = []

    for t in Tin:
        booltemp = []
        for w in win:
            startw = w[0]
            endw = w[1]

            # CONVENTION : VRAI SI DANS L'INTERVALLE
            # QU'IMPORTE SI KEEP OU DEL
            if startw < t and t < endw :
                booltemp.append(True)
            else:
                booltemp.append(False)

        if mode == 'keep' :
            finalbool = np.any(booltemp)
        elif mode == 'del':
            finalbool = not np.any(booltemp)

        boolfinalist.append(finalbool)

    return np.array(boolfinalist)

def time_gap(tsin,marge=2,mode='del'):
    ''' ENTREE une TimeSerie
        SORTIE une window (liste de listes)'''

    inomi = tsin.i_nomi

    T = tsin.to_list()[3]

    Tdiff = np.diff(np.sort(T))

    boollist = []
    # Dans boollist :
    # False si c'est gap
    # True sinon

    for td in Tdiff:
        if td > inomi * marge:
            boollist.append(False)
        else:
            boollist.append(True)

    # Difference suivant le mode
    if mode == 'keep' :
        modbool = True
    elif mode == 'del':
        modbool = False

    indbool = np.where(np.array(boollist) == modbool)[0]
    protowin = reffram.group_consecutives(indbool)

    winout = []

    for e in protowin:
        if len(e) == 1:
            eout = [T[e[0]],T[e[0]+1]]
        elif len(e) == 2:
            eout = [T[e[0]],T[e[1]+1]]
        else:
            "BUG : time_gap"
        winout.append(eout)

    return winout

def compar_elts_in_ts(ts1,ts2):
    ''' ts2 must contains less elts than ts1 (ts2 = cleaned one)'''
    T1 = ts1.to_list()[3]
    T2 = ts1.to_list()[3]

    for t2 in T2:
        T1.remove(t2)

    return T1


def rotate_pt_cls_solo(tsattin,pointin,Rtype='R1', xyzreftuple = ([1, 0, 0], [0, 1, 0], [0, 0, 1]),angtype='deg'):
    ''' ENTREE : tsattin : une TS d'attitude  (N angles)
                 pointin : UN Point en entrée

        SORTIE : une TSpoint de N points '''

    attlisttmp = tsattin.to_list()
    R = attlisttmp[0]
    P = attlisttmp[1]
    Y = attlisttmp[2]
    T = attlisttmp[3]

    A,B,C,_,_ = pointin()

    ptrotlis = XXX.rotate_points(R,P,Y,[np.array([A,B,C])],Rtype,xyzreftuple,angtype)

    ptrotlis = utils.shrink_listoflist(ptrotlis)

    tsout = time_series.TimeSeriePoint()

    for p,t in zip(ptrotlis,T):
        #print p
        tsout.add_point(time_series.Point(p[0],p[1],p[2],t,initype=pointin.initype))

    return tsout


def rotate_points_class(tsattin,ptslin, Rtype='R1', xyzreftuple = ([1, 0, 0], [0, 1, 0], [0, 0, 1]),angtype='deg'):

    listsout = []

    for pt in ptslin:
        ts = rotate_pt_cls_solo(tsattin,pt,Rtype,xyzreftuple,angtype)
        listsout.append(ts)

    return listsout

def add_offset_point(ptin,dA,dB,dC,coortype='ENU'):
    """
    ONLY IMPLEMENTED FOR ENU FOR THE MOMENT
    150415 : remark still necessary ???

    coortype == 'UXYZ' :
        specific case where we correct an Up offset directly in the XYZ coords
        very usefull for an antenna offset in for a moving GPS
        (but works only for the up)
    """
    ptout = copy.copy(ptin)
    if coortype == 'ENU':
        ptout.ENUset(ptin.E + dA,ptin.N + dB,ptin.U + dC,
                     ptin.sE,ptin.sN,ptin.sU)
    elif coortype == 'XYZ':
        ptout.XYZset(ptin.X + dA,ptin.Y + dB,ptin.Z + dC,
                     ptin.sX,ptin.sY,ptin.sZ)
    elif coortype == 'FLH':
        ptout.FLHset(ptin.F + dA,ptin.L + dB,ptin.H + dC,
                     ptin.sF,ptin.sL,ptin.sH)

    elif coortype == 'UXYZ':
        N = conv.normal_vector(ptin.F , ptin.L , ptin.H)
        dA2 , dB2 , dC2 = N * dC
        ptout.XYZset(ptin.X + dA2,ptin.Y + dB2,ptin.Z + dC2,
                     ptin.sX,ptin.sY,ptin.sZ)

    return ptout

def add_offset_ts(tsin,dA,dB,dC,coortype='ENU'):
    """
    return a copy of the tsin
    (tsin won't be affected)

    coortype == 'UXYZ' :
        specific case where we correct an Up offset directly in the XYZ coords
        very usefull for an antenna offset in for a moving GPS
        (but works only for the up)
    """
    tsout = copy.deepcopy(tsin)
    tsout.del_data()
    newpts = [add_offset_point(pt,dA,dB,dC,coortype) for pt in tsin.pts]
    for pt in newpts:
        tsout.add_point(pt)
    return tsout

def add_offset_smart_for_GINS_kine(tsin,tslist_offset_3ple,list_windows,coortype='XYZ'):
    """
    tslist_offset_3ple : list of len N containing (dX,dY,dZ) offsets
    list_windows       : list of len N-1 containing dates of changes
    """

    if type(list_windows[0]) is int:
        list_windows = [conv.posix2dt(e) for e in list_windows]

    list_windows = [dt.datetime(1980,1,1)] + list_windows + [dt.datetime(2099,1,1)]
    tsout = copy.deepcopy(tsin)
    tsout.del_data()
    olddX ,olddY , olddZ = 12345 , 12345 , 12345
    for pt in tsin:
        for i in range(len(list_windows)-1):
            if list_windows[i] <= pt.Tdt < list_windows[i+1]:
                curdX , curdY , curdZ = tslist_offset_3ple[i][0],tslist_offset_3ple[i][1],tslist_offset_3ple[i][2]
                if curdX != olddX or curdY != olddY or curdZ != olddZ:
                    olddX ,olddY , olddZ = curdX , curdY , curdZ
                    print('INFO : change of time window / offset')
                ptout = add_offset_point(pt,curdX,curdY,curdZ,coortype='XYZ')
                tsout.add_point(ptout)
                continue
    return tsout


def find_pts_from_ts_with_time(tin,tstupin,tol=0.001):
    ptsout = []
    for ts in tstupin:
        pt, i = ts.find_point(tin,tol=tol)
        if np.isnan(i):
            print("WARN : find_pts_from_ts_with_time : no pt find")
        ptsout.append(pt)
    return ptsout

def dist_btwn_2pts(ptA,ptB,coortype='XYZ'):

        A = np.array([ptA.X,ptA.Y,ptA.Z])
        B = np.array([ptB.X,ptB.Y,ptB.Z])

        return conv.dist(A,B)

def dist_diff_btwn_2pts(ptA,ptB):
    #dérive la distance entre 2 points A et B
        A = np.array([ptA.X,ptA.Y,ptA.Z])
        B = np.array([ptB.X,ptB.Y,ptB.Z])

        dAB   = A-B
        dist  = np.linalg.norm(dAB)

        diffA =   dAB / dist
        diffB = - dAB / dist

        return diffA, diffB

def helmert_trans(tsin,params='itrf2008_2_etrf2000', invert=False):
    tsout = copy.deepcopy(tsin)
    [pt.helmert_trans(params,invert) for pt in tsout.pts]
    return tsout


def velocity_trans(tsin , vx, vy, vz, epoc_init = 'auto' , epoc_end = 'auto'):
    tsout = copy.deepcopy(tsin)
    [pt.velocity_trans(vx, vy, vz, epoc_init , epoc_end) for pt in tsout.pts]
    return tsout


def interpolator_light(T,X,Y,Z):
    from scipy.interpolate import interp1d

    IX = interp1d(T,X,bounds_error=False,fill_value=np.nan)
    IY = interp1d(T,Y,bounds_error=False,fill_value=np.nan)
    IZ = interp1d(T,Z,bounds_error=False,fill_value=np.nan)

    return IX , IY , IZ

def interpolator_with_extrapolated(T,X,Y,Z):
    from scipy.interpolate import interp1d

    IX = interp1d(T,X,bounds_error=False,fill_value="extrapolate")
    IY = interp1d(T,Y,bounds_error=False,fill_value="extrapolate")
    IZ = interp1d(T,Z,bounds_error=False,fill_value="extrapolate")

    return IX , IY , IZ

def mean_posi_multi(tstup):
    mp = tstup[0].mean_posi()
    for ts in tstup:
        ts.ENUcalc(mp)
    return None

def refENU_for_tslist(tslist_in, tsref_marker = 0):
    """
    tsref_marker : indice of the reference time serie OR the 'all' keyword
    in this case all the time series mean position will be averaged
    """
    if tsref_marker == 'all':
        refENU_stk =  []
        for ts in tslist_in:
            refENUtmp = ts.mean_posi()
            refENU_stk.append(refENUtmp)

        X = np.mean([renu.X for renu in refENU_stk])
        Y = np.mean([renu.Y for renu in refENU_stk])
        Z = np.mean([renu.Z for renu in refENU_stk])

        refENU = time_series.Point(X,Y,Z,0,initype='XYZ')

    else:
        refENU = tslist_in[tsref_marker].mean_posi()


    for ts in tslist_in:
        ts.ENUcalc(refENU)

    return refENU

def time_win_multi(inplis):
    strt = inplis[0].startdate()
    end  = inplis[0].enddate()
    tsoutlis = [inplis[0]]
    for ts in list(inplis)[1:]:
        tsout = time_win(ts,[[strt,end]])
        tsoutlis.append(tsout)
    return tsoutlis

def ts_from_list(A,B,C,T,initype,sA=[],sB=[],sC=[],stat='STAT',name='NoName'):
    tsout = time_series.TimeSeriePoint()
    if len(sA) == 0:
        for a,b,c,t in zip(A,B,C,T):
            pt = time_series.Point(a,b,c,t,initype)
            tsout.add_point(pt)
    else:
        for a,b,c,t,sa,sb,sc in zip(A,B,C,T,sA,sB,sC):
            pt = time_series.Point(a,b,c,t,initype,sa,sb,sc)
            tsout.add_point(pt)
    tsout.stat = stat
    tsout.name = name
    if initype == 'ENU':
        tsout.boolENU = True
    else:
        tsout.boolENU = False
    return tsout


