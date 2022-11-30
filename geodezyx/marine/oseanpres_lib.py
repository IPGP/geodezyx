#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 17:41:32 2022

@author: psakic
"""

from geodezyx import utils
import matplotlib.pyplot as plt
import numpy as np
import itertools
import scipy

#### Import the logger
import logging
log = logging.getLogger(__name__)



psi_per_meter=1.45038




def freq2temp(f):
    """
    frequence du capteur de temperature > temperature

    """
    U0 =  5.766353
    Y1 = -4025.183
    Y2 = -11970.76
    Y3 =  0.
    X = (1/f) * 10**6 # (1/30000) * 10**6
    U = X - U0
    Temp = Y1 * U + Y2 * U**2 + Y3 * U**3
    return Temp

def freq2U(f):
    U0 =  5.766353
    return 1E6/f-U0

def temp2freq(T,out='freq'):
    """
    temperature > frequence du capteur de temperature (output == "freq") 
    OU  coef U (output == "U") OU  periode du capteur U (output == "tau")
    """
    
    #### Handle T as an interable (array, list....)
    if utils.is_iterable(T):
        return np.array([temp2freq(t,out) for t in T])
    #### Handle T as a scalar
    else:
        U0 =  5.766353
        Y1 = -4025.183
        Y2 = -11970.76
        Y3 =  0.
        U = np.roots( [Y3,Y2,Y1,-T]) + U0
        #X = (1/f) * 10**6 # (1/30000) * 10**6
        #U = X - U0
        #Temp = Y1 * U + Y2 * U**2 + Y3 * U**3
        if ( out == 'freq'):
            return [ 10**6/u  for u in U  if ( (10**6/u< 176000) and (10**6/u > 168000 ))][0]
        elif ( out == 'tau'):
            return [ 10**6/u  for u in U  if ( (10**6/u< 176000) and (10**6/u > 168000 ))][0]**-1
        elif ( out == "U"):
            return [ u-U0  for u in U  if ( (u> 10**6/176000) and (u < 10**6/168000 ))][0]
        else:
            return None


def temp2pres_f0(T):
    """
    temperature > composante temp pour le capteur de pression
    (homogène à une frequence)
    """
    U0 =  5.766353
    T1 =  29.89307
    T2 =  0.337614
    T3 =  56.99511
    T4 = 157.6942
    T5 =   0.
    
    F=temp2freq(T)
    U= (1/F)*10**6 - U0
    
    T0 = T1 + T2*U + T3*U**2 + T4*U**3 + T5*U**4
    return 10**6/T0


def freq2pres_old(f,U,F0):
    """
    freqence capteur pression > pression
    OLD version
    """

    C1 = -22682.65
    C2 = -1143.743
    C3 =  70903.62

    D1 =  0.040903
    D2 =  0.0
    
    C = C1 + C2*U + C3*U**2
    D = D1 + D2*U 

    T  = (1/f) * 10**6
    T0 = (1/F0)* 10**6

    P = C * (1 - (T0**2)/(T**2)) * (1 - D*(1 - (T0**2)/(T**2)))
    return (P)



def freq2pres(fpin,tin,out="psi"):
    """
    frequence capteur pression > pression (output == "psi")
    OU profondeur (output == "meter")
    """
    C1 = -22682.65
    C2 = -1143.743
    C3 =  70903.62

    D1 =  0.040903
    D2 =  0.0
    
    U = temp2freq(tin,'U')
    
    C = C1 + C2*U + C3*U**2
    D = D1 + D2*U 

    T  = (1/fpin)              * 10**6
    T0 = (1/temp2pres_f0(tin)) * 10**6

    P = C * (1 - (T0**2)/(T**2)) * (1 - D*(1 - (T0**2)/(T**2)))
    
    if out == "psi":
        return (P)
    elif out == "meter":
        return P/psi_per_meter
    else:
        raise Exception("bad output format in freq2pres")
        
        

    
def pres2freq(pin,tin,inp='psi',
              return_optimize_object=False):
    """
    pression > frequence capteur pression
    """
    
    #### Handle pin as an interable (array, list....)
    if utils.is_iterable(tin):
        log.err("tin as to be a scalar")
        raise Exception
    
    if utils.is_iterable(pin):
        return np.array([pres2freq(pin_i,tin,inp,return_optimize_object) for pin_i in pin])
    
    #### Handle pin as a scalar
    else:
        if inp == "psi":
            puse = pin
        elif inp == "meter":
            puse = pin * psi_per_meter
        
        
        def wrap_zero(fp):
            return freq2pres(fp, tin) - puse
            
        Opti = scipy.optimize.root(wrap_zero,30000)
        
        if return_optimize_object:
            return Opti
        else:
            return Opti.x[0]
        
        

def freq2counter(freq_or_tau_sensor,
                 count_sensor,
                 freq_or_tau_clk,
                 inp="freq",
                 round_fct=np.floor):
    
    if inp=="freq":
        tau_sensor = 1/freq_or_tau_sensor
        tau_clk = 1/freq_or_tau_clk
    else:
        tau_sensor = freq_or_tau_sensor
        tau_clk = freq_or_tau_clk
        
    N = (count_sensor * tau_sensor)/tau_clk
    
    if round_fct:
        return round_fct(N)
    else:
        return N  


def counter2freq(N_counted_by_clk,
                 count_sensor,
                 freq_clk,
                 out="freq"):
    
    
    tau_sensor = N_counted_by_clk/(freq_clk*count_sensor)
    
    if out=="freq":
        return 1/tau_sensor
    else:
        return tau_sensor
        
        
def pres_resolution(val_presin,
                    val_tempin,
                    val_clkin,
                    count_presin,
                    count_tempin,
                    input_val="freq",
                    output_val="meter",
                    relative_delta=True):
    
    if input_val == "freq": 
        if np.any(val_tempin < 1) or np.any(val_presin < 1) or np.any(val_clkin < 1):
            log.warn("some vals are < 1 while input is freq. Please double check")
        tau_temp = 1/val_tempin
        tau_pres = 1/val_presin
        tau_clk  = 1/val_clkin
    elif input_val == "tau":
        if np.any(val_tempin > 1) or np.any(val_presin > 1) or np.any(val_clkin > 1):
            log.warn("some vals are > 1 while input is tau. Please double check")
        tau_temp = val_tempin
        tau_pres = val_presin
        tau_clk  = val_clkin
        
    # internal values are tau, because at least it simplifies the computation 
    # of pres_bias  and avoids the 1/(a+b) = 1/a + 1/b mistake
        
    freq_clk = 1/tau_clk
    
    temp_central = freq2temp(1/tau_temp)
    pres_central = freq2pres(1/tau_pres,temp_central)
    
    tau_pres_bias =  2 / (freq_clk * count_presin)
    tau_temp_bias =  2 / (freq_clk * count_tempin)
        
    PresResArr = np.zeros((3,3))
    
    if relative_delta:
        divider = pres_central
    else:
        divider = 1.   
    
    for i,j in [(1,0),  (1,-1), 
                (0,-1), (-1,-1), 
                (-1,0), (-1,1), 
                (1,1),  (0,1)]:
        
        pres_biased = freq2pres(1/(tau_pres + i*tau_pres_bias),
                                freq2temp(1/(tau_temp + j*tau_temp_bias)))
                        
        pres_res = np.abs(pres_biased - pres_central)/divider
        
        PresResArr[i+1,j+1] = pres_res
    
    if output_val == "psi":
        pass
    elif output_val == "meter":
        pres_central=pres_central/psi_per_meter
        if not relative_delta:
            PresResArr=PresResArr/psi_per_meter ### correction only for absolute values !!!
    else:
        raise Exception("bad output format in pres_resolution")
        
    return pres_central,temp_central,PresResArr.max(),PresResArr


def resolution_grid_compute(Tau_pres,Tau_temp,fe,
                            count_pres,count_temp,
                            output_val="meter",
                            relative_delta=False):

    PresVals = np.zeros((len(Tau_pres),len(Tau_temp)))
    TempVals = np.zeros((len(Tau_pres),len(Tau_temp)))
    ResVals  = np.zeros((len(Tau_pres),len(Tau_temp)))
    FullResVals = np.zeros((len(Tau_pres),len(Tau_temp),3,3))
    
    for itaupres,itautemp in itertools.product(range(len(Tau_pres)),
                                               range(len(Tau_temp))):
        
        taupres,tautemp = Tau_pres[itaupres],Tau_temp[itautemp]
        
        OutRes = pres_resolution(taupres, 
                                 tautemp,
                                 1/fe, 
                                 count_pres, 
                                 count_temp,
                                 input_val="tau",
                                 output_val=output_val,
                                 relative_delta=relative_delta)
    
        PresVals[itaupres,itautemp] = OutRes[0]    
        TempVals[itaupres,itautemp] = OutRes[1]    
        ResVals[itaupres,itautemp]  = OutRes[2]
        FullResVals[itaupres,itautemp,:,:]  = OutRes[3]

    return PresVals,TempVals,ResVals,FullResVals
    

############ PLOT OF THE RESOLUTION VALUES
def resolution_plot_as_gradient_grid(PresVals,
                                     TempVals,
                                     ResVals,
                                     temp_ref,
                                     freq_input=False,
                                     coef_res=1000,
                                     relative_delta=False):
    
    #### coef_res = 10**6 if PPM
    #### coef_res = 10**3 if error in meters
    
    #################### INITIALISATION ##########################
    
    secx_pres2freq = lambda x:pres2freq(x, temp_ref,'meter')
    secx_freq2pres = lambda x:freq2pres(x, temp_ref,'meter')
    secy_temp2freq = lambda x:temp2freq(x,'freq')
    secy_freq2temp = lambda x:freq2temp(x)
    
    if not freq_input:
        secx_fcts = (secx_pres2freq, secx_freq2pres)
        secy_fcts = (secy_temp2freq, secy_freq2temp)
        xlab = 'Depth (m)'
        ylab = 'Temp (°C)'
        secxlab = 'Pres Sensor Freq (Hz)'
        secylab = 'Temp Sensor Freq (Hz)'
    else:
        secx_fcts = (secx_freq2pres, secx_pres2freq)
        secy_fcts = (secy_freq2temp, secy_temp2freq)        
        secxlab = 'Depth (m)'
        secylab = 'Temp (°C)'
        xlab = 'Pres Sensor Freq (Hz)'
        ylab = 'Temp Sensor Freq (Hz)'
        
        
    if not relative_delta:
        ResVals = ResVals*10**3
        title='Error in millimeter max={:.4f}mm'.format(ResVals.max())
        format_colorbar="%.4f"
    else:
        ResVals = ResVals*10**6        
        title='PPM max={:.4f}'.format(ResVals.max())
        format_colorbar="%.4f"
        
    ### change values of ResVal

    #################### INITIALISATION ##########################
    
    
    fig,ax = plt.subplots()    
    ax.ticklabel_format(useOffset=False)
    cm = plt.cm.get_cmap('viridis')
    fig.subplots_adjust(right=0.75,top=.85)
    im1 = ax.contourf(PresVals,TempVals,ResVals,cmap=cm,levels=200)
    cbar_ax1 = fig.add_axes([0.88, 0.15, 0.04, 0.7])
    cbar1 = fig.colorbar(im1,cax=cbar_ax1,format=format_colorbar)

    secx = ax.secondary_xaxis('top', functions=secx_fcts)
    secy = ax.secondary_yaxis('right', functions=secy_fcts)
    secx.ticklabel_format(useOffset=False)
    secy.ticklabel_format(useOffset=False)

    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    secx.set_xlabel(secxlab)
    secy.set_ylabel(secylab)
    ax.set_title(title)
    
    return fig,ax
    
    
    

                          
                          
                          
