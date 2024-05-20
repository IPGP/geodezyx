# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.stats contains functions for least-squares
processing. 

it can be imported directly with:
from geodezyx import stats

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""

import itertools
#### Import the logger
import logging
import multiprocessing as mp

import matplotlib.pyplot as plt
########## BEGIN IMPORT ##########
#### External modules
import numpy as np
import scipy

#### geodeZYX modules
from geodezyx import utils

log = logging.getLogger(__name__)

##########  END IMPORT  ##########



#  _                    _      _____                  _    _ _   _ _
# | |                  | |    / ____|                | |  | | | (_) |
# | |     ___  __ _ ___| |_  | (___   __ _ _ __ ___  | |  | | |_ _| |___
# | |    / _ \/ _` / __| __|  \___ \ / _` | '__/ __| | |  | | __| | / __|
# | |___|  __/ (_| \__ \ |_   ____) | (_| | |  \__ \ | |__| | |_| | \__ \
# |______\___|\__,_|___/\__| |_____/ \__, |_|  |___/  \____/ \__|_|_|___/
#                                       | |
#                                       |_|

def get_accur_coeff(i):
    """
    accuracy coefficients given by
    https://en.wikipedia.org/wiki/Finite_difference_coefficient
    """
    accur_coeff_mat = np.array([[0,0,0,-1/2.,0,1/2.,0,0,0],
    [0,0,1/12. ,-2/3.,0,2/3.,-1/12.,0,0],
    [0,-1/60.,3/20.,-3/4.,0,3/4.,-3/20.,1/60.,0],
    [1/280.,-4/105.,1/5.,-4/5.,0,4/5.,-1/5.,4/105.,-1/280.]])
    if i > 3:
        return accur_coeff_mat[-1]
    else:
        return accur_coeff_mat[i]


def partial_derive(f,var_in,var_out=0,kwargs_f={},args_f=[],h=0,accur=-1):
    """
    This function computes the partial derivatives of a python function
    
    Parameters
    ----------
    f : Python function
        the python function which will be derivated.
        this function must return a scalar.
        the parameter of f suseptibles to be derivated must be scalars.
        i.e. if for isntace you want to derivate a position vector X = [x,y,z]
        f must take as argument f(x,y,z) and not f(X)
    var_in : int or string
        the detrivation is with respect to this variable
        can be a int (starts with 0) or a string describing the name of
        the var in f arguments.
    var_out : int
        the output of f which needs to be considerated as the output
        ** must be an int **
        The default is 0.
    kwargs_f : dict, optional
        dictionary describing the arguments of f. The default is {}.
    args_f : iterable, optional
        tuple/list & dict describing the arguments of f. The default is [].
    h : float, optional
        derivation step, if h == 0 give x * sqrt(epsilon)
        (source : http://en.wikipedia.org/wiki/Numerical_differentiation) .
    accur : int, optional
        accuracy coefficient index. -1 is the best but the slowest. The default is -1.
        https://en.wikipedia.org/wiki/Finite_difference_coefficient


    Returns
    -------
    dout : float
        the derivative of f w.r.t. var_in.

    """

    # tuple => list pour plus d'aisance
    args_f = list(args_f)

    # operational arguments
    args_f_i = list(args_f)
    kwargs_f_i = dict(kwargs_f)

    args_name_list = list(f.__code__.co_varnames)
    # var in is a int
    if type(var_in) is int:
        var_ind = var_in
        var_name = args_name_list[var_ind]
    # var in is a str
    else:
        var_name = var_in
        try:
            var_ind = args_name_list.index(var_name)
        except ValueError:
            log.info("arguments %s",args_name_list)
            raise Exception('wrong var_in name (not in args name list)')

    if var_ind < len(args_f):
        x = args_f[var_ind]
    else:
        x = kwargs_f[var_name]
        
    if utils.is_iterable(x):
        log.error("var_in is not a scalar")
        raise Exception

    if h == 0:
        h = x * np.sqrt(np.finfo(float).eps)
        if h == 0:
            log.warning('h == 0 ! setting @ 10**-6')
            h = 10**-6

    res_stk = []
    accur_coeff = get_accur_coeff(accur)
    for i,k in enumerate(accur_coeff):
        ### if the coeff is nul, no computation
        if k == 0:
            res_stk.append(0.)
        
        else:
            if var_ind < len(args_f):
                args_f_i[var_ind]    = x + h * float(i-4)
            else:
                kwargs_f_i[var_name] = x + h * float(i-4)
            ### i-4 because i=0 is actually the index -4
            ### cf Wikipedia table 
            ### https://en.wikipedia.org/wiki/Finite_difference_coefficient
            
            res = f(*args_f_i,**kwargs_f_i)

            if utils.is_iterable(res):
                res = res[var_out]

            res_stk.append(res)

    dout = np.dot(np.array(res_stk) , accur_coeff) / h

    return dout


def jacobian_line(f,var_in_list,var_out=0,kwargs_f={},args_f=[],h=0,aray=True):
    """ same argument as ``partial_derive`` except var_in becomes var_in_list :
    it's an iterable of all variables the derivation must be performed
    """
    line_out =  []
    for var_in in var_in_list:
        line_out.append(partial_derive(f,var_in,var_out,kwargs_f,args_f,h))
    if aray:
        return np.array(line_out)
    else:
        return line_out


def kwargs_for_jacobian(kwdic_generik,kwdic_variables):
    """ Building a list of kwargs for the jacobian function
        kwdic_generik : parameters which not gonna change
        kwdic_variable : parameters which gonna change, so must be associated
        with iterable
        """
    keys_list = list(kwdic_variables.keys())

    for k,v in kwdic_variables.items():
        if not utils.is_iterable(v):
            log.warning('key %s, val %s is not iterable !!!',k,v)

    values_combined = itertools.product(*list(kwdic_variables.values()))

    kwdic_list_out =[]
    for values in values_combined:
        kwdic_out = dict(kwdic_generik)
        for k,v in zip(keys_list,values):
            kwdic_out[k] = v
        kwdic_list_out.append(kwdic_out)

    return kwdic_list_out

def jacobian(f,var_in_list,var_out,kwargs_f_list=[],h=10**-6,nproc=4):
    """ il n'y a que les kwargs qui sont gérés """
    args_f_list = []
    for i in range(len(kwargs_f_list)):
        args_f_list.append([])


    if len(args_f_list) != len(kwargs_f_list):
        log.info("%s, %s",len(args_f_list) , len(kwargs_f_list))
        raise Exception("Jacobian : len(args_f_list) != len(kwargs_f_list)")
    jacob_temp = []
    args_list = []
    for kwargs_f , args_f in zip(kwargs_f_list,args_f_list):
        args_list.append((f,var_in_list,var_out,kwargs_f,args_f,h))

    pool = mp.Pool(processes=nproc)
    results = [pool.apply(jacobian_line, args=x) for x in args_list]
    pool.close()
    jacob_temp = results

    return np.vstack(jacob_temp)


def nan_cleaner(Ain,Bin):
    """
    remove A & B of their respective NaN

    Args :
        Ain , Bin : lists/arrays

    Returns:
        Aout , Bout : A & B withour NaN
    """
    Ain = np.array(Ain)
    Bin = np.array(Bin)

    notnanA = np.logical_not(np.isnan(Ain))
    notnanB = np.logical_not(np.isnan(Bin))

    notnan = notnanA * notnanB

    Aout = Ain[notnan]
    Bout = Bin[notnan]

    return Aout , Bout


def weight_mat(Sinp,Ninp=[],fuvinp=1,sparsediag=False):
    """
    Args :
        Sinp : liste des Sigmas sig = sqrt(var)
        Ninp : liste de la taille de chaque blocs (obs)
        fuvinp = 1 : facteur unitaire de variance
        inspiré de mat_poids , fct écrite dans la lib resolution de GPShApy

    Returns :
        K : matrice de var-covar
        Q : matrice des cofacteurs
        P : matrice des poids inv(Q)
    """

    if Ninp == []:
        Ninp = [1] * len(Sinp)

    Sinp = np.array(Sinp)
    if np.any(Sinp == 0):
        log.warning("some sigma inputs are == 0 !!!")

    if len(Sinp) != len(Ninp):
        raise Exception("S et N de taille differente")

    Ktemp = []
    for i in range(len(Sinp)):
        try:
            Ktemp = Ktemp + Ninp[i] * [Sinp[i]**2]
        except DeprecationWarning:
            log.warning("Are you sure you didn't invert Sinp <> Ninp?")
    Ktemp = np.array(Ktemp).astype(np.float64)
    Qtemp = (1. /fuvinp) * Ktemp
    Qtemp = Ktemp
    Ptemp = 1. / Qtemp
    
    if sparsediag:
        K = scipy.sparse.diags(Ktemp,0)
        Q = scipy.sparse.diags(Qtemp,0)
        P = scipy.sparse.diags(Ptemp,0)
    else:
        K = np.diag(Ktemp)
        Q = np.diag(Qtemp,0)
        P = np.diag(Ptemp,0)

    #    #K = scipy.linalg.block_diag(*Ktemp)
    #    Q = (1/fuvinp) * K
    #    # normalement P = scipy.linalg.inv(Q)
    #    # mais pour que ce soit + rapide
    #    Qdiag = np.diagonal(Q)
    #    Pdiag = 1/Qdiag
    #    P = scipy.linalg.block_diag(*Pdiag)

    return K , Q , P



def weight_mat_simple(Pinp,Ninp=[],sparsediag=False,
                      return_digaonal_only=False):
    """
    Simple version of weight_mat : takes directly the weights (Pinp)
    and the size for each weigths blocks (Ninp)
    
    Pinp and Ninp have to have the same length
    
    Args :
        Pinp : list of weigths
        Ninp : list of the size of each block (obs number)
        fuvinp = 1 : facteur unitaire de variance

    Returns :
        P : weigth matrix
    """

    if Ninp == []:
        Ninp = [1] * len(Pinp)

    if len(Pinp) != len(Ninp):
        raise Exception("different sizes for Pinp and Ninp !!!")

    Ptemp = []
    for i in range(len(Pinp)):
        try:
            Ptemp = Ptemp + Ninp[i] * [Pinp[i]]
        except DeprecationWarning:
            log.warning("Are you sure you didn't invert Pinp <> Ninp ?")

    if return_digaonal_only:
        P = np.array(Ptemp)
    
    else:
        if sparsediag:
            P = scipy.sparse.diags(Ptemp,0)
        else:
            P = np.diag(Ptemp,0)

    return P


def fuv_calc(V,A,P=1,normafuv=1):
    """
    Args :
        V : residuals vector

        A : Jacobian matrix

        P : weight matrix

        Can manage both standard arrays and sparse array

    Returns :
        fuv : Facteur unitaire de variance (unitary variance factor)

    Notes :
        le fuv dépend de la martice de poids
        mais les sigmas non
        ex :
        poids de 10**-6
        fuv    :  439828.260843
        sigma  :  [ 5.21009306  5.09591568  0.04098106]
        poids de 1
        fuv    :  4.39828260843e-07
        sigma  :  [ 5.21009306  5.09591568  0.04098106]
    """
    V = np.squeeze(np.array(V))

    if not utils.is_iterable(P):
        P = np.ones(len(V)) * P
    elif scipy.sparse.issparse(P):
        P = P.diagonal()
    else:
        log.warning("DEPRECIATION : modification done in fuv_calc, P should be given as Matrix-shaped now")
        P = P.diagonal()

    P = P * (1 / normafuv)

    VPV = np.column_stack((V,P,V))
    numera = np.sum(np.product(VPV,1))

    # nummera is just an adaptation of VT * P * V
    # EDIT 1806 : Je veux bien ... mais c'est une adaptation pourrie !
    # ou alors il faut bien s'assurer que l'on a extrait la diagonale de P
    
    
    if A.ndim == 1:
      A = np.expand_dims(A,0)

    fuv = numera / np.abs(A.shape[0] - A.shape[1])
    return fuv


def sigmas_formal_calc(N,V,A,fuv=None,P=None):
    if not fuv:
        fuv      = fuv_calc(V,A,P)
        
    varcovar = (fuv) * scipy.linalg.inv(N)
    sigmas   = np.sqrt(varcovar.diagonal())
    return sigmas , varcovar
  

def smart_i_giver(subgrp_len_list,i_in_sublis,sublis_id,
                  advanced=False,sublis_id_list=[]):
    """
    eg
    subgrp_len_list = [4201, 4186, 4157, 4041, 4058, 4204, 4204, 4204, 4204]
    i_in_sublis = 2
    sublis_id = 3
    return 4201 + 4186 + 4157 + 2

    advanced = True:
    the sublis_id is not a int but and generic identifier ( str ,int , set ... )
    present in sublis_id_list
    else
    must be an int
    """

    if sublis_id_list == [] and not type(sublis_id) is int:
        log.error("smart_i_giver")
        return None

    if advanced:
        sublis_id = sublis_id_list.index(sublis_id)

    if sublis_id == 0:
        decaleur = 0
    else:
        decaleur = -1

    before = np.sum(subgrp_len_list[:sublis_id])
    return int(before + decaleur + i_in_sublis)


def constraint_improve_N(N,C,trans=False,outsparsetype = 'csc'):
    """ give N normal matrix and C constraints matrix
        returns N compined with C
        trans is a (dirty) way to transpose C
        if made in wrong shape

        convention Ghilani 2011 p424 :

        N    C.T
        C     0
    """

    if trans:
        C = C.T

    O = np.zeros((C.shape[0],C.shape[0]))

    if scipy.sparse.issparse(C):
        L1   = scipy.sparse.hstack((N,C.T))
        L2   = scipy.sparse.hstack((C,O))
        Nout = scipy.sparse.vstack((L1,L2))

        if outsparsetype == 'csc':
            Nout = scipy.sparse.csc_matrix(Nout)
        elif outsparsetype == 'csr':
            Nout = scipy.sparse.csr_matrix(Nout)

    else:
        L1 = np.hstack((N,C.T))
        L2 = np.hstack((C,O))
        Nout = np.vstack((L1,L2))

    return Nout

def triangle_arr2vect(triarrin,k=1):
    ind = np.triu_indices_from(triarrin,k)
    return triarrin[ind]


def bins_middle(bin_edges):
    bin_edges2 = []
    for i in range(len(bin_edges)-1):
        bin_edges2.append((bin_edges[i]+bin_edges[i+1])/2.)
    return bin_edges2


def chi2_test_frontend(dist_inp,nbins=10,ddof=2,debug=0,mode2=False,aaa=1):
    """
    mode1 (par def) : on normalise la theorique et pas la observée
    mode2 : on normalise la distribution observée est pas la théorique
    INCOHERENT AVEC MATLAB => A EVITER

    Enfin bon, la manière dont on fabrique les valeurs théoriques est quand même
    un peu vaseuse ... penser à porter le code matlab chi2gof.m l.185
    Et comprendre aussi pourquoi ils ont un ddof de 2 (quon ajoute ici aussi
    par bete copiage) ...

    en debug mode bin_edges,bin_edges2,hist,gauss,chi2

    """
    hist , bin_edges =  np.histogram(dist_inp,nbins,density=mode2)
    bin_edges2 = bins_middle(bin_edges)
    #    gauss = scipy.stats.norm.pdf(bin_edges2,
    #                         np.mean(dist_inp),np.std(dist_inp))
    #    if not mode2:
    #        koefnorm = scipy.integrate.trapz(hist,bin_edges2)
    #    else:
    #        koefnorm = 1.

    cdf = scipy.stats.norm.cdf(bin_edges[1:-1],np.mean(dist_inp),np.std(dist_inp)*aaa)
    cdf = np.concatenate(([0],cdf,[1]))
    dif = np.diff(cdf)

    gauss = len(dist_inp) * dif

    #    gauss1 = gauss
    #    gauss = gauss * koefnorm
    try:
        chi2 = scipy.stats.chisquare(hist , gauss , ddof)
    except:
        chi2 = (np.nan,np.nan)
    if not debug:
        return chi2
    else:
        plt.plot(bin_edges2,gauss,'+')
        plt.plot(bin_edges2,hist,'x')
        return bin_edges,bin_edges2,hist,gauss,chi2


def chi2_test_lsq(V , A ,  P = None , fuvin = None , risk = 0.05,
                  cleaning_std = False , cleaning_normalized = False,
                  koefP = 1):

    """
    P est uniquement la diagonale de la matrice des poids

    les cleaning sont des astuces pour se rapprocher de 1 (en nettoyant les plus mauvaises valeurs)
    cleaning_normalized est à privilégier (et override cleaning_std)

    koefP est coefficient qu'on donne a P pour trouver une solution viable

    si koefP != 1, le nouveau P est donné en avant dernier argument

    """

    if fuvin is None and P is None:
        log.error("fuvin == None and P == None")
        return None

    ddl = (np.max(A.shape) - np.min(A.shape))

    if cleaning_std:
        Vwork = V
        bb = Vwork < np.std(Vwork) * 3

    if cleaning_normalized:
        Vnorma = V / np.sqrt(1/P)
        bb = Vnorma < 3

    if cleaning_std or cleaning_normalized:
        V = V[bb]
        P = P[bb]

    if koefP != 1:
        P = P * koefP

    if fuvin is not None:
        esfuv = fuvin
    else:
        esfuv = np.sum( V.T * P * V ) / ddl

    pmin = scipy.stats.chi2.ppf(risk/2,ddl) / ddl
    pmax = scipy.stats.chi2.ppf(1 - risk/2,ddl) / ddl

    if pmin < esfuv and esfuv < pmax:
        boolchi2 = True
    else:
        boolchi2 = False

    if koefP == 1:
        return esfuv , pmin , pmax , boolchi2
    else:
        return esfuv , pmin , pmax , P , boolchi2
    
    
def error_ellipse(xm,ym,sigx,sigy,sigxy,nsig=1,ne = 100,scale=1):
    """
    from matlab fct
    http://kom.aau.dk/~borre/matlab/geodesy/errell.m
    It works but don't ask why ...

    (X,Y) orientation convention is inverted  => (Y,X) ...
    so in a practical way you must invert X ,Y
    (it is not important for the axis but it is for the orientation)
    AND
    sigx,sigy,sigxy must be first normalized with the fuv

    sigx, sigy, sigxy :
        so as we can generate a covariance matrix
        cov = np.array([[sigx ** 2,sigxy],[sigxy,sigy ** 2]])

    ne :
        nb of segements for the ellipse

    RETURNS :
        xe,ye,dx2,dy2


    DEBUG:

    si on a
    xe1,ye1,_,_ = stats.error_ellipse(pxp[0],pxp[1], sigxB , sigyB , sigxyB, scale= 10000)
    xe2,ye2,_,_ = stats.error_ellipse(pxp[0],pxp[1], sigyB , sigxB , sigxyB, scale= 10000)
    et PAS les - à D et dxy0 => on a 2 ellipses differentes

    si on a
    xe1,ye1,_,_ = stats.error_ellipse(pxp[0],pxp[1], sigxB , sigyB , sigxyB, scale= 10000)
    xe2,ye2,_,_ = stats.error_ellipse(pxp[0],pxp[1], sigyB , sigxB , sigxyB, scale= 10000)
    et AVEC les - à D et dxy0 => on a 2 ellipses differentes
    au moins une ellipse coincide avec celle de Ghiliani

    A investiguer, en attendant, à éviter
    """

    cov = np.array([[sigx ** 2,sigxy],[sigxy,sigy ** 2]])
    V,D = np.linalg.eig(cov)

    D = -D.T

    std1 = np.sqrt(V[0])
    std2 = np.sqrt(V[1])

    if std1 < std2:
        z1 = np.linspace(-std1,std1,ne)
        z2 = np.sqrt(V[1] * (np.ones(ne) - (z1 / std1) **2))
        dxy = np.dot(D,np.vstack((np.hstack((z1,np.flipud(z1))),np.hstack((z2,-z2)))))
    else:
        z2 = np.linspace(-std2,std2,ne)
        z1 = np.sqrt(V[0] * (np.ones(ne) - (z2 / std2) **2))
        dxy = np.dot(D,np.vstack((np.hstack((z1,-z1)),np.hstack((z2,np.flipud(z2))))))

    dx =  -dxy[0,:]
    dy =  dxy[1,:]

    dx2 =  scale*nsig*dx
    dy2 =  scale*nsig*dy
    xe = xm * np.ones(2*ne) + dx2
    ye = ym * np.ones(2*ne) + dy2

    return xe,ye,dx2,dy2

def error_ellipse_parameters(qxx,qyy,qxy,fuv,out_t=False):
    """
    INPUT :
        qxx,qyy,qxy  : factors as in the varcovar matrix (no normalisation with the fuv or other)
        fuv
    OUTPUT :
        Su/a Sb/b    : semimajor and semiminor axis
        t            : angle that the u/a axis makes with the y axis in clockwise direction
        OR
        phi          : angle that the u/a axis makes with the x axis in trigo direction

    (this one is the natural way, perfect to test those of the Strang & Borre)
    """
    S0   = np.sqrt(fuv)
    t    = 0.5 * np.arctan2( 2*qxy , (qyy - qxx))
    tdeg = np.rad2deg(t)

    quu = qxx * np.sin(t)**2 + 2*qxy * np.cos(t) * np.sin(t) + qyy * np.cos(t)**2
    qvv = qxx * np.cos(t)**2 - 2*qxy * np.cos(t) * np.sin(t) + qyy * np.sin(t)**2
    Su = S0 * np.sqrt(quu)
    Sv = S0 * np.sqrt(qvv)

    a   = Su
    b   = Sv
    phi = 90 - tdeg
    if out_t:
        return a,b,tdeg
    else:
        return a,b,phi


def error_ellipse_parameters_2(sigx,sigy,sigxy,out_deg=True):
    """
    ref : Strang & Borre p 337
    (X,Y) orientation convention is inverted  => (Y,X) ...
    so in a practical way you must invert X ,Y
    (it is not important for the axis but it is for the orientation)
    AND
    sigx,sigy,sigxy must be normalized with the fuv
    """

    cov = np.array([[sigx ** 2,sigxy],[sigxy,sigy ** 2]])

    sigx2 = sigx ** 2
    sigy2 = sigy ** 2
    sigxy2 = sigxy ** 2

    cov = np.array([[sigx ** 2,sigxy],[sigxy,sigy ** 2]])

    V,D   = np.linalg.eig(cov)
    V2,D2 = np.linalg.eig(np.linalg.inv(cov))

    #Dsqrt = np.sqrt(V)

    l1 = 0.5 * (sigx2 + sigy2 + np.sqrt((sigx2 + sigy2)**2 - 4*(sigx2 * sigy2 - sigxy2)))
    l2 = 0.5 * (sigx2 + sigy2 - np.sqrt((sigx2 + sigy2)**2 - 4*(sigx2 * sigy2 - sigxy2)))

    a = np.sqrt(l1)
    b = np.sqrt(l2)
    phi = ( 0.5 * np.arctan2(2*sigxy,sigx2-sigy2))

    if out_deg:
        phi = np.rad2deg(phi)
    return a,b,phi


def ellipse_get_coords(a=0.0, b=0.0, x=0.0, y=0.0, angle=0.0, k=2 ,
                       out_separate_X_Y = True , trigo = True):
    """ Draws an ellipse using (360*k + 1) discrete points; based on pseudo code
    given at http://en.wikipedia.org/wiki/Ellipse
    k = 1 means 361 points (degree by degree)
    a = major axis distance,
    b = minor axis distance,
    x = offset along the x-axis
    y = offset along the y-axis
    angle = trigo/clockwise rotation [in degrees] of the ellipse;
        * angle=0  : the ellipse is aligned with the positive x-axis
        * angle=30 : rotated 30 degrees trigo/clockwise from positive x-axis

    trigo sense is the standard convention

    NB : clockwise is the internal convention, but we prefer trigo
         convention for the Ghiliani ellipses made by error_ellipse_parameters

    source : scipy-central.org/item/23/1/plot-an-ellipse
    """

    if trigo:
        angle = - angle

    pts = np.zeros((360*k+1, 2))

    beta = -angle * np.pi/180.0
    sin_beta = np.sin(beta)
    cos_beta = np.cos(beta)
    alpha = np.radians(np.r_[0.:360.:1j*(360*k+1)])

    sin_alpha = np.sin(alpha)
    cos_alpha = np.cos(alpha)

    pts[:, 0] = x + (a * cos_alpha * cos_beta - b * sin_alpha * sin_beta)
    pts[:, 1] = y + (a * cos_alpha * sin_beta + b * sin_alpha * cos_beta)

    if not out_separate_X_Y:
        return pts
    else:
        return pts[:, 0] , pts[:, 1]

def ellipse_center(a):
    """
    http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
    """
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])


def ellipse_angle_of_rotation( a , outdeg = True ):
    """
    core fct for ellipse_fit
    http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
    """
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    if not outdeg:
        return 0.5*np.arctan(2*b/(a-c))
    else:
        return np.rad2deg(0.5*np.arctan(2*b/(a-c)))


def ellipse_axis_length( a ):
    """
    core fct for ellipse_fit
    http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
    """
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(up/down1)
    res2=np.sqrt(up/down2)
    return np.array([res1, res2])

def fitEllipse_core(x,y):
    """
    core fct for ellipse_fit
    http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
    """
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  np.linalg.eig(np.dot(np.linalg.inv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    return a

def ellipse_fit(x,y):
    """
    find the parameters
    a,b,phi,x0,y0 of an ellipse

    from
    http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
    """
    tmp   = fitEllipse_core(x,y)
    a,b   = ellipse_axis_length( tmp )
    phi   = ellipse_angle_of_rotation(tmp)
    x0,y0 = ellipse_center(tmp)
    return a,b,phi,x0,y0


 #  ______                _   _                _____ _                _                  
 # |  ____|              | | (_)              / ____(_)              | |                 
 # | |__ _   _ _ __   ___| |_ _  ___  _ __   | |     _ _ __ ___   ___| |_ ___ _ __ _   _ 
 # |  __| | | | '_ \ / __| __| |/ _ \| '_ \  | |    | | '_ ` _ \ / _ \ __/ _ \ '__| | | |
 # | |  | |_| | | | | (__| |_| | (_) | | | | | |____| | | | | | |  __/ ||  __/ |  | |_| |
 # |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|  \_____|_|_| |_| |_|\___|\__\___|_|   \__, |
 #                                                                                  __/ |
 #                                                                                 |___/ 



def clean_nan(A,L):
    """
    DISCONTINUED

    A est un array bi dimentionnel
    L est un array mono dimentionnel
    renvoie un A et un L nettoyé mutuellement
    de leurs NaN respectifs
    return np.sqrt(a + a.T - np.diag(a.diagonal()))
    """
    if A.ndim != 2  or L.ndim != 1:
        raise Exception("ERREUR : revoir les dimensions de A (2) et L (1)")

    if A.shape[0] != L.shape[0]:
        raise Exception("ERREUR : A et L ont des longeurs differentes")

    # 11) Detecter les nan de L
    indLnonan = np.where(~np.isnan(L))
    nbnanL = np.isnan(L).sum()
    nbnanAtot = np.isnan(A).sum()
    nbnanA = np.isnan(A).any(1).sum()

    # 12) Virer les nan de L dans A
    A2 = A[indLnonan[0],:]
    # 13) Virer les nan de L dans L
    L2 = L[indLnonan]
    # 21) Detecter les nan de A
    # (même si cetains ont potentiellement déjà été supprimé)
    boolA2arenan = np.isnan(A2)
    nbnanA2 = boolA2arenan.any(1).sum()
    # 22) Virer les nan de A dans A
    Aout = A2[~boolA2arenan.any(1)]

    # 23) Virer les nan de A dans L
    Lout = L2[~boolA2arenan.any(1)]

    if  np.isnan(Lout).sum() != 0 or  np.isnan(Aout).sum() != 0 :
        log.info("pour info, isnan(Lout).sum(), isnan(Aout).sum()")
        log.info(np.isnan(Lout).sum(), np.isnan(Aout).sum())
        raise Exception("ERREUR : A_out et L_out contiennent des NaN")

    if Aout.shape[0] != Lout.shape[0]:
        raise Exception("ERREUR : A_out et L_out ont des longeurs differentes")

    log.info("%i NaN dans le V",nbnanL)
    log.info("%i NaN dans la M vert., APRES suppr. des NaN du V", nbnanA2)
    log.info("%i lignes supprimées dans le V et la M (en théorie)", nbnanA2 + nbnanL)
    log.info("")
    log.info("pour info :")
    log.info("%i NaN dans la M vert., AVANT suppr. des NaN du V", nbnanA)
    log.info("%i NaN dans la M AU TOTAL (vert. et horiz.)", nbnanAtot)

    return Aout , Lout


def fuv_calc_OLD(V,A):
    fuv = np.dot(V.T,V) / np.abs(A.shape[0] - A.shape[1])
    return fuv

def fuv_calc_OLD2(V,A,P=None):
    if P is None:
        P = np.eye(len(V))

    P = np.matrix(P)
    V = np.matrix(V)
    A = np.matrix(A)

    if V.shape[0] == 1:
        V = V.T

    fuv = (V.T * P * V) / np.abs(A.shape[0] - A.shape[1])
    return fuv[0,0]

#def weight_mat_OLD(Sinp,Ninp,fuvinp=1):
#    """ Sinp : liste des Sigmas sig = sqrt(var)
#        Ninp : liste de la taille de chaque blocs (obs)
#        fuvinp = 1 : facteur unitaire de variance
#        inspiré de mat_poids , fct écrite dans la lib resolution de GPShApy
#
#        retourne :
#        K : matrice de var-covar
#        Q : matrice des cofacteurs
#        P : matrice des poids inv(Q)
#    """
#
#    if len(Sinp) != len(Ninp):
#        raise Exception("S et N de taille differente")
#
#    Ktemp = []
#
#    for i in range(len(Sinp)):
#        try:
#            Ktemp.append(np.eye(Ninp[i]) * Sinp[i]**2)
#        except DeprecationWarning:
#            print "weight_mat : Are you sure you don't invert Sinp <> Ninp ?"
#    K = scipy.linalg.block_diag(*Ktemp)
#    Q = (1/fuvinp) * K
#    # normalement P = scipy.linalg.inv(Q)
#    # mais pour que ce soit + rapide
#    Qdiag = np.diagonal(Q)
#    Pdiag = 1/Qdiag
#    P = scipy.linalg.block_diag(*Pdiag)
#
#    return K , Q , P


def partial_derive_old(f,var_in,var_out=0,kwargs_f={},args_f=[],h=0):
    ''' var_in :
            detrivation with respect to this variable
            can be a int (starts with 0) or a string descirbing the name of
            the var in f
        var_out :
            the output of f which needs to be considerated as the output
            ** must be a int **
        args_f & kwargs_f :
            tuple/list & dict describing the arguments of f
        h :
            derivation step, if h == 0 give x * sqrt(epsilon)
            (source : http://en.wikipedia.org/wiki/Numerical_differentiation) '''

    # tuple => list pour plus d'aisance
    args_f = list(args_f)

    # operational arguments
    args_f_m = list(args_f)
    args_f_p = list(args_f)
    kwargs_f_m = dict(kwargs_f)
    kwargs_f_p = dict(kwargs_f)

    args_name_list = list(f.__code__.co_varnames)
    # var in is a int
    if type(var_in) is int:
        var_ind = var_in
        var_name = args_name_list[var_ind]
    # var in is a str
    else:
        var_name = var_in
        try:
            var_ind = args_name_list.index(var_name)
        except ValueError:
            log.error(args_name_list)
            raise Exception('wrong var_in name (not in args name list)')

#    if var_ind < len(args_f):
#        x = args_f[var_ind]
#        if h == 0:
#            h = x * np.sqrt(np.finfo(float).eps)
#            if h == 0:
#                print 'WARN : h == 0 ! setting @ 10**-6 '
#                h = 10**-6
#        args_f_m[var_ind] = x - h
#        args_f_p[var_ind] = x + h
#    else:
#        x = kwargs_f[var_name]
#        if h == 0:
#            h = x * np.sqrt(np.finfo(float).eps)
#            if h == 0:
#                print 'WARN : h == 0 ! setting @ 10**-6 '
#                h = 10**-6
#        kwargs_f_m[var_name] = x - h
#        kwargs_f_p[var_name] = x + h


    if var_ind < len(args_f):
        x = args_f[var_ind]
    else:
        x = kwargs_f[var_name]

    if h == 0:
        h = x * np.sqrt(np.finfo(float).eps)
        if h == 0:
            log.warning('WARN : h == 0 ! setting @ 10**-6 ')
            h = 10**-6

    if var_ind < len(args_f):
        args_f_m[var_ind] = x - h
        args_f_p[var_ind] = x + h
    else:
        kwargs_f_m[var_name] = x - h
        kwargs_f_p[var_name] = x + h

    m = f(*args_f_m,**kwargs_f_m)
    p = f(*args_f_p,**kwargs_f_p)

    if utils.is_iterable(m):
        m = m[var_out]
        p = p[var_out]

    dout = (p - m) / (2. * float(h))

    return dout



#def jacobian_old_bkp(f,var_in_list,var_out,kwargs_f_list=[],args_f_list=[],h=10**-6):
#    if args_f_list == [] and kwargs_f_list == []:
#        raise Exception('jacobian : args_f_list == [] and kwargs_f_list ==  []')
#    elif args_f_list == [] and kwargs_f_list != []:
#        args_f_list =  [ [] * len(kwargs_f_list)]
#    elif args_f_list != [] and kwargs_f_list == []:
#        kwargs_f_list =  [ [] * len(args_f_list)]
#
#    if len(args_f_list) != len(kwargs_f_list):
#        raise Exception("Jacobian : len(args_f_list) != len(kwargs_f_list)")
#    jacob_temp = []
#    for kwargs_f , args_f in zip(kwargs_f_list,args_f_list):
#        line = jacobian_line(f,var_in_list,var_out,kwargs_f,args_f,h)
#        jacob_temp.append(line)
#
#    return np.vstack(*jacob_temp)
