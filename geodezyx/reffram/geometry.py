#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 17:15:57 2019

@author: psakicki
"""

from geodezyx import *
from geodezyx import np,scipy


def BL_from_points(listpointin):
    ''' 
    A partir d'une liste de points,
    retourne les baselines entre ces points dans une matrice 
    '''
    """    
    From a list of 2-D or 3-dD points, returns the a matrix with distance 
    between each points 
    
    Parameters
    ----------
    listpointin : list or numpy.array
        List of N 2D or 3D points [[x1,y1,z1] ... [xn , yn , zn]]
                
    Returns
    -------
    BL : numpy.array
        matrix with distances between each points
        
    """

    N = len(listpointin)
    BL = np.empty((N,N))

    for i,pt1 in enumerate(listpointin):
        for j,pt2 in enumerate(listpointin):

            if i == j:
                BL[i,j] = 0
            else:
                BL[i,j] = np.linalg.norm(pt1 - pt2)

    return BL



def mat_poids(Sinp,Ninp,fuvinp=1):
    """
    discontinued
    """
    # Sinp : liste des Sigmas sig = sqrt(var)
    # Ninp : liste de la taille de chaque blocs (obs)
    # fuvinp = 1 : facteur unitaire de variance

    if len(Sinp) != len(Ninp):
        raise Exception("S et N de taille differente")

    Ktemp = []

    for i in range(len(Sinp)):
        print(Sinp[i])
        Ktemp.append(np.eye(Ninp[i]) * Sinp[i]**2)

    K = scipy.linalg.block_diag(*Ktemp)
    Q = (1/fuvinp) * K
    P = scipy.linalg.inv(Q)

    return K , Q , P

def rotmat2(theta,angtype='deg'):

    if angtype == 'deg':
        theta = np.deg2rad(theta)

    rotmat = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])

    return rotmat


def rotmat3(alpha,beta,gamma,xyzreftuple = ([1, 0, 0], [0, 1, 0], [0, 0, 1]),angtype='deg'):

    xaxis, yaxis, zaxis = xyzreftuple

    if angtype == 'deg':
        alpha = np.deg2rad(alpha)
        beta  = np.deg2rad(beta)
        gamma = np.deg2rad(gamma)

    Rx = trans.rotation_matrix(alpha, xaxis)
    Ry = trans.rotation_matrix(beta, yaxis)
    Rz = trans.rotation_matrix(gamma, zaxis)
    R = trans.concatenate_matrices(Rz, Ry, Rx)[:3,:3]

    return R

def rotate_points(alphal,betal,gammal,pointlin,Rtype='R1',
                  xyzreftuple = ([1, 0, 0], [0, 1, 0], [0, 0, 1]),
                  angtype='deg',fullout = False):
    '''
    R1  = Rz(g) * Ry(b) * Rx(a)
         si les RPY sont donnés dans le NED
         alors les positions résultantes sont dans le NED

    R2  =  matrice RPY2ENU
        si les RPY sont donnés dans le NED
        alors les  résultantes sont DANS LE ENU
        pas besoin de rotation NED2ENU

        Grewal et al. 2007

    Entrée :
        Angles n = A
        liste de listes de P * [ points ]

    Sortie :
        liste de listes [ [ xA ] [ xA ] ... xP [ xA ] ]  '''

    xaxis, yaxis, zaxis = xyzreftuple

    if not utils.is_iterable(alphal):
        alphal = np.array([alphal])
        betal = np.array([betal])
        gammal = np.array([gammal])
        boolnotiterable = True
    else:
        boolnotiterable = False

    pointlout = []
    R_out = []


    for pt in pointlin:

        if not utils.is_iterable(pt) or len(pt) != 3:
            print("ERR : rotate_points : pts != 3 coords")
            return 0

        pointltmp = []

        for a,b,g in zip(alphal,betal,gammal):

            R1 = rotmat3(a,b,g,angtype=angtype,xyzreftuple=xyzreftuple)
            R2 = C_rpy2enu(a,b,g,angtype=angtype)

            if Rtype == 'R1':
                R = R1
            elif Rtype == 'R2':
                R = R2
            R_out.append(R)

            pointltmp.append(np.dot(R,pt))

        pointlout.append(pointltmp)

        if boolnotiterable:
            pointlout = pointltmp

        pointlout = np.array(pointlout)

    if fullout:
        return pointlout , R_out
    else:
        return pointlout


def guess_seq_len(seq):
    #source
    #http://stackoverflow.com/questions/11385718/python-finding-repeating-sequence-in-list-of-integers
    guess = 1

    if len(set(seq)) == 1:
        return 1

    max_len = len(seq) / 2
    for x in range(2, max_len):
        if seq[0:x] == seq[x:2*x] :
            return x

    return guess

def wrapTo2Pi(lon):
    """
     wrapTo2Pi Wrap angle in radians to [0 2*pi]

    lambdaWrapped = wrapTo2Pi(LAMBDA) wraps angles in LAMBDA, in radians,
    to the interval [0 2*pi] such that zero maps to zero and 2*pi maps
    to 2*pi. (In general, positive multiples of 2*pi map to 2*pi and
    negative multiples of 2*pi map to zero.)

    See also wrapToPi, wrapTo180, wrapTo360.

    """
    lon = np.array(lon)
    positiv = lon > 0
    outlon = np.mod(lon , 2*np.pi)
    outlon[np.logical_and(outlon == 0 , positiv)] = 2 * np.pi
    return outlon


def wrapToPi(lon):
    """
    wrapToPi Wrap angle in radians to [-pi pi]

       lambdaWrapped = wrapToPi(LAMBDA) wraps angles in LAMBDA, in radians,
       to the interval [-pi pi] such that pi maps to pi and -pi maps to
       -pi.  (In general, odd, positive multiples of pi map to pi and odd,
       negative multiples of pi map to -pi.)

       See also wrapTo2Pi, wrapTo180, wrapTo360.

    """

    outlon = np.array(lon)
    q =  np.logical_and((outlon < -np.pi) , (np.pi < outlon))
    outlon[q] = wrapTo2Pi(outlon[q] + np.pi) - np.pi
    return outlon


def wrapTo180(lonin):
    """
    wrapTo180 Wrap angle in degrees to [-180 180]

    lonWrapped = wrapTo180(LON) wraps angles in LON, in degrees, to the
    interval [-180 180] such that 180 maps to 180 and -180 maps to -180.
    (In general, odd, positive multiples of 180 map to 180 and odd,
    negative multiples of 180 map to -180.)

    See also wrapTo360, wrapTo2Pi, wrapToPi.
    """
    lon = np.array(lonin)
    q = (lon < -180) and (180 < lon)
    lon[q] = wrapTo360(lon[q] + 180) - 180

    return lon

def wrapTo360(lonin):
    """
    wrapTo360 Wrap angle in degrees to [0 360]

    lonWrapped = wrapTo360(LON) wraps angles in LON, in degrees, to the
    interval [0 360] such that zero maps to zero and 360 maps to 360.
    (In general, positive multiples of 360 map to 360 and negative
    multiples of 360 map to zero.)

    See also wrapTo180, wrapToPi, wrapTo2Pi.
    """

    lon = np.array(lonin)

    positiveInput = (lon > 0)
    lon = np.mod(lon, 360)
    lon[(lon == 0) & positiveInput] = 360
    return lon



# Pas convaincu de son utilité
def unwrap180(anglist,angtype='deg'):

    if angtype == 'deg':
        seuil = 360

    angout = []

    for a in anglist:
        if a > seuil / 2:
            a = a - seuil
        angout.append(a)

    return angout


def wrap360(anglist,angtype='deg'):

    angout = []

    if angtype == 'deg':
        seuil = 360
    elif angtype == 'rad':
        seuil = 2*np.pi

    for a in anglist:
        if a < 0:
            a = a + seuil

        angout.append(a)

    return angout

class interp1d_ang():

    def __init__(self,T,A,angtype='deg',kind='linear',bounds_error=False):

        if angtype == 'deg':
            A = np.deg2rad(A)

        self.A = A
        self.T = T
        self.C = np.cos(A)
        self.S = np.sin(A)

        self.CfT = interpolate.interp1d(T,self.C,kind=kind,bounds_error=bounds_error)
        self.SfT = interpolate.interp1d(T,self.S,kind=kind,bounds_error=bounds_error)


    def __call__(self,T,angtype='deg'):

        I = np.arctan2(self.SfT(T) ,self.CfT(T) )
        I = wrap360(I,angtype='rad')

        if angtype == 'deg':
            return np.rad2deg(I)
        else:
            return I

def group_consecutives(vals, step=1):
    """
    Return list of consecutive lists of numbers from vals (number list).
    """
    run = []
    result = [run]
    expect = None
    for v in vals:
        if (v == expect) or (expect is None):
            run.append(v)
        else:
            run = [v]
            result.append(run)
        expect = v + step

        result2 = []
        for r in result:
            if len(r) > 1:
                result2.append([r[0],r[-1]])
            else:
                result2.append(r)

    return result2

def randomwalk_normal(N=100, d=2 , moy = 0 , sigma = 1):
    """
    d = dimension
    """
    return np.cumsum(moy + np.random.randn(N,d) * sigma)

def randomwalk_uniform(N=100, d=2 , bound = 0.5):
    """
    d = dimension
    bound = contraint of the random walk
    """
    return np.cumsum(np.random.uniform(-bound,bound,(N,d)))


def circle_draw(xc,yc,R,N):
    theta = np.linspace(0,2 * np.pi,N)
    X = np.cos(theta) * R + xc
    Y = np.sin(theta) * R + yc
    return X,Y


def random_walk_in_a_circle(x0 , y0 , xc , yc ,
                            R , N , step_size ,  param = 1 ,
                            polar = True , uniform_or_normal = 'n',
                            rand_seed = -1):

    """
    random : normal ou uniform
    coords : polaire ou cartésien

    param est un paramètre très versatile pour controler le random :
    plage pour le uniform
    sigma pour le normal

    on recommande plutot le polaire normal : on a un pas constant & une derive réaliste sur le cap

    Returns :
        X,Y, Xcircle , Ycircle

    Exemple :

    for un in ('u','n'):
        for pol in range(2):
            X,Y , Xcircle , Ycircle = random_walk_in_a_circle(10,10,0,0,50,10000,polar = pol,uniform_or_normal=un)

            plt.figure()
            plt.plot(Xcircle,Ycircle)
            plt.plot(X,Y)
            plt.axis('equal')
            plt.suptitle(un + str(pol))
    """

    X = [x0]
    Y = [y0]

    if rand_seed > -1:
        RAND = np.random.RandomState(rand_seed)
    else:
        RAND = np.random.RandomState(np.random.randint(10**6))

    Xcircle,Ycircle = circle_draw(xc,yc,R,500)

    for i in range(N-1):
        D = R+1
        iwhil = 0
        while D > R:
            iwhil += 1
            if iwhil > 500:
                print('WARN : infinite loop in random_walk_in_a_circle ...' , iwhil)
            if polar:
                if uniform_or_normal == 'u':
                    dalpha = RAND.uniform(-param,param) * 2 * np.pi
                else:
                    dalpha = RAND.normal(0,param)       * 2 * np.pi
                drho = step_size
                dx = drho * np.cos(dalpha)
                dy = drho * np.sin(dalpha)
                #print dx,dy
            else:
                if uniform_or_normal == 'u':
                    dx = np.random.uniform(-param,param)
                    dy = np.random.uniform(-param,param)
                else:
                    dx = np.random.normal(0,param)
                    dy = np.random.normal(0,param)
                print(dx , dy)
            xtemp = X[-1] + dx
            ytemp = Y[-1] + dy
            D = np.sqrt((xtemp - xc)**2 + (ytemp - yc)**2)

        X.append(xtemp)
        Y.append(ytemp)

    X = np.array(X)
    Y = np.array(Y)

    return X,Y, Xcircle , Ycircle



def randn_bool(N,true_ratio = 0.5,RandGene = None):
    if RandGene is None:
        RandGene = np.random.RandomState()
    if type(RandGene) is int:
        RandGene = np.random.RandomState(RandGene)
    try:
        randlis = RandGene.uniform(size=N)
    except AttributeError:
        "ERR : AttributeError : RandGene  may be an int32/int64, but it's not an authentic int as required ..."
    boolout_lis = []
    for r in randlis:
        if r < true_ratio:
            boolout_lis.append(True)
        else:
            boolout_lis.append(False)

    return boolout_lis


def points_circle_border(Npts,r,r_sigma,az_type_normal=True,
                         main_dir=3.14159,dir_range=3.14159,seed=None):
    if not seed:
        seed = np.random.randint(10000)

    S = np.random.RandomState(seed)

    if not az_type_normal:
        Az = S.rand(Npts)  * 2 * np.pi
    else:
        Az = S.randn(Npts) * dir_range + main_dir


    R = np.array(Npts * [r]) - np.abs(S.randn(Npts) * r_sigma)

    X,Y = polar2cartesian(R,Az,'rad')

    return X , Y


def estimated_autocorrelation(x):
    """
    http://stackoverflow.com/q/14297012/190597
    http://en.wikipedia.org/wiki/Autocorrelation#Estimation
    """
    n = len(x)
    variance = x.var()
    x = x-x.mean()
    r = np.correlate(x, x, mode = 'full')[-n:]
    assert np.allclose(r, np.array([(x[:n-k]*x[-(n-k):]).sum() for k in range(n)]))
    result = r/(variance*(np.arange(n, 0, -1)))
    return result

def savage_buford_formula(Vs,X,d):
    """
    X : distance à la faille , un iterable pour toutes le profil,
    un nombre pour la longeur max

    d : profondeur de la faille
    retourne X , et Vdeform(X)

    X et d doivent être dans la même unité, Vs pas forcément
    """

    if not utils.is_iterable(X):
        X = np.arange(-X,X,1)
    return X , ( Vs / np.pi ) * np.arctan2(X,d)




def R2_calc(y_obs,y_fit,with_r2_bis=False):
    #https://en.wikipedia.org/wiki/Coefficient_of_determination
    ybar = np.mean(y_obs)
    SStot = np.sum((y_obs - ybar)**2)
    SSreg = np.sum((y_fit - ybar)**2)
    SSres = np.sum((y_obs - y_fit)**2)

    r2 = 1. - ( SSres / SStot)
    r2bis = ( SSreg / SStot)

    if not with_r2_bis:
        return r2
    else:
        return r2 , r2bis

def R2_from_a_line_regress(Xobs,Yobs,a,b):
    #https://en.wikipedia.org/wiki/Coefficient_of_determination
    Xfit , Yfit = stats.linear_reg_getvalue(Xobs,a,b)
    r2 = R2_calc(Yobs,Yfit)
    return r2



#  _    _ _       _       _                    _    _____                _      _   _        ______   _       
# | |  | (_)     | |     | |                  | |  / ____|              | |    | | (_)      |  ____| | |      
# | |__| |_  __ _| |__   | |     _____   _____| | | |  __  ___  ___   __| | ___| |_ _  ___  | |__ ___| |_ ___ 
# |  __  | |/ _` | '_ \  | |    / _ \ \ / / _ \ | | | |_ |/ _ \/ _ \ / _` |/ _ \ __| |/ __| |  __/ __| __/ __|
# | |  | | | (_| | | | | | |___|  __/\ V /  __/ | | |__| |  __/ (_) | (_| |  __/ |_| | (__  | | | (__| |_\__ \
# |_|  |_|_|\__, |_| |_| |______\___| \_/ \___|_|  \_____|\___|\___/ \__,_|\___|\__|_|\___| |_|  \___|\__|___/
#            __/ |                                                                                            
#           |___/                                                                         
    
### High level geodetic function

def itrf_speed_calc(x0,y0,z0,t0,vx,vy,vz,t):
    """
    Args :
        x0,y0,z0 (floats) : coordinates at the reference epoch (m)

        t0 (float) : reference epoch (decimal year)

        vx,vy,vz (floats) : speed of the point (m/yr)

        t (float) : output epoch
    Returns :
        xout,yout,zout : coordinates of the point @ the ref. epoch (m)
    """

    xout = x0 + vx * ( t - t0 )
    yout = y0 + vy * ( t - t0 )
    zout = z0 + vz * ( t - t0 )

    return xout,yout,zout


def itrf_psd_fundamuntal_formula(t,A_l,t_l,tau_l,A_e,t_e,tau_e):
    """
    http://itrf.ensg.ign.fr/ITRF_solutions/2014/doc/ITRF2014-PSD-model-eqs-IGN.pdf
    """
    
    dL = A_l * np.log(1 + (t - t_l)/ tau_l) + A_e * (1 + (t - t_e)/ tau_e)
    
    return dL


def calc_pos_speed_itrf(x0,y0,z0,t0,vx,vy,vz,t):
    """
    just a wrapper of itrf_speed_calc
    for legacy reasons
    """
    return itrf_speed_calc(x0,y0,z0,t0,vx,vy,vz,t)


def helmert_trans(Xa,params='itrf2008_2_etrf2000',invert=True,workepoc=2009.):
    """
    NB 1 : http://etrs89.ensg.ign.fr/memo-V8.pdf
    NB 2 :
    Transformation inverse : * -1 pour les paramètres
    cf https://en.wikipedia.org/wiki/Helmert_transformation

    optimisé pour RGF93 => ITRF2008 (d'ou le invert = True & workepoc = 2009.)

    NB3 : Attention losque l'on compare avec
          le convertisseur EUREF avec une vitesse
          parce que elle aussi est modifiée dans la conversion ETRS => ITRS.
          Conclusion : le faire en 2 étapes ETRS = > ITRS dans la même epoc
                                            ITRS epoc 1 => ITRS epoc 2

    if manual
    then params is a tuple
    params = (t1,t2,t3,dab,r1,r2,r3)

    NE MARCHE PAS PARCE BESOIN DU RATE(TAUX) DES PARAMS D'HELMERT !!!!!!
    (160923)
    """
    if invert:
        inver = -1.
    else:
        inver = 1.

    mas2rad = 0.0000000048481368
    mas2rad = 4.8481368111e-6 * 1e-3

    if params == 'itrf2008_2_etrf2000':

        t1rate = .1 * 10**-3
        t2rate = .1 * 10**-3
        t3rate = -1.8 * 10**-3
        dabrate = .08 * 10**-9
        r1rate =  .081 * mas2rad
        r2rate =  .490 * mas2rad
        r3rate = -.792 * mas2rad

        t1rate = .1 * 10**-3
        t2rate = .1 * 10**-3
        t3rate = -1.8 * 10**-3
        dabrate = .08 * 10**-9
        r1rate =  .081 * mas2rad
        r2rate =  .490 * mas2rad
        r3rate = -.792 * mas2rad


        t1  =(52.1  * 10**-3   + t1rate  * ( workepoc - 2000.)) *inver
        t2  =(49.3  * 10**-3   + t2rate  * ( workepoc - 2000.)) *inver
        t3  =(-58.5 * 10**-3   + t3rate  * ( workepoc - 2000.)) *inver
        dab =( 1.34 * 10**-9   + dabrate * ( workepoc - 2000.)) *inver
        r1  =( 0.891 * mas2rad + r1rate  * ( workepoc - 2000.)) *inver
        r2  =( 5.39  * mas2rad + r2rate  * ( workepoc - 2000.)) *inver
        r3  =( -8.712* mas2rad + r3rate  * ( workepoc - 2000.)) *inver


    elif params == 'itrf2000_2_etrf2000':
        t1  =54.0  * 10**-3   *inver
        t2  =51.0  * 10**-3   *inver
        t3  =-48.0 * 10**-3   *inver
        dab = 0.0  * 10**-9   *inver
        r1  = 0.891 * mas2rad *inver
        r2  = 5.390 * mas2rad *inver
        r3  = -8.712* mas2rad *inver




    R = np.matrix([[dab,-r3,r2],
                   [r3,dab,-r1],
                   [-r2,r1,dab]])

    Xb = Xa + np.matrix([t1,t2,t3]) + np.dot(R,Xa)
    Xb = np.squeeze(np.array(Xb))

    return Xb


def helmert_trans_estim_matrixs_maker(X1 , X2):
    """
    internal function for helmert_trans_estim
    """
    x1 , y1 , z1 = X1
    x2 , y2 , z2 = X2
    
    
    block_1 = np.eye(3)
    
    block_2 = np.array([[ 0. , -z1,  y1, x1],
                       [ z1,   0., -x1, y1],
                       [-y1,  x1,   0., z1]])
    
    l = X2 - X1
    A = np.hstack((block_1 , block_2))

    return l,A


def helmert_trans_estim(X1list , X2list, Weights=[]):
    """
    estimates 7 parameters of a 3D Helmert transformation between a set of points
    X1 and a set of points X2 (compute transformation X1 => X2)
    
    Parameters
    ----------
    
    X1list & X2list : list of N (x,y,z) points ,
        or an numpy array of shape (3, N)

    Weights : list of N Wieghtd,
        or an numpy array of shape (3, N)

    
    Returns
    -------
    Res :
        7 Helmert params. : x,y,z translations, x,y,z rotations, scale
    A :
        Design matrix    
    l :
        X2 - X1
        
    Source
    ------
    https://elib.uni-stuttgart.de/bitstream/11682/9661/1/BscThesis_GaoYueqing.pdf
    """

    l_stk = []
    A_stk = []
        
    for X1 , X2 in zip(X1list , X2list):
        lmono , Amono = helmert_trans_estim_matrixs_maker(X1,X2)
            
        l_stk.append(lmono)
        A_stk.append(Amono)
        
    A = np.vstack(A_stk)
    l = np.hstack(l_stk)
    
    if len(Weights) == 0:
        W = np.eye(len(l))
    else:
        W = np.repeat(np.array(Weights),3) * np.eye(len(l))
    
    N    = (A.T).dot(W).dot(A)
    AtWB = (A.T).dot(W).dot(l)
    
    Res = scipy.linalg.inv(N).dot(AtWB)
    
    return Res , A , l


def helmert_trans_apply(Xin,SevenParam_in):
    tx,ty,tz,rx,ry,rz,scal = SevenParam_in 
    
    R = np.array([[1.,rz,-ry],
                  [-rz,1.,rx],
                  [ry,-rx,1.]])
    
    T = np.array([tx,ty,tz])
    S = (1. + scal)
    
    typ=utils.get_type_smart(Xin)
    
    Xout = []
    for X1 in Xin:
        X2 = S * np.dot(R,X1) + T
        Xout.append(X2)
    
    Xout=typ(Xout)
    
    return Xout



def helmert_trans_estim_minimisation(X1in,X2in,HParam_apri=np.zeros(7),
                                     L1norm=True,tol=10**-9,full_output=False):
    
    def minimiz_helmert_fct(HParam_mini_in,X1in,X2in,L1norm_mini=L1norm):
        """
        This fct is the input for the scipy optimization fct
        """
        HParam_mini_wk = HParam_mini_in.copy()
        X12wrk = helmert_trans_apply(X1in,HParam_mini_wk)

        if not L1norm_mini: #return a L2 norm 
            SUM=np.sum((np.sum(np.power(X2in - X12wrk,2),axis=1)))
        else: #Return L1 norm
            SUM=np.sum(np.sum(np.abs(X2in - X12wrk),axis=1))
        return SUM
        
    
    
    RES = scipy.optimize.minimize(minimiz_helmert_fct,HParam_apri,
                                  (X1in,X2in,L1norm),
                                  method="Powell",tol=tol,
                                  options={"maxiter":100,
                                         'xtol':tol,
                                         'ftol':tol})
    
    if not full_output:
        return RES.x
    else:
        return RES
    
    
    

    
def helmert_trans_estim_minimisation_scalar(X1,X2,HParam_opti_apriori,
                                            L1norm=True,itera=2):
    """
    Estimates the Helmert parameters but based on a minimisation approach 
    between X1 and X2 (as suggested in the IGS combination software)
    
    NOT STABLE AVOID THE USE
    """
    
    print("WARN : unstable, avoid !!!!!")
    
    def minimiz_helmert_fct_scalar(hparam_mono_in,hparam_mono_id,
                            HParam_mini_in,X1in,X2in):
        """
        This fct is the input for the scipy optimization fct
        """
        HParam_mini_wk = HParam_mini_in.copy()
        HParam_mini_wk[hparam_mono_id] = hparam_mono_in
        X12wrk = helmert_trans_apply(X1in,HParam_mini_wk)

        if not L1norm: #return a L2 norm 
            return np.sum(np.sqrt(np.sum(np.power(X2in - X12wrk,2),axis=1)))
        else: #Return L1 norm
            return np.sum(np.sum(np.abs(X2in - X12wrk),axis=1))
        
    HParam_opti_wrk = HParam_opti_apriori.copy()

    for j in range(itera): #iter iterations
        for i in range(7): #7 Helmert Parameters:
            RES=scipy.optimize.minimize_scalar(minimiz_helmert_fct,
                                               args=(i,HParam_opti_wrk,X1,X2),
                                               tol=10**-20)
            HParam_opti_wrk[i] = RES.x
            
        if not j:
            HParam_opti_prev = HParam_opti_wrk.copy()
        else:
            print("Helmert Param minimisation iter",j+1,HParam_opti_wrk - HParam_opti_prev)
            HParam_opti_prev = HParam_opti_wrk.copy()
                        
    return HParam_opti_wrk





        
#### ASTRONOMY FUNCTION
def semi_major_axis_from_mean_motion(n):
    """
    source : https://space.stackexchange.com/questions/18289/how-to-get-semi-major-axis-from-tle
    """
    mu = 3.9860044189 * 10**14
    a  = (mu**(1./3.)) / ((2*n*np.pi/86400)**(2./3.))
    return a    

