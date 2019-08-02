#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 17:15:57 2019

@author: psakicki
"""

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

    if not genefun.is_iterable(alphal):
        alphal = np.array([alphal])
        betal = np.array([betal])
        gammal = np.array([gammal])
        boolnotiterable = True
    else:
        boolnotiterable = False

    pointlout = []
    R_out = []


    for pt in pointlin:

        if not genefun.is_iterable(pt) or len(pt) != 3:
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

    if not genefun.is_iterable(X):
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
    Xfit , Yfit = geok.linear_reg_getvalue(Xobs,a,b)
    r2 = R2_calc(Yobs,Yfit)
    return r2
