# ***** BEGIN LICENSE BLOCK *****
# Version: MPL 1.1/GPL 2.0/LGPL 2.1
#
# The contents of this file are subject to the Mozilla Public License Version
# 1.1 (the "License"); you may not use this file except in compliance with
# the License. You may obtain a copy of the License at
# http://www.mozilla.org/MPL/
#
# Software distributed under the License is distributed on an "AS IS" basis,
# WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License
# for the specific language governing rights and limitations under the
# License.
#
# The Original Code is the Python Computer Graphics Kit.
#
# The Initial Developer of the Original Code is Matthias Baas.
# Portions created by the Initial Developer are Copyright (C) 2004
# the Initial Developer. All Rights Reserved.
#
# Contributor(s):
#
# Alternatively, the contents of this file may be used under the terms of
# either the GNU General Public License Version 2 or later (the "GPL"), or
# the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
# in which case the provisions of the GPL or the LGPL are applicable instead
# of those above. If you wish to allow use of your version of this file only
# under the terms of either the GPL or the LGPL, and not to allow others to
# use your version of this file under the terms of the MPL, indicate your
# decision by deleting the provisions above and replace them with the notice
# and other provisions required by the GPL or the LGPL. If you do not delete
# the provisions above, a recipient may use your version of this file under
# the terms of any one of the MPL, the GPL or the LGPL.
#
# ***** END LICENSE BLOCK *****
# -------------------------------------------------------------
# The RenderMan (R) Interface Procedures and Protocol are:
# Copyright 1988, 1989, 2000, Pixar
# All Rights Reserved
#
# RenderMan (R) is a registered trademark of Pixar
# -------------------------------------------------------------
# $Id: sl.py,v 1.2 2006/02/14 19:29:39 mbaas Exp $

"""RenderMan Shading Language functionality.

This module provides some of the functionality from the RenderMan
Shading Language. Those functions that require an actual rendering
context (surfaces, light sources, etc.) are not supported.

Most of the functions can be used just like in the Shading
Language. An exception are those functions whose return type is
dependant on the context as it is the case with random() or
noise(). Here, those functions have to be prepended with the return
type, for example float_random() or point_noise() (that is, the cast
is part of the name).
"""

import math, random, string, re
from math import acos, asin, ceil, cos, exp, floor, pow, sin, sqrt, tan
from . import noise
from .cgtypes import vec3 as _vec3

# Builtin functions
abs = abs
max = max
min = min
round = round

PI = math.pi

def _tovec3(arg):
    try:
        a = len(arg)
    except:
        a = 1
    if a==1:
        return _vec3(arg,0.0,0.0)
    elif a==2:
        x,y=arg
        return _vec3(x,y,0.0)
    elif a==3:
        x,y,z=arg
        return _vec3(x,y,z)
    elif a==4:
        x,y,z,t=arg
        return _vec3(x,y,z)
    else:
        return _vec3()
    

def atan(*args):
    """Returns the arc tangent.

    With one argument it math.atan() is called, with two arguments
    math.atan2() is called.
    """
    
    if len(args)==1:
        return math.atan(args[0])
    elif len(args)==2:
        return math.atan2(args[0],args[1])
    else:
        raise TypeError("atan only takes 1 or 2 arguments.")

def log(*args):
    """Returns the natural logarithm of x or the logarithm to the specified base."""

    if len(args)==1:
        return math.log(args[0])
    elif len(args)==2:
        return math.log(args[0])/math.log(args[1])
    else:
        raise TypeError("log only takes 1 or 2 arguments.")
    

# clamp
def clamp(a, amin, amax):
    """Returns amin if a < amin, amax if a > amax, otherwise a."""
    return min(max(a,amin), amax)

# degrees
def degrees(rad):
    """Convert from radians to degrees."""
    return rad*180.0/PI

# radians
def radians(deg):
    """Convert from degrees to radians."""
    return deg*PI/180.0

def Du(p):
    pass

def Dv(p):
    pass

def Deriv(num, den):
    pass

def filterstep(edge, s1):
    pass

def inversesqrt(x):
    """Returns 1/sqrt(x)."""  
    return 1.0/sqrt(x)

def mix(val0, val1, t):
    """Mix two values.

    For t=0 the value val0 is returned, for t=1 the value val1 is
    returned. For values of t between 0 and 1 a linearly interpolated
    value is returned.
    """
    
    return (1.0-t)*val0 + t*val1

def mod(a,b):
    """Returns a%b. This is just an equivalent for the %-operator."""
    return a%b

def float_noise(*args):
    """Returns a float value which is a (pseudo) random function of its arguments.

    This function is imported from the noise module.
    """
    
    la = len(args)
    if la==1:
        return noise.noise(args[0])
    elif la==2:
        return noise.noise(args[0],args[1])
    elif la==3:
        return noise.noise(args[0],args[1],args[2])
    elif la==4:
        return noise.noise(args[0],args[1],args[2],args[3])
    else:
        raise TypeError("the function takes between 1 and 4 arguments (%s given)"%(la))


def point_noise(*args):
    """Returns a point whose value is a (pseudo) random function of its arguments."""
    
    la = len(args)
    if la==1:
        try:
            a = len(args[0])
        except:
            a = 1
        if a==1:
            return noise.vnoise(args[0],0,0)
        elif a==2:
            return noise.vnoise(args[0],0)
        elif a==3:
            return noise.vnoise(args[0])
        else:
            raise ValueError("arg1: invalid argument length")
    elif la==2:
        try:
            a = len(args[0])
        except:
            a = 1
        if a==1:
            return noise.vnoise(args[0],args[1],0)
        elif a==3:
            return noise.vnoise(args[0],args[1])
        else:
            raise ValueError("arg1: invalid argument length")
    elif la==3:
        return noise.vnoise(args[0],args[1],args[2])
    elif la==4:
        x,y,z,t = noise.vnoise(args[0],args[1],args[2],args[3])
        return _vec3(x,y,z)
    else:
        raise TypeError("the function takes between 1 and 4 arguments (%s given)"%(la))
   
color_noise = point_noise
vector_noise = point_noise

def float_pnoise(*args):
    """Returns a float value which is a periodic (pseudo) random function of its arguments.

    This function is imported from the noise module."""
    
    la = len(args)
    try:
        a = len(args[0])
    except:
        a = 1
        
    if la==2:
        if a==1:
            return noise.pnoise((args[0],),(args[1],))
        else:
            return noise.pnoise(args[0],args[1])
    elif la==4:
        if a==1:
            return noise.pnoise((args[0],args[1]),(args[2],args[3]))
        else:
            return noise.pnoise(args[0],args[1],args[2],args[3])
    else:
        raise TypeError("the function takes between 1 and 4 arguments (%s given)"%(la))

def point_pnoise(*args):
    """Returns a point whose value is a periodic (pseudo) random function of its arguments."""
    
    la = len(args)
    try:
        a = len(args[0])
    except:
        a = 1
        
    if la==2:
        if a==1:
            res = noise.vpnoise((args[0],),(args[1],))
        else:
            res = noise.vpnoise(args[0],args[1])
    elif la==4:
        if a==1:
            res = noise.vpnoise((args[0],args[1]),(args[2],args[3]))
        else:
            res = noise.vpnoise(args[0],args[1],args[2],args[3])
    else:
        raise TypeError("the function takes between 1 and 4 arguments (%s given)"%(la))

    return _tovec3(res)

color_pnoise = point_pnoise
vector_pnoise = point_pnoise

def float_cellnoise(*args):
    """Returns a float value which is a (pseudo) random function of its arguments.

    The return value is constant between integer lattice points. This
    function is imported from the noise module.
    """
    
    la = len(args)
    if la==1:
        return noise.cellnoise(args[0])
    elif la==2:
        return noise.cellnoise(args[0],args[1])
    elif la==3:
        return noise.cellnoise(args[0],args[1],args[2])
    elif la==4:
        return noise.cellnoise(args[0],args[1],args[2],args[3])
    else:
        raise TypeError("the function takes between 1 and 4 arguments (%s given)"%(la))

def point_cellnoise(*args):
    """Returns a point whose value is a (pseudo) random function of its arguments.

    The return value is constant between integer lattice points.
    """
    
    la = len(args)
    if la==1:
        try:
            a = len(args[0])
        except:
            a = 1
        if a==1:
            return noise.vcellnoise(args[0],0,0)
        elif a==2:
            return noise.vcellnoise(args[0],0)
        elif a==3:
            return noise.vcellnoise(args[0])
        else:
            raise ValueError("arg1: invalid argument length")
    elif la==2:
        try:
            a = len(args[0])
        except:
            a = 1
        if a==1:
            return noise.vcellnoise(args[0],args[1],0)
        elif a==3:
            return noise.vcellnoise(args[0],args[1])
        else:
            raise ValueError("arg1: invalid argument length")
    elif la==3:
        return noise.vcellnoise(args[0],args[1],args[2])
    elif la==4:
        x,y,z,t = noise.vcellnoise(args[0],args[1],args[2],args[3])
        return _vec3(x,y,z)
    else:
        raise TypeError("the function takes between 1 and 4 arguments (%s given)"%(la))

color_cellnoise = point_cellnoise
vector_cellnoise = point_cellnoise

def float_random():
    """Return a random number between 0 and 1.

    This call is equivalent to random.random()."""
    
    return random.random()

def color_random():
    """Return a color whose componenets are a random number between 0 and 1.

    The function actually returns a vec3."""
    
    return _vec3(random.random(), random.random(), random.random())

def point_random():
    """Return a point (a vec3) whose componenets are a random number between 0 and 1."""
    return _vec3(random.random(), random.random(), random.random())

def sign(x):
    """Returns -1 with a negative argument, +1 with a positive argument, and 0 if its argument is zero."""
    
    if x<0:
        return -1
    elif x>0:
        return 1
    else:
        return 0

def smoothstep(min, max, x):
    """Returns the value of a smooth step function.

    Returns 0 if x < min, 1 if x > max, and performs a smooth Hermite
    interpolation between 0 and 1 in the interval min to max.
    """
    
    if x<min:
        return 0.0
    if x>max:
        return 1.0
    x = (x-min)/(max-min)
    return x*x*(3.0-2.0*x)

def spline(x, knots):
    """Return the value of a spline function.

    Fits a spline to the control points given and returns the value at
    t which ranges from 0 to 1. At least four control points must
    always be given.
    """

    nknots = len(knots)
    nspans = nknots-3
 
    if nspans<1:
        raise ValueError("spline(): there must be at least 4 control points (%s given)"%nknots)

    x = clamp(x, 0.0, 1.0)*nspans
    span = int(x)
    if span>=nknots-3:
        span = nknots-4
    x    -= span
    knot0, knot1, knot2, knot3 = knots[span:span+4]

    c3 = -0.5*knot0 + 1.5*knot1 - 1.5*knot2 + 0.5*knot3
    c2 =      knot0 - 2.5*knot1 + 2.0*knot2 - 0.5*knot3
    c1 = -0.5*knot0 + 0.5*knot2
    c0 =      knot1

    return ((c3*x + c2)*x + c1)*x + c0


def step(min, x):
    """Returns 0 if x < min, otherwise 1."""
    if x<min:
        return 0.0
    else:
        return 1.0


def distance(P1,P2):
    """Returns the distance between two points.

    The arguments should be of type vec3."""
    return (P2-P1).length()

def ptlined(P0,P1,Q):
    """Returns the distance between a point and a line segment.

    The arguments should be of type vec3.    
    """
    a = P1-P0
    b = Q-P0
    
    x  = a*b
    if x<=0:
        return b.length()
    
    aa = a*a
    if x>=aa:
        return (Q-P1).length()
    
    return sqrt(b*b-(x*x/aa))

def faceforward(N,I,Nref):
    """Flips N so that it faces in the direction opposite to I."""
    return sign(-I*Nref)*N

def length(v):
    """Returns the length of a vector.

    This is equivalent to calling v.length().
    """
    return v.length()

def normalize(v):
    """Returns a unit vector in the direction of v.

    This is equivalent to calling v.normalize().
    """
    return v.normalize()

def reflect(I, N):
    """Returns the reflection vector given an incident direction I and a normal N.

    This is equivalent to calling I.reflect(N).
    """
    return I.reflect(N)

def refract(I, N, eta):
    """Returns the transmitted vector.

    Returns the transmitted vector given an incident direction I, the
    normal vector N and the relative index of refraction eta. This is
    equivalent to calling I.refract(N, eta).
    """
    
    return I.refract(N,eta)

def xcomp(P):
    """Return the x component of p.

    This is equivalent to p.x.
    """
    return P.x

def ycomp(P):
    """Return the y component of p.

    This is equivalent to p.y.
    """
    return P.y

def zcomp(P):
    """Return the z component of p.

    This is equivalent to p.z.
    """
    return P.z

def setxcomp(P,x):
    """Set the x component of p.

    This is equivalent to p.x = x."""
    P.x = x

def setycomp(P,y):
    """Set the y component of p.

    This is equivalent to p.y = y."""
    P.y = y

def setzcomp(P,z):
    """Set the z component of p.

    This is equivalent to p.z = z."""
    P.z = z

def comp(c, index):
    """Get an individual color component.

    This is equivalent to c[index]."""
    return c[index]

def setcomp(c, index, value):
    """Set an individual color component.

    This is equivalent to c[index] = value."""
    c[index]=value

def concat(*args):
    """Returns a concatenated string."""
    return "".join(args)
    
def match(pattern, subject):
    """String pattern matching."""
    return re.search(pattern, subject)!=None

def format(pattern, *args, **keyargs):
    """Returns a formatted string (similar to the C function sprintf())."""
    res=""
    while 1:
        n=pattern.find("%")
        if n==-1:
            res+=pattern
            break
        res+=pattern[0:n]+"%s"
        pattern=pattern[n+2:]

    if "args" in keyargs:
        args=keyargs["args"]

    return res % args

def printf(pattern, *args):
    """Prints the values of the specified variables."""
    print((format(pattern,args=args)))

######################################################################

if __name__=="__main__":

    printf("Hallo s=%f c=%c",0.5,(1,2,3))
    
