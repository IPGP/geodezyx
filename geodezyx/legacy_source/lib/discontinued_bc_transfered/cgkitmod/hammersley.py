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
# $Id: hammersley.py,v 1.1 2006/02/17 18:56:09 mbaas Exp $

"""
  Hammersley & Halton point generation

  This module contains functions to generate points that are uniformly
  distributed and stochastic-looking on either a unit square or a unit
  sphere. The Hammersley point set is more uniform but is
  non-hierarchical, i.e.  for different n arguments you get an
  entirely new sequence. If you need hierarchical behavior you can use
  the Halton point set.

  This is a Python version of the implementation provided in:

  Tien-Tsin Wong, Wai-Shing Luk, Pheng-Ann Heng
  "Sampling with Hammersley and Halton points"
  Journal of Graphics Tools, Vol. 2, No. 2, 1997, pp. 9-24
  http://www.acm.org/jgt/papers/WongLukHeng97/
  http://www.cse.cuhk.edu.hk/~ttwong/papers/udpoint/udpoints.html

  The original C versions of these functions are distributed under
  the following license:
  
  (c) Copyright 1997, Tien-Tsin Wong
  ALL RIGHTS RESERVED
  Permission to use, copy, modify, and distribute this software for
  any purpose and without fee is hereby granted, provided that the above
  copyright notice appear in all copies and that both the copyright notice
  and this permission notice appear in supporting documentation,
 
  THE MATERIAL EMBODIED ON THIS SOFTWARE IS PROVIDED TO YOU "AS-IS"
  AND WITHOUT WARRANTY OF ANY KIND, EXPRESS, IMPLIED OR OTHERWISE,
  INCLUDING WITHOUT LIMITATION, ANY WARRANTY OF MERCHANTABILITY OR
  FITNESS FOR A PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHOR
  BE LIABLE TO YOU OR ANYONE ELSE FOR ANY DIRECT,
  SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY
  KIND, OR ANY DAMAGES WHATSOEVER, INCLUDING WITHOUT LIMITATION,
  LOSS OF PROFIT, LOSS OF USE, SAVINGS OR REVENUE, OR THE CLAIMS OF
  THIRD PARTIES, WHETHER OR NOT THE AUTHOR HAS BEEN
  ADVISED OF THE POSSIBILITY OF SUCH LOSS, HOWEVER CAUSED AND ON
  ANY THEORY OF LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE
  POSSESSION, USE OR PERFORMANCE OF THIS SOFTWARE.
"""

from math import pi, sqrt, cos, sin
from .cgtypes import *

# planeHammersley
def planeHammersley(n):
    """Yields n Hammersley points on the unit square in the xy plane.

    This function yields a sequence of n tuples (x,y) which
    represent a point on the unit square. The sequence of points for
    a particular n is always the same.  When n changes an entirely new
    sequence will be generated.
    
    This function uses a base of 2.
    """
    for k in range(n):
        u = 0
        p=0.5
        kk = k
        while kk>0:
            if kk & 1:
                u += p
            p *= 0.5
            kk >>= 1
        v = (k+0.5)/n
        yield (u, v)

# sphereHammersley
def sphereHammersley(n):
    """Yields n Hammersley points on the unit sphere.

    This function yields n vec3 objects representing points on the
    unit sphere. The sequence of points for a particular n is always
    the same.  When n changes an entirely new sequence will be
    generated.

    This function uses a base of 2.    
    """
    for k in range(n):
        t = 0
        
        p = 0.5
        kk = k
        while kk>0:
            if kk & 1:
                t += p
            p *= 0.5
            kk >>= 1
            
        t = 2.0*t - 1.0
        phi = (k+0.5)/n
        phirad = phi*2.0*pi
        st = sqrt(1.0-t*t)
        yield vec3(st*cos(phirad), st*sin(phirad), t)

# planeHalton
def planeHalton(n=None, p2=3):
    """Yields a sequence of Halton points on the unit square in the xy plane.

    This function yields a sequence of two floats (x,y) which
    represent a point on the unit square. The number of points to
    generate is given by n. If n is set to None, an infinite number of
    points is generated and the caller has to make sure the loop stops
    by checking some other critera.  The sequence of generated points
    is always the same, no matter what n is (i.e. the first n elements
    generated by the sequence planeHalton(n+1) is identical to the
    sequence planeHalton(n)).

    This function uses 2 as its first prime base whereas the second
    prime p2 (which must be a prime number) can be provided by the user.
    """
    
    k = 0
    while k+1!=n:
        u = 0
        p = 0.5
        kk = k
        while kk>0:
            if kk & 1:
                u += p
            p *= 0.5
            kk >>= 1
        v = 0
        ip = 1.0/p2
        p = ip
        kk = k
        while kk>0:
            a = kk % p2
            if a!=0:
                v += a*p
            p *= ip
            kk = int(kk/p2)
        yield (u,v)
        k += 1

# sphereHalton
def sphereHalton(n=None, p2=3):
    """Yields a sequence of Halton points on the unit sphere.

    This function yields a sequence of vec3 objects representing
    points on the unit sphere. The number of points to generate is
    given by n. If n is set to None, an infinite number of points is
    generated and the caller has to make sure the loop stops by
    checking some other critera. The sequence of generated points is
    always the same, no matter what n is (i.e. the first n elements
    generated by the sequence sphereHalton(n+1) is identical to the
    sequence sphereHalton(n)).

    This function uses 2 as its first prime base whereas the second
    base p2 (which must be a prime number) can be provided by the user.
    """
    k = 0
    while k+1!=n:
        t = 0
        p = 0.5
        kk = k
        while kk>0:
            if kk & 1:
                t += p
            p *= 0.5
            kk >>= 1
        t = 2.0*t - 1.0
        st = sqrt(1.0-t*t)
        phi = 0
        ip = 1.0/p2
        p = ip
        kk = k
        while kk>0:
            a = kk % p2
            if a!=0:
                phi += a*p
            p *= ip
            kk = int(kk/p2)
        phirad = phi*4.0*pi
        yield vec3(st*cos(phirad), st*sin(phirad), t)
        k += 1
    
