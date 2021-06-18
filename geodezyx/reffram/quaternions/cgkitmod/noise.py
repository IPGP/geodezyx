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
# $Id: noise.py,v 1.3 2006/02/14 19:29:39 mbaas Exp $

"""This module contains various noise functions.
"""

from ._core import noise, snoise, cellnoise, scellnoise, fBm, turbulence
from ._core import pnoise as _pnoise
from ._core import spnoise as _spnoise
from ._core import vnoise, vsnoise, vcellnoise, vscellnoise, vfBm, vturbulence
from ._core import vpnoise as _vpnoise
from ._core import vspnoise as _vspnoise
from .cgtypes import vec3

# pnoise
def pnoise(*args):
    """pnoise(point, period) / pnoise(point, t, pperiod, tperiod) -> float

    Periodic noise function. Basically this is the same than noise
    but with a periodic return value: pnoise(point) = pnoise(point+period).
    The time value can be either part of the point or it can be
    specified separately. The point and period must always have the
    same dimension. The return value is in the range from 0 to 1.
    """
    
    n = len(args)
    if n==2:
        v, pv = args
        try:
            m = len(v)
        except:
            return _pnoise(v, 0, pv, 1)
        if m==1:
            return _pnoise(v[0], 0, pv[0], 1)
        elif m==2:
            x,y = v
            px,py = pv
            return _pnoise(x,y,px,py)
        elif m==3:
            x,y,z = v
            px,py,pz = pv
            return _pnoise(x,y,z,px,py,pz)
        elif m==4:
            x,y,z,t = v
            px,py,pz,pt = pv
            return _pnoise(x,y,z,t,px,py,pz,pt)
        else:
            raise ValueError("Invalid arguments")
        
    elif n==4:
        v,t,pv,pt = args
        try:
            m = len(v)
        except:
            return _pnoise(v, t, pv, pt)
        if m==1:
            return _pnoise(v[0], t, pv[0], pt)
        elif m==2:
            x,y = v
            px,py = pv
            return _pnoise(x,y,t,px,py,pt)
        elif m==3:
            x,y,z = v
            px,py,pz = pv
            return _pnoise(x,y,z,t,px,py,pz,pt)
        else:
            raise ValueError("Invalid arguments")
    
    else:
        raise TypeError("only 2 or 4 arguments allowed")


# spnoise
def spnoise(*args):
    """spnoise(point, period) / spnoise(point, t, pperiod, tperiod) -> float

    Signed periodic noise function. The return value is in the range
    from -1 to 1. A call to spnoise(args) is equivalent to
    2*pnoise(args)-1.
    """
    
    n = len(args)
    if n==2:
        v, pv = args
        try:
            m = len(v)
        except:
            return _spnoise(v, 0, pv, 1)
        if m==1:
            return _spnoise(v[0], 0, pv[0], 1)
        elif m==2:
            x,y = v
            px,py = pv
            return _spnoise(x,y,px,py)
        elif m==3:
            x,y,z = v
            px,py,pz = pv
            return _spnoise(x,y,z,px,py,pz)
        elif m==4:
            x,y,z,t = v
            px,py,pz,pt = pv
            return _spnoise(x,y,z,t,px,py,pz,pt)
        else:
            raise ValueError("Invalid arguments")
        
    elif n==4:
        v,t,pv,pt = args
        try:
            m = len(v)
        except:
            return _spnoise(v, t, pv, pt)
        if m==1:
            return _spnoise(v[0], t, pv[0], pt)
        elif m==2:
            x,y = v
            px,py = pv
            return _spnoise(x,y,t,px,py,pt)
        elif m==3:
            x,y,z = v
            px,py,pz = pv
            return _spnoise(x,y,z,t,px,py,pz,pt)
        else:
            raise ValueError("Invalid arguments")
    
    else:
        raise TypeError("only 2 or 4 arguments allowed")


# vpnoise
def vpnoise(*args):
    """vpnoise(point, period) / vpnoise(point, t, pperiod, tperiod) -> noiseval

    Periodic noise function. Basically this is the same than vnoise
    but with a periodic return value: vpnoise(point) = vpnoise(point+period).
    The time value can be either part of the point or it can be
    specified separately. The point and period must always have the
    same dimension. The components of the return value are in the range
    from 0 to 1.
    """
    
    n = len(args)
    if n==2:
        v, pv = args
        try:
            m = len(v)
        except:
            return _vpnoise(v, pv)
        if m==1:
            return _vpnoise(v[0], pv[0])
        elif m==2:
            x,y = v
            px,py = pv
            return _vpnoise(x,y,px,py)
        elif m==3:
            x,y,z = v
            px,py,pz = pv
            return _vpnoise(x,y,z,px,py,pz)
        elif m==4:
            x,y,z,t = v
            px,py,pz,pt = pv
            return _vpnoise(x,y,z,t,px,py,pz,pt)
        else:
            raise ValueError("Invalid arguments")
        
    elif n==4:
        v,t,pv,pt = args
        try:
            m = len(v)
        except:
            return _vpnoise(v, t, pv, pt)
        if m==1:
            return _vpnoise(v[0], t, pv[0], pt)
        elif m==2:
            x,y = v
            px,py = pv
            return _vpnoise(x,y,t,px,py,pt)
        elif m==3:
            x,y,z = v
            px,py,pz = pv
            res = _vpnoise(x,y,z,t,px,py,pz,pt)
            return vec3(res.x, res.y, res.z)
        else:
            raise ValueError("Invalid arguments")
    
    else:
        raise TypeError("only 2 or 4 arguments allowed")


# vspnoise
def vspnoise(*args):
    """vspnoise(point, period) / vspnoise(point, t, pperiod, tperiod) -> noiseval

    Periodic noise function. Basically this is the same than vsnoise
    but with a periodic return value: vspnoise(point) = vspnoise(point+period).
    The time value can be either part of the point or it can be
    specified separately. The point and period must always have the
    same dimension. The components of the return value are in the range
    from -1 to 1.
    """
    
    n = len(args)
    if n==2:
        v, pv = args
        try:
            m = len(v)
        except:
            return _vspnoise(v, pv)
        if m==1:
            return _vspnoise(v[0], pv[0])
        elif m==2:
            x,y = v
            px,py = pv
            return _vspnoise(x,y,px,py)
        elif m==3:
            x,y,z = v
            px,py,pz = pv
            return _vspnoise(x,y,z,px,py,pz)
        elif m==4:
            x,y,z,t = v
            px,py,pz,pt = pv
            return _vspnoise(x,y,z,t,px,py,pz,pt)
        else:
            raise ValueError("Invalid arguments")
        
    elif n==4:
        v,t,pv,pt = args
        try:
            m = len(v)
        except:
            return _vspnoise(v, t, pv, pt)
        if m==1:
            return _vspnoise(v[0], t, pv[0], pt)
        elif m==2:
            x,y = v
            px,py = pv
            return _vspnoise(x,y,t,px,py,pt)
        elif m==3:
            x,y,z = v
            px,py,pz = pv
            res = _vspnoise(x,y,z,t,px,py,pz,pt)
            return vec3(res.x, res.y, res.z)
        else:
            raise ValueError("Invalid arguments")
    
    else:
        raise TypeError("only 2 or 4 arguments allowed")


