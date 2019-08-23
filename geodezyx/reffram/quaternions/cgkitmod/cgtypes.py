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
# $Id: cgtypes.py,v 1.2 2006/02/14 19:29:39 mbaas Exp $

## \file cgtypes.py
## Contains the wrapped types vec3, vec4, mat3, mat4 and quat.

"""This module contains the basic types for doing 3D computations.

Classes:

- vec3: 3D vector (x,y,z)
- vec4: 4D vector (x,y,z,w) or (x,y,z,t)
- mat3: 3x3 matrix
- mat4: 4x4 matrix
- quat: Quaternion (w,x,y,z)

Functions:

- slerp(): Spherical linear interpolation between two quaternions.
- squad(): Spherical cubic interpolation.
- getEpsilon(): Return the threshold used for float comparisons
- setEpsilon(): Set a new threshold used for float comparisons
"""

from . import _core
from ._core import slerp, squad

# vec3
class vec3(_core.vec3):
    """A 3-dimensional vector of floats.

    This class is derived from the vec3 class in the _core module and
    adds some missing stuff that's better done in Python.
    """
    
    def __init__(self, *args):
        
        if len(args)==1:
            T = type(args[0])
            # scalar
            if T in [float, int]:
                _core.vec3.__init__(self, float(args[0]))
            # vec3  
            elif isinstance(args[0], vec3) or isinstance(args[0], _core.vec3):
                _core.vec3.__init__(self, *args)
            # Tuple/List
#            elif T in [tuple, list]:
#                _core.vec3.__init__(self, *args[0])
#                if len(args[0])==0:
#                    self.x = self.y = self.z = 0.0
#                elif len(args[0])==1:
#                    self.x = self.y = self.z = args[0][0]
#                elif len(args[0])==2:
#                    self.x, self.y = args[0]
#                    self.z         = 0.0
#                elif len(args[0])==3:
#                    self.x, self.y, self.z = args[0]
#                else:
#                    raise TypeError("vec3() takes at most 3 arguments")
            # String
            elif T is str:
                s=args[0].replace(","," ").replace("  "," ").strip().split(" ")
                if s==[""]:
                    s=[]
                f = [float(x) for x in s]
                _core.vec3.__init__(self, *f)
            # error
            else:
                try:
                    # Try to treat the argument as a sequence of 3 floats...
                    lst = list(args[0])
                    _core.vec3.__init__(self, *lst)
                except:
                    raise TypeError("vec3() arg can't be converted to vec3")

        else:
            # everything else that is not just one argument...
            _core.vec3.__init__(self, *args)

# vec4
class vec4(_core.vec4):
    """A 4-dimensional vector of floats.

    This class is derived from the vec4 class in the _core module and
    adds some missing stuff that's better done in Python.
    """
    
    def __init__(self, *args):
        
        if len(args)==1:
            T = type(args[0])
            # scalar
            if T in [float, int]:
                _core.vec4.__init__(self, float(args[0]))
            # vec4
            elif isinstance(args[0], vec4) or isinstance(args[0], _core.vec4):
                _core.vec4.__init__(self, *args)
            # String
            elif T is str:
                s=args[0].replace(","," ").replace("  "," ").strip().split(" ")
                if s==[""]:
                    s=[]
                f = [float(x) for x in s]
                _core.vec4.__init__(self, *f)
            # error
            else:
                try:
                    # Try to treat the argument as a sequence of floats...
                    lst = list(args[0])
                    _core.vec4.__init__(self, *lst)
                except:
                    raise TypeError("vec4() arg can't be converted to vec4")

        else:
            # everything else that is not just one argument...
            _core.vec4.__init__(self, *args)

# mat3
class mat3(_core.mat3):
    """A 3x3 matrix of floats.

    This class is derived from the mat3 class in the _core module and
    adds some missing stuff that's better done in Python.
    """

    def __init__(self, *args):
        
        if len(args)==1:
            T = type(args[0])
            # scalar
            if T in [float, int]:
                _core.mat3.__init__(self, float(args[0]))
            # mat3
            elif isinstance(args[0], mat3) or isinstance(args[0], _core.mat3):
                _core.mat3.__init__(self, *args)
            # Tuple/List
            elif T in [tuple, list]:
                _core.mat3.__init__(self, mat3(*args[0]))
            # String
            elif T is str:
                s=args[0].replace(","," ").replace("  "," ").strip().split(" ")
                if s==[""]:
                    s=[]
                f = [float(x) for x in s]
                _core.mat3.__init__(self, *f)
            # error
            else:
                raise TypeError("mat3() arg can't be converted to mat3")

        elif len(args)==3:
            try:
                c1,c2,c3 = args
                a,d,g = c1
                b,e,h = c2
                c,f,i = c3
                _core.mat3.__init__(self, a,b,c,d,e,f,g,h,i)
            except:
                raise TypeError("mat3() arg can't be converted to mat3")
        else:
            # everything else that is not just one argument...
            _core.mat3.__init__(self, *args)


    def __getitem__(self, key):
        
        T = type(key)
        if T is int:
            return _core.mat3.__getitem__(self, key)
        elif T is tuple:
            if len(key)!=2:
                raise ValueError("index tuple must be a 2-tuple")
            i,j = key
            if i<0 or i>2 or j<0 or j>2:
                raise IndexError("index out of range")
            return self[j][i]
        else:
            raise TypeError("index must be integer or 2-tuple")

# mat4
class mat4(_core.mat4):
    """A 4x4 matrix of floats.

    This class is derived from the mat4 class in the _core module and
    adds some missing stuff that's better done in Python.
    """

    def __init__(self, *args):
        
        if len(args)==1:
            T = type(args[0])
            # scalar
            if T in [float, int]:
                _core.mat4.__init__(self, float(args[0]))
            # mat4
            elif isinstance(args[0], mat4) or isinstance(args[0], _core.mat4):
                _core.mat4.__init__(self, *args)
            # Tuple/List
            elif T in [tuple, list]:
                _core.mat4.__init__(self, mat4(*args[0]))
            # String
            elif T is str:
                s=args[0].replace(","," ").replace("  "," ").strip().split(" ")
                if s==[""]:
                    s=[]
                f = [float(x) for x in s]
                _core.mat4.__init__(self, mat4(*f))
            # error
            else:
                raise TypeError("mat4() arg can't be converted to mat4")

        elif len(args)==4:
            try:
                c1,c2,c3,c4 = args
                _core.mat4.__init__(self)
                self.setColumn(0,vec4(c1))
                self.setColumn(1,vec4(c2))
                self.setColumn(2,vec4(c3))
                self.setColumn(3,vec4(c4))
            except:
                raise TypeError("mat4() arg can't be converted to mat4")
        elif len(args)==16:
            _core.mat4.__init__(self)
            self.setRow(0, vec4(args[0:4]))
            self.setRow(1, vec4(args[4:8]))
            self.setRow(2, vec4(args[8:12]))
            self.setRow(3, vec4(args[12:16]))
        else:
            # everything else that is not just one argument...
            _core.mat4.__init__(self, *args)

    def __getitem__(self, key):
        
        T = type(key)
        if T is int:
            return _core.mat4.__getitem__(self, key)
        elif T is tuple:
            if len(key)!=2:
                raise ValueError("index tuple must be a 2-tuple")
            i,j = key
            if i<0 or i>3 or j<0 or j>3:
                raise IndexError("index out of range")
            return self[j][i]
        else:
            raise TypeError("index must be integer or 2-tuple")

# quat
class quat(_core.quat):
    """A quaternion class.

    This class is derived from the quat class in the _core module and
    adds some missing stuff that's better done in Python.
    """
    
    def __init__(self, *args):
        
        if len(args)==1:
            T = type(args[0])
            # scalar
            if T in [float, int]:
                _core.quat.__init__(self, float(args[0]))
            # quat
            elif isinstance(args[0], quat) or isinstance(args[0], _core.quat):
                _core.quat.__init__(self, *args)
            # mat3/mat4
            elif isinstance(args[0], mat3) or isinstance(args[0], _core.mat3) or isinstance(args[0], mat4) or isinstance(args[0], _core.mat4):
                _core.quat.__init__(self)
                self.fromMat(args[0])
            # String
            elif T is str:
                s=args[0].replace(","," ").replace("  "," ").strip().split(" ")
                if s==[""]:
                    s=[]
                f = [float(x) for x in s]
                _core.quat.__init__(self, *f)
            # error
            else:
                try:
                    # Try to treat the argument as a sequence of floats...
                    lst = list(args[0])
                    _core.quat.__init__(self, *lst)
                except:
                    raise TypeError("quat() arg can't be converted to quat")

        # 2 arguments (angle & axis)
        elif len(args)==2:
            _core.quat.__init__(self)
            angle, axis = args
            self.fromAngleAxis(angle,axis)

        else:
            # everything else that is not just one argument...
            _core.quat.__init__(self, *args)

    def __pow__(self, other, modulo=None):
        """Return self**q."""
        if modulo!=None:
            raise TypeError("unsupported operation")
        q = quat(other)
        return (q*self.log()).exp()


def getEpsilon():
    """Return the threshold used for float comparisons.
    """
    return _core._getEpsilon()

def setEpsilon(eps):
    """Set the threshold used for float comparisons.
    """
    return _core._setEpsilon(eps)
