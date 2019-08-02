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
# $Id: slots.py,v 1.7 2005/08/28 19:55:40 mbaas Exp $

## \file slots.py
## Contains the standard slot classes.

#from _core import ISlot, DoubleSlot, IntSlot, BoolSlot, PySlot
#from _core import Vec3Slot, Vec4Slot, Mat3Slot, Mat4Slot
from .cgtypes import vec3, vec4, mat3, mat4, quat
from . import _core
from ._core import Dependent, UserSizeConstraint, LinearSizeConstraint
from ._core import ISlot, IArraySlot

# Factory functions for the individual slots:

def DoubleSlot(value=0.0, flags=0):
    return _core.DoubleSlot(value, flags)

def IntSlot(value=0, flags=0):
    return _core.IntSlot(value, flags)

def BoolSlot(value=False, flags=0):
    return _core.BoolSlot(value, flags)

def StrSlot(value="", flags=0):
    return _core.StrSlot(value, flags)

def PySlot(value=None, flags=0):
    return _core.PySlot(value, flags)

def Vec3Slot(value=vec3(), flags=0):
    return _core.Vec3Slot(value, flags)

def Vec4Slot(value=vec4(), flags=0):
    return _core.Vec4Slot(value, flags)

def Mat3Slot(value=mat3(), flags=0):
    return _core.Mat3Slot(value, flags)

def Mat4Slot(value=mat4(), flags=0):
    return _core.Mat4Slot(value, flags)

def QuatSlot(value=quat(), flags=0):
    return _core.QuatSlot(value, flags)



def DoubleArraySlot(multiplicity=1, constraint=None):
    return _core.DoubleArraySlot(multiplicity, constraint)

def IntArraySlot(multiplicity=1, constraint=None):
    return _core.IntArraySlot(multiplicity, constraint)

def BoolArraySlot(multiplicity=1, constraint=None):
    return _core.BoolArraySlot(multiplicity, constraint)

def StrArraySlot(multiplicity=1, constraint=None):
    return _core.StrArraySlot(multiplicity, constraint)

def Vec3ArraySlot(multiplicity=1, constraint=None):
    return _core.Vec3ArraySlot(multiplicity, constraint)

def Vec4ArraySlot(multiplicity=1, constraint=None):
    return _core.Vec4ArraySlot(multiplicity, constraint)

#def Mat3ArraySlot(multiplicity=1, constraint=None):
#    return _core.Mat3ArraySlot(multiplicity, constraint)

def Mat4ArraySlot(multiplicity=1, constraint=None):
    return _core.Mat4ArraySlot(multiplicity, constraint)

#def QuatArraySlot(multiplicity=1, constraint=None):
#    return _core.QuatArraySlot(multiplicity, constraint)


# Test

class ProceduralDoubleSlot(_core.DoubleSlot):
    def __init__(self, proc):
        _core.DoubleSlot.__init__(self, 0.0, 2) # 2 = NO_INPUT_CONNECTIONS
        self._proc = proc

    def computeValue(self):
        self._value = self._proc()

class ProceduralIntSlot(_core.IntSlot):
    def __init__(self, proc):
        _core.IntSlot.__init__(self, 0, 2) # 2 = NO_INPUT_CONNECTIONS
        self._proc = proc

    def computeValue(self):
        self._value = self._proc()

class ProceduralVec3Slot(_core.Vec3Slot):
    def __init__(self, proc):
        _core.Vec3Slot.__init__(self, vec3(), 2) # 2 = NO_INPUT_CONNECTIONS
        self._proc = proc

    def computeValue(self):
        self._value = self._proc()

class ProceduralVec4Slot(_core.Vec4Slot):
    def __init__(self, proc):
        _core.Vec4Slot.__init__(self, vec4(), 2) # 2 = NO_INPUT_CONNECTIONS
        self._proc = proc

    def computeValue(self):
        self._value = self._proc()

class ProceduralMat3Slot(_core.Mat3Slot):
    def __init__(self, proc):
        _core.Mat3Slot.__init__(self, mat3(1), 2) # 2 = NO_INPUT_CONNECTIONS
        self._proc = proc

    def computeValue(self):
        self._value = self._proc()

class ProceduralMat4Slot(_core.Mat4Slot):
    def __init__(self, proc):
        _core.Mat4Slot.__init__(self, mat4(1), 2) # 2 = NO_INPUT_CONNECTIONS
        self._proc = proc

    def computeValue(self):
        self._value = self._proc()

class ProceduralQuatSlot(_core.QuatSlot):
    def __init__(self, proc):
        _core.QuatSlot.__init__(self, quat(), 2) # 2 = NO_INPUT_CONNECTIONS
        self._proc = proc

    def computeValue(self):
        self._value = self._proc()


# NotificationForwarder
class NotificationForwarder(_core.Dependent):
    """Forwards slot notifications.

    This class has the same functionality than the C++ class with the
    same name. However, it's not wrapping the C++ class but is pure
    Python.
    """
    def __init__(self, onvalchanged, onresize=None):
        """Constructor.

        If the forwarder is used with an array slot, the onvalchanged
        callback must take two arguments start and end. If it's used with
        a normal slot, then it does not take any arguments.

        \param onvalchanged A callable object that gets called whenever the value changes
        \param onresize A callable object that gets called whenever the size of an array slot changes
        """
        _core.Dependent.__init__(self)
#        self.callable = callable
        self.onvalchanged = onvalchanged
        self.onresize = onresize
        
    def onValueChanged(self, start=None, end=None):
        # Does the call come from a normal slot?
        if start==None:
            self.onvalchanged()
        # or an array slot?
        else:
            self.onvalchanged(start, end)

    def onResize(self, size):
        if self.onresize!=None:
            self.onresize(size)


# slotPropertyCode
def slotPropertyCode(name, slotname=None):
    """Create the code to add a slot property to a class.

    This is a helper function that creates the code string which
    adds a get() method, a set() method and a property to access
    the value of a slot. The resulting string can then be executed
    inside a class definition via an exec statement. The slot itself
    must be present as a variable in the class.

    The returned string looks like this:
    
    \code
    # Property: <name>
    def _get<Name>(self):
        return self.<slotname>.getValue()

    def _set<Name>(self, val):
        self.<slotname>.setValue(val)

    <name> = property(_get<Name>, _set<Name>, None, "<name>")
    \endcode

    In a %component class (or whatever class) the property is created
    as shown in the following example:

    \code
    class Sphere:
        
        exec slotPropertyCode("radius")

        def __init__(self, ...):
            self.radius_slot = ...
    \endcode

    If no slot name is given then the name is assumed to be '<name>_slot'.
    
    \param name (\c str) Name of the new attribute (property)
    \param slotname (\c str) Name of the corresponding slot
    """
    capname = name.capitalize()
    getmethname = "_get%s"%(capname)
    setmethname = "_set%s"%(capname)
    if slotname==None:
        slotname = "%s_slot"%name

    s="""# Property: %s
def %s(self):
    return self.%s.getValue()

def %s(self, val):
    if isinstance(val, _core.Component):
        val.output_slot.connect(self.%s)
    else:
        self.%s.setValue(val)

%s = property(%s, %s, None, "%s")
"""%(name,
     getmethname,slotname,
     setmethname,slotname,slotname,
     name,getmethname,setmethname,capname)
    return s
