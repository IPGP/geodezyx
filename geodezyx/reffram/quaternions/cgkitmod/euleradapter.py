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
# $Id: euleradapter.py,v 1.1.1.1 2004/12/12 14:30:59 mbaas Exp $

## \file euleradapters.py
## Contains the EulerAdapter class.

from . import protocols
from .Interfaces import *
from .component import *
from . import slots
from .cgtypes import *
from math import pi
from . import _core

# EulerAdapter
class EulerAdapter(Component):
    """Euler angle to mat3, mat4 or quat adapter.

    This class can be used to convert euler angles either to a mat3, a mat4
    or a quat. The input slots are \c anglex_slot, \c angley_slot and
    \c anglez_slot. The output slot is \c output_slot. The type of the output
    can be determined in the constructor.
    """

    protocols.advise(instancesProvide=[ISceneItem])

    def __init__(self,
                 anglex = 0,
                 angley = 0,
                 anglez = 0,
                 radians = False,
                 order = "xyz",
                 outtype = "mat3",
                 name="EulerAdapter",
                 auto_insert=True):
        """Constructor.

        \param anglex (\c float) Initial angle around x axis
        \param angley (\c float) Initial angle around y axis
        \param anglez (\c float) Initial angle around z axis
        \param radians (\c bool) True = Angles are specified in radians instead of degrees
        \param order (\c str) The rotation order ("xyz", "xzy", ...)
        \param outtpye (\c str) Output type ("mat3", "mat4", "quat")
        \param name (\c str) Component name
        \param auto_insert (\c bool) Auto insert flag
        """
        Component.__init__(self, name=name, auto_insert=auto_insert)

        if radians:
            self.factor = 1.0
        else:
            self.factor = pi/180.0

        self.anglex_slot = slots.DoubleSlot(anglex)
        self.angley_slot = slots.DoubleSlot(angley)
        self.anglez_slot = slots.DoubleSlot(anglez)
        self.addSlot("anglex", self.anglex_slot)
        self.addSlot("angley", self.angley_slot)
        self.addSlot("anglez", self.anglez_slot)

        if outtype=="mat3":
            self.output_slot = slots.ProceduralMat3Slot(self.computeMat3)
        elif outtype=="mat4":
            self.output_slot = slots.ProceduralMat4Slot(self.computeMat4)
        elif outtype=="quat":
            self.output_slot = slots.ProceduralQuatSlot(self.computeQuat)
        else:
            raise ValueError("Unknown output type: %s"%outtype)
        
        self.addSlot("output", self.output_slot)
            
        self.anglex_slot.addDependent(self.output_slot)
        self.angley_slot.addDependent(self.output_slot)
        self.anglez_slot.addDependent(self.output_slot)

        # self.fromEuler is the mat3() method that computes the matrix
        # from the euler angles. Which one exactly it is depends on the
        # order
        exec("self.fromEuler = mat3.fromEuler%s"%order.upper())

    def protocols(self):
        return [ISceneItem, IComponent]

    def computeMat3(self):
        """Slot procedure."""
        f = self.factor
        return self.fromEuler(f*self.anglex_slot.getValue(),
                              f*self.angley_slot.getValue(),
                              f*self.anglez_slot.getValue())

    def computeMat4(self):
        """Slot procedure."""
        f = self.factor
        m3 = self.fromEuler(f*self.anglex_slot.getValue(),
                            f*self.angley_slot.getValue(),
                            f*self.anglez_slot.getValue())
        res = mat4(1)
        res.setMat3(m3)
        return res

    def computeQuat(self):
        """Slot procedure."""
        f = self.factor
        m3 = self.fromEuler(f*self.anglex_slot.getValue(),
                            f*self.angley_slot.getValue(),
                            f*self.anglez_slot.getValue())
        res = quat().fromMat(m3)
        return res

    ## protected:
        
    # angle properties...
    exec(slotPropertyCode("anglex"))
    exec(slotPropertyCode("angley"))
    exec(slotPropertyCode("anglez"))

    # "output" property...
    exec(slotPropertyCode("output"))
    

