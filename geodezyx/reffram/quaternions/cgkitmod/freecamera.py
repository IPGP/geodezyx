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
# $Id: freecamera.py,v 1.6 2005/05/12 12:58:29 mbaas Exp $

## \file freecamera.py
## Contains the FreeCamera class.

"""This module contains the FreeCamera class."""

from .Interfaces import *
from . import protocols
from . import slots
from .cgtypes import *
from math import pi
from .worldobject import WorldObject
from .globalscene import getScene
from . import camerabase
import sys
from . import _core
       

# FreeCamera
class FreeCamera(camerabase.CameraBase):
    """A camera object that is free to move and rotate.
    """

    protocols.advise(instancesProvide=[ISceneItem, ICamera])

    def __init__(self,
                 name = "FreeCamera",
                 target = None,
                 fov = 45.0,
                 fstop = 0,
                 focallength = 0,
                 **params):
        
        camerabase.CameraBase.__init__(self, name=name, **params)

        # FOV
        self.fov_slot = slots.DoubleSlot(fov)
        self.addSlot("fov", self.fov_slot)

        # fstop
        self.fstop_slot = slots.DoubleSlot(fstop)
        self.addSlot("fstop", self.fstop_slot)
        
        # focal length
        self.focallength_slot = slots.DoubleSlot(focallength)
        self.addSlot("focallength", self.focallength_slot)

        # Initial targeting
        if target!=None:
            up = getScene().up
            self.transform = mat4().lookAt(self.pos, vec3(target), up)


    def protocols(self):
        return [ISceneItem, IComponent, IWorldObject, ICamera]

    def destroy(self):
#        self.fov_slot.setController(None)
#        self.target_slot.setController(None)
#        self._lookat.pos_slot.setController(None)
#        self._lookat.target_slot.setController(None)
#        self.transform_slot.setController(None)
#        self.pos_slot.setController(None)
#        self.rot_slot.setController(None)
#        self.scale_slot.setController(None)
        del self.fov_slot
        del self.fstop_slot
        del self.focallength_slot

    # projection
    def projection(self, width, height, near, far):
        return mat4().perspective(self.fov, float(width)/height, near, far)

    # viewTransformation
    def viewTransformation(self):
        return self.worldtransform.inverse()


    ## protected:

    exec(slots.slotPropertyCode("fstop"))
    exec(slots.slotPropertyCode("focallength"))

    # "fov" property...
    
    def _getFOV(self):
        """Return the current field of view.

        This method is used for retrieving the \a fov property.

        \return Field of view in angles (\c float)
        """
        return self.fov_slot.getValue()

    def _setFOV(self, fov):
        """Set the field of view.

        This method is used for setting the \a fov property.

        \param fov (\c float) Field of view in angles (0-180)
        """
        fov = float(fov)
        if fov<0:
            fov = 0.0
        if fov>180:
            fov = 180.0
        self.fov_slot.setValue(fov)

    fov = property(_getFOV, _setFOV, None, "Field of view (in angles)")
        
