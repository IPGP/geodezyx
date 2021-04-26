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
# $Id: targetcamera.py,v 1.8 2005/07/21 15:07:11 mbaas Exp $

## \file targetcamera.py
## Contains the TargetCamera class.

"""This module contains the TargetCamera class."""

from .Interfaces import *
from . import protocols
from . import slots
from .cgtypes import *
from math import pi
from .worldobject import WorldObject
from .globalscene import getScene
from . import camerabase, lookat
import sys
from . import _core

# TargetCamera
class TargetCamera(camerabase.CameraBase):
    """A camera object that always looks at a particular target.
    """

    protocols.advise(instancesProvide=[ISceneItem, ICamera])

    def __init__(self,
                 name = "TargetCamera",
                 fov = 45.0,
                 target = vec3(0,0,0),
                 roll = 0.0,
                 up = None,
                 fstop = 0,
                 focallength = 0,
                 **params):
        
        camerabase.CameraBase.__init__(self, name=name, **params)

        target = vec3(target)

        # FOV
        self.fov_slot = slots.DoubleSlot(fov)
        self.addSlot("fov", self.fov_slot)

        # Target
        self.target_slot = slots.Vec3Slot(target)
        self.addSlot("target", self.target_slot)

        # Roll
        self.roll_slot = slots.DoubleSlot(roll)
        self.addSlot("roll", self.roll_slot)

        # Up
        self.up_slot = slots.Vec3Slot()
        self.addSlot("up", self.up_slot)
        if up==None:
            self.up_slot.setValue(getScene().up)
        else:
            self.up_slot.setValue(vec3(up))

        self._lookat = lookat.LookAt()
        self._lookat.name = "TargetCamera_LookAt"
        self._lookat.output_slot.connect(self.rot_slot)
        self.pos_slot.connect(self._lookat.pos_slot)
        self.target_slot.connect(self._lookat.target_slot)
        self.roll_slot.connect(self._lookat.roll_slot)
        self.up_slot.connect(self._lookat.up_slot)

        # fstop
        self.fstop_slot = slots.DoubleSlot(fstop)
        self.addSlot("fstop", self.fstop_slot)
        
        # focal length
        self.focallength_slot = slots.DoubleSlot(focallength)
        self.addSlot("focallength", self.focallength_slot)


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
        del self.target_slot
        del self._lookat.output_slot
        del self._lookat.pos_slot
        del self._lookat.target_slot
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
    exec(slots.slotPropertyCode("roll"))
    exec(slots.slotPropertyCode("up"))

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
        
    # "target" property...
    
    def _getTarget(self):
        """Return the current target position.

        This method is used for retrieving the \a target property.

        \return Target position (\c vec3)
        """
        return self.target_slot.getValue()

    def _setTarget(self, pos):
        """Set a new target position.

        This method is used for setting the \a target property.

        \param pos (\c vec3) Target position
        """
        pos = vec3(pos)
        self.target_slot.setValue(pos)

    target = property(_getTarget, _setTarget, None, "Target position")
