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
# $Id: autocam.py,v 1.1.1.1 2004/12/12 14:30:54 mbaas Exp $

## \file autocam.py
## Contains the AutoCam class.

"""This module contains the AutoCam class."""

from .Interfaces import *
from . import protocols
from . import slots
from . import _core
from .eventmanager import eventManager
from .globalscene import getScene
from . import events
from .cgtypes import *
from .component import *
       

# AutoCam
class AutoCam(Component):
    """AutoCam
    """

    protocols.advise(instancesProvide=[ISceneItem])

    def __init__(self, name="AutoCam"):
        Component.__init__(self, name=name)

        # Input
        self.input_slot = slots.Vec3Slot()
        self.addSlot("input", self.input_slot)

        # Output
        self.output_slot = slots.Vec3Slot()
        self.addSlot("output", self.output_slot)

        self.accel_factor = 0.5
        self.out_offset = vec3(0,0,0)

        self.out = vec3(0,0,0)
        self.out_velocity = vec3(0,0,0)

        eventManager().connect(events.STEP_FRAME, self)

    def protocols(self):
        return [ISceneItem, IComponent, IWorldObject, ICamera]

#    def destroy(self):
#        del self.fov_slot
#        del self.target_slot

    def onStepFrame(self):
        dt = getScene().timer().timestep
        diff = self.input-self.out
        if (diff.length()<1.0):
            a = -0.75*self.out_velocity
        else:
            a = self.accel_factor*(self.input-self.out)
        self.out_velocity += dt*a
        self.out = self.out + dt*self.out_velocity
        self.output = self.out+self.out_offset
     
    ## protected:
        
    # "input" property...
    exec(slotPropertyCode("input"))

    # "output" property...
    exec(slotPropertyCode("output"))

        

