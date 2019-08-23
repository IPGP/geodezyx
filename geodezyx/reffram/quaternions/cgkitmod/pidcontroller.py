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
# $Id: pidcontroller.py,v 1.1.1.1 2004/12/12 14:31:11 mbaas Exp $

## \file pidcontroller.py
## Contains the PIDController class.

from . import component
from . import eventmanager, events
from .globalscene import getScene
from .slots import *
from .cgtypes import *
from . import _core

# PIDController
class PIDController(component.Component):
    """PID controller.

    """

    exec(slotPropertyCode("input"))
    exec(slotPropertyCode("output"))
    exec(slotPropertyCode("setpoint"))
    exec(slotPropertyCode("maxout"))
    exec(slotPropertyCode("minout"))
    exec(slotPropertyCode("Kp"))
    exec(slotPropertyCode("Ki"))
    exec(slotPropertyCode("Kd"))

    def __init__(self,
                 name = "PIDController",
                 setpoint = 0.0,
                 Kp = 0.0,
                 Ki = 0.0,
                 Kd = 0.0,
                 maxout = 999999,
                 minout = -999999,
                 auto_insert = True):
        """Constructor.
        """
        
        component.Component.__init__(self, name, auto_insert)

        self.input_slot = DoubleSlot()
        self.setpoint_slot = DoubleSlot(setpoint)
        self.maxout_slot = DoubleSlot(maxout)
        self.minout_slot = DoubleSlot(minout)
        self.Kp_slot = DoubleSlot(Kp)
        self.Ki_slot = DoubleSlot(Ki)
        self.Kd_slot = DoubleSlot(Kd)

        self.output_slot = ProceduralDoubleSlot(self.computeOutput)

        self.addSlot("input", self.input_slot)
        self.addSlot("setpoint", self.setpoint_slot)
        self.addSlot("output", self.output_slot)
        self.addSlot("maxout", self.maxout_slot)
        self.addSlot("minout", self.minout_slot)
        self.addSlot("Kp", self.Kp_slot)
        self.addSlot("Ki", self.Ki_slot)
        self.addSlot("Kd", self.Kd_slot)

        self.input_slot.addDependent(self.output_slot)
        self.setpoint_slot.addDependent(self.output_slot)
        self.maxout_slot.addDependent(self.output_slot)
        self.minout_slot.addDependent(self.output_slot)
        self.Kp_slot.addDependent(self.output_slot)
        self.Ki_slot.addDependent(self.output_slot)
        self.Kd_slot.addDependent(self.output_slot)

        self._integral = 0.0
        self._prev_err = 0.0

        eventmanager.eventManager().connect(events.STEP_FRAME, self)
        eventmanager.eventManager().connect(events.RESET, self)
        

    def onStepFrame(self):
        err = self.setpoint-self.input
        dt = getScene().timer().timestep
        self._integral += dt*err

    def onReset(self):
        self._integral = 0.0
        self._prev_err = 0.0

    def computeOutput(self):
        err = self.setpoint-self.input
        dt = getScene().timer().timestep
        I = self._integral
        D = (err-self._prev_err)/dt
#        print "D:",D
        res = self.Kp*err + self.Ki*I + self.Kd*D
        
        self._prev_err = err
        
        maxout = self.maxout
        minout = self.minout
        if res>maxout:
            res = maxout
        elif res<minout:
            res = minout
            
        return res
