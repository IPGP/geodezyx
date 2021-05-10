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
# $Id: flockofbirds.py,v 1.2 2005/07/07 08:24:57 mbaas Exp $

## \file flockofbirds.py
## Contains the FlockOfBirds class.

import types
from . import component
from . import eventmanager, events
from .slots import *
from .cgtypes import *
from . import fob

# FlockOfBirds
class FlockOfBirds(component.Component):
    """Receives tracker values from a Flock of Birds motion tracker.

    For each bird the component creates the following four slots:

    - \c pos<n>_slot (\c Vec3Slot) - Position
    - \c angle<n>_slot (\c Vec3Slot) - Euler angles
    - \c matrix<n>_slot (\c Mat3Slot) - Rotation matrix
    - \c quat<n>_slot (\c QuatSlot) - Quaternion

    where <n> is the number of the bird (1-based). Not all of the slots
    are active at the same time. You can select which slot should carry
    the corresponding value via the bird mode (see constructor).

    """

    def __init__(self,
                 name = "FlockOfBirds",
                 com_port = 0,
                 baud_rate = 115200,
                 timeout = 2.0,
                 num_birds = 1,
                 bird_mode = "p",
                 hemisphere = "forward",
                 auto_insert = True):
        """Constructor.

        The bird mode is one of "p" (position), "a" (angles), "m" (matrix),
        "q" (quaternion) or a combination like "pa", "pm" and "pq".
        The argument \a bird_mode can either be a mode string that will
        be used for all birds, or it can be a list of mode strings for
        every bird.

        \param name (\c str) Component name
        \param com_port (\c int) COM port to use for communicating with the flock (0,1,...)
        \param baud_rate (\c int) Baud rate
        \param timeout (\c float) Time out value in seconds for RS232 operations
        \param num_birds (\c int) Number of birds to use (including the ERC)
        \param bird_mode (\c str or \c list) Tracker mode (selects what kind of data the birds will send)
        \param hemisphere (\c str) Hemisphere setting for the transmitter ("forward", "aft", "upper", "lower", "left", "right")
        \param auto_insert (\c bool) True if the component should be inserted into the scene automatically
        """
        
        component.Component.__init__(self, name, auto_insert)

        self.com_port = com_port
        self.baud_rate = baud_rate
        self.timeout = timeout
        self.num_birds = num_birds
        self.hemisphere = hemisphere

        if isinstance(bird_mode, str):
            self.bird_mode = self.num_birds*[bird_mode]
        else:
            self.bird_mode = bird_mode

        self.bird_mode = [x.upper() for x in self.bird_mode]

        if len(self.bird_mode)!=self.num_birds:
            raise ValueError("%d bird modes expected (got %d)"%(self.num_birds, len(self.bird_mode)))

        # Create slots
        self._slots = []
        for i in range(self.num_birds):
            p = Vec3Slot()
            a = Vec3Slot()
            m = Mat3Slot(mat3(1))
            q = QuatSlot(quat(1))
            self._slots.append((p,a,m,q))
            setattr(self, "pos%d_slot"%(i+1), p)
            setattr(self, "angle%d_slot"%(i+1), a)
            setattr(self, "matrix%d_slot"%(i+1), m)
            setattr(self, "quat%d_slot"%(i+1), q)
            self.addSlot("pos%d"%(i+1), p)
            self.addSlot("angle%d"%(i+1), a)
            self.addSlot("matrix%d"%(i+1), m)
            self.addSlot("quat%d"%(i+1), q)

        # Initialize the FOB...
        
        self._fob = fob.FOB()
        self._fob.open(self.com_port, self.baud_rate, self.timeout)
        
        self._fob.fbbAutoConfig(numbirds=self.num_birds)
        self._fob.fbbReset()
        hs = { "forward":fob.FORWARD,
               "aft":fob.AFT,
               "upper":fob.UPPER,
               "lower":fob.LOWER,
               "left":fob.LEFT,
               "right":fob.RIGHT }
        for i in range(1, self.num_birds):
            self._fob.hemisphere(hs[self.hemisphere.lower()], addr=i+1)
            self.setTrackerMode(self.bird_mode[i], i+1)
        self._fob.run()
        
        eventmanager.eventManager().connect(events.STEP_FRAME, self)
        

    ## protected:

    def onStepFrame(self):
        # Poll values
        for i in range(1, self.num_birds):
            values = self._fob.point(i+1)
            if values==None:
                print(("No values received from bird %d."%(i+1)))                
            else:
                # 1. Convert values into proper types
                #    (position, angles, matrix, quat)
                # 2. Set the values on the corresponding slots
#                print values
                ps,angs,ms,qs = self._slots[i]
                m = self.bird_mode[i]
                if m[0]=="P":
                    p = vec3(values[:3])
                    ps.setValue(p)
                c = m[-1]                    
                if c=="A":
                    a = vec3(values[-3:])
                    angs.setValue(a)
                elif c=="M":
                    m = mat3(values[-9:])
                    ms.setValue(m)
                elif c=="Q":
                    q = quat(values[-4:])
                    qs.setValue(q)
                    


    # setTrackerMode
    def setTrackerMode(self, mode, addr):
        """Issue the appropriate fob command to set the tracker mode.

        \param mode (\c str) Tracker mode ("P", "A", "M", "Q", "PA", "PM", "PQ")
        \param addr (\c int) Bird address
        """
        if mode=="P":
            self._fob.position(addr)
        elif mode=="A":
            self._fob.angles(addr)
        elif mode=="M":
            self._fob.matrix(addr)
        elif mode=="Q":
            self._fob.quaternion(addr)
        elif mode=="PA":
            self._fob.position_angles(addr)            
        elif mode=="PM":
            self._fob.position_matrix(addr)
        elif mode=="PQ":
            self._fob.position_quaternion(addr)
        else:
            raise ValueError("Invalid tracker mode: %s"%mode)
        
