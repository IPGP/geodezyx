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
# $Id: timer.py,v 1.1.1.1 2004/12/12 14:31:27 mbaas Exp $

## \file timer.py
## Contains the Timer class.

import time
from .component import Component
from .eventmanager import eventManager
from . import slots, events

# Timer
class Timer(Component):
    """%Timer class that manages the virtual time.

    Attributes:

    - time (Slot) - Current virtual time
    - timestep - Delta time step
    - duration - Total duration of the animation/simulation
    
    - frame - Current frame number (as float)
    - fps - Frame rate (as float)
    - framecount - Total number of frames (as float)
    
    - clock - Contains the real time (this is a real stop watch)

    The actual slot object of each slot attribute can be accessed by appending
    "_slot" the attribute name.
    """

    def __init__(self, name="Timer", auto_insert=True):
        Component.__init__(self, name=name, auto_insert=auto_insert)

        # The time() value when the clock was started
        self._clock_start  = 0.0
        # This offset is added to the returned clock value. It stores the
        # clock value at the time the clock was stopped
        self._clock_offset = 0.0
        # This flag indicates if the clock is running or not
        self._clock_running = False

        self._duration = 1.0
        self._framecount = 25
        self._fps = 25
        self._timestep = 0.04
        
        # Virtual time...
        self.time_slot = slots.DoubleSlot()
        self.addSlot("time", self.time_slot)

        eventManager().connect(events.RESET, self)
        

    def startClock(self):
        self._clock_start = time.time()
        self._clock_running = True        

    def stopClock(self):
        self._clock_offset = self.clock
        self._clock_running = False

    def step(self):
        t = self.time_slot.getValue() + self._timestep
        self.time_slot.setValue(t)
        eventManager().event(events.STEP_FRAME)

    ## protected:

    def onReset(self):
        """Reset callback."""
        self.time = 0.0
        
    # "time" property...
    
    def _getTime(self):
        """Return the current virtual time.

        This method is used for retrieving the \a time property.

        \return Virtual time (\c float)
        """
        return self.time_slot.getValue()

    def _setTime(self, t):
        """Set the virtual time.

        This method is used for setting the \a time property.

        \param t (\c float) Time value
        """
        self.time_slot.setValue(t)

    time = property(_getTime, _setTime, None, "Virtual time")

    # "frame" property...
    
    def _getFrame(self):
        """Return the current frame number.

        This method is used for retrieving the \a frame property.

        \return Frame number (\c float)
        """
        return self.time_slot.getValue()/self._timestep

    def _setFrame(self, f):
        """Set the frame number.

        This method is used for setting the \a frame property.

        \param f (\c float) Frame number
        """
        self.time_slot.setValue(f*self._timestep)

    frame = property(_getFrame, _setFrame, None, "Frame number")

    # "clock" property...

    def _getClock(self):
        """Return the current clock value (real time).

        This method is used for retrieving the \a clock property.

        \return Real time (\c float)
        """
        if self._clock_running:
            return self._clock_offset + time.time()-self._clock_start
        else:
            return self._clock_offset

    def _setClock(self, t):
        """Set the clock value.

        This method is used for setting the \a clock property.

        \param t (\c float) Time value
        """
        self._clock_offset = float(t)
        self._clock_start = time.time()

    clock = property(_getClock, _setClock, None, "Real time")

    # "timestep" property...
    
    def _getTimeStep(self):
        """Return the current time step.

        This method is used for retrieving the \a timestep property.

        \return Time step (\c float)
        """
        return self._timestep

    def _setTimeStep(self, dt):
        """Set the time step.

        This method is used for setting the \a timestep property.

        \param t (\c float) Time value
        """
        self._timestep = float(dt)

    timestep = property(_getTimeStep, _setTimeStep, None, "Time step")

    # "fps" property...
    
    def _getFPS(self):
        """Return the current frame rate.

        This method is used for retrieving the \a fps property.

        \return Frame rate (\c float)
        """
        return 1.0/self._timestep

    def _setFPS(self, fps):
        """Set the frame rate.

        This method is used for setting the \a fps property.

        \param fps (\c float) New frame rate
        """
        self._timestep = 1.0/float(fps)

    fps = property(_getFPS, _setFPS, None, "Frames per second")

    # "duration" property...
    
    def _getDuration(self):
        """Return the duration of the animation/simulation.

        This method is used for retrieving the \a duration property.

        \return Duration (\c float)
        """
        return self._duration

    def _setDuration(self, d):
        """Set the duration.

        This method is used for setting the \a duration property.

        \param d (\c float) New duration
        """
        self._duration = float(d)

    duration = property(_getDuration, _setDuration, None, "Total duration")

    # "framecount" property...
    
    def _getFrameCount(self):
        """Return the total number of frames.

        This method is used for retrieving the \a framecount property.

        \return Total Number of frames (\c float)
        """
        return self._duration/self._timestep

    def _setFrameCount(self, n):
        """Set the total number of frames.

        This method is used for setting the \a framecount property.

        \param n (\c float) Total number of frames
        """
        self._duration = float(n)*self._timestep

    framecount = property(_getFrameCount, _setFrameCount, None, "Total number of frames")
