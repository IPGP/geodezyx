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
# $Id: events.py,v 1.3 2005/09/20 13:35:14 mbaas Exp $

## \file events.py
## Contains the standard event names and some predefined event classes.

import string
from . import keydefs

# STEP_FRAME is called whenever the timer is stepped forward one frame
# using timer.step().
# The event takes no arguments.
STEP_FRAME = "StepFrame"

# Reset the simulation/animation
RESET = "Reset"

# A key was pressed on the keyboard.
# The argument to the event is a KeyEvent object.
KEY_PRESS = "KeyPress"
# A key was released on the keyboard
# The argument to the event is a KeyEvent object.
KEY_RELEASE = "KeyRelease"

LEFT_DOWN = "LeftDown"
LEFT_UP = "LeftUp"
MIDDLE_DOWN = "MiddleDown"
MIDDLE_UP = "MiddleUp"
RIGHT_DOWN = "RightDown"
RIGHT_UP = "RightUp"
MOUSE_BUTTON_DOWN = "MouseButtonDown"
MOUSE_BUTTON_UP = "MouseButtonUp"
MOUSE_MOVE = "MouseMove"
MOUSE_WHEEL = "MouseWheel"

JOYSTICK_AXIS = "JoystickAxis"
JOYSTICK_BALL = "JoystickBall"
JOYSTICK_HAT = "JoystickHat"
JOYSTICK_BUTTON_DOWN = "JoystickButtonDown"
JOYSTICK_BUTTON_UP = "JoystickButtonUp"

SPACE_MOTION = "SpaceMotion"
SPACE_BUTTON_DOWN = "SpaceButtonDown"
SPACE_BUTTON_UP = "SpaceButtonUp"
SPACE_ZERO = "SpaceZero"

TABLET = "Tablet"

# KeyEvent
class KeyEvent:
    """Keyboard event (key press or release).

    This event is sent as argument to the KEY_PRESS and KEY_RELEASE events.
    The data is stored in the attributes key, keycode and mods.
    """
    
    def __init__(self, key, keycode, mods=0):
        """Constructor.

        \param key (\c unicode) Unicode key
        \param keycode Key (\c int) Key code (untranslated)
        \param mods (\c int) Modifier flags
        """
        # Unicode key
        self.key = key
        # Key code
        self.keycode = keycode
        # Modifier flags
        self.mods = mods

    def __str__(self):
        c = repr(self.key)
#        if self.key in string.printable:
#           c = self.key
#        else:
#            c = "."
        return "<KeyEvent key:%s (%d) mods:%d>"%(c,self.keycode,self.mods)

    # shiftKey
    def shiftKey(self):
        """Return True if the key is a Shift key."""
        return self.keycode==keydefs.KEY_SHIFT_LEFT or \
               self.keycode==keydefs.KEY_SHIFT_RIGHT

    # controlKey
    def controlKey(self):
        """Return True if the key is a Control key."""
        return self.keycode==keydefs.KEY_CONTROL_LEFT or \
               self.keycode==keydefs.KEY_CONTROL_RIGHT

    # altKey
    def altKey(self):
        """Return True if the key is an Alt key."""
        return self.keycode==keydefs.KEY_ALT_LEFT or \
               self.keycode==keydefs.KEY_ALT_RIGHT

# MouseButtonEvent
class MouseButtonEvent:
    """Mouse button event.
    """

    def __init__(self, button, x, y, x0, y0):
        """Constructor.

        \param button (\c int) Button number (1=Left / 2=Middle / 3=Right)
        \param x (\c int) Mouse position (pixel coordinate)
        \param y (\c int) Mouse position (pixel coordinate)
        \param x0 (\c float) Mouse position (normalized)
        \param y0 (\c float) Mouse position (normalized)
        """
        self.button = button
        self.x = x
        self.y = y
        self.x0 = x0
        self.y0 = y0

    def __str__(self):
        return "<MouseButtonEvent button:%d x:%d y:%d x0:%1.2f y0:%1.2f>"%(self.button, self.x, self.y, self.x0, self.y0)

# MouseWheelEvent
class MouseWheelEvent:
    """Mouse wheel event.
    """

    def __init__(self, delta, x, y, x0, y0):
        """Constructor.

        \param delta (\c int) Wheel delta
        \param x (\c int) Mouse position (pixel coordinate)
        \param y (\c int) Mouse position (pixel coordinate)
        \param x0 (\c float) Mouse position (normalized)
        \param y0 (\c float) Mouse position (normalized)
        """
        self.delta = delta
        self.x = x
        self.y = y
        self.x0 = x0
        self.y0 = y0

    def __str__(self):
        return "<MouseWheelEvent delta:%d x:%d y:%d x0:%1.2f y0:%1.2f>"%(self.delta, self.x, self.y, self.x0, self.y0)


# MouseMoveEvent
class MouseMoveEvent:
    """Mouse move event.
    """

    def __init__(self, x, y, dx, dy, x0, y0, dx0, dy0, buttons):
        """Constructor.

        \param x (\c int) Mouse position (pixel coordinate)
        \param y (\c int) Mouse position (pixel coordinate)
        \param dx (\c int) Mouse delta (pixel)
        \param dy (\c int) Mouse delta (pixel)
        \param x0 (\c float) Mouse position (normalized)
        \param y0 (\c float) Mouse position (normalized)
        \param dx0 (\c float) Mouse delta (normalized)
        \param dy0 (\c float) Mouse delta (normalized)
        \param buttons (\c int) Mouse buttons (each bit is a mouse button)
        """
        self.x = x
        self.y = y
        self.dx = dx
        self.dy = dy
        self.x0 = x0
        self.y0 = y0
        self.dx0 = dx0
        self.dy0 = dy0
        self.buttons = buttons

    def __str__(self):
        return "<MouseMoveEvent %d/%d (%d/%d) %1.2f/%1.2f (%1.3f/%1.3f) btns:%s>"%(self.x, self.y, self.dx, self.dy, self.x0, self.y0, self.dx0, self.dy0, hex(self.buttons))
#        return "<MouseMoveEvent x:%d y:%d dx:%d dy:%d x0:%1.2f y0:%1.2f dx:%1.2f dy:%1.2f buttons:%s>"%(self.x, self.y, self.dx, self.dy, self.x0, self.y0, self.dx0, self.dy0, hex(self.buttons))


# JoystickAxisEvent
class JoystickAxisEvent:
    """Joystick axis event.
    """

    def __init__(self, joystick, axis, value):
        """Constructor.

        \param joystick (\c int) Joystick ID (0,1,...)
        \param axis (\c int) Axis ID (0,1,...)
        \param value (\c float) The current value of the specified axis
        """
        self.joystick = joystick
        self.axis = axis
        self.value = value

    def __str__(self):
        return "<JoystickAxisEvent joystick:#%d axis:#%d value:%f>"%(self.joystick, self.axis, self.value)

# JoystickHatEvent
class JoystickHatEvent:
    """Joystick hat event.
    """

    def __init__(self, joystick, hat, x, y):
        """Constructor.

        \param joystick (\c int) Joystick ID (0,1,...)
        \param hat (\c int) Hat ID (0,1,...)
        \param x (\c int) The current x value of the specified hat
        \param y (\c int) The current x value of the specified hat
        """
        self.joystick = joystick
        self.hat = hat
        self.x = x
        self.y = y

    def __str__(self):
        return "<JoystickHatEvent joystick:#%d hat:#%d x:%d y:%d>"%(self.joystick, self.hat, self.x, self.y)

# JoystickBallEvent
class JoystickBallEvent:
    """Joystick ball event.
    """

    def __init__(self, joystick, ball, value):
        """Constructor.

        \param joystick (\c int) Joystick ID (0,1,...)
        \param ball (\c int) Ball ID (0,1,...)
        \param value (\c float) The current value of the specified ball
        """
        self.joystick = joystick
        self.ball = ball
        self.value = value

    def __str__(self):
        return "<JoystickBallEvent joystick:#%d ball:#%d value:%f>"%(self.joystick, self.ball, self.value)

# JoystickButtonEvent
class JoystickButtonEvent:
    """Joystick button event.

    This event is sent as argument to the JOYSTICK_BUTTON_DOWN and
    JOYSTICK_BUTTON_UP events.
    """

    def __init__(self, joystick, button):
        """Constructor.

        \param joystick (\c int) Joystick ID (0,1,...)
        \param button (\c int) Button ID (0,1,...)
        """
        self.joystick = joystick
        self.button = button

    def __str__(self):
        return "<JoystickButtonEvent joystick:#%d button:#%d>"%(self.joystick, self.button)


# SpaceMotion
class SpaceMotionEvent:
    """SpaceMotion event.

    This event is created when a SpaceMouse or SpaceBall is moved or rotated.
    """

    def __init__(self, translation, rotation, period):
        """Constructor.

        \param translation (\c vec3) Translation vector
        \param rotation (\c vec3) Rotation vector
        \param period (\c int) Time in milliseconds since the last event
        """
        self.translation = translation
        self.rotation = rotation
        self.period = period

    def __str__(self):
        return "<SpaceMotion t:%s r:%s period:%d>"%(self.translation, self.rotation, self.period)

# SpaceButton
class SpaceButtonEvent:
    """SpaceButton event.

    This event is generated when a SpaceMouse or SpaceBall button was
    pressed or released.
    """

    def __init__(self, button):
        """Constructor.
        
        \param button (\c int) Button number (1-29)
        """
        self.button = button

    def __str__(self):
        return "<SpaceButton button:%d>"%(self.button)
