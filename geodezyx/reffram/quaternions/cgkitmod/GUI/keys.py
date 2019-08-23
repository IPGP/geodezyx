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

## \file keys.py
## Contains the Keys class.

import wx
import string

# Exceptions:
class InvalidKeyDescription(Exception):
    """Exception class."""
    pass

class KeyNotBound(Exception):
    """Exception class."""
    pass

# Keys
class Keys(object):
    """This class manages key strokes.

    Key bindings correspond to special key attributes. The attribute
    name is the name of the key and the value is a callable object
    that takes no arguments.

    Example:

    \code
    >>> keys = Keys()

    # Bind functions to keys...
    >>> keys.ctrl_c = onExit
    keys.shift_ctrl_tab = onFoo

    # Remove a key binding
    >>> del keys.ctrl_c
    
    \endcode
    """

    def __init__(self):

        # Command dictionary. The keys are the tuples that contain the
        # modifier flags and the key code. The value is a callable object.
        self._commands = {}
        
        # A dictionary that translates wx key codes into a string description
        # The dictionary is built from the WXK_xyz constants
        self._wxkeystrs = {}
        wxkeys = [x for x in dir(wx) if x[0:4]=="WXK_"]
        for k in wxkeys:
            code = getattr(wx, k)
            self._wxkeystrs[code] = k[4:].capitalize()

        # The Window object that uses this key manager
        self._window = None


    def attach(self, window):
        self._window = window
        wx.EVT_KEY_DOWN(window, self._onKeyDown)
#        wx.EVT_KEY_UP(window, self._onKeyUp)
#        wx.EVT_CHAR(window, self._onChar)

    def detach(self):
        pass

    def findCommandKeys(self, func):
        """Search the key table for a command and return all associated keys.

        \param func (\c callable) The bound function
        \return A list of readable key strings.
        """

        res = []
        # Compare all commands...
        for key in self._commands.keys():
            cmd = self._commands[key]
            if func==cmd:
                res.append(self._cmdkey2text(key))

        return res

    ######################################################################
    ## protected:

#    def __setitem__(self, name, value):
#        print "setitem"
#        key = self._str2cmdkey(name)
#        print "Setting",key,"to",value
#        self._commands[key] = value        

    def _onKeyDown(self, event):
        """Handle a KeyDown event.

        Note: This event handler is \b not called if the corresponding
        description text is present in the menu. In that case, the menu
        handles the events itself (a wx feature). But it seems this only
        works with English modifier names.
        """
        print(("KeyDown",))
        print(("KeyCode:",event.GetKeyCode(),"RawKeyCode:",event.GetRawKeyCode()))

        # Create the key tuple
        key = self._event2cmdkey(event)
        # Check if the key was bound
        if key in self._commands:
            # Call the bound function
            func = self._commands[key]
            func()
        else:
            event.Skip()
    
    def _onKeyUp(self, event):
#        print "KeyUp",
#        print "KeyCode:",event.GetKeyCode(),"RawKeyCode:",event.GetRawKeyCode()
        event.Skip()

    def __setattr__(self, name, value):
        if name=="":
            return
        # If the name starts with an underscore then it's not a key description
        if name[0]=="_":
            object.__setattr__(self, name, value)
            return

        # Bind key...

        # Create the key tuple
        key = self._str2cmdkey(name)
        # Remove an existing binding
        if key in self._commands:
            self.__delattr__(name)
        # Store the key binding
        self._commands[key] = value

        # Update menu items
        self._updateMenu(value)

    def __delattr__(self, name):
        # If the name starts with an underscore then it's not a key description
        if name[0]=="_":
            object.__delattr__(self, name)
            return

        # Create the key tuple
        key = self._str2cmdkey(name)
        # Check if a key binding exists
        if key in self._commands:
            func = self._commands[key]
            del self._commands[key]
            # Update menu items
            self._updateMenu(func)
        else:
            raise KeyNotBound("Key '%s' is not bound to a function."%name)


    def _updateMenu(self, func):
        """Update all menu items that are bound to func."""
        # Check if any menu items have to be updated
        menu = getattr(self._window, "menu", None)
        if menu==None:
            return
        nodes = menu.findCommandNodes(func)
        for n in nodes:
            n.update()        
        

    def _event2cmdkey(self, event):
        """Convert a key event into a tuple which can be used as key."""
        return (bool(event.ShiftDown()), bool(event.ControlDown()),
                bool(event.AltDown()), bool(event.MetaDown()),
                event.GetKeyCode())

    def _str2cmdkey(self, s):
        """Convert a key description string into the key tuple."""

        f = s.upper().split("_")

        # Check if the modifiers are correct
        for m in f[:-1]:
            if m not in ["SHIFT", "CTRL", "ALT", "META"]:
                raise InvalidKeyDescription("Key '%s' contains invalid modifiers."%s)

        # Convert the last argument into a key code
        skey = f[-1]
        if len(skey)==1:
            key = ord(skey)
        else:
            key = getattr(wx, "WXK_%s"%skey, None)
            if key==None:
                raise InvalidKeyDescription("Key '%s' does not exist."%s)
        
        return ("SHIFT" in f, "CTRL" in f, "ALT" in f, "META" in f, key)

    def _cmdkey2text(self, key):
        shift, ctrl, alt, meta, key = key
        mods = []
        if shift:
            mods.append("Shift")
#            mods.append("Umschalt")
        if ctrl:
            mods.append("Ctrl")
#            mods.append("Strg")
        if alt:
            mods.append("Alt")
        if meta:
            mods.append("Meta")

        ks = self._wxkeystrs.get(key, None)
        if ks==None:
            ks=chr(key)
        mods.append(ks)
            
        return string.join(mods, "+")
        
            
        
    
        
