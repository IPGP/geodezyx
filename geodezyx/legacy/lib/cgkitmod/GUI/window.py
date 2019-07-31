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

import wx
from . import idmanager
from . import keys

class Window(object):
    """Frame class.

    Attributes:

    - wx
    - keys
    - menu
    """


    def __init__(self,
                 parent = None,
                 id = -1,
                 title = 'Title',
                 pos = wx.DefaultPosition,
                 size = wx.DefaultSize,
                 style = wx.DEFAULT_FRAME_STYLE | wx.CLIP_CHILDREN):
        """Create a Frame instance.
        """

        # Create the actual window (a wx Frame)
        self._wx = wx.Frame(parent, id, title, pos, size, style)

        self._keys = keys.Keys()
        self._keys.attach(self)

        # Set up the wx main menu bar
        self._menu_id_manager = idmanager.IDManager(10)
        menubar = wx.MenuBar()
        self._wx.SetMenuBar(menubar)
        self._menu = None

        self.CreateStatusBar()

        wx.EVT_LEFT_DOWN(self, self.onLeftDown)
        wx.EVT_LEFT_UP(self, self.onLeftUp)
        wx.EVT_MIDDLE_DOWN(self, self.onMiddleDown)
        wx.EVT_MIDDLE_UP(self, self.onMiddleUp)
        wx.EVT_RIGHT_DOWN(self, self.onRightDown)
        wx.EVT_RIGHT_UP(self, self.onRightUp)

    def popupMenu(self, menu, pos):
        """Open a menu as a popup menu at the given position.

        \param menu (\c Menu)  Menu tree
        \param pos (\c 2-sequence) Position of the popup menu
        """
        menu._attach(self._menu_id_manager, self, None)
        x,y = pos
        self.PopupMenuXY(menu.wx, x, y)
        menu._detach()

    def onLeftDown(self, event):
        event.Skip()

    def onLeftUp(self, event):
        event.Skip()

    def onMiddleDown(self, event):
        event.Skip()

    def onMiddleUp(self, event):
        event.Skip()

    def onRightDown(self, event):
        event.Skip()

    def onRightUp(self, event):
        event.Skip()


    ######################################################################
    ## protected:

    def __getattr__(self, name):
        """Forward attribute requests to the wxFrame object."""
        
        val = getattr(self._wx, name, None)
        if val==None:
            raise AttributeError("'Window' object has no attribute '%s'."%name)
        else:
            return val

        
    # "wx" property...
    
    def _getWx(self):
        """Return the corresponding wxFrame object.

        This method is used for retrieving the \a wx property.

        \return Frame (\c wxFrame).
        """
        return self._wx

    wx = property(_getWx, None, None, "wxFrame object")


    # "keys" property...
    
    def _getKeys(self):
        """Return the key manager.

        This method is used for retrieving the \a keys property.

        \return Key manager (\c Keys).
        """
        return self._keys

    keys = property(_getKeys, None, None, "Key manager")


    # "menu" property...
    
    def _getMenu(self):
        """Return the menu tree.

        This method is used for retrieving the \a menu property.

        \return Root menu node (\c MenuNode).
        """
        return self._menu

    def _setMenu(self, menu):
        """Set a menu.

        This method is used for setting the \a menu property.
        """

        if self._menu!=None:
            self._menu._detach()

        self._menu = menu
        self._menu._attach(self._menu_id_manager, self, self._wx.GetMenuBar())


    menu = property(_getMenu, _setMenu, None, "Menu tree")


