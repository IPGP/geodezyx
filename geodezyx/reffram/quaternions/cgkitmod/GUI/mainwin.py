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
# $Id: mainwin.py,v 1.2 2006/01/27 07:52:40 mbaas Exp $

import wx
from . import window
#from cgkit import getApp
from .menu import *
from .panels import Panels, LayoutNode, PanelNode, PanelWidget, HORIZONTAL, VERTICAL
from .toolbars import ToolBarPalette

# MainWin
class MainWindow(window.Window):
    def __init__(self):
        window.Window.__init__(self)

        # Create menu
        file = Menu("&File", name="file", items=
                    [("&New"),
                     ("&Open..."),
                     ("&Save"),
                     ("Save &as...", None),
                     ("&Close"),
                     MenuItem("---",name="sep1"),
                     ("&Exit", None)
                     ])

        self.menu = Menu(items=[file])


        # ToolBar palette
        self.toolbars = ToolBarPalette(self.wx, self._menu_id_manager)

        # Create panels
        self.panels = Panels(self.wx)
        self.panels.layout = PanelNode(name="Dummy")


        self.mainlayout = wx.BoxSizer(wx.VERTICAL)
        self.mainlayout.AddSizer(self.toolbars.sizer, 0, wx.EXPAND, 0)
        self.mainlayout.Add(self.panels.wx, 1, wx.EXPAND, 0)
        self.wx.SetSizer(self.mainlayout)


#        views = LayoutNode(splittertype=VERTICAL, name="splitter1")

#        self.btn = wx.Button(parent, -1, "Button", wx.Point(0,0))
#        front = PanelNode(name="front", activatable=True,
#                          widget=PanelWidget(wx=self.btn))

#        dict = globals()
#        dict["app"]=getApp()
#        self.shell = wx.py.shell.Shell(parent, -1, locals=dict)
#        shell = PanelNode(name="shell", widget=PanelWidget(wx=self.shell))
        
#        views.setChildren((front, shell))
#        self.panels.layout = views
#        front.makeCurrent()

#        self.panels.updateLayout()

        if "mainwin.geometry" in getApp().prefs:
            x,y,w,h = getApp().prefs["mainwin.geometry"]
            print(("set",x,y,w,h))
            self.SetDimensions(x,y,w,h)
            if getApp().prefs["mainwin.maximized"]:
                self.Maximize(True)
           
        wx.EVT_SIZE(self, self.onResize)
        wx.EVT_MOVE(self, self.onMove)
        return

    def onMove(self, event):
        app = getApp()
        if "mainwin.geometry" in app.prefs:
            x,y,w,h = app.prefs["mainwin.geometry"]
        else:
            w,h = self.GetSizeTuple()
        x,y = self.GetPosition()
        if not self.IsMaximized():
            app.prefs["mainwin.geometry"] = [x,y,w,h]

    def onResize(self, event):
        app = getApp()
        if "mainwin.geometry" in app.prefs:
            x,y,w,h = app.prefs["mainwin.geometry"]
        else:
            x,y = self.GetPositionTuple()
        w,h = self.GetSize()
        if self.IsMaximized():
            app.prefs["mainwin.maximized"] = True
        else:
            app.prefs["mainwin.maximized"] = False
            app.prefs["mainwin.geometry"] = [x,y,w,h]
        event.Skip()

