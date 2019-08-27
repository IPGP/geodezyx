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
# $Id: panels.py,v 1.2 2006/01/27 07:52:40 mbaas Exp $

import wx
import cgkit.cmds
#from cgkit import getApp, pluginmanager
from cgkit import pluginmanager
import re
from . import panelicons

HORIZONTAL = 0x01
VERTICAL   = 0x02

# Exceptions:

class DuplicateNames(Exception):
    """Exception.

    Raised when there's a name clash between two nodes."""
    pass

class LayoutError(Exception):
    """Exception.

    Raises when a widget is added to a layout that already manages this
    widget."""
    pass

######################################################################

# LayoutRepository
class LayoutRepository(object):
    def __init__(self):
        self._layouts = []


# PanelWidgetRepository
class PanelWidgetRepository(object):
    """Stores all panel widgets that have been created.
    """
    
    def __init__(self):
        """Constructor."""
        self._widgets = []

    def __iter__(self):
        """Return an iterator that iterates over all widgets.
        """
        return iter(self._widgets)

    def insert(self, widget):
        """Inserts a widget into the repository.

        Nothing happens if the widget is already present.
        
        \param widget (\c PanelWidget) Widget that should be inserted into the repository
        """

        if widget not in self._widgets:
            self._widgets.append(widget)

    

# PanelWidget
class PanelWidget(object):
    """Container for a wx widget that serves as a panel widget.

    It's not allowed to have several %PanelWidget object share the same
    wx object.

    Attributes:

    - wx (\c wx \c widget)
    - name (\c str) The name of the widget (this mirrors the name of the wx widget)
    """
    def __init__(self, wx=None):
        """Constructor.

        \param wx (\c wx \c widget) The managed wx widget
        """
        self._wx = wx
        self._active = False
        self._usedby = []

    def __str__(self):
        return '<%s "%s">'%(self.__class__.__name__, self.name)

    # isActive
    def isActive(self):
        """Check if the widget is currently active.

        \return True if active
        """
        return self._active

    # activate
    def activate(self, root):
        """Activate the node.

        This method marks the widget as active.

        Derived classes may use this method to create the wx widget. In
        this case they have to use root.wx as wx parent.        

        \param root (\c Panels) The panels object where the widget is located.
        """
        self._active = True

    # deactivate
    def deactivate(self):
        """Deactivate the node.

        This method marks the widget as inactive.
        """
        self._active = False

    # isInUse
    def isInUse(self):
        """Returns whether the widget is used by any layout or not.

        \return True if the widget is used by a layout (be it active or not).
        """
        return self._usedby!=[]

    # isUsedBy
    def isUsedBy(self, layoutroot):
        """Checks if the widget is used by the given layout.

        \param layoutroot (\c ILayoutNode) The root node of the layout in question
        """
        for n in self._usedby:
            if layoutroot==n._layoutRoot():
                return True
        return False

    # acquire
    def acquire(self, panelnode):
        """Mark the widget as being used by the given node.

        \param panelnode (\c PanelNode) The panel node that controls the position and size of the widget
        """
        self._usedby.append(panelnode)

    # release
    def release(self, panelnode):
        """Mark the widget as being no longer in use by the given node.

        \param panelnode (\c PanelNode) The panel node that was controlling the position and size of the widget
        """        
        self._usedby.remove(panelnode)

    # hide
    def hide(self):
        if self._wx!=None:
            self._wx.Hide()

    # show
    def show(self):
        if self._wx!=None:
            self._wx.Show()

    def setGeom(self, x, y, w, h):
        if self._wx!=None:
            self._wx.SetDimensions(x,y,w,h)
        

    ######################################################################
    ## protected:
    
    # "wx" property...
    
    def _getWx(self):
        """Return the encapsulated wx widget.

        This method is used for retrieving the \a wx property.

        \return wxPython object.
        """
        return self._wx

    def _setWx(self, widget):
        """Set the wx widget.

        This method is used for setting the \a wx property.

        \param widget (\c wx object) wx Widget
        """
        self._wx = widget

    wx = property(_getWx, _setWx, None, "Encapsulated wx widget")

    # "name" property...
    
    def _getName(self):
        """Return the name of the widget.

        This method is used for retrieving the \a name property.

        \return Name (\c str)
        """
        if self._wx==None:
            return ""
        else:
            return str(self._wx.GetName())

    def _setName(self, name):
        """Set the name of the widget.

        This method is used for setting the \a name property.

        \param name (\c str) New name
        """
        if self._wx!=None:
            self._wx.SetName(name)

    name = property(_getName, _setName, None, "Widget name")

######################################################################
    
# Panels
class Panels(object):
    """%Panels area.

    Attributes:

    - \a layout (\c ILayoutNode): Panel layout tree
    - \a wx (\c wx widget): Corresponding wx widget
    - \a repository (\c PanelWidgetRepository) Widget repository

    Panel nodes can be accessed in one of the following ways:

    - The attribute \a mousepanel contains the panel that's currently
      underneath the mouse pointer (or None if the mouse is outside)
    - The attribute \a activepanel contains the currently active panel
      (i.e. the one that's activatable and that received the last mouse
      click)
    - Via attribute access you can address the panels by their name
    - You can iterate over all panel nodes
    

    """
        
    def __init__(self, parent, id=-1, pos=wx.DefaultPosition, size=wx.DefaultSize, style=wx.CLIP_CHILDREN):
        """Constructor.
        """
#        wx.Window.__init__(self, parent, id, pos, size, style)

        self._wx = wx.Window(parent, id, pos, size, style)

        self._repository = PanelWidgetRepository()

        self._layout = LayoutNode(root=self)

        self.mousepanel = None
        self.activepanel = None

        self._xclick = None
        self._yclick = None
        self._xsplitter = None
        self._ysplitter = None
        self._splitter_flags = None
        self._layout_node = None

        self._current_cursor = None
        self._hcrsr = wx.StockCursor(wx.CURSOR_SIZEWE)
        self._vcrsr = wx.StockCursor(wx.CURSOR_SIZENS)
        self._hvcrsr = wx.StockCursor(wx.CURSOR_SIZING)

        wx.EVT_SIZE(self._wx, self.onSize)
        wx.EVT_LEFT_DOWN(self._wx, self.onLeftDown)
        wx.EVT_LEFT_UP(self._wx, self.onLeftUp)
        wx.EVT_MOTION(self._wx, self.onMouseMove)
        wx.EVT_PAINT(self._wx, self.onPaint)

#    def __repr__(self):
#        res = "<Panellayout>"
#        return res

    def __str__(self):
        return str(self._layout)

    def updateLayout(self):
        w,h = self._wx.GetClientSizeTuple()
        self._layout.setGeom(0,0,w,h)
        self._wx.Refresh(False)

    ######################################################################
    ## protected:
        
    def __iter__(self):
        if self._layout!=None:
            return iter(self._layout.panelNodes())
        else:
            return iter([])

    def __getitem__(self, key):
        if self._layout!=None:
            return self._layout.findNodeByName(key)
        else:
            return None

    def __getattr__(self, name):
        if self._layout!=None:
            return self._layout.findNodeByName(name)
        else:
            return None
        
    def onSize(self, event):
        w,h = event.GetSize()
        self._layout.setGeom(0,0,w,h)

    def onEnterPanel(self, panel):
        self.mousepanel = panel
        self._wx.SetCursor(wx.NullCursor)
#        print "enter", panel.name

    def onLeavePanel(self, panel):
        pass
#        print "leave", panel.name

    def onClickPanel(self, panel):
        self.activepanel = panel
        self._wx.Refresh(False)

    def onLeftDown(self, event):
        x = event.GetX()
        y = event.GetY()
        flags, node = self._layout.findPanel(x,y)
        if flags!=0:
            self._splitter_flags = flags
            self._layout_node = node
            self._xclick = x
            self._yclick = y
            self._xsplitter, self._ysplitter = node.getPixelPos()
            self._wx.CaptureMouse()
            
    def onLeftUp(self, event):
        if self._xclick!=None:
            self._xclick = None
            self._yclick = None
            self._xsplitter = None
            self._ysplitter = None
            self._splitter_flags = None
            self._layout_node = None
            self._wx.ReleaseMouse()
            self._wx.SetCursor(wx.NullCursor)
#            self.Refresh(False)
        
    def onMouseMove(self, event):
        x = event.GetX()
        y = event.GetY()
#        print x,y
        if self._xclick!=None:
            dx = x-self._xclick
            dy = y-self._yclick
            if self._splitter_flags & HORIZONTAL==0:
                dx = 0
            if self._splitter_flags & VERTICAL==0:
                dy = 0
            self._layout_node.setPixelPos(self._xsplitter+dx, self._ysplitter+dy, True)
            self._wx.Refresh(False)
        else:
            flags, node = self._layout.findPanel(x,y)
            if flags==HORIZONTAL:
                self._wx.SetCursor(self._hcrsr)
            elif flags==VERTICAL:
                self._wx.SetCursor(self._vcrsr)
            elif flags==HORIZONTAL|VERTICAL:
                self._wx.SetCursor(self._hvcrsr)
            else:
                self._wx.SetCursor(wx.NullCursor)
            
    def onPaint(self, event):
        dc = wx.PaintDC(self._wx)

        dc.BeginDrawing()

        # Paint empty panels black
        dc.SetPen(wx.TRANSPARENT_PEN)
        dc.SetBrush(wx.BLACK_BRUSH)
        for n in self:
            if n.widget==None:
                dc.DrawRectangle(n._x0, n._y0, n._width, n._height)

        # Draw a border around the current panel
        activepen = wx.Pen(wx.Color(0,0,255), 3, wx.SOLID)
        dc.SetPen(activepen)
        panel = self.activepanel
        if panel!=None:
            dc.DrawRectangle(panel._x0, panel._y0, panel._width, panel._height)

        dc.EndDrawing()

    # "layout" property...
    
    def _getLayout(self):
        """Return the layout tree.

        This method is used for retrieving the \a layout property.

        \return Layout tree (\c ILayoutNode).
        """
        return self._layout

    def _setLayout(self, layout):
        """Set the layout tree.

        This method is used for setting the \a layout property.

        \param layout (\c ILayoutNode) Layout tree.
        """
        if self._layout!=None:
            self._layout.deactivate()
        self._layout = layout
        self.updateLayout()
        self._layout.activate(self)

    layout = property(_getLayout, _setLayout, None, "Panel layout tree")

    # "wx" property...
    
    def _getWx(self):
        """Return the encapsulated wx widget.

        This method is used for retrieving the \a wx property.

        \return wxPython object.
        """
        return self._wx

    wx = property(_getWx, None, None, "Encapsulated wx widget")

    # "repository" property...
    
    def _getRepository(self):
        """Return the widget repository.

        This method is used for retrieving the \a repository property.

        \return Widget repository (\c PanelWidgetRepository)
        """
        return self._repository

    repository = property(_getRepository, None, None, "Panel widget repository")
        


######################################################################
######################################################################

# ILayoutNode
class ILayoutNode(object):
    """This is the base class for a node in the layout tree.

    This class maintains the two attributes \em _parent (parent node)
    and \em _root (corresponding Panels object). The root attribute of
    every node in a tree must always point to the same %Panels object.
    If the layout is not active, then the %Panels object is None.

    The root object is set when the layout gets activated (i.e attached
    to a %Panels object). If a new node is created inside an active layout,
    the root has to be provided in the constructor.

    Attributes:

    - \b name (\c str): Node name
    """
    
    def __init__(self, root=None, name=None):
        """Constructor.

        \param root (\c Panels) Panels object which belongs to this layout or None.
        """
        self._parent = None
        self._root = root

        # Node name
        self._name = name
        

    # activate
    def activate(self, root):
        """Activate the node.

        Attaches the layout to the Panels object \a root. If the node has
        children this method has to call their %activate() method as well.
        If there are already PanelWidgets attached to this layout this
        method has to make them visible and do the layout.

        This method is called by the %Panels object whenever the layout
        is switched.

        \pre \a root must not have another layout activated.
        \pre self must be inactive
        \param root (\c Panels) The panels object to which the layout is attached.
        \see deactivate()
        """
        pass

    # deactivate
    def deactivate(self):
        """Deactivates this node.

        Detaches a layout from its root Panels object. All the widgets
        have to be hidden during deactivation.

        This method is called by the %Panels object whenever the layout
        is switched.

        \pre The layout has to be previously activated.
        \see activate()
        """
        pass

    # isActive
    def isActive(self):
        """Checks if this layout is active or not.

        \return True if the layout is active
        """
        return self._root!=None

    # setRoot
    def setRoot(self, root):
        """Set a new root.

        If the node contains children this method has to be overwritten
        and it has to set the new root on the children as well.

        This method is for internal use to propagate a new root through
        the entire tree when the layout is activated or deactivated.

        \param root (\c Panels) New root object.
        """
        self._root = root

    # setGeom
    def setGeom(self, x, y, width, height):
        """Set the position and size of the managed area.

        \param x (\c int) X position in pixel
        \param y (\c int) Y position in pixel
        \param width (\c int) Width in pixel
        \param height (\c int) Height in pixel
        """
        pass

    # applyConstraint
    def applyConstraint(self, width, height, interactive=False):
        """Check if any active constraints accept the given width and height.

        This method checks if \a width and \a height are acceptable for
        this node. If they are, the method has to return the width and
        height unmodified, otherwise an adjusted value has to be returned
        that's as close as possible to the input values.

        The parameter \a interactive specifies if the resizing was
        directly initiated by the user (i.e. a splitter was dragged).
        The usual policy is to allow resizing in this case. Otherwise
        the size is kept fixed (for example, if you only resize the entire
        application window).

        \param width (\c int)  Width to check
        \param height (\c int)  Height to check
        \param interactive (\c bool)  Flag that specifies if the resizing
               stems from a user interaction (the user drags a splitter)
        \return Tuple (width, height) with adjusted size values
        """
        pass

    # isResizable
    def isResizable(self, direction):
        """Check if the size of the node may be changed by the user.

        \param direction (\c int) A combination of HORIZONTAL and VERTICAL
        \return True if resizable.
        """
        pass

    # panelNodes
    def panelNodes(self):
        """Return a list of all panel nodes in this subtree.

        This method only has to return the panel nodes, i.e. the leafes
        of the tree. LayoutNode objects are skipped.

        This method enables the Panels object to iterate over all
        panel nodes.

        \return A list of PanelNode objects.
        """
        pass

    # findNodeByName
    def findNodeByName(self, name):
        """Returns the node in this subtree with the given name.

        \param name (\c str) Node name
        \return Node or None.
        """
        pass

    # makeNameUnique
    def makeNameUnique(self, name):
        """Makes a given name unique by appending an appropriate number.

        The argument \a name is checked if it's already used somewhere
        in the tree. If it's unused it is returned unchanged.
        Otherwise a number is appended that makes the name unique (if
        the original name aready contained a number, this number is
        increased).

        \param name (\c str) Input name
        """
        while self.findNodeByName(name)!=None:
            m = re.search("[0-9]+$", name)
            if m==None:
                name = name+"1"
            else:
                n = int(name[m.start():])
                name = name[:m.start()]+str(n+1)
        return name
            

    # findPanel
    def findPanel(self, x, y):
        return (0, self)

    # showConfigPanel
    def showConfigPanel(self):
        """Show the configuration panel.

        A LayoutNode passes this method call forward to its children
        which finally create the panel.
        """
        pass

    # hideConfigPanel
    def hideConfigPanel(self):
        """Hide the configuration panel.

        A LayoutNode passes this method call forward to its children
        which finally remove the panel.
        """
        pass

    ######################################################################
    ## protected:

    def _layoutRoot(self):
        """Return the root layout node.

        \param Root layout node (\c ILayoutNode).
        """
        p = self
        while p._parent!=None:
            p = p._parent
        return p

    # "name" property...
    
    def _getName(self):
        """Return the node name.

        This method is used for retrieving the \a name property.

        \return Name (\c str)
        """
        return self._name

    def _setName(self, name):
        """Set the node name.

        This method is used for setting the \a name property.

        \param name (\c str) New name
        """
        # Check if the new name exists already....
        prevname = self._name
        self._name = None
        if self._layoutRoot().findNodeByName(name)!=None:
            self._name = prevname
            raise DuplicateNames('There is already a node with the name "%s"'%name)

        self._name = name
        
    name = property(_getName, _setName, None, "Name")


# LayoutNodeCoords
class LayoutNodeCoords:
    """Helper class for the LayoutNode class.
    """
    def __init__(self,x0=0,x1=0,x2=0,x3=0,y0=0,y1=0,y2=0,y3=0):
        self.x0 = x0
        self.y0 = y0
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2
        self.x3 = x3
        self.y3 = y3


# LayoutNode
class LayoutNode(ILayoutNode):
    """Layout node.

    This class represents a layout scheme that can manage areas
    that are laid out as depicted in the following image:
    
    \image html "panel_layout_node.png" "Layout"

    This class is also used if the region should only be split in
    horizontal or vertical direction. If the region is split horizontally
    then x1 = x2 = x3, if it's split vertically then y1 = y2 = y3.
    The following constraints are always enforced on the coordinates:
    x0 <= x1 <= x2 <= x3 and y0 <= y1 <= y2 <= y3. If a vertical splitter
    is present then x1 < x2, if a horizontal splitter is present then
    y1 < y2. All coordinates of all layout nodes are always relative to
    the widget that's associated with the layout hierarchy.

    There are two coordinates associated with each splitter. One is
    the \em logical coordinate that lies within [0,1] (0=left/top,
    1=right/bottom) and the other is the true \em pixel \em coordinate
    relative to the parent widget. In the presence of size constraints
    among the children panels/layouts the logical coordinates remain
    unchanged whereas the true pixel position changes. The pixel coordinates
    always try to represent the logical coordinates as close as possible.
    """
    def  __init__(self, x=0, y=0, width=0, height=0, root=None, name=None,
                  splittertype=HORIZONTAL|VERTICAL, x0=0.5, y0=0.5,
                  children=None):
        """Constructor.

        The arguments should be given as keyword arguments.

        \param x (\c int) Initial X position
        \param y (\c int) Initial Y position
        \param width (\c int) Initial width
        \param height (\c int) Initial height
        \param root (\c Panels) Root object
        \param name (\c str) Node name
        \param splittertype (\c int) A combination of HORIZONTAL and VERTICAL
        \param x0 (\c float) Initial logical splitter x position
        \param y0 (\c float) Initial logical splitter y position
        \param children  A tuple with children nodes (either LayoutNode or PanelNode). The number of children is determined by the splitter type.
        """
        
        ILayoutNode.__init__(self, root=root, name=name)

        # This is the region that's managed by this node
        # This coordinates are relative to the root widget
        # x1/y1 is the true pixel coordinate of the splitter 
        self._x0 = 0
        self._y0 = 0
        self._x1 = 0
        self._y1 = 0
        self._x2 = 0
        self._y2 = 0
        self._x3 = 0
        self._y3 = 0

        # Logical splitter coordinate (0-1)
        self._splitter_x = x0
        self._splitter_y = y0

        # Splitter width in pixel
        self._splitter_width = 6

        # Splitter type
        self._splitter_type = splittertype

        # Children nodes (may not be None)
        # Upper left
        self._child00 = None
        # Upper right
        self._child01 = None
        # Lower left
        self._child10 = None
        # Lower right
        self._child11 = None
        
        if children==None:
            if splittertype==HORIZONTAL|VERTICAL:
                children = (None,None,None,None)
            else:
                children = (None, None)
        self.setChildren(children)

        self.setGeom(x, y, width, height)

    def __str__(self):
        return self._toStr(self, "", 0)

    def _toStr(self, node, res, d):
        """Helper method for the __str__ operator."""

        if isinstance(node, LayoutNode):
            s="?"
            if node._splitter_type==HORIZONTAL:
                s = "horizontal"
            elif node._splitter_type==VERTICAL:
                s = "vertical"
            elif node._splitter_type==HORIZONTAL | VERTICAL:
                s = "horizontal & vertical"
            if node.name!=None:
                name = '"%s"'%node.name
            else:
                name = "<None>"
            s = 'LayoutNode %s active:%d type:%s'%(name, node.isActive(), s)
            res += d*" "+s+"\n"
            childs = node.getChildren()
            for c in childs:
                res = self._toStr(c, res, d+2)
            return res
        elif isinstance(node, PanelNode):
            return res + d*" "+str(node)+"\n"
        else:
            return res + d*" "+"???\n"

    def activate(self, root):
        self.setRoot(root)
        self._child00.activate(root)
        self._child01.activate(root)
        self._child10.activate(root)
        self._child11.activate(root)

    def deactivate(self):
        self.setRoot(None)
        self._child00.deactivate()
        self._child01.deactivate()
        self._child10.deactivate()
        self._child11.deactivate()

    def setRoot(self, root):
        ILayoutNode.setRoot(self, root)
        self._child00.setRoot(root)
        self._child01.setRoot(root)
        self._child10.setRoot(root)
        self._child11.setRoot(root)

    def getSplitterType(self):
        """Return the type of layout node.

        Returns the direction that are split by a splitter. A horizontal
        splitter splits in X direction, a vertical splitter in Y direction.

        \return A combination of HORIZONTAL and VERTICAL
        \see setSplitterType()
        """
        return self._splitter_type

    def setSplitterType(self, stype):
        """Set the layout node type.

        \param stype (\c int) A combination of HORIZONTAL and VERTICAL
        \see getSplitterType()
        """
        self._splitter_type = stype
        

    # setGeom
    def setGeom(self, x, y, width, height):
        # Initialize coordinates...
        self._x0 = int(x)
        self._y0 = int(y)
        self._x3 = int(x+width)
        self._y3 = int(y+height)
        # The following call updates the pixel position of the splitter
        self.setLogicalPos(self._splitter_x, self._splitter_y)

    def findNodeByName(self, name):
        if self._name==name:
            return self
        
        res = self._child00.findNodeByName(name)
        if res==None:
            res = self._child01.findNodeByName(name)
        if res==None:
            res = self._child10.findNodeByName(name)
        if res==None:
            res = self._child11.findNodeByName(name)
        return res

    def showConfigPanel(self):
        self._child00.showConfigPanel()
        if self._splitter_type==HORIZONTAL:
            self._child01.showConfigPanel()
        elif self._splitter_type==VERTICAL:
            self._child10.showConfigPanel()
        else:
            self._child01.showConfigPanel()
            self._child10.showConfigPanel()
            self._child11.showConfigPanel()

    def hideConfigPanel(self):
        self._child00.hideConfigPanel()
        if self._splitter_type==HORIZONTAL:
            self._child01.hideConfigPanel()
        elif self._splitter_type==VERTICAL:
            self._child10.hideConfigPanel()
        else:
            self._child01.hideConfigPanel()
            self._child10.hideConfigPanel()
            self._child11.hideConfigPanel()    
        
    def getChildren(self):
        """Return the children nodes.

        The number of children returned depends on the splitter type
        (either 2 or 4).

        \return A tuple with two or four children nodes (\c ILayoutNode).
        """
        if self._splitter_type==HORIZONTAL:
            return self._child00, self._child01
        elif self._splitter_type==VERTICAL:
            return self._child00, self._child10
        else:
            return self._child00, self._child10, self._child01, self._child11


    def setChildren(self, childs):
        """Set the children of this node.

        \a childs is a tuple of nodes which should be set as children.
        A node may also be None in which case an empty PanelNode is set.
        The children nodes must not be part of another tree (if they are
        an exception is thrown).
        
        \param childs A sequence containing the children. The number of children depends on the splitter type (2 or 4).
        """

        if self._splitter_type==HORIZONTAL:
            self.setChild(childs[0], 0)
            self.setChild(childs[1], 1)
            self.setChild(None, 2)
            self.setChild(None, 3)
        elif self._splitter_type==VERTICAL:
            self.setChild(childs[0], 0)
            self.setChild(None, 1)
            self.setChild(childs[1], 2)
            self.setChild(None, 3)
        else:
            self.setChild(childs[0], 0)
            self.setChild(childs[1], 1)
            self.setChild(childs[2], 2)
            self.setChild(childs[3], 3)


    def childIndex(self, child):
        """Return the index of the children node.

        \a child must actually be a child of self, otherwise a ValueError
        exception is thrown.

        \param child (\c ILayoutNode) Children node
        \return Index (0-3)
        """
        
        lst = [self._child00, self._child01, self._child10, self._child11]
        return lst.index(child)

    def getChild(self, idx):
        """Return a single children node.

        \param idx (\c int) Children index (0-3)
        \return Children node (\c ILayoutNode).
        """
        if idx<0 or idx>3:
            raise IndexError("Children index out of range (%d)"%idx)

        return [self._child00, self._child01, self._child10, self._child11][idx]

    def setChild(self, child, idx):
        """Replace a single children node.

        \param child (\c ILayoutNode) New children or None
        \param idx (\c int) Children index (0-3)
        \todo Pruefen, ob Child tree ok ist (keine doppelten Widgets oder Namen)
        """
        if idx<0 or idx>3:
            raise IndexError("Children index out of range (%d)"%idx)
        if child==None:
#            if self._root!=None:
#                black = wx.Window(self._root.wx, -1)
#                black.SetBackgroundColour(wx.Colour())
#            else:
#                black = None
#            child = PanelNode(widget=black)
            child = PanelNode()
        if child._parent!=None:
            raise ValueError('The panel layout node is already part of a layout tree.')
            
        if idx==0:
            self.removeChild(self._child00)
            self._child00 = child
        elif idx==1:
            self.removeChild(self._child01)
            self._child01 = child
        elif idx==2:
            self.removeChild(self._child10)
            self._child10 = child
        elif idx==3:
            self.removeChild(self._child11)
            self._child11 = child

        child._parent = self
        if self.isActive():
            child.activate(self._root)


    def removeChild(self, child):
        """Remove a child node.

        \a child may be None in which case the method returns immediately.

        \param child (\c ILayoutNode) Children which is to be removed or None
        """
        if child==None:
            return
        
        if child==self._child00:
            self._child00 = None
        elif child==self._child01:
            self._child01 = None
        elif child==self._child10:
            self._child10 = None
        elif child==self._child11:
            self._child11 = None
        else:
            raise ValueError("Layout node is not a children node.")

        child._parent = None
        if self.isActive():
            child.deactivate()

        self._fillEmptyChildren()

    def remove(self, idx=0):
        """Remove this layout node and replace it with a children node.
        """

        child = self.getChild(idx)
        self.setChildren((None,None,None,None))

        root = self._root

        parent = self._parent
        if parent!=None:
            i = parent.childIndex(self)
            parent.removeChild(self)
            parent.setChild(child, i)
        else:
            if root!=None:
                root.layout = child

        if root!=None:
            root.updateLayout()
        
        
    # isResizable
    def isResizable(self, direction):
        res = False
        if direction&HORIZONTAL:
            res |= ( (self._child00.isResizable(HORIZONTAL) and
                      self._child10.isResizable(HORIZONTAL) )
                     or
                     (self._child01.isResizable(HORIZONTAL) and
                      self._child11.isResizable(HORIZONTAL)) )
        if direction&VERTICAL:
            res |= ( (self._child00.isResizable(VERTICAL) and
                      self._child01.isResizable(VERTICAL) )
                     or
                     (self._child10.isResizable(VERTICAL) and
                      self._child11.isResizable(VERTICAL)) )

        return res

    # applyConstraint
    def applyConstraint(self, width, height, interactive=False):
        coords = LayoutNodeCoords(self._x0, 0, 0, self._x0+width,
                                  self._y0, 0, 0, self._y0+height)
        self._calcSplitterCoords(self._splitter_x, self._splitter_y,
                                 coords, interactive)
        return (coords.x3-coords.x0, coords.y3-coords.y0)

    # panelNodes
    def panelNodes(self):
        res = []
        for c in self.getChildren():
            res += c.panelNodes()
        return res

    # isInside
    def isInside(self, x, y):
        """Check if the mouse coordinate (x,y) is inside the managed region.

        The coordinate must be given relative to the root widget.

        \param x (\c int) X coordinate
        \param y (\c int) Y coordinate
        """
        return (x>=self._x0 and x<self._x3 and y>=self._y0 and y<self._y3)

    # findPanel
    def findPanel(self, x, y):
        """Find a panel by mouse position.

        The method assumes that x,y lies somehwere in the managed region.

        The return value contains a flag that has the following values:

        - 0: A panel was hit
        - 0x1: A horizontal splitter was hit
        - 0x2: A vertical splitter was hit
        - 0x3: The middle of a splitter was hit

        \return Tuple (Flag, LayoutNode)
        """
        # Left column?
        if x<self._x1:
            # Upper left panel?
            if y<self._y1:
                return self._child00.findPanel(x,y)
            # Vertical splitter?
            elif y<self._y2:
                return (VERTICAL, self)
            # Lower left panel
            else:
                return self._child10.findPanel(x,y)
        # Horizontal splitter? (or middle)
        elif x<self._x2:
            if y>=self._y1 and y<self._y2:
                return (VERTICAL | HORIZONTAL, self)
            else:
                return (HORIZONTAL, self)
        # Right column
        else:
            # Upper right panel?
            if y<self._y1:
                return self._child01.findPanel(x,y)
            # Vertical splitter?
            elif y<self._y2:
                return (VERTICAL, self)
            # Lower right panel
            else:
                return self._child11.findPanel(x,y)

    # setLogicalPos
    def setLogicalPos(self, x0=0.5, y0=0.5, interactive=False):
        """Set the logical position of the splitters.

        The pixel position is updated as well.
        
        \param x0 (\c float) Logical X position 
        \param y0 (\c float) Logical Y position 
        """

        self._splitter_x = x0
        self._splitter_y = y0

        coords = LayoutNodeCoords(self._x0, 0, 0, self._x3,
                                  self._y0, 0, 0, self._y3)
        self._calcSplitterCoords(x0, y0, coords, interactive)
        self._x1 = coords.x1
        self._x2 = coords.x2
        self._y1 = coords.y1
        self._y2 = coords.y2
        self._update()

    # getLogicalPos
    def getLogicalPos(self):
        """Return the logical position of the splitter.

        \return Tuple (x0, y0) containing the logical position.
        """
        return self._splitter_x, self._splitter_y

    # setPixelPos
    def setPixelPos(self, x, y, interactive=False):
        """Set the pixel position of the splitter.

        This method doesn't set the pixel position directly as the
        position might break some size constraints. Instead the
        position is converted to logical coordinates and setLogicalPos()
        is called instead (so the logical position is updated as well).

        \param x (\c int) X coordinate of the destination pixel position
        \param y (\c int) Y coordinate of the destination pixel position
        """
        x0, y0 = self._pixel2logical(x, y)
        self.setLogicalPos(x0, y0, interactive)

    # getPixelPos
    def getPixelPos(self):
        """Return the pixel position of the splitter.

        \param Tuple (x, y) containing the pixel position.
        """
        return self._x1, self._y1
        

    ######################################################################
    ## protected:

    def _logical2pixel(self, x0, y0):
        """Convert logical coordinates to pixel coordinates.

        \param x0 (\c float) Logical X coordinate
        \param y0 (\c float) Logical Y coordinate
        \return Tuple (x,y) with pixel coordinates (\c int, \c int)
        """
        xmax = self._x3 - self._splitter_width
        ymax = self._y3 - self._splitter_width
        return (int((1.0-x0)*self._x0 + x0*xmax),
                int((1.0-y0)*self._y0 + y0*ymax))

    def _pixel2logical(self, x, y):
        """Convert pixel coordinates to logical coordinates.

        \param x (\c int) X coordinate
        \param y (\c int) Y coordinate
        \return Tuple (x0, y0) with logical coordinates (\c float, \c float)
        """
        w = self._x3 - self._x0 - self._splitter_width
        h = self._y3 - self._y0 - self._splitter_width
        
        if w==0:
            x0 = 0.0
        else:
            x0 = float(x-self._x0)/w
            
        if h==0:
            y0 = 0.0
        else:
            y0 = float(y-self._y0)/h
            
        return x0,y0

    def _update(self):
        """Update the children areas.

        \pre x0-x3/y0-y3 contain valid values
        """
        x0,y0 = self._x0, self._y0
        x1,y1 = self._x1, self._y1
        x2,y2 = self._x2, self._y2
        x3,y3 = self._x3, self._y3
        
        self._child00.setGeom(x0, y0, x1-x0, y1-y0)
        self._child01.setGeom(x2, y0, x3-x2, y1-y0)
        self._child10.setGeom(x0, y2, x1-x0, y3-y2)
        self._child11.setGeom(x2, y2, x3-x2, y3-y2)        

    def _calcSplitterCoords(self, x0, y0, coords, interactive):
        """Calculate new splitter positions.
      
        \param x0 (\c float) Logical X coordinate of the splitter
        \param y0 (\c float) Logical Y coordinate of the splitter
        \param coords (\c LayoutNodeCoords) Coordinates container
        \param interactive (\c bool) Interactive flag
        \pre x0,x3,y0,y3 are set in coords and x0 <= x3 and y0 <= y3
        """
        
        x,y = self._logical2pixel(x0,y0)
        coords.x1 = x
        coords.x2 = x + self._splitter_width
        coords.y1 = y
        coords.y2 = y + self._splitter_width
        self._enforceSplitterConstraints(coords)
        self._enforceSizeConstraints(coords, interactive)
        return coords
        

    def _enforceSizeConstraints(self, coords, interactive):
        """Change the splitter pixel coordinates to enforce size constraints.

        This method changes pixel coordinates of the splitter so that
        the size constraints on the children will be met. If there are
        no constraints then x0-x3/y0-y3 won't be modified.
        The logical coordinates remain unchanged in any case.
        """

        # Is the left column resizable?
        if (self._child00.isResizable(HORIZONTAL) and self._child10.isResizable(HORIZONTAL)):
            # Right column has priority:
            # Determine suggested width of the right column by asking
            # child01 and child11 in both orders and taking the maximum
            w = coords.x3-coords.x2
            w1,h = self._child01.applyConstraint(w,1, interactive)
            w1,h = self._child11.applyConstraint(w1,1, interactive)
            w2,h = self._child11.applyConstraint(w,1, interactive)
            w2,h = self._child01.applyConstraint(w2,1, interactive)
            w = max(w1,w2)
            # Set the final splitter x position 
            coords.x2 = coords.x3 - w
            coords.x1 = coords.x2 - self._splitter_width
            self._enforceSplitterConstraints(coords)
        else:
            # Left column has priority:
            # Determine suggested width of the left column by asking
            # child00 and child10 in both orders and taking the maximum
            w = coords.x1-coords.x0
            w1,h = self._child00.applyConstraint(w,1, interactive)
            w1,h = self._child10.applyConstraint(w1,1, interactive)
            w2,h = self._child10.applyConstraint(w,1, interactive)
            w2,h = self._child00.applyConstraint(w2,1, interactive)
            w = max(w1,w2)
            # Set the final splitter x position 
            coords.x1 = coords.x0 + w
            coords.x2 = coords.x1 + self._splitter_width
            self._enforceSplitterConstraints(coords)

        # Is the top row resizable?
        if (self._child00.isResizable(VERTICAL) and self._child01.isResizable(VERTICAL)):
            # Bottom row has priority:
            # Determine suggested width of the bottom row by asking
            # child10 and child11 in both orders and taking the maximum
            h = coords.y3-coords.y2
            w,h1 = self._child10.applyConstraint(1,h, interactive)
            w,h1 = self._child11.applyConstraint(1,h1, interactive)
            w,h2 = self._child11.applyConstraint(1,h, interactive)
            w,h2 = self._child10.applyConstraint(1,h2, interactive)
            h = max(h1,h2)
            # Set the final splitter y position 
            coords.y2 = coords.y3 - h
            coords.y1 = coords.y2 - self._splitter_width
            self._enforceSplitterConstraints(coords)
        else:
            # Top row has priority:
            # Determine suggested height of the top row by asking
            # child00 and child01 in both orders and taking the maximum
            h = coords.y1-coords.y0
            w,h1 = self._child00.applyConstraint(1,h, interactive)
            w,h1 = self._child01.applyConstraint(1,h1, interactive)
            w,h2 = self._child01.applyConstraint(1,h, interactive)
            w,h2 = self._child00.applyConstraint(1,h2, interactive)
            h = max(h1,h2)
            # Set the final splitter y position 
            coords.y1 = coords.y0 + h
            coords.y2 = coords.y1 + self._splitter_width
            self._enforceSplitterConstraints(coords)
        

    def _enforceSplitterConstraints(self, coords):
        """Enforces the pixel coordinate constraints on both of the splitter.

        Enforces the following constraints x0 <= x1 <= x2 <= x3 and
        y0 <= y1 <= y2 <= y3 provided that x0 <= x3 and y0 <= y3 already
        hold. The distance x2-x1 resp. y2-y1 is kept intact, unless x2<x1
        or y2<y1 in which case both values will be the same. If a splitter
        is not present in a particular direction, the corresponding splitter
        coordinates will be fixed at x3 resp. y3.

        The logical coordinate remains unchanged.

        This method is called whenever the splitter position changes.

        \pre x0 <= x3 and y0 <= y3
        """

        if self._splitter_type & HORIZONTAL == 0:
            coords.x1 = coords.x2 = coords.x3
        if self._splitter_type & VERTICAL == 0:
            coords.y1 = coords.y2 = coords.y3
            
        # d is the current width of the splitter
        d = coords.x2-coords.x1
        if d<0:
            coords.x2 = coords.x1
            d = 0
            
        if coords.x1<coords.x0:
            coords.x1 = coords.x0
            coords.x2 = coords.x1 + d
        if coords.x2>coords.x3:
            coords.x2 = coords.x3
            coords.x1 = max(coords.x0, coords.x2-d)

        # d is the current width of the splitter
        d = coords.y2-coords.y1
        if d<0:
            coords.y2 = coords.y1
            d = 0
            
        if coords.y1<coords.y0:
            coords.y1 = coords.y0
            coords.y2 = coords.y1 + d
        if coords.y2>coords.y3:
            coords.y2 = coords.y3
            coords.y1 = max(coords.y0, coords.y2-d)

    def _fillEmptyChildren(self):
        if self._child00==None:
            self.setChild(PanelNode(), 0)
        if self._child01==None:
            self.setChild(PanelNode(), 1)
        if self._child10==None:
            self.setChild(PanelNode(), 2)
        if self._child11==None:
            self.setChild(PanelNode(), 3)

# PanelNode
class PanelNode(ILayoutNode):
    """This class represents a leaf node in the layout tree.

    Each node in the layout tree must be a %PanelNode object. This object
    represents one panel is associated with a wx widget. The class
    has the following attributes:

    - \a constraint (\c SizeConstraint): A constraint on the size of the panel
    - \a widget (\c PanelWidget): The associated widget
    - [\a wx (\c wx \c widget): The associated wx widget]
    """
    
    def __init__(self, name="", activatable=False, constraint=None, widget=None):
        """Constructor.

        \param name (\c str) Node name
        \param activatable (\c bool) A flag that determines if the panel can
               be activated or not (only 3D views should be made activatable)
        \param constraint (\c SizeConstraint) Size constraint
        \param widget (\c PanelWidget) The panel widget or a sequence of widgets or None
        """
        ILayoutNode.__init__(self, name=name)

        # Position and size of the panel
        self._x0 = 0
        self._y0 = 0
        self._width = 0
        self._height = 0

        # Activatable flag
        self._activatable = activatable

        # Panel widgets. This is actually a stack where the last widget
        # is on top (all other widgets are hidden).
        self._widgets = []

        # Hide all widgets
        for w in self._widgets:
            w.hide()

        # Push the widgets...
        if widget!=None:
            try:
                wlst = list(widget)
            except TypeError:
                wlst = [widget]
            
            for w in wlst:
                self.pushWidget(w)

        self.showConfigPanel()

        # Constraint
        if constraint==None:
            constraint = NoConstraint()
        self._constraint = constraint

    def __str__(self):
        if self.name!=None:
            name = '"%s"'%self.name
        else:
            name = "<None>"
        return 'PanelNode %s active:%d widgets:%s'%(name, self.isActive(), [str(x) for x in self._widgets])

    def setGeom(self, x, y, width, height):
#        print "Panel '%s', setGeom(%d, %d, %d, %d)"%(self.name, x, y, width, height)
        self._x0 = x
        self._y0 = y
        self._width = width
        self._height = height

        x+=1
        y+=1
        width-=2
        height-=2
        self._constraint.setSize(width, height)
        if len(self._widgets)>0:
            self._widgets[-1].setGeom(x,y,width,height)

    def activate(self, root):
        # Store root (Panels object)
        self.setRoot(root)

        # Mark all widgets as activated
        for w in self._widgets:
            if not isinstance(w, PanelNodeDlg):
                self._root.repository.insert(w)
            w.activate(root)

        # Set the position and size of the panel widget
        self.setGeom(self._x0, self._y0, self._width, self._height)
        # Connect and show the top widget
        if len(self._widgets)>0:
            w = self._widgets[-1]
            self._connectPanelEvents(w.wx)
            w.show()

    def deactivate(self):
        # Disconnect and hide the top widget
        if len(self._widgets)>0:
            w = self._widgets[-1]
            w.hide()
            self._disconnectPanelEvents(w.wx)
        
        # Mark all widgets as inactivated
        for w in self._widgets:
            w.deactivate()

        self.setRoot(None)


    def pushWidget(self, widget):
        """Add a widget on top of the stack.

        Adds a widget on top of the stack. If the layout is active the
        new widget will be displayed.

        \param widget (\c PanelWidget) The new panel widget
        """
        if widget.isUsedBy(self._layoutRoot()):
            raise LayoutError('The widget "%s" is already managed by this layout.'%widget.name)

        # Hide and disconnect the previous top widget...
        if self.isActive():
            if len(self._widgets)>0:
                w = self._widgets[-1]
                w.hide()
                self._disconnectPanelEvents(w.wx)

        # Store the widget...
        widget.acquire(self)
        self._widgets.append(widget)

        if self.isActive():
            if not isinstance(widget, PanelNodeDlg):
                self._root.repository.insert(widget)
            # activate the widget...
            widget.activate(self._root)
            # Set the position and size of the panel widget
            self.setGeom(self._x0, self._y0, self._width, self._height)
            # Connect and show the top widget
            if len(self._widgets)>0:
                w = self._widgets[-1]
                self._connectPanelEvents(w.wx)
                w.show()

    def popWidget(self):
        """Remove the top widget from the stack.

        \pre The stack must not be empty
        """
        w = self._widgets[-1]
        self._widgets.pop()
        w.release(self)
        # Disconnect and hide the top widget and show and connect the
        # previous one
        if self.isActive():
            w.deactivate()
            w.hide()
            self._disconnectPanelEvents(w.wx)
            if len(self._widgets)>0:
                w2 = self._widgets[-1]
                self._connectPanelEvents(w2.wx)
                w2.show()
                self.setGeom(self._x0, self._y0, self._width, self._height)
        

    def findNodeByName(self, name):
        if self._name==name:
            return self
        else:
            return None

    def showConfigPanel(self):
        # Is a config panel already on top?
        if len(self._widgets)>0 and isinstance(self._widgets[-1], PanelNodeDlg):
            return
        # Remove any config panel that's somewhere inside the stack
        self._widgets = [w for w in self._widgets if not isinstance(w, PanelNodeDlg)]
        # Create a new config panel and push it onto the stack
        dlg = PanelNodeDlg(self)
        self.pushWidget(dlg)

    def hideConfigPanel(self):
        # If there's only the config panel present then leave it there
        if len(self._widgets)==1 and isinstance(self._widgets[0], PanelNodeDlg):
            return
        
        # Is a config panel on top?
        if len(self._widgets)>0 and isinstance(self._widgets[-1], PanelNodeDlg):
            # pop the panel
            self.popWidget()
        else:
            # Remove any config panel that's somewhere inside the stack
            self._widgets = [w for w in self._widgets if not isinstance(w, PanelNodeDlg)]
            
    def isResizable(self, direction):
        return self._constraint.isResizable(direction)

    def applyConstraint(self, width, height, interactive=False):
        w,h = self._constraint(width-2, height-2, interactive)
        w+=2
        h+=2
        return w,h

    def panelNodes(self):
        return [self]

#    def findPanel(self, x, y):
#        return (0, self)

    def makeCurrent(self):
        """Make this panel the current panel.
        """
        if self._root!=None and self._activatable:
            self._root.onClickPanel(self)

    def split(self, x=None, y=None):
        """Split this panel in two or four new panels.

        \param x (\c float) Split coordinate in x direction (0-1) or None
        \param y (\c float) Split coordinate in y direction (0-1) or None
        \return Newly created layout node (\c LayoutNode) that contains
               the original panel and the new empty ones.
        """
        if x==None and y==None:
            return
        elif x!=None and y==None:
            stype = HORIZONTAL
            childs = (self, None)
        elif x==None and y!=None:
            stype = VERTICAL
            childs = (self, None)
        else:
            stype = HORIZONTAL | VERTICAL
            childs = (self, None, None, None)

        node = LayoutNode(splittertype=stype)

        parent = self._parent
        if parent!=None:
            idx = parent.childIndex(self)
            parent.removeChild(self)
            parent.setChild(node, idx)
        else:
            if self._root!=None:
                self._root.layout = node

        node.setChildren(childs)
        if self._root!=None:
            self._root.updateLayout()

        return node


    ######################################################################
    ## protected:

    def _connectPanelEvents(self, wxwin):
        """Connect the given wx window to the Panels object.

        This method sets the connections so that the Panels object will
        get notified when the mouse enters or leaves this panel and
        when it gets clicked.

        \param wxwin (\c wx \c window) A wx.Window instance
        """
#        print "connect",wxwin
        wx.EVT_ENTER_WINDOW(wxwin, self._onEnter)
        wx.EVT_LEAVE_WINDOW(wxwin, self._onLeave)
        wx.EVT_LEFT_DOWN(wxwin, self._onClick)
        wx.EVT_MIDDLE_DOWN(wxwin, self._onClick)
        wx.EVT_RIGHT_DOWN(wxwin, self._onClick)

    def _disconnectPanelEvents(self, wxwin):
        """Disconnect the given wx window from the Panels object.

        \param wxwin (\c wx \c window) A wx.Window instance
        """
#        print "disconnect",wxwin
        wx.EVT_ENTER_WINDOW(wxwin, None)
        wx.EVT_LEAVE_WINDOW(wxwin, None)
        wx.EVT_LEFT_DOWN(wxwin, None)
        wx.EVT_MIDDLE_DOWN(wxwin, None)
        wx.EVT_RIGHT_DOWN(wxwin, None)

    def _onEnter(self, event):
        self._root.onEnterPanel(self)
#        print "onEnter", self.name, self._root
        event.Skip()

    def _onLeave(self, event):
        self._root.onLeavePanel(self)
#        print "onLeave", self.name, self._root
        event.Skip()

    def _onClick(self, event):
        if self._activatable:
            self._root.onClickPanel(self)
        event.Skip()

    # "widget" property...
    
    def _getWidget(self):
        """Return the top panel widget.

        This method is used for retrieving the \a widget property.

        \return PanelWidget object or None.
        """
        if len(self._widgets)==0:
            return None
        return self._widgets[-1]

    def _setWidget(self, widget):
        """Set the top panel widget.

        This method is used for setting the \a widget property.

        \param widget (\c PanelWidget) Widget
        """
        if len(self._widgets)>0:
            self.popWidget()
        self.pushWidget(widget)

    widget = property(_getWidget, _setWidget, None, "Associated panel widget")

    # "constraint" property...
    
    def _getConstraint(self):
        """Return the constraint object.

        This method is used for retrieving the \a constraint property.

        \return Constraint object (\c SizeConstraint).
        """
        return self._constraint

    def _setConstraint(self, constraint):
        """Set the constraint object.

        This method is used for setting the \a constraint property.

        \param constraint (\c SizeConstraint) Constraint object
        """
        self._constraint = constraint

    constraint = property(_getConstraint, _setConstraint, None, "Size constraint")
    
        
        
# NoConstraint
class NoConstraint:
    """An unconstrained size constraint class."""
    
    def __call__(self, width, height, interactive=False):
        return width, height

    def setSize(self, width, height):
        pass

    def isResizable(self, direction):
        return True


# SizeConstraint
class SizeConstraint:
    """This class is used to enforce size constraints on widgets/panels.

    This class computes acceptable width/height values that conform
    to the current constraints. The constraints can either keep the
    width or height constant or at a multiple of a base width/height.
    The constraint class is a callable object that can be called
    with an input width and height and it returns a new width/height
    tuple that conforms to the constraint and is close to the input
    values.

    Examples:

    \code
    # A constraint that keeps the width fixed at 270 pixel, the height is arbitrary
    >>> c = SizeConstraint(width=270, wfixed=True)
    >>> c(100,100)
    (270, 100)

    # Now keep the height at a multiple of 32
    >>> c = SizeConstraint(height=32)
    >>> c(100,100)
    (100, 96)    
    \endcode
    
    """
    
    def __init__(self, width=None, height=None):
        """Constructor.

        \param width  A tuple (width, step, resizable) describing the width constraint
        \param height  A tuple (height, step, resizable) describing the height constraint
        """
        if width==None:
            width = None, None, False
        if height==None:
            height = None, None, False
            
        self.width, self.wstep, self.wresizable = width
        self.height, self.hstep, self.hresizable = height

    def __call__(self, width, height, interactive=False):
        """Check if the constraint accepts the given size.

        This method has to check if the given width and height
        are acceptable for the constraint. If they are the method
        just has to return those values again, if they are not
        it has to return a suggested size that would be acceptable.

        \param width (\c int) A width value
        \param height (\c int) A height value
        \param interactive (\c bool) If True the constraints also accepts multiples of its width (default: False)
        \return Tuple (width, height) specifying an acceptable size that's
               as close as possible to the input size
        """
        if self.width!=None:
            w = self.width
        else:
            w = width
            
        if self.height!=None:
            h = self.height
        else:
            h = height
        
        if self.wresizable and interactive:
            dw = width-self.width
            if dw%self.wstep==0:
                w = width
            else:
                margin = int(0.15*self.wstep)
                w = self.width + int((dw+margin)/self.wstep)*self.wstep

        if self.hresizable and interactive:
            dh = height-self.height
            if dh%self.hstep==0:
                h = height
            else:
                margin = int(0.15*self.hstep)
                h = self.height + int((dh+margin)/self.hstep)*self.hstep

        return w,h

    def setSize(self, width, height):
        """Set a new width and height.

        The new size is not set directly but run through the
        constraint (in interactive mode). So the actual size that
        gets set is a size that's compliant to the constraint.
        """
        w, h = self(width, height, True)
        if self.width!=None:
            self.width = w
        if self.height!=None:
            self.height = h

    def isResizable(self, direction):
        """Check if the size of the constraint may be changed.

        \param direction (\int) A combination of HORIZONTAL and VERTICAL
        \return True if resizable.
        """
        res = False
        if direction&HORIZONTAL:
            res |= self.wresizable
        if direction&VERTICAL:
            res |= self.hresizable
    
######################################################################


class PanelNodeDlg(PanelWidget):
    def __init__(self, panelnode):
        """Constructor.

        \param parentwin (\c MainWindow) Parent window (must have a "panels.wx" attribute)
        \param name (\c str) Widget name
        """
        PanelWidget.__init__(self)
        self.panelnode = panelnode
    
    def activate(self, root):
        PanelWidget.activate(self, root)
        self.wx = wx.Panel(root.wx, -1, style=wx.SIMPLE_BORDER)
        
        self.editlabel = wx.StaticText(self.wx, -1, "Name:")
        self.nameedit = wx.TextCtrl(self.wx, -1, value=self.panelnode.name)

        self.splitbox = wx.StaticBox(self.wx, -1, "Split")
        self.hsplitbtn = wx.BitmapButton(self.wx, -1, panelicons.getHSplitBitmap())
        self.hsplitbtn.SetToolTip(wx.ToolTip("Split this panel horizontally"))
        self.vsplitbtn = wx.BitmapButton(self.wx, -1, panelicons.getVSplitBitmap())
        self.vsplitbtn.SetToolTip(wx.ToolTip("Split this panel vertically"))
        self.hvsplitbtn = wx.BitmapButton(self.wx, -1, panelicons.getHVSplitBitmap())
        self.hvsplitbtn.SetToolTip(wx.ToolTip("Split this panel in both directions"))

        self.toplayout = wx.BoxSizer(wx.VERTICAL)

        self.s1 = wx.BoxSizer(wx.HORIZONTAL)
        self.s1.Add(self.editlabel, 0, wx.ALL | wx.ALIGN_CENTRE, 4)
        self.s1.Add(self.nameedit, 1, wx.ALL | wx.ALIGN_CENTRE, 4)

        self.s2 = wx.StaticBoxSizer(self.splitbox, wx.HORIZONTAL)
        self.s2.Add(self.hsplitbtn, 0, wx.ALL, 2)
        self.s2.Add(self.vsplitbtn, 0, wx.ALL, 2)
        self.s2.Add(self.hvsplitbtn, 0, wx.ALL, 2)

        self.toplayout.Add((0,1), 1, 0,0)
        self.toplayout.Add(self.s1, 0, wx.ALIGN_CENTRE, 0)
        self.toplayout.Add(self.s2, 0, wx.ALIGN_CENTRE, 0)
        self.toplayout.Add((0,1), 1, 0,0)

        self.wx.SetSizer(self.toplayout)

        wx.EVT_TEXT_ENTER(self.wx, self.nameedit.GetId(), self.onNameEntered)
        wx.EVT_KILL_FOCUS(self.nameedit, self.onNameEntered)
        wx.EVT_BUTTON(self.wx, self.hsplitbtn.GetId(), self.onHSplit)
        wx.EVT_BUTTON(self.wx, self.vsplitbtn.GetId(), self.onVSplit)
        wx.EVT_BUTTON(self.wx, self.hvsplitbtn.GetId(), self.onHVSplit)

    def onNameEntered(self, event):
        newname = self.nameedit.GetValue()
        try:
            self.panelnode.name = newname
        except DuplicateNames:
#            dlg = wx.MessageDialog(self.wx,
#                                   'There is already a layout node called "%s".\nPlease choose another name.'%newname, "Duplicate names", wx.OK)
#            dlg.ShowModal()
            newname = self.panelnode._layoutRoot().makeNameUnique(newname)
            self.panelnode.name = newname
            self.nameedit.SetValue(newname)
            self.nameedit.SetSelection(-1,-1)
 #           self.nameedit.SetFocus()
        
    def onHSplit(self, event):
        self.panelnode.split(x=0.5)

    def onVSplit(self, event):
        self.panelnode.split(y=0.5)

    def onHVSplit(self, event):
        self.panelnode.split(x=0.5, y=0.5)


######################################################################

class WidgetDlg(wx.Frame):
    def __init__(self, parent, id):
        wx.Frame.__init__(self, parent, id, "Panel widgets", style=wx.DEFAULT_FRAME_STYLE | wx.FRAME_FLOAT_ON_PARENT)

        self.lst = wx.ListCtrl(self, -1, style=wx.LC_REPORT | wx.LC_SINGLE_SEL | wx.LC_HRULES | wx.LC_VRULES)
        self.lst.InsertColumn(0, "Panel")
        self.lst.InsertColumn(1, "Name")
        self.lst.InsertColumn(2, "Ref.count")
        self.lst.InsertStringItem(0, "Shell")
        self.lst.SetStringItem(0,1,"shell")
        self.lst.SetStringItem(0,2,"1")
        self.lst.InsertStringItem(1, "Shader Editor")
        self.lst.SetStringItem(1,1,"shadereditor")
        self.lst.SetStringItem(1,2,"1")
        self.lst.InsertStringItem(2, "3D Viewport")
        self.lst.SetStringItem(2,1,"front")
        self.lst.SetStringItem(2,2,"1")
        self.lst.InsertStringItem(3, "3D Viewport")
        self.lst.SetStringItem(3,1,"perspective")
        self.lst.SetStringItem(3,2,"1")

def widgetDlg():
    w = WidgetDlg(getApp().window.wx, -1)
    w.Show()

cgkit.cmds.widgetDlg = widgetDlg

######################################################################

# createPanelWidget
def createPanelWidget(name, modname=None):
    """Looks for a panel widget plugin and creates an instance of it.

    \param name (\c str) Widget name
    \param modname (\c str) Module name or None (which matches every module)
    """

    for odesc in pluginmanager.iterProtoObjects(PanelWidget):
        if odesc.name==name and (modname==None or odesc.plugindesc.modulename==modname):
            ClassName = odesc.object
            widget = ClassName()
            return widget
        
    return None

cgkit.cmds.createPanelWidget = createPanelWidget

######################################################################

if __name__=="__main__":

    c = SizeConstraint(5,5, wfixed=True)
    print((c(172,22)))
