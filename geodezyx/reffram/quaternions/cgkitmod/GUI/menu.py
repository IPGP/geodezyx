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

## \file menu.py
## Contains the menu tree classes.

import wx
import types


class MenuNodeNotFound(Exception):
    """Exception class."""
    pass

class DuplicateNames(Exception):
    """Exception class."""
    pass

class InvalidMenuRoot(Exception):
    """Exception class."""
    pass

class NoWxMenuObjectFound(Exception):
    """Exception class."""
    pass

class MenuAlreadyInUse(Exception):
    """Exception class."""
    pass

class NoMenuRoot(Exception):
    """Exception class."""
    pass

NODE_NORMAL_ITEM  = 0x01
NODE_CHECK_ITEM   = 0x02
NODE_COMMAND_ITEM = NODE_NORMAL_ITEM | NODE_CHECK_ITEM
NODE_SEPARATOR    = 0x04
NODE_MENU         = 0x08
NODE_MARKER       = 0x10
NODE_ALL_NODES    = 0xff


# MenuNode
class MenuNode(object):
    """Base class for a menu node.

    This is the base class for all menu nodes (inner node or leaf).
    A %MenuNode class can be used as a sequence where the items are
    the children nodes. Leaf nodes are always empty.
    
    Each node has the following attributes:

    - \b name (\c str): Node name
    - \b absname (\c str): Absolute node name
    - \b parent (\c %MenuNode): Parent node
    - \b root (\c %MenuNode): Root node
    - \b text (\c str): Text that's visible to the user
    - \b enabled (\c bool): If this flag is false the menu or menu item is grayed out.
    - \b visible (\c bool): If this flag is false the menu or menu item isn't displayed at all.
    - \b wx (\c wxMenu, \c wxMenuBar or \c wxMenuItem): Corresponding wxPython object (or None if the menu isn't attached).
    
    - flags

    A menu node can be either in an \em attached or \em detached state
    which tells if a node is attached to a wxPython menu or not. An
    attached menu is "active" and usually has a corresponding wxPython
    menu object associated with it. When you create a menu node it's
    not attached to a menu. A node automatically becomes attached when
    it's added to another node that is already attached.

    Any class that is derived from this base class has to implement the
    following methods:

    - setEnabled()
    - setVisibility()
    - _wx_menu_count()
    - _create_wx()
    
    """
    
    def __init__(self, name=None, parent=None, text=None,
                 enabled=True, visible=True):
        """Constructor.

        If no name is given the constructor tries to create a name
        from the \a text argument (using _create_auto_name()).

        \param name (\c str) Node name
        \param parent (\c %MenuNode) Parent node
        \param text (\c str) %Menu title or item text
        \param enabled (\c bool) Initial enabled state
        \param visible (\c bool) Initial visibility state
        """
        # Node name
        self._name = name
        # Parent node
        self._parent = parent
        # Text (menu title or item text)
        self._text = text
        # NODE_xxx flags
        self._flags = 0
        # Enabled state
        self._enabled = enabled
        # Visbility state
        self._visible = visible
        # The corresponding wxPython object (or None)
        self._wx = None
        # The wxPython ID of the object
        self._id = None
        # Activation flag
        self._attached = False

        if self._name==None:
            self._create_auto_name(self._text)


    # isAttached
    def isAttached(self):
        return self._attached

    def setEnabled(self, flag=True):
        """Set the enabled flag of this node.

        If a node is disabled all of its children are considered disabled
        as well.

        This method should be overridden by derived classes. Those
        implementations should call the inherited method (where the
        internal flag _enabled is set). If the node is attached to a
        menu this method should update their corresponding menu items
        to reflect the state change.

        This method is called automatically when the \em enabled property
        receives a new value.

        \param flag (\c bool) New enabled state.
        """
        self._enabled = flag

    def setVisibility(self, flag=True):
        """Set the visibility flag of this node.

        If a node is invisible all of its children are considered invisible
        as well.

        This method should be overridden by derived classes. Those
        implementations should call the inherited method (where the
        internal flag _visible is set). If the node is attached to a
        menu this method should remove their corresponding menu items
        to reflect the state change.

        This method is called automatically when the \em visible property
        receives a new value.

        \param flag (\c bool) New visibility state.
        """
        
        self._visible = flag

        if flag:
            if not self.isAttached():
                self._create_wx()
        else:
            if self.isAttached():
                self._destroy_wx()


    # remove
    def remove(self, node=""):
        """Remove a menu node.

        Removes the node with the given name. The name may contain an
        entire (relative) path to a menu node or may be left blank if
        self should be removed (which is the default behavior). The node
        is removed from the visible wx menu and from the menu tree and
        is returned to the caller which may ignore it or use it to add
        it later or to another menu.

        \param node (\c str) Name of the node which should be removed (default: "").
        \return The removed menu node object (\c MenuNode).
        """

        # Find the node n which is to be removed
        n = self.findNode(node)

        # Remove the corresponding wx objects from the visible menu
        if n.isAttached():
            n._destroy_wx()

        parent = n.parent
        if parent!=None:
            del parent[n]

        return n

    def findCommandNodes(self, func):
        """Return all nodes that are bound to func.

        \return A list of all nodes that are bound to func.
        """

        res = []
        
        # Check if this node is bound to func
        f = getattr(self, "_command", None)
        if f!=None:
            if f==func:
                res.append(self)

        # Now check the children
        for n in self:
            res+=n.findCommandNodes(func)

        return res                
        

    # findNode
    def findNode(self, path, exc=True):
        """Search for a node with a specific relative name.

        The name can include an entire path. Each node name is separated
        by a dot. If path is empty then self is returned.

        \param path (\c str) Name of the item
        \param exc (\c bool) If True then an exception is generated when the
                             node is not found.
        \return Node or None/Exception.
        """

        if path=="":
            return self

        f = path.split(".")
        name = f[0]
        for n in self:
            if name==n.name:
                if len(f)==1:
                    return n
                else:
                    return n.findNode(".".join(f[1:]))

        if exc:
            errpath = name
            s = self.absname
            if s!="":
                errpath = ".".join([s,errpath])
            raise MenuNodeNotFound('Menu node "%s" not found.'%(errpath))
        else:
            return None
            

    ######################################################################
    ## protected:

    def __len__(self):
        return 0

    def __getitem__(self, key):
        if isinstance(key, int):
            raise IndexError("sequence index out of range")
        else:
            raise TypeError("sequence indices must be integers")

    # _create_auto_name
    def _create_auto_name(self, text):
        """Create a name from the given text.

        The argument \a text is the text that's visible to the user
        in the menu. The method tries to generate a node name out
        of the given text. If the text only consists of dashes or is None
        then no name is generated (None), otherwise blanks and tabs are
        collapsed and replaced by underscores, dots and ampersands are
        deleted and the text is made lower case to create the name.

        The result of the automatic name generation is stored in the
        name property (_name).
        
        \param text (\c str) Menu text.
        """
        if text==None:
            self._name = None
            return
        text = text.strip()
        # Collapse dash sequences
        while text.find("--")!=-1:
            text = text.replace("--", "-")
        # Was the text that of a separator? then no name is generated
        if text=="-":
            self._name = None
            return
        # Replace tabs with spaces
        text = text.replace("\t", " ")
        # Collapse spaces
        while text.find("  ")!=-1:
            text = text.replace("  ", " ")
        # Replace/delete special characters
        text = text.replace("&","")
        text = text.replace(".","")
        text = text.replace(" ","_")
        self._name = text.lower()

    def _wx_menu_count(self):
        """Return the number of menu "slots" that are occupied by this node.

        The return value is the number of attached menu items (or
        menus) that are beneath this node. Invisble nodes must not be
        counted. This value is used to determine the position of a
        node in the wx menu.

        By default this method returns 0. Usually, a derived class should
        override this method.

        \return Number of visible items occupied by this node.
        \see _wx_menu_pos()
        """
        return 0

    def _wx_menu_pos(self):
        """Return the position of the node inside the parent wx menu.

        The return value is the starting index where self is located
        (or would be located) in the corresponding wx menu. This value
        is used to insert the node into a menu or menubar.

        The node has to be added to the menu tree but doesn't have to
        be attached.
        
        \pre self is in the children list of parent.
        \pre The parent nodes must already be attached.
        \return Position (\c int) or None if the node doesn't have a parent.
        \see _wx_menu_count()
        """

        parent = self.parent
        if parent==None:
            return None

        # Determine the index of self within parent
        idx = 0
        for n in parent:
            if n==self:
                break
            idx += n._wx_menu_count()

        if parent._wx==None:
            idx += parent._wx_menu_pos()

        return idx

        

    def _attach(self, IDmanager, window, wxmenubar):
        """Attach the menu to a window.
        
        \param IDmanager ID generator
        \param window (\c Window) wxPython window
        \param wxmenubar (\c wxMenuBar) wxPython menu bar or None
        """
        root = self.root
        if not isinstance(root, Menu):
            raise InvalidMenuRoot('The root of a menu has to be a Menu object.')
        # The node must not be attached to another menu
        if self.isAttached():
            raise MenuAlreadyInUse("The menu node is already in use.")

        root._IDmanager = IDmanager
        root._window = window
        if wxmenubar!=None:
            root._wx_menubar = wxmenubar
        root._menu_id_lut = {}
        root._create_wx()

    def _detach(self):
        """Detach the menu from a window."""
        
        if self.root!=self:
            raise NoMenuRoot("Only root nodes can be detached.")

        self._destroy_wx()
        del self._IDmanager
        del self._window
        if hasattr(self, "_wx_menubar"):
            del self._wx_menubar
        del self._menu_id_lut
        

    def _create_wx(self):
        """Create the corresponding wxPython object.

        This method has to create the appropriate wxPython object which
        is represented by the node and insert it into the parent menu,
        submenu or menubar at the correct position as returned by
        _wx_menu_pos(). 
        
        The wxPython object has to be stored in the attribute self._wx
        and its ID in self._id.
        The method may only be called if the parent wxPython object
        of the node exists or if the node represents the entire menu bar.

        A derived class must call this inherited method!
        """
        self._attached = True

    def _destroy_wx(self):
        """Destroy the associated wxPython object.

        This method has to remove its corresponding wxPython object
        from the visible wx menu. The method may only be called in
        attached state!

        A derived class must call this inherited method!
        """
        self._attached = False
        self._wx = None
        self._id = -1        

    def _wx_parent(self):
        """Return the wxPython parent object.

        The object returned is either a wx.Menu, a wx.MenuBar or None.
        If a node has no direct wx parent, the wx parent of the
        parent node is returned.

        This method must not be called when the node is part of an
        invisible tree.

        \return wxPython object or None.
        """
        p = self.parent

        while p!=None:
            if p._wx!=None:
                if isinstance(p._wx, wx.MenuItem):
                    return p._wx.GetSubMenu()
                else:
                    return p._wx
            p = p.parent

        return None


    def _is_menubar(self):
        return False
    
    def _is_embedded_menu(self):
        return False


    # "name" property...
    
    def _getName(self):
        """Return the node name.

        This method is used for retrieving the \a name property.

        \return Node name (\c str).
        """
        return self._name

    name = property(_getName, None, None, "Node name")

    # "absname" property...

    def _getAbsName(self):
        """Return the absolute name of the node.

        Returns the absolute name of the node. This name
        uniquely identifies the node within its tree and is composed
        of all the names from the root to the node separated by a dot.
        The root node is always unnamed.

        This method is used for retrieving the \a absname property.

        \return Absolute name of the node (\c str).
        """
        if self._parent==None:
            return ""
        else:
            pre = self._parent.getFullName()
            if self._name==None:
                name = "<unnamed>"
            else:
                name = self._name
            if pre=="":
                return name
            else:
                return pre+"."+name

    absname = property(_getAbsName, None, None, "Absolute node name")

    # "parent" property...
    
    def _getParent(self):
        """Return the parent node.

        This method is used for retrieving the \a parent property.

        \return Parent node (\c MenuNode).
        """
        return self._parent

    parent = property(_getParent, None, None, "Parent node")

    # "root" property...

    def _getRoot(self):
        """Return the root node.

        This method is used for retrieving the \a root property.

        \return Root node (\c MenuNode).
        """
        p = self
        while 1:
            if p._parent==None:
                return p
            p = p._parent

    root = property(_getRoot, None, None, "Root node")

    # "text" property...
    
    def _getText(self):
        """Return the menu text which is displayed to the user.

        This method is used for retrieving the \a text property.

        \return Text (\c str).
        """
        return self._text

    text = property(_getText, None, None, "Menu text")

    # "enabled" property...
    
    def _getEnabled(self):
        """Return the enabled state.

        This method is used for retrieving the \a enabled property.

        \return True if the node is enabled.
        """
        return self._enabled

    def _setEnabled(self, flag=True):
        """Set the enabled state.

        This method is used for setting the \a enabled property.

        \param flag (\c bool) New enabled state.
        """
        self.setEnabled(flag)

    enabled = property(_getEnabled, _setEnabled, None, "Enabled state")

    # "visible" property...
    
    def _getVisible(self):
        """Return the visible state.

        This method is used for retrieving the \a visible property.

        \return True if the node is visible.
        """
        return self._visible

    def _setVisible(self, flag=True):
        """Set the visible state.

        This method is used for setting the \a visible property.

        \param flag (\c bool) New visibility state.
        """
        self.setVisibility(flag)

    visible = property(_getVisible, _setVisible, None, "Visibility state")

    # "wx" property...
    
    def _getWx(self):
        """Return the corresponding wxPython object.

        This method is used for retrieving the \a wx property.

        \return wxPython object.
        """
        return self._wx

    wx = property(_getWx, None, None, "Corresponding wxPython object")



######################################################################
######################################################################
        
# MenuMarker
class MenuMarker(MenuNode):
    """This class represents a position in a menu.
    """
    def __init__(self, name=None, parent=None):
        MenuNode.__init__(self, name=name, parent=parent)
        self._flags = NODE_MARKER

    def __repr__(self):
        return "<MenuMarker name=%s>"%(self.name)

######################################################################
######################################################################

# MenuItem
class MenuItem(MenuNode):
    """An individual menu item.

    This class represents one individual menu item which can be a
    normal item, a checkable item or a separator. A menu item object
    has the following attributes in addition to the ones inherited by
    MenuNode:

    - \b command (\c callable): The command to execute (may be None)
    - \b checked (\c bool): The state of a checkable menu item
    """
    
    def __init__(self, itemdef, name=None, enabled=True, visible=True):
        """Constructor.

        The definition \a itemdef of the item is generally composed of
        two items, the text that appears in the menu and the
        corresponding command which is to be executed whenever the
        user selects the menu item. The command must be a callable
        that takes no arguments.
        
        The text and command is given as a 2-tuple (text, command).
        It's also possible to just provide a string which is equivalent
        to (text, None).

        If the text starts with "[x]" then a checkable menu item is
        generated.

        If the text consists of 3 or more dashes ("---"), then a
        separator is generated.

        \param itemdef (\c 2-tuple or \c str) Item definition
        \param name (\c str) The internal name of the item
        \param enabled (\c bool) Initial enabled state
        \param visible (\c bool) Initial visibility state
        """

        flags = NODE_NORMAL_ITEM

        self._checkitem = False
        self._separator = False
        self._checked   = False

        # Split the item definition into the text and command component
        if isinstance(itemdef, (str,)):
            itemstr = itemdef
            itemcmd = None
        else:
            itemstr, itemcmd = itemdef

        # Is the menu item a checkable item?
        if itemstr[:3]=="[x]" or itemstr[:3]=="[X]":
            flags = NODE_CHECK_ITEM
            self._checkitem = True
            if itemstr[1]=="X":
                self._checked=True
            itemstr = itemstr[3:]

        # Is the item a separator?
        if len(itemstr)>2 and itemstr.count("-")==len(itemstr):
            flags = NODE_SEPARATOR
            self._separator = True
            itemcmd = None

        MenuNode.__init__(self, name=name, text=itemstr, enabled=enabled, visible=visible)
        self._flags = flags

        self._command   = itemcmd

    def __repr__(self):
        s = 'MenuItem name=%s text="%s"'%(self.name, self.text)
        if self._checkitem:
            s+=" checkable"
        return "<%s>"%s

    def isSeparator(self):
        return self._separator

    def isCheckItem(self):
        return self._checkitem


    def setEnabled(self, flag=True):
        MenuNode.setEnabled(self, flag)
        if self._wx!=None:
            self._wx.Enable(flag)

    def update(self):
        if self._wx!=None:
            root = self.root
            text = self.text
            keys = root._window.keys.findCommandKeys(self._command)
            if len(keys)>0:
                text+="\t"+keys[0]
            self._wx.SetText(text)


    ######################################################################

    ## protected:

    def _wx_menu_count(self):
        if self._wx==None:
            return 0
        else:
            return 1

    def _create_wx(self):

        MenuNode._create_wx(self)

        root = self.root

        # Determine the type of MenuItem to create...
        if self.isCheckItem():
            kind = wx.ITEM_CHECK
            self._id = root._IDmanager.allocID()
        elif self.isSeparator():
            kind = wx.ITEM_SEPARATOR
            self._id = wx.ID_SEPARATOR
        else:
            kind = wx.ITEM_NORMAL
            self._id = root._IDmanager.allocID()

        # Get the parent which must be a wxMenu object
        parent = self._wx_parent()
        # Create the MenuItem...
        text = self.text
        keys = root._window.keys.findCommandKeys(self._command)
        if len(keys)>0:
            text+="\t"+keys[0]
        self._wx = wx.MenuItem(parent, self._id, text, self._getHelpString(), kind)
        # ...and insert it into the parent menu
        pos = self._wx_menu_pos()
#        print 'Inserting item "%s" at position %d'%(self._text, pos)
        parent.InsertItem(pos, self._wx)
        # Set initial attributes
        self._wx.Enable(self.enabled)
        if self.isCheckItem():
            self._wx.Check(self.checked)

        # Connect the new wxMenuItem to the menu event handler
        if not self.isSeparator():
            root._connect(self)


    def _destroy_wx(self):

        root = self.root
        
        # Disconnect the menu item from the event handler
        if not self.isSeparator():
            root._disconnect(self)

        # Remove the wxMenuItem from the parent menu
        parent = self._wx_parent()
        parent.Remove(self._id)
        # Make the ID available again
        if self._id!=-1:
            root._IDmanager.freeID(self._id)

        MenuNode._destroy_wx(self)

    def _getHelpString(self):
        doc = getattr(self._command, "__doc__", None)
        if doc==None:
            return ""
        else:
            return doc.split("\n")[0]

    # "command" property...

    def _getCommand(self):
        """Return the associated command (callable object).

        This method is used for retrieving the \a command property.

        \return Command (\c callable).
        """
        return self._command

    def _setCommand(self, command):
        """Set a new command.

        This method is used for setting the \a command property.

        \param command (\c callable) The new command or None
        """
        self._command = command

    command = property(_getCommand, _setCommand, None, "Menu item command")

    # "checked" property...

    def _getChecked(self):
        """Get the checked flag.

        This method is used for retrieving the \a checked property.

        \return Checked flag (\c bool).
        """
        return self._checked
    
    def _setChecked(self, flag=True):
        """Get the checked flag.

        This method is used for setting the \a checked property.

        \param flag (\c bool) New state.
        """
        self._checked = flag
        if self._wx!=None:
            self._wx.Check(flag)

    checked = property(_getChecked, _setChecked, None, "Checked flag")


######################################################################
######################################################################
       
# Menu
class Menu(MenuNode):
    """A menu node containing other menus, markers or menu items.

    This class is used as a container for MenuNode classes and it
    represents a node in the entire menu tree.
    """
    
    def __init__(self, title=None, name=None, enabled=True, visible=True,
                 items=[]):
        """Constructor.

        \param title (\c str) Menu title which will be display in the GUI
        \param name (\c str) Node name
        \param enabled (\c bool) Initial enabled state
        \param visible (\c bool) Initial visibility state
        \param items (\c sequence) A sequence of children node descriptions
        """

        MenuNode.__init__(self, name=name, text=title, enabled=enabled, visible=visible)
        self._flags = NODE_MENU

        # Children nodes
        self._children = [MenuMarker("START", self), MenuMarker("END", self)]

        for item in items:
            if isinstance(item, MenuNode):
                mi = item
            else:
                mi = MenuItem(item)
            self.append(mi)


    def insertAfter(self, newnode, treenodename):
        """Insert a node after another existing node.

        \param newnode  The node that should be inserted
        \param treenodename The node name after which newnode is to be inserted
        """

        # The new node must not be attached to another menu
        if newnode.isAttached() or newnode.parent!=None:
            raise MenuAlreadyInUse("The menu node is already in use.")

        treenode = self.findNode(treenodename)
        # If treenode is not a direct child then delegate the insert
        # operation to the parent of treenode
        if (treenode.parent!=self):
            treenode.parent.insertAfter(newnode, treenodename.split(".")[-1])
            return

        # Check if the name of the new node collides with an existing node
        newnodename = newnode.name
        if newnodename!=None and self.findNode(newnodename, False)!=None:
            raise DuplicateNames('Menu node "%s" already exists.'%newnodename)

        # Add the new node into the children list
        idx = self._children.index(treenode)
        self._children.insert(idx+1, newnode)
        newnode._parent = self

        if self.isAttached():
            newnode._create_wx()


    def insertBefore(self, newnode, treenodename):
        """Insert a node before another existing node.

        \param newnode  The node that should be inserted
        \param treenodename The node name before which newnode is to be inserted
        """
        
        # The new node must not be attached to another menu
        if newnode.isAttached() or newnode.parent!=None:
            raise MenuAlreadyInUse("The menu node is already in use.")

        treenode = self.findNode(treenodename)
        # If treenode is not a direct child then delegate the insert
        # operation to the parent of treenode
        if (treenode.parent!=self):
            treenode.parent.insertBefore(newnode, treenodename.split(".")[-1])
            return

        # Check if the name of the new node collides with an existing node
        newnodename = newnode.name
        if newnodename!=None and self.findNode(newnodename, False)!=None:
            raise DuplicateNames('Menu node "%s" already exists.'%newnodename)

        # Add the new node into the children list
        idx = self._children.index(treenode)
        self._children.insert(idx, newnode)
        newnode._parent = self

        if self.isAttached():
            newnode._create_wx()

    # append
    def append(self, node):
        """Append a node at the end of the children list.

        The node is added right before the END marker.

        \param node (\c MenuNode) The node the be appended to the menu.
        """
        self.insertBefore(node, "END")

       

    def setEnabled(self, flag=True, node=""):
        if node!="":
            # Find the node n whose enabled state should be changed
            n = self.findNode(node)
            n.setEnabled(flag)
            return

        MenuNode.setEnabled(self, flag)
        if not self.isAttached():
            return

        # 1. Node = MenuBar
        if self._is_menubar():
            for n in self:
                n.setEnabled(flag)
        # 2. Node = Toplevel menu
        elif self.parent==None or self.parent._is_menubar():
            parent = self._wx_parent()
            # parent may be None if the menu is used for a popup menu
            if parent!=None:
                parent.EnableTop(self._wx_menu_pos(), flag)
        # 3. Node = Sub menu
        elif self.text!=None:
            parent = self._wx_parent()
            parent.Enable(self._wx.GetId(), flag)
        # 3. Node = Embedded menu
        else:
            for n in self:
                n.setEnabled(flag)              


    ######################################################################
    ## protected:

    def _wx_menu_count(self):
        if self._is_embedded_menu():
            res = 0
            for n in self:
                res += n._wx_menu_count()
            return res
        
        if self._wx==None:
            return 0
        else:
            return 1

    def __repr__(self):
        return '<Menu name=%s title="%s">'%(self.name, self.text)

    def __len__(self):
        return len(self._children)

    def __getitem__(self, key):
        return self._children[key]

    def __delitem__(self, key):
        # Remove the node from the tree
        self._children.remove(key)
        key._parent = None

    def __getattr__(self, name):
        items = [x for x in self._children if name==x.name]
        if items==[]:
            raise AttributeError('Menu node "%s" not found.'%(name))
        return items[0]

    def __setattr__(self, name, value):
        if name=="":
            return
        if name[0]=="_":
            MenuNode.__setattr__(self, name, value)
        elif isinstance(value, Menu):
            self.append(value)
        elif isinstance(value, MenuItem):
            self.append(value)
        elif (isinstance(value, (str,)) or
              (isinstance(value, tuple) and len(value)==2 and
              isinstance(value[0], (str,)) and callable(value[1]))):
            mi = MenuItem(value, name=name)
            self.append(mi)
        else:
            MenuNode.__setattr__(self, name, value)
            

##    def _getIndex(self, name):
##        """Return the index of the children with the given name.

##        \return Index (0-based).
##        """
        
##        items = filter(lambda x: name==x.name, self._children)
##        if items==[]:
##            raise MenuNodeNotFound, 'Menu node "%s" not found.'%(name)
##        return self._children.index(items[0])

    def _create_wx(self):

        MenuNode._create_wx(self)

        # 1. Node = MenuBar
        if self._is_menubar():
#            print "Menubar"
            self._create_wx_menubar()    
        # 2. Node = Toplevel menu
        elif self.parent==None or self.parent._is_menubar():
#            print "Toplevel"
            self._create_wx_toplevel_menu()
        # 3. Node = Sub menu
        elif self.text!=None:
#            print "Submenu"
            self._create_wx_submenu()
        # 3. Node = Embedded menu
        else:
#            print "Embedded"
            self._create_wx_embedded_menu()

        self.setEnabled(self.enabled)


    def _create_wx_menubar(self):
        """Implements _create_wx() for the root node.

        This method just calls _create_wx() on all its children.
        (The wxMenuBar object must already be set!)

        \see _create_wx()
        """
        # _wx_menubar must already be initialized with a wx.MenuBar object

        self._wx = self._wx_menubar
        
        # Create wx objects for all menus...
        for node in self:
            node._create_wx()


    def _create_wx_toplevel_menu(self):
        """Implements _create_wx() for toplevel menus.

        This method creates a wxMenu object and fills it with its
        children.

        \see _create_wx()
        """

        self._wx = wx.Menu()

        # Create wx objects for all items...
        for node in self:
            node._create_wx()

        # Parent is a wx.MenuBar object or None (popup menus)
        parent = self._wx_parent()
        if parent!=None:
            pos = self._wx_menu_pos()
            parent.Insert(pos, self._wx, self.text)


    def _create_wx_submenu(self):
        """Implements _create_wx() for submenus.

        This method creates a wxMenu object and fills it with its
        children. Then the wxMenu object is wrapped in a wxMenuItem
        object which serves as the actual menu item.

        \see _create_wx()
        """

        # Store the menu in _wx temporarily, so that the children have a
        # valid parent
        self._wx = wx.Menu()

        # Create wx objects for all items...
        for node in self:
            node._create_wx()

        # Parent is a wx.Menu object
        parent   = self._wx_parent()
        self._id = self.root._IDmanager.allocID()
        self._wx = wx.MenuItem(parent, self._id, self.text, "",
                               wx.ITEM_NORMAL, self._wx)
        pos = self._wx_menu_pos()
        parent.InsertItem(pos, self._wx)
       

    def _create_wx_embedded_menu(self):
        """Implements _create_wx() for embedded menus.

        This method does not create a wx menu object but only inserts
        its children into the parent menu.        

        \see _create_wx()
        """

        self._wx = None
        self._id = None
        wxmenu = self._wx_parent()
        # Create wx objects for all items...
        for node in self:
            node._create_wx()


    def _destroy_wx(self):

        # Destroy the wx objects for all items...
        for node in self:
            node._destroy_wx()

        # 1. Node = MenuBar
        if self._is_menubar():
            pass
        # 2. Node = Toplevel menu
        elif self.parent==None or self.parent._is_menubar():
            # Remove the menu from the menu bar (if it's not a popup menu)
            parent = self._wx_parent()
            if parent!=None:
                pos = self._wx_menu_pos()
                parent.Remove(pos)
        # 3. Node = Sub menu
        elif self.text!=None:
            # Remove the wxMenuItem from the parent menu
            parent = self._wx_parent()
            parent.Remove(self._id)
            # Make the ID available again
            self.root._IDmanager.freeID(self._id)
        # 3. Node = Embedded menu
#        else: pass

        MenuNode._destroy_wx(self)


    def _is_menubar(self):
        """Returns True if the node represents a menu bar."""
        return hasattr(self, "_wx_menubar")

    def _is_embedded_menu(self):
        """Returns True if the node represents an embedded menu."""
        return self.text==None and self.parent!=None and (not self._is_menubar())

    def _get_embedded_items(self):
        """
        May only be called on embedded menus.
        """

        res = []
        for cn in self._children:
            if cn._wx!=None:
                res.append(cn)
            elif isinstance(cn, Menu):
                res += cn._get_embedded_items()
        return res
        

    def _connect(self, node):
        """Connect a menu event with the handler.

        May only be called on the root node.
        """

        id = node._wx.GetId()
        wx.EVT_MENU(self._window, id, self._onMenuSelection)
        self._menu_id_lut[id] = node

    def _disconnect(self, node):
        """Disconnect a menu event from the handler.

        May only be called on the root node.
        """
        
        id = node._wx.GetId()
        print(("Disconnecting",node.name,"ID:",id))
        wx.EVT_MENU(self._window, id, None)
        del self._menu_id_lut[id]


    def _onMenuSelection(self, event):
        print(("Item",event.GetId(),"selected"))
        item = self._menu_id_lut.get(event.GetId(), None)
        if item==None:
            print ("Menu item ID not found!")
            return

        cmd = item.command
        if cmd!=None:
            cmd()
        

    def _dump(self, depth=0):
        print(("%s%s"%(2*depth*" ",self)))
        for n in self._children:
            if hasattr(n, "_dump"):
                n._dump(depth+1)
            else:
                print(("%s%s"%(2*(depth+1)*" ",n)))

