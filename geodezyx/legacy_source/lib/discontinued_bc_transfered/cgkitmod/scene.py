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
# $Id: scene.py,v 1.2 2005/04/14 07:58:55 mbaas Exp $

## \file scene.py
## Contains the Scene class.

from .cgtypes import vec3
from .Interfaces import ISceneItem
from . import timer
from . import protocols
from .boundingbox import BoundingBox
from . import joystick
from . import _core
from . import globalscene

# Scene
class Scene(object):
    """%Scene class."""
    
    def __init__(self):
        object.__init__(self)
        
        # Handedness ('r' or 'l')
        self._handedness = "r"
        # Normalized up direction
        self._up = vec3(0,0,1)
        # Unit   (1 unit = unitscale [Unit])
        self._unit = "m"
        # Unit scale
        self._unitscale = 1.0

        # Global options
        self._globals = {}
        
        self.items = []
#        self.items_by_name = {}

        self._worldroot = _core.WorldObject("WorldRoot")
        self._timer = timer.Timer(auto_insert=False)

        # Key: ID  Value:Joystick object
        self._joysticks = {}
        
        self.clear()

    def __len__(self):
        return len(self.items)

    def __iter__(self):
        return iter(self.items)

    def hasGlobal(self, name):
        return name in self._globals
            
    def getGlobal(self, name, default=None):
        return self._globals.get(name, default)

    def setGlobal(self, name, value):
        self._globals[name] = value

    # createDefaultItems
#    def createDefaultItems(self):
#        self._worldroot = WorldObject("root", auto_insert=False)

    # worldRoot
    def worldRoot(self):
        """Return the root object of the world.
        """
        return self._worldroot

    # walkWorld
    def walkWorld(self, root=None):
        """Walk the world tree and yield each object.

        This method can be used to iterator over the entire world tree or
        a subtree thereof. The argument root specifies the root of the
        tree which is to traverse (the root itself will not be returned).
        """
        if root==None:
            root = self._worldroot
        return self._walkWorld(root)

    # timer
    def timer(self):
        """Return the global timer object."""
        return self._timer

    # clear
    def clear(self):
        """Clear the entire scene."""

        # Remove all children of the world root...
        for obj in list(self._worldroot.iterChilds()):
            self._worldroot.removeChild(obj)
            
        self.items = [self._timer, self._worldroot]
#        self.items_by_name = {}

    # insert
    def insert(self, item):
        """Insert an item into the scene."""
#        protocols.adapt(item, ISceneItem)
#        if isinstance(item, worldobject.WorldObject):
        if isinstance(item, _core.WorldObject):
            self._worldroot.addChild(item)
        else:
            self.items.append(item)

    # item
    def item(self, name):
        """Return the item with the specified name.

        \param name (\c str) Item name
        \return Item or None if no item was found
        """
        for item in self.items:
            if name==item.name:
                return item
        return None

    # worldObject
    def worldObject(self, name):
        """Return the world object with the specified name.

        You can use '|' as a path separator.

        \param name (\c str) Object name
        \return WorldObject
        """
        obj = self._worldroot
        if len(name)>0:
            names = name.split("|")
        else:
            names = []
        for n in names:
            obj = obj.child(n)
        return obj

    # boundingBox
    def boundingBox(self):
        """Return the bounding box of the entire scene.

        \todo Use the world transform instead of the local transform
        """
        return self._worldroot.boundingBox()
#        return self._boundingBox(self._worldroot)
            
#    def _boundingBox(self, node):
#        res = BoundingBox()
#        for obj in node.iterChilds():
#            L = obj.localTransform()
#            g = obj.geom
#            if g!=None:
#                bb = g.boundingBox()
#                if not bb.isEmpty():
#                    res.addBoundingBox(bb.transform(L))

#            res.addBoundingBox(self._boundingBox(obj).transform(L))

#        return res

    # setJoystick
    def setJoystick(self, joystick):
        """Set a joystick object.
        """
        self._joysticks[joystick.id] = joystick

    # getJoystick
    def getJoystick(self, id):
        """Get a joystick object.

        A dummy joystick object is returned if there is no joystick
        with the specified id.
        """
        if id in self._joysticks:
            return self._joysticks[id]
        else:
            return joystick.Joystick(id=id, name="Dummy-Joystick")

    ## protected:


    def _walkWorld(self, obj):
        """Helper method for the walkWorld() method."""
        for child in obj.iterChilds():
            yield child
            for c in self._walkWorld(child):
                yield c

    
    # "handedness" property...
    
    def _getHandedness(self):
        """Return the handedness.

        This method is used for retrieving the \a handedness property.

        \return Handedness (\c str) which is either 'l' or 'r'.
        """
        return self._handedness

    def _setHandedness(self, h):
        """Set the handedness.

        This method is used for setting the \a handedness property.

        \param h (\c str) Handedness ('l' or 'r')
        """
        h = h.lower()
        if h!="l" and h!="r":
            raise ValueError("Handedness must be either 'l' or 'r'.")
        self._handedness = h

    handedness = property(_getHandedness, _setHandedness, None, "Handedness")

    # "up" property...
    
    def _getUp(self):
        """Return the up vector.

        This method is used for retrieving the \a up property.

        \return Up vector (\c vec3)
        """
        return self._up

    def _setUp(self, up):
        """Set the up vector.

        This method is used for setting the \a up property.

        \param up (\c vec3) Up vector
        """
        self._up = vec3(up).normalize()

    up = property(_getUp, _setUp, None, "Up vector")
    

######################################################################

# getScene
def getScene():
    """Return the global %scene instance."""
    return globalscene.getScene()

