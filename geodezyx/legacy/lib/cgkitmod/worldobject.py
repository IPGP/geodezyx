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

## \file worldobject.py
## Contains the WorldObject class.

#from _core import WorldObject as _WorldObject
#from _core import _WorldObjectChildIterator
from . import _core
from .Interfaces import ISceneItem, ISceneItemContainer
from . import protocols
from . import globalscene
from .cgtypes import *

# WorldObject
class WorldObject(_core.WorldObject):
    """The base class for a world object.

    A world object is a scene item that has a position and orientation in
    space and that can usually be represented by 3D geometry.

    Attributes:

    - name
    - transform
    - pos
    - rot
    - scale
    - pivot
    - mass

    If you access an attribute that's actually stored in the geom then
    the attribute access is forwarded to the geom. This means you can
    access the geom slots through the world object (for example, if you
    attach a sphere geom to a world object, then the world object "has"
    an attribute 'radius').
    """

    protocols.advise(instancesProvide=[ISceneItem, ISceneItemContainer])
    
    def __init__(self,
                 name="object",
                 transform = None,
                 pos = None, rot = None, scale = None,
                 pivot = None,
                 offsetTransform = None,
                 parent = None,
                 mass = None,
                 material = None,
                 visible = True,
                 linearvel = None,
                 angularvel = None,
                 auto_insert = True):
        """Constructor.

        \param name (\c str) Object name
        \param transform (\c mat4) Initial transform
        \param pos (\c vec3) Initial position
        \param rot (\c mat3) Initial rotation
        \param scale (\c vec3) Initial scaling
        \param pivot (\c vec3) Initial pivot point (takes precedence over offsetTransform)
        \param offsetTransform (\c mat4) Initial offset transform
        \param parent (\c WorldObject or \c str) Parent object or None
        \param mass (\c float) Total mass
        \param material (\c Material) Material class (or a sequence of materials)
        \param visible (\c Bool) Visibility flag
        \param linearvel (\c vec3) Linear velocity
        \param angularvel (\c vec3) Angular velocity
        \param auto_insert (\c bool) If True, the object is inserted into the
                        scene automatically
        """

        _initWorldObject(self, baseClass=_core.WorldObject,
                         name=name,
                         transform=transform, pos=pos, rot=rot,
                         scale=scale, pivot=pivot,
                         offsetTransform=offsetTransform,
                         parent=parent,
                         mass=mass, material=material,
                         visible=visible,
                         linearvel=linearvel, angularvel=angularvel,
                         auto_insert=auto_insert)
        
#        if offsetTransform!=None:
#            self.setOffsetTransform(offsetTransform)
#        if pivot!=None:
#            self.pivot = pivot
#        if transform!=None:
#            self.transform = transform
#        if pos!=None:
#            self.pos = pos
#        if rot!=None:
#            self.rot = rot
#        if scale!=None:
#            self.scale = scale
#        if mass!=None:
#            self.mass = mass
#        if material!=None:
#            self.material = material

#        if auto_insert:
#            scene.getScene().insert(self)

    def protocols(self):
        return [ISceneItem, ISceneItemContainer]

    def __getattr__(self, name):
        if self.geom!=None and name[:2]!="__":
            if hasattr(self.geom, name):
                return getattr(self.geom, name)
            elif self.geom.hasSlot(name):
                return self.geom.slot(name).getValue()
#        if self.geom!=None and self.geom.hasSlot(name):
#            exec "res=self.geom.%s"%name
#            return res
        raise AttributeError('Object "%s" has no attribute "%s"'%(self.name, name))
    
    def __setattr__(self, name, val):
        if self.geom!=None and self.geom.hasSlot(name):
#            print "self.geom.%s=%s"%(name, val)
            exec("self.geom.%s=%s"%(name, val))
        else:
            _core.WorldObject.__setattr__(self, name, val)
        

# Common WorldObject initializations
def _initWorldObject(self, baseClass,
                     name, parent, transform=None,
                     pos=None, rot=None, scale=None,
                     pivot=None, offsetTransform=None,
                     mass=None, material=None, visible=True,
                     linearvel=None, angularvel=None,
                     auto_insert=True):
    """Helper function for usage in constructors.
    """
    if auto_insert:
        if parent is None:
            parent = globalscene.getScene().worldRoot()
        else:
            if type(parent) is str:
                parent = globalscene.getScene().worldObject(parent)
        name = parent.makeChildNameUnique(name)
    
    baseClass.__init__(self, name)
    
    if offsetTransform!=None:
        self.setOffsetTransform(offsetTransform)
    if pivot!=None:
        self.pivot = vec3(pivot)
    if transform!=None:
        self.transform = transform
    if pos!=None:
        self.pos = vec3(pos)
    if rot!=None:
        self.rot = mat3(rot)
    if scale!=None:
        self.scale = vec3(scale)
    if mass!=None:
        self.mass = mass
    if material!=None:
        try:
            # Check if material is a sequence or not. If it is not a
            # sequence the following line will raise an exception.
            len(material)
        except:
            material = [material]
        self.setNumMaterials(len(material))
        for i,mat in enumerate(material):
            self.setMaterial(mat, i)
    if linearvel!=None:
        self.linearvel = vec3(linearvel)
    if angularvel!=None:
        self.angularvel = vec3(angularvel)

    self.visible = visible

    if auto_insert:
        parent.addChild(self)
#        scene.getScene().insert(self)

