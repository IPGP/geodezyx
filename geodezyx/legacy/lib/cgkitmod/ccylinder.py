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
# $Id: ccylinder.py,v 1.2 2005/04/19 08:45:56 mbaas Exp $

## \file ccylinder.py
## Contains the capped cylinder class.

from .cgtypes import vec3
from .Interfaces import *
from .worldobject import WorldObject
from .ccylindergeom import CCylinderGeom
from .slots import *
from . import protocols


# CCylinder
class CCylinder(WorldObject):

    protocols.advise(instancesProvide=[ISceneItem, IRigidBody])

    def __init__(self,
                 name="CCylinder",
                 dynamics=True,
                 static=False,
                 radius=1.0, length=1.0,
                 segmentsu = 16,
                 segmentsvl = 1,
                 segmentsvr = 4,
                 **params):
        WorldObject.__init__(self, name=name, **params)

        self.geom = CCylinderGeom(radius, length, segmentsu, segmentsvl, segmentsvr)

        self.dynamics_slot = BoolSlot(dynamics)
        self.static_slot = BoolSlot(static)
        self.addSlot("dynamics", self.dynamics_slot)
        self.addSlot("static", self.static_slot)

    exec(slotPropertyCode("dynamics"))
    exec(slotPropertyCode("static"))
        
    def protocols(self):
        return [ISceneItem, IRigidBody]
