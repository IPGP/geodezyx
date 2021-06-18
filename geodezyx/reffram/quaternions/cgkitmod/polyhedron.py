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
# $Id: polyhedron.py,v 1.1 2005/02/20 18:30:41 mbaas Exp $

## \file polyhedron.py
## Contains the Polyhedron class.

from .cgtypes import *
from .Interfaces import *
from .worldobject import WorldObject
from .polyhedrongeom import PolyhedronGeom
from .slots import *
from . import protocols
from . import _core


# Polyhedron
class Polyhedron(WorldObject):

    protocols.advise(instancesProvide=[ISceneItem])

    def __init__(self,
                 name = "Polyhedron",
#                 dynamics = False,
#                 static = False,
                 verts = [],
                 polys = [],
                 **params):
        WorldObject.__init__(self, name=name, **params)

        self.geom = PolyhedronGeom()

#        self.dynamics = dynamics
#        self.static_slot = BoolSlot(static)

        ph = self.geom
        
        if len(verts)>0:
            ph.verts.resize(len(verts))
            i = 0
            for v in verts:
                ph.verts.setValue(i, vec3(v))
                i+=1

        if len(polys)>0:
            ph.setNumPolys(len(polys))
            i = 0
            for poly in polys:
                ph.setPoly(i, poly)
                i+=1
        
    exec(slotPropertyCode("static"))

