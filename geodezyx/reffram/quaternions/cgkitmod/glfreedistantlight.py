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

## \file glfreedistantlight.py
## Contains the GLFreeDistantLight class.

from . import protocols
from .Interfaces import *
from .slots import *
from .cgtypes import vec3
from .worldobject import _initWorldObject
from . import cmds
from . import _core

# GLFreeDistantLight
class GLFreeDistantLight(_core.GLDistantLight):
    """This class represents an OpenGL distant light.

    The direction of the light is the local positive Z axis.
    """

    protocols.advise(instancesProvide=[ISceneItem])

    def __init__(self,
                 name="GLDistantLight",
                 parent=None,
                 enabled=True,
                 intensity=1.0,
                 ambient=None,
                 diffuse=None,
                 specular=None,
                 cast_shadow = False,
                 auto_insert=True,
                 **params
                 ):
        _initWorldObject(self, baseClass=_core.GLDistantLight,
                         name=name, parent=parent,
                         auto_insert=auto_insert,
                         **params)

        self.enabled = enabled
        self.intensity = intensity
        if ambient!=None:
            self.ambient = vec3(ambient)
        if diffuse!=None:
            self.diffuse = vec3(diffuse)
        if specular!=None:
            self.specular = vec3(specular)

        self.cast_shadow = cast_shadow

    def protocols(self):
        return [ISceneItem]


        
