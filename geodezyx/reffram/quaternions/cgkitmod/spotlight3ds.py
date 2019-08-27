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
# $Id: spotlight3ds.py,v 1.2 2005/07/03 09:43:36 mbaas Exp $

## \file spotlight3ds.py
## Contains the SpotLight3DS class.

from . import protocols
from .Interfaces import *
from .slots import *
from . import lookat
from .lightsource import LightSource

# SpotLight3DS
class SpotLight3DS(LightSource):
    """This class represents a spotlight as it appears in 3DS files.

    The direction of the light is the local positive Z axis.
    """

    protocols.advise(instancesProvide=[ISceneItem])

    def __init__(self,
                 name="SpotLight3DS",
                 enabled = True,   # off
                 intensity = 1.0,  # multiplier

                 see_cone = False,
                 roll = 0.0, 
                 outer_range = 0,
                 inner_range = 0,
                 attenuation = 0,
                 rectangular_spot = 0,
                 shadowed = False,
                 shadow_bias = 0,
                 shadow_filter = 4.0,
                 shadow_size = 256,
                 spot_aspect = 0,
                 use_projector = False,
                 projector = 0,
                 overshoot = False,
                 ray_shadows = False,
                 ray_bias = False,
                 hotspot = 43,
                 falloff = 45,
                 color = (1,1,1),
                 target=vec3(0,0,0),  # spot
#                 transform = None,
#                 pos = None, rot = None, scale = None,
#                 pivot= None,
#                 offsetTransform = None,
#                 auto_insert=True,
                 **params
                 ):

        LightSource.__init__(self, name=name, **params)

        target = vec3(target)

        # Target
        self.target_slot = Vec3Slot(target)
        self.addSlot("target", self.target_slot)

        self.enabled = enabled
        self.intensity = intensity
        self.attenuation = attenuation
        self.attenuation = attenuation
        self.inner_range = inner_range
        self.outer_range = outer_range
        self.falloff = falloff
        self.hotspot = hotspot
        self.overshoot = overshoot
        self.color = color
        self.shadowed = shadowed
        self.cast_shadow = shadowed
        self.shadow_filter = shadow_filter
        self.shadow_bias = shadow_bias
        self.shadow_size = shadow_size

        # Create the internal LookAt component
        self._lookat = lookat.LookAt()
        self._lookat.name = "SpotLight3DS_LookAt"
        self._lookat.output_slot.connect(self.rot_slot)
        self.pos_slot.connect(self._lookat.pos_slot)
        self.target_slot.connect(self._lookat.target_slot)
        

    # Create the "target" property
    exec(slotPropertyCode("target"))

    def protocols(self):
        return [ISceneItem]

    

        
