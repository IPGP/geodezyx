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
# $Id$

import os.path
from . import protocols
from .Interfaces import *
from .slots import *
from . import lookat
from .lightsource import LightSource
from . import ribexport
from math import *

# MayaSpotLight
class MayaSpotLight(LightSource):
    """This class represents a spotlight that approximates a Maya spot light.

    The direction of the light is the local negative Z axis.
    """

    protocols.advise(instancesProvide=[ISceneItem, ribexport.ILightSource])

    def __init__(self,
                 name="MayaSpotLight",
                 enabled = True,
                 color = (1,1,1),
                 intensity = 1.0,
                 decayRate = 0,
                 coneAngle = 40,
                 penumbraAngle = 0,
                 dropoff = 0,

                 useDepthMapShadows = False,
                 dmapResolution = 512,
                 useMidDistDmap = True,
                 useDmapAutoFocus = True,
                 dmapFocus = 90.0,
                 dmapFilterSize = 1,
                 dmapBias = 0.001,
                 
                 target=vec3(0,0,0),  # spot
                 **params
                 ):

        LightSource.__init__(self, name=name, **params)

#        target = vec3(target)

        # Target
#        self.target_slot = Vec3Slot(target)
#        self.addSlot("target", self.target_slot)

        self.enabled = enabled
        self.color = vec3(color)
        self.intensity = float(intensity)
        self.decayRate = int(decayRate)
        self.coneAngle = float(coneAngle)
        self.penumbraAngle = float(penumbraAngle)
        self.dropoff = float(dropoff)

        self.useDepthMapShadows = bool(useDepthMapShadows)
        self.dmapResolution = int(dmapResolution)
        self.useMidDistDmap = bool(useMidDistDmap)
        self.useDmapAutoFocus = bool(useDmapAutoFocus)
        self.dmapFocus = float(dmapFocus)
        self.dmapFilterSize = float(dmapFilterSize)
        self.dmapBias = float(dmapBias)
        
    def protocols(self):
        return [ISceneItem]


    def createPasses(self):
        """Returns a list of RenderPass objects."""
       
        if self.useDepthMapShadows:
            if self.useDmapAutoFocus:
                fov = max(self.coneAngle, self.coneAngle+2*self.penumbraAngle)
            else:
                fov = self.dmapFocus
            shadpass = ribexport.ShadowPass(
                                  output = [("shadow.z", "zfile", "z", {})],
                                  light = self,
                                  fov = fov,
                                  resolution = self.dmapResolution,
                                  orientoffset = mat4(1).rotation(pi, vec3(0,1,0))
)
            return [shadpass]

        return []

    def shaderName(self):
        """Returns the name of the corresponding light shader or None.
        """
        return "mayaspotlight"

    def shaderSource(self):
        """Returns surface shader source code as a string or None.
        """
        return """// MayaSpotLight shader
        
light $SHADERNAME(
         uniform color lightcolor = color "rgb" (1, 1, 1);
         uniform float intensity = 1.0;
         uniform float decayRate = 0;
         uniform float coneAngle = 40.0;
         uniform float penumbraAngle = 0.0;
         uniform float dropoff = 0.0;

         uniform string dmapName = "";
         uniform float dmapFilterSize = 1;
         uniform float dmapBias = 0.001;
         )
{
  uniform vector axis = normalize(vector "shader" (0,0,-1));
  uniform float coneAngle2 = radians(0.5*coneAngle);
  uniform float penumbra_start = coneAngle2;
  uniform float penumbra_end = coneAngle2;
  if (penumbraAngle<0)
    penumbra_start += radians(penumbraAngle);
  else
    penumbra_end += radians(penumbraAngle);
  uniform float penumbra_range = penumbra_end-penumbra_start;

  illuminate(point "shader" (0,0,0), axis, penumbra_end)
  {
    vector L0 = normalize(L);
    float dist = length(vtransform("world", L));
    float att = 1.0;

    // Penumbra
    float angle = acos(L0.axis);
    // Linear decay inside the penumbra
//    att *= clamp(1.0-(angle-penumbra_start)/penumbra_range, 0.0, 1.0);
    // Smooth decay inside the penumbra
    att *= 1.0-smoothstep(penumbra_start, penumbra_end, angle);

    // Dropoff
//    att *= clamp(2.0-pow(exp(angle),dropoff/255), 0.0, 1.0);

    // Decay
    if (dist>1.0)
      att /= pow(dist, decayRate);

    Cl = att * intensity * lightcolor;

    if (dmapName!="")
      Cl *= 1.0 - shadow(dmapName, Ps, "width", dmapFilterSize, "bias", dmapBias, "samples", 32);
    
  }
}        
        """

    def shaderParams(self, passes):
        """Return a dictionary with shader parameters and their values."""
        params = {"intensity":self.intensity,
                  "uniform color lightcolor":self.color,
                  "uniform float decayRate":self.decayRate,
                  "uniform float coneAngle":self.coneAngle,
                  "uniform float penumbraAngle":self.penumbraAngle,
                  "uniform float dropoff":self.dropoff,
                  }
        if self.useDepthMapShadows and passes[0].done():
            zfile = passes[0].realFilename(passes[0].output[0][0])
            mapname = os.path.splitext(zfile)[0]+".map"
            params["uniform string dmapName"] = mapname
            params["uniform float dmapFilterSize"] = self.dmapFilterSize
            params["uniform float dmapBias"] = self.dmapBias
        return params
    

        
