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
# $Id: objmaterial.py,v 1.7 2005/08/28 19:42:43 mbaas Exp $

## \file objmaterial.py
## Contains the OBJMaterial class.

import os.path, sys
from . import protocols
from .Interfaces import *
from .slots import *
from .material import Material
from . import ribexport
from . import _Image as Image
from .expression import Expression
from . import _core
from .cgtypes import *

# OBJTextureMap
class OBJTextureMap:
    """A texture map definition as it appears in OBJ/MTL files.
    """

    def __init__(self,
                 filename,
                 offset = (0,0,0),
                 scale = (1,1,1),
                 turb = (0,0,0),
                 mm = (0.0, 1.0),
                 clamp = False,
                 blendu = True,
                 blendv = True,
                 bumpsize = None,
                 refltype = None):
        
        self.filename = filename
        self.offset = offset
        self.scale = scale
        self.turb = turb
        self.mm = mm
        self.clamp = clamp
        self.blendu = blendu
        self.blendv = blendv
        if bumpsize!=None:
            self.bumpsize = bumpsize
        if refltype!=None:
            self.refltype = refltype
        
        
# OBJMaterial
#class OBJMaterial(Material):
class OBJMaterial(_core.GLMaterial):
    """This class represents a material as it appears in OBJ (or rather MTL) files.

    """

    protocols.advise(instancesProvide=[ribexport.IMaterial])
    
    def __init__(self,
                 name="OBJMaterial",
                 illum = 2,
                 Ka = (0.2, 0.2, 0.2),
                 Kd = (0.8, 0.8, 0.8),
                 Ks = (0.0, 0.0, 0.0),
                 Ke = (0.0, 0.0, 0.0),
                 Ns = 0.0,
                 Ni = 1.0,
                 d = 1.0,
                 Tr = 1.0,
                 Tf = (1.0, 1.0, 1.0),
                 sharpness = 0.0,
                 map_Ka = None,
                 map_Kd = None,
                 map_Ks = None,
                 map_Ke = None,
                 map_Ns = None,
                 map_d = None,
                 map_Bump = None,
                 refl = None,  # one single map or a list of maps
                 density = 1.0,
                 ):

#        Material.__init__(self, name=name, density=density)
        _core.GLMaterial.__init__(self, name, density)

        # Create a Kd slot and connect it with the GL diffuse slot
        self.Kd_slot = Vec3Slot(vec3(Kd))
        self.addSlot("Kd", self.Kd_slot)
        self.diffuse_adapter = Expression("(col.x,col.y,col.z,1)", col=vec3())
        self.Kd_slot.connect(self.diffuse_adapter.col_slot)
        self.diffuse_adapter.output_slot.connect(self.diffuse_slot)

        self.illum = illum
        self.Ka = Ka
#        self.Kd = Kd
        self.Ks = Ks
        self.Ke = Ke
        self.Ns = Ns
        self.Ni = Ni
        self.d = d
        self.Tr = Tr
        self.Tf = vec3(Tf)
        self.sharpness = sharpness
        self.map_Ka = map_Ka
        self.map_Kd = map_Kd
        self.map_Ks = map_Ks
        self.map_Ke = map_Ke
        self.map_Ns = map_Ns
        self.map_d = map_d
        self.map_Bump = map_Bump

        if refl==None:
            refl = []
        else:
            # Turn refl into a list
            try:
                refl = list(refl)
            except:
                refl = [refl]
        self.refl = refl

        # "Hide" the GLMaterial slots
        self.removeSlot("ambient")
        self.removeSlot("diffuse")
        self.removeSlot("specular")
        self.removeSlot("shininess")
        self.removeSlot("emission")

    exec(slotPropertyCode("Kd"))

    def mtlDefinition(self):
        s = """illum %d
Ka %f %f %f
Kd %f %f %f
Ks %f %f %f
Ke %f %f %f
Ns %f
Ni %f
"""%(self.illum,
     self.Ka[0], self.Ka[1], self.Ka[2],
     self.Kd[0], self.Kd[1], self.Kd[2],
     self.Ks[0], self.Ks[1], self.Ks[2],
     self.Ke[0], self.Ke[1], self.Ke[2],
     self.Ns,
     self.Ni)
        if self.d!=1.0:
            s += "d %f\n"%self.d
        if self.Tr!=1.0:
            s += "Tr %f\n"%self.d
        if self.Tf[0]!=1.0 or self.Tf[1]!=1.0 or self.Tf[2]!=1.0:
            s += "Tf %f %f %f\n"%(self.Tf[0], self.Tf[1], self.Tf[2])
        s += self.mapDefinition("map_Ka", self.map_Ka)
        s += self.mapDefinition("map_Kd", self.map_Kd)
        s += self.mapDefinition("map_Ks", self.map_Ks)
        s += self.mapDefinition("map_Ke", self.map_Ke)
        s += self.mapDefinition("map_Ns", self.map_Ns)
        s += self.mapDefinition("map_d", self.map_d)
        s += self.mapDefinition("map_Bump", self.map_Bump)
        for map in self.refl:
            s += self.mapDefinition("refl", map)
        return s

    # mapDefinition
    def mapDefinition(self, mapname, map):
        """Return the map defintion line or an empty string.
        """
        if map==None:
            return ""

        onoff = {True:"on", False:"off"}

        res = mapname
        # -o
        u,v,w = map.offset
        if u!=0 or v!=0 or w!=0:
            res += " -o %f %f %f"%(u,v,w)
        # -s
        u,v,w = map.scale
        if u!=1 or v!=1 or w!=1:
            res += " -s %f %f %f"%(u,v,w)
        # -t
        u,v,w = map.turb
        if u!=0 or v!=0 or w!=0:
            res += " -t %f %f %f"%(u,v,w)
        # -mm
        base, gain = map.mm
        if base!=0 or gain!=1:
            res += " -mm %f %f"%(base, gain)
        # -clamp
        if map.clamp:
            res += " -clamp %s"%onoff[map.clamp]
        # -blendu
        if not map.blendu:
            res += " -blendu %s"%onoff[map.blendu]
        # -blendv
        if not map.blendv:
            res += " -blendv %s"%onoff[map.blendv]
        # -bm
        if mapname=="map_Bump" and getattr(map, "bumpsize", 1)!=1:
            res += " -bm %f"%map.bumpsize
        # -type
        if mapname=="refl":
            res += " -type %s"%getattr(map, "refltype", "sphere")

        return res+" %s\n"%map.filename        
        
        
    def createPasses(self):
        """Returns a list of RenderPass objects."""
        texdefs = []
        if self.map_Ka!=None:
            texdefs.append((self.map_Ka.filename, "periodic", "periodic", "gaussian", 1.0, 1.0, {}))
        if self.map_Kd!=None:
            texdefs.append((self.map_Kd.filename, "periodic", "periodic", "gaussian", 1.0, 1.0, {}))
        if self.map_Ks!=None:
            texdefs.append((self.map_Ks.filename, "periodic", "periodic", "gaussian", 1.0, 1.0, {}))
        if self.map_d!=None:
            texdefs.append((self.map_d.filename, "periodic", "periodic", "gaussian", 1.0, 1.0, {}))
        if self.map_Bump!=None:
            texdefs.append((self.map_Bump.filename, "periodic", "periodic", "gaussian", 1.0, 1.0, {}))
        
        if texdefs!=[]:
            return [ribexport.TexPass(maps=texdefs)]
        else:
            return []

    def preProcess(self, exporter):
        """Preprocessing method.

        This method is called before the image is rendered and can be used
        to create or copy image maps.
        """
        pass
            

    def color(self):
        """Return the color for the RiColor() call or None.
        """
        return None
#        return self.Kd

    def opacity(self):
        """Return the opacity for the RiOpacity() call or None.
        """
        return None
#        a = self.Tr*self.d
#        return (a,a,a)

    def surfaceShaderName(self):
        """Returns the name of the corresponding surface shader or None.
        """
        return "matobj"

    def surfaceShaderSource(self):
        """Returns surface shader source code as a string or None.

        If the return value is None, then shaderName() must return
        the name of the shader to use.
        """
        return """// OBJ/MTL material shader
        
surface $SHADERNAME(float illum = 2;
           color Ka = color "rgb" (0.2, 0.2, 0.2);
           color Kd = color "rgb" (0.8, 0.8, 0.8);
           color Ks = color "rgb" (0.8, 0.8, 0.8);
           color Ke = color "rgb" (0, 0, 0);
           float Ni = 0.0;
           float Ns = 0.0;
           float d = 1.0;
           float Tr = 1.0;
           color Tf = color "rgb" (1.0, 1.0, 1.0);
           
           string map_Ka = "";
           float map_Ka_offset[3] = {0, 0, 0};
           float map_Ka_scale[3] = {1, 1, 1};
           
           string map_Kd = "";
           float map_Kd_offset[3] = {0, 0, 0};
           float map_Kd_scale[3] = {1, 1, 1};
           
           string map_Ks = "";
           float map_Ks_offset[3] = {0, 0, 0};
           float map_Ks_scale[3] = {1, 1, 1};

           string map_Ke = "";
           float map_Ke_offset[3] = {0, 0, 0};
           float map_Ke_scale[3] = {1, 1, 1};

           string map_Ns = "";
           float map_Ns_offset[3] = {0, 0, 0};
           float map_Ns_scale[3] = {1, 1, 1};

           string map_d = "";
           float map_d_offset[3] = {0, 0, 0};
           float map_d_scale[3] = {1, 1, 1};

           varying point Pref = point(0,0,0);
           )
{
  BAKE_BEGIN
  normal Nf = BAKE_NORMAL(N);
  color _Ka = Ka;
  color _Kd = Kd;
  color _Ks = Ks;
  color _Ke = Ke;
  float _Ns = Ns;
  float _d = d;

  if (map_Ka!="")
  {
    _Ka = texture(map_Ka,
                  s*map_Ka_scale[0]+map_Ka_offset[0],
                  (1-t)*map_Ka_scale[1]+map_Ka_offset[1]);
  }

  if (map_Kd!="")
  {
    _Kd = texture(map_Kd,
                  s*map_Kd_scale[0]+map_Kd_offset[0],
                  (1-t)*map_Kd_scale[1]+map_Kd_offset[1]);
  }

  if (map_Ks!="")
  {
    _Ks = texture(map_Ks,
                  s*map_Ks_scale[0]+map_Ks_offset[0],
                  (1-t)*map_Ks_scale[1]+map_Ks_offset[1]);
  }

  if (map_Ke!="")
  {
    _Ke = texture(map_Ke,
                  s*map_Ke_scale[0]+map_Ke_offset[0],
                  (1-t)*map_Ke_scale[1]+map_Ke_offset[1]);
  }

  if (map_Ns!="")
  {
    _Ns = texture(map_Ns,
                  s*map_Ns_scale[0]+map_Ns_offset[0],
                  (1-t)*map_Ns_scale[1]+map_Ns_offset[1]);
  }

  if (map_d!="")
  {
    _d = texture(map_d,
                 s*map_d_scale[0]+map_d_offset[0],
                 (1-t)*map_d_scale[1]+map_d_offset[1]);
  }

  if (illum==0)
  {
    Ci = _Kd;
  }
  else if (illum==1)
  {
    Ci = _Ka + _Kd*diffuse(Nf);
  }
  else
  {
    // Ns==0? then disable specular part
    if (_Ns<1E-6)
    {
      _Ks = color "rgb" (0,0,0);
    }
    Ci = _Ka*ambient() + _Kd*diffuse(Nf) + _Ks*specular(Nf,-normalize(I),1/_Ns) + _Ke;
  }  

  Oi = _d*Tr*Tf;
  Ci *= Oi;
  BAKE_END
}        
        """

    def surfaceShaderParams(self, passes):
        """Return a dictionary with shader parameters and their values."""
        res = {"uniform float illum" : self.illum,
               "uniform color Ka" : self.Ka,
               "uniform color Kd" : self.Kd,
               "uniform color Ks" : self.Ks,
               "uniform color Ke" : self.Ke,
               "uniform float Ns" : self.Ns,
               "uniform float Ni" : self.Ni,
               "uniform float d" : self.d,
               "uniform float Tr" : self.Tr,
               "uniform color Tf" : self.Tf
               }

        if self.map_Ka!=None:
            texname = os.path.basename(self.map_Ka.filename)
            name, ext = os.path.splitext(texname)
            res["uniform string map_Ka"] = name+".tex"       
            res["uniform float[3] map_Ka_offset"] = self.map_Ka.offset
            res["uniform float[3] map_Ka_scale"] = self.map_Ka.scale

        if self.map_Kd!=None:
            texname = os.path.basename(self.map_Kd.filename)
            name, ext = os.path.splitext(texname)
            res["uniform string map_Kd"] = name+".tex"       
            res["uniform float[3] map_Kd_offset"] = self.map_Kd.offset
            res["uniform float[3] map_Kd_scale"] = self.map_Kd.scale

        if self.map_Ks!=None:
            texname = os.path.basename(self.map_Ks.filename)
            name, ext = os.path.splitext(texname)
            res["uniform string map_Ks"] = name+".tex"       
            res["uniform float[3] map_Ks_offset"] = self.map_Ks.offset
            res["uniform float[3] map_Ks_scale"] = self.map_Ks.scale

        if self.map_Ke!=None:
            texname = os.path.basename(self.map_Ke.filename)
            name, ext = os.path.splitext(texname)
            res["uniform string map_Ke"] = name+".tex"  
            res["uniform float[3] map_Ke_offset"] = self.map_Ke.offset
            res["uniform float[3] map_Ke_scale"] = self.map_Ke.scale

        if self.map_Ns!=None:
            texname = os.path.basename(self.map_Ns.filename)
            name, ext = os.path.splitext(texname)
            res["uniform string map_Ns"] = name+".tex"
            res["uniform float[3] map_Ns_offset"] = self.map_Ns.offset
            res["uniform float[3] map_Ns_scale"] = self.map_Ns.scale

        if self.map_d!=None:
            texname = os.path.basename(self.map_d.filename)
            name, ext = os.path.splitext(texname)
            res["uniform string map_d"] = name+".tex"       
            res["uniform float[3] map_d_offset"] = self.map_d.offset
            res["uniform float[3] map_d_scale"] = self.map_d.scale

        return res
        
    def surfaceShaderTransform(self):
        return mat4(1)

    def displacementShaderName(self):
        if self.map_Bump==None:
            return None
        
        return "matobj_bump"

   
    def displacementShaderSource(self):
        if self.map_Bump==None:
            return None
        
        return """// OBJ/MTL material shader (bump)
        
displacement $SHADERNAME(
           string map_Bump = "";
           float bump_size = 1.0;
           float map_Bump_offset[3] = {0, 0, 0};
           float map_Bump_scale[3] = {1, 1, 1};
           )
{
  float amount = 0;
  
  if (map_Bump!="")
  {
    color bc = texture(map_Bump,
                       s*map_Bump_scale[0]+map_Bump_offset[0],
                      (1-t)*map_Bump_scale[1]+map_Bump_offset[1]);
    amount = (comp(bc,0)+comp(bc,1)+comp(bc,2))/3;
  }

  P = P + bump_size*amount*normalize(N);
  N = calculatenormal(P);
}
        """
    
    def displacementShaderParams(self, passes):
        if self.map_Bump==None:
            return {}

        res = {}
        texname = os.path.basename(self.map_Bump.filename)
        name, ext = os.path.splitext(texname)
        res["uniform string map_Bump"] = name+".tex"
        res["uniform float bump_size"] = getattr(self.map_Bump, "bumpsize", 1.0)
        res["uniform float[3] map_Bump_offset"] = self.map_Bump.offset
        res["uniform float[3] map_Bump_scale"] = self.map_Bump.scale
        return res

    def displacementShaderTransform(self):
        return mat4(1)

    def displacementBound(self):
        if self.map_Bump==None:
            return "current", 0
        else:
            return "current", abs(getattr(self.map_Bump, "bumpsize", 1.0))
	
    def interiorShaderName(self):
        return None
    
    def interiorShaderSource(self):
        return None
    
    def interiorShaderParams(self, passes):
        return {}

    def interiorShaderTransform(self):
        return mat4(1)
      
        
