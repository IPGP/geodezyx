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
# $Id: material3ds.py,v 1.4 2005/08/28 19:42:43 mbaas Exp $

## \file material3ds.py
## Contains the Material3DS class.

import os.path, sys, shutil
import math
from . import protocols
from ._OpenGL.GL import *
from .Interfaces import *
from .slots import *
from . import lookat
from . import sl
from .material import Material
from .glmaterial import GLMaterial, GLTexture
from . import ribexport
from . import _Image as Image

# TextureMap3DS
class TextureMap3DS:
    """A texture map definition as it appears in 3DS files.
    """
    
    def __init__(self,
                 name,
                 flags = 0,
                 percent = 0.0,
                 blur = 0.0,
                 scale = (1.0, 1.0),
                 offset = (0,0),
                 rotation = 0,
                 tint1 = (0,0,0),
                 tint2 = (0,0,0),
                 tintr = (0,0,0),
                 tintg = (0,0,0),
                 tintb = (0,0,0)
                 ):
        self.name = name
        self.flags = flags
        self.percent = percent
        self.blur = blur
        self.scale = scale
        self.offset = offset
        self.rotation = rotation
        self.tint1 = tint1
        self.tint2 = tint2
        self.tintr = tintr
        self.tintg = tintg
        self.tintb = tintb

    def __str__(self):
        s = "TextureMap3DS:\n"
        s += '  Name    : "%s"\n'%self.name
        s += '  Flags   : %d\n'%self.flags
        s += '  Percent : %s\n'%self.percent
        s += '  Blur    : %s\n'%self.blur
        s += '  Scale   : %s, %s\n'%self.scale
        s += '  Offset  : %s, %s\n'%self.offset
        s += '  Rotation: %s\n'%self.rotation
        s += '  Tint1   : %s\n'%(self.tint1,)
        s += '  Tint2   : %s\n'%(self.tint2,)
        s += '  TintR   : %s\n'%(self.tintr,)
        s += '  TintG   : %s\n'%(self.tintg,)
        s += '  TintB   : %s'%(self.tintb,)
        return s


# Material3DS
class Material3DS(Material):
    """This class represents a material as it appears in 3DS files.

    """

    protocols.advise(instancesProvide=[ribexport.IMaterial])
    
#    protocols.advise(instancesProvide=[ISceneItem])

    def __init__(self,
                 name="Material3DS",
                 ambient = (0,0,0,0),
                 diffuse = (1.0, 1.0, 1.0, 1.0),
                 specular = (1.0, 1.0, 1.0, 1.0),
                 shininess = 1,
                 shin_strength = 0,
                 use_blur = 0,
                 transparency = 0.0,
                 falloff = 0,
                 additive = 0,
                 use_falloff = 0,
                 self_illum = False,
                 self_ilpct = 0.0,
                 shading = 0,
                 soften = 0,
                 face_map = 0,
                 two_sided = 0,
                 map_decal = 0,
                 use_wire = 0,
                 use_wire_abs = 0,
                 wire_size = 0,
                 density = 1.0,
                 texture1_map = None,
                 texture1_mask = None,
                 texture2_map = None,
                 texture2_mask = None,
                 opacity_map = None,
                 opacity_mask = None,
                 bump_map = None,
                 bump_mask = None,
                 specular_map = None,
                 specular_mask = None,
                 shininess_map = None,
                 shininess_mask = None,
                 self_illum_map = None,
                 self_illum_mask = None,
                 reflection_map = None,
                 reflection_mask = None,

                 bump_size = 1.0   # Extra parameter to control the bump map
                 ):

        Material.__init__(self, name=name, density=density)

        self.ambient = ambient
        self.diffuse = diffuse
        self.specular = specular
        self.shininess = shininess
        self.shin_strength = shin_strength
        self.transparency = transparency
        self.self_illum = self_illum
        self.self_ilpct = self_ilpct
        self.texture1_map = texture1_map
        self.texture1_mask = texture1_mask
        self.texture2_map = texture2_map
        self.texture2_mask = texture2_mask
        self.opacity_map = opacity_map
        self.opacity_mask = opacity_mask
        self.bump_map = bump_map
        self.bump_mask = bump_mask
        self.specular_map = specular_map
        self.specular_mask = specular_mask
        self.shininess_map = shininess_map
        self.shininess_mask = shininess_mask
        self.self_illum_map = self_illum_map
        self.self_illum_mask = self_illum_mask
        self.reflection_map = reflection_map
        self.reflection_mask = reflection_mask
        self.bump_size = bump_size

        if texture1_map==None:
            self._gltexture = None
            glambient = ambient
            gldiffuse = diffuse
        else:
            map = texture1_map
            T = mat4(1,0,0,-map.offset[0]-0.5, 0,-1,0,0.5-map.offset[1], 0,0,1,0, 0,0,0,1)
            a = sl.radians(map.rotation)
            ca = math.cos(a)
            sa = math.sin(a)
            R = mat4(ca,-sa,0,0, sa,ca,0,0, 0,0,1,0, 0,0,0,1)
            S = mat4(map.scale[0],0,0,0.5, 0,map.scale[1],0,0.5, 0,0,1,0, 0,0,0,1)
            self._gltexture = GLTexture(imagename = map.name,
                                        mode = GL_MODULATE,
#                                        transform = mat4().scaling(vec3(1,-1,1))
                                        transform = S*R*T
                                        )
            glambient = map.percent*vec4(1,1,1,1) + (1-map.percent)*vec4(ambient)
            gldiffuse = map.percent*vec4(1,1,1,1) + (1-map.percent)*vec4(diffuse)

        self._glmaterial = GLMaterial(ambient = glambient,
                                      diffuse = gldiffuse,
                                      specular = shin_strength*specular,
                                      shininess = 25*shininess,
                                      texture = self._gltexture
                                      )

##        print >>sys.stderr, "---",name,"---"
##        print >>sys.stderr, "ambient",self.ambient
##        print >>sys.stderr, "diffuse",self.diffuse
##        print >>sys.stderr, "Map1:"
##        print >>sys.stderr, self.texture1_map
##        print >>sys.stderr, "Map2:"
##        print >>sys.stderr, self.texture2_map
##        print >>sys.stderr, "Opacity:"
##        print >>sys.stderr, self.opacity_map
##        print >>sys.stderr, "Bump:"
##        print >>sys.stderr, self.bump_map
##        print >>sys.stderr, "Specular:"
##        print >>sys.stderr, self.specular_map
##        print >>sys.stderr, "Shininess:"
##        print >>sys.stderr, self.shininess_map
##        print >>sys.stderr, "Self illum:"
##        print >>sys.stderr, self.self_illum_map
##        print >>sys.stderr, "Reflection:"
##        print >>sys.stderr, self.reflection_map
##        print >>sys.stderr, "shininess",shininess
##        print >>sys.stderr, "shin_strength",shin_strength
##        print >>sys.stderr, "transparency",transparency
##        print >>sys.stderr, "falloff",falloff
##        print >>sys.stderr, "use_falloff",use_falloff
##        print >>sys.stderr, "use_blur", use_blur
##        print >>sys.stderr, "additive",additive
##        print >>sys.stderr, "shading",shading
##        print >>sys.stderr, "soften",soften
##        print >>sys.stderr, "self_illum", self_illum
##        print >>sys.stderr, "face_map", face_map
##        print >>sys.stderr, "two_sided", two_sided
##        print >>sys.stderr, "map_decal", map_decal
##        print >>sys.stderr, "use_wire", use_wire
##        print >>sys.stderr, "use_wire_abs", use_wire_abs
##        print >>sys.stderr, "wire_size", wire_size

    def applyGL(self):
        self._glmaterial.applyGL()

    def usesBlending(self):
        return self._glmaterial.usesBlending()

        
    def createPasses(self):
        """Returns a list of RenderPass objects."""
        texdefs = []
        texdefs += self._createTexDef(self.texture1_map)
        texdefs += self._createTexDef(self.opacity_map)
        texdefs += self._createTexDef(self.bump_map)
        texdefs += self._createTexDef(self.specular_map)
        texdefs += self._createTexDef(self.shininess_map)
        texdefs += self._createTexDef(self.self_illum_map)
        texdefs += self._createTexDef(self.reflection_map)

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
#        self._copyImageMap(self.texture1_map, exporter)                 
#        self._copyImageMap(self.opacity_map, exporter)
#        self._copyImageMap(self.bump_map, exporter)
#        self._copyImageMap(self.specular_map, exporter)
#        self._copyImageMap(self.shininess_map, exporter)
#        self._copyImageMap(self.self_illum_map, exporter)
#        self._copyImageMap(self.reflection_map, exporter)
            

    def color(self):
        """Return the color for the RiColor() call or None.
        """
        return self.diffuse

    def opacity(self):
        """Return the opacity for the RiOpacity() call or None.
        """
        a = 1.0-self.transparency
        return (a,a,a)

    def surfaceShaderName(self):
        """Returns the name of the corresponding surface shader or None.
        """
        return "mat3ds"

    def surfaceShaderSource(self):
        """Returns surface shader source code as a string or None.

        If the return value is None, then shaderName() must return
        the name of the shader to use.
        """
        return """// 3DS material shader
        
surface $SHADERNAME(color ambient_col = color "rgb" (0, 0, 0);
           color diffuse_col = color "rgb" (1, 1, 1);
           color specular_col = color "rgb" (1, 1, 1);
           float shininess = 0.0;
           float shin_strength = 0.0;
           float self_ilpct = 0.0;
           // Texture map 1
           string texture1_map = "";
           float t1_uscale = 1.0;
           float t1_vscale = 1.0;
           float t1_uoffset = 0;
           float t1_voffset = 0;
           float t1_rotation = 0;
           float t1_blur = 0.0;
           // Specular map
           string specular_map = "";
           float sp_uscale = 1.0;
           float sp_vscale = 1.0;
           float sp_uoffset = 0;
           float sp_voffset = 0;
           float sp_rotation = 0;
           float sp_blur = 0.0;
           // Shininess map
           string shininess_map = "";
           float sh_uscale = 1.0;
           float sh_vscale = 1.0;
           float sh_uoffset = 0;
           float sh_voffset = 0;
           float sh_rotation = 0;
           float sh_blur = 0.0;
           // Opacity map
           string opacity_map = "";
           float op_uscale = 1.0;
           float op_vscale = 1.0;
           float op_uoffset = 0;
           float op_voffset = 0;
           float op_rotation = 0;
           float op_blur = 0.0;
           // Self illumination map
           string self_illum_map = "";
           float si_uscale = 1.0;
           float si_vscale = 1.0;
           float si_uoffset = 0;
           float si_voffset = 0;
           float si_rotation = 0;
           float si_blur = 0.0;
           // Reflection map
           string reflection_map = "";
           float refl_uscale = 1.0;
           float refl_vscale = 1.0;
           float refl_uoffset = 0;
           float refl_voffset = 0;
           float refl_rotation = 0;
           float refl_percent = 0.0;
           float refl_blur = 0.0;
           varying point Pref = point(0,0,0);
           )
{
  BAKE_BEGIN
  normal Nf = BAKE_NORMAL(N);
  color C_diffuse = diffuse_col;
  color C_specular = specular_col;
  color C_opacity = Os;
  color C_refl = 0;
  float final_shininess = shininess;
  float final_self_ilpct = self_ilpct;
  float s0, t0, ss, tt, a, ca, sa;
  color col;

  if (texture1_map!="")
  {
    a = radians(t1_rotation);
    ca = cos(a);
    sa = sin(a);
    s0 = s-t1_uoffset-0.5;
    t0 = 0.5-t-t1_voffset;
    ss = ca*s0 - sa*t0;
    tt = sa*s0 + ca*t0;
    ss = t1_uscale*ss+0.5;
    tt = t1_vscale*tt+0.5;
    C_diffuse = texture(texture1_map, ss, tt, "blur", t1_blur);
  }

  if (specular_map!="")
  {
    a = radians(sp_rotation);
    ca = cos(a);
    sa = sin(a);
    s0 = s-sp_uoffset-0.5;
    t0 = 0.5-t-sp_voffset;
    ss = ca*s0 - sa*t0;
    tt = sa*s0 + ca*t0;
    ss = sp_uscale*ss+0.5;
    tt = sp_vscale*tt+0.5;
    C_specular = texture(specular_map, ss, tt, "blur", sp_blur);
  }

  if (shininess_map!="")
  {
    a = radians(sh_rotation);
    ca = cos(a);
    sa = sin(a);
    s0 = s-sh_uoffset-0.5;
    t0 = 0.5-t-sh_voffset;
    ss = ca*s0 - sa*t0;
    tt = sa*s0 + ca*t0;
    ss = sh_uscale*ss+0.5;
    tt = sh_vscale*tt+0.5;
    col = texture(shininess_map, ss, tt, "blur", sh_blur);
    final_shininess = (comp(col,0)+comp(col,1)+comp(col,2))/3;
  }

  if (opacity_map!="")
  {
    a = radians(op_rotation);
    ca = cos(a);
    sa = sin(a);
    s0 = s-op_uoffset-0.5;
    t0 = 0.5-t-op_voffset;
    ss = ca*s0 - sa*t0;
    tt = sa*s0 + ca*t0;
    ss = op_uscale*ss+0.5;
    tt = op_vscale*tt+0.5;
    float op = texture(opacity_map, ss, tt, "blur", op_blur);
    C_opacity = color "rgb" (op, op, op);
  }

  if (self_illum_map!="")
  {
    a = radians(si_rotation);
    ca = cos(a);
    sa = sin(a);
    s0 = s-si_uoffset-0.5;
    t0 = 0.5-t-si_voffset;
    ss = ca*s0 - sa*t0;
    tt = sa*s0 + ca*t0;
    ss = si_uscale*ss+0.5;
    tt = si_vscale*tt+0.5;
    col = texture(self_illum_map, ss, tt, "blur", si_blur);
    final_self_ilpct = (comp(col,0)+comp(col,1)+comp(col,2))/3;
  }

  if (reflection_map!="")
  {
    vector R = normalize(vtransform("world", reflect(I, Nf)));
    float lat = acos(zcomp(R));
    float long = acos(ycomp(R)/sin(lat));
    lat = lat/PI;
    long = long/(2*PI);
    if (xcomp(R)>0)
      long = -long;
    C_refl = texture(reflection_map, long, lat, "blur", refl_blur);
  }

  Ci = C_diffuse*diffuse(Nf) +
       shin_strength*C_specular*phong(Nf, -normalize(I), 25*final_shininess) +
       refl_percent*C_refl;
  Ci = mix(Ci, C_diffuse, final_self_ilpct);
  Oi = C_opacity;
  Ci *= Oi;
  BAKE_END
}        
        """

    def surfaceShaderParams(self, passes):
        """Return a dictionary with shader parameters and their values."""
        a = self.ambient
        d = self.diffuse
        s = self.specular
        res = {"uniform color ambient_col" : (a[0], a[1], a[2]),
                "uniform color diffuse_col" : (d[0], d[1], d[2]),
                "uniform color specular_col" : (s[0], s[1], s[2]),
                "uniform float shininess" : self.shininess,
                "uniform float shin_strength" : self.shin_strength,
                "uniform float self_ilpct" : self.self_ilpct
                }

        if self.texture1_map!=None:
            texname = os.path.basename(self.texture1_map.name)
            name, ext = os.path.splitext(texname)
            res["uniform string texture1_map"] = name+".tex"
            res["uniform float t1_uscale"] = self.texture1_map.scale[0]
            res["uniform float t1_vscale"] = self.texture1_map.scale[1]
            res["uniform float t1_uoffset"] = self.texture1_map.offset[0]
            res["uniform float t1_voffset"] = self.texture1_map.offset[1]
            res["uniform float t1_rotation"] = self.texture1_map.rotation

        if self.specular_map!=None:
            texname = os.path.basename(self.specular_map.name)
            name, ext = os.path.splitext(texname)
            res["uniform string specular_map"] = name+".tex"
            res["uniform float sp_uscale"] = self.specular_map.scale[0]
            res["uniform float sp_vscale"] = self.specular_map.scale[1]
            res["uniform float sp_uoffset"] = self.specular_map.offset[0]
            res["uniform float sp_voffset"] = self.specular_map.offset[1]
            res["uniform float sp_rotation"] = self.specular_map.rotation

        if self.shininess_map!=None:
            texname = os.path.basename(self.shininess_map.name)
            name, ext = os.path.splitext(texname)
            res["uniform string shininess_map"] = name+".tex"
            res["uniform float sh_uscale"] = self.shininess_map.scale[0]
            res["uniform float sh_vscale"] = self.shininess_map.scale[1]
            res["uniform float sh_uoffset"] = self.shininess_map.offset[0]
            res["uniform float sh_voffset"] = self.shininess_map.offset[1]
            res["uniform float sh_rotation"] = self.shininess_map.rotation

        if self.opacity_map!=None:
            texname = os.path.basename(self.opacity_map.name)
            name, ext = os.path.splitext(texname)
            res["uniform string opacity_map"] = name+".tex"
            res["uniform float op_uscale"] = self.opacity_map.scale[0]
            res["uniform float op_vscale"] = self.opacity_map.scale[1]
            res["uniform float op_uoffset"] = self.opacity_map.offset[0]
            res["uniform float op_voffset"] = self.opacity_map.offset[1]
            res["uniform float op_rotation"] = self.opacity_map.rotation

        if self.self_illum_map!=None:
            texname = os.path.basename(self.self_illum_map.name)
            name, ext = os.path.splitext(texname)
            res["uniform string self_illum_map"] = name+".tex"
            res["uniform float si_uscale"] = self.self_illum_map.scale[0]
            res["uniform float si_vscale"] = self.self_illum_map.scale[1]
            res["uniform float si_uoffset"] = self.self_illum_map.offset[0]
            res["uniform float si_voffset"] = self.self_illum_map.offset[1]
            res["uniform float si_rotation"] = self.self_illum_map.rotation

        if self.reflection_map!=None:
            texname = os.path.basename(self.reflection_map.name)
            name, ext = os.path.splitext(texname)
            res["uniform string reflection_map"] = name+".tex"
            res["uniform float refl_uscale"] = self.reflection_map.scale[0]
            res["uniform float refl_vscale"] = self.reflection_map.scale[1]
            res["uniform float refl_uoffset"] = self.reflection_map.offset[0]
            res["uniform float refl_voffset"] = self.reflection_map.offset[1]
            res["uniform float refl_rotation"] = self.reflection_map.rotation
            res["uniform float refl_percent"] = self.reflection_map.percent
            res["uniform float refl_blur"] = self.reflection_map.blur

        return res

    def surfaceShaderTransform(self):
        return mat4(1)

    def displacementShaderName(self):
        if self.bump_map==None:
            return None
        
        return "mat3ds_bump"

    
    def displacementShaderSource(self):
        if self.bump_map==None:
            return None
        
        return """// 3DS material shader (bump)
        
displacement $SHADERNAME(
           // Bump map
           string bump_map = "";
           float uscale = 1.0;
           float vscale = 1.0;
           float uoffset = 0;
           float voffset = 0;
           float rotation = 0;
           float blur = 0;
           float bump_size = 1.0;
           )
{
  float s0, t0, ss, tt, a, ca, sa;
  color bc;
  float amount = 0.0;

  if (bump_map!="")
  {
    a = radians(rotation);
    ca = cos(a);
    sa = sin(a);
    s0 = s-uoffset-0.5;
    t0 = 0.5-t-voffset;
    ss = ca*s0 - sa*t0;
    tt = sa*s0 + ca*t0;
    ss = uscale*ss+0.5;
    tt = vscale*tt+0.5;
    bc = texture(bump_map, ss, tt, "blur", blur);
    amount = (comp(bc,0)+comp(bc,1)+comp(bc,2))/3;
  }

  P = P - bump_size*amount*normalize(N);
  N = calculatenormal(P);
}
        """
    
    def displacementShaderParams(self, passes):
        if self.bump_map==None:
            return {}

        res = {}
        texname = os.path.basename(self.bump_map.name)
        name, ext = os.path.splitext(texname)
        res["uniform string bump_map"] = name+".tex"
        res["uniform float uscale"] = self.bump_map.scale[0]
        res["uniform float vscale"] = self.bump_map.scale[1]
        res["uniform float uoffset"] = self.bump_map.offset[0]
        res["uniform float voffset"] = self.bump_map.offset[1]
        res["uniform float rotation"] = self.bump_map.rotation
        res["uniform float blur"] = self.bump_map.blur
        res["uniform float bump_size"] = self.bump_size
        return res

    def displacementBound(self):
        if self.bump_map==None:
            return "current", 0
        else:
            return "current", self.bump_size

    def displacementShaderTransform(self):
        return mat4(1)
	
    def interiorShaderName(self):
        return None
    
    def interiorShaderSource(self):
        return None
    
    def interiorShaderParams(self, passes):
        return {}

    def interiorShaderTransform(self):
        return mat4(1)


##    def _copyImageMap(self, texmap, exporter):
##        """Copy the texture map image into the map folder.

##        \param texmap (\c TextureMap3DS) Texture map object or None
##        \param exporter Exporter instance
##        """
##        return
##        if texmap==None:
##            return

##        exporter.checkMapPath()
##        texname = os.path.basename(texmap.name)
##        name, ext = os.path.splitext(texname)
##        if ext.lower()!=".tif":
##            print 'Converting "%s"'%texmap.name
##            tifname = os.path.join(exporter.map_path, name+".tif")
##            # Read original map
##            try:
##                img = Image.open(texmap.name)
##            except IOError, e:
##                print e
##                return
##            # Save map as TIF file
##            img.save(tifname)
##        else:
##            print 'Copying "%s"'%texmap.name
##            shutil.copyfile(texmap.name, os.path.join(exporter.map_path, texmap.name))
         

    def _createTexDef(self, texmap):
        """Create texture definition for the TexPass.

        The method returns a list that is either empty or contains
        one texture definition tuple.

        \param texmap (\c TextureMap3DS) Texture map object or None
        \return An empty list or a list containing one texture def
        """
        if texmap==None:
            return []
        
#        texname = os.path.basename(texmap.name)
        texname = texmap.name
#        texname, ext = os.path.splitext(texname)
        blur = texmap.blur+1.0
        texdef = (texname, #+".tif",
                  "periodic", "periodic",
                  "gaussian", blur, blur, {})
        return [texdef]
        
