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
# $Id: dddsimport.py,v 1.5 2005/06/07 07:43:03 mbaas Exp $

import os.path
from . import _core
from .cgtypes import *
from .globalscene import getScene
from .worldobject import WorldObject
from .trimesh import TriMesh
from .trimeshgeom import TriMeshGeom
from .targetcamera import TargetCamera
from .glpointlight import GLPointLight
from .gltargetspotlight import GLTargetSpotLight
from .spotlight3ds import SpotLight3DS
from .glmaterial import GLMaterial
from .material3ds import Material3DS, TextureMap3DS
from . import pluginmanager
from .sl import *

TYPE_UNKNOWN = 0
TYPE_AMBIENT = 1
TYPE_OBJECT = 2
TYPE_CAMERA = 3
TYPE_TARGET = 4
TYPE_LIGHT = 5
TYPE_SPOT = 6

GEOM_INIT_NORMALS   = 0x01
GEOM_INIT_MATIDS    = 0x02
GEOM_INIT_TEXELS    = 0x04
GEOM_INIT_FLAGS     = 0x08
GEOM_INIT_SMOOTHING = 0x10
GEOM_INIT_ALL       = 0xff


# DDDSImporter
class DDDSImporter:

    _protocols = ["Import"]

    # extension
    def extension():
        """Return the file extensions for this format."""
        return ["3ds"]
    extension = staticmethod(extension)

    # description
    def description(self):
        """Return a short description for the file dialog."""
        return "3D Studio"
    description = staticmethod(description)

    # importFile
    def importFile(self, filename, flags=GEOM_INIT_ALL, parent=None):
        """Import a 3DS file."""

        self.filename = filename
        self.ddds = _core.File3ds()
        self.ddds.load(filename)
#        f = self.ddds.current_frame
        f = getScene().timer().frame
        if f<self.ddds.segment_from:
            f = self.ddds.segment_from
        if f>self.ddds.segment_to:
            f = self.ddds.segment_to
        self.ddds.eval(f)

        # Create a dictionary containing all available meshes.
        # Key is the mesh name. This is used to check if all meshes have
        # been processed (i.e. if there were corresponding nodes).
        self.meshes = {}
        m = self.ddds.meshes()
        while m!=None:
            if m.name in self.meshes:
                print("Warning: Duplicate mesh names in 3ds file")
            self.meshes[m.name] = m
            m = next(m)

        # Create objects...
        self.materials = {}
        self.createObject(self.ddds.nodes(), parent=parent, flags=flags)
        del self.materials

        # Create TriMeshes for all left over meshes...
        for n in self.meshes:
            mesh = self.meshes[n]
            tm = TriMeshGeom()
            mesh.initGeom(tm, flags)
            worldobj = TriMesh(name = mesh.name, parent=parent)
            worldobj.geom = tm

        self.ddds.free()

        
    def createObject(self, node, parent=None, flags=GEOM_INIT_ALL):
        """
        \param node (\c Lib3dsNode) Current node
        \param parent (\c WorldObject) Parent object or None
        \param flags (\c int) Flags for mesh creation

        \todo worldtransform als Slot
        \todo Pivot so setzen wie in 3DS
        """
        if node==None:
            return

        worldobj = None
        auto_insert = True
        if parent!=None:
            auto_insert = False

        # Object?
        if node.type==TYPE_OBJECT:
            if node.name!="$$$DUMMY" and node.name!="_Quader01":
                data = node.object_data
                mesh = self.ddds.meshByName(node.name)
                if mesh.name in self.meshes:
                    del self.meshes[mesh.name]
#                print "###",node.name,"###"
#                print "Node matrix:"
#                print node.matrix
#                print "Pivot:",data.pivot
#                print "Mesh matrix:"
#                print mesh.matrix
                tm = TriMeshGeom()
                mesh.initGeom(tm, flags)

                if parent==None:
                    PT = mat4().translation(-data.pivot)
                    m = node.matrix*PT*mesh.matrix.inverse()
                else:
                    PT = mat4().translation(-data.pivot)
                    m = node.matrix*PT*mesh.matrix.inverse()
                    # worldtransform von Parent bestimmen...
                    pworldtransform = parent.worldtransform
#                    pworldtransform = parent.localTransform()
#                    p = parent.parent
#                    while p!=None:
#                        pworldtransform = p.localTransform()*pworldtransform
#                        p = p.parent
                    m = pworldtransform.inverse()*m
                worldobj = TriMesh(name=node.name,
#                                   pivot=data.pivot,
                                   transform=m,
                                   auto_insert=auto_insert)
                worldobj.geom = tm

                # Set the materials
                matnames = mesh.materialNames()
                worldobj.setNumMaterials(len(matnames))
                for i in range(len(matnames)):
                    if matnames[i]=="":
                        material = GLMaterial()  
                    # Was the material already instantiated?
                    elif matnames[i] in self.materials:
                        material = self.materials[matnames[i]]
                    else:
                        mat = self.ddds.materialByName(matnames[i])
#                        print "Material:",matnames[i]
#                        print "  self_illum:",mat.self_illum
#                        print "  self_ilpct:",mat.self_ilpct
                        material = Material3DS(name = matnames[i],
                                               ambient = mat.ambient,
                                               diffuse = mat.diffuse,
                                               specular = mat.specular,
                                               shininess = mat.shininess,
                                               shin_strength = mat.shin_strength,
                                               use_blur = mat.use_blur,
                                               transparency = mat.transparency,
                                               falloff = mat.falloff,
                                               additive = mat.additive,
                                               use_falloff = mat.use_falloff,
                                               self_illum = mat.self_illum,
                                               self_ilpct = getattr(mat, "self_ilpct", 0),
                                               shading = mat.shading,
                                               soften = mat.soften,
                                               face_map = mat.face_map,
                                               two_sided = mat.two_sided,
                                               map_decal = mat.map_decal,
                                               use_wire = mat.use_wire,
                                               use_wire_abs = mat.use_wire_abs,
                                               wire_size = mat.wire_size,
                                               texture1_map = self._createTexMap(mat.texture1_map),
                                               texture2_map = self._createTexMap(mat.texture2_map),
                                               opacity_map = self._createTexMap(mat.opacity_map),
                                               bump_map = self._createTexMap(mat.bump_map),
                                               specular_map = self._createTexMap(mat.specular_map),
                                               shininess_map = self._createTexMap(mat.shininess_map),
                                               self_illum_map = self._createTexMap(mat.self_illum_map),
                                               reflection_map = self._createTexMap(mat.reflection_map)
                                               )
                        self.materials[matnames[i]] = material
                    worldobj.setMaterial(material, i)
                
        # Camera?
        elif node.type==TYPE_CAMERA:
            cam = self.ddds.cameraByName(node.name)
#            print "Camera:",node.name
#            print "  Roll:",cam.roll
#            print node.matrix
            # Convert the FOV from horizontal to vertical direction
            # (Todo: What aspect ratio to use?)
            fov = degrees(atan(480/640.0*tan(radians(cam.fov/2.0))))*2.0
            worldobj = TargetCamera(name = node.name,
                                    pos = cam.position,
                                    target = cam.target,
                                    fov = fov,
                                    auto_insert = auto_insert)
        # Light?
        elif node.type==TYPE_LIGHT:
            lgt = self.ddds.lightByName(node.name)
            worldobj = self._createLight(node, lgt, auto_insert)

        if worldobj!=None and parent!=None:
            parent.addChild(worldobj)

        if worldobj!=None:
            self.createObject(node.childs(), worldobj, flags=flags)
        else:
            self.createObject(node.childs(), parent, flags=flags)
        self.createObject(next(node), parent, flags=flags)
        

    ## protected:

    def _createLight(self, node, lgt, auto_insert):
        """Create a light source.

        Helper method for the createObject() method.
        """
        shadow_size = lgt.shadow_size
        shadow_bias = lgt.shadow_bias
        shadow_filter = lgt.shadow_filter
        # Check if a shadow map is used but the parameters are all zero
        # (this happens when the Use Global Settings checkbox in MAX
        # is set). If it happens then use the same default values as MAX
        if lgt.shadowed and shadow_size==0 and shadow_bias==0 and shadow_filter==0:
            shadow_size = 256
            shadow_bias = 3.0
            shadow_filter = 5.0
            
        if lgt.spot_light:
            worldobj = SpotLight3DS(name = node.name,
                                    pos = lgt.position,
                                    target = lgt.spot,
                                    enabled = not lgt.off,
                                    intensity = lgt.multiplier,
                                    see_cone = lgt.see_cone,
                                    roll = lgt.roll, 
                                    outer_range = lgt.outer_range,
                                    inner_range = lgt.inner_range,
                                    attenuation = lgt.attenuation,
                                    rectangular_spot = lgt.rectangular_spot,
                                    shadowed = lgt.shadowed,
                                    shadow_bias = shadow_bias,
                                    shadow_filter = shadow_filter,
                                    shadow_size = shadow_size,
                                    spot_aspect = lgt.spot_aspect,
                                    use_projector = lgt.use_projector,
                                    projector = lgt.projector,
                                    overshoot = lgt.spot_overshoot,
                                    ray_shadows = lgt.ray_shadows,
                                    ray_bias = lgt.ray_bias,
                                    hotspot = lgt.hotspot,
                                    falloff = lgt.falloff,
                                    color = lgt.color,
                                    auto_insert = auto_insert)
#            worldobj = GLTargetSpotLight(name = node.name,
#                                         pos = lgt.position,
#                                         target = lgt.spot,
#                                         cutoff = lgt.falloff/2,
#                                         intensity = lgt.multiplier,
#                                         diffuse = lgt.color,
#                                         auto_insert = auto_insert)
#            if lgt.shadowed:
#                worldobj.shadowmap = (lgt.shadow_size, 1.0, 0.001)
        else:
            worldobj = SpotLight3DS(name = node.name,
                                    pos = lgt.position,
                                    target = lgt.spot,
                                    enabled = not lgt.off,
                                    intensity = lgt.multiplier,
                                    see_cone = lgt.see_cone,
                                    roll = lgt.roll, 
                                    outer_range = lgt.outer_range,
                                    inner_range = lgt.inner_range,
                                    attenuation = lgt.attenuation,
                                    rectangular_spot = lgt.rectangular_spot,
                                    shadowed = False,
                                    shadow_bias = shadow_bias,
                                    shadow_filter = shadow_filter,
                                    shadow_size = shadow_size,
                                    spot_aspect = lgt.spot_aspect,
                                    use_projector = lgt.use_projector,
                                    projector = lgt.projector,
                                    overshoot = True,
                                    ray_shadows = lgt.ray_shadows,
                                    ray_bias = lgt.ray_bias,
                                    hotspot = lgt.hotspot,
                                    falloff = lgt.falloff,
                                    color = lgt.color,
                                    auto_insert = auto_insert)
#            worldobj = GLPointLight(name = node.name,
#                                    pos = lgt.position,
#                                    intensity = lgt.multiplier,
#                                    diffuse = lgt.color,
#                                    auto_insert = auto_insert)
            
        return worldobj
        
    def _createTexMap(self, texmap):
        """Create a TextureMap3DS object (or None).

        \param texmap (\c Lib3dsTextureMap) Texture map
        """
        if texmap.name=="":
            return None

        dir = os.path.dirname(self.filename)
        dir = os.path.abspath(dir)
        return TextureMap3DS(name=os.path.join(dir, texmap.name),
                             flags = texmap.flags,
                             percent = texmap.percent,
                             blur = texmap.blur,
                             scale = texmap.scale,
                             offset = texmap.offset,
                             rotation = texmap.rotation,
                             tint1 = texmap.tint_1,
                             tint2 = texmap.tint_2,
                             tintr = texmap.tint_r,
                             tintg = texmap.tint_g,
                             tintb = texmap.tint_b)

        


######################################################################

# Register the Importer class as a plugin class
if hasattr(_core, "File3ds"):
    pluginmanager.register(DDDSImporter)
