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
# $Id: lwobimport.py,v 1.1 2006/03/12 22:41:42 mbaas Exp $

import os.path, sys
from .cgtypes import *
from .worldobject import WorldObject
from .trimesh import TriMesh
from .trimeshgeom import TriMeshGeom
from .polyhedron import Polyhedron
from .polyhedrongeom import PolyhedronGeom
from . import pluginmanager
from . import cmds
from cgkit.all import UNIFORM, INT, GLMaterial
from . import lwob


class _LWOBReader(lwob.LWOBReader):
    """Read a Lightwave Object file (*.lwo).
    """
    
    def __init__(self, parent=None):
        lwob.LWOBReader.__init__(self)

        # Parent node for the Lightwave object
        self.parent = parent

        # The TriMeshGeom that receives the triangle mesh
        self.trimeshgeom = TriMeshGeom()
        # The PolyhedronGeom that receives the mesh (if it is no triangle mesh)
        self.polyhedrongeom = PolyhedronGeom()
        # The generated WorldObject
        self.worldobj = None
        # The number of surfaces in the file
        self.numsurfaces = 0
        # A mapping from surface name to material id
        # Key: Surface name / Value: Material id (0-based)
        self.surface_ids = {}

        # Message flags so that warning messages are only output once
        self.crv_msg = False
        self.patch_msg = False

    def handlePNTS(self, points):
        """Handle the points chunk.

        Stores the points in the TriMeshGeom.
        """
        verts = self.trimeshgeom.verts
        verts.resize(len(points))
        for i,p in enumerate(points):
            verts[i] = p
    
    def handleSRFS(self, names):
        """Handle the surface names chunk.
        """
        self.numsurfaces = len(names)
        for i,n in enumerate(names):
            self.surface_ids[n] = i
    
    def handlePOLS(self, polys):
        """Handle the polygons.

        This method creates the actual object. It is assumed that the
        points have been read before and are stored in the TriMeshGeom.
        It is also assumed that a SRFS chunk was present and numsurfaces
        is initialized.
        """
        # Assume the mesh is a triangle mesh and initialize the TriMeshGeom
        # first. If this fails, use a PolyhedronGeom instead...
        if self._initTriMesh(polys):
            geom = self.trimeshgeom
        else:
            # Copy the vertices into the polyhedron geom...
            numverts = self.trimeshgeom.verts.size()
            self.polyhedrongeom.verts.resize(numverts)
            self.trimeshgeom.verts.copyValues(0, numverts, self.polyhedrongeom.verts, 0)
            del self.trimeshgeom
            # Initialize the polys...
            self._initPolyhedron(polys)
            geom = self.polyhedrongeom

        w = WorldObject(name="lwob", parent=self.parent)
        w.setNumMaterials(self.numsurfaces)
        w.geom = geom
        self.worldobj = w

    def _initPolyhedron(self, polys):
        """Initialize the polys of the PolyhedronGeom.

        Sets the faces a poly mesh and adds a matid slot with the
        material indices.
        """
        geom = self.polyhedrongeom
        geom.setNumPolys(len(polys))
        geom.newVariable("matid", UNIFORM, INT)
        matid = geom.slot("matid")
        for i,(verts,surfid) in enumerate(polys):
            geom.setLoop(i, 0, verts)
            matid[i] = max(0, surfid-1)

    def _initTriMesh(self, polys):
        """Initialize the faces of the TriMeshGeom.

        Sets the faces of a triangle mesh and adds a matid slot with
        the material indices.
        If the mesh contains faces with more than 3 vertices the
        method aborts and returns False.
        """
        faces = self.trimeshgeom.faces
        faces.resize(len(polys))
        self.trimeshgeom.newVariable("matid", UNIFORM, INT)
        matid = self.trimeshgeom.slot("matid")
        for i,(verts,surfid) in enumerate(polys):
            if len(verts)!=3:
                self.trimeshgeom.deleteVariable("matid")
                return False
            faces[i] = verts
            matid[i] = max(0, surfid-1)
            
        return True
        

    def handleCRVS(self, curves):
        if not self.crv_msg:
            print("Curves are not yet supported.")
            self.crv_msg = True
    
    def handlePCHS(self, patches):
        if not self.patch_msg:
            print("Patches are not yet supported.")
            self.patch_msg = True
    
    def handleSURF(self, surface):
        """Handle a surface chunk.

        Currently this just creates a GLMaterial with the base color of
        the surface. Everything else is ignored so far.
        """
        if surface.name not in self.surface_ids:
            raise lwob.LWOBError('Invalid surface name "%s" (name not available in SRFS chunk)'%surface.name)
        
        id = self.surface_ids[surface.name]
        
        col = surface.color
        if col==None:
            col = (255,255,255)
            
        mat = GLMaterial(diffuse=vec3(col)/255)
        self.worldobj.setMaterial(mat, id)


# LWOBImporter
class LWOBImporter:

    _protocols = ["Import"]

    # extension
    def extension():
        """Return the file extensions for this format."""
        return ["lwo"]
    extension = staticmethod(extension)

    # description
    def description(self):
        """Return a short description for the file dialog."""
        return "Lightwave object file"
    description = staticmethod(description)

    # importFile
    def importFile(self, filename, parent=None):
        """Import an LWOB file."""

        f = file(filename, "rb")
        reader = _LWOBReader(parent=parent)
        reader.read(f)
        f.close()

######################################################################

# Register the Importer class as a plugin class
pluginmanager.register(LWOBImporter)
