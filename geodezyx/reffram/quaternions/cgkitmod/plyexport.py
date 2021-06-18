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
# $Id: plyexport.py,v 1.2 2005/05/08 22:17:42 mbaas Exp $

import os.path, sys
from .cgtypes import *
from .globalscene import getScene
from .geomobject import *
from .trimeshgeom import TriMeshGeom
from .polyhedrongeom import PolyhedronGeom
from . import pluginmanager
from . import cmds
from . import _core

# PLYExporter
class PLYExporter:

    _protocols = ["Export"]

    # extension
    def extension():
        """Return the file extensions for this format."""
        return ["ply"]
    extension = staticmethod(extension)

    # description
    def description(self):
        """Return a short description for the file dialog."""
        return "Polygon (PLY)"
    description = staticmethod(description)

    # exportFile
    def exportFile(self, filename, object=None, mode="ascii"):
        """Export a PLY file.

        object is the object to export. If it is None, it will be taken
        from the scene. If there is more than one object in the scene,
        an exception is thrown.
        mode specifies whether the output file will be ascii or binary.
        The values can be "ascii", "little_endian", "big_endian".
        """

        if object==None:
            # Get a list of all objects that have a geom
            objs = list(getScene().walkWorld())
            objs = [obj for obj in objs if obj.geom!=None]
            if len(objs)==0:
                raise ValueError("No object to export.")
            elif len(objs)>1:
                raise ValueError("Only a single object can be exported.")
            object = objs[0]
            
        object = cmds.worldObject(object)
        if object.geom==None:
            raise ValueError("No geometry attached to object %s"%object.name)
        geom = self.convertObject(object)
        if geom==None:
            raise ValueError("Cannot export geometry of type %s as a PLY file"%(object.geom.__class__.__name__))

        # Open the file...
        ply = _core.PLYWriter()
        try:
            mode = eval ("_core.PlyStorageMode.%s"%mode.upper())
        except:
            raise ValueError("Invalid mode: %s"%mode)
        ply.create(filename, mode)

        # Set comment
        var = geom.findVariable("comment")
        if var!=None and var[1]==CONSTANT and var[2]==STRING and var[3]==1:
            slot = geom.slot("comment")
            for s in slot[0].split("\n"):
                ply.addComment(s)
        # Set obj_info
        var = geom.findVariable("obj_info")
        if var!=None and var[1]==CONSTANT and var[2]==STRING and var[3]==1:
            slot = geom.slot("obj_info")
            for s in slot[0].split("\n"):
                ply.addObjInfo(s)

        # Write the model
        ply.write(geom, object.worldtransform)
        ply.close()

    # convertObject
    def convertObject(self, obj):
        """Converts an object into a polyhedron or trimesh if necessary.

        The return value is a GeomObject (TriMeshGeom or PolyhedronGeom)
        or None.
        """
        geom = obj.geom
        if isinstance(geom, TriMeshGeom):
            return geom
        
        if not isinstance(geom, PolyhedronGeom):
            # Try to convert into a polyhedron...
            pg = PolyhedronGeom()
            try:
                geom.convert(pg)
                geom = pg
            except:
                pass

        # Is it a PolyhedronGeom that has no polys with holes? then return
        # the geom...
        if isinstance(geom, PolyhedronGeom) and not geom.hasPolysWithHoles():
            return geom

        # Try to convert into a triangle mesh...
        tm = TriMeshGeom()
        try:
            geom.convert(tm)
            return tm
        except:
            pass

        return None
        

######################################################################

# Register the exporter class as a plugin class
pluginmanager.register(PLYExporter)
