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
# $Id: ifsimport.py,v 1.1.1.1 2004/12/12 14:31:04 mbaas Exp $

import os.path, struct
from .cgtypes import *
from .worldobject import WorldObject
from .trimesh import TriMesh
from .trimeshgeom import TriMeshGeom
from . import pluginmanager

# IfsImporter
class IfsImporter:
    """IFS importer.

    This class imports models from the Brown Mesh Set library which are
    stored in the Indexed Face Set (IFS) format.
    """

    _protocols = ["Import"]

    # extension
    def extension():
        """Return the file extensions for this format."""
        return ["ifs"]
    extension = staticmethod(extension)

    # description
    def description(self):
        """Return a short description for the file dialog."""
        return "Indexed face set"
    description = staticmethod(description)

    # importFile
    def importFile(self, filename):
        """Import an IFS file."""

        f = file(filename, "rb")

        # Check the header...
        s = self.readString(f)
        if s!="IFS":
            raise ValueError('The file "%s" is is not a IFS file.'%filename)

        # Read (and ignore) the version number
        s = f.read(4)
        ver = struct.unpack("<f", s)[0]

        # Read the model name
        modelname = self.readString(f)

        # Create the mesh geom
        tm = TriMeshGeom()

        # Read vertices...
        s = self.readString(f)
        if s!="VERTICES":
            raise ValueError("Vertices chunk expected, got '%s' instead."%s)

        s = f.read(4)
        numverts = int(struct.unpack("<I", s)[0])
        tm.verts.resize(numverts)

        for i in range(numverts):
            s = f.read(12)
            x,y,z = struct.unpack("<fff", s)
            tm.verts[i] = vec3(x,y,z)

        # Read faces...
        s = self.readString(f)
        if s!="TRIANGLES":
            raise ValueError("Triangle chunk expected, got '%s' instead."%s)
            
        s = f.read(4)
        numfaces = int(struct.unpack("<I", s)[0])
        tm.faces.resize(numfaces)

        for i in range(numfaces):
            s = f.read(12)
            a,b,c = struct.unpack("<III", s)
            tm.faces[i] = (int(a), int(b), int(c))

        # Create a world object
        obj = TriMesh(name=modelname)
        obj.geom = tm


    def readString(self, fhandle):
        """Read a string.

        \param fhandle Open file handle
        \return String
        """
        s = fhandle.read(4)
        w = int(struct.unpack("<I", s)[0])
        s = fhandle.read(w)
        # Return the string without the trailing \000
        return s[:-1]


######################################################################

# Register the IfsImporter class as a plugin class
pluginmanager.register(IfsImporter)
