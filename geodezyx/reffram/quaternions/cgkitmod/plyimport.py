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
# $Id: plyimport.py,v 1.4 2005/05/09 14:34:24 mbaas Exp $

import os.path
from . import _core
from .cgtypes import *
from .geomobject import *
from .polyhedron import *
from . import pluginmanager


# PLYImporter
class PLYImporter:

    _protocols = ["Import"]

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

    # importFile
    def importFile(self, filename, includevar=None, excludevar=None, invertfaces=False):
        """Import a PLY file.

        includevar is a list of primitive variable names that should be
        imported. excludevar is a list of primitive variable names that
        should not be imported.
        invertfaces specifies if the face orientations should be inverted
        or not.
        """

        self.includevar = includevar
        self.excludevar = excludevar

        imp = _core.PLYReader()
        imp.open(filename)
        # Obtain the header information
        header = imp.readHeader()

        # Create the variable declaration tuples for importing the data...
        inttypes = ["int8", "uint8", "int16", "uint16", "int32", "uint32",
                    "char", "uchar", "short", "ushort", "int", "uint"]
        vardecl = []
        vec3vars = {}
        elements, comment, objinfo = header
        for elname, ninstances, properties in elements:
            for propname, type, len_type, val_type in properties:
                # Vertex indices? (they are handled by default)
                if elname=="vertex" and propname in ["x", "y", "z"]:
                    continue
                # Vertex indices? (they are handled by default)
                if elname=="face" and propname in ["vertex_indices"]:
                    continue
                # Normals?
                if propname in ["nx", "ny", "nz"]:
                    if "N" not in vec3vars and self.isVarAccepted("N"):
                        vardecl.append(("N", NORMAL, elname, ("nx", "ny", "nz")))
                        vec3vars["N"] = True
                    continue

                if type=="list":
                    type = val_type
                if type in inttypes:
                    t = INT
                else:
                    t = FLOAT
                if self.isVarAccepted(propname):
                    vardecl.append((propname, t, elname, (propname,)))

#        print vardecl

        # Create a polyhedron
        name = os.path.splitext(os.path.basename(filename))[0]
        p = Polyhedron(name=name)
        # Set the comment and objinfo
        if comment!="":
            p.geom.newVariable("comment", CONSTANT, STRING)
            s = p.geom.slot("comment")
            s[0] = comment
        if objinfo!="":
            p.geom.newVariable("obj_info", CONSTANT, STRING)
            s = p.geom.slot("obj_info")
            s[0] = objinfo

        # Read the model
        imp.read(p.geom, vardecl, invertfaces)

        imp.close()

    # isVarAccepted
    def isVarAccepted(self, name):
        """Return True if the variable should be imported.

        name is the name of the primitive variable (which is not necessarily
        the ply property!).
        """
        if self.includevar!=None:
            if name not in self.includevar:
                return False
        if self.excludevar!=None:
            if name in self.excludevar:
                return False
        return True

######################################################################

# Register the Importer class as a plugin class
pluginmanager.register(PLYImporter)
