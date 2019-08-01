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
# $Id: x3dimport.py,v 1.2 2005/04/15 13:30:35 mbaas Exp $

import os.path
from . import _core
from .cgtypes import *
from .globalscene import getScene
from .worldobject import WorldObject
from .trimesh import TriMesh
from .trimeshgeom import TriMeshGeom
from .targetcamera import TargetCamera
from .glpointlight import GLPointLight
from .glfreespotlight import GLFreeSpotLight
from .glfreedistantlight import GLFreeDistantLight
from .glmaterial import GLMaterial, GLTexture
from .box import Box
from .quadrics import Sphere
from .group import Group
from . import pluginmanager
from .sl import *


# VRMLImporter
class VRMLImporter:

    _protocols = ["Import"]

    # extension
    def extension():
        """Return the file extensions for this format."""
        return ["wrl"]
    extension = staticmethod(extension)

    # description
    def description(self):
        """Return a short description for the file dialog."""
        return "Virtual Reality Modeling Language (VRML)"
    description = staticmethod(description)

    # importFile
    def importFile(self, filename):
        """Import a VRML file."""

        imp = _core.X3DReader(False)
        imp.read(filename)
        
#        _core.cyber(filename)
#        return

#        sg = _core.X3DSceneGraph()
#        if not sg.load(filename):
#            print "%s (%d): %s"%(sg.getParserErrorMessage(),
#                                 sg.getParserErrorLineNumber(),
#                                 sg.getParserErrorToken())
#            return

#        sg.foo("b = TriMeshGeom()")

# X3DImporter
class X3DImporter:

    _protocols = ["Import"]

    # extension
    def extension():
        """Return the file extensions for this format."""
        return ["x3d"]
    extension = staticmethod(extension)

    # description
    def description(self):
        """Return a short description for the file dialog."""
        return "Extensible 3D (X3D)"
    description = staticmethod(description)

    # importFile
    def importFile(self, filename):
        """Import a X3D file."""

        imp = _core.X3DReader()
        imp.read(filename)

     
######################################################################

# Register the Importer class as a plugin class
if hasattr(_core, "cyber"):
    pluginmanager.register(VRMLImporter)
    pluginmanager.register(X3DImporter)
