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
# $Id: pyimport.py,v 1.4 2005/08/30 09:59:07 mbaas Exp $

import sys, os.path, copy
from . import pluginmanager
from . import globalscene
from . import eventmanager

# PyImporter
class PyImporter:

    _protocols = ["Import"]

    # extension
    def extension():
        return ["py"]
    extension = staticmethod(extension)

    # description
    def description(self):
        return "Python file"
    description = staticmethod(description)

    # importFile
    def importFile(self, filename):
        """Import a Python source file."""

        file_globals = {}
        # Use all globals from the cgkit package
        # Commented out the following line as is imported from cgkit.all now
#        file_globals = copy.copy(cgkit.__dict__)
        # Add some modules...
        exec("from cgkit.all import *", file_globals)
        exec("from cgkit.sl import *", file_globals)
        exec("from math import *", file_globals)
        # ...and some special global names...
        scene = globalscene.getScene()
        file_globals["scene"] = scene
        file_globals["timer"] = scene.timer()
        file_globals["worldroot"] = scene.worldRoot()
        file_globals["eventmanager"] = eventmanager.eventManager()
        
        paths = sys.path
        # Add the directory of the input file to the module search paths
        # so that local imports do work
        sys.path = [os.path.abspath(os.path.dirname(filename))] + sys.path
        # Import the file
        exec(compile(open(filename).read(), filename, 'exec'), file_globals)
        # Restore the old search paths
        sys.path = paths


######################################################################

# Register the OffImporter class as a plugin class
pluginmanager.register(PyImporter)
