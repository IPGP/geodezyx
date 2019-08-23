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
# Portions created by the Initial Developer are Copyright (C) 2009
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
# $Id: riutil.py,v 1.1.1.1 2004/12/12 14:31:21 mbaas Exp $

import sys, os, os.path
import ctypes.util

def resolveRManLib(libName, renderer=None):
    """Resolve the given library name.

    If the name is an absolute file name, it is just returned unmodified.
    Otherwise the method tries to resolve the name and return an absolute
    path to the library. If no library file could be found, the name
    is returned unmodified.
    
    renderer may be the name of a renderer package (aqsis, pixie, 3delight,
    prman) which serves as a hint to which package the library belongs to.
    If not given, the renderer is determined from the library name.
    """
    if os.path.isabs(libName):
        return libName
    
    # Try to figure out the location of the lib
    lib = ctypes.util.find_library(libName)
    if lib is not None:
        return lib

    # A list of library search paths...
    searchPaths = []

    # Is there a renderer-specific search path?
    if renderer is None:
        renderer = rendererFromLib(libName)
    libDir = rendererLibDir(renderer)
    if libDir is not None:
        searchPaths.append(libDir)

    # Also examine LD_LIBRARY_PATH if we are on Linux
    if sys.platform.startswith("linux"):
        libPaths = os.getenv("LD_LIBRARY_PATH")
        if libPaths is not None:
            searchPaths.extend(libPaths.split(":"))

    # Check the search paths...
    libFileName = libraryFileName(libName)
    for path in searchPaths:
        lib = os.path.join(path, libFileName)
        if os.path.exists(lib):
            return lib

    # Nothing found, then just return the original name
    return libName

def rendererFromLib(libName):
    """Return the renderer package name given a library name.
    
    libName is the name of a RenderMan library. The function tries to
    determine from which renderer it is and returns the name of the
    render package. None is returned if the library name is unknown.
    """
    name = os.path.basename(libName)
    name = os.path.splitext(name)[0]
    if name.startswith("lib"):
        name = name[3:]
        
    if name in ["aqsislib", "ri2rib", "slxargs"]:
        return "aqsis"
    elif libName in ["sdr", "ri"]:
        return "pixie"
    elif libName in ["3delight"]:
        return "3delight"
    elif libName in ["prman"]:
        return "prman"
    
    return None

def rendererLibDir(renderer):
    """Return a renderer-specific library path.
    
    The return path is based on a renderer-specific environment variable.
    None is returned when no path could be determined.
    renderer may be None, "aqsis", "pixie", "3delight" or "prman" (case-insensitive).
    """
    if renderer is None:
        return None
    
    envVarDict = {"aqsis":"AQSISHOME",
                  "pixie":"PIXIEHOME",
                  "3delight":"DELIGHT",
                  "prman":"RMANTREE"}
    
    envVar = envVarDict.get(renderer.lower())
    if envVar is not None:
        base = os.getenv(envVar)
        if base is not None:
            return os.path.join(base, "lib")
    return None

def libraryFileName(libName):
    """Extend a base library name to a file name.

    Example:  "foo" -> "libfoo.so"    (Linux)
                    -> "foo.dll"      (Windows)
                    -> "libfoo.dylib" (OSX)
    """
    if sys.platform.startswith("linux"):
        return "lib%s.so"%libName
    elif sys.platform=="darwin":
        return "lib%s.dylib"%libName
    elif sys.platform.startswith("win"):
        return "%s.dll"%libName
    return libName

#####################################################################
if __name__=='__main__':
    pass
