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
# Portions created by the Initial Developer are Copyright (C) 2010
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

import sys, os.path, ctypes

def findFfmpegLib(name):
    """Find a ffmpeg library.
    
    name is the name of the library that uses a dot to separate name from
    version number (e.g. "avutil.49").
    The function tries to find the library in some common location and it
    also tries the name where a dot has been replaced by a dash (e.g. "avutil-49").
    
    If the library can't be found, an exception is thrown, otherwise a ctypes
    handle to the library is returned.
    """
    
    # A list of library names that are tried to load the lib
    tries = [name]
    
    if not os.path.isabs(name):
        paths = []
        
        # OSX?
        if sys.platform=="darwin":
            paths = ["/usr/local/lib", "/opt/local/lib", "/sw/lib"]
            for path in paths:
                tries.append(os.path.join(path, "lib%s.dylib"%name))
        # Linux?
        elif sys.platform=="linux2":
            f = name.split(".")
            if len(f)==2:
                name = "lib%s.so.%s"%(f[0], f[1])
                tries.append(name)
            paths = ["/usr/local/lib"]
            for path in paths:
                tries.append(os.path.join(path, name))
        # Windows?
        elif sys.platform=="win32":
            winName = name.replace(".", "-")
            tries.append(winName)
            paths = []

    # Try all the alternatives we have...    
    exc = None
    for libName in tries:
        try:
            lib = ctypes.cdll[libName]
            return lib
        except OSError:
            if exc is None:
                exc = sys.exc_info()[1]
    
    raise exc
