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

# Shadow module for module "Image"

try:

    # Try to import the original module...
    from Image import *
    
except ImportError:

    # Create dummy symbols...

    ADAPTIVE = 1
    AFFINE = 0
    ANTIALIAS = 1
    BICUBIC = 3
    BILINEAR = 2
    CONTAINER = 2
    CUBIC = 3
    DEBUG = 0
    EXTENSION = {}
    EXTENT = 1
    FLIP_LEFT_RIGHT = 0
    FLIP_TOP_BOTTOM = 1
    FLOYDSTEINBERG = 3
    ID = []
    LINEAR = 2
    MESH = 4
    MIME = {}
    MODES = ['1', 'CMYK', 'F', 'I', 'L', 'P', 'RGB', 'RGBA', 'RGBX', 'YCbCr']
    NEAREST = 0
    NONE = 0
    NORMAL = 0
    OPEN = {}
    ORDERED = 1
    PERSPECTIVE = 2
    QUAD = 3
    RASTERIZE = 2
    ROTATE_180 = 3
    ROTATE_270 = 4
    ROTATE_90 = 2
    SAVE = {}
    SEQUENCE = 1
    VERSION = '1.1.4'
    WEB = 0

    def Image(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def IntType(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def StringType(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def TupleType(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def UnicodeStringType(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def blend(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def composite(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def eval(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def frombuffer(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def fromstring(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def getmodebands(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def getmodebase(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def getmodetype(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def init(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def isDirectory(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def isImageType(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def isNumberType(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def isSequenceType(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def isStringType(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def isTupleType(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def merge(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def new(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def open(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def preinit(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def register_extension(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def register_mime(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def register_open(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

    def register_save(*args, **keyargs):
        raise ImportError("No module named Image. Please install PIL (http://www.pythonware.com/products/pil/index.htm).")

