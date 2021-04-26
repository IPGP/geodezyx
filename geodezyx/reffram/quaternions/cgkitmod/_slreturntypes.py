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
# Portions created by the Initial Developer are Copyright (C) 2008
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

"""This module contains the return types used by the slparams and sloargs module.
"""

class _ShaderInfo:
    """Holds information about a shader.
    
    For backwards compatibility this class behaves like a tuple
    (type, name, params). The meta data has to be accessed via attribute
    access.
    """
    
    def __init__(self, type=None, name=None, params=None, meta=None):
        """Constructor.
        """
        if params is None:
            params = []
        if meta is None:
            meta = {}
        
        # The shader type (surface, displacement, ...)
        self.type = type
        # The shader name
        self.name = name
        # The shader parameters (a list of _ShaderParam objects)
        self.params = params
        # The meta data attached to the shader (a dictionary)
        self.meta = meta
        
    def __str__(self):
        return "(%r, %r, %s)"%(self.type, self.name, self.params)
    
    __repr__ = __str__
        
    def __len__(self):
        return 3
    
    def __iter__(self):
        yield self.type
        yield self.name
        yield self.params
        
    def __getitem__(self, idx):
        return (self.type, self.name, self.params)[idx]


class _ShaderParam:
    """Holds information about a shader parameter.
    
    For backwards compatibility this class behaves like a 7-tuple
    (outSpec, storage, type, size, name, space, default). 
    """
    
    def __init__(self, outputSpec=None, storage=None, type=None, size=None, name=None, space=None, default=None):
        """Constructor.
        """
        
        # The output specifier ("output" or empty string)
        self.outputSpec = outputSpec
        # The storage class ("uniform", "varying")
        self.storage = storage
        # The parameter type
        self.type = type
        # The array length or None if the param is not an array
        self.size = size
        # The parameter name
        self.name = name
        # The space (or list of spaces in case of an array) or None
        self.space = space
        # The default value
        self.default = default
        
    def __str__(self):
        return "(%r, %r, %r, %r, %r, %r, %r)"%(self.outputSpec, self.storage, self.type, self.size, self.name, self.space, self.default)
    
    __repr__ = __str__

    def __len__(self):
        return 7
    
    def __iter__(self):
        yield self.outputSpec
        yield self.storage
        yield self.type
        yield self.size
        yield self.name
        yield self.space
        yield self.default
        
    def __getitem__(self, idx):
        return (self.outputSpec, self.storage, self.type, self.size, self.name, self.space, self.default)[idx]

