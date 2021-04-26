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

"""Dummy module that is used as a replacement for PyProtocols for the time being.
"""

import inspect

# Key:Class name - Value:Dict (key:interface obj, value:None or adapter class)
_interfaces = {}

def adapt(obj, interface):
    """Dummy implementation.
    """
    global _interfaces
    
    className = obj.__class__.__name__
    interfaces = _interfaces.get(className, {})
    
#    print "ADAPT",obj.__class__.__name__, interface, interfaces.has_key(interface)
    if interface in interfaces:
        return obj
    else:
        raise NotImplementedError("interface not implemented")

def advise(instancesProvide, asAdapterForTypes=None):
    """Dummy implementation.
    """
    global _interfaces
    
    frame = inspect.currentframe()
    outer = inspect.getouterframes(frame)
    currentClassName = outer[1][3]

#    print "ADVISE",currentClassName, asAdapterForTypes
    
    for interface in instancesProvide:
        if currentClassName not in _interfaces:
            _interfaces[currentClassName] = {}
        _interfaces[currentClassName][interface] = 1
#        print "    %s implements %s"%(currentClassName, interface)

class Interface:
    """Dummy interface base class.
    """
    pass
