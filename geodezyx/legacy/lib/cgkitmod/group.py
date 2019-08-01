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
# $Id: group.py,v 1.3 2005/04/19 08:45:56 mbaas Exp $

## \file group.py
## Contains the Group class.

from .Interfaces import *
from . import protocols
from .slots import *
from .worldobject import WorldObject
from .globalscene import getScene

# Group
class Group(WorldObject):
    """Group class.

    This class is almost identical to the WorldObject class, except
    that its constructor accepts additional arguments.
    """

    protocols.advise(instancesProvide=[ISceneItem, IRigidBody])

    def __init__(self,
                 name="group",
                 dynamics=True,
                 static=False,
                 childs=[],
                 **params):
        """Constructor.

        \param name (\c str) Object name
        \param childs (\c list) A list of WorldObjects or object names.
        """
        WorldObject.__init__(self, name=name, **params)

        self.dynamics_slot = BoolSlot(dynamics)
        self.static_slot = BoolSlot(static)
        self.addSlot("dynamics", self.dynamics_slot)
        self.addSlot("static", self.static_slot)
        
        for c in childs:
            c = _worldObject(c)
            if c.parent!=None:
                c.parent.removeChild(c)
            self.addChild(c)

    exec(slotPropertyCode("dynamics"))
    exec(slotPropertyCode("static"))



def _worldObject(obj):
    """Return a world object.

    (this has been pulled in from the cmds module, so that we don't need
    to import the cmds module (producing a cycle). This needs to be fixed
    in a future version)
    """
    
    if isinstance(obj, str):
        return getScene().worldObject(obj)
    else:
        return obj
