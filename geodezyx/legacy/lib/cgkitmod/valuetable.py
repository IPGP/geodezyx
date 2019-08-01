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
# $Id: valuetable.py,v 1.1.1.1 2004/12/12 14:31:30 mbaas Exp $

## \file valuetable.py
## Contains the ValueTable component.

from .globalscene import getScene
from . import component
import bisect
from .cgtypes import *
from .slots import *

# TableEntry
class _TableEntry:
    """Stores a time/value pair for the ValueTable component.

    This is an internal class.
    """
    def __init__(self, t, v):
        self.t = t
        self.v = v

    def __str__(self):
        return "%f : %s"%(self.t, self.v)

    __repr__ = __str__

    def __cmp__(self, other):
        if other==None:
            return 1
        else:
            if self.t<other.t:
                return -1
            elif self.t>other.t:
                return 1
            else:
                return 0
            
# ValueTable
class ValueTable(component.Component):
    """ValueTable component.

    This class stores time/value pairs and has an output slot that
    holds the appropriate value for the current time. The type of
    the value can be specified in the constructor. The name of the
    output slot is always \c output_slot.
    """
    
    def __init__(self,
                 name = "ValueTable",
                 type = "vec3",
                 values = [],
                 modulo = None,
                 tscale = 1.0,
                 auto_insert = True):
        """Constructor.

        \param name (\c str) Component name
        \param type (\c str) Value type
        \param values A list of tuples (time, value).
        \param modulo (\c float) Loop duration (None = no loop)
        \param tscale (\c float) Scaling factor for the time. A value of less than 1.0 makes the animation slower.
        """
        
        component.Component.__init__(self, name=name, auto_insert=auto_insert)

        # Value list. Contains a sorted list of TableEntry objects
        self.values = []
        # Time modulo value (or None)
        self.modulo = modulo
        # Scale factor for the time
        self.tscale = 1.0
        # Type of the value slot
        self.type = type
        
        self.time_slot = DoubleSlot()
        self.addSlot("time", self.time_slot)
        typ = type.lower()
        exec("self.output_slot = Procedural%sSlot(self.computeValue)"%typ.capitalize())
        self.addSlot("output", self.output_slot)
        pytypes = {"double":"float"}
        exec("self.default_value = %s()"%pytypes.get(typ, typ))

        self.time_slot.addDependent(self.output_slot)
        getScene().timer().time_slot.connect(self.time_slot)

        for t,v in values:
            self.add(t,v)

    def __iter__(self):
        return self.iterValues()

    def __call__(self, time):
        if len(self.values)==0:
            return self.default_value

        time *= self.tscale
        if self.modulo!=None:
            time = time % self.modulo

        idx = bisect.bisect_left(self.values, _TableEntry(time,self.default_value))
        if idx>=len(self.values):
            idx = len(self.values)-1
        e = self.values[idx]
        if time<e.t:
            if idx>0:
                return self.values[idx-1].v
        return e.v       

    def __getitem__(self, time):
        return self(time)

    def __setitem__(self, time, value):
        self.add(time, value)

    # iterValues
    def iterValues(self):
        """Iterate over all time/value pairs.
        """
        for e in self.values:
            yield e.t, e.v

    # add
    def add(self, t, v):
        """Add a value to the table.

        \param t (\c float) Time
        \param v Value
        """
        entry = _TableEntry(t, v)
        idx = bisect.bisect_left(self.values, entry)
        # Check if times are identical and the previous value has
        # to be replaced
        if idx<len(self.values):
            if self.values[idx].t==t:
                # Replace the value
                self.values[idx] = entry
                return

        # Insert the value
        self.values.insert(idx, entry)

    ## protected:
        
    def computeValue(self):
        """Computes a new output value."""
        return self(self.time_slot.getValue())

    exec(slotPropertyCode("output"))
