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
# $Id: expression.py,v 1.3 2005/10/18 07:02:39 mbaas Exp $

## \file expression.py
## Contains the Expression class.

"""This module contains the Expression class."""

from .Interfaces import *
from . import protocols
from . import slots
from .eventmanager import eventManager
from .globalscene import getScene
from . import events
from .cgtypes import *
from .component import *
from .sl import *
from math import *
from . import _core
       

# Expression
class Expression(Component):
    """Compute an output value using a user defined expression.

    This component outputs a single value that is driven by a user defined
    expression. The expression is specified by a string and can use an
    arbitrary number of parameters. The parameters and their default values
    have to be provided to the constructor via keyword arguments. An exception
    is the special variable "t" which will always hold the current time
    (unless you declare it explicitly).
    For each parameter a slot is created (<name>_slot), so it is also
    possible to animate the parameters. The output value can be accessed via
    the "output" and "output_slot" attributes.

    Example:

    \code
    s = Sphere()
    e = Expression("1.0 + amp*sin(freq*t)", amp=0.2, freq=2.0)
    e.output_slot.connect(s.radius_slot)
    \endcode
    
    """

    protocols.advise(instancesProvide=[ISceneItem])

    def __init__(self,
                 expr = "",
                 exprtype = None,
                 name = "Expression",
                 **keyargs
                 ):
        """Constructor.

        If no expression type is given the component tries to determine
        the type itself by executing the expression and inspecting the
        return type.

        \param expr (\c str) Expression
        \param exprtype (\c str) Output type or None
        \param name (\c str) Component name
        \param keyargs Parameters used in the expression
        """
        Component.__init__(self, name=name)

        self.expr = expr
        self.exprtype = exprtype

        # Create a parameter slot for every extra key arg...
        for k in keyargs:
            T = type(keyargs[k])
            if T==float or T==int:
                typ = "double"
                valstr = str(keyargs[k])
            elif T==vec3:
                typ = "vec3"
                x, y, z = keyargs[k]
                valstr = "vec3(%s, %s, %s)"%(x,y,z)
            elif T==vec4:
                typ = "vec4"
                x, y, z, w = keyargs[k]
                valstr = "vec4(%s, %s, %s, %s)"%(x,y,z,w)
            elif T==mat3:
                typ = "mat3"
                valstr = "mat3(%s)"%keyargs[k].toList(rowmajor=True)
            elif T==mat4:
                typ = "mat4"
                valstr = "mat4(%s)"%keyargs[k].toList(rowmajor=True)
            elif T==quat:
                typ = "quat"
                w, x, y, z = keyargs[k]
                valstr = "quat(%s, %s, %s, %s)"%(w,x,y,z)
            else:
                typ = "py"
                valstr = "keyargs[k]"
#                raise ValueError("Unsupported type: %s"%T)
            # Create slot
            exec("self.%s_slot = %sSlot(%s)"%(k, typ.capitalize(), valstr))
            exec("self.addSlot(k, self.%s_slot)"%k)

        # If t was not explicitly given use the timer...
        if "t" not in keyargs:
            self.t_slot = DoubleSlot()
            getScene().timer().time_slot.connect(self.t_slot)
            self.addSlot("t", self.t_slot)

        # Store a list of all parameter names
        self.vars = list(keyargs.keys())
        if "t" not in self.vars:
            self.vars.append("t")

        if self.exprtype==None:
            self.exprtype = self._determineReturnType()

        # Create the output slot
        e = self.exprtype
        if e.lower()=="float":
            e = "double"
        exec("self.output_slot = Procedural%sSlot(self.outProc)"%e.capitalize())
        self.addSlot("output", self.output_slot)

        # Create dependencies
        for v in self.vars:
            exec("self.%s_slot.addDependent(self.output_slot)"%(v))

            

    def protocols(self):
        return [ISceneItem, IComponent]

    def outProc(self):
        for _v in self.vars:
            exec("%s = self.%s_slot.getValue()"%(_v, _v))
        return eval("%s(%s)"%(self.exprtype, self.expr))
     
    ## protected:
        
    # "output" property...
    exec(slotPropertyCode("output"))

    def _determineReturnType(self):
        """Try to execute the stored expression and return the output type.
        """
        for _v in self.vars:
            exec("%s = self.%s_slot.getValue()"%(_v, _v))
        out = eval(self.expr)
        T = type(out)
        if T==float or T==int:
            return "float"
        if isinstance(out, _core.vec3):
            return "vec3"
        if isinstance(out, _core.vec4):
            return "vec4"
        if isinstance(out, _core.mat3):
            return "mat3"
        if isinstance(out, _core.mat4):
            return "mat4"
        if isinstance(out, _core.quat):
            return "quat"

        if T==tuple or T==list:
            if len(out)==3:
                return "vec3"
            if len(out)==4:
                return "vec4"
            if len(out)==9:
                return "mat3"
            if len(out)==16:
                return "mat4"
            raise ValueError("Unsupported sequence size: %d"%len(out))

        raise ValueError("Unknown expression type: %s"%T)
        
        
        

        

