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
# $Id: component.py,v 1.3 2005/08/15 15:47:56 mbaas Exp $

## \file component.py
## Contains the Component base class.

from ._core import Component as _Component
from .Interfaces import ISceneItem
import copy, inspect, re
from . import protocols
from .slots import *
from .globalscene import getScene

# Component
class Component(_Component):
    """Base component class.

    Attributes:

    - name (\c str): The name of the component
    """

    protocols.advise(instancesProvide=[ISceneItem])
    
    def __init__(self, name="", auto_insert=True):
        _Component.__init__(self, name)
        
        if auto_insert:
            getScene().insert(self)

    def protocols(self):
        return [ISceneItem]

    def connect(self, srcslot, dstslot):
        dstslot.setController(srcslot)

############

def _parseFuncDecl(funcdecl):
    """Parse a function declaration string.

    The syntax of the string is like a C function declaration:

    <result_type> <func_name>(<parameters>)

    Example: "vec3 f(double x, double x, double z)"
    
    The function name may be left out

    The function returns a tuple (result_type, func_name, params)
    where params is a list of tuples (name, type, default). The name of a
    parameter might also be empty.
    """
    funcdecl = funcdecl.strip()

    # Parse the result type...
    m = re.match("[a-zA-Z0-9]+ ", funcdecl)
    if m==None:
        raise SyntaxError("No valid return type found")
        
    restype = funcdecl[m.start():m.end()].strip()

    # Parse the function name...
    s = funcdecl[m.end():]
    i = s.find("(")
    if i==-1:
        raise SyntaxError("Parameter block is missing")
    funcname = s[:i].strip()

    # Parse the arguments...
    args = s[i+1:-1].split(",")
    inputs = []
    for s in args:
        s = s.strip()
        f = s.split()
        if len(f)==1:
            type = f[0]
            name = ""
        elif len(f)==2:
            type,name = f
        else:
            raise SyntaxError("Invalid parameter declaration: '%s'"%s)
        inputs.append((name,type,None))

    return restype, funcname, inputs
    
def _inspectFuncDecl(func):
    """Determine the function declaration from the function itself.

    The function returns a tuple (result_type, func_name, params)
    where params is a list of tuples (name, type, default).

    result_type and params may also be None if they cannot be determined
    (i.e. if not every argument to the function has a default value).
    """
    # Extract function name...
    funcname = func.__name__

    args,va,vkw,defaults = inspect.getargspec(func)
    # Only proceed if every argument has a default value
    if defaults==None or len(args)!=len(defaults):
        return None, funcname, None

    # Extract the argument types (=the types of the default values)
    inputs = []
    for name, value in zip(args, defaults):
        type = value.__class__.__name__
        if type=="float":
            type = "double"
        inputs.append((name, type, value))

    # Determine the result type
    v = func()
    restype = v.__class__.__name__
    if restype=="float":
        restype = "double"
    
    return restype, funcname, inputs
        

def _defaultValueLiteral(type, default):
    """Convert a value into a string so that it can be executed again.
    """
    if default==None:
        return ""
    
    if type in ["vec3", "vec4"]:
        return "%s(%s)"%(type,default)
    else:
        return "%s"%default


def createFunctionComponentSource(clsname, restype, funcname, inputs):
    """Create a string that creates a component class.

    """

    res = """from cgkit import _core
from cgkit.globalscene import getScene
    
class %s(Component):
    
    def __init__(self):
        Component.__init__(self, func.__name__)
"""%(clsname)

    # Create output slot
    outname = "output"
    slotclassname = "Procedural%sSlot"%restype.capitalize()
    slotname = "%s_slot"%outname
    res += "        self.%s = %s(self.compute%s)\n"%(slotname, slotclassname, outname.capitalize())
    res += '        self.addSlot("%s", self.%s)\n'%(outname, slotname)

    # Create input slots
    for name, type, default in inputs:
        slotclassname = "%sSlot"%type.capitalize()
        slotname = "%s_slot"%name
        res += "        self.%s = %s(%s)\n"%(slotname, slotclassname, _defaultValueLiteral(type, default))
        res += '        self.addSlot("%s", self.%s)\n'%(name, slotname)
        res += "        self.%s.addDependent(self.%s_slot)\n"%(slotname, outname)
        if name=="time" and type=="double":
            res += "        getScene().timer().time_slot.connect(self.%s)\n"%slotname

#    res += "        self._func_obj = %s\n"%funcname
    res += "        self._func_obj = func\n"

    res += '\n    exec(slotPropertyCode("%s"))\n'%(outname)
    for name,type,default in inputs:
        res += '    exec(slotPropertyCode("%s"))\n'%(name)

    res += """\n    def compute%s(self):
        return self._func_obj("""%(outname.capitalize())

    inames = ["%s=self.%s"%(x[0],x[0]) for x in inputs]
    res += ", ".join(inames)
    res += ")\n"
    
    return res


# createFunctionComponent
def createFunctionComponent(func, funcdeclaration=None):

    if funcdeclaration==None:
        restype, funcname, inputs = _inspectFuncDecl(func)
        if restype==None:
            raise ValueError("Cannot determine the types of the arguments and the return value.")
    else:
        restype, funcname, inputs = _parseFuncDecl(funcdeclaration)

    if funcname==None or funcname=="":
        funcname = func.__name__
    
    clsname = funcname.capitalize()
    s = createFunctionComponentSource(clsname, restype, funcname, inputs)
    ns = copy.copy(globals())
    ns["func"]=func
    exec(s, ns)
    return ns[clsname]
    
