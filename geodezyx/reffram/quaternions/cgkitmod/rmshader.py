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
# $Id: rmshader.py,v 1.9 2006/05/26 21:33:29 mbaas Exp $

import sys, os, os.path, re, types, io, copy
from . import protocols
from .Interfaces import *
from . import ri
from . import material, lightsource
from . import ribexport
from . import slparams
from .slots import *
from .cgtypes import *

# RMShader
class RMShader(object):
    """RenderMan shader source file.

    This class encapsulates an external RenderMan shader source file.
    The shader is parsed in the constructor and for each parameter
    a corresponding slot is created. The parameter values themselves
    are also made accessible via normal attribute access.

    When the shader is to be instantiated, you can get the shader name
    by calling shaderName() and the parameter dictionary for the
    current time by calling params(). The type of the shader is returned
    by the method shaderType().

    The source code for the material class can be obtained by calling
    loadShaderSource().

    The class can also be used to just store a shader name. In this
    case you only pass the shader name instead of the file name to
    the constructor. 

    This class has to be used in conjunction with the RMMaterial or
    RMLightSource class.

    \see RMMaterial, RMLightSource
    """
    
    def __init__(self,
                 shader = None,
                 transform = mat4(1),
                 cpp = None,
                 cpperrstream = sys.stderr,
                 includedirs = None,
                 defines = None,
                 params = None,
                 **keyargs):
        """Constructor.

        The first argument is the shader name or file name. cpp determines
        the preprocessor that should be used when extracting parameters.
        cpperrstream is used to output errors from the preprocessor. 
        includedirs is a list of strings that contain directories where to
        look for include files. defines is a list of tuples (name, value)
        that specify the predefined symbols to use (see the function 
        slparams.slparams() for details).
        params can be used to declare parameters if the shader source
        is not available. The value must be a dictionary that contains
        token/value pairs. The token may contain an inline declaration.
        Any additional keyword argument is also considered to be a shader
        parameter. However, this parameter cannot have an inline declaration,
        so it is recommended to declare the parameter afterwards using
        the declare() method.
        
        \param name (\c str) Shader file name or shader name
        """
      
        # Shader file or None
        self.filename = None
        # Shader name (without path and extension)
        self.shadername = None
        # Shader type (surface, displacement, ...)
        self.shadertype = None
        # Shader transformation
        self.transform = transform
        # Shader parameters as dictionary.
        # Key:Parameter name   Value:Declaration
        self.shaderparams = {}

        # String parameters are not stored as slots because they can either
        # contain a string or a RenderPass object
        # str_params: Key:Parameter name /  Value:Value
        self.str_params = {}

        # Params for which there is no declaration
        # Key: Param name   Value:Value
        self.undeclared = keyargs

        if params!=None:
            for param in params:
                f = param.split()
                # Was the parameter an empty string? then ignore
                if len(f)==0:
                    continue
                # Add the param to the 'undeclared' dict
                # (even if there is a declaration this will ensure that the
                # value is taken as default value during the declaration)
                paramname = f[-1]
                if paramname not in self.undeclared:
                    self.undeclared[paramname] = params[param]
                if len(f)>1:
                    self.declare(param)
            
        
        # If there is no extension then name is just the shader name
        # and not the shader file
        if shader!=None:
            if os.path.splitext(shader)[1]=="":
                self.shadername = shader
            else:
                self.filename = shader

        # Read the shader parameters from the shader source file...
        if self.filename!=None:
            slinfo = slparams.slparams(shader, cpp=cpp, cpperrstream=cpperrstream, includedirs=includedirs, defines=defines)
            if len(slinfo)==0:
                raise ValueError("no shader found in %s"%shader)
            if len(slinfo)>1:
                print(("WARNING: There is more than one shader in %s"%shader))

            # Declare the variables...
            self.shadertype, self.shadername, params = slinfo[0]
            for p in params:
                cls = p[1]
                type = p[2]
                arraysize = p[3]
                name = p[4]
                default = p[6]
                self.declare(name, type, cls, arraysize, default)

    # getattr
    def __getattr__(self, name):
        # Shader parameter?
        slot = self.__dict__.get("%s_slot"%name, None)
        if slot!=None:
            # Array slot or normal slot?
            if isinstance(slot, IArraySlot):
                return list(slot)
            else:
                return slot.getValue()

        # String parameter?
        str_params = self.__dict__.get("str_params", {})
        if name in str_params:
            return str_params[name]
    
        raise AttributeError('shader "%s" has no attribute "%s"'%(self.shadername, name))

    # setattr
    def __setattr__(self, name, val):
        # Shader parameter?
        slot = self.__dict__.get("%s_slot"%name, None)
        if slot!=None:
            # Array slot or normal slot?
            if isinstance(slot, IArraySlot):
                for i,v in enumerate(val):
                    slot[i] = v
            else:
                slot.setValue(val)
            return

        # String parameter?
        str_params = self.__dict__.get("str_params", {})
        if name in str_params:
            str_params[name] = val
            return
            
        object.__setattr__(self, name, val)


    # shaderName
    def shaderName(self):
        """Return the shader name.

        \return Shader name
        """
        return self.shadername

    # shaderType
    def shaderType(self):
        """Return the shader type.

        None is returned if only the shader name was specified.

        \return Shader type ("surface", "light", ...)
        """
        return self.shadertype

    # params
    def params(self, passes=None):
        """Return the parameter dictionary for the current time.

        If available the parameters will contain inline declarations.

        \return Dictionary containing all parameters.
        """
        
        res = copy.copy(self.undeclared)

        for name in self.shaderparams:
            decl = self.shaderparams[name]
            p = "%s %s"%(decl, name)
            val = getattr(self, name)
            if isinstance(val, ribexport.RenderPass):
                if val.done():
                    mapname = val.realFilename(val.output[0][0])
#                    mapname = os.path.splitext(mapname)[0]+".map"
                    val = mapname
                else:
                    val = ""
            res[p] = val

        return res

    # getPasses
    def getPasses(self):
        """Return a list containing the passes required for this shader.
        """
        res = []
        for val in list(self.str_params.values()):
            if isinstance(val, ribexport.RenderPass):
                res.append(val)
        return res

    # declare
    def declare(self, name, type=None, cls=None, arraysize=None, default=None):
        """Declare a shader parameter.

        name is the parameter name. name may also contain the entire
        declaration in SL syntax. In this case, all other arguments
        are ignored, otherwise they provide the missing information.
        type is the only parameter that is mandatory if name does not
        contain the entire declaration.

        When a parameter is declared it is added to the list of known
        parameters and a corresponding slot (<name>_slot) is created.

        Examples:

        shader.declare('uniform float Ka=0.5')
        shader.declare('uniform float Ka')
        shader.declare('float Ka')
        shader.declare('Ka', type='float')
        """

#        print 'declare("%s", type=%s, cls=%s, arraysize=%s, default=%s)'%(name, type, cls, arraysize, default)
        
        # Get a slparams-stype params tuple from either name alone or from
        # all the args
        params = self._declare_getParams(name, type, cls, arraysize, default)

        typelut = {"float":"double",
                   "string":"str",
                   "color":"vec3",
                   "point":"vec3",
                   "vector":"vec3",
                   "normal":"vec3",
                   "matrix":"mat4"}

        # Iterate over all parameters (as there could be several when name
        # contains something like "float Ka; float Kd")...
        for p in params:
            ptype = p[2]
            parraylen = p[3]
            pname = p[4]
            pdefault = slparams.convertdefault(p)
            slottype = typelut[ptype]
            if parraylen is None:
                decl = "%s %s"%(p[1], ptype)
            else:
                decl = "%s %s[%d]"%(p[1], ptype, parraylen)
            # Check if this variable was already declared...
            # (it's ok if the new and old declarations are identical)
            if pname in self.shaderparams:
                if decl!=self.shaderparams[pname]:
                    raise ValueError('"%s" is already declared as "%s"'%(pname, self.shaderparams[pname]))
                continue
            # Check if the parameter was specified in the constructor.
            # If so, use the value to initialize the slot            
            if pname in self.undeclared:
                pytype = slottype
                if pytype=="double":
                    pytype = "float"

                # Don't cast strings because they might contain a RenderPass
                if pytype=="str":
                    pdefault = self.undeclared[pname]
                else:
                    # Scalar?
                    if parraylen is None:
                        pdefault = eval("%s(%s)"%(pytype, repr(self.undeclared[pname])))
                    # Array
                    else:
                        userDefault = self.undeclared[pname]
                        try:
                            if len(userDefault)!=parraylen:
                                raise ValueError('Invalid default value for shader parameter "%s". Expected %s values, but got %s'%(pname, parraylen, len(userDefault)))
                        except TypeError:
                            raise TypeError('Invalid default value for shader parameter "%s". Expected a sequence of %s values.'%(pname, parraylen))
                        pdefault = [eval("%s(%s)"%(pytype, repr(v))) for v in self.undeclared[pname]]
                del self.undeclared[pname]
                
            # Create the slot and add the variable to the params dictionary
            if ptype=="string":
                self.str_params[pname] = pdefault
            else:
                self.createSlot(pname, slottype, parraylen, pdefault)
            self.shaderparams[pname] = decl
    
    def _declare_getParams(self, name, type, cls, arraysize, default):
        """Helper method for declare().
        
        Turns the arguments into a params "tuple". See declare() for a description 
        of the arguments.
        The return value is either a params object as returned by slparams()
        or an old-style params tuple.
        """
        # Create a "dummy shader" which will be passed to slparams to parse
        # the declaration in name
        shd = "surface spam(%s) {}"%name
        try:
            # Force a syntax error when name contains no declaration
            if " " not in name:
                raise slparams.SyntaxError()
            slinfo = slparams.slparams(io.StringIO(shd))
            shdtype, shdname, params = slinfo[0]
        except slparams.SyntaxError as e:
            # Check if name is only a single name or if there was an attempt
            # to specify the entire declaration
            invalid = " []():;'\"'"
            for inv in invalid:
                if inv in name:
                    raise ValueError('Invalid declaration: "%s"'%name)
            # It's probably really just the name, so use the remaining
            # arguments to create a parameter tuple...
            if cls is None:
                cls = "uniform"
            if type is None:
                raise ValueError('No type for parameter "%s" specified'%name)
            if type not in ["float", "string", "color", "point", "vector",
                            "normal", "matrix"]:
                raise ValueError('Invalid type for parameter "%s": %s'%(name, type))
            params = [("", cls, type, arraysize, name, "", str(default))]
            
        return params

        
    # loadShaderSource
    def loadShaderSource(self):
        """Load shader source and replace the shader name.

        This method loads a shader file, replaces the shader name
        with "$SHADERNAME" and returns the source code as a string.
        None is returned if \a filename is None.

        \return Shader source code or None
        """

        if self.filename is None:
            return None
        
        f = file(self.filename)
        src = f.read()
        # Search for <shader type> + one or more white space + <shadername>..
        match = re.search("%s\s+%s"%(self.shadertype, self.shadername), src)
        if match is not None:
            s, e = match.start(), match.end()
            src = src[:e-len(self.shadername)] + "$SHADERNAME" + src[e:]
        else:
            print(('Shader name "%s" not found in %s'%(self.shadername, self.filename)))

        return src
        

    ## protected:

    # createSlot
    def createSlot(self, name, type, arraylen, default):
        """Create a slot for a shader parameter.

        \param name (\c str) Parameter name
        \param type (\c str) Slot type ("double", "vec3", ...)
        \param arraylen (\c int) Array length or None if the parameter is not an array
        \param default Default value or None
        """
        
        if arraylen is None:
            locals = {}
            exec("slot = %sSlot()"%type.capitalize(), globals(), locals)
            slot = locals["slot"]
            if default is not None:
                slot.setValue(default)
        else:
            locals = {}
            exec("slot = %sArraySlot()"%type.capitalize(), globals(), locals)
            slot = locals["slot"]
            slot.resize(arraylen)
            if default is not None:
                for i,v in enumerate(default):
                    slot[i] = v
                
        setattr(self, "%s_slot"%name, slot)

        

# RMMaterial
class RMMaterial(material.Material):
    """RenderMan material that takes shader source files as input.

    Use this material class if you want to write the RenderMan shaders
    yourself in external source files or if you want to use external shaders
    that you will compile manually.

    The material may consist of a surface shader, a displacement shader
    and an interior shader.

    The shader source files (or only the shader names) are passed via
    RMShader instances as arguments to the constructor. If the RMShader
    instance points to a file, the material object will take care of the
    compilation of the file. Otherwise, it is up to you to compile the
    shader and make sure that the renderer can find it.

    The parameters of the shaders are made available as attributes of
    the material objects. The corresponding slots can be obtained
    by adding the suffix \c _slot to the name. Attribute names in the
    surface shader have priority over the attributes in the displacement
    shader which in turn has priority over the interior shader. This means,
    if there are identical parameter names in all shaders you will access
    the parameter of the surface shader. You can also access the attributes
    of each shader via the \c surface, \c displacement and \c interior
    attributes which contain the corresponding RMShader instances.

    Example:
    
    \code
    mat = RMMaterial(surface = RMShader("mysurface.sl"),
                     displacement = RMShader("c:\\shaders\\dented.sl"),
                     displacementbound = ("current", 1.0),
                     color = (1,0.5,0.8)
                     )
    ...
    Sphere(pos=(1,2,3), radius=0.5, material=mat)
    \endcode 

    In this example, the material uses the surface shader \c mysurface.sl
    and the displacement shader \c c:\shaders\dented.sl. The shaders will
    be compiled automatically because the shader source files are given
    (instead of just the shader names).

    \see RMShader
    """

    protocols.advise(instancesProvide=[ribexport.IMaterial])

    def __init__(self,
                 name = "RMMaterial",
                 surface = None,
                 displacement = None,
                 displacementbound = ("current", 0.0),
                 interior = None,
                 color = (1,1,1),
                 opacity = (1,1,1)):
        
        material.Material.__init__(self, name, 0)

        if isinstance(surface, str):
            surface = RMShader(surface)
        if isinstance(displacement, str):
            displacement = RMShader(displacement)
        if isinstance(interior, str):
            interior = RMShader(interior)

        self.surface = surface
        self.displacement = displacement
        self.displacementbound = displacementbound
        self.interior = interior

        self._color = color
        self._opacity = opacity

        if self.surface==None:
            self.surface = RMShader()
        if self.displacement==None:
            self.displacement = RMShader()
        if self.interior==None:
            self.interior = RMShader()

    def __getattr__(self, name):
        if name[-5:]=="_slot":
            slotname = name
        else:
            slotname = "%s_slot"%name
        if hasattr(self.surface, slotname):
            return getattr(self.surface, name)
        elif hasattr(self.displacement, slotname):
            return getattr(self.displacement, name)
        elif hasattr(self.interior, slotname):
            return getattr(self.interior, name)            
        else:
            raise AttributeError('material "%s" has no attribute "%s"'%(self.name, name))


    def __setattr__(self, name, val):
        slotname = "%s_slot"%name
        srf = self.__dict__.get("surface", None)
        dsp = self.__dict__.get("displacement", None)
        itr = self.__dict__.get("itr", None)
        if hasattr(srf, slotname):
            setattr(srf, name, val)
        elif hasattr(dsp, slotname):
            setattr(dsp, name, val)
        elif hasattr(itr, slotname):
            setattr(itr, name, val)
        else:
            material.Material.__setattr__(self, name, val)


    def createPasses(self):
        """Returns a list of RenderPass objects."""
        res = self.surface.getPasses()
        res += self.displacement.getPasses()
        res += self.interior.getPasses()
        return res

    def preProcess(self, exporter):
        """Preprocessing method."""
        pass

    def color(self):
        """Return the color for the RiColor() call or None.
        """
        return self._color

    def opacity(self):
        """Return the opacity for the RiOpacity() call or None.
        """
        return self._opacity

    # surfaceShaderName
    def surfaceShaderName(self):
        """Returns the name of the corresponding surface shader or None.
        """
        return self.surface.shaderName()

    def surfaceShaderSource(self):
        return self.surface.loadShaderSource()

    def surfaceShaderParams(self, passes):
        """Return a dictionary with shader parameters and their values."""
        return self.surface.params(passes)

    def surfaceShaderTransform(self):
        return self.surface.transform

    # displacementShaderName
    def displacementShaderName(self):
        return self.displacement.shaderName()
    
    def displacementShaderSource(self):
        return self.displacement.loadShaderSource()
    
    def displacementShaderParams(self, passes):
        return self.displacement.params(passes)

    def displacementBound(self):
        return self.displacementbound

    def displacementShaderTransform(self):
        return self.displacement.transform

    # interiorShaderName
    def interiorShaderName(self):
        return self.interior.shaderName()
    
    def interiorShaderSource(self):
        return self.interior.loadShaderSource()
    
    def interiorShaderParams(self, passes):
        return self.interior.params(passes)

    def interiorShaderTransform(self):
        return self.interior.transform

# RMLightSource
class RMLightSource(lightsource.LightSource):
    """RenderMan light source.

    Use this light source class if you want to write the RenderMan light shader
    yourself in an external source file or if you want to use an external
    shader that you will compile manually.

    The shader source file (or only the shader name) is passed via a
    RMShader instances as argument to the constructor. If the RMShader
    instance points to a file, the light source will take care of the
    compilation of the file. Otherwise, it is up to you to compile the
    shader and make sure that the renderer can find it.

    The parameters of the shader are made available as attributes of
    the light source object. The corresponding slots can be obtained
    by adding the suffix \c _slot to the name. 

    \see RMShader
    """

    protocols.advise(instancesProvide=[ISceneItem, ribexport.ILightSource])

    def __init__(self,
                 name = "RMLightSource",
                 shader = None,
                 **params):
        
        lightsource.LightSource.__init__(self, name=name, **params)

        if isinstance(shader, str):
            shader = RMShader(shader)

        self.shader = shader

        if self.shader==None:
            self.shader = RMShader()

    def protocols(self):
        return [ISceneItem, ribexport.ILightSource]

    def __getattr__(self, name):
#        if name in self.__dict__:
#            return self.__dict__.get(name)
        
        if name[-5:]=="_slot":
            slotname = name
        else:
            slotname = "%s_slot"%name
        if hasattr(self.shader, slotname):
            return getattr(self.shader, name)
        else:
            raise AttributeError('light source "%s" has no attribute "%s"'%(self.name, name))


    def __setattr__(self, name, val):
        slotname = "%s_slot"%name
        shd = self.__dict__.get("shader", None)
        if hasattr(shd, slotname):
            setattr(srf, name, val)
        else:
            material.Material.__setattr__(self, name, val)


    def createPasses(self):
        """Returns a list of RenderPass objects."""
        return self.shader.getPasses()

    # shaderName
    def shaderName(self):
        """Returns the name of the corresponding shader or None.
        """
        return self.shader.shaderName()

    # surfaceShaderSource
    def shaderSource(self):
        return self.shader.loadShaderSource()

    # surfaceShaderParams
    def shaderParams(self, passes):
        """Return a dictionary with shader parameters and their values."""
        return self.shader.params(passes)


        
