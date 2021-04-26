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
# -------------------------------------------------------------
# The RenderMan (R) Interface Procedures and Protocol are:
# Copyright 1988, 1989, 2000, Pixar
# All Rights Reserved
#
# RenderMan (R) is a registered trademark of Pixar
# -------------------------------------------------------------
# $Id: slparams.py,v 1.4 2006/02/14 19:29:39 mbaas Exp $

"""Extract shader parameters from a RenderMan shader source file.

Functions:

- slparams()
- convertdefault()
"""

import os, os.path, sys, string, io, types, math
from . import sltokenize
from . import cgtypes, sl, simplecpp
from . import _slparser
from ._slreturntypes import _ShaderInfo, _ShaderParam
import collections
try:
    from . import sloargs
    _has_sloargs = True
except ImportError as exc:
    _has_sloargs = False
    _sloargs_importerror = exc

class SLParamsError(Exception):
    pass

class PreprocessorNotFound(SLParamsError):
    pass

class SyntaxError(SLParamsError):
    pass

class NoMoreTokens(SLParamsError):
    pass


# Parser class (subclassed from the generated Yapps parser)
class _SLParser(_slparser._SLParserBase):
    """SL parser class.

    This class is derived from the generated parser and implements
    some missing methods.
    """
    def __init__(self, scanner):
        _slparser._SLParserBase.__init__(self, scanner)

        # Current filename
        self.filename = "?"
        # Offset which has to be subtracted from the line number to
        # get the true line number within the file.
        self.linenroffset = 0

        # Parameter attributes...
        self.output   = ""
        self.detail   = ""
        self.type     = ""
        self.arraylen = None
        self.name     = ""
        self.space    = None
        self.default  = ""
        self.spaces   = []

        # Shader parameters
        self.params = []

    def newParams(self):
        """Start a new shader.

        The parameter list is cleared.
        """
        self.params = []

    def newType(self):
        """Clear all type parameters.

        This is called when a new type is declared (which is not equivalent
        to a new parameter, e.g. "float a,b,c")
        """
        self.output   = ""
        self.detail   = "uniform"
        self.type     = ""
        self.arraylen = None
        self.name     = ""
        self.space    = None
        self.default  = ""
        self.spaces   = []

    def defaultSpace(self, typ):
        """Return the default space for typ."""
        if typ in ["point", "vector", "normal", "matrix"]:
            return "current"
        elif typ=="color":
            return "rgb"
        else:
            return None
        
    def appendSpace(self):
        """Append self.space to self.spaces."""
        if self.space!=None:
            self.spaces.append(self.space)
        else:
            self.spaces.append(self.defaultSpace(self.type))
                
        self.space = None

    def storeParam(self):
        """Store the current set of attributes as a new parameter.

        The attributes are reset so that a new parameter can begin.
        """
        if self.arraylen==None:
            if self.space==None:
                self.space = self.defaultSpace(self.type)
            self.params.append(_ShaderParam(self.output, self.detail, self.type, None,
                            self.name, self.space, self.default))
        else:
            spaces = self.spaces
            if self.defaultSpace(self.type)==None:
                spaces = None
            self.params.append(_ShaderParam(self.output, self.detail, self.type,
                               self.arraylen, self.name, spaces, self.default))
        self.arraylen = None
        self.name = ""
        self.space = None
        self.default = ""
        self.spaces   = []
        

    def switchFile(self, cppline):
        """Switch to another file.

        This method is called when a preprocessor line (starting with #)
        is read. This line contains the current file name and the line number.
        """
        f = cppline.strip().split(" ")
        linenr = int(f[1])
        filename = f[2][1:-1]
        self.filename = filename
        self.linenroffset = self._scanner.get_line_number()-linenr+1

        
# _SLfilter
class _SLfilter:
    """Used by the sltokenizer to filter the Shading Language source.

    Only the shader and function definitions remain, the bodies will
    be dropped. The filtered result will be used as input for the
    actual parser.
    """
    def __init__(self):
        # Current {}-depth (0 = outside any {..})
        self.curly_depth = 0
        # Current ()-depth
        self.bracket_depth = 0
        # Receives the source code that'll get passed to the parser
        self.SLsource = ""
        self.stream_enabled = True

        self.current_filename = None
        
    def eater(self, type, tok, start, end, line, filename):
        """Record only tokens with depth 0."""
        if filename!=self.current_filename:
            self.current_filename = filename
            if self.SLsource!="" and self.SLsource[-1]!="\n":
                self.SLsource+="\n"
            self.SLsource+='# %d "%s"\n'%(start[0],filename)
        
        if tok=="}":
            self.curly_depth-=1
            if self.curly_depth==0 and self.bracket_depth==0:
                self.stream_enabled = True
        elif tok==")":
            self.bracket_depth-=1

        # Always record newlines so that line numbers won't get messed up
        if self.stream_enabled or tok=="\n":
            self.SLsource+=tok
            
        if tok=="{":
            if self.curly_depth==0 and self.bracket_depth==0:
                self.stream_enabled = False
            self.curly_depth+=1
        elif tok=="(":
            self.bracket_depth+=1


# slparams
def slparams(slfile=None, cpp=None, cpperrstream=sys.stderr, slname=None, includedirs=None, defines=None):
    """Extracts the shader parameters from a RenderMan Shader file.

    The argument *slfile* is either the name of a compiled shader, the name of
    the shader source file (``*.sl``) or a file-like object that provides the
    shader sources.
    
    *cpp* determines how the shader source is preprocessed.
    It can either be a string containing the name of an external
    preprocessor tool (such as ``cpp``) that must take the file name as
    parameter and dump the preprocessed output to stdout or it can be
    a callable that takes *slfile* and *cpperrstream* as input and returns
    the preprocessed sources as a string. If the external
    preprocessor does not produce any data a :exc:`PreprocessorNotFound`
    exception is thrown.
    The error stream of the preprocessor is written to the object
    specified by *cpperrstream* which must have a :meth:`write()`
    method. If *cpperrstream* is ``None``, the error stream is ignored.

    If *cpp* is ``None`` a simple internal preprocessor based on the
    :mod:`simplecpp` module is used.
    The *slname* argument is an alias for *slfile*, it is only available
    for backwards compatibility.
    
    *includedirs* is a list of strings that contain directories where to
    look for include files. *defines* is a list of tuples (*name*, *value*)
    that specify the predefined symbols to use.
    
    The function returns a list of shader info objects. These objects have
    four attributes:
    
    - ``type``: The type of the shader (surface, displacement, etc.)
    - ``name``: The name of the shader
    - ``params``: The shader parameters (see below)
    - ``meta``: The shader meta data. How exactly meta data is specified depends
      on the renderer you are using.
     
    The parameters are given as a list of shader parameter objects
    describing each parameter. A shader parameter object has the
    following attributes:

    - ``outputSpec``: The output specifier (either ``"output"`` or an empty string)
    - ``storage``: The storage class (``"uniform"`` or ``"varying"``)
    - ``type``: The parameter type
    - ``size``: The array length or ``None`` if the parameter is not an array
    - ``name``: The name of the parameter
    - ``space``: The space in which a point-like type was defined
    - ``default``: The "raw" default value. If the input was a shader source file,
      this will always be a string containing an expression. If the input was
      a compiled shader this will already be an appropriate Python value.
      You should never use this value directly, but always use :func:`convertdefault()`
      to obtain a value which can be further processed. This way, your code
      will work for both, compiled shaders and shader source files.
    
    For backwards compatibility, the shader info object behaves like a
    3-tuple (*type*, *name*, *params*). The meta data can only be accessed via
    name though. The shader parameter objects can also be used like 7-tuples
    containing the above data (in the order given above).
    
    Example (output slightly reformatted for better readability)::

      >>> from cgkit import slparams
      >>> shaders = lparams.slparams("plastic.sl")
      >>> print shaders
      [('surface', 'plastic', 
        [('', 'uniform', 'float', None, 'Ka', None, '1'),
         ('', 'uniform', 'float', None, 'Kd', None, '0.5'),
         ('', 'uniform', 'float', None, 'Ks', None, '0.5'),
         ('', 'uniform', 'float', None, 'roughness', None, '0.1'),
         ('', 'uniform', 'color', None, 'specularcolor', 'rgb', '1')])]
      >>> shaders[0].type
      'surface'
      >>> shaders[0].name
      'plastic'
      >>> for param in shaders[0].params: print param.name
      ... 
      Ka
      Kd
      Ks
      roughness
      specularcolor
      >>> shaders[0].meta
      {}

    The parser used inside this function was generated using the parser generator
    `Yapps <http://theory.stanford.edu/~amitp/Yapps/>`_  by Amit Patel.
    """

    if slname is not None:
        slfile = slname
    if slfile is None:
        return []
    
    # Check if the input file is a string referring to a compiled shader
    # (suffix != .sl). If so, use the sloargs module to get the shader information
    if isinstance(slfile, str):
        if os.path.splitext(slfile)[1].lower()!=".sl":
            if (_has_sloargs):
                return sloargs.slparams(slfile)
            else:
                raise _sloargs_importerror

    # Run the preprocessor on the input file...
    
    slsrc = preprocess(cpp, slfile, cpperrstream=cpperrstream, defines=defines, includedirs=includedirs)
    f = io.StringIO(slsrc)
    
    # ...and filter it, so that only the shader and function
    # definitions remain...
    filter = _SLfilter()
    sltokenize.tokenize(f.readline, filter.eater)
    f.close()

#    print filter.SLsource

    # Parse the filtered source code...
    scanner = _slparser._SLParserBaseScanner(filter.SLsource)
    parser = _SLParser(scanner)
    
#    return wrap_error_reporter(parser, "definitions")

    try:
        lst = getattr(parser, "definitions")()
        # Turn the 3-tuples into _ShaderInfo objects
        return [_ShaderInfo(*tup) for tup in lst]
    except _slparser.NoMoreTokens as err:
        raise NoMoreTokens("No more tokens")
    except _slparser.SyntaxError as err:
        scanner = parser._scanner
        input = scanner.input
        cpplineno = scanner.get_line_number()
        lineno = cpplineno-parser.linenroffset
        # '+1' at the end to exclude the newline. If no newline is found
        # rfind() returns -1, so by adding +1 we get 0 which is what we want.
        start = input.rfind("\n", 0, err.charpos)+1
        end   = input.find("\n", err.charpos, -1)
        origline = input[start:end].replace("\t", " ")
        line = " "+origline.lstrip()
        errpos = err.charpos-start-(len(origline)-len(line))
#        print "%s^"%((errpos)*" ")
        msg = 'Syntax error in "%s", line %d:\n'%(parser.filename, lineno)
        msg += '\n%s\n%s^'%(line, errpos*" ")
        exc = SyntaxError(msg)
        exc.filename = parser.filename
        exc.lineno = lineno
        exc.line = line
        exc.errpos = errpos
        raise exc

# preprocess
def preprocess(cpp, file, cpperrstream=sys.stderr, defines=None, includedirs=None):
    """Preprocess an input file.

    cpp is either a string containing the name of an external preprocessor
    command (e.g. 'cpp'), a callable function that takes file and
    cpperrstream as input and returns the preprocessed source code as
    a string or None in which case the built-in preprocessor is used.
    If the external preprocessor produces no data it is assumed that
    it was not found a a PreprocessorNotFound exception is thrown.

    file is either a string containing the input file name or it's a
    file-like object.
    If the input file doesn't exist, an IOError is thrown.

    cpperrstream is a stream that receives any error messages from
    the preprocessor. If it's None, error messages are ignored.

    defines is a list of tuples (name, value) that specify the predefined
    symbols. includedirs is a list of strings with include paths.

    The function returns the preprocessed sources as a string.
    """

    if cpp==None:
        cpp = simplecpp.PreProcessor(defines=defines, includedirs=includedirs)
        return cpp(file, cpperrstream)
    # Is cpp a callable then invoke it and return the result
    elif isinstance(cpp, collections.Callable):
        return cpp(file, cpperrstream)

    # no callable, so cpp contains the name of an external preprocessor tool

    # Variables:
    #
    # cmd: The command that is used to execute the preprocessor.
    #      It either contains only the cpp name or the cpp name + input file
    # slsrc: (string) Unprocessed SL source if it should be passed via stdin
    #      otherwise it's None.
    # ppslsrc: Preprocessed SL source (this is the result).

    # Is file a string? Then it's a file name...
    if isinstance(file, str):
        cmdtoks = [cpp]
        if defines!=None:
            for name,value in defines:
                if value==None:
                    cmdtoks.append("-D%s"%name)
                else:
                    cmdtoks.append("-D%s=%s"%(name,value))
        if includedirs!=None:
            cmdtoks.extend(["-I%s"%dir for dir in includedirs])
        cmdtoks.append(file)
        cmd = " ".join(cmdtoks)
        print (cmd)
        slsrc = None
        # Check if the file exists and can be accessed (by trying to open it)
        # If the file doesn't exist, an exception is thrown.
        dummy = open(file)
        dummy.close()
    # ...otherwise it's a file-like object
    else:
        cmd = "%s"%cpp
        slsrc = file.read()

    # The preprocessed data will be in slsrc.
    stdin, stdout, stderr = os.popen3(cmd)
    if slsrc!=None:
        stdin.write(slsrc)
    stdin.close()
    ppslsrc = stdout.read()
    stdout.close()
    if cpperrstream!=None:
        errs = stderr.read()
        cpperrstream.write(errs)
    # No data? Then it's assumed that the preprocessor couldn't be found
    if len(ppslsrc)==0:
        raise PreprocessorNotFound("Calling '%s' didn't produce any data."%cmd)

    return ppslsrc


# Setup local namespace for convertdefault()
_local_namespace = {}
exec ("from cgkit.sl import *", _local_namespace)
   
def convertdefault(paramtuple):
    """Converts the default value of a shader parameter into a Python type.

    *paramtuple* must be a 7-tuple (or parameter object) as returned by 
    :func:`slparams()`.
    The function returns a Python object that corresponds to the default
    value of the parameter. If the default value can't be converted
    then ``None`` is returned. Only the functions that are present in the
    :mod:`sl<cgkit.sl>` module are evaluated. If a default value calls a user defined
    function then ``None`` is returned.

    The SL types will be converted into the following Python types:

        +------------+-------------+
        | SL type    | Python type |
        +============+=============+
        | ``float``  | ``float``   |
        +------------+-------------+
        | ``string`` | ``string``  |
        +------------+-------------+
        | ``color``  | ``vec3``    |
        +------------+-------------+
        | ``point``  | ``vec3``    |
        +------------+-------------+
        | ``vector`` | ``vec3``    |
        +------------+-------------+
        | ``normal`` | ``vec3``    |
        +------------+-------------+
        | ``matrix`` | ``mat4``    |
        +------------+-------------+

    Arrays will be converted into lists of the corresponding type.
    """
    global _local_namespace
    
    typ = paramtuple[2]
    arraylen = paramtuple[3]
    defstr = paramtuple[6]
    
    # The default value is not a string? Then it already contains the
    # converted default value (this is the case when the value was
    # extracted from a compiled shader).
    if not isinstance(defstr, str):
        # Make sure that point-like types are returned as vec3 and matrix types
        # are returned as mat4.
        if typ in ["color", "point", "vector", "normal"]:
            retType = cgtypes.vec3
        elif typ=="matrix":
            retType = cgtypes.mat4
        else:
            # No vec3/mat4 type, then just return the value
            return defstr
        # Cast the value...
        if arraylen is None:
            return retType(defstr)
        else:
            return [retType(v) for v in defstr]

    # Replace {} with [] so that SL arrays look like Python lists
    if arraylen is not None:
        defstr = defstr.replace("{","[").replace("}","]")
    # If the parameter is not an array, then create an array with one
    # element (to unify further processing). It will be unwrapped in the end
    if arraylen is None:
        defstr = "[%s]"%defstr
    # Evaluate the string to create "raw" Python types (lists and tuples)
    try:
        rawres = eval(defstr, globals(), _local_namespace)
    except:
        return None

    # Convert into the appropriate type...
    if typ=="float":
        try:
            res = [float(x) for x in rawres]
        except:
            return None
    elif typ=="color" or typ=="point" or typ=="vector" or typ=="normal":
        try:
            res = [cgtypes.vec3(x) for x in rawres]
        except:
            return None
    elif typ=="matrix":
        try:
            res = [cgtypes.mat4(x) for x in rawres]
        except:
            return None
    elif typ=="string":
        try:
            res = [str(x) for x in rawres]
        except:
            return None

    if arraylen is None:
        if len(res)==0:
            return None
        else:
            res = res[0]

    return res

######################################################################

if __name__=="__main__":
    pass
