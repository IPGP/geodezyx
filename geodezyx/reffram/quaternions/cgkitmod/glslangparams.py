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
# $Id: glslangparams.py,v 1.2 2006/04/12 11:52:09 mbaas Exp $

"""Extract shader parameters from an OpenGL shader source file.
"""

import sys, copy
import io
from . import glslangtokenize
from . import simplecpp
from . import slparams
from .glslangtokenize import WHITESPACE, NAME, NUMBER, STRING, NEWLINE, OPERATOR, CHARACTER, TYPE, QUALIFIER

class GLSLangParseError(Exception):
    pass

# _GLSLangParser
class _GLSLangParser:
    """Extract variables from a glslang shader.
    """
    
    def __init__(self, file, structs=None):
        """Constructor.

        fils is a file-like object.
        structs is a dictionary containing the 'predefined' structs.
        This is used internally when the contents of a struct is
        recursively parsed.
        """

        # The following lists receive the result:

        # Contains 2-tuples (type, identifier)
        self.attribute = []
        # Contains 3-tuples (type, identifier, arraysize)
        self.varying = []
        # Contains 5-tuples (type, identifier, arraysize, structname, struct)
        self.uniform = []
        # Contains 5-tuples (type, identifier, arraysize, structname, struct)
        self.const = []
        # Contains 5-tuples (type, identifier, arraysize, structname, struct)
        self.other = []

        # file
        self.file = file
        # Current state
        self.state = self.initialState
        
        # Current qualifier
        self.qualifier = None
        # Current type
        self.type = None
        # Identifier
        self.identifier = None
        # String with the array size
        self.arraysize = None
        # Struct name
        self.structname = None
        self.struct = None

        self.named_structs = {}
        if structs!=None:
            self.named_structs = copy.deepcopy(structs)

        # Number of open parentheses
        self.openparen = 0

        # The contents of a struct as a string
        # (for calling the parser again on this)
        self.structcontents = ""
        

    def read(self):
        """Read the file.
        """
        glslangtokenize.tokenize(self.file.readline, self.tokeater)

    def tokeater(self, type, s, start, end, line, filename):
        """Token eater.

        This method skips whitespace and newlines and calls the
        current state method.
        """
        if type in [WHITESPACE, NEWLINE]:
            return

        self.state(type, s, start, end, line, filename)

    def switchState(self, state):
        """Switch to a different state.

        state is a callable.
        """
        self.state = state

    def variable(self):
        """Store a result variable.

        This method is called by the states whenever all information for
        a variable has been read.
        """
        if self.qualifier=="attribute":
            self.attribute.append((self.type, self.identifier))
        elif self.qualifier=="varying":
            self.varying.append((self.type, self.identifier, self.arraysize))
        elif self.qualifier=="uniform":
            self.uniform.append((self.type, self.identifier, self.arraysize, self.structname, self.struct))
        elif self.qualifier=="const":
            self.const.append((self.type, self.identifier, self.arraysize, self.structname, self.struct))
        else:
            self.other.append((self.type, self.identifier, self.arraysize, self.structname, self.struct))

    # States:
    # Each state takes the same arguments as a token reader

    def initialState(self, type, s, start, end, line, filename):
        """Initial state.
        """
        self.qualifier = None
        self.type = None
        self.identifier = None
        self.arraysize = None
        self.structname = None
        self.struct = None
        if type==QUALIFIER:
            self.qualifier = s
            self.switchState(self.qualifierState)
        elif type==TYPE and s!="struct":
            self.type = s
            self.switchState(self.typeState)
        elif s=="struct":
            self.type = s
            self.switchState(self.structState)
        else:
            # Is the type a previously defined struct?
            if s in self.named_structs:
                self.type = "struct"
                self.structname = s
                self.struct = self.named_structs[s]
                self.switchState(self.typeState)
            else:
                raise GLSLangParseError("%s, line %d: Syntax error: %s"%(filename, start[0], s))

    def qualifierState(self, type, s, start, end, line, filename):
        """A qualifier has been read.
        """
        if type==TYPE and s!="struct":
            self.type = s
            if self.qualifier in ["attribute", "varying"] and s not in ["float", "vec2", "vec3", "vec4", "mat2", "mat3", "mat4"]:
                raise GLSLangParseError("%s, line %d: Invalid type for an %s variable: %s"%(filename, start[0], self.qualifier, s))
            self.switchState(self.typeState)
        elif s=="struct":
            self.type = s
            self.switchState(self.structState)
        else:
            # Is the type a previously defined struct?
            if s in self.named_structs:
                self.type = "struct"
                self.structname = s
                self.struct = self.named_structs[s]
                if self.qualifier in ["attribute", "varying"]:
                    raise GLSLangParseError("%s, line %d: %s variables cannot be declared as structs"%(filename, start[0], self.qualifier))
                self.switchState(self.typeState)
            else:
                raise GLSLangParseError("%s, line %d: Syntax error: %s"%(filename, start[0], s))

    def typeState(self, type, s, start, end, line, filename):
        """The type of a varibale has been read.
        """
        if type==NAME:
            self.identifier = s
            self.switchState(self.nameState)
        elif s==";":
            self.switchState(self.initialState)
        else:
            raise GLSLangParseError("%s, line %d: Syntax error: %s"%(filename, start[0], s))

    def nameState(self, type, s, start, end, line, filename):
        """The name of a variable/function has been read.
        """
        if s==",":
            self.variable()
            self.switchState(self.typeState)
        elif s==";":
            self.variable()
            self.switchState(self.initialState)
        elif s=="[":
            self.arraysize = ""
            if self.qualifier=="attribute":
                raise GLSLangParseError("%s, line %d: attribute variables cannot be declared as arrays"%(filename, start[0]))
            self.switchState(self.arrayState)
        elif s=="(":
            self.openparen = 1
            self.switchState(self.functionState)
        elif s=="=":
            self.switchState(self.initState)
        else:
            raise GLSLangParseError("%s, line %d: Syntax error: %s"%(filename, start[0], s))

    def initState(self, type, s, start, end, line, filename):
        """A '=' has been encountered.

        Skip the initializer.
        """
        if s==";":
            self.variable()
            self.switchState(self.initialState)

    def arrayState(self, type, s, start, end, line, filename):
        """A '[' has been encountered.
        """
        if s=="]":
            self.switchState(self.arrayState2)
        else:
            self.arraysize += s

    def arrayState2(self, type, s, start, end, line, filename):
        if s==";":
            self.variable()
            self.switchState(self.initialState)
        elif s==",":
            self.variable()
            self.arraysize = None
            self.switchState(self.typeState)
        else:
            raise GLSLangParseError("%s, line %d: Syntax error: %s"%(filename, start[0], s))

    def functionState(self, type, s, start, end, line, filename):
        """Skip function args.
        """
        if s=="(":
            self.openparen += 1
        elif s==")":
            self.openparen -= 1
            if self.openparen==0:
                self.switchState(self.functionState2)
            
    def functionState2(self, type, s, start, end, line, filename):
        """Skip function body.
        """

        # If this was only a function prototype there will be an immediate
        # semicolon instead of the function body
        if self.openparen==0 and s==";":
            self.switchState(self.initialState)
        elif s=="{":
            self.openparen += 1
        elif s=="}":
            if self.openparen==0:
                raise GLSLangParseError("%s, line %d: '{' expected, got '}'"%(filename, start[0]))
            self.openparen -= 1
            if self.openparen==0:
                self.switchState(self.initialState)
        elif self.openparen==0:
            raise GLSLangParseError("%s, line %d: Syntax error: %s"%(filename, start[0], s))

    def structState(self, type, s, start, end, line, filename):
        """The keyword 'struct' has been read.
        """
        self.structcontents = ""
        if type==NAME:
            self.structname = s
            self.switchState(self.structState2)
        elif s=="{":
            self.switchState(self.structState3)
        else:
            raise GLSLangParseError("%s, line %d: Syntax error: %s"%(filename, start[0], s))

    def structState2(self, type, s, start, end, line, filename):
        """The struct name has been read.
        """
        if s=="{":
            self.switchState(self.structState3)
        else:
            raise GLSLangParseError("%s, line %d: Syntax error: %s"%(filename, start[0], s))

    def structState3(self, type, s, start, end, line, filename):
        """Collect the struct contents.
        """
        if s=="}":
            f = io.StringIO(self.structcontents)
            p = _GLSLangParser(f, structs=self.named_structs)
            p.read()
            self.struct = p.other
            if self.structname!="":
                self.named_structs[self.structname] = self.struct
            self.switchState(self.typeState)
        else:
            self.structcontents += s+" "
        


# glslangparams
def glslangparams(shader=None, cpp=None, cpperrstream=sys.stderr):
    """Extracts the shader parameters from an OpenGL 2 shader source file.

    The argument *shader* is either the name of the shader source file
    or a file-like object that provides the shader sources.
    *cpp* determines how the shader source is preprocessed. It
    can either be a string containing the name of an external
    preprocessor tool (such as ``cpp``) that must take the file name as
    parameter and dump the preprocessed output to stdout or it can be
    a callable that takes *shader* and *cpperrstream* as input and returns
    the preprocessed sources as a string. If the external
    preprocessor does not produce any data a :exc:`PreprocessorNotFound`
    exception is thrown.
    The error stream of the preprocessor is written to the object
    that is specified by *cpperrstream* which must have a :meth:`write()`
    method. If *cpperrstream* is ``None``, the error stream is ignored.

    If *cpp* is ``None`` a simple internal preprocessor based on the
    :mod:`simplecpp` module is used.

    The function returns three lists (*uniform*, *attribute*, *varying*)
    that contain the variables with the corresponding qualifier.

    A uniform variable is a 5-tuple (*type*, *identifier*, *arraysize*,
    *structname*, *struct*). *arraysize* is a string containing the
    expression for the length of the array (i.e. the value between
    the square brackets). If the variable is no array, *arraysize* is ``None``.
    When the variable is a struct, *type* has the value ``'struct'``. In this
    case, the struct is given in *struct* (which is itself a list of
    variables as 5-tuples). If the struct has a name, this name is
    given in *structname*, otherwise *structname* is ``None``.

    An attribute variable is a 2-tuple (*type*, *identifier*) and a
    varying variable is a 3-tuple (*type*, *identifier*, *arraysize*) where
    *arraysize* is defined as in the uniform case.
    """

    # Run the preprocessor on the input file...
    
    if cpp==None:
        cpp = simplecpp.PreProcessor()
        
    glslangsrc = slparams.preprocess(cpp, shader, cpperrstream)
    f = io.StringIO(glslangsrc)

    # Extract the variables...
    parser = _GLSLangParser(f)
    parser.read()
    return parser.uniform, parser.attribute, parser.varying


