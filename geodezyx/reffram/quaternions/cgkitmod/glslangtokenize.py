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
# $Id: glslangtokenize.py,v 1.1 2006/03/19 20:47:33 mbaas Exp $

"""OpenGL Shading Language Tokenizer."""

import re

WHITESPACE = 0
NAME       = 1
NUMBER     = 2
STRING     = 3
NEWLINE    = 4
OPERATOR   = 5
CHARACTER  = 6
TYPE       = 7
QUALIFIER  = 8

# tokenize
def tokenize(readline, tokeater):
    """Reads a Shading Language input stream and creates tokens.

    The first parameter, *readline*, must be a callable object which
    provides the same interface as the :meth:`readline()` method of built-in
    file objects. Each call to the function should return one line of
    input as a string.

    The second parameter, *tokeater*, must also be a callable object.
    It is called with six parameters: the token type, the token
    string, a tuple (*srow*, *scol*) specifying the row and column where
    the token begins in the source, a tuple (*erow*, *ecol*) giving the
    ending position of the token, the line on which the token was
    found and the filename of the current file.

    The token type can be one of

    * ``WHITESPACE``: This is a series of blanks and/or tabs.
    * ``NAME``: A valid identifier name or keyword.
    * ``NUMBER``: An integer or float.
    * ``STRING``: A string enclosed in ``'"'``.
    * ``NEWLINE``: A newline character.
    * ``OPERATOR``: An operator such as ``'+', '-', '!', '==', '!='``, etc.
    * ``CHARACTER``: A single character that doesn't fit anything else.
    * ``TYPE``: A language type (``void``, ``bool``, ``float``, ``int``, ``vec2``,
      ``vec3``, ``vec4``, etc.)
    * ``QUALIFIER``: A type qualifier (``const``, ``attribute``, ``uniform``,
      ``varying``, ``in``, ``out``, ``inout``)

    By default, the filename argument is an empty string. It will only
    be the actual filename if you provide a preprocessed file stream
    as input (so you should first run ``cpp`` on any shader). The
    tokenizer actually expects preprocessed data as it doesn't handle
    comments.
    """

    types = ["void", "bool", "float", "int", "bvec2", "bvec3", "bvec4",
             "ivec2", "ivec3", "ivec4", "vec2", "vec3", "vec4",
             "mat2", "mat3", "mat4", "struct",
             "sampler1D", "sampler2D", "sampler3D", "samplerCube",
             "sampler1DShadow", "sampler2DShadow"]
    
    qualifiers = ["const", "attribute", "uniform", "varying", "in", "out", "inout"]
                      
    regs =  ( (WHITESPACE, re.compile(r"[ \t]+")),
              (NAME,       re.compile(r"[A-Za-z_][A-Za-z_0-9]*")),
              (NUMBER,     re.compile(r"[0-9]+(\.[0-9]+)?(E(\+|-)?[0-9]+)?")),
              (STRING,     re.compile(r"\"[^\"]*\"")),
              (OPERATOR,   re.compile(r"\(|\)|\[|\]|\+\+|\-\-|\+|-|!|\.|~|\*|/|\^|%|&|\||<<|>>|<|>|<=|>=|==|!=|&&|\|\||\^\^|\?:|=|\+=|-=|\*=|/=|%=|<<=|>>=|&=|\^=|\|=")),
              (NEWLINE,    re.compile(r"\n"))
            )

    linenr   = 0
    filename = ""
    while 1:
        # Read next line
        line = readline()
        # No more lines? then finish
        if line=="":
            break

        linenr+=1
        # Base for starting column...
        scolbase = 0

        # Process preprocessor lines...
        if line[0]=="#":
            try:
                f = line.strip().split(" ")
                linenr = int(f[1])-1
                filename = f[2][1:-1]
            except:
                pass
            continue

        s = line

        # Create tokens...
        while s!="":
            unmatched=1
            # Check all regular expressions...
            for r in regs:
                m=r[1].match(s)
                # Does it match? then the token is found
                if m!=None:
                    scol = m.start()
                    ecol = m.end()
                    tok = s[scol:ecol]
                    s   = s[ecol:]
                    typ = r[0]
                    if typ==NAME:
                        if tok in types:
                            typ = TYPE
                        elif tok in qualifiers:
                            typ = QUALIFIER
                    tokeater(typ, tok, (linenr, scolbase+scol), (linenr, scolbase+ecol), line, filename)
                    scolbase += ecol
                    unmatched=0
                    continue

            # No match? then report a single character...
            if unmatched:
                tok = s[0]
                tokeater(CHARACTER, tok, (linenr, scolbase), (linenr, scolbase+1), line, filename)
                s = s[1:]
                scolbase += 1
            
            

def _tokeater(type, s, start, end, line, filename):
    if type==WHITESPACE or type==NEWLINE:
        return
    types = { 0:"WHITESPACE", 1:"NAME", 2:"NUMBER", 3:"STRING", 4:"NEWLINE",
              5:"OPERATOR", 6:"CHARACTER", 7:"TYPE", 8:"QUALIFIER" }
    print(("%s: %-10s: %-20s %s %s"%(filename, types[type], s, start, end)))

######################################################################

if __name__=="__main__":

    f=open("test.shader")
    tokenize(f.readline, _tokeater)
    
