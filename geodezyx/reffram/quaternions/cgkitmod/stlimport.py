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
# $Id: stlimport.py,v 1.1 2005/01/09 20:13:06 mbaas Exp $

import os.path, sys, struct
from .cgtypes import *
from .trimesh import TriMesh
from . import pluginmanager

# STLReader
class STLReader:
    
    def __init__(self, filename):
        self.filename = filename

    # read
    def read(self):
        if self.isASCII():
            self.read_ascii()
        else:
            self.read_bin()

    # read_bin
    def read_bin(self):
        """Read a binary STL file.
        """
        f = file(self.filename, "rb")

        s = f.read(80)
        objname = s.split("\000")[0]
        s = f.read(4)
        numfaces = struct.unpack("<i", s)[0]

        self.begin(objname)

        for i in range(numfaces):
            s = f.read(50)
            normal = vec3(struct.unpack("<fff", s[:12]))
            v = struct.unpack("<fffffffff", s[12:48])
            verts = [vec3(v[0:3]), vec3(v[3:6]), vec3(v[6:9])]
            self.triangle(normal, verts)

        self.end(objname)

    # read_ascii
    def read_ascii(self):
        """Read a ASCII STL file.
        """
        
        f = file(self.filename)
        
        objname = "unnamed"
        normal = vec3()
        verts = []

        state = 0
        linenr = 0
        for s in f:
            linenr+=1
            a = s.split()
            if len(a)==0:
                continue
            
            # Begin of file
            if state==0:
                if a[0]=="solid":
                    if len(a)>1:
                        objname = a[1]
                    self.begin(objname)
                    state = 1
                else:
                    raise SyntaxError('Keyword "solid" expected in line %d'%linenr)
            # Begin of a triangle ("facet") or end of file ("endsolid")
            elif state==1:
                if a[0]=="facet":
                    if len(a)!=5 or a[1]!="normal":
                        raise SyntaxError('Syntax error in line %d'%linenr)
                    normal = vec3(float(a[2]), float(a[3]), float(a[4]))
                    verts = []
                    state = 2
                elif a[0]=="endsolid":
                    self.end(objname)
                    state = 0
                else:
                    raise SyntaxError('Syntax error in line %d. Facet or end of solid expected.'%linenr)
            # "outer loop"
            elif state==2:
                if len(a)!=2 or a[0]!="outer" or a[1]!="loop":
                    raise SyntaxError('Keyword "outer loop" expected in line %d.'%linenr)
                state = 3
            # A vertex definition
            elif state==3:
                if len(a)!=4 or a[0]!="vertex":
                    raise SyntaxError('Syntax error in line %d: Vertex expected'%linenr)
                v = vec3(float(a[1]), float(a[2]), float(a[3]))
                verts.append(v)
                if len(verts)==3:
                    state = 4
            # "endloop"
            elif state==4:
                if a[0]!="endloop":
                    raise SyntaxError('Keyword "endloop" expected in line %d.'%linenr)
                state = 5
            # "endfacet"
            elif state==5:
                if a[0]!="endfacet":
                    raise SyntaxError('Keyword "endfacet" expected in line %d.'%linenr)
                self.triangle(normal, verts)
                state = 1

    # isASCII
    def isASCII(self):
        f = file(self.filename)
        s1 = f.readline().strip()
        s2 = f.readline().strip()
        f.close()
        if s1[:5]!="solid":
            return False
        if s2[:5]!="facet":
            return False
        return True

    # begin
    def begin(self, name):
        """Solid begin callback.

        This callback is called before the triangles are read.
        """
        pass

    # end
    def end(self, name):
        """Solid begin callback.

        This callback is called after the triangles were read.
        """
        pass

    # triangle
    def triangle(self, normal, verts):
        """Triangle callback.

        normal is the normal vector as vec3 and verts is a list with the
        three vertices (vec3s).
        """
        pass


# STLImport
class STLImport(STLReader):
    
    def __init__(self, filename):
        STLReader.__init__(self, filename)
        self.verts = []
        self.numfaces = 0

    # begin
    def begin(self, name):
        pass

    # end
    def end(self, name):
        faces = []
        for i in range(self.numfaces):
            faces.append(list(range(i*3, i*3+3)))
        TriMesh(name=name, verts=self.verts, faces=faces)

    # triangle
    def triangle(self, normal, verts):
        """Triangle callback.

        normal is the normal vector as vec3 and verts is a list with the
        three vertices (vec3s).
        """
        self.numfaces += 1
        self.verts += verts

# STLImporter
class STLImporter:

    _protocols = ["Import"]

    # extension
    def extension():
        """Return the file extensions for this format."""
        return ["stl"]
    extension = staticmethod(extension)

    # description
    def description(self):
        """Return a short description for the file dialog."""
        return "StereoLithography"
    description = staticmethod(description)

    # importFile
    def importFile(self, filename):
        """Import a STL file."""

        reader = STLImport(filename)
        reader.read()


######################################################################

# Register the Importer class as a plugin class
pluginmanager.register(STLImporter)
