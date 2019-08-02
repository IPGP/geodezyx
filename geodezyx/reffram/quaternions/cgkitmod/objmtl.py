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
# $Id: objmtl.py,v 1.3 2006/03/20 19:32:12 mbaas Exp $

from .cgtypes import *

# WavefrontReaderBase
class WavefrontReaderBase:
    """Wavefront OBJ/MTL reader base class.

    This is the base class for the OBJ and the MTL reader.
    """
    
    def __init__(self):
        """Constructor.
        """
        self.linenr = 0
        self.line = ""

    # read
    def read(self, f):
        """Read the content of a file.

        f is a file like object that can be used to read the content
        of the file.
        The file is read and for each keyword a handle_<keyword>() method
        is called with the data as argument (the number of arguments depends
        on the keyword). Each argument is of type str.
        A syntax error is generated if such a handler method is not available.
        The task of these handler methods is to preprocess the data, check
        for errors and invoke the final handler methods which are just called
        after the keyword.
        Before the file is read, the begin() method is called. At the end,
        the end() method is called.
        """

        self.linenr = 0
        self.begin()
        for self.line in f:
            self.linenr+=1
            line2 = self.line.strip()
            # Ignore empty lines and comments
            if line2=="" or line2[0] in ['#', '$', '!', '@']:
                continue

            a = line2.split()
            cmd = a[0]
            args = a[1:]
            handler = getattr(self, "handle_%s"%cmd, None)
            if handler!=None:
                handler(*args)
            else:
                self.handleUnknown(cmd, args)

        self.end()

    # begin
    def begin(self):
        """Begin reading a file.

        This method is called before the file is read.
        """
        pass

    # end
    def end(self):
        """End of reading.

        This method is called after the file was read.
        """
        pass

    # handleUnknown
    def handleUnknown(self, cmd, arglist):
        """Handle unknown keywords.

        cmd is the command keyword (the 1st argument in the current line)
        and arglist is the data (the remaining arguments).

        The default implementation raises a SyntaxError exception.
        """
        raise SyntaxError("Unknown statement in line %d: %s"%(self.linenr, self.line))


# MTLReader
class MTLReader(WavefrontReaderBase):
    """Wavefront MTL reader.

    This class can be used as base class to read a MTL file.
    """
    
    def __init__(self):
        """Constructor.
        """
        WavefrontReaderBase.__init__(self)

#    def handleUnknown(self, cmd, arglist):
#        pass
#        print "Unknown:",cmd,arglist

    # parseMapArgs
    def parseMapArgs(self, args):
        """Parse the arguments from a map_xyz line into options and names.

        args must be a list of strings containing the arguments (already
        split).
        Returns a 2-tuple (options, names) where options is a dictionary
        that contains the options that were present in args (key is
        the option string). names is a list of strings that were no option.
        """
        options = {}
        names = []
        while len(args)>0:
            a = args[0]
            # Is this an option?
            if a[0]=="-":
                if a=="-s":
                    val = tuple([float(x) for x in args[1:4]])
                    options[a] = val
                    args = args[4:]
                elif a=="-t":
                    val = tuple([float(x) for x in args[1:4]])
                    options[a] = val
                    args = args[4:]
                elif a=="-o":
                    val = tuple([float(x) for x in args[1:4]])
                    options[a] = val
                    args = args[4:]
                elif a=="-mm":
                    val = tuple([float(x) for x in args[1:3]])
                    options[a] = val
                    args = args[3:]
                elif a=="-bm":
                    options[a] = float(args[1])
                    args = args[2:]
                elif a=="-type":
                    options[a] = args[1]
                    args = args[2:]
                elif a=="-clamp" or a=="-blendu" or a=="-blendv":
                    options[a] = args[1].lower()=="on"
                    args = args[2:]
                else:
                    raise SyntaxError("Unknown map option in line %d: %s"%(self.linenr, a))
            else:  # no option
                names.append(a)
                args = args[1:]

        return options, names

    # Pre handler methods (they must be called "handle_<keyword>")

    def handle_newmtl(self, *name):
        """New material command.

        name is the name of the material (all names are concatenated).
        """
        self.newmtl(" ".join(name))

    def handle_illum(self, model):
        """Illumination model."""
        self.illum(int(model))

    def handle_Ns(self, shininess):
        """Shininess coefficient.
        """
        self.Ns(float(shininess))

    def handle_Ka(self, r, g, b):
        """Ambient color.
        """
        self.Ka(vec3(float(r), float(g), float(b)))

    def handle_Kd(self, r, g, b):
        """Diffuse color.
        """
        self.Kd(vec3(float(r), float(g), float(b)))

    def handle_Ks(self, r, g, b):
        """Specular color.
        """
        self.Ks(vec3(float(r), float(g), float(b)))

    def handle_Ke(self, r, g, b):
        """Emissive color.
        """
        self.Ke(vec3(float(r), float(g), float(b)))

    def handle_Tr(self, alpha):
        """Transparency."""
        self.Tr(float(alpha))

    def handle_d(self, alpha):
        """Dissolve factor (transparency)."""
        self.d(float(alpha))

    def handle_Tf(self, *args):
        """Transparency (as color)."""
        if len(args)==1:
            v = vec3(float(args[0]))
        else:
            r,g,b = args
            v = vec3(float(r), float(g), float(b))
        
        self.Tf(v)

    def handle_Ni(self, ref):
        """Refraction index."""
        self.Ni(float(ref))

    def handle_sharpness(self, v):
        """Sharpness value."""
        self.sharpness(float(v))

    def handle_map_Ka(self, *mapargs):
        """Ambient texture file."""
        opts, args = self.parseMapArgs(mapargs)
        if len(args)>0:
            self.map_Ka(args[0], opts)

    def handle_map_Kd(self, *mapargs):
        """Diffuse texture file."""
        opts, args = self.parseMapArgs(mapargs)
        if len(args)>0:
            self.map_Kd(args[0], opts)

    def handle_map_Ks(self, *mapargs):
        """Specular texture file."""
        opts, args = self.parseMapArgs(mapargs)
        if len(args)>0:
            self.map_Ks(args[0], opts)

    def handle_map_Ke(self, *mapargs):
        """Emission texture file."""
        opts, args = self.parseMapArgs(mapargs)
        if len(args)>0:
            self.map_Ke(args[0], opts)

    def handle_map_Ns(self, *mapargs):
        """Shininess texture file."""
        opts, args = self.parseMapArgs(mapargs)
        if len(args)>0:
            self.map_Ns(args[0], opts)

    def handle_map_d(self, *mapargs):
        """Opacity texture file."""
        opts, args = self.parseMapArgs(mapargs)
        if len(args)>0:
            self.map_d(args[0], opts)

    def handle_map_Bump(self, *mapargs):
        """Bump map texture file."""
        opts, args = self.parseMapArgs(mapargs)
        if len(args)>0:
            self.map_Bump(args[0], opts)

    def handle_refl(self, *mapargs):
        """Reflection map."""
        opts, args = self.parseMapArgs(mapargs)
        if len(args)>0:
            self.refl(args[0], opts)

    # Handler methods

    def newmtl(self, name):
        pass

    def illum(self, model):
        """Illumination model."""
        pass

    def Ns(self, shininess):
        """Shininess."""
        pass

    def Ka(self, c):
        """Ambient color.

        c is a vec3 containing the color.
        """
        pass

    def Kd(self, c):
        """Diffuse color.

        c is a vec3 containing the color.
        """
        pass

    def Ks(self, c):
        """Specular color.

        c is a vec3 containing the color.
        """
        pass

    def Ke(self, c):
        """Emissive color.

        c is a vec3 containing the color.
        """
        pass

    def Tr(self, alpha):
        """Transparency."""
        pass

    def d(self, alpha):
        """Dissolve factor (transparency)."""
        pass

    def Tf(self, c):
        """Transparency.

        c is a vec3 containing a color.
        """
        pass

    def Ni(self, ref):
        """Refraction index."""
        pass

    def sharpness(self, v):
        """Sharpness value."""
        pass

    def map_Ka(self, mapname, options):
        pass

    def map_Kd(self, mapname, options):
        pass

    def map_Ks(self, mapname, options):
        pass

    def map_Ke(self, mapname, options):
        pass

    def map_Ns(self, mapname, options):
        pass

    def map_d(self, mapname, options):
        pass

    def map_Bump(self, mapname, options):
        pass

    def refl(self, mapname, options):
        pass


# OBJReader
class OBJReader(WavefrontReaderBase):
    """Wavefront OBJ reader.

    This class can be used as base class to read an OBJ file.
    """
    
    def __init__(self):
        """Constructor.
        """
        WavefrontReaderBase.__init__(self)
        self.v_count = 0
        self.vp_count = 0
        self.vt_count = 0
        self.vn_count = 0

    # read
    def read(self, f):
        self.v_count = 0
        self.vp_count = 0
        self.vt_count = 0
        self.vn_count = 0
        WavefrontReaderBase.read(self, f)

    # Pre handler methods (they must be called "handle_<keyword>")

    def handle_mtllib(self, *files):
        """Material library command.

        files contains the filenames.
        """
        if len(files)==0:
            raise SyntaxError("No material library given in line %d: %s"%(self.linenr, self.line))
        self.mtllib(*files)

    def handle_usemtl(self, *name):
        """Material name.

        (all names are concatenated)
        """
        self.usemtl(" ".join(name))

    def handle_g(self, *groups):
        """Group command.

        groups are the group names that the following geometry belongs to.
        """
        if len(groups)==0:
            groups = ("default",)
#            raise SyntaxError("No group name given in line %d: %s"%(self.linenr, self.line))
        self.g(*groups)

    def handle_s(self, group_number):
        """Smoothing group command.

        group_number is a string that contains the smoothing group number
        or "off".
        """
        group_number = group_number.lower()
        if group_number=="off":
            group_number = 0
        try:
            group_number = int(group_number)
        except:
            raise SyntaxError("Invalid smoothing group number in line %d: %s"%(self.linenr, self.line))
        self.s(group_number)

    def handle_v(self, x, y, z, w=1):
        """Vertex definition.
        """
        self.v_count += 1
        self.v(vec4(float(x), float(y), float(z), float(w)))

    def handle_vp(self, u, v=0, w=1):
        """Vertex in parameter space.
        """
        self.vp_count += 1
        self.vp(vec3(float(u), float(v), float(w)))

    def handle_vn(self, x, y, z):
        """Normal.
        """
        self.vn_count += 1
        self.vn(vec3(float(x), float(y), float(z)))

    def handle_vt(self, u, v=0, w=0):
        """Texture vertex.
        """
        self.vt_count += 1
        self.vt(vec3(float(u), float(v), float(w)))

    def handle_p(self, *verts):
        """Points."""
        
        if len(verts)==0:
            raise SyntaxError("At least 1 vertex required in line %d: %s"%(self.linenr, self.line))

        vlist = []
        for s in verts:
            vert = int(s)
            if vert<0:
                vert = self.v_count+vert+1
            if vert==0:
                raise ValueError("0-index in line %d: %s"%(self.linenr, self.line))
            vlist.append(vert)
            
        self.p(*vlist)     

    def handle_l(self, *verts):
        """Line."""

        if len(verts)<2:
            raise SyntaxError("At least 2 vertices required in line %d: %s"%(self.linenr, self.line))

        vlist = []
        for s in verts:
            a = s.split("/")
            if len(a)==0 or len(a)>2:
                raise SyntaxError("Syntax error in line %d: %s"%(self.linenr, self.line))
            vert = int(a[0])
            if vert<0:
                vert = self.v_count+vert+1
            tvert = None
            if len(a)>1:
                if a[1]!="":
                    tvert = int(a[1])
                    if tvert<0:
                        tvert = self.vt_count+tvert+1
            if vert==0 or tvert==0:
                raise ValueError("0-index in line %d: %s"%(self.linenr, self.line))
            vlist.append((vert, tvert))
        self.l(*vlist)

    def handle_f(self, *verts):
        """Polygonal face.
        """
        if len(verts)<3:
            raise SyntaxError("At least 3 vertices required in line %d: %s"%(self.linenr, self.line))

        vlist = []
        for s in verts:
            a = s.split("/")
            if len(a)==0 or len(a)>3:
                raise SyntaxError("Syntax error in line %d: %s"%(self.linenr, self.line))
            vert = int(a[0])
            if vert<0:
                vert = self.v_count+vert+1
            tvert = None
            normal = None
            if len(a)>1:
                if a[1]!="":
                    tvert = int(a[1])
                    if tvert<0:
                        tvert = self.vt_count+tvert+1
                if len(a)>2 and a[2]!="":
                    normal = int(a[2])
                    if normal<0:
                        normal = self.vn_count+normal+1
            if vert==0 or tvert==0 or normal==0:
                raise ValueError("0-index in line %d: %s"%(self.linenr, self.line))
            vlist.append((vert, tvert, normal))
        self.f(*vlist)

    def handle_o(self, name):
        """Object name.
        """
        self.o(name)

    def handle_bevel(self, on_off):
        """Bevel interpolation on/off.
        """
        self.bevel(on_off)

    def handle_c_interp(self, on_off):
        """Color interpolation on/off.
        """
        self.c_interp(on_off)

    def handle_d_interp(self, on_off):
        """Dissolve interpolation on/off.
        """
        self.d_interp(on_off)

    def handle_lod(self, level):
        """Level of Detail.
        """
        self.lod(int(level))

    def handle_shadow_obj(self, filename):
        """Shadow object.
        """
        self.shadow_obj(filename)

    def handle_trace_obj(self, filename):
        """Shadow object.
        """
        self.trace_obj(filename)

    # Handler methods

    def call(self, filename, *args):
        pass

    def csh(self, cmd):
        pass

    def mtllib(self, *files):
        """Specification of material libraries.

        files is a sequence of file names that contain material definitions.
        """
        pass

    def usemtl(self, name):
        """Material name.

        name is a string containing the name of the material to use for
        the following elements.
        """
        pass

    def g(self, *groups):
        """Grouping statement.

        groups is a sequence of group names that the following geometry
        belongs to.
        """
        pass

    def s(self, group_number):
        """Smoothing group.

        group_number is an integer containing the smoothing group number.
        Smoothing groups should be turned off if the group number is 0.
        """
        pass

    def v(self, vert):
        """Geometric vertex.

        vert is always a vec4."""
        pass

    def vp(self, vert):
        """A point in parameter space.

        This vertex is used for free-form curves or surfaces.
        vert is always a vec3.
        """
        pass

    def vn(self, normal):
        """Normal vector.
        
        normal is a vec3.
        """
        pass

    def vt(self, tvert):
        """Texture vertex.
        
        tvert is always a vec3.
        """
        pass

    def p(self, *verts):
        """Points.

        verts is a list of vertex indices. The indices are always >0
        (negative indices are automatically converted to their
        corresponding positive indices).
        All indices are 1-based. If an index in the file was 0, an exception
        was already thrown.
        """
        print(("POINT",verts))

    def l(self, *verts):
        """Line.

        verts contains 2-tuples (vert, tvert) which contains the indices
        to the vertexa and the texture vertex. tvert may be None.
        The indices are always >0 (negative indices are automatically
        converted to their corresponding positive indices).
        All indices are 1-based. If an index in the file was 0, an exception
        was already thrown.
        """
        pass

    def f(self, *verts):
        """Polygonal face.
        
        verts contains 3-tuples (vert, tvert, normal) which contains
        the indices to the vertex, the texture vertex and the normal.
        tvert and normal may be None, otherwise the values are always >0
        (negative indices are automatically converted to their corresponding
        positive indices).
        All indices are 1-based. If an index in the file was 0, an exception
        was already thrown.
        """
        pass

    def o(self, name):
        """Optional object name.

        name is a string containing the specified name for the elements
        following this statement.
        """
        pass

    def bevel(self, on_off):
        pass

    def c_interp(self, on_off):
        pass

    def d_interp(self, on_off):
        pass

    def lod(self, level):
        """
        level is an integer.
        """
        pass

    def shadow_obj(self, filename):
        pass

    def trace_obj(self, filename):
        pass

