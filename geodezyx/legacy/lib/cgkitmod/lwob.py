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
# $Id: lwob.py,v 1.1 2006/03/12 22:41:42 mbaas Exp $

import struct
from .cgtypes import vec3

class LWOBError(Exception):
    pass

# Texture
class Texture:
    """Texture description.

    A value remains None if it wasn't set in the LWOB file.
    """
    def __init__(self):
        self.type = None        # str
        self.flags = None       # int
        self.size = None        # vec3
        self.center = None      # vec3
        self.falloff = None     # vec3
        self.velocity = None    # vec3
        self.color = None       # 3-tuple
        self.value = None       # int
        self.bumpamplitude = None # float
        self.floatparams = 10*[None]  # list of floats
        self.intparams = 10*[None]    # list of ints
        self.imagename = None   # str
        self.alphaname = None   # str
        self.widthwrap = None   # int
        self.heightwrap = None  # int
        self.antialiasingstrength = None # float
        self.opacity = None     # float
        self.shader = []        # list of str
        self.shaderdata = []    # list of str
        self.seq_offset = None  # int
        self.seq_flags = None   # int
        self.seq_looplength = None # int
        self.flyer_begin = None # int
        self.flyer_end = None   # int
        self.cycle_speed = None # int
        self.cycle_low = None   # int
        self.cycle_high = None  # int

# Surface
class Surface:
    """Surface description.

    A value remains None if it wasn't set in the LWOB file.    
    """
    
    def __init__(self, name):
        self.name = name
        
        self.color = None          # 3-tuple
        self.flags = None          # int
        self.luminosity = None     # float
        self.diffuse = None        # float
        self.specular = None       # float
        self.reflection = None     # float
        self.transparency = None   # float
        self.glossiness = None     # int
        self.reflectionmode = None # int
        self.refmap = None         # str
        self.refmap_seamangle = None # float
        self.refractiveindex = None  # float
        self.edgetransp = None     # float
        self.maxsmoothangle = None # float

        self.tex_color = []        # list of Texture objects
        self.tex_diffuse = []
        self.tex_specular = []
        self.tex_reflection = []
        self.tex_transparency = []
        self.tex_luminosity = []
        self.tex_bump = []
        

# LWOBReader
class LWOBReader:
    """Lightwave Object reader.

    This class reads Lightwave Object (*.lwo) files and calls handler
    methods for each chunk found. The following handler methods can be
    implemented in derived classes:

    - handlePNTS(points):  points is a list of vec3 objects containing
                           all vertices.
    - handleSRFS(names):   names is a list of strings containing the names
                           of all surface definitions.
    - handlePOLS(polys):   polys is a list of polygons. Each polygon is a
                           tuple (verts, surface) where verts is a list of
                           vertex indices and surface the index of the surface
                           that is used by this polygon.
    - handleCRVS(curves):  curves is a list of curves. Each curve is a
                           tuple (verts, surface, flags).
    - handlePCHS(patches): patches is a list of patches. Each patch is a
                           tuple (verts, surface).
    - handleSURF(surface): surface is a Surface object that contains the
                           surface attributes.
    """

    def __init__(self):
        """Constructor.
        """
        # The file handle
        self._file = None
        # The number of bytes that are left in the file.
        # After the FORM has been read, it is initialized with the
        # FORM length and is used to determine whether there are still
        # more chunks left or not.
        self._bytes_left = 12

        # Current Surface object when reading a SURF chunk
        self._current_surface = None
        # Current Texture object
        self._current_tex = None

    def read(self, file):
        """Read the lwo file.

        file is a file-like object that provides the read() and seek()
        methods.
        """
        self._file = file
        self._bytes_left = 12
        
        # Read the FORM header...
        tag, length = self.readChunkHeader()
        self.onChunk(tag, length, 0)
        self._bytes_left = length
        if tag!="FORM":
            raise LWOBError("Not a Lightwave object file")

        # Read the LWOB tag...
        lwob = self.readNBytes(4)
        if lwob!="LWOB":
            raise LWOBError("Not a Lightwave object file")

        # Read the contained chunks and call the 'reader' method if available.
        # If no reader is defined for a chunk this chunk is skipped.
        while not self._eofReached():
            tag,length = self.readChunkHeader()
            self.onChunk(tag, length, 1)
#            print tag, length
            handlername = "read%s"%tag
            handler = getattr(self, handlername, None)
            if handler==None:
                self.skipChunk(length)
            else:
                handler(length)


    # Handler methods. These methods can be overridden in derived classes.

    def onChunk(self, tag, length, level): pass

    def handlePNTS(self, points): pass
    def handleSRFS(self, names): pass
    def handlePOLS(self, polys): pass
    def handleCRVS(self, curves): pass
    def handlePCHS(self, patches): pass
    def handleSURF(self, surface): pass


    ### readXXXX() reader methods (XXXX is the chunk tag)
    ### A reader takes the chunk length as argument and has to read
    ### its data from the file and call an appropriate handler method
    ### (handleXXXX(...))

    def readPNTS(self, length):
        """Read a PNTS chunk and call the PNTS handler.
        """
        if length%12!=0:
            raise LWOBError("Invalid PNTS chunk size (%d)"%length)
        
        numpnts = int(length/12)
        # Read the chunk data...
        data = self.readNBytes(length)
        # Unpack the coordinates...
        coords = struct.unpack(">%s"%(numpnts*"fff"), data)
        # Initialize vec3s...
        pnts = []
        for i in range(numpnts):
            pnts.append(vec3(coords[i*3:(i+1)*3]))
            
        self.handlePNTS(pnts)

    def readSRFS(self, length):
        """Read a SRFS chunk and call the SRFS handler.
        """
        data = self.readNBytes(length)
        names = []
        while len(data)>0:
            n = data.find("\000")
            if n==-1:
                n = len(data)
            names.append(data[:n])
            if n%2==0:
                n+=1
            data = data[n+1:]

        self.handleSRFS(names)

    def readPOLS(self, length):
        """Read a POLS chunk and call the POLS handler.
        """
        data = self.readNBytes(length)
        idx = 0
        unpack = struct.unpack
        polys = []
        # The number of detail polygons that follow
        # (detail polygons are ignored)
        details_left = 0
        while idx<len(data):
            numverts = unpack(">H", data[idx:idx+2])[0]
            idx += 2
            verts = unpack(">%s"%(numverts*"H"), data[idx:idx+numverts*2])
            idx += numverts*2
            surf = unpack(">H", data[idx:idx+2])[0]
            idx += 2
            if details_left==0:
                polys.append((verts, surf))
            else:
                details_left -= 1

            # Does the polygon have detail polygons?
            if surf<0:
                details_left = unpack(">H", data[idx:idx+2])[0]
                idx += 2

        self.handlePOLS(polys)

    def readCRVS(self, length):
        """Read a CRVS chunk and call the CRVS handler.
        """
        data = self.readNBytes(length)
        idx = 0
        unpack = struct.unpack
        curves = []
        while idx<len(data):
            numverts = unpack(">H", data[idx:idx+2])[0]
            idx += 2
            verts = unpack(">%s"%(numverts*"H"), data[idx:idx+numverts*2])
            idx += numverts*2
            surf,flags = unpack(">HH", data[idx:idx+4])
            idx += 4
            curves.append((verts, surf, flags))

        self.handleCRVS(curves)
        
    def readPCHS(self, length):
        """Read a PCHS chunk and call the PCHS handler.
        """
        data = self.readNBytes(length)
        idx = 0
        unpack = struct.unpack
        patches = []
        while idx<len(data):
            numverts = unpack(">H", data[idx:idx+2])[0]
            idx += 2
            verts = unpack(">%s"%(numverts*"H"), data[idx:idx+numverts*2])
            idx += numverts*2
            surf = unpack(">H", data[idx:idx+2])[0]
            idx += 2
            patched.append((verts, surf))

        self.handlePCHS(patches)

    def readSURF(self, length):
        """Read a SURF chunk and call the SURF handler.
        """
        # Read the surface name byte by byte...
        name = ""
        while 1:
            c = self.readNBytes(1)
            length -= 1
            if ord(c)==0:
                if len(name)%2==0:
                    # Skip the pad byte
                    self.readNBytes(1)
                    length -= 1
                break
            name += c

        self._current_surface = Surface(name)
        self._current_tex = None
        
#        print "SURFACE",name
        while length>0:
            tag,sublength = self.readSubChunkHeader()
            self.onChunk(tag, sublength, 2)
            length -= 6+sublength
#            print " ",tag,sublength
            handlername = "readSURF_%s"%tag
            handler = getattr(self, handlername, None)
            if handler==None:
                self.skipChunk(sublength)
            else:
                handler(sublength)

        self.handleSURF(self._current_surface)
        self._current_tex = None
        self._current_surface = None

    def readSURF_COLR(self, length):
        """COLR chunk."""
        if length!=4:
            raise LWOBError("Invalid COLR size (%d instead of 4)"%length)
        data = self.readNBytes(4)
        color = struct.unpack(">BBB", data[:3])
        self._current_surface.color = color

    def readSURF_FLAG(self, length):
        """FLAG chunk."""
        if length!=2:
            raise LWOBError("Invalid FLAG size (%d instead of 2)"%length)
        data = self.readNBytes(2)
        flags = struct.unpack(">H", data)[0]
        self._current_surface.flags = flags
        
    def readSURF_LUMI(self, length):
        """LUMI chunk."""
        if length!=2:
            raise LWOBError("Invalid LUMI size (%d instead of 2)"%length)
        data = self.readNBytes(2)
        lumi = struct.unpack(">H", data)[0]
        if self._current_surface.luminosity==None:
            self._current_surface.luminosity = lumi/256.0

    def readSURF_VLUM(self, length):
        """VLUM chunk."""
        if length!=4:
            raise LWOBError("Invalid VLUM size (%d instead of 4)"%length)
        data = self.readNBytes(4)
        lum = struct.unpack(">f", data)[0]
        self._current_surface.luminosity = lum

    def readSURF_DIFF(self, length):
        """DIFF chunk."""
        if length!=2:
            raise LWOBError("Invalid DIFF size (%d instead of 2)"%length)
        data = self.readNBytes(2)
        diff = struct.unpack(">H", data)[0]
        if self._current_surface.diffuse==None:
            self._current_surface.diffuse = diff/256.0

    def readSURF_VDIF(self, length):
        """VDIF chunk."""
        if length!=4:
            raise LWOBError("Invalid VDIF size (%d instead of 4)"%length)
        data = self.readNBytes(4)
        diff = struct.unpack(">f", data)[0]
        self._current_surface.diffuse = diff

    def readSURF_SPEC(self, length):
        """SPEC chunk."""
        if length!=2:
            raise LWOBError("Invalid SPEC size (%d instead of 2)"%length)
        data = self.readNBytes(2)
        spec = struct.unpack(">H", data)[0]
        if self._current_surface.specular==None:
            self._current_surface.specular = spec/256.0

    def readSURF_VSPC(self, length):
        """VSPC chunk."""
        if length!=4:
            raise LWOBError("Invalid VSPC size (%d instead of 4)"%length)
        data = self.readNBytes(4)
        spec = struct.unpack(">f", data)[0]
        self._current_surface.specular = spec

    def readSURF_REFL(self, length):
        """REFL chunk."""
        if length!=2:
            raise LWOBError("Invalid REFL size (%d instead of 2)"%length)
        data = self.readNBytes(2)
        refl = struct.unpack(">H", data)[0]
        if self._current_surface.reflection==None:
            self._current_surface.reflection = refl/256.0

    def readSURF_VRFL(self, length):
        """VRFL chunk."""
        if length!=4:
            raise LWOBError("Invalid VRFL size (%d instead of 4)"%length)
        data = self.readNBytes(4)
        refl = struct.unpack(">f", data)[0]
        self._current_surface.reflection = refl

    def readSURF_TRAN(self, length):
        """TRAN chunk."""
        if length!=2:
            raise LWOBError("Invalid TRAN size (%d instead of 2)"%length)
        data = self.readNBytes(2)
        tran = struct.unpack(">H", data)[0]
        if self._current_surface.transparency==None:
            self._current_surface.transparency = tran/256.0

    def readSURF_VTRN(self, length):
        """VTRN chunk."""
        if length!=4:
            raise LWOBError("Invalid VTRN size (%d instead of 4)"%length)
        data = self.readNBytes(4)
        tran = struct.unpack(">f", data)[0]
        self._current_surface.transparency = tran

    def readSURF_GLOS(self, length):
        """GLOS chunk."""
        if length!=2:
            raise LWOBError("Invalid GLOS size (%d instead of 2)"%length)
        data = self.readNBytes(2)
        glos = struct.unpack(">H", data)[0]
        self._current_surface.glossiness = glos

    def readSURF_RFLT(self, length):
        """RFLT chunk."""
        if length!=2:
            raise LWOBError("Invalid RFLT size (%d instead of 2)"%length)
        data = self.readNBytes(2)
        mode = struct.unpack(">H", data)[0]
        self._current_surface.reflectionmode = mode

    def readSURF_RIMG(self, length):
        """RIMG chunk."""
        data = self.readNBytes(length)
        n = data.find("\000")
        if n!=-1:
            data = data[:n]
        self._current_surface.refmap = data

    def readSURF_RSAN(self, length):
        """RSAN chunk."""
        if length!=4:
            raise LWOBError("Invalid RSAN size (%d instead of 4)"%length)
        data = self.readNBytes(4)
        deg = struct.unpack(">f", data)[0]
        self._current_surface.refmap_seamangle = deg

    def readSURF_RIND(self, length):
        """RIND chunk."""
        if length!=4:
            raise LWOBError("Invalid RIND size (%d instead of 4)"%length)
        data = self.readNBytes(4)
        ior = struct.unpack(">f", data)[0]
        self._current_surface.refractiveindex = ior

    def readSURF_EDGE(self, length):
        """EDGE chunk."""
        if length!=4:
            raise LWOBError("Invalid EDGE size (%d instead of 4)"%length)
        data = self.readNBytes(4)
        edge = struct.unpack(">f", data)[0]
        self._current_surface.edgetransp = edge

    def readSURF_SMAN(self, length):
        """SMAN chunk."""
        if length!=4:
            raise LWOBError("Invalid SMAN size (%d instead of 4)"%length)
        data = self.readNBytes(4)
        sman = struct.unpack(">f", data)[0]
        self._current_surface.maxsmoothangle = sman

    def readSURF_CTEX(self, length):
        """CTEX chunk."""
        data = self.readNBytes(length)
        n = data.find("\000")
        if n!=-1:
            data = data[:n]

        tex = Texture()
        tex.type = data
        self._current_tex = tex
        self._current_surface.tex_color.append(tex)

    def readSURF_DTEX(self, length):
        """DTEX chunk."""
        data = self.readNBytes(length)
        n = data.find("\000")
        if n!=-1:
            data = data[:n]

        tex = Texture()
        tex.type = data
        self._current_tex = tex
        self._current_surface.tex_diffuse.append(tex)

    def readSURF_STEX(self, length):
        """STEX chunk."""
        data = self.readNBytes(length)
        n = data.find("\000")
        if n!=-1:
            data = data[:n]

        tex = Texture()
        tex.type = data
        self._current_tex = tex
        self._current_surface.tex_specular.append(tex)

    def readSURF_RTEX(self, length):
        """RTEX chunk."""
        data = self.readNBytes(length)
        n = data.find("\000")
        if n!=-1:
            data = data[:n]

        tex = Texture()
        tex.type = data
        self._current_tex = tex
        self._current_surface.tex_reflection.append(tex)

    def readSURF_TTEX(self, length):
        """TTEX chunk."""
        data = self.readNBytes(length)
        n = data.find("\000")
        if n!=-1:
            data = data[:n]

        tex = Texture()
        tex.type = data
        self._current_tex = tex
        self._current_surface.tex_transparency.append(tex)

    def readSURF_LTEX(self, length):
        """LTEX chunk."""
        data = self.readNBytes(length)
        n = data.find("\000")
        if n!=-1:
            data = data[:n]

        tex = Texture()
        tex.type = data
        self._current_tex = tex
        self._current_surface.tex_luminosity.append(tex)

    def readSURF_BTEX(self, length):
        """BTEX chunk."""
        data = self.readNBytes(length)
        n = data.find("\000")
        if n!=-1:
            data = data[:n]

        tex = Texture()
        tex.type = data
        self._current_tex = tex
        self._current_surface.tex_bump.append(tex)

    def readSURF_TFLG(self, length):
        """TFLG chunk."""
        if length!=2:
            raise LWOBError("Invalid TFLG size (%d instead of 2)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TFLG chunk")
        data = self.readNBytes(2)
        flags = struct.unpack(">H", data)[0]
        self._current_tex.flags = flags

    def readSURF_TSIZ(self, length):
        """TSIZ chunk."""
        if length!=12:
            raise LWOBError("Invalid TSIZ size (%d instead of 12)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TSIZ chunk")
        data = self.readNBytes(12)
        v = vec3(struct.unpack(">fff", data))
        self._current_tex.size = v

    def readSURF_TCTR(self, length):
        """TCTR chunk."""
        if length!=12:
            raise LWOBError("Invalid TCTR size (%d instead of 12)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TCTR chunk")
        data = self.readNBytes(12)
        v = vec3(struct.unpack(">fff", data))
        self._current_tex.center = v

    def readSURF_TFAL(self, length):
        """TFAL chunk."""
        if length!=12:
            raise LWOBError("Invalid TFAL size (%d instead of 12)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TFAL chunk")
        data = self.readNBytes(12)
        v = vec3(struct.unpack(">fff", data))
        self._current_tex.falloff = v

    def readSURF_TVEL(self, length):
        """TVEL chunk."""
        if length!=12:
            raise LWOBError("Invalid TVEL size (%d instead of 12)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TVEL chunk")
        data = self.readNBytes(12)
        v = vec3(struct.unpack(">fff", data))
        self._current_tex.velocity = v

    def readSURF_TCLR(self, length):
        """TCLR chunk."""
        if length!=4:
            raise LWOBError("Invalid TCLR size (%d instead of 4)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TCLR chunk")
        data = self.readNBytes(4)
        col = struct.unpack(">BBB", data[:3])
        self._current_tex.color = v

    def readSURF_TVAL(self, length):
        """TVAL chunk."""
        if length!=2:
            raise LWOBError("Invalid TVAL size (%d instead of 2)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TVAL chunk")
        data = self.readNBytes(2)
        val = struct.unpack(">H", data)[0]
        self._current_tex.value = val

    def readSURF_TAMP(self, length):
        """TAMP chunk."""
        if length!=4:
            raise LWOBError("Invalid TAMP size (%d instead of 4)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TAMP chunk")
        data = self.readNBytes(4)
        amp = struct.unpack(">f", data)[0]
        self._current_tex.bumpamplitude = amp

    def readSURF_TFP0(self, length):
        """TFP0 chunk."""
        if length!=4:
            raise LWOBError("Invalid TFP0 size (%d instead of 4)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TFP0 chunk")
        data = self.readNBytes(4)
        f = struct.unpack(">f", data)[0]
        self._current_tex.floatparams[0] = f

    def readSURF_TFP1(self, length):
        """TFP1 chunk."""
        if length!=4:
            raise LWOBError("Invalid TFP1 size (%d instead of 4)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TFP1 chunk")
        data = self.readNBytes(4)
        f = struct.unpack(">f", data)[0]
        self._current_tex.floatparams[1] = f

    def readSURF_TFP2(self, length):
        """TFP2 chunk."""
        if length!=4:
            raise LWOBError("Invalid TFP2 size (%d instead of 4)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TFP2 chunk")
        data = self.readNBytes(4)
        f = struct.unpack(">f", data)[0]
        self._current_tex.floatparams[2] = f

    def readSURF_TFP3(self, length):
        """TFP3 chunk."""
        if length!=4:
            raise LWOBError("Invalid TFP3 size (%d instead of 4)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TFP3 chunk")
        data = self.readNBytes(4)
        f = struct.unpack(">f", data)[0]
        self._current_tex.floatparams[3] = f

    def readSURF_TFP4(self, length):
        """TFP4 chunk."""
        if length!=4:
            raise LWOBError("Invalid TFP4 size (%d instead of 4)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TFP4 chunk")
        data = self.readNBytes(4)
        f = struct.unpack(">f", data)[0]
        self._current_tex.floatparams[4] = f

    def readSURF_TFP5(self, length):
        """TFP5 chunk."""
        if length!=4:
            raise LWOBError("Invalid TFP5 size (%d instead of 4)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TFP5 chunk")
        data = self.readNBytes(4)
        f = struct.unpack(">f", data)[0]
        self._current_tex.floatparams[5] = f

    def readSURF_TFP6(self, length):
        """TFP6 chunk."""
        if length!=4:
            raise LWOBError("Invalid TFP6 size (%d instead of 4)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TFP6 chunk")
        data = self.readNBytes(4)
        f = struct.unpack(">f", data)[0]
        self._current_tex.floatparams[6] = f

    def readSURF_TFP7(self, length):
        """TFP7 chunk."""
        if length!=4:
            raise LWOBError("Invalid TFP7 size (%d instead of 4)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TFP7 chunk")
        data = self.readNBytes(4)
        f = struct.unpack(">f", data)[0]
        self._current_tex.floatparams[7] = f

    def readSURF_TFP8(self, length):
        """TFP8 chunk."""
        if length!=4:
            raise LWOBError("Invalid TFP8 size (%d instead of 4)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TFP8 chunk")
        data = self.readNBytes(4)
        f = struct.unpack(">f", data)[0]
        self._current_tex.floatparams[8] = f

    def readSURF_TFP9(self, length):
        """TFP9 chunk."""
        if length!=4:
            raise LWOBError("Invalid TFP9 size (%d instead of 4)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TFP9 chunk")
        data = self.readNBytes(4)
        f = struct.unpack(">f", data)[0]
        self._current_tex.floatparams[9] = f

    def readSURF_TIP0(self, length):
        """TIP0 chunk."""
        if length!=2:
            raise LWOBError("Invalid TIP0 size (%d instead of 2)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TIP0 chunk")
        data = self.readNBytes(2)
        i = struct.unpack(">h", data)[0]
        self._current_tex.intparams[0] = i

    def readSURF_TIP1(self, length):
        """TIP1 chunk."""
        if length!=2:
            raise LWOBError("Invalid TIP1 size (%d instead of 2)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TIP1 chunk")
        data = self.readNBytes(2)
        i = struct.unpack(">h", data)[0]
        self._current_tex.intparams[1] = i

    def readSURF_TIP2(self, length):
        """TIP2 chunk."""
        if length!=2:
            raise LWOBError("Invalid TIP2 size (%d instead of 2)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TIP2 chunk")
        data = self.readNBytes(2)
        i = struct.unpack(">h", data)[0]
        self._current_tex.intparams[2] = i

    def readSURF_TIP3(self, length):
        """TIP3 chunk."""
        if length!=2:
            raise LWOBError("Invalid TIP3 size (%d instead of 2)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TIP3 chunk")
        data = self.readNBytes(2)
        i = struct.unpack(">h", data)[0]
        self._current_tex.intparams[3] = i

    def readSURF_TIP4(self, length):
        """TIP4 chunk."""
        if length!=2:
            raise LWOBError("Invalid TIP4 size (%d instead of 2)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TIP4 chunk")
        data = self.readNBytes(2)
        i = struct.unpack(">h", data)[0]
        self._current_tex.intparams[4] = i

    def readSURF_TIP5(self, length):
        """TIP5 chunk."""
        if length!=2:
            raise LWOBError("Invalid TIP5 size (%d instead of 2)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TIP5 chunk")
        data = self.readNBytes(2)
        i = struct.unpack(">h", data)[0]
        self._current_tex.intparams[5] = i

    def readSURF_TIP6(self, length):
        """TIP6 chunk."""
        if length!=2:
            raise LWOBError("Invalid TIP6 size (%d instead of 2)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TIP6 chunk")
        data = self.readNBytes(2)
        i = struct.unpack(">h", data)[0]
        self._current_tex.intparams[6] = i

    def readSURF_TIP7(self, length):
        """TIP7 chunk."""
        if length!=2:
            raise LWOBError("Invalid TIP7 size (%d instead of 2)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TIP7 chunk")
        data = self.readNBytes(2)
        i = struct.unpack(">h", data)[0]
        self._current_tex.intparams[7] = i

    def readSURF_TIP8(self, length):
        """TIP8 chunk."""
        if length!=2:
            raise LWOBError("Invalid TIP8 size (%d instead of 2)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TIP8 chunk")
        data = self.readNBytes(2)
        i = struct.unpack(">h", data)[0]
        self._current_tex.intparams[8] = i

    def readSURF_TIP9(self, length):
        """TIP9 chunk."""
        if length!=2:
            raise LWOBError("Invalid TIP9 size (%d instead of 2)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TIP9 chunk")
        data = self.readNBytes(2)
        i = struct.unpack(">h", data)[0]
        self._current_tex.intparams[9] = i

    def readSURF_TIMG(self, length):
        """TIMG chunk."""
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TIMG chunk")
        data = self.readNBytes(length)
        n = data.find("\000")
        if n!=-1:
            data = data[:n]

        self._current_tex.imagename = data

    def readSURF_TALP(self, length):
        """TALP chunk."""
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TALP chunk")
        data = self.readNBytes(length)
        n = data.find("\000")
        if n!=-1:
            data = data[:n]

        self._current_tex.alphaname = data

    def readSURF_TWRP(self, length):
        """TWRP chunk."""
        if length!=4:
            raise LWOBError("Invalid TWRP size (%d instead of 4)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TWRP chunk")
        data = self.readNBytes(4)
        w,h = struct.unpack(">HH", data)
        self._current_tex.widthwrap = w
        self._current_tex.heightwrap = h

    def readSURF_TAAS(self, length):
        """TAAS chunk."""
        if length!=4:
            raise LWOBError("Invalid TAAS size (%d instead of 4)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TAAS chunk")
        data = self.readNBytes(4)
        aa = struct.unpack(">f", data)[0]
        self._current_tex.antialiasingstrength = aa

    def readSURF_TOPC(self, length):
        """TOPC chunk."""
        if length!=4:
            raise LWOBError("Invalid TOPC size (%d instead of 4)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the TOPC chunk")
        data = self.readNBytes(4)
        op = struct.unpack(">f", data)[0]
        self._current_tex.opacity = op

    def readSURF_SHDR(self, length):
        """SHDR chunk."""
        if self._current_tex==None:
            raise LWOBError("Invalid position of the SHDR chunk")
        name = self.readNBytes(length)
        n = name.find("\000")
        if n!=-1:
            name = name[:n]

        self._current_tex.shader.append(name)

    def readSURF_SDAT(self, length):
        """SDAT chunk."""
        if self._current_tex==None:
            raise LWOBError("Invalid position of the SDAT chunk")
        data = self.readNBytes(length)
        self._current_tex.shaderdata.append(data)

    def readSURF_IMSQ(self, length):
        """IMSQ chunk."""
        if length!=6:
            raise LWOBError("Invalid IMSQ size (%d instead of 6)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the IMSQ chunk")
        data = self.readNBytes(6)
        offset,flags,looplen = struct.unpack(">HHH", data)
        self._current_tex.seq_offset = offset
        self._current_tex.seq_flags = flags
        self._current_tex.seq_looplength = looplen

    def readSURF_FLYR(self, length):
        """FLYR chunk."""
        if length!=8:
            raise LWOBError("Invalid FLYR size (%d instead of 8)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the FLYR chunk")
        data = self.readNBytes(8)
        b,e = struct.unpack(">II", data)
        self._current_tex.flyer_begin = b
        self._current_tex.flyer_end = e

    def readSURF_IMCC(self, length):
        """IMCC chunk."""
        if length!=6:
            raise LWOBError("Invalid IMCC size (%d instead of 6)"%length)
        if self._current_tex==None:
            raise LWOBError("Invalid position of the IMCC chunk")
        data = self.readNBytes(6)
        speed,low,high = struct.unpack(">HHH", data)
        self._current_tex.cycle_speed = speed
        self._current_tex.cycle_low = low
        self._current_tex.cycle_high = high


    #####

    def skipChunk(self, length):
        """Skip a chunk/sub chunk in the file (without reading it).

        Sets the file's position to the begin of the next chunk.
        length is the chunk length as it is stored in the file.
        The method takes care of adding the pad byte if length is odd.
        """
        if length%2==1:
            length += 1
        self._file.seek(length, 1)
        self._bytes_left -= length

    def readChunkHeader(self):
        """Read the head of the next chunk.

        Returns a tuple (tag, length).
        """
        s = self.readNBytes(8)
        tag = s[:4]
        length = struct.unpack(">I", s[4:])[0]
        return tag, length

    def readSubChunkHeader(self):
        """Read the head of the next sub chunk.

        Returns a tuple (tag, length).
        """
        s = self.readNBytes(6)
        tag = s[:4]
        length = struct.unpack(">H", s[4:])[0]
        return tag, length

    def readNBytes(self, n):
        """Read n bytes.

        Throws an exception if less than n bytes were read.
        """
        s = self._file.read(n)
        self._bytes_left -= len(s)
        if len(s)!=n:
            raise LWOBError("premature end of file")
        return s

    def _eofReached(self):
        """Return True when the end of the file has been reached.
        """
        return self._bytes_left<=1
            
