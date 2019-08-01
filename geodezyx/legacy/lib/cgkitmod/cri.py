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
# Portions created by the Initial Developer are Copyright (C) 2008
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
# $Id: ri.py,v 1.4 2006/02/14 19:29:39 mbaas Exp $

"""High level cri module that can be used pretty much like the cgkit.ri module.
"""

import sys, re
from . import ri
from .cgtypes import vec3
try:
    import ctypes
    _has_ctypes = True
except ImportError:
    # If importing ctypes fails (which should only happen on older Python versions)
    # the loadRI() function can still be used to load the standard ri module
    # (and get the shortened function names without the "Ri" prefix).
    _has_ctypes = False
    
if _has_ctypes:
    from . import _cri
    from ._cri import importRINames

try:
    import numpy
    from numpy.ctypeslib import ndpointer
    _has_numpy = True
except ImportError:
    _has_numpy = False
    

def loadRI(libName):
    """Load a RenderMan library and return a module-like handle to it.
    
    libName is the name of a shared library that implements the RenderMan
    interface. The name can either be an absolute file name or just the
    name of the library (without suffix or "lib" prefix) in which case 
    the function tries to find the library file itself.
    The return value is the library handle as returned by the ctypes 
    LoadLibrary() function. This handle has already been prepared so that
    it can be used like a module that contains the RenderMan interface.
    
    It is also possible to pass in None or an empty string in which case
    the built-in ri module will be returned.
    """
    if libName is None or libName=="":
        _ri = ri
    else: 
        if not _has_ctypes:
            raise ImportError("The ctypes module is not available. Please upgrade to a newer Python version (2.5 or newer) or install ctypes separately.")
        _ri = _cri.loadRI(libName)
        _ri = _RenderManAPI(_ri)
    
    # Create an alias for every Ri function and RI_ constant that has the prefix removed...
    for name in dir(_ri):
        if name.startswith("Ri"):
            setattr(_ri, name[2:], getattr(_ri, name))
        elif name.startswith("RI_"):
            setattr(_ri, name[3:], getattr(_ri, name))

    return _ri

class _ProceduralWrapper:
    """Helper class for RiProcedural.
    
    The main purpose of this class is to ensure that all relevant Python
    objects (the ctypes functions) are kept alive for as long as the
    procedural is active. And while we have an instance of this class
    around, it also stores the data value which is passed on the user
    functions. This means, the C data value is not used at all (as every
    procedural will have its own _ProceduralWrapper instance).
    
    The subdiv() and free() methods of this class are the actual functions
    that get passed to the C RiProcedural function.
    """
    def __init__(self, RMAPI, data, subdivFunc, freeFunc):
        # A reference to the "ri" object that has spawned this wrapper instance
        self.RMAPI = RMAPI
        # The wrapped up data value
        self.data = data
        # The user subdiv function
        self.subdivFunc = subdivFunc
        # The user free function (may be None)
        self.freeFunc = freeFunc
        
        # The ctypes function wrappers that get passed to the C function
        self.cSubdivideFunc = RMAPI.RtProcSubdivFunc(self.subdivide)
        self.cFreeFunc = RMAPI.RtProcFreeFunc(self.free)
        
    def subdivide(self, dataPtr, detail):
        """Subdivide function wrapper.
        
        Unwraps the data value and calls the user subdiv function.
        """
        self.subdivFunc(self.data, detail)
        
    def free(self, dataPtr):
        """Free function wrapper.
        
        Unwraps the data value, calls the user free function and
        removes this wrapper instance from the ri object (because
        the associated Python objects need no longer be kept alive).
        """
        if self.freeFunc is not None:
            self.freeFunc(self.data)
        try:
            del self.RMAPI._procWrappers[self]
        except KeyError:
            pass
    

class _RenderManAPI:
    """RenderMan interface.
    
    An instance of this class can be used in a "module-like" fashion,
    i.e. it defines all the functions and constants from the RenderMan
    standard and behaves much like the "ri" module.
    """
    
    def __init__(self, rimod):
        """Constructor.
        
        rimod is the underlying ctypes library handle as returned by _cri.loadRI().
        """
        self._ri = rimod
        
        # Regular expression to parse declarations
        self._declRe = re.compile(r"^\s*(?:(constant|uniform|varying|vertex) )?\s*"
                                   "(float|integer|string|color|point|vector|normal|matrix|hpoint)\s*"
                                   "(?:\[\s*([0-9]+)\s*\])?\s*(\w+)$")
        
        # Copy the RI_ and Rt attributes...
        for name in dir(rimod):
            if name[:2] in ["Rt", "RI"]:
                setattr(self, name, getattr(rimod, name))
                
        # "Import" the error handlers, filters and basis functions...
        for name in ["RiErrorAbort", "RiErrorIgnore", "RiErrorPrint",
                     "RiBoxFilter", "RiTriangleFilter", "RiCatmullRomFilter",
                     "RiGaussianFilter", "RiSincFilter",
                     "RiBezierBasis", "RiBSplineBasis", "RiCatmullRomBasis",
                     "RiHermiteBasis", "RiPowerBasis",
                     "RiProcDelayedReadArchive", "RiProcRunProgram", 
                     "RiProcDynamicLoad", "RiProcFree"]:
            if hasattr(rimod, name):
                setattr(self, name, getattr(rimod, name))

        # Variable declarations.
        # Key: Variable name - Value: (class, type, n)
        self._declarations = {}
        self._standardDeclarations(self._ri)
        
        # Stores the _ProceduralWrapper objects
        self._procWrappers = {}
        
        # Lookup table ctypes type -> numpy type
        self._numpyTypes = {}
        if _has_numpy:
            # Fill the _numpyTypes dict
            RtFloat = self._ri.RtFloat
            if RtFloat==ctypes.c_float:
                self._numpyTypes[RtFloat] = numpy.float32
            elif RtFloat==ctypes.c_double:
                self._numpyTypes[RtFloat] = numpy.float64
            else:
                raise TypeError("RtFloat is of an unknown type")
            
            # numpy.int_ is 64bit on a 64bit system (so it's rather long instead of int)
            # That's why we use an int32 here which should also match the
            # C int type on 64bit systems.
            self._numpyTypes[self._ri.RtInt] = numpy.int32

    def _getRiLastError(self):
        return self._ri.RiLastError
    
    def _setRiLastError(self, val):
        self._ri.RiLastError = val
        
    RiLastError = property(_getRiLastError, _setRiLastError)
    LastError = property(_getRiLastError, _setRiLastError)
        
    def RiArchiveBegin(self, archivename, *paramlist, **keyparams):
        """Begin an inline archive.
        
        Example: RiArchiveBegin("myarchive")
                 ...
                 RiArchiveEnd()
                 RiReadArchive("myarchive")
        """
        return self._ri.RiArchiveBegin(archivename, *self._createCParamList(paramlist, keyparams))
    
    def RiArchiveEnd(self):
        """Terminate an inline archive.
        
        Example: RiArchiveBegin("myarchive")
                 ...
                 RiArchiveEnd()
                 RiReadArchive("myarchive")
        """
        self._ri.RiArchiveEnd()

    def RiArchiveRecord(self, type, format, *args):
        """Output a user data record.
    
        type is one of RI_COMMENT, RI_STRUCTURE or RI_VERBATIM.
    
        Example: RiArchiveRecord(RI_COMMENT, "Frame %d", 2)
        """
        self._ri.RiArchiveRecord(type, format, *args)

    def RiAreaLightSource(self, name, *paramlist, **keyparams):
        """Start the definition of an area light and return the light handle.
    
        Example: RiAttributeBegin()
                 area1 = RiAreaLightSource("arealight", intensity=0.75)
                 ....
                 RiAttributeEnd()
                 RiIlluminate(area1, RI_TRUE)
        """
        self._ri.RiAreaLightSource(name, *self._createCParamList(paramlist, keyparams))

    def RiAtmosphere(self, name, *paramlist, **keyparams):
        """Set the current atmosphere shader.
    
        If name is RI_NULL then no atmosphere shader is used.
    
        Example: RiAtmosphere("fog")
        """
        self._ri.RiAtmosphere(name, *self._createCParamList(paramlist, keyparams))

    def RiAttribute(self, name, *paramlist, **keyparams):
        """Set an implementation-specific attribute.
    
        Example: RiAttribute("displacementbound", "sphere", 0.5)
        """
        a = self._createCParamList(paramlist, keyparams)
        self._ri.RiAttribute(name, *self._createCParamList(paramlist, keyparams))
        
    def RiAttributeBegin(self):
        """Push the current set of attributes onto the attribute stack.
    
        Example: RiAttributeBegin()
                 ...
                 RiAttributeEnd()
        """
        self._ri.RiAttributeBegin()

    def RiAttributeEnd(self):
        """Pops the current set of attributes from the attribute stack."""
        self._ri.RiAttributeEnd()

    def RiBasis(self, ubasis, ustep, vbasis, vstep):
        """Set the current basis for the u and v direction.
    
        ubasis/vbasis can either be one of the predefined basis matrices
        RiHermiteBasis, RiCatmullRomBasis, RiBezierBasis, RiBSplineBasis,
        RiPowerBasis or it can be a user defined matrix.
    
        For the predefined matrices there are also predefined variables
        which can be used for the step parameters:
        RI_HERMITESTEP, RI_CATMULLROMSTEP, RI_BEZIERSTEP, RI_BSPLINESTEP,
        RI_POWERSTEP.
    
        Example: RiBasis(RiBezierBasis, RI_BEZIERSTEP,
                         RiHermiteBasis, RI_HERMITESTEP)
        """
        ubasis = self._toCArray(self._ri.RtFloat, ubasis)
        vbasis = self._toCArray(self._ri.RtFloat, vbasis)
        self._ri.RiBasis(ubasis, ustep, vbasis, vstep)

    def RiBegin(self, name):
        """Starts the main block using a particular rendering method.
        
        The default renderer is selected by passing RI_NULL as name.
    
        Example: RiBegin(RI_NULL)
                 ...
                 RiEnd()
        """
        self._ri.RiBegin(name)
        
    def RiBlobby(self, nleaf, code, floats, strings, *paramlist, **keyparams):
        """Create a blobby surface.
    
        Number of array elements for primitive variables:
        -------------------------------------------------
        constant: 1              varying: nleaf
        uniform:  1              vertex:  nleaf
    
        Example: RiBlobby(2, [1001,0, 1003,0,16, 0,2,0,1],
                          [1.5,0,0,0, 0,1.5,0,0, 0,0,1.5,0, 0,0,-.1,1,
                          0.4, 0.01,0.3, 0.08], ["flat.zfile"])
        """
        code = self._toCArray(self._ri.RtInt, code)
        floats = self._toCArray(self._ri.RtFloat, floats)
        strings = self._toCArray(self._ri.RtString, strings)
        self._ri.RiBlobby(nleaf, len(code), code, len(floats), floats, len(strings), strings, *self._createCParamList(paramlist, keyparams))

    def RiBound(self, bound):
        """Set the bounding box for subsequent primitives.
    
        bound must be a sequence of six floating point values specifying
        the extent of the box along each coordinate direction:
        bound = [xmin, xmax, ymin, ymax, zmin, zmax]
    
        Example: RiBound([-1,1, 0,1, 0.5,0.75])
        """
        bound = self._toCArray(self._ri.RtFloat, bound)
        self._ri.RiBound(bound)

    def RiCamera(self, name, *paramlist, **keyparams):
        """Mark the current camera description.
    
        Example: RiCamera("rightcamera")
        """
        self._ri.RiCamera(name, *self._createCParamList(paramlist, keyparams))

    def RiClipping(self, near, far):
        """Sets the near and the far clipping plane along the direction of view.
    
        near and far must be positive values in the range from RI_EPSILON to
        RI_INFINITY.
    
        Example: RiClipping(0.1, 100)
        """
        self._ri.RiClipping(near, far)
        
    def RiClippingPlane(self, x, y, z, nx, ny, nz):
        """Adds a new clipping plane, defined by a surface point and its normal.
    
        All the geometry in positive normal direction is clipped.
    
        Example: RiClippingPlane(0,0,0, 0,0,-1) clips everything below the XY plane
        """
        self._ri.RiClippingPlane(x, y, z, nx, ny, nz)
    
    def RiColor(self, Cs):
        """Set the current color.
    
        Cs must be a sequence of at least N values where N is the number of
        color samples (set by RiColorSamples(), default is 3).
    
        Example: RiColor([0.2,0.5,0.2])
        """
        Cs = self._toCArray(self._ri.RtFloat, Cs)
        self._ri.RiColor(Cs)
        
    def RiColorSamples(self, nRGB, RGBn):
        """Redefine the number of color components to be used for specifying colors.
    
        nRGB is a n x 3 matrix that can be used to transform the n component color
        to a RGB color (n -> RGB).
        RGBn is just the opposite, its a 3 x n matrix that's used to transform
        a RGB color to a n component color (RGB -> n).
        Thus, the new number of color components is len(matrix)/3 (matrix is
        either nRGB or RGBn).
    
        Example: RiColorSamples([0.3,0.3,0.3], [1,1,1])
        """
        nRGB = self._toCArray(self._ri.RtFloat, nRGB)
        RGBn = self._toCArray(self._ri.RtFloat, RGBn)
        if len(nRGB)!=len(RGBn):
            raise ValueError("The conversion matrices must have the same number of elements.")
        if len(nRGB)%3!=0:
            raise ValueError("Invalid number of elements in the conversion matrices.")
        n = len(nRGB)//3
        self._ri.RiColorSamples(n, nRGB, RGBn)

    def RiConcatTransform(self, transform):
        """Concatenate a transformation onto the current transformation.
    
        transform must be a sequence that evaluates to 16 floating point
        values (4x4 matrix).
    
        Example: RiConcatTransform([2,0,0,0, 0,2,0,0, 0,0,2,0, 0,0,0,1])
                 RiConcatTransform([[2,0,0,0], [0,2,0,0], [0,0,2,0], [0,0,0,1]])
        """
        transform = self._toCArray(self._ri.RtFloat, transform)
        self._ri.RiConcatTransform(transform)
        
    def RiCone(self, height, radius, thetamax, *paramlist, **keyparams):
        """Create a cone (along the z axis).
    
        Number of array elements for primitive variables:
        -------------------------------------------------
        constant: 1              varying: 4
        uniform:  1              vertex:  4
    
        Example: RiCone(1.5, 0.7, 360)
        """
        self._ri.RiCone(height, radius, thetamax, *self._createCParamList(paramlist, keyparams))

    def RiContext(self, handle):
        """Set the current active rendering context.
    
        Example: ctx1 = RiGetContext()
                 ...
                 RiContext(ctx1)
        """
        self._ri.RiContext(handle)
    
    def RiCoordinateSystem(self, spacename):
        """Mark the current coordinate system with a name.
    
        Example: RiCoordinateSystem("lamptop")
        """
        self._ri.RiCoordinateSystem(spacename)
        
    def RiCoordSysTransform(self, spacename):
        """Replace the current transformation matrix with spacename.
    
        Example: RiCoordSysTransform("lamptop")
        """
        self._ri.RiCoordSysTransform(spacename)
    
    def RiCropWindow(self, left, right, bottom, top):
        """Specify a subwindow to render.
    
        The values each lie between 0 and 1.
    
        Example: RiCropWindow(0.0, 1.0 , 0.0, 1.0)  (renders the entire frame)
                 RiCropWindow(0.5, 1.0 , 0.0, 0.5)  (renders the top right quarter)
        """
        self._ri.RiCropWindow(left, right, bottom, top)
        
    def RiCurves(self, type, nvertices, wrap, *paramlist, **keyparams):
        """Create a number of curve primitives.
    
        type is either RI_LINEAR or RI_CUBIC.
        nvertices is an array with the number of vertices in each curve.
        wrap is either RI_PERIODIC or RI_NONPERIODIC.
        The width of the curves can be specified with the parameter
        RI_WIDTH (varying float) or RI_CONSTANTWIDTH (constant float).
    
        Number of array elements for primitive variables:
        -------------------------------------------------
        constant: 1              varying: #segments (depends on type and wrap)
        uniform:  #curves        vertex:  #points
    
        Example: RiCurves(RI_CUBIC, [4], RI_NONPERIODIC,
                          P=[0,0,0, -1,-0.5,1, 2,0.5,1, 1,0,-1],
                          width=[0.1, 0.04])
        """
        nvertices = self._toCArray(self._ri.RtInt, nvertices)
        self._ri.RiCurves(type, len(nvertices), nvertices, wrap, *self._createCParamList(paramlist, keyparams))

    def RiCylinder(self, radius,zmin,zmax,thetamax,*paramlist, **keyparams):
        """Create a cylinder (along the z axis).
    
        Number of array elements for primitive variables:
        -------------------------------------------------
        constant: 1              varying: 4
        uniform:  1              vertex:  4
    
        Example: RiCylinder(1.5, 0.0, 1.0, 360)
        """
        self._ri.RiCylinder(radius, zmin, zmax, thetamax, *self._createCParamList(paramlist, keyparams))
   
    def RiDeclare(self, name, declaration):
        """Declare the name and type of a variable.
    
        The syntax of the declaration is:  [class] [type] ['['n']']
    
        class ::= constant | uniform | varying | vertex
        type  ::= float | integer | string | color | point | vector | normal |
                  matrix | hpoint
        
        Example: RiDeclare("foo","uniform float")
                 RiDeclare("bar","constant integer [4]")
                 RiDeclare("mycolor", "varying color")
        """
        # Process the declaration internally so that parameter lists can be
        # constructed properly...
        decl = "%s %s"%(declaration, name)
        m = self._declRe.match(decl)
        if m is None:
            raise ValueError("Invalid declaration: %s"%decl)
        # Groups is a 4-tuple (class, type, n, name)
        grps = m.groups()
        self._declarations[grps[3]] = grps[:3]
        
        # Forward the call...
        return self._ri.RiDeclare(name, declaration)

    def RiDepthOfField(self, fstop, focallength, focaldistance):
        """Set depth of field parameters.
    
        If fstop is RI_INFINITY depth of field is turned off.
    
        Example: RiDepthOfField(22,45,1200)
        """
        self._ri.RiDepthOfField(fstop, focallength, focaldistance)
        
    def RiDetail(self, bound):
        """Set the current bounding box.
    
        bound must be a sequence of six floating point values specifying
        the extent of the box along each coordinate direction:
        bound = [xmin, xmax, ymin, ymax, zmin, zmax]
    
        Example: RiDetail([10,20,40,70,0,1])
        """
        bound = self._toCArray(self._ri.RtFloat, bound)
        self._ri.RiDetail(bound)
    
    def RiDetailRange(self, minvisible, lowertransition, uppertransition, maxvisible):
        """Set the current detail range.
    
        The values of the parameters must satisfy the following ordering:
        minvisible <= lowertransition <= uppertransition <= maxvisible
    
                    lowertransition  uppertransition
     visibility         |____________________|       
        ^               /                    \\
        |              /                      \\
        |_____________/                        \\_____________\\ Level of detail
                     |                          |            /
                 minvisible                      maxvisible
    
        Example: RiDetailRange(0,0,10,20)
        """
        self._ri.RiDetailRange(minvisible, lowertransition, uppertransition, maxvisible)

    def RiDisplayChannel(self, channel, *paramlist, **keyparams):
        """Defines a new display channel.
    
        Example: RiDisplayChannel("color aovCi", "string opacity", "aovOi")
        """
        self._ri.RiDisplayChannel(channel, *self._createCParamList(paramlist, keyparams))

    def RiDisk(self, height, radius, thetamax, *paramlist, **keyparams):
        """Create a disk (parallel to the XY plane).
    
        Number of array elements for primitive variables:
        -------------------------------------------------
        constant: 1              varying: 4
        uniform:  1              vertex:  4
    
        Example: RiDisk(0.0, 1.0, 360)"""
        self._ri.RiDisk(height, radius, thetamax, *self._createCParamList(paramlist, keyparams))

    def RiDisplacement(self, name, *paramlist, **keyparams):
        """Set the current displacement shader.
    
        Example: RiDisplacement("dented", km=1.5)
        """
        self._ri.RiDisplacement(name, *self._createCParamList(paramlist, keyparams))
        
    def RiDisplay(self, name, type, mode, *paramlist, **keyparams):
        """Specify the destination and type of the output.
    
        Example: RiDisplay("frame0001.tif", RI_FILE, RI_RGB)
                 RiDisplay("myimage.tif", RI_FRAMEBUFFER, RI_RGB)
        """
        self._ri.RiDisplay(name, type, mode, *self._createCParamList(paramlist, keyparams))
    
    def RiElse(self):
        """Add an else block to a conditional block.
        """
        
        self._ri.RiElse()

    def RiElseIf(self, expression, *paramlist, **keyparams):
        """Add an else-if block to a conditional block.
        """
        self._ri.RiElseIf(expression, *self._createCParamList(paramlist, keyparams))

    def RiEnd(self):
        """Terminates the main block.
        """
        self._ri.RiEnd()
        
    def RiErrorHandler(self, handler):
        try:
            self._ri.RiErrorHandler(handler)
        except ctypes.ArgumentError:
            self._ri.RiErrorHandler(self._ri.RtErrorHandler(handler))

    def RiExposure(self, gain, gamma):
        """Sets the parameters for the output color transformation.
    
        The transformation is color_out = (color_in*gain)^(1/gamma)
    
        Example: RiExposure(1.3, 2.2)
        """
        self._ri.RiExposure(gain, gamma)

    def RiExterior(self, name, *paramlist, **keyparams):
        """Set the current exterior volume shader.
    
        Example: RiExterior("fog")
        """
        self._ri.RiExterior(name, *self._createCParamList(paramlist, keyparams))
    
    def RiFormat(self, xres, yres, aspect):
        """Set the resolution of the output image and the aspect ratio of a pixel.
    
        Example: RiFormat(720,576,1)"""
        self._ri.RiFormat(xres, yres, aspect)
        
    def RiFrameAspectRatio(self, frameratio):
        """Set the ratio between width and height of the image.
    
        Example: RiFrameAspectRatio(4.0/3)
        """
        self._ri.RiFrameAspectRatio(frameratio)

    def RiFrameBegin(self, number):
        """Start a new frame.
    
        Example: RiFrameBegin(1)
                 ...
                 RiFrameEnd()
        """
        self._ri.RiFrameBegin(number)

    def RiFrameEnd(self):
        """Terminates a frame."""
        self._ri.RiFrameEnd()

    def RiGeneralPolygon(self, nverts, *paramlist, **keyparams):
        """Create a general planar concave polygon with holes.
    
        Number of array elements for primitive variables:
        -------------------------------------------------
        constant: 1              varying: #total vertices
        uniform:  1              vertex:  #total vertices
    
        Example: RiGeneralPolygon([4,3], P=[0,0,0, 0,1,0, 0,1,1, 0,0,1, \\
                                            0,0.25,0.5, 0,0.75,0.75, 0,0.75,0.25])
        """
        nverts = self._toCArray(self._ri.RtInt, nverts)
        self._ri.RiGeneralPolygon(len(nverts), nverts, *self._createCParamList(paramlist, keyparams))

    def RiGeometricApproximation(self, type, value):
        """Sets parameters for approximating surfaces.
    
        Example: RiGeometricApproximation(RI_FLATNESS, 0.5)
        """
        self._ri.RiGeometricApproximation(type, value)
        
    def RiGeometry(self, type, *paramlist, **keyparams):
        """Create an implementation-specific geometric primitive.
    
        Example: RiGeometry("teapot")
        """
        self._ri.RiGeometry(type, *self._createCParamList(paramlist, keyparams))

    def RiGetContext(self):
        """Set the current active rendering context.
    
        Example: ctx1 = RiGetContext()
                 ...
                 RiContext(ctx1)
        """
        return self._ri.RiGetContext()
    
    def RiHider(self, type, *paramlist, **keyparams):
        """Choose a hidden-surface elimination technique.
    
        Example: RiHider(RI_HIDDEN)  (default)
        """
        self._ri.RiHider(type, *self._createCParamList(paramlist, keyparams))

    def RiHyperboloid(self, point1, point2, thetamax, *paramlist, **keyparams):
        """Create a hyperboloid (with the z axis as symmetry axis).
    
        Example: RiHyperboloid([1,0,0],[1,1,1],360)
        """
        point1 =self._toCArray(self._ri.RtFloat, point1)
        point2 =self._toCArray(self._ri.RtFloat, point2)
        self._ri.RiHyperboloid(point1, point2, thetamax, *self._createCParamList(paramlist, keyparams))

    def RiIdentity(self):
        """Set the current transformation to the identity.
    
        Example: RiIdentity()
        """
        self._ri.RiIdentity()
        
    def RiIfBegin(self, expression, *paramlist, **keyparams):
        """Begin a conditional block.
        """
        self._ri.RiIfBegin(expression, *self._createCParamList(paramlist, keyparams))

    def RiIfEnd(self):
        """Terminate a conditional block.
        """
        self._ri.RiIfEnd()

    def RiIlluminate(self, light, onoff):
        """Activate or deactive a light source.
    
        Example: RiIlluminate(lgt, RI_TRUE)
        """
        self._ri.RiIlluminate(light, onoff)
        
    def RiImager(self, name, *paramlist, **keyparams):
        """Set an imager shader.
    
        if name is RI_NULL, no imager shader is used.
    
        Example: RiImager("background", "color bgcolor", [0.3,0.3,0.9])
        """
        self._ri.RiImager(name, *self._createCParamList(paramlist, keyparams))
   
    def RiInterior(self, name, *paramlist, **keyparams):
        """Set the current interior volume shader.
    
        Example: RiInterior("water")
        """
        self._ri.RiInterior(name, *self._createCParamList(paramlist, keyparams))

    def RiLightSource(self, name, *paramlist, **keyparams):
        """Add another light source and return its light handle.
    
        name is the name of the light source shader. 
    
        Example: light1 = RiLightSource("distantlight", intensity=1.5)
        """
        return self._ri.RiLightSource(name, *self._createCParamList(paramlist, keyparams))

    def RiMakeBrickMap(self, ptcnames, bkmname, *paramlist, **keyparams):
        """Create a brick map file from a list of point cloud file names.
    
        Example: RiMakeBrickMap(["sphere.ptc", "box.ptc"], "spherebox.bkm", "float maxerror", 0.002)
        """
        n = len(ptcnames)
        names = (n*ctypes.c_char_p)(*ptcnames)
        self._ri.RiMakeBrickMap(n, names, bkmname, *self._createCParamList(paramlist, keyparams))

    def RiMakeCubeFaceEnvironment(self, px,nx,py,ny,pz,nz, texname, fov, filterfunc, swidth, twidth, *paramlist, **keyparams):
        """Convert six image files into an environment map.
    
        The px/nx images are the views in positive/negative x direction.
        fov is the field of view that was used to generate the individual images.
        filterfunc is either one of the predefined filter (RiGaussianFilter,
        RiBoxFilter, RiTriangleFilter, RiSincFilter, RiCatmullRomFilter) or
        a callable that takes four float arguments x, y, width, height. 
        swidth and twidth define the support of the filter.
    
        Example: RiMakeCubeFaceEnvironment("px.tif","nx.tif","py.tif","ny.tif",
                                           "pz.tif","nz.tif", "tex.tif", 92.0,
                                            RiGaussianFilter, 2,2)
        """
        swidth = self._ri.RtFloat(swidth)
        twidth = self._ri.RtFloat(twidth)
        # Try passing the function in directly first (in case it's a standard
        # filter). If this fails, function must be a custom filter.
        try:
            self._ri.RiMakeCubeFaceEnvironment(px, nx, py, ny, pz, nz, texname, fov, filterfunc, swidth, twidth, *self._createCParamList(paramlist, keyparams))
        except ctypes.ArgumentError:
            filterfunc = self._ri.RtFilterFunc(filterfunc)
            self._ri.RiMakeCubeFaceEnvironment(px, nx, py, ny, pz, nz, texname, fov, filterfunc, swidth, twidth, *self._createCParamList(paramlist, keyparams))

    def RiMakeLatLongEnvironment(self, picname, texname, filterfunc, swidth, twidth, *paramlist, **keyparams):
        """Convert an image file into an environment map.
    
        filterfunc is either one of the predefined filter (RiGaussianFilter,
        RiBoxFilter, RiTriangleFilter, RiSincFilter, RiCatmullRomFilter) or
        a callable that takes four float arguments x, y, width, height. 
        swidth and twidth define the support of the filter.
    
        Example: RiMakeLatLongEnvironment("img.tif", "tex.tif",
                                          RiGaussianFilter, 2,2)
        """
        swidth = self._ri.RtFloat(swidth)
        twidth = self._ri.RtFloat(twidth)
        # Try passing the function in directly first (in case it's a standard
        # filter). If this fails, function must be a custom filter.
        try:
            self._ri.RiMakeLatLongEnvironment(picname, texname, filterfunc, swidth, twidth, *self._createCParamList(paramlist, keyparams))
        except ctypes.ArgumentError:
            filterfunc = self._ri.RtFilterFunc(filterfunc)
            self._ri.RiMakeLatLongEnvironment(picname, texname, filterfunc, swidth, twidth, *self._createCParamList(paramlist, keyparams))
        
    def RiMakeShadow(self, picname, shadowname, *paramlist, **keyparams):
        """Transform a depth image into a shadow map.
    
        Example: RiMakeShadow("depthimg.tif", "shadow.tif")
        """
        self._ri.RiMakeShadow(picname, shadowname, *self._createCParamList(paramlist, keyparams))

    def RiMakeTexture(self, picname, texname, swrap, twrap, filterfunc, swidth, twidth, *paramlist, **keyparams):
        """Convert an image file into a texture file.
    
        swrap and twrap are one of RI_PERIODIC, RI_CLAMP or RI_BLACK.
        filterfunc is either one of the predefined filter (RiGaussianFilter,
        RiBoxFilter, RiTriangleFilter, RiSincFilter, RiCatmullRomFilter) or
        a callable that takes four float arguments x, y, width, height. 
        swidth and twidth define the support of the filter.
    
        Example: RiMakeTexture("img.tif", "tex.tif", RI_PERIODIC, RI_CLAMP, \\
                               RiGaussianFilter, 2,2)
        """
        swidth = self._ri.RtFloat(swidth)
        twidth = self._ri.RtFloat(twidth)
        # Try passing the function in directly first (in case it's a standard
        # filter). If this fails, function must be a custom filter.
        try:
            self._ri.RiMakeTexture(picname, texname, swrap, twrap, filterfunc, swidth, twidth, *self._createCParamList(paramlist, keyparams))
        except ctypes.ArgumentError:
            filterfunc = self._ri.RtFilterFunc(filterfunc)
            self._ri.RiMakeTexture(picname, texname, swrap, twrap, filterfunc, swidth, twidth, *self._createCParamList(paramlist, keyparams))            

    def RiMatte(self, onoff):
        """Indicates whether subsequent primitives are matte objects.
    
        Example: RiMatte(RI_TRUE)
        """
        self._ri.RiMatte(onoff)

    def RiMotionBegin(self, *times):
        """Start the definition of a moving primitive.
    
        You can specify the time values directly or inside a sequence,
        for example, RiMotionBegin(0,1) or RiMotionBegin([0,1]).
    
        Example: RiMotionBegin(0.0, 1.0)
                 RiTranslate(1.0, 0.0, 0.0)
                 RiTranslate(1.0, 2.0, 0.0)
                 RiMotionEnd()
        """
        # For some reason the time values must be doubles...
        times = tuple([ctypes.c_double(v) for v in self._flatten(times)])
        self._ri.RiMotionBegin(len(times), *times)
        
    def RiMotionEnd(self):
        "Terminates the definition of a moving primitive."
        self._ri.RiMotionEnd()

    def RiNuPatch(self, nu, uorder, uknot, umin, umax, nv, vorder, vknot, vmin, vmax, *paramlist, **keyparams):
        """Create a NURBS patch.
    
        Number of array elements for primitive variables:
        -------------------------------------------------
        constant: 1              varying: #segment corners
        uniform:  #segments      vertex:  nu*nv
        """
        uknot = self._toCArray(self._ri.RtFloat, uknot)
        vknot = self._toCArray(self._ri.RtFloat, vknot)
        self._ri.RiNuPatch(nu, uorder, uknot, umin, umax, nv, vorder, vknot, vmin, vmax, *self._createCParamList(paramlist, keyparams))

    def RiObjectBegin(self, *paramlist, **keyparams):
        """Start the definition of a retained model and return the object handle.
    
        Example: obj1 = RiObjectBegin()
                 ...
                 RiObjectEnd()
        """
        return self._ri.RiObjectBegin(*self._createCParamList(paramlist, keyparams))
        
    def RiObjectEnd(self):
        """Terminate the definition of a retained model."""
        self._ri.RiObjectEnd()
    
    def RiObjectInstance(self, handle):
        """Create an instance of a previously defined model.
    
        Example: RiObjectInstance(obj1)
        """
        self._ri.RiObjectInstance(handle)

    def RiOpacity(self, Os):
        """Set the current opacity.
    
        Os must be a sequence of at least N values where N is the number
        of color samples (set by RiColorSamples(), default is 3). The
        opacity values must lie in the range from 0 to 1 (where 0 means
        completely transparent and 1 means completely opaque).
    
        Example: RiOpacity([0,0,1])
        """
        Os = self._toCArray(self._ri.RtFloat, Os)
        self._ri.RiOpacity(Os)

    def RiOption(self, name, *paramlist, **keyparams):
        """Set an implementation-specific option.
    
        Example: RiOption("searchpath", "shader","~/shaders:&")
        """
        self._ri.RiOption(name, *self._createCParamList(paramlist, keyparams))

    def RiOrientation(self, orientation):
        """Set the orientation of subsequent surfaces.
    
        orientation is either RI_OUTSIDE, RI_INSIDE, RI_LH (left handed)
        or RI_RH (right handed).
        """
        self._ri.RiOrientation(orientation)

    def RiParaboloid(self, rmax, zmin, zmax, thetamax, *paramlist, **keyparams):
        """Create a paraboloid (with the z axis as symmetry axis).
    
        Number of array elements for primitive variables:
        -------------------------------------------------
        constant: 1              varying: 4
        uniform:  1              vertex:  4
    
        Example: RiParaboloid(1.0, 0.0, 1.0, 360)
        """
        self._ri.RiParaboloid(rmax, zmin, zmax, thetamax, *self._createCParamList(paramlist, keyparams))

    def RiPatch(self, type, *paramlist, **keyparams):
        """RiPatch(type, paramlist)
    
        type is one of RI_BILINEAR (4 vertices) or RI_BICUBIC (16 vertices).
    
        Number of array elements for primitive variables:
        -------------------------------------------------
        constant: 1              varying: 4
        uniform:  1              vertex:  4/16 (depends on type)
    
        Example: RiPatch(RI_BILINEAR, P=[0,0,0, 1,0,0, 0,1,0, 1,1,0])
        """
        self._ri.RiPatch(type, *self._createCParamList(paramlist, keyparams))

    def RiPatchMesh(self, type, nu, uwrap, nv, vwrap, *paramlist, **keyparams):
        """Create a mesh made of patches.
    
        type is one of RI_BILINEAR or RI_BICUBIC.
        uwrap/vwrap can be RI_PERIODIC or RI_NONPERIODIC.
        The number of control points is nu*nv.
    
        Number of array elements for primitive variables:
        -------------------------------------------------
        constant: 1              varying: #patch corners (depends on uwrap/vwrap)
        uniform:  #patches       vertex:  nu*nv (same as "P")
    
        """
        self._ri.RiPatchMesh(type, nu, uwrap, nv, vwrap, *self._createCParamList(paramlist, keyparams))

    def RiPerspective(self, fov):
        """Concatenate a perspective transformation onto the current transformation.
        
        Example: RiPerspective(45)"""
        self._ri.RiPerspective(fov)

    def RiPixelFilter(self, function, xwidth, ywidth):
        """Set a pixel filter function and its width in pixels.

        function is either one of the predefined filter (RiGaussianFilter,
        RiBoxFilter, RiTriangleFilter, RiSincFilter, RiCatmullRomFilter) or
        a callable that takes four float arguments x, y, width, height. 
        xwidth and ywidth define the support of the filter.
    
        Example: RiPixelFilter(RiGaussianFilter, 2.0, 1.0)
        """
        xwidth = self._ri.RtFloat(xwidth)
        ywidth = self._ri.RtFloat(ywidth)
        # Try passing the function in directly first (in case it's a standard
        # filter). If this fails, function must be a custom filter.
        try:
            self._ri.RiPixelFilter(function, xwidth, ywidth)
        except ctypes.ArgumentError:
            self._ri.RiPixelFilter(self._ri.RtFilterFunc(function), xwidth, ywidth)
    
    def RiPixelSamples(self, xsamples, ysamples):
        """Set the sampling rate in horizontal and vertical direction.
    
        Example: RiPixelSamples(2,2)"""
        self._ri.RiPixelSamples(xsamples, ysamples)

    def RiPixelVariance(self, variance):
        """Limit the acceptable variance in the output value of pixels.
    
        Example: RiPixelVariance(0.01)"""
        self._ri.RiPixelVariance(variance)
    
    def RiPoints(self, *paramlist, **keyparams):
        """Create individual points.
    
        The size of the points can be either set with the primitive variable
        RI_WIDTH (one float per point) or RI_CONSTANTWIDTH (one float for all
        points).
    
        Number of array elements for primitive variables:
        -------------------------------------------------
        constant: 1              varying: #points
        uniform:  1              vertex:  #points    
        """
        params = self._createCParamList(paramlist, keyparams)
        for i in range(0, len(params), 2):
            if params[i]==b"P":
                n = len(params[i+1])
                if n%3!=0:
                    raise ValueError('Invalid number of floats in the "P" parameter.')
                n //= 3
                break
        else:
            raise ValueError('Parameter "P" is missing.')
        self._ri.RiPoints(n, *params)

    def RiPointsGeneralPolygons(self, nloops, nverts, vertids, *paramlist, **keyparams):
        """Create a polyhedron made of general planar concave polygons.
    
        nloops:  The number of loops for each polygon
        nverts:  The number of vertices in each loop
        vertids: The vertex indices of the loop vertices (0-based)
        The vertices themselves are stored in the parameter list (parameter "P").
    
        Number of array elements for primitive variables:
        -------------------------------------------------
        constant: 1              varying: #vertices (*)
        uniform:  #polygons      vertex:  #vertices (*)
    
        (*) max(vertids)+1
        """
        nloops = self._toCArray(self._ri.RtInt, nloops)
        nverts = self._toCArray(self._ri.RtInt, nverts)
        vertids = self._toCArray(self._ri.RtInt, vertids)
        self._ri.RiPointsGeneralPolygons(len(nloops), nloops, nverts, vertids, *self._createCParamList(paramlist, keyparams))

    def RiPointsPolygons(self, nverts, vertids, *paramlist, **keyparams):
        """Create a polyhedron made of planar convex polygons that share vertices.
    
        nverts:  An array with the number of vertices in each polygon
        vertids: The vertex indices of the polygon vertices (0-based)
        The vertices themselves are stored in the parameter list (parameter "P").
    
        Number of array elements for primitive variables:
        -------------------------------------------------
        constant: 1              varying: #vertices (*)
        uniform:  #polygons      vertex:  #vertices (*)
    
        (*) max(vertids)+1
        """
        nverts = self._toCArray(self._ri.RtInt, nverts)
        vertids = self._toCArray(self._ri.RtInt, vertids)
        self._ri.RiPointsPolygons(len(nverts), nverts, vertids, *self._createCParamList(paramlist, keyparams))

    def RiPolygon(self, *paramlist, **keyparams):
        """Create a planar and convex polygon.
    
        The parameter list must include at least position ("P") information.
    
        Number of array elements for primitive variables:
        -------------------------------------------------
        constant: 1              varying: #vertices
        uniform:  1              vertex:  #vertices
    
        Example: RiPolygon(P=[0,1,0, 0,1,1, 0,0,1, 0,0,0])
        """
        params = self._createCParamList(paramlist, keyparams)
        for i in range(0, len(params), 2):
            if params[i]==b"P":
                n = len(params[i+1])
                if n%3!=0:
                    raise ValueError('Invalid number of floats in the "P" parameter.')
                n //= 3  
                break
        else:
            raise ValueError('Parameter "P" is missing.')
        self._ri.RiPolygon(n, *params)

    def RiProcedural(self, data, bound, subdividefunc, freefunc=None):
        """Declare a procedural model.
    
        subdividefunc and freefunc may either be the standard RenderMan
        procedurals (RiProcDelayedReadArchive, RiProcRunProgram,
        RiProcDynamicLoad and RiProcFree) or Python callables.
        In the former case, data must be a sequence of strings or a single
        string containing the data for the functions. In the latter case,
        data may be any Python object which is just passed on to the
        functions.
        freefunc is optional and defaults to None.
    
        Example: RiProcedural("mymodel.rib", [-1,1,-1,1,-1,1], \\
                              RiProcDelayedReadArchive, RI_NULL)
                              
                 RiProcedural(["python teapot.py",""],[0,1,0,1,0,1], \\
                              RiProcRunProgram, RI_NULL)
                              
                 RiProcedural(["teapot.so",""],[0,1,0,1,0,1], \\
                              RiProcDynamicLoad, RI_NULL)
        """
        bound = self._toCArray(self._ri.RtFloat, bound)
            
        callMode = 0
        # Try to convert the data into a string array. If this fails,
        # don't even try to do a call without function wrapper...
        try:
            if data is None:
                data2 = ""
            else:
                data2 = data
            # Put everything into a list and append an empty string. This is
            # done so that you can also just pass one single string in all cases.
            data2 = [data2]+[""]
            stringArray = self._toCArray(self._ri.RtString, data2)
        except:
            callMode = 1
        
        # Try to use the functions directly (this will only work for builtin functions)...
        if callMode==0:
            try:
                self._ri.RiProcedural(ctypes.byref(stringArray), bound, subdividefunc, freefunc)
                callMode = 2
            except ctypes.ArgumentError:
                callMode = 1
        
        # Use a wrapper function (this must be used for Python callables)...
        if callMode==1:
            procWrapper = _ProceduralWrapper(RMAPI=self, data=data, subdivFunc=subdividefunc, freeFunc=freefunc)
            # Keep a reference to the instance so that it is kept alive...
            self._procWrappers[procWrapper] = 1
            # Call the procedural function (the data value is now the wrapped
            # one from within the wrapper object)...
            self._ri.RiProcedural(0, bound, procWrapper.cSubdivideFunc, procWrapper.cFreeFunc)
        

    def RiProjection(self, name, *paramlist, **keyparams):
        """Specify a projection method.
    
        The standard projections are RI_PERSPECTIVE and RI_ORTHOGRAPHIC.
        The perspective projection takes one optional parameter, RI_FOV.
    
        Example: RiProjection(RI_PERSPECTIVE, fov=45)
        """
        self._ri.RiProjection(name, *self._createCParamList(paramlist, keyparams))

    def RiQuantize(self, type, one, min, max, ditheramplitude):
        """Set the quantization parameters for colors and depth.
    
        Example: RiQuantize(RI_RGBA, 2048, -1024, 3071, 1.0)
        """
        self._ri.RiQuantize(type, one, min, max, ditheramplitude)

    def RiReadArchive(self, filename, callback=None, *paramlist, **keyparams):
        """Include an archive file.
    
        RiExample: RiReadArchive("teapot.rib")
        """
        if hasattr(callback, "__call__"):
            callback = self._ri.RtArchiveCallback(callback)
        self._ri.RiReadArchive(filename, callback, *self._createCParamList(paramlist, keyparams))

    def RiRelativeDetail(self, relativedetail):
        """Set the factor for all level of detail calculations.
    
        Example: RiRelativeDetail(0.7)"""
        self._ri.RiRelativeDetail(relativedetail)

    def RiResource(self, handle, type, *paramlist, **keyparams):
        """Create or operate on a named resource of a particular type.
        """
        self._ri.RiResource(handle, type, *self._createCParamList(paramlist, keyparams))
    
    def RiResourceBegin(self):
        """Push the current set of resources. 
        """
        self._ri.RiResourceBegin()
    
    def RiResourceEnd(self):
        """Pop the current set of resources. 
        """
        self._ri.RiResourceEnd()

    def RiReverseOrientation(self):
        """Causes the current orientation to be toggled.
    
        Example: RiReverseOrientation()
        """
        self._ri.RiReverseOrientation()

    def RiRotate(self, angle, *axis):
        """Rotate about angle degrees about the given axis.
    
        The axis is either given as 3 scalars or a sequence of 3 scalars.
    
        Example: RiRotate(90, 1,0,0)
        """
        axis = self._toCArray(self._ri.RtFloat, axis)
        self._ri.RiRotate(angle, *tuple(axis))

    def RiScale(self, *scale):
        """Concatenate a scaling onto the current transformation.
    
        The scaling is either given as 3 scalars or a sequence of 3 scalars.
    
        Example: RiScale(2,2,2)"""
        scale = self._toCArray(self._ri.RtFloat, scale)
        self._ri.RiScale(*tuple(scale))

    def RiScopedCoordinateSystem(self, spacename):
        """Mark the current coordinate system with a name but store it on a separate stack.
    
        Example: RiScopedCoordinateSystem("lamptop")
        """
        self._ri.RiScopedCoordinateSystem(spacename)

    def RiScreenWindow(self, left, right, bottom, top):
        """Specify the extents of the output image on the image plane.
    
        Example: RiScreenWindow(-1,1,-1,1)
        """
        self._ri.RiScreenWindow(left, right, bottom, top)

    def RiShader(self, name, handle, *paramlist, **keyparams):
        """Set the current coshader.
    
        Example: RiShader("plastic", "plastic_layer", Kd=0.7, Ks=0.3)
        """
        self._ri.RiShader(name, handle, *self._createCParamList(paramlist, keyparams))

    def RiShadingInterpolation(self, type):
        """Specify how shading samples are interpolated.
    
        type can be RI_CONSTANT or RI_SMOOTH.
    
        Example: RiShadingInterpolation(RI_SMOOTH)"""
        self._ri.RiShadingInterpolation(type)

    def RiShadingRate(self, size):
        """Set the current shading rate to an area of size pixels.
    
        Example: RiShadingRate(1.0)
        """
        self._ri.RiShadingRate(size)

    def RiShutter(self, opentime, closetime):
        """Set the times at which the shutter opens and closes.
    
        Example: RiShutter(0.1, 0.9)
        """
        self._ri.RiShutter(opentime, closetime)

    def RiSides(self, nsides):
        """Specify the number of visible sides of subsequent surfaces.
    
        Example: RiSides(1)"""
        self._ri.RiSides(nsides)

    def RiSkew(self, angle, *vecs):
        """Concatenate a skew onto the current transformation.
    
        angle is given in degrees.
        The two vectors are each given as 3 scalars or a sequence
        of 3 scalars.
    
        Example: RiSkew(45, 0,1,0, 1,0,0)
        """
        vecs = self._flatten(vecs)
        self._ri.RiSkew(angle, *vecs)
        
    def RiSolidBegin(self, type):
        """Start the definition of a solid object.
    
        type is one of RI_PRIMITIVE, RI_UNION, RI_DIFFERENCE and RI_INTERSECTION.
    
        Example: RiSolidBegin(RI_INTERSECTION)
                 RiSolidBegin(RI_PRIMITIVE)
                 ...
                 RiSolidEnd()
                 RiSolidBegin(RI_PRIMITIVE)
                 ...
                 RiSolidEnd()
                 RiSolidEnd()
        """
        self._ri.RiSolidBegin(type)

    def RiSolidEnd(self):
        """Terminate the definition of a solid object."""
        self._ri.RiSolidEnd()

    def RiSphere(self, radius,zmin,zmax,thetamax,*paramlist, **keyparams):
        """Create a sphere.
    
        Number of array elements for primitive variables:
        -------------------------------------------------
        constant: 1              varying: 4
        uniform:  1              vertex:  4
    
        Example: RiSphere(1.0, -1.0, 1.0, 360)
        """
        self._ri.RiSphere(radius, zmin, zmax, thetamax, *self._createCParamList(paramlist, keyparams))

    def RiSubdivisionMesh(self, scheme, nverts, vertids, tags, nargs, intargs, floatargs, *paramlist, **keyparams):
        """Create a subdivision surface.
    
        The only standard scheme is currently "catmull-clark".
        nverts:  The number of vertices in each face
        vertids: The vertex indices of the face vertices (0-based)
        tags: A string array of tag names.
        nargs: The number of int and float args for each tag.
        intargs: The integer arguments.
        floatargs: The float arguments.
        The vertices themselves are stored in the parameter list (parameter "P").
        
    
        Number of array elements for primitive variables:
        -------------------------------------------------
        constant: 1              varying: #vertices (*)
        uniform:  #faces         vertex:  #vertices (*)
    
        (*) max(vertids)+1
        """
        nverts = self._toCArray(self._ri.RtInt, nverts)
        vertids = self._toCArray(self._ri.RtInt, vertids)
        tags = self._toCArray(self._ri.RtToken, tags)
        nargs = self._toCArray(self._ri.RtInt, nargs)
        intargs = self._toCArray(self._ri.RtInt, intargs)
        floatargs = self._toCArray(self._ri.RtFloat, floatargs)
        self._ri.RiSubdivisionMesh(scheme, len(nverts), nverts, vertids, len(tags), tags, nargs, intargs, floatargs, *self._createCParamList(paramlist, keyparams))

    def RiSurface(self, name, *paramlist, **keyparams):
        """Set the current surface shader.
    
        Example: RiSurface("plastic", Kd=0.7, Ks=0.3)"""
        self._ri.RiSurface(name, *self._createCParamList(paramlist, keyparams))

    def RiSystem(self, cmd):
        """Execute an arbitrary command in the same environment as the current rendering pass.
        """
        self._ri.RiSystem(cmd)

    def RiTextureCoordinates(self, s1, t1, s2, t2, s3, t3, s4, t4):
        """Set the current set of texture coordinates.
    
        Declares a projection from the unit square [(0,0), (1,0), (0,1), (1,1)]
        in parameter space to quadrilateral [(s1,t1), (s2,t2), (s3,t3), (s4,t4)]
        in texture space.
    
        Example: RiTextureCoordinates(0.0, 0.0, 2.0, -0.5, -0.5, 1.75, 3.0, 3.0)"""
        self._ri.RiTextureCoordinates(s1, t1, s2, t2, s3, t3, s4, t4)

    def RiTorus(self, major, minor, phimin, phimax, thetamax, *paramlist, **keyparams):
        """Create a torus (with the z axis as symmetry axis).
    
        Number of array elements for primitive variables:
        -------------------------------------------------
        constant: 1              varying: 4
        uniform:  1              vertex:  4
    
        Example: RiTorus(1.5, 0.1, 0, 360, 360)
        """
        self._ri.RiTorus(major, minor, phimin, phimax, thetamax, *self._createCParamList(paramlist, keyparams))
        
    def RiTransform(self, transform):
        """Set the current transformation.
    
        transform must be a sequence that evaluates to 16 floating point
        values (4x4 matrix).
    
        Example: RiTransform([2,0,0,0, 0,2,0,0, 0,0,2,0, 0,0,0,1])
                 RiTransform([[2,0,0,0], [0,2,0,0], [0,0,2,0], [0,0,0,1]])
        """
        transform = self._toCArray(self._ri.RtFloat, transform)
        self._ri.RiTransform(transform)

    def RiTransformBegin(self):
        """Push the current transformation on the transformation stack.
    
        Example: RiTransformBegin()
                 ...
                 RiTransformEnd()
        """
        self._ri.RiTransformBegin()

    def RiTransformEnd(self):
        """Pop the current transformation from the stack."""
        self._ri.RiTransformEnd()

    def RiTransformPoints(self, fromspace, tospace, points):
        """Transform a set of points from one space to another.
        
        points must be a sequence of points where each point is another
        sequence with 3 floats.
        The return value is a sequence of vec3 objects containing the
        transformed points. None is returned in case of an error.
        """
        points = (len(points)*self._ri.RtPoint)(*[tuple(p) for p in points])
        n = len(points)
        res = self._ri.RiTransformPoints(fromspace, tospace, n, points)
        if res==0:
            return None
        
        return list([vec3(*p) for p in points])
    
    def RiTranslate(self, *translate):
        """Concatenate a translation onto the current transformation.
    
        The translation is either given as 3 scalars or a sequence of
        3 scalars.
    
        Example: RiTranslate(1.2, 4.3, -0.5)
                 or
                 RiTranslate( (1.2, 4.3, -0.5) )
        """
        translate = self._toCArray(self._ri.RtFloat, translate)
        self._ri.RiTranslate(*tuple(translate))
    
    def RiTrimCurve(self, ncurves, order, knot, min, max, n, u, v, w):
        """Set the current trim curve.
        """
        ncurves = self._toCArray(self._ri.RtInt, ncurves)
        order = self._toCArray(self._ri.RtInt, order)
        knot = self._toCArray(self._ri.RtFloat, knot)
        min = self._toCArray(self._ri.RtFloat, min)
        max = self._toCArray(self._ri.RtFloat, max)
        n = self._toCArray(self._ri.RtInt, n)
        u = self._toCArray(self._ri.RtFloat, u)
        v = self._toCArray(self._ri.RtFloat, v)
        w = self._toCArray(self._ri.RtFloat, w)
        self._ri.RiTrimCurve(len(ncurves), ncurves, order, knot, min, max, n, u, v, w)

    def RiWorldBegin(self):
        """Start the world block.
    
        Example: RiWorldBegin()
                 ...
                 RiWorldEnd()
        """
        self._ri.RiWorldBegin()

    def RiWorldEnd(self):
        """Terminates the world block."""
        self._ri.RiWorldEnd()

        
    def _toCArray(self, ctype, seq):
        """Convert and flatten a sequence into a ctypes array.
        
        ctype is the base type of the array and seq is the sequence that is
        to be flattened and converted. The return value is a ctypes
        array type with ctype as element type.
        If seq is already a ctypes sequence or a numpy array then no
        data conversion is done and the input is used directly.
        """
        # Is the value already a ctypes array? Then there's nothing to do
        if isinstance(seq, ctypes.Array):
            if seq._type_!=ctype:
                raise TypeError("ctypes array should be of type %s instead of %s"%(getattr(ctype, "__name__", "?"), getattr(seq._type_, "__name__", "?")))
            return seq
        
        # Is the sequence a numpy array?
        if hasattr(seq, "ctypes"):
            try:
                self._assertNumPyArrayType(seq, ctype)
            except TypeError:
                raise TypeError("numpy error: %s"%sys.exc_info()[1])
            n = seq.size
            Cls = (n*ctype)
            seq = Cls.from_address(seq.ctypes.data)
            return seq
            
        # Convert the sequence into a ctypes array...
        seq = self._flatten(seq)
        return (len(seq)*ctype)(*seq)
    
    def _assertNumPyArrayType(self, seq, ctype):
        """Make sure a numpy array contains elements of the correct type.
        
        seq must be a numpy array and ctype is the ctypes type that the
        array elements must match.
        Raises a TypeError exception when seq is of the wrong type.
        """
        dtype = self._numpyTypes.get(ctype, None)
        if dtype is None:
            raise TypeError("numpy arrays cannot be used for this parameter type (%s)"%getattr(ctype, "__name__", "?"))
        np = ndpointer(dtype=dtype, flags='CONTIGUOUS')
        # Just call from_param() to test the type. The actual return value is not needed.
        np.from_param(seq)
       
    def _flatten(self, seq):
        """Return a list of the individual items in a (possibly nested) sequence.
        """
        res = []
        ScalarTypes = [int, float, str]
        for v in seq:
            vtype = type(v)
            # v=scalar?
            if vtype in ScalarTypes:
                res.append(v)
            # no standard scalar or string. Then it might be a sequence..
            else:
                # Check if it is really a sequence...
                try:
                    n = len(v)
                except:
                    res.append(v)
                    continue
                res += self._flatten(v)
        return res
        
    def _createCParamList(self, paramlist, keyparams):    
        """
        Combine the keyparams with the paramlist into one paramlist
        and convert sequence values.
        Appends None (RI_NULL) to the parameter list.
        The tokens are converted to bytes objects.
        """
        # Combine the paramlist and keyparams values...
        res = ri._merge_paramlist(paramlist, keyparams)
    
        # Check if the values need conversion...
        for i in range(0, len(res), 2):
            token = res[i].strip()
            res[i] = bytes(token, "ascii")
            
            ##### Determine the type of the variable #####
            
            # Try to process any inline declaration...
            m = self._declRe.match(token)
            # No inline declaration or invalid declaration...
            if m is None:
                # If the token doesn't contain a space then it must have been 
                # just a name without inline declaration. In this case, try
                # to get the declaration from the previously declared tokens...
                if token.find(" ")==-1:
                    decl = self._declarations.get(token, None)
                    if decl is None:
                        raise ValueError('Token "%s" is not declared.'%token)
                else:
                    raise ValueError('Invalid inline declaration: %s'%token)
            else:
                decl = m.groups()[:3]
                
            cls,vtype,n = decl
            
            # The actual size is not checked here, so we only determine
            # the "base" type...
            if vtype=="integer":
                ctype = self._ri.RtInt
            elif vtype=="string":
                ctype = self._ri.RtString
            else:
                ctype = self._ri.RtFloat 

            # Check if the value is a sequence. If not, turn it into a list...
            if type(res[i+1]) is str:
                res[i+1] = [res[i+1]]
            try:
                n = len(res[i+1])
            except:
                res[i+1] = [res[i+1]]
            
            # Convert the value(s) into a ctypes array (even single values
            # must be converted)...
            try:
                res[i+1] = self._toCArray(ctype, res[i+1])
            except TypeError:
                raise TypeError('Parameter "%s": %s'%(res[i], sys.exc_info()[1]))

        res.append(None)
        return res

    def _standardDeclarations(self, ri):
        """Predeclare standard parameters.
        """
        decls = {ri.RI_AMPLITUDE:"uniform float",
                 ri.RI_BACKGROUND:"color",
                 ri.RI_BEAMDISTRIBUTION:"float",
                 ri.RI_CONEANGLE:"float",
                 ri.RI_CONEDELTAANGLE:"float",
                 ri.RI_CONSTANTWIDTH:"constant float",
                 ri.RI_CS:"varying color",
                 ri.RI_DISTANCE:"float",
                 ri.RI_FOV:"float",
                 ri.RI_FROM:"point",
                 ri.RI_HANDLEID:"string",
                 ri.RI_INTENSITY:"float",
                 ri.RI_KA:"uniform float",
                 ri.RI_KD:"uniform float",
                 ri.RI_KR:"uniform float",
                 ri.RI_KS:"uniform float",
                 ri.RI_LIGHTCOLOR:"color",
                 ri.RI_MAXDISTANCE:"float",
                 ri.RI_MINDISTANCE:"float",
                 ri.RI_N:"varying normal",
                 ri.RI_NP:"uniform normal",
                 ri.RI_ORIGIN:"integer[2]",
                 ri.RI_OS:"varying color",
                 ri.RI_P:"vertex point", 
                 ri.RI_PW:"vertex hpoint",
                 ri.RI_PZ:"vertex point",
                 ri.RI_ROUGHNESS:"uniform float",
                 ri.RI_S:"varying float", 
                 ri.RI_SPECULARCOLOR:"uniform color",
                 ri.RI_ST:"varying float[2]",
                 ri.RI_T:"varying float",
                 ri.RI_TEXTURENAME:"string",
                 ri.RI_TO:"point",
                 ri.RI_WIDTH:"varying float",
                 "shader":"string",
                 "archive":"string",
                 "texture":"string",
                 "procedural":"string",
                 "endofframe":"integer",
                 "sphere":"float",
                 "coordinatesystem":"string",
                 "name":"string",
                 "sense":"string"
                }
        
        # Fill the declarations dict...
        for name,decl in list(decls.items()):
            m = self._declRe.match("%s %s"%(decl, name))
            grps = m.groups()
            self._declarations[grps[3]] = grps[:3]
