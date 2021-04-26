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
# $Id: ri.py,v 1.4 2006/02/14 19:29:39 mbaas Exp $

"""A RenderMan(R) binding for Python.

This module contains a complete binding of Pixar's RenderMan API. The
binding was written to be compliant to v3.2 of Pixar's RenderMan
Interface specification. However, it also supports some newer features
such as string handles for light sources or object instances.

It is safe to import the module using::

  from cgkit.ri import *

All the functions that get imported start with the prefix Ri, all
constants start with RI_ or RIE_, so you probably won't get into a
naming conflict.

After importing the module this way you can use the functions just as
you are used to from the C API::

  from cgkit.ri import *

  RiBegin(RI_NULL)
  RiWorldBegin()
  RiSurface("plastic")
  RiSphere(1,-1,1,360)
  RiWorldEnd()
  RiEnd()

For more details on using the module see the cgkit manual at
http://cgkit.sourceforge.net/
"""

import sys, time, os, os.path, getpass, inspect, gzip
try:
    from ._core import vec3 as _vec3
except:
    from .cgtypes import vec3 as _vec3

########################### Constants #############################

RI_NULL         = None

RI_TRUE         = 1
RI_FALSE        = 0

RI_HIDDEN       = "hidden"
RI_PAINT        = "paint"

RI_EPSILON      = 1.0e-10
RI_INFINITY     = 1.0e38

RI_FILE         = "file"
RI_FRAMEBUFFER  = "framebuffer"
RI_RGB          = "rgb"
RI_RGBA         = "rgba"
RI_RGBZ         = "rgbZ"
RI_RGBAZ        = "rgbaz"
RI_A            = "a"
RI_Z            = "z"
RI_AZ           = "az"

RI_ORIGIN       = "origin"

RI_PERSPECTIVE  = "perspective"
RI_ORTHOGRAPHIC = "orthographic"
RI_FOV          = "fov"

RI_LH           = "lh"
RI_RH           = "rh"
RI_INSIDE       = "inside"
RI_OUTSIDE      = "outside"

RI_BILINEAR     = "bilinear"
RI_BICUBIC      = "bicubic"

RI_LINEAR       = "linear"
RI_CUBIC        = "cubic"

RI_CONSTANT     = "constant"
RI_SMOOTH       = "smooth"

RI_P            = "P"
RI_PW           = "Pw"
RI_PZ           = "Pz"
RI_N            = "N"
RI_NP           = "Np"
RI_NG           = "Ng"
RI_CI           = "Ci"
RI_OI           = "Oi"
RI_CS           = "Cs"
RI_OS           = "Os"
RI_S            = "s"
RI_T            = "t"
RI_ST           = "st"

RI_COMMENT      = "comment"
RI_STRUCTURE    = "structure"
RI_VERBATIM     = "verbatim"

RI_HERMITESTEP    = 2
RI_CATMULLROMSTEP = 1
RI_BEZIERSTEP     = 3
RI_BSPLINESTEP    = 1
RI_POWERSTEP      = 4

RI_PERIODIC     = "periodic"
RI_NONPERIODIC  = "nonperiodic"
RI_CLAMP        = "clamp"
RI_BLACK        = "black"

RI_FLATNESS     = "flatness"

RI_PRIMITIVE    = "primitive"
RI_UNION        = "union"
RI_DIFFERENCE   = "difference"
RI_INTERSECTION = "intersection"

RI_WIDTH        = "width"
RI_CONSTANTWIDTH = "constantwidth"

RI_HOLE         = "hole"
RI_CREASE       = "crease"
RI_CORNER       = "corner"
RI_INTERPOLATEBOUNDARY = "interpolateboundary"

RI_AMBIENTLIGHT = "ambientlight"
RI_POINTLIGHT   = "pointlight"
RI_DISTANTLIGHT = "distantlight"
RI_SPOTLIGHT    = "spotlight"

RI_INTENSITY    = "intensity"
RI_LIGHTCOLOR   = "lightcolor"
RI_FROM         = "from"
RI_TO           = "to"
RI_CONEANGLE    = "coneangle"
RI_CONEDELTAANGLE = "conedeltaangle"
RI_BEAMDISTRIBUTION = "beamdistribution"

RI_MATTE        = "matte"
RI_METAL        = "metal"
RI_SHINYMETAL   = "shinymetal"
RI_PLASTIC      = "plastic"
RI_PAINTEDPLASTIC = "paintedplastic"

RI_KA           = "Ka"
RI_KD           = "Kd"
RI_KS           = "Ks"
RI_ROUGHNESS    = "roughness"
RI_KR           = "Kr"
RI_TEXTURENAME  = "texturename"
RI_SPECULARCOLOR = "specularcolor"

RI_DEPTHCUE     = "depthcue"
RI_FOG          = "fog"
RI_BUMPY        = "bumpy"

RI_MINDISTANCE  = "mindistance"
RI_MAXDISTANCE  = "maxdistance"
RI_BACKGROUND   = "background"
RI_DISTANCE     = "distance"
RI_AMPLITUDE    = "amplitude"

RI_RASTER       = "raster"
RI_SCREEN       = "screen"
RI_CAMERA       = "camera"
RI_WORLD        = "world"
RI_OBJECT       = "object"

RI_IDENTIFIER   = "identifier"
RI_NAME         = "name"
RI_SHADINGGROUP = "shadinggroup"

RI_IGNORE       = "ignore"
RI_PRINT        = "print"
RI_ABORT        = "abort"
RI_HANDLER      = "handler"

RI_HANDLEID     = "__handleid"

# Tokens specific to the cgkit binding...
RI_RIBOUTPUT    = "_riboutput"
RI_VERSION      = "_version"
RI_ROUND_NDIGITS = "roundndigits"
RI_NUM_SIGNIFICANT_DIGITS = "numsignificantdigits"
RI_FLOAT_FMT_STRING = "floatfmtstring"

# Error handling: severity levels
RIE_INFO        = 0
RIE_WARNING     = 1
RIE_ERROR       = 2
RIE_SEVERE      = 3

RIE_INCAPABLE   = 11

RIE_NOTOPTIONS  = 25       # Invalid state for options
RIE_NOTATTRIBS  = 26       # Invalid state for attributes
RIE_NOTPRIMS    = 27       # Invalid state for primitives
RIE_ILLSTATE    = 28       # Other invalid state 
RIE_RANGE       = 42       # Parameter out of range 
RIE_CONSISTENCY = 43       # Parameters inconsistent

RIE_INVALIDSEQLEN = 80     # A sequence hasn't had the required length
RIE_UNDECLARED    = 81     # An undeclared parameter is used

RiLastError     = 0

############################ Types ###################################

RtBoolean = bool
RtInt = int
RtFloat = float
RtString = str
RtToken = str
RtVoid = None
RtPointer = lambda x: x

RtColor = tuple
RtPoint = tuple
RtVector = tuple
RtNormal = tuple
RtHpoint = tuple
RtMatrix = tuple
RtBasis = tuple
RtBound = tuple

RtObjectHandle = lambda x: x
RtLightHandle = lambda x: x
RtContextHandle = lambda x: x

RtFilterFunc = lambda x: x
RtErrorHandler = lambda x: x
RtProcSubdivFunc = lambda x: x
RtProcFreeFunc = lambda x: x
RtArchiveCallback = lambda x: x


######################################################################

class RIBStream:
    """This class encapsulates the output stream.

    The version number is automatically placed into the stream before
    any "real" Ri calls are made. Output from RiArchiveRecord() will
    be placed before the version number. (Note: The version line is disabled
    for now).
    """
    
    def __init__(self, outstream):
        self.out = outstream
        self.output_version = 1

    def close(self):
        """Close the stream, unless it's stdout."""
        if self.out!=sys.stdout:
            self.out.close()

    def flush(self):
        """Flush the internal buffer."""
        self.out.flush()

    def write(self, data):
        """Write data into the stream."""
        if self.output_version:
            # The binding contains newer calls, so this version number
            # might not be accurate anyway.
#            self.out.write('version 3.03\n')
            self.output_version = 0
        self.out.write(data)

    def writeArchiveRecord(self, data):
        """Same as write() but suppresses the version number.

        This method is used by RiArchiveRecord(), everyone else uses
        write().
        """
        self.out.write(data)
        

###################### Standard error handlers #######################

def RiErrorIgnore(code, severity, message):
    """Standard error handler.

    Ignores error messages."""
    
    pass

def RiErrorPrint(code, severity, message):
    """Standard error handler.

    Prints the message to stderr."""

    if severity==RIE_WARNING:
        print("WARNING:", end=' ', file=sys.stderr)
    elif severity==RIE_ERROR or severity==RIE_SEVERE:
        print("ERROR (%d):"%(code), end=' ', file=sys.stderr)        
    print(message, file=sys.stderr)

def RiErrorAbort(code, severity, message):
    """Standard error handler.

    Prints the message to stderr and aborts if it was an error."""

    RiErrorPrint(code, severity, message)
    if severity>=RIE_ERROR:
        sys.exit(1)


class RIException(Exception):
    """RenderMan Interface exception

    This exception is thrown by the error handler RiErrorException()."""
    pass

def RiErrorException(code, severity, message):
    """This error handler raises an exception when an error occurs.

    If the "error" is only an info or warning message the message is
    printed to stderr, otherwise the exception RIException is thrown.
    The actual error message is given as an argument to the constructor of
    RIException (the line with the file name, line number and offending
    Ri call is removed. You will have that information in the Traceback).
    """
    if severity<RIE_ERROR:
        RiErrorPrint(code, severity, message)
    else:
        if message[:7]=="In file":
            n=message.find("\n")
            message=message[n+1:]
        raise RIException(message)

########################### Functions ################################

# RiErrorHandler
def RiErrorHandler(handler):
    """Install a new error handler.

    The handler takes three arguments: code, severity, message.
    Besides the three standard error handler RiErrorIgnore, RiErrorPrint
    and RiErrorAbort there's an additional error handler available called
    RiErrorException. Whenever an error occurs RiErrorException raises
    the exception RIException.

    If you use one of the standard error handlers the corresponding RIB
    request is written to the output. If you supply RiErrorException or
    your own handler then the handler is installed but no output is
    written to the output stream.

    The last error code is always stored in the variable RiLastError.
    Note: If you import the module with "from ri import *" you have to
    import it with "import ri" as well and you must access RiLastError
    via "ri.RiLastError" otherwise the variable will always be 0.

    Example: RiErrorHandler(RiErrorAbort)
    """

    global _errorhandler

    _errorhandler = handler

    if handler==RiErrorIgnore:
        _ribout.write('ErrorHandler "ignore"\n')
    elif handler==RiErrorPrint:
        _ribout.write('ErrorHandler "print"\n')
    elif handler==RiErrorAbort:
        _ribout.write('ErrorHandler "abort"\n')


# RiBegin
def RiBegin(name):
    """Starts the main block using a particular rendering method.
    
    The default renderer is selected by passing RI_NULL as name.
    Here this means the output is written to stdout.
    If the name has the extension ".rib" then the output is written into
    a file with that name. Otherwise the name is supposed to be an
    external renderer (e.g. "rendrib" (BMRT), "rgl" (BMRT), "aqsis" (Aqsis),
    "renderdl" (3Delight),...) which is started and fed with the data.

    Example: RiBegin(RI_NULL)
             ...
             RiEnd()
    """
    global _ribout, _colorsamples, _lighthandle, _errorhandler
    global _insideframe, _insideworld, _insideobject, _insidesolid
    global _insidemotion, _declarations

    _create_new_context()

    # Determine where the output should be directed to...
    if name==RI_NULL or name=="":
        # -> stdout
        outstream = sys.stdout
    else:
        root, ext = os.path.splitext(name)
        ext=ext.lower()
        if ext==".rib":
            # -> file (rib)
            outstream = open(name,"w")
        elif ext==".gz":
            outstream = gzip.open(name,"wb")
        else:
            # -> pipe
            outstream = os.popen(name,"w")

    _ribout = RIBStream(outstream)

    # Initialize internal variables
    _colorsamples = 3
    _lighthandle  = 0
    _errorhandler = RiErrorPrint
    _insideframe  = 0
    _insideworld  = 0
    _insideobject = 0
    _insidesolid  = 0
    _insidemotion = 0
    _init_declarations()


# RiEnd
def RiEnd():
    """Terminates the main block.
    """
    global _ribout

    _ribout.flush()
    
    if _ribout != sys.stdout:
        _ribout.close()
        _ribout = sys.stdout

    _destroy_context()

# RiWorldBegin
def RiWorldBegin():
    """Start the world block.

    Example: RiWorldBegin()
             ...
             RiWorldEnd()
    """

    global _insideworld

    if _insideworld:
        _error(RIE_ILLSTATE, RIE_ERROR, "World blocks cannot be nested.")
    
    _ribout.write("WorldBegin\n")
    _insideworld = 1

# RiWorldEnd
def RiWorldEnd():
    """Terminates the world block."""

    global _insideworld
    
    _ribout.write("WorldEnd\n")
    _insideworld = 0

# RiOption
def RiOption(name, *paramlist, **keyparams):
    """Set an implementation-specific option.

    Example: RiOption("searchpath", "shader","~/shaders:&")
    """
    global _ribout
    global _round_ndigits
    global _float_conversion_string

    # cgkit specific options?
    if name==RI_RIBOUTPUT:
        keyparams = _paramlist2lut(paramlist, keyparams)
        if keyparams.get(RI_VERSION, None)==0:
            # Disable the "version" call in the RIB stream...
            if hasattr(_ribout, "output_version"):
                _ribout.output_version = 0
        
        _round_ndigits = keyparams.get(RI_ROUND_NDIGITS, _round_ndigits)
        numDigits = keyparams.get(RI_NUM_SIGNIFICANT_DIGITS, None)
        if numDigits is not None:
            _float_conversion_string = "%%1.%dg"%int(numDigits)
        _float_conversion_string = keyparams.get(RI_FLOAT_FMT_STRING, _float_conversion_string)
        return
            
    _ribout.write('Option "'+name+'"'+_paramlist2string(paramlist, keyparams)+"\n")

# RiAttribute
def RiAttribute(name, *paramlist, **keyparams):
    """Set an implementation-specific attribute.

    Example: RiAttribute("displacementbound", "sphere", 0.5)
    """
    
    _ribout.write('Attribute "'+name+'"'+_paramlist2string(paramlist, keyparams)+"\n")

# RiAttributeBegin
def RiAttributeBegin():
    """Push the current set of attributes onto the attribute stack.

    Example: RiAttributeBegin()
             ...
             RiAttributeEnd()
    """
    
    _ribout.write("AttributeBegin\n")

# RiAttributeEnd
def RiAttributeEnd():
    """Pops the current set of attributes from the attribute stack."""

    _ribout.write("AttributeEnd\n")

# RiTransformBegin
def RiTransformBegin():
    """Push the current transformation on the transformation stack.

    Example: RiTransformBegin()
             ...
             RiTransformEnd()
    """
    
    _ribout.write("TransformBegin\n")

# RiTransformEnd
def RiTransformEnd():
    """Pop the current transformation from the stack."""

    _ribout.write("TransformEnd\n")

# RiFrameBegin
def RiFrameBegin(number):
    """Start a new frame.

    Example: RiFrameBegin(1)
             ...
             RiFrameEnd()
    """

    global _insideframe

    if _insideframe:
        _error(RIE_ILLSTATE, RIE_ERROR, "Frame blocks cannot be nested.")
            
    _ribout.write("FrameBegin %d\n"%number)
    _insideframe = 1
    

# RiFrameEnd
def RiFrameEnd():
    """Terminates a frame."""
    
    global _insideframe
    
    _ribout.write("FrameEnd\n")
    _insideframe = 0

# RiHider
def RiHider(type, *paramlist, **keyparams):
    """Choose a hidden-surface elimination technique.

    Example: RiHider(RI_HIDDEN)  (default)
    """
    
    if type==RI_NULL: type="null"
    _ribout.write('Hider "'+type+'"'+_paramlist2string(paramlist, keyparams)+"\n")

# RiSphere
def RiSphere(radius,zmin,zmax,thetamax,*paramlist, **keyparams):
    """Create a sphere.

    Number of array elements for primitive variables:
    -------------------------------------------------
    constant: 1              varying: 4
    uniform:  1              vertex:  4

    Example: RiSphere(1.0, -1.0, 1.0, 360)
    """

    _ribout.write('Sphere %s %s %s %s'%(radius, zmin, zmax, thetamax)+ \
                 _paramlist2string(paramlist, keyparams)+"\n")

# RiCone
def RiCone(height, radius, thetamax, *paramlist, **keyparams):
    """Create a cone (along the z axis).

    Number of array elements for primitive variables:
    -------------------------------------------------
    constant: 1              varying: 4
    uniform:  1              vertex:  4

    Example: RiCone(1.5, 0.7, 360)
    """

    _ribout.write('Cone %s %s %s'%(height, radius, thetamax)+ \
                 _paramlist2string(paramlist, keyparams)+"\n")

# RiDisk
def RiDisk(height, radius, thetamax, *paramlist, **keyparams):
    """Create a disk (parallel to the XY plane).

    Number of array elements for primitive variables:
    -------------------------------------------------
    constant: 1              varying: 4
    uniform:  1              vertex:  4

    Example: RiDisk(0.0, 1.0, 360)"""

    _ribout.write('Disk %s %s %s'%(height, radius, thetamax)+ \
                 _paramlist2string(paramlist, keyparams)+"\n")

# RiCylinder
def RiCylinder(radius,zmin,zmax,thetamax,*paramlist, **keyparams):
    """Create a cylinder (along the z axis).

    Number of array elements for primitive variables:
    -------------------------------------------------
    constant: 1              varying: 4
    uniform:  1              vertex:  4

    Example: RiCylinder(1.5, 0.0, 1.0, 360)
    """

    _ribout.write('Cylinder %s %s %s %s'%(radius, zmin, zmax, thetamax)+ \
                 _paramlist2string(paramlist, keyparams)+"\n")

# RiTorus
def RiTorus(major, minor, phimin, phimax, thetamax, *paramlist, **keyparams):
    """Create a torus (with the z axis as symmetry axis).

    Number of array elements for primitive variables:
    -------------------------------------------------
    constant: 1              varying: 4
    uniform:  1              vertex:  4

    Example: RiTorus(1.5, 0.1, 0, 360, 360)
    """

    _ribout.write('Torus %s %s %s %s %s'%(major, minor, phimin, phimax, thetamax)+ \
                 _paramlist2string(paramlist, keyparams)+"\n")

# RiHyperboloid
def RiHyperboloid(point1, point2, thetamax, *paramlist, **keyparams):
    """Create a hyperboloid (with the z axis as symmetry axis).

    Example: RiHyperboloid([1,0,0],[1,1,1],360)
    """

    p1 = _seq2list(point1, 3)
    p2 = _seq2list(point2, 3)
    _ribout.write('Hyperboloid %s %s %s%s\n'%(p1[1:-1], p2[1:-1], thetamax, _paramlist2string(paramlist, keyparams)))

# RiParaboloid
def RiParaboloid(rmax, zmin, zmax, thetamax, *paramlist, **keyparams):
    """Create a paraboloid (with the z axis as symmetry axis).

    Number of array elements for primitive variables:
    -------------------------------------------------
    constant: 1              varying: 4
    uniform:  1              vertex:  4

    Example: RiParaboloid(1.0, 0.0, 1.0, 360)
    """

    _ribout.write('Paraboloid %s %s %s %s'%(rmax, zmin, zmax, thetamax)+ \
                 _paramlist2string(paramlist, keyparams)+"\n")

# RiPolygon
def RiPolygon(*paramlist, **keyparams):
    """Create a planar and convex polygon.

    The parameter list must include at least position ("P") information.

    Number of array elements for primitive variables:
    -------------------------------------------------
    constant: 1              varying: #vertices
    uniform:  1              vertex:  #vertices

    Example: RiPolygon(P=[0,1,0, 0,1,1, 0,0,1, 0,0,0])
    """

    _ribout.write('Polygon'+_paramlist2string(paramlist, keyparams)+"\n")

# RiGeneralPolygon
def RiGeneralPolygon(nverts, *paramlist, **keyparams):
    """Create a general planar concave polygon with holes.

    Number of array elements for primitive variables:
    -------------------------------------------------
    constant: 1              varying: #total vertices
    uniform:  1              vertex:  #total vertices

    Example: RiGeneralPolygon([4,3], P=[0,0,0, 0,1,0, 0,1,1, 0,0,1, \\
                                        0,0.25,0.5, 0,0.75,0.75, 0,0.75,0.25])
    """

    _ribout.write('GeneralPolygon '+_seq2list(nverts)+ \
                 _paramlist2string(paramlist, keyparams)+"\n")

# RiPointsPolygons
def RiPointsPolygons(nverts, vertids, *paramlist, **keyparams):
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

    _ribout.write('PointsPolygons '+_seq2list(nverts)+' '+ \
                 _seq2list(vertids)+ \
                 _paramlist2string(paramlist, keyparams)+"\n")

# RiPointsGeneralPolygons
def RiPointsGeneralPolygons(nloops, nverts, vertids, *paramlist, **keyparams):
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

    _ribout.write('PointsGeneralPolygons '+_seq2list(nloops)+' '+ \
                 _seq2list(nverts)+' '+_seq2list(vertids)+ \
                 _paramlist2string(paramlist, keyparams)+"\n")


# Predefined basis matrices
RiHermiteBasis = "hermite"
RiCatmullRomBasis = "catmull-rom"
RiBezierBasis = "bezier"
RiBSplineBasis = "b-spline"
RiPowerBasis = "power"

# RiBasis
def RiBasis(ubasis, ustep, vbasis, vstep):
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

    if type(ubasis) is str:
        ubasis = '"'+ubasis+'"'
    else:
        ubasis = _seq2list(ubasis, 16)
        
    if type(vbasis) is str:
        vbasis = '"'+vbasis+'"'
    else:
        vbasis = _seq2list(vbasis, 16)
        
    _ribout.write('Basis '+ubasis+' '+str(ustep)+' '+vbasis+' '+str(vstep)+"\n")

# RiPatch
def RiPatch(type, *paramlist, **keyparams):
    """RiPatch(type, paramlist)

    type is one of RI_BILINEAR (4 vertices) or RI_BICUBIC (16 vertices).

    Number of array elements for primitive variables:
    -------------------------------------------------
    constant: 1              varying: 4
    uniform:  1              vertex:  4/16 (depends on type)

    Example: RiPatch(RI_BILINEAR, P=[0,0,0, 1,0,0, 0,1,0, 1,1,0])
    """

    _ribout.write('Patch "'+type+'"'+_paramlist2string(paramlist, keyparams)+"\n")

# RiPatchMesh
def RiPatchMesh(type, nu, uwrap, nv, vwrap, *paramlist, **keyparams):
    """Create a mesh made of patches.

    type is one of RI_BILINEAR or RI_BICUBIC.
    uwrap/vwrap can be RI_PERIODIC or RI_NONPERIODIC.
    The number of control points is nu*nv.

    Number of array elements for primitive variables:
    -------------------------------------------------
    constant: 1              varying: #patch corners (depends on uwrap/vwrap)
    uniform:  #patches       vertex:  nu*nv (same as "P")

    """

    _ribout.write('PatchMesh "'+type+'" '+str(nu)+' "'+uwrap+'" '+\
                 str(nv)+' "'+vwrap+'"'+\
                 _paramlist2string(paramlist, keyparams)+"\n")
    

# RiNuPatch
def RiNuPatch(nu, uorder, uknot, umin, umax, nv, vorder, vknot, vmin, vmax, *paramlist, **keyparams):
    """Create a NURBS patch.

    Number of array elements for primitive variables:
    -------------------------------------------------
    constant: 1              varying: #segment corners
    uniform:  #segments      vertex:  nu*nv
    """

    _ribout.write('NuPatch '+str(nu)+" "+str(uorder)+' '+_seq2list(uknot)+" "+ \
                 str(umin)+" "+str(umax)+" "+ \
                str(nv)+" "+str(vorder)+' '+_seq2list(vknot)+" "+ \
                 str(vmin)+" "+str(vmax)+_paramlist2string(paramlist, keyparams)+"\n")

# RiTrimCurve
def RiTrimCurve(ncurves, order, knot, min, max, n, u, v, w):
    """Set the current trim curve.
    """

    _ribout.write('TrimCurve '+_seq2list(ncurves)+' '+\
                 _seq2list(order)+' '+_seq2list(knot)+' '+\
                 _seq2list(min)+' '+_seq2list(max)+' '+_seq2list(n)+' '+ \
                 _seq2list(u)+' '+ \
                 _seq2list(v)+' '+ \
                 _seq2list(w)+'\n')

# RiPoints
def RiPoints(*paramlist, **keyparams):
    """Create individual points.

    The size of the points can be either set with the primitive variable
    RI_WIDTH (one float per point) or RI_CONSTANTWIDTH (one float for all
    points).

    Number of array elements for primitive variables:
    -------------------------------------------------
    constant: 1              varying: #points
    uniform:  1              vertex:  #points    
    """

    _ribout.write('Points'+_paramlist2string(paramlist, keyparams)+"\n")

# RiCurves
def RiCurves(type, nvertices, wrap, *paramlist, **keyparams):
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

    _ribout.write('Curves "'+type+'" '+_seq2list(nvertices)+' "'+wrap+'"'+
                  _paramlist2string(paramlist, keyparams)+'\n')

# RiSubdivisionMesh
def RiSubdivisionMesh(scheme, nverts, vertids, tags, nargs, intargs, floatargs, *paramlist, **keyparams):
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

    if len(tags)==0:
        _ribout.write('SubdivisionMesh "'+scheme+'" '+_seq2list(nverts)+' '+ \
                 _seq2list(vertids)+' '+ \
                 _paramlist2string(paramlist, keyparams)+"\n")
    else:
        _ribout.write('SubdivisionMesh "'+scheme+'" '+_seq2list(nverts)+' '+ \
                 _seq2list(vertids)+' '+_seq2list(tags)+' '+ \
                 _seq2list(nargs)+' '+_seq2list(intargs)+' '+ \
                 _seq2list(floatargs)+' '+ \
                 _paramlist2string(paramlist, keyparams)+"\n")

# RiBlobby
def RiBlobby(nleaf, code, floats, strings, *paramlist, **keyparams):
    """Create a blobby surface.

    Number of array elements for primitive variables:
    -------------------------------------------------
    constant: 1              varying: nleaf
    uniform:  1              vertex:  nleaf

    Example: RiBlobby(2, [1001,0, 1003,0,16, 0,2,0,1],
                      [1.5,0,0,0, 0,1.5,0,0, 0,0,1.5,0, 0,0,-.1,1,
                      0.4, 0.01,0.3, 0.08], ["flat.zfile"])
    """

    _ribout.write('Blobby '+str(nleaf)+' '+_seq2list(code)+' '+_seq2list(floats)+
                  ' '+_seq2list(strings)+
                  _paramlist2string(paramlist, keyparams)+'\n')

# RiColorSamples
def RiColorSamples(nRGB, RGBn):
    """Redefine the number of color components to be used for specifying colors.

    nRGB is a n x 3 matrix that can be used to transform the n component color
    to a RGB color (n -> RGB).
    RGBn is just the opposite, its a 3 x n matrix that's used to transform
    a RGB color to a n component color (RGB -> n).
    Thus, the new number of color components is len(matrix)/3 (matrix is
    either nRGB or RGBn).

    Example: RiColorSamples([0.3,0.3,0.3], [1,1,1])
    """
    global _colorsamples

    if len(nRGB)!=len(RGBn):
        _error(RIE_CONSISTENCY, RIE_ERROR,
               "The color transformation matrices must have the same number of values.")

    if len(nRGB)%3!=0 or len(nRGB)==0:
        _error(RIE_CONSISTENCY, RIE_ERROR,
               "The number of values in the transformation matrices must be a multiple of 3.")
        
    _colorsamples = len(_flatten(nRGB))//3
    _ribout.write('ColorSamples '+_seq2list(nRGB)+' '+_seq2list(RGBn)+'\n')

# RiColor
def RiColor(Cs):
    """Set the current color.

    Cs must be a sequence of at least N values where N is the number of
    color samples (set by RiColorSamples(), default is 3).

    Example: RiColor([0.2,0.5,0.2])
    """

    col=_seq2col(Cs)
    _ribout.write("Color "+col+"\n")

# RiOpacity
def RiOpacity(Os):
    """Set the current opacity.

    Os must be a sequence of at least N values where N is the number
    of color samples (set by RiColorSamples(), default is 3). The
    opacity values must lie in the range from 0 to 1 (where 0 means
    completely transparent and 1 means completely opaque).

    Example: RiOpacity([0,0,1])
    """

    col=_seq2col(Os)
    _ribout.write("Opacity "+col+"\n")

# RiShadingRate
def RiShadingRate(size):
    """Set the current shading rate to an area of size pixels.

    Example: RiShadingRate(1.0)
    """

    _ribout.write("ShadingRate %s\n"%size)

# RiShadingInterpolation
def RiShadingInterpolation(type):
    """Specify how shading samples are interpolated.

    type can be RI_CONSTANT or RI_SMOOTH.

    Example: RiShadingInterpolation(RI_SMOOTH)"""

    _ribout.write('ShadingInterpolation "'+type+'"\n')

# RiShader
def RiShader(name, handle, *paramlist, **keyparams):
    """Set the current coshader.

    Example: RiShader("plastic", "plastic_layer", Kd=0.7, Ks=0.3)"""

    _ribout.write('Shader "'+name+'"'+' "'+handle+'"'+ \
                 _paramlist2string(paramlist, keyparams)+"\n")
    
# RiSurface
def RiSurface(name, *paramlist, **keyparams):
    """Set the current surface shader.

    Example: RiSurface("plastic", Kd=0.7, Ks=0.3)"""

    _ribout.write('Surface "'+name+'"'+ \
                 _paramlist2string(paramlist, keyparams)+"\n")

# RiInterior
def RiInterior(name, *paramlist, **keyparams):
    """Set the current interior volume shader.

    Example: RiInterior("water")
    """

    _ribout.write('Interior "'+name+'"'+_paramlist2string(paramlist, keyparams)+"\n")

# RiExterior
def RiExterior(name, *paramlist, **keyparams):
    """Set the current exterior volume shader.

    Example: RiExterior("fog")
    """

    _ribout.write('Exterior "'+name+'"'+_paramlist2string(paramlist, keyparams)+"\n")

# RiAtmosphere
def RiAtmosphere(name, *paramlist, **keyparams):
    """Set the current atmosphere shader.

    If name is RI_NULL then no atmosphere shader is used.

    Example: RiAtmosphere("fog")
    """

    if name==RI_NULL:
        _ribout.write('Atmosphere\n')
    else:
        _ribout.write('Atmosphere "'+name+'"'+ \
                      _paramlist2string(paramlist, keyparams)+"\n")

# RiDisplacement
def RiDisplacement(name, *paramlist, **keyparams):
    """Set the current displacement shader.

    Example: RiDisplacement("dented", km=1.5)
    """

    _ribout.write('Displacement "'+name+'"'+_paramlist2string(paramlist, keyparams)+"\n")

# RiImager
def RiImager(name, *paramlist, **keyparams):
    """Set an imager shader.

    if name is RI_NULL, no imager shader is used.

    Example: RiImager("background", "color bgcolor", [0.3,0.3,0.9])
    """

    if name==RI_NULL:
        _ribout.write('Imager\n')
    else:
        _ribout.write('Imager "'+name+'"'+_paramlist2string(paramlist, keyparams)+"\n")

# RiClipping
def RiClipping(near, far):
    """Sets the near and the far clipping plane along the direction of view.

    near and far must be positive values in the range from RI_EPSILON to
    RI_INFINITY.

    Example: RiClipping(0.1, 100)
    """

    _ribout.write("Clipping %s %s\n"%(near,far))

# RiClippingPlane
def RiClippingPlane(x, y, z, nx, ny, nz):
    """Adds a new clipping plane, defined by a surface point and its normal.

    All the geometry in positive normal direction is clipped.

    Example: RiClippingPlane(0,0,0, 0,0,-1) clips everything below the XY plane
    """

    _ribout.write('ClippingPlane %s %s %s %s %s %s\n'%(x,y,z,nx,ny,nz))

# RiDisplay
def RiDisplay(name,type,mode, *paramlist, **keyparams):
    """Specify the destination and type of the output.

    Example: RiDisplay("frame0001.tif", RI_FILE, RI_RGB)
             RiDisplay("myimage.tif", RI_FRAMEBUFFER, RI_RGB)
    """

    _ribout.write('Display "'+name+'" "'+type+'" "'+mode+'"'+ \
                 _paramlist2string(paramlist, keyparams)+"\n")

# RiDisplayChannel
def RiDisplayChannel(channel, *paramlist, **keyparams):
    """Defines a new display channel.

    Example: RiDisplayChannel("color aovCi", "string opacity", "aovOi")
    """
    _ribout.write('DisplayChannel "%s"%s\n'%(channel, _paramlist2string(paramlist, keyparams)))

# RiFormat
def RiFormat(xres, yres, aspect):
    """Set the resolution of the output image and the aspect ratio of a pixel.

    Example: RiFormat(720,576,1)"""

    _ribout.write('Format %s %s %s\n'%(xres, yres, aspect))

# RiFrameAspectRatio
def RiFrameAspectRatio(frameratio):
    """Set the ratio between width and height of the image.

    Example: RiFrameAspectRatio(4.0/3)
    """

    _ribout.write('FrameAspectRatio %s\n'%frameratio)

# RiGeometricApproximation
def RiGeometricApproximation(type, value):
    """Sets parameters for approximating surfaces.

    Example: RiGeometricApproximation(RI_FLATNESS, 0.5)  (default value)
    """

    _ribout.write('GeometricApproximation "'+type+'" '+str(value)+"\n")


# RiProjection
def RiProjection(name, *paramlist, **keyparams):
    """Specify a projection method.

    The standard projections are RI_PERSPECTIVE and RI_ORTHOGRAPHIC.
    The perspective projection takes one optional parameter, RI_FOV.

    Example: RiProjection(RI_PERSPECTIVE, fov=45)
    """

    if name==RI_NULL:
        _ribout.write('Projection\n')
    else:
        _ribout.write('Projection "'+name+'"'+ \
                      _paramlist2string(paramlist, keyparams)+"\n")

# RiCamera
def RiCamera(name, *paramlist, **keyparams):
    """Mark the current camera description.

    Example: RiCamera("rightcamera")
    """
    _ribout.write('Camera "%s"%s\n'%(name, _paramlist2string(paramlist, keyparams)))

# RiScreenWindow
def RiScreenWindow(left, right, bottom, top):
    """Specify the extents of the output image on the image plane.

    Example: RiScreenWindow(-1,1,-1,1)
    """

    _ribout.write('ScreenWindow %s %s %s %s\n'%(left, right, bottom, top))

# RiCropWindow
def RiCropWindow(left, right, bottom, top):
    """Specify a subwindow to render.

    The values each lie between 0 and 1.

    Example: RiCropWindow(0.0, 1.0 , 0.0, 1.0)  (renders the entire frame)
             RiCropWindow(0.5, 1.0 , 0.0, 0.5)  (renders the top right quarter)
    """

    _ribout.write('CropWindow %s %s %s %s\n'%(left, right, bottom, top))

# RiPixelSamples
def RiPixelSamples(xsamples, ysamples):
    """Set the sampling rate in horizontal and vertical direction.

    Example: RiPixelSamples(2,2)"""

    _ribout.write("PixelSamples %s %s\n"%(max(1,xsamples), max(1,ysamples)))

# RiPixelVariance
def RiPixelVariance(variance):
    """Limit the acceptable variance in the output value of pixels.

    Example: RiPixelVariance(0.01)"""

    _ribout.write('PixelVariance %s\n'%variance)

# Predefined filter functions:
RiGaussianFilter   = "gaussian"
RiBoxFilter        = "box"
RiTriangleFilter   = "triangle"
RiSincFilter       = "sinc"
RiCatmullRomFilter = "catmull-rom"

# RiPixelFilter
def RiPixelFilter(function, xwidth, ywidth):
    """Set a pixel filter function and its width in pixels.

    Note: Here you can only use one of the predefined filter functions:
    RiGaussianFilter, RiBoxFilter, RiTriangleFilter, RiSincFilter and
    RiCatmullRomFilter.
    Passing a callable as filter function will issue a warning and ignore
    the call.

    Example: RiPixelFilter(RiGaussianFilter, 2.0, 1.0)"""
    
    if hasattr(function, "__call__"):
        _error(RIE_INCAPABLE, RIE_WARNING, "Only the standard filters can be stored in a RIB stream.")
        return

    _ribout.write('PixelFilter "'+function+'" '+str(xwidth)+' '+str(ywidth)+'\n')

# RiExposure
def RiExposure(gain, gamma):
    """Sets the parameters for the output color transformation.

    The transformation is color_out = (color_in*gain)^(1/gamma)

    Example: RiExposure(1.3, 2.2)
    """

    _ribout.write("Exposure %s %s\n"%(gain, gamma))

# RiQuantize
def RiQuantize(type, one, min, max, ditheramplitude):
    """Set the quantization parameters for colors and depth.

    Example: RiQuantize(RI_RGBA, 2048, -1024, 3071, 1.0)
    """

    _ribout.write('Quantize "%s" %s %s %s %s\n'%(type, one, min, max, ditheramplitude))
    

# RiDepthOfField
def RiDepthOfField(fstop, focallength, focaldistance):
    """Set depth of field parameters.

    If fstop is RI_INFINITY depth of field is turned off.

    Example: RiDepthOfField(22,45,1200)
    """

    if fstop==RI_INFINITY:
        _ribout.write('DepthOfField\n')
    else:
        _ribout.write('DepthOfField %s %s %s\n'%(fstop, focallength, focaldistance))

# RiMotionBegin
def RiMotionBegin(*times):
    """Start the definition of a moving primitive.

    You can specify the time values directly or inside a sequence,
    for example, RiMotionBegin(0,1) or RiMotionBegin([0,1]).

    Example: RiMotionBegin(0.0, 1.0)
             RiTranslate(1.0, 0.0, 0.0)
             RiTranslate(1.0, 2.0, 0.0)
             RiMotionEnd()
    """

    global _insidemotion

    if _insidemotion:
        _error(RIE_ILLSTATE, RIE_ERROR, "Motion blocks cannot be nested.")

    _ribout.write('MotionBegin '+_seq2list(times)+'\n')
    _insidemotion = 1

# RiMotionEnd
def RiMotionEnd():
    "Terminates the definition of a moving primitive."

    global _insidemotion

    _ribout.write('MotionEnd\n')
    _insidemotion = 0

# RiShutter
def RiShutter(opentime, closetime):
    """Set the times at which the shutter opens and closes.

    Example: RiShutter(0.1, 0.9)
    """

    _ribout.write('Shutter %s %s\n'%(opentime, closetime))

# RiTranslate
def RiTranslate(*translation):
    """Concatenate a translation onto the current transformation.

    The translation is either given as 3 scalars or a sequence of
    3 scalars.

    Example: RiTranslate(1.2, 4.3, -0.5)
             or
             RiTranslate( (1.2, 4.3, -0.5) )
    """

    # Argument = sequence?
    if len(translation)==1:
        s=_seq2list(translation,3)
        _ribout.write('Translate '+s[1:-1]+"\n")
    # Argument = 3 scalars?
    elif len(translation)==3:
        dx,dy,dz=translation
        _ribout.write('Translate %s %s %s\n'%(dx,dy,dz))
    # Invalid argument size
    else:
        raise TypeError("RiTranslate() only takes a sequence or three scalars as arguments")

# RiRotate
def RiRotate(angle, *axis):
    """Rotate about angle degrees about the given axis.

    The axis is either given as 3 scalars or a sequence of 3 scalars.

    Example: RiRotate(90, 1,0,0)
    """

    # Argument = sequence?
    if len(axis)==1:
        s=_seq2list(axis,3)
        _ribout.write('Rotate %s %s\n'%(angle, s[1:-1]))
    # Argument = 3 scalars?
    elif len(axis)==3:
        ax,ay,az=axis
        _ribout.write('Rotate %s %s %s %s\n'%(angle, ax, ay, az))
    # Invalid argument size
    else:
        raise TypeError("RiRotate() only takes 2 or 4 arguments (%s given)"%(len(axis)+1))


# RiScale
def RiScale(*scaling):
    """Concatenate a scaling onto the current transformation.

    The scaling is either given as 3 scalars or a sequence of 3 scalars.

    Example: RiScale(2,2,2)"""

    # Argument = sequence?
    if len(scaling)==1:
        s=_seq2list(scaling,3)
        _ribout.write('Scale '+s[1:-1]+"\n")
    # Argument = 3 scalars?
    elif len(scaling)==3:
        sx,sy,sz=scaling
        _ribout.write('Scale %s %s %s\n'%(sx, sy, sz))
    # Invalid argument size
    else:
        raise TypeError("RiScale() only takes a sequence or three scalars as arguments")


# RiSkew
def RiSkew(angle, *vecs):
    """Concatenate a skew onto the current transformation.

    angle is given in degrees.
    The two vectors are each given as 3 scalars or a sequence
    of 3 scalars.

    Example: RiSkew(45, 0,1,0, 1,0,0)
    """

    # Argument = two sequences?
    if len(vecs)==2:
        s1=_seq2list(vecs[0],3)
        s2=_seq2list(vecs[1],3)
        _ribout.write('Skew '+str(angle)+" "+s1[1:-1]+" "+s2[1:-1]+"\n")
    # Argument = 6 scalars?
    elif len(vecs)==6:
        dx1,dy1,dz1,dx2,dy2,dz2=vecs
        _ribout.write('Skew %s %s %s %s %s %s %s\n'%(angle, dx1, dy1, dz1, dx2, dy2, dz2))
    # Invalid argument size
    else:
        raise TypeError("RiSkew() only takes 3 or 7 arguments (%s given)"%(len(vecs)+1))


# RiPerspective
def RiPerspective(fov):
    """Concatenate a perspective transformation onto the current transformation.
    
    Example: RiPerspective(45)"""

    _ribout.write('Perspective %s\n'%fov)

# RiIdentity
def RiIdentity():
    """Set the current transformation to the identity.

    Example: RiIdentity()
    """

    _ribout.write('Identity\n')

# RiConcatTransform
def RiConcatTransform(transform):
    """Concatenate a transformation onto the current transformation.

    transform must be a sequence that evaluates to 16 floating point
    values (4x4 matrix).

    Example: RiConcatTransform([2,0,0,0, 0,2,0,0, 0,0,2,0, 0,0,0,1])
             RiConcatTransform([[2,0,0,0], [0,2,0,0], [0,0,2,0], [0,0,0,1]])
    """

    _ribout.write('ConcatTransform '+_seq2list(transform,16)+"\n")

# RiTransform
def RiTransform(transform):
    """Set the current transformation.

    transform must be a sequence that evaluates to 16 floating point
    values (4x4 matrix).

    Example: RiTransform([2,0,0,0, 0,2,0,0, 0,0,2,0, 0,0,0,1])
             RiTransform([[2,0,0,0], [0,2,0,0], [0,0,2,0], [0,0,0,1]])
    """

    _ribout.write('Transform '+_seq2list(transform,16)+"\n")

# RiSides
def RiSides(nsides):
    """Specify the number of visible sides of subsequent surfaces.

    Example: RiSides(1)"""

    if nsides!=1 and nsides!=2:
        _error(RIE_RANGE, RIE_ERROR, "The number of sides (nsides) must be either 1 or 2.")

    _ribout.write('Sides %s\n'%nsides)

# RiOrientation
def RiOrientation(orientation):
    """Set the orientation of subsequent surfaces.

    orientation is either RI_OUTSIDE, RI_INSIDE, RI_LH (left handed)
    or RI_RH (right handed).
    """

    _ribout.write('Orientation "'+orientation+'"\n')

# RiReverseOrientation
def RiReverseOrientation():
    """Causes the current orientation to be toggled.

    Example: RiReverseOrientation()
    """

    _ribout.write('ReverseOrientation\n')

# RiMatte
def RiMatte(onoff):
    """Indicates whether subsequent primitives are matte objects.

    Example: RiMatte(RI_TRUE)
    """
    if onoff:
        _ribout.write('Matte 1\n')
    else:
        _ribout.write('Matte 0\n')

# RiLightSource
def RiLightSource(name, *paramlist, **keyparams):
    """Add another light source and return its light handle.

    name is the name of the light source shader. You can set a user defined
    string handle via the RI_HANDLEID parameter.

    Example: light1 = RiLightSource("distantlight", intensity=1.5)
    """
    global _lighthandle

    paramlist = _merge_paramlist(paramlist, keyparams)
    lshandle = None
    for i in range(0, len(paramlist), 2):
        token = paramlist[i]
        if token==RI_HANDLEID:
            lshandle = str(paramlist[i+1])
            paramlist = paramlist[:i]+paramlist[i+2:]
            break

    # Check if the user provided a handle id...
    if lshandle is None:
        _lighthandle += 1
        lshandle = _lighthandle
        _ribout.write('LightSource "%s" %d%s\n'%(name, lshandle, _paramlist2string(paramlist, {})))
    else:
        _ribout.write('LightSource "%s" "%s"%s\n'%(name, lshandle, _paramlist2string(paramlist, {})))
        
    return lshandle

# RiIlluminate
def RiIlluminate(light, onoff):
    """Activate or deactive a light source.

    Example: RiIlluminate(lgt, RI_TRUE)
    """

    if type(light) is str:
        _ribout.write('Illuminate "%s" %d\n'%(light, onoff))
    else:
        _ribout.write('Illuminate %d %d\n'%(light, onoff))
        

# RiAreaLightSource
def RiAreaLightSource(name, *paramlist, **keyparams):
    """Start the definition of an area light and return the light handle.

    You can set a user defined string handle via the RI_HANDLEID parameter.
    
    Example: RiAttributeBegin()
             area1 = RiAreaLightSource("arealight", intensity=0.75)
             ....
             RiAttributeEnd()
             RiIlluminate(area1, RI_TRUE)
    """
    global _lighthandle

    # Check if the user provided a handle id...
    if RI_HANDLEID in keyparams:
        lshandle = str(keyparams[RI_HANDLEID])
        del keyparams[RI_HANDLEID]
        _ribout.write('AreaLightSource "%s" "%s"%s\n'%(name, lshandle, _paramlist2string((), keyparams)))
    else:
        _lighthandle+=1
        lshandle = _lighthandle
        _ribout.write('AreaLightSource "%s" %d%s\n'%(name, lshandle, _paramlist2string((), keyparams)))

    return lshandle

# RiDeclare
def RiDeclare(name, declaration):
    """Declare the name and type of a variable.

    The syntax of the declaration is:  [class] [type] ['['n']']

    class ::= constant | uniform | varying | vertex
    type  ::= float | integer | string | color | point | vector | normal |
              matrix | hpoint
    
    Example: RiDeclare("foo","uniform float")
             RiDeclare("bar","constant integer [4]")
             RiDeclare("mycolor", "varying color")
    """

    global _declarations

    if declaration==RI_NULL:
        declaration=""
        
    _ribout.write('Declare "'+name+'" "'+declaration+'"\n')
    _declarations[name]=declaration
    return name

# RiArchiveRecord
def RiArchiveRecord(type, format, *args):
    """Output a user data record.

    type is one of RI_COMMENT, RI_STRUCTURE or RI_VERBATIM.

    For comments and structural hints you can use the special variables
    $CREATOR, $DATE and $USER which will be replaced by their appropriate
    value (program name, date string, user name).

    Example: RiArchiveRecord(RI_COMMENT, "Frame %d", 2)
             RiArchiveRecord(RI_STRUCTURE, "CreationDate $DATE")
    """

    if type!=RI_VERBATIM and format.find("$")!=-1:
        format = format.replace("$DATE", time.ctime())
        format = format.replace("$CREATOR", sys.argv[0])
        try:
            user = getpass.getuser()
        except:
            user = "<unknown>"
        format = format.replace("$USER", user)

    if type==RI_COMMENT:
        outstr = "# "+format%args+"\n"
    elif type==RI_STRUCTURE:
        outstr = "##"+format%args+"\n"
    elif type==RI_VERBATIM:
        outstr = format%args
    else:
        return

    # Use the writeArchiveRecord() if there is any, otherwise use write()
    # (the latter case happens when RiBegin() wasn't called. _ribout is
    # then set to stdout)
    if hasattr(_ribout, "writeArchiveRecord"):
        _ribout.writeArchiveRecord(outstr)
    else:
        _ribout.write(outstr)
    
    
# RiReadArchive
def RiReadArchive(filename, callback=None, *ignore):
    """Include an archive file.

    In this implementation the callback function is not used and can
    be left out.

    Example: RiReadArchive("teapot.rib")"""

    _ribout.write('ReadArchive "'+filename+'"\n')

# RiArchiveBegin
def RiArchiveBegin(archivename, *paramlist, **keyparams):
    """Begin an inline archive.
    
    Example: RiArchiveBegin("myarchive")
             ...
             RiArchiveEnd()
             RiReadArchive("myarchive")
    """
    _ribout.write('ArchiveBegin "%s"%s\n'%(archivename, _paramlist2string((), keyparams)))
    return archivename

# RiArchiveEnd
def RiArchiveEnd():
    """Terminate an inline archive.
    
    Example: RiArchiveBegin("myarchive")
             ...
             RiArchiveEnd()
             RiReadArchive("myarchive")
    """
    _ribout.write('ArchiveEnd\n')


def RiProcDelayedReadArchive(): return "DelayedReadArchive"
def RiProcRunProgram(): return "RunProgram"
def RiProcDynamicLoad(): return "DynamicLoad"
def RiProcFree(data): pass

# RiProcedural
def RiProcedural(data, bound, subdividefunc, freefunc=None):
    """Declare a procedural model.

    subdividefunc and freefunc may either be the standard RenderMan
    procedurals (RiProcDelayedReadArchive, RiProcRunProgram,
    RiProcDynamicLoad and RiProcFree) or Python callables.
    In the former case, data must be a sequence of strings or a single
    string containing the data for the functions. In the latter case,
    data may be any Python object which is just passed on to the
    functions. 
    freefunc is optional and defaults to None.

    Because this module can only produce RIB, a custom subdivide function is
    simply called with a detail value of RI_INFINITY to generate all the
    data at once.
    
    Example: RiProcedural("mymodel.rib", [-1,1,-1,1,-1,1], \\
                          RiProcDelayedReadArchive, RI_NULL)
                          
             RiProcedural(["python teapot.py",""],[0,1,0,1,0,1], \\
                          RiProcRunProgram, RI_NULL)
                          
             RiProcedural(["teapot.so",""],[0,1,0,1,0,1], \\
                          RiProcDynamicLoad, RI_NULL)
    """
    if subdividefunc in [RiProcDelayedReadArchive, RiProcRunProgram, RiProcDynamicLoad]:
        if type(data) is str:
            data=[data]
        _ribout.write('Procedural "'+subdividefunc()+'" '+_seq2list(data)+ \
                     ' '+_seq2list(bound,6)+"\n")
    else:
        # Call the custom procedure to generate all the data...
        subdividefunc(data, RI_INFINITY)
        if freefunc is not None:
            freefunc(data)

# RiGeometry
def RiGeometry(type, *paramlist, **keyparams):
    """Create an implementation-specific geometric primitive.

    Example: RiGeometry("teapot")
    """

    _ribout.write('Geometry "'+type+'"'+_paramlist2string(paramlist, keyparams)+"\n")

# RiBound
def RiBound(bound):
    """Set the bounding box for subsequent primitives.

    bound must be a sequence of six floating point values specifying
    the extent of the box along each coordinate direction:
    bound = [xmin, xmax, ymin, ymax, zmin, zmax]

    Example: RiBound([-1,1, 0,1, 0.5,0.75])
    """

    _ribout.write('Bound '+_seq2list(bound,6)+'\n')
    
# RiSolidBegin
def RiSolidBegin(type):
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

    _ribout.write('SolidBegin "'+type+'"\n')

# RiSolidEnd
def RiSolidEnd():
    """Terminate the definition of a solid object."""

    _ribout.write('SolidEnd\n')

# RiObjectBegin
def RiObjectBegin(*paramlist, **keyparams):
    """Start the definition of a retained model and return the object handle.

    You can pass a user defined string handle via the RI_HANDLEID parameter.

    Example: obj1 = RiObjectBegin()
             ...
             RiObjectEnd()
    """
    global _objecthandle, _insideobject

    if _insideobject:
        _error(RIE_ILLSTATE, RIE_ERROR, "Object blocks cannot be nested.")

    keyparams = _paramlist2dict(paramlist, keyparams)

    # Check if the user provided a handle id...
    if RI_HANDLEID in keyparams:
        objhandle = str(keyparams[RI_HANDLEID])
        del keyparams[RI_HANDLEID]
        _ribout.write('ObjectBegin "%s"\n'%objhandle)
    else:
        _objecthandle+=1
        objhandle = _objecthandle
        _ribout.write('ObjectBegin %d\n'%objhandle)
    
    _insideobject=1

    return objhandle

# RiObjectEct
def RiObjectEnd():
    """Terminate the definition of a retained model."""

    global _insideobject
    
    _ribout.write('ObjectEnd\n')
    _insideobject=0

# RiObjectInstance
def RiObjectInstance(handle):
    """Create an instance of a previously defined model.

    Example: RiObjectInstance(obj1)
    """

    if type(handle) is str:
        _ribout.write('ObjectInstance "%s"\n'%handle)
    else:
        _ribout.write('ObjectInstance %d\n'%handle)

# RiTextureCoordinates
def RiTextureCoordinates(s1, t1, s2, t2, s3, t3, s4, t4):
    """Set the current set of texture coordinates.

    Declares a projection from the unit square [(0,0), (1,0), (0,1), (1,1)]
    in parameter space to quadrilateral [(s1,t1), (s2,t2), (s3,t3), (s4,t4)]
    in texture space.

    Example: RiTextureCoordinates(0.0, 0.0, 2.0, -0.5, -0.5, 1.75, 3.0, 3.0)"""

    _ribout.write('TextureCoordinates %s %s %s %s %s %s %s %s\n'%(s1,t1,s2,t2,s3,t3,s4,t4))

# RiMakeTexture
def RiMakeTexture(picname, texname, swrap, twrap, filterfunc, swidth, twidth, *paramlist, **keyparams):
    """Convert an image file into a texture file.

    swrap and twrap are one of RI_PERIODIC, RI_CLAMP or RI_BLACK.
    filterfunc has to be one of the predefined filter functions:
    RiGaussianFilter, RiBoxFilter, RiTriangleFilter, RiSincFilter or
    RiCatmullRomFilter (otherwise a warning is issued and the call is ignored).
    swidth and twidth define the support of the filter.

    Example: RiMakeTexture("img.tif", "tex.tif", RI_PERIODIC, RI_CLAMP, \\
                           RiGaussianFilter, 2,2)
    """
    
    if hasattr(filterfunc, "__call__"):
        _error(RIE_INCAPABLE, RIE_WARNING, "Only the standard filters can be stored in a RIB stream.")
        return

    _ribout.write('MakeTexture "'+picname+'" "'+texname+'" "'+swrap+'" "'+
                  twrap+'" "'+filterfunc+'" '+str(swidth)+' '+str(twidth)+
                  _paramlist2string(paramlist, keyparams)+'\n')

# RiMakeLatLongEnvironment
def RiMakeLatLongEnvironment(picname, texname, filterfunc, swidth, twidth, *paramlist, **keyparams):
    """Convert an image file into an environment map.

    filterfunc has to be one of the predefined filter functions:
    RiGaussianFilter, RiBoxFilter, RiTriangleFilter, RiSincFilter or
    RiCatmullRomFilter (otherwise a warning is issued and the call is ignored).
    swidth and twidth define the support of the filter.

    Example: RiMakeLatLongEnvironment("img.tif", "tex.tif",
                                      RiGaussianFilter, 2,2)
    """
    
    if hasattr(filterfunc, "__call__"):
        _error(RIE_INCAPABLE, RIE_WARNING, "Only the standard filters can be stored in a RIB stream.")
        return

    _ribout.write('MakeLatLongEnvironment "'+picname+'" "'+texname+'" "'+
                  filterfunc+'" '+str(swidth)+' '+str(twidth)+
                  _paramlist2string(paramlist, keyparams)+'\n')

# RiMakeCubeFaceEnvironment
def RiMakeCubeFaceEnvironment(px,nx,py,ny,pz,nz, texname, fov, filterfunc, swidth, twidth, *paramlist, **keyparams):
    """Convert six image files into an environment map.

    The px/nx images are the views in positive/negative x direction.
    fov is the field of view that was used to generate the individual images.
    filterfunc has to be one of the predefined filter functions:
    RiGaussianFilter, RiBoxFilter, RiTriangleFilter, RiSincFilter or
    RiCatmullRomFilter (otherwise a warning is issued and the call is ignored).
    swidth and twidth define the support of the filter.

    Example: RiMakeCubeFaceEnvironment("px.tif","nx.tif","py.tif","ny.tif",
                                       "pz.tif","nz.tif", "tex.tif", 92.0,
                                        RiGaussianFilter, 2,2)
    """
    
    if hasattr(filterfunc, "__call__"):
        _error(RIE_INCAPABLE, RIE_WARNING, "Only the standard filters can be stored in a RIB stream.")
        return

    _ribout.write('MakeCubeFaceEnvironment "'+px+'" "'+nx+'" "'+
                  py+'" "'+ny+'" "'+pz+'" "'+nz+'" "'+ texname+'" '+
                  str(fov)+' "'+filterfunc+'" '+str(swidth)+' '+str(twidth)+
                  _paramlist2string(paramlist, keyparams)+'\n')

# RiMakeShadow
def RiMakeShadow(picname, shadowname, *paramlist, **keyparams):
    """Transform a depth image into a shadow map.

    Example: RiMakeShadow("depthimg.tif", "shadow.tif")
    """
    
    _ribout.write('MakeShadow "'+picname+'" "'+shadowname+'"'+_paramlist2string(paramlist, keyparams)+'\n')

# RiMakeBrickMap
def RiMakeBrickMap(ptcnames, bkmname, *paramlist, **keyparams):
    """Create a brick map file from a list of point cloud file names.

    Example: RiMakeBrickMap(["sphere.ptc", "box.ptc"], "spherebox.bkm", "float maxerror", 0.002)
    """
    names = " ".join(['"%s"'%name for name in ptcnames])
    _ribout.write('MakeBrickMap [%s] "%s"%s\n'%(names, bkmname, _paramlist2string(paramlist, keyparams)))

# RiDetail
def RiDetail(bound):
    """Set the current bounding box.

    bound must be a sequence of six floating point values specifying
    the extent of the box along each coordinate direction:
    bound = [xmin, xmax, ymin, ymax, zmin, zmax]

    Example: RiDetail([10,20,40,70,0,1])
    """

    _ribout.write('Detail '+_seq2list(bound,6)+'\n')

# RiRelativeDetail
def RiRelativeDetail(relativedetail):
    """Set the factor for all level of detail calculations.

    Example: RiRelativeDetail(0.7)"""

    _ribout.write('RelativeDetail %s\n'%relativedetail)

# RiDetailRange
def RiDetailRange(minvisible, lowertransition, uppertransition, maxvisible):
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

    _ribout.write('DetailRange %s %s %s %s\n'%(minvisible, lowertransition, uppertransition, maxvisible))

# RiCoordinateSystem
def RiCoordinateSystem(spacename):
    """Mark the current coordinate system with a name.

    Example: RiCoordinateSystem("lamptop")
    """
    _ribout.write('CoordinateSystem "'+spacename+'"\n')

# RiScopedCoordinateSystem
def RiScopedCoordinateSystem(spacename):
    """Mark the current coordinate system with a name but store it on a separate stack.

    Example: RiScopedCoordinateSystem("lamptop")
    """
    _ribout.write('ScopedCoordinateSystem "'+spacename+'"\n')

# RiTransformPoints
def RiTransformPoints(fromspace, tospace, points):
    """Transform a set of points from one space to another.

    This function is not implemented and always returns None.
    """

    return None

# RiCoordSysTransform
def RiCoordSysTransform(spacename):
    """Replace the current transformation matrix with spacename.

    Example: RiCoordSysTransform("lamptop")
    """

    _ribout.write('CoordSysTransform "'+spacename+'"\n')

# RiContext
def RiContext(handle):
    """Set the current active rendering context.

    Example: ctx1 = RiGetContext()
             ...
             RiContext(ctx1)
    """

    _switch_context(handle)

# RiGetContext
def RiGetContext():
    """Get a handle for the current active rendering context.

    Example: ctx1 = RiGetContext()
             ...
             RiContext(ctx1)"""

    global _current_context

    return _current_context

# RiSystem
def RiSystem(cmd):
    """Execute an arbitrary command in the same environment as the current rendering pass.
    """
    # Escape quotes
    cmd = cmd.replace('"', r'\"')
    _ribout.write('System "%s"\n'%cmd)

# RiIfBegin
def RiIfBegin(expression, *paramlist, **keyparams):
    """Begin a conditional block.
    """
    _ribout.write('IfBegin "%s"%s\n'%(expression, _paramlist2string(paramlist, keyparams)))

# RiElseIf
def RiElseIf(expression, *paramlist, **keyparams):
    """Add an else-if block to a conditional block.
    """
    
    _ribout.write('ElseIf "%s"%s\n'%(expression, _paramlist2string(paramlist, keyparams)))

# RiElse
def RiElse():
    """Add an else block to a conditional block.
    """
    
    _ribout.write('Else\n')

# RiIfEnd
def RiIfEnd():
    """Terminate a conditional block.
    """
    
    _ribout.write('IfEnd\n')
    
# RiResource
def RiResource(handle, type, *paramlist, **keyparams):
    """Create or operate on a named resource of a particular type.
    """
    _ribout.write('Resource "%s" "%s"%s\n'%(handle, type, _paramlist2string(paramlist, keyparams)))

# RiResourceBegin
def RiResourceBegin():
    """Push the current set of resources. 
    """
    _ribout.write('ResourceBegin\n')

# RiResourceEnd
def RiResourceEnd():
    """Pop the current set of resources. 
    """
    _ribout.write('ResourceEnd\n')

##################### Global variabels (internal) ####################

_contexts     = {}
_current_context = None

# The number of digits that floats are rounded to (may be negative).
# This is the second argument to round(). The parameter is global and not
# part of a Ri context.
_round_ndigits = 10
# The conversion format string for floats. Floats are first rounded and
# then converted using this string.
_float_conversion_string = "%1.6g"

####

# Initially the output stream is stdout (and not an instance of RIBStream)
# In interactive sessions this prevents the version number to be written.
_ribout       = sys.stdout
_colorsamples = 3
_lighthandle  = 0
_objecthandle = 0
_errorhandler = RiErrorPrint
_declarations = {}

_insideframe  = 0
_insideworld  = 0
_insideobject = 0
_insidesolid  = 0
_insidemotion = 0

# If you're adding new global variables then make sure that they're
# saved and loaded from the context handling functions and initialized
# in RiBegin() and during the module initialization.

##################### Internal helper functions #######################

def _save_context(handle):
    "Save a context."
    ctx = (_ribout, _colorsamples, _lighthandle, _objecthandle,
           _errorhandler, _declarations,
           _insideframe, _insideworld, _insideobject, _insidesolid,
           _insidemotion)
    _contexts[handle]=ctx

def _load_context(handle):
    "Load a context."
    global _ribout, _colorsamples, _lighthandle, _objecthandle
    global _errorhandler, _declarations, _insideframe, _insideworld
    global _insideobject, _insidesolid, _insidemotion

    _ribout, _colorsamples, _lighthandle, _objecthandle, \
    _errorhandler, _declarations, \
    _insideframe, _insideworld, _insideobject, _insidesolid, \
    _insidemotion = _contexts[handle]

def _create_new_context():
    "Create a new context and make it the active one."
    global _current_context

    keys = list(_contexts.keys())
    if len(keys)>0:
        handle = max(keys)+1
    else:
        handle = 1
    _contexts[handle]=()
    
    if _current_context!=None:
        _save_context(_current_context)

    _current_context = handle

def _switch_context(handle):
    "Save the current context and make another context the active one."
    global _current_context
    
    if _current_context!=None:
        _save_context(_current_context)
    _current_context = handle
    _load_context(handle)

def _destroy_context():
    "Destroy the current active context"
    global _contexts, _current_context

    handle = _current_context
    del _contexts[handle]
    _current_context = None

def _init_declarations():
    global _declarations
    _declarations = {RI_P:"vertex point", RI_PZ:"vertex point",
                     RI_PW:"vertex hpoint",
                     RI_N:"varying normal", RI_NP:"uniform normal",
                     RI_CS:"varying color", RI_OS:"varying color",
                     RI_S:"varying float", RI_T:"varying float",
                     RI_ST:"varying float[2]",
                     RI_ORIGIN:"integer[2]",
                     RI_KA:"uniform float",
                     RI_KD:"uniform float",
                     RI_KS:"uniform float",
                     RI_ROUGHNESS:"uniform float",
                     RI_KR:"uniform float",
                     RI_TEXTURENAME:"string",
                     RI_SPECULARCOLOR:"uniform color",
                     RI_INTENSITY:"float",
                     RI_LIGHTCOLOR:"color",
                     RI_FROM:"point",
                     RI_TO:"point",
                     RI_CONEANGLE:"float",
                     RI_CONEDELTAANGLE:"float",
                     RI_BEAMDISTRIBUTION:"float",
                     RI_AMPLITUDE:"uniform float",
                     RI_MINDISTANCE:"float",
                     RI_MAXDISTANCE:"float",
                     RI_BACKGROUND:"color",
                     RI_DISTANCE:"float",
                     RI_FOV:"float",
                     RI_WIDTH:"varying float",
                     RI_CONSTANTWIDTH:"constant float",
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

def _error(code, severity, message):
    global _errorhandler, RiLastError

    RiLastError = code

    st = inspect.stack(1)
    # Search the offending Ri function in the stack...
    j=None
    for i in range(len(st)):
        if st[i][3][:2]=="Ri":
            j=i
    # No function beginning with "Ri" found? That's weird. Maybe someone
    # messed with the function names.
    if j==None:
        where=""
    else:
        # name of the Ri function
        call   = inspect.stack(0)[j][3]
        # filename and line number where the offending Ri call occured
        file   = inspect.stack(1)[j+1][1]
        line   = inspect.stack(1)[j+1][2]
        if file==None: file="?"
        where = 'In file "%s", line %s - %s():\n'%(file, line, call)

    _errorhandler(code,severity,where+message)

def _seq2col(seq):
    """Convert a sequence containing a color into a string."""
    if len(seq)<_colorsamples:
        _error(RIE_INVALIDSEQLEN, RIE_ERROR, "Invalid sequence length (%s instead of %s)"%(len(seq), _colorsamples))
    colseq = tuple(seq)
    return '[%s]'%(" ".join( [str(x) for x in colseq[:_colorsamples]] ))

def _flatten(seq):
    """Return a list of the individual items in a (possibly nested) sequence.

    Returns a list with all items as strings.
    If an item was already a string it's enclosed in apostrophes.
    
    Example: _flatten( [(1,2,3), (4,5,6)] ) -> ["1","2","3","4","5","6"]
             _flatten( ("str1","str2") )    -> ['"str1"','"str2"']
    """
    global _round_ndigits
    global _float_conversion_string
    
    ndigits = _round_ndigits
    floatFmtStr = _float_conversion_string
    
    res = []
    for v in seq:
        vtype = type(v)
        # v=scalar?
        if vtype is float:
            res.append(floatFmtStr%round(v,ndigits))
        elif vtype is int:
            res.append(str(v))
        # vec3?
        elif isinstance(v, _vec3):
            res.extend([floatFmtStr%round(v.x,ndigits),
                        floatFmtStr%round(v.y,ndigits),
                        floatFmtStr%round(v.z,ndigits)])
        # v=string?
        elif type(v) is str:
            res.append('"%s"'%v)
        # no scalar or string. Then it might be a sequence...
        else:
            # Check if it is really a sequence...
            try:
                n = len(v)
            except:
                res.append(str(v))
                continue
            res += _flatten(v)
    return res

def _seq2list(seq, count=None):
    """Convert a sequence into a string.

    The function checks if the sequence contains count elements (unless
    count is None). If it doesn't an error is generated.
    The return value is a string containing the sequence. The string can
    be used as parameter value to RIB commands.
    """

    f = _flatten(seq)
    # Has the sequence an incorrect length? then generate an error
    if count!=None and len(f)!=count:
        _error(RIE_INVALIDSEQLEN, RIE_ERROR, "Invalid sequence length (%s instead of %s)"%(len(f), count))
        
    return '[%s]'%" ".join(f)

def _paramlist2dict(paramlist, keyparams):
    """Combine the paramlists (tuple & dict) into one dict.
    
    paramlist is a tuple with function arguments (token/value pairs or
    a dictionary). keyparams is a dictionary with keyword arguments.
    The dictionary keyparams will be modified and returned.
    """

    if len(paramlist)==1 and type(paramlist[0]) is dict:
        keyparams = paramlist[0]
        paramlist = ()
    
    # Add the paramlist tuple to the keyword argument dict
    for i in range(len(paramlist)//2):
        token = paramlist[i*2]
        value = paramlist[i*2+1]
        keyparams[token]=value

    return keyparams

def _paramlist2lut(paramlist, keyparams):
    """Combine the paramlists into one dict without inline declarations.

    paramlist is a tuple with function arguments. keyparams is a
    dictionary with keyword arguments. The dictionary keyparams will
    be modified and returned.

    The resulting dictionary can be used to look up the value of tokens.
    """
    # Add the paramlist tuple to the keyword argument dict
    for i in range(len(paramlist)//2):
        token = paramlist[i*2]
        value = paramlist[i*2+1]
        # Extract the name of the token (without inline declaration
        # if there is one)
        f = token.split(" ")
        tokname = f[-1]
        keyparams[tokname]=value

    return keyparams
    
def _merge_paramlist(paramlist, keyparams):
    """Merge a paramlist tuple and keyparams dict into one single list.
    """
    if len(paramlist)==1 and type(paramlist[0]) is dict:
        keyparams = paramlist[0]
        paramlist = ()

    res = list(paramlist)
    # Check if the number of values is uneven (this is only allowed when
    # the last value is None (RI_NULL) in which case this last value is ignored)
    if (len(res)%2==1):
        if res[-1] is None:
            res = res[:-1]
        else:
            raise ValueError("The parameter list must contain an even number of values")

    # Append the params from the keyparams dict to the parameter list
    for param,value in sorted(keyparams.items()):
        res.extend([param,value])
    return res
    

def _paramlist2string(paramlist, keyparams={}):
    """Convert the paramlist into a string representation.

    paramlist is a tuple with function arguments (token/value pairs or
    a dictionary). keyparams is a dictionary with keyword arguments.
    Each token has to be a string, the value can be of any type. If the
    value is a string, it's enclosed in apostrophes. A trailing token
    without a value is ignored, which also means that a trailing RI_NULL
    can be passed.
    The resulting string contains a leading space, unless there are no
    token/value pairs.
    """

    global _declarations
    global _round_ndigits
    global _float_conversion_string

    paramlist = _merge_paramlist(paramlist, keyparams)

    res=""
    for i in range(0, len(paramlist), 2):
        token = paramlist[i].strip()
        value = paramlist[i+1]
        # Extract the name of the token (without inline declaration
        # if there is one)
        f = token.split(" ")
        tokname = f[-1:][0]
        inline  = f[:-1]

        if not (tokname in _declarations or inline!=[]):
            _error(RIE_UNDECLARED,RIE_ERROR,'Parameter "'+tokname+
                   '" is not declared.')
        
        # Check if the value is a sequence (if it returns an iterator)
        isseq=0
        try:
            isseq = (iter(value)!=None)
        except:
            pass
        # Convert value into the appropriate string representation
        if type(value) is str:
            value='["'+value+'"]'
#        elif type(value)==types.ListType or type(value)==types.TupleType:
        elif isseq:
            value = _seq2list(value)
        else:
            if type(value) is float:
                value = '[%s]'%(_float_conversion_string%round(value, _round_ndigits))
            else:
                value='[%s]'%value
        res+=' "%s" %s'%(token, value)

    if (res==" "): res=""
    return res


############################################################

_init_declarations()

if __name__=='__main__':

    RiBegin(RI_NULL)
    RiErrorHandler(RiErrorAbort)

    RiWorldBegin();

    RiWorldEnd()

    RiSkew(45,0,1,0,1,0,0)
    RiSkew(45,[0,1,0],[1,0,0])
    
    RiEnd()
    
