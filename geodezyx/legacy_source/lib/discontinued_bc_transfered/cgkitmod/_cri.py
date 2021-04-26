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

"""This is the "lower-level" _cri module which is used by the cri module.

This module provides the function loadRI() which loads a shared RenderMan
library and prepares the returned library handle, so that all RenderMan
constants, tokens and functions are properly defined.

The ri library returned by loadRI() has to be used like a C library and
can crash the application if not used properly (like forgetting to terminate
parameter lists with RI_NULL). This is why this low-level ri library is not
passed back to the user but instead the cri module provides another wrapper
layer that allows the library to be used just like the pure Python ri module.
"""

import os.path
from ctypes import *
from . import rmanlibutil

def loadRI(libName):
    """Load a RenderMan library and return a module-like handle to it.
    
    libName is the name of a shared library that implements the RenderMan
    interface. The name can either be an absolute file name or just the
    name of the library (without suffix or "lib" prefix) in which case 
    the function tries to find the library file itself.
    The return value is the library handle as returned by the ctypes 
    LoadLibrary() function. This handle has already been prepared so that
    it can be used like a module that contains the RenderMan interface.
    """
    
    # Try to figure out the location of the lib if the name is not an absolute path...
    libName = rmanlibutil.resolveRManLib(libName)
        
    # Load the library...
    ri = cdll.LoadLibrary(libName)
    
    _createRiTypes(ri)
    _createRiTokens(ri)
    _createRiConstants(ri)
    _createRiFunctions(ri)

    # Create an alias for every Ri function and RI_ constant that has the prefix removed...
    for name in dir(ri):
        if name.startswith("Ri"):
            setattr(ri, name[2:], getattr(ri, name))
        elif name.startswith("RI_"):
            setattr(ri, name[3:], getattr(ri, name))

    ri.__class__.RiLastError = property(_getLastError, _setLastError)
    ri.__class__.LastError = property(_getLastError, _setLastError)

    return ri

def importRINames(ri, ns):
    """Import the RenderMan names into the given namespace.
    """
    for name in dir(ri):
        if name.startswith("Ri"):
            localName = name
            ns[localName] = getattr(ri, name)
            
        if name[:2] in ["Rt", "RI"]:
            ns[name] = getattr(ri, name)


def _createRiTypes(ri):
    """Create the RenderMan types.
    
    ri must be the open ctypes library handle.
    The RenderMan types are added as attributes to the ri object. All names
    begin with "Rt" (RtInt, RtFloat, ...).
    
    This is a helper function for loadRI().
    """
    
    # Base types that are not composed of other RenderMan types...
    baseTypes = dict(RtBoolean = c_short,
                     RtInt = c_int,
                     RtFloat = c_float,
                     RtString = c_char_p,
                     RtToken = c_char_p,
                     RtVoid = None,
                     RtPointer = c_void_p)
    
    for name,ctype in list(baseTypes.items()):
        setattr(ri, name, ctype)
        
    ri.RtColor = 3*ri.RtFloat
    ri.RtPoint = 3*ri.RtFloat
    ri.RtVector = 3*ri.RtFloat
    ri.RtNormal = 3*ri.RtFloat
    ri.RtHpoint = 4*ri.RtFloat
    ri.RtMatrix = 16*ri.RtFloat
    ri.RtBasis = 16*ri.RtFloat
    ri.RtBound = 6*ri.RtFloat
    
    ri.RtObjectHandle = ri.RtPointer
    ri.RtLightHandle = ri.RtPointer
    ri.RtContextHandle = ri.RtPointer
    
    ri.RtFilterFunc = CFUNCTYPE(ri.RtFloat,  ri.RtFloat, ri.RtFloat, ri.RtFloat, ri.RtFloat)
    ri.RtErrorHandler = CFUNCTYPE(ri.RtVoid,  ri.RtInt, ri.RtInt, c_char_p)
    ri.RtProcSubdivFunc = CFUNCTYPE(ri.RtVoid,  ri.RtPointer, ri.RtFloat)
    ri.RtProcFreeFunc = CFUNCTYPE(ri.RtVoid,  ri.RtPointer)
    ri.RtArchiveCallback = CFUNCTYPE(ri.RtVoid,  ri.RtToken, c_char_p)   # var args are missing
    
    
def _createRiConstants(ri):
    """Create the RenderMan constants.
    
    Add the RenderMan constants to the ri object which must be an open
    ctypes library handle.
    The types must already be available on ri (i.e. _createRiTypes() must
    already have been called).
    
    This is a helper function for loadRI().
    """
    
    ri.RI_NULL         = None
    ri.RI_TRUE         = 1
    ri.RI_FALSE        = 0
    ri.RI_EPSILON      = 1.0e-10
    ri.RI_INFINITY     = ri.RtFloat(1.0e38).value

    ri.RI_HERMITESTEP    = 2
    ri.RI_CATMULLROMSTEP = 1
    ri.RI_BEZIERSTEP     = 3
    ri.RI_BSPLINESTEP    = 1
    ri.RI_POWERSTEP      = 4

    bases = ["RiBezierBasis", "RiBSplineBasis", "RiCatmullRomBasis",
             "RiHermiteBasis", "RiPowerBasis"]
    
    for basis in bases:
        try:
            value = ri.RtBasis.in_dll(ri, basis)
        except ValueError:
            raise ValueError('The RenderMan implementation "%s" does not define the standard basis "%s"'%(getattr(ri, "_name", "?"), basis))
        setattr(ri, basis, value)
    
def _createRiTokens(ri):
    """Create the RenderMan tokens.
    
    The RenderMan constants are added as attributes to the ri object which
    must be an open ctypes library handle. All names begin with "RI_".
    
    This is a helper function for loadRI().
    """

    # All the following constants will be added to the ri object.
    # The items are (Token Name, Default). If a token does not exist in the
    # library, the default value is used.
    tokens = [("RI_A", "a"),
              ("RI_ABORT", "abort"),
              ("RI_AMBIENTLIGHT", "ambientlight"),
              ("RI_AMPLITUDE", "amplitude"),
              ("RI_AZ", "az"),
              ("RI_BACKGROUND", "background"),
              ("RI_BEAMDISTRIBUTION", "beamdistribution"),
              ("RI_BICUBIC", "bicubic"),
              ("RI_BILINEAR", "bilinear"),
              ("RI_BLACK", "black"),
              ("RI_BUMPY", "bumpy"),
              ("RI_CAMERA", "camera"),
              ("RI_CI", "Ci"),
              ("RI_CLAMP", "clamp"),
              ("RI_COMMENT", "comment"),
              ("RI_CONEANGLE", "coneangle"),
              ("RI_CONEDELTAANGLE", "conedeltaangle"),
              ("RI_CONSTANT", "constant"),
              ("RI_CONSTANTWIDTH", "constantwidth"),
              ("RI_CS", "Cs"),
              ("RI_CUBIC", "cubic"),
              ("RI_DEPTHCUE", "depthcue"),
              ("RI_DIFFERENCE", "difference"),
              ("RI_DISTANCE", "distance"),
              ("RI_DISTANTLIGHT", "distantlight"),
              ("RI_FILE", "file"),
              ("RI_FLATNESS", "flatness"),
              ("RI_FOG", "fog"),
              ("RI_FOV", "fov"),
              ("RI_FRAMEBUFFER", "framebuffer"),
              ("RI_FROM", "from"),
              ("RI_HANDLEID", "__handleid"),
              ("RI_HANDLER", "handler"),
              ("RI_HIDDEN", "hidden"),
              ("RI_IDENTIFIER", "identifier"),
              ("RI_IGNORE", "ignore"),
              ("RI_INSIDE", "inside"),
              ("RI_INTENSITY", "intensity"),
              ("RI_INTERSECTION", "intersection"),
              ("RI_KA", "Ka"),
              ("RI_KD", "Kd"),
              ("RI_KR", "Kr"),
              ("RI_KS", "Ks"),
              ("RI_LH", "lh"),
              ("RI_LIGHTCOLOR", "lightcolor"),
              ("RI_LINEAR", "linear"),
              ("RI_MATTE", "matte"),
              ("RI_MAXDISTANCE", "maxdistance"),
              ("RI_METAL", "metal"),
              ("RI_MINDISTANCE", "mindistance"),
              ("RI_N", "N"),
              ("RI_NAME", "name"),
              ("RI_NG", "Ng"),
              ("RI_NONPERIODIC", "nonperiodic"),
              ("RI_NP", "Np"),
              ("RI_OBJECT", "object"),
              ("RI_OI", "Oi"),
              ("RI_ORIGIN", "origin"),
              ("RI_ORTHOGRAPHIC", "orthographic"),
              ("RI_OS", "Os"),
              ("RI_OUTSIDE", "outside"),
              ("RI_P", "P"),
              ("RI_PAINT", "paint"),
              ("RI_PAINTEDPLASTIC", "paintedplastic"),
              ("RI_PERIODIC", "periodic"),
              ("RI_PERSPECTIVE", "perspective"),
              ("RI_PLASTIC", "plastic"),
              ("RI_POINTLIGHT", "pointlight"),
              ("RI_PRIMITIVE", "primitive"),
              ("RI_PRINT", "print"),
              ("RI_PW", "Pw"),
              ("RI_PZ", "Pz"),
              ("RI_RASTER", "raster"),
              ("RI_RGB", "rgb"),
              ("RI_RGBA", "rgba"),
              ("RI_RGBAZ", "rgbaz"),
              ("RI_RGBZ", "rgbz"),
              ("RI_RH", "rh"),
              ("RI_ROUGHNESS", "roughness"),
              ("RI_S", "s"),
              ("RI_SCREEN", "screen"),
              ("RI_SHADINGGROUP", "shadinggroup"),
              ("RI_SHINYMETAL", "shinymetal"),
              ("RI_SMOOTH", "smooth"),
              ("RI_SPECULARCOLOR", "specularcolor"),
              ("RI_SPOTLIGHT", "spotlight"),
              ("RI_ST", "st"),
              ("RI_STRUCTURE", "structure"),
              ("RI_T", "t"),
              ("RI_TEXTURENAME", "texturename"),
              ("RI_TO", "to"),
              ("RI_UNION", "union"),
              ("RI_VERBATIM", "verbatim"),
              ("RI_WIDTH", "width"),
              ("RI_WORLD", "world"),
              ("RI_Z", "z"),
              ]
    
    for tokSpec in tokens:
        if type(tokSpec) is str:
            name = tokSpec
            default = None
        else:
            try:
                name,default = tokSpec
                if type(name) is not str:
                    raise ValueError()
            except:
                raise ValueError("Invalid token spec: %s"%(tokSpec,))
            
        try:
            value = c_char_p.in_dll(ri, name).value
            # Turn into a unicode string
            value = value.decode("ascii")
        except ValueError:
            value = default
            
        setattr(ri, name, value)

        
def _createRiFunctions(ri):
    """Declare the RenderMan functions.
    
    ri must be an open ctypes library handle.
    
    This is a helper function for loadRI().
    """
    RtToken = ri.RtToken
    RtBoolean = ri.RtBoolean
    RtInt = ri.RtInt
    RtFloat = ri.RtFloat
    RtString = ri.RtString
    RtMatrix = ri.RtMatrix
    RtPoint = ri.RtPoint
    RtColor = ri.RtColor
    RtBasis = ri.RtBasis
    RtBound = ri.RtBound
    RtContextHandle = ri.RtContextHandle
    RtLightHandle = ri.RtLightHandle
    RtObjectHandle = ri.RtObjectHandle
    RtPointer = ri.RtPointer
    
    ri.RiArchiveRecord.argtypes = [RtToken, c_char_p]
    ri.RiAreaLightSource.argtypes = [RtToken]
    ri.RiAreaLightSource.restype = RtLightHandle 
    ri.RiAtmosphere.argtypes = [RtToken]
    ri.RiAttribute.argtypes = [RtToken]
    ri.RiAttributeBegin.argtypes = []
    ri.RiAttributeEnd.argtypes = []
    ri.RiBasis.argtypes = [RtBasis, RtInt, RtBasis, RtInt]
    ri.RiBegin.argtypes = [RtToken]
    ri.RiBlobby.argtypes = [RtInt, RtInt, POINTER(RtInt), RtInt, POINTER(RtFloat), RtInt, POINTER(RtString)]
    ri.RiBound.argtypes = [RtBound]
    ri.RiClipping.argtypes = [RtFloat, RtFloat]
    ri.RiClippingPlane.argtypes = [RtFloat, RtFloat, RtFloat, RtFloat, RtFloat, RtFloat]
    ri.RiColor.argtypes = [RtColor]
    ri.RiColorSamples.argtypes = [RtInt, POINTER(RtFloat), POINTER(RtFloat)]
    ri.RiConcatTransform.argtypes = [RtMatrix]
    ri.RiCone.argtypes = [RtFloat, RtFloat, RtFloat]
    ri.RiContext.argtypes = [RtContextHandle]
    ri.RiCoordinateSystem.argtypes = [RtToken]
    ri.RiCoordSysTransform.argtypes = [RtToken]
    ri.RiCropWindow.argtypes = [RtFloat, RtFloat, RtFloat, RtFloat]
    ri.RiCurves.argtypes = [RtToken, RtInt]
    ri.RiCylinder.argtypes = [RtFloat, RtFloat, RtFloat, RtFloat]
    ri.RiDeclare.argtypes = [c_char_p, c_char_p]
    ri.RiDeclare.restype = RtToken
    ri.RiDepthOfField.argtypes = [RtFloat, RtFloat, RtFloat]
    ri.RiDetail.argtypes = [RtBound]
    ri.RiDetailRange.argtypes = [RtFloat, RtFloat, RtFloat, RtFloat]
    ri.RiDisk.argtypes = [RtFloat, RtFloat, RtFloat]
    ri.RiDisplacement.argtypes = [RtToken]
    ri.RiDisplay.argtypes = [RtToken, RtToken, RtToken]
    ri.RiEnd.argtypes = []
    # See RiPixelFilter for an explanation why RiErrorHandler is commented out
#    ri.RiErrorHandler.argtypes = [RtErrorHandler]
    ri.RiExposure.argtypes = [RtFloat, RtFloat]
    ri.RiExterior.argtypes = [RtToken]
    ri.RiFormat.argtypes = [RtInt, RtInt, RtFloat]
    ri.RiFrameAspectRatio.argtypes = [RtFloat]
    ri.RiFrameBegin.argtypes = [RtInt]
    ri.RiFrameEnd.argtypes = []
    ri.RiGeneralPolygon.argtypes = [RtInt, POINTER(RtInt)]
    ri.RiGeometricApproximation.argtypes = [RtToken, RtFloat]
    ri.RiGeometry.argtypes = [RtToken]
    ri.RiGetContext.argtypes = []
    ri.RiGetContext.restype = RtContextHandle
    ri.RiHider.argtypes = [RtToken]
    ri.RiHyperboloid.argtypes = [RtPoint, RtPoint, RtFloat]
    ri.RiIdentity.argtypes = []
    ri.RiIlluminate.argtypes = [RtLightHandle, RtBoolean]
    ri.RiImager.argtypes = [RtToken]
    ri.RiInterior.argtypes = [RtToken]
    ri.RiLightSource.argtypes = [RtToken]
    ri.RiLightSource.restype = RtLightHandle
    # In the following texture calls the declaration of the filter function is removed
    # (see RiPixelFilter for more infos)
#    ri.RiMakeCubeFaceEnvironment.argtypes = [c_char_p, c_char_p, c_char_p, c_char_p, c_char_p, c_char_p, c_char_p, RtFloat, RtFilterFunc, RtFloat, RtFloat]
    ri.RiMakeCubeFaceEnvironment.argtypes = [c_char_p, c_char_p, c_char_p, c_char_p, c_char_p, c_char_p, c_char_p, RtFloat]
#    ri.RiMakeLatLongEnvironment.argtypes = [c_char_p, c_char_p, RtFilterFunc, RtFloat, RtFloat]
    ri.RiMakeLatLongEnvironment.argtypes = [c_char_p, c_char_p]
    ri.RiMakeShadow.argtypes = [c_char_p, c_char_p]
#    ri.RiMakeTexture.argtypes = [c_char_p, c_char_p, RtToken, RtToken, RtFilterFunc, RtFloat, RtFloat]
    ri.RiMakeTexture.argtypes = [c_char_p, c_char_p, RtToken, RtToken]
    ri.RiMatte.argtypes = [RtBoolean]
    ri.RiMotionBegin.argtypes = [RtInt]
    ri.RiMotionEnd.argtypes = []
    ri.RiNuPatch.argtypes = [RtInt, RtInt, POINTER(RtFloat), RtFloat, RtFloat, RtInt, RtInt, POINTER(RtFloat), RtFloat, RtFloat]
    ri.RiObjectBegin.argtypes = []
    ri.RiObjectBegin.restype = RtObjectHandle
    ri.RiObjectEnd.argtypes = []
    ri.RiObjectInstance.argtypes = [RtObjectHandle]
    ri.RiOpacity.argtypes = [RtColor]
    ri.RiOption.argtypes = [RtToken]
    ri.RiOrientation.argtypes = [RtToken]
    ri.RiParaboloid.argtypes = [RtFloat, RtFloat, RtFloat, RtFloat]
    ri.RiPatch.argtypes = [RtToken]
    ri.RiPatchMesh.argtypes = [RtToken, RtInt, RtToken, RtInt, RtToken]
    ri.RiPerspective.argtypes = [RtFloat]
    # The following line is commented out because declaring the first argument
    # as RtFilterFunc will prevent the standard filters from being passed in directly
    # so that the renderer recognizes them as standard filters (which is necessary
    # during RIB generation so that the renderer can generate the appropriate
    # filter name).
#    ri.RiPixelFilter.argtypes = [RtFilterFunc, RtFloat, RtFloat]
    ri.RiPixelSamples.argtypes = [RtFloat, RtFloat]
    ri.RiPixelVariance.argtypes = [RtFloat]
    ri.RiPoints.argtypes = [RtInt]
    ri.RiPointsGeneralPolygons.argtypes = [RtInt, POINTER(RtInt), POINTER(RtInt), POINTER(RtInt)]
    ri.RiPointsPolygons.argtypes = [RtInt, POINTER(RtInt), POINTER(RtInt)]
    ri.RiPolygon.argtypes = [RtInt]
#    ri.RiProcedural.argtypes = [RtPointer, RtBound, RtProcSubdivFunc, RtProcFreeFunc]
    ri.RiProcedural.argtypes = [RtPointer, RtBound]
    ri.RiProjection.argtypes = [RtToken]
    ri.RiQuantize.argtypes = [RtToken, RtInt, RtInt, RtInt, RtFloat]
    # When the second argument is set to RtArchiveCallback then it is not possible to pass None
#    ri.RiReadArchive.argtypes = [RtToken, RtArchiveCallback]
    ri.RiReadArchive.argtypes = [RtToken]
    ri.RiRelativeDetail.argtypes = [RtFloat]
    ri.RiReverseOrientation.argtypes = []
    ri.RiRotate.argtypes = [RtFloat, RtFloat, RtFloat, RtFloat]
    ri.RiScale.argtypes = [RtFloat, RtFloat, RtFloat]
    ri.RiScreenWindow.argtypes = [RtFloat, RtFloat, RtFloat, RtFloat]
    ri.RiShadingInterpolation.argtypes = [RtToken]
    ri.RiShadingRate.argtypes = [RtFloat]
    ri.RiShutter.argtypes = [RtFloat, RtFloat]
    ri.RiSides.argtypes = [RtInt]
    ri.RiSkew.argtypes = [RtFloat, RtFloat, RtFloat, RtFloat, RtFloat, RtFloat, RtFloat]
    ri.RiSolidBegin.argtypes = [RtToken]
    ri.RiSolidEnd.argtypes = []
    ri.RiSphere.argtypes = [RtFloat, RtFloat, RtFloat, RtFloat]
    ri.RiSubdivisionMesh.argtypes = [RtToken, RtInt, POINTER(RtInt), POINTER(RtInt), RtInt, POINTER(RtToken), POINTER(RtInt), POINTER(RtInt), POINTER(RtFloat)]
    ri.RiSurface.argtypes = [RtToken]
    ri.RiTextureCoordinates.argtypes = [RtFloat, RtFloat, RtFloat, RtFloat, RtFloat, RtFloat, RtFloat, RtFloat]
    ri.RiTorus.argtypes = [RtFloat, RtFloat, RtFloat, RtFloat, RtFloat]
    ri.RiTransform.argtypes = [RtMatrix]
    ri.RiTransformBegin.argtypes = []
    ri.RiTransformEnd.argtypes = []
    ri.RiTransformPoints.argtypes = [RtToken, RtToken, RtInt, POINTER(RtPoint)]
    ri.RiTranslate.argtypes = [RtFloat, RtFloat, RtFloat]
    ri.RiTrimCurve.argtypes = [RtInt, POINTER(RtInt), POINTER(RtInt), POINTER(RtFloat), POINTER(RtFloat), POINTER(RtFloat), POINTER(RtInt), POINTER(RtFloat), POINTER(RtFloat), POINTER(RtFloat)]
    ri.RiWorldBegin.argtypes = []
    ri.RiWorldEnd.argtypes = []
    
    ri.RiBoxFilter.argtypes = [c_float, c_float, c_float, c_float]
    ri.RiBoxFilter.restype = c_float
    ri.RiTriangleFilter.argtypes = [c_float, c_float, c_float, c_float]
    ri.RiTriangleFilter.restype = c_float
    ri.RiCatmullRomFilter.argtypes = [c_float, c_float, c_float, c_float]
    ri.RiCatmullRomFilter.restype = c_float
    ri.RiGaussianFilter.argtypes = [c_float, c_float, c_float, c_float]
    ri.RiGaussianFilter.restype = c_float
    ri.RiSincFilter.argtypes = [c_float, c_float, c_float, c_float]
    ri.RiSincFilter.restype = c_float
    
    # Optional calls from the 3.4 API:
    
    if hasattr(ri, "RiArchiveBegin"):
        ri.RiArchiveBegin.argtypes = [RtToken]
        ri.RiArchiveBegin.restype = c_void_p
    else:
        ri.RiArchiveBegin = _undefinedFunctionSubstitute

    if hasattr(ri, "RiArchiveEnd"):
        ri.RiArchiveEnd.argtypes = []
    else:
        ri.RiArchiveEnd = _undefinedFunctionSubstitute

    if hasattr(ri, "RiCamera"):
        ri.RiCamera.argtypes = [RtToken]
    else:
        ri.RiCamera = _undefinedFunctionSubstitute

    if hasattr(ri, "RiDisplayChannel"):
        ri.RiDisplayChannel.argtypes = [RtToken]
    else:
        ri.RiDisplayChannel = _undefinedFunctionSubstitute

    if hasattr(ri, "RiElse"):
        ri.RiElse.argtypes = []
    else:
        ri.RiElse = _undefinedFunctionSubstitute

    if hasattr(ri, "RiElseIf"):
        ri.RiElseIf.argtypes = [RtToken]
    else:
        ri.RiElseIf = _undefinedFunctionSubstitute
        
    if hasattr(ri, "RiIfBegin"):
        ri.RiIfBegin.argtypes = [RtToken]
    else:
        ri.RiIfBegin = _undefinedFunctionSubstitute

    if hasattr(ri, "RiIfEnd"):
        ri.RiIfEnd.argtypes = []
    else:
        ri.RiIfEnd = _undefinedFunctionSubstitute

    if hasattr(ri, "RiMakeBrickMap"):
        ri.RiMakeBrickMap.argtypes = [RtInt, POINTER(c_char_p), c_char_p]
    else:
        ri.RiMakeBrickMap = _undefinedFunctionSubstitute

    if hasattr(ri, "RiResource"):
        ri.RiResource.argtypes = [RtToken, RtToken]
    else:
        ri.RiResource = _undefinedFunctionSubstitute

    if hasattr(ri, "RiResourceBegin"):
        ri.RiResourceBegin.argtypes = []
    else:
        ri.RiResourceBegin = _undefinedFunctionSubstitute

    if hasattr(ri, "RiResourceEnd"):
        ri.RiResourceEnd.argtypes = []
    else:
        ri.RiResourceEnd = _undefinedFunctionSubstitute

    if hasattr(ri, "RiScopedCoordinateSystem"):
        ri.RiScopedCoordinateSystem.argtypes = [RtToken]
    else:
        ri.RiScopedCoordinateSystem = _undefinedFunctionSubstitute

    if hasattr(ri, "RiShader"):
        ri.RiShader.argtypes = [RtToken, RtToken]
    else:
        ri.RiShader = _undefinedFunctionSubstitute

    if hasattr(ri, "RiSystem"):
        ri.RiSystem.argtypes = [c_char_p]
    else:
        ri.RiSystem = _undefinedFunctionSubstitute
    
def _undefinedFunctionSubstitute(*args, **kwargs):
    raise NotImplementedError("This function call is not implemented by this implementation of the RenderMan interface")
    
def _getLastError(ri):
    """Getter function for the RiLastError variable.
    """
    try:
        return c_int.in_dll(ri, "RiLastError").value
    except ValueError:
        return None    

def _setLastError(ri, value):
    """Setter function for the RiLastError variable.
    """
    try:
        c_int.in_dll(ri, "RiLastError").value = value
    except ValueError:
        pass
    
   
