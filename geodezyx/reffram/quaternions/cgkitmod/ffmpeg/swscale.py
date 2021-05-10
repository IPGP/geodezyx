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
# Portions created by the Initial Developer are Copyright (C) 2010
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

from ctypes import *
from . import findlib
from . import decls
from .decls import SwsContext

def swscale_version():
    """Return the libswscale library version.
    
    Returns a tuple (major,minor,micro) containing the three parts of the
    version number. 
    """
    v = _lib().swscale_version()
    major = v>>16
    minor = (v>>8)&0xff
    micro = v&0xff
    return (major,minor,micro)

def sws_getContext(srcW, srcH, srcFormat, dstW, dstH, dstFormat, flags):
    """Allocate a SwsContext object.
    
    Must be deallocated using sws_freeContext().
    Returns a :class:`SwsContext` object.
    """
    func = _lib().sws_getContext
    func.restype = POINTER(SwsContext)
    res = func(srcW, srcH, srcFormat, dstW, dstH, dstFormat, flags, 0, 0, 0)
    if res:
        return res.contents
    else:
        return None

def sws_freeContext(swsCtx):
    """Free a SwsContext object.
    """
    _lib().sws_freeContext(byref(swsCtx))   
    
def sws_scale(swsCtx, src, srcStride, srcSliceY, srcSliceH, dst, dstStride):
    """Scale/convert an image.
    
    Scales the image slice in srcSlice and puts the resulting scaled
    slice in the image in dst. A slice is a sequence of consecutive
    rows in an image.

    Slices have to be provided in sequential order, either in
    top-bottom or bottom-top order. If slices are provided in
    non-sequential order the behavior of the function is undefined.

    *swsCtx* is the scaling context as returned by sws_getContext().
    *src* is the array containing the pointers to the planes of the source slice.
    *srcStride* is the array containing the strides for each plane of the
    source image
    *srcSliceY* is the position in the source image of the slice to process,
    that is the number (counted starting from zero) in the image of the first
    row of the slice.
    *srcSliceH* is the height of the source slice, that is the number
    of rows in the slice.
    *dst* is the array containing the pointers to the planes of the destination
    image.
    *dstStride* is the array containing the strides for each plane of
    the destination image
    Returns the height of the output slice.
    """
    ret = _lib().sws_scale(byref(swsCtx), src, srcStride, srcSliceY, srcSliceH, dst, dstStride)
    if ret<0:
        raise RuntimeError("Error: %s"%ret)
    return ret


# By default, try to load the swscale library that has the same major version
# than the one that was used for creating the cppdefs and decls module.
_libname = "swscale.%s"%(decls.LIBSWSCALE_VERSION_MAJOR)
_libswscale = None

def _lib():
    """Return the swscale shared library.
    """
    global _libswscale
    global _libname
    if _libswscale is None:
        _libswscale = findlib.findFfmpegLib(_libname)
    return _libswscale

    