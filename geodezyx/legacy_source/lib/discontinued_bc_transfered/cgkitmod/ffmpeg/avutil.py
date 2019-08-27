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

import ctypes
from . import findlib
from . import decls
from .decls import AVRational

class AVError(Exception):
    """Base AV error class.
    """
    def __init__(self, msgOrErrnum):
        """Constructor.
        
        msgOrErrnum is either a string containing an error message or an integer
        containing the ffmpeg error number. If a number is passed, the number
        is converted into an error message.
        The error number is stored in the attribute errnum. If a string has
        been passed, that attribute will be None.
        """
        self.errnum = None
        # If the input argument is an int, then try to convert it into a message
        if type(msgOrErrnum) is int:
            self.errnum = msgOrErrnum
            msg = av_strerror(self.errnum)
            if msg is None:
                msg = "Error %s"%(msgOrErrnum)
        else:
            msg = msgOrErrnum
        
        Exception.__init__(self, msg)


def avutil_version():
    """Return the libavutil library version.
    
    Returns a tuple (major,minor,micro) containing the three parts of the
    version number. 
    """
    v = _lib().avutil_version()
    major = v>>16
    minor = (v>>8)&0xff
    micro = v&0xff
    return (major,minor,micro)

def av_free(obj):
    """Free memory which has been allocated with av_malloc(z)() or av_realloc().
    
    *obj* may be a AVFrame object that was allocated using avcodec_alloc_frame()
    (note: *obj* must be the object itself, not a pointer to it).
    """
    if obj is None:
        return
    _lib().av_free(ctypes.byref(obj))

def av_d2q(d, max):
    """Convert a floating point number into a rational.
    
    *d* is the floating point number and *max* the maximum allowed numerator
    and denominator
    Returns a :class:`AVRational` object representing *d*.
    """
    func = _lib().av_d2q
    func.args = [ctypes.c_double, ctypes.c_int]
    func.restype = AVRational
    return func(ctypes.c_double(d), int(max))

def av_strerror(errnum):
    """Return a description of an AVERROR code.
    
    Returns ``None`` if no description could be found. 
    """
    buf = ctypes.create_string_buffer(200)
    func = _lib().av_strerror
    res = func(int(errnum), buf, ctypes.sizeof(buf))
    if res==0:
        return buf.value
    else:
        return None

def av_get_sample_fmt_name(sampleFmt):
    """Return the name of the sample format or ``None`` if an unknown value is passed.
    
    *sampleFmt* is an ``AV_SAMPLE_FMT_*`` value.
    """
    func = _lib().av_get_sample_fmt_name
    func.restype = ctypes.c_char_p
    name = func(sampleFmt)
    return name

def av_get_sample_fmt(name):
    """Return the sample format corresponding to the given name.
    
    Returns ``AV_SAMPLE_FMT_NONE`` if the name is not recognized.
    """
    return _lib().av_get_sample_fmt(name)

def av_get_bytes_per_sample(sampleFmt):
    """Return the number of bytes per sample.
    
    *sampleFmt* is an ``AV_SAMPLE_FMT_*`` value. Returns 0 if sampleFmt refers to
    an unknown sample format.
    """
    return _lib().av_get_bytes_per_sample(sampleFmt)

def av_get_pix_fmt_name(pixFmt):
    """Return the name of the pixel format or ``None`` if an unknown value is passed.
    
    *pixFmt* is a ``PIX_FMT_*`` value.
    """
    func = _lib().av_get_pix_fmt_name
    func.restype = ctypes.c_char_p
    name = func(pixFmt)
    return name
    
def av_get_pix_fmt(name):
    """Return the pixel format corresponding to the given name.
    
    Returns ``PIX_FMT_NONE`` if the name is not recognized.
    """
    return _lib().av_get_pix_fmt(name)


# By default, try to load the avutil library that has the same major version
# than the one that was used for creating the cppdefs and decls module.
_libname = "avutil.%s"%(decls.LIBAVUTIL_VERSION_MAJOR)
_libavutil = None

def _lib():
    """Return the avutil shared library.
    """
    global _libavutil
    global _libname
    
    if _libavutil is None:
        _libavutil = findlib.findFfmpegLib(_libname)
    return _libavutil

