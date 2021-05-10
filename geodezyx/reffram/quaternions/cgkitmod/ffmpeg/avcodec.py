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
from ctypes import *
from . import findlib
from . import avutil
from . import decls
from .decls import AVCodec, AVFrame


class AVCodecError(avutil.AVError):
    pass


def avcodec_version():
    """Return the libavcodec library version.
    
    Returns a tuple (major,minor,micro) containing the three parts of the
    version number. 
    """
    v = _lib().avcodec_version()
    major = v>>16
    minor = (v>>8)&0xff
    micro = v&0xff
    return (major,minor,micro)

def avcodec_find_decoder(id):
    """Finds a decoder with a matching codec ID.
    
    id is an integer containing the codec id.
    Returns a AVCodec object or None if no decoder was found.
    """
    # The returned pointer is owned by the library, so we don't have to
    # do anything to get rid of it.
    func = _lib().avcodec_find_decoder
    func.restype = POINTER(AVCodec)
    res = func(id)
    if res:
        return res.contents
    else:
        return None

def avcodec_find_decoder_by_name(name):
    """Finds a registered decoder with the specified name.
    
    Returns a AVCodec object or None if no decoder was found.
    """
    # The returned pointer is owned by the library, so we don't have to
    # do anything to get rid of it.
    func = _lib().avcodec_find_decoder_by_name
    func.args = [ctypes.c_char_p]
    func.restype = POINTER(AVCodec)
    res = func(name)
    if res:
        return res.contents
    else:
        return None

def avcodec_find_encoder(id):
    """Finds a registered encoder with a matching codec ID.
    
    id is an integer containing the codec id.
    Returns an AVCodec object or None if no encoder was found.
    """
    # The returned pointer is owned by the library, so we don't have to
    # do anything to get rid of it.
    func = _lib().avcodec_find_encoder
    func.restype = POINTER(AVCodec)
    res = func(id)
    if res:
        return res.contents
    else:
        return None

def avcodec_find_encoder_by_name(name):
    """Finds a registered encoder with the specified name.
    
    Returns an AVCodec object or None if no encoder was found.
    """
    # The returned pointer is owned by the library, so we don't have to
    # do anything to get rid of it.
    func = _lib().avcodec_find_encoder_by_name
    func.args = [ctypes.c_char_p]
    func.restype = POINTER(AVCodec)
    res = func(name)
    if res:
        return res.contents
    else:
        return None

def avcodec_open(avctx, codec):
    """Initializes the AVCodecContext to use the given AVCodec. 

    avctx is a AVCodecContext object and codec a AVCodec object.
    Raises an exception when an error occurs.

    Prior to using this function the context has to be allocated.

    The functions avcodec_find_decoder_by_name(), avcodec_find_encoder_by_name(),
    avcodec_find_decoder() and avcodec_find_encoder() provide an easy way for 
    retrieving a codec.
    
    Warning: This function is not thread safe!
    """
    ret = _lib().avcodec_open(byref(avctx), byref(codec))
    if ret<0:
        # It can happen that the file was actually opened, so try to close it again.
        # The line is commented out because closing it here seems to cause more harm
        # (sometimes it leads to the lib not being able to open anything at all anymore)
        #avcodec_close(avctx)
        raise AVCodecError(ret)
    return ret    
    
def avcodec_close(codecCtx):
    """Close the codec.
    """
    return _lib().avcodec_close(byref(codecCtx))
    
def avcodec_alloc_frame():
    """Allocates an AVFrame and sets its fields to default values. 

    The resulting struct can be deallocated by simply calling avutil.av_free().
    Returns an AVFrame object.
    """
    func = _lib().avcodec_alloc_frame
    func.restype = POINTER(AVFrame)
    res = func()
    if res:
        return res.contents
    else:
        return None

def av_init_packet(pkt):
    """Initialize optional fields of a packet with default values.
    
    *pkt* is a :class:`AVPacket` object.
    Don't call this function on a packet that already contains data as this
    would lead to a memory leak.
    """
    _lib().av_init_packet(byref(pkt))
        
def av_free_packet(pkt):
    """Free a packet.
    
    *pkt* is a :class:`AVPacket` object.
    Frees the data associated with the packet (if there is any data set).
    """
    _lib().av_free_packet(byref(pkt))

#def avcodec_decode_video(codecCtx, picture, buf, bufsize):
#    """Decodes a video frame from buf into picture. 
#
#    The avcodec_decode_video() function decodes a video frame from the input
#    buffer buf of size buf_size. To decode it, it makes use of the video
#    codec which was coupled with avctx using avcodec_open(). The resulting 
#    decoded frame is stored in picture.    
#    Returns a tuple (got_picture, bytes) where got_picture is a boolean that
#    indicates whether a frame was decoded (True) or not and bytes is the number
#    of bytes used.
#    """
#    got_picture = c_int()
#    ret = _lib().avcodec_decode_video(byref(codecCtx), byref(picture), byref(got_picture), buf, bufsize)
#    if ret<0:
#        raise AVCodecError("Error: %s"%ret)
#    return bool(got_picture.value), ret

def avcodec_decode_video2(codecCtx, picture, pkt):
    """Decodes a video frame from pkt into picture. 

    Some decoders may support multiple frames in a single AVPacket, such
    decoders would then just decode the first frame.
    
    Some codecs have a delay between input and output, these need to be fed
    with avpkt.data=NULL, avpkt.size=0 at the end to return the remaining frames.

    *codecCtx* is an :class:`AVCodecCtx` object, *picture* is an :class:`AVFrame`
    object which will receive the decoded video frame. *pkt* is a :class:`AVPacket`
    object containing the encoded data.
    
    Returns a tuple (*got_picture*, *bytes*) where *got_picture* is a boolean that
    indicates whether a frame was decoded (True) or not and *bytes* is the number
    of bytes used or 0 if no frame could be decompressed.
    """
    got_picture = c_int()
    ret = _lib().avcodec_decode_video2(byref(codecCtx), byref(picture), byref(got_picture), byref(pkt))
    if ret<0:
        raise AVCodecError(ret)
    return bool(got_picture.value), ret

#def avcodec_decode_audio2(codecCtx, sampleBuf, buf, bufsize):
#    """Decode an audio frame.
#    
#    sampleBuf is a ctypes short array.
#    buf is the encoded data and bufsize the size of the encoded data.
#    Returns the frameSize (size in bytes of the decoded frame) and the number
#    of bytes that have been used from buf.
#    """
#    frameSize = c_int(ctypes.sizeof(sampleBuf))
#    ret = _lib().avcodec_decode_audio2(byref(codecCtx), sampleBuf, byref(frameSize), buf, bufsize)
#    if ret<0:
#        raise AVCodecError(ret)
#    return frameSize.value, ret

def avcodec_decode_audio3(codecCtx, sampleBuf, pkt):
    """Decode an audio frame.
    
    Some decoders may support multiple frames in a single AVPacket, such
    decoders would then just decode the first frame. In this case,
    avcodec_decode_audio3 has to be called again with an AVPacket that contains
    the remaining data in order to decode the second frame etc.
    
    *codecCtx* is an :class:`AVCodecCtx` object, *sampleBuf* is a ctypes short
    array that will receive the decoded audio data.
    *pkt* is a :class:`AVPacket` object containing the encoded data.
    
    Returns the frameSize (size in bytes of the decoded frame) and the number
    of bytes that have been used from the input buffer in *pkt*.
    If there are no more frames, both values are 0.
    """
    # The size of the output buffer in bytes
    frameSize = c_int(ctypes.sizeof(sampleBuf))
    ret = _lib().avcodec_decode_audio3(byref(codecCtx), sampleBuf, byref(frameSize), byref(pkt))
    if ret<0:
        raise AVCodecError(ret)
    return frameSize.value, ret

def avcodec_encode_video(codecCtx, buf, bufsize, picture):
    """Encodes a video frame from picture into buf.

    *codecCtx* is a :class:`AVCodecContext` object, *buf* is a pointer to
    the output buffer, *bufsize* is the size in bytes of the buffer and
    *picture* is a :class:`AVFrame` object containing the picture to encode.
    
    The input picture should be stored using a specific format, namely
    avctx.pix_fmt.
    Returns 0 or the number of bytes used from the output buffer.
    """
    ret = _lib().avcodec_encode_video(byref(codecCtx), buf, bufsize, byref(picture))
    if ret<0:
        raise AVCodecError(ret)
    return ret

def avpicture_alloc(picture, pix_fmt, width, height):
    """Allocate memory for a picture. 

    Call avpicture_free() to free it.
    
    picture     the AVPicture obejct to be filled in 
    pix_fmt     the format of the picture 
    width     the width of the picture 
    height     the height of the picture
    """
    ret = _lib().avpicture_alloc(byref(picture), pix_fmt, width, height)
    if ret<0:
        raise RuntimeError("Error: %s"%ret)
    
def avpicture_free(picture):
    """Free a picture previously allocated by avpicture_alloc().
    """
    _lib().avpicture_free(byref(picture))
    
def avpicture_get_size(pix_fmt, width, height):
    """Calculate image size.
    
    Calculate the size in bytes that a picture of the given 
    width and height would occupy if stored in the given picture format.
    
    pix_fmt     the given picture format 
    width     the width of the image 
    height     the height of the image
    
    Returns the size in bytes.
    """
    ret = _lib().avpicture_get_size(pix_fmt, width, height)
    if ret<0:
        raise AVCodecError(ret)
    return ret

def avpicture_fill(picture, ptr, pix_fmt, width, height):
    """Fill in the AVPicture fields. 

    The fields of the given AVPicture are filled in by using the 'ptr'
    address which points to the image data buffer. Depending on the
    specified picture format, one or multiple image data pointers and
    line sizes will be set. If a planar format is specified, several
    pointers will be set pointing to the different picture planes and
    the line sizes of the different planes will be stored in the
    lines_sizes array.
    
    picture     AVPicture whose fields are to be filled in 
    ptr     Buffer which will contain or contains the actual image data 
    pix_fmt     The format in which the picture data is stored. 
    width     the width of the image in pixels 
    height     the height of the image in pixels  

    Returns the size of the image data in bytes
    """
    return _lib().avpicture_fill(byref(picture), ptr, pix_fmt, width, height)

def avcodec_flush_buffers(codecCtx):
    """Flush buffers.

    Should be called when seeking or when switching to a different stream.
    """
    _lib().avcodec_flush_buffers(byref(codecCtx))

def avcodec_get_context_defaults2(codecCtx, codecType):
    """Set context defaults.
    
    In avcodec.h this function is marked as not being public, but what
    is the official way to do the same thing? (we are using it anyway
    as it's also used by the ffmpeg command line tool which we use as
    a reference)
    
    codecType is CODEC_TYPE_AUDIO, CODEC_TYPE_VIDEO, etc.
    """
    _lib().avcodec_get_context_defaults2(byref(codecCtx), codecType);

######################################################################

# By default, try to load the avcodec library that has the same major version
# than the one that was used for creating the cppdefs and decls module.
_libname = "avcodec.%s"%(decls.LIBAVCODEC_VERSION_MAJOR)
_libavcodec = None

def _lib():
    """Return the avcodec shared library.
    """
    global _libavcodec
    global _libname
    if _libavcodec is None:
        _libavcodec = findlib.findFfmpegLib(_libname)
    return _libavcodec

