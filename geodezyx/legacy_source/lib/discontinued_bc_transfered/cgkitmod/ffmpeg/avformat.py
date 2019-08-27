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

import sys
from ctypes import *
from . import findlib
from . import avutil
from . import avcodec
from . import decls
from .decls import AVOutputFormat, AVFormatContext, AVStream


class AVFormatError(avutil.AVError):
    pass

######################################################################
# Functions
######################################################################

def avformat_version():
    """Return the libavformat library version.
    
    Returns a tuple (major,minor,micro) containing the three parts of the
    version number. 
    """
    v = _lib().avformat_version()
    major = v>>16
    minor = (v>>8)&0xff
    micro = v&0xff
    return (major,minor,micro)

def av_register_all():
    """Initialize libavformat and register all the (de)muxers and protocols.
    
    This must be called at the beginning before any file is opened.
    """
    _lib().av_register_all()

def av_oformat_next(fmt=None):
    """Return the first/next registered output format.
    
    If *fmt* is ``None``, the first registered output format is returned
    as a :class:`AVOutputFormat` object. If that object is passed back in
    as input, the next registered format is returned and so on. If there
    are no more formats available, ``None`` is returned.
    """
    func = _lib().av_oformat_next
    func.restype = POINTER(AVOutputFormat)
    if fmt is None:
        ptr = None
    else:
        ptr = byref(fmt)
    res = func(ptr)
    if res:
        res = res.contents
    else:
        res = None
    return res

def av_iformat_next(fmt=None):
    """Return the first/next registered input format.
    
    If *fmt* is ``None``, the first registered input format is returned
    as a :class:`AVInputFormat` object. If that object is passed back in
    as input, the next registered format is returned and so on. If there
    are no more formats available, ``None`` is returned.
    """
    func = _lib().av_iformat_next
    func.restype = POINTER(AVOutputFormat)
    if fmt is None:
        ptr = None
    else:
        ptr = byref(fmt)
    res = func(ptr)
    if res:
        res = res.contents
    else:
        res = None
    return res

def av_open_input_file(fileName, format=None, buf_size=0, params=None):
    """Open a media file as input. 

    *fileName* is a string containing the file to open.
    *format* is an optional :class:`AVInputFormat` object that can be specified
    to force a particular file format (AVInputFormat is not yet exposed).
    *buf_size* is an optional buffer size (or 0/``None`` if the default size is ok) 
    *params* is an optional :class:`AVFormatParameters` object if extra parameters
    have to be passed.
    
    The return value is a :class:`AVFormatContext` object that contains information
    about the format and that serves as a handle to the open file.
    In case of an error, an :exc:`AVFormatError` exception is thrown.
    """
    # Check input parameters
    if not isinstance(fileName, str):
        raise ValueError("fileName must be a string")
    if buf_size is None:
        buf_size = 0
    elif type(buf_size)!=int:
        raise ValueError("buf_size must be an int or None")
    if format is not None:
        raise ValueError("format parameter is not yet supported")
    if params is not None:
        raise ValueError("params parameter is not yet supported")
    
    formatCtxPtr = POINTER(AVFormatContext)()
    ret = _lib().av_open_input_file(byref(formatCtxPtr), fileName.encode("utf-8"), format, buf_size, params)
    if ret!=0:
        raise AVFormatError(ret)
    return formatCtxPtr.contents

def av_close_input_file(formatCtx):
    """Close a media file (but not its codecs).

    *formatCtx* is the format context as returned by :func:`av_open_input_file()`
    If it is ``None``, the function returns immediately.
    """
    if formatCtx is None:
        return
    _lib().av_close_input_file(byref(formatCtx))

def av_find_stream_info(formatCtx):
    """Read packets of a media file to get stream information.

    This is useful for file formats with no headers such as MPEG. This
    function also computes the real frame rate in case of mpeg2 repeat
    frame mode. The logical file position is not changed by this
    function; examined packets may be buffered for later processing.

    *formatCtx* is the media file handle as returned by :func:`av_open_input_file()`.
    The return value is the integer that was returned by the underlying
    C function (it is always >=0). In the case of an error, an :exc:`AVFormatError`
    exception is thrown.
    """
    ret = _lib().av_find_stream_info(byref(formatCtx))
    if ret<0:
        raise AVFormatError(ret)
    return ret

def av_read_frame(formatCtx, pkt):
    """Return the next frame of a stream. 

    *formatCtx* is the media file handle as returned by :func:`av_open_input_file()`.
    *pkt* is a :class:`AVPacket` object that will be filled with the packet data.

    The returned packet is valid until the next :func:`av_read_frame()` or until
    :func:`av_close_input_file()` and must be freed with :func:`av_free_packet()`.
    For video, the packet contains exactly one frame. For audio, it contains an
    integer number of frames if each frame has a known fixed size (e.g. PCM
    or ADPCM data). If the audio frames have a variable size (e.g. MPEG audio),
    then it contains one frame.

    ``pkt.pts``, ``pkt.dts`` and ``pkt.duration`` are always set to correct values
    in ``AVStream.timebase`` units (and guessed if the format cannot provided them).
    ``pkt.pts`` can be ``AV_NOPTS_VALUE`` if the video format has B frames, so it is
    better to rely on ``pkt.dts`` if you do not decompress the payload.
    
    Returns False when EOF was reached.
    """
    ret = _lib().av_read_frame(byref(formatCtx), byref(pkt))
    if ret<0:
        if ret==decls.AVERROR_EOF:
            return False
        if formatCtx.pb:
            if formatCtx.pb.contents.eof_reached:
                return False
        raise AVFormatError(ret)
    return True

def av_seek_frame(formatCtx, stream_index, timestamp, flags):
    """Seek to the key frame at timestamp. 

    'timestamp' in 'stream_index'. 

    stream_index     If stream_index is (-1), a default stream is selected, and timestamp is automatically converted from AV_TIME_BASE units to the stream specific time_base. 
    timestamp     timestamp in AVStream.time_base units or if there is no stream specified then in AV_TIME_BASE units 
    flags     flags which select direction and seeking mode (AVSEEK_FLAG_*)
    """
    func = _lib().av_seek_frame
    func.args = [POINTER(AVFormatContext), c_int, c_longlong, c_int]
    func.restype = c_int
    ret = func(byref(formatCtx), stream_index, timestamp, flags)
    if ret<0:
        raise AVFormatError(ret)
    return ret

def av_write_header(formatCtx):
    """Allocate stream private data and write the stream header.
    """
    ret = _lib().av_write_header(byref(formatCtx))
    if ret!=0:
        raise AVFormatError(ret)

def av_write_trailer(formatCtx):
    """Writes the stream trailer to an output media file and frees the file private data.
    
    May only be called after a successful call to av_write_header().
    """
    ret = _lib().av_write_trailer(byref(formatCtx))
    if ret!=0:
        raise AVFormatError(ret)
 
def av_write_frame(formatCtx, pkt):
    """Writes a packet to an output media file.

    The packet shall contain one audio or video frame.
    The packet must be correctly interleaved according to the container
    specification, if not then av_interleaved_write_frame() must be used.
    Returns 1 if end of stream is wanted, 0 otherwise.
    """
    ret = _lib().av_write_frame(byref(formatCtx), byref(pkt))
    if ret<0:
        raise AVFormatError(ret)
    return ret

def av_interleaved_write_frame(formatCtx, pkt):
    """Writes a packet to an output media file ensuring correct interleaving.

    The packet must contain one audio or video frame.
    If the packets are already correctly interleaved, the application should
    call av_write_frame() instead as it is slightly faster. It is also important
    to keep in mind that completely non-interleaved input will need huge amounts
    of memory to interleave with this, so it is preferable to interleave at the
    demuxer level.
    Returns 1 if end of stream is wanted, 0 otherwise.
    """
    ret = _lib().av_interleaved_write_frame(byref(formatCtx), byref(pkt))
    if ret<0:
        raise AVFormatError(ret)
    return ret

def avformat_alloc_context():
    """Allocates and initializes an :class:`AVFormatContext` object.

    Returns an initialized :class:`AVFormatContext` object that can be used to
    write a new file.
    
    The allocated structure must be freed with :func:`avutil.av_free()`
    (everything that has been explicitly allocated must also be freed explicitly).
    """
    func = _lib().avformat_alloc_context
    func.restype = POINTER(AVFormatContext)
    ctx = func()
    if ctx:
        return ctx.contents
    else:
        raise MemoryError("Failed to allocate AVFormatContext struct")

def av_guess_format(shortName=None, fileName=None, mimeType=None):
    """Search for an output format.
    
    Returns the output format in the list of registered output formats
    which best matches the provided parameters, or returns NULL if
    there is no match.

    shortName (if not None) checks if shortName matches with the
    names of the registered formats.
    fileName (if not None) checks if file name terminates with the
    extensions of the registered formats.
    mimeType (if not None) checks if mime type matches with the
    MIME type of the registered formats.
    
    Returns an :class:`AVOutputFormat` object or ``None``. The returned object
    is owned by ffmpeg and must not be freed.
    """
    # in v52 the function is still called guess_format instead of av_guess_format
    func = _lib().guess_format
    func.restype = POINTER(AVOutputFormat)
    fmt = func(shortName, fileName, mimeType)
    if fmt:
        return fmt.contents
    else:
        return None

def av_guess_codec(outputFormat, codecType, shortName=None, fileName=None, mimeType=None):
    """Guesses the codec ID based upon muxer and filename.
    
    *outputFormat* is an :class:`AVOutputFormat` object.
    *codecType* determines what kind of codec is required, it can be either
    ``CODEC_TYPE_VIDEO`` or ``CODEC_TYPE_AUDIO``.
    Returns a codec id (int) or ``None`` if no codec could be found.
    """
    id = _lib().av_guess_codec(byref(outputFormat), shortName, fileName, mimeType, codecType)
    if id==decls.CODEC_ID_NONE:
        return None
    else:
        return id

def av_new_stream(formatCtx, id):
    """Adds a new stream to a media file.
    
    Allocates a new stream and adds it to formatCtx.
    *formatCtx* is an :class:`AVFormatContext` object which will receive the
    new stream. *id* is an integer containing a file-format-dependent stream id.
    Returns an :class:`AVStream` object (which will also be available in
    ``formatCtx.streams``).
    
    If the stream couldn't be allocated, an :exc:`AVFormatError` exception is
    thrown.
    """
    func = _lib().av_new_stream
    func.restype = POINTER(AVStream)
    stream = func(byref(formatCtx), id)
    if stream:
        return stream.contents
    else:
        raise AVFormatError("Failed to allocate new media stream.")

def av_init_packet(pkt):
    """Initialize optional fields of a packet with default values.
    
    *pkt* must be a :class:`AVPacket` object.
    """
    _lib().av_init_packet(byref(pkt));


def url_fopen(s, fileName, flags):
    """

    avio.h: int url_fopen(ByteIOContext **s, const char *filename, int flags);
    
    flags:
    #define URL_RDONLY 0
    #define URL_WRONLY 1
    #define URL_RDWR   2
    """
    ret = _lib().url_fopen(byref(s), fileName, flags)
    if ret<0:
        raise AVFormatError('Error opening file "%s" for writing (error %s)'%(fileName, ret))
    return ret

def url_fclose(s):
    """
    avio.h: int url_fclose(ByteIOContext *s);
    """
    ret = _lib().url_fclose(s)
    return ret

def dump_format(formatCtx, index, url, is_output):
    _lib().dump_format(byref(formatCtx), index, url, is_output) 

#####################################################

# By default, try to load the avformat library that has the same major version
# than the one that was used for creating the cppdefs and decls module.
_libname = "avformat.%s"%(decls.LIBAVFORMAT_VERSION_MAJOR)
_libavformat = None

def _lib():
    """Return the avformat shared library.
    """
    global _libavformat
    global _libname
    if _libavformat is None:
        _libavformat = findlib.findFfmpegLib(_libname)
    return _libavformat

