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
import fractions
from .ffmpeg import decls, cppdefs, avutil, swscale, avformat, avcodec
import ctypes

try:
    import Image
    _pilImportException = None
except ImportError as exc:
    _pilImportException = exc
    
try:
    import numpy
    _numpyImportException = None
except ImportError as exc:
    _numpyImportException = exc

# Pixel formats
RGB = 0
BGR = 1
RGBA = 2
ARGB = 3
BGRA = 4
ABGR = 5
GRAY = 6
 
# Pixel access
WIDTH_HEIGHT = 100
HEIGHT_WIDTH = 101

# Color access
SEPARATE_CHANNELS = 200
COMBINED_CHANNELS = 201

class MediaFileError(Exception):
    pass

class VideoData:
    """Decoded image data.
    
    This class serves as a container for the data that is passed to the
    video callbacks.
    """
    def __init__(self, stream, size, srcPixFmt, frame=None, pts=None):
        """Constructor.
        """
        # The parent stream object (VideoStream_ffmpeg)
        self.stream = stream
        # Presentation timestamp (float)
        self.pts = pts
        # Image size (width, height)
        self.size = size
        
        # The source pixel format of the video stream.
        # This is one of the ffmpeg PIX_FMT_* integer values and it's used to
        # convert the video frame into the target format.
        self._srcPixFmt = srcPixFmt
        
        # AVFrame object (containing the decoded frame (but not as RGB))
        # The contents of this object change on every frame iteration, so the
        # same object is used for every frame.
        self._frame = frame
        
        # If the frame is queried as a numpy array, this will be the array that
        # is kept between calls.
        self._numpyArray = None
        # Stores a tuple (pixelFormat, pixelAccess, colorAccess) that contains
        # the values that were used to initialize the numpy array. Whenever
        # one of these values changes, the array needs to be re-initialized.
        self._numpyBufferParams = None
        # The sws context that is used to convert the image into a RGB numpy array
        self._numpySwsCtx = None
        # The uint*[4] array containing the data pointers into the numpy array
        # (actually only, the first item is used, the others are NULL). This
        # is passed into the sws_scale() function to fill the numpy buffer.
        self._dataPtrs = None
        # The int[4] array of line sizes (only the first one is used). This
        # is passed into the sws_scale() function to fill the numpy buffer.
        self._lineSizes = None
        
        # AVPicture object that contains the target buffer for conversion to PIL.
        # This is just a buffer that can hold the converted image (just like
        # a numpy array).
        self._picture = None
        # The size in bytes of a full (RGB) picture
        self._pictureSize = None
        # The sws context for conversion into a PIL image
        self._pilSwsCtx = None
    
    def _free(self):
        """Free any resources that have been allocated by the object.
        
        This must be called by the stream after all frames have been processed.
        """
        if self._numpySwsCtx is not None:
            swscale.sws_freeContext(self._numpySwsCtx)
            self._numpySwsCtx = None
            
        if self._pilSwsCtx is not None:
            swscale.sws_freeContext(self._pilSwsCtx)
            self._pilSwsCtx = None
        
        if self._picture is not None:
            avcodec.avpicture_free(self._picture)
            self._picture = None

        
    def isKeyFrame(self):
        """Check whether the current frame is a key frame or not.
        """
        return self._frame.key_frame==1
    
    def numpyArray(self, pixelFormat=RGB, pixelAccess=WIDTH_HEIGHT, colorAccess=SEPARATE_CHANNELS):
        """Return the current frame as a numpy array.

        pixelFormat specifies the number and order of the channels in the
        returned buffer. The source video frames are converted into the
        specified format.
        
        pixelAccess defines the order of the x,y indices when accessing
        a pixel in the returned buffer. The value can either be WIDTH_HEIGHT
        to specify the pixel position as (x,y) or HEIGHT_WIDTH if (y,x)
        is required.
        
        colorAccess specifies whether the color channels should be accessed
        as a third index (SEPARATE_CHANNELS) or if the entire pixel should
        be stored as a single integer value (COMBINED_CHANNELS). The latter
        is only allowed if the pixel format contains 1 or 4 channels.
        
        Note that the pixelAccess and colorAccess parameters do not affect
        the underlying memory layout of the buffer, they only affect the way
        that individual pixels are accessed via the numpy array.
        The memory layout is always so that rows are stored one after another
        without any gaps and the channel values of an individual pixel are
        always stored next to each other.  
        
        The underlying numpy array is reused for every frame, so if you need
        to keep the array you have to copy it.
        The returned array has a shape of (width,height,3) and a dtype of uint8.
        """
        # Make sure numpy has been imported successfully
        global _numpyImportException
        if _numpyImportException is not None:
            raise _numpyImportException
        
        # Do we need to allocate a new numpy buffer?
        if self._numpyArray is None or (pixelFormat, pixelAccess, colorAccess)!=self._numpyBufferParams:
            self._initNumpyProcessing(pixelFormat, pixelAccess, colorAccess)

        # Convert the frame into an rgb image. The result is put directly into
        # the numpy array.
        height = self.size[1]
        swscale.sws_scale(self._numpySwsCtx, self._frame.data, self._frame.linesize, 0, height, self._dataPtrs, self._lineSizes)
        
        return self._numpyArray
    
    def _initNumpyProcessing(self, pixelFormat, pixelAccess, colorAccess):
        """Initialize decoding into a numpy array.
        
        This has to be called at least once before a frame is converted into
        a numpy array.
        If the layout of the numpy array changes, the method has to be called
        again.
        
        The method creates the numpy array, the swscale context to convert
        the image data and the _lineSizes and _dataPtrs arrays for the sws_scale()
        call.
        """
        # Check the pixelAccess value...
        if pixelAccess not in [WIDTH_HEIGHT, HEIGHT_WIDTH]:
            raise ValueError("Invalid pixelAcess value")
        
        # Check the colorAccess value...
        if colorAccess not in [SEPARATE_CHANNELS, COMBINED_CHANNELS]:
            raise ValueError("Invalid colorAccess value")
        
        # Check the pixel format...
        numChannels = None
        dstPixFmt = None
        if pixelFormat==RGB:
            numChannels = 3
            dstPixFmt = decls.PIX_FMT_RGB24
        elif pixelFormat==BGR:
            numChannels = 3
            dstPixFmt = decls.PIX_FMT_BGR24
        elif pixelFormat==RGBA:
            numChannels = 4
            dstPixFmt = decls.PIX_FMT_RGBA
        elif pixelFormat==ARGB:
            numChannels = 4
            dstPixFmt = decls.PIX_FMT_ARGB
        elif pixelFormat==ABGR:
            numChannels = 4
            dstPixFmt = decls.PIX_FMT_ABGR
        elif pixelFormat==BGRA:
            numChannels = 4
            dstPixFmt = decls.PIX_FMT_BGRA
        elif pixelFormat==GRAY:
            numChannels = 1
            dstPixFmt = decls.PIX_FMT_GRAY8
        else:
            raise ValueError("Invalid pixelFormat value")

        width,height = self.size
        # Allocate the numpy array that will hold the converted image.
        # The memory layout of the buffer is always so that the image is stored
        # in rows and all the channels are stored together with a pixel.
        # Example: Row 0: RGB-RGB-RGB-RGB...
        #          Row 1: RGB-RGB-RGB-RGB...
        #          ...
        # To get from one channel to the next channel, you always have to add 1 byte
        # (for int8 channels). To get to the next x position, you have to add
        # numChannels bytes and to get to the next y position you have to add
        # width*numChannels bytes.
        # The pixelAccess and colorAccess parameters don't affect this memory
        # layout, they only affect how a pixel is accessed via numpy.
        if colorAccess==SEPARATE_CHANNELS:
            # Allocate a RGB buffer...
            if pixelAccess==WIDTH_HEIGHT:
                self._numpyArray = numpy.empty((width,height,numChannels), dtype=numpy.uint8)
                # Adjust the strides so that image rows are consecutive
                self._numpyArray.strides = (numChannels, numChannels*width, 1)
            else:
                self._numpyArray = numpy.empty((height,width,numChannels), dtype=numpy.uint8)
        else:
            if numChannels not in [1,4]:
                raise ValueError("COMBINED_CHANNELS pixel access can only be used with 1-channel or 4-channel pixel formats")
            
            if numChannels==1:
                dtype = numpy.uint8
            else:
                dtype = numpy.uint32
            
            # Allocate a RGBA buffer (where a RGBA value is stored as a uint32)
            if pixelAccess==WIDTH_HEIGHT:
                self._numpyArray = numpy.empty((width,height), dtype=dtype, order="F")
            else:
                self._numpyArray = numpy.empty((height,width), dtype=dtype)

        # Free any previously allocated context 
        if self._numpySwsCtx is not None:
            swscale.sws_freeContext(self._numpySwsCtx)
            self._numpySwsCtx = None

        # Allocate the swscale context
        self._numpySwsCtx = swscale.sws_getContext(width, height, self._srcPixFmt, 
                                                   width, height, dstPixFmt, 1)

        # Initialize the lineSizes and dataPtrs array that are passed into sws_scale().
        # The first data pointer points to the beginning of the numpy buffer.
        DataPtrType = ctypes.POINTER(ctypes.c_uint8)
        self._lineSizes = (4*ctypes.c_int)(numChannels*width,0,0,0)
        self._dataPtrs = (4*DataPtrType)(self._numpyArray.ctypes.data_as(DataPtrType), None, None, None)
        
        # Store the parameters so that we can detect parameter changes
        self._numpyBufferParams = (pixelFormat, pixelAccess, colorAccess)
    
    def pilImage(self):
        """Return the current frame as a PIL image.
        """
        global _pilImportException
        if _pilImportException is not None:
            raise _pilImportException

        # Do we need to initialize the PIL conversion?
        if self._pilSwsCtx is None:
            self._initPilProcessing()

        # Convert the image into RGB format (the convered image is stored
        # in the AVPicture buffer)...
        height = self.size[1]
        swscale.sws_scale(self._pilSwsCtx, self._frame.data, self._frame.linesize, 0, height, self._picture.data, self._picture.linesize)
            
        # Obtain a Python string containing the RGB image data...
        dataStr = ctypes.string_at(self._picture.data[0], self._pictureSize)
        
        # Convert to PIL image
        return Image.fromstring("RGB", self.size, dataStr)

    def _initPilProcessing(self):
        """Initialize decoding into a PIL image.
        
        Initializes an intermediate buffer (AVPicture) that will receive
        the converted image and a swscale context to do the conversion. 
        """
        width,height = self.size
        
        dstPixFmt = decls.PIX_FMT_RGB24
        
        # In case, this is called several times, make sure we don't leak memory.
        self._free()
        
        # Allocate a buffer that can hold the converted image
        self._picture = decls.AVPicture()
        avcodec.avpicture_alloc(self._picture, dstPixFmt, width, height)
        # Get the size of the destination image buffer
        self._pictureSize = avcodec.avpicture_get_size(dstPixFmt, width, height)
        
        # Allocate the swscale context
        self._pilSwsCtx = swscale.sws_getContext(width, height, self._srcPixFmt, 
                                                 width, height, dstPixFmt, 1)


class AudioData:
    """Decoded audio data.
    """
    def __init__(self, stream, pts=None, channels=None, framerate=None, samples=None, sampleSize=None):
        # The parent stream
        self.stream = stream
        # Presentation timestamp as an integer in stream time base units
        self.pts = pts
        # The number of channels (int)
        self.channels = channels
        # The framerate in Hz
        self.framerate = framerate
        # The decoded samples (this is a ctypes array of shorts)
        self.samples = samples
        # The number of bytes in the samples buffer
        self.sampleSize = sampleSize


class StreamBase_ffmpeg(object):
    """Base class for the audio/video streams.
    """
    def __init__(self, stream):
        """Constructor.
        
        stream is an AVStream object.
        """
        object.__init__(self)

        # Data callback
        self._dataCallback = None
        
        # AVStream object
        self._stream = stream

        # Get the codec context of the stream (AVCodecContext)
        self._codecCtx = stream.codec.contents

        # Obtain an appropriate decoder (AVCodec)...
        self._codec = avcodec.avcodec_find_decoder(self._codecCtx.codec_id)

        if self._codec is not None:
            # Open the codec
            avcodec.avcodec_open(self._codecCtx, self._codec)

    def close(self):
        """Close the stream/codec.
        """
        if self._codecCtx is not None:
            avcodec.avcodec_close(self._codecCtx)
            self._codecCtx = None
            self._codec = None
    
    def setDataCallback(self, callback):
        """Set a callback that gets called with decoded data.
        
        *callback* is a callable that must take a data object as input.
        It may also be ``None`` to remove any previously set callback.
        """
        self._dataCallback = callback
    
    def getDataCallback(self):
        """Return the currently set data callback.
        """
        return self._dataCallback
    
    def decodeBegin(self):
        """This is called right before decoding begins.
        
        The method returns a boolean indicating whether the stream should
        be enabled (True) or not (False). A stream might be enabled if there
        is no callback (i.e. noone is interested in the decoded data).
        """
        pass
    
    def decodeEnd(self):
        """This is called after decoding has finished.
        
        This is only called when decodeBegin() didn't return False or raise
        an error.
        """
        pass
    
    def handlePacket(self, pkt):
        """This is called during decoding to handle an encoded file packet.
        
        This is only called when decodeBegin() didn't return False or raise
        an error.
        Yields the decoded data (in case of audio data, a packet may contain
        more than one frame).
        """
        pass

    @property
    def index(self):
        """Return the stream index.
        """
        return self._stream.index
        
    @property
    def fourCC(self):
        v = self._codecCtx.codec_tag
        return "%s%s%s%s"%(chr(v&0xff), chr((v>>8)&0xff), chr((v>>16)&0xff), chr((v>>24)&0xff))

    @property
    def codecName(self):
        if self._codec is None:
            return None
        else:
            return self._codec.name

    @property
    def codecLongName(self):
        if self._codec is None:
            return None
        else:
            return self._codec.long_name

    @property
    def bitRate(self):
        """Return the average bit rate as an integer (bits per second).
        """
        return self._codecCtx.bit_rate

    @property
    def timeBase(self):
        """Return the time base as a Fraction object.
        
        This is the fundamental unit of time (in seconds) in terms of which
        frame timestamps are represented
        """
        tb = self._stream.time_base
        num = tb.num
        den = tb.den
        # Make sure we don't divide by 0
        if den==0:
            num = 0
            den = 1
        return fractions.Fraction(num, den)
    
    @property
    def duration(self):
        """Return the duration in timeBase units.
        """
        return self._stream.duration
    


class AudioStream_ffmpeg(StreamBase_ffmpeg):
    """Audio stream class.
    """
    def __init__(self, stream):
        StreamBase_ffmpeg.__init__(self, stream)
    
    def decodeBegin(self):

        # A dummy packet that is used to pass the output buffer to avcodec_decode_audio3
        self._dummyPkt = decls.AVPacket()
        avcodec.av_init_packet(self._dummyPkt)

        # Allocate a buffer for the decoded audio frame
        size = 192000
        self._sampleBuf = (size*ctypes.c_short)()

        self._audioData = AudioData(self, channels=self._codecCtx.channels, framerate=self._codecCtx.sample_rate, samples=self._sampleBuf)
        return True

    def decodeEnd(self):
        pass
    
    def handlePacket(self, pkt):
        # Decode the frame...
        codecCtx = self._codecCtx
        
        # Set the data pointer and size of the dummy packet to that of packet
        # (the dummy packet is used because there may be more than one frame
        # and we have to adjust the data pointer and size to decode all of them)
        dummyPkt = self._dummyPkt
        dummyPkt.data = pkt.data
        dummyPkt.size = pkt.size
        
        bytesUsed = 1
        while bytesUsed>0 and bool(dummyPkt.data):
            # Decode the next frame
            frameSize,bytesUsed = avcodec.avcodec_decode_audio3(codecCtx, self._sampleBuf, dummyPkt)
            # frameSize is the number of bytes that have been written into the sample buffer
            # bytesUsed is the number of bytes that have been used in the input buffer

            # The following two lines are equivalent to dummyPkt.data += bytesUsed
            # (the first line converts the pointer value to an int, the second
            # line adds bytesUsed and converts the result back to a ctypes pointer)
            addr = ctypes.addressof(dummyPkt.data.contents)
            dummyPkt.data = ctypes.cast(addr+bytesUsed, type(dummyPkt.data))
            dummyPkt.size -= bytesUsed
            
            # Report the frame if we have one
            if frameSize>0:
                self._audioData.pts = pkt.pts
                self._audioData.sampleSize = frameSize
                
                yield self._audioData

    @property
    def numChannels(self):
        return self._codecCtx.channels

    @property
    def sampleRate(self):
        """Return the sample rate as an integer (samples per second).
        """
        return self._codecCtx.sample_rate


class VideoStream_ffmpeg(StreamBase_ffmpeg):
    """Video stream class.
    """

    def __init__(self, stream):
        """Constructor.
        
        stream is an AVStream object containing the video stream.
        """
        StreamBase_ffmpeg.__init__(self, stream)
        self._frame = None
        self._videoData = None
    
    def decodeBegin(self):
        
        # Allocate a video frame that will store the decoded frames and
        # create the VideoData object that will be passed to the callers
        # (the actual conversion to RGB (or whatever) is done by the
        # VideoData object)
        try:
            self._frame = avcodec.avcodec_alloc_frame()
            if self._frame is None:
                raise MemoryError("Failed to allocate AVFrame object")
            
            width = self._codecCtx.width
            height = self._codecCtx.height
            self._videoData = VideoData(self, size=(width, height), frame=self._frame, srcPixFmt=self._codecCtx.pix_fmt)
            
            tb = self._stream.time_base
            self._timebase = float(tb.num)/tb.den
        except:
            self.decodeEnd()
            raise

        return True

    def decodeEnd(self):
        if self._videoData is not None:
            self._videoData._free()
            self._videoData = None
        if self._frame is not None:
            avutil.av_free(self._frame)
            self._frame = None

    def handlePacket(self, pkt):
        # Decode the frame...
        codecCtx = self._codecCtx
        hasFrame,bytes = avcodec.avcodec_decode_video2(codecCtx, self._frame, pkt)

        if hasFrame:
            # Conversion into RGB or creation of a PIL image is done by the
            # video data object on demand
            
            #print "PTS:%s  DTS:%s  Duration:%s"%(pkt.pts, pkt.dts, pkt.duration)
            self._videoData.pts = pkt.dts*self._timebase
        
            yield self._videoData

    @property
    def size(self):
        """Return the width and height of a video frame in pixels.
        """
        return (self._codecCtx.width, self._codecCtx.height)
        
    @property
    def width(self):
        """Return the width of a video frame in pixels.
        """
        return self._codecCtx.width

    @property
    def height(self):
        """Return the height of a video frame in pixels.
        """
        return self._codecCtx.height
    
    @property
    def frameRate(self):
        """Return the frame rate as a Fraction object.
        """
        fr = self._stream.r_frame_rate
        num = fr.num
        den = fr.den
        # Make sure we don't divide by 0
        if den==0:
            num = 0
            den = 1
        return fractions.Fraction(num, den)
    
    @property
    def numFrames(self):
        """Return the number of frames if known.
        
        Returns an integer or None if the number of frames is not known.
        """
        n = self._stream.nb_frames
        if n==0:
            n = None
        return n
        
    @property
    def pixelAspect(self):
        """Return the pixel aspect ratio if the value is known.
        
        Returns a Fraction object or None if the aspect ratio is not known.
        """
        a = self._stream.sample_aspect_ratio
        num = a.num
        den = a.den
        if den!=0 and num!=0:
            return fractions.Fraction(num, den)
        else:
            return None


class Media_Read_ffmpeg(object):
    """Media file reader.
    """
    
    def __init__(self, fileName):
        """Constructor.
        
        fileName is the name of the media file.
        """
        object.__init__(self)
        
        # Audio streams (AudioStream_ffmpeg objects)
        self.audioStreams = []
        # Video streams (VideoStream_ffmpeg objects)
        self.videoStreams = []
    
        # All streams in file order
        self._streams = []
    
        # The AVFormatContext object for the open file
        self._formatCtx = None
        
        # Open the video file
        self._formatCtx = avformat.av_open_input_file(fileName, None, 0, None)
        # Fill the 'streams' fields...
        avformat.av_find_stream_info(self._formatCtx)

        # Create the stream wrapper objects
        for i in range(self._formatCtx.nb_streams):
            stream = self._formatCtx.streams[i].contents
            codec = stream.codec.contents
            if codec.codec_type==decls.CODEC_TYPE_VIDEO:
                stream = VideoStream_ffmpeg(stream)
                self._streams.append(stream)
                self.videoStreams.append(stream)
            elif codec.codec_type==decls.CODEC_TYPE_AUDIO:
                stream = AudioStream_ffmpeg(stream)
                self._streams.append(stream)
                self.audioStreams.append(stream)
            else:
                self._streams.append(None)

    def __enter__(self):
        return self
    
    def __exit__(self, errorType, errorValue, traceback):
        self.close()
        return False
    
    def close(self):
        """Close the file.
        """
        for stream in self.audioStreams:
            stream.close()
        for stream in self.videoStreams:
            stream.close()
        
        self.audioStreams = []
        self.videoStreams = []
        
        if self._formatCtx is not None:
            avformat.av_close_input_file(self._formatCtx)
            self._formatCtx = None
    
    def numAudioStreams(self):
        """Return the number of audio streams defined in the file.
        """
        return len(self.audioStreams)
    
    def numVideoStreams(self):
        """Return the number of video streams defined in the file.
        """
        return len(self.videoStreams)

    def iterData(self, streams=None):
        """Iterate over the frames of the first video stream.
        """
        # Use the first video stream by default...
        if streams is None:
            if len(self.videoStreams)==0:
                return
            videoStream = self.videoStreams[0]
            streams = [videoStream]

        if len(streams)==0:
            return

        # Initialize the streams...
        # streamDict: Key:Stream index - Value:Stream object
        streamDict = {}
        for stream in streams:
            # Check the type of the stream
            if not isinstance(stream, StreamBase_ffmpeg):
                raise TypeError("Stream object expected")
            # Is this a stream from someone else?
            if stream not in self._streams:
                raise ValueError("Invalid stream (stream is from a different source)")
            
            stream.decodeBegin()
            streamDict[stream.index] = stream

        try:
            # Iterate over the packets in the file and decode the frames...
            for pkt in self.iterPackets():
                stream = streamDict.get(pkt.stream_index, None)
                if stream is not None:
                    for data in stream.handlePacket(pkt):
                        yield data
            
            # Some codecs have a delay between input and output, so go on
            # feeding empty packets until we receive no more frames.            
            pkt = decls.AVPacket()
            avcodec.av_init_packet(pkt)
            for stream in list(streamDict.values()):
                hasData = True
                while hasData: 
                    hasData = False
                    for data in stream.handlePacket(pkt):
                        hasData = True
                        yield data
        finally:
            for stream in list(streamDict.values()):
                stream.decodeEnd()
    
    def decode(self):
        """Decode stream data and pass it to the callbacks stored in the streams.
        """
        # Determine the streams that have a callback
        streams = []
        for stream in self._streams:
            cb = stream.getDataCallback()
            if cb is not None:
                streams.append(stream)

        # Iterate over the data...
        for data in self.iterData(streams):
            callback = data.stream.getDataCallback()
            callback(data)

    def iterPackets(self):
        """Iterate over all raw packets in the file.
        
        Yields AVPacket objects with the current packet. Actually, it always
        yields the same object that just has a new packet in it, so it's not
        valid to store the packet for later use.
        """
        formatCtx = self._formatCtx
        if formatCtx is None:
            return

        pkt = decls.AVPacket()
        eof = False
        while not eof:
            # Read the next frame packet...
            try:
                eof = not avformat.av_read_frame(formatCtx, pkt)
                if not eof:
                    yield pkt
            finally:
                avcodec.av_free_packet(pkt)

    @property
    def title(self):
        if self._formatCtx is None:
            return ""
        else:
            return self._formatCtx.title

    @property
    def author(self):
        if self._formatCtx is None:
            return ""
        else:
            return self._formatCtx.author

    @property
    def copyright(self):
        if self._formatCtx is None:
            return ""
        else:
            return self._formatCtx.copyright

    @property
    def album(self):
        if self._formatCtx is None:
            return ""
        else:
            return self._formatCtx.album

    @property
    def comment(self):
        if self._formatCtx is None:
            return ""
        else:
            return self._formatCtx.comment

    @property
    def year(self):
        if self._formatCtx is None:
            return 0
        else:
            return self._formatCtx.year

    @property
    def track(self):
        if self._formatCtx is None:
            return 0
        else:
            return self._formatCtx.track

    @property
    def genre(self):
        if self._formatCtx is None:
            return ""
        else:
            return self._formatCtx.genre



class VideoStream_Write_ffmpeg:
    def __init__(self, formatCtx, size, frameRate, pixelAspect):
        """Constructor.
        
        formatCtx is a AVFormatContext object.
        size is a tuple (width, height) containing the resolution of the
        video images.
        frameRate is the target framerate. It can be a single int or float
        or a tuple (num,den).
        pixelAspect is a float containing the pixel aspect ratio.
        """
        self._formatCtx = formatCtx
        self._size = size
        self._pixelAspect = pixelAspect
        
        try:
            # Check if frameRate is a 2-tuple
            num,den = frameRate
        except:
            # Turn the framerate into a rational...
            r = avutil.av_d2q(frameRate, 255)
            frameRate = (r.num, r.den)
        
        # The framerate as a tuple (num, den).
        self._frameRate = frameRate

        # Create an AVStream
        self._stream = self._createStream()
        
        # Get the codec context of the stream
        self._codecCtx = self._stream.codec.contents
        
        width,height = self._size
        # Set up the buffer that can store the encoded frame
        self._bufSize = 6*width*height+200
        self._buffer = ctypes.create_string_buffer(self._bufSize)

        self._frame = avcodec.avcodec_alloc_frame()
        if self._frame is None:
            raise MemoryError("Failed to allocate AVFrame object")
        avcodec.avpicture_alloc(self._frame, self._codecCtx.pix_fmt, width, height)

        # Allocate a picture for the converted image...
        self._picture = decls.AVPicture()
        avcodec.avpicture_alloc(self._picture, self._codecCtx.pix_fmt, width, height)
        # Get the size of the destination image buffer
        self._pictureSize = avcodec.avpicture_get_size(self._codecCtx.pix_fmt, width, height)
        
        self._pkt = decls.AVPacket()
        avformat.av_init_packet(self._pkt)
        #pkt.stream_index = outputstream.index
        
        self._currentPts = 0


    def close(self):
        pass
    
    def writeFrame(self, frame):
        buf = ctypes.addressof(self._buffer)
        self._frame.pts = self._currentPts
        self._currentPts += 1
        bytes = avcodec.avcodec_encode_video(self._codecCtx, buf, self._bufSize, self._frame)
#        print "pts", self._frame.pts
        print(("bytes %s"%bytes))
        if bytes>0:
            self._pkt.data = ctypes.cast(self._buffer, ctypes.POINTER(ctypes.c_uint8))
            self._pkt.size = bytes
            self._pkt.pts = self._frame.pts
#            if self._codecCtx.coded_frame.contents.key_frame:
#                self._pkt.flags |= cppdefs.PKT_FLAG_KEY
                
            avformat.av_interleaved_write_frame(self._formatCtx, self._pkt)

    def _createStream(self):
        """Create the AVStream and open the codec.
        """
        
        # Allocate a new stream
        stream = avformat.av_new_stream(self._formatCtx, self._formatCtx.nb_streams)
        # Initialize it as a video stream
        avcodec.avcodec_get_context_defaults2(stream.codec.contents, decls.CODEC_TYPE_VIDEO)

        outputFormat = self._formatCtx.oformat.contents
        codecCtx = stream.codec.contents
        
        if outputFormat.flags & decls.AVFMT_GLOBALHEADER:
            stream.codec.contents.flags |= decls.CODEC_FLAG_GLOBAL_HEADER
        
        # Prepare the codec...
        codecId = avformat.av_guess_codec(outputFormat, decls.CODEC_TYPE_VIDEO, fileName=self._formatCtx.filename)
        if codecId is None:
            raise ValueError("Could not determine video codec to use for encoding the video.")
        codec = avcodec.avcodec_find_encoder(codecId)
        if codec is None:
            raise ValueError("Could not find encoder.")

        if codec.supported_framerates:
            print ("mediafile: supported framerates available")
        else:
            print ("mediafile: supported framerates not available")

        codecCtx.codec_id = codecId
        num,den = self._frameRate
        codecCtx.time_base.den = num
        codecCtx.time_base.num = den
        width,height = self._size
        codecCtx.width = width
        codecCtx.height = height
        codecCtx.sample_aspect_ratio = avutil.av_d2q(self._pixelAspect, 255)
        codecCtx.pix_fmt = decls.PIX_FMT_YUV420P
        stream.sample_aspect_ratio = codecCtx.sample_aspect_ratio
        
        # Check if the codec supports the pixel format (if not, switch the format)
        if codec.pix_fmts:
            i = 0
            while codec.pix_fmts[i]!=-1:
                if codecCtx.pix_fmt==codec.pix_fmts[i]:
                    break
                i += 1
            else:
                codecCtx.pix_fmt = codec.pix_fmts[0]
        
        # Open the codec
        avcodec.avcodec_open(codecCtx, codec)
        
        return stream
    

class Media_Write_ffmpeg(object):
    """Media file writer.
    """
    
    def __init__(self, fileName):
        object.__init__(self)
        
        self._formatCtx = None
        self._headerWritten = False
        
        # A list of VideoStream_Write_ffmpeg objects
        self._videoStreams = []
        
        # Create and initialize the AVFormatContext structure...
        self._formatCtx = self._createFormatContext(fileName)
        
        # Open the file
        status = avformat.url_fopen(self._formatCtx.pb, fileName, decls.URL_WRONLY)
        if status<0:
            self._freeFormatContext()
            raise MediaFileError('Could not write file "%s"'%fileName)

    def close(self):
        """Close the file.
        """
        if self._formatCtx is not None:
            for stream in self._videoStreams:
                stream.close()
            self._videoStreams = []

            if self._formatCtx.pb:
                if self._headerWritten:
                    avformat.av_write_trailer(self._formatCtx)
                avformat.url_fclose(self._formatCtx.pb)
            self._freeFormatContext()

    def createVideoStream(self, size, frameRate=30, pixelAspect=1.0):
        """Add a new video stream to the file.
        """
        try:
            stream = VideoStream_Write_ffmpeg(self._formatCtx,
                                              size=size,
                                              frameRate=frameRate,
                                              pixelAspect=pixelAspect)
            self._videoStreams.append(stream)
        except:
            self.close()
            raise

    def writeFrame(self, frame):
        if not self._headerWritten:
            self._writeHeader()
        
        self._videoStreams[0].writeFrame(frame)
    
    def _writeHeader(self):
        avformat.av_write_header(self._formatCtx)
        self._headerWritten = True
    
    def _createFormatContext(self, fileName):
        """Create a new AVFormatContext and initialize it.
        
        fileName is the name of the output file (this will determine the format
        of the file).
        Returns a AVFormatContext object.
        """
        # Determine the output format (AVOutputFormat object) from the file name...
        fmt = avformat.av_guess_format(fileName=fileName)
        if fmt is None:
            raise MediaFileError('Cannot determine output format for file "%s"'%fileName)
        
        # Create a new AVFormatContext struct
        formatCtx = avformat.avformat_alloc_context()
        
        # Initialize some fields...
        formatCtx.oformat = ctypes.pointer(fmt)
        formatCtx.filename = fileName
        formatCtx.timestamp = 0          # Where should this timestamp come from?

        mux_preload = 0.5
        mux_max_delay = 0.7
        formatCtx.preload = int(mux_preload*cppdefs.AV_TIME_BASE)
        formatCtx.max_delay= (int)(mux_max_delay*cppdefs.AV_TIME_BASE)
        formatCtx.loop_output = cppdefs.AVFMT_NOOUTPUTLOOP
        formatCtx.flags |= cppdefs.AVFMT_FLAG_NONBLOCK
        
        # Do we have to call av_set_parameters() here? (this is done in ffmpeg.c)
        
        return formatCtx
    
    def _freeFormatContext(self):
        """Free the internal AVFormatContext struct.
        """
        if self._formatCtx is not None:
            avutil.av_free(self._formatCtx)
            self._formatCtx = None

    def writeFrame__dummy(self):
        # Create a packet struct and initialize it
        pkt = AVPacket()
        avcodec.av_init_packet(pkt)
        #pkt.stream_index = outputstream.index
        
        #big_picture.pts= ost->sync_opts;
        
        # Encode the frame
        bytes = avcodec.avcodec_encode_video(codecCtx, buf, bufsize, picture)
        if bytes>0:
            pkt.data = buf
            pkt.size = bytes
#            pkt.pts = ...
            if codecCtx.coded_frame.key_frame:
                pkt.flags |= PKT_FLAG_KEY
                
            avformat.av_interleaved_write_frame(formatCtx, pkt)


######################################################################

_av_is_registered = False

def open(name, mode="r", **keyargs):
    """Open an audio/video file.

    *name* is the name of the media file to read or write. *mode* determines
    whether an existing file will be read or a new file will be created. Valid
    values are ``"r"`` for reading a file and ``"w"`` for writing a file.
     
    Returns a :class:`Media_Read` or :class:`Media_Write` object representing
    the open file.
    """
    global _av_is_registered
    if not _av_is_registered:
        avformat.av_register_all()
        _av_is_registered = True

    if mode=="r":
        return Media_Read_ffmpeg(name, **keyargs)
    elif mode=="w":
        return Media_Write_ffmpeg(name, **keyargs)
    else:
        raise ValueError("Unknown mode '%s' (must be 'r' or 'w')"%mode)
        
