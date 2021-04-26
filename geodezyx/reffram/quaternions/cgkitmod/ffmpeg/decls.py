# ***** BEGIN LICENSE BLOCK *****
# Version: MPL 1.1/GPL 2.0/LGPL 2.1
#
# The contents of this file are subject to the Mozilla Public License Version
# 1.1 (the License); you may not use this file except in compliance with
# the License. You may obtain a copy of the License at
# http://www.mozilla.org/MPL/
#
# Software distributed under the License is distributed on an AS IS basis,
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
# either the GNU General Public License Version 2 or later (the GPL), or
# the GNU Lesser General Public License Version 2.1 or later (the LGPL),
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

from .cppdefs import *
from ctypes import *

STRING = c_char_p


class AVPicture(Structure):
    pass
uint8_t = c_uint8
AVPicture._fields_ = [
    ('data', POINTER(uint8_t) * 4),
    ('linesize', c_int * 4),
]
class AVFormatContext(Structure):
    pass
class AVClass(Structure):
    pass
class AVOption(Structure):
    pass

# values for enumeration 'AVOptionType'
FF_OPT_TYPE_FLAGS = 0
FF_OPT_TYPE_INT = 1
FF_OPT_TYPE_INT64 = 2
FF_OPT_TYPE_DOUBLE = 3
FF_OPT_TYPE_FLOAT = 4
FF_OPT_TYPE_STRING = 5
FF_OPT_TYPE_RATIONAL = 6
FF_OPT_TYPE_BINARY = 7
FF_OPT_TYPE_CONST = 128
AVOptionType = c_int # enum
class N8AVOption4DOLLAR_19E(Union):
    pass
int64_t = c_int64
class AVRational(Structure):
    pass
AVRational._fields_ = [
    ('num', c_int),
    ('den', c_int),
]
N8AVOption4DOLLAR_19E._pack_ = 4
N8AVOption4DOLLAR_19E._fields_ = [
    ('dbl', c_double),
    ('str', STRING),
    ('i64', int64_t),
    ('q', AVRational),
]
AVOption._pack_ = 4
AVOption._fields_ = [
    ('name', STRING),
    ('help', STRING),
    ('offset', c_int),
    ('type', AVOptionType),
    ('default_val', N8AVOption4DOLLAR_19E),
    ('min', c_double),
    ('max', c_double),
    ('flags', c_int),
    ('unit', STRING),
]
AVClass._fields_ = [
    ('class_name', STRING),
    ('item_name', CFUNCTYPE(STRING, c_void_p)),
    ('option', POINTER(AVOption)),
    ('version', c_int),
    ('log_level_offset_offset', c_int),
    ('parent_log_context_offset', c_int),
    ('opt_find', CFUNCTYPE(POINTER(AVOption), c_void_p, STRING, STRING, c_int, c_int)),
]
class AVInputFormat(Structure):
    pass
class AVOutputFormat(Structure):
    pass
class AVIOContext(Structure):
    pass
class AVStream(Structure):
    pass
class AVPacketList(Structure):
    pass
class AVPacket(Structure):
    pass
AVPacket._pack_ = 4
AVPacket._fields_ = [
    ('pts', int64_t),
    ('dts', int64_t),
    ('data', POINTER(uint8_t)),
    ('size', c_int),
    ('stream_index', c_int),
    ('flags', c_int),
    ('duration', c_int),
    ('destruct', CFUNCTYPE(None, POINTER(AVPacket))),
    ('priv', c_void_p),
    ('pos', int64_t),
    ('convergence_duration', int64_t),
]
class AVProgram(Structure):
    pass

# values for enumeration 'CodecID'
CODEC_ID_NONE = 0
CODEC_ID_MPEG1VIDEO = 1
CODEC_ID_MPEG2VIDEO = 2
CODEC_ID_MPEG2VIDEO_XVMC = 3
CODEC_ID_H261 = 4
CODEC_ID_H263 = 5
CODEC_ID_RV10 = 6
CODEC_ID_RV20 = 7
CODEC_ID_MJPEG = 8
CODEC_ID_MJPEGB = 9
CODEC_ID_LJPEG = 10
CODEC_ID_SP5X = 11
CODEC_ID_JPEGLS = 12
CODEC_ID_MPEG4 = 13
CODEC_ID_RAWVIDEO = 14
CODEC_ID_MSMPEG4V1 = 15
CODEC_ID_MSMPEG4V2 = 16
CODEC_ID_MSMPEG4V3 = 17
CODEC_ID_WMV1 = 18
CODEC_ID_WMV2 = 19
CODEC_ID_H263P = 20
CODEC_ID_H263I = 21
CODEC_ID_FLV1 = 22
CODEC_ID_SVQ1 = 23
CODEC_ID_SVQ3 = 24
CODEC_ID_DVVIDEO = 25
CODEC_ID_HUFFYUV = 26
CODEC_ID_CYUV = 27
CODEC_ID_H264 = 28
CODEC_ID_INDEO3 = 29
CODEC_ID_VP3 = 30
CODEC_ID_THEORA = 31
CODEC_ID_ASV1 = 32
CODEC_ID_ASV2 = 33
CODEC_ID_FFV1 = 34
CODEC_ID_4XM = 35
CODEC_ID_VCR1 = 36
CODEC_ID_CLJR = 37
CODEC_ID_MDEC = 38
CODEC_ID_ROQ = 39
CODEC_ID_INTERPLAY_VIDEO = 40
CODEC_ID_XAN_WC3 = 41
CODEC_ID_XAN_WC4 = 42
CODEC_ID_RPZA = 43
CODEC_ID_CINEPAK = 44
CODEC_ID_WS_VQA = 45
CODEC_ID_MSRLE = 46
CODEC_ID_MSVIDEO1 = 47
CODEC_ID_IDCIN = 48
CODEC_ID_8BPS = 49
CODEC_ID_SMC = 50
CODEC_ID_FLIC = 51
CODEC_ID_TRUEMOTION1 = 52
CODEC_ID_VMDVIDEO = 53
CODEC_ID_MSZH = 54
CODEC_ID_ZLIB = 55
CODEC_ID_QTRLE = 56
CODEC_ID_SNOW = 57
CODEC_ID_TSCC = 58
CODEC_ID_ULTI = 59
CODEC_ID_QDRAW = 60
CODEC_ID_VIXL = 61
CODEC_ID_QPEG = 62
CODEC_ID_XVID = 63
CODEC_ID_PNG = 64
CODEC_ID_PPM = 65
CODEC_ID_PBM = 66
CODEC_ID_PGM = 67
CODEC_ID_PGMYUV = 68
CODEC_ID_PAM = 69
CODEC_ID_FFVHUFF = 70
CODEC_ID_RV30 = 71
CODEC_ID_RV40 = 72
CODEC_ID_VC1 = 73
CODEC_ID_WMV3 = 74
CODEC_ID_LOCO = 75
CODEC_ID_WNV1 = 76
CODEC_ID_AASC = 77
CODEC_ID_INDEO2 = 78
CODEC_ID_FRAPS = 79
CODEC_ID_TRUEMOTION2 = 80
CODEC_ID_BMP = 81
CODEC_ID_CSCD = 82
CODEC_ID_MMVIDEO = 83
CODEC_ID_ZMBV = 84
CODEC_ID_AVS = 85
CODEC_ID_SMACKVIDEO = 86
CODEC_ID_NUV = 87
CODEC_ID_KMVC = 88
CODEC_ID_FLASHSV = 89
CODEC_ID_CAVS = 90
CODEC_ID_JPEG2000 = 91
CODEC_ID_VMNC = 92
CODEC_ID_VP5 = 93
CODEC_ID_VP6 = 94
CODEC_ID_VP6F = 95
CODEC_ID_TARGA = 96
CODEC_ID_DSICINVIDEO = 97
CODEC_ID_TIERTEXSEQVIDEO = 98
CODEC_ID_TIFF = 99
CODEC_ID_GIF = 100
CODEC_ID_FFH264 = 101
CODEC_ID_DXA = 102
CODEC_ID_DNXHD = 103
CODEC_ID_THP = 104
CODEC_ID_SGI = 105
CODEC_ID_C93 = 106
CODEC_ID_BETHSOFTVID = 107
CODEC_ID_PTX = 108
CODEC_ID_TXD = 109
CODEC_ID_VP6A = 110
CODEC_ID_AMV = 111
CODEC_ID_VB = 112
CODEC_ID_PCX = 113
CODEC_ID_SUNRAST = 114
CODEC_ID_INDEO4 = 115
CODEC_ID_INDEO5 = 116
CODEC_ID_MIMIC = 117
CODEC_ID_RL2 = 118
CODEC_ID_8SVX_EXP = 119
CODEC_ID_8SVX_FIB = 120
CODEC_ID_ESCAPE124 = 121
CODEC_ID_DIRAC = 122
CODEC_ID_BFI = 123
CODEC_ID_CMV = 124
CODEC_ID_MOTIONPIXELS = 125
CODEC_ID_TGV = 126
CODEC_ID_TGQ = 127
CODEC_ID_TQI = 128
CODEC_ID_AURA = 129
CODEC_ID_AURA2 = 130
CODEC_ID_V210X = 131
CODEC_ID_TMV = 132
CODEC_ID_V210 = 133
CODEC_ID_DPX = 134
CODEC_ID_MAD = 135
CODEC_ID_FRWU = 136
CODEC_ID_FLASHSV2 = 137
CODEC_ID_CDGRAPHICS = 138
CODEC_ID_R210 = 139
CODEC_ID_ANM = 140
CODEC_ID_BINKVIDEO = 141
CODEC_ID_IFF_ILBM = 142
CODEC_ID_IFF_BYTERUN1 = 143
CODEC_ID_KGV1 = 144
CODEC_ID_YOP = 145
CODEC_ID_VP8 = 146
CODEC_ID_PICTOR = 147
CODEC_ID_ANSI = 148
CODEC_ID_A64_MULTI = 149
CODEC_ID_A64_MULTI5 = 150
CODEC_ID_R10K = 151
CODEC_ID_MXPEG = 152
CODEC_ID_LAGARITH = 153
CODEC_ID_PRORES = 154
CODEC_ID_JV = 155
CODEC_ID_DFA = 156
CODEC_ID_8SVX_RAW = 157
CODEC_ID_PCM_S16LE = 65536
CODEC_ID_PCM_S16BE = 65537
CODEC_ID_PCM_U16LE = 65538
CODEC_ID_PCM_U16BE = 65539
CODEC_ID_PCM_S8 = 65540
CODEC_ID_PCM_U8 = 65541
CODEC_ID_PCM_MULAW = 65542
CODEC_ID_PCM_ALAW = 65543
CODEC_ID_PCM_S32LE = 65544
CODEC_ID_PCM_S32BE = 65545
CODEC_ID_PCM_U32LE = 65546
CODEC_ID_PCM_U32BE = 65547
CODEC_ID_PCM_S24LE = 65548
CODEC_ID_PCM_S24BE = 65549
CODEC_ID_PCM_U24LE = 65550
CODEC_ID_PCM_U24BE = 65551
CODEC_ID_PCM_S24DAUD = 65552
CODEC_ID_PCM_ZORK = 65553
CODEC_ID_PCM_S16LE_PLANAR = 65554
CODEC_ID_PCM_DVD = 65555
CODEC_ID_PCM_F32BE = 65556
CODEC_ID_PCM_F32LE = 65557
CODEC_ID_PCM_F64BE = 65558
CODEC_ID_PCM_F64LE = 65559
CODEC_ID_PCM_BLURAY = 65560
CODEC_ID_PCM_LXF = 65561
CODEC_ID_S302M = 65562
CODEC_ID_ADPCM_IMA_QT = 69632
CODEC_ID_ADPCM_IMA_WAV = 69633
CODEC_ID_ADPCM_IMA_DK3 = 69634
CODEC_ID_ADPCM_IMA_DK4 = 69635
CODEC_ID_ADPCM_IMA_WS = 69636
CODEC_ID_ADPCM_IMA_SMJPEG = 69637
CODEC_ID_ADPCM_MS = 69638
CODEC_ID_ADPCM_4XM = 69639
CODEC_ID_ADPCM_XA = 69640
CODEC_ID_ADPCM_ADX = 69641
CODEC_ID_ADPCM_EA = 69642
CODEC_ID_ADPCM_G726 = 69643
CODEC_ID_ADPCM_CT = 69644
CODEC_ID_ADPCM_SWF = 69645
CODEC_ID_ADPCM_YAMAHA = 69646
CODEC_ID_ADPCM_SBPRO_4 = 69647
CODEC_ID_ADPCM_SBPRO_3 = 69648
CODEC_ID_ADPCM_SBPRO_2 = 69649
CODEC_ID_ADPCM_THP = 69650
CODEC_ID_ADPCM_IMA_AMV = 69651
CODEC_ID_ADPCM_EA_R1 = 69652
CODEC_ID_ADPCM_EA_R3 = 69653
CODEC_ID_ADPCM_EA_R2 = 69654
CODEC_ID_ADPCM_IMA_EA_SEAD = 69655
CODEC_ID_ADPCM_IMA_EA_EACS = 69656
CODEC_ID_ADPCM_EA_XAS = 69657
CODEC_ID_ADPCM_EA_MAXIS_XA = 69658
CODEC_ID_ADPCM_IMA_ISS = 69659
CODEC_ID_ADPCM_G722 = 69660
CODEC_ID_AMR_NB = 73728
CODEC_ID_AMR_WB = 73729
CODEC_ID_RA_144 = 77824
CODEC_ID_RA_288 = 77825
CODEC_ID_ROQ_DPCM = 81920
CODEC_ID_INTERPLAY_DPCM = 81921
CODEC_ID_XAN_DPCM = 81922
CODEC_ID_SOL_DPCM = 81923
CODEC_ID_MP2 = 86016
CODEC_ID_MP3 = 86017
CODEC_ID_AAC = 86018
CODEC_ID_AC3 = 86019
CODEC_ID_DTS = 86020
CODEC_ID_VORBIS = 86021
CODEC_ID_DVAUDIO = 86022
CODEC_ID_WMAV1 = 86023
CODEC_ID_WMAV2 = 86024
CODEC_ID_MACE3 = 86025
CODEC_ID_MACE6 = 86026
CODEC_ID_VMDAUDIO = 86027
CODEC_ID_SONIC = 86028
CODEC_ID_SONIC_LS = 86029
CODEC_ID_FLAC = 86030
CODEC_ID_MP3ADU = 86031
CODEC_ID_MP3ON4 = 86032
CODEC_ID_SHORTEN = 86033
CODEC_ID_ALAC = 86034
CODEC_ID_WESTWOOD_SND1 = 86035
CODEC_ID_GSM = 86036
CODEC_ID_QDM2 = 86037
CODEC_ID_COOK = 86038
CODEC_ID_TRUESPEECH = 86039
CODEC_ID_TTA = 86040
CODEC_ID_SMACKAUDIO = 86041
CODEC_ID_QCELP = 86042
CODEC_ID_WAVPACK = 86043
CODEC_ID_DSICINAUDIO = 86044
CODEC_ID_IMC = 86045
CODEC_ID_MUSEPACK7 = 86046
CODEC_ID_MLP = 86047
CODEC_ID_GSM_MS = 86048
CODEC_ID_ATRAC3 = 86049
CODEC_ID_VOXWARE = 86050
CODEC_ID_APE = 86051
CODEC_ID_NELLYMOSER = 86052
CODEC_ID_MUSEPACK8 = 86053
CODEC_ID_SPEEX = 86054
CODEC_ID_WMAVOICE = 86055
CODEC_ID_WMAPRO = 86056
CODEC_ID_WMALOSSLESS = 86057
CODEC_ID_ATRAC3P = 86058
CODEC_ID_EAC3 = 86059
CODEC_ID_SIPR = 86060
CODEC_ID_MP1 = 86061
CODEC_ID_TWINVQ = 86062
CODEC_ID_TRUEHD = 86063
CODEC_ID_MP4ALS = 86064
CODEC_ID_ATRAC1 = 86065
CODEC_ID_BINKAUDIO_RDFT = 86066
CODEC_ID_BINKAUDIO_DCT = 86067
CODEC_ID_AAC_LATM = 86068
CODEC_ID_QDMC = 86069
CODEC_ID_CELT = 86070
CODEC_ID_DVD_SUBTITLE = 94208
CODEC_ID_DVB_SUBTITLE = 94209
CODEC_ID_TEXT = 94210
CODEC_ID_XSUB = 94211
CODEC_ID_SSA = 94212
CODEC_ID_MOV_TEXT = 94213
CODEC_ID_HDMV_PGS_SUBTITLE = 94214
CODEC_ID_DVB_TELETEXT = 94215
CODEC_ID_SRT = 94216
CODEC_ID_MICRODVD = 94217
CODEC_ID_TTF = 98304
CODEC_ID_PROBE = 102400
CODEC_ID_MPEG2TS = 131072
CODEC_ID_FFMETADATA = 135168
CodecID = c_int # enum
class AVChapter(Structure):
    pass
class AVDictionary(Structure):
    pass
AVFormatContext._pack_ = 4
AVFormatContext._fields_ = [
    ('av_class', POINTER(AVClass)),
    ('iformat', POINTER(AVInputFormat)),
    ('oformat', POINTER(AVOutputFormat)),
    ('priv_data', c_void_p),
    ('pb', POINTER(AVIOContext)),
    ('nb_streams', c_uint),
    ('streams', POINTER(AVStream) * 20),
    ('filename', c_char * 1024),
    ('timestamp', int64_t),
    ('title', c_char * 512),
    ('author', c_char * 512),
    ('copyright', c_char * 512),
    ('comment', c_char * 512),
    ('album', c_char * 512),
    ('year', c_int),
    ('track', c_int),
    ('genre', c_char * 32),
    ('ctx_flags', c_int),
    ('packet_buffer', POINTER(AVPacketList)),
    ('start_time', int64_t),
    ('duration', int64_t),
    ('file_size', int64_t),
    ('bit_rate', c_int),
    ('cur_st', POINTER(AVStream)),
    ('cur_ptr_deprecated', POINTER(uint8_t)),
    ('cur_len_deprecated', c_int),
    ('cur_pkt_deprecated', AVPacket),
    ('data_offset', int64_t),
    ('index_built', c_int),
    ('mux_rate', c_int),
    ('packet_size', c_uint),
    ('preload', c_int),
    ('max_delay', c_int),
    ('loop_output', c_int),
    ('flags', c_int),
    ('loop_input', c_int),
    ('probesize', c_uint),
    ('max_analyze_duration', c_int),
    ('key', POINTER(uint8_t)),
    ('keylen', c_int),
    ('nb_programs', c_uint),
    ('programs', POINTER(POINTER(AVProgram))),
    ('video_codec_id', CodecID),
    ('audio_codec_id', CodecID),
    ('subtitle_codec_id', CodecID),
    ('max_index_size', c_uint),
    ('max_picture_buffer', c_uint),
    ('nb_chapters', c_uint),
    ('chapters', POINTER(POINTER(AVChapter))),
    ('debug', c_int),
    ('raw_packet_buffer', POINTER(AVPacketList)),
    ('raw_packet_buffer_end', POINTER(AVPacketList)),
    ('packet_buffer_end', POINTER(AVPacketList)),
    ('metadata', POINTER(AVDictionary)),
    ('raw_packet_buffer_remaining_size', c_int),
    ('start_time_realtime', int64_t),
    ('fps_probe_size', c_int),
    ('ts_id', c_int),
]
class SwsContext(Structure):
    pass
SwsContext._fields_ = [
]
class AVCodecTag(Structure):
    pass
AVCodecTag._fields_ = [
]
class AVMetadataConv(Structure):
    pass
AVOutputFormat._fields_ = [
    ('name', STRING),
    ('long_name', STRING),
    ('mime_type', STRING),
    ('extensions', STRING),
    ('priv_data_size', c_int),
    ('audio_codec', CodecID),
    ('video_codec', CodecID),
    ('write_header', CFUNCTYPE(c_int, POINTER(AVFormatContext))),
    ('write_packet', CFUNCTYPE(c_int, POINTER(AVFormatContext), POINTER(AVPacket))),
    ('write_trailer', CFUNCTYPE(c_int, POINTER(AVFormatContext))),
    ('flags', c_int),
    ('dummy', c_void_p),
    ('interleave_packet', CFUNCTYPE(c_int, POINTER(AVFormatContext), POINTER(AVPacket), POINTER(AVPacket), c_int)),
    ('codec_tag', POINTER(POINTER(AVCodecTag))),
    ('subtitle_codec', CodecID),
    ('metadata_conv', POINTER(AVMetadataConv)),
    ('priv_class', POINTER(AVClass)),
    ('next', POINTER(AVOutputFormat)),
]
class AVProbeData(Structure):
    pass
class AVFormatParameters(Structure):
    pass
AVInputFormat._fields_ = [
    ('name', STRING),
    ('long_name', STRING),
    ('priv_data_size', c_int),
    ('read_probe', CFUNCTYPE(c_int, POINTER(AVProbeData))),
    ('read_header', CFUNCTYPE(c_int, POINTER(AVFormatContext), POINTER(AVFormatParameters))),
    ('read_packet', CFUNCTYPE(c_int, POINTER(AVFormatContext), POINTER(AVPacket))),
    ('read_close', CFUNCTYPE(c_int, POINTER(AVFormatContext))),
    ('read_seek', CFUNCTYPE(c_int, POINTER(AVFormatContext), c_int, int64_t, c_int)),
    ('read_timestamp', CFUNCTYPE(int64_t, POINTER(AVFormatContext), c_int, POINTER(int64_t), int64_t)),
    ('flags', c_int),
    ('extensions', STRING),
    ('value', c_int),
    ('read_play', CFUNCTYPE(c_int, POINTER(AVFormatContext))),
    ('read_pause', CFUNCTYPE(c_int, POINTER(AVFormatContext))),
    ('codec_tag', POINTER(POINTER(AVCodecTag))),
    ('read_seek2', CFUNCTYPE(c_int, POINTER(AVFormatContext), c_int, int64_t, int64_t, int64_t, c_int)),
    ('metadata_conv', POINTER(AVMetadataConv)),
    ('priv_class', POINTER(AVClass)),
    ('next', POINTER(AVInputFormat)),
]
class AVCodecContext(Structure):
    pass
class AVFrac(Structure):
    pass
AVFrac._pack_ = 4
AVFrac._fields_ = [
    ('val', int64_t),
    ('num', int64_t),
    ('den', int64_t),
]

# values for enumeration 'AVDiscard'
AVDISCARD_NONE = -16
AVDISCARD_DEFAULT = 0
AVDISCARD_NONREF = 8
AVDISCARD_BIDIR = 16
AVDISCARD_NONKEY = 32
AVDISCARD_ALL = 48
AVDiscard = c_int # enum

# values for enumeration 'AVStreamParseType'
AVSTREAM_PARSE_NONE = 0
AVSTREAM_PARSE_FULL = 1
AVSTREAM_PARSE_HEADERS = 2
AVSTREAM_PARSE_TIMESTAMPS = 3
AVSTREAM_PARSE_FULL_ONCE = 4
AVStreamParseType = c_int # enum
class AVCodecParserContext(Structure):
    pass
class AVIndexEntry(Structure):
    pass
AVProbeData._fields_ = [
    ('filename', STRING),
    ('buf', POINTER(c_ubyte)),
    ('buf_size', c_int),
]
class N8AVStream4DOLLAR_21E(Structure):
    pass
AVStream._pack_ = 4
AVStream._fields_ = [
    ('index', c_int),
    ('id', c_int),
    ('codec', POINTER(AVCodecContext)),
    ('r_frame_rate', AVRational),
    ('priv_data', c_void_p),
    ('first_dts', int64_t),
    ('pts', AVFrac),
    ('time_base', AVRational),
    ('pts_wrap_bits', c_int),
    ('stream_copy', c_int),
    ('discard', AVDiscard),
    ('quality', c_float),
    ('start_time', int64_t),
    ('duration', int64_t),
    ('language', c_char * 4),
    ('need_parsing', AVStreamParseType),
    ('parser', POINTER(AVCodecParserContext)),
    ('cur_dts', int64_t),
    ('last_IP_duration', c_int),
    ('last_IP_pts', int64_t),
    ('index_entries', POINTER(AVIndexEntry)),
    ('nb_index_entries', c_int),
    ('index_entries_allocated_size', c_uint),
    ('nb_frames', int64_t),
    ('unused', int64_t * 5),
    ('filename', STRING),
    ('disposition', c_int),
    ('probe_data', AVProbeData),
    ('pts_buffer', int64_t * 17),
    ('sample_aspect_ratio', AVRational),
    ('metadata', POINTER(AVDictionary)),
    ('cur_ptr', POINTER(uint8_t)),
    ('cur_len', c_int),
    ('cur_pkt', AVPacket),
    ('reference_dts', int64_t),
    ('probe_packets', c_int),
    ('last_in_packet_buffer', POINTER(AVPacketList)),
    ('avg_frame_rate', AVRational),
    ('codec_info_nb_frames', c_int),
    ('stream_identifier', c_int),
    ('info', POINTER(N8AVStream4DOLLAR_21E)),
    ('request_probe', c_int),
]
AVProgram._fields_ = [
    ('id', c_int),
    ('provider_name', STRING),
    ('name', STRING),
    ('flags', c_int),
    ('discard', AVDiscard),
    ('stream_index', POINTER(c_uint)),
    ('nb_stream_indexes', c_uint),
    ('metadata', POINTER(AVDictionary)),
    ('program_num', c_int),
    ('pmt_pid', c_int),
    ('pcr_pid', c_int),
]
AVChapter._pack_ = 4
AVChapter._fields_ = [
    ('id', c_int),
    ('time_base', AVRational),
    ('start', int64_t),
    ('end', int64_t),
    ('title', STRING),
    ('metadata', POINTER(AVDictionary)),
]
AVPacketList._fields_ = [
    ('pkt', AVPacket),
    ('next', POINTER(AVPacketList)),
]
AVIOContext._pack_ = 4
AVIOContext._fields_ = [
    ('buffer', POINTER(c_ubyte)),
    ('buffer_size', c_int),
    ('buf_ptr', POINTER(c_ubyte)),
    ('buf_end', POINTER(c_ubyte)),
    ('opaque', c_void_p),
    ('read_packet', CFUNCTYPE(c_int, c_void_p, POINTER(uint8_t), c_int)),
    ('write_packet', CFUNCTYPE(c_int, c_void_p, POINTER(uint8_t), c_int)),
    ('seek', CFUNCTYPE(int64_t, c_void_p, int64_t, c_int)),
    ('pos', int64_t),
    ('must_flush', c_int),
    ('eof_reached', c_int),
    ('write_flag', c_int),
    ('is_streamed', c_int),
    ('max_packet_size', c_int),
    ('checksum', c_ulong),
    ('checksum_ptr', POINTER(c_ubyte)),
    ('update_checksum', CFUNCTYPE(c_ulong, c_ulong, POINTER(uint8_t), c_uint)),
    ('error', c_int),
    ('read_pause', CFUNCTYPE(c_int, c_void_p, c_int)),
    ('read_seek', CFUNCTYPE(int64_t, c_void_p, c_int, int64_t, c_int)),
    ('seekable', c_int),
]
AVDictionary._fields_ = [
]

# values for enumeration 'PixelFormat'
PIX_FMT_NONE = -1
PIX_FMT_YUV420P = 0
PIX_FMT_YUYV422 = 1
PIX_FMT_RGB24 = 2
PIX_FMT_BGR24 = 3
PIX_FMT_YUV422P = 4
PIX_FMT_YUV444P = 5
PIX_FMT_YUV410P = 6
PIX_FMT_YUV411P = 7
PIX_FMT_GRAY8 = 8
PIX_FMT_MONOWHITE = 9
PIX_FMT_MONOBLACK = 10
PIX_FMT_PAL8 = 11
PIX_FMT_YUVJ420P = 12
PIX_FMT_YUVJ422P = 13
PIX_FMT_YUVJ444P = 14
PIX_FMT_XVMC_MPEG2_MC = 15
PIX_FMT_XVMC_MPEG2_IDCT = 16
PIX_FMT_UYVY422 = 17
PIX_FMT_UYYVYY411 = 18
PIX_FMT_BGR8 = 19
PIX_FMT_BGR4 = 20
PIX_FMT_BGR4_BYTE = 21
PIX_FMT_RGB8 = 22
PIX_FMT_RGB4 = 23
PIX_FMT_RGB4_BYTE = 24
PIX_FMT_NV12 = 25
PIX_FMT_NV21 = 26
PIX_FMT_ARGB = 27
PIX_FMT_RGBA = 28
PIX_FMT_ABGR = 29
PIX_FMT_BGRA = 30
PIX_FMT_GRAY16BE = 31
PIX_FMT_GRAY16LE = 32
PIX_FMT_YUV440P = 33
PIX_FMT_YUVJ440P = 34
PIX_FMT_YUVA420P = 35
PIX_FMT_VDPAU_H264 = 36
PIX_FMT_VDPAU_MPEG1 = 37
PIX_FMT_VDPAU_MPEG2 = 38
PIX_FMT_VDPAU_WMV3 = 39
PIX_FMT_VDPAU_VC1 = 40
PIX_FMT_RGB48BE = 41
PIX_FMT_RGB48LE = 42
PIX_FMT_RGB565BE = 43
PIX_FMT_RGB565LE = 44
PIX_FMT_RGB555BE = 45
PIX_FMT_RGB555LE = 46
PIX_FMT_BGR565BE = 47
PIX_FMT_BGR565LE = 48
PIX_FMT_BGR555BE = 49
PIX_FMT_BGR555LE = 50
PIX_FMT_VAAPI_MOCO = 51
PIX_FMT_VAAPI_IDCT = 52
PIX_FMT_VAAPI_VLD = 53
PIX_FMT_YUV420P16LE = 54
PIX_FMT_YUV420P16BE = 55
PIX_FMT_YUV422P16LE = 56
PIX_FMT_YUV422P16BE = 57
PIX_FMT_YUV444P16LE = 58
PIX_FMT_YUV444P16BE = 59
PIX_FMT_VDPAU_MPEG4 = 60
PIX_FMT_DXVA2_VLD = 61
PIX_FMT_RGB444LE = 62
PIX_FMT_RGB444BE = 63
PIX_FMT_BGR444LE = 64
PIX_FMT_BGR444BE = 65
PIX_FMT_GRAY8A = 66
PIX_FMT_BGR48BE = 67
PIX_FMT_BGR48LE = 68
PIX_FMT_YUV420P9BE = 69
PIX_FMT_YUV420P9LE = 70
PIX_FMT_YUV420P10BE = 71
PIX_FMT_YUV420P10LE = 72
PIX_FMT_YUV422P10BE = 73
PIX_FMT_YUV422P10LE = 74
PIX_FMT_YUV444P9BE = 75
PIX_FMT_YUV444P9LE = 76
PIX_FMT_YUV444P10BE = 77
PIX_FMT_YUV444P10LE = 78
PIX_FMT_NB = 79
PixelFormat = c_int # enum
class AVFrame(Structure):
    pass

# values for enumeration 'AVSampleFormat'
AV_SAMPLE_FMT_NONE = -1
AV_SAMPLE_FMT_U8 = 0
AV_SAMPLE_FMT_S16 = 1
AV_SAMPLE_FMT_S32 = 2
AV_SAMPLE_FMT_FLT = 3
AV_SAMPLE_FMT_DBL = 4
AV_SAMPLE_FMT_NB = 5
AVSampleFormat = c_int # enum
class AVCodec(Structure):
    pass

# values for enumeration 'AVMediaType'
AVMEDIA_TYPE_UNKNOWN = -1
AVMEDIA_TYPE_VIDEO = 0
AVMEDIA_TYPE_AUDIO = 1
AVMEDIA_TYPE_DATA = 2
AVMEDIA_TYPE_SUBTITLE = 3
AVMEDIA_TYPE_ATTACHMENT = 4
AVMEDIA_TYPE_NB = 5
AVMediaType = c_int # enum
class RcOverride(Structure):
    pass
uint64_t = c_uint64
uint16_t = c_uint16
class AVPaletteControl(Structure):
    pass
class AVHWAccel(Structure):
    pass

# values for enumeration 'AVColorPrimaries'
AVCOL_PRI_BT709 = 1
AVCOL_PRI_UNSPECIFIED = 2
AVCOL_PRI_BT470M = 4
AVCOL_PRI_BT470BG = 5
AVCOL_PRI_SMPTE170M = 6
AVCOL_PRI_SMPTE240M = 7
AVCOL_PRI_FILM = 8
AVCOL_PRI_NB = 9
AVColorPrimaries = c_int # enum

# values for enumeration 'AVColorTransferCharacteristic'
AVCOL_TRC_BT709 = 1
AVCOL_TRC_UNSPECIFIED = 2
AVCOL_TRC_GAMMA22 = 4
AVCOL_TRC_GAMMA28 = 5
AVCOL_TRC_NB = 6
AVColorTransferCharacteristic = c_int # enum

# values for enumeration 'AVColorSpace'
AVCOL_SPC_RGB = 0
AVCOL_SPC_BT709 = 1
AVCOL_SPC_UNSPECIFIED = 2
AVCOL_SPC_FCC = 4
AVCOL_SPC_BT470BG = 5
AVCOL_SPC_SMPTE170M = 6
AVCOL_SPC_SMPTE240M = 7
AVCOL_SPC_NB = 8
AVColorSpace = c_int # enum

# values for enumeration 'AVColorRange'
AVCOL_RANGE_UNSPECIFIED = 0
AVCOL_RANGE_MPEG = 1
AVCOL_RANGE_JPEG = 2
AVCOL_RANGE_NB = 3
AVColorRange = c_int # enum

# values for enumeration 'AVChromaLocation'
AVCHROMA_LOC_UNSPECIFIED = 0
AVCHROMA_LOC_LEFT = 1
AVCHROMA_LOC_CENTER = 2
AVCHROMA_LOC_TOPLEFT = 3
AVCHROMA_LOC_TOP = 4
AVCHROMA_LOC_BOTTOMLEFT = 5
AVCHROMA_LOC_BOTTOM = 6
AVCHROMA_LOC_NB = 7
AVChromaLocation = c_int # enum

# values for enumeration 'AVLPCType'
AV_LPC_TYPE_DEFAULT = -1
AV_LPC_TYPE_NONE = 0
AV_LPC_TYPE_FIXED = 1
AV_LPC_TYPE_LEVINSON = 2
AV_LPC_TYPE_CHOLESKY = 3
AV_LPC_TYPE_NB = 4
AVLPCType = c_int # enum

# values for enumeration 'AVAudioServiceType'
AV_AUDIO_SERVICE_TYPE_MAIN = 0
AV_AUDIO_SERVICE_TYPE_EFFECTS = 1
AV_AUDIO_SERVICE_TYPE_VISUALLY_IMPAIRED = 2
AV_AUDIO_SERVICE_TYPE_HEARING_IMPAIRED = 3
AV_AUDIO_SERVICE_TYPE_DIALOGUE = 4
AV_AUDIO_SERVICE_TYPE_COMMENTARY = 5
AV_AUDIO_SERVICE_TYPE_EMERGENCY = 6
AV_AUDIO_SERVICE_TYPE_VOICE_OVER = 7
AV_AUDIO_SERVICE_TYPE_KARAOKE = 8
AV_AUDIO_SERVICE_TYPE_NB = 9
AVAudioServiceType = c_int # enum
AVCodecContext._pack_ = 4
AVCodecContext._fields_ = [
    ('av_class', POINTER(AVClass)),
    ('bit_rate', c_int),
    ('bit_rate_tolerance', c_int),
    ('flags', c_int),
    ('sub_id', c_int),
    ('me_method', c_int),
    ('extradata', POINTER(uint8_t)),
    ('extradata_size', c_int),
    ('time_base', AVRational),
    ('width', c_int),
    ('height', c_int),
    ('gop_size', c_int),
    ('pix_fmt', PixelFormat),
    ('rate_emu', c_int),
    ('draw_horiz_band', CFUNCTYPE(None, POINTER(AVCodecContext), POINTER(AVFrame), POINTER(c_int), c_int, c_int, c_int)),
    ('sample_rate', c_int),
    ('channels', c_int),
    ('sample_fmt', AVSampleFormat),
    ('frame_size', c_int),
    ('frame_number', c_int),
    ('real_pict_num', c_int),
    ('delay', c_int),
    ('qcompress', c_float),
    ('qblur', c_float),
    ('qmin', c_int),
    ('qmax', c_int),
    ('max_qdiff', c_int),
    ('max_b_frames', c_int),
    ('b_quant_factor', c_float),
    ('rc_strategy', c_int),
    ('b_frame_strategy', c_int),
    ('hurry_up', c_int),
    ('codec', POINTER(AVCodec)),
    ('priv_data', c_void_p),
    ('rtp_payload_size', c_int),
    ('rtp_callback', CFUNCTYPE(None, POINTER(AVCodecContext), c_void_p, c_int, c_int)),
    ('mv_bits', c_int),
    ('header_bits', c_int),
    ('i_tex_bits', c_int),
    ('p_tex_bits', c_int),
    ('i_count', c_int),
    ('p_count', c_int),
    ('skip_count', c_int),
    ('misc_bits', c_int),
    ('frame_bits', c_int),
    ('opaque', c_void_p),
    ('codec_name', c_char * 32),
    ('codec_type', AVMediaType),
    ('codec_id', CodecID),
    ('codec_tag', c_uint),
    ('workaround_bugs', c_int),
    ('luma_elim_threshold', c_int),
    ('chroma_elim_threshold', c_int),
    ('strict_std_compliance', c_int),
    ('b_quant_offset', c_float),
    ('error_recognition', c_int),
    ('get_buffer', CFUNCTYPE(c_int, POINTER(AVCodecContext), POINTER(AVFrame))),
    ('release_buffer', CFUNCTYPE(None, POINTER(AVCodecContext), POINTER(AVFrame))),
    ('has_b_frames', c_int),
    ('block_align', c_int),
    ('parse_only', c_int),
    ('mpeg_quant', c_int),
    ('stats_out', STRING),
    ('stats_in', STRING),
    ('rc_qsquish', c_float),
    ('rc_qmod_amp', c_float),
    ('rc_qmod_freq', c_int),
    ('rc_override', POINTER(RcOverride)),
    ('rc_override_count', c_int),
    ('rc_eq', STRING),
    ('rc_max_rate', c_int),
    ('rc_min_rate', c_int),
    ('rc_buffer_size', c_int),
    ('rc_buffer_aggressivity', c_float),
    ('i_quant_factor', c_float),
    ('i_quant_offset', c_float),
    ('rc_initial_cplx', c_float),
    ('dct_algo', c_int),
    ('lumi_masking', c_float),
    ('temporal_cplx_masking', c_float),
    ('spatial_cplx_masking', c_float),
    ('p_masking', c_float),
    ('dark_masking', c_float),
    ('idct_algo', c_int),
    ('slice_count', c_int),
    ('slice_offset', POINTER(c_int)),
    ('error_concealment', c_int),
    ('dsp_mask', c_uint),
    ('bits_per_coded_sample', c_int),
    ('prediction_method', c_int),
    ('sample_aspect_ratio', AVRational),
    ('coded_frame', POINTER(AVFrame)),
    ('debug', c_int),
    ('debug_mv', c_int),
    ('error', uint64_t * 4),
    ('mb_qmin', c_int),
    ('mb_qmax', c_int),
    ('me_cmp', c_int),
    ('me_sub_cmp', c_int),
    ('mb_cmp', c_int),
    ('ildct_cmp', c_int),
    ('dia_size', c_int),
    ('last_predictor_count', c_int),
    ('pre_me', c_int),
    ('me_pre_cmp', c_int),
    ('pre_dia_size', c_int),
    ('me_subpel_quality', c_int),
    ('get_format', CFUNCTYPE(PixelFormat, POINTER(AVCodecContext), POINTER(PixelFormat))),
    ('dtg_active_format', c_int),
    ('me_range', c_int),
    ('intra_quant_bias', c_int),
    ('inter_quant_bias', c_int),
    ('color_table_id', c_int),
    ('internal_buffer_count', c_int),
    ('internal_buffer', c_void_p),
    ('global_quality', c_int),
    ('coder_type', c_int),
    ('context_model', c_int),
    ('slice_flags', c_int),
    ('xvmc_acceleration', c_int),
    ('mb_decision', c_int),
    ('intra_matrix', POINTER(uint16_t)),
    ('inter_matrix', POINTER(uint16_t)),
    ('stream_codec_tag', c_uint),
    ('scenechange_threshold', c_int),
    ('lmin', c_int),
    ('lmax', c_int),
    ('palctrl', POINTER(AVPaletteControl)),
    ('noise_reduction', c_int),
    ('reget_buffer', CFUNCTYPE(c_int, POINTER(AVCodecContext), POINTER(AVFrame))),
    ('rc_initial_buffer_occupancy', c_int),
    ('inter_threshold', c_int),
    ('flags2', c_int),
    ('error_rate', c_int),
    ('antialias_algo', c_int),
    ('quantizer_noise_shaping', c_int),
    ('thread_count', c_int),
    ('execute', CFUNCTYPE(c_int, POINTER(AVCodecContext), CFUNCTYPE(c_int, POINTER(AVCodecContext), c_void_p), c_void_p, POINTER(c_int), c_int, c_int)),
    ('thread_opaque', c_void_p),
    ('me_threshold', c_int),
    ('mb_threshold', c_int),
    ('intra_dc_precision', c_int),
    ('nsse_weight', c_int),
    ('skip_top', c_int),
    ('skip_bottom', c_int),
    ('profile', c_int),
    ('level', c_int),
    ('lowres', c_int),
    ('coded_width', c_int),
    ('coded_height', c_int),
    ('frame_skip_threshold', c_int),
    ('frame_skip_factor', c_int),
    ('frame_skip_exp', c_int),
    ('frame_skip_cmp', c_int),
    ('border_masking', c_float),
    ('mb_lmin', c_int),
    ('mb_lmax', c_int),
    ('me_penalty_compensation', c_int),
    ('skip_loop_filter', AVDiscard),
    ('skip_idct', AVDiscard),
    ('skip_frame', AVDiscard),
    ('bidir_refine', c_int),
    ('brd_scale', c_int),
    ('crf', c_float),
    ('cqp', c_int),
    ('keyint_min', c_int),
    ('refs', c_int),
    ('chromaoffset', c_int),
    ('bframebias', c_int),
    ('trellis', c_int),
    ('complexityblur', c_float),
    ('deblockalpha', c_int),
    ('deblockbeta', c_int),
    ('partitions', c_int),
    ('directpred', c_int),
    ('cutoff', c_int),
    ('scenechange_factor', c_int),
    ('mv0_threshold', c_int),
    ('b_sensitivity', c_int),
    ('compression_level', c_int),
    ('use_lpc', c_int),
    ('lpc_coeff_precision', c_int),
    ('min_prediction_order', c_int),
    ('max_prediction_order', c_int),
    ('prediction_order_method', c_int),
    ('min_partition_order', c_int),
    ('max_partition_order', c_int),
    ('timecode_frame_start', int64_t),
    ('request_channels', c_int),
    ('drc_scale', c_float),
    ('reordered_opaque', int64_t),
    ('bits_per_raw_sample', c_int),
    ('channel_layout', int64_t),
    ('request_channel_layout', int64_t),
    ('rc_max_available_vbv_use', c_float),
    ('rc_min_vbv_overflow_use', c_float),
    ('hwaccel', POINTER(AVHWAccel)),
    ('ticks_per_frame', c_int),
    ('hwaccel_context', c_void_p),
    ('color_primaries', AVColorPrimaries),
    ('color_trc', AVColorTransferCharacteristic),
    ('colorspace', AVColorSpace),
    ('color_range', AVColorRange),
    ('chroma_sample_location', AVChromaLocation),
    ('execute2', CFUNCTYPE(c_int, POINTER(AVCodecContext), CFUNCTYPE(c_int, POINTER(AVCodecContext), c_void_p, c_int, c_int), c_void_p, POINTER(c_int), c_int)),
    ('weighted_p_pred', c_int),
    ('aq_mode', c_int),
    ('aq_strength', c_float),
    ('psy_rd', c_float),
    ('psy_trellis', c_float),
    ('rc_lookahead', c_int),
    ('crf_max', c_float),
    ('log_level_offset', c_int),
    ('lpc_type', AVLPCType),
    ('lpc_passes', c_int),
    ('slices', c_int),
    ('subtitle_header', POINTER(uint8_t)),
    ('subtitle_header_size', c_int),
    ('pkt', POINTER(AVPacket)),
    ('is_copy', c_int),
    ('thread_type', c_int),
    ('active_thread_type', c_int),
    ('thread_safe_callbacks', c_int),
    ('vbv_delay', uint64_t),
    ('audio_service_type', AVAudioServiceType),
    ('request_sample_fmt', AVSampleFormat),
    ('pts_correction_num_faulty_pts', int64_t),
    ('pts_correction_num_faulty_dts', int64_t),
    ('pts_correction_last_pts', int64_t),
    ('pts_correction_last_dts', int64_t),
]
class AVCodecParser(Structure):
    pass
AVCodecParserContext._pack_ = 4
AVCodecParserContext._fields_ = [
    ('priv_data', c_void_p),
    ('parser', POINTER(AVCodecParser)),
    ('frame_offset', int64_t),
    ('cur_offset', int64_t),
    ('next_frame_offset', int64_t),
    ('pict_type', c_int),
    ('repeat_pict', c_int),
    ('pts', int64_t),
    ('dts', int64_t),
    ('last_pts', int64_t),
    ('last_dts', int64_t),
    ('fetch_timestamp', c_int),
    ('cur_frame_start_index', c_int),
    ('cur_frame_offset', int64_t * 4),
    ('cur_frame_pts', int64_t * 4),
    ('cur_frame_dts', int64_t * 4),
    ('flags', c_int),
    ('offset', int64_t),
    ('cur_frame_end', int64_t * 4),
    ('key_frame', c_int),
    ('convergence_duration', int64_t),
    ('dts_sync_point', c_int),
    ('dts_ref_dts_delta', c_int),
    ('pts_dts_delta', c_int),
    ('cur_frame_pos', int64_t * 4),
    ('pos', int64_t),
    ('last_pos', int64_t),
]
AVMetadataConv._fields_ = [
]
AVFormatParameters._fields_ = [
    ('time_base', AVRational),
    ('sample_rate', c_int),
    ('channels', c_int),
    ('width', c_int),
    ('height', c_int),
    ('pix_fmt', PixelFormat),
    ('channel', c_int),
    ('standard', STRING),
    ('mpeg2ts_raw', c_uint, 1),
    ('mpeg2ts_compute_pcr', c_uint, 1),
    ('initial_pause', c_uint, 1),
    ('prealloced_context', c_uint, 1),
    ('video_codec_id', CodecID),
    ('audio_codec_id', CodecID),
]
AVIndexEntry._fields_ = [
    ('pos', int64_t),
    ('timestamp', int64_t),
    ('flags', c_int, 2),
    ('size', c_int, 30),
    ('min_distance', c_int),
]
N8AVStream4DOLLAR_21E._pack_ = 4
N8AVStream4DOLLAR_21E._fields_ = [
    ('last_dts', int64_t),
    ('duration_gcd', int64_t),
    ('duration_count', c_int),
    ('duration_error', c_double * 725),
    ('codec_info_duration', int64_t),
]
RcOverride._fields_ = [
    ('start_frame', c_int),
    ('end_frame', c_int),
    ('qscale', c_int),
    ('quality_factor', c_float),
]

# values for enumeration 'AVPictureType'
AV_PICTURE_TYPE_NONE = 0
AV_PICTURE_TYPE_I = 1
AV_PICTURE_TYPE_P = 2
AV_PICTURE_TYPE_B = 3
AV_PICTURE_TYPE_S = 4
AV_PICTURE_TYPE_SI = 5
AV_PICTURE_TYPE_SP = 6
AV_PICTURE_TYPE_BI = 7
AVPictureType = c_int # enum
int8_t = c_int8
int16_t = c_int16
uint32_t = c_uint32
class AVPanScan(Structure):
    pass
AVFrame._pack_ = 4
AVFrame._fields_ = [
    ('data', POINTER(uint8_t) * 4),
    ('linesize', c_int * 4),
    ('base', POINTER(uint8_t) * 4),
    ('key_frame', c_int),
    ('pict_type', AVPictureType),
    ('pts', int64_t),
    ('coded_picture_number', c_int),
    ('display_picture_number', c_int),
    ('quality', c_int),
    ('age', c_int),
    ('reference', c_int),
    ('qscale_table', POINTER(int8_t)),
    ('qstride', c_int),
    ('mbskip_table', POINTER(uint8_t)),
    ('motion_val', POINTER(int16_t * 2) * 2),
    ('mb_type', POINTER(uint32_t)),
    ('motion_subsample_log2', uint8_t),
    ('opaque', c_void_p),
    ('error', uint64_t * 4),
    ('type', c_int),
    ('repeat_pict', c_int),
    ('qscale_type', c_int),
    ('interlaced_frame', c_int),
    ('top_field_first', c_int),
    ('pan_scan', POINTER(AVPanScan)),
    ('palette_has_changed', c_int),
    ('buffer_hints', c_int),
    ('dct_coeff', POINTER(c_short)),
    ('ref_index', POINTER(int8_t) * 2),
    ('reordered_opaque', int64_t),
    ('hwaccel_picture_private', c_void_p),
    ('pkt_pts', int64_t),
    ('pkt_dts', int64_t),
    ('owner', POINTER(AVCodecContext)),
    ('thread_opaque', c_void_p),
    ('best_effort_timestamp', int64_t),
    ('pkt_pos', int64_t),
    ('sample_aspect_ratio', AVRational),
    ('width', c_int),
    ('height', c_int),
    ('format', c_int),
]
class AVProfile(Structure):
    pass
AVCodec._fields_ = [
    ('name', STRING),
    ('type', AVMediaType),
    ('id', CodecID),
    ('priv_data_size', c_int),
    ('init', CFUNCTYPE(c_int, POINTER(AVCodecContext))),
    ('encode', CFUNCTYPE(c_int, POINTER(AVCodecContext), POINTER(uint8_t), c_int, c_void_p)),
    ('close', CFUNCTYPE(c_int, POINTER(AVCodecContext))),
    ('decode', CFUNCTYPE(c_int, POINTER(AVCodecContext), c_void_p, POINTER(c_int), POINTER(AVPacket))),
    ('capabilities', c_int),
    ('next', POINTER(AVCodec)),
    ('flush', CFUNCTYPE(None, POINTER(AVCodecContext))),
    ('supported_framerates', POINTER(AVRational)),
    ('pix_fmts', POINTER(PixelFormat)),
    ('long_name', STRING),
    ('supported_samplerates', POINTER(c_int)),
    ('sample_fmts', POINTER(AVSampleFormat)),
    ('channel_layouts', POINTER(int64_t)),
    ('max_lowres', uint8_t),
    ('priv_class', POINTER(AVClass)),
    ('profiles', POINTER(AVProfile)),
    ('init_thread_copy', CFUNCTYPE(c_int, POINTER(AVCodecContext))),
    ('update_thread_context', CFUNCTYPE(c_int, POINTER(AVCodecContext), POINTER(AVCodecContext))),
]
AVHWAccel._fields_ = [
    ('name', STRING),
    ('type', AVMediaType),
    ('id', CodecID),
    ('pix_fmt', PixelFormat),
    ('capabilities', c_int),
    ('next', POINTER(AVHWAccel)),
    ('start_frame', CFUNCTYPE(c_int, POINTER(AVCodecContext), POINTER(uint8_t), uint32_t)),
    ('decode_slice', CFUNCTYPE(c_int, POINTER(AVCodecContext), POINTER(uint8_t), uint32_t)),
    ('end_frame', CFUNCTYPE(c_int, POINTER(AVCodecContext))),
    ('priv_data_size', c_int),
]
AVPaletteControl._fields_ = [
    ('palette_changed', c_int),
    ('palette', c_uint * 256),
]
AVCodecParser._fields_ = [
    ('codec_ids', c_int * 5),
    ('priv_data_size', c_int),
    ('parser_init', CFUNCTYPE(c_int, POINTER(AVCodecParserContext))),
    ('parser_parse', CFUNCTYPE(c_int, POINTER(AVCodecParserContext), POINTER(AVCodecContext), POINTER(POINTER(uint8_t)), POINTER(c_int), POINTER(uint8_t), c_int)),
    ('parser_close', CFUNCTYPE(None, POINTER(AVCodecParserContext))),
    ('split', CFUNCTYPE(c_int, POINTER(AVCodecContext), POINTER(uint8_t), c_int)),
    ('next', POINTER(AVCodecParser)),
]
AVPanScan._fields_ = [
    ('id', c_int),
    ('width', c_int),
    ('height', c_int),
    ('position', int16_t * 2 * 3),
]
AVProfile._fields_ = [
    ('profile', c_int),
    ('name', STRING),
]
__all__ = ['CODEC_ID_ADPCM_IMA_AMV', 'CODEC_ID_SMC', 'CODEC_ID_AAC',
           'PIX_FMT_YUV422P', 'CODEC_ID_ADPCM_IMA_ISS',
           'AVMEDIA_TYPE_VIDEO', 'AVOutputFormat', 'uint8_t',
           'CODEC_ID_ANM', 'CODEC_ID_SONIC', 'CODEC_ID_SVQ1',
           'CODEC_ID_TRUEHD', 'CODEC_ID_SVQ3', 'CODEC_ID_TSCC',
           'CODEC_ID_BINKVIDEO', 'CODEC_ID_ASV2', 'CODEC_ID_SPEEX',
           'CODEC_ID_ASV1', 'CODEC_ID_ADPCM_G726', 'CODEC_ID_THP',
           'CODEC_ID_KGV1', 'CODEC_ID_VIXL', 'CODEC_ID_TIFF',
           'CODEC_ID_H264', 'CODEC_ID_H261', 'CODEC_ID_MIMIC',
           'CODEC_ID_H263', 'AVPanScan',
           'AV_AUDIO_SERVICE_TYPE_VOICE_OVER', 'AVCodec',
           'CODEC_ID_ALAC', 'AVCOL_PRI_FILM', 'CODEC_ID_PCM_ALAW',
           'PIX_FMT_VDPAU_H264', 'CODEC_ID_ADPCM_IMA_WAV',
           'AVCOL_SPC_NB', 'CODEC_ID_SMACKAUDIO',
           'AVCOL_TRC_UNSPECIFIED', 'AVDISCARD_BIDIR',
           'CODEC_ID_LOCO', 'AVPicture', 'CODEC_ID_INTERPLAY_VIDEO',
           'AV_SAMPLE_FMT_NB', 'PIX_FMT_YUVJ440P', 'CODEC_ID_FLAC',
           'CODEC_ID_PCM_U24LE', 'FF_OPT_TYPE_BINARY',
           'CODEC_ID_MPEG2TS', 'CODEC_ID_PCM_DVD', 'CODEC_ID_THEORA',
           'CODEC_ID_PBM', 'CODEC_ID_FFMETADATA', 'AVDictionary',
           'CODEC_ID_BMP', 'AV_SAMPLE_FMT_S32', 'AVDISCARD_NONE',
           'PIX_FMT_BGR555BE', 'PIX_FMT_DXVA2_VLD',
           'CODEC_ID_BETHSOFTVID', 'CODEC_ID_CAVS',
           'PIX_FMT_YUV444P10BE', 'CODEC_ID_AASC',
           'AV_PICTURE_TYPE_SP', 'CODEC_ID_AURA',
           'AV_AUDIO_SERVICE_TYPE_EFFECTS',
           'AV_AUDIO_SERVICE_TYPE_NB', 'CODEC_ID_TGV', 'CODEC_ID_TGQ',
           'PIX_FMT_YUVJ444P', 'CODEC_ID_WMAPRO', 'CODEC_ID_NONE',
           'CODEC_ID_MICRODVD', 'AVMEDIA_TYPE_DATA',
           'CODEC_ID_PCM_S16BE', 'CODEC_ID_PCM_F32BE', 'AVLPCType',
           'CODEC_ID_ADPCM_IMA_EA_SEAD', 'AV_SAMPLE_FMT_FLT',
           'CODEC_ID_QDRAW', 'AVInputFormat', 'CODEC_ID_WAVPACK',
           'CODEC_ID_HDMV_PGS_SUBTITLE', 'PIX_FMT_BGR24',
           'PIX_FMT_VDPAU_WMV3', 'AVHWAccel', 'CODEC_ID_WMV2',
           'CODEC_ID_WMV3', 'PIX_FMT_YUVJ420P', 'CODEC_ID_WMV1',
           'CODEC_ID_FLIC', 'CODEC_ID_PAM', 'CODEC_ID_AVS',
           'PIX_FMT_VAAPI_MOCO', 'CODEC_ID_PCM_S32LE',
           'CODEC_ID_PCM_S16LE', 'PIX_FMT_RGB48LE',
           'CODEC_ID_VOXWARE', 'AVCodecContext', 'PIX_FMT_RGB444LE',
           'CODEC_ID_ADPCM_IMA_WS', 'CODEC_ID_ADPCM_G722',
           'CODEC_ID_MP3ON4', 'CODEC_ID_PGMYUV', 'AVMediaType',
           'CODEC_ID_DVB_SUBTITLE', 'CODEC_ID_PCM_S24LE',
           'AV_AUDIO_SERVICE_TYPE_KARAOKE', 'CODEC_ID_NELLYMOSER',
           'AVColorPrimaries', 'AVCOL_RANGE_UNSPECIFIED',
           'PIX_FMT_BGR4_BYTE', 'PIX_FMT_UYVY422', 'PIX_FMT_BGRA',
           'AVCOL_PRI_BT470BG', 'AVDISCARD_NONREF', 'CODEC_ID_XVID',
           'AVCHROMA_LOC_TOP', 'CODEC_ID_AMR_NB', 'CODEC_ID_VORBIS',
           'AV_SAMPLE_FMT_DBL', 'AVPaletteControl',
           'CODEC_ID_CDGRAPHICS', 'CODEC_ID_AMR_WB',
           'CODEC_ID_PCM_F64BE', 'PIX_FMT_MONOWHITE',
           'CODEC_ID_FFH264', 'AVOption', 'CODEC_ID_RV40',
           'AV_PICTURE_TYPE_BI', 'CODEC_ID_TRUEMOTION1',
           'CODEC_ID_TRUEMOTION2', 'AVMEDIA_TYPE_ATTACHMENT',
           'CODEC_ID_RL2', 'CODEC_ID_MP4ALS', 'PIX_FMT_BGR8',
           'PIX_FMT_VDPAU_MPEG1', 'PIX_FMT_BGR4',
           'PIX_FMT_VDPAU_MPEG2', 'CODEC_ID_SONIC_LS',
           'PIX_FMT_VDPAU_MPEG4', 'CODEC_ID_VMDAUDIO',
           'PIX_FMT_RGB444BE', 'CODEC_ID_FLV1', 'CODEC_ID_FRAPS',
           'CODEC_ID_INDEO5', 'CODEC_ID_INDEO4', 'CODEC_ID_INDEO3',
           'CODEC_ID_INDEO2', 'FF_OPT_TYPE_DOUBLE', 'AVDISCARD_ALL',
           'CODEC_ID_ADPCM_EA_XAS', 'AV_LPC_TYPE_FIXED',
           'CODEC_ID_PCM_MULAW', 'CODEC_ID_APE', 'PIX_FMT_YUYV422',
           'CODEC_ID_RV30', 'AVSTREAM_PARSE_NONE',
           'PIX_FMT_UYYVYY411', 'CODEC_ID_FFV1', 'CODEC_ID_TTA',
           'CODEC_ID_TTF', 'AVDISCARD_NONKEY', 'CODEC_ID_MPEG4',
           'CODEC_ID_8SVX_FIB', 'CODEC_ID_IFF_ILBM',
           'AVCOL_RANGE_MPEG', 'CODEC_ID_DIRAC', 'CodecID',
           'CODEC_ID_EAC3', 'AVFormatContext', 'CODEC_ID_PROBE',
           'CODEC_ID_PCM_U24BE', 'CODEC_ID_PCM_S24BE',
           'AV_PICTURE_TYPE_NONE', 'AVCHROMA_LOC_BOTTOM', 'uint16_t',
           'CODEC_ID_RPZA', 'PIX_FMT_YUV420P9LE', 'CODEC_ID_FLASHSV2',
           'CODEC_ID_NUV', 'CODEC_ID_SNOW', 'CODEC_ID_VP8',
           'CODEC_ID_PCM_S8', 'CODEC_ID_SMACKVIDEO', 'CODEC_ID_MSZH',
           'PIX_FMT_RGB555BE', 'FF_OPT_TYPE_INT', 'CODEC_ID_S302M',
           'CODEC_ID_4XM', 'CODEC_ID_SGI', 'CODEC_ID_ESCAPE124',
           'PIX_FMT_YUV420P10BE', 'CODEC_ID_A64_MULTI',
           'PIX_FMT_NONE', 'CODEC_ID_ROQ', 'CODEC_ID_ATRAC3P',
           'CODEC_ID_PCM_F32LE', 'CODEC_ID_TARGA', 'PIX_FMT_ARGB',
           'CODEC_ID_MPEG2VIDEO', 'CODEC_ID_ADPCM_SWF',
           'CODEC_ID_WESTWOOD_SND1', 'CODEC_ID_8SVX_EXP',
           'CODEC_ID_DPX', 'AVCOL_RANGE_NB', 'CODEC_ID_IDCIN',
           'CODEC_ID_SUNRAST', 'CODEC_ID_VCR1', 'CODEC_ID_XAN_WC3',
           'AVPictureType', 'AV_AUDIO_SERVICE_TYPE_DIALOGUE',
           'CODEC_ID_ADPCM_IMA_EA_EACS', 'uint64_t',
           'PIX_FMT_YUVA420P', 'CODEC_ID_SP5X', 'CODEC_ID_ATRAC3',
           'CODEC_ID_ATRAC1', 'AV_PICTURE_TYPE_P',
           'AVCHROMA_LOC_LEFT', 'AVMEDIA_TYPE_AUDIO', 'RcOverride',
           'CODEC_ID_TXD', 'PIX_FMT_YUV420P', 'CODEC_ID_MSMPEG4V2',
           'CODEC_ID_MSMPEG4V3', 'CODEC_ID_MSMPEG4V1',
           'CODEC_ID_SHORTEN', 'CODEC_ID_H263I',
           'AVSTREAM_PARSE_TIMESTAMPS', 'CODEC_ID_RAWVIDEO',
           'CODEC_ID_V210X', 'CODEC_ID_H263P', 'CODEC_ID_RA_288',
           'AV_SAMPLE_FMT_U8', 'AVCOL_TRC_BT709', 'AVCOL_TRC_GAMMA28',
           'CODEC_ID_BINKAUDIO_RDFT', 'AVSTREAM_PARSE_HEADERS',
           'CODEC_ID_MJPEG', 'CODEC_ID_CLJR', 'AVCOL_SPC_RGB',
           'CODEC_ID_C93', 'CODEC_ID_RA_144', 'AVCOL_SPC_SMPTE170M',
           'AVOptionType', 'AV_PICTURE_TYPE_B', 'AVCHROMA_LOC_NB',
           'PIX_FMT_YUV422P10BE', 'PIX_FMT_YUV410P',
           'CODEC_ID_GSM_MS', 'CODEC_ID_PCM_LXF', 'CODEC_ID_VB',
           'AV_AUDIO_SERVICE_TYPE_EMERGENCY', 'AV_PICTURE_TYPE_I',
           'AVSTREAM_PARSE_FULL', 'CODEC_ID_ADPCM_XA', 'CODEC_ID_VC1',
           'CODEC_ID_BFI', 'AVCHROMA_LOC_CENTER',
           'PIX_FMT_YUV420P9BE', 'PIX_FMT_RGB8', 'CODEC_ID_DFA',
           'PIX_FMT_RGB4', 'PIX_FMT_YUV422P16BE', 'CODEC_ID_DTS',
           'FF_OPT_TYPE_INT64', 'CODEC_ID_ADPCM_4XM', 'PIX_FMT_RGBA',
           'PIX_FMT_BGR565LE', 'AVCOL_PRI_BT709', 'AVStream',
           'FF_OPT_TYPE_RATIONAL', 'PIX_FMT_BGR555LE',
           'AV_AUDIO_SERVICE_TYPE_HEARING_IMPAIRED', 'CODEC_ID_ULTI',
           'PIX_FMT_YUV444P16BE', 'AVMEDIA_TYPE_SUBTITLE',
           'CODEC_ID_TWINVQ', 'CODEC_ID_CINEPAK', 'CODEC_ID_FLASHSV',
           'int16_t', 'PIX_FMT_YUV422P16LE', 'AVProbeData',
           'CODEC_ID_TRUESPEECH', 'CODEC_ID_ADPCM_EA', 'PIX_FMT_NV21',
           'AVClass', 'CODEC_ID_RV20', 'PIX_FMT_YUV422P10LE',
           'CODEC_ID_GIF', 'CODEC_ID_FRWU', 'CODEC_ID_MLP',
           'PIX_FMT_VDPAU_VC1', 'AVProfile', 'PIX_FMT_GRAY8',
           'PIX_FMT_YUV444P', 'CODEC_ID_PCM_BLURAY', 'CODEC_ID_YOP',
           'CODEC_ID_A64_MULTI5', 'AVCOL_PRI_SMPTE240M',
           'PIX_FMT_YUV420P16BE', 'AVIOContext', 'CODEC_ID_SRT',
           'CODEC_ID_PGM', 'CODEC_ID_DSICINVIDEO',
           'AVCHROMA_LOC_TOPLEFT', 'CODEC_ID_DSICINAUDIO',
           'CODEC_ID_AC3', 'CODEC_ID_HUFFYUV', 'PIX_FMT_RGB48BE',
           'CODEC_ID_8BPS', 'CODEC_ID_PNG', 'AVCOL_SPC_UNSPECIFIED',
           'AV_SAMPLE_FMT_NONE', 'PIX_FMT_BGR565BE', 'CODEC_ID_QDMC',
           'CODEC_ID_JPEG2000', 'PIX_FMT_GRAY8A',
           'CODEC_ID_ADPCM_EA_MAXIS_XA', 'AV_PICTURE_TYPE_SI',
           'AVRational', 'AVColorSpace', 'CODEC_ID_ADPCM_IMA_SMJPEG',
           'CODEC_ID_AAC_LATM', 'CODEC_ID_DXA', 'CODEC_ID_V210',
           'CODEC_ID_PCM_F64LE', 'AVStreamParseType', 'CODEC_ID_VP6A',
           'CODEC_ID_ADPCM_IMA_DK4', 'CODEC_ID_ADPCM_IMA_DK3',
           'CODEC_ID_VP6F', 'FF_OPT_TYPE_CONST', 'CODEC_ID_MAD',
           'AVCodecTag', 'PIX_FMT_YUVJ422P', 'CODEC_ID_VMDVIDEO',
           'AVMetadataConv', 'PIX_FMT_YUV411P', 'AVIndexEntry',
           'AV_LPC_TYPE_CHOLESKY', 'CODEC_ID_CYUV', 'CODEC_ID_ZLIB',
           'PIX_FMT_BGR48LE', 'CODEC_ID_MP3ADU',
           'CODEC_ID_PCM_S16LE_PLANAR', 'PIX_FMT_ABGR',
           'PIX_FMT_PAL8', 'CODEC_ID_DVD_SUBTITLE',
           'CODEC_ID_XAN_WC4', 'CODEC_ID_ADPCM_EA_R1',
           'CODEC_ID_ADPCM_EA_R2', 'CODEC_ID_ADPCM_EA_R3',
           'CODEC_ID_MACE3', 'PIX_FMT_YUV420P10LE', 'CODEC_ID_MACE6',
           'AVCOL_PRI_NB', 'uint32_t', 'CODEC_ID_JPEGLS', 'AVProgram',
           'CODEC_ID_8SVX_RAW', 'PIX_FMT_MONOBLACK', 'AVCOL_SPC_FCC',
           'CODEC_ID_IFF_BYTERUN1', 'CODEC_ID_PCX',
           'FF_OPT_TYPE_STRING', 'CODEC_ID_PRORES', 'CODEC_ID_AURA2',
           'AVFrame', 'CODEC_ID_PCM_S32BE', 'PIX_FMT_VAAPI_IDCT',
           'AVAudioServiceType', 'CODEC_ID_MJPEGB', 'CODEC_ID_MSRLE',
           'CODEC_ID_PCM_S24DAUD', 'CODEC_ID_DVB_TELETEXT',
           'CODEC_ID_WMALOSSLESS', 'CODEC_ID_ADPCM_YAMAHA',
           'PIX_FMT_NV12', 'PIX_FMT_NB', 'CODEC_ID_ANSI',
           'CODEC_ID_CSCD', 'CODEC_ID_ADPCM_IMA_QT',
           'CODEC_ID_PICTOR', 'AVCOL_TRC_GAMMA22', 'PIX_FMT_YUV440P',
           'CODEC_ID_XSUB', 'AVCOL_PRI_UNSPECIFIED',
           'PIX_FMT_RGB565LE', 'AV_LPC_TYPE_NONE', 'AVPacketList',
           'CODEC_ID_CMV', 'CODEC_ID_PCM_U8', 'CODEC_ID_QTRLE',
           'PIX_FMT_XVMC_MPEG2_MC', 'PIX_FMT_RGB555LE',
           'AVCOL_SPC_SMPTE240M', 'CODEC_ID_MOTIONPIXELS',
           'AVCHROMA_LOC_BOTTOMLEFT', 'PIX_FMT_YUV444P9BE',
           'AVCOL_TRC_NB', 'AVDISCARD_DEFAULT', 'CODEC_ID_ADPCM_THP',
           'AVCodecParser', 'AVCodecParserContext', 'CODEC_ID_TMV',
           'CODEC_ID_PTX', 'CODEC_ID_MMVIDEO',
           'CODEC_ID_ADPCM_SBPRO_2', 'CODEC_ID_DVAUDIO',
           'CODEC_ID_DVVIDEO', 'CODEC_ID_WMAVOICE', 'CODEC_ID_MXPEG',
           'AVFrac', 'CODEC_ID_ROQ_DPCM', 'CODEC_ID_MOV_TEXT',
           'PIX_FMT_GRAY16BE', 'N8AVOption4DOLLAR_19E',
           'CODEC_ID_TIERTEXSEQVIDEO', 'PIX_FMT_BGR444BE',
           'PIX_FMT_YUV444P9LE', 'AV_LPC_TYPE_DEFAULT',
           'CODEC_ID_FFVHUFF', 'PIX_FMT_YUV444P16LE', 'CODEC_ID_VMNC',
           'AVCOL_RANGE_JPEG', 'AV_LPC_TYPE_LEVINSON',
           'CODEC_ID_ZMBV', 'CODEC_ID_RV10', 'CODEC_ID_QCELP',
           'PIX_FMT_YUV420P16LE', 'AVSampleFormat', 'AV_LPC_TYPE_NB',
           'AVCOL_SPC_BT709', 'CODEC_ID_XAN_DPCM',
           'CODEC_ID_MPEG1VIDEO', 'AVCOL_SPC_BT470BG',
           'CODEC_ID_INTERPLAY_DPCM', 'CODEC_ID_GSM',
           'CODEC_ID_BINKAUDIO_DCT', 'CODEC_ID_MSVIDEO1',
           'PIX_FMT_XVMC_MPEG2_IDCT', 'CODEC_ID_ADPCM_CT',
           'AVFormatParameters', 'CODEC_ID_VP5', 'CODEC_ID_VP6',
           'CODEC_ID_ADPCM_ADX', 'CODEC_ID_PCM_U32BE', 'CODEC_ID_VP3',
           'PIX_FMT_BGR48BE', 'CODEC_ID_WNV1', 'PIX_FMT_RGB565BE',
           'CODEC_ID_SIPR', 'CODEC_ID_R10K', 'AVPacket',
           'AV_PICTURE_TYPE_S', 'FF_OPT_TYPE_FLOAT',
           'CODEC_ID_MUSEPACK8', 'CODEC_ID_MDEC', 'CODEC_ID_ADPCM_MS',
           'CODEC_ID_MUSEPACK7', 'int8_t', 'AVChapter',
           'CODEC_ID_IMC', 'CODEC_ID_LJPEG', 'CODEC_ID_CELT',
           'AVCOL_PRI_SMPTE170M', 'CODEC_ID_LAGARITH',
           'AV_SAMPLE_FMT_S16',
           'AV_AUDIO_SERVICE_TYPE_VISUALLY_IMPAIRED', 'SwsContext',
           'PIX_FMT_RGB24', 'CODEC_ID_MPEG2VIDEO_XVMC',
           'CODEC_ID_COOK', 'CODEC_ID_AMV', 'CODEC_ID_SOL_DPCM',
           'CODEC_ID_ADPCM_SBPRO_3', 'CODEC_ID_ADPCM_SBPRO_4',
           'CODEC_ID_PCM_U16LE', 'AVColorRange',
           'AVSTREAM_PARSE_FULL_ONCE', 'CODEC_ID_PCM_U32LE',
           'CODEC_ID_WMAV2', 'AVColorTransferCharacteristic',
           'CODEC_ID_WMAV1', 'AV_AUDIO_SERVICE_TYPE_MAIN', 'int64_t',
           'CODEC_ID_PPM', 'CODEC_ID_R210', 'PIX_FMT_BGR444LE',
           'PixelFormat', 'CODEC_ID_KMVC', 'CODEC_ID_PCM_ZORK',
           'FF_OPT_TYPE_FLAGS', 'AVMEDIA_TYPE_NB', 'CODEC_ID_TQI',
           'AVDiscard', 'PIX_FMT_RGB4_BYTE', 'PIX_FMT_GRAY16LE',
           'CODEC_ID_PCM_U16BE', 'CODEC_ID_SSA',
           'AVMEDIA_TYPE_UNKNOWN', 'CODEC_ID_QPEG',
           'AVCOL_PRI_BT470M', 'PIX_FMT_YUV444P10LE',
           'CODEC_ID_WS_VQA', 'N8AVStream4DOLLAR_21E',
           'PIX_FMT_VAAPI_VLD', 'AVChromaLocation', 'CODEC_ID_JV',
           'CODEC_ID_DNXHD', 'AV_AUDIO_SERVICE_TYPE_COMMENTARY',
           'CODEC_ID_TEXT', 'AVCHROMA_LOC_UNSPECIFIED',
           'CODEC_ID_QDM2', 'CODEC_ID_MP2', 'CODEC_ID_MP3',
           'CODEC_ID_MP1']
