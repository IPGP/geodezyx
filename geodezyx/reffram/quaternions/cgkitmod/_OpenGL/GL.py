# Shadow module for module "OpenGL.GL"

try:

    # Try to import the original module...
    from OpenGL.GL import *
    
except:

    # Create dummy symbols...

    GL_2D = 1536
    GL_2_BYTES = 5127
    GL_3D = 1537
    GL_3D_COLOR = 1538
    GL_3D_COLOR_TEXTURE = 1539
    GL_3_BYTES = 5128
    GL_4D_COLOR_TEXTURE = 1540
    GL_4_BYTES = 5129
    GL_ACCUM = 256
    GL_ACCUM_ALPHA_BITS = 3419
    GL_ACCUM_BLUE_BITS = 3418
    GL_ACCUM_BUFFER_BIT = 512
    GL_ACCUM_CLEAR_VALUE = 2944
    GL_ACCUM_GREEN_BITS = 3417
    GL_ACCUM_RED_BITS = 3416
    GL_ADD = 260
    GL_ALL_ATTRIB_BITS = 1048575
    GL_ALPHA = 6406
    GL_ALPHA12 = 32829
    GL_ALPHA16 = 32830
    GL_ALPHA4 = 32827
    GL_ALPHA8 = 32828
    GL_ALPHA_BIAS = 3357
    GL_ALPHA_BITS = 3413
    GL_ALPHA_SCALE = 3356
    GL_ALPHA_TEST = 3008
    GL_ALPHA_TEST_FUNC = 3009
    GL_ALPHA_TEST_REF = 3010
    GL_ALWAYS = 519
    GL_AMBIENT = 4608
    GL_AMBIENT_AND_DIFFUSE = 5634
    GL_AND = 5377
    GL_AND_INVERTED = 5380
    GL_AND_REVERSE = 5378
    GL_ATTRIB_STACK_DEPTH = 2992
    GL_AUTO_NORMAL = 3456
    GL_AUX0 = 1033
    GL_AUX1 = 1034
    GL_AUX2 = 1035
    GL_AUX3 = 1036
    GL_AUX_BUFFERS = 3072
    GL_BACK = 1029
    GL_BACK_LEFT = 1026
    GL_BACK_RIGHT = 1027
    GL_BITMAP = 6656
    GL_BITMAP_TOKEN = 1796
    GL_BLEND = 3042
    GL_BLEND_DST = 3040
    GL_BLEND_SRC = 3041
    GL_BLUE = 6405
    GL_BLUE_BIAS = 3355
    GL_BLUE_BITS = 3412
    GL_BLUE_SCALE = 3354
    GL_BYTE = 5120
    GL_C3F_V3F = 10788
    GL_C4F_N3F_V3F = 10790
    GL_C4UB_V2F = 10786
    GL_C4UB_V3F = 10787
    GL_CCW = 2305
    GL_CLAMP = 10496
    GL_CLEAR = 5376
    GL_CLIENT_ALL_ATTRIB_BITS = -1
    GL_CLIENT_ATTRIB_STACK_DEPTH = 2993
    GL_CLIENT_PIXEL_STORE_BIT = 1
    GL_CLIENT_VERTEX_ARRAY_BIT = 2
    GL_CLIP_PLANE0 = 12288
    GL_CLIP_PLANE1 = 12289
    GL_CLIP_PLANE2 = 12290
    GL_CLIP_PLANE3 = 12291
    GL_CLIP_PLANE4 = 12292
    GL_CLIP_PLANE5 = 12293
    GL_COEFF = 2560
    GL_COLOR = 6144
    GL_COLOR_ARRAY = 32886
    GL_COLOR_ARRAY_POINTER = 32912
    GL_COLOR_ARRAY_SIZE = 32897
    GL_COLOR_ARRAY_STRIDE = 32899
    GL_COLOR_ARRAY_TYPE = 32898
    GL_COLOR_BUFFER_BIT = 16384
    GL_COLOR_CLEAR_VALUE = 3106
    GL_COLOR_INDEX = 6400
    GL_COLOR_INDEXES = 5635
    GL_COLOR_LOGIC_OP = 3058
    GL_COLOR_MATERIAL = 2903
    GL_COLOR_MATERIAL_FACE = 2901
    GL_COLOR_MATERIAL_PARAMETER = 2902
    GL_COLOR_WRITEMASK = 3107
    GL_COMPILE = 4864
    GL_COMPILE_AND_EXECUTE = 4865
    GL_CONSTANT_ATTENUATION = 4615
    GL_COPY = 5379
    GL_COPY_INVERTED = 5388
    GL_COPY_PIXEL_TOKEN = 1798
    GL_CULL_FACE = 2884
    GL_CULL_FACE_MODE = 2885
    GL_CURRENT_BIT = 1
    GL_CURRENT_COLOR = 2816
    GL_CURRENT_INDEX = 2817
    GL_CURRENT_NORMAL = 2818
    GL_CURRENT_RASTER_COLOR = 2820
    GL_CURRENT_RASTER_DISTANCE = 2825
    GL_CURRENT_RASTER_INDEX = 2821
    GL_CURRENT_RASTER_POSITION = 2823
    GL_CURRENT_RASTER_POSITION_VALID = 2824
    GL_CURRENT_RASTER_TEXTURE_COORDS = 2822
    GL_CURRENT_TEXTURE_COORDS = 2819
    GL_CW = 2304
    GL_DECAL = 8449
    GL_DECR = 7683
    GL_DEPTH = 6145
    GL_DEPTH_BIAS = 3359
    GL_DEPTH_BITS = 3414
    GL_DEPTH_BUFFER_BIT = 256
    GL_DEPTH_CLEAR_VALUE = 2931
    GL_DEPTH_COMPONENT = 6402
    GL_DEPTH_FUNC = 2932
    GL_DEPTH_RANGE = 2928
    GL_DEPTH_SCALE = 3358
    GL_DEPTH_TEST = 2929
    GL_DEPTH_WRITEMASK = 2930
    GL_DIFFUSE = 4609
    GL_DITHER = 3024
    GL_DOMAIN = 2562
    GL_DONT_CARE = 4352
    GL_DOUBLE = 5130
    GL_DOUBLEBUFFER = 3122
    GL_DRAW_BUFFER = 3073
    GL_DRAW_PIXEL_TOKEN = 1797
    GL_DST_ALPHA = 772
    GL_DST_COLOR = 774
    GL_EDGE_FLAG = 2883
    GL_EDGE_FLAG_ARRAY = 32889
    GL_EDGE_FLAG_ARRAY_POINTER = 32915
    GL_EDGE_FLAG_ARRAY_STRIDE = 32908
    GL_EMISSION = 5632
    GL_ENABLE_BIT = 8192
    GL_EQUAL = 514
    GL_EQUIV = 5385
    GL_EVAL_BIT = 65536
    GL_EXP = 2048
    GL_EXP2 = 2049
    GL_EXTENSIONS = 7939
    GL_EYE_LINEAR = 9216
    GL_EYE_PLANE = 9474
    GL_FALSE = 0
    GL_FASTEST = 4353
    GL_FEEDBACK = 7169
    GL_FEEDBACK_BUFFER_POINTER = 3568
    GL_FEEDBACK_BUFFER_SIZE = 3569
    GL_FEEDBACK_BUFFER_TYPE = 3570
    GL_FILL = 6914
    GL_FLAT = 7424
    GL_FLOAT = 5126
    GL_FOG = 2912
    GL_FOG_BIT = 128
    GL_FOG_COLOR = 2918
    GL_FOG_DENSITY = 2914
    GL_FOG_END = 2916
    GL_FOG_HINT = 3156
    GL_FOG_INDEX = 2913
    GL_FOG_MODE = 2917
    GL_FOG_START = 2915
    GL_FRONT = 1028
    GL_FRONT_AND_BACK = 1032
    GL_FRONT_FACE = 2886
    GL_FRONT_LEFT = 1024
    GL_FRONT_RIGHT = 1025
    GL_GEQUAL = 518
    GL_GREATER = 516
    GL_GREEN = 6404
    GL_GREEN_BIAS = 3353
    GL_GREEN_BITS = 3411
    GL_GREEN_SCALE = 3352
    GL_HINT_BIT = 32768
    GL_INCR = 7682
    GL_INDEX_ARRAY = 32887
    GL_INDEX_ARRAY_POINTER = 32913
    GL_INDEX_ARRAY_STRIDE = 32902
    GL_INDEX_ARRAY_TYPE = 32901
    GL_INDEX_BITS = 3409
    GL_INDEX_CLEAR_VALUE = 3104
    GL_INDEX_LOGIC_OP = 3057
    GL_INDEX_MODE = 3120
    GL_INDEX_OFFSET = 3347
    GL_INDEX_SHIFT = 3346
    GL_INDEX_WRITEMASK = 3105
    GL_INT = 5124
    GL_INTENSITY = 32841
    GL_INTENSITY12 = 32844
    GL_INTENSITY16 = 32845
    GL_INTENSITY4 = 32842
    GL_INTENSITY8 = 32843
    GL_INVALID_ENUM = 1280
    GL_INVALID_OPERATION = 1282
    GL_INVALID_VALUE = 1281
    GL_INVERT = 5386
    GL_KEEP = 7680
    GL_LEFT = 1030
    GL_LEQUAL = 515
    GL_LESS = 513
    GL_LIGHT0 = 16384
    GL_LIGHT1 = 16385
    GL_LIGHT2 = 16386
    GL_LIGHT3 = 16387
    GL_LIGHT4 = 16388
    GL_LIGHT5 = 16389
    GL_LIGHT6 = 16390
    GL_LIGHT7 = 16391
    GL_LIGHTING = 2896
    GL_LIGHTING_BIT = 64
    GL_LIGHT_MODEL_AMBIENT = 2899
    GL_LIGHT_MODEL_LOCAL_VIEWER = 2897
    GL_LIGHT_MODEL_TWO_SIDE = 2898
    GL_LINE = 6913
    GL_LINEAR = 9729
    GL_LINEAR_ATTENUATION = 4616
    GL_LINEAR_MIPMAP_LINEAR = 9987
    GL_LINEAR_MIPMAP_NEAREST = 9985
    GL_LINES = 1
    GL_LINE_BIT = 4
    GL_LINE_LOOP = 2
    GL_LINE_RESET_TOKEN = 1799
    GL_LINE_SMOOTH = 2848
    GL_LINE_SMOOTH_HINT = 3154
    GL_LINE_STIPPLE = 2852
    GL_LINE_STIPPLE_PATTERN = 2853
    GL_LINE_STIPPLE_REPEAT = 2854
    GL_LINE_STRIP = 3
    GL_LINE_TOKEN = 1794
    GL_LINE_WIDTH = 2849
    GL_LINE_WIDTH_GRANULARITY = 2851
    GL_LINE_WIDTH_RANGE = 2850
    GL_LIST_BASE = 2866
    GL_LIST_BIT = 131072
    GL_LIST_INDEX = 2867
    GL_LIST_MODE = 2864
    GL_LOAD = 257
    GL_LOGIC_OP = 3057
    GL_LOGIC_OP_MODE = 3056
    GL_LUMINANCE = 6409
    GL_LUMINANCE12 = 32833
    GL_LUMINANCE12_ALPHA12 = 32839
    GL_LUMINANCE12_ALPHA4 = 32838
    GL_LUMINANCE16 = 32834
    GL_LUMINANCE16_ALPHA16 = 32840
    GL_LUMINANCE4 = 32831
    GL_LUMINANCE4_ALPHA4 = 32835
    GL_LUMINANCE6_ALPHA2 = 32836
    GL_LUMINANCE8 = 32832
    GL_LUMINANCE8_ALPHA8 = 32837
    GL_LUMINANCE_ALPHA = 6410
    GL_MAP1_COLOR_4 = 3472
    GL_MAP1_GRID_DOMAIN = 3536
    GL_MAP1_GRID_SEGMENTS = 3537
    GL_MAP1_INDEX = 3473
    GL_MAP1_NORMAL = 3474
    GL_MAP1_TEXTURE_COORD_1 = 3475
    GL_MAP1_TEXTURE_COORD_2 = 3476
    GL_MAP1_TEXTURE_COORD_3 = 3477
    GL_MAP1_TEXTURE_COORD_4 = 3478
    GL_MAP1_VERTEX_3 = 3479
    GL_MAP1_VERTEX_4 = 3480
    GL_MAP2_COLOR_4 = 3504
    GL_MAP2_GRID_DOMAIN = 3538
    GL_MAP2_GRID_SEGMENTS = 3539
    GL_MAP2_INDEX = 3505
    GL_MAP2_NORMAL = 3506
    GL_MAP2_TEXTURE_COORD_1 = 3507
    GL_MAP2_TEXTURE_COORD_2 = 3508
    GL_MAP2_TEXTURE_COORD_3 = 3509
    GL_MAP2_TEXTURE_COORD_4 = 3510
    GL_MAP2_VERTEX_3 = 3511
    GL_MAP2_VERTEX_4 = 3512
    GL_MAP_COLOR = 3344
    GL_MAP_STENCIL = 3345
    GL_MATRIX_MODE = 2976
    GL_MAX_ATTRIB_STACK_DEPTH = 3381
    GL_MAX_CLIENT_ATTRIB_STACK_DEPTH = 3387
    GL_MAX_CLIP_PLANES = 3378
    GL_MAX_EVAL_ORDER = 3376
    GL_MAX_LIGHTS = 3377
    GL_MAX_LIST_NESTING = 2865
    GL_MAX_MODELVIEW_STACK_DEPTH = 3382
    GL_MAX_NAME_STACK_DEPTH = 3383
    GL_MAX_PIXEL_MAP_TABLE = 3380
    GL_MAX_PROJECTION_STACK_DEPTH = 3384
    GL_MAX_TEXTURE_SIZE = 3379
    GL_MAX_TEXTURE_STACK_DEPTH = 3385
    GL_MAX_VIEWPORT_DIMS = 3386
    GL_MODELVIEW = 5888
    GL_MODELVIEW_MATRIX = 2982
    GL_MODELVIEW_STACK_DEPTH = 2979
    GL_MODULATE = 8448
    GL_MULT = 259
    GL_N3F_V3F = 10789
    GL_NAME_STACK_DEPTH = 3440
    GL_NAND = 5390
    GL_NEAREST = 9728
    GL_NEAREST_MIPMAP_LINEAR = 9986
    GL_NEAREST_MIPMAP_NEAREST = 9984
    GL_NEVER = 512
    GL_NICEST = 4354
    GL_NONE = 0
    GL_NOOP = 5381
    GL_NOR = 5384
    GL_NORMALIZE = 2977
    GL_NORMAL_ARRAY = 32885
    GL_NORMAL_ARRAY_POINTER = 32911
    GL_NORMAL_ARRAY_STRIDE = 32895
    GL_NORMAL_ARRAY_TYPE = 32894
    GL_NOTEQUAL = 517
    GL_NO_ERROR = 0
    GL_OBJECT_LINEAR = 9217
    GL_OBJECT_PLANE = 9473
    GL_ONE = 1
    GL_ONE_MINUS_DST_ALPHA = 773
    GL_ONE_MINUS_DST_COLOR = 775
    GL_ONE_MINUS_SRC_ALPHA = 771
    GL_ONE_MINUS_SRC_COLOR = 769
    GL_OR = 5383
    GL_ORDER = 2561
    GL_OR_INVERTED = 5389
    GL_OR_REVERSE = 5387
    GL_OUT_OF_MEMORY = 1285
    GL_PACK_ALIGNMENT = 3333
    GL_PACK_LSB_FIRST = 3329
    GL_PACK_ROW_LENGTH = 3330
    GL_PACK_SKIP_PIXELS = 3332
    GL_PACK_SKIP_ROWS = 3331
    GL_PACK_SWAP_BYTES = 3328
    GL_PASS_THROUGH_TOKEN = 1792
    GL_PERSPECTIVE_CORRECTION_HINT = 3152
    GL_PIXEL_MAP_A_TO_A = 3193
    GL_PIXEL_MAP_A_TO_A_SIZE = 3257
    GL_PIXEL_MAP_B_TO_B = 3192
    GL_PIXEL_MAP_B_TO_B_SIZE = 3256
    GL_PIXEL_MAP_G_TO_G = 3191
    GL_PIXEL_MAP_G_TO_G_SIZE = 3255
    GL_PIXEL_MAP_I_TO_A = 3189
    GL_PIXEL_MAP_I_TO_A_SIZE = 3253
    GL_PIXEL_MAP_I_TO_B = 3188
    GL_PIXEL_MAP_I_TO_B_SIZE = 3252
    GL_PIXEL_MAP_I_TO_G = 3187
    GL_PIXEL_MAP_I_TO_G_SIZE = 3251
    GL_PIXEL_MAP_I_TO_I = 3184
    GL_PIXEL_MAP_I_TO_I_SIZE = 3248
    GL_PIXEL_MAP_I_TO_R = 3186
    GL_PIXEL_MAP_I_TO_R_SIZE = 3250
    GL_PIXEL_MAP_R_TO_R = 3190
    GL_PIXEL_MAP_R_TO_R_SIZE = 3254
    GL_PIXEL_MAP_S_TO_S = 3185
    GL_PIXEL_MAP_S_TO_S_SIZE = 3249
    GL_PIXEL_MODE_BIT = 32
    GL_POINT = 6912
    GL_POINTS = 0
    GL_POINT_BIT = 2
    GL_POINT_SIZE = 2833
    GL_POINT_SIZE_GRANULARITY = 2835
    GL_POINT_SIZE_RANGE = 2834
    GL_POINT_SMOOTH = 2832
    GL_POINT_SMOOTH_HINT = 3153
    GL_POINT_TOKEN = 1793
    GL_POLYGON = 9
    GL_POLYGON_BIT = 8
    GL_POLYGON_MODE = 2880
    GL_POLYGON_OFFSET_FACTOR = 32824
    GL_POLYGON_OFFSET_FILL = 32823
    GL_POLYGON_OFFSET_LINE = 10754
    GL_POLYGON_OFFSET_POINT = 10753
    GL_POLYGON_OFFSET_UNITS = 10752
    GL_POLYGON_SMOOTH = 2881
    GL_POLYGON_SMOOTH_HINT = 3155
    GL_POLYGON_STIPPLE = 2882
    GL_POLYGON_STIPPLE_BIT = 16
    GL_POLYGON_TOKEN = 1795
    GL_POSITION = 4611
    GL_PROJECTION = 5889
    GL_PROJECTION_MATRIX = 2983
    GL_PROJECTION_STACK_DEPTH = 2980
    GL_PROXY_TEXTURE_1D = 32867
    GL_PROXY_TEXTURE_2D = 32868
    GL_Q = 8195
    GL_QUADRATIC_ATTENUATION = 4617
    GL_QUADS = 7
    GL_QUAD_STRIP = 8
    GL_R = 8194
    GL_R3_G3_B2 = 10768
    GL_READ_BUFFER = 3074
    GL_RED = 6403
    GL_RED_BIAS = 3349
    GL_RED_BITS = 3410
    GL_RED_SCALE = 3348
    GL_RENDER = 7168
    GL_RENDERER = 7937
    GL_RENDER_MODE = 3136
    GL_REPEAT = 10497
    GL_REPLACE = 7681
    GL_RETURN = 258
    GL_RGB = 6407
    GL_RGB10 = 32850
    GL_RGB10_A2 = 32857
    GL_RGB12 = 32851
    GL_RGB16 = 32852
    GL_RGB4 = 32847
    GL_RGB5 = 32848
    GL_RGB5_A1 = 32855
    GL_RGB8 = 32849
    GL_RGBA = 6408
    GL_RGBA12 = 32858
    GL_RGBA16 = 32859
    GL_RGBA2 = 32853
    GL_RGBA4 = 32854
    GL_RGBA8 = 32856
    GL_RGBA_MODE = 3121
    GL_RIGHT = 1031
    GL_S = 8192
    GL_SCISSOR_BIT = 524288
    GL_SCISSOR_BOX = 3088
    GL_SCISSOR_TEST = 3089
    GL_SELECT = 7170
    GL_SELECTION_BUFFER_POINTER = 3571
    GL_SELECTION_BUFFER_SIZE = 3572
    GL_SET = 5391
    GL_SHADE_MODEL = 2900
    GL_SHININESS = 5633
    GL_SHORT = 5122
    GL_SMOOTH = 7425
    GL_SPECULAR = 4610
    GL_SPHERE_MAP = 9218
    GL_SPOT_CUTOFF = 4614
    GL_SPOT_DIRECTION = 4612
    GL_SPOT_EXPONENT = 4613
    GL_SRC_ALPHA = 770
    GL_SRC_ALPHA_SATURATE = 776
    GL_SRC_COLOR = 768
    GL_STACK_OVERFLOW = 1283
    GL_STACK_UNDERFLOW = 1284
    GL_STENCIL = 6146
    GL_STENCIL_BITS = 3415
    GL_STENCIL_BUFFER_BIT = 1024
    GL_STENCIL_CLEAR_VALUE = 2961
    GL_STENCIL_FAIL = 2964
    GL_STENCIL_FUNC = 2962
    GL_STENCIL_INDEX = 6401
    GL_STENCIL_PASS_DEPTH_FAIL = 2965
    GL_STENCIL_PASS_DEPTH_PASS = 2966
    GL_STENCIL_REF = 2967
    GL_STENCIL_TEST = 2960
    GL_STENCIL_VALUE_MASK = 2963
    GL_STENCIL_WRITEMASK = 2968
    GL_STEREO = 3123
    GL_SUBPIXEL_BITS = 3408
    GL_T = 8193
    GL_T2F_C3F_V3F = 10794
    GL_T2F_C4F_N3F_V3F = 10796
    GL_T2F_C4UB_V3F = 10793
    GL_T2F_N3F_V3F = 10795
    GL_T2F_V3F = 10791
    GL_T4F_C4F_N3F_V4F = 10797
    GL_T4F_V4F = 10792
    GL_TEXTURE = 5890
    GL_TEXTURE_1D = 3552
    GL_TEXTURE_2D = 3553
    GL_TEXTURE_ALPHA_SIZE = 32863
    GL_TEXTURE_BINDING_1D = 32872
    GL_TEXTURE_BINDING_2D = 32873
    GL_TEXTURE_BIT = 262144
    GL_TEXTURE_BLUE_SIZE = 32862
    GL_TEXTURE_BORDER = 4101
    GL_TEXTURE_BORDER_COLOR = 4100
    GL_TEXTURE_COMPONENTS = 4099
    GL_TEXTURE_COORD_ARRAY = 32888
    GL_TEXTURE_COORD_ARRAY_POINTER = 32914
    GL_TEXTURE_COORD_ARRAY_SIZE = 32904
    GL_TEXTURE_COORD_ARRAY_STRIDE = 32906
    GL_TEXTURE_COORD_ARRAY_TYPE = 32905
    GL_TEXTURE_ENV = 8960
    GL_TEXTURE_ENV_COLOR = 8705
    GL_TEXTURE_ENV_MODE = 8704
    GL_TEXTURE_GEN_MODE = 9472
    GL_TEXTURE_GEN_Q = 3171
    GL_TEXTURE_GEN_R = 3170
    GL_TEXTURE_GEN_S = 3168
    GL_TEXTURE_GEN_T = 3169
    GL_TEXTURE_GREEN_SIZE = 32861
    GL_TEXTURE_HEIGHT = 4097
    GL_TEXTURE_INTENSITY_SIZE = 32865
    GL_TEXTURE_INTERNAL_FORMAT = 4099
    GL_TEXTURE_LUMINANCE_SIZE = 32864
    GL_TEXTURE_MAG_FILTER = 10240
    GL_TEXTURE_MATRIX = 2984
    GL_TEXTURE_MIN_FILTER = 10241
    GL_TEXTURE_PRIORITY = 32870
    GL_TEXTURE_RED_SIZE = 32860
    GL_TEXTURE_RESIDENT = 32871
    GL_TEXTURE_STACK_DEPTH = 2981
    GL_TEXTURE_WIDTH = 4096
    GL_TEXTURE_WRAP_S = 10242
    GL_TEXTURE_WRAP_T = 10243
    GL_TRANSFORM_BIT = 4096
    GL_TRIANGLES = 4
    GL_TRIANGLE_FAN = 6
    GL_TRIANGLE_STRIP = 5
    GL_TRUE = 1
    GL_UNPACK_ALIGNMENT = 3317
    GL_UNPACK_LSB_FIRST = 3313
    GL_UNPACK_ROW_LENGTH = 3314
    GL_UNPACK_SKIP_PIXELS = 3316
    GL_UNPACK_SKIP_ROWS = 3315
    GL_UNPACK_SWAP_BYTES = 3312
    GL_UNSIGNED_BYTE = 5121
    GL_UNSIGNED_INT = 5125
    GL_UNSIGNED_SHORT = 5123
    GL_V2F = 10784
    GL_V3F = 10785
    GL_VENDOR = 7936
    GL_VERSION = 7938
    GL_VERSION_1_1 = 1
    GL_VERTEX_ARRAY = 32884
    GL_VERTEX_ARRAY_POINTER = 32910
    GL_VERTEX_ARRAY_SIZE = 32890
    GL_VERTEX_ARRAY_STRIDE = 32892
    GL_VERTEX_ARRAY_TYPE = 32891
    GL_VIEWPORT = 2978
    GL_VIEWPORT_BIT = 2048
    GL_XOR = 5382
    GL_ZERO = 0
    GL_ZOOM_X = 3350
    GL_ZOOM_Y = 3351

    def GLerror(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def contiguous(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glAccum(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glAlphaFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glAreTexturesResident(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glArrayElement(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glBegin(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glBindTexture(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glBitmap(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glBlendFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glCallList(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glCallLists(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glClear(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glClearAccum(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glClearColor(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glClearDepth(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glClearIndex(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glClearStencil(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glClipPlane(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor3(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor3b(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor3bv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor3d(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor3dv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor3f(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor3fv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor3i(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor3iv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor3s(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor3sv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor3ub(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor3ubv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor3ui(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor3uiv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor3us(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor3usv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor4(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor4b(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor4bv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor4d(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor4dv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor4f(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor4fv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor4i(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor4iv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor4s(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor4sv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor4ub(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor4ubv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor4ui(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor4uiv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor4us(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColor4usv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColorMask(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColorMaterial(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColorPointer(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColorPointerb(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColorPointerd(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColorPointerf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColorPointeri(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColorPointers(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColorPointerub(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColorPointerui(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColorPointerus(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColorb(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColord(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColorf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColori(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColors(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColorub(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColorui(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glColorus(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glCopyPixels(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glCopyTexImage1D(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glCopyTexImage2D(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glCopyTexSubImage1D(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glCopyTexSubImage2D(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glCullFace(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glDeleteLists(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glDeleteTextures(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glDepthFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glDepthMask(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glDepthRange(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glDisable(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glDisableClientState(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glDrawArrays(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glDrawBuffer(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glDrawElements(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glDrawElementsub(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glDrawElementsui(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glDrawElementsus(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glDrawPixels(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glDrawPixelsb(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glDrawPixelsf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glDrawPixelsi(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glDrawPixelss(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glDrawPixelsub(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glDrawPixelsui(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glDrawPixelsus(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEdgeFlag(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEdgeFlagPointer(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEdgeFlagPointerb(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEdgeFlagv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEnable(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEnableClientState(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEnd(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEndList(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEvalCoord(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEvalCoord1(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEvalCoord1d(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEvalCoord1dv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEvalCoord1f(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEvalCoord1fv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEvalCoord2(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEvalCoord2d(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEvalCoord2dv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEvalCoord2f(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEvalCoord2fv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEvalCoordd(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEvalCoordf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEvalMesh1(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEvalMesh2(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEvalPoint(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEvalPoint1(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glEvalPoint2(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glFeedbackBuffer(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glFinish(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glFlush(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glFog(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glFogf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glFogfv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glFogi(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glFogiv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glFrontFace(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glFrustum(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGenLists(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGenTextures(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetBoolean(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetBooleanv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetClipPlane(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetDouble(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetDoublev(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetFloat(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetFloatv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetInteger(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetIntegerv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetLightfv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetLightiv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetMapdv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetMapfv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetMapiv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetMaterialfv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetMaterialiv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetPixelMapfv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetPixelMapuiv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetPixelMapusv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetPolygonStipple(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetPolygonStippleub(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetString(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetTexEnvfv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetTexEnviv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetTexGendv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetTexGenfv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetTexGeniv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetTexImage(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetTexImageb(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetTexImaged(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetTexImagef(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetTexImagei(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetTexImages(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetTexImageub(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetTexImageui(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetTexImageus(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetTexLevelParameterfv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetTexLevelParameteriv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetTexParameterfv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glGetTexParameteriv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glHint(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glIndex(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glIndexMask(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glIndexPointer(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glIndexPointerb(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glIndexPointerd(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glIndexPointerf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glIndexPointeri(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glIndexPointers(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glIndexPointerub(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glIndexd(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glIndexdv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glIndexf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glIndexfv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glIndexi(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glIndexiv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glIndexs(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glIndexsv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glIndexub(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glIndexubv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glInitNames(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glInterleavedArrays(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glIsEnabled(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glIsList(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glIsTexture(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glLight(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glLightModel(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glLightModelf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glLightModelfv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glLightModeli(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glLightModeliv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glLightf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glLightfv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glLighti(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glLightiv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glLineStipple(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glLineWidth(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glListBase(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glLoadIdentity(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glLoadMatrixd(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glLoadMatrixf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glLoadName(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glLogicOp(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glMap1d(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glMap1f(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glMap2d(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glMap2f(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glMapGrid1d(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glMapGrid1f(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glMapGrid2d(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glMapGrid2f(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glMaterial(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glMaterialf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glMaterialfv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glMateriali(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glMaterialiv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glMatrixMode(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glMultMatrixd(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glMultMatrixf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNewList(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormal(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormal3(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormal3b(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormal3bv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormal3d(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormal3dv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormal3f(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormal3fv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormal3i(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormal3iv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormal3s(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormal3sv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormal4(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormalPointer(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormalPointerb(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormalPointerd(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormalPointerf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormalPointeri(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormalPointers(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormalb(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormald(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormalf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormali(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glNormals(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glOrtho(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPassThrough(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPixelMapfv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPixelMapuiv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPixelMapusv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPixelStoref(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPixelStorei(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPixelTransferf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPixelTransferi(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPixelZoom(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPointSize(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPolygonMode(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPolygonOffset(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPolygonStipple(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPolygonStippleub(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPopAttrib(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPopClientAttrib(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPopMatrix(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPopName(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPrioritizeTextures(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPushAttrib(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPushClientAttrib(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPushMatrix(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glPushName(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos2(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos2d(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos2dv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos2f(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos2fv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos2i(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos2iv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos2s(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos2sv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos3(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos3d(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos3dv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos3f(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos3fv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos3i(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos3iv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos3s(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos3sv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos4(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos4d(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos4dv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos4f(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos4fv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos4i(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos4iv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos4s(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPos4sv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPosd(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPosf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPosi(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRasterPoss(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glReadBuffer(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glReadPixels(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glReadPixelsb(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glReadPixelsd(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glReadPixelsf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glReadPixelsi(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glReadPixelss(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glReadPixelsub(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glReadPixelsui(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glReadPixelsus(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRectd(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRectdv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRectf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRectfv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRecti(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRectiv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRects(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRectsv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRenderMode(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRotate(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRotated(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glRotatef(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glScale(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glScaled(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glScalef(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glScissor(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glSelectBuffer(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glSelectWithCallback(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glShadeModel(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glStencilFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glStencilMask(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glStencilOp(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord1(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord1d(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord1dv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord1f(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord1fv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord1i(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord1iv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord1s(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord1sv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord2(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord2d(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord2dv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord2f(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord2fv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord2i(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord2iv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord2s(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord2sv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord3(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord3d(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord3dv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord3f(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord3fv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord3i(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord3iv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord3s(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord3sv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord4(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord4d(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord4dv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord4f(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord4fv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord4i(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord4iv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord4s(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoord4sv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoordPointer(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoordPointerb(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoordPointerd(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoordPointerf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoordPointeri(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoordPointers(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoordb(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoordd(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoordf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoordi(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexCoords(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexEnvf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexEnvfv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexEnvi(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexEnviv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexGen(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexGend(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexGendv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexGenf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexGenfv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexGeni(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexGeniv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexImage1D(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexImage1Db(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexImage1Df(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexImage1Di(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexImage1Ds(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexImage1Dub(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexImage1Dui(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexImage1Dus(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexImage2D(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexImage2Db(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexImage2Df(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexImage2Di(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexImage2Ds(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexImage2Dub(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexImage2Dui(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexImage2Dus(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexParameter(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexParameterf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexParameterfv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexParameteri(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexParameteriv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexSubImage1D(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexSubImage1Db(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexSubImage1Df(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexSubImage1Di(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexSubImage1Ds(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexSubImage1Dub(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexSubImage1Dui(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexSubImage1Dus(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexSubImage2D(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexSubImage2Db(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexSubImage2Df(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexSubImage2Di(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexSubImage2Ds(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexSubImage2Dub(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexSubImage2Dui(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTexSubImage2Dus(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTranslate(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTranslated(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glTranslatef(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex2d(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex2dv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex2f(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex2fv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex2i(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex2iv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex2s(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex2sv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex3d(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex3dv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex3f(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex3fv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex3i(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex3iv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex3s(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex3sv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex4d(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex4dv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex4f(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex4fv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex4i(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex4iv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex4s(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertex4sv(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertexPointer(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertexPointerb(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertexPointerd(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertexPointerf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertexPointeri(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertexPointers(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertexb(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertexd(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertexf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertexi(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glVertexs(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glViewport(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def isSequenceType(*args, **keyargs):
        raise ImportError("No module named OpenGL.GL. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

