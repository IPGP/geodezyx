# Shadow module for module "OpenGL.GLU"

try:

    # Try to import the original module...
    from OpenGL.GLU import *
    
except:

    # Create dummy symbols...

    GLU_AUTO_LOAD_MATRIX = 100200
    GLU_BEGIN = 100100
    GLU_CCW = 100121
    GLU_CULLING = 100201
    GLU_CW = 100120
    GLU_DISPLAY_MODE = 100204
    GLU_DOMAIN_DISTANCE = 100217
    GLU_EDGE_FLAG = 100104
    GLU_END = 100102
    GLU_ERROR = 100103
    GLU_EXTENSIONS = 100801
    GLU_EXTERIOR = 100123
    GLU_FILL = 100012
    GLU_FLAT = 100001
    GLU_INCOMPATIBLE_GL_VERSION = 100903
    GLU_INSIDE = 100021
    GLU_INTERIOR = 100122
    GLU_INVALID_ENUM = 100900
    GLU_INVALID_VALUE = 100901
    GLU_LINE = 100011
    GLU_MAP1_TRIM_2 = 100210
    GLU_MAP1_TRIM_3 = 100211
    GLU_NONE = 100002
    GLU_NURBS_ERROR1 = 100251
    GLU_NURBS_ERROR10 = 100260
    GLU_NURBS_ERROR11 = 100261
    GLU_NURBS_ERROR12 = 100262
    GLU_NURBS_ERROR13 = 100263
    GLU_NURBS_ERROR14 = 100264
    GLU_NURBS_ERROR15 = 100265
    GLU_NURBS_ERROR16 = 100266
    GLU_NURBS_ERROR17 = 100267
    GLU_NURBS_ERROR18 = 100268
    GLU_NURBS_ERROR19 = 100269
    GLU_NURBS_ERROR2 = 100252
    GLU_NURBS_ERROR20 = 100270
    GLU_NURBS_ERROR21 = 100271
    GLU_NURBS_ERROR22 = 100272
    GLU_NURBS_ERROR23 = 100273
    GLU_NURBS_ERROR24 = 100274
    GLU_NURBS_ERROR25 = 100275
    GLU_NURBS_ERROR26 = 100276
    GLU_NURBS_ERROR27 = 100277
    GLU_NURBS_ERROR28 = 100278
    GLU_NURBS_ERROR29 = 100279
    GLU_NURBS_ERROR3 = 100253
    GLU_NURBS_ERROR30 = 100280
    GLU_NURBS_ERROR31 = 100281
    GLU_NURBS_ERROR32 = 100282
    GLU_NURBS_ERROR33 = 100283
    GLU_NURBS_ERROR34 = 100284
    GLU_NURBS_ERROR35 = 100285
    GLU_NURBS_ERROR36 = 100286
    GLU_NURBS_ERROR37 = 100287
    GLU_NURBS_ERROR4 = 100254
    GLU_NURBS_ERROR5 = 100255
    GLU_NURBS_ERROR6 = 100256
    GLU_NURBS_ERROR7 = 100257
    GLU_NURBS_ERROR8 = 100258
    GLU_NURBS_ERROR9 = 100259
    GLU_OUTLINE_PATCH = 100241
    GLU_OUTLINE_POLYGON = 100240
    GLU_OUTSIDE = 100020
    GLU_OUT_OF_MEMORY = 100902
    GLU_PARAMETRIC_ERROR = 100216
    GLU_PARAMETRIC_TOLERANCE = 100202
    GLU_PATH_LENGTH = 100215
    GLU_POINT = 100010
    GLU_SAMPLING_METHOD = 100205
    GLU_SAMPLING_TOLERANCE = 100203
    GLU_SILHOUETTE = 100013
    GLU_SMOOTH = 100000
    GLU_TESS_BEGIN = 100100
    GLU_TESS_BEGIN_DATA = 100106
    GLU_TESS_BOUNDARY_ONLY = 100141
    GLU_TESS_COMBINE = 100105
    GLU_TESS_COMBINE_DATA = 100111
    GLU_TESS_COORD_TOO_LARGE = 100155
    GLU_TESS_EDGE_FLAG = 100104
    GLU_TESS_EDGE_FLAG_DATA = 100110
    GLU_TESS_END = 100102
    GLU_TESS_END_DATA = 100108
    GLU_TESS_ERROR = 100103
    GLU_TESS_ERROR1 = 100151
    GLU_TESS_ERROR2 = 100152
    GLU_TESS_ERROR3 = 100153
    GLU_TESS_ERROR4 = 100154
    GLU_TESS_ERROR5 = 100155
    GLU_TESS_ERROR6 = 100156
    GLU_TESS_ERROR7 = 100157
    GLU_TESS_ERROR8 = 100158
    GLU_TESS_ERROR_DATA = 100109
    GLU_TESS_MAX_COORD = 9.9999999999999998e+149
    GLU_TESS_MISSING_BEGIN_CONTOUR = 100152
    GLU_TESS_MISSING_BEGIN_POLYGON = 100151
    GLU_TESS_MISSING_END_CONTOUR = 100154
    GLU_TESS_MISSING_END_POLYGON = 100153
    GLU_TESS_NEED_COMBINE_CALLBACK = 100156
    GLU_TESS_TOLERANCE = 100142
    GLU_TESS_VERTEX = 100101
    GLU_TESS_VERTEX_DATA = 100107
    GLU_TESS_WINDING_ABS_GEQ_TWO = 100134
    GLU_TESS_WINDING_NEGATIVE = 100133
    GLU_TESS_WINDING_NONZERO = 100131
    GLU_TESS_WINDING_ODD = 100130
    GLU_TESS_WINDING_POSITIVE = 100132
    GLU_TESS_WINDING_RULE = 100140
    GLU_UNKNOWN = 100124
    GLU_U_STEP = 100206
    GLU_VERSION = 100800
    GLU_VERSION_1_1 = 1
    GLU_VERSION_1_2 = 1
    GLU_VERTEX = 100101
    GLU_V_STEP = 100207

    def GLUerror(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluBeginCurve(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluBeginPolygon(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluBeginSurface(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluBeginTrim(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluBuild1DMipmaps(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluBuild1DMipmapsb(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluBuild1DMipmapsf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluBuild1DMipmapsi(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluBuild1DMipmapss(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluBuild1DMipmapsub(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluBuild1DMipmapsui(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluBuild1DMipmapsus(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluBuild2DMipmaps(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluBuild2DMipmapsb(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluBuild2DMipmapsf(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluBuild2DMipmapsi(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluBuild2DMipmapss(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluBuild2DMipmapsub(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluBuild2DMipmapsui(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluBuild2DMipmapsus(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluCylinder(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluDeleteNurbsRenderer(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluDeleteQuadric(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluDeleteTess(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluDisk(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluEndCurve(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluEndPolygon(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluEndSurface(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluEndTrim(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluErrorString(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluGetNurbsProperty(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluGetString(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluGetTessProperty(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluLoadSamplingMatrices(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluLookAt(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluNewNurbsRenderer(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluNewQuadric(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluNewTess(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluNextContour(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluNurbsCallback(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluNurbsCurve(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluNurbsProperty(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluNurbsSurface(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluOrtho2D(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluPartialDisk(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluPerspective(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluPickMatrix(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluProject(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluPwlCurve(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluQuadricCallback(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluQuadricDrawStyle(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluQuadricNormals(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluQuadricOrientation(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluQuadricTexture(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluScaleImage(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluScaleImageb(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluScaleImagef(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluScaleImagei(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluScaleImages(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluScaleImageub(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluScaleImageui(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluScaleImageus(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluSphere(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluTessBeginContour(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluTessBeginPolygon(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluTessCallback(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluTessEndContour(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluTessEndPolygon(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluTessNormal(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluTessProperty(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluTessVertex(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def gluUnProject(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLU. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

