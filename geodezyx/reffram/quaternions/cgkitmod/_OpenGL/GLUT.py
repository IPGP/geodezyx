# Shadow module for module "OpenGL.GLUT"

try:

    # Try to import the original module...
    from OpenGL.GLUT import *
    
except:

    # Create dummy symbols...

    GLUT_ACCUM = 4
    GLUT_ACTIVE_ALT = 4
    GLUT_ACTIVE_CTRL = 2
    GLUT_ACTIVE_SHIFT = 1
    GLUT_ALPHA = 8
    GLUT_API_VERSION = 3
    GLUT_BITMAP_8_BY_13 = 3
    GLUT_BITMAP_9_BY_15 = 2
    GLUT_BITMAP_HELVETICA_10 = 6
    GLUT_BITMAP_HELVETICA_12 = 7
    GLUT_BITMAP_HELVETICA_18 = 8
    GLUT_BITMAP_TIMES_ROMAN_10 = 4
    GLUT_BITMAP_TIMES_ROMAN_24 = 5
    GLUT_BLUE = 2
    GLUT_CURSOR_BOTTOM_LEFT_CORNER = 19
    GLUT_CURSOR_BOTTOM_RIGHT_CORNER = 18
    GLUT_CURSOR_BOTTOM_SIDE = 13
    GLUT_CURSOR_CROSSHAIR = 9
    GLUT_CURSOR_CYCLE = 5
    GLUT_CURSOR_DESTROY = 3
    GLUT_CURSOR_FULL_CROSSHAIR = 102
    GLUT_CURSOR_HELP = 4
    GLUT_CURSOR_INFO = 2
    GLUT_CURSOR_INHERIT = 100
    GLUT_CURSOR_LEFT_ARROW = 1
    GLUT_CURSOR_LEFT_RIGHT = 11
    GLUT_CURSOR_LEFT_SIDE = 14
    GLUT_CURSOR_NONE = 101
    GLUT_CURSOR_RIGHT_ARROW = 0
    GLUT_CURSOR_RIGHT_SIDE = 15
    GLUT_CURSOR_SPRAY = 6
    GLUT_CURSOR_TEXT = 8
    GLUT_CURSOR_TOP_LEFT_CORNER = 16
    GLUT_CURSOR_TOP_RIGHT_CORNER = 17
    GLUT_CURSOR_TOP_SIDE = 12
    GLUT_CURSOR_UP_DOWN = 10
    GLUT_CURSOR_WAIT = 7
    GLUT_DEPTH = 16
    GLUT_DEVICE_IGNORE_KEY_REPEAT = 610
    GLUT_DEVICE_KEY_REPEAT = 611
    GLUT_DISPLAY_MODE_POSSIBLE = 400
    GLUT_DOUBLE = 2
    GLUT_DOWN = 0
    GLUT_ELAPSED_TIME = 700
    GLUT_ENTERED = 1
    GLUT_FULLY_COVERED = 3
    GLUT_FULLY_RETAINED = 1
    GLUT_GAME_MODE_ACTIVE = 0
    GLUT_GAME_MODE_DISPLAY_CHANGED = 6
    GLUT_GAME_MODE_HEIGHT = 3
    GLUT_GAME_MODE_PIXEL_DEPTH = 4
    GLUT_GAME_MODE_POSSIBLE = 1
    GLUT_GAME_MODE_REFRESH_RATE = 5
    GLUT_GAME_MODE_WIDTH = 2
    GLUT_GREEN = 1
    GLUT_HAS_DIAL_AND_BUTTON_BOX = 603
    GLUT_HAS_JOYSTICK = 612
    GLUT_HAS_KEYBOARD = 600
    GLUT_HAS_MOUSE = 601
    GLUT_HAS_OVERLAY = 802
    GLUT_HAS_SPACEBALL = 602
    GLUT_HAS_TABLET = 604
    GLUT_HIDDEN = 0
    GLUT_INDEX = 1
    GLUT_INIT_DISPLAY_MODE = 504
    GLUT_INIT_WINDOW_HEIGHT = 503
    GLUT_INIT_WINDOW_WIDTH = 502
    GLUT_INIT_WINDOW_X = 500
    GLUT_INIT_WINDOW_Y = 501
    GLUT_JOYSTICK_AXES = 615
    GLUT_JOYSTICK_BUTTONS = 614
    GLUT_JOYSTICK_BUTTON_A = 1
    GLUT_JOYSTICK_BUTTON_B = 2
    GLUT_JOYSTICK_BUTTON_C = 4
    GLUT_JOYSTICK_BUTTON_D = 8
    GLUT_JOYSTICK_POLL_RATE = 616
    GLUT_KEY_DOWN = 103
    GLUT_KEY_END = 107
    GLUT_KEY_F1 = 1
    GLUT_KEY_F10 = 10
    GLUT_KEY_F11 = 11
    GLUT_KEY_F12 = 12
    GLUT_KEY_F2 = 2
    GLUT_KEY_F3 = 3
    GLUT_KEY_F4 = 4
    GLUT_KEY_F5 = 5
    GLUT_KEY_F6 = 6
    GLUT_KEY_F7 = 7
    GLUT_KEY_F8 = 8
    GLUT_KEY_F9 = 9
    GLUT_KEY_HOME = 106
    GLUT_KEY_INSERT = 108
    GLUT_KEY_LEFT = 100
    GLUT_KEY_PAGE_DOWN = 105
    GLUT_KEY_PAGE_UP = 104
    GLUT_KEY_REPEAT_DEFAULT = 2
    GLUT_KEY_REPEAT_OFF = 0
    GLUT_KEY_REPEAT_ON = 1
    GLUT_KEY_RIGHT = 102
    GLUT_KEY_UP = 101
    GLUT_LAYER_IN_USE = 801
    GLUT_LEFT = 0
    GLUT_LEFT_BUTTON = 0
    GLUT_LUMINANCE = 512
    GLUT_MENU_IN_USE = 1
    GLUT_MENU_NOT_IN_USE = 0
    GLUT_MENU_NUM_ITEMS = 300
    GLUT_MIDDLE_BUTTON = 1
    GLUT_MULTISAMPLE = 128
    GLUT_NORMAL = 0
    GLUT_NORMAL_DAMAGED = 804
    GLUT_NOT_VISIBLE = 0
    GLUT_NUM_BUTTON_BOX_BUTTONS = 607
    GLUT_NUM_DIALS = 608
    GLUT_NUM_MOUSE_BUTTONS = 605
    GLUT_NUM_SPACEBALL_BUTTONS = 606
    GLUT_NUM_TABLET_BUTTONS = 609
    GLUT_OVERLAY = 1
    GLUT_OVERLAY_DAMAGED = 805
    GLUT_OVERLAY_POSSIBLE = 800
    GLUT_OWNS_JOYSTICK = 613
    GLUT_PARTIALLY_RETAINED = 2
    GLUT_RED = 0
    GLUT_RGB = 0
    GLUT_RGBA = 0
    GLUT_RIGHT_BUTTON = 2
    GLUT_SCREEN_HEIGHT = 201
    GLUT_SCREEN_HEIGHT_MM = 203
    GLUT_SCREEN_WIDTH = 200
    GLUT_SCREEN_WIDTH_MM = 202
    GLUT_SINGLE = 0
    GLUT_STENCIL = 32
    GLUT_STEREO = 256
    GLUT_STROKE_MONO_ROMAN = 1
    GLUT_STROKE_ROMAN = 0
    GLUT_TRANSPARENT_INDEX = 803
    GLUT_UP = 1
    GLUT_VIDEO_RESIZE_HEIGHT = 909
    GLUT_VIDEO_RESIZE_HEIGHT_DELTA = 905
    GLUT_VIDEO_RESIZE_IN_USE = 901
    GLUT_VIDEO_RESIZE_POSSIBLE = 900
    GLUT_VIDEO_RESIZE_WIDTH = 908
    GLUT_VIDEO_RESIZE_WIDTH_DELTA = 904
    GLUT_VIDEO_RESIZE_X = 906
    GLUT_VIDEO_RESIZE_X_DELTA = 902
    GLUT_VIDEO_RESIZE_Y = 907
    GLUT_VIDEO_RESIZE_Y_DELTA = 903
    GLUT_VISIBLE = 1
    GLUT_WINDOW_ACCUM_ALPHA_SIZE = 114
    GLUT_WINDOW_ACCUM_BLUE_SIZE = 113
    GLUT_WINDOW_ACCUM_GREEN_SIZE = 112
    GLUT_WINDOW_ACCUM_RED_SIZE = 111
    GLUT_WINDOW_ALPHA_SIZE = 110
    GLUT_WINDOW_BLUE_SIZE = 109
    GLUT_WINDOW_BUFFER_SIZE = 104
    GLUT_WINDOW_COLORMAP_SIZE = 119
    GLUT_WINDOW_CURSOR = 122
    GLUT_WINDOW_DEPTH_SIZE = 106
    GLUT_WINDOW_DOUBLEBUFFER = 115
    GLUT_WINDOW_FORMAT_ID = 123
    GLUT_WINDOW_GREEN_SIZE = 108
    GLUT_WINDOW_HEIGHT = 103
    GLUT_WINDOW_NUM_CHILDREN = 118
    GLUT_WINDOW_NUM_SAMPLES = 120
    GLUT_WINDOW_PARENT = 117
    GLUT_WINDOW_RED_SIZE = 107
    GLUT_WINDOW_RGBA = 116
    GLUT_WINDOW_STENCIL_SIZE = 105
    GLUT_WINDOW_STEREO = 121
    GLUT_WINDOW_WIDTH = 102
    GLUT_WINDOW_X = 100
    GLUT_WINDOW_Y = 101
    GLUT_XLIB_IMPLEMENTATION = 13

    def glutAddMenuEntry(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutAddSubMenu(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutAttachMenu(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutBitmapCharacter(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutBitmapLength(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutBitmapWidth(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutButtonBoxFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutChangeToMenuEntry(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutChangeToSubMenu(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutCopyColormap(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutCreateMenu(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutCreateSubWindow(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutCreateWindow(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutDestroyMenu(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutDestroyWindow(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutDetachMenu(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutDeviceGet(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutDialsFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutDisplayFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutEnterGameMode(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutEntryFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutEstablishOverlay(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutExtensionSupported(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutForceJoystickFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutFullScreen(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutGameModeGet(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutGameModeString(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutGet(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutGetColor(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutGetMenu(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutGetModifiers(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutGetWindow(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutHideOverlay(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutHideWindow(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutIconifyWindow(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutIdleFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutIgnoreKeyRepeat(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutInit(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutInitDisplayMode(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutInitDisplayString(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutInitWindowPosition(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutInitWindowSize(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutJoystickFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutKeyboardFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutKeyboardUpFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutLayerGet(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutLeaveGameMode(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutMainLoop(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutMenuStateFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutMenuStatusFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutMotionFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutMouseFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutOverlayDisplayFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutPassiveMotionFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutPopWindow(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutPositionWindow(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutPostOverlayRedisplay(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutPostRedisplay(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutPostWindowOverlayRedisplay(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutPostWindowRedisplay(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutPushWindow(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutRemoveMenuItem(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutRemoveOverlay(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutReportErrors(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutReshapeFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutReshapeWindow(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSetColor(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSetCursor(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSetIconTitle(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSetKeyRepeat(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSetMenu(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSetWindow(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSetWindowTitle(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSetupVideoResizing(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutShowOverlay(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutShowWindow(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSolidCone(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSolidCube(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSolidDodecahedron(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSolidIcosahedron(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSolidOctahedron(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSolidSphere(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSolidTeapot(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSolidTetrahedron(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSolidTorus(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSpaceballButtonFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSpaceballMotionFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSpaceballRotateFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSpecialFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSpecialUpFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutStopVideoResizing(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutStrokeCharacter(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutStrokeLength(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutStrokeWidth(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutSwapBuffers(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutTabletButtonFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutTabletMotionFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutTimerFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutUseLayer(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutVideoPan(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutVideoResize(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutVideoResizeGet(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutVisibilityFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutWarpPointer(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutWindowStatusFunc(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutWireCone(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutWireCube(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutWireDodecahedron(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutWireIcosahedron(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutWireOctahedron(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutWireSphere(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutWireTeapot(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutWireTetrahedron(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

    def glutWireTorus(*args, **keyargs):
        raise ImportError("No module named OpenGL.GLUT. Please install PyOpenGL (http://pyopengl.sourceforge.net/).")

