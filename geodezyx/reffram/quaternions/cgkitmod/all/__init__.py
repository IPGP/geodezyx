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
# Portions created by the Initial Developer are Copyright (C) 2004
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
# $Id: __init__.py,v 1.8 2006/05/26 21:32:01 mbaas Exp $

"""The 'all' sub package imports all names from cgkit.

You can import from this package if you simply want to make available
(almost) all the objects, functions, etc. defined in cgkit (e.g. for
interactive sessions or simple scripts):

from cgkit.all import *
"""

from cgkit import _core

from cgkit.eventmanager import eventManager
from cgkit.events import *
from cgkit.keydefs import *
from cgkit.application import getApp

from cgkit.cgtypes import vec3, vec4, mat3, mat4, quat, getEpsilon, setEpsilon, slerp, squad
from cgkit.scene import Scene, getScene
from cgkit.sceneglobals import Globals
from cgkit.worldobject import WorldObject
from cgkit.material import Material
from cgkit.glmaterial import GLMaterial, GLTexture, GLShader, GLSLANG_VERTEX, GLSLANG_FRAGMENT, GL_DECAL, GL_REPLACE, GL_BLEND, GL_MODULATE, GL_NEAREST, GL_LINEAR, GL_NEAREST_MIPMAP_NEAREST, GL_NEAREST_MIPMAP_LINEAR, GL_LINEAR_MIPMAP_NEAREST, GL_LINEAR_MIPMAP_LINEAR, GL_CLAMP, GL_REPEAT, GL_RGB, GL_RGBA
from cgkit.lightsource import LightSource
from cgkit.component import Component, createFunctionComponent
from cgkit.slots import DoubleSlot, BoolSlot, IntSlot, Vec3Slot, Vec4Slot, Mat3Slot, Mat4Slot, QuatSlot, PySlot, slotPropertyCode, ProceduralIntSlot, ProceduralDoubleSlot, ProceduralVec3Slot, ProceduralVec4Slot, ProceduralMat3Slot, ProceduralMat4Slot, ProceduralQuatSlot, NotificationForwarder, UserSizeConstraint, LinearSizeConstraint
from cgkit.slots import Dependent
from cgkit.boundingbox import BoundingBox

### Geom objects:
from cgkit.spheregeom import SphereGeom
from cgkit.ccylindergeom import CCylinderGeom
from cgkit.torusgeom import TorusGeom
from cgkit.boxgeom import BoxGeom
from cgkit.planegeom import PlaneGeom
from cgkit.trimeshgeom import TriMeshGeom
from cgkit.polyhedrongeom import PolyhedronGeom
from cgkit.drawgeom import DrawGeom
from cgkit.beziercurvegeom import BezierCurveGeom, BezierPoint

### Geometry world objects:
from cgkit.quadrics import Sphere
from cgkit.ccylinder import CCylinder
from cgkit.torus import Torus
from cgkit.box import Box
from cgkit.plane import Plane
from cgkit.trimesh import TriMesh
from cgkit.polyhedron import Polyhedron
from cgkit.draw import Draw
from cgkit.ribarchive import RIBArchive
from cgkit.beziercurve import BezierCurve

from cgkit.joint import Joint

### Dynamics world objects
from cgkit.odedynamics import ODEDynamics, ODEContactProperties, ODEBallJoint, ODEHingeJoint, ODESliderJoint, ODEHinge2Joint, ODE_COLLISION
from cgkit.joints import HingeJoint

### Camera/light
from cgkit.targetcamera import TargetCamera
from cgkit.freecamera import FreeCamera
from cgkit.lookat import LookAt
from cgkit.glpointlight import GLPointLight
from cgkit.glfreespotlight import GLFreeSpotLight
from cgkit.gltargetspotlight import GLTargetSpotLight
from cgkit.glfreedistantlight import GLFreeDistantLight
from cgkit.gltargetdistantlight import GLTargetDistantLight

from cgkit.spotlight3ds import SpotLight3DS
from cgkit.material3ds import Material3DS, TextureMap3DS
from cgkit.objmaterial import OBJMaterial, OBJTextureMap
from cgkit.mayaspotlight import MayaSpotLight

from cgkit.camcontrol import CameraControl

from cgkit.group import Group
from cgkit.tunnel import Tunnel
from cgkit.flockofbirds import FlockOfBirds
from cgkit.valuetable import ValueTable
from cgkit.expression import Expression
from cgkit.euleradapter import EulerAdapter
from cgkit.pidcontroller import PIDController
from cgkit.gnuplotter import GnuPlotter
from cgkit.slideshow import SlideShow, Slide, XFade, XCube
from cgkit.motionpath import MotionPath

from cgkit.glrenderer import GLRenderInstance

from cgkit.joystick import Joystick

CONSTANT = _core.VarStorage.CONSTANT
UNIFORM = _core.VarStorage.UNIFORM
VARYING = _core.VarStorage.VARYING
VERTEX = _core.VarStorage.VERTEX
FACEVARYING = _core.VarStorage.FACEVARYING
FACEVERTEX = _core.VarStorage.FACEVERTEX
USER = _core.VarStorage.USER

INT = _core.VarType.INT
FLOAT = _core.VarType.FLOAT
STRING = _core.VarType.STRING
COLOR = _core.VarType.COLOR
POINT = _core.VarType.POINT
VECTOR = _core.VarType.VECTOR
NORMAL = _core.VarType.NORMAL
MATRIX = _core.VarType.MATRIX
HPOINT = _core.VarType.HPOINT

### Importer
import cgkit.offimport
import cgkit.pyimport
import cgkit.dddsimport
import cgkit.x3dimport
import cgkit.ifsimport
import cgkit.objimport
import cgkit.stlimport
import cgkit.asfamcimport
import cgkit.bvhimport
import cgkit.maimport
import cgkit.plyimport
import cgkit.lwobimport
### Exporter
import cgkit.ribexport
import cgkit.offexport
import cgkit.objexport
import cgkit.plyexport

from cgkit.rmshader import RMMaterial, RMLightSource, RMShader
from cgkit.ribexport import ShadowPass, FlatReflectionPass

from cgkit.cmds import *
