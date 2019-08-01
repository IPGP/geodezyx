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
# $Id: camerabase.py,v 1.1 2005/05/12 12:58:28 mbaas Exp $

## \file camerabase.py
## Contains the CameraBase class.

"""This module contains the CameraBase class."""

from .Interfaces import *
from . import protocols
from . import slots
from .cgtypes import *
from math import pi
from .worldobject import WorldObject
from .globalscene import getScene
from . import _core

# CameraBase
class CameraBase(WorldObject):
    """Base class for camera objects.
    """

    def __init__(self,
                 auto_nearfar = True,
                 nearplane = 0.1,
                 farplane = 1000,
                 **params):
        WorldObject.__init__(self, **params)

        self.nearplane_slot = slots.DoubleSlot(nearplane)
        self.farplane_slot = slots.DoubleSlot(farplane)
        self.autonearfar_slot = slots.BoolSlot(auto_nearfar)
        self.addSlot("nearplane", self.nearplane_slot)
        self.addSlot("farplane", self.farplane_slot)
        self.addSlot("autonearfar", self.autonearfar_slot)

    # "output" property...
    exec(slots.slotPropertyCode("nearplane"))
    exec(slots.slotPropertyCode("farplane"))
    exec(slots.slotPropertyCode("autonearfar"))
    

    def protocols(self):
        return [ISceneItem, IComponent, IWorldObject, ICamera]

    # eyeRay
    def eyeRay(self, x0, y0, width, height):
        """Return a ray from the eye position through an image point.

        This method returns a ray whose origin is at the eye position
        and that goes through a given point on the image plane. The
        point on the plane is given by (x0, y0) which each ranges from
        0 to 1. (0,0) is at the upper left and (1,1) at the lower right.
        The arguments width and height determine the ratio of the image
        plane (the absolute values of width and height are irrelevant).
        The return value is a 2-tuple (p,u) where p is the ray origin
        and u the normalized direction. Both vectors are given in world
        space.
        """
        V = self.viewTransformation()
        P = self.projection(width, height, 1, 10)
        R = mat4().rotation(pi, vec3(0,1,0))
        if getScene().handedness=='l':            
            S = mat4().scaling(vec3(-1,1,1))
            I = (P*S*R*V).inverse()
        else:
            I = (P*R*V).inverse()
        x = 2.0*x0-1.0
        y = 1.0-2.0*y0
        q = I*vec3(x,y,0)
        p = self.worldtransform[3]
        p = vec3(p.x, p.y, p.z)
        return (p, (q-p).normalize())

    # getNearFar
    def getNearFar(self):
        """Return the distances to the near and far clipping plane.

        If auto_nearfar is True, the near/far values are computed from
        the scene extent, otherwise the stored values are used.
        
        Compute near and far clipping plane distances from the bounding
        box of the scene. The scene bounding box is converted to a
        bounding sphere and the near and far clipping planes are set
        as tangent planes to the bounding sphere.
        """

        if not self.autonearfar:
            return self.nearplane, self.farplane

        # Get the bounding box of the entire scene
        bbox = getScene().boundingBox()

        # Determine bounding sphere
        bmin,bmax = bbox.getBounds()
        if bmin!=None and bmin!=bmax:
            # Box center (resp. sphere center)
            c = 0.5*(bmin+bmax)
            # Radius of the bounding sphere
            r = (bmax-c).length()
        else:
            c = vec3(0,0,0)
            r = 10

        # Transformation World->Camera
        VT = self.viewTransformation()

#        minnear = (bmax-bmin).length()/1000
        minnear = self.nearplane
        minfar = minnear+1

        # cz: Depth of the center point
        cz = (VT*c).z
        near = max(cz-r, minnear)
        far  = max(cz+r, minfar)

        if (far-near)<0.01:
            far+=1

        return (near,far)
