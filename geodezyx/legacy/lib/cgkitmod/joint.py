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
# $Id: joint.py,v 1.5 2005/08/28 19:42:43 mbaas Exp $

## \file joint.py
## Contains the Joint class.

from ._OpenGL.GL import *
from .cgtypes import *
from .component import *
from .worldobject import WorldObject
from .geomobject import GeomObject
from .spheregeom import SphereGeom
from .boundingbox import BoundingBox
from .euleradapter import EulerAdapter
from . import _core
from . import cmds
from .sl import *

# JointGeom
class JointGeom(GeomObject):
    """The geom object class for Joint objects.

    The geometry draws a wire frame sphere at the joint position
    and connections to all children joints (visualizing the bones). 
    """

    def __init__(self, joint, radius=0.05):
        """Constructor.

        \param joint (\c Joint) The joint that uses this geom
        \param radius (\c float) Radius of the wire sphere at the joint position
        """
        GeomObject.__init__(self)
        self.joint = joint
        self.radius = radius

        self.spheregeom = SphereGeom(radius=self.radius, segmentsu=8, segmentsv=4)

    def uniformCount(self):
        return 0

    def varyingCount(self):
        return 0

    def vertexCount(self):
        return 0

    def boundingBox(self):
        r = vec3(self.radius, self.radius, self.radius)
        return BoundingBox(-r, r)

    def drawGL(self):
        # Draw sphere
        glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT | GL_POLYGON_BIT)
        glDisable(GL_LIGHTING)
        glPushMatrix()

        # The joint is located at the pivot point
        p = self.joint.getOffsetTransform()[3]
        pivot = vec3(p.x, p.y, p.z)
        glTranslate(pivot.x, pivot.y, pivot.z)
        
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
        self.spheregeom.drawGL()

        # Draw bone
        glBegin(GL_LINES)
        for child in self.joint.iterChilds():
            if not isinstance(child, Joint):
                continue
            p = child.pos-pivot
            try:
                b1 = p.ortho().normalize()
                b2 = p.cross(b1).normalize()
            except:
                continue
            b1 *= 0.75*self.radius
            b2 *= 0.75*self.radius
            
            glVertex3d(b1.x, b1.y, b1.z)
            glVertex3d(p.x, p.y, p.z)
            
            glVertex3d(-b1.x, -b1.y, -b1.z)
            glVertex3d(p.x, p.y, p.z)

            glVertex3d(b2.x, b2.y, b2.z)
            glVertex3d(p.x, p.y, p.z)

            glVertex3d(-b2.x, -b2.y, -b2.z)
            glVertex3d(p.x, p.y, p.z)
        glEnd()

        # Pivot coordinate system...
        P = self.joint.getOffsetTransform()
        r = 1.5*self.radius
        glBegin(GL_LINES)
        # X axis
        b = r*P[0]
        glColor3f(1,0,0)
        glVertex3f(0,0,0)
        glVertex3f(b.x,b.y,b.z)
        # Y axis
        b = r*P[1]
        glColor3f(0,1,0)
        glVertex3f(0,0,0)
        glVertex3f(b.x,b.y,b.z)
        # Z axis
        b = r*P[2]
        glColor3f(0,0,1)
        glVertex3f(0,0,0)
        glVertex3f(b.x,b.y,b.z)
        glEnd()

        glPopMatrix()
        glPopAttrib()


class Mult(Component):

    def __init__(self, name="Mult", auto_insert=True):
        Component.__init__(self, name=name, auto_insert=auto_insert)

        # Create the input slots
        self.op1_slot = Mat3Slot()
        self.op2_slot = Mat3Slot()
        # Create the output slot
        self.output_slot = ProceduralMat3Slot(self.computeOutput)

        # Add the slots to the component
        self.addSlot("op1", self.op1_slot)
        self.addSlot("op2", self.op2_slot)
        self.addSlot("output", self.output_slot)

        # Set up slot dependencies
        self.op1_slot.addDependent(self.output_slot)
        self.op2_slot.addDependent(self.output_slot)

    def computeOutput(self):
        return self.op1*self.op2

    # Create value attributes
    exec(slotPropertyCode("op1"))
    exec(slotPropertyCode("op2"))
    exec(slotPropertyCode("output"))

# Joint
class Joint(WorldObject):
    """Joint class.
    """

    def __init__(self,
                 name = "Joint",
                 radius = 0.05,
                 rotationorder = "xyz",
                 **params):
        
        WorldObject.__init__(self, name=name, **params)

        self.geom = JointGeom(self, radius=radius)

        self.euleradapter = EulerAdapter(order=rotationorder, outtype="mat3", auto_insert=False)
        self.nullrot = Mult(auto_insert=False)
        
#        self.euleradapter.output_slot.connect(self.rot_slot)
        self.nullrot.op1 = mat3(1)
        self.euleradapter.output_slot.connect(self.nullrot.op2_slot)
        self.nullrot.output_slot.connect(self.rot_slot)

        self.anglex_slot = self.euleradapter.anglex_slot
        self.angley_slot = self.euleradapter.angley_slot
        self.anglez_slot = self.euleradapter.anglez_slot

#    def setEuler(self, x,y,z):
#        m = mat3().fromEulerXYZ(radians(x), radians(y), radians(z))
#        self.rot = m

    ## protected:

    # angle properties...
    exec(slotPropertyCode("anglex"))
    exec(slotPropertyCode("angley"))
    exec(slotPropertyCode("anglez"))

    # freezePivot
    def freezePivot(self):
        """Make the current pivot coordinate system the default pose.

        After calling this method, the current rotation of the pivot
        coordinate system will define the default pose. This means,
        rotations are now defined around the local pivot axes.
        """
        self.nullrot.op1 = self.getOffsetTransform().getMat3()
        
