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
# $Id: bvhimport.py,v 1.1 2005/02/11 14:32:34 mbaas Exp $

from .cgtypes import *
from .joint import Joint
from .valuetable import ValueTable
from . import bvh
from . import pluginmanager
from .sl import *

# BVHReader
class BVHReader(bvh.BVHReader):
    """Specialized BVH reader class.

    This class creates a hierarchy of joints and applies the motion to it.
    """
    
    def __init__(self, filename):
        bvh.BVHReader.__init__(self, filename)
      

    def onHierarchy(self, root):
        self.createSkeleton(root)
        self.root = root

    def onMotion(self, frames, dt):
        self.frames = frames
        self.dt = dt
        self.currentframe = 0

    def onFrame(self, values):
        self.applyMotion(self.root, values)
        self.currentframe += 1

    def applyMotion(self, node, values):
        """Apply a motion sample to the skeleton.

        node is the current joint and values the joint angles for the
        entire skeleton.
        The method returns the remaining joint angles.
        """

        t = self.currentframe*self.dt
        
        nc = len(node.channels)
        vals = values[:nc]
        pos = vec3()
        pos_flag = False
        for ch,v in zip(node.channels, vals):
            if ch=="Xrotation":
                node.vtx.add(t, v)
            elif ch=="Yrotation":
                node.vty.add(t, v)
            elif ch=="Zrotation":
                node.vtz.add(t, v)
            elif ch=="Xposition":
                pos.x = v
                pos_flag = True
            elif ch=="Yposition":
                pos.y = v
                pos_flag = True
            elif ch=="Zposition":
                pos.z = v
                pos_flag = True

        if pos_flag:
            node.vtpos.add(t, pos)
            
        values = values[nc:]
        for c in node.children:
            values = self.applyMotion(c, values)
        return values

    # createSkeleton
    def createSkeleton(self, node, parent=None):
        """Create the skeleton hierarchy.

        This method creates the skeleton recursively. Each invocation
        creates one joint.
        """
        order = self.rotationOrder(node.channels)
        # Create a new Joint object
        j = Joint(name = node.name,
                  pos = vec3(node.offset),
                  rotationorder = order,
                  parent = parent)
        # Store the joint in the node so that later the motion can be applied
        node.joint = j
        
        vtx = ValueTable(type="double")
        vty = ValueTable(type="double")
        vtz = ValueTable(type="double")
        vtx.output_slot.connect(j.anglex_slot)
        vty.output_slot.connect(j.angley_slot)
        vtz.output_slot.connect(j.anglez_slot)
        node.vtx = vtx
        node.vty = vty
        node.vtz = vtz
        if node.isRoot():
            vtpos = ValueTable(type="vec3")
            vtpos.output_slot.connect(j.pos_slot)
            node.vtpos = vtpos
            
        for c in node.children:
            self.createSkeleton(c, j)

    # rotationOrder
    def rotationOrder(self, channels):
        """Determine rotation order string from the channel names.
        """
        res = ""
        for c in channels:
            if c[-8:]=="rotation":
                res += c[0]

        # Complete the order string if it doesn't already contain
        # all three axes
        m = { "":"XYZ",
              "X":"XYZ", "Y":"YXZ", "Z":"ZXY",
              "XY":"XYZ", "XZ":"XZY",
              "YX":"YXZ", "YZ":"YZX",
              "ZX":"ZXY", "ZY":"ZYX" }
        if res in m:
            res = m[res]
        return res

######################################################################

# BVHImporter
class BVHImporter:

    _protocols = ["Import"]

    # extension
    def extension():
        """Return the file extensions for this format."""
        return ["bvh"]
    extension = staticmethod(extension)

    # description
    def description(self):
        """Return a short description for the file dialog."""
        return "Biovision Hierarchical"
    description = staticmethod(description)

    # importFile
    def importFile(self, filename):
        """Import a BVH file."""
        
        bvh = BVHReader(filename)
        bvh.read()


######################################################################

# Register the Importer class as a plugin class
pluginmanager.register(BVHImporter)

