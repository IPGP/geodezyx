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
# $Id: glrenderer.py,v 1.3 2005/08/28 19:42:43 mbaas Exp $

## \file glrenderer.py
## Contains the GLRenderer class.

from ._OpenGL.GL import *
from ._OpenGL.GLU import *
from .globalscene import getScene
from . import _core

# GLRenderInstance
class GLRenderInstance(_core.GLRenderInstance):
    """GLRenderInstance.

    """

#    protocols.advise(instancesProvide=[ISceneItem])
    
    def __init__(self):
        _core.GLRenderInstance.__init__(self)


    # pick
    def pick(self, x, y, dx=2, dy=2, root=None):
        """Do an OpenGL picking pass at the specified cursor position.

        \param x (\c float) X window coordinate where to pick
        \param y (\c float) Y window coordinate where to pick
        \param dx (\c float) Width of the picking region
        \param dy (\c float) Height of the picking region
        \param root (\c WorldObject) Only check this subtree (None=entire world).
        \return Returns a list of hits. Each hit entry is a 3-tuple (zmin, zmax, object).
        """
        
        scene = getScene()
        vpx,vpy,width,height = self.getViewport()
        if self.stereo_mode==1:
            width /= 2

        glSelectBuffer(50)
        glRenderMode(GL_SELECT)
        glInitNames()
        glPushName(0)

        glDisable(GL_DEPTH_TEST)
        glDisable(GL_LIGHTING)
        glDisable(GL_BLEND)
        glDisable(GL_NORMALIZE)

        # Projection
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPickMatrix(x, height-y, dx, dy, [vpx,vpy,width,height])
        glMultMatrixd(self.getProjection().toList())

        # Viewing transformation
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        if scene.handedness=='l':
            glScaled(-1,1,1)
        glRotated(180,0,1,0);
        glMultMatrixd(self.getViewTransformation().toList())

        # 'Render' scene
        lut = {}
        if root==None:
            root = scene.worldRoot()
        glPushMatrix()
        # Change to the coordinate system of the parent of the root
        # so that we can continue with the local system
        if root.parent!=None:
            glMultMatrixd(root.parent.worldtransform.toList())
        self._pick_node(root, 0, lut)
        glPopMatrix()

        # Get results
        buffer = glRenderMode(GL_RENDER)
        res = []
        for zmin,zmax,names in buffer:
            idx = names[0]
            res.append((zmin,zmax,lut[idx]))

        # Sort so that nearest object is first
        res.sort()
        return res

    def _pick_node(self, obj, idx, lut):
        """Helper method for the pick() method.

        Returns the current index to use for the object 'names'.
        """
        if not obj.visible:
            return idx

        glPushMatrix()
        glMultMatrixd(obj.localTransform().toList())

        geom = obj.geom
        if geom!=None:
            glLoadName(idx)
            lut[idx] = obj
            idx+=1
            geom.drawGL()

        for child in obj.iterChilds():
            idx = self._pick_node(child, idx, lut)

        glPopMatrix()
        return idx
