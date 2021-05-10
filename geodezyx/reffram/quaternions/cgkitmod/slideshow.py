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

from ._OpenGL.GL import *
from . import _Image as Image
from . import _ImageDraw as ImageDraw
import glob, types
from . import component
from .eventmanager import eventManager
from .events import *
from .globalscene import getScene
from .targetcamera import TargetCamera
from .plane import Plane
from .quadrics import Sphere
from .glmaterial import GLMaterial, GLTexture
from .sl import *
from .slots import *
from .cgtypes import *
from math import *
from . import cmds
from . import _core

# Cube transition
class XCube:
    def __init__(self, duration=2.0):
        self.duration = duration
    
    # preTransition
    def preTransition(self, slideshow):
        """Prepare transition.

        This method is called at the beginning of a transition.
        """
        slideshow.backplate.rot = mat3().rotation(pi/2, vec3(0,1,0))
        slideshow.backplate.pos = vec3(2,0,-2)
        slideshow.helper.pos = vec3(0,0,-2)
        slideshow.helper.rot = mat3(1)
        
        cmds.link(slideshow.frontplate, slideshow.helper)
        cmds.link(slideshow.backplate, slideshow.helper)

    # doTransition
    def doTransition(self, slideshow, t):
        """Do the transition.

        This method is called every frame during the transition.
        t is the current transition value (0-1).
        At t=0 the frontplate shows and the backplate is not visible,
        at t=1 the frontplate should not be visible anymore and the
        backplate has come to sight
        """
        w = smoothstep(0,1,smoothstep(0, 1, t))
        slideshow.helper.rot = mat3().rotation(-w*pi/2, vec3(0,1,0))

    # postTransition
    def postTransition(self, slideshow):
        """Post transition.

        This method is called at the end of a transition.
        post: frontplate must be visible
        """
        slideshow.helper.rot = mat3(1)
        cmds.link(slideshow.frontplate, None)
        cmds.link(slideshow.backplate, None)


# Crossfade
class XFade:
    def __init__(self, duration=2.0, zmove=0.0):
        self.duration = duration
        self.zmove = zmove

    ############ XFade transition ####################

    def preTransition(self, slideshow):
        slideshow.backplate.pos = slideshow.frontplate.pos+vec3(0,0,-0.001)
        mat = slideshow.frontplate.getMaterial()
        mat.blend_sfactor = GL_SRC_ALPHA
        mat.blend_dfactor = GL_ONE_MINUS_SRC_ALPHA

    def doTransition(self, slideshow, t):
        w = smoothstep(0,1,t)
        mat = slideshow.frontplate.getMaterial()
        mat.diffuse = vec4(1,1,1,1-w)
        slideshow.frontplate.pos = vec3(0,0,self.zmove*t*t)

    def postTransition(self, slideshow):
        mat = slideshow.frontplate.getMaterial()
        mat.diffuse = vec4(1,1,1,1)
        mat.blend_sfactor = -1
        mat.blend_dfactor = -1
        mat = slideshow.backplate.getMaterial()
        mat.diffuse = vec4(1,1,1,1)

######################################################################

# Slide
class Slide:
    def __init__(self, filepattern, transition=XCube()):
        self.files = glob.glob(filepattern)
        self.transition = transition

######################################################################

# SlideShow
class SlideShow(component.Component):
    """Slide show component.
    """
    
    def __init__(self, name="SlideShow", slides=[], auto_insert=True):
        
        component.Component.__init__(self, name, auto_insert)

        if isinstance(slides, str):
            slides = Slide(slides)

        try:
            len(slides)
        except:
            slides = [slides]

        self.images = []
        for s in slides:
            for f in s.files:
                self.images.append((f, s.transition))

        # Currently displayed image index
        self.imageidx = 0

        self.setup()

        em = eventManager()
        em.connect(STEP_FRAME, self)
        em.connect(LEFT_DOWN, self)
        em.connect(KEY_PRESS, self)

        self.transition = self.images[0][1]

        # Global time when the next transition starts
#        self.transition_start = 3.0
        self.transition_start = 999999.0
        # Slide show state
        self.state = 0


    # setup
    def setup(self):
        """Setup the scene.
        """
        scene = getScene()

        # Set up direction
        scene.up = (0,1,0)

        # Camera
        self.cam = TargetCamera(
            name = "SlideShow_Cam",
            pos = (0,0,5.6),
            fov = 30
        )

        # Set background...
        gradient = Image.new("RGB", (2,128))
        draw = ImageDraw.Draw(gradient)
        colbottom = vec3(0.1, 0.1, 0.3)
        colmid = vec3(0.95,0.95,1)
        coltop = vec3(0.2, 0.2, 0.4)
        colormap = [colbottom, colbottom, colmid, coltop, coltop]
        for y in range(128):
            t = float(y)/128
            c = spline(t, colormap)
            draw.line([0,y, 1,y], fill=(int(c.x*255), int(c.y*255), int(c.z*255)))

        back = Plane(
            lx=4,
            ly=3,
            pos=(0,0,-5),
            scale=2.4,
            material = GLMaterial( diffuse = (0,1,1),
                                   texture = GLTexture(image = gradient,
                                                       mode  = GL_REPLACE))
        )


        # Create plate materials...
        initial = Image.new("RGB", (4,4))
        invY = mat4(1).translate(vec3(0,1,0)).scale(vec3(1,-1,1))
        
        self.mat1 = GLMaterial(
           diffuse = (1,1,1,1),
           texture = GLTexture( image = initial,
                                mode = GL_REPLACE,
#                                size = (512,512),
                                transform = invY,
                                wrap_s = GL_CLAMP,
                                wrap_t = GL_CLAMP
                               )
        )
        self.backmaterial = self.mat1
        self.setBackImage(self.images[0][0])

        self.mat2 = GLMaterial(
           diffuse = (1,1,1,1),
           texture = GLTexture( image = initial,
                                mode = GL_REPLACE,
#                                size = (512,512),
                                transform = invY,
                                wrap_s = GL_CLAMP,
                                wrap_t = GL_CLAMP
                               )
        )
        self.backmaterial = self.mat2
        self.setBackImage(self.images[1][0])

        self.frontmaterial = self.mat1
        self.backmaterial = self.mat2

        # Create plates...
        self.frontplate = Plane(lx = 4, ly = 3, material = self.frontmaterial)

        self.backplate = Plane(lx = 4, ly = 3,
                               pos = (0,0,-0.1),
                               material = self.backmaterial
                               )

        self.helper = Sphere(radius=0.1, pos=(0,0,-1))


    # onLeftDown
    def onLeftDown(self, e):
        if self.state==0:
            self.transition_start = getScene().timer().time

    def onKeyPress(self, e):
        # Page down or Enter?
        if e.keycode==281 or e.keycode==13:
            self.transition_start = getScene().timer().time
        # Page up
        elif e.keycode==280:
            pass
        # q (reset cam)
        elif e.key=="q":
            getScene().up = vec3(0,1,0)
            self.cam.pos = vec3(0,0,5.6)
            self.cam.target = vec3(0,0,0)
            self.cam.fov = 30
            cmds.link(self.cam)
        
    # onStepFrame
    def onStepFrame(self):
        timer = getScene().timer()

        ### State 0: Waiting for the transition
        if self.state==0:
            if timer.time>=self.transition_start:
                self.state = 1

        ### State 1: Begin of transition
        if self.state==1:
            self.transition.preTransition(self)
            eventManager().event("SlidePreTransition", self.images[self.imageidx][0])
            self.state = 2
            
        ### State 2: Transition
        if self.state==2:
            if self.transition.duration>0:
                t = (timer.time-self.transition_start)/self.transition.duration
            else:
                t = 1.0
            if t>1.0:
                t = 1.0
            self.transition.doTransition(self, t)
            eventManager().event("SlideDoTransition", t)
            if t==1.0:
                self.state = 3
                
        ### State 3: Transition is over
        elif self.state==3:
            eventManager().event("SlidePostTransition")
            self.transition.postTransition(self)
            self.switchSlide()

            self.frontplate.transform = mat4(1)
            self.backplate.transform = mat4(1)
            self.backplate.pos = vec3(0,0,-1)
                       
#            self.transition_start += 5.0
            self.transition_start = 999999.0
            self.state = 0

    # switchSlide
    def switchSlide(self):
        """Prepare everything for the next slide.
        """
        # Swap front and back material...
        dummy = self.frontmaterial
        self.frontmaterial = self.backmaterial
        self.backmaterial = dummy

        # Update the materials of the plates
        self.frontplate.setMaterial(self.frontmaterial)
        self.backplate.setMaterial(self.backmaterial)

        # Set the next image on the back material
        idx = (self.imageidx+2)%len(self.images)
        self.setBackImage(self.images[idx][0])
        self.imageidx = (self.imageidx+1)%len(self.images)

        # Set the new transition object
        self.transition = self.images[self.imageidx][1]
        

    # setBackImage
    def setBackImage(self, name):
        """Set a new image on the back material."""
        self.backmaterial.texture.imagename = name
        img = Image.open(name)
        if img.mode!="RGB" and img.mode!="RGBA":
            img = img.convert("RGB")
        self.backmaterial.texture.image = img
        w,h = img.size
        f = (4.0*h)/(3.0*w)
        if f>=1.0:
            S = mat4(1).scale(vec3(f,1,1))
            S.translate(vec3(-0.5+0.5/f,0,0))
            draw = ImageDraw.Draw(img)
            draw.line([(0,0), (0,h)], fill=(0,0,0))
            draw.line([(1,0), (1,h)], fill=(0,0,0))
            draw.line([(w-1,0), (w-1,h)], fill=(0,0,0))
            draw.line([(w-2,0), (w-2,h)], fill=(0,0,0))
        else:
            f = (3.0*w)/(4.0*h)
            S = mat4(1).scale(vec3(1,f,1))
            S.translate(vec3(0,-0.5+0.5/f,0))
            draw = ImageDraw.Draw(img)
            draw.line([(0,0), (w,0)], fill=(0,0,0))
            draw.line([(0,1), (w,1)], fill=(0,0,0))
            draw.line([(0,h-1), (w,h-1)], fill=(0,0,0))
            draw.line([(0,h-2), (w,h-2)], fill=(0,0,0))
        invY = mat4(1).translate(vec3(0,1,0)).scale(vec3(1,-1,1))
        self.backmaterial.texture.transform = invY*S
               
            

######################################################################

#sl = SlideShow()
