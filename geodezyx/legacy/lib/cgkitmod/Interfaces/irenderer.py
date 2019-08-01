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
# $Id: irenderer.py,v 1.1.1.1 2004/12/12 14:31:42 mbaas Exp $

class IRenderer:
    """The basic renderer interface.
    
    The renderer is responsible for displaying the scene on the screen
    in real time (i.e. frame rates that allow working with the scene
    interactively). 
    
    There are actually two classes involved in rendering, the Renderer
    (this class) and a RenderInstance class (IRenderInstance). The
    Renderer class represents a particlar rendering algorithm or
    render engine while the RenderInstance is used for rendering on
    a particular widget. All RenderInstances that were spawned from
    the same Renderer can share resources like display lists, textures, etc.
        
    Before a scene can be rendered a renderer object has to be created
    that is connected to one particular scene (which is passed to the 
    constructor). Then, for each view you create a RenderInstance object
    by calling createRenderInstance().
    On a render instance you can set various parameters like render flags, 
    camera and viewport parameters via 
    \link IRenderInstance::setRenderFlags setRenderFlags() \endlink, 
    \link IRenderInstance::setCamera setCamera() \endlink and
    \link IRenderInstance::setViewport setViewport() \endlink. Once 
    everything is set up, the scene can be
    rendered by calling \link IRenderInstance::paint paint() \endlink for 
    every frame. A picking operation can be performed using the 
    \link IRenderInstance::pick pick() \endlink method.
    
    Of course, an actual implementation of the interface may actually
    contain more methods to activate special features supported by
    this particular renderer (e.g. shadows, shaders, ...).
    
    There are two "flavors" of renderers, one that expects the application
    to provide the rendering widget and one that creates its own window.
    The interface is the same, only the semantics of some methods change
    slightly. If a renderer creates its own window it has to do so in
    the createRenderInstance() method and close the window again in 
    \link IRenderInstance::close RenderInstance.close() \endlink. 
    An embedded renderer can assume that its 
    associated rendering context was previously made current by the 
    application for most of the calls (as noted in the
    method documentations).
    
    \see IRenderInstance
    """
        
    def __init__(self, scene):
        """Constructor.
        
        \param scene (\c IScene) The scene interface of the scene that
            is to be rendered.
        """
        
    def createRenderInstance(self):
        """Open & initialize a render instance.
        
        The renderer will set default values for the flags, camera
        and viewport.
        
        If the renderer is embedded this method may only be called if
        the associated rendering context (OpenGL context) is made current.        

        \return Render instance object (\c RenderInstance)
        """
        
        
        
class IRenderInstance:
    """The render instance.
    
    \see IRenderer
    """

    def close(self):
        """Close the renderer interface.
        
        Frees any resources that were previously allocated.
        """

    def setRenderFlags(self, flags):
        """Set the flags that determine how the scene is rendered.
        
        The \a flags parameter is a combination of the following flags:
            
        - \c SOLID: Draw the scene using shaded solid faces
        - \c WIREFRAME: Draw the scene using unlit wireframes (this flag can 
             be combined with SOLID)
        - \c SMOOTH: Use smooth shading, otherwise flat (has only an effect
             when SOLID is set)
        - \c TEXTURE: Use texture mapping (has only an effect when SOLID
             is set)
        
        \param flags (\c int) Flags
        """
        
    def getRenderFlags(self):
        """Return the current render flags.
        
        \return Flags (\c int).
        \see setRenderFlags()
        """
        
    def setCamera(self, cam):
        """Set the camera to use for viewing the scene.
        
        The camera interface \a cam is used to retrieve the location
        and orientation of the camera (= inverse view transformation)
        and the projection matrix (which determines if the camera
        is an orthographic or perspective camera).
        
        \param cam (\c IViewCamera) Camera interface
        """
        
    def getCamera(self):
        """Return the currently active camera.
        
        \return Camera (\c IViewCamera)
        \see setCamera()
        """
        
    def setViewport(self, viewport):
        """Set the viewport settings.
        
        \a viewport contains the position and size (x,y,width,height) of 
        the viewing area in pixels.
        
        If the renderer is embedded this method may only be called if
        the associated rendering context (OpenGL context) is made current.
        
        \param viewport (\c IViewportSettings) Viewport settings
        """
        
    def getViewport(self):
        """Return the viewport settings.
        
        \return Viewport settings (\c IViewportSettings).
        \see setViewport()
        """
    
    def paint(self):
        """Renders the scene.
        
        The scene that was given to the open() call is rendered
        using the current camera, viewport settings and render flags.
        
        If the renderer is embedded this method may only be called if
        the associated rendering context (OpenGL context) is made current.
        
        \see setRenderFlags(), setCamera(), setViewport()
        """
        
    def pick(self, x, y):
        """Performs a picking operation and returns the picked objects.
        
        Performs a picking operation at pixel position (\a x, \a y)
        and returns all objects that were hit in depth sorted order (nearest
        first).
        
        If the renderer is embedded this method may only be called if
        the associated rendering context (OpenGL context) is made current.
                
        \param x (\c int) Pick position in pixels (X)
        \param y (\c int) Pick position in pixels (Y)
        """
        
