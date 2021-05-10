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
# $Id: odedynamics.py,v 1.8 2006/05/26 17:08:04 mbaas Exp $

## \file odedynamics.py
## Contains the ODEDynamics class.

"""This module contains a Dynamics component using the ODE rigid body
dynamics package."""

import weakref
from . import protocols
from .Interfaces import *
from .component import Component
from .globalscene import getScene
from .eventmanager import eventManager
from .cgtypes import *
from .worldobject import WorldObject
from .spheregeom import SphereGeom
from .ccylindergeom import CCylinderGeom
from .boxgeom import BoxGeom
from .trimeshgeom import TriMeshGeom
from .planegeom import PlaneGeom
from .joints import *
from .events import *
from .slots import *
from . import cmds
from . import _core
try:
    import ode
    has_ode = True
except:
    has_ode = False

#import numarray
#import numarray.linear_algebra

ODE_COLLISION = "ODECollision"

# ODECollisionEvent
class ODECollisionEvent:
    """Event object that is used as argument for collision events.
    """
    
    def __init__(self, obj1, obj2, contacts, contactproperties):
        self.obj1 = obj1
        self.obj2 = obj2
        self.contacts = contacts
        self.contactproperties = contactproperties

    def __str__(self):
        if self.obj1==None:
            n1 = "None"
        else:
            n1 = '"%s"'%self.obj1.name
        if self.obj2==None:
            n2 = "None"
        else:
            n2 = '"%s"'%self.obj2.name
        return '<ODECollisionEvent: obj1:%s obj2:%s #Contacts:%d>'%(n1, n2, len(self.contacts))

    # averageContactGeom
    def averageContactGeom(self):
        """Return the average contact position, normal and depth.
        """
        pos = vec3(0)
        normal = vec3(0)
        depth = 0.0
        for c in self.contacts:
            p,n,d,g1,g2 = c.getContactGeomParams()
            pos += vec3(p)
            normal += vec3(n)
            depth += d

        # n is not 0 (the event wouldn't have been generated if it was)
        n = len(self.contacts)
        pos /= n
        normal /= n
        depth /= n

        return pos,normal,depth
    

# ODEDynamics
class ODEDynamics(Component):
    """Dynamics component using the ODE rigid body dynamics package.
    """

    def __init__(self,
                 name="ODEDynamics",
                 gravity = 9.81,
                 substeps = 1,
                 cfm = None,
                 erp = None,
                 defaultcontactproperties = None,
                 enabled = True,
                 show_contacts = False,
                 contactmarkersize = 0.1,
                 contactnormalsize = 1.0,
                 collision_events = False,
                 use_quick_step = True,
                 auto_add = False,
                 auto_insert = True):
        """Constructor.

        \param name (\c str) Component name
        \param use_quick_step (\c bool) if True, use QuickStep method for ODE
               stepping (default); if False, use Step (slower, but more accurate)
        \param auto_add (\c bool) Automatically add the world objects to the simulation
        \param auto_insert (\c bool) Automatically add the component to the scene
        """
        Component.__init__(self, name=name, auto_insert=auto_insert)

        scene = getScene()

        self.use_quick_step = use_quick_step
        
        self.substeps = substeps
        self.collision_events = collision_events

        self.world = ode.World()
        g = -gravity*scene.up
        self.world.setGravity(g)
        if cfm!=None:
            self.world.setCFM(cfm)
        if erp!=None:
            self.world.setERP(erp)
        
#        self.world.setAutoDisableFlag(True)

        self.enabled = enabled

        self.eventmanager = eventManager()
        self.eventmanager.connect(STEP_FRAME, self)
        self.eventmanager.connect(RESET, self.reset)

        # Space object
        self.space = ode.Space()

        # Joint group for the contact joints
        self.contactgroup = ode.JointGroup()

        # A dictionary with WorldObjects as keys and ODEBodies as values
        self.body_dict = {}
        # A list of all bodies
        self.bodies = []
        # A list of all joints
        self.joints = []

        # A dictionary with contact properties
        # Key is a 2-tuple (material1, material2). Value is a
        # ODEContactProperties object.
        # The value of (mat1, mat2) should always be the same than
        # the value of (mat2, mat1).
        self.contactprops = {}

        # Default contact properties
        if defaultcontactproperties==None:
            defaultcontactproperties = ODEContactProperties()
        self.defaultcontactprops = defaultcontactproperties

        self.show_contacts = show_contacts
        self.contactmarkersize = contactmarkersize
        self.contactnormalsize = contactnormalsize

        # A list of weakrefs to body manipulators
        self.body_manips = []

        # Debug statistics (the number of contacts per simulation step)
        self.numcontacts = 0

        # Automatically add world objects
        if auto_add:
            # Add all rigid bodies first...
            for obj in scene.worldRoot().iterChilds():
                try:
                    obj = protocols.adapt(obj, IRigidBody)
                    self.add(obj)
                except NotImplementedError:
                    pass

            # Then add all joints...
            for obj in scene.worldRoot().iterChilds():
                if isinstance(obj, ODEJointBase):
                    self.add(obj)

    # remove
    def remove(self, objs):
        """Remove one or more world objects from the simulation.

        \param objs  World objects to be removed (given as names or objects)
        """
        objs = cmds.worldObjects(objs)
        for obj in objs:
            self._remove(obj)
    
    def _remove(self, obj):
        for body in self.bodies:
            if body.obj == obj:
                # after deleting all the references to the body, 
                # PyODE will remove it from simulation (in the destructor)
                self.bodies.remove(body)
                del(self.body_dict[obj])
                del(body.odegeoms)
                del(body.odebody)
                del(obj.manip)
                break
        

    # add
    def add(self, objs, categorybits=None, collidebits=None):
        """Add a world object to the simulation.

        \param objs World objects to be simulated (given as names or objects)
        """
        objs = cmds.worldObjects(objs)
        for obj in objs:
            self._add(obj, categorybits, collidebits)

    def _add(self, obj, categorybits, collidebits):
        
        if isinstance(obj, ODEJointBase):
            obj.activate(self)
            self.joints.append(obj)
            return

        try:
            obj = protocols.adapt(obj, IRigidBody)
        except NotImplementedError:
            print(('Object "%s" is not a rigid body.'%obj.name))
            return

        # Should the object be ignored?
        if not obj.dynamics:
            return

#       print "#####################",obj.name,"#########################"
#       print "Old center of gravity:",obj.cog
#       print "Old inertiatensor:"
#        print obj.inertiatensor
#       print "Old Offset transform:"
        P = obj.getOffsetTransform()
#       print P
        obj.pivot = P*obj.cog
#       print "Setting pivot to", obj.pivot
#        print "New center of gravity:", obj.cog
#        print "Intermediate inertia tensor:"
#        print obj.inertiatensor
#        I = obj.inertiatensor
#        a = numarray.array([I.getRow(0), I.getRow(1), I.getRow(2)], type=numarray.Float64)
#        vals, vecs = numarray.linear_algebra.eigenvectors(a)
#        print "Eigenvalues:",vals
        # Eigenvektoren sind ungenau! (bei cube_base_cog grosse Abweichungen)
#        b1 = vec3(vecs[0]).normalize()
#        b2 = vec3(vecs[1]).normalize()
#        b3 = vec3(vecs[2]).normalize()
#        P = obj.getOffsetTransform()        
#        Pnull = P*vec3(0,0,0)
#        b1 = P*b1 - Pnull
#        b2 = P*b2 - Pnull
#        b3 = P*b3 - Pnull        
#        I2 = mat3(b1,b2,b3)
#        print I2
#        if I2.determinant()<0:
#            I2.setColumn(0, -I2.getColumn(0))
#        print "Det of new basis",I2.determinant()
#        P = obj.getOffsetTransform()
#        P.setMat3(I2)
#        obj.setOffsetTransform(P)
#        print "New offset transform:"
#        print P
#        print "New center of gravity:",obj.cog
#        print "New inertia tensor:"
#        print obj.inertiatensor


        # If the mass of a rigid body is <= 0, ODE will crash
        # Using a more reasonable default: mass = 1
        if obj.mass <= 0:
            print(("Using default mass=1 for " + str(obj)))
            obj.mass = 1
        body = ODEBody(obj, self, categorybits=categorybits, collidebits=collidebits)
        self.bodies.append(body)
        self.body_dict[obj] = body

        obj.manip = self.createBodyManipulator(obj) # Quick access to manipulator

    # reset
    def reset(self):
        """Reset the body states.

        All bodies are reset to their position and velocity at the time
        they were added.
        
        This method is also called when the RESET event is issued.
        """
        for b in self.bodies:
            b.reset()
            b.updateObj()

        for j in self.joints:
            j.reset()

    # setContactProperties
    def setContactProperties(self, mat_tuple, props):
        """Set the contact properties of a pair of materials.

        The contact properties \a props are applied whenever an object
        with material \a mat1 collides with an object with material \a mat2.

        \param mat1 (\c Material) Material 1
        \param mat2 (\c Material) Material 2
        \param props (\c ODEContactProperties) Contact properties
        """
        mat1,mat2 = mat_tuple
        self.contactprops[(mat1,mat2)] = props  
        # Collision events are forced to appear in (mat1,mat2) order
     
    # getContactProperties
    def getContactProperties(self, matpair):
        """Return the contact properties for a material pair.

        \param matpair A 2-tuple of two Material objects
        \return Contact properties (\c ODEContactProperties)
        """
        cp = self.contactprops.get(matpair)
        if cp==None:
#            print 'ODEDynamics: Warning: No contact properties defined for material "%s" and "%s"'%(matpair[0].name, matpair[1].name)
            cp = self.defaultcontactprops
        return cp

    # createBodyManipulator
    def createBodyManipulator(self, obj):
#        return self.body_dict[obj].odebody
        bm = ODEBodyManipulator(self.body_dict[obj])
        self.body_manips.append(weakref.ref(bm))
        return bm

    # nearCallback
    def nearCallback(self, args, geom1, geom2):        
        try:            
            # Force collision event to appear in (mat1, mat2) order
            # (i.e. the order in which the material pair was defined in contactprops)
            if not (geom1.material, geom2.material) in self.contactprops:
                (geom1, geom2) = (geom2, geom1)
        except:
            pass

        # Check if the objects do collide
        contacts = ode.collide(geom1, geom2)

        # No contacts? then return immediately
        if len(contacts)==0:
            return
        
#        print len(contacts),"contacts"
        self.numcontacts += len(contacts)

        # Get the contact properties
        cp = self.getContactProperties((geom1.material, geom2.material))

        # Create a collision event?
        if self.collision_events:
            obj1 = getattr(geom1, "worldobj", None)
            obj2 = getattr(geom2, "worldobj", None)
            evt = ODECollisionEvent(obj1, obj2, contacts, cp)
            self.eventmanager.event(ODE_COLLISION, evt)

        # Create contact joints
        for c in contacts:
            if self.show_contacts:
                p,n,d,g1,g2 = c.getContactGeomParams()
                cmds.drawMarker(p, size=self.contactmarkersize, color=(1,0,0))
                cmds.drawLine(p, vec3(p)+self.contactnormalsize*vec3(n), color=(1,0.5,0.5))
#                print p,n,d
            
            # Set the properties
            cp.apply(c)
            # Create the contact joint
            j = ode.ContactJoint(self.world, self.contactgroup, c)
            b1 = geom1.getBody()
            b2 = geom2.getBody()
#            if b1==None:
#                b1=ode.environment
#            if b2==None:
#                b2=ode.environment
            j.attach(b1, b2)

    # onStepFrame
    def onStepFrame(self):
        """Callback for the StepFrame event.
        """
        if self.substeps==0 or not self.enabled:
            return

        if self.show_contacts:
            cmds.drawClear()

        # Remove dead body manipulators...
        self.body_manips = [x for x in self.body_manips if x() is not None]

        # Sim loop...
        subdt = getScene().timer().timestep/self.substeps
        for i in range(self.substeps):
            self.numcontacts = 0
            
            # Apply body manipulators
            for bmref in self.body_manips:
                bm = bmref()
                if bm is not None:
                    bm._apply()
            
            # Detect collisions and create contact joints
            self.space.collide(None, self.nearCallback)
#            print "#Contacts:",self.numcontacts

            # Simulation step
            if self.use_quick_step:
                self.world.quickStep(subdt)
            else:
                self.world.step(subdt)

            # Remove all contact joints
            self.contactgroup.empty()
            
        # Update the world objects
        for body in self.bodies:
            body.updateObj()

        # Reset body manipulators
        for bmref in self.body_manips:
            bm = bmref()
            if bm is not None:
                bm._reset()


######################################################################

# ODEContactProperties
class ODEContactProperties:

    """Contact properties.

    This class stores contact properties that are used whenever
    two objects collide. The attributes are used to initialize
    ODE contact joints.
    """    
    
    def __init__(self,
                 mode = 0,
                 mu = 0.3,
                 mu2 = None,
                 bounce = None,
                 bounce_vel = None,
                 soft_erp = None,
                 soft_cfm = None,
                 motion1 = None,
                 motion2 = None,
                 slip1 = None,
                 slip2 = None,
                 fdir1 = None):
        """Constructor.

        See the ODE manual for an explanation of the parameters.

        The flags for the mode parameter are automatically set.
        However, you can initially set the ContactApprox flags.
        """

        if mu2!=None:
            mode |= ode.ContactMu2
        else:
            mu2 = 0.0

        if bounce!=None:
            mode |= ode.ContactBounce
        else:
            bounce = 0.0

        if bounce_vel==None:
            bounce_vel = 0.0

        if soft_erp!=None:
            mode |= ode.ContactSoftERP
        else:
            soft_erp = 0.0

        if soft_cfm!=None:
            mode |= ode.ContactSoftCFM
        else:
            soft_cfm = 0.0

        if motion1!=None:
            mode |= ode.ContactMotion1
        else:
            motion1 = 0.0

        if motion2!=None:
            mode |= ode.ContactMotion2
        else:
            motion2 = 0.0

        if slip1!=None:
            mode |= ode.ContactSlip1
        else:
            slip1 = 0.0

        if slip2!=None:
            mode |= ode.ContactSlip2
        else:
            slip2 = 0.0

        if fdir1!=None:
            mode |= ode.ContactFDir1
        else:
            fdir1 = (0,0,0)

        self.mode = mode
#        print "ODEContactProps: mode =",mode
        self.mu = mu
        self.mu2 = mu2
        self.bounce = bounce
        self.bounce_vel = bounce_vel
        self.soft_erp = soft_erp
        self.soft_cfm = soft_cfm
        self.motion1 = motion1
        self.motion2 = motion2
        self.slip1 = slip1
        self.slip2 = slip2
        self.fdir1 = fdir1

    # apply
    def apply(self, contact):
        """Apply the contact properties to a contact joint."""
        contact.setMode(self.mode)
        contact.setMu(self.mu)
        contact.setMu2(self.mu2)
        contact.setBounce(self.bounce)
        contact.setBounceVel(self.bounce_vel)
        contact.setMotion1(self.motion1)
        contact.setMotion2(self.motion2)
        contact.setSlip1(self.slip1)
        contact.setSlip2(self.slip2)
        contact.setSoftERP(self.soft_erp)
        contact.setSoftCFM(self.soft_cfm)
        contact.setFDir1(self.fdir1)
                         

######################################################################

# ODEBodyManipulator
class ODEBodyManipulator(object):
    """Body manipulator class.

    This class can be used to apply external forces/torques to a body
    or to manipulate the position/orientation/velocities directly.

    You should not create instances of this class yourself. Instead, use
    the createBodyManipulator() method of the ODEDynamics component.
    """
    
    def __init__(self, body):
        """Constructor.

        \param body (\c ODEBody) The internal body object representing the
               rigid body to manipulate.
        """
        self._body = body
        self._odebody = body.odebody
        self._reset()

    # setPos
    def setPos(self, pos):
        """Set the position of the body.

        pos must be a 3-sequence of floats.
        """
        self._odebody.enable()
        self._odebody.setPosition(pos)

    # setRot
    def setRot(self, rot):
        """Set the rotation of the body.

        rot must be a mat3 containing the orientation or a list of 9 values
        in row-major order.
        """
        self._odebody.enable()
        self._odebody.setRotation(mat3(rot).toList(True))  # Now setRot really accepts a mat3 (and also a list with 9 elements)

    # setLinearVel
    def setLinearVel(self, vel):
        """Set the linear velocity.

        vel must be a 3-sequence of floats.
        """
        self._odebody.enable()
        self._odebody.setLinearVel(vel)

    # setAngularVel
    def setAngularVel(self, vel):
        """Set the angular velocity.

        vel must be a 3-sequence of floats.
        """
        self._odebody.enable()
        self._odebody.setAngularVel(vel)

    # setInitialPos
    def setInitialPos(self, pos):
        """Set the initial position.

        pos must be a 3-sequence of floats.
        """
        self._body.initial_pos = vec3(pos)

    # setInitialRot
    def setInitialRot(self, rot):
        """Set the initial orientation.

        rot must be a mat3.
        """
        self._body.initial_rot = mat3(rot)

    # setInitialLinearVel
    def setInitialLinearVel(self, vel):
        """Set the initial linear velocity.

        vel must be a 3-sequence of floats.
        """
        self._odebody.initial_linearvel = vec3(vel)

    # setInitialAngularVel
    def setInitialAngularVel(self, vel):
        """Set the initial angular velocity.

        vel must be a 3-sequence of floats.
        """
        self._odebody.initial_angularvel = vec3(vel)

    # addForce
    def addForce(self, force, relforce=False, pos=None, relpos=False):
        """Apply an external force to a rigid body.
        
        Add an external force to the current force vector. force is a
        vector containing the force to apply. If relforce is True the
        force is interpreted in local object space, otherwise it is
        assumed to be given in global world space. By default, the
        force is applied at the center of gravity. You can also pass a
        different position in the pos argument which must describe a
        point in space. relpos determines if the point is given in
        object space or world space (default).
        """
        
        R = None
        force = vec3(force)
        if relforce:
            R = mat3(self._odebody.getRotation())
            force = R*force
        # Is a position given? then add the corresponding torque
        if pos!=None:
            pos = vec3(pos)
            bodyorig = vec3(self._odebody.getPosition())
            if relpos:
                if R==None:
                    R = mat3(self._odebody.getRotation())
                pos = R*pos + bodyorig

            self._torque += (pos-bodyorig).cross(force)
            self._torque_flag = True
            
        self._force += vec3(force)
        self._force_flag = True
                    
    # addTorque
    def addTorque(self, torque, reltorque=False):
        """Apply an external torque to a rigid body.
        
        Add an external torque to the current torque vector. torque is
        a vector containing the torque to apply. reltorque determines
        if the torque vector is given in object space or world space
        (default).
        """
        
        torque = vec3(torque)
        if reltorque:
            R = mat3(self._odebody.getRotation())
            torque = R*torque
            
        self._torque += torque
        self._torque_flag = True

    # _apply
    def _apply(self):
        """Apply the stored force/torque.

        This method is called by the ODEDynamics object during the simulation
        step (once for every sub step).
        """
        
        if self._force_flag or self._torque_flag:
            self._odebody.enable()
            
        if self._force_flag:
            self._odebody.addForce(self._force)
        if self._torque_flag:
            self._odebody.addTorque(self._torque)

    # _reset
    def _reset(self):
        """Reset the stored force/torque.

        This method is called by the ODEDynamics object at the end of one
        entire simulation step.
        """
        self._force = vec3(0)
        self._torque = vec3(0)
        
        self._force_flag = False
        self._torque_flag = False

    ## protected:
        
    # "body" property...
    
    def _getBody(self):
        """Return the current body (WorldObject).

        This method is used for retrieving the \a body property.

        \return Rigid body
        """
        return self._body.obj

    body = property(_getBody, None, None, "Rigid body")

    # "odebody" property...
    # Someone may think this method would return an ODEBody, not an ode.Body
    # That's why I changed the case.
    def _get_odeBody(self):
        """Return the current ODE body (type ode.Body from PyODE). 

        This method is used for retrieving the \a odebody property.

        \return ODE body
        """
        return self._odebody

    odebody = property(_get_odeBody, None, None, "ODE body")

    # odegeoms property
    
    def _get_odeGeoms(self):
        """Return the current ODE geom list.

        This method is used for retrieving the \a odegeoms property.

        \return ODE geom list
        """
        return self._body.odegeoms

    odegeoms = property(_get_odeGeoms, None, None, "ODE geoms")

    
    
######################################################################

# ODEBody
class ODEBody:
    """This class is the interface between ODE bodies and world objects.

    It encapsulates the ODE body and the ODE geom objects.

    """
    def __init__(self, obj, dyncomp, categorybits=None, collidebits=None):
        """Constructor.

        \param obj (\c WorldObject) World object
        \param dyncomp (\c ODEDynamics) Dynamics component
        """
        
        # Store the world object
        self.obj = obj

        self.dyncomp = dyncomp
        world = self.dyncomp.world
        space = self.dyncomp.space

        # Store initial position/orientation
        self.initial_pos = obj.pos
        self.initial_rot = obj.rot
        self.initial_linearvel = obj.linearvel #vec3(0,0,0)
        self.initial_angularvel = obj.angularvel #vec3(0,0,0)

        # The ODE body
        self.odebody = None
        # A list of ODE geoms 
        self.odegeoms = None
    
        # Create ODE collision geoms
        self.odegeoms = self.createGeoms(obj, space)

        # Set the category and collide bits
        for geom in self.odegeoms:
            if categorybits!=None:
                geom.setCategoryBits(categorybits)
            if collidebits!=None:
                geom.setCollideBits(collidebits)

        # Create an ODE body
        if not obj.static:
            self.odebody = self.createBody(obj, world)
            for g in self.odegeoms:
                g.setBody(self.odebody)

        # Store the material in the ODE geoms so that the contact
        # properties can be determined in the near callback.
        # Also store the corresponding WorldObject
        for g in self.odegeoms:
            g.material = obj.getMaterial()
            g.worldobj = obj

        # Initialize the ODE body/geoms
        self.reset()

        # Create the notification forwarders and make connections
        self.initialization = True
        # Forward changes to the "static" slot
        self._static_forwarder = NotificationForwarder(self.onStaticChanged)
        self.obj.static_slot.addDependent(self._static_forwarder)
        # Forward changes to the "mass" slot
        self._mass_forwarder = NotificationForwarder(self.onMassChanged)
        self.obj.mass_slot.addDependent(self._mass_forwarder)
        self.initialization = False

    # updateObj
    def updateObj(self):
        """Update the transformation of the world object.

        The transformation from the ODE body is copied to the corresponding
        world object. The linearvel and angularvel attributes are
        updated as well.
        """
        if self.obj.static:
            return

#        if type(self.odegeom)==ode.GeomTriMesh:
#            self.odegeom.clearTCCache()

        # Get the ODE transformation
        pos = self.odebody.getPosition()
        rot = self.odebody.getRotation()
        # Create the transformation matrix for the world object
        M = mat4().translation(vec3(pos))
        M.setMat3(mat3(rot))
        # Set the transformation on the world object
        self.obj.transform = M

        self.obj.linearvel = vec3(self.odebody.getLinearVel())
        self.obj.angularvel = vec3(self.odebody.getAngularVel())

    ## protected:

    # onStaticChanged
    def onStaticChanged(self):
        if self.initialization:
            return
        
        if self.obj.static:
            if self.odebody!=None:
                for g in self.odegeoms:
                    g.setBody(None)
            self.odebody = None
        else:
            self.odebody = self.createBody(self.obj, self.dyncomp.world)
            for g in self.odegeoms:
                g.setBody(self.odebody)
            self.odebody.setPosition(self.obj.pos)
            self.odebody.setRotation(self.obj.rot.toList(rowmajor=True))
            self.odebody.setLinearVel((0,0,0))
            self.odebody.setAngularVel((0,0,0))            
            
    # onMassChanged
    def onMassChanged(self):
        if self.initialization:
            return
        
        print(("Mass changed to",self.obj.mass))

    # reset
    def reset(self):
        """Restore the initial state of the body."""
        if self.obj.static:
            for g in self.odegeoms:
                if g.placeable():
                    g.setPosition(self.initial_pos)
                    g.setRotation(self.initial_rot.toList(rowmajor=True))
        else:
            self.odebody.setPosition(self.initial_pos)
            self.odebody.setRotation(self.initial_rot.toList(rowmajor=True))
            self.odebody.setLinearVel(self.initial_linearvel)
            self.odebody.setAngularVel(self.initial_angularvel)

    # createBody
    def createBody(self, obj, world):
        """Create an ODE body.

        Creates an ODE body and initializes the mass properties.
        """
        # Create an ODE Mass object
        M = ode.Mass()
        m = obj.totalmass
        if m==0:
            print(('WARNING: Object "%s" has a mass of zero.'%obj.name))
        cog = obj.cog
        I = obj.inertiatensor
#        print '---Rigid body "%s"--------'%obj.name
#        print 'cog:',cog
#        print I
        I = mat3(I)
        M.setParameters(m, cog.x, cog.y, cog.z, I[0,0], I[1,1], I[2,2], I[0,1], I[0,2], I[1,2])

        # Create an ODE body
        body = ode.Body(world)
        body.setMass(M)

        return body

    # createGeoms
    def createGeoms(self, obj, space):
        """Create all the geoms for one world object including the children.

        The generated geoms will be encapsulated in GeomTransforms
        so that all of them are defined with respect to the pivot coordinate
        system of \a obj.

        \param obj (\c WorldObject) Top level world object
        \param space (\c ode.Space) ODE Space object
        \return List of ODE geom objects
        """

        # Plane? This is a special case because ODE planes are not placeable
        if isinstance(obj.geom, PlaneGeom):
            if not obj.static:
                raise RuntimeError("Planes can only be used as static bodies")
            L = obj.localTransform()
            n = L*vec3(0,0,1) - L*vec3(0)
            n = n.normalize()
            d = n*obj.pos
            odegeom = ode.GeomPlane(space, n, d)
            return [odegeom]

        res = []

        # Create all ODE geoms in the subtree starting at obj
        # Each geom will have a position/rotation that moves it to
        # the local coordinate system of obj.
        geoms = self._createLGeoms(obj, mat4(1))

        # Apply the offset transform and encapsulate the geoms inside
        # geom transforms...
        P = obj.getOffsetTransform()
        Pinv = P.inverse()
        pos,rot,scale = Pinv.decompose()
        if scale!=vec3(1,1,1):
            print('WARNING: ODEDynamics: Scaled geometries are not supported')
        res = []
        for g in geoms:
            M = mat4(1).translation(vec3(g.getPosition()))
            M.setMat3(mat3(g.getRotation()))
            M *= Pinv;
            p4 = M.getColumn(3)
            g.setPosition((p4.x, p4.y, p4.z))
            # row major or column major?
            g.setRotation(M.getMat3().toList(rowmajor=True))
        
            gt = ode.GeomTransform(space)
            gt.setGeom(g)
            res.append(gt)

        return res
        
    # _createLGeoms
    def _createLGeoms(self, obj, M):
        """Create ODE geom objects for a world object.

        Create the ODE geoms for \a obj (and its children).
        The transformation \a M is applied to each geom.

        The returned geoms are not assigned to a Body yet.

        \param obj (\c WorldObject) The world object
        \param M (\c mat4) Transformation
        \return List of ODE geoms
        """
        res = []
        # Create an ODE geom for the object itself
        geom = obj.geom
        if geom!=None:
            res.append(self._createGeom(geom, M))

        # Create ODE geoms for the children
        for child in obj.iterChilds():
            if not getattr(child, "dynamics", False):
                continue
            L = child.localTransform()
            res += self._createLGeoms(child, M*L)

        return res
        
    # _createGeom
    def _createGeom(self, geom, M):
        """Create an ODE geom object for the given cgkit geom.

        Create an ODE geom of the appropriate type (depening on \a geom)
        and apply the transformation \a M to it (currently, M may only
        contain a translation and a rotation, otherwise a warning is
        printed).

        The returned geom is not assigned to a Body yet.

        \param geom (\c GeomObject) cgkit geom object
        \param M (\c mat4) Transformation that's applied to the geom
        \return ODE geom
        """

        # Create the raw geom object that's defined in the local
        # coordinate system L of the world object this geom belongs to.

        # Sphere geom?
        if isinstance(geom, SphereGeom):
            odegeom = ode.GeomSphere(None, geom.radius)
        # CCylinder geom?
        elif isinstance(geom, CCylinderGeom):
            odegeom = ode.GeomCCylinder(None, geom.radius, geom.length)
        # Box geom?
        elif isinstance(geom, BoxGeom):
            odegeom = ode.GeomBox(None, (geom.lx, geom.ly, geom.lz))
        # TriMesh geom?
        elif isinstance(geom, TriMeshGeom):
            verts = []
            for i in range(geom.verts.size()):
                verts.append(geom.verts.getValue(i))
#                print "V",i,verts[i]
            faces = []
            for i in range(geom.faces.size()):
                f = geom.faces.getValue(i)
                faces.append(f)
#                print "F",i,faces[i]
                ia,ib,ic = faces[i]
                a = verts[ia]
                b = verts[ib]
                c = verts[ic]
                sa = ((b-a).cross(c-a)).length()
#                if sa<0.0001:
#                    print "*****KLEIN*****",sa

            tmd = ode.TriMeshData()
            tmd.build(verts, faces)            
            odegeom = ode.GeomTriMesh(tmd, None)
        # Plane geom?
        elif isinstance(geom, PlaneGeom):
            L = obj.localTransform()
            n = L*vec3(0,0,1) - L*vec3(0)
            n = n.normalize()
            d = n*obj.pos
            odegeom = ode.GeomPlane(space, n, d)
        # Unknown geometry
        else:
            raise ValueError('WARNING: ODEDynamics: Cannot determine collision geometry of object "%s".'%geom.name)
#            print 'WARNING: ODEDynamics: Cannot determine collision geometry of object "%s". Using bounding box instead.'%obj.name
#            bmin, bmax = obj.boundingBox().getBounds()
#            s = bmax-bmin
#            odegeom = ode.GeomBox(None, s)
#            pos,rot,scale = P.inverse().decompose()
#            odegeom.setPosition(pos + 0.5*(bmax+bmin))
#            odegeomtransform = ode.GeomTransform(space)
#            odegeomtransform.setGeom(odegeom)
#            return odegeomtransform

        # Displace the geom by M
        pos,rot,scale = M.decompose()
        if scale!=vec3(1,1,1):
            print('WARNING: ODEDynamics: Scaled geometries are not supported')
        odegeom.setPosition(pos)
        # row major or column major?
        odegeom.setRotation(rot.getMat3().toList(rowmajor=True))
        
        return odegeom

######################################################################

# ODEJointBase
class ODEJointBase(WorldObject):

    protocols.advise(instancesProvide=[ISceneItem])

    exec(slotPropertyCode("lostop"))
    exec(slotPropertyCode("histop"))
    exec(slotPropertyCode("motorvel"))
    exec(slotPropertyCode("motorfmax"))
    exec(slotPropertyCode("fudgefactor"))
    exec(slotPropertyCode("bounce"))
    exec(slotPropertyCode("cfm"))
    exec(slotPropertyCode("stoperp"))
    exec(slotPropertyCode("stopcfm"))

    def __init__(self,
                 name = "",
                 body1 = None,
                 body2 = None,
                 **params):
        WorldObject.__init__(self, name=name, **params)
        
        # ODEDynamics component
        self.odedynamics = None

        # The corresponding ODE joint
        self.odejoint = None

        self.body1 = body1
        self.body2 = body2

    # attach
    def attach(self, body1, body2=None):
        self.body1 = body1
        self.body2 = body2
        if self.odejoint!=None:
            if body1==None:
                b1 = ode.environment
            else:
                b1 = self.odedynamics.body_dict[body1].odebody
            
            if body2==None:
                b2 = ode.environment
            else:
                b2 = self.odedynamics.body_dict[body2].odebody

            self.odejoint.attach(b1, b2)

    # activate
    def activate(self, odedynamics):
        """
        this method is called by the ODEDynamics componenent
        """
        if odedynamics==self.odedynamics:
            return
        if self.odedynamics!=None:
            print(('Warning: Joint "%s" is already in use'%self.name))
        self.odedynamics = odedynamics
        self._createODEjoint()
        self._initODEjoint()

    # reset
    def reset(self):
        """Reset method.

        This may only be called after the bodies have been reset.
        """
        pass

    def _createSlots(self):
        
        self._createSlot("lostop", None, "-ode.Infinity", "onLoStopChanged")
        self._createSlot("histop", None, "ode.Infinity", "onHiStopChanged")
        self._createSlot("motorvel", None, "0.0", "onMotorVelChanged")
        self._createSlot("motorfmax", None, "0.0", "onMotorFMaxChanged")
        self._createSlot("fudgefactor", None, "1.0", "onFudgeFactorChanged")
        self._createSlot("bounce", None, "0.1", "onBounceChanged")
        self._createSlot("cfm", None, 1E-5, "onCFMChanged")
        self._createSlot("stoperp", None, 0.2, "onStopERPChanged")
        self._createSlot("stopcfm", None, 1E-5, "onStopCFMChanged")
        
    def _createSlot(self, name, nr, default, callback):
        """Create one joint param slot.

        This method creates a DoubleSlot called \a name (with additional
        number, if \a nr is not None). \a default is the initial default
        value (as string) and \a callback the name of the callback method.

        \param name (\c str) Parameter name
        \param nr (\c int) Parameter set number or None
        \param default (\c str) Default value as string
        \param callback (\c str) Name of the callback method
        """
        if nr!=None:
            name += str(nr)
        # Create the slot
        s = "self.%s_slot = DoubleSlot(%s)"%(name, default)
        exec(s)
        # Add the slot to the component
        s = "self.addSlot('%s', self.%s_slot)"%(name, name)
        exec(s)
        # Create the forwarder
        s = "self._%s_forwarder = NotificationForwarder(self.%s)"%(name, callback)
        exec(s)
        s = "self.%s_slot.addDependent(self._%s_forwarder)"%(name, name)
        exec(s)
       

    # _createODEjoint
    def _createODEjoint(self):
        """Create the corresponding ODE joint.

        This may only be called if the ODEDynamics component has been set.

        This method has to create an ODE joint object and store it
        in self.odejoint.
        """
        pass

    # _initODEjoint
    def _initODEjoint(self):
        """Initialize the ODE joint.
        """
        # This will call the attach() method of the ODE joint
        self.attach(self.body1, self.body2)
        
        try:
            print(("***",self.odejoint.getParam(ode.ParamStopCFM)))
            
            self.odejoint.setParam(ode.ParamFMax, self.motorfmax)
            self.odejoint.setParam(ode.ParamVel, self.motorvel)
            self.odejoint.setParam(ode.ParamLoStop, self.lostop)
            self.odejoint.setParam(ode.ParamHiStop, self.histop)
            self.odejoint.setParam(ode.ParamFudgeFactor, self.fudgefactor)
            self.odejoint.setParam(ode.ParamBounce, self.bounce)
            self.odejoint.setParam(ode.ParamCFM, self.cfm)
            self.odejoint.setParam(ode.ParamStopERP, self.stoperp)
            self.odejoint.setParam(ode.ParamStopCFM, self.stopcfm)
        except AttributeError:
            pass  # not all joints have these attributes

    def onLoStopChanged(self):
#        print "Lostop has been changed to",self.lostop_slot.getValue()
        if self.odejoint!=None:
            self.odejoint.setParam(ode.ParamLoStop, self.lostop_slot.getValue())

    def onHiStopChanged(self):
#        print "Histop has been changed to",self.histop_slot.getValue()
        if self.odejoint!=None:
            self.odejoint.setParam(ode.ParamHiStop, self.histop_slot.getValue())

    def onMotorVelChanged(self):
#        print "Motor velocity has been changed to",self.motorvel_slot.getValue()
        if self.odejoint!=None:
            self.odejoint.setParam(ode.ParamVel, self.motorvel)


    def onMotorFMaxChanged(self):
#        print "Motor FMax has been changed to",self.motorfmax_slot.getValue()
        if self.odejoint!=None:
            self.odejoint.setParam(ode.ParamFMax, self.motorfmax)

    def onFudgeFactorChanged(self):
#        print "Fudgefactor has been changed to",self.fudgefactor_slot.getValue()
        if self.odejoint!=None:
            self.odejoint.setParam(ode.ParamFudgeFactor, self.fudgefactor_slot.getValue())

    def onBounceChanged(self):
#        print "Bounce has been changed to",self.bounce_slot.getValue()
        if self.odejoint!=None:
            self.odejoint.setParam(ode.ParamBounce, self.bounce_slot.getValue())


    def onMotorVel2Changed(self):
#        print "Motor2 velocity has been changed to",self.motorvel2_slot.getValue()
        if self.odejoint!=None:
            self.odejoint.setParam(ode.ParamVel2, self.motorvel2)


    def onMotorFMax2Changed(self):
#        print "Motor2 FMax has been changed to",self.motorfmax2_slot.getValue()
        if self.odejoint!=None:
            self.odejoint.setParam(ode.ParamFMax2, self.motorfmax2)

    def onCFMChanged(self):
#        print "CFM has been changed to",self.cfm_slot.getValue()
        if self.odejoint!=None:
            self.odejoint.setParam(ode.ParamCFM, self.cfm)

    def onStopCFMChanged(self):
#        print "Stop CFM has been changed to",self.stopcfm_slot.getValue()
        if self.odejoint!=None:
            self.odejoint.setParam(ode.ParamStopCFM, self.stopcfm)

    def onStopERPChanged(self):
#        print "Stop ERP has been changed to",self.stoperp_slot.getValue()
        if self.odejoint!=None:
            self.odejoint.setParam(ode.ParamStopERP, self.stoperp)


# BallJoint
class ODEBallJoint(ODEJointBase):

    def __init__(self,
                 name = "ODEBallJoint",
                 body1 = None,
                 body2 = None,
                 **params):
        ODEJointBase.__init__(self,
                              name=name, body1=body1, body2=body2,
                              **params)

        self._createSlots()

    # _createODEjoint
    def _createODEjoint(self):
        # Create the ODE joint
        self.odejoint = ode.BallJoint(self.odedynamics.world)

    # _initODEjoint
    def _initODEjoint(self):
        ODEJointBase._initODEjoint(self)

        W = self.worldtransform
        p = W[3]
        self.odejoint.setAnchor((p.x, p.y, p.z))


# HingeJoint
class ODEHingeJoint(ODEJointBase):
    """
    the rotation axis is the local x axis
    """

    exec(slotPropertyCode("angle"))

    def __init__(self,
                 name = "ODEHingeJoint",
                 body1 = None,
                 body2 = None,
                 **params):
        ODEJointBase.__init__(self, name=name, body1=body1, body2=body2,
                              **params)

        self._createSlots()

        self.angle_slot = ProceduralDoubleSlot(self.computeAngle)
        self.addSlot("angle", self.angle_slot)
        getScene().timer().time_slot.addDependent(self.angle_slot)

#???
#    def reset(self):
#        if self.odejoint!=None:
#            print "ANGLE:",self.odejoint.getAngle(), self.angle
#            self.odejoint.setAnchor(self.odejoint.getAnchor())
        
    # computeAngle
    def computeAngle(self):
        """Return the current angle.

        This method is used as procedure for the angle_slot.
        """
        if self.odejoint==None:
            return 0.0
        else:
            return self.odejoint.getAngle()            

    # _createODEjoint
    def _createODEjoint(self):
        # Create the ODE joint
        self.odejoint = ode.HingeJoint(self.odedynamics.world)

    # _initODEjoint
    def _initODEjoint(self):
        ODEJointBase._initODEjoint(self)

        W = self.worldtransform
        p = W[3]
        self.odejoint.setAnchor((p.x, p.y, p.z))
        a = -W[0]
#        print "AXIS:",a,self.name
        self.odejoint.setAxis((a.x, a.y, a.z))
        self.odejoint.setParam(ode.ParamFMax, self.motorfmax)
        self.odejoint.setParam(ode.ParamVel, self.motorvel)



# FixedJoint
class ODEFixedJoint(ODEJointBase):
    """
    Fixed Joint: Glues two bodies together.
    Not recommended by ODE manual, but useful when a solid body has different contact properties on different sides.
    """

    def __init__(self,
                 name = "ODEFixedJoint",
                 body1 = None,
                 body2 = None,
                 **params):
        ODEJointBase.__init__(self, name=name, body1=body1, body2=body2,
                              **params)

        self._createSlots()

    # _createODEjoint
    def _createODEjoint(self):
        # Create the ODE joint
        self.odejoint = ode.FixedJoint(self.odedynamics.world)

    # _initODEjoint
    def _initODEjoint(self):
        ODEJointBase._initODEjoint(self)
        self.odejoint.setFixed()




# SliderJoint
class ODESliderJoint(ODEJointBase):
    """
    the slider axis is the local x axis
    """

    exec(slotPropertyCode("position"))

    def __init__(self,
                 name = "ODEHingeJoint",
                 body1 = None,
                 body2 = None,
                 **params):
        ODEJointBase.__init__(self, name=name, body1=body1, body2=body2,
                              **params)

        self._createSlots()

        self.position_slot = ProceduralDoubleSlot(self.computePosition)
        self.addSlot("position", self.position_slot)
        getScene().timer().time_slot.addDependent(self.position_slot)

    # computePosition
    def computePosition(self):
        """Return the current slider position.

        This method is used as procedure for the position_slot.
        """
        if self.odejoint==None:
            return 0.0
        else:
            return self.odejoint.getPosition()            

    # _createODEjoint
    def _createODEjoint(self):
        # Create the ODE joint
        self.odejoint = ode.SliderJoint(self.odedynamics.world)

    # _initODEjoint
    def _initODEjoint(self):
        ODEJointBase._initODEjoint(self)

        W = self.worldtransform
        a = W[0]
        self.odejoint.setAxis((a.x, a.y, a.z))
        self.odejoint.setParam(ode.ParamFMax, self.motorfmax)
        self.odejoint.setParam(ode.ParamVel, self.motorvel)


# Hinge2Joint
class ODEHinge2Joint(ODEJointBase):
    """
    axis1 is the local z axis
    """

    exec(slotPropertyCode("motorfmax"))
    exec(slotPropertyCode("motorvel"))
    exec(slotPropertyCode("suspensionerp"))
    exec(slotPropertyCode("suspensioncfm"))
    exec(slotPropertyCode("motorfmax2"))
    exec(slotPropertyCode("motorvel2"))

    def __init__(self,
                 name = "ODEHingeJoint",
                 body1 = None,
                 body2 = None,
                 **params):
        ODEJointBase.__init__(self, name=name, body1=body1, body2=body2,
                              **params)

        self._createSlots()
        self._createSlot("suspensionerp", None, "0.2", "onSuspensionERPChanged")
        self._createSlot("suspensioncfm", None, "1E-7", "onSuspensionCFMChanged")

        self._createSlot("motorvel2", None, "0.0", "onMotorVel2Changed")
        self._createSlot("motorfmax2", None, "0.0", "onMotorFMax2Changed")


        # default axis2 is the local y axis
        self.axis2_local = vec4(0,1,0,0)

    # _createODEjoint
    def _createODEjoint(self):
        # Create the ODE joint
        self.odejoint = ode.Hinge2Joint(self.odedynamics.world)

    # _initODEjoint
    def _initODEjoint(self):
        ODEJointBase._initODEjoint(self)

        W = self.worldtransform
        p = W[3]
        self.odejoint.setAnchor((p.x, p.y, p.z))
        a = W[2]
        self.odejoint.setAxis1((a.x, a.y, a.z))
        a = W*self.axis2_local
        self.odejoint.setAxis2((a.x, a.y, a.z))
        
        self.odejoint.setParam(ode.ParamFMax, self.motorfmax)
        self.odejoint.setParam(ode.ParamVel, self.motorvel)

    def onSuspensionERPChanged(self):
        print(("susp. erp has been changed to",self.suspensionerp_slot.getValue()))
        if self.odejoint!=None:
            self.odejoint.setParam(ode.ParamSuspensionERP,
                                   self.suspensionerp_slot.getValue())

    def onSuspensionCFMChanged(self):
        print(("susp. cfm has been changed to",self.suspensioncfm_slot.getValue()))
        if self.odejoint!=None:
            self.odejoint.setParam(ode.ParamSuspensionCFM,
                                   self.suspensioncfm_slot.getValue())


# UniversalJoint
class ODEUniversalJoint(ODEJointBase):
    """

    """

    exec(slotPropertyCode("angle"))

    def __init__(self,
                 name = "ODEUniversalJoint",
                 body1 = None,
                 body2 = None,
                 **params):
        ODEJointBase.__init__(self, name=name, body1=body1, body2=body2,
                              **params)

        self._createSlots()

#        self.angle_slot = ProceduralDoubleSlot(self.computeAngle)
#        self.addSlot("angle", self.angle_slot)
#        getScene().timer().time_slot.addDependent(self.angle_slot)


    # _createODEjoint
    def _createODEjoint(self):
        # Create the ODE joint
        self.odejoint = ode.UniversalJoint(self.odedynamics.world)

    # _initODEjoint
    def _initODEjoint(self):
        ODEJointBase._initODEjoint(self)

        W = self.worldtransform
        p = W[3]
        self.odejoint.setAnchor((p.x, p.y, p.z))
        a = W[0]
        self.odejoint.setAxis1((a.x, a.y, a.z))
        a = W[1]
        self.odejoint.setAxis2((a.x, a.y, a.z))
        self.odejoint.setParam(ode.ParamFMax, self.motorfmax)
        self.odejoint.setParam(ode.ParamVel, self.motorvel)
