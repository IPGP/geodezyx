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
# $Id: tunnel.py,v 1.2 2005/08/08 18:05:34 mbaas Exp $

## \file tunnel.py
## Contains the Tunnel class.

from .cgtypes import *
from . import component
from .slots import *
from . import _core
import socket, threading, struct

# _ImmediateForwarder
class _ImmediateForwarder(NotificationForwarder):
    """Internal class for the Tunnel component.

    This forwarder immediately sends the new value to the server part
    of the tunnel.
    """
    def __init__(self, tunnel, slot, msgid, idx):
        NotificationForwarder.__init__(self, self.sendValue)
        self.tunnel = tunnel
        self.slot = slot
        self.msgid = msgid
        self.idx = idx

    def sendValue(self):
        t = self.tunnel
        msg = t.encodeValueMessage(self.msgid, self.idx, self.slot)
        t.send(msg)

# _GatherForwarder
class _GatherForwarder(NotificationForwarder):
    """Internal class for the Tunnel component.

    This forwarder marks the slot as changed and sends a combined
    message if it's the sender slot (the last slot in the tunnel).
    """
    def __init__(self, tunnel, idx, sender):
        NotificationForwarder.__init__(self, self.markSlot)
        self.tunnel = tunnel
        self.idx = idx
        self.sender = sender

    def markSlot(self):
        self.tunnel.markAsChanged(self.idx)
        if self.sender:
            self.tunnel.sendBulk()

                 

# Tunnel
class Tunnel(component.Component):
    """%Tunnel component which can connect slots on different machines.

    A tunnel is a component that has an arbitrary number of input slots
    and a corresponding output slot for each input slot.
    The value of the input slot is simply passed to the output slot.
    The speciality of the tunnel is that both ends may reside on two
    different machines.

    \image html tunnel.png

    An instance of the tunnel component either represents the input part
    (client) or the output part (server). Value changes on the client are
    then propagated to the server. When creating either side of the tunnel
    you have to define the name and type of the slots that should be created.
    The slot types of the client and the server should always match.

    Note: ArraySlots are currently not supported.

    Example:

    \code
    # On machine A
    t = Tunnel(
       slots = [("pos", "Vec3Slot()"),
                ("rot", "Mat3Slot()")],
       host  = "<name or IP of machine B>",
    )

    # Do connections between other components and t.pos_slot/t.rot_slot
    # or set the values directly on the tunnel
    \endcode

    \code
    # On machine B
    t = Tunnel(
       server = True,
       slots = [("pos", "Vec3Slot()"),
                ("rot", "Mat3Slot()")],
    )
    
    # Make connections between t.pos_slot/t.rot_slot and other components
    # or retrieve the values directly from the tunnel
    \endcode

    \par Protocol

    The client sends its messages as UDP datagrams to the server.
    One packet may contain several individual messages that are just
    concatenated. The first two bytes of each message is the message
    ID. The remaining part depends on the ID. The byte order is always
    in big endian.

    - ID 0: Init message (must be the only message in a packet)
    - ID 1: Binary int value
    - ID 2: Binary double value
    - ID 3: Binary bool value (actually short int)
    - ID 4: Binary vec3 value (3 doubles)
    - ID 5: Binary vec4 value (4 doubles)
    - ID 6: Binary mat3 value (9 doubles)
    - ID 7: Binary mat4 value (16 doubles)
    - ID 8: Binary quat value (4 doubles)

    All the value messages are composed of the slot index number (short int),
    followed by the binary value.    
    """

    def __init__(self,
                 name = "Tunnel",
                 server = False,
                 slots = None,
                 port = 64738,
                 host = "localhost",
                 gather_messages = False,
                 verbose = False,
                 auto_insert = True):
        """Constructor.

        The argument \a slot specifies the slots to create on the tunnel.
        It must be a list of 2-tuples containing two strings. The first
        string is the name of the attribute and the second contains the
        type and initializer of the slot (in Python syntax). The order and
        types of the slots on the client and server should always match.

        The arguments \a host and \a gather_messages are only meaningful
        on the client side. If \a gather_messages is True, then value changes
        on a slot are not immediately sent but only when the last slot
        receives a new value. You can use this option if you know that
        all slot values will change at the same time. In this case it's more
        efficient to change all values and then send only one message
        containing all new values instead of sending individual messages for
        every slot.

        \param name (\c str) Component name
        \param server (\c bool) True if this is the server part
        \param slots (\c list) A list of slot definitions. Each definition is a tuple (attribute name, slot class and initial arguments).
        \param port (\c int) Port number to use for the data transfer
        \param host (\c str) Host where the tunnel server is running (name or IP address). Client only.
        \param gather_messages (\c bool) True if the messages should be combined. Client only.
        \param verbose (\c bool) True if debug messages should be printed
        \param auto_insert (\c bool) True if the component should be inserted into the scene automatically
        """
        component.Component.__init__(self, name=name, auto_insert=auto_insert)

        if slots==None:
            slots = []

        self.server = server
        self.host = host
        self.port = port
        self.sock = None
        self.slots = slots
        self.gather_messages = gather_messages
        self.verbose = verbose

        if self.server:
            self.host = socket.gethostname()
            
        # A list with tuples (slot, id)
        self.slotobjs = []
        # A dictionary with slot indices as keys. If an index is present
        # in a dictionary then the corresponding slot was modified and the
        # value has to be sent to the server
        self.changed = {}

        self.forwarders = []

        # Initialize client/server...
        if server:
            self.initServer()
        else:
            self.initClient()

    def __str__(self):
        if self.server:
            return 'Tunnel server "%s" listening at %s:%d. %d slots.'%(self.name, self.host, self.port, len(self.slots))
        else:
            return 'Tunnel client "%s". %d slots. Server at %s:%d.'%(self.name, len(self.slots), self.host, self.port)

    ## protected:

    # initClient
    def initClient(self):
        """Initialization method for the client.
        """
        # Create socket object for sending UDP datagrams
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

        # Create slot attributes...
        self.createClientSlotAttribs()

        # Send init message
        s = ["%s:%s"%x for x in self.slots]
        s = ";".join(s)
        self.send(struct.pack(">H", 0)+s)

    # initServer
    def initServer(self):
        """Initialization method for the server.
        """
        # Create slot attributes...
        self.createServerSlotAttribs()
        
        # Open UDP socket
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        ipaddr = socket.gethostbyname(self.host)
        if self.verbose:
            print(("Open UDP port %d on %s (%s)"%(self.port, ipaddr, self.host)))
        self.sock.bind((ipaddr, self.port))

        # Split thread with the server loop
        self.serverthread = threading.Thread(name="Tunnel-Server",
                                             target=self.serverLoop)
        self.serverthread.setDaemon(True)
        self.serverthread.start()    

    # serverLoop
    def serverLoop(self):
        """Server loop.

        This method reads messages from the previously created socket
        and process them.
        The method is run in its own thread (because reading from the
        socket is a blocking operation).
        """

        if self.verbose:
            print(("Tunnel server running at port",self.port))
        while 1:
            # Read a packet (blocking)
            rawdata, addr = self.sock.recvfrom(5000)
            if self.verbose:
                print("---Tunnel messages--")
            # Process all messages...
            while rawdata!="":
                # Get the id of the first message...
                try:
                    id = struct.unpack(">H", rawdata[:2])[0]
                except:
                    print("Error: Invalid UDP packet (no valid message id).")
                    rawdata = ""
                    continue

                # Process the message...
                msgdata = rawdata[2:]
                if id==0:
                    if self.verbose:
                        print("Init")
                    f = msgdata.split(";")
                    localslots = [x[1] for x in list(self.iterSlotDefs(self.slots))]
                    remoteslots = [tuple(x.split(":")) for x in f]
                    remoteslots = [x[1] for x in list(self.iterSlotDefs(remoteslots))]
                    mismatch = True
                    if len(localslots)==len(remoteslots):
                        mismatch = localslots!=remoteslots
                    if mismatch:
                        print(("%s: The types of the remote and local slots don't match."%self.name))
                    rawdata = ""
                else:
                    # Decode the slot index and the value to set...
                    try:
                        n,idx,v = self.decodeValueMessage(id, msgdata)
                    except struct.error as e:
                        print(e)
                        rawdata = ""
                        continue
                    except ValueError as e:
                        print(e)
                        rawdata = ""
                        continue
                    # Prepare the next message...
                    rawdata = msgdata[n:]
                    if self.verbose:
                        print(("Message id: %d - Slot:%d Value:%s"%(id,idx,v)))
                    # Get the slot that will receive the new value...
                    try:
                        slot,id = self.slotobjs[idx]
                    except IndexError as e:
                        print(("Error: Invalid slot index (%d)"%idx))
                    # Set the new value...
                    try:
                        slot.setValue(v)
                    except:
                        print(("Error: Could not assign value",v,"to slot of type",slot.typeName()))

        self.sock.close()
        if self.verbose:
            print(("Tunnel server stopped at port",self.port))

    # markAsChanged
    def markAsChanged(self, idx):
        """Mark a slot that has changed its value.
        """
        self.changed[idx] = 1

    # send
    def send(self, msg):
        """Send one or more messages to the server.

        \param msg (\c str) Message(s)
        """
        self.sock.sendto(msg, (self.host, self.port))

    # sendBuld
    def sendBulk(self):
        """Send the values of all slots that have changed since the last call.
        """
        msgs = ""
        for idx in list(self.changed.keys()):
            slot,id = self.slotobjs[idx]
            msgs += self.encodeValueMessage(id, idx, slot)

        self.send(msgs)
        self.changed = {}

    # createClientSlotAttribs
    def createClientSlotAttribs(self):
        slots = self.slots
        if slots==None:
            return

        self.forwarders = []
        slotnames = ["<none>", "IntSlot", "DoubleSlot", "BoolSlot", "Vec3Slot",
                     "Vec4Slot", "Mat3Slot", "Mat4Slot", "QuatSlot"]

        for varname, slotname, slotparams in self.iterSlotDefs(slots):
            # Create the slot object
            exec("slot = %s%s"%(slotname, slotparams))
            # Add the slot as attribute
            setattr(self, "%s_slot"%varname, slot)
            exec("self.addSlot('%s', self.%s_slot)"%(varname, varname))

            id = slotnames.index(slotname)
            idx = len(self.forwarders)
            
            self.slotobjs.append((slot,id))

            # Create the forwarder object
            if self.gather_messages:
                f = _GatherForwarder(self, idx, idx==len(slots)-1)
            else:
                f = _ImmediateForwarder(self, slot, id, idx)
            slot.addDependent(f)
            self.forwarders.append(f)

    # createServerSlotAttribs
    def createServerSlotAttribs(self):
        slots = self.slots
        if slots==None:
            return

        slotnames = ["<none>", "IntSlot", "DoubleSlot", "BoolSlot", "Vec3Slot",
                     "Vec4Slot", "Mat3Slot", "Mat4Slot", "QuatSlot"]

        for varname, slotname, slotparams in self.iterSlotDefs(slots):
            # Create the slot object
            exec("slot = %s%s"%(slotname, slotparams))
            # Add the slot as attribute
            setattr(self, "%s_slot"%varname, slot)
            exec("self.addSlot('%s', self.%s_slot)"%(varname, varname))

            id = slotnames.index(slotname)
            idx = len(self.forwarders)
            
            self.slotobjs.append((slot,id))


    # deleteSlotAttribs
    def deleteSlotAttribs(self):
        if self.slots==None:
            return
        
        for name,slot in self.slots:
            slotname = "%s_slot"%name
            exec("del self.%s"%slotname)

        self.slots = None

    # encodeValueMessage
    def encodeValueMessage(self, id, idx, slot):
        """Encode a binary value message.

        \param id (\c int) Message ID
        \param idx (\c int) Slot index
        \param slot (\c Slot) Corresponding slot
        """

        # int
        if id==1:
            return struct.pack(">HHi", id, idx, slot.getValue())
        # double
        elif id==2:
            return struct.pack(">HHd", id, idx, slot.getValue())
        # bool
        elif id==3:
            return struct.pack(">HHh", id, idx, slot.getValue())
        # vec3
        elif id==4:
            x,y,z = slot.getValue()
            return struct.pack(">HHddd", id, idx, x, y, z)
        # vec4
        elif id==5:
            x,y,z,w = slot.getValue()
            return struct.pack(">HHdddd", id, idx, x, y, z,w)
        # mat3
        elif id==6:
            M = slot.getValue()
            return struct.pack(">HHddddddddd", id, idx, *M.toList(rowmajor=True))
        # mat4
        elif id==7:
            M = slot.getValue()
            return struct.pack(">HHdddddddddddddddd", id, idx, *M.toList(rowmajor=True))
        # quat
        elif id==8:
            q = slot.getValue()
            return struct.pack(">HHdddd", id, idx, q.w, q.x, q.y, q.z)

    # decodeValueMessage
    def decodeValueMessage(self, id, msg):
        """Decode a binary value message.

        Raises a ValueError exception if the id is invalid. If the
        size of msg is too short, a struct.error exception is thrown.

        \param id (\c int) Message ID
        \param msg (\c str) Data part of the message (may contain additional messages)
        \return Tuple (Message size, slot index, value)
        """

        # int
        if id==1:
            id, v = struct.unpack(">Hi", msg[:6])
            return 6,id,v
        # double
        elif id==2:
            id, v = struct.unpack(">Hd", msg[0:10])
            return 10,id,v
        # bool
        elif id==3:
            id, v = struct.unpack(">Hh", msg[:4])
            return 4,id,v
        # vec3
        elif id==4:
            idx, x,y,z = struct.unpack(">Hddd", msg[:26])
            return (26, idx, vec3(x,y,z))
        # vec4
        elif id==5:
            idx, x,y,z,w = struct.unpack(">Hdddd", msg[:34])
            return (34, idx, vec4(x,y,z,w))
        # mat3
        elif id==6:
            idx, m11,m12,m13,m21,m22,m23,m31,m32,m33 = struct.unpack(">Hddddddddd", msg[:74])
            return (74, idx, mat3(m11,m12,m13,m21,m22,m23,m31,m32,m33))
        # mat4
        elif id==7:
            idx, m11,m12,m13,m14,m21,m22,m23,m24,m31,m32,m33,m34,m41,m42,m43,m44 = struct.unpack(">Hdddddddddddddddd", msg[:130])
            return (130, idx, mat4(m11,m12,m13,m14,m21,m22,m23,m24,m31,m32,m33,m34,m41,m42,m43,m44))
        # quat
        elif id==8:
            idx, w,x,y,z = struct.unpack(">Hdddd", msg[:34])
            return (34, idx, quat(w,x,y,z))
        else:
            raise ValueError("Unknown message id (%d)"%id)
            

    # iterSlotDefs
    def iterSlotDefs(self, slots):
        """Generates the variable name, slot name and slot parameters.

        This generator method takes a slot definition list and generates
        3-tuples (variable name, slot name, slot parameters).

        For example, the slot definition list [("xpos", "IntSlot(2)"),
        ("ypos", "IntSlot(12)")] will generate two tuples ("xpos", "IntSlot",
        "(2)") and ("ypos", "IntSlot", "(12)").
        """
        
        for varname, slotdef in slots:
            n = slotdef.find("(")
            if n!=-1:
                slotname = slotdef[:n]
                slotparams = slotdef[n:]
            else:
                slotname = slotdef
                slotparams = "()"
            yield varname, slotname, slotparams
        
