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
# $Id: eventmanager.py,v 1.4 2005/03/08 16:13:42 mbaas Exp $

## \file eventmanager.py
## Contains the EventManager class.

import sys, types, bisect
from copy import copy

# Receiver
class _Receiver:
    """Stores a receiver (=callable) and a priority.

    This class defines a comparison operator (that's why this class is
    used instead of a tuple (priority, receiver)).
    """
    def __init__(self, receiver, priority=1):
        self.receiver = receiver
        self.priority = priority


    def __str__(self):
        # If the receiver is an instance method, then obtain the corresponding
        # class name
        cls = getattr(self.receiver, "im_class", None)
        if cls!=None:
            s = getattr(cls, "__name__", "?")+"."
        else:
            s = ""
        s += getattr(self.receiver, "__name__", "<unnamed>")
        return "(%d, %s)"%(self.priority, s)

    __repr__ = __str__

    def __lt__(self, other):
        return self.priority<other.priority
    

# EventManager
class EventManager:
    """Manages any kind of events.

    System wide event receivers must indicate if the event should
    be consumed (return value = True) or passed to the scene wide
    receivers.
    """

    def __init__(self):
        """Constructor."""

        # System wide connections.
        # Key: Event name - Value: Sorted list of _Receivers
        self.system_connections = {}
        # Scene wide connections.
        # Key: Event name - Value: Sorted list of _Receivers
        self.scene_connections = {}

    def __str__(self):
        s = 70*"-"+"\n"
        s += "System events\n"
        s += 70*"-"+"\n"
        for event in self.system_connections:
            s+='Event: "%s"\n'%event
            for rec in self.system_connections[event]:
                s+="  %s\n"%rec
        s+="\n"
                
        s += 70*"-"+"\n"
        s += "Scene events\n"
        s += 70*"-"+"\n"
        for event in self.scene_connections:
            s+='Event: "%s"\n'%event
            for rec in self.scene_connections[event]:
                s+="  %s\n"%rec
        return s
        

    # event
    def event(self, name, *params, **keyargs):
        """Signal an event.

        When an event handler returns True the notification is interrupted
        (i.e. the event is consumed and not available anymore for other
        handlers).

        \param name (\c str) Name of the event.
        \return True, if any event handler returned True.

        \todo Aufruf der Empfaenger-Methode kann Exception verursachen.
        """

        # Process system wide connections
        receivers = self.system_connections.get(name, [])
        for rec in copy(receivers):
            if rec.receiver(*params, **keyargs):
                return True

        # Process scene wide connections
        receivers = self.scene_connections.get(name, [])
        for rec in copy(receivers):
            if rec.receiver(*params, **keyargs):
                return True

        return False

    # connect
    def connect(self, name, receiver, priority=10, system=False):
        """Connect a function or method to an event.

        The priority determines the calling order when the corresponding
        event is emitted. Receivers with a smaller priority are invoked first.

        \param name (\c str) Name of the event.
        \param receiver The receiving function or method, or an instance
                        of a class that must implement an on<Event>() method.
        \param priority (\c int) Priority of the receiver
        \param system (\c bool) Specifies if the connection is system wide or
                      not (default: \c False).
        """

        receiver = self._determine_receiver(name, receiver)

        # Disconnect the receiver if it was already connected
        try:
            self.disconnect(name,receiver,system)
        except KeyError:
            pass

        if system:
            connections = self.system_connections
        else:
            connections = self.scene_connections

        rec = _Receiver(receiver, priority)
        # Has the event already any connections? then add the new receiver
        if name in connections:
            bisect.insort(connections[name], rec)
#            connections[name].append(receiver)
        # otherwise create a new list
        else:
            connections[name] = [rec]

        return (name,receiver)

    # disconnect
    def disconnect(self, name, receiver=None, system=False):
        """Disconnect a function or method from an event.

        \param name (\c str) Name of the event.
        \param receiver The receiving function or method, or an instance
                        of a class that must implement an on<Event>() method.
        \param system (\c bool) Specifies if the connection is system wide or
                      not (default: \c False).
        """

        if isinstance(name,tuple):
            name,receiver = name

        if system:
            connections = self.system_connections
        else:
            connections = self.scene_connections

        # Are all connections to be removed?
        if receiver==None:
            if name in connections:
                del connections[name]
                return

        receiver = self._determine_receiver(name, receiver)

        # Has the event any connections at all?
        if name not in connections:
            raise KeyError('Receiver is not connected to event "%s"'%name)

        # Try to remove the connection
#        try:
#            connections[name].remove(receiver)
#        except ValueError:
#            raise KeyError('Receiver is not connected to event "%s"'%name)
        
        for i,rec in enumerate(connections[name]):
            if rec.receiver==receiver:
                break
        else:
            raise KeyError('Receiver is not connected to event "%s"'%name)

        del connections[name][i]

    # disconnectAll
    def disconnectAll(self, system=False):
        """Remove all connections.

        \param system (\c bool) Specifies if the system wide connections or
                      the scene wide connections shall be removed
                      (default: \c False).
        """
        if system:
            self.system_connections = {}
        else:
            self.scene_connections = {}


    ## private:
    def _determine_receiver(self, name, receiver):
        """Returns the receiver object.

        The receiver is the callable that gets called when an event
        occurs.

        \param name (\c str) Name of the event.
        \param receiver The receiving function or method, or an instance
                        of a class that must implement an on<Event>() method.
        \return Callable object.
        """

        methodname = "on"+name
        if hasattr(receiver, methodname):
            return getattr(receiver, methodname)
        
        if hasattr(receiver, "__call__"):
            return receiver
            
        raise ValueError("Receiver argument must be a callable or an object with an %s() method"%methodname)
        
######################################################################

_global_event_manager = EventManager()

# eventManager
def eventManager():
    """Returns the global event manager.

    \return Event manager (EventManager)
    """
    global _global_event_manager
    return _global_event_manager
