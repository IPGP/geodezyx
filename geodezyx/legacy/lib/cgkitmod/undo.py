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
# $Id: undo.py,v 1.1.1.1 2004/12/12 14:31:30 mbaas Exp $

## \file undo.py
## \brief Provides the undo/redo framework.

class UndoError(Exception):
    pass

# _ReverseIterator
class _ReverseIterator:
    """Iterator that iterates over a list in reverse order.
    """
    def __init__(self, list):
        self._idx = len(list)-1
        self._list = list

    def __iter__(self):
        return self

    def __next__(self):
        if self._idx==-1:
            raise StopIteration()
        res = self._list[self._idx]
        self._idx -= 1
        return res


# UndoManager
class UndoManager:
    """This class manages undo objects.

    This class can be used to undo/redo operations. Each operation is
    represented by an undo object (UndoObject) that knows how to undo
    and redo this operation. Whenever an operation is invoked for the
    first time it has to create an undo object and push() it onto the
    undo stack in the manager. The undo() method then pops the undo
    objects from the stack performing their undo action and pushes
    them on the redo stack. The redo() method just does the opposite.
    Whenever push() is called the redo stack is emptied.

    You can iterate over all undo objects in the undo stack by using
    the manager as a sequence or by the iterUndo() method. The iterRedo()
    method iterates over all the undo objects in the redo stack.
    """
    
    def __init__(self, maxundoops=None):
        """Constructor.

        \param maxundoops (\c int) Maximum number of undo operations to maintain or None for an unlimited number.
        """
        # The last element is the top element
        self._undo_stack = []
        self._redo_stack = []
        # Maximum number of undo objects or None (=unlimited)
        self._max_undo_ops = maxundoops

    def __len__(self):
        return len(self._undo_stack)

    def __iter__(self):
        return _ReverseIterator(self._undo_stack)

    # undoCount
    def undoCount(self):
        """Return the number of operations on the undo stack.

        \return Number of operation son the undo stack (\c int).
        """
        return len(self._undo_stack)
    
    # iterUndo
    def iterUndo(self):
        """Iterate over all undo objects in the undo stack (from top to bottom)."""
        return _ReverseIterator(self._undo_stack)

    # redoCount
    def redoCount(self):
        """Return the number of operations on the redo stack.

        \return Number of operation son the redo stack (\c int).
        """
        return len(self._redo_stack)

    # iterRedo
    def iterRedo(self):
        """Iterate over all undo objects in the redo stack (from top to bottom)."""
        return _ReverseIterator(self._redo_stack)

    # clear
    def clear(self):
        """Clear the undo/redo stack.

        All undo objects are removed from both stacks.
        """
        del self._undo_stack[:]
        del self._redo_stack[:]

    def undoBegin(self, desc):
        """Start an undo block.

        All following undo operations are combined into one single
        undo operation.

        \param desc (\c str) Description text.
        """
        pass

    def undoEnd(self):
        pass

    # undo
    def undo(self):
        """Performs an undo operation.

        The last operation is undone and pushed on the redo stack.
        If the undo stack is empty an UndoError exception is thrown.

        If the undo operation throws an exception, then both stacks
        are discarded and the exception is propagated to the caller.
        """
        if len(self._undo_stack)==0:
            raise UndoError("There is no operation to undo.")
        u = self._undo_stack.pop()
        try:
            u.undo()
        except:
            self.clear()
            raise
        self._redo_stack.append(u)

    # redo
    def redo(self):
        """Performs a redo operation.

        The last undo operation is redone. If the redo stack is empty
        an UndoError exception is thrown.

        If the redo operation throws an exception, then both stacks
        are discarded and the exception is propagated to the caller.
        """
        if len(self._redo_stack)==0:
            raise UndoError("There is no operation to redo.")
        u = self._redo_stack.pop()
        try:
            u.redo()
        except:
            self.clear()
            raise
        self._undo_stack.append(u)

    # push
    def push(self, undoobj):
        """Push an undo object on the stack.

        \param undoobj (\c UndoObject) Undo object
        """
        self._undo_stack.append(undoobj)
        del self._redo_stack[:]
        # Make sure there are no more items on the stack than the
        # specified maximum number of operations...
        if self._max_undo_ops!=None:
            if len(self._undo_stack)>self._max_undo_ops:
                del self._undo_stack[0]


# UndoObject
class UndoObject:
    """Base undo object which represents an undoable operation.

    This class has the following attributes:

    - desc (\c str): A short description describing the operation.
         This description might be shown in the undo menu.
    """
    
    def __init__(self, desc):
        self.description = desc

    def __str__(self):
        return "<Undo object '%s'>"%self.desc

    # undo
    def undo(self):
        """Performs an undo operation."""
        raise UndoError("No undo operation implemented.")

    # redo
    def redo(self):
        """Performs a redo operation.

        This method may only be called if undo() was called previously.
        """
        raise UndoError("No redo operation implemented.")

######################################################################

_undo_manager = {None:UndoManager()}

def undoCount(stackid=None):
    global _undo_manager
    return _undo_manager[stackid].undoCount()

def redoCount(stackid=None):
    global _undo_manager
    return _undo_manager[stackid].redoCount()

def iterUndo(stackid=None):
    global _undo_manager
    return _undo_manager[stackid].iterUndo()

def iterRedo(stackid=None):
    global _undo_manager
    return _undo_manager[stackid].iterRedo()

def clear(stackid=None):
    global _undo_manager
    _undo_manager[stackid].clear()

def undo(stackid=None):
    global _undo_manager
    _undo_manager[stackid].undo()

def redo(stackid=None):
    global _undo_manager
    _undo_manager[stackid].redo()

def push(undoobj, stackid=None):
    global _undo_manager
    _undo_manager[stackid].push(undoobj)
