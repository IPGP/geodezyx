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

## \file idmanager.py
## \brief Contains the IDManager class.

# Exceptions...
class IDNotAllocated(Exception):
    """Exception class."""
    pass

# IDManager
class IDManager:
    """This class manages IDs.

    This class is used to generate IDs for GUI elements such as menu items.
    You can call allocID() to get an unused ID and later on freeID() to
    make the ID available again.
    """
    def __init__(self, startid=1):
        """Constructor.

        \param startid (\c int) The starting value for the IDs (default: 1)
        """
        # Each interval is given by its start and its end (=last item+1)
        self._intervals = [[startid,None]]
        self._startid = startid

    def __call__(self):
        return self.allocID()

    # allocID
    def allocID(self):
        """Request an ID.

        \return An unused ID (\c int).
        """
        id = self._intervals[0][0]
        self._intervals[0][0]+=1
        # Is the interval eaten up? then remove it
        if self._intervals[0][0]==self._intervals[0][1]:
            del self._intervals[0]
        return id

    # freeID
    def freeID(self, id):
        """Releases an ID.

        \param id (\c int) A previously allocated ID
        """
        if id<self._startid:
            raise IDNotAllocated("ID %d hasn't been allocated."%id)
        
        idx = 0
        for istart,iend in self._intervals:
            if id<istart:
                break
            # Check if id is in the current interval (=attempt to free an ID
            # which wasn't allocated)
            if iend==None:
                iend = id+1
            if id>=istart and id<iend:
                raise IDNotAllocated("ID %d hasn't been allocated."%id)
            idx+=1

        # Is the freed id id directly before the interval
        if id==self._intervals[idx][0]-1:
            self._intervals[idx][0]-=1
            # Check if the interval can be merged with the previous interval
            if idx!=0 and self._intervals[idx-1][1]==self._intervals[idx][0]:
                self._intervals[idx-1][1]=self._intervals[idx][1]
                del self._intervals[idx]

        # Is the freed id directly behind the previous interval? then enlarge
        elif idx!=0 and id==self._intervals[idx-1][1]:
            self._intervals[idx-1][1]+=1
            # Check if the enlarged interval can be merged with its successor
            if self._intervals[idx-1][1]==self._intervals[idx][0]:
                self._intervals[idx-1][1]=self._intervals[idx][1]
                del self._intervals[idx]
        else:
            # Create a new interval with one element
            self._intervals.insert(idx,[id,id+1])

