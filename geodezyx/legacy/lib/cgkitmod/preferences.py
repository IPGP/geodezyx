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
# $Id: preferences.py,v 1.1.1.1 2004/12/12 14:31:13 mbaas Exp $

## \file preferences.py
## \brief Contains the Preferences class.

import sys, os, os.path, pickle

# Preferences
class Preferences(object):
    """This class stores arbitrary values persistently.

    This class can be used just like a dictionary.
    """

    def __init__(self, filename):
        """Constructor."""

        # If filename is no absolute path make it relative to the default
        # config file path
        if not os.path.isabs(filename):
            filename = os.path.join(configPath(), filename)
            
        self._prefs = {}
        self._filename = filename

    def __len__(self):
        return len(self._prefs)

    def __iter__(self):
        return iter(list(self._prefs.items()))

    def __getitem__(self, key):
        if key in self._prefs:
            return self._prefs[key]
        else:
            return None

    def __setitem__(self, key, value):
        self._prefs[key]=value

    def __delitem__(self, key):
        if key in self._prefs:
            del self._prefs[key]

    def __contains__(self, item):
        return item in self._prefs

    def get(self, key, default=None):
        return self._prefs.get(key, default)

    # load
    def load(self):
        """Load the preferences.

        """

        f = open(self._filename)
        id = pickle.load(f)
        if id!=1:
            raise Exception("Unknown config file format")
        self._prefs = pickle.load(f)
        f.close()


    # save
    def save(self):
        """Save the preferences.

        """
        
        self._preparePath(os.path.dirname(self._filename))
        
        f = open(self._filename, "w")
        pickle.dump(1, f)
        pickle.dump(self._prefs, f)
        f.close()


    ######################################################################
    ## protected:

    def _preparePath(self, path):
        """Prepare a path so that every directory on the path exists.

        Checks if the path exists and creates it if it does not exist.
        """
        if not os.path.exists(path):
            parent = os.path.dirname(path)
            if parent!="":
                self._preparePath(parent)
            os.mkdir(path)        

    # "filename" property

    def _getFilename(self):
        """Return the filename.

        This method is used for retrieving the \a filename property.

        \return Filename (\c str).
        """
        return self._filename

    filename = property(_getFilename, None, None, "File name")
    

######################################################################

def configPath(appname="gaia"):
    """Return the full path where config files are located.

    \todo Change "gaia" as default application name
    """
    
    # Windows? (XP)
    if sys.platform=="win32":
        if "APPDATA" in os.environ:
            return os.path.abspath(os.path.join(os.environ["APPDATA"], appname.capitalize()))
        else:
            return None
    # Other
    else:
        if "HOME" in os.environ:
            return os.path.abspath(os.environ["HOME"])
        else:
            return None


######################################################################

if __name__=="__main__":

    print((configPath()))
            
            

    
