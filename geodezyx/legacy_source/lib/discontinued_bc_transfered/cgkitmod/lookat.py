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
# $Id: lookat.py,v 1.2 2005/07/21 15:07:11 mbaas Exp $

## \file lookat.py
## Contains the LookAt component.

"""This module contains the LookAt component which has a position
slot (pos) and a target slot (target) which are both of type vec3.
The output (output) is a mat3 that contains a rotation so that when positioned
at pos the z axis points to target."""

from .component import createFunctionComponent
from .cgtypes import *
from .sl import radians

def _lookat(pos=vec3(0), target=vec3(0), up=vec3(0,0,1), roll=0.0):
    try:
        M = mat4().lookAt(pos, target, up)
    except:
        M = mat4(1)

    return M.getMat3()*mat3().rotation(radians(roll), vec3(0,0,1))

LookAt = createFunctionComponent(_lookat) # "mat3 (vec3 pos, vec3 target, vec3 up, double roll)")
