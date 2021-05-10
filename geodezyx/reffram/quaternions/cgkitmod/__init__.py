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
# $Id: __init__.py,v 1.19 2006/02/14 19:29:39 mbaas Exp $

"""Python Computer Graphics Kit
============================

The Python Computer Graphics Kit is a generic 3D package that can be
useful in any domain where you have to deal with 3D data of any kind,
be it for visualization, creating photorealistic images, Virtual
Reality or even games.

At its lowest level, the package provides the basic functionality that
is useful for writing your own tools that process 3D data. For
example, the cgtypes module defines the fundamental types for computer
graphics such as vectors and matrices, the ri module contains the
complete RenderMan API to create RIB files, and so on.

Using these relatively low level modules the second level provides the
functionality to store a 3D scene in memory. This is the major part of
the package and this actually turns the general purpose language
Python into a specialized scripting language like MEL or MaxScript,
for example.

Eventually, the package provides small tools that let you actually see
your 3D scene. The two standard tools are for interactive rendering
using OpenGL and for offline rendering using a RenderMan renderer.

See the manual and tutorials at: http://cgkit.sourceforge.net/
"""

import os.path
from . import cgkitinfo

# Is this cgkit light? Then modify the subpackage searchpath so that the
# special 'light' modules are found instead of the normal ones.
if cgkitinfo.cgkit_light:
    # Package's main folder
    dirname = __path__[0]
    # Put the 'light' subpackage in front so that its contents shadows
    # the normal extension modules...
    __path__.insert(0, os.path.join(dirname, "light"))
