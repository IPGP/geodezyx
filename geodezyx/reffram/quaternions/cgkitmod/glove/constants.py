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
# $Id: constants.py,v 1.1 2005/01/17 21:52:55 mbaas Exp $

# Glove handedness
FD_HAND_LEFT = 0
FD_HAND_RIGHT = 1

HANDEDNESS = { FD_HAND_LEFT : "Left",
               FD_HAND_RIGHT : "Right" }

# Glove type
FD_GLOVENONE = 0
FD_GLOVE5U = 1       # DG5 Ultra serial
FD_GLOVE5UW = 2      # DG5 Ultra serial, wireless
FD_GLOVE5U_USB = 3   # DG5 Ultra USB
FD_GLOVE7 = 4        # 7-sensor
FD_GLOVE7W = 5       # 7-sensor, wireless
FD_GLOVE16 = 6       # 16-sensor
FD_GLOVE16W = 7      # 16-sensor, wireless
FD_GLOVE14U = 8      # DG14 Ultra serial
FD_GLOVE14UW = 9     # DG14 Ultra serial, wireless
FD_GLOVE14U_USB	= 10 # DG14 Ultra USB

GLOVETYPE = { FD_GLOVENONE    : "none",
              FD_GLOVE5U      : "Data Glove 5 Ultra serial",
              FD_GLOVE5UW     : "Data Glove 5 Ultra serial, wireless",
              FD_GLOVE5U_USB  : "Data Glove 5 Ultra USB",
              FD_GLOVE7       : "Data Glove 5",
              FD_GLOVE7W      : "Data Glove 5, wireless",
              FD_GLOVE16      : "Data Glove 16",
              FD_GLOVE16W     : "Data Glove 16, wireless",
              FD_GLOVE14U     : "Data Glove 14 Ultra serial",
              FD_GLOVE14UW    : "Data Glove 14 Ultra serial, wireless",
              FD_GLOVE14U_USB : "Data Glove 14 Ultra USB" }

# the following are the glove type values of the 1.x SDK
#FD_GLOVENONE = 0 # no glove
#FD_GLOVE7 = 1    # 7-sensor
#FD_GLOVE7W = 2   # 7-sensor, wireless
#FD_GLOVE16 = 3   # 16-sensor
#FD_GLOVE16W = 4  # 16-sensor, wireless

# Sensors
FD_THUMBNEAR = 0
FD_THUMBFAR = 1
FD_THUMBINDEX = 2
FD_INDEXNEAR = 3
FD_INDEXFAR = 4
FD_INDEXMIDDLE = 5
FD_MIDDLENEAR = 6
FD_MIDDLEFAR = 7
FD_MIDDLERING = 8
FD_RINGNEAR = 9
FD_RINGFAR = 10
FD_RINGLITTLE = 11
FD_LITTLENEAR = 12
FD_LITTLEFAR = 13
FD_THUMBPALM = 14
FD_WRISTBEND = 15
FD_PITCH = 16
FD_ROLL = 17

# Product IDs for USB Gloves
DG14U_R = 0x00
DG14U_L	= 0x01
DG5U_R = 0x10
DG5U_L = 0x11
