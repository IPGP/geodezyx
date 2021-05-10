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
# $Id: constants.py,v 1.1 2005/01/09 20:13:06 mbaas Exp $

# Information categories
WTI_INTERFACE = 1
WTI_STATUS = 2
WTI_DEFCONTEXT = 3
WTI_DEFSYSCTX = 4
WTI_DDCTXS = 400
WTI_DSCTXS = 500
WTI_DEVICES = 100
WTI_CURSORS = 200
WTI_EXTENSIONS = 300

# Hardware capabilities
HWC_INTEGRATED = 0x0001
HWC_TOUCH = 0x0002
HWC_HARDPROX = 0x0004
HWC_PHYSID_CURSORS = 0x0008

# Unit specifiers
TU_NONE = 0
TU_INCHES = 1
TU_CENTIMETERS = 2
TU_CIRCLE = 3

# Cursor capabilities
CRC_MULTIMODE = 0x0001
CRC_AGGREGATE = 0x0002
CRC_INVERT = 0x0004

# System button assignment values 
SBN_NONE = 0x00
SBN_LCLICK = 0x01
SBN_LDBLCLICK = 0x02
SBN_LDRAG = 0x03
SBN_RCLICK = 0x04
SBN_RDBLCLICK = 0x05
SBN_RDRAG = 0x06
SBN_MCLICK = 0x07
SBN_MDBLCLICK = 0x08
SBN_MDRAG = 0x09
# for Pen Windows
SBN_PTCLICK = 0x10
SBN_PTDBLCLICK = 0x20
SBN_PTDRAG = 0x30
SBN_PNCLICK = 0x40
SBN_PNDBLCLICK = 0x50
SBN_PNDRAG = 0x60
SBN_P1CLICK = 0x70
SBN_P1DBLCLICK = 0x80
SBN_P1DRAG = 0x90
SBN_P2CLICK = 0xA0
SBN_P2DBLCLICK = 0xB0
SBN_P2DRAG = 0xC0
SBN_P3CLICK = 0xD0
SBN_P3DBLCLICK = 0xE0
SBN_P3DRAG = 0xF0

# Context option values
CXO_SYSTEM = 0x0001
CXO_PEN = 0x0002
CXO_MESSAGES = 0x0004
CXO_MARGIN = 0x8000
CXO_MGNINSIDE = 0x4000
CXO_CSRMESSAGES = 0x0008

# Context status values
CXS_DISABLED = 0x0001
CXS_OBSCURED = 0x0002
CXS_ONTOP = 0x0004

# Context lock values 
CXL_INSIZE = 0x0001
CXL_INASPECT = 0x0002
CXL_SENSITIVITY = 0x0004
CXL_MARGIN = 0x0008
CXL_SYSOUT = 0x0010

# WTPKT bits
PK_CONTEXT = 0x0001          # reporting context
PK_STATUS = 0x0002           # status bits 
PK_TIME = 0x0004             # time stamp 
PK_CHANGED = 0x0008          # change bit vector 
PK_SERIAL_NUMBER = 0x0010    # packet serial number
PK_CURSOR = 0x0020           # reporting cursor 
PK_BUTTONS = 0x0040          # button information
PK_X = 0x0080                # x axis 
PK_Y = 0x0100                # y axis
PK_Z = 0x0200                # z axis
PK_NORMAL_PRESSURE = 0x0400  # normal or tip pressure 
PK_TANGENT_PRESSURE = 0x0800 # tangential or barrel pressure
PK_ORIENTATION = 0x1000      # orientation info: tilts */
PK_ROTATION = 0x2000         # rotation info; 1.1 

# Packet status values
TPS_PROXIMITY =	0x0001
TPS_QUEUE_ERR =	0x0002
TPS_MARGIN = 0x0004
TPS_GRAB = 0x0008
TPS_INVERT = 0x0010

# Relative buttons
TBN_NONE = 0
TBN_UP = 1
TBN_DOWN = 2
