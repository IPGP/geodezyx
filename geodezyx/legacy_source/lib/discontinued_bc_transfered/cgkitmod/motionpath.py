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
# $Id: motionpath.py,v 1.2 2006/01/27 07:51:02 mbaas Exp $

from cgkit.component import Component
from cgkit.slots import *
from math import *

# MotionPath
class MotionPath(Component):
    """Motion path.
    """
    
    def __init__(self,
                 name = "MotionPath",
                 curve = None,
                 begintime = 0.0,
                 endtime = 1.0,
                 loop = False,
                 follow = False,
                 bank = False,
                 bankamplitude = 0.1,
                 auto_insert = True):
        """Constructor.
        """
        Component.__init__(self, name=name, auto_insert=auto_insert)
        
        self.curve = curve

        self.begintime = begintime
        self.endtime = endtime
        self.loop = loop
        self.follow = follow
        self.bank = bank
        self.bankamplitude = bankamplitude

        self.transform_slot = ProceduralMat4Slot(self.computeTransform)
        self.time_slot = DoubleSlot()
        self.time_slot.addDependent(self.transform_slot)
        self.addSlot("transform", self.transform_slot)
        self.addSlot("time", self.time_slot)

    # computeTransform
    def computeTransform(self):
        """Procedural for the transform slot.
        """

        # Determine the native curve parameter t...
        time = self.time_slot.getValue()
        len = self.curve.length()
        s = (time-self.begintime)/(self.endtime-self.begintime)
        if self.loop:
            s %= 1.0
        s *= len
        if s<0:
            s = 0.0
        if s>len:
            s = len
        t = self.arcLenToCurveParam(s)

        # Evaluate the curve
        p, dp, d2p = self.curve.evalFrame(t)

        up = vec3(0,0,1)
        if self.follow:
            try:
                T = mat4().lookAt(p, p+dp, up)
            except:
                T = mat4().translation(p)
        else:
            T = mat4().translation(p)

        if self.bank:
            dp = dp.normalize()

            dt = 0.001
            dp_prev = self.curve.deriv(t-dt).normalize()
            d2p = (dp-dp_prev)/dt
            # Make the 2nd derivative orthogonal to the tangent...
            len = d2p.length()
            u = d2p.cross(dp)
            d2p = dp.cross(u)
            d2p = len*d2p.normalize()            
            if (dp.cross(up)*d2p<0):
                len = -len
            T = T*mat4().rotation(self.bankamplitude*len, vec3(0,0,1))
        return T

    # arcLenToCurveParam
    def arcLenToCurveParam(self, s, eps=0.0001, maxiter=100):
        """Determine the native curve parameter for a given arc length.
        """
        tmin, tmax = self.curve.paraminterval
        totallen = self.curve.arcLen(tmax)
        # Initial "guess"...
        t = tmin+(s/totallen)*(tmax-tmin)

        while maxiter>0:
            F = self.curve.arcLen(t)-s
            if abs(F)<=eps:
                return t

            dF = self.curve.deriv(t).length()
            if abs(dF)<1E-12:
                # TODO: do something
                print("MotionPath: ************ dF==0 !!!")
                dF = 1 
            t -= F/dF
            maxiter -= 1

        print(("MotionPath: Maximum number of iterations reached (s=%f)"%s))
        return t
        
        
