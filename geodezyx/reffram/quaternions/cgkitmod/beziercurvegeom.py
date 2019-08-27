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
# $Id: beziercurvegeom.py,v 1.2 2005/09/02 12:58:02 mbaas Exp $

from ._OpenGL.GL import *
import bisect
from .geomobject import *
from .slots import *
from .cgtypes import *
from .boundingbox import BoundingBox
from . import protocols
from .ribexport import IGeometry
from .ri import *

# BezierPoint
class BezierPoint:
    """A single Bezier point with tangents.

    This object is used during construction of a BezierCurveGeom object.
    """
    def __init__(self, pos, intangent=vec3(0), outtangent=vec3(0)):
        self.pos = vec3(pos)
        self.intangent = vec3(intangent)
        self.outtangent = vec3(outtangent)

_uniform_size_constraint = UserSizeConstraint(1)
_null_size_constraint = UserSizeConstraint(0)

# BezierCurveGeom
class BezierCurveGeom(GeomObject):
    """Cubic Bezier curve.
    """

    protocols.advise(instancesProvide=[IGeometry])
    
    def __init__(self,
                 pnts = None,
                 closed = False,
                 epsilon = 0.01,
                 subdiv = 4,
                 show_tangents = False):
        """Constructor.

        pnts is a list of BezierPoint objects.
        """
        GeomObject.__init__(self)

        # Epsilon value for segLength()
        self.epsilon = epsilon
        # Number of subdivisions for drawGL()
        self.subdiv = subdiv
        # Show tangents or not?
        self.show_tangents = show_tangents

        # Create the slots...
        
        # The end points of the Bezier segments, i.e. the curve passes
        # through these points.
        self.pnts_slot = Vec3ArraySlot()
        self._var_sizeconstraint = LinearSizeConstraint(self.pnts_slot,1,0)
        if closed:
            b = 0
        else:
            b = -2
        self._vtx_sizeconstraint = LinearSizeConstraint(self.pnts_slot,3,b)
        # The incoming tangent for point i
        self.intangents_slot = Vec3ArraySlot(constraint=self._var_sizeconstraint)
        # The outgoing tangent for point i
        self.outtangents_slot = Vec3ArraySlot(constraint=self._var_sizeconstraint)
        
        # Create aliases for the slot "value" attribute
        self.pnts = self.pnts_slot
        self.intangents = self.intangents_slot
        self.outtangents = self.outtangents_slot

        # Add the slots to the component
        self.addSlot("pnts", self.pnts_slot)
        self.addSlot("intangents", self.intangents_slot)
        self.addSlot("outtangents", self.outtangents_slot)

        self._closed = closed

        if pnts!=None:
            # Initialize the Bezier points...
            self.pnts.resize(len(pnts))
            for i,bp in enumerate(pnts):
                self.pnts[i] = bp.pos
                self.intangents[i] = bp.intangent
                self.outtangents[i] = bp.outtangent

        # A list of lengths (the lengths are accumulated, so it's
        # the length of the entire curve until segment i)
        self._seg_lengths = []
        # Flag that indicates if all internal precomputed attributes have
        # to be recomputed or if they are still valid.
        self._internal_data_invalid = True

        # Have the slots call the onCurveChanged() method whenever
        # a point or a tangent is modified.
        self._nf = NotificationForwarder(self.onCurveChanged, self.onCurveResized)
        self.pnts_slot.addDependent(self._nf)
        self.intangents_slot.addDependent(self._nf)
        self.outtangents_slot.addDependent(self._nf)


    def _updateInternalData(self):
        """Precalculate the lengths of the Bezier segments.
        """
        self._seg_lengths = []
        sl = 0.0
        if self._closed:
            k = 0
        else:
            k = 1
        for i in range(self.pnts.size()-k):
            sl += self.segLength(self.segCtrlPoints(i), self.epsilon)
            self._seg_lengths.append(sl)

        self._internal_data_invalid = False
        

#    def uniformCount(self):
#        return 1

#    def varyingCount(self):
#        # One value per segment boundary, i.e.
#        # self.numsegs   for periodic curves and
#        # self.numsegs+1 for non-periodic curves
#        # which can be reduced to self.pnts.size() for both cases
#        # (as the curve points just mark segment boundaries)
#        return self.pnts.size()

#    def vertexCount(self):
#        # One value per control points,
#        n = 3*self.numsegs
#        if not self._closed:
#            n += 1
#        return n

    # slotSizeConstraint
    def slotSizeConstraint(self, storage):
        global _uniform_size_constraint
        
        if storage==UNIFORM:
            return _uniform_size_constraint
        elif storage==VARYING:
            return self._var_sizeconstraint
        elif storage==VERTEX:
            return self._vtx_sizeconstraint
        else:
            return _null_size_constraint

    # boundingBox
    def boundingBox(self):
        """Return the bounding box of the control polygon.
        """
        bb = BoundingBox()
        for i,p in enumerate(self.pnts):
            bb.addPoint(p)
            bb.addPoint(p+self.intangents[i])
            bb.addPoint(p+self.outtangents[i])
        return bb

    # drawGL
    def drawGL(self):
        if self.pnts.size()==0:
            return

        glPushAttrib(GL_LIGHTING_BIT)
        glDisable(GL_LIGHTING)
        glColor3f(1,1,1)

        # Draw the curve...
        p0 = self.pnts[0]
        in0 = self.intangents[0]
        out0 = self.outtangents[0]
        for i in range(1, self.pnts.size()):
            p1 = self.pnts[i]
            in1 = self.intangents[i]
            out1 = self.outtangents[i]
            self.segDrawGL([p0, p0+out0, p1+in1, p1], self.subdiv)
            p0 = p1
            in0 = in1
            out0 = out1

        # Draw an additional segment if the curve is closed...
        if self._closed:
            p0 = self.pnts[-1]
            out0 = self.outtangents[-1]
            p1 = self.pnts[0]
            in1 = self.intangents[0]
            self.segDrawGL([p0, p0+out0, p1+in1, p1], self.subdiv)
            

        # Draw the in and out tangents...
        if self.show_tangents:
            for i in range(self.pnts.size()):
                p = self.pnts[i]
                pin = self.intangents[i]
                pout = self.outtangents[i]
                glBegin(GL_LINE_STRIP)
                glColor3f(0,1,0)
                v = p+pin
                glVertex3d(v.x, v.y, v.z)
                glVertex3d(p.x, p.y, p.z)
                glColor3f(0,0,1)
                glVertex3d(p.x, p.y, p.z)
                v = p+pout
                glVertex3d(v.x, v.y, v.z)
                glEnd()            

        glPopAttrib()

    # render
    def render(self, matid):
        if matid==0:
            # Set Bezier basis
            RiBasis(RiBezierBasis, RI_BEZIERSTEP, RiBezierBasis, RI_BEZIERSTEP)

            # Create the list of curve vertices...
            pnts = []
            for i in range(self.numsegs):
                ps = self.segCtrlPoints(i)
                pnts += ps[0:3]

            if self._closed:
                wrap = RI_PERIODIC
            else:
                wrap = RI_NONPERIODIC
                pnts.append(self.pnts[-1])

            # Create parameter list...
            params = {"P":pnts}
            clss = ["constant", "uniform", "varying", "vertex",
                    "facevarying", "facevertex", "user"]
            typs = ["integer", "float", "string", "color", "point",
                    "vector", "normal", "matrix", "hpoint"]
            for name, storage, type, multiplicity in self.iterVariables():
                cls = clss[abs(storage)]
                typ = typs[abs(type)]
                if cls=="user":
                    continue
                s = self.slot(name)
                params[name] = list(s)
                if multiplicity==1:
                    decl = "%s %s"%(cls, typ)
                else:
                    decl = "%s %s[%d]"%(cls, typ, multiplicity)
                RiDeclare(name, decl)

            # Output a curve primitive
            RiCurves(RI_CUBIC, [len(pnts)], wrap, params)

    # eval
    def eval(self, t):
        """Evaluate the curve at parameter t and return the curve point.
        """
        # Segment number
        seg = int(t)
        # Segment parameter
        t = t%1.0
#        if seg>=self.pnts.size()-1:
#            seg = self.pnts.size()-2
#            t = 1.0
        return self.segEval(t, self.segCtrlPoints(seg))        

    # evalFrame
    def evalFrame(self, t):
        """Evaluate the curve at parameter t and return a coordinate frame.
        """
        # Segment number
        seg = int(t)
        # Segment parameter
        t = t%1.0
#        if seg>=self.pnts.size()-1:
#            seg = self.pnts.size()-2
#            t = 1.0
        return self.segEvalFrame(t, self.segCtrlPoints(seg))

    # deriv
    def deriv(self, t):
        """Return the derivative at parameter t.
        """
        # Segment number
        seg = int(t)
        # Segment parameter
        t = t%1.0
#        if seg>=self.pnts.size()-1:
#            seg = self.pnts.size()-2
#            t = 1.0
        return self.segDeriv(t, self.segCtrlPoints(seg))

    # arcLen
    def arcLen(self, t):
        """Return the arc length up to t.
        """
        if self._internal_data_invalid:
            self._updateInternalData()
            
        seg = int(t)
        t = t%1.0
        if seg>=len(self._seg_lengths):
            return self._seg_lengths[-1]

        res = 0.0
        if seg>0:
            res = self._seg_lengths[seg-1]
        b,c = self.segSplit(self.segCtrlPoints(seg), t)
        return res + self.segLength(b, self.epsilon)

    # length
    def length(self):
        """Return the length of the entire curve.
        """
        if self._internal_data_invalid:
            self._updateInternalData()
        if len(self._seg_lengths)==0:
            return 0.0
        else:
            return self._seg_lengths[-1]


    def eval0(self, t):
        if self._internal_data_invalid:
            self._updateInternalData()

        # Determine the Bezier segment index we're currently in
        seg = bisect.bisect_right(self._seg_lengths, t)
        if seg==self.pnts.size()-1:
            return self.pnts[-1], vec3(), vec3()

        if seg>0:
            t -= self._seg_lengths[seg-1]

        p,t = self.segEval0(t, self.segCtrlPoints(seg))
        if p!=None:
            return p, vec3(), vec3()
        else:
            return self.pnts[-1], vec3(), vec3()

    # onCurveChanged
    def onCurveChanged(self, start, end):
        """This method is called whenever a point or tangent is modified.
        """
        self._internal_data_invalid = True

    # onCurveResized
    def onCurveResized(self, size):
        """This method is called whenever the number of points is modified.
        """
        self._internal_data_invalid = True
        
    ## protected:

    def segSmooth(self, seg):
        b0,b1,b2,b3 = self.segCtrlPoints(seg-1)
        d1 = b2-b1
        d2 = b3-b2
        self.outtangents[seg] = d2
        self.intangents[seg+1] = self.pnts[seg]+3*d2-d1 - self.pnts[seg+1]

    def segCtrlPoints(self, i):
        """Return the Bezier control points of segment i.

        i may be an arbitrary segment number (the numbering is cyclic).
        """
        s = self.pnts.size()
        i = i%s
        j = (i+1)%s
        b0 = self.pnts[i]
        b3 = self.pnts[j]
        b1 = b0+self.outtangents[i]
        b2 = b3+self.intangents[j]
        return [b0, b1, b2, b3]

    # segEval0
    def segEval0(self, t, ctrlpnts, eps=0.001):
        """

        Return value: (pnt, t2)
        pnt is None if t>length. t2 is t-length.
        """
        l1 = self.segHullLength(ctrlpnts)
        b,c = self.segSplit(ctrlpnts)
        bl = self.segHullLength(b)
        cl = self.segHullLength(c)
        l2 = bl+cl
        # Is the result accurate enough? then return the current result,
        # otherwise subdivide further
        if l1-l2<=eps:
            p,t = self.segHullEval0(t, b)
            if p!=None:
                return p,t
            p,t = self.segHullEval0(t, c)
            if p!=None:
                return p,t
            return None, t
        else:
            p,t = self.segEval0(t, b, eps=eps)
            if p!=None:
                return p,t
            p,t = self.segEval0(t, c, eps=eps)
            if p!=None:
                return p,t
            return None, t

    def segHullEval0(self, t, ctrlpnts):
        b0,b1,b2,b3 = ctrlpnts
        # 1st segment
        d = (b1-b0).length()
        if t<d:
            a = t/d
            return (1-a)*b0 + a*b1, t
        t-=d
        # 2nd segment
        d = (b2-b1).length()
        if t<d:
            a = t/d
            return (1-a)*b1 + a*b2, t
        t-=d
        # 3rd segment
        d = (b3-b2).length()
        if t<d:
            a = t/d
            return (1-a)*b2 + a*b3, t
        t-=d
        return None, t
        
    # segEval
    def segEval(self, t, ctrlpnts):
        """Evaluate the curve at parameter t using the de Casteljau scheme.

        t ranges from 0 to 1.
        """
        _t = 1.0-t
        b0,b1,b2,b3 = ctrlpnts
        c0 = _t*b0 + t*b1
        c1 = _t*b1 + t*b2
        c2 = _t*b2 + t*b3
        d0 = _t*c0 + t*c1
        d1 = _t*c1 + t*c2
        return _t*d0+t*d1

    # segEvalFrame
    def segEvalFrame(self, t, ctrlpnts):
        """Evaluate the curve at parameter t using the de Casteljau scheme.

        t ranges from 0 to 1.
        """
        _t = 1.0-t
        b0,b1,b2,b3 = ctrlpnts
        c0 = _t*b0 + t*b1
        c1 = _t*b1 + t*b2
        c2 = _t*b2 + t*b3
        d0 = _t*c0 + t*c1
        d1 = _t*c1 + t*c2
        e0 = _t*d0+t*d1
        return e0, 3*(d1-d0), 6*(c2-2*c1+c0)

    # segDeriv
    def segDeriv(self, t, ctrlpnts):
        """Return the derivation at t.

        t ranges from 0 to 1.
        """
        _t = 1.0-t
        b0,b1,b2,b3 = ctrlpnts
        c0 = _t*b0 + t*b1
        c1 = _t*b1 + t*b2
        c2 = _t*b2 + t*b3
        d0 = _t*c0 + t*c1
        d1 = _t*c1 + t*c2
        return 3*(d1-d0)


    # segDrawGL
    def segDrawGL(self, ctrlpnts, subdiv=4):
        """Draw one Bezier segment.

        Actually the control polygon is drawn after a number of
        subdivisions.
        """
        if subdiv==0:
            b0,b1,b2,b3 = ctrlpnts
            glBegin(GL_LINE_STRIP)
            glVertex3d(b0.x, b0.y, b0.z)
            glVertex3d(b1.x, b1.y, b1.z)
            glVertex3d(b2.x, b2.y, b2.z)
            glVertex3d(b3.x, b3.y, b3.z)
            glEnd()
        else:
            b,c = self.segSplit(ctrlpnts)
            self.segDrawGL(b, subdiv=subdiv-1)
            self.segDrawGL(c, subdiv=subdiv-1)
            
    # segLength
    def segLength(self, ctrlpnts, eps=0.001):
        """Return the length of one Bezier segment.

        ctrlpnts is a list of four vec3s.
        """
        l1 = self.segHullLength(ctrlpnts)
        b,c = self.segSplit(ctrlpnts)
        bl = self.segHullLength(b)
        cl = self.segHullLength(c)
        l2 = bl+cl
        # Is the result accurate enough? then return the current result,
        # otherwise subdivide further
        if l1-l2<=eps:
            return l2
        else:
            return self.segLength(b,eps)+self.segLength(c,eps)

    # segHullLength
    def segHullLength(self, ctrlpnts):
        """Return the length of the control "polygon".
        """
        b0,b1,b2,b3 = ctrlpnts
        return (b1-b0).length() + (b2-b1).length() + (b3-b2).length()

    # segSplit
    def segSplit(self, ctrlpnts, t=0.5):
        """Split a segment into two segments.
        """
        b0,b1,b2,b3 = ctrlpnts
        _t = 1.0-t
        c0 = _t*b0 + t*b1
        c1 = _t*b1 + t*b2
        c2 = _t*b2 + t*b3
        d0 = _t*c0 + t*c1
        d1 = _t*c1 + t*c2
        e0 = _t*d0 + t*d1
        return [b0,c0,d0,e0], [e0,d1,c2,b3]

    # "paraminterval" property...
    
    def _getParamInterval(self):
        """Return the current parameter interval.

        This method is used for retrieving the \a paraminterval property.

        \return (t_min, t_max)
        """
        if self._closed:
            return 0, self.pnts.size()
        else:
            return 0, self.pnts.size()-1

    paraminterval = property(_getParamInterval, None, None, "Parameter range")

    # "closed" property...
    
    def _getClosed(self):
        return self._closed

    def _setClosed(self, c):
        if c==self._closed:
            return
        self._closed = c
        self._internal_data_invalid = True
        # Update the size constraint for vertex variables
        if self._closed:
            self._vtx_sizeconstraint.setCoeffs(3,0)
        else:
            self._vtx_sizeconstraint.setCoeffs(3,-2)

    closed = property(_getClosed, _setClosed, None, "Closed")

    # "numsegs" property...
    
    def _getNumSegs(self):
        if self._closed:
            return self.pnts.size()
        else:
            return self.pnts.size()-1

    numsegs = property(_getNumSegs, None, None, "Number of Bezier segments")


