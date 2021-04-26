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
# Portions created by the Initial Developer are Copyright (C) 2009
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
# $Id: riutil.py,v 1.1.1.1 2004/12/12 14:31:21 mbaas Exp $

import os.path
import ctypes
import functools
from . import cgtypes
from . import rmanlibutil
from . import _pointcloud
try:
    import numpy
    _numpy_available = True
except ImportError:
    _numpy_available = False

def _arrayPointer(a, n):
    """Check an array and return a pointer to its start.
    
    a is an array object and n the minimum number of values the array
    must contain. An exception is thrown when the array is too small
    or its type is not float (4 bytes) or when a is not a recognized
    array at all.
    """
    # Is the array a ctypes array?
    if isinstance(a, ctypes.Array):
        if a._type_!=ctypes.c_float:
            raise TypeError("Float array expected")
        if len(a)<n:
            raise TypeError("Array is not large enough")
        return ctypes.addressof(a)
    # Does the array support the array interface? (e.g. a numpy array)
    elif hasattr(a, "__array_interface__"):
        strides = a.__array_interface__.get("strides", None)
        typestr = a.__array_interface__.get("typestr", None)
        shape = a.__array_interface__.get("shape", None)
        data = a.__array_interface__.get("data", None)
        if strides is not None:
            raise TypeError("Unsupported array type (strides are not supported)")
        if not (isinstance(typestr, str) and len(typestr)==3 and typestr[1:3]=="f4"):
            raise TypeError("Unsupported array type (the array must contain 4-byte floats)")
        if type(shape) is not tuple:
            raise TypeError("Unsupported array type (unexpected shape value)")
        numFloats = functools.reduce(lambda x,y: x*y, shape)
        if numFloats<n:
            raise TypeError("Array is not large enough")
        if type(data) is not tuple or len(data)!=2:
            raise TypeError("Unsupported array type (unsupported data value)")
        if data[1]==True:
            raise TypeError("Array is read-only")
        if type(data[0]) not in [int, int]:
            raise TypeError("Unsupported array type (data pointer is not an int)")
        return data[0]
    # Unknown array
    else:
        raise TypeError("Unknown array type")


class PtcReader(object):
    """Point cloud reader class.
    
    An instance of this class is returned by the open() function.
    """
    
    def __init__(self, fileName, libName):
        """Constructor.
        
        fileName is the name of the point cloud file and libName the name
        of the shared library that implements the point cloud API.
        """
        object.__init__(self)
        
        self._handle = None
        
        ptclib = self._loadPtcLib(libName)
        
        # Store the functions pointers that are still required
        self._PtcReadDataPoint = ptclib.PtcReadDataPoint
        self._PtcClosePointCloudFile = ptclib.PtcClosePointCloudFile

        self._PtcGetPointCloudInfo = ptclib.PtcGetPointCloudInfo

        # Set 64 as default (which is the maximum in PRMan (when using this API call))
        nvars = ctypes.c_int(64)

        fileNameEncoded = fileName.encode("ascii")

        # Just open the file to find out the number of variables in the file...
        # (3Delight only)
        if "3delight" in os.path.basename(libName):
            handle = ptclib.PtcOpenPointCloudFile(fileNameEncoded, ctypes.byref(nvars), None, None)
            if handle is None:
                raise IOError("Cannot open point cloud file %s"%fileName)
            ptclib.PtcClosePointCloudFile(handle)

        # Now prepare storage for the variable names and types and open the file for real...
        numVars = nvars.value
        types = (numVars*ctypes.c_char_p)()
        names = (numVars*ctypes.c_char_p)()
        handle = ptclib.PtcOpenPointCloudFile(fileNameEncoded, ctypes.byref(nvars), types, names)
        if handle is None:
            raise IOError("Cannot open point cloud file %s"%fileName)
        
        self._handle = handle

        self.name = fileName
        # The dictionary containing the point cloud attributes.
        # Access to the attributes is provided via properties.
        self._ptcAttrs = {}
        
        code = ""
        idx = 0
        vars = []
        for i in range(nvars.value):
            name = names[i].decode("ascii")
            type = types[i].decode("ascii")
            vars.append((type, name))
            if type=="float":
                code += "dataDict['%s'] = data[%s]\n"%(name,idx)
                idx += 1
            elif type in ["vector", "point", "normal", "color"]:
                code += "dataDict['%s'] = data[%s:%s]\n"%(name,idx,idx+3)
                idx += 3
            elif type=="matrix":
                code += "dataDict['%s'] = data[%s:%s]\n"%(name,idx,idx+16)
                idx += 16
            else:
                raise RuntimeError("Unknown variable type in point cloud file: %s"%type)
        # Compile the code that will pick the correct data components and put them into a dict
        self._dataCollectionCode = compile(code, "<string>", "exec")
        
        self._ptcAttrs["variables"] = vars
        
        # Get the npoints attribute...
        n = ctypes.c_int()
        ptclib.PtcGetPointCloudInfo.argtypes = [ptclib.PtcPointCloud, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int)]
        if ptclib.PtcGetPointCloudInfo(handle, b"npoints", ctypes.byref(n))==1:
            self._ptcAttrs["npoints"] = n.value
        # Get the datasize attribute
        if ptclib.PtcGetPointCloudInfo(handle, b"datasize", ctypes.byref(n))==1:
            self._ptcAttrs["datasize"] = n.value
        # Get the bbox attribute
        ptclib.PtcGetPointCloudInfo.argtypes = [ptclib.PtcPointCloud, ctypes.c_char_p, ctypes.POINTER(ctypes.c_float)]
        b = (6*ctypes.c_float)()
        if ptclib.PtcGetPointCloudInfo(handle, b"bbox", b)==1:
            self._ptcAttrs["bbox"] = list(b)
        # Get the format attribute
        f = (3*ctypes.c_float)()
        if ptclib.PtcGetPointCloudInfo(handle, b"format", f)==1:
            self._ptcAttrs["format"] = tuple(f)
        # Get the world2eye attribute
        m = (16*ctypes.c_float)()
        if ptclib.PtcGetPointCloudInfo(handle, b"world2eye", m)==1:
            self._ptcAttrs["world2eye"] = list(m)
        # Get the world2ndc attribute
        if ptclib.PtcGetPointCloudInfo(handle, b"world2ndc", m)==1:
            self._ptcAttrs["world2ndc"] = list(m)

        if self.npoints is None:
            raise IOError("Could not obtain the number of points in point cloud file %s."%fileName)
        if self.npoints<0:
            raise ValueError("The number of points in the point cloud file is negative (%s)."%self.npoints)
        if self.datasize is None:
            raise IOError("Could not obtain datasize value from point cloud file %s."%fileName)

        # The number of points that can still be read before eof is hit
        self._numPointsLeft = self.npoints

        # Set up storage for reading individual data points
        self._pos = (3*ctypes.c_float)()
        self._normal = (3*ctypes.c_float)()
        self._radius = ctypes.c_float()
        self._data = (self.datasize*ctypes.c_float)()

    def __del__(self):
        """Destructor.
        
        Closes the file.
        """
        self.close()

    @property
    def variables(self):
        """Return the variables attribute."""
        return self._ptcAttrs.get("variables", None)
    
    @property
    def npoints(self):
        """Return the npoints attribute."""
        return self._ptcAttrs.get("npoints", None)

    @property
    def datasize(self):
        """Return the datasize attribute."""
        return self._ptcAttrs.get("datasize", None)

    @property
    def bbox(self):
        """Return the bbox attribute."""
        return self._ptcAttrs.get("bbox", None)

    @property
    def format(self):
        """Return the format attribute."""
        return self._ptcAttrs.get("format", None)

    @property
    def world2eye(self):
        """Return the world2eye attribute."""
        return self._ptcAttrs.get("world2eye", None)

    @property
    def world2ndc(self):
        """Return the world2ndc attribute."""
        return self._ptcAttrs.get("world2ndc", None)
    
    def iterAttrs(self):
        """Iterate over all attributes defined in the file.
        
        Yields tuples (name,value).
        """
        return list(self._ptcAttrs.items())
    
    def close(self):
        """Close the point cloud file.
        
        This method is also called from the destructor.
        """
        if self._handle is not None:
            self._PtcClosePointCloudFile(self._handle)
            self._handle = None
    
    def readDataPoint(self):
        """Read the next data point.
        
        Returns a tuple (pos, normal, radius, dataDict) where pos and normal
        are 3-tuples of floats, radius is a single float and dataDict a
        dictionary with the extra variables that are attached to the point.
        If no more point is available an EOFError exception is thrown.
        An IOErrror exception is thrown when an error occurs during reading or
        when the file has already been closed.
        """
        if self._handle is None:
            raise IOError("The point cloud file has already been closed (%s)"%self.name)
        if self._numPointsLeft==0:
            raise EOFError("There are no more points left to read from point cloud file %s"%self.name)
        
        self._numPointsLeft -= 1
        res = self._PtcReadDataPoint(self._handle, self._pos, self._normal, self._radius, self._data)
        if res==0:
            raise IOError("Error while reading data point from point cloud file %s"%self.name)
        else:
            dataDict = {}
            data = self._data
            exec(self._dataCollectionCode)
            return tuple(self._pos), tuple(self._normal), self._radius.value, dataDict

    def readDataPoints(self, numPoints, buffer):
        """Read a sequence of data points.
        
        numPoints is the number of points to read. buffer is either a single
        buffer that will receive all values or a tuple (pointbuf, normalbuf,
        radiusbuf, databuf) that contains the individual buffers for the
        respective values. A buffer must always be large enough to hold
        numPoints values.
        The function accepts ctypes arrays as buffers or any sequence object
        that supports the array interface (such as numpy arrays).
        
        The return value is the number of points that have actually
        been read (additional items in the buffers remain at their previous
        value). When 0 is returned, the end of the file has been reached.
        """
        if numPoints<=0:
            return 0
        
        # Are there 4 individual buffers?
        if type(buffer) is tuple:
            if len(buffer)!=4:
                raise ValueError("Expected four individual buffers, but got %s"%len(buffer))
            pbuf,nbuf,rbuf,dbuf = buffer
            pntPtr = _arrayPointer(pbuf, 3*numPoints)
            normPtr = _arrayPointer(nbuf, 3*numPoints)
            radPtr = _arrayPointer(rbuf, numPoints)
            dataPtr = _arrayPointer(dbuf, self.datasize*numPoints)
            pntStride = 3
            normStride = 3
            radStride = 1
            dataStride = self.datasize
        # There is only one single buffer for all values
        else:
            sizeOfFloat = 4
            pntStride = 7+self.datasize
            normStride = pntStride
            radStride = pntStride
            dataStride = pntStride
            pntPtr = _arrayPointer(buffer, pntStride*numPoints)
            normPtr = pntPtr+3*sizeOfFloat
            radPtr = normPtr+3*sizeOfFloat
            dataPtr = radPtr+sizeOfFloat
        
        num = min(numPoints, self._numPointsLeft)
        
        # Read the points
        self._numPointsLeft -= num
        _pointcloud.readDataPoints(ctypes.addressof(self._PtcReadDataPoint), self._handle, num,
                                   pntPtr, pntStride, normPtr, normStride, radPtr, radStride, dataPtr, dataStride)
        return num

    def iterPoints(self):
        """Iterate over all the points in the file.
        
        Yields tuples (point,normal,radius,data) for every point in the file.
        This is equivalent to calling readDataPoint() repeatedly.
        """
        while self._numPointsLeft>0:
            yield self.readDataPoint()
            
    def iterBatches(self, batchSize=1000, combinedBuffer=False, numpyArray=False):
        """Iterate over point batches.
        
        Reads batchSize points at once and yields one or more buffers
        containing the data.
        combinedBuffer determines whether all data is written into one single
        buffer or if there is an individual buffer for the point, normal, radius
        and data.
        numpyArray determines whether the buffers are created as numpy arrays
        or ctypes arrays.
        """
        global _numpy_available
        
        if batchSize<1:
            raise ValueError("Invalid point batch size: %s"%batchSize)
        if numpyArray and not _numpy_available:
            raise ImportError("numpy is not available") 
        
        num = self._numPointsLeft
        buffer = None
        bufLen = 0
        while num>0:
            n = min(batchSize, num)
            if bufLen!=n:
                if combinedBuffer:
                    if numpyArray:
                        buffer = numpy.zeros(shape=(n,7+self.datasize), dtype=numpy.float32)
                    else:
                        buffer = (((7+self.datasize)*n)*ctypes.c_float)()
                else:
                    if numpyArray:
                        ps = numpy.zeros(shape=(n,3), dtype=numpy.float32)
                        ns = numpy.zeros(shape=(n,3), dtype=numpy.float32)
                        rs = numpy.zeros(shape=(n,), dtype=numpy.float32)
                        ds = numpy.zeros(shape=(n,self.datasize), dtype=numpy.float32)
                    else:
                        ps = (3*n*ctypes.c_float)()
                        ns = (3*n*ctypes.c_float)()
                        rs = (n*ctypes.c_float)()
                        ds = (self.datasize*n*ctypes.c_float)()
                    buffer = (ps,ns,rs,ds)
                bufLen = n
            self.readDataPoints(n, buffer)
            yield buffer
            num -= n

    def _loadPtcLib(self, libName):
        """Load a RenderMan library providing the point cloud API.
        """
        resolvedLibName = rmanlibutil.resolveRManLib(libName)
        ptclib = ctypes.cdll.LoadLibrary(resolvedLibName)
        
        ptclib.PtcPointCloud = ctypes.c_void_p
        
        # Reading point cloud files
        ptclib.PtcOpenPointCloudFile.argtypes = [ctypes.c_char_p,
                                                 ctypes.POINTER(ctypes.c_int),
                                                 ctypes.POINTER(ctypes.c_char_p),
                                                 ctypes.POINTER(ctypes.c_char_p)]
        ptclib.PtcOpenPointCloudFile.restype = ptclib.PtcPointCloud
        
        ptclib.PtcGetPointCloudInfo.argtypes = [ptclib.PtcPointCloud, ctypes.c_char_p, ctypes.c_char_p]
        ptclib.PtcGetPointCloudInfo.restype = ctypes.c_int
        
        ptclib.PtcReadDataPoint.argtypes = [ptclib.PtcPointCloud,
                                            ctypes.POINTER(ctypes.c_float),
                                            ctypes.POINTER(ctypes.c_float),
                                            ctypes.POINTER(ctypes.c_float),
                                            ctypes.POINTER(ctypes.c_float)]
        ptclib.PtcReadDataPoint.restype = ctypes.c_int
        
        ptclib.PtcClosePointCloudFile.argtypes = [ptclib.PtcPointCloud]
        ptclib.PtcClosePointCloudFile.restype = None
        
        return ptclib
    
    
class PtcWriter:
    """Point cloud writer class.
    
    An instance of this class is returned by the open() function.
    """
    
    def __init__(self, fileName, libName, vars, world2eye, world2ndc, format):
        """Constructor.
        
        fileName is the name of the point cloud file and libName the name
        of the shared library that implements the point cloud API.
        vars is a list of tuples (type, name) that specifies what additional
        variables to write into the file. world2eye and world2ndc are 4x4
        matrices and format a tuple
        """
        
        self._handle = None
        
        ptclib = self._loadPtcLib(libName)
        
        # Store the functions pointers that are still required
        self._PtcWriteDataPoint = ptclib.PtcWriteDataPoint
        self._PtcFinishPointCloudFile = ptclib.PtcFinishPointCloudFile

        self.name = fileName
        
        xres,yres,aspect = format
        
        w2e = self._matrixToCTypes(world2eye)
        w2n = self._matrixToCTypes(world2ndc)
        
        n = len(vars)
        cvartypes = (n*ctypes.c_char_p)()
        cvarnames = (n*ctypes.c_char_p)()
        idx = 0
        code = ""
        for i in range(n):
            type,name = vars[i]
            cvartypes[i] = type.encode("ascii")
            cvarnames[i] = name.encode("ascii")
            if type=="float":
                code += "cdata[%s] = data.get('%s', 0.0)\n"%(idx, name)
                idx += 1
            elif type in ["vector", "point", "normal", "color"]:
                code += "cdata[%s:%s] = data.get('%s', (0.0,0.0,0.0))\n"%(idx,idx+3,name)
                idx += 3
            elif type=="matrix":
                code += "cdata[%s:%s] = data.get('%s'), (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))\n"%(idx,idx+16,name)
                idx += 16
            else:
                raise RuntimeError("Unknown point cloud variable type: %s"%type)

        self.datasize = idx
        self._dataInitCode = code
        
        cformat = (3*ctypes.c_float)(float(xres), float(yres), float(aspect))
        
        fileNameEncoded = fileName.encode("ascii")
        handle = ptclib.PtcCreatePointCloudFile(fileNameEncoded, n, cvartypes, cvarnames, w2e, w2n, cformat)
        if handle is None:
            raise IOError("Cannot open point cloud file %s for writing"%fileName)

        self._handle = handle
        
    def _matrixToCTypes(self, m):
        """Convert a matrix into a ctypes array.
        
        m can be any object that contains 16 floats (the values may be nested).
        """
        if isinstance(m, cgtypes.mat4):
            values = m.toList(rowmajor=True)
        else:
            values = []
            for v in m:
                try:
                    values.extend(list(v))
                except:
                    values.append(v)
            if len(values)!=16:
                raise ValueError("Matrix must be composed of 16 values, got %s instead."%len(values))
        return (16*ctypes.c_float)(*values)
        
    def __del__(self):
        self.close()
    
    def close(self):
        """Close the point cloud file.
        
        This method is also called from the destructor.
        """
        if self._handle is not None:
            self._PtcFinishPointCloudFile(self._handle)
            self._handle = None

    def writeDataPoint(self, point, normal, radius, data):
        """Write a point into the point cloud file.
        
        point and normal are vectors (any 3-sequence of floats) and radius
        a float. data is a dict that contains the extra variables that
        must have been declared in the constructor. Undeclared values are
        ignored, missing declared values are set to 0.
        """
        if self._handle is None:
            raise IOError("The point cloud file has already been closed.")

        p = (3*ctypes.c_float)(*tuple(point))
        n = (3*ctypes.c_float)(*tuple(normal))
        cdata = (self.datasize*ctypes.c_float)()
        exec(self._dataInitCode)
        res = self._PtcWriteDataPoint(self._handle, p, n, radius, cdata)
        if res==0:
            raise IOError("Failed to write point cloud data point")

    def writeDataPoints(self, numPoints, buffer):
        """Write a sequence of data points.
        
        numPoints is the number of points to write. buffer is either a single
        buffer that contains all values or a tuple (pointbuf, normalbuf,
        radiusbuf, databuf) that each contains the respective value.
        The buffers must contain at least numPoints items.
        The function accepts ctypes arrays as buffers or any sequence object
        that supports the array interface (such as numpy arrays).
        """
        # Are there 4 individual buffers?
        if type(buffer) is tuple:
            if len(buffer)!=4:
                raise ValueError("Expected four individual buffers, but got %s"%len(buffer))
            pbuf,nbuf,rbuf,dbuf = buffer
            pntPtr = _arrayPointer(pbuf, 3*numPoints)
            normPtr = _arrayPointer(nbuf, 3*numPoints)
            radPtr = _arrayPointer(rbuf, numPoints)
            dataPtr = _arrayPointer(dbuf, self.datasize*numPoints)
            pntStride = 3
            normStride = 3
            radStride = 1
            dataStride = self.datasize
        # There is only one single buffer for all values
        else:
            sizeOfFloat = 4
            pntStride = 7+self.datasize
            normStride = pntStride
            radStride = pntStride
            dataStride = pntStride
            pntPtr = _arrayPointer(buffer, pntStride*numPoints)
            normPtr = pntPtr+3*sizeOfFloat
            radPtr = normPtr+3*sizeOfFloat
            dataPtr = radPtr+sizeOfFloat
        
        # Write the points
        _pointcloud.writeDataPoints(ctypes.addressof(self._PtcWriteDataPoint), self._handle, numPoints,
                                    pntPtr, pntStride, normPtr, normStride, radPtr, radStride, dataPtr, dataStride)

    def _loadPtcLib(self, libName):
        """Load a RenderMan library providing the point cloud API.
        """
        resolvedLibName = rmanlibutil.resolveRManLib(libName)
        ptclib = ctypes.cdll.LoadLibrary(resolvedLibName)
        
        ptclib.PtcPointCloud = ctypes.c_void_p

        # Writing point cloud files
        ptclib.PtcCreatePointCloudFile.argtypes = [ctypes.c_char_p,
                                                   ctypes.c_int,
                                                   ctypes.POINTER(ctypes.c_char_p),
                                                   ctypes.POINTER(ctypes.c_char_p),
                                                   ctypes.POINTER(ctypes.c_float),
                                                   ctypes.POINTER(ctypes.c_float),
                                                   ctypes.POINTER(ctypes.c_float)]
        ptclib.PtcCreatePointCloudFile.restype = ptclib.PtcPointCloud
        
        ptclib.PtcWriteDataPoint.argtypes = [ptclib.PtcPointCloud,
                                             ctypes.POINTER(ctypes.c_float),
                                             ctypes.POINTER(ctypes.c_float),
                                             ctypes.c_float,
                                             ctypes.POINTER(ctypes.c_float)]
        ptclib.PtcWriteDataPoint.restype = ctypes.c_int
        
        ptclib.PtcFinishPointCloudFile.argtypes = [ptclib.PtcPointCloud]
        ptclib.PtcFinishPointCloudFile.restype = None
        
        return ptclib


def open(fileName, mode="r", libName=None, **kwargs):
    """Open a point cloud file for reading or writing.
    
    *fileName* is the name of the point cloud file. *mode* is either ``"r"``
    for reading a file or ``"w"`` for writing a new point cloud file.
    *libName* is the library name that implements the point cloud API.
    When mode is ``"w"``, the following additional keyword arguments must
    be present:
    
    - ``vars``: A list of tuples (*type*, *name*) that defines what additional variables to write
    - ``world2eye``: The world2eye matrix 
    - ``world2ndc``: The world2ndc matrix
    - ``format``: A tuple (*xres*, *yres*, *aspect*)
    
    Depending on the mode, the function either returns a :class:`PtcReader` or
    :class:`PtcWriter` object.
    """
    if mode=="r":
        return PtcReader(fileName, libName, **kwargs)
    elif mode=="w":
        return PtcWriter(fileName, libName=libName, **kwargs)
    else:
        raise ValueError('Invalid file mode: "%s" (expected "r" or "w")'%mode)

###################################################################

if __name__=="__main__":
    from cgkit.cgtypes import *
    import random
    
    if 1:
        print("---WRITE---")
        ptc = open("test.ptc", "w", "3delight", vars=[("float", "spam"), ("vector", "dir")], world2eye=None, world2ndc=None, format=(320,240,1.333))
        for i in range(100):
            x = random.random()
            y = random.random()
            z = random.random()
            ptc.writeDataPoint((x,y,z), (0,1,0), 0.2, {"spam":0.5})
        ptc.close()
    
    print("----READ----")
    rd = open("test.ptc", "r", "3delight")
    print((rd.variables))
    print(("npoints:",rd.npoints))
    print(("datasize",rd.datasize))
    print(("bbox",rd.bbox))
    print(("format",rd.format))
    print(("world2eye",rd.world2eye))
    print(("world2ndc",rd.world2ndc))
    p,n,r,data = rd.readDataPoint()
    print((p,n,r,data))
    print((rd.readDataPoint()))
