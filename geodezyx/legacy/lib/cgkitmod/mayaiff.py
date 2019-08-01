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

"""This module contains the IFFReader base class to parse Maya IFF files.
"""

import struct, os.path

class Chunk:
    """Chunk class.
    
    This class stores information about a chunk. It is passed to the handler
    methods which can also use this class to read the actual chunk data.
    
    The class has the following attributes:
    
    - tag: The chunk name (four characters)
    - size: The size in bytes of the data part of the chunk
    - pos: The absolute position of the data part within the input file
    - parent: The GroupChunk object of the parent chunk
    - depth: The depth of the node (i.e. how deep it is nested). The root
             has a depth of 0.

    The binary chunk data can be read by using the read() method. You
    can specify the number of bytes to read or read the entire data
    at once. It is not possible to read data that lies outside this
    chunk.
    """
    
    def __init__(self, file, parent, tag, size, pos, depth):
        """Constructor.
        """
        # The file handle that currently points to the start of the chunk data
        self.file = file
        # The parent chunk object
        self.parent = parent
        # The chunk name (four characters)
        self.tag = tag.decode("ascii")
        # The chunk size in bytes (only the data part)
        self.size = size
        # The absolute position (of the data) within the input file
        self.pos = pos
        # The depth of the node (i.e. how deep it is nested)
        self.depth = depth
        
        # The number of bytes read so far
        self._bytesRead = 0

    def __str__(self):
        return "<Chunk %s at pos %d (%d bytes)>"%(self.tag, self.pos, self.size)

    def chunkPath(self):
        """Return the full path to this chunk.

        The result is a concatenation of all chunk names that lead to this chunk.
        """
        name = "%s"%(self.tag)
        if self.parent is None:
            return name
        else:
            return self.parent.chunkPath()+".%s"%name

    def read(self, bytes=-1):
        """Read the specified number of bytes from the chunk.
        
        If bytes is -1 the entire chunk data is read.
        """
        if self.file is None:
            raise RuntimeError("This chunk is not active anymore")

        maxbytes = self.size-self._bytesRead
        if bytes<0:
            bytes = maxbytes
        else:
            bytes = min(bytes, maxbytes)
        self._bytesRead += bytes
        return self.file.read(bytes)

  
class GroupChunk(Chunk):
    """Specialized group chunk class.

    In addition to the Chunk class this class has an attribute "type"
    that contains the group type (the first four characters of the data part).
    """
    
    def __init__(self, file, parent, tag, size, pos, type, depth):
        Chunk.__init__(self, file=file, parent=parent, tag=tag, size=size, pos=pos, depth=depth)

        # The group type
        self.type = type.decode("ascii")

        # Start with 4 because the type was already read
        self._bytesRead = 4
        
    def __str__(self):
        return "<GroupChunk %s (%s) at pos %d (%d bytes)>"%(self.tag, self.type, self.pos, self.size)
 
    def chunkPath(self):
        """Return the full path to this chunk.
        """
        name = "%s[%s]"%(self.tag, self.type)
        if self.parent is None:
            return name
        else:
            return self.parent.chunkPath()+".%s"%name
   

class IFFReader:
    """Read Maya IFF files and call an appropriate handler for each chunk.
    
    This is the base class for any Maya IFF file reader class. Derived classes
    can implement the chunk handler methods onBeginGroup(), onEndGroup() and
    onDataChunk(). These handlers receive a Chunk (or GroupChunk) object as
    input that contains information about the current chunk and that can also
    be used to read the actual chunk data.
    """
    
    def __init__(self, iffType=None):
        """Constructor.
        
        iffType may be a four-character string defining the type of IFF files
        the class is supposed to handle (it's the type of the first group chunk).
        If None is passed, any IFF file will be processed. For example, when set
        to "Maya", the read() method will raise an error if the input file is
        not a Maya binary file.
        """
        # If the type is given as a string, turn it into a bytes object
        if type(iffType) is str:
            iffType = iffType.encode("utf-8")
        self._iffType = iffType

    def abort(self):
        """Aborts reading the current file.
        
        This method can be called in handler functions to abort reading
        the file.
        """
        self._abortFlag = True

    def onBeginGroup(self, chunk):
        """Callback that is called whenever a new group tag begins.
        
        chunk is a GroupChunk object containing information about the group chunk.
        """
        pass
        
    def onEndGroup(self, chunk):
        """Callback that is called whenever a group goes out of scope.
        
        chunk is a GroupChunk object containing information about the group chunk
        (it is the same instance that was passed to onBeginGroup()).
        """
        pass
        
    def onDataChunk(self, chunk):
        """Callback that is called for each data chunk.

        chunk is a Chunk object that contains information about the chunk
        and that can be used to read the actual chunk data.
        """
        pass

    def read(self, file):
        """Read the IFF file.
        
        This method reads all chunks sequentially and invokes appropriate
        callback methods. 
        file is a file-like object or the name of a file.
        """
        if type(file) is str:
            self.filename = file
            file = open(file, "rb")
            closeFile = True
        else:
            self.filename = getattr(file, "name", "?")
            closeFile = False
        
        self._file = file
        self._abortFlag = False
       
        try:
            self._read(file)
        finally:
            if closeFile:
                file.close()

    def _read(self, file):
        """Actual read() implementation.
        
        file must be a file handle.
        """ 
        # Check if this actually is a Maya IFF file
        # (and that it starts with a group tag)
        header = file.read(12)
        file.seek(0)
        if len(header)<8 or header[0:4]!=b"FOR4":# or header[8:12]!="Maya":
            raise ValueError('The file "%s" is not a Maya IFF file.'%self.filename)
        if self._iffType is not None and header[8:12]!=self._iffType:
            raise ValueError('The file "%s" is not a %s file.'%(self.filename, self._iffType))
        
        # The current byte position inside the file
        pos = 0
        # A stack with alignment values. Each group tag pushes a new value
        # which is popped when the group goes out of scope
        alignments = []
        # A stack with the currently open group chunks. The items are
        # 2-tuples (endpos, groupchunk).
        pendingGroups = []
        # The current depth of the chunks
        depth = 0
        while not self._abortFlag:
            tag,size = self.readChunkHeader()
            if tag is None:
                break
            pos += 8
            if self.isGroupChunk(tag):
                # Group chunk...
                chtype = file.read(4)
                if len(pendingGroups)==0:
                    parent = None
                else:
                    parent = pendingGroups[-1][1]
                chunk = GroupChunk(file=file, parent=parent, tag=tag, size=size, pos=pos, type=chtype, depth=depth) 
                self.onBeginGroup(chunk)
                av = self.alignmentValue(tag)
                alignments.append(av)
                end = pos+self.paddedSize(size, av)
                if len(pendingGroups)>0 and end>pendingGroups[-1][0]:
                    raise ValueError('Chunk %s at position %s in file "%s" has an invalid size (%d) that goes beyond its contained group chunk.'%(tag,pos-8,os.path.basename(self.filename),size))
                pendingGroups.append((end, chunk))
                pos += 4
                depth += 1
            else:
                # Data chunk...
                chunk = Chunk(file=file, parent=pendingGroups[-1][1], tag=tag, size=size, pos=pos, depth=depth)
                self.onDataChunk(chunk)
                pos += self.paddedSize(size, alignments[-1])

            # Check which groups are to be closed...
            while len(pendingGroups)>0 and pos>=pendingGroups[-1][0]:
                end,chunk = pendingGroups.pop()
                self.onEndGroup(chunk)
                depth -= 1
                
            # Seek to the next chunk position. This is done here (even though it
            # wouldn't be necessary in some cases) so that the callbacks have
            # no chance to mess with the file handle and bring the reader
            # out of sync.
            file.seek(pos)

    def readChunkHeader(self):
        """Read the tag and size of the next chunk.
        
        Returns a tuple (tag, size) where tag is the four character
        chunk name and size is an integer containing the size of the
        data part of the chunk.
        Returns None,None if the end of the file has been reached.
        Throws an exception when an incomplete tag/size was read.
        """
        header = self._file.read(8)
        if len(header)==0:
            return None,None
        if len(header)!=8:
            raise ValueError('Premature end of file "%s" (chunk tag & size expected)'%os.path.basename(self.filename))
        return (header[:4], struct.unpack(">L", header[4:])[0])
        
    def isGroupChunk(self, tag):
        """Check if the given tag refers to a group chunk.

        tag is the chunk name. Returns True when tag is the name
        of a group chunk.
        """
        return tag in [b"FORM", b"CAT ", b"LIST", b"PROP",
                       b"FOR4", b"CAT4", b"LIS4", b"PRO4",
                       b"FOR8", b"CAT8", b"LIS8", b"PRO8"]
    
    def alignmentValue(self, tag):
        """Return the alignment value for a group chunk.
        
        Returns 2, 4 or 8.
        """
        if tag in [b"FORM", b"CAT ", b"LIST", b"PROP"]:
            return 2
        elif tag in [b"FOR4", b"CAT4", b"LIS4", b"PRO4"]:
            return 4
        elif tag in [b"FOR8", b"CAT8", b"LIS8", b"PRO8"]:
            return 8
        else:
            return 2
        
    def paddedSize(self, size, alignment):
        """Return the padded size that is aligned to the given value.
        
        size is an arbitrary chunk size and alignment an integer
        containing the alignment value (2,4,8). If size is already
        aligned it is returned unchanged, otherwise an appropriate
        number of padding bytes is added and the aligned size is
        returned.
        """
        # Padding required?
        if size%alignment!=0:
            padding = alignment-size%alignment
            size += padding
        return size
      
    @staticmethod
    def dump(buf):
        """Helper method to do a hex dump of chunk data.
        
        buf is a string containing the data to dump.
        """
        offset = 0
        while len(buf)>0:
            data = buf[:16]
            buf = buf[16:]
            s = "%04x: "%offset
            s += " ".join(["%02x"%ord(x) for x in data])
            s += (55-len(s))*' '
            for c in data:
                if ord(c)<32:
                    c = '.'
                s += c
            print(s)
            offset += 16
            
            
if __name__=="__main__":
    
    import sys
    
    class IFFDumper(IFFReader):
        def onBeginGroup(self, chunk):
            indent = chunk.depth*"\t"
            print(("%s%s (%s bytes)"%(indent, chunk.type, chunk.size)))
        
        def onEndGroup(self, chunk):
            pass
            #print ("GRP END %s"%chunk)
        
        def onDataChunk(self, chunk):
            indent = chunk.depth*"\t"
            print(("%s%s (%s bytes)"%(indent, chunk.tag, chunk.size)))
#            if chunk.size<=64:
#                data = chunk.read()
#                self.dump(data)
            
    rd = IFFDumper()
    rd.read(open(sys.argv[1], "rb"))
 