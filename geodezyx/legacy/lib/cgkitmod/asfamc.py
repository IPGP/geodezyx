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
# $Id: asfamc.py,v 1.3 2005/02/06 22:24:27 mbaas Exp $

## \file asfamc.py
## Contains the ASFReader and AMCReader class.

import string

# ReaderBase
class ReaderBase:
    """Base class for the ASFReader and AMCReader classes.

    A derived class must implement readXyz() method where Xyz is
    the section name (the keyword in the file is ':xyz').
    If such a read method is not implemented, the corresponding
    section is ignored.
    """

    def __init__(self, filename):

        self.filename = filename
        # A list of unprocessed tokens (strings)
        self.tokenlist = []
        # The last returned token (string)
        self.lasttok = None
        self.lastline = None
        self.undoline = None
        # The current line number
        self.linenr = 0

    # read
    def read(self):
        """Read the entire file.
        """
        self.fhandle = file(self.filename)

        while 1:
            try:
                tok = self.token()
            except StopIteration:
                break

            # Check if the keyword starts with ':' (a new section)
            if tok[0]!=":":
                raise SyntaxError("Syntax error in line %d: Invalid keyword '%s'"%(self.linenr, tok))

            # Read the current section...

            # Determine the method name (only letters and digits allowed,
            # everything else will be removed)
            mn = tok[1:].capitalize()
            methodname = ""
            for c in mn:
                if c not in string.ascii_letters and c not in string.digits:
                    continue
                methodname += c

            handler = getattr(self, "read"+methodname, None)
            if handler!=None:
                handler()
            else:
                try:
                    self.skipSection()
                except:
                    break        

    # intToken
    def intToken(self):
        """Return the next token which must be an int.
        """

        tok = self.token()
        try:
            return int(tok)
        except ValueError:
            raise SyntaxError("Syntax error in line %d: Integer expected, got '%s' instead"%(self.linenr, tok))

    # floatToken
    def floatToken(self):
        """Return the next token which must be a float.
        """

        tok = self.token()
        try:
            return float(tok)
        except ValueError:
            raise SyntaxError("Syntax error in line %d: Float expected, got '%s' instead"%(self.linenr, tok))

    # token
    def token(self):
        """Return the next token."""

        # Are there still some tokens left? then just return the next one
        if self.tokenlist!=[]:
            tok = self.tokenlist[0]
            self.lasttok = tok
            self.tokenlist = self.tokenlist[1:]
            return tok

        # Read a new line
        s = self.readLine()
        self.createTokens(s)
        return self.token()

    # undoToken
    def undoToken(self):
        """Put back the last token.
        """
        self.tokenlist = [self.lasttok]+self.tokenlist
        self.lasttok = None

    # undoLine
    def undoLine(self):
        """Put back the last line.
        """
        self.undoline = self.lastline
        self.lastline = None
        self.tokenlist = []

    # skipSection
    def skipSection(self):
        """Read and ignore data until the next section is reached.

        The next call of token() will result in a ':' keyword.
        """
        while 1:
            s = self.readLine().strip()
            if s[0]==':':
                self.createTokens(s)
                return
        

    # readLine
    def readLine(self):
        """Return the next line containing data.

        Empty lines and comments are skipped. If the end of the file
        has been reached, a StopIteration exception is thrown.
        The return value is the next line containing data (this will
        never be an empty string).
        """
        # Discard any remaining tokens
        self.tokenlist = []

        # Is there a line that was undone?
        if self.undoline!=None:
            self.lastline = self.undoline
            self.undoline = None
            return self.lastline
        
        # Read the next line
        while 1:
            s = self.fhandle.readline()
            self.linenr += 1
            if s=="":
                raise StopIteration()
            if s[0]!="#":
                self.lastline = s
                return s

    # createTokens
    def createTokens(self, s):
        """Populate the token list from the content of s.
        """
        s = s.strip()
        s = s.replace("(", " ")
        s = s.replace(")", " ")
        s = s.replace(",", " ")
        a = s.split()
        self.tokenlist = a

######################################################################

# ASFReader
class ASFReader(ReaderBase):
    """Read Acclaim Skeleton File Definition (ASF) files.

    This class serves as base class for reading ASF files.
    Derived classes may implement the following methods:

    - onVersion(ver)
    - onName(name)
    - onUnits(units)
    - onDocumentation(doc)
    - onRoot(data)
    - onBonedata(bones)
    - onHierarchy(links)
    """

    def __init__(self, filename):
        ReaderBase.__init__(self, filename)

    def onVersion(self, ver): pass
    def onName(self, name): pass
    def onUnits(self, units): pass
    def onDocumentation(self, doc): pass
    def onRoot(self, data): pass
    def onBonedata(self, bones): pass
    def onHierarchy(self, links): pass

    # readVersion
    def readVersion(self):
        ver = self.floatToken()
        self.onVersion(ver)

    # readName
    def readName(self):
        name = self.token()
        self.onName(name)

    # readUnits
    def readUnits(self):
        units = {}
        while 1:
            s = self.readLine()
            if s[0]==':':
                self.undoLine()
                self.onUnits(units)
                return
            a = s.strip().split()
            # Get the argument
            if len(a)>1:
                val = a[1]
            # Try to convert to float. If it doesn't work, use the string value
            try:
                val = float(val)
            except ValueError:
                pass
            units[a[0]] = val

    # readDocumentation
    def readDocumentation(self):
        lines = []
        while 1:
            s = self.readLine()
            if s[0]==':':
                self.undoLine()
                doc = "\n".join(lines)
                self.onDocumentation(doc)
                return
            lines += [s[:-1]]

    # readRoot
    def readRoot(self):
        data = {}
        while 1:
            s = self.readLine()
            if s[0]==':':
                self.undoLine()
                self.onRoot(data)
                return
            a = s.strip().split()
            val = tuple(a[1:])
            data[a[0]] = val
    
    # readBonedata
    def readBonedata(self):
        bones = []
        while 1:
            stop = False
            try:
                # Get the next token
                begintok = self.token()
            except StopIteration:
                stop = True

            # Is it already the begin of the next section?
            if begintok[0]==':':
                self.undoToken()
                stop = True

            # End of file or next section? -> no more bone data
            if stop:
                self.onBonedata(bones)
                return
            
            if begintok!="begin":
                raise SyntaxError("Syntax error in line %d: 'begin' expected, got '%s'"%(self.linenr, begintok))

            data = {}

            while 1:
                s = self.readLine().strip()
                if s=='end':
                    break
                a = s.strip().split()
                # Check for 'limits' which is a special case as it might
                # span several lines.
                if a[0]=="limits":
                    if "dof" not in data:
                        raise ValueError("Line %d: 'dof' settings must appear before the 'limits' settings"%(self.linenr))
                    dof = data["dof"]
                    # Put back the line and use the token mechanism
                    self.undoLine()
                    # Read the 'limits' token
                    tok = self.token()
                    limits = []
                    for i in range(len(dof)):
                        minval = self.token()
                        maxval = self.token()
                        if minval.lower()!="-inf":
                            try:
                                minval = float(minval)
                            except ValueError:
                                raise SyntaxError("Syntax error in line %d: Float or '-inf' expected, got '%s'"%(self.linenr, minval))
                        if maxval.lower()!="inf":
                            try:
                                maxval = float(maxval)
                            except ValueError:
                                raise SyntaxError("Syntax error in line %d: Float or 'inf' expected, got '%s'"%(self.linenr, maxval))
                        limits.append((minval,maxval))
                    data["limits"] = limits
                else:
                    val = tuple(a[1:])
                    data[a[0]] = val

            bones.append(data)

    # readHierarchy
    def readHierarchy(self):
        links = []
        while 1:
            stop = False
            try:
                # Get the next token
                begintok = self.token()
            except StopIteration:
                stop = True
            
            # Is it already the begin of the next section?
            if begintok[0]==':':
                self.undoToken()
                stop = True

            # End of file or next section? -> no more bone data
            if stop:
                self.onHierarchy(links)
                return

            if begintok!="begin":
                raise SyntaxError("Syntax error in line %d: 'begin' expected, got '%s'"%(self.linenr, begintok))

            while 1:
                s = self.readLine().strip()
                if s=='end':
                    break
                a = s.strip().split()
                links.append((a[0], a[1:]))

# AMCReader
class AMCReader(ReaderBase):
    """Read Acclaim Motion Capture Data (AMC) files.

    This class serves as base class for reading AMC files.
    Derived classes have to implement the following method:

    - onFrame(framenr, data)

    """

    def __init__(self, filename):
        ReaderBase.__init__(self, filename)

    def onFrame(self, framenr, data):
        pass
#        print framenr, data

    # read
    def read(self):
        self.fhandle = file(self.filename)

        while 1:
            try:
                tok = self.token()
            except StopIteration:
                return

            # Check if the keyword starts with ':'. If it does, ignore it
            if tok[0]!=":":
                self.undoToken()
                break

        # Read the motion data...
        while 1:
            try:
                framenr = self.intToken()
            except StopIteration:
                return

            data = []
            while 1:
                try:
                    s = self.readLine()
                except StopIteration:
                    self.onFrame(framenr, data)
                    return
                
                a = s.strip().split()
                # Check if it was a frame number
                try:
                    int(a[0])
                    self.undoLine()
                    self.onFrame(framenr, data)
                    break
                except ValueError:
                    pass
                bone = a[0]
                values = [float(x) for x in a[1:]]
                data.append((bone, values))


