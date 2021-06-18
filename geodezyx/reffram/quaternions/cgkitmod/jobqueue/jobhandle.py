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

import os, os.path, glob


class JobHandle(object):
    """Job object (queued).
    
    Objects of this class represent jobs that have already been queued.
    
    Every job has:
    
    - A type name
    - An instance name
    - A job description
    - An associated directory
    """
    
    def __init__(self, location, rootLocation):
        """Constructor.
        
        location is the path to the job directory (this may also point to a
        link file).
        rootLocation is the path to the job root directory (this is used
        for resolving job relative relative links).
        """
        object.__init__(self)
        # The job location (if useSymLinks is False, this may also point to a link file)
        self._location = str(location)

        # The job root directory
        self._rootLocation = rootLocation
        
        # The resolved location. This is always a directory (or technically, it
        # could be a non-existent path). If useSymLinks is True, this is
        # always the same as the above location.
        self._realLocation = self._resolveTxtLinks(self._location)
        
        
    def _resolveTxtLinks(self, location):
        """Resolve a path so that it points to a directory and not to a link file.
        
        If location already points to a directory, it is returned unmodified.
        If it points to a link file, the file is followed and the process
        repeats until a directory is encountered or the targegt does not exist.
        """
        visited = {}
        while os.path.exists(location) and os.path.isfile(location):
            visited[location] = True
            location = self._getTxtLinkTarget(location)
            if location in visited:
                raise ValueError("Cyclic directory links detected")
        return location
    
    def _getTxtLinkTarget(self, linkFileName):
        """Return the target that a link file points to.
        """
        f = open(linkFileName, "rt")
        # Read the first line which must be "[link]"
        s = f.readline()
        if s!="[link]\n":
            raise ValueError("%s is not a link file"%linkFileName)
        # Read the line containing the target
        targetLine = f.readline()
        f.close()
        # Check the target line
        if not targetLine.startswith("target="):
            raise ValueError("Error in link file %s"%linkFileName)
        target = targetLine[7:].strip()
        # The link must always be relative to the job root
        if not target.startswith("$JOBROOT"):
            raise ValueError("Invalid link target: %s"%target)
        target = self._rootLocation+target[8:]
        return target
        
    def __str__(self):
        if self.isFinished():
            state = "finished"
        elif self.isRunning():
            state = "running"
        elif self.isWaiting():
            state = "waiting"
        else:
            state = "?"
            
        res = "%s (%s)"%(os.path.basename(self._location), state)
        return res
        
    @property
    def location(self):
        """Return the job directory.
        """
        return self._location
    
    @property
    def label(self):
        """Return a short job label to quickly identify the job.
        """
        try:
            try:
                f = open(self.labelFile, "rt")
                label = f.read().strip()
            finally:
                f.close()
        except:
            label = ""

        if label=="":
            label = "Job %s"%self.number
        
        return label

    @property
    def number(self):
        """Return the job number.
        """
        nrStr = os.path.basename(self._location)[3:]
        try:
            return int(nrStr)
        except:
            return -1
        
    @property
    def submitTime(self):
        """Return the submission time (in seconds).
        
        May return None if the job doesn't exist or is broken.
        """
        try:
            s = os.stat(self.procDefFile)
            return s.st_mtime
        except OSError:
            return None

    @property
    def startTime(self):
        """Return the time the job was started (in seconds).
        
        Returns None if the job hasn't been started yet (or there is a problem
        with the job directory).
        """
        try:
            s = os.stat(self.pidFile)
            return s.st_mtime
        except OSError:
            return None

    @property
    def endTime(self):
        """Return the time the job was finished (in seconds).
        
        Returns None if the job hasn't been finished yet (or there is a problem
        with the job directory).
        """
        try:
            s = os.stat(self.finishedDir)
            return s.st_mtime
        except OSError:
            return None
        
    @property
    def progress(self):
        """Return the progress percentage value as an int.
        
        Returns None if the value couldn't be obtained for some reason.
        """
        try:
            progressFile = self.progressFile
            if os.path.exists(progressFile):
                res = os.path.getsize(self.progressFile)
                if res>100:
                    res = 100
            else:
                res = 0
            return res
        except:
            return None
    
    @property
    def statusLine(self):
        """Return a string containing the current status line for the GUI.
        
        The status line indicates what a running job is currently doing.
        """
        try:
            try:
                f = open(self.statusLineFile, "rt")
                line = f.read().strip()
            finally:
                f.close()
        except:
            line = ""

        return line
        

    @property
    def finishedDir(self):
        return os.path.join(self._realLocation, ".finished")
    
    @property
    def runningDir(self):
        return os.path.join(self._realLocation, ".running")

    @property
    def procDefFile(self):
        return os.path.join(self._realLocation, ".proc_def")

    @property
    def labelFile(self):
        """Return the location of the file that contains the job label.
        """
        return os.path.join(self._realLocation, ".label")

    @property
    def pidFile(self):
        """Return the location of the PID file.
        
        This file contains the PID of the process that is/was running the job.
        """
        return os.path.join(self.runningDir, ".pid")

    @property
    def hostFile(self):
        """Return the location of the host file.
        
        The host file contains the name of the host that is/was running the job.
        """
        return os.path.join(self.runningDir, ".host")

    @property
    def progressFile(self):
        """Return the location of the progress indicator file.
        """
        return os.path.join(self.runningDir, ".progress")

    @property
    def statusLineFile(self):
        """Return the location of the status line file.
        """
        return os.path.join(self.runningDir, ".statusline")

    @property
    def errorMarkerFile(self):
        """Return the location of the error marker file.
        
        The presence of this file indicates that running the job resulted
        in an error.
        """
        return os.path.join(self.runningDir, ".error_marker")

    @property
    def procTracebackFile(self):
        """Return the location of the file that contains the job procedure traceback.
        """
        return os.path.join(self.runningDir, ".proc_traceback")

    @property
    def stdoutFile(self):
        """Return the location of the file that contains stdout.
        """
        return os.path.join(self.runningDir, ".stdout")

    @property
    def stderrFile(self):
        """Return the location of the file that contains stderr.
        """
        return os.path.join(self.runningDir, ".stderr")
    
    def isWaiting(self):
        """Check if this job is currently in 'waiting' state.
        """
        if os.path.exists(self.finishedDir):
            return False
        if os.path.exists(self.runningDir):
            return False
        return os.path.exists(self.procDefFile)
    
    def isRunning(self):
        """Check if this job is currently in 'running' state.
        """
        if self.isFinished():
            return False
        else:
            return os.path.exists(self.runningDir)
    
    def isFinished(self):
        """Check if this job is currently in 'finished' state.
        """
        return os.path.exists(self.finishedDir)
    
    def listSubJobs(self):
        """Return a list of sub-jobs.
        
        The return value is a list of JobHandle objects.
        """
        subDirs = glob.glob(os.path.join(self._realLocation, "job*"))
        #subDirs = filter(lambda p: os.path.isdir(p), subDirs)
        jobs = [(jobDir, int(os.path.basename(jobDir)[3:])) for jobDir in subDirs]
        jobs.sort(key=lambda a: a[1])
        return [JobHandle(jobDir, self._rootLocation) for jobDir,nr in jobs]
    
    def hasError(self, recursive=False):
        """Check if this job produced an error.
        
        If recursive is False, only this job is considered, not the children
        jobs. In this case, the return value is only meaningful when the job
        is in finished state.
        If recursive is True, the result will be True if any job in the entire
        sub-hierarchy has an error. This can be an expensive operation because
        all sub-directories have to be checked.
        """
        err = os.path.exists(self.errorMarkerFile)
        if err:
            return True
        
        if recursive:
            subJobs = self.listSubJobs()
            for j in subJobs:
                if j.hasError(recursive=True):
                    return True
            
        return False
    
    def setError(self):
        """Set the error flag.
        
        This method may only be called by the process that is currently
        running the job (the .running directory must exist).
        """
        errFile = self.errorMarkerFile
        if not os.path.exists(errFile):
            f = open(errFile, "wb")
            f.close()
