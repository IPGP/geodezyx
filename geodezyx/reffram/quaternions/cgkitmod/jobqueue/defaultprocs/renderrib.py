# Render a RIB file

import os, subprocess
import cgkit.jobqueue

class renderrib(cgkit.jobqueue.JobProc):
    
    def __init__(self, rib, renderer="aqsis"):
        cgkit.jobqueue.JobProc.__init__(self)
        self._rib = rib
        self._renderer = renderer.lower()
        
    def run(self):
        args = [self._rendererExecutable(), self._rib]
        cmd = " ".join(args)
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                shell=True)
        outdata,errdata = proc.communicate()
        outdata = outdata.rstrip()
        errdata = errdata.rstrip()
        print("OUT")
        print(outdata)
        print("ERR")
        print(errdata)
        print("RET")
        print((proc.returncode))
        
        if proc.returncode!=0 or errdata!="":
            self.setError()
        
        
    def _rendererExecutable(self):
        """Return the name of the renderer executable.
        """
        if self._renderer=="aqsis":
            return "aqsis"
        elif self._renderer=="pixie":
            return "rndr"
        elif self._renderer=="3delight":
            return "renderdl"
        elif self._renderer=="prman":
            return "prman"
        else:
            raise ValueError("Unknown renderer: %s"%self._renderer)
