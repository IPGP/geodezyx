# Test job to dump the environment to stdout

import os, socket, getpass
import cgkit.jobqueue

class dumpenv(cgkit.jobqueue.JobProc):
    """Render a Blender file.
    """
    
    def __init__(self):
        cgkit.jobqueue.JobProc.__init__(self, label="Dump Environment")
                 
    def run(self):
        print(("Host name   : %s"%socket.gethostname()))
        print(("User name   : %s"%getpass.getuser()))
        print(("Current dir : %s"%os.getcwd()))
        print ("\nEnvironment variables:\n")
        vars = list(os.environ.keys())
        for var in sorted(vars):
            print(("  %s = %s"%(var, os.environ.get(var))))
