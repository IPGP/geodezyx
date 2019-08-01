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
# $Id: call.py,v 1.1.1.1 2004/12/12 14:31:33 mbaas Exp $

## \file call.py
## Contains the Call class.

import string, types


# Call
class Call:
    """Call adapter.

    This class can either be used to wrap a string or a callable
    object with a fixed set of arguments (positional and keyword
    arguments).  When the %Call object is called (without any
    arguments) it executes the string or calls the wrapped function
    with the stored arguments.

    Objects of this class can be compared using == or != and they're
    hashable (so can be used in dictionaries, etc.)

    The purpose of this class is to provide an adapter for GUI commands
    which take no arguments.

    Example:

    \code

    def foo(a, b, c=0):
        print "a:%s  b:%s  c:%s"%(a,b,c)

    >>> f = Call(foo, 1, 2, c=3)
    
    # Call foo with the specified arguments
    >>> f()
    a:1  b:2  c:3
    
    # Print a readable description of the command object
    >>> print f
    foo(1, 2, c=3)

    \endcode

    """
    
    def __init__(self, func, *args, **keyargs):
        """Constructor.

        You must provide the callable object as first argument. All
        additional arguments will be passed to the callable object
        when its called.

        If \a func is a string then there may be one additional argument
        that specifies a dictionary which should be used as global
        (and local) namespace. The dictionary may be specified as positional
        argument following the string or as keyword argument "globals".

        The constructor sets the variable __doc__ to a reasonable value.
        It'll be either the function doc string if there is one or the
        command string itself.

        \param func (\c callable or \c str) The function that should be called
        \param *args  The positional arguments (default = no args)
        \param **keyargs  The keyword arguments (default = no keyword args)
        """
        if isinstance(func, (str,)):
            self._func = func
            self._args = None
            self._keyargs = None
            compilestr = func
            if "globals" in keyargs:
                self._globals = keyargs["globals"]
            elif len(args)==1:
                self._globals = args[0]
            self.__doc__ = None
        elif callable(func):
            self._func = func
            self._args = args
            self._keyargs = keyargs
            compilestr = "apply(self._func, self._args, self._keyargs)"
            self.__doc__ = getattr(func, "__doc__", None)
        else:
            raise TypeError("The 'func' argument must be a string or a callable object.")

        self._code = compile(compilestr, "<command: '%s'>"%str(self), "exec")
        if self.__doc__==None or self.__doc__=="":
            self.__doc__ = str(self)

    def __str__(self):
        if callable(self._func):
            res = getattr(self._func, "__name__", "<unknown>")
            args = [str(x) for x in self._args]
            args += ["%s=%s"%(x,self._keyargs[x]) for x in self._keyargs]
            res = "%s(%s)"%(res,string.join(args,", "))
        else:
            res = self._func
        return res

    def __repr__(self):
        return "<Call cmd='%s'>"%str(self)

    def __call__(self):
        if hasattr(self, "_globals"):
            exec (self._code, self._globals)
        else:
            exec (self._code)

    def __eq__(self, other):
        if not isinstance(other, Call):
            return False

        return (self._func==other._func and
                self._args==other._args and
                self._keyargs==other._keyargs)

    def __ne__(self, other):
        if not isinstance(other, Call):
            return True

        return (self._func!=other._func or
                self._args!=other._args or
                self._keyargs!=other._keyargs)

    def __hash__(self):
        return hash(self._code)
#        return hash(self._func) ^ hash(self._args)


