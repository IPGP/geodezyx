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
# $Id: pluginmanager.py,v 1.2 2005/07/12 16:19:40 mbaas Exp $

## \file pluginmanager.py
## Contains the PluginManager class.

# (this line separates the \file comment from the module docs)

"""This module contains the plugin manager and associated functions.

The PluginManager class contains all the functionality to load plugins
and retrieve the objects inside the plugins.

The module already initializes a global plugin manager that can be
accessed using the functions with the same name than the manager methods.
"""

import os, os.path, sys, copy, types, glob, imp, inspect, traceback
try:
    import logging
except:
    # Define a dummy class that does nothing
    class DummyLogger:
        def debug(self, *args, **kwargs): pass
        def info(self, *args, **kwargs): pass
        def warning(self, *args, **kwargs): pass
        def error(self, *args, **kwargs): pass
        def critical(self, *args, **kwargs): pass
        def exception(self, *args, **kwargs): pass
    logging = DummyLogger()

# Exceptions:
class DuplicateObject(Exception):
    """Exception."""
    pass

class PluginAlreadyLoaded(Exception):
    """Exception."""
    pass

class UnknownFileType(Exception):
    """Exception."""
    pass

class ProtocolSpecsMissing(Exception):
    """Exception."""
    pass

class MissingName(Exception):
    """Exception."""
    pass


# Status flags
STATUS_OK         = 0x00
STATUS_EXCEPTION  = 0x01
STATUS_DUPLICATE  = 0x02

# PluginDescriptor
class PluginDescriptor:
    """Status descriptor of a plugin file.

    This class is just a container for a couple of info attributes about
    a plugin. The attributes are:

    - \b filename: The absolute file name of the plugin
    - \b filedesc: A file descriptor tuple (\em suffix, \em mode, \em type)
                 as returned by the get_suffixes() function in the \c imp
                 module (\em type is one of imp.PY_SOURCE, imp.PY_COMPILED
                 or imp.C_EXTENSION)
    - \b objcount: Number of imported objects (usually this should be
                 #classes + #functions)
    - \b classcount: Number of imported classes
    - \b funccount: Number of imported functions
    - \b status: Status flags. If everything is ok the status is STATUS_OK,
                otherwise its a combination of STATUS_EXCEPTION and
                STATUS_DUPLICATE.
    - \b traceback: A string containing the traceback message if the 
                STATUS_EXCEPTION flag is set.
    - \b objdescs: A list of object descriptors of all imported objects
    - \b module: The module object of the imported plugin
    - \b modulename: The name of the plugin module

    \see PluginManager, PluginObjectDescriptor
    """
    def __init__(self, filename="<unknown>", filedesc=None,
                 status=STATUS_OK, traceback=None,
                 objcount=0, classcount=0, funccount=0, module=None):
        self.filename   = filename
        self.filedesc   = filedesc
        self.status     = status
        self.traceback  = traceback
        self.objcount   = objcount
        self.classcount = classcount
        self.funccount  = funccount
        self.module     = module
        self.objdescs   = []
        self.modulename = None


    def __repr__(self):
        res = '<PluginDescriptor file:"%s" status:%s #objs:%s #cls:%s #func:%s>'%(self.filename, self.status, self.objcount, self.classcount, self.funccount)
        return res

    def __str__(self):
        res = 'Plugin "%s":\n'%self.filename
        s = "?"
        if self.filedesc!=None:
            t = self.filedesc[2]
            if t==imp.PY_SOURCE:
                s = "Python source"
            elif t==imp.PY_COMPILED:
                s = "Compiled Python code"
            elif t==imp.C_EXTENSION:
                s = "C extension"
            else:
                s += " (%s)"%t
        res += "  File type  : %s\n"%s
        res += "  Status     : %s\n"%self.status
        res += "  #Objects   : %s\n"%self.objcount
        res += "  #Classes   : %s\n"%self.classcount
        res += "  #Functions : %s\n"%self.funccount
        res += "  Objects    : %s\n"%self.objdescs
        res += "  Module     : %s"%getattr(self.module, "__name__", None)
        if self.status==STATUS_EXCEPTION:
            res += "\n  Traceback:\n%s"%self.traceback
        return res

    def setStatusFlag(self, flags):
        self.status |= flags

# PluginObjectDescriptor
class PluginObjectDescriptor:
    """Descriptor for a plugin object.

    - \b object: The actual plugin object
    - \b name: Object name
    - \b plugindesc: Associated plugin
    - \b status: Status flags. If everything is ok the status is STATUS_OK,
                otherwise its STATUS_DUPLICATE.

    \see PluginManager, PluginDescriptor
    """
    def __init__(self, object=None, name=None,
                 plugindesc=None, status=STATUS_OK):
        self.object = object
        self.name = name
        self.plugindesc = plugindesc
        self.status = status

    def __repr__(self):
        res = "<Plugin object '%s', status=%s>"%(self.objectIdentifier(), self.status)
        return res

    def moduleName(self):
        if self.plugindesc==None:
            return None
        else:
            return self.plugindesc.modulename

    def objectIdentifier(self):
        res = self.name
        modname = self.moduleName()
        if modname!=None:
            res = modname+"."+res
        return res

    def setStatusFlag(self, flags):
        self.status |= flags


# PluginManager
class PluginManager:
    """Loads and manages plugins.

    This class imports and manages plugins. A plugin is a file that
    gets imported just like a regular Python module. Each plugin can
    be uniquely identified either by its absolute filename or its
    module name. The module name is only available once the plugin is
    loaded. Any objects (usually classes) that are specially marked
    will be made available by the plugin manager for later retrieval.

    A class gets imported if it has an attribute \c _protocols which
    is a list of supported protocols (a protocol is identified by an
    arbitrary hashable object). Each class can be uniquely identified
    by its module name and class name.
    
    Plugins are loaded via calls to importPlugin() or importPlugins().
    You can also register any class as plugin class any time via
    the register() function. When a plugin is loaded or an object
    is registered the plugin/object is represented by a descriptor
    class (PluginDescriptor resp. PluginObjectDescriptor). These
    descriptor objects contain information about the plugin/object.
    The object descriptor also holds a reference to the actual plugin
    object.

    You can iterate over all currently available protocols via the
    iterProtocols() method. The iterProtoObjects() method can be used
    to iterate over all objects that support a particular protocol
    and iterObjects() can be used to iterate over all objects.
    The default iterator iterates over all plugin descriptors.

    Example plugin class:
    
    \code
    # file: plugin.py
    
    class PluginClass:

        _protocols = ["MyProtocol"]

        def __init__(self):
            ...
        ...
    \endcode

    And here's how to use it:

    \code
    pm = PluginManager()
    
    # Load the plugin
    pdesc = pm.importPlugin("plugin.py")
    if pdesc.status!=STATUS_OK:
      # there was an error...

    # Search a plugin class
    objdesc = pm.findObject("plugin.PluginClass")
    # Create an instance of the plugin class...
    PluginClass = objdesc.object
    instance = PluginClass()
    
    \endcode

    \see PluginDescriptor, PluginObjectDescriptor
    """

    def __init__(self):
        """Constructor."""
        
        # A list with plugin descriptors of all imported plugins
        self._plugins = []

        # A dictionary containing all objects.
        # Key is the object identifier string, value is the descriptor
        self._objects = {}

        # A dictionary with all available protocols. The key is the protocol
        # identifier (an arbitrary hashable object), the value is a list
        # of object descriptors.
        self._protocols = {}
        


    ######################################################################
    # Neue Schnittstelle:

    def __iter__(self):
        """Return an iterator that iterates over all plugin descriptors."""
        return iter(self._plugins)

    def iterProtocols(self):
        """Return an iterator that iterates over all protocols."""
        return iter(self._protocols)

    def iterProtoObjects(self, proto):
        """Return an iterator that iterates over all object descriptors supporting a particular protocol."""
        if proto in self._protocols:
            return iter(self._protocols[proto])
        else:
            return iter([])

    def iterObjects(self):
        """Return an iterator that iterates over all object descriptors.
        """
        return iter(list(self._objects.values()))

    # removeAll
    def removeAll(self):
        """Remove all plugins and plugin objects.
        """
        
        while len(self._objects)>0:
            key = list(self._objects.keys())[0]
            self.remove(self._objects[key])
            
        self._plugins = []
        self._objects = {}
        self._protocols = {}

    # importPlugin
    def importPlugin(self, filename, overwrite=False):
        """Load a plugin file.

        The given file is executed, i.e. it is imported like a module.
        Any object (usually class or function) that has an attribute
        "_protocols" is stored in the plugin manager. This attribute
        is a list of protocol specifiers (which are just arbitrary
        hashable objects).

        The file must have a suffix that's recognized as a module by
        the Python interpreter (one of the suffixes returned by
        imp.get_suffixes()), otherwise an UnknownFileType exception is
        thrown.

        It is possible to load a plugin several times if \a overwrite
        is set to True. The new definitions override the previous ones.
        However, if references or instances to old objects are kept somewhere
        they still refer to the old definition. When writing a plugin you
        should always bear in mind that the file could be executed several
        times and write your initialization code accordingly.

        The function returns a PluginDescriptor object which contains
        information about the imported plugin. 
        
        \param filename (\c str) Name of the plugin file (including path)
        \param overwrite (\c bool) Flag that indicates if objects are allowed
                  to overwrite existing plugin objects
        \return Plugin descriptor object (\c PluginDescriptor)
        \todo Lokalen Modul-Import innerhalb des Plugins testen
        """
        
        # Make the filename absolute (and normalized)
        filename = os.path.abspath(filename)
        # path: Absolute path to the file
        path, name = os.path.split(filename)
        # ext: Extension of the file
        name, ext  = os.path.splitext(filename)
        # modname: Module name (only the file name, without path and extension)
        modname = os.path.basename(name)

        # Check if the file is already imported
        oldpdesc = self.findPlugin(filename)
        if oldpdesc!=None:
            if overwrite:
                self.removePlugin(oldpdesc)
            else:
                raise PluginAlreadyLoaded("Plugin '%s' is already loaded."%filename)

        # Find the file description tuple  (suffix, mode, type)
        filedesc = [x for x in imp.get_suffixes() if x[0]==ext]
        if filedesc==[]:
            raise UnknownFileType('File "%s" is of unknown type."'%filename)
        filedesc = filedesc[0]

        # Create a plugin descriptor class where all sorts of
        # info is gathered
        pdesc = PluginDescriptor(filename=filename, filedesc=filedesc)

        # Open the file
        f = file(filename)

        # Add the plugin to the list of imported plugins
        self._plugins.append(pdesc)

        # Import the file as a module
        # The path to the file is added to the search path for modules,
        # so the file can import local modules
        try:
            oldpath = copy.copy(sys.path)
            sys.path.append(path)
            try:
                mod = imp.load_module(modname, f, filename, filedesc)
                pdesc.module = mod
                pdesc.modulename = mod.__name__
            finally:
                f.close()
                sys.path = oldpath
        except:
            # Extract the traceback message...
            lst = traceback.extract_tb(sys.exc_info()[2])
            flst = traceback.format_list(lst[1:])
            pdesc.traceback = "".join(flst)
            # Add the actual exception to the message...
            t,v = sys.exc_info()[:2]
            elst = traceback.format_exception_only(t,v)
            pdesc.traceback += "".join(elst)
            
            pdesc.setStatusFlag(STATUS_EXCEPTION)
            return pdesc

        # Examine the module namespace and register all objects that
        # have a "_protocols" attribute...
        for objname in mod.__dict__:
            obj = mod.__dict__[objname]
            if hasattr(obj, "_protocols"):
                self.register(obj, pdesc=pdesc)
                
        return pdesc

    # importPlugins
    def importPlugins(self, plugins, out = sys.stderr):
        """Import several plugins at once.

        \a plugins can be a single file/directory name or a sequence
        of file/directory names. Directories are recursively descended.

        \return A list of plugin descriptors.
        """

        if isinstance(plugins, str):
            plugins = [plugins]
        
        res = []
        for plugin in plugins:
            # Check if the file/dir is there...
            if not os.path.exists(plugin):
                out.write('ERROR: Plugin file or directory "%s" does not exist.\n'%plugin)
                continue

            # Get a list of files (either the single file or the contents
            # of the directory)
            if os.path.isfile(plugin):
                files = [plugin]
            else:
                files = glob.glob(os.path.join(plugin,"*"))

            # Import each file...
            for f in files:
                # Directory? then import recursively
                if os.path.isdir(f):
                    res += self.importPlugins(f, out)
                # File...
                else:
                    name, ext = os.path.splitext(f)
                    ext = ext.lower()
                    # Only try to import if the file is actually a
                    # Python module (other than .pyc)...
                    filedesc = [x for x in imp.get_suffixes() if x[0]==ext]
                    if filedesc!=[] and ext!=".pyc":
                        out.write('Loading plugin "'+f+'"...\n')
                        d = self.importPlugin(f)
                        res.append(d)
        return res

    # register
    def register(self, obj, name=None, pdesc=None, overwrite=False):
        """Register a plugin object.

        If \a name is None then the name of the object will be used.      

        \param obj  The plugin object
        \param name (\c str) The name of the object or None
        \param pdesc (\c PluginDescriptor) The corresponding plugin or None
        \param overwrite (\c bool) If True the object will overwrite any existing objects.
        \return Object descriptor (\c PluginObjectDescriptor)
        """

        # Fill in name automatically?
        if name==None:
            name = getattr(obj, "__name__", None)
            if name==None:
                raise MissingName("Cannot determine the name of the object.")
        if not hasattr(obj, "_protocols"):
            raise ProtocolSpecsMissing("Attribut '_protocols' in object '%s' is missing."%name)

        # Create an object descriptor
        desc = PluginObjectDescriptor(object=obj, name=name, plugindesc=pdesc)

        # Check if the object already exists...
        id = desc.objectIdentifier()
        if id in self._objects:
            if overwrite:
                olddesc = self._objects[id]
                self._removeObjectDesc(olddesc)
            else:                
                raise DuplicateObject("Plugin object '%s' already exists."%id)

        # Add the descriptor to the protocol lists...
        self._insertObjectDesc(desc)

        return desc

    # remove
    def remove(self, objdesc):
        """Remove a plugin object.

        \param objdesc (\c PluginObjectDescriptor) Object descriptor
        """
        # Check if the object is really managed by the plugin manager
        id = objdesc.objectIdentifier()
        if objdesc.objectIdentifier() not in self._objects:
            raise ValueError("The object '%s' does not exist."%id)

        self._removeObjectDesc(objdesc)

    # removePlugin
    def removePlugin(self, pdesc):
        """Remove a plugin and its imported objects from the manager.

        Note: If some classes are still used somewhere those classes
        remain active. It's just that you won't be able to retrieve
        those classes via the plugin manager.

        \param pdesc (\c PluginDescriptor) Descriptor of the plugin that should be removed
        """
        # Remove all imported objects...
        while len(pdesc.objdescs)>0:
            self.remove(pdesc.objdescs[0])

        # Remove the plugin
        self._plugins.remove(pdesc)


    # findPlugin
    def findPlugin(self, filename):
        """Return the plugin descriptor for a given filename.

        \param filename (\c str) Plugin filename
        \return Plugin descriptor (\c PluginDescriptor) or None
        """
        filename = os.path.abspath(filename)
        pl = [d for d in self._plugins if d.filename==filename]
        if pl==[]:
            return None
        else:
            return pl[0]

    # findObject
    def findObject(self, name, modname=None):
        """Return the object descriptor for the specified object.

        \param name (\c str) Object name
        \param modname (\c str) Module name
        \return Object descriptor or None
        """
        if modname!=None:
            name = modname+"."+name
        return self._objects.get(name, None)

        
    ######################################################################
    ## protected:

    def _insertObjectDesc(self, desc):
        """Insert a plugin object.

        \pre The object isn't already stored in the manager
        \pre The field desc.name contains a string
        \pre The object has a valid _protocols attribute (a sequence)
        \param desc (\c PluginObjectDescriptor) Object descriptor
        """

        # Add the object descriptor to the manager
        self._objects[desc.objectIdentifier()] = desc

        # Insert the object descriptor into the corresponding protocol dicts
        for prot in desc.object._protocols:
            if prot not in self._protocols:
                self._protocols[prot] = []
            self._protocols[prot].append(desc)

        # Update the plugin descriptor
        pdesc = desc.plugindesc
        if pdesc!=None:
            pdesc.objdescs.append(desc)
            pdesc.objcount += 1
            if inspect.isclass(desc.object):
                pdesc.classcount += 1
            if inspect.isfunction(desc.object):
                pdesc.funccount += 1


    def _removeObjectDesc(self, desc):
        """Remove a plugin object from the manager.

        \pre The object is stored in the manager
        \param desc (\c PluginObjectDescriptor) Object descriptor
        """

        # Update the plugin descriptor
        pdesc = desc.plugindesc
        if pdesc!=None:
            pdesc.objdescs.remove(desc)
            pdesc.objcount -= 1
            if inspect.isclass(desc.object):
                pdesc.classcount -= 1
            if inspect.isfunction(desc.object):
                pdesc.funccount -= 1
        
        # Remove the object from the protocol dicts
        for prot in desc.object._protocols:
            self._protocols[prot].remove(desc)
                
        # Remove the object descriptor
        del self._objects[desc.objectIdentifier()]
        

######################################################################

# Global default plugin manager
_plugin_manager = PluginManager()

def iterPlugins():
    """Iterate over all global plugins.

    This is equivalent to calling \c iter(pm) (with \c pm being the global
    plugin manager).

    \see PluginManager::__iter__()
    """
    global _plugin_manager
    return iter(_plugin_manager)

def iterProtocols():
    """Global %iterPrototols() function.

    \see PluginManager::iterProtocols()
    """
    global _plugin_manager
    return _plugin_manager.iterProtocols()

def iterProtoObjects(proto):
    """Global %iterProtoObjects() function.

    \see PluginManager::iterProtoObjects()
    """
    global _plugin_manager
    return _plugin_manager.iterProtoObjects(proto)

def iterObjects(proto):
    """Global %iterObjects() function.

    \see PluginManager::iterObjects()
    """
    global _plugin_manager
    return _plugin_manager.iterObjects()

def removeAll(obj):
    """Global %removeAll() function.

    \see PluginManager::removeAll()
    """
    global _plugin_manager
    _plugin_manager.removeAll()

def importPlugin(filename, overwrite=False):
    """Global %importPlugin() function.

    \see PluginManager::importPlugin()
    """
    global _plugin_manager
    return _plugin_manager.importPlugin(filename, overwrite)

def importPlugins(dir, out=sys.stdout):
    """Global %importPlugins() function.

    \see PluginManager::importPlugins()
    """
    global _plugin_manager
    return _plugin_manager.importPlugins(dir, out)

def register(obj, name=None, pdesc=None, overwrite=False):
    """Global %register() function.

    \see PluginManager::register()
    """
    global _plugin_manager
    return _plugin_manager.register(obj, name, pdesc, overwrite)

def remove(objdesc):
    """Global %remove() function.

    \see PluginManager::remove()
    """
    global _plugin_manager
    _plugin_manager.remove(objdesc)

def removePlugin(pdesc):
    """Global %removePlugin() function.

    \see PluginManager::removePlugin()
    """
    global _plugin_manager
    _plugin_manager.removePlugin(pdesc)

def findPlugin(filename):
    """Global %findPlugin() function.

    \see PluginManager::findPlugin()
    """
    global _plugin_manager
    return _plugin_manager.findPlugin(filename)

def findObject(name, modname=None):
    """Global %findObject() function.

    \see PluginManager::findObject()
    """
    global _plugin_manager
    return _plugin_manager.findObject(name, modname)



######################################################################

if __name__=="__main__":

    class TestClass:
        _protocols = ["myproto"]

    pm = PluginManager()

#    desc = pm.register(TestClass)                      
#    print desc.objectIdentifier()

    pm.importPlugins(".")
#    pm.importPlugin("testplugin.py")
#    pm.importPlugin("testplugin.py", overwrite=True)

    for pdesc in pm:
        print (pdesc)

    for prot in pm.iterProtocols():
        print(("Protokoll: %s"%prot))
        for odesc in pm.iterProtoObjects(prot):
            print (odesc)

    print ("Objects:")
    for desc in pm.iterObjects():
        print (desc)

#    print "Plugins:",pm._plugins
#    print "Objects:",pm._objects
#    print "Protocols:",pm._protocols

