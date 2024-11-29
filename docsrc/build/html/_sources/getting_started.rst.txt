.. _getting_started:

===============
Getting started
===============

.. _install: 

------------
Installation
------------

with PyPi and pip
------------------

Standard installation (for non-admin users)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We recommend to use ``pip`` to do a proper installation:  
``pip install geodezyx``

To get the latest working version you can install directly the GitHub-hosted version:  
``pip install git+https://github.com/IPGP/geodezyx``

Full installation (for sudoers/advanced users)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This version installs facultative external modules, which require system libraries only installable when you have the sudo/admin rights, namely ``netCDF4``, ``kepler.py``, ``ncompress``.

We recommend to use ``pip`` to do a propper installation:  
``pip install geodezyx[full]``

To get the latest working version you can install directly the GitHub-hosted version:  
``pip install 'git+https://github.com/IPGP/geodezyx#egg=geodezyx[full]'``


clone and manually install from GitHub
--------------------------------------

You can manually fork and/or clone the GitHub repository (https://github.com/IPGP/geodezyx/) using your favorite flavor:

* SSH: ``git clone git@github.com:IPGP/geodezyx.git``
* https: ``git clone https://github.com/IPGP/geodezyx.git``

and install the Toolbox you downloaded with ``python setup.py install``

Alternatively, you can also add the ``geodezyx`` folder in your ``PYTHONPATH`` (for experimented users)

---------------
Minimal exemple
---------------

To test if the `geodezyx` toolbox is well installed, import:
::

    #### Import
    import geodezyx           # Import the geodezyx modules

If the module are well imported (without errors in the console), fine! you have installed the `geodezyx` toolbox!
If not, check again the potential errors during the installation.

It might be better to import the modules separatelly, e.g.:

::

    import geodezyx.conv as conv                  # Import the conversion module

See the description of each module in the following section

--------------------------
Description of the Modules
--------------------------

* :py:mod:`geodezyx.athmo`: functions implementing troposhere/ionosphere models
* :py:mod:`geodezyx.conv`: time and coordinates conversion
* :py:mod:`geodezyx.files_rw`: read, import & write geodetic files
* :py:mod:`geodezyx.geodyn`: Euler pole calculator and plot velocity field maps
* :py:mod:`geodezyx.marine`: contains a GEBCO bathymetry extrapolator
* :py:mod:`geodezyx.operational`: Download GNSS data/products from various servers / Read and Import sitelogs
* :py:mod:`geodezyx.reffram`: Reference Frame & high-level coordinates conversion
* :py:mod:`geodezyx.stats`: Low-level statistics functions & outlier detection functions / Low-level least-square functions
* :py:mod:`geodezyx.time_series`: module to handle Geodetic time-series
* :py:mod:`geodezyx.utils`: Shell-like functions (grep, find in folder ...) and functions to optimize list management

