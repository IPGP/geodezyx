.. _getting_started:

Getting started
===============

Installation
------------

See the ``README.md`` file


Minimal exemple
---------------

To test if the `GeodeZYX Toolbox` is well installed, import:
::

    #### Import star style
    from geodezyx import *                   # Import the GeodeZYX modules
    from geodezyx.externlib import *         # Import the external modules
    from geodezyx.megalib.megalib import *   # Import the legacy modules names

If the modules are well imported (without errors in the console), fine! you have the `GeodeZYX Toolbox`!
If not, check again the installation README

Nevertheless, star imports are not recomended. It is better to import the modules separatelly, e.g.: 

::

    import geodezyx.conv as conv                  # Import the conversion module

See the description of each module in the following section

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

