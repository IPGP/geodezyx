.. _getting_started:

===============
Getting started
===============

.. _install: 

------------
Installation
------------

Standard installation (for non-admin users)
-------------------------------------------

We recommend to use ``pip`` to install the last stable `PyPi-hosted <https://pypi.org/project/geodezyx/>`_ version:

``pip install geodezyx``

To get the **latest developpement version** you can install directly the `GitHub-hosted <https://github.com/IPGP/geodezyx>`_ version:

``pip install git+https://github.com/IPGP/geodezyx``

If you **reinstall** a new version the package, you should uninstall the existing one and/or use the ``--upgrade`` option:

``pip uninstall geodezyx``

``pip install --upgrade geodezyx``


Full installation (for sudoers/advanced users)
----------------------------------------------

This version installs facultative external modules, which require system libraries only installable when you have the sudo/admin rights, namely ``netCDF4``, ``kepler.py``, ``ncompress``, ``seawater``, ``gsw``.

We recommend to use ``pip`` to install the last stable `PyPi-hosted <https://pypi.org/project/geodezyx/>`_ version:

``pip install geodezyx[full]``

To get the **latest developpement version** you can install directly the `GitHub-hosted <https://github.com/IPGP/geodezyx>`_ version:

``pip install "geodezyx[full] @ git+https://github.com/IPGP/geodezyx"``

Clone and manually install from GitHub (for devloppers)
-------------------------------------------------------

Manually fork and/or clone the GitHub repository (https://github.com/IPGP/geodezyx) using your favorite flavor:

* https: ``git clone https://github.com/IPGP/geodezyx.git`` (recommended)
* SSH: ``git clone git@github.com:IPGP/geodezyx.git``

Go to the folder where you cloned the repository:
* ``cd geodezyx``

And then install the toolbox you downloaded. Three solutions are possible:

1. ``pip install -e .`` (`editable mode <https://setuptools.pypa.io/en/latest/userguide/development_mode.html>`_, recommended if you want to edit the source code)
2. ``python setup.py install`` (standard mode)
3. Add the ``geodezyx`` folder in your ``PYTHONPATH``, for experimented (and old-fashoned) users


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

