<img src="./geodezyx_toolbox_logo.png" width="300">

# GeodeZYX Toolbox

**Version 4.4.1 / 2023-02-08**, README Revision: 2024-02-08


**Authors:** Pierre Sakic (IPGP, Paris, France) & Gustavo Mansur (GFZ, Potsdam, Germany)

**Documentation:** [https://ipgp.github.io/geodezyx-toolbox](https://ipgp.github.io/geodezyx-toolbox/)

**GitHub repository:** [https://github.com/IPGP/geodezyx-toolbox](https://github.com/IPGP/geodezyx-toolbox) 

**PyPi project:** [https://pypi.org/project/geodezyx](https://pypi.org/project/geodezyx/)

**Contact e-mail:** sakic@ipgp.fr

**Citation:** Sakic, Pierre; Mansur, Gustavo; Chaiyaporn, Kitpracha; Ballu, Valérie (2019):
*The geodeZYX toolbox: a versatile Python 3 toolbox for geodetic-oriented purposes*. 
V. 4.0. GFZ Data Services. [http://doi.org/10.5880/GFZ.1.1.2019.002](http://doi.org/10.5880/GFZ.1.1.2019.002)

**Licence:** GNU LGPL v3 (see below) 

**Contributors:**
* Kitpracha "Na" Chaiyaporn (GFZ, Potsdam, Germany)
* Valérie Ballu (CNRS/La Rochelle University, France)


## Introduction

The purpose of the GeodeZYX toolbox (pronounced *geode-**zeecks***) is to provide all the functions which
can be useful for Geodesy and Geophysics. 

It includes low-level functions, file management functions,
time and space-coordinates conversion functions, 
data (especially GNSS observations and orbits) retrieve functions, 
plots and visual selection functions ...

It is designed for Python 3 on a LINUX Ubuntu-like system.
Also tested with Anaconda

## Documentation

See the following link:
[https://geodezyx.github.io/geodezyx-toolbox/](https://geodezyx.github.io/geodezyx-toolbox/)

## Installation 

See the following link:
[https://geodezyx.github.io/geodezyx-toolbox/getting_started.html#installation](https://geodezyx.github.io/geodezyx-toolbox/getting_started.html#installation)

## Installation (detailed legacy instructions)

### with PyPI and pip

We recommend using `pip` to do a proper installation.  
To get the latest working version, you can install directly the GitHub-hosted version (recommended):  

``pip install git+https://github.com/IPGP/geodezyx-toolbox``

You can also install the version hosted on PyPI (but you will not get the latest changes)

``pip install geodezyx``

### clone and manually install from GitHub

You can manually fork and clone the GitHub repository
[https://github.com/GeodeZYX/geodezyx-toolbox/](https://github.com/GeodeZYX/geodezyx-toolbox/)

and install the Toolbox you downloaded with ``python setup.py install``

Alternatively, you can also add the ``geodezyx`` folder in your ``PYTHONPATH`` (for experimented users)

## Changelog

### v4.4.1, 2023-02-08
  * The GitHub repository has now been moved under the IPGP organization.
    * It must be transparent for your clones but updating them is recommended \
      https://docs.github.com/en/repositories/creating-and-managing-repositories/renaming-a-repository  
  * The numbering goes without a starting zero anymore. The GeodeZYX toolbox is a grown-up project now!
  * Bugs corrected for `read_rnx2_obs` and `OrbDF` manipulation function.

### v0.4.4.0, 2023-11-24
  * The toolbox turns to the _GNU Lesser General Public License version 3_
  * Module docstring has been updated
  * Angle conversion functions have been refactored
  * sp3/clk DataFrame column names are renamed to better fit the data content:
    * `sat` > `prn`, `const` > `sys`, `sv` > `prni`, `AC` > `ac`
      
### v0.4.3.6, 2023-06-30
  * a routine version update (bug corrections...)
  * colors in the logger
  * NB: v0.4.3.5, v0.4.3.4 & v0.4.3.3 are cancelled releases

### v0.4.3.2, 2023-03-08
  * 1000th commit on GitHub 🥳
  * GitHub repository is renamed for simplification 
    * `GeodeZYX-Toolbox_v4` becomes `geodezyx-toolbox` (lowercase and without version)
    * It must be transparent for your clones but updating them is recommended \
      https://docs.github.com/en/repositories/creating-and-managing-repositories/renaming-a-repository  
  * beta for GROOPS automatized run functions
  * Dropbox download function
  * UTM coordinates conversion function (IGN algorithm) 
  * implementation of UTM coordinates in the TimeSeries class
  * logger improvement (shorter timestamps & message ranks)

### v0.4.3.1, 2022-12-09
  * Routine update: multiple new features/functions and bug corrections

### v0.4.3.0, 2022-02-10
  * enhanced logger replaces basic prints.
  * GeodeZYX is now "virtualenv-ready", i.e. stand-alone based on the ``setup.py`` required modules.
  * a ``full`` version is set, for advanced installation.

## Licence

GNU Lesser General Public License, Version 3, 29 June 2007

Copyright © 2019 Helmholtz Centre Potsdam GFZ 
German Research Centre for Geosciences, Potsdam, Germany 
(Pierre Sakic, Gustavo Mansur, and Kitpracha "Na" Chaiyaporn, Valérie Ballu)

The geodeZYX toolbox is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published by the 
Free Software Foundation, either version 3 of the License, or 
(at your option) any later version. The geodeZYX toolbox is distributed 
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
See the GNU Lesser General Public License for more details. 

You should have received a copy of the GNU Lesser General Public License 
along with this program. If not, see http://www.gnu.org/licenses/.

```
       _____                _      ________     ____   __  _______          _ _
      / ____|              | |    |___  /\ \   / /\ \ / / |__   __|        | | |
     | |  __  ___  ___   __| | ___   / /  \ \_/ /  \ V /     | | ___   ___ | | |__   _____  __
     | | |_ |/ _ \/ _ \ / _` |/ _ \ / /    \   /    > <      | |/ _ \ / _ \| | '_ \ / _ \ \/ /
     | |__| |  __/ (_) | (_| |  __// /__    | |    / . \     | | (_) | (_) | | |_) | (_) >  <
      \_____|\___|\___/ \__,_|\___/_____|   |_|   /_/ \_\    |_|\___/ \___/|_|_.__/ \___/_/\_\
```

_I'd rather be burned alive than to program in Perl again._  
 G. Mansur, 28th IUGG General Assembly, Berlin, Germany, 2023

