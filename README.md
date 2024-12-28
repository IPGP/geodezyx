<img src="./geodezyx_toolbox_logo.png" width="300">

# geodezyx (aka ___The GeodeZYX Toolbox___) 

**Version: 4.5.2**  
**Date: 2024-12-28**  
**README Revision:** 2024-11-28  

**Authors:** Pierre Sakic (IPGP, Paris, France) & Gustavo Mansur (GFZ, Potsdam, Germany) 

**Documentation:** [https://ipgp.github.io/geodezyx](https://ipgp.github.io/geodezyx/)

**GitHub repository:** [https://github.com/IPGP/geodezyx](https://github.com/IPGP/geodezyx) 

**PyPi project:** [https://pypi.org/project/geodezyx](https://pypi.org/project/geodezyx/)

**Contact e-mail:** sakic@ipgp.fr

**Citation:** Sakic, Pierre; Mansur, Gustavo; Chaiyaporn, Kitpracha; Ballu, Valérie (2019):
*The geodeZYX toolbox: a versatile Python 3 toolbox for geodetic-oriented purposes*. 
V. 4.0. GFZ Data Services. [http://doi.org/10.5880/GFZ.1.1.2019.002](http://doi.org/10.5880/GFZ.1.1.2019.002)

**Licence:** GNU LGPL v3 (see below) 

**Contributors:**
* Chaiyaporn "Na" Kitpracha (GFZ, Potsdam, Germany)
* Valérie Ballu (CNRS/La Rochelle University, France)

## Introduction

The purpose of _geodezyx_ (pronounced *geode-**zeecks***), also known as _the GeodeZYX toolbox_, is to provide all the functions which
can be useful for Geodesy and Geophysics. 

It includes low-level functions, file management functions,
time and space-coordinates conversion functions, 
data (especially GNSS observations and orbits) retrieve functions, 
plots and visual selection functions ...

It is designed for Python 3 on a LINUX Ubuntu-like system.
Also tested with Anaconda

## Documentation

See the following link:
[https://ipgp.github.io/geodezyx/](https://ipgp.github.io/geodezyx/)

## Installation 

See the following link:
[https://ipgp.github.io/geodezyx/getting_started.html#installation](https://ipgp.github.io/geodezyx/getting_started.html#installation)

## Installation (detailed legacy instructions)

### with PyPI and pip

We recommend using `pip` to do a proper installation.  
To get the latest working version, you can install directly the GitHub-hosted version (recommended):  

``pip install git+https://github.com/IPGP/geodezyx``

You can also install the version hosted on PyPI (but you will not get the latest changes)

``pip install geodezyx``

### clone and manually install from GitHub

You can manually fork and/or clone the GitHub repository
([https://github.com/IPGP/geodezyx/](https://github.com/IPGP/geodezyx/))
using your favorite flavor:
* SSH: ``git clone git@github.com:IPGP/geodezyx.git``
* https: ``git clone https://github.com/IPGP/geodezyx.git``

and install the Toolbox you downloaded with ``python setup.py install``

Alternatively, you can also add the ``geodezyx`` folder in your ``PYTHONPATH`` (for experimented users)

## Changelog

### v4.5.2, 2024-12-28
  * First version using CI/CD with GitHub Actions
  * Correct bug for mono-GNSS RINEX in `read_rnx2_obs` and `read_rnx3_obs`.
  * Correct python 3.12's annoying regex's raw string SyntaxWarning.

### v4.5.1, 2024-11-28
  * Make `geodezyx` compatible with Python 3.12
  * Improvement of PRIDE-related execution functions
  * New leap second functions, to cope with new Ubuntu 24.04 LTS
  * `seawater` and `gsw` modules are installed in `full` mode & imported on demand only
  * Correct angle conversion functions
  * New version of `conv.numpydt2dt` function
  * Misc routine improvements
  * v4.5.0 is cancelled, because `requirements.txt` is missing

### v4.4.3, 2024-05-20
  * Speed execution optimization for orbit/clock-related functions
  * Misc routine improvements

### v4.4.2, 2024-04-17
  * The GitHub repository, and the project in general, has been renamed as ``geodezyx`` (in lower case)
    to uniformize its multiple spellings and then avoid confusion.
    * It must be transparent for your clones but updating them is recommended \
      https://docs.github.com/en/repositories/creating-and-managing-repositories/renaming-a-repository
  * misc routine improvements
  * Refactoring of the GNSS data/products dowmload functions

### v4.4.1, 2024-02-08
  * The GitHub repository has now been moved under the IPGP organization.
    * It must be transparent for your clones but updating them is recommended \
      https://docs.github.com/en/repositories/creating-and-managing-repositories/renaming-a-repository  
  * The version numbering goes without a starting zero from now on.  
    Initially, this first zero was kept as a "perpetual beta" marker,  
    but The GeodeZYX toolbox is a grown-up project now!
  * Bugs corrected for `read_rnx2_obs` and `OrbDF` manipulation functions.

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

## Credits
`geodezyx` implements several methods, algorithms... initially developed by others authors.
We thank them for their work and their indirect contribution to this toolbox. Namely:
  * Goudarzi, M.A., Cocard, M. & Santerre, R. EPC: Matlab software to estimate Euler pole parameters. 
GPS Solut 18, 153–162 (2014). https://doi.org/10.1007/s10291-013-0354-4
  * Caroline Geisert, ENSTA Bretagne, for the plate motion model functions, based on Goudarzi et al.'s work.
  * Yann-Treden Tranchant, LIENSs La Rochelle, for the ocean circulation model and OBP processing functions.
  * Jean-Mathieu Nocquet and its PYACS toolbox for the plate motion model functions. https://github.com/JMNocquet/pyacs36

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

