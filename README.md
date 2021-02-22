       _____                _      ________     ____   __  _______          _ _
      / ____|              | |    |___  /\ \   / /\ \ / / |__   __|        | | |
     | |  __  ___  ___   __| | ___   / /  \ \_/ /  \ V /     | | ___   ___ | | |__   _____  __
     | | |_ |/ _ \/ _ \ / _` |/ _ \ / /    \   /    > <      | |/ _ \ / _ \| | '_ \ / _ \ \/ /
     | |__| |  __/ (_) | (_| |  __// /__    | |    / . \     | | (_) | (_) | | |_) | (_) >  <
      \_____|\___|\___/ \__,_|\___/_____|   |_|   /_/ \_\    |_|\___/ \___/|_|_.__/ \___/_/\_\


# GeodeZYX Toolbox

**Version 0.4.1.0 / 2021-02-22**

README Revision: 2021-02-22


**Authors: **Pierre Sakic, Gustavo Mansur, and Kitpracha "Na" Chaiyaporn
GFZ, Potsdam, Germany, with contributions from Valérie Ballu (CNRS/La Rochelle University)

**Contact e-mail:** pierre.sakic@gfz-potsdam.de

**Citation: ***Sakic, Pierre; Mansur, Gustavo; Chaiyaporn, Kitpracha; Ballu, Valérie (2019):
The geodeZYX toolbox: a versatile Python 3 toolbox for geodetic-oriented purposes. 
V. 4.0. GFZ Data Services. http://doi.org/10.5880/GFZ.1.1.2019.002*

## LICENCE

GNU General Public License, Version 3, 29 June 2007

Copyright © 2019 Helmholtz Centre Potsdam GFZ 
German Research Centre for Geosciences, Potsdam, Germany 
(Pierre Sakic, Gustavo Mansur, and Kitpracha "Na" Chaiyaporn, Valérie Ballu)

The geodeZYX toolbox is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the 
Free Software Foundation, either version 3 of the License, or 
(at your option) any later version. The geodeZYX toolbox is distributed 
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
See the GNU General Public License for more details. 

You should have received a copy of the GNU General Public License 
along with this program. If not, see http://www.gnu.org/licenses/.

## INTRODUCTION

The purpose of this GeodeZYX toolbox is to provide all the functions which
can be useful for Geodesy and Geophysics. 

It includes low level functions, file management functions,
time and space-coordinates conversion functions, 
data (especially GNSS observations and orbits) retrieve functions, 
plots and visual selection functions ...

It is designed for Python 3 on a LINUX Ubuntu-like system.
Also tested with Anaconda

## COMPONENTS

### Main modules

  * `conv`           : time and coordinates conversion
  * `utils`          : Shell-like functions (grep, find in folder ...)
                 and functions to optimize list management
  * `time_series`    : module to handle Geodetic time-series
  * `stats`          : Low-level statistics functions & outlier detection functions
  * `reffram`        : Reference Frame & high-level coordinates conversion
  * `operational`    : Download GNSS data/products from various servers
                 Read sitelogs 
  * `files_rw`       : read, import & write geodetic files
  * `geodyn`         : Euler pole calculator and plot velocity field maps

### Other modules
TBC


## INSTALLATION





