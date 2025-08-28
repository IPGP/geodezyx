<img src="./geodezyx_toolbox_logo.png" width="300">

# geodezyx (aka ___The GeodeZYX Toolbox___) 

**Version: 5.0.0**  
**Date: 2025-08-28**  

**Authors:** Pierre Sakic (IPGP, Paris, France), Gustavo Mansur (GFZ, Potsdam, Germany), Samuel Nahmani (IPGP/IGN, Paris, France), 

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

## Toolbox's highlights
* read RINEX2 and RINEX3/4 as Pandas' DataFrame
* read SINEX, SP3 & Clock RINEX as Pandas' DataFrame
* convert most time representations and time scales used in geodesy to/from Python's `datetime`
* convert coordinates in different frames (geographic, geocentric, topocentric)
* perform easily Helmert's Transformation
* And many more!

## Installation
See the following link:  
[https://ipgp.github.io/geodezyx/getting_started.html#installation](https://ipgp.github.io/geodezyx/getting_started.html#installation)

## Documentation
See the following link:  
[https://ipgp.github.io/geodezyx/](https://ipgp.github.io/geodezyx/)

## Credits
`geodezyx` implements several methods, algorithms... initially developed by others authors.
We thank them for their work and their indirect contribution to this toolbox. Namely:
  * Goudarzi, M.A., Cocard, M. & Santerre, R. EPC: Matlab software to estimate Euler pole parameters. 
GPS Solut 18, 153–162 (2014). https://doi.org/10.1007/s10291-013-0354-4
  * Caroline Geisert, ENSTA Bretagne, for the plate motion model functions, based on Goudarzi et al.'s work.
  * Yann-Treden Tranchant, LIENSs La Rochelle, for the ocean circulation model and OBP processing functions.
  * Jean-Mathieu Nocquet and its PYACS toolbox for the plate motion model functions. https://github.com/JMNocquet/pyacs36
  * Médéric Gravelle, CNRS/LIENSs La Rochelle, who initiate the SPOTGINS import/export functions, now grouped in `spytgins` module.

## Licence
GNU Lesser General Public License, Version 3, 29 June 2007

_geodezyx_ is free software: you can redistribute it and/or modify it
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

