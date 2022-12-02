       _____                _      ________     ____   __  _______          _ _
      / ____|              | |    |___  /\ \   / /\ \ / / |__   __|        | | |
     | |  __  ___  ___   __| | ___   / /  \ \_/ /  \ V /     | | ___   ___ | | |__   _____  __
     | | |_ |/ _ \/ _ \ / _` |/ _ \ / /    \   /    > <      | |/ _ \ / _ \| | '_ \ / _ \ \/ /
     | |__| |  __/ (_) | (_| |  __// /__    | |    / . \     | | (_) | (_) | | |_) | (_) >  <
      \_____|\___|\___/ \__,_|\___/_____|   |_|   /_/ \_\    |_|\___/ \___/|_|_.__/ \___/_/\_\


# GeodeZYX Toolbox

**Version 0.4.1.0 / 2021-02-22**

README Revision: 2021-06-22


**Authors:** Pierre Sakic, Gustavo Mansur, and Kitpracha "Na" Chaiyaporn
(GFZ, Potsdam, Germany), with contributions from Valérie Ballu (CNRS/La Rochelle University, France)

**Contact e-mail:** pierre.sakic@gfz-potsdam.de

**Citation:*** Sakic, Pierre; Mansur, Gustavo; Chaiyaporn, Kitpracha; Ballu, Valérie (2019):
The geodeZYX toolbox: a versatile Python 3 toolbox for geodetic-oriented purposes. 
V. 4.0. GFZ Data Services. http://doi.org/10.5880/GFZ.1.1.2019.002*

## Licence

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

## Introduction

The purpose of this GeodeZYX toolbox is to provide all the functions which
can be useful for Geodesy and Geophysics. 

It includes low level functions, file management functions,
time and space-coordinates conversion functions, 
data (especially GNSS observations and orbits) retrieve functions, 
plots and visual selection functions ...

It is designed for Python 3 on a LINUX Ubuntu-like system.
Also tested with Anaconda

## Installation 

See the following link:
[https://geodezyx.github.io/GeodeZYX-Toolbox_v4/getting_started.html#installation](https://geodezyx.github.io/GeodeZYX-Toolbox_v4/getting_started.html#installation)

## Installation (legacy instruction)

### Automatic installation from PyPi

Should be implemented soon !

### Download and install the toolbox from GitHub

#### Github in a nutshell
If you are not familiar with Git and Github have a look in 
https://guides.github.com/

We decide to adopt the forking+pull request workflow
https://guides.github.com/activities/forking/

##### the proper basic steps
* Create an account on GitHub
* Go on the page of the project :
  https://github.com/GeodeZYX/GeodeZYX-Toolbox_v4
* Fork it, i.e. click on the "Fork" button
* Clone your fork on your computer ("Clone" button) 

#### Install the toolbox

### Old style installation (should be avoided)

If you do not want to use the “pythonic” installation way (with `setup.py`) yet. Here are the two alternative options to install the toolbox from GitHub.

#### Mandatory external librairies as prerequisite

We suppose that the regular Python libraires (matplotlib, numpy, scipy ...)
are already installed and up-to-date.

You have to install some annex Python libraries which are not installed per default(tested with Anaconda and Ubuntu). Namely :

tabulate
collections
natsort
pyorbital
sympy
bs4 

##### With the standard Python Interpreter of Ubuntu (most cases)
In a terminal use pip to install the libraries :

    sudo pip3 -I -U install <name of the library above>

##### With Anaconda
update Anaconda first

    conda update --all

Install the package in Anaconda
    conda install <name of the library above>


#### Installation for your whole environnement

* In a terminal, edit your `.bashrc` using e.g. `nano ~/.bashrc`
* add a line:   
`export PYTHONPATH=$PYTHONPATH:<path to the toolbox>/GeodeZYX-Toolbox_v4/geodezyx`

* The toolbox is installed !

#### Spyder user

For the ones who use the Graphical User Interface Spyder (which I recommend):

* click on the “PYTHONPATH manager” button (the Python logo-shaped icon between the spanner icon and the left arrow icon on the main toolbar)
*  add the folder
<path to the toolbox>/GeodeZYX-Toolbox_v4/geodezyx

* The toolbox is installed !

### check the installation

in a python terminal run

    from geodezyx import *

If no error appears, then the toolbox is well installed!
Else, it may miss a library. Follow the *prerequisite* section




