## Changelog

### v5.0.0, 2025-07-10
**Breaking changes**:
  * (Much) faster import. The main `geodezyx` doesn't import all modules per default anymore. Thus, favor the following import methods:
    * `from geodezyx import <module>`
    * `import geodezyx.<module>`
  * We recommend `gzyx_<module>` as import alias, e.g.:
    * `from geodezyx import files_rw as gzyx_files_rw`
    * `import geodezyx.files_rw as gzyx_files_rw`
  * Two exceptions: widely-used `conv` and `utils` submodules remain loaded by the main `geodezyx` module.   
    Thus, `import geodezyx` will only load `conv` and `utils` modules.
  * NB: `from geodezyx import *` still load __all modules__. But its use is __highly discouraged__ (very slow).
  * PEP8 compliance: `conv` module functions are now PEP8 compliant, e.g. upper case function names are now lower case.
    * Exemples:
      * `dt2MJD` > `dt2_mjd`, `MJD2dt` > `mjd2_dt`
      * `XYZ2ENU` > `xyz2enu`, `GEO2XYZ` > `geo2xyz`
    * Deprecated functions are still available with a warning, but will be removed in the future.
  
  * Lot of minor imporvements and bug corrections

### v4.6.0, 2025-02-04
  * Faster RINEX3 reading (~ x5)
  * Frontend `read_rinex_obs` fonction
  * For `find_recursive` function: no more ambiguous  `case_sensitive` option, replaced with `regex`

### v4.5.3, 2024-12-28
  * First version using CI/CD with GitHub Actions
  * Correct bug for mono-GNSS RINEX in `read_rnx2_obs` and `read_rnx3_obs`.
  * Correct python 3.12's annoying regex's raw string SyntaxWarning.
  * v4.5.2 is cancelled, because `requirements.txt` is too version-restrictive. 

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