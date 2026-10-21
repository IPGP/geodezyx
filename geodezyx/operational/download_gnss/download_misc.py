#### Import the logger
import logging
import os
import shutil

########## BEGIN IMPORT ##########
#### External modules
from geodezyx import conv

log = logging.getLogger("geodezyx")


##########  END IMPORT  ##########


def orbclk_long2short_name(
    longname_filepath_in,
    rm_longname_file=False,
    center_id_last_letter=None,
    center_manual_short_name=None,
    force=False,
    dryrun=False,
    output_dirname=None,
):
    """
    Rename a long naming new convention IGS product file to the short old
    convention.

    Naming will be done automatically based on the 3 first characters of the
    long AC id. e.g. CODE => cod, GRGS => grg, NOAA => noa ...

    Parameters
    ----------
    longname_filepath_in : str
        Full path of the long name product file.
    rm_longname_file : bool
        Remove the original long name product file.
    center_id_last_letter : str
        Replace the last letter of the short AC id by another letter
        (see note below).
    center_manual_short_name : str
        Replace completely the long name with this one.
        Overrides center_id_last_letter.
    force : bool
        If False, skip if the file already exists.
    dryrun : bool
        If True, don't rename effectively, just output the new name.
    output_dirname : str
        Directory where the output shortname will be created.
        If None, will be created in the same folder as the input longname.

    Returns
    -------
    shortname_filepath : str
        Path of the short old-named product file.

    Note
    ----
    If you rename MGEX orbits, we advise to set
    center_id_last_letter="m".
    The AC code name will be changed to keep a MGEX convention
    (but any other character can be used too).

    e.g. for Bern's products, the long id is CODE

    if center_id_last_letter=None, it will become cod,
    if center_id_last_letter=m, it will become com
    """

    log.info("will rename " + longname_filepath_in)

    longname_basename = os.path.basename(longname_filepath_in)
    longname_dirname = os.path.dirname(longname_filepath_in)

    if not output_dirname:
        output_dirname = longname_dirname

    center = longname_basename[:3]

    if center_manual_short_name:
        center = center_manual_short_name
    elif center_id_last_letter:
        center_as_list = list(center)
        center_as_list[-1] = center_id_last_letter
        center = "".join(center_as_list)

    yyyy = int(longname_basename.split("_")[1][:4])
    doy = int(longname_basename.split("_")[1][4:7])

    day_dt = conv.doy2dt(yyyy, doy)

    wwww, dow = conv.dt2gpstime(day_dt)

    shortname_prefix = center.lower() + str(wwww).zfill(4) + str(dow)

    ### Type handling
    shortname = shortname_prefix  # default
    if "SP3" in longname_basename:
        shortname = shortname_prefix + ".sp3"
    elif "CLK" in longname_basename:
        shortname = shortname_prefix + ".clk"
    elif "ERP" in longname_basename:
        shortname = shortname_prefix + ".erp"
    elif "BIA" in longname_basename:
        shortname = shortname_prefix + ".bia"
    elif "SNX" in longname_basename:
        shortname = shortname_prefix + ".snx"
    else:
        log.error("filetype not found for " + longname_basename)

    ### Compression handling
    if longname_basename[-3:] == ".gz":
        shortname = shortname + ".gz"
    elif longname_basename[-2:] == ".Z":
        shortname = shortname + ".Z"

    shortname_filepath = os.path.join(output_dirname, shortname)

    if not force and os.path.isfile(shortname_filepath):
        log.info("skip " + longname_filepath_in)
        log.info(shortname_filepath + " already exists")
        return shortname_filepath

    if not dryrun:
        log.info("renaming " + longname_filepath_in + " => " + shortname_filepath)
        shutil.copy2(longname_filepath_in, shortname_filepath)

    if rm_longname_file and not dryrun:
        log.info("remove " + longname_filepath_in)
        os.remove(longname_filepath_in)

    return shortname_filepath
