import collections

########## BEGIN IMPORT ##########
#### External modules
import concurrent.futures
import datetime as dt

#### Import the logger
import logging
import os
import subprocess
import numpy as np
import shutil as shutils
from threading import Lock

#### geodeZYX modules
from geodezyx import files_rw
from geodezyx import operational
from geodezyx import utils
from geodezyx import conv
from geodezyx.operational.soft_frontend.rtklib_frontend import _prods2tmp, _read_cfg, _write_cfg

log = logging.getLogger("geodezyx")
##########  END IMPORT  ##########


def rtklib_run_from_rinex(
    rnx_rover,
    rnx_base,
    generik_conf,
    out_dir,
    tmp_dir=None,
    prod_dir=None,
    experience_prefix="",
    rover_auto_conf=False,
    base_auto_conf=True,
    xyz_rover=[0, 0, 0],
    xyz_base=[0, 0, 0],
    outtype="auto",
    calc_center="IGS0OPSFIN",
    force=False,
    clean_tmp=False,
    exe_path="/home/psakicki/SOFTWARE/RTKLIB_explorer/RTKLIB/app/consapp/rnx2rtkp/gcc/rnx2rtkp",
):
    r"""
    Run RTKLIB `rnx2rtkp` from rover/base RINEX observations using a generic RTKLIB
    configuration file, optionally overriding antenna/receiver metadata from RINEX
    headers and downloading required GNSS products (SP3 precise orbits and BRDC nav).

    The function:
      - creates output and temporary directory structure
      - uncompresses compressed RINEX inputs if needed
      - reads start/end times and sampling interval from the rover/base RINEX
      - builds a run-specific RTKLIB config file from `generik_conf`
      - optionally fills rover/base antenna positions and eccentricities from RINEX headers
      - downloads required products into `prod_dir`:
        - SP3 precise orbit from `calc_center`
        - BRDC navigation RINEX
      - calls the external executable `rnx2rtkp` with assembled arguments
      - writes outputs to the `out_dir` directory

    Parameters
    ----------
    rnx_rover : str | os.PathLike
        Path to rover observation RINEX. May be compressed (e.g. `.gz`, `.Z` and
        other formats supported by `geodezyx.operational.check_if_compressed_rinex`).
    rnx_base : str | os.PathLike
        Path to base observation RINEX. May be compressed.
    generik_conf : str | os.PathLike
        Path to a generic RTKLIB configuration file (key=value format). This file is
        parsed by `read_conf_file` and then overridden according to function options.
    out_dir : str | os.PathLike
        Directory where results are saved. This parameter is mandatory.
    tmp_dir : str | os.PathLike | None, default=None
        Temporary directory for intermediate files. Optional. If not provided,
        defaults to `out_dir/TMP`.
    prod_dir : str | os.PathLike | None, default=None
        Directory where to search for orbits, clocks, and BRDC files. Optional.
        If not provided, defaults to `tmp_dir`.
    experience_prefix : str, default=""
        Prefix added to the output file stem used to name the generated `.conf`
        and `.out` files.
    rover_auto_conf : bool, default=False
        If True, reads rover RINEX header and injects rover antenna type, XYZ position,
        and antenna eccentricities into the produced config (keys `ant1-\*`).
        If False, rover settings remain those defined in `generik_conf`.
    base_auto_conf : bool, default=True
        If True, reads base RINEX header and injects base antenna type, XYZ position,
        and antenna eccentricities into the produced config (keys `ant2-\*`).
    xyz_rover : list[float], default=[0, 0, 0]
        Rover station ECEF XYZ coordinates in meters \([X, Y, Z]\). If `xyz_rover[0] != 0`,
        these coordinates override the rover RINEX header coordinates; otherwise the rover
        header coordinates are used.
    xyz_base : list[float], default=[0, 0, 0]
        Base station ECEF XYZ coordinates in meters \([X, Y, Z]\). If `xyz_base[0] != 0`,
        these coordinates override the base RINEX header coordinates; otherwise the base
        header coordinates are used.
        Note: this is a simple sentinel check; a real X coordinate of 0 would be treated
        as "not provided".
    outtype : str, default="auto"
        Output solution format. If not `"auto"`, sets RTKLIB `out-solformat` to one of:
        `"dms"`, `"deg"`, `"xyz"`, `"enu"` (case-insensitive).
        If `"auto"`, uses whatever is defined in `generik_conf`.
    calc_center : str, default="IGS0OPSFIN"
        Analysis center/product identifier used when downloading GNSS precise products
        via `geodezyx.operational.download_gnss_products`. Must match what that
        downloader expects.
    force : bool, default=False
        If True, forces reprocessing even if output file already exists.
    clean_tmp : bool, default=False
        If True, deletes all contents of `tmp_dir` at the beginning of the
        function execution.
    exe_path : str, default=...
        Filesystem path to the RTKLIB `rnx2rtkp` executable.

    Returns
    -------
    out_result_fil : str
        Path to the RTKLIB output solution file (`.out`).

    Raises
    ------
    FileNotFoundError
        If no SP3 orbit file or no BRDC navigation file can be found locally or downloaded.

    Notes
    -----
    - The generated config file is written to `out_dir` and named using:
      `<experience_prefix>\_<rover4>\_<base4>\_<YYYY\_DOY>.conf`.
    - The RTKLIB command is executed via `subprocess.call(..., shell=True)` using
      `/bin/bash`, with a single concatenated command string.
    - A warning is emitted if the base RINEX time span does not fully cover the rover span.
    """
    # Directory structure setup
    # out_dir: mandatory, where results are saved
    out_dir = utils.create_dir(out_dir)

    # tmp_dir: optional, defaults to out_dir/TMP if not provided
    if tmp_dir is None:
        tmp_dir = os.path.join(out_dir, "TMP")
    tmp_dir = utils.create_dir(tmp_dir)

    if clean_tmp:
        shutils.rmtree(tmp_dir + "/*")

    # prod_dir: optional, defaults to tmp_dir if not provided
    # This is where we search for orbits/clocks/BRDC files
    if prod_dir is None:
        prod_dir = tmp_dir
    else:
        prod_dir = utils.create_dir(prod_dir)

    # RINEX START & END - FAST
    rov_srt_fast = conv.rinexname2dt(rnx_rover)
    bas_srt_fast = conv.rinexname2dt(rnx_base)

    # RINEX NAMES
    rov_name = os.path.basename(rnx_rover)[0:4]
    bas_name = os.path.basename(rnx_base)[0:4]

    # EXPERIENCE / OUTPUT NAMES
    srt_str = rov_srt_fast.strftime("%Y_%j_%H%M")
    exp_full_name = "_".join((experience_prefix, rov_name, bas_name, srt_str))

    out_conf_fil = os.path.join(out_dir, exp_full_name + ".conf")
    out_result_fil = os.path.join(out_dir, exp_full_name + ".out")

    if (
        not force
        and os.path.isfile(out_result_fil)
        and os.path.getsize(out_result_fil) > 2000
    ):
        log.info(f"RTKLIB output file {out_result_fil} already exists. Skipping...")
        return out_result_fil

    #### START HEAVY PROCESSING ####
    log.info("Starting RTKLIB RUN for {}".format(exp_full_name))
    # uncompressing rinex if compressed
    if operational.check_if_compressed_rinex(rnx_rover):
        rnx_rover = operational.crz2rnx(rnx_rover, tmp_dir)
    if operational.check_if_compressed_rinex(rnx_base):
        rnx_base = operational.crz2rnx(rnx_base, tmp_dir)

    # RINEX START & END - PRECISE
    rov_srt, rov_end, rov_itv = operational.rinex_start_end(rnx_rover, True)
    bas_srt, bas_end, bas_itv = operational.rinex_start_end(rnx_base, True)

    log.info(
        "Rover: {} | Base: {} | Start: {} | End: {} | Interval: {} s".format(
            rov_name, bas_name, rov_srt, rov_end, rov_itv
        )
    )

    # READ GENERIC CONF FILE
    dicoconf = _read_cfg(generik_conf)

    if not outtype.lower() == "auto":
        dicoconf["out-solformat"] = outtype.lower()
        log.info(f"out-solformat was 'auto', set to: {dicoconf['out-solformat']}")

    def _edit_dicoconf(rnx_inp, xyz_inp, ant_n=1):
        antobj, recobj, siteobj, locobj = files_rw.read_rinex_2_dataobjts(rnx_inp)
        n = str(ant_n)
        dicoconf[f"ant{n}-postype"] = "xyz"
        dicoconf[f"ant{n}-anttype"] = antobj.Antenna_Type
        if xyz_inp[0] != 0:
            dicoconf[f"ant{n}-pos1"] = xyz_inp[0]
            dicoconf[f"ant{n}-pos2"] = xyz_inp[1]
            dicoconf[f"ant{n}-pos3"] = xyz_inp[2]
        else:
            dicoconf[f"ant{n}-pos1"] = locobj.X_coordinate_m
            dicoconf[f"ant{n}-pos2"] = locobj.Y_coordinate_m
            dicoconf[f"ant{n}-pos3"] = locobj.Z_coordinate_m
        dicoconf[f"ant{n}-antdelu"] = antobj.Up_Ecc
        dicoconf[f"ant{n}-antdeln"] = antobj.North_Ecc
        dicoconf[f"ant{n}-antdele"] = antobj.East_Ecc

    # Edit conf file dic
    if rover_auto_conf:
        _edit_dicoconf(rnx_rover, xyz_rover, ant_n=1)
    if base_auto_conf:
        _edit_dicoconf(rnx_base, xyz_base, ant_n=2)

    if not (bas_srt <= rov_srt <= rov_end <= bas_end):
        log.warning(
            "rover/base epoch inconsistency: not bas_srt <= rov_srt <= rov_end <= bas_end !!!"
        )

    # write conf file
    _write_cfg(dicoconf, out_conf_fil)

    ##### ORBITS
    # Function to decompress files into tmp_dir
    ## SP3
    # if "FIN" in calc_center or "RAP" in calc_center:
    #     orb_srt = conv.round_dt(bas_srt,"1D", mode="floor")
    #     orb_end = conv.round_dt(bas_end,"1D", mode="floor")
    # elif "ULT" in calc_center:
    #     orb_srt = conv.round_dt(bas_srt,"6H", mode="floor")
    #     orb_end = conv.round_dt(bas_end,"6H", mode="floor")
    # else:
    #     orb_srt = bas_srt
    #     orb_end = bas_end
    #
    # prodlis = operational.download_gnss_products(
    #     archive_dir=prod_dir,
    #     startdate=orb_srt,
    #     enddate=orb_end,
    #     archtype="year/doy",
    #     AC_names=(calc_center , ),
    #     archive_center="ign"
    # )

    prodlis = operational.dl_orbclk(
        prod_dir,
        (bas_srt, bas_end),
        calc_center,
        prod_types=("sp3", "clk"),
        data_centers=("ign",),
    )

    if not prodlis or not prodlis[0]:
        log.error("No SP3/CLK product files found remotely nor locally")
        log.error(f"Is analysis center/latency correct?: {calc_center}")
        log.warning(f"We continue without precise orbits/only with broadcast ones")
        prodlis_ok = []
    else:
        prodlis_ok = [_prods2tmp(orb, tmp_dir) for orb in prodlis]

    ### BRDC
    # statdic = dict()
    # statdic["nav"] = ["BRDC"]
    # nav_srt = conv.round_dt(bas_srt, "1D", mode="floor")
    # nav_end = conv.round_dt(bas_end, "1D", mode="floor")
    brdclis = operational.dl_brdc(prod_dir, (bas_srt, bas_end))

    # brdclis = operational.download_gnss_rinex(
    #    statdic, prod_dir, nav_srt, nav_end, archtype="year/doy"
    # )

    # brdc_path_lis , brdc_bool_lis = zip(*brdclis)
    if len(brdclis) == 0:  # or sum(brdc_bool_lis) == 0:
        log.error("No BRDC nav file found remotely nor locally")
        raise FileNotFoundError("No BRDC nav file found remotely nor locally")
    else:
        # brdclis_ok = [prods2tmp(n,tmp_dir) for n,b in brdclis if b]
        brdclis_ok = [_prods2tmp(n, tmp_dir) for n in brdclis]

    # Command
    arg_config = "-k " + out_conf_fil
    arg_interval = "-ti " + str(np.round(rov_itv, 6))
    arg_mode = ""
    # arg_mode="-p 4"
    arg_resultfile = "-o " + out_result_fil
    # com_combinsol="-c"

    bigcomand = " ".join(
        (
            exe_path,
            arg_config,
            arg_interval,
            arg_mode,
            arg_resultfile,
            rnx_rover,
            rnx_base,
            " ".join(brdclis_ok),
            " ".join(prodlis_ok),
        )
    )
    log.info(
        "Rover: {} | Base: {} | Start: {} | End: {} | Interval: {} s".format(
            rov_name, bas_name, rov_srt, rov_end, rov_itv
        )
    )
    log.info(bigcomand)
    subprocess.call([bigcomand], executable="/bin/bash", shell=True)
    log.info("RTKLIB RUN FINISHED for {}".format(exp_full_name))

    return out_result_fil
