#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.files_rw contains functions to
write misc. geodetic data in dedicated files.

it can be imported directly with:
from geodezyx import files_rw

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""

########## BEGIN IMPORT ##########
#### External modules
import datetime as dt

#### Import the logger
import logging
import os
import re

import numpy as np
import pandas as pd

#### geodeZYX modules
from geodezyx import conv
from geodezyx import files_rw
from geodezyx import reffram
from geodezyx import utils
from geodezyx import operational

# from geodezyx.megalib.megalib import *   # Import the legacy modules names
log = logging.getLogger("geodezyx")


#### Import star style
# from geodezyx import *                   # Import the GeodeZYX modules
# from geodezyx.externlib import *         # Import the external modules
# from geodezyx.megalib.megalib import *   # Import the legacy modules names
##########  END IMPORT  ##########


def write_sp3(
    orb_df_in,
    outpath,
    outname=None,
    prefix="orb",
    skip_null_epoch=True,
    force_format_c=False,
):
    """
    Write a SP3 file from an Orbit DataFrame

    Parameters
    ----------
    orb_df_in : DataFrame
        Input Orbit DataFrame.
    outpath : str
        The output path of the file (see also outname).
    outname : None or str, optional
        None = outpath is the full path (directory + filename) of the output.
        A string = a manual name for the file.
        'auto_old_cnv' = automatically generate the filename (old convention)
        'auto_new_cnv' = automatically generate the filename (new convention)
        The default is None.
    prefix : str, optional
        the output 3-char. name of the AC. The default is 'orb'.
    skip_null_epoch : bool, optional
        Do not write an epoch if all sats are null (filtering).
        The default is True.
    force_format_c : bool, optional
        force SP3's format c. The default is False.

    Returns
    -------
    The string containing the formatted SP3 data.
    """

    ################## MAIN DATA
    lines_stk = []

    sp3_df_wrk = orb_df_in.sort_values(["epoch", "prn"])

    epoch_raw_list = sp3_df_wrk["epoch"].unique()
    sat_list = sorted(sp3_df_wrk["prn"].unique())
    sat_list = list(sorted(sp3_df_wrk["prn"].unique()))
    ## sat_list    = list(reversed(sat_list))
    #### PS 210721
    #### reversed bc the sats are sorted ascending=False, but why???
    #### to have G before E ??
    sat_list_set = set(sat_list)
    epoch_used_list = []

    if not "clk" in sp3_df_wrk.columns:
        sp3_df_wrk["clk"] = 999999.999999

    for epoc in epoch_raw_list:
        sp3epoc = pd.DataFrame(sp3_df_wrk[sp3_df_wrk["epoch"] == epoc])

        ######## if keep_missing_sat_in_epoch:
        ## manage missing sats for the current epoc
        missing_sats = sat_list_set.difference(set(sp3epoc["prn"]))

        for miss_sat in missing_sats:
            miss_line = sp3epoc.iloc[0].copy()
            miss_line["prn"] = miss_sat
            miss_line["sys"] = miss_sat[0]
            ### check the sp3 doc
            # bad position = 0.000000
            # bad clock    = 999999.9999999
            miss_line["x"] = 0.000000
            miss_line["y"] = 0.000000
            miss_line["z"] = 0.000000
            miss_line["clk"] = 999999.999999

            sp3epoc = sp3epoc.append(miss_line)
        #### end of missing sat bloc

        sp3epoc.sort_values("prn", inplace=True, ascending=True)
        timestamp = conv.dt2sp3_timestamp(conv.numpy_dt2dt(epoc)) + "\n"

        linefmt = "p{:}{:14.6f}{:14.6f}{:14.6f}{:14.6f}\n"

        LinesStkEpoch = []
        sum_val_epoch = 0
        for ilin, lin in sp3epoc.iterrows():
            if not "clk" in lin.index:  # manage case if no clk in columns
                lin["clk"] = 999999.999999
            line_out = linefmt.format(
                lin["prn"], lin["x"], lin["y"], lin["z"], lin["clk"]
            )

            sum_val_epoch += lin["x"] + lin["y"] + lin["z"]

            LinesStkEpoch.append(line_out)

        ### if skip_null_epoch activated, print only if valid epoch
        if not (np.isclose(sum_val_epoch, 0) and skip_null_epoch):
            lines_stk.append(timestamp)  # stack the timestamp
            lines_stk = lines_stk + LinesStkEpoch  # stack the values
            epoch_used_list.append(epoc)  # stack the epoc as dt

    ################## HEADER
    ######### SATELLITE LIST

    satline_stk = []
    sigmaline_stk = []

    if force_format_c:
        nlines = 5
    else:
        div, mod = np.divmod(len(sat_list), 17)

        if div < 5:
            nlines = 5
        else:
            nlines = div

            if mod != 0:
                nlines += 1

    for i in range(nlines):
        sat_line = sat_list[17 * i : 17 * (i + 1)]
        sat_line_sigma = len(sat_line) * " 01"

        if len(sat_line) < 17:
            complem = " 00" * (17 - len(sat_line))
        else:
            complem = ""

        if i == 0:
            nbsat4line = len(sat_list)
        else:
            nbsat4line = ""

        satline = "+  {:3}   ".format(nbsat4line) + "".join(sat_line) + complem + "\n"
        sigmaline = "++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n"
        sigmaline = "++       " + sat_line_sigma + complem + "\n"

        satline_stk.append(satline)
        sigmaline_stk.append(sigmaline)

    ######### 2 First LINES
    start_dt = conv.numpy_dt2dt(np.min(epoch_used_list))

    header_line1 = (
        "#cP"
        + conv.dt2sp3_timestamp(start_dt, False)
        + "     {:3}".format(len(epoch_used_list))
        + "   u+U IGSXX FIT  XXX\n"
    )

    delta_epoch = int(utils.most_common(np.diff(epoch_used_list) * 10**-9))
    mjd = conv.dt2mjd(start_dt)
    mjd_int = int(np.floor(mjd))
    mjd_dec = mjd - mjd_int
    gps_wwww, gps_sec = conv.dt2gpstime(start_dt, True, "gps")

    header_line2 = "## {:4} {:15.8f} {:14.8f} {:5} {:15.13f}\n".format(
        gps_wwww, gps_sec, delta_epoch, mjd_int, mjd_dec
    )

    ######### HEADER BOTTOM
    header_bottom = """%c M  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
%c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
%f  1.2500000  1.025000000  0.00000000000  0.000000000000000
%f  0.0000000  0.000000000  0.00000000000  0.000000000000000
%i    0    0    0    0      0      0      0      0         0
%i    0    0    0    0      0      0      0      0         0
/* PCV:IGSXX_XXXX OL/AL:FESXXXX  NONE     YN CLK:CoN ORB:CoN
/*     GeodeZYX Toolbox Output
/*
/*
"""

    ################## FINAL STACK

    final_lines_stk = []

    final_lines_stk.append(header_line1)
    final_lines_stk.append(header_line2)
    final_lines_stk = final_lines_stk + satline_stk + sigmaline_stk
    final_lines_stk.append(header_bottom)
    final_lines_stk = final_lines_stk + lines_stk + ["EOF"]

    final_str = "".join(final_lines_stk)

    ### Manage the file path
    prefix_opera = prefix

    if not outname:
        outpath_opera = outpath
    elif outname == "auto_old_cnv":
        week, dow = conv.dt2gpstime(start_dt)
        filename = prefix_opera + str(week) + str(dow) + ".sp3"
        outpath_opera = os.path.join(outpath, filename)

    elif outname == "auto_new_cnv":
        log.error("not implemented yet !!!!!")
        raise Exception

    f = open(outpath_opera, "w+")
    f.write(final_str)


def write_clk(
    df_clk_in, outpath, outname=None, prefix="orb", header="", output_std_values=False
):
    """
    Write a SP3 Clock file from an Clock DataFrame

    Parameters
    ----------
    df_clk_in : DataFrame
        Input Clock DataFrame.
    outpath : str
        The output path of the file (see also outname).
    outname : None or str, optional
        None = outpath is the full path (directory + filename) of the output.
        A string = a manual name for the file.
        'auto_old_cnv' = automatically generate the filename (old convention)
        'auto_new_cnv' = automatically generate the filename (new convention)
        The default is None.
    prefix : str, optional
        the output 3-char. name of the AC. The default is 'orb'.
    header : str, optional
        A string describing the clk file header. The default is "".
    output_std_values : bool, optional
        Add observation sigmas as the last column. The default is False.

    Returns
    -------
    The string containing the formatted clock data.
    """

    HEAD = header
    row_str_stk = []

    if output_std_values:
        row_str_proto = "{:2} {:4} {:4d} {:02d} {:02d} {:02d} {:02d} {:9.6f} {:2d}   {:19.12e} {:19.12e}"
    else:
        row_str_proto = (
            "{:2} {:4} {:4d} {:02d} {:02d} {:02d} {:02d} {:9.6f} {:2d}   {:19.12e}"
        )

    for irow, row in df_clk_in.iterrows():

        if output_std_values:
            one_or_two = 2
            row_str = row_str_proto.format(
                row["type"],
                row["name"],
                int(row["year"]),
                int(row["month"]),
                int(row["day"]),
                int(row["hour"]),
                int(row["minute"]),
                int(row["second"]),
                one_or_two,
                row["bias"],
                row["sigma"],
            )
        else:
            one_or_two = 1
            row_str = row_str_proto.format(
                row["type"],
                row["name"],
                int(row["year"]),
                int(row["month"]),
                int(row["day"]),
                int(row["hour"]),
                int(row["minute"]),
                int(row["second"]),
                one_or_two,
                row["bias"],
            )
        row_str_stk.append(row_str)

    ## Add EOF
    row_str_stk.append("EOF")

    corpse = "\n".join(row_str_stk)

    out = (
        HEAD
        + "                                                            END OF HEADER\n"
        + corpse
    )

    ### Manage the file path
    prefix_opera = prefix

    if not outname:
        outpath_opera = outpath
    elif outname == "auto_old_cnv":
        start_dt = dt.datetime(
            int(df_clk_in.iloc[0]["year"]),
            int(df_clk_in.iloc[0]["month"]),
            int(df_clk_in.iloc[0]["day"]),
        )
        week, dow = conv.dt2gpstime(start_dt)
        filename = prefix_opera + str(week) + str(dow) + ".clk"
        outpath_opera = os.path.join(outpath, filename)

    elif outname == "auto_new_cnv":
        log.error("not implemented yet !!!!!")
        raise Exception

    else:
        outpath_opera = os.path.join(outpath, outname)

    out = out

    with open(outpath_opera, "w+") as Fout:
        Fout.write(out)
        Fout.close()

    return out


def ine_block_mono(sat, dt_in, extra_intrvl_strt=0.1, extra_intrvl_end=0.4, step=300):
    """
    Write an EPOS INE block
    """

    fields = [
        "orb____1",
        "orb____2",
        "orb____3",
        "orb____4",
        "orb____5",
        "orb____6",
        "orb___db",
        "orb_s2db",
        "orb_c2db",
        "orb_s4db",
        "orb_c4db",
        "orb___yb",
        "orb___xb",
        "orb_sixb",
        "orb_coxb",
        "orb___cr",
    ]

    mjd = np.floor(conv.dt2mjd(dt_in))
    mjd_strt = mjd - extra_intrvl_strt
    mjd_end = mjd + extra_intrvl_end + 1

    lines = []

    l1 = " sat_nr  : " + sat + "\n"
    l2 = " stepsize: {:3}  {:6.2f}\n".format(sat, step)

    lines.append(l1)
    lines.append(l2)

    for field in fields:
        line = " {:}: {:3}  0.000000000000000E+00 {:11.5f} {:11.5f}\n".format(
            field, sat, mjd_strt, mjd_end
        )
        lines.append(line)

    lines.append(" end_sat\n")

    str_out = "".join(lines)

    return str_out


def write_ine_dummy_file(
    Sat_list,
    dt_in,
    extra_intrvl_strt=0.1,
    extra_intrvl_end=0.4,
    step=300,
    out_file_path=None,
):
    """
    Write an EPOS INE dummy (empty values) file
    """

    Lines = []

    mjd = np.floor(conv.dt2mjd(dt_in))
    mjd_strt = mjd - extra_intrvl_strt
    mjd_end = mjd + extra_intrvl_end + 1

    datestr = conv.dt2str(dt.datetime.now(), str_format="%Y/%m/%d %H:%M:%S")

    mjd_strt_deci = mjd_strt - np.floor(mjd_strt)

    head_proto = """%=INE 1.00 {:} NEWSE=INE+ORBCOR                                                                                 
+global
 day_info: 
 epoch   :                            {:5}  {:16.14f}
 interval:                            {:11.5f} {:11.5f}
 stepsize:      {:6.2f}
-global
+initial_orbit
"""
    head = head_proto.format(datestr, int(mjd), 0, mjd_strt, mjd_end, step)

    Lines.append(head)

    for sat in Sat_list:
        Lines.append(
            "******************************************************************\n"
        )
        sat_str = ine_block_mono(sat, dt_in, extra_intrvl_strt, extra_intrvl_end, step)
        Lines.append(sat_str)
        Lines.append(
            "******************************************************************\n"
        )

    str_end = """-initial_orbit
%ENDINE
"""

    Lines.append(str_end)

    str_out = "".join(Lines)

    if out_file_path:
        with open(out_file_path, "w") as f:
            f.write(str_out)
            f.close()

    return str_out


def sp3_overlap_creator(
    ac_list,
    dir_in,
    dir_out,
    suffix_out_input=None,
    overlap_size=7200,
    force=False,
    manage_missing_sats="exclude_missing_epoch",
    eliminate_null_sat=True,
    severe=False,
    separated_systems_export=False,
    first_date=None,
    end_date=None,
    # new_naming = False,
    exclude_bad_epoch=True,
    sys=None,
    new_name=False,
):
    """
    Generate an SP3 Orbit file with overlap based on the SP3s of the
    days before and after

    Parameters
    ----------
    ac_list : list
        3-character codes of the ACs.
    dir_in : str
        where the input sp3 are.
    dir_out : str
         where the output sp3 will be outputed.
    suffix_out_input : str, optional
        last char of the 3-char. code. if None, then it is the same as input.
    overlap_size : int, optional
        Overlapsize. The default is 7200.
    force : True, optional
        force overwrite. The default is False.
    manage_missing_sats : str, optional
        'exclude_missing_day' : generate a file with only the common sat
        between the 3 days. Thus, exclude the missing sats for a complete day\n
        'exclude_missing_epoch' : generate a file with only sat with full epochs\n
        'extrapolate' : extrapolate the missing sats based on the first/last epoch\n
        The default is 'exclude_missing_epoch'.
    eliminate_null_sat : bool, optional
        eliminate null sat. The default is True.
    severe : bool, optional
        raise an exception if problem. The default is False.
    separated_systems_export : bool, optional
        export different sp3 for different system. The default is False.
    first_date : datetime, optional
        exclude SP3 before this epoch
    end_date :datetime, optional
        exclude SP3 after this epoch
    exclude_bad_epoch : bool, optional
        remove bad epoch (usually filled with 99999999.9999999 or 0.00000000)
    sys : string, optional
        if just to keep one system (e.g.: G)

    Returns
    -------
    None.

    Note
    ----
    start/end date are not implemented
    the force option skips existing files
    not implemented for new names
    """

    Dict_Lfiles_ac = dict()

    for ac in ac_list:
        Dict_Lfiles_ac[ac] = []
        lfile = Dict_Lfiles_ac[ac]

        if new_name:
            lfile = operational.find_igs_products_files(
                dir_in,
                ["sp3"],
                ac_list,
                first_date,
                date_end=end_date,
                severe=False,
                recursive_search=True,
                regex_old_naming=True,
                regex_new_naming=True,
                regex_igs_tfcc_naming=False,
                compressed="incl",
            )

        else:
            # extlist = ["sp3","SP3","sp3.gz","SP3.gz"]
            extlist = ["sp3", "sp3.gz", "eph", "SP3", "SP3.gz"]
            for ext in extlist:
                lfile = lfile + utils.find_recursive(dir_in, "*" + ac + "*" + ext)

        log.info("Nb of SP3 found for %s %s", ac, len(lfile))

        if not suffix_out_input:
            suffix_out = ac
        else:
            suffix_out = (
                ac[:2] + suffix_out_input
            )  ## GM  2021-09-06 this is not good for the new namning convention

        d = []
        wwwwd = []

        new_name_list = []

        regex_suffix = "[0-9]{4}_[0-9]{2}[A-Z]_[0-9]{2}[A-Z]_ORB.SP3"

        for sp3 in lfile:
            # wwwwd_str = os.path.basename(sp3)[3:8]
            # d.append(conv.gpstime2dt(int(wwwwd_str[:4]),int(wwwwd_str[4:])))
            extension = sp3[-3:]
            dat = conv.sp3name2dt(sp3)
            d.append(dat)

            if re.search(regex_suffix, sp3):
                new_name_list.append(True)
            else:
                new_name_list.append(False)

        #        ## just a check for the files list. May be excluded in the future.
        #        for item in lfile:
        #            if 'COD0MGXFIN_20200050000_01D_05M_ORB.SP3' in item:
        #                print(lfile.index(item))

        for dat, newnamebool in zip(
            d[1:-1], new_name_list[1:-1]
        ):  ####if selection manuel, zip > 2lists !!!
            #        for dat in d[875:876]: ####if selection manuel, zip > 2lists !!!

            dummy = 1
            try:
                log.info("*********** %s %s", ac, dat)

                if first_date:
                    if dat < first_date:
                        log.info("SKIP date", dat)
                        continue
                if end_date:
                    if dat > end_date:
                        log.info("SKIP date", dat)
                        continue

                wwwwd_str = conv.dt_2_sp3_datestr(dat).zfill(5)

                dat_bef = dat - dt.timedelta(days=1)
                dat_aft = dat + dt.timedelta(days=1)

                wwwwd_str_bef = utils.join_improved(
                    "", *conv.dt2gpstime(dat_bef)
                ).zfill(5)
                wwwwd_str_aft = utils.join_improved(
                    "", *conv.dt2gpstime(dat_aft)
                ).zfill(5)

                ###### check if exists
                dir_out_wk = os.path.join(dir_out, "wk" + str(wwwwd_str)[:4])
                utils.create_dir(dir_out_wk)
                fil_out = dir_out_wk + "/" + suffix_out + wwwwd_str + ".sp3"

                if not force and os.path.isfile(fil_out):
                    log.info("0)) " + fil_out + " exists, skipping...")
                    continue

                ### *************** STEP 1 ***************
                log.info("1)) Search for the days before/after")
                log.info("1)) %s %s", dat_bef, dat_aft)

                ## GM 2021-09-06 in case of the new naming conventation this is step needs to be done:
                ## not the best way to handle this, since it is hardcoded and the step needs to be given in a "if"
                ## maybe try to find a better way with regex!
                # if new_naming:
                #     if ac in ['WUM','GRG','SHA']:
                #         step_str = '15'
                #     else:
                #         step_str = '05'

                #     year,day = conv.dt_to_doy(conv.gpstime2dt(int(wwwwd_str[:4]),int(wwwwd_str[-1])))
                #     p1    = utils.find_regex_in_list(str(year)+str(day).zfill(3)  + "0000_01D_"+step_str+"M_ORB.SP3",lfile,True)

                #     year_bef,day_bef = conv.dt_to_doy(conv.gpstime2dt(int(wwwwd_str_bef[:4]),int(wwwwd_str_bef[-1])))
                #     p_bef = utils.find_regex_in_list(str(year_bef)+str(day_bef).zfill(3)  + "0000_01D_"+step_str+"M_ORB.SP3",lfile,True)

                #     year_aft,day_aft = conv.dt_to_doy(conv.gpstime2dt(int(wwwwd_str_aft[:4]),int(wwwwd_str_aft[-1])))
                #     p_aft = utils.find_regex_in_list(str(year_aft)+str(day_aft).zfill(3)  + "0000_01D_"+step_str+"M_ORB.SP3",lfile,True)

                # if re.search(regex_suffix,)
                if newnamebool:
                    day, year = conv.dt2doy_year(dat)
                    regex_prefix = str(year) + str(day).zfill(3)
                    p1 = utils.find_regex_in_list(
                        regex_prefix + regex_suffix, lfile, True
                    )

                    day_bef, year_bef = conv.dt2doy_year(dat_bef)
                    regex_prefix_bef = str(year_bef) + str(day_bef).zfill(3)
                    p_bef = utils.find_regex_in_list(
                        regex_prefix_bef + regex_suffix, lfile, True
                    )

                    day_aft, year_aft = conv.dt2doy_year(dat_aft)
                    regex_prefix_aft = str(year_aft) + str(day_aft).zfill(3)
                    p_aft = utils.find_regex_in_list(
                        regex_prefix_aft + regex_suffix, lfile, True
                    )

                else:
                    if extension == "eph":
                        p1 = utils.find_regex_in_list(wwwwd_str + ".eph", lfile, True)
                        p_bef = utils.find_regex_in_list(
                            wwwwd_str_bef + ".eph", lfile, True
                        )
                        p_aft = utils.find_regex_in_list(
                            wwwwd_str_aft + ".eph", lfile, True
                        )
                    else:
                        p1 = utils.find_regex_in_list(wwwwd_str + ".sp3", lfile, True)
                        p_bef = utils.find_regex_in_list(
                            wwwwd_str_bef + ".sp3", lfile, True
                        )
                        p_aft = utils.find_regex_in_list(
                            wwwwd_str_aft + ".sp3", lfile, True
                        )

                log.info("1)) Files found for the days before/after")
                log.info("0b) %s", p_bef)
                log.info("01) %s", p1)
                log.info("0a) %s", p_aft)

                if not p1 or not p_bef or not p_aft:
                    log.error("with day %s", dat)
                    continue

                sp3 = files_rw.read_sp3(p1)
                sp3_bef = files_rw.read_sp3(p_bef)
                sp3_aft = files_rw.read_sp3(p_aft)

                ### Filtering to keep p only
                sp3 = sp3[sp3.type == "p"]
                sp3_bef = sp3_bef[sp3_bef.type == "p"]
                sp3_aft = sp3_aft[sp3_aft.type == "p"]

                if sys:
                    sp3 = sp3[sp3.sys == sys]
                    sp3_bef = sp3_bef[sp3_bef.sys == sys]
                    sp3_aft = sp3_aft[sp3_aft.sys == sys]

                sp3_bef = sp3_bef[sp3_bef["epoch"] < sp3["epoch"].min()]
                sp3_aft = sp3_aft[sp3_aft["epoch"] > sp3["epoch"].max()]

                sp3concat = pd.concat((sp3_bef, sp3, sp3_aft))

                dat_filter_bef = dat - dt.timedelta(seconds=overlap_size)
                dat_filter_aft = (
                    dat + dt.timedelta(seconds=overlap_size) + dt.timedelta(days=1)
                )

                ### *************** STEP 2 ***************
                log.info("2)) dates of the overlap period before/after")
                log.info("2)) %s %s", dat_filter_bef, dat_filter_aft)

                ### *************** STEP 3 ***************
                log.info("3)) Dates of: SP3 concatenated, before, current, after")
                log.info(
                    "3)) %s %s", sp3concat["epoch"].min(), sp3concat["epoch"].max()
                )
                log.info("3b) %s %s", sp3_bef["epoch"].min(), sp3_bef["epoch"].max())
                log.info("31) %s %s", sp3["epoch"].min(), sp3["epoch"].max())
                log.info("3a) %s %s", sp3_aft["epoch"].min(), sp3_aft["epoch"].max())

                sp3concat = sp3concat[
                    (sp3concat["epoch"] >= dat_filter_bef)
                    & (sp3concat["epoch"] <= dat_filter_aft)
                ]

                if exclude_bad_epoch:
                    good_epochs_bool_999 = sp3concat[["x", "y", "z"]] < 999999.0
                    good_epochs_bool_000 = np.logical_not(
                        np.isclose(sp3concat[["x", "y", "z"]], 0.0)
                    )

                    good_epochs_bool = np.logical_and(
                        good_epochs_bool_999, good_epochs_bool_000
                    )
                    good_epochs_bool = np.all(good_epochs_bool, axis=1)

                    sp3concat = sp3concat[good_epochs_bool]

                ########## HERE WE MANAGE THE MISSING SATS
                if manage_missing_sats == "exclude_missing_day":
                    log.info("4))", "remove missing sats -- day")
                    common_sats = (
                        set(sp3_bef["prn"])
                        .intersection(set(sp3["prn"]))
                        .intersection(set(sp3_aft["prn"]))
                    )
                    sp3concat = sp3concat[sp3concat["prn"].isin(common_sats)]

                elif manage_missing_sats == "exclude_missing_epoch":
                    log.info("4))", "remove missing sats -- epoch")
                    nepoc = len(sp3concat["epoch"].unique())
                    sp3concat_satgrp = sp3concat.groupby("prn")

                    all_sats = sp3concat["prn"].unique()
                    good_sats = sp3concat_satgrp.count() == nepoc

                    ###### good_sats = good_sats.reset_index()["sat"]
                    ## we get the good sats based one column containing a boolean
                    ## because of the test just before (abitrarily epoch column)
                    ## and after get the corresponding good sats names
                    good_sats = good_sats[good_sats["epoch"]].reset_index()["prn"]

                    bad_sats = list(set(all_sats) - set(good_sats))
                    log.info("excluded bad sats: %s", bad_sats)

                    sp3concat = sp3concat[sp3concat["prn"].isin(good_sats)]

                elif manage_missing_sats == "extrapolate":
                    log.info("4)) extrapolate missing sats ")
                    for iovl, SP3_ovl in enumerate((sp3_bef, sp3_aft)):
                        if iovl == 0:
                            backward = True
                            forward = False
                            backfor = "backward"
                        elif iovl == 1:
                            backward = False
                            forward = True
                            backfor = "forward"

                        sats = set(sp3["prn"])
                        sats_ovl = set(SP3_ovl["prn"])

                        sats_miss = sats.difference(sats_ovl)
                        if not sats_miss:
                            continue
                        log.info(
                            "4a) extrapolate missing sats %s %s", backfor, sats_miss
                        )

                        sp3extrapo_in = sp3concat[sp3concat["prn"].isin(sats_miss)]

                        # step = utils.most_common(sp3concat["epoch"].diff().dropna())
                        # step = step.astype('timedelta64[s]').astype(np.int32)
                        step = 900
                        # print(step)

                        # print("sp3extrapo_in",sp3extrapo_in)

                        sp3extrapo = reffram.extrapolate_sp3_data_frame(
                            sp3extrapo_in,
                            step=step,
                            n_step=int(overlap_size / step),
                            backward=backward,
                            forward=forward,
                            until_backward=dat_filter_bef,
                            until_forward=dat_filter_aft,
                            return_all=False,
                        )

                        sp3concat = pd.concat((sp3concat, sp3extrapo))
                        log.info(sp3extrapo)

                else:
                    log.error("check manage_missing_sats value")
                    raise Exception

                if eliminate_null_sat:
                    good_sats = []
                    for sat in sp3concat["prn"].unique():
                        xyzvals = sp3concat[sp3concat["prn"] == sat][
                            ["x", "y", "z"]
                        ].sum(axis=1)

                        v = np.sum(np.isclose(xyzvals, 0)) / len(xyzvals)

                        if v < 0.50:
                            good_sats.append(sat)
                        else:
                            log.info("6) eliminate because null position %s", sat)

                    sp3concat = sp3concat[sp3concat["prn"].isin(good_sats)]

                ### *************** STEP 7 ***************
                log.info("7)) Start/End Epoch of the concatenated file ")
                log.info(
                    "7)) %s %s", sp3concat["epoch"].min(), sp3concat["epoch"].max()
                )

                #### All systems
                log.info("8)) outputed file")
                log.info(fil_out)
                write_sp3(sp3concat, fil_out)

                #### system separated
                if False:
                    for sys in sp3concat["sys"].unique():
                        try:
                            sp3concat_sys = sp3concat[sp3concat["sys"] == sys]
                            fil_out_sys = (
                                dir_out_wk
                                + "/"
                                + suffix_out[:2]
                                + sys.lower()
                                + wwwwd_str.zfill(5)
                                + ".sp3"
                            )
                            log.info("9)) outputed file")
                            log.info(fil_out_sys)
                            write_sp3(sp3concat_sys, fil_out_sys)
                        except:
                            continue

            except KeyboardInterrupt:
                raise KeyboardInterrupt

            except Exception as e:
                if severe:
                    log.error(e)
                    raise e
                else:
                    log.warning("Error %s but no severe mode, continue...", e)
                    continue

    """
    sort_wrt="site" or "site_num"
    
    soln_in_df
    use soln AND pt information in the input DataFrame
    """


def write_epos_sta_coords(
    df_inp,
    file_out,
    sort_wrt="site",
    no_time_limit_for_first_period=True,
    no_time_limit_for_last_period=True,
    soln_in_df=True,
    trf_name="xTRFnn",
):
    """
    Write an EPOS coordinate file

    Parameters
    ----------
    df_inp : DataFrame
        Input Orbit DataFrame.
    file_out : str
        The output path of the file.
    sort_wrt : bool, optional
        Sort the values with respect to a DF column.
        The default is "site".
    no_time_limit_for_first_period : bool, optional
        No time limit for the first period.
        The default is True.
    no_time_limit_for_last_period : bool, optional
        No time limit for the last period.
        The default is True.
    soln_in_df : bool, optional
        Soln in DF.
        The default is True.

    Returns
    -------
    None.

    """

    df_work = df_inp.sort_values([sort_wrt, "mjd_start"])

    stat_lines_blk_stk = []

    generic_header = """+info
 FLATTENING                  298.2550
 MAJOR_AXIS              6378140.0000
 REFERENCE_FRAME                IGS14
 NUMBER_OF_STATIONS             {:5d}
 REF_MJD                        {:5d}
-info
"""

    generic_header = generic_header.format(
        len(df_work["site_num"].unique()), int(utils.most_common(df_work["MJD_ref"]))
    )

    stat_lines_blk_stk.append(generic_header)

    stat_lines_blk_stk.append("+station_coordinates")

    for site in df_work[sort_wrt].unique():

        stat_lines_blk_stk.append(
            "*------------------------- ---- ----- -beg- -end- -**- ------------------------------------------------\n*"
        )

        df_site_block = df_work[df_work[sort_wrt] == site]

        df_site_block.reset_index(inplace=True)

        for i_l, (_, l) in enumerate(df_site_block.iterrows()):

            if soln_in_df:
                iope = int(l["soln"])
                pt = l["pt"]
            else:
                iope = i_l + 1
                pt = "A"

            if no_time_limit_for_first_period and i_l == 0:
                mjd_start = 0
            else:
                mjd_start = l["mjd_start"]

            if no_time_limit_for_last_period and (i_l + 1) == len(df_site_block):
                mjd_end = 0
            else:
                mjd_end = l["mjd_end"]

            line_site_fmt = " SITE            m {:4d}  {:1d} {:} {:5d} {:5d} {:5d} {:}   {:}  {:1d}      LOG_CAR       LOG_CAR"
            line_site_fmt = " SITE            m {:4d}  {:1d} {:} {:5d} {:5d} {:5d} {:}   {:}  {:1d}{:>13}       {:<13}"
            line_valu_fmt = " POS_VEL:XYZ     m {:4d}  {:1d} {:+15.4f} {:+15.4f} {:+15.4f}      {:+6.4f} {:+6.4f} {:+6.4f}"
            line_sigm_fmt = " SIG_PV_XYZ      m {:4d}  {:1d} {:+15.4f} {:+15.4f} {:+15.4f}      {:+6.4f} {:+6.4f} {:+6.4f}"

            line_site = line_site_fmt.format(
                int(l["site_num"]),
                int(iope),
                l["tecto_plate"].upper(),
                int(l["MJD_ref"]),
                int(mjd_start),
                int(mjd_end),
                l["site"],
                pt,
                int(iope),
                trf_name,
                trf_name,
            )

            line_valu = line_valu_fmt.format(
                int(l["site_num"]),
                int(iope),
                l["x"],
                l["y"],
                l["z"],
                l["Vx"],
                l["Vy"],
                l["Vz"],
            )

            line_sigm = line_sigm_fmt.format(
                int(l["site_num"]),
                int(iope),
                l["sx"],
                l["sy"],
                l["sz"],
                l["sVx"],
                l["sVy"],
                l["sVz"],
            )

            stat_lines_blk_stk.append(line_site)
            stat_lines_blk_stk.append(line_valu)
            stat_lines_blk_stk.append(line_sigm)
            stat_lines_blk_stk.append("*")

    stat_lines_blk_stk.append("-station_coordinates")

    final_str = "\n".join(stat_lines_blk_stk)

    with open(file_out, "w+") as f:
        f.write(final_str)

    return final_str


# def write_sndy_light_dat(ts_in,outdir,outprefix):
#     """ Not properly implemented """
#     fil = open(os.path.join(outdir,outprefix),'w+')
#     if isinstance(ts_in,TimeSeriePoint):
#         if ts_in.initype() == 'FLH':
#             for pt in ts_in.pts:
#                 lin = ' '.join([str(e) for e in [pt.F , pt.L , pt.H , pt.T , pt.sF , pt.sL , pt.sH ]])
#                 fil.write(lin + '\n')
#     elif isinstance(ts_in,TimeSerieObs):
#         if ts_in.typeobs == 'RPY':
#             for att in ts_in.obs:
#                 lin = ' '.join([str(e) for e in [att.R , att.P , att.Y , att.T , att.Q.w , att.Q.x , att.Q.y , att.Q.z ]])
#                 fil.write(lin + '\n')
#     fil.close()
