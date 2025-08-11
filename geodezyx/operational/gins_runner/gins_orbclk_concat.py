#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 10/08/2025 11:38:34

@author: psakic
"""

from geodezyx import utils, conv
import geodezyx.operational.gins_runner.gins_common as gynscmn
import os
import subprocess
import datetime as dt
import multiprocessing as mp
import logging
# Import the logger
log = logging.getLogger("geodezyx")

def concat_orb_clk(date_srt, date_end, nprocs=1, prod="G20", verbose=True):
    """
    Concatenate orbit and clock files for a given date range using multiprocessing.

    Parameters
    ----------
    date_srt : datetime
        Start date for the concatenation process.
    date_end : datetime
        End date for the concatenation process.
    nprocs : int
        Number of processes to use for multiprocessing.
    prod : str, optional
        Product type (e.g., "G20"). Default is "G20".
    verbose : bool, optional
        If True, enable verbose logging. Default is True.

    Returns
    -------
    None
    """

    def _chk_cat_orbclk(orbclk_out, silent=False):
        """
        Check if the concatenated orbit/clock file exists.

        Parameters
        ----------
        orbclk_out : str
            Path to the output orbit/clock file.
        silent : bool, optional
            If True, suppress logging. Default is False.

        Returns
        -------
        bool
            True if the file exists, False otherwise.
        """
        orbclk_out = orbclk_out + ".gz"
        if os.path.isfile(orbclk_out):
            log.info("%s here :)", orbclk_out)
            return True
        elif not silent:
            log.error("%s not here :(", orbclk_out)
            return False
        else:
            return False

    def _run_cat_orbclk(cmd, out_fil):
        """
        Run the command to generate the concatenated orbit/clock file.

        Parameters
        ----------
        cmd : str
            Command to execute.
        out_fil : str
            Path to the output file.

        Returns
        -------
        None
        """
        if _chk_cat_orbclk(out_fil, silent=True):
            return None

        if verbose:
            log.info(cmd)
        subprocess.run(
            cmd,
            executable="/bin/bash",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True,
        )
        utils.gzip_compress(out_fil, rm_inp=True)
        _chk_cat_orbclk(out_fil)
        return None

    global cat_orbclk_wrap

    def cat_orbclk_wrap(date_inp):
        """
        Wrapper function to process a single date for orbit/clock concatenation.

        Parameters
        ----------
        date_inp : datetime
            Date to process.

        Returns
        -------
        tuple
            Paths to the concatenated orbit and clock files.
        """
        jjul_srt = str(conv.dt2jjul_cnes(date_inp - dt.timedelta(days=0)))
        jjul_end = str(conv.dt2jjul_cnes(date_inp + dt.timedelta(days=2)))
        gs_user = gynscmn.get_gin_path(extended=False)
        gins_data = gynscmn.get_gin_path(extended=True) + "/data"

        orb_out = os.path.join(
            gs_user, "GPSDATA", "_".join((prod + "ORB", "AUTOM", jjul_srt, jjul_end))
        )
        cmd_orb = " ".join(
            ["rapat_orb_gnss.sh", jjul_srt, jjul_end, "3", gins_data, orb_out, "0"]
        )

        clk_out = os.path.join(
            gs_user, "GPSDATA", "_".join((prod, "AUTOM", jjul_srt, jjul_end))
        )
        cmd_clk = " ".join(["get_hor_hautes", jjul_srt, jjul_end, prod, clk_out])

        for cmd, out_fil in [(cmd_orb, orb_out), (cmd_clk, clk_out)]:
            _run_cat_orbclk(cmd, out_fil)

        return orb_out, clk_out

    # Generate a list of dates to process
    date_lis = conv.dt_range(date_srt, date_end)

    # Use multiprocessing to process the dates
    pool = mp.Pool(processes=nprocs)
    try:
        _ = pool.map(cat_orbclk_wrap, date_lis, chunksize=1)
    except Exception as e:
        log.error("error in the pool.map : %s", e)
    pool.close()

    return None


def main():
    """
    Command-line interface for concatenating orbit and clock files.
    """
    import argparse
    from datetime import datetime

    parser = argparse.ArgumentParser(
        description="Concatenate orbit and clock files for a given date range using multiprocessing."
    )
    parser.add_argument('-s', '--date_srt', type=str, required=True,
                        help='Start date for the concatenation process (format: YYYY-MM-DD)')
    parser.add_argument('-e', '--date_end', type=str, required=True,
                        help='End date for the concatenation process (format: YYYY-MM-DD)')
    parser.add_argument('-n', '--nprocs', type=int, default=1,
                        help='Number of processes to use for multiprocessing (optional). Default is 1')
    parser.add_argument('-p', '--prod', type=str, default="G20",
                        help='Product type (optional). Default is "G20"')
    parser.add_argument('-v', '--verbose', action='store_true', default=True,
                        help='Enable verbose logging (optional). Default is True')

    args = parser.parse_args()

    # Convert string dates to datetime objects
    try:
        date_srt = datetime.strptime(args.date_srt, '%Y-%m-%d')
        date_end = datetime.strptime(args.date_end, '%Y-%m-%d')
    except ValueError as e:
        print(f"Error parsing dates: {e}")
        print("Please use YYYY-MM-DD format for dates")
        return

    # Call the function
    concat_orb_clk(
        date_srt=date_srt,
        date_end=date_end,
        nprocs=args.nprocs,
        prod=args.prod,
        verbose=args.verbose
    )

if __name__ == "__main__":
    main()