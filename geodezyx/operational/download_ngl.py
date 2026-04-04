#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NGL (Nevada Geodetic Laboratory) tropospheric time-series downloader.

This sub-module of geodezyx.operational contains functions to download,
extract, and read GNSS tropospheric solutions from the UNR/NGL server:
    https://geodesy.unr.edu/gps_timeseries/

URL pattern for tropospheric files (one zip per station per year):
    https://geodesy.unr.edu/gps_timeseries/<frame>/trop/<STAT>/<STAT>.<YYYY>.trop.zip

Each zip contains daily SINEX tropospheric files compressed as .gz:
    <STAT>.<YYYY>.<DOY>.trop.gz

The data is read into pandas DataFrames using
geodezyx.files_rw.read.read_atmo.read_snx_trop (version=1).

The geodezyx toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License
Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""

########## BEGIN IMPORT ##########
import logging
import os
import zipfile

import pandas as pd

from geodezyx.files_rw import unzip_gz_z
from geodezyx.files_rw.read.read_atmo import read_snx_trop
from geodezyx.operational.download_gnss.download_utils import download_http

##########  END IMPORT  ##########

log = logging.getLogger("geodezyx")

# Base URL for the NGL GPS tropospheric time series
_NGL_BASE_URL = "https://geodesy.unr.edu/gps_timeseries"



# ============================================================
#  Low-level helpers
# ============================================================


def _ngl_trop_url(station: str, year: int, frame: str = "IGS20") -> str:
    """
    Build the URL for a NGL yearly tropospheric zip archive.

    Parameters
    ----------
    station : str
        4-character station name (case-insensitive; uppercased automatically).
    year : int
        Four-digit year.
    frame : str, optional
        Reference frame label used in the URL path.  Default is ``"IGS20"``.

    Returns
    -------
    str
        Full HTTPS URL pointing to the zip archive.

    Examples
    --------
    >>> _ngl_trop_url("HOUZ", 2000)
    'https://geodesy.unr.edu/gps_timeseries/IGS20/trop/HOUZ/HOUZ.2000.trop.zip'
    """
    stat = station.upper()
    return f"{_NGL_BASE_URL}/{frame}/trop/{stat}/{stat}.{year}.trop.zip"


def _extract_zip(zip_path: str, extract_dir: str) -> list:
    """
    Extract a zip archive and return the list of extracted file paths.

    Parameters
    ----------
    zip_path : str
        Path to the ``.zip`` file.
    extract_dir : str
        Directory into which the contents will be extracted.

    Returns
    -------
    list of str
        Absolute paths of all extracted files.
    """
    os.makedirs(extract_dir, exist_ok=True)
    extracted = []
    with zipfile.ZipFile(zip_path, "r") as zf:
        for member in zf.namelist():
            zf.extract(member, extract_dir)
            extracted.append(os.path.join(extract_dir, member))
    log.info("Extracted %d file(s) from %s", len(extracted), os.path.basename(zip_path))
    return extracted


# ============================================================
#  Main public function
# ============================================================


def download_ngl_trop(
    stations,
    year_start: int,
    year_end: int,
    output_dir: str,
    mode: str = "download",
    frame: str = "IGS20",
    keep_zip: bool = True,
    keep_gz: bool = True,
    timeout: int = 120,
    max_try: int = 4,
):
    """
    Download (and optionally process) NGL tropospheric SINEX solutions.

    Data source
    -----------
    Nevada Geodetic Laboratory (UNR):
    ``https://geodesy.unr.edu/gps_timeseries/<frame>/trop/<STAT>/<STAT>.<YYYY>.trop.zip``

    Each yearly ``.zip`` archive contains one ``.gz`` file per day of year,
    e.g. ``HOUZ.2000.136.trop.gz``.  Each ``.gz`` holds a standard SINEX
    tropospheric solution that is read with
    :func:`geodezyx.files_rw.read_atmo.read_snx_trop` (``version=1``).

    Parameters
    ----------
    stations : str or list of str
        One station name or a list of 4-character station names
        (e.g. ``"HOUZ"`` or ``["HOUZ", "REYK"]``).
    year_start : int
        First year to download (inclusive).
    year_end : int
        Last year to download (inclusive).
    output_dir : str
        Root directory for all downloaded / extracted / parquet files.
        Sub-directories are created automatically following the layout::

            <output_dir>/
                zips/           # yearly .zip archives
                trop/           # extracted .gz files (per station/year)
                decompressed/   # decompressed SINEX files
                parquet/        # daily parquet files (mode "parquet" only)

    mode : {"download", "decompress", "parquet"}, optional
        Processing mode:

        ``"download"``
            Only download the yearly ``.zip`` archives.
        ``"decompress"``
            Download **and** extract the ``.zip`` archives, then decompress
            each daily ``.gz`` file into a plain SINEX ``.trop`` file.
        ``"parquet"``
            Download, extract, decompress, **and** read each daily SINEX file
            into a :class:`pandas.DataFrame`, then save it as a ``.parquet``
            file.

        Default is ``"download"``.
    frame : str, optional
        Reference frame label in the URL path.  Default is ``"IGS20"``.
    keep_zip : bool, optional
        Whether to keep the ``.zip`` archives after extraction.
        Only relevant when *mode* is ``"decompress"`` or ``"parquet"``.
        Default is ``True``.
    keep_gz : bool, optional
        Whether to keep the ``.gz`` files after decompression.
        Only relevant when *mode* is ``"decompress"`` or ``"parquet"``.
        Default is ``True``.
    timeout : int, optional
        HTTP timeout in seconds.  Default is 120.
    max_try : int, optional
        Maximum download retry attempts.  Default is 4.

    Returns
    -------
    results : dict
        A dictionary summarising the operation.  Keys depend on *mode*:

        * ``"download"`` → ``{"zip_paths": [...]}``.
        * ``"decompress"`` → ``{"zip_paths": [...], "trop_paths": [...]}``.
        * ``"parquet"`` → ``{"zip_paths": [...], "trop_paths": [...],
          "parquet_paths": [...]}``.

    Examples
    --------
    Download-only for one station over two years:

    >>> download_ngl_trop("HOUZ", 2000, 2001, "/tmp/ngl_trop")

    Download, extract, and convert to parquet for multiple stations:

    >>> download_ngl_trop(
    ...     ["HOUZ", "REYK"],
    ...     year_start=2020,
    ...     year_end=2022,
    ...     output_dir="/data/ngl_trop",
    ...     mode="parquet",
    ... )
    """
    # --- Normalise inputs ---------------------------------------------------
    if isinstance(stations, str):
        stations = [stations]
    stations = [s.upper() for s in stations]

    valid_modes = ("download", "decompress", "parquet")
    if mode not in valid_modes:
        raise ValueError(f"mode must be one of {valid_modes}, got '{mode}'.")

    years = range(int(year_start), int(year_end) + 1)

    # --- Prepare output directories ----------------------------------------
    dir_zips = os.path.join(output_dir, "zips")
    dir_trop = os.path.join(output_dir, "trop")
    dir_decomp = os.path.join(output_dir, "decompressed")
    dir_parquet = os.path.join(output_dir, "parquet")

    for d in (dir_zips, dir_trop, dir_decomp, dir_parquet):
        os.makedirs(d, exist_ok=True)

    # --- Accumulate results -------------------------------------------------
    zip_paths = []
    trop_paths = []
    parquet_paths = []

    for station in stations:
        for year in years:
            # ------------------------------------------------------------------
            # STEP 1 – Download yearly zip
            # ------------------------------------------------------------------
            url = _ngl_trop_url(station, year, frame=frame)
            zip_save_dir = os.path.join(dir_zips, station)
            zip_path, _ = download_http(
                url,
                output_dir=zip_save_dir,
                timeout=timeout,
                max_try=max_try,
            )

            if not os.path.isfile(zip_path):
                log.warning("Could not download %s – skipping year %d.", url, year)
                continue

            zip_paths.append(zip_path)

            if mode == "download":
                continue

            # ------------------------------------------------------------------
            # STEP 2 – Extract zip → .gz files
            # ------------------------------------------------------------------
            extract_dir = os.path.join(dir_trop, station, str(year))
            gz_files = _extract_zip(zip_path, extract_dir)

            if not keep_zip:
                os.remove(zip_path)
                log.debug("Removed zip: %s", zip_path)

            # ------------------------------------------------------------------
            # STEP 3 – Decompress .gz files → plain SINEX .trop files
            # ------------------------------------------------------------------
            decomp_dir = os.path.join(dir_decomp, station, str(year))
            os.makedirs(decomp_dir, exist_ok=True)  # required by unzip_gz_z
            for gz_path in sorted(gz_files):
                if not gz_path.endswith(".gz"):
                    log.debug("Skipping non-gz file: %s", gz_path)
                    continue

                trop_path = unzip_gz_z(
                    gz_path,
                    out_gzip_dir=decomp_dir,
                    remove_inp=not keep_gz,
                )
                trop_paths.append(trop_path)

                if mode == "decompress":
                    continue

                # --------------------------------------------------------------
                # STEP 4 – Read SINEX → DataFrame → parquet
                # --------------------------------------------------------------
                # Derive DOY from filename: <STAT>.<YYYY>.<DOY>.trop
                basename = os.path.basename(trop_path)
                parts = basename.split(".")
                # Expected pattern: STAT.YYYY.DOY.trop
                if len(parts) >= 3:
                    try:
                        doy_str = parts[2]
                        doy = int(doy_str)
                    except (ValueError, IndexError):
                        doy_str = "000"
                        doy = 0
                else:
                    doy_str = "000"
                    doy = 0

                parquet_fname = f"{station}.{year}.{doy_str:>03}.trop.parquet"
                parquet_dir_stat = os.path.join(dir_parquet, station, str(year))
                os.makedirs(parquet_dir_stat, exist_ok=True)
                parquet_path = os.path.join(parquet_dir_stat, parquet_fname)

                if os.path.isfile(parquet_path):
                    log.debug("Parquet already exists: %s", parquet_path)
                    parquet_paths.append(parquet_path)
                    continue

                try:
                    df = read_snx_trop(trop_path, version=1)
                except Exception as exc:
                    log.warning(
                        "Could not read SINEX file %s: %s", trop_path, exc
                    )
                    continue

                if df is None or df.empty:
                    log.warning("Empty DataFrame for %s – skipping parquet.", trop_path)
                    continue

                df.to_parquet(parquet_path, index=False)
                log.info("Saved parquet: %s", parquet_path)
                parquet_paths.append(parquet_path)

    # --- Build summary dict -------------------------------------------------
    results = {"zip_paths": zip_paths}
    if mode in ("decompress", "parquet"):
        results["trop_paths"] = trop_paths
    if mode == "parquet":
        results["parquet_paths"] = parquet_paths

    return results


