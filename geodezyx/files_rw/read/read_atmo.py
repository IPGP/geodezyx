#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 10:52:09 2021

@author: psakic
"""

import re
import os

########## BEGIN IMPORT ##########
#### External modules
import numpy as np
import pandas as pd

#### geodeZYX modules
from geodezyx import conv
from geodezyx import utils

##########  END IMPORT  ##########


############### Reading Tropospheric ################################################

#  _______                              _                     ______ _ _
# |__   __|                            | |                   |  ____(_) |
#    | |_ __ ___  _ __   ___  ___ _ __ | |__   ___ _ __ ___  | |__   _| | ___  ___
#    | | '__/ _ \| '_ \ / _ \/ __| '_ \| '_ \ / _ \ '__/ _ \ |  __| | | |/ _ \/ __|
#    | | | | (_) | |_) | (_) \__ \ |_) | | | |  __/ | |  __/ | |    | | |  __/\__ \
#    |_|_|  \___/| .__/ \___/|___/ .__/|_| |_|\___|_|  \___| |_|    |_|_|\___||___/
#                | |             | |
#                |_|             |_|


def _parse_trop_description(snxfile):
    """
    Parse the ``TROP/DESCRIPTION`` block of a SINEX troposphere file to extract
    the ordered list of solution column names.

    Handles both format generations:

    * **SINEX_TRO v0.01** – keywords ``SOLUTION_FIELDS_1`` and (optionally)
      ``SOLUTION_FIELDS_2`` list the column names.  ``SOLUTION_FIELDS_2`` is a
      continuation line used when more than 7 fields are needed.
    * **SINEX_TRO v2.00** – keyword ``TROPO PARAMETER NAMES`` (three tokens)
      lists the column names.

    Post-processing applied to the raw column list:

    * Every ``STDDEV`` token is renamed to ``{previous_column}_sig``, following
      the SINEX_TRO convention that ``STDDEV`` is the standard deviation of the
      immediately preceding parameter.
    * All names are lower-cased.
    * Leading ``#`` characters (e.g. ``#ACTAK`` → ``actak``) are stripped.

    Parameters
    ----------
    snxfile : str
        Path to the SINEX troposphere file.

    Returns
    -------
    cols : list of str
        Ordered list of lowercase column names found in ``TROP/DESCRIPTION``,
        or an empty list if the block or the relevant keyword is absent.
    """
    cols_v2 = []   # SINEX_TRO v2.00  – TROPO PARAMETER NAMES
    cols_v1 = []   # SINEX_TRO v0.01  – SOLUTION_FIELDS_1 / _2
    in_block = False

    for line in open(snxfile, "r", encoding="ISO-8859-1"):
        if re.search(r"\+TROP/DESCRIPTION", line):
            in_block = True
            continue
        if re.search(r"-TROP/DESCRIPTION", line):
            break
        if not in_block:
            continue
        # Only process data lines (first character is a space)
        if not line.startswith(" "):
            continue
        parts = line.split()
        if not parts:
            continue

        # v2.00: "TROPO PARAMETER NAMES   col1 col2 ..."
        if len(parts) >= 3 and " ".join(parts[:3]).upper() == "TROPO PARAMETER NAMES":
            cols_v2.extend(parts[3:])

        # v0.01: "SOLUTION_FIELDS_1   col1 col2 ..." (SOLUTION_FIELDS_2 is continuation)
        elif parts[0].upper().startswith("SOLUTION_FIELDS"):
            cols_v1.extend(parts[1:])

    # Prefer v2.00 format when available
    raw_cols = cols_v2 if cols_v2 else cols_v1

    if not raw_cols:
        return []

    # Post-process: rename STDDEV → {prev_col}_sig, lowercase, strip leading '#'
    result = []
    for col in raw_cols:
        col_clean = col.lstrip("#").lower()
        if col_clean == "stddev" and result:
            result.append(result[-1] + "_sig")
        else:
            result.append(col_clean)

    return result


# ---------------------------------------------------------------------------
# Fallback column names used when TROP/DESCRIPTION is absent.
# Keys are the number of *value* fields (i.e. total fields minus STAT and EPOCH).
# ---------------------------------------------------------------------------
_TROP_FALLBACK_COLS = {
    # ZTD-only (e.g. simple AC submissions)
    2: ["trotot", "trotot_sig"],
    # Standard IGS 8-field layout: ZTD + grad North + grad East
    6: ["trotot", "trotot_sig", "tgntot", "tgntot_sig", "tgetot", "tgetot_sig"],
    # NGL extended 12-field layout (STAT EPOCH + 10 values)
    10: [
        "trotot", "trotot_sig",
        "trowet",
        "tgetot", "tgetot_sig",
        "tgntot", "tgntot_sig",
        "wvapor", "wvapor_sig",
        "mtemp",
    ],
}


def read_snx_trop(snxfile, dataframe_output=True, auto_desc_cols=False):
    """
    Read troposphere solutions from a Troposphere SINEX file.

    Dynamically parses the ``TROP/DESCRIPTION`` block to discover the column
    layout before reading the ``TROP/SOLUTION`` block.  Both
    **SINEX_TRO v0.01** (``SOLUTION_FIELDS_1``/``_2`` keywords, 2-digit years)
    and **SINEX_TRO v2.00** (``TROPO PARAMETER NAMES`` keyword, 4-digit years,
    long station names up to 9 characters) are supported.

    Column naming convention (v2.00 Table 1-3):

    * ``trotot`` / ``trotot_sig`` – zenith total delay and its std dev
    * ``trowet`` / ``trowet_sig`` – zenith wet delay
    * ``trodry`` / ``trodry_sig`` – zenith dry (hydrostatic) delay
    * ``tgntot`` / ``tgntot_sig`` – total north gradient (wet + dry)
    * ``tgetot`` / ``tgetot_sig`` – total east gradient (wet + dry)
    * ``tgnwet``, ``tgnwet_sig``, ``tgedry``, ``tgedry_sig`` –
      wet/dry components of the horizontal gradients
    * ``iwv``    – integrated water vapour
    * ``press``  – surface pressure
    * ``epress`` – partial water vapour pressure
    * ``temdry`` – dry temperature
    * ``humrel`` – relative humidity
    * ``nsat``   – number of satellites used
    * ``gdop``   – geometric dilution of precision
    * ``wmtemp`` / ``wmtemp_sig`` – weighted mean temperature
    * ``acok``, ``acdl`` – combination quality counters (combined products)
    * ``dstax``, ``dstay``, ``dstaz`` – station coordinate differences
      (combined products, v0.01)
    * ``pwv`` / ``pwv_sig`` – precipitable water vapour (v0.01)

    Every ``STDDEV`` token in ``TROP/DESCRIPTION`` is automatically renamed to
    ``{preceding_parameter}_sig``.

    When no ``TROP/DESCRIPTION`` block is found the function falls back to
    heuristics based on the number of value fields per data line (2, 6, or 10
    values; see ``_TROP_FALLBACK_COLS``).  For any other count, generic names
    ``col1``, ``col2``, … are assigned.

    The epoch year is auto-detected from the string length of the year token:

    * 4-character token → direct 4-digit year (SINEX_TRO v2.00 / modern files).
    * 2-character token → standard SINEX 2-digit convention:
      ``YY`` 80–99 → 1980–1999, ``YY`` 00–79 → 2000–2079.

    Parameters
    ----------
    snxfile : str
        Path to the SINEX troposphere solution file (plain text or
        ISO-8859-1 encoded).
    dataframe_output : bool, optional
        If ``True`` (default), return a :class:`pandas.DataFrame`.
        If ``False``, return a sorted list of dicts (one dict per epoch/station
        record) with keys ``STAT``, ``epoc``, and one key per data column.
    auto_desc_cols : bool, optional
        If ``True``, column names are extracted from the ``TROP/DESCRIPTION`` block and
        used as-is (after post-processing).
        If ``False`` (default), the function falls back to predefined column names.
        Set to ``True`` to preserve original column names from the file;
        set to ``False`` to enforce a consistent naming convention based on field count.

    Returns
    -------
    df : pandas.DataFrame or list of dict
        If ``dataframe_output=True``: DataFrame with columns ``STAT``,
        ``epoc`` (datetime), plus one column per parameter found in
        ``TROP/DESCRIPTION`` (or heuristic fallback names).
        All parameter columns are cast to ``float``.

        If ``dataframe_output=False``: list of dicts sorted by
        ``(STAT, epoc)``.

    See Also
    --------
    troposinex2df : Convert the list-of-dicts output to a DataFrame.
    _parse_trop_description : Parse the TROP/DESCRIPTION block.

    References
    ----------
    * SINEX_TRO v0.01 – https://files.igs.org/pub/data/format/sinex_tropo.txt
    * SINEX_TRO v2.00 – https://files.igs.org/pub/data/format/sinex_tro_v2.00.pdf
    """
    # --- Step 1: discover column names from TROP/DESCRIPTION ---
    desc_cols = _parse_trop_description(snxfile)

    # --- Step 2: parse TROP/SOLUTION ---
    records = []
    in_trop = False

    def _fnan(x):
        """Replace asterisk-flagged values with NaN."""
        return np.nan if isinstance(x, str) and "*" in x else x

    for line in open(snxfile, "r", encoding="ISO-8859-1"):
        if re.compile("TROP/SOLUTION").search(line):
            in_trop = not in_trop
            continue

        if not in_trop:
            continue

        # Skip comment lines, block markers, etc.
        if not line.startswith(" "):
            continue

        fields = line.split()
        if len(fields) < 3:
            continue

        stat = fields[0].upper()
        epoch_str = fields[1]

        # --- Epoch parsing (auto-detect 2- vs 4-digit year) ---
        if ":" not in epoch_str:
            epoch = conv.year_decimal2dt(epoch_str)
        else:
            date_elts = epoch_str.split(":")
            year_str = date_elts[0]
            yy_raw = int(year_str)
            if len(year_str) == 4:
                # SINEX_TRO v2.00 – 4-digit year
                yy = yy_raw
            else:
                # SINEX_TRO v0.01 – 2-digit year, standard SINEX convention
                yy = (1900 + yy_raw) if yy_raw >= 80 else (2000 + yy_raw)
            doy = int(date_elts[1])
            sec = int(date_elts[2])
            epoch = conv.doy2dt(yy, doy, seconds=sec)

        values = fields[2:]
        n = len(values)
        rec = {"site": stat, "epoch": epoch}

        if desc_cols and auto_desc_cols:
            # Dynamic path: use column names from TROP/DESCRIPTION
            for i, col in enumerate(desc_cols):
                rec[col] = _fnan(values[i]) if i < n else np.nan
        else:
            # Fallback path: heuristic based on value-field count
            col_names = _TROP_FALLBACK_COLS.get(
                n, [f"col{i + 1}" for i in range(n)]
            )
            for i, col in enumerate(col_names):
                rec[col] = _fnan(values[i]) if i < n else np.nan

        records.append(rec)

    # Sort by station then epoch
    records.sort(key=lambda r: (r["site"], r["epoch"]))

    if dataframe_output:
        return troposinex2df(records)
    else:
        return records

def troposinex2df(records):
    """
    Convert a list of SINEX troposphere records to a pandas DataFrame.

    Parameters
    ----------
    records : list of dict
        List of per-record dicts as returned by :func:`read_snx_trop` with
        ``dataframe_output=False``.  Each dict must contain at least the keys
        ``STAT`` (str) and ``epoc`` (datetime); remaining keys are
        tropospheric parameter columns.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame with ``STAT`` and ``epoc`` as the first two columns,
        followed by one column per parameter.  All parameter columns are
        cast to ``float``.
    """
    if not records:
        return pd.DataFrame(columns=["site", "epoch"])
    df = pd.DataFrame(records)
    fixed = ["site", "epoch"]
    param_cols = [c for c in df.columns if c not in fixed]
    df = df[fixed + param_cols]
    df[param_cols] = df[param_cols].apply(pd.to_numeric, errors="coerce")
    return df


def read_gfz_trop(trpfile):
    """
    Read GFZ troposphere SINEX solution into pandas DataFrame.

    Parses GFZ SINEX format and extracts station name, epoch, and ZTD with
    associated uncertainties plus horizontal gradients.

    Parameters
    ----------
    trpfile : str
        Path to GFZ troposphere SINEX file.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame with columns: STAT, epoc, year, doy, secofday, ztd_est,
        ztd_est_std, num_sat, tgn_est, tgn_est_std, tge_est, tge_est_std.
        Numeric columns are converted to float type.
    """
    fields = []
    flagtrop = False

    for line in open(trpfile, "r", encoding="ISO-8859-1"):
        if re.compile("TROP/SOLUTION").search(line):
            flagtrop = not flagtrop
            continue

        if flagtrop == True and line[0] == " ":
            field = line.split()
            fields.append(field)
        else:
            continue

    DF = pd.DataFrame(fields)
    DF.drop(columns=[0, 8, 9, 10, 11, 12, 18, 19, 20], inplace=True)
    DF.columns = [
        "STAT",
        "epoc",
        "year",
        "doy",
        "secofday",
        "ztd_est",
        "ztd_est_std",
        "num_sat",
        "tgn_est",
        "tgn_est_std",
        "tge_est",
        "tge_est_std",
    ]
    cols_numeric = [
        "epoc",
        "ztd_est",
        "ztd_est_std",
        "num_sat",
        "tgn_est",
        "tgn_est_std",
        "tge_est",
        "tge_est_std",
    ]
    DF[cols_numeric] = DF[cols_numeric].apply(pd.to_numeric, errors="coerce")
    DF["epoc"] = conv.mjd2dt(DF["epoc"].values)
    DF["epoc"] = DF["epoc"].dt.floor("H")
    return DF



def read_bernese_trp(trpfile):
    """
    Read tropospheric solution in TRP format from Bernese GNSS software.

    Parses Bernese TRP format and returns troposphere parameters including
    ZTD, north-south and east-west gradients with associated uncertainties.

    Parameters
    ----------
    trpfile : str
        Path to TRP file from Bernese GNSS software.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame with columns: STAT, year, month, day, hour, minute, second,
        MOD_U, CORR_U, SIGMA_U, TOTAL_U, CORR_N, SIGMA_N, CORR_E, SIGMA_E,
        and datetime column 'dt'. Numeric columns are converted to float type.

    Notes
    -----
    Written by Chaiyaporn Kitpracha.
    """
    flagtrop = False
    field = []
    for line in open(trpfile, "r", encoding="ISO-8859-1"):
        if re.compile("STATION NAME").search(line):
            headers = line.split()
            headers.remove("YYYY")
            headers.remove("MM")
            headers.remove("DD")
            headers.remove("HH")
            headers.remove("MM")
            headers.remove("SS")
            headers[3] = "year"
            headers[4] = "month"
            headers[5] = "day"
            headers[6] = "hour"
            headers[7] = "minute"
            headers[8] = "second"
            flagtrop = True
            continue

        if flagtrop and not line == "\n":
            fields = line.split()
            field.append(fields)
        else:
            continue

    DF = pd.DataFrame(field, columns=headers)
    DF["dt"] = pd.to_datetime(DF[["year", "month", "day", "hour", "minute", "second"]])
    cols_num = [
        "MOD_U",
        "CORR_U",
        "SIGMA_U",
        "TOTAL_U",
        "CORR_N",
        "SIGMA_N",
        "CORR_E",
        "SIGMA_E",
    ]
    DF[cols_num] = DF[cols_num].apply(pd.to_numeric, errors="coerce")
    DF.drop(["year", "month", "day", "hour", "minute", "second"], axis=1, inplace=True)
    return DF


def read_rinex_met(metfile):
    """
    Read RINEX meteorological files and convert to pandas DataFrame.

    Handles single file or multiple files (as list or iterable). Concatenates
    multiple files into a single DataFrame with proper time indexing.

    Parameters
    ----------
    metfile : str or list of str
        Path(s) to RINEX meteorological file(s). Can be a single filename string
        or a list/iterable of filenames (e.g., from glob).

    Returns
    -------
    df : pandas.DataFrame
        Meteorological data from RINEX file(s). Index is set to epoch (datetime).
        Includes columns for temperature, pressure, humidity, and associated
        uncertainties if available.

    Notes
    -----
    Written by Chaiyaporn Kitpracha.
    """
    if utils.is_iterable(metfile):
        merge_df = pd.DataFrame()
        for metfile_m in metfile:
            met_df = read_rinex_met_2(str(metfile_m))
            merge_df = pd.concat([merge_df, met_df])
        return merge_df
    else:
        met_df = read_rinex_met_2(metfile)
        return met_df


def read_rinex_met_2(metfile):
    """
    Read a single RINEX meteorological file (internal function).

    Worker function for read_rinex_met. Parses header for sensor metadata and
    converts meteorological observations to DataFrame with proper types and
    datetime indexing.

    Parameters
    ----------
    metfile : str
        Path to a single RINEX meteorological file.

    Returns
    -------
    df : pandas.DataFrame
        Meteorological data with datetime index. Includes temperature, pressure,
        humidity observations and sensor uncertainties. Station code (STA column)
        is added from RINEX header.

    Notes
    -----
    Written by Chaiyaporn Kitpracha.
    """
    ln = 0
    for line in open(metfile, "r", encoding="ISO-8859-1"):
        if re.compile("MARKER NAME").search(line):
            marker = line.split()[0]
            marker = marker.upper()
        if re.compile("# / TYPES OF OBSERV").search(line):
            tmp = line.split()
            headers = tmp[1 : int(tmp[0]) + 1]
        if re.compile("TD SENSOR MOD/TYPE/ACC").search(line):
            tmp = line.split()
            temp_unc = float(tmp[-4])
        if re.compile("PR SENSOR MOD/TYPE/ACC").search(line):
            tmp = line.split()
            press_unc = float(tmp[-4])
        if re.compile("HR SENSOR MOD/TYPE/ACC").search(line):
            tmp = line.split()
            humrel_unc = float(tmp[-4])
        if re.compile("END OF HEADER").search(line):
            break
        ln = ln + 1

    df = pd.read_csv(
        metfile,
        skiprows=range(0, ln + 1),
        delim_whitespace=True,
        names=["year", "month", "day", "hour", "minute", "second"] + headers,
    )
    df["year"] = (
        df["year"] + 2000 if df["year"].any() <= 79 else df["year"].any() + 1900
    )
    df["STA"] = marker
    df["epoch"] = pd.to_datetime(
        df[["year", "month", "day", "hour", "minute", "second"]], errors="coerce"
    )
    df.drop(["year", "month", "day", "hour", "minute", "second"], axis=1, inplace=True)
    if press_unc is not None:
        df["PR_std"] = press_unc
    if temp_unc is not None:
        df["TD_std"] = temp_unc
    if humrel_unc is not None:
        df["HR_std"] = humrel_unc
    df.set_index("epoch", inplace=True)
    return df


# ---------------------------------------------------------------------------
# Column mapping from SPOTGINS native names to SINEX_TRO convention
# (as used by read_snx_trop).  Applied when sinex_columns=True.
# ---------------------------------------------------------------------------
_SPOTGINS_TO_SINEX_COLS = {
    # Station / epoch
    "site":            "site",
    "epoch":           "epoch",
    # Auxiliary
    "MJD":             "mjd",
    "DecimalYear":     "decimal_year",
    "Const":           "const",
    "Dateofexe":       "dateofexe",
    "GinsVersion":     "ginsversion",
    "PrairieVersion":  "prairieversion",
    # ZTD parameters
    "TROTOT":          "trotot",
    "TRODRY":          "trodry",
    "TROWET":          "trowet",
    "STDWET":          "trowet_sig",
    "STDTOT":          "trotot_sig",
    "STDDRY":          "trodry_sig",
    # Gradient parameters
    "TGNTOT":          "tgntot",
    "STDTGN":          "tgntot_sig",
    "TGETOT":          "tgetot",
    "STDTGE":          "tgetot_sig",
    "TGNWET":          "tgnwet",
    "STDTGNWET":       "tgnwet_sig",
    "TGEWET":          "tgewet",
    "STDTGEWET":       "tgewet_sig",
}


def read_spotgins_tropo(filepath, sinex_columns=False):
    """
    Read a SPOTGINS tropospheric time series file (ZTD or GRAD format).

    Parses SPOTGINS ZTD (zenith total delay) or GRAD (horizontal gradients)
    format and returns tropospheric parameters with associated metadata.
    Column names are extracted from file headers and kept in their original
    case (e.g. ``TROTOT``, ``MJD``, ``Const``).  The datetime column
    (``yyyy-mm-ddTHH:MM:SS``) is always renamed to ``epoch``.
    Header metadata is dynamically inferred from the file (no hardcoded fields).

    Parameters
    ----------
    filepath : str
        Path to a SPOTGINS ``.ztd`` or ``.grad`` file (can be gzip compressed).
    sinex_columns : bool, optional
        If ``False`` (default), column names keep the original case from the
        file (e.g. ``TROTOT``, ``TRODRY``, ``Const``), with ``site`` and
        ``epoch`` added by the reader.

        If ``True``, column names are renamed to the SINEX_TRO convention
        used by :func:`read_snx_trop` (lowercase):

        * ``site``  → ``site``
        * ``epoch`` → ``epoch``
        * ``MJD``   → ``mjd``
        * ``TROTOT`` / ``TRODRY`` / ``TROWET`` → ``trotot`` / ``trodry`` / ``trowet``
        * ``STDWET`` → ``trowet_sig``
        * ``TGNTOT`` / ``TGETOT`` → ``tgntot`` / ``tgetot``
        * ``STDTGN`` / ``STDTGE`` → ``tgntot_sig`` / ``tgetot_sig``
        * Auxiliary columns (``Const``, ``Dateofexe``, ``GinsVersion``,
          ``PrairieVersion``, ``DecimalYear``) are lowercased.

        The full mapping is available in the module-level dict
        ``_SPOTGINS_TO_SINEX_COLS``.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame with one row per epoch. Columns use original file case
        (or SINEX_TRO names when ``sinex_columns=True``). No index is set.
    meta : dict
        Header metadata extracted from comment block. Keys are converted to
        lowercase with spaces replaced by underscores. Includes station
        information, constellation, coordinates, and other metadata.
    """
    meta = {}
    data_lines = []
    cols = None
    datetime_col = None

    with open(filepath, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                # The last comment line is the column header
                cols = line.lstrip("#").strip().split()
                # Find the datetime column (contains timestamp format)
                datetime_col = next((col for col in cols if any(x in col.lower() for x in ["hh", "mm", "ss"])), None)
                
                # Extract metadata from key: value lines
                m = re.match(r"^#\s*([^:]+?)\s*:\s*(.+)", line)
                if m:
                    key = m.group(1).strip()
                    value = m.group(2).strip()
                    # Convert key to lowercase with underscores for spaces
                    meta_key = key.lower().replace(" ", "_")
                    try:
                        meta[meta_key] = float(value)
                    except ValueError:
                        meta[meta_key] = value
            else:
                stripped = line.strip()
                if stripped:
                    data_lines.append(stripped)

    if cols is None:
        raise ValueError("Could not find column header in file")

    if not data_lines:
        empty_df = pd.DataFrame(columns=cols)
        if sinex_columns:
            empty_df = empty_df.rename(columns=_SPOTGINS_TO_SINEX_COLS)
        return empty_df, meta

    from io import StringIO

    df = pd.read_csv(
        StringIO("\n".join(data_lines)),
        sep=r"\s+",
        header=None,
        names=cols,
    )

    # Rename datetime column to epoch (keep original case for all other columns)
    if datetime_col and datetime_col in df.columns:
        df = df.rename(columns={datetime_col: "epoch"})

    # Convert numeric columns to float (exclude string/metadata columns)
    # Match known string columns case-insensitively to preserve original case
    _known_string_cols = {"epoch", "const", "dateofexe", "ginsversion", "prairieversion"}
    string_cols = {c for c in df.columns if c.lower() in _known_string_cols}
    float_cols = [c for c in df.columns if c not in string_cols]
    df[float_cols] = df[float_cols].apply(pd.to_numeric, errors="coerce")
    df["epoch"] = pd.to_datetime(df["epoch"], errors="coerce")

    # Add station from metadata as first column named 'site'
    if "station" in meta:
        df.insert(0, "site", meta["station"])

    # Optionally rename to SINEX_TRO convention (as used by read_snx_trop)
    if sinex_columns:
        df = df.rename(columns=_SPOTGINS_TO_SINEX_COLS)

    return df, meta
