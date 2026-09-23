########## BEGIN IMPORT ##########
#### Import the logger
import logging
import os
from os import PathLike
import shutil
from datetime import datetime, timedelta


import pandas as pd
import pyarrow as pa
import  pyarrow.parquet as pq
from tqdm import tqdm

# from threading import Lock

#### geodeZYX modules
from geodezyx import files_rw
from geodezyx import utils
from geodezyx import conv


log = logging.getLogger("geodezyx")


##########  END IMPORT  ##########

def rtklib_out2prq(resdir, pattern="*out", force=False, sample=None):
    """
    Convert RTKLIB output files to Parquet format.

    Parameters
    ----------
    resdir : str
        Results directory containing RTKLIB output files.
    pattern : str, optional
        File pattern to search for (default: "*out").
        gzipped files will be considered automatically.
    force : bool, optional
        Force conversion even if parquet file already exists (default: False).
    sample : str, optional
        Resampling interval for position data (default: None, no resampling).
        If provided, uses _resample_df to resample to the specified interval.
        Examples: "1min", "15min", "1H" (1 hour), "1D" (1 day).

    Returns
    -------
    list
        List of created/updated parquet files.
    """
    l_out0 = utils.find_recursive(resdir, pattern)
    l_out_gz = utils.find_recursive(resdir, pattern  + ".gz")

    if len(l_out_gz) > 0:
        unzip_dir = resdir + "/tmp_unzip"
        utils.create_dir(unzip_dir)
        l_out1 = [utils.uncompress(f, dirout=unzip_dir, opts="") for f in l_out_gz]
    else:
        l_out1 = []

    l_out = l_out0 + l_out1

    f_prq_lis = []
    for f in l_out:
        f_prq = f.replace(".out", ".parquet")
        if not os.path.isfile(f_prq) or force:
            df_out2prq = files_rw.read_rtklib(f, return_df=True)
            if sample:
                df_out2prq = _resample_df(df_out2prq, sample)
            df_out2prq.to_parquet(f_prq, engine="auto", compression="zstd")
            f_prq_lis.append(f_prq)
            log.info(f"Created parquet file: {f_prq}")

    if len(l_out_gz) > 0:
        shutil.rmtree(unzip_dir)

    return f_prq_lis

def _get_year_doy_dirs(base_dir, start_date, end_date):
    """
    Generate list of year/doy subdirectory paths within a date range.

    Parameters
    ----------
    base_dir : str
        Base directory path.
    start_date : datetime
        Start date (inclusive) for generating paths.
    end_date : datetime
        End date (inclusive) for generating paths.

    Returns
    -------
    list of str
        List of paths like <base_dir>/<year>/<doy> within the date range.
    """
    year_doy_dirs = []
    current_date = start_date.date()
    end_date_only = end_date.date()

    while current_date <= end_date_only:
        year = current_date.year
        doy = current_date.timetuple().tm_yday
        year_doy_path = os.path.join(base_dir, str(year), str(doy))
        year_doy_dirs.append(year_doy_path)
        current_date += timedelta(days=1)

    log.info(f"Date range: {start_date.date()} to {end_date.date()} -> {len(year_doy_dirs)} directories")
    return year_doy_dirs

def _drop_pandas_meta(tbl):
    """Drop the 'pandas' metadata key so all tables share the same schema."""
    meta = {k: v for k, v in tbl.schema.metadata.items() if k != b"pandas"}
    return tbl.replace_schema_metadata(meta)

def rtklib_merge_prq(
    parquet_inp,
    exp_prefix="",
    fast_merge=False,
    rtklib_out_files=None,
    sample=None,
    start_date=None,
    end_date=None,
    days=None,
    output_dir=None,
):
    """
    Merge individual RTKLIB parquet files into a single consolidated parquet file.

    Parameters
    ----------
    parquet_inp : str or os.PathLike or list of str
        Either a directory path (all ``*.parquet`` files inside are collected
        recursively) **or** an explicit list of parquet file paths.
        When date filtering is used (start_date, end_date), parquet_inp should be
        a directory structured as <directory>/<year>/<doy>/*.parquet
        The merged output file is written to the directory (or, for a list,
        to the directory of the first file in the list).
    exp_prefix : str, default=""
        Prefix used to name the merged output file (<exp_prefix>_all.parquet).
    fast_merge : bool, default=False
        If True, only merges the parquet files corresponding to
        ``rtklib_out_files`` (or those in the explicit list) and appends them
        to an already-existing ``<exp_prefix>_all.parquet`` file.
        If False, scans the whole directory recursively for parquet files.
    rtklib_out_files : list of str, optional
        List of ``.out`` file paths produced by a previous RTKLIB run.
        Only used when ``fast_merge=True`` and ``parquet_inp`` is a directory,
        to avoid a full recursive scan.
    sample : str, optional
        Resampling interval for position data (default: None, no resampling).
        If provided, uses _resample_df to resample each table to the specified interval
        before merging.
        Warning: resampling slows down the process.
        Examples: "1min", "15min", "1H" (1 hour), "1D" (1 day).
    start_date : str or datetime, optional
        Start date for filtering parquet files (inclusive).
        Can be a string (any format accepted by conv.date_pattern2dt) or datetime object.
        If provided with days or end_date, filters files based on directory structure <year>/<doy>.
    end_date : str or datetime, optional
        End date for filtering parquet files (inclusive).
        Can be a string (any format accepted by conv.date_pattern2dt) or datetime object.
        If provided with days or start_date, filters files based on directory structure <year>/<doy>.
    days : int, optional
        Number of days to process. Only used if start_date XOR end_date is provided.
        If start_date is given, processes N days starting from start_date.
        If end_date is given, processes N days ending at end_date.
    output_dir : str or os.PathLike, optional
        Output directory where the merged parquet file will be saved.
        If not provided, the output is saved to the input directory (or the directory
        of the first file if an explicit list is provided).

    Returns
    -------
    str
        Path to the merged parquet file.
    """
    # --- Parse and validate date parameters ---
    filter_by_date = start_date is not None or end_date is not None

    if filter_by_date:
        # Convert string dates to datetime if needed
        if isinstance(start_date, str):
            start_date = conv.date_pattern2dt(start_date)
        if isinstance(end_date, str):
            end_date = conv.date_pattern2dt(end_date)

        # Compute the actual start and end dates based on days parameter
        if start_date is not None and end_date is None and days is not None:
            end_date = start_date + timedelta(days=days - 1)
        elif end_date is not None and start_date is None and days is not None:
            start_date = end_date - timedelta(days=days - 1)
        elif start_date is None or end_date is None:
            # Need both start_date and end_date for filtering
            raise ValueError(
                "When using date filtering, both start_date and end_date must be provided "
                "(or one of them with days parameter)"
            )

        log.info(f"Date range for filtering: {start_date.date()} to {end_date.date()}")

    # --- resolve source files and output directory ---
    if isinstance(parquet_inp, (str, os.PathLike)) and os.path.isdir(parquet_inp):
        # Use provided output_dir if available, otherwise use input directory
        prq_out_dir = str(output_dir) if output_dir else str(parquet_inp)
        if fast_merge and rtklib_out_files:
            l_prq = [f.replace(".out", ".parquet") for f in rtklib_out_files]
            l_prq = [f for f in l_prq if os.path.exists(f)]
        else:
            # If date filtering is enabled, scan only year/doy directories within the range
            if filter_by_date:
                year_doy_dirs = _get_year_doy_dirs(str(parquet_inp), start_date, end_date)
                l_prq = []
                for year_doy_dir in year_doy_dirs:
                    if os.path.isdir(year_doy_dir):
                        l_prq.extend(utils.find_recursive(year_doy_dir, "*parquet"))
            else:
                l_prq = utils.find_recursive(str(parquet_inp), "*parquet")
    else:
        # parquet_inp is an explicit list of parquet files
        l_prq = list(parquet_inp)
        # Use provided output_dir if available, otherwise use directory of first file
        prq_out_dir = str(output_dir) if output_dir else (os.path.dirname(os.path.abspath(l_prq[0])) if l_prq else ".")

    prq_path_out = os.path.join(prq_out_dir, exp_prefix + "_all.parquet")
    prq_path_tmp = prq_path_out + ".tmp"

    # Exclude the output file itself and stray temp files from the source list
    l_prq = [
        f for f in l_prq if not f.endswith("_all.parquet") and not f.endswith(".tmp")
    ]

    # When fast-merging, prepend the existing merged file so it is streamed
    # first; write to a temp path to avoid reading and writing the same file.
    if fast_merge and os.path.exists(prq_path_out):
        l_prq_merge = [prq_path_out] + l_prq
        prq_path_wrk = prq_path_tmp
    else:
        l_prq_merge = l_prq
        prq_path_wrk = prq_path_out

    # Stream each source table directly through a ParquetWriter —
    # no pandas conversion, no in-memory concat.
    writer = None
    try:
        pbar = tqdm(l_prq_merge, desc="Merging parquet", unit="file")
        for f in pbar:
            pbar.set_postfix_str(os.path.basename(f), refresh=False)
            if sample:
                # Read as pandas, resample, then convert back to pyarrow
                try:
                    df = pd.read_parquet(f)
                    df = _resample_df(df, sample)
                    tbl = pa.Table.from_pandas(df)
                    tbl = _drop_pandas_meta(tbl)
                    corrupt_tbl = False
                except:
                    corrupt_tbl = True
            else:
                tbl = _drop_pandas_meta(pa.parquet.read_table(f))
                corrupt_tbl = True if tbl.num_columns == 0 else False

            if corrupt_tbl:
                log.warning(f"Skipping empty/corrupt parquet file: {f}")
                continue
            if writer is None:
                writer = pa.parquet.ParquetWriter(prq_path_wrk, tbl.schema, compression="zstd")
            writer.write_table(tbl)
    finally:
        if writer:
            writer.close()

    # Atomically replace the previous merged file when using a temp path
    if prq_path_wrk == prq_path_tmp and os.path.exists(prq_path_tmp):
        os.replace(prq_path_tmp, prq_path_out)

    log.info(f"Merged parquet saved to {prq_path_out}")
    return prq_path_out


def _resample_df(df_inp: pd.DataFrame, sample: str = "15min"):
    """
    Resample a DataFrame with GNSS position data to a specified time interval.

    This helper function resamples DataFrame containing GNSS solution positions (x, y, z)
    to a coarser time resolution using median aggregation. It removes any duplicate
    entries resulting from the resampling operation.

    Parameters
    ----------
    df_inp : pandas.DataFrame
        Input DataFrame with an 'epoch' column (datetime or datetime-like) and
        position columns 'x', 'y', 'z' (numeric).
    sample : str, default="15min"
        Resampling interval as a pandas time offset string.
        Examples: "1min", "15min", "1H" (1 hour), "1D" (1 day).

    Returns
    -------
    pandas.DataFrame
        Resampled DataFrame with:
        - 'epoch' column (reset from index)
        - 'x', 'y', 'z' columns containing median values over each resampling interval
        - No duplicate rows

    Notes
    -----
    - Uses median aggregation to provide robust resampling (resistant to outliers and NaN values)
    - Expects the input DataFrame to have an 'epoch' column with datetime values
    - The original epoch index is reset in the output
    """
    # ...existing code...
    df_epo = df_inp.set_index("epoch")

    # Resample position columns with median aggregation
    df_out = df_epo[["x", "y", "z"]].resample(sample).median()

    # Resample quality and uncertainty columns if they exist
    for col in ["Q", "ns", "sdx", "sdy", "sdz", "sdxy", "sdyz", "sdxz", "age", "ratio"]:
        if col in df_epo.columns:
            df_out[col] = df_epo[[col]].resample(sample).median()
    for col in ["rover", "base"]:
        if col in df_epo.columns:
            df_out[col] = df_epo[[col]].resample(sample).first()

    # Reset index to convert epoch back to a regular column
    df_out = df_out.reset_index(inplace=False)
    df_out = df_out.drop_duplicates(inplace=False)
    return df_out


def parquet2csv(
    prq_inp: str | PathLike, out_dir: str | PathLike, sample: str = "15min"
):
    """
    Convert merged RTKLIB parquet file to CSV format, processed by rover/base pairs.

    Reads a merged parquet file containing GNSS solutions from multiple rover/base
    station pairs, resamples each pair's data independently, and exports to CSV files.
    This function uses PyArrow filters to minimize memory usage by reading only the
    data relevant to each rover/base pair.

    Parameters
    ----------
    prq_inp : str or os.PathLike
        Path to the merged parquet file containing rover/base pair GNSS solutions.
        The file must contain columns: 'epoch', 'rover', 'base', 'x', 'y', 'z'.
    out_dir : str or os.PathLike
        Output directory where CSV files will be saved. Created if it doesn't exist.
    sample : str, default="15min"
        Resampling interval for position data. Passed to _resample_df().
        Examples: "1min", "15min", "1H", "1D"

    Returns
    -------
    None

    Output Files
    ------------
    CSV files in `out_dir` with naming pattern:
        {rover}_{base}_{sample}.csv

    Each CSV contains columns:
        - epoch: datetime of the resampled position
        - x, y, z: median position coordinates for the resampling interval

    Notes
    -----
    - Uses PyArrow filter expressions to read only necessary data from the parquet file
    - Efficiently handles large parquet files by filtering at read time (minimal RAM usage)
    - Processes each rover/base pair sequentially
    - Prints rover/base pair names to console as they are processed
    """
    # Read the merged parquet file to get unique rover/base combinations
    log.info(f"Identify rover/base pairs in: {prq_inp}")
    df_rovbas = pd.read_parquet(
        prq_inp, engine="auto", columns=["rover", "base"]
    ).drop_duplicates()

    # Process each unique rover/base pair
    for irow, (rov, bas) in df_rovbas.iterrows():
        log.info("loading rover/base: %s/%s", rov, bas)

        # Define PyArrow filters to read only this rover/base pair
        # Filters minimize memory usage by selecting data at read time
        filters = [
            ("rover", "==", rov),
            ("base", "==", bas),
        ]

        # Read only the filtered data for this pair
        df_grp = pd.read_parquet(prq_inp, engine="auto", filters=filters)

        # Resample to the specified time interval
        df_out = _resample_df(df_grp, sample)

        df_out["year"] = df_out["epoch"].dt.year
        df_out["doy"] = df_out["epoch"].apply(conv.dt2doy, args=(int,))
        df_out["mjd"] = df_out["epoch"].apply(conv.dt2mjd)

        # Export to CSV file with naming convention: rover_base_sample.csv
        out_path = f"{out_dir}/{rov}_{bas}_{sample}.csv"
        log.info("saving %s resampled data to CSV: %s", sample, out_path)
        df_out.to_csv(out_path, index=False)

    return None


def rtklib_prq2out(
    prq_inp: str | PathLike, out_dir: str | PathLike, sample: str | None = None, rover: str | None = None, base: str | None = None
):
    """
    Convert merged RTKLIB parquet file to .out format, processed by rover/base pairs.

    Reads a merged parquet file containing GNSS solutions from multiple rover/base
    station pairs and exports to RTKLIB .out format. This function uses PyArrow filters
    to minimize memory usage by reading only the data relevant to each rover/base pair.
    Optionally resamples data to a specified time interval.

    Parameters
    ----------
    prq_inp : str or os.PathLike
        Path to the merged parquet file containing rover/base pair GNSS solutions.
        The file must contain columns: 'epoch', 'rover', 'base', 'x', 'y', 'z',
        'Q', 'ns', 'sdx', 'sdy', 'sdz', 'sdxy', 'sdyz', 'sdxz', 'age', 'ratio'.
    out_dir : str or os.PathLike
        Output directory where .out files will be saved. Created if it doesn't exist.
    sample : str, default=None
        Optional resampling interval for position data. If provided, uses _resample_df()
        to resample to the specified interval before writing.
        Examples: "1min", "15min", "1H" (1 hour), "1D" (1 day).
    rover : str, default=None
        Optional filter to process only a specific rover site. If None, processes all rovers.
    base : str, default=None
        Optional filter to process only a specific base site. If None, processes all bases.

    Returns
    -------
    list
        List of paths to generated .out files.

    Output Files
    ------------
    .out files in `out_dir` with naming pattern:
        {rover}_{base}.out  (if no sample)
        {rover}_{base}_{sample}.out  (if sample is specified)

    Each .out file contains:
        - RTKLIB header with metadata
        - Column headers: GPST, x-ecef(m), y-ecef(m), z-ecef(m), Q, ns, sdx(m), sdy(m), sdz(m), sdxy(m), sdyz(m), sdzx(m), age(s), ratio
        - Data rows with position and quality information

    Notes
    -----
    - Uses PyArrow filter expressions to read only necessary data from the parquet file
    - Efficiently handles large parquet files by filtering at read time (minimal RAM usage)
    - Processes each rover/base pair sequentially
    - RTKLIB .out format uses XYZ ECEF coordinates and GPST time scale
    - Quality codes: 1=fix, 2=float, 3=sbas, 4=dgps, 5=single, 6=ppp
    """
    # Create output directory if it doesn't exist
    out_dir = utils.create_dir(str(out_dir))

    # Read the merged parquet file to get unique rover/base combinations
    log.info(f"Identify rover/base pairs in: {prq_inp}")
    df_rovbas = pd.read_parquet(
        prq_inp, engine="auto", columns=["rover", "base"]
    ).drop_duplicates()

    # Apply optional filters
    if rover is not None:
        df_rovbas = df_rovbas[df_rovbas["rover"] == rover]
    if base is not None:
        df_rovbas = df_rovbas[df_rovbas["base"] == base]

    out_files = []

    # Process each unique rover/base pair
    for irow, (rov, bas) in df_rovbas.iterrows():
        log.info("Converting rover/base: %s/%s to .out format", rov, bas)

        # Define PyArrow filters to read only this rover/base pair
        filters = [
            ("rover", "==", rov),
            ("base", "==", bas),
        ]

        # Read only the filtered data for this pair
        df_grp = pd.read_parquet(prq_inp, engine="auto", filters=filters)

        # Resample to the specified time interval if requested
        if sample:
            df_grp = _resample_df(df_grp, sample)
            out_filename = f"{rov}_{bas}_{sample}.out"
        else:
            out_filename = f"{rov}_{bas}.out"

        out_path = os.path.join(out_dir, out_filename)

        # Write RTKLIB .out format
        _write_rtklib_out(df_grp, out_path, rover=rov, base=bas)
        out_files.append(out_path)
        log.info("Saved .out file: %s", out_path)

    return out_files


def _write_rtklib_out(df: pd.DataFrame, out_path: str, rover: str = "XXXX00XXX", base: str = "XXXX00XXX"):
    """
    Write a DataFrame to RTKLIB .out format file.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing GNSS solution data with columns:
        epoch, x, y, z, Q, ns, sdx, sdy, sdz, sdxy, sdyz, sdxz, age, ratio
    out_path : str
        Output file path for the .out file
    rover : str, default="XXXX00XXX"
        Rover site name for header metadata
    base : str, default="XXXX00XXX"
        Base site name for header metadata
    """
    if df.empty:
        log.warning(f"DataFrame is empty, skipping: {out_path}")
        return

    # Sort by epoch to ensure chronological order
    df = df.sort_values("epoch").reset_index(drop=True)

    # Extract time range from data
    epoch_start = df["epoch"].min()
    epoch_end = df["epoch"].max()

    # Write file
    with open(out_path, "w") as f:
        # Write RTKLIB header
        f.write("% program   : geodezyx parquet2rtklib_out converter\n")
        f.write("% (converted from parquet format)\n")
        f.write(f"% rover site: {rover}\n")
        f.write(f"% base site : {base}\n")
        f.write(f"% obs start : {epoch_start.strftime('%Y/%m/%d %H:%M:%S.%f')[:-3]} GPST\n")
        f.write(f"% obs end   : {epoch_end.strftime('%Y/%m/%d %H:%M:%S.%f')[:-3]} GPST\n")
        f.write("% pos mode  : (unknown)\n")
        f.write("% freqs     : (unknown)\n")
        f.write("% ephemeris : (unknown)\n")
        f.write("% solution  : (unknown)\n")
        f.write("%\n")
        f.write("% (x/y/z-ecef=WGS84,Q=1:fix,2:float,3:sbas,4:dgps,5:single,6:ppp,ns=# of satellites)\n")
        f.write("%  GPST                      x-ecef(m)      y-ecef(m)      z-ecef(m)   Q  ns   sdx(m)   sdy(m)   sdz(m)  sdxy(m)  sdyz(m)  sdzx(m) age(s)  ratio\n")

        # Write data rows
        for idx, row in df.iterrows():
            epoch = row["epoch"]
            date_str = epoch.strftime("%Y/%m/%d %H:%M:%S.%f")[:-3]

            # Extract values with safe defaults
            x = float(row["x"])
            y = float(row["y"])
            z = float(row["z"])
            Q = int(row["Q"]) if "Q" in row and pd.notna(row["Q"]) else 0
            ns = int(row["ns"]) if "ns" in row and pd.notna(row["ns"]) else 0
            sdx = float(row["sdx"]) if "sdx" in row and pd.notna(row["sdx"]) else 0.0
            sdy = float(row["sdy"]) if "sdy" in row and pd.notna(row["sdy"]) else 0.0
            sdz = float(row["sdz"]) if "sdz" in row and pd.notna(row["sdz"]) else 0.0
            sdxy = float(row["sdxy"]) if "sdxy" in row and pd.notna(row["sdxy"]) else 0.0
            sdyz = float(row["sdyz"]) if "sdyz" in row and pd.notna(row["sdyz"]) else 0.0
            sdxz = float(row["sdxz"]) if "sdxz" in row and pd.notna(row["sdxz"]) else 0.0
            age = float(row["age"]) if "age" in row and pd.notna(row["age"]) else 0.0
            ratio = float(row["ratio"]) if "ratio" in row and pd.notna(row["ratio"]) else 0.0

            # Format output line (matching RTKLIB .out format exactly)
            # Format: date time, then %14.4f x, %14.4f y, %14.4f z, %3d Q, %3d ns, %8.4f sd*, %6.2f age, %6.1f ratio
            line = f"{date_str} {x:14.4f} {y:14.4f} {z:14.4f} {Q:3d} {ns:3d} {sdx:8.4f} {sdy:8.4f} {sdz:8.4f} {sdxy:8.4f} {sdyz:8.4f} {sdxz:8.4f} {age:6.2f} {ratio:6.1f}\n"
            f.write(line)
