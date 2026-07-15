#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Shared I/O safeguards for GNSS educational scripts."""

import datetime as dt
import os
from pathlib import Path

import hatanaka

from geodezyx import files_rw
from geodezyx import operational


def ensure_matplotlib_cache() -> Path:
    """Point matplotlib to a writable cache directory."""
    mplconfigdir = Path("/tmp") / "geodezyx_matplotlib"
    mplconfigdir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mplconfigdir))
    return mplconfigdir


def extract_download_path(entry):
    """Return the local file path from a geodezyx download entry."""
    if isinstance(entry, tuple):
        return entry[0]
    return entry


def expected_rinex_path(station_name: str, date: dt.datetime, work_dir: Path) -> Path:
    """Return the expected local path of a teaching RINEX file."""
    yy = f"{date.year % 100:02d}"
    doy = date.timetuple().tm_yday
    filename = f"{station_name.lower()}{doy:03d}0.{yy}d.Z"
    return work_dir / f"{date.year}" / f"{doy:03d}" / filename


def expected_product_directory(date: dt.datetime, work_dir: Path) -> Path:
    """Return the local directory where daily precise products are stored."""
    doy = date.timetuple().tm_yday
    return work_dir / f"{date.year}" / f"{doy:03d}"


def is_valid_local_rinex(rinex_path: Path) -> bool:
    """Check that a local compressed RINEX can be decompressed."""
    if not rinex_path.exists():
        return False

    try:
        hatanaka.decompress(rinex_path)
    except Exception as exc:
        print(f"Invalid local RINEX detected: {rinex_path}")
        print(f"Reason: {exc}")
        return False

    return True


def is_valid_local_sp3(sp3_path: Path) -> bool:
    """Check that a local SP3 file can be read."""
    if not sp3_path.exists():
        return False

    try:
        files_rw.read_sp3(str(sp3_path), returns_pandas=True, new_col_names=True)
    except Exception as exc:
        print(f"Invalid local SP3 detected: {sp3_path}")
        print(f"Reason: {exc}")
        return False

    return True


def resolve_rinex_files(
    statdico: dict[str, list[str]],
    date: dt.datetime,
    work_dir: Path,
):
    """Reuse valid local RINEX files when available, otherwise download them."""
    expected_paths = []
    missing_paths = []
    invalid_paths = []

    for station_names in statdico.values():
        for station_name in station_names:
            rinex_path = expected_rinex_path(station_name, date, work_dir)
            expected_paths.append(rinex_path)
            if not rinex_path.exists():
                missing_paths.append(rinex_path)
            elif not is_valid_local_rinex(rinex_path):
                invalid_paths.append(rinex_path)

    if not missing_paths and not invalid_paths:
        print("Using local RINEX files already present in the working directory.")
        return [str(path) for path in expected_paths]

    if missing_paths:
        print("Some RINEX files are missing locally, starting download...")
    if invalid_paths:
        print("Some local RINEX files are corrupted, forcing a fresh download...")

    try:
        download_output = operational.download_gnss_rinex(
            statdico=statdico,
            output_dir=str(work_dir),
            startdate=date,
            enddate=date,
            parallel_download=1,
            force=bool(invalid_paths),
        )
    except Exception as exc:
        raise RuntimeError(
            "Unable to obtain valid teaching RINEX files. "
            "Please check your network access or delete the corrupted local copies in "
            f"{work_dir / f'{date.year}' / f'{date.timetuple().tm_yday:03d}'}."
        ) from exc

    print("Download output:")
    print(download_output)
    return download_output


def resolve_sp3_product(
    work_dir: Path,
    processing_date: dt.datetime,
):
    """Reuse a valid local SP3 file when available, otherwise download it."""
    product_dir = expected_product_directory(processing_date, work_dir)
    local_sp3_candidates = sorted(product_dir.glob("*.sp3*"))
    valid_local_sp3 = [path for path in local_sp3_candidates if is_valid_local_sp3(path)]

    if valid_local_sp3:
        print("Using local SP3 product already present in the working directory.")
        return str(valid_local_sp3[0])

    if local_sp3_candidates:
        print("Some local SP3 products are corrupted, forcing a fresh download...")
    else:
        print("No local SP3 product found, starting download...")

    try:
        download_products = operational.download_gnss_products(
            archive_dir=str(work_dir),
            startdate=processing_date,
            enddate=processing_date,
            archtype="year/doy",
            ac_names=("IGS",),
            prod_types=("sp3",),
            repro=0,
            archive_center="ign",
            parallel_download=1,
        )
    except Exception as exc:
        raise RuntimeError(
            "Unable to obtain a valid SP3 product. "
            "Please check your network access or delete the corrupted local SP3 files in "
            f"{product_dir}."
        ) from exc

    print("Precise product download output:")
    print(download_products)
    return extract_download_path(download_products[0])
