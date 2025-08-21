#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Author: Samuel Nahmani
# Date: 2025-08-14
# Description: Generate SPOTGINS tropospheric time series (.tropo) from listing files.
#              Reads gzip files in reverse order to keep only the last hourly estimation.
#              Produces an easy-to-read, pandas-friendly ISO 8601 date format.
#


import argparse
import gzip
import re
import datetime
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description="Create SPOTGINS tropospheric time series (.tropo) from listing files."
    )
    p.add_argument("-name", required=True, help="Station name (e.g., ABMF00GLP)")
    p.add_argument("-base_dir", required=True,
                   help="Base directory containing <STATION>/020_listings")
    p.add_argument("-output_dir", required=True,
                   help="Directory where the .tropo file will be written")
    p.add_argument("-stations_file",
                   default="/root/spotgins/metadata/stations/station_master_file.dat",
                   help="Path to station_master_file.dat (default: %(default)s)")
    return p.parse_args()


def read_master_file(station_name: str, master_file_path: Path):
    """Retourne un dict avec X,Y,Z,Lon,Lat,H ou None si fichier absent / station non trouvée."""
    if not master_file_path.is_file():
        print(f"Warning: master file not found at {master_file_path}")
        return None
    try:
        with master_file_path.open("r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                # on cherche la station avec un peu de tolérance sur les espaces
                if f" {station_name} " in line or line.strip().startswith(station_name) or station_name in line:
                    fields = re.split(r"\s+", line.strip())
                    if len(fields) >= 9:
                        return {
                            "X": float(fields[6]),
                            "Y": float(fields[7]),
                            "Z": float(fields[8]),
                            "Lon": float(fields[3]),
                            "Lat": float(fields[4]),
                            "H": float(fields[5]),
                        }
    except Exception as e:
        print(f"Error reading master file: {e}")
    return None


def extract_data_from_gz(gz_path: Path):
    """
    Lit le .gz à l'envers, ne garde que la dernière estimation par horodatage.
    Retourne une liste triée: [iso_date, ZTD, ZTD_STD, ZHD, ZWD, TGN, STD_TGN, TGE, STD_TGE].
    """
    records_dict = {}
    try:
        with gzip.open(gz_path, "rt", encoding="utf-8", errors="ignore") as f:
            for line in f.readlines()[::-1]:  # lecture inversée
                if "SINEX_TT2:" not in line:
                    continue
                line = line.split("SINEX_TT2:", 1)[-1].strip()
                fields = re.split(r"\s+", line)
                # format attendu: [STATION, YY:DOY:SEC, ZTD, ZTD_STD, ZHD, ZWD, TGN, STD_TGN, TGE, STD_TGE]
                if len(fields) != 10:
                    continue
                try:
                    y, doy, sec = map(int, fields[1].split(":"))
                    y += 2000 if y < 100 else 0
                    dt = datetime.datetime.strptime(f"{y}-{doy:03d}", "%Y-%j") + datetime.timedelta(seconds=sec)
                    iso_date = dt.isoformat()
                except Exception:
                    continue

                # si on a déjà cet horodatage, on s'arrête pour ce fichier (on a la dernière estimation)
                if iso_date in records_dict:
                    break

                record = [
                    iso_date,
                    float(fields[2]),  # ZTD
                    float(fields[3]),  # ZTD_STD
                    float(fields[4]),  # ZHD
                    float(fields[5]),  # ZWD
                    float(fields[6]),  # TGN
                    float(fields[7]),  # STD_TGN
                    float(fields[8]),  # TGE
                    float(fields[9]),  # STD_TGE
                ]
                records_dict[iso_date] = record
    except Exception as e:
        print(f"Error reading {gz_path.name}: {e}")

    # tri par date ISO
    return list(sorted(records_dict.values(), key=lambda r: r[0]))


def generate_header(station_name: str, position: dict | None):
    """Construit l'entête SPOTGINS. Si position est None, utilise des 0.0."""
    if position is None:
        position = {"X": 0.0, "Y": 0.0, "Z": 0.0, "Lon": 0.0, "Lat": 0.0, "H": 0.0}
        center_note = "SPOTGINS (position unavailable)"
    else:
        center_note = "SPOTGINS - IPGP"

    now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d at %H:%M:%S (UTC)")
    return (
        "# SPOTGINS SOLUTION [TROPO] v1\n"
        f"# Creation {now}\n"
        "#--------------------------------------\n"
        f"# STATION          : {station_name}\n"
        f"# ANALYSIS_CENTRE  : {center_note}\n"
        "# HOW_TO_CITE      : Santamaría-Gómez et al. (2025) https://doi.org/10.5194/essd-2025-223\n"
        "# STRATEGY_SUMMARY : https://en.poleterresolide.fr/geodesy-plotter-en/#/solution/SPOTGINS\n"
        "#--------------------------------------\n"
        "# REF_FRAME        : IGS20\n"
        "# PRODUCTS         : G20/GRG\n"
        "# CONSTELLATION    : GPS only until 2018-10, GPS+Galileo thereafter when available\n"
        "# UNITS            : meters\n"
        "# ELLIPSOID        : GRS80\n"
        "#--------------------------------------\n"
        f"# X_pos            : {position['X']:.6f}\n"
        f"# Y_pos            : {position['Y']:.6f}\n"
        f"# Z_pos            : {position['Z']:.6f}\n"
        f"# Longitude        : {position['Lon']:.6f}\n"
        f"# Latitude         : {position['Lat']:.6f}\n"
        f"# Height           : {position['H']:.6f}\n"
        "#--------------------------------------\n"
        "# Date                 ZTD        ZTD_STD    ZHD        ZWD        TGN        STD_TGN    TGE        STD_TGE\n"
        "# -------------------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------\n"
    )


def main():
    args = parse_args()
    station = args.name
    base_dir = Path(args.base_dir)
    output_dir = Path(args.output_dir)
    master_file = Path(args.stations_file)

    listings_dir = base_dir / station / "020_listings"
    output_file = output_dir / f"SPOTGINS_{station}.tropo"

    print(f"master file SPOTGINS : {master_file}")
    print(output_file)
    print(listings_dir)

    position = read_master_file(station, master_file)
    if not listings_dir.is_dir():
        print(f"Missing directory: {listings_dir}")
        return

    records = []
    for gz_file in sorted(listings_dir.glob("*.gz")):
        records.extend(extract_data_from_gz(gz_file))

    if not records:
        print(f"No valid records found for {station}.")
        return

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", encoding="utf-8") as f:
        f.write(generate_header(station, position))
        for r in records:
            # largeur fixe, supporte les signes négatifs
            f.write(f"{r[0]:<20} " + " ".join(f"{v:>10.6f}" for v in r[1:]) + "\n")


if __name__ == "__main__":
    main()

