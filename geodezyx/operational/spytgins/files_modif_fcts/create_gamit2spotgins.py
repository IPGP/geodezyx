#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 28/03/2025 18:43:35

@author: psakic
"""

from spytgins import StationFile, OffsetsFile
import os

def create_gamit2spotgins(station_info_path, lfile_path=None, ninecharfile=None, write_amplitude=True, output_dir=None):
        """
        Create station and offsets files from GAMIT station information.

        Parameters
        ----------
        station_info_path : str
            Path to the GAMIT station information file.
        lfile_path : str, optional
            Path to the lfile. Default is None.
        ninecharfile : str, optional
            Path to the nine-character file. Default is None.
        write_amplitude : bool, optional
            Flag to indicate whether to write amplitude information. Default is True.
        output_dir : str, optional
            Directory where the output files will be saved.
            If not provided, it will be the directory of station_info_path.
            Default is None.

        """
        if not output_dir:
            output_dir = os.path.dirname(station_info_path)

        stafil = StationFile.from_gamit(station_info_path, lfile_path=lfile_path, ninecharfile=ninecharfile)
        stafil_outpath = os.path.join(output_dir, 'station_file_from_gamit.dat')
        stafil.write(stafil_outpath)

        offfil = OffsetsFile.from_station_file(stafil)
        offfil_outpath = os.path.join(output_dir, 'offsets_file_from_gamit.dat')
        offfil.write(offfil_outpath, write_amplitude=write_amplitude)

        print(f"Station file created at: {stafil_outpath}")
        print(f"Offsets file created at: {offfil_outpath}")

        return stafil_outpath , offfil_outpath , stafil, offfil


def main():
    parser = argparse.ArgumentParser(
        description="Create station and offsets files from GAMIT station information."
    )
    parser.add_argument('-s', '--station_info_path', required=True,
                        help='Path to the GAMIT station information file.')
    parser.add_argument('-l', '--lfile_path', default=None,
                        help='Path to the lfile. Default is None.')
    parser.add_argument('-n', '--ninecharfile', default=None,
                        help='Path to the nine-character file. Default is None.')
    parser.add_argument('-wa', '--write_amplitude', action='store_true',
                        help='Flag to indicate whether to write amplitude information. Default is True.')
    parser.add_argument('-o', '--output_dir', default=None,
                        help='Directory where the output files will be saved. Default is the directory of station_info_path.')

    args = parser.parse_args()

    stafil_outpath, offfil_outpath, stafil, offfil = create_gamit2spotgins(
        station_info_path=args.station_info_path,
        lfile_path=args.lfile_path,
        ninecharfile=args.ninecharfile,
        write_amplitude=args.write_amplitude,
        output_dir=args.output_dir
    )

if __name__ == "__main__":
    main()