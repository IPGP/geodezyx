#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 23:31:09 2025

@author: psakic
"""
import numpy as np
import os


def read_antex(filepath, ant_type=None):
    """
    Parse an ANTEX file and extract antenna information.

    Parameters
    ----------
    filepath : str
        Path to the ANTEX file to be read.
    ant_type : str, optional
        Specific antenna type to filter. If provided, only data for this antenna type will be included.

    Returns
    -------
    dict
        A dictionary containing parsed antenna data. The structure includes:
        - HEADER: List of header lines (if `header_in_output` is True).
        - Antenna type as keys, with corresponding metadata and frequency-specific data as values.

    Note
    ----
    Structure of the returned dictionary:

    dict
    ├── HEADER (list of str)
    ├── ANTS (dict)
    │   ├── <Antenna Type> (str)
    │   │   ├── TYPE (str)
    │   │   ├── SERIAL (str)
    │   │   ├── METH (str)
    │   │   ├── DAZI (float)
    │   │   ├── ZEN (list of float)
    │   │   ├── NFREQ (int)
    │   │   ├── SINEX_CODE (str)
    │   │   ├── COMMENT (list of str)
    │   │   ├── FREQS (dict)
    │   │   │   ├── <Frequency> (str)
    │   │   │   │   ├── PCO (list of float)
    │   │   │   │   ├── NOAZI (numpy array)
    │   │   │   │   ├── AZI (numpy array)
    """
    atx_dic = {"ANTS": {}, "HEADER": []}
    ant_dic = None
    freq_dic = None
    freq = None
    header = []
    skip_to_next_ant = False

    # Read the file content
    with open(filepath, "r") as f:
        lines = f.readlines()

    in_header = True
    for line in lines:
        label = line[60:].strip()

        # Process header lines
        if in_header:
            atx_dic["HEADER"].append(line)
            if label == "END OF HEADER":
                in_header = False
            continue

        # Start a new antenna block
        if label == "START OF ANTENNA":
            skip_to_next_ant = False
            ant_dic = {
                "COMMENT": [],
                "FREQS": {},
                "METH": "",
                "TYPE": "",
                "SERIAL": "",
                "SVN": "",
                "COSPAR": "",
                "SINEX_CODE": "",
                "NFREQ": 0,
                "ZEN": 0,
            }

        # Skip processing if the antenna type does not match
        if skip_to_next_ant:
            continue

        # Parse antenna type and serial number
        elif label == "TYPE / SERIAL NO":
            if ant_type and not line[:20].strip().startswith(ant_type):
                skip_to_next_ant = True

            ant_dic["TYPE"] = line[0:20].strip()
            ant_dic["SERIAL"] = line[20:40].strip()
            ant_dic["SVN"] = line[40:50].strip()
            ant_dic["COSPAR"] = line[50:60].strip()


        # Parse azimuth increment
        elif label == "DAZI":
            ant_dic["DAZI"] = float(line[:6])

        # Parse zenith angle range and increment
        elif label == "ZEN1 / ZEN2 / DZEN":
            ant_dic["ZEN"] = list(map(float, line[:60].split()))

        # Parse comments
        elif label == "COMMENT":
            ant_dic["COMMENT"].append(line[:60])

        # Parse the number of frequencies
        elif label == "# OF FREQUENCIES":
            ant_dic["NFREQ"] = int(line[:20].strip())

        # Parse SINEX code
        elif label == "SINEX CODE":
            ant_dic["SINEX_CODE"] = line[:20].strip()

        # Parse method, author, and date
        elif label == "METH / BY / # / DATE":
            ant_dic["METH"] = line[0:20].strip()
            ant_dic["BY"] = line[20:40].strip()
            ant_dic["NCAL"] = line[40:46].strip()
            ant_dic["DATE"] = line[50:60].strip()

        # Start a new frequency block
        elif label == "START OF FREQUENCY":
            freq = line[:20].strip()
            ant_dic["FREQS"][freq] = {"AZI": [], "NOAZI": "", "PCO": []}
            freq_dic = ant_dic["FREQS"][freq]

        # End the current frequency block
        elif label == "END OF FREQUENCY":
            freq = None

        # Parse phase center offset (PCO) values
        elif freq and label == "NORTH / EAST / UP":
            freq_dic["PCO"] = list(map(float, line[:60].split()))

        # Parse elevation-dependent phase center variations (PCV) for NOAZI
        elif freq and line.strip().startswith("NOAZI"):
            values = np.array(list(map(float, line[:].split()[1:])))
            freq_dic["NOAZI"] = values

        # Parse azimuth-dependent PCV values
        elif freq and not line.strip().startswith("NOAZI"):
            values = np.array(list(map(float, line[:].split())))
            freq_dic["AZI"].append(values)

        # End the current antenna block
        elif label == "END OF ANTENNA":
            for f in ant_dic["FREQS"].keys():
                if ant_dic["FREQS"][f]["AZI"]:
                    ant_dic["FREQS"][f]["AZI"] = np.vstack(ant_dic["FREQS"][f]["AZI"])
            atx_dic["ANTS"][ant_dic["TYPE"]] = ant_dic
            ant_dic = None

    return atx_dic


def l_atx_std(val_inp, label_inp, newline=True):
    """
    Format a standard ANTEX line with a value and label.

    Parameters
    ----------
    val_inp : str
        The value to be written in the line.
    label_inp : str
        The label to be written in the line.

    Returns
    -------
    str
        Formatted ANTEX line.
    """
    nl = "\n" if newline else  ""
    return f"{val_inp:60}{label_inp:20}" + nl

def l_atx_cal(val_inp, noazi=True, newline=True):
    """
    Format calibration data for ANTEX.

    Parameters
    ----------
    val_inp : list
        Calibration values to be written.
    noazi : bool, optional
        If True, format as NOAZI calibration. Default is True.

    Returns
    -------
    str
        Formatted calibration data.
    """
    nl = "\n" if newline else  ""
    if noazi:
        n = len(val_inp)
        lin = "   NOAZI" + n * "{:+8.2f}" + nl
        return lin.format(*val_inp)
    else:
        lin_stk = []
        for lmat in val_inp:
            n = len(lmat)
            lin = "{:8.1f}" + (n - 1) * "{:+8.2f}"
            lin_stk.append(lin.format(*lmat))
        return "\n".join(lin_stk) + nl


def write_antex(atx_dic_inp, dir_out, fname_out="out.atx",erase_header=False):
    """
    Write antenna calibration information to an ANTEX file.

    Parameters
    ----------
    atx_dic_inp : dict
        Dictionary containing antenna calibration
        data to be written to the ANTEX file.
        See `read_antex` for dictionnary structure
    dir_out : str
        Directory where the output ANTEX file will be saved.
    fname_out : str, optional
        Name of the output ANTEX file. Default is 'out.atx'.
    erase_header : bool, optional
        If True, the header in the input dictionary will be replaced with a default header. Default is False.

    Notes
    -----
    - The function formats the antenna data into the ANTEX file format and writes it to the specified directory.
    - If `erase_header` is True, a default header is used instead of the one provided in the input dictionary.

    Returns
    -------
    str
        Path to the generated ANTEX file.
    """

    lout = []

    # Add header if present in the input dictionary
    if "HEADER" in atx_dic_inp and not erase_header:
        lhead = atx_dic_inp["HEADER"]
    else:
        lhead = [l_atx_std("     1.4            M", "ANTEX VERSION / SYST"),
                 l_atx_std("A", "PCV TYPE / REFANT"),
                 l_atx_std("", "END OF HEADER")]

    lhead = lhead[:-1] + [l_atx_std("HANDELED WITH GEODEZYX TOOLBOX - github.com/IPGP/geodezyx", "COMMENT")] + lhead[-1:]
    lout.append("".join(lhead))

    # Iterate over antennas in the input dictionary
    for kant, ant in atx_dic_inp["ANTS"].items():

        # Write antenna start label
        label = "START OF ANTENNA"
        lout.append(l_atx_std("", label))

        # Write antenna type and serial number
        label = "TYPE / SERIAL NO"
        val = f"{ant['TYPE']:20}{ant['SERIAL']:40}{ant['SVN']:20}{ant['COSPAR']:40}"
        lout.append(l_atx_std(val, label))

        # Write method, date, and other metadata
        label = "METH / BY / # / DATE"
        val = f"{ant['METH']:20}{ant['BY']:20}{ant['NCAL']:6}    {ant['DATE']:10}"
        lout.append(l_atx_std(val, label))

        label = "DAZI"
        val = f"  {ant['DAZI']:6.1f}"
        lout.append(l_atx_std(val, label))

        label = "ZEN1 / ZEN2 / DZEN"
        val = "  {:6.1f}{:6.1f}{:6.1f}".format(*ant["ZEN"])
        lout.append(l_atx_std(val, label))

        label = "# OF FREQUENCIES"
        nfreq = len(ant["FREQS"])
        if nfreq != ant['NFREQ']:
            print("# OF FREQUENCIES will be updated!")
            ant['NFREQ'] = nfreq
        val = f"{ant['NFREQ']:6}"
        lout.append(l_atx_std(val, label))

        # Facultative labels
        if "SINEX_CODE" in ant.keys():
            label = "SINEX CODE"
            val = f"{ant['SINEX_CODE']:60}"
            lout.append(l_atx_std(val, label))

        # Write comments
        if "COMMENT" in ant.keys():
            for com in ant["COMMENT"]:
                label = "COMMENT"
                val = f"{com:60}"
                lout.append(l_atx_std(val, label))

        # Write frequency-specific data
        for kfreq, freq in ant["FREQS"].items():

            label = "START OF FREQUENCY"
            val = f"   {kfreq:57}"
            lout.append(l_atx_std(val, label))

            label = "NORTH / EAST / UP"
            val = "{:10.2f}{:10.2f}{:10.2f}".format(*freq["PCO"])
            lout.append(l_atx_std(val, label))

            if len(freq["NOAZI"]) > 0:
                lout.append(l_atx_cal(freq["NOAZI"], noazi=True))
            if len(freq["AZI"]) > 0:
                lout.append(l_atx_cal(freq["AZI"], noazi=False))

            label = "END OF FREQUENCY"
            val = f"   {kfreq:57}"
            lout.append(l_atx_std(val, label))

        # Write antenna end label
        label = "END OF ANTENNA"
        lout.append(l_atx_std("", label))

    # Write the output to the specified file
    path_out = os.path.join(dir_out, fname_out)
    f_out = open(path_out, "w+")
    for l in lout:
        f_out.write(l)

    return path_out

pigs = "/home/psakicki/Downloads/igs20(2).atx"


p = "/home/psakicki/Downloads/AS-ANT3BCAL_NONE.atx"
antname = "AS-ANT3BCAL     NONE"

p = "/home/psakicki/Downloads/igs20(2).atx"
antname = "TPSPG_A1+GP     NONE"

atx_igs_dic = read_antex(pigs)

atx_dic = read_antex(p, antname)

atx_dic_use = atx_dic.copy()
#atx_dic_use["HEADER"] = atx_dic_use["HEADER"][:2]
freq_use = atx_dic_use["ANTS"][antname]["FREQS"]
atx_dic_use["ANTS"][antname]["COMMENT"].append("*** PSK (IPGP-OVS) 2025-04:")
atx_dic_use["ANTS"][antname]["COMMENT"].append("*** MANUAL DUPLICATION OF MULTI-GNSS SIGNALS")
freq_use["E01"] = freq_use["G01"]
freq_use["C01"] = freq_use["G01"]
freq_use["J01"] = freq_use["G01"]
freq_use["C02"] = freq_use["G02"]
freq_use["G05"] = freq_use["G02"]
freq_use["E05"] = freq_use["G02"]
freq_use["C05"] = freq_use["G02"]
freq_use["E06"] = freq_use["G02"]
freq_use["C06"] = freq_use["G02"]
freq_use["R06"] = freq_use["G02"]
freq_use["J06"] = freq_use["G02"]
freq_use["E07"] = freq_use["G02"]
freq_use["C07"] = freq_use["G02"]
freq_use["E08"] = freq_use["G02"]
freq_use["C08"] = freq_use["G02"]

# f = open("/home/psakicki/Downloads/igs20(2).atx")

write_antex(atx_dic_use, "/home/psakicki/aaa_FOURBI",
            antname + ".atx",
            erase_header=True)

p1 = "/home/psakicki/aaa_FOURBI/AS-ANT3BCAL     NONE.atx"
p2 = "/home/psakicki/aaa_FOURBI/TPSPG_A1+GP     NONE.atx"

atx_dic_perso1 = read_antex(p1)
atx_dic_perso2 = read_antex(p2)

import mergedeep
merged_dic = dict(mergedeep.merge({},atx_igs_dic,atx_dic_perso1,atx_dic_perso2))
merged_dic['HEADER'] = atx_igs_dic['HEADER']
write_antex(merged_dic, "/home/psakicki/aaa_FOURBI",
            "merged" + ".atx",
            erase_header=False)