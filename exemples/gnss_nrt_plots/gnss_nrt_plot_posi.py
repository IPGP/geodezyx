#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 07/01/2026 11:42:12

@author: psakic
"""

import pandas as pd
import matplotlib.pyplot as plt

from geodezyx import utils, conv, stats, utils_xtra
import numpy as np
import matplotlib as mpl


xyz_dic = {
    "BOMG": [3352122.640732463, 4915723.753034734, -2297076.907396815],
    "BORG": [3352179.468244755, 4915389.35801885, -2297785.106621451],
    "C98G": [3349647.071016342, 4918290.539762681, -2292685.261361883],
    "CASG": [3342033.023264375, 4922773.708655806, -2290074.971696334],
    "CFNG": [3349537.184260634, 4916753.630176714, -2296761.787795528],
    "CRAG": [3347500.97944078, 4917914.352085327, -2294826.506682635],
    "DERG": [3351077.500251723, 4916251.558666458, -2297303.247570688],
    "DSRG": [3351348.575051963, 4915747.889372488, -2298064.221329384],
    "ENCG": [3354597.116872568, 4913638.028335221, -2297376.880814932],
    "ENOG": [3353268.038184023, 4914797.26627217, -2296549.933107008],
    "FEUG": [3340804.842314329, 4916311.646266831, -2305336.814085809],
    "FJAG": [3350995.968928102, 4916475.044280337, -2295856.459849963],
    "FOAG": [3350528.01283789, 4915139.525874629, -2299179.008076069],
    "FREG": [3353758.5735742, 4915576.557677006, -2292887.943071418],
    "GBNG": [3343297.13748501, 4919838.611709842, -2294437.898300031],
    "GBSG": [3344238.331729389, 4917088.44551912, -2299848.399194048],
    "GITG": [3354343.97136548, 4914925.741250496, -2294766.595995225],
    "GPNG": [3347812.884805933, 4917375.372526649, -2296671.256797257],
    "GPSG": [3346186.407024792, 4916676.421682362, -2299378.496052825],
    "HDLG": [3343663.622590312, 4918310.037725286, -2297365.60844327],
    "MAIG": [3383479.378135485, 4901531.448891763, -2280412.059167621],
    "PBRG": [3357258.662478508, 4913169.130239947, -2294728.333650563],
    "PRAG": [3350809.181810517, 4913568.597968383, -2302266.991206395],
    "PVDG": [3348379.438174965, 4916189.740796518, -2299156.356559711],
    "RVAG": [3352538.751430386, 4914422.556831216, -2298369.016684549],
    "SNEG": [3351359.653965392, 4916209.867089082, -2296999.327082343],
    "TRCG": [3342109.774250217, 4917359.715981821, -2301525.627671301],
}


if utils.get_computer_name() == 'HPEB8a':
    tot_prq_path = "/home/sakic/aaa_FOURBI/OVPF_static-start_ULT_01min_03days_all.parquet"
    outdir_plots = "/home/sakic/IPGP_WORK/OVS/GNSS_OVS/2601_OVPF_erruption_all/060_plots"
    mpl_cfg = "/home/sakic/CODES/geodezyx_toolbox_PS_perso_scripts/MISC/matplotlib_env/2602_new_on_HPEB8a/matplotlibrc_PSmpl02a.rc"
elif utils.get_computer_name() == 'volcalcgnss':
    tot_prq_path = "/backuped/calcgnss/rtklib_results/OVPF_static-start_ULT/07days_01min/OVPF_static-start_ULT_all.parquet"
    outdir_plots = "/backuped/calcgnss/rtklib_results/OVPF_static-start_ULT/07days_01min/plots"
    mpl_cfg = "/opt/softs_gnss/geodezyx/misc/matplotlib_env/matplotlibrc_PSmpl02a.rc"
    
mpl.rc_file(mpl_cfg)

out_exts = (".png",".eps",)

# read the merged parquet file
df_raw = pd.read_parquet(tot_prq_path, engine="auto")
df_all = df_raw

df_all = df_all[df_all["base"] == "GITG"]


### decimation at minute
# df_all = df_all[df_all["epoch"].dt.second == 0]

#START = dt.datetime(2026,3,8)
#df_all = df_all[(pd.Timestamp(START) < df_all["epoch"])]

rovers = df_all["rover"].unique()

shift = 2 # cm


figall, axall = plt.subplots(3,1)
        
for irov , (rov, df_rovbas_ori) in enumerate(df_all.groupby("rover")):
    fig, ax = plt.subplots(3,1)

    df_rovbas = df_rovbas_ori.set_index("epoch")
    # resample to 1 minute and take median
    df_rovbas = df_rovbas.resample("1min").median(numeric_only=True)

    df_rovbas[["e","n","u"]] = conv.xyz2enu_vector(df_rovbas[["x","y","z"]],
                                                    xyz_dic[rov])

    df_rovbas_cln, boolbad = stats.outlier_mad_df(df_rovbas,
                                            ["e","n","u"],
                                            3.5,
                                            np.logical_and)

    for i,enu in enumerate(["e","n","u"]):
        t = df_rovbas_cln.index

        y = df_rovbas_cln[enu] - df_rovbas_cln[enu].median()
        y = y * 100 # m > cm
        
        ax[i].plot(t,y,
                   marker=".",
                   )
        
        axall[i].plot(t,y + shift * irov,
                   marker=".",
                   label=rov
                   )
        
        ax[i].set_ylabel(enu.upper() + " (cm)")
        axall[i].set_ylabel(enu.upper() + " (cm)")

        ylims = .015
        ylims_margin = ylims + 0.01

        if np.min(y) > -ylims_margin and np.max(y) < ylims_margin:
            ax[i].set_ylim((-ylims,ylims))
        
        axall[i].legend()

    last_epoc = df_all["epoch"].max()
    last_epoc_str = conv.dt2str_iso(last_epoc)
    now_str = conv.dt2str_iso(conv.now('utc'))
        
    ax[0].set_title(f"generated: {now_str}, last epoch: {last_epoc_str}")
    fig.suptitle("Position variation " + rov)
    
    axall[0].set_title(f"generated: {now_str}, last epoch: {last_epoc_str}")
    figall.suptitle("Position variation - Summary")
    figall.tight_layout()


    utils_xtra.plot_utils.figure_saver(
        fig,
        outdir_plots,
        "PDF_GNSS_NRT_posi_" + rov,
        dpi=400, 
        outtype=out_exts,
        formt="A4",
        tight_layout=True
        
    )
    
    
utils_xtra.plot_utils.figure_saver(
    figall,
    outdir_plots, 
    "PDF_GNSS_NRT_posi_all",
    dpi=400, 
    outtype=out_exts,
    formt="A4",
    tight_layout=True
)




    
