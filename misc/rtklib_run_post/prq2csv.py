import pandas as pd

#inp_prq = "/home/sakic/IPGP_WORK/OVS/GNSS_OVS/2601_OVPF_erruption_all/040_pos_from_volcalc/OVPF_static-start_ULT_all.parquet"
#out_dir = "/home/sakic/IPGP_WORK/OVS/GNSS_OVS/2601_OVPF_erruption_all/050_csv/"

inp_prq = "/backuped/calcgnss/rtklib_process/OVPF_static-start_ULT/OVPF_static-start_ULT_all.parquet"
out_dir = "/backuped/calcgnss/rtklib_process/OVPF_static-start_ULT/CSV/"

# read the merged parquet file
df_raw = pd.read_parquet(inp_prq, engine="auto")

sample = "15min"

for (rov, bas), df_grp in df_raw.groupby(['rover', 'base']):
    print(rov, bas)
    df_epo = df_grp.set_index("epoch")
    df_med = df_epo[['x', 'y', 'z']].resample(sample).median()
    df_med.reset_index(inplace=True)
    df_med.drop_duplicates(inplace=True)
    df_med.to_csv(f"{out_dir}/{rov}_{bas}_{sample}.csv", index=False)
