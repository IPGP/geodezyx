import geodezyx
import glob
import datetime as dt

p = "/home/ovsgnss/090_TEMP_STUFFS/HOUE00GLP_S_20240050000_01D_30S_MO.rnx" 
p = "/scratch/calcgnss/data_spotgins_labeli/METS00FIN/*"
p = "/scratch/calcgnss/data_spotgins_labeli/ALBH00CAN/ALBH00CAN_R_2023*_01D_30S_MO.crx.gz"
p = "/scratch/calcgnss/data_spotgins_labeli/ALBH00CAN/albh*1?d*"
p = "/scratch/calcgnss/data_spotgins_labeli/"
p = "/root/020_BDRNX/data_ovs_glass/GL"
archive_folder = "/root/030_RESULTS/spotgins_singu_2503a" 

specific_sites=['HOUZ00GLP']
specific_sites=['ALBH00CAN',  'CASC00PRT', 'HOB200AUS',  'METS00FIN', 'STJO00CAN', 'ZIMM00CHE']

l = glob.glob(p)

l = geodezyx.operational.rinex_finder(p,
                                      specific_sites=specific_sites,
                                      start_epoch=dt.datetime(2000,5,3),
                                      end_epoch=dt.datetime(2024,12,31),
)


print("# OF RINEXS",len(l))

#prefix = "SPOTGINS_TEST"

#stafil = "/home/ovsgnss/010_SOFTS/spotgins/metadata/stations/station_file.dat"
#ocloadfil = "/home/ovsgnss/010_SOFTS/spotgins/metadata/oceanloading/load_fes2022_cf.spotgins"
#ocloadfil = "/home/ovsgnss/010_SOFTS/spotgins/metadata/oceanloading/load_fes2014b_cf.spotgins"
#optprafil = "/home/ovsgnss/010_SOFTS/spotgins/metadata/directeur/options_prairie_static"

#geodezyx.operational.spotgins_run(l, results_folder_inp=arcfld,force=False, nprocs=8, stations_file_inp=stafil)
geodezyx.operational.spotgins_run(rnxs_path_inp=l,
                                  results_folder_inp=archive_folder,
                                  nprocs=32,
                                  force=False,
                                  updatebd_login='sakic')
