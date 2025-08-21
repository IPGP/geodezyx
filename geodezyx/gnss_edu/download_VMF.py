# -*- coding: utf-8 -*-
"""
"""
import wget



def download_VMF(year=None,month=None,day=None): 
# =============================================================================
#  This subroutine downloads VMF files from ^1e VMF ^^rver. Here, only operational VMF1 data downloads. 
#  To download VMF1_FC, 'VMF_OP' should be changed to 'VMF1_FC'. 

#  File created by Zohreh Adavi zohreh.adavi@tuwien.ac.at, 2023-02-20 
# =============================================================================  
   LISTt=['00','06','12','18']
   for lstnum in range(0,len(LISTt)):
         try:
            filename='VMFG_'+str(year)+str(month)+str(day).zfill(2)+'.H'+LISTt[lstnum] 
            url="https://vmf.geo.tuwien.ac.at/trop_products/GRID/2.5x2/VMF1/VMF1_OP/"+str(year)+'/'+filename
            wget.download(url)
         except:
             pass