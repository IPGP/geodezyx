#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 17:15:06 2019

@author: psakicki
"""

def renamedic_fast_4_pandas(*inpnames):
    """
    EXEMPLE :
    rnamedic = gf.renamedic_fast_4_pandas(*["zmax","ang","zsmooth","smoothtype","xgrad","ygrad",
                                       'r_eiko','z_eiko','pt_eiko_x','pt_eiko_y',"t_eiko",
                                       'r_sd',  'z_sd',  'pt_sd_x'  ,'pt_sd_y'  ,'t_sd',
                                       'diff_x','diff_y','diff','diff_t'])

    pda = pda.rename(columns = rnamedic)
    """
    renamedic = dict()

    for i,nam in enumerate(inpnames) :
        renamedic[i] = nam

    return renamedic

def pandas_column_rename_dic(*inpnames):
    """
    wrapper of renamedic_fast_4_pandas

    EXEMPLE :
    rnamedic = gf.renamedic_fast_4_pandas(*["zmax","ang","zsmooth","smoothtype","xgrad","ygrad",
                                       'r_eiko','z_eiko','pt_eiko_x','pt_eiko_y',"t_eiko",
                                       'r_sd',  'z_sd',  'pt_sd_x'  ,'pt_sd_y'  ,'t_sd',
                                       'diff_x','diff_y','diff','diff_t'])

    pda = pda.rename(columns = rnamedic)
    """
    return renamedic_fast_4_pandas(*inpnames)

def pandas_DF_2_tuple_serie(DFin,columns_name_list,reset_index_first=False):
    """
    This function is made to solve the multiple columns selection 
    problem
    the idea is :
        S1 = pandas_DF_2_tuple_serie(DF1 , columns_name_list)
        S2 = pandas_DF_2_tuple_serie(DF2 , columns_name_list)
        BOOL = S1.isin(S2)
        DF1[BOOL]
        
    Source :
        https://stackoverflow.com/questions/53432043/pandas-dataframe-selection-of-multiple-elements-in-several-columns
    """
    if reset_index_first:
        DF = DFin.reset_index(level=0, inplace=False)
    else:
        DF = DFin
        
    Sout = pd.Series(list(map(tuple, DF[columns_name_list].values.tolist())),index=DF.index)
    return Sout