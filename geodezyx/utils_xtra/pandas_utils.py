#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.utils contains functions for operations 
related to Python's Pandas object manipulations. 

it can be imported directly with:
from geodezyx import utils

The geodezyx toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/IPGP/geodezyx
"""


########## BEGIN IMPORT ##########
#### External modules
import pandas as pd
import numpy as np


#### geodeZYX modules


##########  END IMPORT  ##########



def renamedic_fast_4_pandas(*inpnames):
    """
    Create rename dictionary for Pandas columns.

    Parameters
    ----------
    *inpnames : str
        Column names

    Returns
    -------
    dict
        Dictionary mapping column indices to new names

    Notes
    -----
    Example::

        rnamedic = utils.renamedic_fast_4_pandas(*["zmax","ang","zsmooth","smoothtype","xgrad","ygrad",
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
    Wrapper of renamedic_fast_4_pandas for Pandas.

    Parameters
    ----------
    *inpnames : str
        Column names to rename

    Returns
    -------
    dict
        Dictionary mapping column indices to new names

    Notes
    -----
    Example::

        rnamedic = utils.renamedic_fast_4_pandas(*["zmax","ang","zsmooth","smoothtype","xgrad","ygrad",
                                           'r_eiko','z_eiko','pt_eiko_x','pt_eiko_y',"t_eiko",
                                           'r_sd',  'z_sd',  'pt_sd_x'  ,'pt_sd_y'  ,'t_sd',
                                           'diff_x','diff_y','diff','diff_t'])

        pda = pda.rename(columns = rnamedic)
    """
    return renamedic_fast_4_pandas(*inpnames)

def pandas_DF_2_tuple_serie(DFin,columns_name_list,reset_index_first=False):
    """
    Solve the multiple columns selection problem.

    Parameters
    ----------
    DFin : pandas.DataFrame
        Input DataFrame
    columns_name_list : list
        List of column names to select
    reset_index_first : bool, optional
        Reset index before conversion (default False)

    Returns
    -------
    pandas.Series
        Series of tuples

    Notes
    -----
    The idea is::

        S1 = pandas_DF_2_tuple_serie(DF1, columns_name_list)
        S2 = pandas_DF_2_tuple_serie(DF2, columns_name_list)
        BOOL = S1.isin(S2)
        DF1[BOOL]

    References
    ----------
    https://stackoverflow.com/questions/53432043/pandas-dataframe-selection-of-multiple-elements-in-several-columns
    """
    if reset_index_first:
        DF = DFin.reset_index(level=0, inplace=False)
    else:
        DF = DFin
        
    Sout = pd.Series(list(map(tuple, DF[columns_name_list].values.tolist())),index=DF.index)
    return Sout


def weighted_average(df,data_col,weight_col,by_col):
    """
    Source
    ------
    https://stackoverflow.com/questions/31521027/groupby-weighted-average-and-sum-in-pandas-dataframe
    """
    df['_data_times_weight'] = df[data_col]*df[weight_col]
    df['_weight_where_notnull'] = df[weight_col]*pd.notnull(df[data_col])
    g = df.groupby(by_col)
    result = g['_data_times_weight'].sum() / g['_weight_where_notnull'].sum()
    del df['_data_times_weight'], df['_weight_where_notnull']
    return result


def diff_pandas(df, col_name, use_np_diff=False):
    """
    Differentiate a Pandas DataFrame, if index is time.

    This function calculates the difference between consecutive elements in a specified column of a DataFrame.
    The difference is divided by the difference in time (seconds) between the corresponding indices.
    This is essentially a derivative operation, assuming the index represents time.

    Parameters
    ----------
    df : pandas.DataFrame
        The input DataFrame. The index should represent time.

    col_name : str
        The name of the column in the DataFrame that you want to differentiate.

    use_np_diff : bool, optional
        If True, use Numpy's diff.
        Default is False.
        This option has a (much) faster execution speed.

    Returns
    -------
    pandas.DataFrame or numpy.array
        The differentiated column of the input DataFrame. The type of the return value depends on the
        'return_array' parameter. If 'return_array' is False (default), a DataFrame is returned. If
        'return_array' is True, a numpy array is returned.

    """
    if not use_np_diff:
        out = df[col_name].diff() / df[col_name].index.to_series().diff().dt.total_seconds()
    else:
        dif = np.diff(df[col_name].values) / np.diff(df.index).astype(np.float32) * 10 ** -9 ## because it is in nanosec per def
        out = pd.Series(np.insert(dif,0,np.nan),  ## add NaN as the 1st value
                        index=df[col_name].index,
                        name=col_name)
    return out

        

def pandas_DF_print(DFin):
    string = DFin.to_string()
    print(string)
    return string

