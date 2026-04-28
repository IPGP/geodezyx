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


#### Import the logger
import logging
log = logging.getLogger("geodezyx")



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

def pandas_df_2_tuple_serie(df_inp, columns_name_list, reset_index_first=False):
    """
    Solve the multiple columns selection problem.

    Parameters
    ----------
    df_inp : pandas.DataFrame
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
        df_inp = df_inp.reset_index(level=0, inplace=False)
    else:
        df_inp = df_inp
        
    sout = pd.Series(list(map(tuple, df_inp[columns_name_list].values.tolist())),index=df_inp.index)
    return sout


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
        dif = np.diff(df[col_name].values) / (np.diff(df.index).astype(np.float32) * 10 ** -9) ## because it is in nanosec per def
        out = pd.Series(np.insert(dif,0,np.nan),  ## add NaN as the 1st value
                        index=df[col_name].index,
                        name=col_name)
    return out

        

def df_print(df_inp):
    string = df_inp.to_string()
    print(string)
    return string


def df_reg_2_multidx(df_inp, index_order):
    """
    Convert a regular DataFrame to a multi-index DataFrame.

    General function for creating multi-index DataFrames from regular columns.
    Can be used for any DataFrame type (orbits, clocks, stations, etc.).

    Parameters
    ----------
    df_inp : DataFrame
        Input DataFrame.
    index_order : list of str
        List of column names to set as multi-index, in order.

    Returns
    -------
    DataFrame
        Multi-index DataFrame.

    Examples
    --------
    >>> df = pd.DataFrame({'prn': ['G01', 'G01'], 'epoch': [...], 'x': [1, 2]})
    >>> df_multi = df_reg_2_multidx(df, ['prn', 'epoch'])
    """
    df_wrk = df_inp.reset_index()
    df_wrk = df_wrk.sort_values(index_order)
    df_wrk = df_wrk.set_index(index_order, inplace=False)
    return df_wrk


def df_common_finder(
    df_a_inp,
    df_b_inp,
    return_index=False,
    supplementary_sort=False,
    order=None,
    skip_reg2multidx_df_a=False,
    skip_reg2multidx_df_b=False,
):
    """
    Find common index entries in two DataFrames and output the corresponding filtered DataFrames.

    This is a general-purpose function that finds common rows based on a multi-index.
    It can be used for any DataFrame type (orbits, clocks, stations, etc.).

    Parameters
    ----------
    df_a_inp : DataFrame
        The first input DataFrame.
    df_b_inp : DataFrame
        The second input DataFrame.
    return_index : bool, optional
        If True, the function also returns the common index. Default is False.
    supplementary_sort : bool, optional
        If True, an additional sort is performed. This is useful for multi GNSS where the output DataFrame may not be well sorted. Default is False.
    order : list of str, optional
        The order of the index for the multi-index DataFrame.
        If None, the function will infer it from the input DataFrames.
        Default is None.
    skip_reg2multidx_df_a : bool, optional
        If True, skips the conversion of the first input DataFrame to a multi-index DataFrame. Default is False.
        The inputs are assumed to be already in multi-index format to optimize execution speed.
        (For advanced use only)
    skip_reg2multidx_df_b : bool, optional
        If True, skips the conversion of the second input DataFrame to a multi-index DataFrame. Default is False.
        The inputs are assumed to be already in multi-index format to optimize execution speed.
        (For advanced use only)

    Returns
    -------
    df_a_out : DataFrame
        The first output DataFrame with common index entries.
    df_b_out : DataFrame
        The second output DataFrame with common index entries.
    iinter : Index, optional
        The common index. Only returned if return_index is True.

    Notes
    -----
    This function is designed to be flexible enough for various use cases:
    - Orbits/SP3 files with satellite and epoch as index
    - Clocks with station/satellite name and epoch as index
    - SNX files with station and epoch as index
    - Any other DataFrame type with a natural multi-index structure

    Examples
    --------
    >>> df_a, df_b = df_common_finder(df_orb_a, df_orb_b, order=['prn', 'epoch'])
    >>> df_a, df_b, idx = df_common_finder(df_a, df_b, order=['col1', 'col2'], return_index=True)
    """
    # Auto-detect index order if not provided
    if order is None:
        # Use the existing multi-index if available
        if isinstance(df_a_inp.index, pd.MultiIndex):
            order = list(df_a_inp.index.names)
        else:
            raise ValueError("order parameter must be provided if DataFrame is not already multi-indexed")

    # Convert to multi-index if needed
    if not skip_reg2multidx_df_a:
        df_a = df_reg_2_multidx(df_a_inp, index_order=order)
    else:
        df_a = df_a_inp

    if not skip_reg2multidx_df_b:
        df_b = df_reg_2_multidx(df_b_inp, index_order=order)
    else:
        df_b = df_b_inp

    # Find common index entries
    i1 = df_a.index
    i2 = df_b.index

    iinter = i1.intersection(i2)
    iinter = iinter.sort_values()

    # Filter DataFrames to keep only common entries
    df_a_out = df_a.loc[iinter]
    df_b_out = df_b.loc[iinter]

    # Additional sorting if requested
    if supplementary_sort:
        df_a_out = df_a_out.sort_values(order)
        df_b_out = df_b_out.sort_values(order)

    # Warn if length mismatch
    if len(df_a_out) != len(df_b_out):
        log.warning("len(df_a_out) != len(df_b_out)")
        log.warning("TIPS: df_a_in and/or df_b_in might contain duplicates")

    if return_index:
        return df_a_out, df_b_out, iinter
    else:
        return df_a_out, df_b_out

