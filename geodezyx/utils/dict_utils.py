#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.utils contains functions for operations 
related to Python's dictionary manipulations. 

it can be imported directly with:
from geodezyx import utils

The geodezyx toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/IPGP/geodezyx
"""

#### Import the logger
import logging
log = logging.getLogger('geodezyx')

def dicts_merge(*dict_args):
    """
    Merge multiple dictionaries into a single dictionary.

    Performs a shallow copy and merge of any number of dictionaries with
    precedence going to key-value pairs in later dictionaries.

    Parameters
    ----------
    *dict_args : dict
        Variable number of dictionaries to merge.

    Returns
    -------
    dict
        A new merged dictionary.

    Warnings
    --------
    First values will be erased if the same key is present in following
    dictionaries. Later dictionaries override earlier ones.

    Notes
    -----
    See https://stackoverflow.com/questions/38987/how-can-i-merge-two-python-dictionaries-in-a-single-expression
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result

def dicts_of_list_merge_mono(dol1, dol2):
    """
    Merge two dictionaries of lists by combining list values for common keys.

    Parameters
    ----------
    dol1 : dict
        First dictionary with lists as values.
    dol2 : dict
        Second dictionary with lists as values.

    Returns
    -------
    dict
        Merged dictionary where lists from both input dictionaries are combined.

    Notes
    -----
    See https://stackoverflow.com/questions/1495510/combining-dictionaries-of-lists-in-python
    """
    keys = set(dol1).union(dol2)
    no = []
    return dict((k, dol1.get(k, no) + dol2.get(k, no)) for k in keys)


def dicts_of_list_merge(*dict_args):
    """
    Merge multiple dictionaries of lists into a single dictionary.

    Parameters
    ----------
    *dict_args : dict
        Variable number of dictionaries with lists as values.

    Returns
    -------
    dict
        Merged dictionary where lists from all input dictionaries are combined.

    See Also
    --------
    dicts_of_list_merge_mono : Merge two dictionaries of lists.
    """
    result = dict()
    for dictionary in dict_args:
        result = dicts_of_list_merge_mono(result,dictionary)
    return result


def dic_key_for_vals_list_finder(dic_in , value_in):
    """
    Find the key in a dictionary of lists that contains a given value.

    Parameters
    ----------
    dic_in : dict
        Dictionary with lists as values, e.g.:

            dic_in[key1] = [val1a, val1b]
            dic_in[key2] = [val2a, val2b, val2c]

    value_in : object
        Value to search for in the lists.

    Returns
    -------
    key or None
        The key associated with the list containing value_in, or None if not found.

    Warnings
    --------
    This function returns the first key found. The input dictionary should be
    injective (no duplicate values across lists) for predictable behavior.

    Notes
    -----
    Example: if value_in = val2b, the function returns key2.

    Uses log.warning() to report when no key is found for the given value.
    """
    for k , v in dic_in.items():
        if value_in in v:
            return k
    
    log.warning("no key found for value %s... None returned",value_in)
    return None
