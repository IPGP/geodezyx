#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 14:27:56 2019

@author: psakicki
"""

from geodezyx.extern import * 

def is_iterable(inp,consider_str_as_iterable=False):
    """
    Test if the input is an iterable like a list or a numpy array or not

    Parameters
    ----------
    inp : list, numpy.array, ...

    consider_str_as_iterable : bool
        string are considered as iterable by Python per default
        This boolean will avoid True as return if you test a string
        
    Returns
    -------
    out : bool
        True if inp is iterable, False either
    """
    if not consider_str_as_iterable and type(inp) is str:
        return False

    try:
        iter(inp)
    except TypeError:
        out = False
    else:
        out = True
    return out

def is_not_iterable(inp,consider_str_as_iterable=False):
    """
    Simple negation of is_iterable()
    """
    if not consider_str_as_iterable and type(inp) is str:
        return False

    try:
        iter(inp)
    except TypeError:
        out = True
    else:
        out = False
    return out
