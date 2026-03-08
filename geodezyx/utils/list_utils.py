#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.utils contains functions for operations
related to Python's list manipulations.

it can be imported directly with:
from geodezyx import utils

The geodezyx toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/IPGP/geodezyx
"""


########## BEGIN IMPORT ##########
#### Import the logger
import logging
import re
import os

import numpy as np

#### geodeZYX modules
log = logging.getLogger("geodezyx")


##########  END IMPORT  ##########


def is_listoflist(inp):
    """
    Check if inp is a list of list.

    Parameters
    ----------
    inp : iterable
        Input object to check.

    Returns
    -------
    bool
        True if inp contains at least one list or numpy array element,
        False otherwise.

    Examples
    --------
    >>> is_listoflist([[1, 2], [3, 4]])
    True
    >>> is_listoflist([1, 2, 3])
    False
    """
    # On fait l'hypothèse que c'est un iterable qui va dedans de toute façon
    return any(isinstance(el, (list, np.ndarray)) for el in inp)


def shrink_listoflist(lin):
    """
    Shrink a list of list if it contains only one sublist.

    If ``lin`` is a list of list and contains only one element,
    returns the inner sublist, e.g. ``[[a, b, c]]`` => ``[a, b, c]``.

    Parameters
    ----------
    lin : list
        Input list, potentially a list of lists.

    Returns
    -------
    list
        The single sublist if ``lin`` is a one-element list of lists,
        otherwise ``lin`` unchanged.

    Examples
    --------
    >>> shrink_listoflist([[1, 2, 3]])
    [1, 2, 3]
    >>> shrink_listoflist([[1, 2], [3, 4]])
    [[1, 2], [3, 4]]
    """
    if is_listoflist(lin) and len(lin) == 1:
        return lin[0]
    else:
        return lin


def uniqify_list(seq, idfun=None):
    """
    Remove duplicate elements from a sequence while preserving order.

    Parameters
    ----------
    seq : iterable
        Input sequence to deduplicate.
    idfun : callable, optional
        Function to extract the identifier for uniqueness comparison.
        If None, the elements themselves are used as identifiers.

    Returns
    -------
    list
        Deduplicated sequence with order preserved.

    Notes
    -----
    Based on: https://www.peterbe.com/plog/uniqifiers-benchmark
    """
    if idfun is None:

        def idfun(x):
            return x

    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        # in old Python versions:
        # if seen.has_key(marker)
        # but in new ones:
        if marker in seen:
            continue
        seen[marker] = 1
        result.append(item)
    return result


def uniqify_list_of_lists(l):
    """
    Remove duplicate sublists while preserving uniqueness.

    Parameters
    ----------
    l : list of lists
        Input list containing sublists.

    Returns
    -------
    list
        List of unique sublists.

    Notes
    -----
    Source: http://stackoverflow.com/questions/3724551/python-uniqueness-for-list-of-lists
    """
    return [list(x) for x in set(tuple(x) for x in l)]


def find_common_elts(*lists):
    """
    Find common elements across multiple lists.

    Parameters
    ----------
    *lists : list
        Variable number of input lists.

    Returns
    -------
    numpy.ndarray
        Sorted array of elements common to all input lists.
    """
    sett = set(lists[0])

    for l in lists[1:]:
        sett = sett.intersection(l)

    return np.array(sorted(list(sett)))


def uniq_and_sort(l, natural_sort=True):
    """
    Remove duplicates from a list and sort it.

    Parameters
    ----------
    l : list
        Input list to deduplicate and sort.
    natural_sort : bool, optional
        If True, use natural sorting (default). If False, use standard sorting.

    Returns
    -------
    list
        Sorted list with duplicate elements removed.
    """
    import natsort
    if natural_sort:
        return natsort.natsorted(list(set(l)))
    else:
        return sorted(list(set(l)))


def df_sel_val_in_col(df, col_name, col_val):
    """
    Select rows from a DataFrame where a column matches a specific value.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame.
    col_name : str
        Name of the column to filter by.
    col_val : scalar
        Value to match in the specified column.

    Returns
    -------
    pandas.DataFrame
        Filtered DataFrame containing only rows where col_name == col_val.
    """
    return df[df[col_name] == col_val]


def uniquetol(a, tol):
    """
    Find unique elements in an array within a tolerance threshold.

    Parameters
    ----------
    a : array-like
        Input array.
    tol : float
        Absolute tolerance for uniqueness comparison.

    Returns
    -------
    numpy.ndarray
        Array of unique elements within the specified tolerance.

    Notes
    -----
    Source: http://stackoverflow.com/questions/37847053/uniquify-an-array-list-with-a-tolerance-in-python-uniquetol-equivalent
    """
    a = np.array(a)
    return a[~(np.triu(np.abs(a[:, None] - a) <= tol, 1)).any(0)]


def uniquetol2(a, tol=10 ** -6):
    """
    Find unique elements in an array within a tolerance threshold (optimized version).

    Parameters
    ----------
    a : array-like
        Input array.
    tol : float, optional
        Tolerance for rounding before finding unique elements. Default is 10**-6.

    Returns
    -------
    numpy.ndarray
        Array of unique elements.

    Notes
    -----
    This is a faster alternative to uniquetol.
    Source: https://stackoverflow.com/questions/5426908/find-unique-elements-of-floating-point-array-in-numpy-with-comparison-using-a-d
    """
    return np.unique(a.round(decimals=4))


def groups_near_central_values(a, tol, b=None):
    """
    group elements of an array by proximity to unique central values.

    Parameters
    ----------
    a : array-like
        Input array to group.
    tol : float
        Absolute tolerance for grouping elements near central values.
    b : array-like, optional
        Auxiliary array corresponding to elements in A. Default is None.

    Returns
    -------
    list or tuple
        If B is None, returns a list of lists where each sublist contains
        elements from A grouped around a central value.
        If B is provided, returns a tuple (groups_A, groups_B) containing
        grouped elements from both arrays.

    Notes
    -----
    This function is in beta status and may have bugs if tolerance is poorly chosen.
    """
    a = np.array(a)

    auniq = uniquetol(a, tol)
    grouplis = []

    if b is None:
        bbool = False
    else:
        bbool = True

    if bbool:
        group_b_lis = []

    for iauniq, auniq in enumerate(auniq):
        group = []
        grouplis.append(group)

        if bbool:
            group_b = []
            group_b_lis.append(group_b)

        for ia, a in enumerate(a):
            count = 0
            if np.abs(a - auniq) <= tol:
                group.append(a)
                count += 1

                if bbool:
                    group_b.append(b[ia])

        if count > 1:
            log.warning("%s in %s groups", a, count)

    if bbool:
        return grouplis, group_b_lis
    else:
        return grouplis


def uniq_set_list(setlis, frozen=True):
    """
    Remove duplicate sets from a list of sets.

    Parameters
    ----------
    setlis : list of sets
        Input list containing sets or set-like iterables.
    frozen : bool, optional
        If True, returns frozensets (hashable and immutable). Default is True.

    Returns
    -------
    list
        List of unique sets (either sets or frozensets depending on frozen parameter).
    """
    if not frozen:
        return [set(e) for e in list(set([frozenset(e) for e in setlis]))]
    else:
        return [frozenset(e) for e in list(set([frozenset(e) for e in setlis]))]


def sort_binom_list(x, y, array_out=False):
    """
    Sort Y according to X and sort X.

    Parameters
    ----------
    x : list or array-like
        Reference values to sort by.
    y : list or array-like
        Values to sort according to X ordering.
    array_out : bool, optional
        If True, return numpy arrays. If False, return lists. Default is False.

    Returns
    -------
    tuple
        (xnew, ynew) where both are sorted according to X.
        Type depends on array_out parameter.

    Raises
    ------
    Warning
        If len(X) != len(Y), a warning is logged.
    """
    if len(x) != len(y):
        log.warning("len(X) != len(Y) !!!")

    ynew = [y for (x, y) in sorted(zip(x, y))]
    xnew = sorted(x)
    if not array_out:
        return xnew, ynew
    else:
        return np.array(xnew), np.array(ynew)


def sort_multinom_list(x, *y):
    """
    Sort multiple Y sequences according to X and sort X.

    Parameters
    ----------
    x : list or array-like
        Reference values to sort by.
    *y : list or array-like
        Variable number of sequences to sort according to X ordering.

    Returns
    -------
    tuple
        (xnew, Ynew_1, Ynew_2, ...) where all sequences are sorted according to X.
        X is returned as a numpy array, while Y sequences are returned as lists.
    """
    ynew_stk = []
    for YY in y:
        ynew = [y for (x, y) in sorted(zip(x, YY))]
        ynew_stk.append(ynew)
    xnew = sorted(x)
    fintup = (np.array(xnew),) + tuple(ynew_stk)
    return fintup


def sort_basename(file_paths):
    """
    Sort a list of file paths by their basenames.

    Parameters
    ----------
    file_paths : list
        List of file paths to be sorted.

    Returns
    -------
    list
        Sorted list of file paths by their basenames.
    """
    return list(sorted(file_paths, key=os.path.basename))


def most_common(lst):
    """
    Find the most frequently occurring element in a list.

    Parameters
    ----------
    lst : list or iterable
        Input sequence to analyze.

    Returns
    -------
    scalar
        The element with the highest frequency in the list.

    Notes
    -----
    Source: http://stackoverflow.com/questions/1518522/python-most-common-element-in-a-list
    """
    if type(lst) is not list:
        lst = list(lst)
    return max(set(lst), key=lst.count)


def median_improved(l):
    """
    Calculate the median of a list, handling even-length lists differently.

    Parameters
    ----------
    l : list or array-like
        Input sequence.

    Returns
    -------
    scalar
        The median value. For even-length lists, returns the nearest value
        in the list to the actual median instead of interpolating.

    Notes
    -----
    For even-length lists, does not return the mean of the two middle values
    but instead returns the nearest value from the input list.
    """
    l = np.array(l)
    if len(l) % 2 == 0:
        return find_nearest(l, np.median(l))[0]
    else:
        return np.median(l)


def trio_lists_2_tab(xlis, ylis, vlis):
    """
    Convert three lists into a 2D table structure.

    Parameters
    ----------
    xlis : list
        List of X values that will become column headers.
    ylis : list
        List of Y values that will become row headers.
    vlis : list
        List of data values corresponding to (X, Y) pairs.

    Returns
    -------
    list
        A 2D table structure where the first row contains X values,
        and subsequent rows contain Y value and corresponding V values.
        Compatible with the tabulate module.

    Notes
    -----
    The data lookup is performed with a brute-force approach (nested loops),
    which is not optimized for large datasets.

    Examples
    --------
    >>> trio_lists_2_tab([1, 2, 1, 2], [1, 1, 2, 2], [10, 20, 30, 40])
    [[1, 2], [1, 10, 20], [2, 30, 40]]
    """

    xlis_uniq = sorted(np.unique(xlis))
    ylis_uniq = sorted(np.unique(ylis))

    finalis = []
    finalis.append(xlis_uniq)

    for y in ylis_uniq:
        curlis = [y]
        finalis.append(curlis)
        for x in xlis_uniq:
            for xx, yy, vv in zip(xlis, ylis, vlis):
                if xx == x and yy == y:
                    curlis.append(vv)
    return finalis


def minmax(l):
    """
    Find the minimum and maximum values in a list.

    Parameters
    ----------
    l : list or array-like
        Input sequence.

    Returns
    -------
    tuple
        (min_value, max_value) of the input sequence.
    """
    return np.min(l), np.max(l)


def middle(linp):
    """
    Calculate the midpoints between consecutive elements of a list.

    Parameters
    ----------
    linp : list or array-like
        Input sequence with at least 2 elements.

    Returns
    -------
    list
        List of midpoint values between consecutive elements.
    """
    lout = []
    for i in range(len(linp) - 1):
        lout.append((linp[i] + linp[i + 1]) / 2)
    return lout


def second_smallest(numbers):
    """
    Find the second smallest value in a sequence.

    Parameters
    ----------
    numbers : iterable
        Input sequence to analyze.

    Returns
    -------
    scalar
        The second smallest element in the sequence.

    Notes
    -----
    Returns infinity if there are fewer than 2 elements.
    """
    m1, m2 = float("inf"), float("inf")
    for x in numbers:
        if x <= m1:
            m1, m2 = x, m1
        elif x < m2:
            m2 = x
    return m2


def find_index_multi_occurences(l, elt):
    """
    Find all indices where an element occurs in a list.

    Parameters
    ----------
    l : list
        Input list to search.
    elt : scalar
        Element to find.

    Returns
    -------
    list
        List of indices where elt appears in L.
    """
    return [i for i, x in enumerate(l) if x == elt]


def find_surrounding(L, v):
    """
    Find the two nearest values surrounding a target value.

    Parameters
    ----------
    L : iterable
        Input list/array to search.
    v : scalar
        Target value to find surrounding values for.

    Returns
    -------
    tuple
        (surrounding_values, surrounding_indices) where:
        - surrounding_values is a tuple of the two nearest values
        - surrounding_indices is a tuple of their indices in L
    """
    a = np.array(L)
    b = v

    surrounding_values = np.sort(
        [b + i for i in sorted(np.subtract(a, b), key=lambda x: abs(x))[:2]]
    )
    surrounding_index = (
        np.where(surrounding_values[0] == a)[0][0],
        np.where(surrounding_values[1] == a)[0][0],
    )

    return tuple(surrounding_values), surrounding_index


##################
### LIST IT FCTs
##################


def chunkIt(seq, num):
    """
    Divide a list into approximately num equal sublists.

    Parameters
    ----------
    seq : list or array-like
        Input sequence to divide.
    num : int
        Desired number of sublists.

    Returns
    -------
    list
        List of sublists, each roughly equal in size.

    Notes
    -----
    The sublists may vary in size by 1 element if the sequence length
    is not evenly divisible by num.

    Source: http://stackoverflow.com/questions/2130016/splitting-a-list-of-arbitrary-size-into-only-roughly-n-equal-parts
    """
    # http://stackoverflow.com/questions/2130016/splitting-a-list-of-arbitrary-size-into-only-roughly-n-equal-parts
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last) : int(last + avg)])
        last += avg
    return out


def sliceIt(seq, num):
    """
    Divide a list into sublists of fixed size.

    Parameters
    ----------
    seq : list or array-like
        Input sequence to divide.
    num : int
        Size of each sublist.

    Returns
    -------
    list
        List of sublists, each containing num elements (last sublist may be shorter).

    Notes
    -----
    Source: http://stackoverflow.com/questions/4501636/creating-sublists
    """
    # http://stackoverflow.com/questions/4501636/creating-sublists
    return [seq[i : i + num] for i in range(0, len(seq), num)]


def sublistsIt(seq, lenofsublis_lis, output_array=False):
    """
    Divide a sequence into sublists of specified sizes.

    Parameters
    ----------
    seq : list or array-like
        Input sequence to divide.
    lenofsublis_lis : list
        List of integers specifying the size of each sublist.
        Example: [2, 3, 4, 2] creates 4 sublists of sizes 2, 3, 4, and 2.
    output_array : bool, optional
        If True, return list of numpy arrays. If False, return list of lists.
        Default is False.

    Returns
    -------
    list
        List of sublists (or arrays if output_array=True).

    Raises
    ------
    Exception
        If sum(lenofsublis_lis) != len(seq).
    """

    if np.sum(lenofsublis_lis) != len(seq):
        raise Exception("sublistsIt : sum(lenofsublis_lis) != len(seq) ")
    seq = list(seq)
    sublis_lis = []
    start = 0
    for l in lenofsublis_lis:
        end = start + l
        sublis = seq[start:end]
        sublis_lis.append(sublis)
        start = end
    if output_array:
        return [np.array(e) for e in sublis_lis]
    else:
        return sublis_lis


def identical_consecutive_eltsIt(linp):
    """
    Group consecutive identical elements together.

    Parameters
    ----------
    linp : list or iterable
        Input sequence with potentially repeated consecutive elements.

    Returns
    -------
    list
        List of lists, where each inner list contains consecutive identical elements.
    """
    lout_big = []
    linter = [linp[0]]
    lout_big.append(linter)

    for e in linp[1:]:
        if e == linter[-1]:
            linter.append(e)
        else:
            linter = [e]
            lout_big.append(linter)

    return lout_big


def find_nearest(listin, value):
    """
    Find the nearest value in a list to a target value.

    Parameters
    ----------
    listin : list or array-like
        Input list/array to search.
    value : scalar
        Target value to find nearest element to.

    Returns
    -------
    tuple
        (nearest_value, index_of_nearest) where nearest_value is the element
        in listin closest to value, and index_of_nearest is its index.
    """
    array = np.array(listin)
    idx = (np.abs(array - value)).argmin()

    return listin[idx], idx


def find_interval_bound(listin, val, outindexes=True):
    """
    Find the bounding values/indices of an interval around a target value.

    Parameters
    ----------
    listin : list or array-like
        Input list/array (assumed to be sorted).
    val : scalar
        Target value to find bounds for.
    outindexes : bool, optional
        If True, return indices of bounds. If False, return the bounding values.
        Default is True.

    Returns
    -------
    tuple
        If outindexes is True, returns (lower_index, upper_index).
        If outindexes is False, returns (lower_value, upper_value).
    """
    import bisect
    i = bisect.bisect(listin, val)
    if outindexes:
        return i - 1, i
    else:
        return listin[i - 1], listin[i]


def occurence(l, tolerence=None, pretty_output=False):
    """
    Count occurrences of elements in a list.

    Parameters
    ----------
    l : list
        Input list
    tolerence : float, optional
        Tolerance to find close elements of L
        if no tolerance is given then a set() is used
    pretty_output : bool
        if False, return a list of 2-tuples:

            (element of the list, number of occurrence of this element in the list)

        if True, return tuple with sorted occurrences and values

    Returns
    -------
    output : list or tuple
        See pretty_output parameter

    Notes
    -----
    pretty_output is implemented because the first mode is not really useful (180612)
    the equal test is also replaced by is close
    """
    l = np.array(l)

    if tolerence:
        lset = uniquetol2(l, tol=tolerence)
    else:
        lset = set(l)

    outlisoccur = []
    for l in lset:
        outlisoccur.append((l, np.sum(l == l)))

    if pretty_output:
        mtmp = np.vstack(outlisoccur)
        vals = mtmp[:, 0]
        occurs = mtmp[:, 1]
        occurs, vals = sort_binom_list(occurs, vals)

        output = (occurs, vals)
    else:
        output = outlisoccur

    return output


def decimateIt(listinp, n):
    """
    Decimate a list by selecting every n-th element.

    Parameters
    ----------
    listinp : list or array-like
        Input sequence to decimate.
    n : int
        Decimation factor. Elements at indices where i % n == 0 are selected.

    Returns
    -------
    list
        Decimated list containing every n-th element.
    """
    outlist = []
    for i in range(len(listinp)):
        if np.mod(i, n) == 0:
            outlist.append(listinp[i])
    return outlist


def consecutive_groupIt(data, only_start_end=False):
    """
    Identify groups of continuous numbers in a list.

    Parameters
    ----------
    data : list or array-like
        Input sequence of numbers.
    only_start_end : bool, optional
        If True, return only (start, end) tuples for each group.
        If False, return full lists of elements in each group. Default is False.

    Returns
    -------
    list
        List of groups. Each group is either a list of consecutive elements
        (if only_start_end=False) or a tuple (start, end) (if only_start_end=True).

    Notes
    -----
    Useful for time periods with a prior conversion to MJD.

    Source :
    https://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
    """
    from operator import itemgetter
    import itertools

    groups_stk = []
    for k, g in itertools.groupby(enumerate(data), lambda A: A[0] - A[1]):
        groups_stk.append(list(map(itemgetter(1), g)))

    if only_start_end:
        return [(g[0], g[-1]) for g in groups_stk]
    else:
        return groups_stk


def identical_groupIt(data):
    """
    Group consecutive identical elements together.

    Parameters
    ----------
    data : list or iterable
        Input sequence with potentially repeated consecutive elements.

    Returns
    -------
    list
        List of lists, where each inner list contains consecutive identical elements.

    Notes
    -----
    Source :
    https://stackoverflow.com/questions/30293071/python-find-same-values-in-a-list-and-group-together-a-new-list
    """
    import itertools
    return [list(j) for i, j in itertools.groupby(data)]


def get_interval(start, end, delta):
    """
    Generate a list of values at regular intervals between start and end.

    Parameters
    ----------
    start : numeric
        Starting value (inclusive).
    end : numeric
        Ending value (exclusive).
    delta : numeric
        Step size between consecutive values.

    Returns
    -------
    list
        List of values from start to end with step delta.

    Notes
    -----
    Source: http://stackoverflow.com/questions/10688006/generate-a-list-of-datetimes-between-an-interval-in-python
    """
    # after http://stackoverflow.com/questions/10688006/generate-a-list-of-datetimes-between-an-interval-in-python
    outlist = []
    curr = start
    while curr < end:
        outlist.append(curr)
        curr += delta
    return outlist


def duplicates_finder(seq):
    """
    Find all duplicate elements in a sequence.

    Parameters
    ----------
    seq : iterable
        Input sequence to search for duplicates.

    Returns
    -------
    list
        List of elements that appear more than once in the input sequence.

    Notes
    -----
    Source: http://stackoverflow.com/questions/9835762/find-and-list-duplicates-in-python-list
    """
    seen = set()
    seen_add = seen.add
    # adds all elements it doesn't know yet to seen and all other to seen_twice
    seen_twice = set(x for x in seq if x in seen or seen_add(x))
    # turn the set into a list (as requested)
    return list(seen_twice)


def sort_table(table, col):
    """
    Sort a table by multiple columns.

    Parameters
    ----------
    table : list of lists or tuple of tuples
        where each inner list represents a row
    col : int
        column number to sort by

    Returns
    -------
    outtable : list
        sorted table
    """
    refcol = list(table[col])
    outtable = []
    for col in table:
        outtable.append([x for (y, x) in sorted(zip(refcol, col))])
    return outtable


def dicofdic(mat, names):
    """
    Create a 2D dictionary from a matrix and corresponding names.

    Parameters
    ----------
    mat : array-like
        N x N matrix of values.
    names : list
        List of N names to use as keys for both dimensions.

    Returns
    -------
    dict
        A nested dictionary where dic[name1][name2] = mat[i, j]
        where i and j are the indices corresponding to name1 and name2.

    Notes
    -----
    Source: http://stackoverflow.com/questions/13326042/2d-dictionary-with-multiple-keys-per-value
    """

    # a partir d'un matrice de N x N
    # et d'une liste de N noms
    # Fabrique un dictionnaire 2D dic[nom1][nom2] = M[n1,n2]
    # http://stackoverflow.com/questions/13326042/2d-dictionary-with-multiple-keys-per-value
    import collections
    d2_dict = collections.defaultdict(dict)

    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            d2_dict[names[i]][names[j]] = mat[i, j]

    return d2_dict


def find_regex_in_list(regex, L, only_first_occurence=False, line_number=False):
    """
    Find elements in a list matching a regular expression pattern.

    Parameters
    ----------
    regex : str
        Regular expression pattern to search for.
    L : list
        List of strings to search.
    only_first_occurence : bool, optional
        If True, return only the first match. If False, return all matches.
        Default is False.
    line_number : bool, optional
        If True, return tuples of (index, element). If False, return just elements.
        Default is False.

    Returns
    -------
    list or scalar
        If only_first_occurence=True, returns a single match (element or tuple).
        If only_first_occurence=False, returns a list of matches.
        Format depends on line_number parameter.
    """
    Lout = []
    for i, e in enumerate(L):
        if re.search(regex, e):
            if not line_number:
                found = e
            else:
                found = (i, e)

            if only_first_occurence:
                return found
            else:
                Lout.append(found)
    return Lout
