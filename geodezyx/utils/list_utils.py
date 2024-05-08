#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.utils contains functions for operations 
related to Python's list manipulations. 

it can be imported directly with:
from geodezyx import utils

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU LGPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (IPGP, sakic@ipgp.fr)
GitHub repository :
https://github.com/GeodeZYX/geodezyx-toolbox
"""


########## BEGIN IMPORT ##########
#### External modules
import bisect
import collections
import itertools
#### Import the logger
import logging
import re

import natsort
import numpy as np

#### geodeZYX modules
log = logging.getLogger(__name__)


##########  END IMPORT  ##########


def is_listoflist(inp):
    """
    check if inp is a list of list [[...] , [...] , ... ,[...]]
    """
    # On fait l'hypothèse que c'est un iterable qui va dedans de toute façon
    return any(isinstance(el, (list,np.ndarray)) for el in inp)

def shrink_listoflist(lin):
    """
    if lin (list in) is a list of list and contains only one element
    then returns the sublist, e.g. : [[a,b,c]] => [a,b,c]
    """
    if is_listoflist(lin) and len(lin) == 1:
        return lin[0]
    else:
        return lin
    
def uniqify_list(seq, idfun=None):
   """
   order preserving uniq function
   based on
   https://www.peterbe.com/plog/uniqifiers-benchmark
   """
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in seq:
       marker = idfun(item)
       # in old Python versions:
       # if seen.has_key(marker)
       # but in new ones:
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result


def uniqify_list_of_lists(L):
    """
    source
    http://stackoverflow.com/questions/3724551/python-uniqueness-for-list-of-lists
    """
    return [list(x) for x in set(tuple(x) for x in L)]


def find_common_elts(*lists):
    """
    Find common elements in different lists
    """
    sett = set(lists[0])

    for l in lists[1:]:
        sett  = sett.intersection(l)

    return np.array(sorted(list(sett)))


def uniq_and_sort(L,natural_sort=True):
    """
    In a list, remove duplicates and sort the list
    """
    if natural_sort:
        return natsort.natsorted(list(set(L)))
    else:
        return sorted(list(set(L)))


def df_sel_val_in_col(DF,col_name,col_val):
    """
    Return a selected value of a column in a DataFrame
    i.e.
    return DF[DF[col_name] == col_val]
    """
    return DF[DF[col_name] == col_val]
              
              
def uniquetol(A,tol):
    """
    tol : absolute tolerance
    source
    http://stackoverflow.com/questions/37847053/uniquify-an-array-list-with-a-tolerance-in-python-uniquetol-equivalent
    """
    A = np.array(A)
    return A[~(np.triu(np.abs(A[:,None] - A) <= tol,1)).any(0)]


def uniquetol2(A,tol=10**-6):
    """
    This one is faster
    https://stackoverflow.com/questions/5426908/find-unique-elements-of-floating-point-array-in-numpy-with-comparison-using-a-d
    """
    return np.unique(A.round(decimals=4))

def groups_near_central_values(A,tol,B=None):
    """
    beta ... bug if tol is bad
    tol : absolute tolerance
    170514 : B is an annex list
    """
    A = np.array(A)

    Auniq = uniquetol(A,tol)
    Grouplis = []

    if B is None:
        Bbool = False
    else:
        Bbool = True

    if Bbool:
        GrouplisB = []

    for iauniq , auniq in enumerate(Auniq):
        Group = []
        Grouplis.append(Group)

        if Bbool:
            GroupB = []
            GrouplisB.append(GroupB)



        for ia , a in enumerate(A):
            count = 0
            if np.abs(a - auniq) <= tol:
                Group.append(a)
                count += 1

                if Bbool:
                    GroupB.append(B[ia])

        if count > 1:
            log.warning("%s in %s groups",a,count)

    if Bbool:
        return Grouplis , GrouplisB
    else:
        return Grouplis
    
def uniq_set_list(setlis,frozen=True):
    """ uniqify a list of sets"""
    if not frozen:
        return [set(e) for e in list(set([frozenset(e) for e in setlis]))]
    else:
        return [frozenset(e) for e in list(set([frozenset(e) for e in setlis]))]
    
    
def sort_binom_list(X,Y,array_out=False):
    """
    sort Y according to X
    and sort X
    """
    if len(X) != len(Y):
        log.warning('len(X) != len(Y) !!!')

    Ynew = [y for (x,y) in sorted(zip(X,Y))]
    Xnew = sorted(X)
    if not array_out:
        return Xnew,Ynew
    else:
        return np.array(Xnew),np.array(Ynew)

def sort_multinom_list(X,*Y):
    """
    sort Y according to X
    and sort X
    """
    Ynew_stk = []
    for YY in Y:
        Ynew = [y for (x,y) in sorted(zip(X,YY))]
        Ynew_stk.append(Ynew)
    Xnew = sorted(X)
    fintup = (np.array(Xnew),) + tuple(Ynew_stk)
    return fintup


def most_common(lst):
    """
    http://stackoverflow.com/questions/1518522/python-most-common-element-in-a-list
    """
    if type(lst) is not list:
        lst = list(lst)
    return max(set(lst), key=lst.count)

def median_improved(L):
    """
    manage the case where len(L) is even
    in this case, doesn't return the mean of the 2 medians, but the nearest value
    """
    L = np.array(L)
    if len(L) % 2 == 0:
        return find_nearest(L,np.median(L))[0]
    else:
        return np.median(L)

def trio_lists_2_tab(Xlis , Ylis , Vlis):
    """
    From a trio of lists
    Xlis => parameter you want in columns
    Ylis => parameter you want in rows
    Vlis => data(X,Y)

    make a tab compatible with mabular module
    tabulate(Finalis, headers="firstrow")

    NB : very dirty research of the elements ...
    """


    Xlis_uniq = sorted(np.unique(Xlis))
    Ylis_uniq = sorted(np.unique(Ylis))

    Finalis = []
    Finalis.append(Xlis_uniq)

    for y in Ylis_uniq:
        curlis = [y]
        Finalis.append(curlis)
        for x in Xlis_uniq:
            for xx , yy , vv in zip(Xlis , Ylis , Vlis):
                if xx == x and yy == y:
                    curlis.append(vv)
    return Finalis


def minmax(L):
    return np.min(L) , np.max(L)

def middle(Lin):
    Lout = []
    for i in range(len(Lin)-1):
        Lout.append((Lin[i] + Lin[i+1]) / 2)
    return Lout

def second_smallest(numbers):
    m1, m2 = float('inf'), float('inf')
    for x in numbers:
        if x <= m1:
            m1, m2 = x, m1
        elif x < m2:
            m2 = x
    return m2


def find_index_multi_occurences(L,elt):
    return [i for i, x in enumerate(L) if x == elt]


def find_surrounding(L,v):
    """
    find the surrounding values

    Parameters
    ----------
    L : Iterable
        Input list/array.

    Returns
    -------
    surrounding_values : tuple
        surounding values.
    surrounding_index : tuple
        surounding indexes.
        
        
    """
    a = np.array(L)
    b = v
    
    surrounding_values = np.sort([b + i for i in sorted(np.subtract(a,b),key=lambda x: abs(x))[:2]])
    surrounding_index = (np.where(surrounding_values[0] == a)[0][0],np.where(surrounding_values[1] == a)[0][0])
    
    return tuple(surrounding_values),surrounding_index

##################
### LIST IT FCTs
##################

def chunkIt(seq, num):
    """ make num sublists of a list """
    #http://stackoverflow.com/questions/2130016/splitting-a-list-of-arbitrary-size-into-only-roughly-n-equal-parts
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out

def sliceIt(seq,num):
    """ make sublist of num elts of a list """
    # http://stackoverflow.com/questions/4501636/creating-sublists
    return [seq[i:i+num] for i in range(0, len(seq), num)]


def sublistsIt(seq,lenofsublis_lis,output_array=False):
    """
    make sublists of seq , accoding to len of the sublist in lenofsublis_lis
    ex lenofsublis_lis = [2,3,4,2]

    if output_array : output list of array , else list of list

    (fct perso)
    """

    if np.sum(lenofsublis_lis) != len(seq):
        raise Exception('sublistsIt : sum(lenofsublis_lis) != len(seq) ')
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
    
    
def identical_consecutive_eltsIt(Lin):
    Lout_big = []
    Linter   = [Lin[0]]
    Lout_big.append(Linter)
            
    for e in Lin[1:]:
        if (e == Linter[-1]):
            Linter.append(e)
        else:
            Linter = [e]
            Lout_big.append(Linter)
                
    return Lout_big
                



def find_nearest(listin,value):
    """
    Returns:
        value of the nearest , index of the nearest
    """
    array = np.array(listin)
    idx = (np.abs(array-value)).argmin()

    return listin[idx] , idx


def find_interval_bound(listin,val,outindexes = True):
    """
    trouve les bornes d'un intervalle
    (listin est supposé triée)
    """
    i = bisect.bisect(listin,val)
    if outindexes:
        return i-1,i
    else:
        return listin[i-1] , listin[i]

def occurence(L,tolerence=None,pretty_output=False):
    """
    Input :
        tolerence : tolerence to find close elements of L
                    if no tolerence is given then a set() is used
    Returns :
        if pretty_output = False :
            return a list of 2-tuples : (element of the list, number of occurence of this element in the list)
        if pretty_output 0 True :
            
    Nota : pretty_output is implemented because the first mode is not really useful (180612)
           the equal test is also replaced by is close
    """
    if tolerence:
        Lset = uniquetol2(L , tol=tolerence)
    else:
        Lset = set(L)
        
    L = np.array(L)
    outlisoccur = []
    for l in Lset:
        outlisoccur.append((l,np.sum(L == l)))
        
    if pretty_output:
        Mtmp = np.vstack(outlisoccur)
        Vals   = Mtmp[:,0]
        Occurs = Mtmp[:,1] 
        Occurs , Vals = sort_binom_list(Occurs , Vals)
    
        output = (Occurs , Vals)
    else:
        output = outlisoccur
       
        
    return output
def decimateIt(listin,n):
    """
    n so as np.mod(i,n) == 0
    """
    outlist = []
    for i in range(len(listin)):
        if np.mod(i,n) == 0:
            outlist.append(listin[i])
    return outlist


def consecutive_groupIt(data,only_start_end=False):
    """
    Identify groups of continuous numbers in a list

    Useful for time period with a prealable conversion to MJD

    Source :
    https://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
    """
    from operator import itemgetter
    from itertools import groupby
    groups_stk = []
    for k, g in groupby(enumerate(data), lambda A : A[0] - A[1]):
        groups_stk.append(list(map(itemgetter(1), g)))

    if only_start_end:
        return [(g[0],g[-1]) for g in groups_stk]
    else:
        return groups_stk


def identical_groupIt(data):
    """
    Source :
        https://stackoverflow.com/questions/30293071/python-find-same-values-in-a-list-and-group-together-a-new-list
    """
    return [list(j) for i, j in itertools.groupby(data)]


def get_interval(start, end, delta):
    # after http://stackoverflow.com/questions/10688006/generate-a-list-of-datetimes-between-an-interval-in-python
    outlist =[]
    curr = start
    while curr < end:
        outlist.append(curr)
        curr += delta
    return outlist


def duplicates_finder(seq):
    """
    from
    http://stackoverflow.com/questions/9835762/find-and-list-duplicates-in-python-list
    moooeeeep solution
    """
    seen = set()
    seen_add = seen.add
    # adds all elements it doesn't know yet to seen and all other to seen_twice
    seen_twice = set( x for x in seq if x in seen or seen_add(x) )
    # turn the set into a list (as requested)
    return list( seen_twice )

def sort_table(table, col):
    """ sort a table by multiple columns
        table: a list of lists (or tuple of tuples) where each inner list
               represents a row
        cols: specifying the column numbers to sort by
    """
    refcol = list(table[col])
    outtable = []
    for col in table:
        outtable.append([x for (y,x) in sorted(zip(refcol,col))])
    return outtable


def dicofdic(mat,names):

    # a partir d'un matrice de N x N
    # et d'une liste de N noms
    # Fabrique un dictionnaire 2D dic[nom1][nom2] = M[n1,n2]
    # http://stackoverflow.com/questions/13326042/2d-dictionary-with-multiple-keys-per-value

    d2_dict = collections.defaultdict(dict)

    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            d2_dict[names[i]][names[j]] = mat[i,j]

    return d2_dict


def find_regex_in_list(regex,L,only_first_occurence=False,
                       line_number=False):
    Lout = []
    for i,e in enumerate(L):
        if re.search(regex,e):
            if not line_number:
                found = e
            else:
                found = (i,e)
                
            if only_first_occurence:
                return found
            else:    
                Lout.append(found)
    return Lout
            
