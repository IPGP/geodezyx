# -*- coding: utf-8 -*-
"""
Created on Thu Jun 19 16:28:08 2014

@author: psakicki

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics

Copyright (C) 2019 Pierre Sakic (GFZ, pierre.sakic@gfz-postdam.de)
GitHub repository :
https://github.com/PierreS1/GeodeZYX-Toolbox-Lite

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
"""
import os
import re
import numpy as np
from collections import defaultdict
import operator
from tempfile import mkstemp
from shutil import move, copy2
# http://stackoverflow.com/questions/123198/how-do-i-copy-a-file-in-python
from os import remove, close
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib
import sys
import subprocess
import shutil
import bisect
from collections import Counter
import pickle
import inspect
import uuid
import time
import tempfile
import scipy
from natsort import natsorted, ns
import fnmatch
import itertools
import pandas as pd

sys.dont_write_bytecode = True

#def clear():
#    '''Clears the shell of the spider application.
#    Use either clear() or cls()'''
#    os.system('cls')
#    return None
#cls = clear

#def geodezyx_checker():
    


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
    if not consider_str_as_iterable and type(inp) is str:
        return False

    try:
        iter(inp)
    except TypeError:
        out = True
    else:
        out = False
    return out

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

def clear_all():
    '''Clears all the variables from the workspace of the spyder
    application.'''
    #cls()
    gl = globals().copy()
    for var in gl:
        if var[0] == '_': continue
        if 'func' in str(globals()[var]): continue
        if 'module' in str(globals()[var]): continue

        del globals()[var]

def spyder_run_check():
    """ Check if the code is run inside Spyder IDE"""
    if any('SPYDER' in name for name in os.environ):
        return True
    else:
        return False

def indice_printer(i,print_every = 10,text_before=''):
    if np.mod(i,print_every) == 0:
        print(text_before ,  i)
    return None

def regex2filelist(dossier,regex,outtype='file'):

    ''' a partir d'un chemin de dossier et d'une regex, donne les éléments
    du dossier correspondant à la regex

    OUTTYPE :
    file : juste les fichiers'''

    templist = []

    for filename in os.listdir(dossier):
        if re.compile(regex).search(filename):
            templist.append(os.path.join(dossier,filename))

    #  Contient fichiers + dossiers
    #  => eliminer les dossiers

    if outtype == 'file':
        outlist = [ f for f in templist if os.path.isfile(f) ]
    else:
        outlist = templist
    outlist.sort()

    return outlist

def check_regex(filein,regex):

    ''' verfie si un fichier contient une regex
        retourne un booleen '''

    outbool = False

    for line in open(filein):
        if re.compile(regex).search(line):
            outbool = True
            break

    return outbool


def patterns_in_string_checker(string,*patterns):
    """
    from
    http://stackoverflow.com/questions/3389574/check-if-multiple-strings-exist-in-another-string
    """
    print(patterns)
    L = [x in string for x in patterns]
    return bool(any( L ))


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
            print("WARN : " , a , " in " , count , " groups")

    if Bbool:
        return Grouplis , GrouplisB
    else:
        return Grouplis


def fileprint(output,outfile):
    print(output)
    with open(outfile, "a") as f:
        f.write("{}\n".format(output))
    return None

def write_in_file(string_to_write,outdir,outname,ext='.txt'):
    outpath = os.path.join(outdir,outname + ext)
    F = open(outpath,'w+')
    try:    
        F.write(string_to_write.encode('utf8'))
        # astuce de http://stackoverflow.com/questions/6048085/writing-unicode-text-to-a-text-file
    except TypeError:
        print("INFO : write_in_file : alternative reading following a TypeError")
        F.write(string_to_write)
        
    F.close()
    return outpath

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

def trunc(f, n):
    ''' troncature un flottant f a la nieme decimale '''
    '''Truncates/pads a float f to n decimal places without rounding'''
    slen = len('%.*f' % (n, f))
    return float(str(f)[:slen])

def multidot(tupin):

    out = np.eye(tupin[0].shape[0])

    for e in tupin:
        out = np.dot(out,e)

    return out

def mdot(*args):
    ret = args[0]
    for a in args[1:]:
        ret = np.dot(ret,a)
    return ret

def mdotr(*args):
    ret = args[-1]
    for a in reversed(args[:-1]):
        ret = np.dot(a,ret)
    return ret

def dicofdic(mat,names):

    # a partir d'un matrice de N x N
    # et d'une liste de N noms
    # Fabrique un dictionnaire 2D dic[nom1][nom2] = M[n1,n2]
    # http://stackoverflow.com/questions/13326042/2d-dictionary-with-multiple-keys-per-value

    d2_dict = defaultdict(dict)

    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            d2_dict[names[i]][names[j]] = mat[i,j]

    return d2_dict


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


def replace(file_path, pattern, subst):
    """ from http://stackoverflow.com/questions/39086/search-and-replace-a-line-in-a-file-in-python """
    #Create temp file
    fh, abs_path = mkstemp()
    new_file = open(abs_path,'w')
    old_file = open(file_path)
    for line in old_file:
        new_file.write(line.replace(pattern, subst))
    #close temp file
    new_file.close()
    close(fh)
    old_file.close()
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)

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

def find_index_multi_occurences(L,elt):
    return [i for i, x in enumerate(L) if x == elt]

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

#################
### SHELL LIKE FCTS
#################

def tail(filename, count=1, offset=1024):
    """
    A more efficent way of getting the last few lines of a file.
    Depending on the length of your lines, you will want to modify offset
    to get better performance.
    """
    f_size = os.stat(filename).st_size
    if f_size == 0:
        return []
    with open(filename, 'r') as f:
        if f_size <= offset:
            offset = int(f_size / 2)
        while True:
            seek_to = min(f_size - offset, 0)
            f.seek(seek_to)
            lines = f.readlines()
            # Empty file
            if seek_to <= 0 and len(lines) == 0:
                return []
            # count is larger than lines in file
            if seek_to == 0 and len(lines) < count:
                return lines
            # Standard case
            if len(lines) >= (count + 1):
                return lines[count * -1:]

def head(filename, count=1):
    """
    This one is fairly trivial to implement but it is here for completeness.
    """
    with open(filename, 'r') as f:
        lines = [f.readline() for line in range(1, count+1)]
        return list(filter(len, lines))

def grep(file_in,search_string,only_first_occur=False,
         invert=False,regex=False,line_number=False,col=(None,None)):
    """
    if nothing is found returns a empty string : ""
    (and NOT a singleton list with an empty string inside)
    
    Args :
        col : Define the columns where the grep is executed 
              (Not delimited columns , one character = a new column)
              from 1st col. / until last col. : use None as index

    search_string can be a list (150721 change)
    """
    if type(search_string) is str:
        search_string = [search_string]

    matching_line_list = []
    line_number_list   = []
    trigger = False
    for iline , line in enumerate(open(file_in,encoding = "ISO-8859-1")):
        trigger = False
        for seastr in search_string:
            if regex:
                if re.search(seastr,line[col[0]:col[1]]):
                    trigger = True
            else:
                if seastr in line[col[0]:col[1]]:
                    trigger = True
        if invert:
            trigger = not trigger
        if trigger:
            matching_line_list.append(line)
            line_number_list.append(iline)
            if only_first_occur:
                break
    if  line_number and len(line_number_list) == 1:
        return line_number_list[0] , matching_line_list[0]
    elif line_number:
        return line_number_list    , matching_line_list
    elif len(matching_line_list) == 1:
        return matching_line_list[0]
    elif len(matching_line_list) == 0:
        return ''
    elif only_first_occur:
        return matching_line_list[0]
    else:
        return matching_line_list

def egrep_big_string(regex,bigstring,only_first_occur=False):
    """
    perform a regex grep on a big string sepatated with \n

    NB : must be improved with regular pattern matching, wo regex
    """

    matching_line_list = []

    for l in bigstring.split("\n"):
        result = re.search(regex, l)
        if result:
            matching_line_list.append(l)
            if only_first_occur:
                break

    if len(matching_line_list) == 1:
        return matching_line_list[0]
    elif len(matching_line_list) == 0:
        return ''
    elif only_first_occur:
        return matching_line_list[0]
    else:
        return matching_line_list

def grep_boolean(file_in,search_string):
    for line in open(file_in):
        if search_string in line:
            return True
    return False

def regex_OR_from_list(listin):
    return "(" + join_improved("|" , *listin) +  ")"

def cat(outfilename, *infilenames):
    """
    Is for concatenating files ...
    For just a print, use cat_print !
    http://stackoverflow.com/questions/11532980/reproduce-the-unix-cat-command-in-python
        kindall response
    """
    with open(outfilename, 'w') as outfile:
        for infilename in infilenames:
            with open(infilename , 'r+') as infile:
                for line in infile:
                    if line.strip():
                        outfile.write(line)
    return outfilename


def cat_remove_header(infilepath,outfilepath,header='',
                      header_included = False):

    bool_out = False
    F = open(infilepath,'r+')

    with open(outfilepath, 'w') as outfile:
        for line in F:
            if header in line:
                bool_out = True
                if not header_included:
                    continue
            if bool_out:
                outfile.write(line)

    return outfilepath

def cat_print(inpfile):
    fil = open(inpfile)
    for l in fil:
        if l[-1] == '\n':
            print(l[:-1])
        else:
            print(l)
    return None


def empty_file_check(fpath):
    """
    Check if a file is empty or not. 
    NB : Also check if the file exists
    
    Parameters
    ----------
    fpath : str
        the file path
                
    Returns
    -------  
    True : 
        if the file is empty
    
    False : 
        if the file is not empty
         
    Source
    ------
    http://stackoverflow.com/questions/2507808/python-how-to-check-file-empty-or-not
    """
    return not (os.path.isfile(fpath) and os.path.getsize(fpath) > 0)


def find_recursive(parent_folder , pattern, 
                   sort_results = True, case_sensitive = True):
    """
    Find files in a folder and his sub-folders in a recursive way

    Parameters
    ----------
    parent_folder : str
        the parent folder path

    pattern : str
        the researched files pattern name (can manage wildcard or regex)
        * wildcard (only * and ?) for case_sensitive = False
        * regex for case_sensitive = True
        
    sort_results : bool
        Sort results
        
    case_sensitive : bool
        Case sensitve or not
                
    Returns
    -------
    matches : list
        Found files        

    Source
    ------
    https://stackoverflow.com/questions/2186525/use-a-glob-to-find-files-recursively-in-python
    
    https://stackoverflow.com/questions/15652594/how-to-find-files-with-specific-case-insensitive-extension-names-in-python
    (for the case unsensitive case)
    """
    matches = []
    if case_sensitive:
        for root, dirnames, filenames in os.walk(parent_folder):
            for filename in fnmatch.filter(filenames, pattern):
                matches.append(os.path.join(root, filename))
    else:
        for root, dirnames, filenames in os.walk(parent_folder):
            for filename in filenames:  
                try:
                    bool_match = re.search(pattern, filename, re.IGNORECASE)
                except Exception as e:
                    print("ERR : if case_sensitive = True, pattern have to be a REGEX (and not only a simple wildcard)")
                    raise Exception
                    
                if bool_match:
                    matches.append(os.path.join(root, filename))

    if sort_results:
        matches = sorted(matches)

    return matches

def glob_smart(dir_path,file_pattern=None,verbose=True):
    if file_pattern:
        dir_path_ok = os.join.path(dir_path,file_pattern)
    else:
        dir_path_ok = dir_path
        
    outlist = glob.glob(dir_path_ok)
    
    if verbose:
        if not outlist:
            print("WARN : no file(s) found as" , dir_path_ok)
        else:
            print("INFO :" , len(outlist),"file(s) found as", dir_path_ok)
            
    return outlist


def insert_lines_in_file(file_path,text_values,lines_ids):
    
    if not is_iterable(text_values):
        text_values = [text_values]
    
    if not is_iterable(lines_ids):
        lines_ids = [lines_ids]
    
    f = open(file_path, "r")
    contents = f.readlines()
    f.close()

    for txt , lin in zip(text_values,lines_ids):
        contents.insert(lin, txt)       
    
    f = open(file_path, "w")
    contents = "".join(contents)
    f.write(contents)
    f.close()

    return file_path

def insert_str_in_file_if_line_contains(file_path,str_to_insert,
                                    line_pattern_tup,
                                    position=None,
                                    only_first_occur=False):
    """
    NB : position is not implemented
    """
    f = open(file_path, "r")
    contents = f.readlines()
    
    f.close()

    pattern_found = False
    
    for ilin , lin in enumerate(contents):
        for pattern in line_pattern_tup:
            if re.search(pattern, lin):
                contents[ilin] = str_to_insert + lin
                pattern_found = True
                if only_first_occur:
                    break
                
        if only_first_occur and pattern_found:
            print("break at line " , ilin)
            break
                
    f = open(file_path, "w")
    contents = "".join(contents)
    f.write(contents)
    f.close()

    return file_path
    


#########################
### PLOT
#########################

def color_list(L , colormap='jet'):
    cm     = plt.get_cmap(colormap)
    NCOL   = len(np.unique(L))
    colist = [cm(1.*i/NCOL) for i in range(NCOL)]
    return colist

def symbols_list(L):

    Lsym = ["o",
    "v",
    "^",
    "<",
    ">",
    ".",
    ",",
    "1",
    "2",
    "3",
    "4",
    "8",
    "s",
    "p",
    "P",
    "*",
    "h",
    "H",
    "+",
    "x",
    "X",
    "D",
    "d",
    "|",
    "_"]

    return Lsym[:len(L)]



def colors_from_colormap_getter(ncolors , colormap = 'viridis'):
    import matplotlib.pyplot as plt
    cm = plt.get_cmap(colormap)
    return [cm(1.*i/ncolors) for i in range(ncolors)]

def ylim_easy(Lin,delta = .1, min_null_if_neg = False):
    minn = np.min(Lin)
    maxx = np.max(Lin)
    rangee = np.abs(maxx - minn)
    if  min_null_if_neg and (minn - delta * rangee) < 0. :
        return (0. , maxx + delta * rangee)
    else:
        return (minn - delta * rangee , maxx + delta * rangee)

def get_figure(figin = 0):
    # Un autre exemple comme défini dans
    # export_ts_figure_pdf
    #    if type(fig) is int:
    #        f = plt.figure(fig)
    #    elif type(fig) is Figure:
    #        f = fig
    if isinstance(figin,matplotlib.figure.Figure):
        figout = figin
    elif figin == 0:
        figout = plt.figure()
    else:
        figout = plt.figure(figin)
    # be sure the fig have a axe (necessary ?)
    if figout.get_axes() == []:
        figout.add_subplot(111)
    plt.figure(figout.number)
    return figout


def figure_saver(figobjt_in , outdir , outname , outtype = '.png' , formt = 'a4' ):
    if not is_iterable(outtype):
         outtype = (outtype,) 
         
    outpath_stk = []
    for outtype_iter in outtype:
        outpath = os.path.join(outdir,outname+outtype_iter)
    #    if formt == 'a4':
    #    elif
        figobjt_in.savefig(outpath)
        outpath_stk.append(outpath)
    if len(outpath_stk) == 1:
        outpath_stk = outpath_stk[0]

    return outpath_stk


def axis_data_coords_sys_transform(axis_obj_in,xin,yin,inverse=False):
    """ inverse = False : Axis => Data
                = True  : Data => Axis
    """
    xlim = axis_obj_in.get_xlim()
    ylim = axis_obj_in.get_ylim()

    xdelta = xlim[1] - xlim[0]
    ydelta = ylim[1] - ylim[0]
    if not inverse:
        xout =  xlim[0] + xin * xdelta
        yout =  ylim[0] + yin * ydelta
    else:
        xdelta2 = xin - xlim[0]
        ydelta2 = yin - ylim[0]
        xout = xdelta2 / xdelta
        yout = ydelta2 / ydelta
    return xout,yout

def id2val(value_lis,id_lis,idin):
    """ from a value list and a id pointer list
        return the good val from the good id
        replace dico bc. set is not supproted as key"""
    return value_lis[id_lis.index(idin)]

################
### "GETTERs" : wrappers to get easily some fcts
################

def get_timestamp(outstring = True):
    """frontend to get easily a timestamp"""
    if outstring:
        return dt.datetime.now().strftime('%Y%m%d_%H%M%S')
    else:
        return dt.datetime.now()

def get_function_name():
    return inspect.stack()[0][3]

def get_computer_name():
    import platform
    return platform.node()


def get_username():
    import pwd
    return pwd.getpwuid( os.getuid() )[ 0 ]
  

def vectorialize(array_in):
    """
    redondant avec .flat ???
    """
    vector_out = array_in.reshape((array_in.size,))
    return vector_out

def get_specific_locals(prefix):
    """ get locals params with 'prefix' in the name
        can actually be a suffix """
    # MARCHE PAS EN L'ETAT
    loctemp = dict(globals())
    print(loctemp)
    outlis = []
    for k,v in loctemp.items():
        if prefix in k:
            outlis.append(str(k[4:]) + ' : ' + str(v))
    outlis.sort()
    return outlis

class Tee(object):
    # based on
    # http://stackoverflow.com/questions/11325019/output-on-the-console-and-file-using-python
    # Secondary links
    # http://stackoverflow.com/questions/616645/how-do-i-duplicate-sys-stdout-to-a-log-file-in-python
    # http://stackoverflow.com/questions/2996887/how-to-replicate-tee-behavior-in-python-when-using-subprocess
    def __init__(self, *files):
        self.original = sys.stdout
        self.files = [sys.stdout] + list(files)
        sys.stdout = self

    def write(self, obj):
        for f in self.files:
            f.write(obj)
    def flush(self):
        for f in self.files:
            f.flush()
    def stop(self):
        for fil in self.files:
            if fil != self.original:
                fil.close()
        sys.stdout = self.original
        self.files = [self.original]

    def pause(self):
        sys.stdout = self.original
        self.files_saved = list(self.files)
        self.files = [self.original]

    def restart(self):
        self.files = self.files_saved
        sys.stdout = self

def Tee_frontend(pathin,prefix,suffix='',ext='log',print_timestamp=True):
    #if suffix != '':
    #    suffix = suffix + '_'
    if print_timestamp:
        ts = '_' + get_timestamp()
    else:
        ts = ''
    f = open(os.path.join(pathin,prefix +'_'+suffix+ts+'.'+ext), 'w')
    f_tee = Tee(f)
    #sys.stdout = f_tee
    return f_tee


def uncompress(pathin,dirout = '', opts='-f'):
    if not os.path.isfile(pathin):
        print('ERR : uncompress : ', pathin , ' dont exist !!!')
        return None
    komand = 'uncompress ' + opts + ' ' + pathin
    subprocess.call([komand], shell=True)
    pathout_temp = '.'.join(pathin.split('.')[:-1])

    if dirout == '':
        pathout = pathout_temp
    else:
        pathout = os.path.join(dirout,os.path.basename(pathout_temp))
        shutil.move(pathout_temp,pathout)
    return pathout


def most_common(lst):
    #http://stackoverflow.com/questions/1518522/python-most-common-element-in-a-list
    if type(lst) is not list:
        lst = list(lst)
    return max(set(lst), key=lst.count)

def line_count(filein):
    num_lines = sum(1 for line in open(filein))
    return num_lines


def str2float_smart(str_in):
    try:
        out = float(str_in)
    except ValueError:
        out = np.nan
    return out


def str2int_smart(str_in):
    try:
        out = int(str_in)
    except ValueError:
        out = np.nan
    return out

def array_from_lists(*listsin):
    """ fonction pour arreter de galerer avec les conversions
    de lists => matrices """
    out = np.vstack(listsin)
    out = out.T
    return out

def stringizer(tupin,separ=' ',eol=True):
    """
    transform elts of a tuple in a string line, ready for write in a file
    """
    if eol:
        end = '\n'
    else:
        end = ''
    lineout =  ' '.join([str(e) for e in tupin]) + end
    return lineout

def eval_a_dict(dictin,where,verbose=True):
    """
    where is most of time globals()
    WARN : doesnt work in a function !!!
    use instead :
    for k,v in booldic.items():
        globals()[k] = v
        locals()[k] = v
    """
    for k,v in dictin.items():
        where[k] = v
        if verbose:
            print('INFO : eval_a_dict : the 2 values must be equals ', where[k] , v , ", variable",k )
    return None

def boolean_dict(list_of_keywords):
    outdict = dict()
    for k in list_of_keywords:
        outdict[k] = (True , False)
    return outdict

 
def line_in_file_checker(file_path,string):
    #Open the file
    fobj = open(file_path)
    
    trigger = False
    for line in fobj:
        if string in line:
            trigger = True
            break
        
    fobj.close()
    return trigger
        

def extract_text_between_elements(file_path,elt_start,elt_end):
    """
    source :
        https://stackoverflow.com/questions/9222106/how-to-extract-information-between-two-unique-words-in-a-large-text-file
    """

    f = open(file_path,'r')
    text = f.read()
    text_extract = text.split(elt_start,1)[-1].split(elt_end)[0]

    return text_extract

def extract_text_between_elements_2(file_path , elt_start , elt_end,
                                           return_string = False,
                                           nth_occur_elt_start=0,
                                           nth_occur_elt_end=0):
    """
    This function is based on REGEX (elt_start , elt_end are REGEX)
    and can manage several blocks in the same file

    return_string = True  : returns a string of the matched lines
    return_string = False : returns a list of the matched lines
    
    NB : in SINEX context, with "+MARKER", use backslash i.e.
         "\+MARKER"    
         
    NB2 : think about StingIO for a Pandas DataFrame Handeling
    https://docs.python.org/2/library/stringio.html
    """
    
    try:
        F = open(file_path,"r",encoding = "ISO-8859-1")
    except:
        F = file_path
    
    out_lines_list = []
    
    trigger = False
    last_triggered_line = False
    
    for line in F:

        if re.search(elt_start , line) and nth_occur_elt_start == 0:
            trigger = True
        elif re.search(elt_start , line) and nth_occur_elt_start > 0:
            nth_occur_elt_start = nth_occur_elt_start - 1
        
        if re.search(elt_end , line) and nth_occur_elt_end == 0:
            last_triggered_line = True
        elif trigger and re.search(elt_end , line) and nth_occur_elt_end > 0:
            nth_occur_elt_end = nth_occur_elt_end - 1
            
        if trigger:
            out_lines_list.append(line)
            
        if last_triggered_line:
            trigger = False
            last_triggered_line = False

    if return_string:
        return "".join(out_lines_list)
    else:
        return out_lines_list
    
    
def str_2_float_line(line , sep=" ",out_type=float):
    """
    convert a line of number (in str)
    to a list of float (or other out_type)
    """
    fields = line.strip().split(sep)
    fields = [e for e in fields if len(e) > 0]
    try:
        return [out_type(e) for e in fields]
    except:
        return fields
    
########################
### WRAPPER OF I/O FCTs
########################

def pickle_saver(datain , outdir = None , outname = None , ext='.pik' ,
                 timestamp = False,full_path=None):
    """
    if full_path is given, override outdir and outname
    RETURN :
        outpath : the output path of the pickle
    """
    if not full_path and not outdir and not outname:
        print('ERR : pickle_saver : not full_path and not outdir and not outname are given !')
        raise Exception

    if timestamp:
        ts = get_timestamp() + '_'
    else:
        ts = ''
    if full_path:
        outpath = full_path
    else:
        outpath = os.path.join(outdir,ts + outname + ext)


    # For Python 2
    # pickle.dump( datain , open( outpath , "w+" ) )
    # For Python 3
    pickle.dump( datain , open( outpath , "wb" ) )
    return outpath

def pickle_loader(pathin):
    try:
        #outdata = pickle.load( open( pathin , "r" ) )
        outdata = pickle.load( open( pathin ,'rb' )  , encoding='latin1')
    except UnicodeDecodeError:
        print('INFO : pickle loader : alternative reading following a UnicodeDecodeError')
        #outdata = pickle.load( open( pathin ,'rb' )  , encoding='latin1')
        outdata = pickle.load( open( pathin , "r" ) )
    except TypeError:
        print('INFO : pickle loader : alternative reading following a TypeError (is it a Py2 pickle ?)')
        # source :
        # https://stackoverflow.com/questions/28218466/unpickling-a-python-2-object-with-python-3
        outdata = pickle.load( open( pathin ,'rb' )  , encoding='latin1')

    return outdata

def memmap_from_array(arrin):
    nam = str(np.random.randint(99999)) + '.mmp.tmp'
    path = os.path.join(tempfile.mkdtemp(), nam)
    mmp = np.memmap(path, dtype='float64', mode='w+',
                    shape=arrin.shape)
    mmp[:] = arrin[:]
    return mmp

def mmpa(arrin):
    # just a light wrapper
    return memmap_from_array(arrin)

def read_mat_file(pathin,full=False):
    """ low level reader of a MATLAB mat file """
    import scipy.io as sio
    import os
    fullmat = sio.loadmat(pathin)
    data = fullmat[os.path.basename(pathin).split('.')[0]]
    if full:
        return data , fullmat
    else:
        return data

def save_obj_as_file(objin,pathin,prefix,ext='.exp',suffix=''):
    """
    OLD proto-version of pickle saver DISCONTINUED
    """
    if suffix != '':
        suffix = suffix + '_'
    ts = get_timestamp()
    fullpath = os.path.join(pathin,prefix +'_'+suffix+ts+ext)


    # Python 2
    #pickle.dump( objin, open( fullpath , "w+" ) )
    # Python 3
    pickle.dump( objin , open( fullpath , "wb" ) )

    return fullpath

def uniq_set_list(setlis,frozen=True):
    """ uniqify a list of sets"""
    if not frozen:
        return [set(e) for e in list(set([frozenset(e) for e in setlis]))]
    else:
        return [frozenset(e) for e in list(set([frozenset(e) for e in setlis]))]

def create_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory

def remove_dir(directory):
    if os.path.exists(directory):
        shutil.rmtree(directory)
    else:
        print('WARN : remove_dir : ' , directory , ' dont exists ... skip ...')
    return None

def walk_dir(parent_dir):
    """
    from a main parent_dir
    returns files_list & dirs_list all the files and all the dirs in the
    parent_dir

    https://www.tutorialspoint.com/python/os_walk.htm
    """
    files_list , dirs_list = [] , []
    for root, dirs, files in os.walk(parent_dir, topdown=False):
        for name in files:
            files_list.append(os.path.join(root, name))
        for name in dirs:
            dirs_list.append(os.path.join(root, name))

    return files_list , dirs_list



def save_array_fast(arrin,outname='',
                    outdir='/home/psakicki/aaa_FOURBI',
                    txt=True):
    if outname == '':
        outname = str(uuid.uuid4())[:8]
    if scipy.sparse.issparse(arrin):
        arrin = arrin.toarray()
    outpath = os.path.join(outdir,outname)
    print(outpath)
    if txt:
        np.savetxt(outpath,arrin)
    else:
        np.save(outpath,arrin)
    print(outpath)
    return outpath

def join_improved(strseparat,*varsin):
    return strseparat.join([str(v) for v in varsin])

def split_improved(strin,sep_left,sep_right):
    return strin.split(sep_left)[1].split(sep_right)[0]

def diagonalize(x,n=10):
    if is_iterable(x):
        x = list(x) * n
        x = np.expand_dims(np.array(x),1)
        return np.diag(x[:,0])
    else:
        return np.diag((np.ones((n,n)) * x)[:,0])


def read_comments(filein,comment='#'):
    outlist = []
    for l in open(filein):
        if l[0] == comment:
            outlist.append(l[1:-1])
    return outlist

def sort_binom_list(X,Y,array_out=False):
    """
    sort Y according to X
    and sort X
    """
    if len(X) != len(Y):
        print('WARN : sort_binom_list : len(X) != len(Y) !!!')

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
    

            
            
def dicts_merge(*dict_args):
    '''
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.

    WARN : first values will be erased if the same key is present in following dicts !!!

    http://stackoverflow.com/questions/38987/how-can-i-merge-two-python-dictionaries-in-a-single-expression
    '''
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


def dicts_of_list_merge_mono(dol1, dol2):
    """
    https://stackoverflow.com/questions/1495510/combining-dictionaries-of-lists-in-python
    """
    keys = set(dol1).union(dol2)
    no = []
    return dict((k, dol1.get(k, no) + dol2.get(k, no)) for k in keys)


def dicts_of_list_merge(*dict_args):
    result = dict()
    for dictionary in dict_args:
        result = dicts_of_list_merge_mono(result,dictionary)
    return result


def dic_key_for_vals_list_finder(dic_in , value_in):
    """
    dic_in is a dict like :
        dic_in[key1] = [val1a , val1b]
        dic_in[key2] = [val2a , val2b , val2c]
        
    E.g. if value_in = val2b then the function returns key2
    
    NB : the function returns the first value found, then
         dic_in has to be injective !!
    """
    for k , v in dic_in.items():
        if value_in in v:
            return k
    
    print("WARN : no key found for value",value_in,"... None returned")
    return None


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

def alphabet(num=None):
    if not num:
        import string
        return list(string.ascii_lowercase)
    else:
        return alphabet()[num]

def dday():
    D = (dt.datetime(2016,10,14) - dt.datetime.now()).days
    print('J -',  D , 'avant la quille')
    return D

def Aformat(A,landscape=True):
    LA = [(841 , 1189),
    (594 , 841),
    (420 , 594),
    (297 , 420),
    (210 , 297),
    (148 , 210),
    (105 , 148),
    (74  , 105),
    (52  , 74)]

    out = np.array(LA[A])/25.4
    if landscape:
        out = np.flipud(out)
    return tuple(out)

def timeout(func, args=(), kwargs={}, timeout_duration=1, default=None):
    '''
    This function will spwan a thread and run the given function
    using the args, kwargs and return the given default
    value if the timeout_duration is exceeded

    http://stackoverflow.com/questions/366682/how-to-limit-execution-time-of-a-function-call-in-python
    '''
    import threading
    class InterruptableThread(threading.Thread):
        def __init__(self):
            threading.Thread.__init__(self)
            self.result = default
        def run(self):
            try:
                self.result = func(*args, **kwargs)
            except:
                self.result = default
    start = time.time()
    it = InterruptableThread()
    it.start()
    it.join(timeout_duration)
    if it.isAlive():
        print('ERR : a timout is triggered !!!')
        stop = time.time() - start
        print("     runtime : " , stop, 'sec.')
        return it.result
    else:
        stop = time.time() - start
        print("INFO : runtime : " , stop, 'sec.')
        return it.result


def find_common_elts(*lists):
    sett = set(lists[0])

    for l in lists[1:]:
        sett  = sett.intersection(l)

    return np.array(sorted(list(sett)))


def uniq_and_sort(L,natural_sort=True):
    if natural_sort:
        return natsorted(list(set(L)))
    else:
        return sorted(list(set(L)))


def df_sel_val_in_col(DF,col_name,col_val):
    """
    Return a selected value of a column in a DataFrame
    i.e.
    return DF[DF[col_name] == col_val]
    """
    return DF[DF[col_name] == col_val]



def docstring_generic():
    """
    prints and returns an prototype generic docstring. Based on Numpy docstring
    convention
    
    Source
    ------
    https://numpydoc.readthedocs.io/en/latest/format.html
    
    """
    
    docstr_out = """    
    General description

    Parameters
    ----------
    param1 : float or int or str or dict or n-tuple or bool or list or numpy.array
        Description param1

    param2 : float or int or str or dict or n-tuple or bool or list or numpy.array
        Description param2
        
    param3 : float or int or str or dict or n-tuple or bool or list or numpy.array
        Description param3
                
    Returns
    -------
    out1 : float or int or str or dict or n-tuple or bool or list or numpy.array
        Description out1
    
    out2 : float or int or str or dict or n-tuple or bool or list or numpy.array
        Description out2
        
    Note
    ----
    Misc. Notes

    Source
    ------
    www.forum-source.com
    
    Examples
    --------
    >>> answer
    42    
    """
    
    print(docstr_out)
    
    return docstr_out

