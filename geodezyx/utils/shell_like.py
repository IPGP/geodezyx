#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.utils contains functions for Shell-like
 operations in Python. 

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
import fnmatch
import glob
import gzip
#### Import the logger
import logging
import os
import re
import shutil
import subprocess

#### geodeZYX modules
from geodezyx import utils

log = logging.getLogger(__name__)

##########  END IMPORT  ##########



#################
### SHELL LIKE FCTS
#################



def subprocess_frontend(cmd_in,
                        save_log=False,
                        log_dir=None,
                        log_name_out="out.log",
                        log_name_err="err.log",
                        logname_timestamp=False):
    
    
    now = utils.get_timestamp()
    
    process1 = subprocess.run([cmd_in],
                              capture_output=True,
                              text=True,
                              shell=True)
    
    process1_stdout = process1.stdout.strip()
    process1_stderr = process1.stderr.strip()
    
    if save_log:
        
        if not log_dir:
            log_dir = os.getcwd()
            
        if logname_timestamp:
            prefix = now + "_" 
        else:
            prefix = ""
        
        out_file = open(log_dir + "/" + prefix + log_name_out, "a+")
        err_file = open(log_dir + "/" + prefix + log_name_err, "a+")
    
        out_file.write(process1_stdout)
        err_file.write(process1_stderr)
    
        out_file.close()
        err_file.close()
    
    return process1,process1_stdout,process1_stderr



def tail(filename, count=1, offset=1024):
    """
    A more efficent way of getting the last few lines of a file.
    Depending on the length of your lines, you will want to modify offset
    to get better performance.
    """
    
    f_size = os.stat(filename).st_size
    if f_size == 0:
        return []
    with open(filename, 'r',errors='ignore') as f:
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
        lines = []
        for iline in range(1, count+1):
            try:
                lines.append(f.readline())
            except Exception as e:
                log.warn(e)
                continue
        return list(filter(len, lines))



def grep(file_in,search_string,only_first_occur=False,
         invert=False,regex=False,line_number=False,col=(None,None),
         force_list_output=False):
    """
    if nothing is found returns a empty string : ""
    (and NOT a singleton list with an empty string inside)
    
    Args :
        col : Define the columns where the grep is executed 
              (Not delimited columns , one character = a new column)
              from 1st col. / until last col. : use None as index
              
              force_list_output : if the output is an unique element, 
              it will be returned in a list anyway

    search_string can be a list (150721 change)
    """
    if type(search_string) is str:
        search_string = [search_string]

    matching_line_list = []
    line_number_list   = []
    trigger = False
    
    
    if file_in[-2:] in ("gz","GZ"): ### cand handle gziped files 
        File = gzip.open(file_in, "r+")
        Lines = [e.decode('utf-8') for e in File]
        File.close()
    else:
        File = open(file_in,encoding = "ISO-8859-1")
        Lines = File.readlines()
        File.close()
        
    for iline , line in enumerate(Lines):
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
    elif len(matching_line_list) == 1 and not force_list_output:
        return matching_line_list[0]
    elif len(matching_line_list) == 1 and force_list_output:
        return matching_line_list
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
        log.warning("nothing found... check the input args order: 1-regex, 2-bigstring")
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
    return "(" + utils.join_improved("|" , *listin) +  ")"

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
        the file is empty or does not exist
    
    False : 
        the file exists and is not empty
         
    Source
    ------
    http://stackoverflow.com/questions/2507808/python-how-to-check-file-empty-or-not
    """
    return not (os.path.isfile(fpath) and os.path.getsize(fpath) > 0)


def find_recursive(parent_folder , pattern, 
                   sort_results = True, case_sensitive = True,
                   extended_file_stats=False,
                   warn_if_empty=True):
    """
    Find files in a folder and his sub-folders in a recursive way

    Parameters
    ----------
    parent_folder : str
        the parent folder path

    pattern : str
        the researched files pattern name (can manage wildcard or regex)
        * wildcard (only * and ?) for case_sensitive = True
        * regex for case_sensitive = False
        
    sort_results : bool
        Sort results
        
    case_sensitive : bool
        Case sensitive or not. If False, the pattern *must* be a regex
        
    extended_file_stats : bool
        if True, returns the stats of the files
        the outputed matches list will be a list of tuples
        (file_path,stat_object), where stat_object has the following attributes
        
        - st_mode - protection bits,
        - st_ino - inode number,
        - st_dev - device,
        - st_nlink - number of hard links,
        - st_uid - user id of owner,
        - st_gid - group id of owner,
        - st_size - size of file, in bytes,
        - st_atime - time of most recent access,
        - st_mtime - time of most recent content modification,
        - st_ctime - platform dependent; time of most recent metadata 
                     change on Unix, or the time of creation on Windows)
                
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
            #for filename in fnmatch.fnmatch(filenames, pattern):
                matches.append(os.path.join(root, filename))
    else: # not case sensitive, use a regex
        for root, dirnames, filenames in os.walk(parent_folder):
            for filename in filenames:  
                try:
                    bool_match = re.search(pattern, filename, re.IGNORECASE)
                except Exception as e:
                    log.error("if case_sensitive = False, pattern have to be a REGEX (and not only a simple wildcard)")
                    raise e
                    
                if bool_match:
                    matches.append(os.path.join(root, filename))

    if sort_results:
        matches = sorted(matches)
        
    
    if extended_file_stats:
        matches_ext = []
        for f in matches:
            try:
                stat = os.stat(f)
            except FileNotFoundError:
                log.warning("file not found %s",f)
                continue
                
            matches_ext.append((f,stat))
        matches = matches_ext
        
    if warn_if_empty and len(matches) == 0:
        log.warning("no files found! check parent folder and pattern")
        log.info("Parent folder: %s",parent_folder)
        log.info("Pattern      : %s",pattern)
                
    return matches

def glob_smart(dir_path,file_pattern=None,verbose=True):
    if file_pattern:
        dir_path_ok = os.join.path(dir_path,file_pattern)
    else:
        dir_path_ok = dir_path
        
    outlist = glob.glob(dir_path_ok)
    
    if verbose:
        if not outlist:
            log.warning("no file(s) found as %s" , dir_path_ok)
        else:
            log.info("%d file(s) found as %s", len(outlist), dir_path_ok)
            
    return outlist


def insert_lines_in_file(file_path,text_values,lines_ids):
    
    if not utils.is_iterable(text_values):
        text_values = [text_values]
    
    if not utils.is_iterable(lines_ids):
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
            log.info("break at line " , ilin)
            break
                
    f = open(file_path, "w")
    contents = "".join(contents)
    f.write(contents)
    f.close()

    return file_path
    
def uncompress(pathin,dirout = '', opts='-f'):
    log.warn("function discontinued, use files_rw.unzip_gz_Z() instead")
    if not os.path.isfile(pathin):
        log.error('uncompress : %s doesnt exist !!!', pathin)
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


def create_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory

def remove_dir(directory):
    if os.path.exists(directory):
        shutil.rmtree(directory)
    else:
        log.warning('%s dont exists ... skip ...',directory)
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

def fileprint(output,outfile):
    log.debug(output)
    with open(outfile, "a") as f:
        f.write("{}\n".format(output))
    return None




def write_in_file(string_to_write,outdir_or_outpath,
                  outname="",ext='.txt',encoding='utf8',append=False):
    """
    encoding : utf8, latin_1
    https://docs.python.org/3/library/codecs.html#standard-encodings
    
    check the following commented old version if troubles
    """
    
    if outname:
        outpath = os.path.join(outdir_or_outpath,outname + ext)
    else:
        outpath = outdir_or_outpath
        
    if append:
        aw = "a+"
        newline = "\n"
    else:
        aw = "w+"
        newline = ""
        
    F = open(outpath,aw,encoding=encoding)
    F.write(string_to_write + newline)
    F.close()
    return outpath

# def write_in_file(string_to_write,outdir,outname,ext='.txt',encoding='utf8'):
#     """

#     encoding : utf8, latin_1
#     https://docs.python.org/3/library/codecs.html#standard-encodings
#     """
#     outpath = os.path.join(outdir,outname + ext)
#     F = open(outpath,'w+', encoding=encoding)
#     try:    
#         F.write(string_to_write.encode(encode))
#         # astuce de http://stackoverflow.com/questions/6048085/writing-unicode-text-to-a-text-file
#     except UnicodeEncodeError as e:
#         print(string_to_write)
#         raise(e)
#     except TypeError as e:
#         print(e)
#         print("INFO : write_in_file : alternative write following a TypeError")
#         F.write(string_to_write)
        
#     F.close()
#     return outpath

def replace(file_path, pattern, subst):
    """
    Replace a string in a file with a substitute

    Parameters
    ----------
    file_path : str
        path of the file.
    pattern : str
        string to be replaced.
    subst : str
        string which will be substituted.

    Note
    ----
    http://stackoverflow.com/questions/39086/search-and-replace-a-line-in-a-file-in-python

    Returns
    -------
    None.

    """
    #Create temp file
    from tempfile import mkstemp
    from os import close
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
    os.remove(file_path)
    #Move new file
    shutil.move(abs_path, file_path)
    
    
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



def is_exe(fpath):
    """
    Chech if a file is executable

    Parameters
    ----------
    fpath : str
        file path.

    Returns
    -------
    bool
        is executable.

    """
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

