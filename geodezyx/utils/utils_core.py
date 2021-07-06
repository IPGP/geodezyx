#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: psakic

This sub-module of geodezyx.utils contains functions for misc. 
low level function 

it can be imported directly with:
from geodezyx import utils

The GeodeZYX Toolbox is a software for simple but useful
functions for Geodesy and Geophysics under the GNU GPL v3 License

Copyright (C) 2019 Pierre Sakic et al. (GFZ, pierre.sakic@gfz-postdam.de)
GitHub repository :
https://github.com/GeodeZYX/GeodeZYX-Toolbox_v4
"""


########## BEGIN IMPORT ##########
#### External modules
import datetime as dt
import gzip
import inspect
import numpy as np
import os
import pandas as pd
import pickle
import re
import scipy
import sys
import tempfile
import time
import uuid
#### geodeZYX modules


##########  END IMPORT  ##########



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
    """
    Check if the code is run inside Spyder IDE
    """
    if any('SPYDER' in name for name in os.environ):
        return True
    else:
        return False

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


def transpose_vector_array(X):
    """
    transpose a Nx3 array to a 3xN array if necessary
    (necessary for some usages)
    
    Parameters
    ----------
    X : iterable
    
    Returns
    -------
    X : iterable
        X transposed if necessary.

    """
    
    if isinstance(X,pd.core.frame.DataFrame):
        X = X.values
        
    if isinstance(X,np.ndarray):
        if not X.shape[0] == 3:
            X = X.T
    return X


def is_lambda(v):
    """
    Check if v is lambda

    Source 
    -------
    https://stackoverflow.com/questions/3655842/how-can-i-test-whether-a-variable-holds-a-lambda
    """
    
    LAMBDA = lambda:0
    return isinstance(v, type(LAMBDA)) and v.__name__ == LAMBDA.__name__



def contains_word(s, w):
    return f' {w} ' in f' {s} '


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


def indice_printer(i,print_every = 10,text_before=''):
    """
    print an index every N iteration
    """
    if np.mod(i,print_every) == 0:
        print(text_before ,  i)
    return None



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



def globals_filtered():
    """
    Filter globals() varirables with only compatible variables for pickle.
    
    https://stackoverflow.com/questions/2960864/how-to-save-all-the-variables-in-the-current-python-session


    Returns
    -------
    data_out : dict
        filtered globals() variables.

    """
    from spyder_kernels.utils.nsview import globalsfilter,get_supported_types

    data = globals()
    
    #settings = VariableExplorer.get_settings()


    data_out = globalsfilter(globals(),                   
                         check_all=True,
                         filters=tuple(get_supported_types()['picklable']),
                         exclude_private=True,
                         exclude_uppercase=False,
                         exclude_capitalized=True,
                         exclude_unsupported=True,
                         excluded_names=[],
                         exclude_callables_and_modules=True)
    
    return data_out



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

def str2int_float_autodetect(str_list_in):
    intflt_list_out = []
    for str_in in str_list_in:
        flt_tmp = float(str_in)
        if np.isclose(flt_tmp - np.floor(flt_tmp),0):
            val = int(flt_tmp)
        else:
            val = flt_tmp
        intflt_list_out.append(val)
    return(intflt_list_out)
            

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
        

def line_count(filein):
    num_lines = sum(1 for line in open(filein))
    return num_lines


def patterns_in_string_checker(string,*patterns):
    """
    recipe to the famous problem of pattern in string
    from
    http://stackoverflow.com/questions/3389574/check-if-multiple-strings-exist-in-another-string
    """
    L = [x in string for x in patterns]
    return bool(any( L ))


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
                                           nth_occur_elt_end=0,
                                           invert=False,
                                           verbose=False):
    """
    This function is based on REGEX (elt_start , elt_end are REGEX)
    and can manage several blocks in the same file

    return_string = True  : returns a string of the matched lines
    return_string = False : returns a list of the matched lines
    invert : exclude text between the pattern
    
    NB : in SINEX context, with "+MARKER", use backslash i.e.
         "\\+MARKER"    
         
    NB2 : think about StingIO for a Pandas DataFrame Handeling
    https://docs.python.org/2/library/stringio.html
    """
    
    
    ## 3 possibilities
    # file_path is a path, uncompressed
    # file_path in a path, gzip compressed
    # file_path is already a list of lines
    
    if type(file_path) is str: ### case 1 : path compressed 
        if file_path[-2:] in ("gz","GZ"):
            F = gzip.open(file_path, "r+")
            F = [e.decode('utf-8') for e in F]
        else:                  ### case 2 : path uncompressed
            try:
                F = open(file_path,"r",encoding = "ISO-8859-1")
            except:
                F = open(file_path,"r")
    else:                      ### case 3 : already a list of lines
        F = file_path
    
    out_lines_list = []
    
    trigger = False
    if invert:
        trigger = not trigger
        
    last_triggered_line = False
    
    for line in F:

        if re.search(elt_start , line) and nth_occur_elt_start == 0:
            trigger = True
            if invert:
                trigger = not trigger
        elif re.search(elt_start , line) and nth_occur_elt_start > 0:
            nth_occur_elt_start = nth_occur_elt_start - 1
        
        if re.search(elt_end , line) and nth_occur_elt_end == 0:
            last_triggered_line = True
        elif trigger and re.search(elt_end , line) and nth_occur_elt_end > 0:
            nth_occur_elt_end = nth_occur_elt_end - 1
                        
        if trigger:
            if verbose:
                print(line)
            out_lines_list.append(line)
            
        if last_triggered_line:
            last_triggered_line = False
            trigger = False
            if invert:
                trigger = not trigger
            

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


def replace_in_file(file_in,str_before,str_after):
    """
    https://stackoverflow.com/questions/17140886/how-to-search-and-replace-text-in-a-file
    """
    with open(file_in, 'r') as file :
        filedata = file.read()
        
    # Replace the target string
    filedata = filedata.replace(str_before,str_after)
        
    # Write the file out again
    with open(file_in + "", 'w') as file:
        file.write(filedata)





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


def get_type_smart(obj_in):
    """
    get type of an object, to convert easily another one to this type
    for instance type(np.array(A)) doesn't return a constructor 
    """
    typ=type(obj_in)
    if typ is np.ndarray:
        return np.array
    else:
        return typ
    

class Tee(object):
    """
    Internal class for Tee_frontend
    
    Source
    ------
        based on
        http://stackoverflow.com/questions/11325019/output-on-the-console-and-file-using-python
        Secondary links
        http://stackoverflow.com/questions/616645/how-do-i-duplicate-sys-stdout-to-a-log-file-in-python
        http://stackoverflow.com/questions/2996887/how-to-replicate-tee-behavior-in-python-when-using-subprocess
    """
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

def Tee_frontend(dir_in,logname_in,suffix='',ext='log',print_timestamp=True):
    """
    Write in a file the console output

    Parameters
    ----------
    dir_in : str
        directory path.
    logname_in : str
        logfile name.
    suffix : str, optional
        An optional suffix. The default is ''.
    ext : str, optional
        file extension. The default is 'log'.
    print_timestamp : bool, optional
        print a timestamp in the filename. The default is True.

    Returns
    -------
    F_tee : F_tee object
        Object controling the output
        
    Note
    ----
        It is recommended to stop the writing at the end of the script
        with F_tee.stop()

    """
    #if suffix != '':
    #    suffix = suffix + '_'
    if print_timestamp:
        ts = '_' + get_timestamp()
    else:
        ts = ''
    f = open(os.path.join(dir_in,logname_in +'_'+suffix+ts+'.'+ext), 'w')
    F_tee = Tee(f)
    #sys.stdout = f_tee
    return F_tee

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

def greek_alphabet(num=None,maj=False):
    if not num:
        a = ['\u03B1',
        '\u03B2',
        '\u03B3',
        '\u03B4',
        '\u03B5',
        '\u03B6',
        '\u03B7',
        '\u03B8',
        '\u03B9',
        '\u03BA',
        '\u03BB',
        '\u03BC',
        '\u03BD',
        '\u03BE',
        '\u03BF',
        '\u03C0',
        '\u03C1',
        '\u03C3',
        '\u03C4',
        '\u03C5',
        '\u03C6',
        '\u03C7',
        '\u03C8',
        '\u03C9']

        A=['\u0391',
        '\u0392',
        '\u0393',
        '\u0394',
        '\u0395',
        '\u0396',
        '\u0397',
        '\u0398',
        '\u0399',
        '\u039A',
        '\u039B',
        '\u039C',
        '\u039D',
        '\u039E',
        '\u039F',
        '\u03A0',
        '\u03A1',
        '\u03A3',
        '\u03A4',
        '\u03A5',
        '\u03A6',
        '\u03A7',
        '\u03A8',
        '\u03A9']

        if maj:
            return A
        else:
            return a
    else:
        return greek_alphabet()[num]

    
def trunc(f, n):
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
