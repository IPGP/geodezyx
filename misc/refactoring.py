#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 09:32:31 2019

@author: psakicki
"""

#### Import the logger
import logging
import os
import re

from geodezyx import utils

log = logging.getLogger(__name__)

################ GENERATION OF THE DICT
p="/home/psakicki/CODES/GeodeZYX-Toolbox_v4/geodezyx"

Lmodul = utils.find_recursive(p,"*py")

Lmodul_exclu = ["legacy_source","quaternions","__init__","legacy","ropeproject"]

Lmodul_clean = [e for e in Lmodul if not utils.patterns_in_string_checker(e,*Lmodul_exclu)]

FctDict = dict()

for mod in Lmodul_clean:
    mod_name = os.path.dirname(mod).split("/")[-1]
    Fcts_in_mod = utils.grep(mod,"def .+\(.*\):",regex=True)
    
    if not Fcts_in_mod:
        continue
    
    if not utils.is_iterable(Fcts_in_mod):
        Fcts_in_mod = [Fcts_in_mod]
        
        
    for fct in Fcts_in_mod:
        if not "self" in fct: #exclude Methods
            fct_name = utils.split_improved(fct,"def ","(")
        
        FctDict[fct_name] = mod_name 

# Manual Cleaning of the Dict
FctDict.pop(".+\\",None)
FctDict.pop("__getitem__",None)
FctDict.pop("__repr__",None)
FctDict.pop("__call__",None)
FctDict.pop("__init__",None)
FctDict.pop("write",None)
FctDict.pop("trunc",None)
FctDict.pop("flush",None)
FctDict.pop("replace",None)
    
############### FIND CALLs
    

p="/home/psakicki/CODES/GeodeZYX-Toolbox_v4/geodezyx"

Lmodul_replac = utils.find_recursive(p,"*py")

Lmodul_exclu_replac = ["legacy_source","quaternions","__init__","legacy","ropeproject"]

Lmodul_clean_replac = [e for e in Lmodul_replac if not utils.patterns_in_string_checker(e,*Lmodul_exclu_replac)]    

#class ParseCall(ast.NodeVisitor):
#    def __init__(self):
#        self.ls = []
#    def visit_Attribute(self, node):
#        ast.NodeVisitor.generic_visit(self, node)
#        self.ls.append(node.attr)
#    def visit_Name(self, node):
#        self.ls.append(node.id)
#
#class FindFuncs(ast.NodeVisitor):
#    def visit_Call(self, node):
#        p = ParseCall()
#        p.visit(node.func)
#        print(".".join(p.ls))
#        ast.NodeVisitor.generic_visit(self, node)

ExcludedLines = ["pt.helmert_trans(params,invert) for pt in tsout.pts",
                 "geok.helmert_trans(np.array([self.X,self.Y,self.Z]),params,invert)",
                 "return Counter(L).most_common(1)[0][0]"]

for fct,mod_good in FctDict.items():
    log.info("FUNCTION",fct)
    for mod_replac in Lmodul_clean_replac:
        for l in open(mod_replac):
            
            match = False
            
#            pat1="(\w|_|-)*\."+fct+"\("
#            repl1=""+mod_good+"."+fct + "("
#
#            pat2="=(\w|_|-)*\."+fct+"\("
#            repl2="="+mod_good+"."+fct + "("
#
#            pat3="\((\w|_|-)*\."+fct+"\("
#            repl3="("+mod_good+"."+fct + "("


            pat1="(\S)*\."+fct+"\("
            repl1=""+mod_good+"."+fct + "("

            pat2="=(\S)*\."+fct+"\("
            repl2="="+mod_good+"."+fct + "("

            pat3="\((\S)*\."+fct+"\("
            repl3="("+mod_good+"."+fct + "("

            
            if utils.patterns_in_string_checker(l,*ExcludedLines):
                continue
            
            elif re.search(pat2,l):
                lnew = re.sub(pat2,repl2,l)
                match = True
            elif re.search(pat3,l):
                lnew = re.sub(pat3,repl3,l)
                match = True
            elif re.search(pat1,l):
                lnew = re.sub(pat1,repl1,l)
                match = True
                


            if match:
                if l == lnew:
                    continue
                log.info("BEF:",l)
                log.info("AFT:",lnew)
                #utils.replace(mod_replac,l,lnew)


            
        
    
