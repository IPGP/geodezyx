#!/bin/bash

hosnam=`hostname`
dat=`date '+%Y-%m-%d %H:%M:%S'`
rm -fv lib/*pyc
git commit -a -m "$hosnam $dat"
git push origin master
