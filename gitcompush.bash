#!/bin/bash

hosnam=`hostname`
user="$USER"
dat=`date '+%Y-%m-%d %H:%M:%S'`
rm -fv lib/*pyc
git commit -a -m "$user $hosnam $dat"
git push origin master
