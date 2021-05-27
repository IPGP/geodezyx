#!/bin/bash

hosnam=`hostname`
user="$USER"
if [ $USER eq 'psakicki' ]
    then user="psakic"
fi
dat=`date '+%Y-%m-%d %H:%M:%S'`
#rm -fv lib/*pyc
#find . -name "__pycache__" -exec rm -fvr "{}" \;
#find . -name "*.pyc" -exec rm -fv "{}" \;
find . -name "__pycache__" | xargs rm -fvr
find . -name "*.pyc"       | xargs rm -fv
git commit -a -m "$user $hosnam $dat"
git push origin master
