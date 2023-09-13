#!/bin/bash
# start in the sweep and then goes into each and parses each log file
# bash ../../parseLammpsLogMdextManyfolders.sh


for f in */; do
    cd $f
    bash ../../../../parseLammpsLogMdext.sh
    cd ..
done