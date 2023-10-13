#!/bin/bash

for dir in */;do
    cd $dir
    python /home/kamron/mdext_KF/examples/plot1DLJmanyPressures.py
    cd ..
done