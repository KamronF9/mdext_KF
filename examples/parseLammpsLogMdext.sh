#!/bin/bash
# start in the sweep
# parseLammpsLogMdext.sh


for f in *.log; do
    # let ct=$ct+1
    # f=$(echo "$fgz" | sed 's/.gz//g')  # remove .gz

    # gzip -dc $fgz > ${f}.jdftxout

    # mv ${f}.jdftxout ../jdftxouts/

    python ../../../parseLammpsLogMdext.py -i $f
done