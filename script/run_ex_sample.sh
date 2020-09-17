#!/bin/bash 

./build/bin/kdtree_sample \
    --input data.PhysioNet/chfdb/chf01.txt \
    --input-format multirecord \
    --input-type double \
    -n 300000 -r 0.3 -m 5 \
    -q -u \
    --random \
    --quasi-type halton \
    --sample-size 2000 \
    --sample-num 50 \
    --output-level 0
