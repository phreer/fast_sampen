#!/bin/bash 

./build/bin/Release/kdtree_sample \
    --input data/chf01.txt \
    --input-format multirecord \
    --input-type double \
    -n 800000 -r 0.1 -m 3 \
    -q -u \
    --quasi-type halton \
    --sample-size 5000 \
    --sample-num 50 \
    --output-level 0
