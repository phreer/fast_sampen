#!/bin/bash 

set -o noclobber
input_files=(pink/pink_noise-2000000.txt\
             gaussian/gaussian_noise-2000000.txt\
             surrogate-data-with-correlations-trends-and-nonstationarities-1.0.0/tns/d2h4pd050918s_2.txt)

r=0.1
m=3
output_file=result/offset_time_r01_m${m}_surrogate.txt.debug
for file in ${input_files[@]}; do 
    ./build/bin/kdtree_sample \
        --input data.PhysioNet/chfdb/chf01.txt \
        --input-format multirecord \
        --input-type double \
        -n 100000 -r 0.1 -m 3 \
        -q -u \
        --quasi-type halton \
        --sample-size 4000 \
        --sample-num 50 \
        --random \
        --output-level 0
done