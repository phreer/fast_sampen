#!/bin/bash

set -o noclobber

DoExperimentTime() 
{
    local filename=$1 
    local m=$2 
    local r=$3 
    local line_offset=$4 
    local output_file=$5 
    local sample_size=$6 
    local sample_num=$7
    date >> $output_file 
    for i in {0..5}; do 
        n=$(( $(python -c "print(2 ** $i)") * 50000 ))
        ./build/bin/kdtree_sample \
            --input $filename \
            --input-format multirecord \
            --input-type double \
            --line-offset $line_offset \
            -m $m -r $r -n $n \
            --sample-size $sample_size \
            --sample-num $sample_num \
            -q -u --random --quasi-type sobol \
            --output-level 0 >> $output_file
    done 
}

input_files=(chfdb/chf01.txt\
             ltafdb/00.txt\
             ltstdb/s20011.txt\
             mghdb/mgh001.txt\
             pink/pink_noise-2000000.txt\
             gaussian/gaussian_noise-2000000.txt\
             surrogate-data-with-correlations-trends-and-nonstationarities-1.0.0/tns/d2h4pd050918s_2.txt
             mit-bih-long-term-ecg-database-1.0.0/14046.txt)
m=5
r=0.3
sample_size=4000
sample_num=50
for f in ${input_files[@]}; do
    input_file='./data.PhysioNet/'$f
    database=${input_file%/*}
    database=${input_file##*/}
    output_file=result/time_r${r}_m${m}_${database}_$(date +%Y-%m-%d).txt 
    DoExperimentTime $input_file 3 0.1 100000 $output_file $sample_size $sample_num &
done 
