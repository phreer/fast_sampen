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
    for i in {1..10}; do 
        n=$(( $i * 200000 ))
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

if [ $# != 4 ]; then 
    echo 'Usage: $0 M R SAMPLE_SIZE SAMPLE_NUM' >&2
    exit -1
fi
r=$1
m=$2
sample_size=$3
sample_num=$4
l=100000
for f in ${input_files[@]}; do
    input_file='./data.PhysioNet/'$f
    database=${input_file%/*}
    database=${input_file##*/}
    output_file=result/time-r${r}-m${m}-n0_${sample_size}-n1_${sample_num}-l_${l}-${database}_$(date +%Y-%m-%d).txt 
    DoExperimentTime $input_file ${m} ${r} $l $output_file $sample_size $sample_num &
done 