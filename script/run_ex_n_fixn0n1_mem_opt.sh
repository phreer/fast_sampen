#!/bin/bash

set -o noclobber

DoExperimentConvergenceN()
{
    local filename=$1
    local m=$2
    local r=$3
    local line_offset=$4
    local output_file=$5
    local sample_size=2000
    local sample_num=150
    n=2048
    for i in `seq 0 9`; do
        ./build/bin/fast_sampen \
            --input $filename \
            --input-format multirecord \
            --input-type double \
            --line-offset $line_offset \
            --sample-size $sample_size \
            --sample-num $sample_num \
            -n $n -m $m -r $r \
            --swr --random \
            --variance --n-computation 1 \
            --output-level 1 >> $output_file
        n=$(python -c "print(int($n * 2))")
    done
}


CONFIG=script/experiment_config.sh
if [ ! -e $CONFIG ]; then
    echo "Configuration file does not exist." >&2
    exit -1
fi
source $CONFIG
input_files=(ltafdb/00.txt\
             chbmit/chb07_01.txt\
             mit-bih-long-term-ecg-database-1.0.0/14046.txt\
             pink/pink-1m.txt)

subdir=n_fixedn0n1_mem_opt_m${m}_r${r}_n02kn1150_220329

if [ ! -e result/$subdir ]; then
    mkdir -p result/$subdir
fi
date >> $output_file

for f in ${input_files[@]}; do
    input_file="$INPUT_DIR"/$f
    echo "input_file: $input_file"
    database=${input_file%/*}
    database=${input_file##*/}
    output_file=result/${subdir}/convergence_r${r}_m${m}_${database}_$(date +%Y-%m-%d).txt 
    DoExperimentConvergenceN $input_file $m $r $line_offset $output_file &
done 
