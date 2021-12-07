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
    for i in `seq 0 10`; do
        ./build/bin/kdtree_sample \
            --input $filename \
            --input-format multirecord \
            --input-type double \
            --line-offset $line_offset \
            --sample-size $sample_size \
            --sample-num $sample_num \
            -n $n -m $m -r $r \
            --swr --random \
            --variance --n-computation 50 \
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
input_files=(RRAF/AF_fa002.rr_multirecord.txt\
             RRCHF/CHF_Filt-time-9643.rr_multirecord.txt\
             RRHealth/Health_Filt-time-1583.01.rr_multirecord.txt)
subdir=convergence_n_fixedn0n1_final_m${m}_r${r}_n02kn1150_211108
line_offset=0

if [ ! -e result/$subdir ]; then
    mkdir -p result/$subdir
fi
date >> $output_file
echo $comment > result/${subdir}/convergence.txt
for f in ${input_files[@]}; do
    input_file="$INPUT_DIR"/$f
    echo "input_file: $input_file"
    database=${input_file%/*}
    database=${input_file##*/}
    output_file=result/${subdir}/convergence_r${r}_m${m}_${database}_$(date +%Y-%m-%d).txt 
    DoExperimentConvergenceN $input_file $m $r $line_offset $output_file &
done 
