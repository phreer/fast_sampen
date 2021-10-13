#!/bin/bash

set -o noclobber

DoExperimentConvergenceSampleNum()
{
    local filename=$1
    local m=$2
    local r=$3
    local line_offset=$4
    local output_file=$5
    local n=$6
    local sample_size=$7
    date >> $output_file
    for i in {1..20}; do 
        local sample_num=$(( $i * 20 ))
        ./build/bin/kdtree_sample \
            --input $filename \
            --input-format multirecord \
            --input-type double \
            --line-offset $line_offset \
            -m $m -r $r -n $n \
            --sample-size $sample_size \
            --sample-num $sample_num \
            --fast-direct \
            --swr --grid -q -u --random --quasi-type sobol --variance \
            --output-level 0 >> $output_file
    done 
}

CONFIG=script/experiment_config.sh
if [ ! -e $CONFIG ]; then
    echo "Configuration file does not exist." >&2
    exit -1
fi
source $CONFIG

n=1000000
sample_size=2000
subdir=swr
for f in ${input_files[@]}; do
    input_file='./data.PhysioNet/'$f
    database=${input_file%/*}
    database=${input_file##*/}
    output_file=result/${subdir}/convergence_n1/n0${sample_size}_r${r}_m${m}_${database}_$(date +%Y-%m-%d).txt 
    DoExperimentConvergenceSampleNum $input_file $m $r 100000 $output_file $n $sample_size &
done 