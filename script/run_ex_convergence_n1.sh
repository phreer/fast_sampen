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
    for i in {1..25}; do 
        local sample_num=$(( $i * 10 ))
        ./build/bin/kdtree_sample \
            --input $filename \
            --input-format multirecord \
            --input-type double \
            --line-offset $line_offset \
            -m $m -r $r -n $n \
            --sample-size $sample_size \
            --sample-num $sample_num \
            --fast-direct \
            --swr \
            --output-level 0 >> $output_file
            # --grid -q -u --random --quasi-type sobol --variance \
    done 
}

CONFIG=script/experiment_config.sh
if [ ! -e $CONFIG ]; then
    echo "Configuration file does not exist." >&2
    exit -1
fi
source $CONFIG

n=1000010
sample_size=2000 # N0
subdir=convergence_n1_m${m}_r${r}_211031
mkdir -p result/$subdir 2>/dev/null
for f in ${input_files[@]}; do
    input_file=${INPUT_DIR}/$f
    database=${input_file%/*}
    database=${input_file##*/}
    output_file=result/${subdir}/n0${sample_size}_r${r}_m${m}_${database}_$(date +%Y-%m-%d).txt 
    DoExperimentConvergenceSampleNum $input_file $m $r $line_offset $output_file $n $sample_size &
done 