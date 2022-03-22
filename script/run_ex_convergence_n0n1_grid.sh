#!/bin/bash

set -o noclobber

DoExperimentConvergenceSampleSize()
{
    local filename=$1
    local n=$2
    local m=$3
    local r=$4
    local line_offset=$5
    local output_file=$6
    date >> $output_file
    ./build/bin/experiment_n0n1 \
        --input $filename \
        --input-format multirecord \
        --input-type double \
        --line-offset $line_offset \
        -n $n -m $m -r $r \
        --swr -u -q \
        --variance --random \
        --output-level 1 >> $output_file
}


CONFIG=script/experiment_config.sh
if [ ! -e $CONFIG ]; then
    echo "Configuration file does not exist." >&2
    exit -1
fi
source $CONFIG
set -x
n=500000
r=0.15
m=5
subdir=convergence_n0n1_grid/m${m}_r${r}_220308

if [ ! -e result/$subdir ]; then
    mkdir -p result/$subdir
fi
comment="Experiment on the convergence of sample size (N_0) and sample number (N_1)"
echo $comment > result/${subdir}/convergence.txt
for f in ${input_files[@]}; do
    input_file="$INPUT_DIR"/$f
    echo "input_file: $input_file"
    database=${input_file%/*}
    database=${input_file##*/}
    output_file=result/${subdir}/convergence_r${r}_m${m}_${database}.txt 
    DoExperimentConvergenceSampleSize $input_file $n $m $r $line_offset $output_file &
done 