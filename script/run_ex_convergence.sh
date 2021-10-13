#!/bin/bash

set -o noclobber

DoExperimentConvergenceSampleSize()
{
    local filename=$1
    local m=$2
    local r=$3
    local line_offset=$4
    local output_file=$5
    date >> $output_file
    ./build/bin/experiment_n0n1 \
        --input $filename \
        --input-format multirecord \
        --input-type double \
        --line-offset $line_offset \
        -n $n -m $m -r $r \
        --swr --random \
        --variance \
        --output-level 0 >> $output_file
}


CONFIG=script/experiment_config.sh
if [ ! -e $CONFIG ]; then
    echo "Configuration file does not exist." >&2
    exit -1
fi
source $CONFIG
line_offset=100000

if [ ! -e result/$subdir ]; then
    mkdir -p result/$subdir
fi
comment="QMC using lattice sequence."
echo $comment > result/${subdir}/convergence.txt
for f in ${input_files[@]}; do
    input_file='./data.PhysioNet/'$f
    database=${input_file%/*}
    database=${input_file##*/}
    output_file=result/${subdir}/convergence_r${r}_m${m}_${database}_$(date +%Y-%m-%d).txt 
    DoExperimentConvergenceSampleSize $input_file $m $r $line_offset $output_file &
done 