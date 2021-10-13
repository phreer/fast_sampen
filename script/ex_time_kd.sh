#!/bin/bash

set -o noclobber

DoExperimentTime() 
{
    local filename=$1 
    local m=$2 
    local r=$3 
    local line_offset=$4 
    local output_file=$5 
    date >> $output_file 
    for i in {0..10}; do 
        n=$(( $(python -c "print(2 ** $i)") * 1000 ))
        ./build/bin/kdtree_sample \
            --input $filename \
            --input-format multirecord \
            --input-type double \
            --line-offset $line_offset \
            -m $m -r $r -n $n \
            --fast-direct \
            --output-level 0 >> $output_file
    done 
}

if [ $# != 4 ]; then 
    echo 'Usage: $0 M R SAMPLE_SIZE SAMPLE_NUM' >&2
    exit -1
fi

CONFIG=script/experiment_config.sh
if [ ! -e $CONFIG ]; then
    echo "Configuration file does not exist." >&2
    exit -1
fi
source $CONFIG
r=$1
m=$2
l=100000
subdir=${subdir}/time_kd_r${r}-m${m}

if [ ! -e result/${subdir} ]; then
    mkdir -p result/$subdir
fi
for f in ${input_files[@]}; do
    input_file='./data.PhysioNet/'$f
    database=${input_file%/*}
    database=${input_file##*/}
    output_file=result/${subdir}/l_${l}-${database}_$(date +%Y-%m-%d).txt 
    DoExperimentTime $input_file ${m} ${r} $l $output_file &
done 