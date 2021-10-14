#!/bin/bash

set -o noclobber

DoExperimentConvergenceN()
{
    local filename=$1
    local m=$2
    local r=$3
    local line_offset=$4
    local output_file=$5
    for i in `seq 0 10`; do
        local n=$(python -c "print((2 ** $i) * 1024)")
        local sample_size=$(python -c "import math; print(int(math.ceil(math.sqrt($n))))")
        local sample_num=$(python -c "import math; print(int(math.ceil(math.log($n, 2))))")
        echo "sample_size: $sample_size"
        echo "sample_num: $sample_num"
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
    done
}


CONFIG=script/experiment_config.sh
if [ ! -e $CONFIG ]; then
    echo "Configuration file does not exist." >&2
    exit -1
fi
source $CONFIG

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