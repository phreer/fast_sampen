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
        local n=$(python3 -c "print((2 ** $i) * 4096)")
        local n0n1_str=$(python3 script/get_n0n1_policy2.py $n)
        IFS=" " read -r -a n0n1_array <<< $n0n1_str
        local sample_size=${n0n1_array[0]}
        local sample_num=${n0n1_array[1]}
        # echo "n: $n"
        # echo "sample_size: $sample_size"
        # echo "sample_num: $sample_num"
        ./build/bin/fast_sampen \
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
    echo "Configuration file does not exist. You need to change your working directory to the root of the fast sample entropy project." >&2
    exit -1
fi
source $CONFIG
input_files=(ltafdb/00.txt\
             chbmit/chb07_01.txt\
             mit-bih-long-term-ecg-database-1.0.0/14046.txt\
             pink/pink-1m.txt)
subdir=n_policy2_mem_opt_m${m}_r${r}_220329
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
