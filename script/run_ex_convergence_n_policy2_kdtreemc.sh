#!/bin/bash

set -o noclobber

DoExperimentConvergenceN()
{
    local filename=$1
    local m=$2
    local r=$3
    local line_offset=$4
    local output_file=$5
    for i in `seq 0 8`; do
        local n=$(python3 -c "print((2 ** $i) * 4096)")
        local sample_size=$(python3 script/get_n0_policy2_kdtreemc.py $n)
        # echo "n: $n"
        # echo "sample_size: $sample_size"
        # echo "sample_num: $sample_num"
        ./build/bin/fast_sampen \
            --input $filename \
            --input-format multirecord \
            --input-type double \
            --line-offset $line_offset \
            --sample-size $sample_size \
            --sample-num 1 \
            -n $n -m $m -r $r \
            --kdtree-sample --random \
            -rkd \
            --variance --n-computation 50 \
            --output-level 0 >> $output_file
    done
}


m=6
r=0.15
INPUT_DIR=./data.PhysioNet
input_files=(chfdb/chf01.txt\
             ltafdb/00.txt\
             ltstdb/s20011.txt\
             mghdb/mgh001.txt\
             chbmit/chb07_01.txt\
             mit-bih-long-term-ecg-database-1.0.0/14046.txt\
             RRAF/AF_fa002.rr_multirecord.txt\
             RRHealth/Health_Filt-time-1583.01.rr_multirecord.txt\
             RRCHF/CHF_Filt-time-9643.rr_multirecord.txt\
             gaussian/gaussian-1m.txt\
             pink/pink-1m.txt\
             uniform/uniform-1m.txt)

subdir=convergence_n_100log2n_kdtreemc_220316/r"$r"_m"$m"

if [ ! -e result/$subdir ]; then
    mkdir -p result/$subdir
fi

date >> $output_file

for f in ${input_files[@]}; do
    input_file="$INPUT_DIR"/$f
    if [[ "$f" == "RR"* ]]; then
        line_offset=0
    else
        line_offset=100000
    fi
    echo "input_file: $input_file"
    echo "line_offset: $line_offset"
    database=${input_file%/*}
    database=${input_file##*/}
    output_file=result/${subdir}/convergence_r${r}_m${m}_${database}_$(date +%Y-%m-%d).txt 
    DoExperimentConvergenceN $input_file $m $r $line_offset $output_file
done 
