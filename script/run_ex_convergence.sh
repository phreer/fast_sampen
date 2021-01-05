#!/bin/bash

set -o noclobber

DoExperimentConvergenceSampleSize()
{
    local filename=$1
    local m=$2
    local r=$3
    local line_offset=$4
    local output_file=$5
    local n=$6
    local sample_num=$7
    date >> $output_file
    for i in {1..20}; do 
        sample_size=$(( $i * 500 ))
        ./build/bin/kdtree_sample \
            --input $filename \
            --input-format multirecord \
            --input-type double \
            --line-offset $line_offset \
            -m $m -r $r -n $n \
            --sample-size $sample_size \
            --sample-num $sample_num \
            --swr --grid -q -u --random --quasi-type sobol --presort \
            --variance \
            --output-level 0 >> $output_file
    done 
}

input_files=(chfdb/chf01.txt\
             ltafdb/00.txt\
             ltstdb/s20011.txt\
             mghdb/mgh001.txt\
             pink/pink_noise-2000000.txt\
             gaussian/gaussian_noise-2000000.txt\
             surrogate-data-with-correlations-trends-and-nonstationarities-1.0.0/tns/d2h4pd050918s_2.txt
             mit-bih-long-term-ecg-database-1.0.0/14046.txt)

CONFIG=script/experiment_config.sh
if [ ! -e $CONFIG ]; then
    echo "Configuration file does not exist." >&2
    exit -1
fi
source $CONFIG
n=300000
sample_num=50
line_offset=100000

if [ ! -e result/$subdir ]; then
    mkdir -p result/$subdir
fi
comment="QMC using lattice sequence."
echo comment > result/${subdir}/convergence.txt
for f in ${input_files[@]}; do
    input_file='./data.PhysioNet/'$f
    database=${input_file%/*}
    database=${input_file##*/}
    output_file=result/${subdir}/convergence_r${r}_m${m}_${database}_$(date +%Y-%m-%d).txt 
    DoExperimentConvergenceSampleSize $input_file $m $r $line_offset $output_file $n $sample_num &
done 