#!/bin/bash

set -o noclobber

DoExperimentTimeKDTree() 
{
    local filename=$1 
    local m=$2 
    local r=$3 
    local line_offset=$4 
    local output_file=$5 
    date >> $output_file 
    for i in {0..11}; do
        n=$(( $(python -c "print(2 ** $i)") * 2000 ))
        ./build/bin/fast_sampen \
            --input $filename \
            --input-format multirecord \
            --input-type double \
            --line-offset $line_offset \
            -m $m -r $r -n $n \
            --fast-direct -rkd --simple-kdtree \
            --output-level 0 >> $output_file
    done 
}

if [ $# != 3 ]; then 
    echo 'Usage: $0 <r> <m> <tag>' >&2
    exit -1
fi


r=$1
m=$2
tag=$3
l=100000
subdir=time_kdtree/${tag}_r${r}-m${m}

if [ ! -e result/${subdir} ]; then
    mkdir -p result/$subdir
fi

input_files=(chfdb/chf01.txt\
             ltafdb/00.txt\
             ltstdb/s20011.txt\
             mghdb/mgh001.txt\
             chbmit/chb07_01.txt\
             mit-bih-long-term-ecg-database-1.0.0/14046.txt\
             pink/pink-1m.txt\
             gaussian/gaussian-1m.txt\
             uniform/uniform-1m.txt\
             chbmit/chb07_01.txt\
             RRCHF/CHF_Filt-time-9643.rr_multirecord.txt\
             RRHealth/Health_Filt-time-1583.01.rr_multirecord.txt\
             RRAF/AF_fa002.rr_multirecord.txt)

for f in ${input_files[@]}; do
    input_file='./data.PhysioNet/'$f
    database=${input_file%/*}
    database=${input_file##*/}
    output_file=result/${subdir}/${database}_$(date +%Y-%m-%d).txt 
    DoExperimentTimeKDTree $input_file ${m} ${r} $l $output_file
done 