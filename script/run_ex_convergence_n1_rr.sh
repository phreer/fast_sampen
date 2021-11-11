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
    ./build/bin/experiment_n0n1 \
        --input $filename \
        --input-format multirecord \
        --input-type double \
        --line-offset $line_offset \
        -m $m -r $r -n $n \
        --sample-size-array $sample_size \
        --fast-direct \
        --swr --random --variance \
        --output-level 1 >> $output_file
        # --grid -q -u --random --quasi-type sobol --variance \
}

CONFIG=script/experiment_config.sh
if [ ! -e $CONFIG ]; then
    echo "Configuration file does not exist." >&2
    exit -1
fi
source $CONFIG
input_files=(chbmit/chb07_01.txt)
n=100000
line_offset=0
sample_size=4000 # N0
subdir=convergence_n1_m${m}_r${r}_211104
mkdir -p result/$subdir 2>/dev/null

# AF
input_file_af="RRAF/AF_fa002.rr_multirecord.txt"
input_file=${INPUT_DIR}/${input_file_af}
database=${input_file%/*}
database=${input_file##*/}
output_file=result/${subdir}/n0${sample_size}_r${r}_m${m}_${database}_$(date +%Y-%m-%d).txt 
DoExperimentConvergenceSampleNum $input_file $m $r $line_offset $output_file $n $sample_size &
# Health
input_file_health="RRHealth/Health_Filt-time-1583.01.rr_multirecord.txt"
input_file=${INPUT_DIR}/${input_file_health}
database=${input_file%/*}
database=${input_file##*/}
output_file=result/${subdir}/n0${sample_size}_r${r}_m${m}_${database}_$(date +%Y-%m-%d).txt 
DoExperimentConvergenceSampleNum $input_file $m $r $line_offset $output_file $n $sample_size &
# CHF
n=70000 # Insufficient data length.
input_file_chf="RRCHF/CHF_Filt-time-9643.rr_multirecord.txt"
input_file=${INPUT_DIR}/${input_file_chf}
database=${input_file%/*}
database=${input_file##*/}
output_file=result/${subdir}/n0${sample_size}_r${r}_m${m}_${database}_$(date +%Y-%m-%d).txt 
DoExperimentConvergenceSampleNum $input_file $m $r $line_offset $output_file $n $sample_size &