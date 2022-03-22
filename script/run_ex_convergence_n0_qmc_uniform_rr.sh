#!/bin/bash
# 方法: QMCSampEn, MCSampEn (with replacement)
# 目标: 误差与 N_0 的变化关系, N_0 = [200, 400, 600, 800,..., 4000]
# 参数: N_1 = 250, r = 0.15, m = 4, n = 1000010
# 数据集: RR/Health, RR/CHF, RR/AF

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
    ./build/bin/experiment_n0n1 \
            --input $filename \
            --input-format multirecord \
            --input-type double \
            --line-offset $line_offset \
            -m $m -r $r -n $n \
            --sample-num-array $sample_num \
            --fast-direct \
            --swr -q -u --quasi-type sobol --random --variance \
            --output-level 1 >> $output_file
}

m=4
r=0.15
line_offset=0
INPUT_DIR=./data.PhysioNet
sample_num=250 # N1

subdir=convergence_n0_m${m}_r${r}_220308
mkdir -p result/$subdir 2>/dev/null

n=100000
# AF
input_file_af="RRAF/AF_fa002.rr_multirecord.txt"
input_file=${INPUT_DIR}/${input_file_af}
database=${input_file%/*}
database=${input_file##*/}
output_file=result/${subdir}/n1${sample_num}_r${r}_m${m}_${database}_$(date +%Y-%m-%d).txt
DoExperimentConvergenceSampleSize $input_file $m $r $line_offset $output_file $n $sample_num &
# Health
input_file_health="RRHealth/Health_Filt-time-1583.01.rr_multirecord.txt"
input_file=${INPUT_DIR}/${input_file_health}
database=${input_file%/*}
database=${input_file##*/}
output_file=result/${subdir}/n1${sample_num}_r${r}_m${m}_${database}_$(date +%Y-%m-%d).txt
DoExperimentConvergenceSampleSize $input_file $m $r $line_offset $output_file $n $sample_num &
# CHF
n=70000 # Insufficient data length.
input_file_chf="RRCHF/CHF_Filt-time-9643.rr_multirecord.txt"
input_file=${INPUT_DIR}/${input_file_chf}
database=${input_file%/*}
database=${input_file##*/}
output_file=result/${subdir}/n1${sample_num}_r${r}_m${m}_${database}_$(date +%Y-%m-%d).txt
DoExperimentConvergenceSampleSize $input_file $m $r $line_offset $output_file $n $sample_num &
