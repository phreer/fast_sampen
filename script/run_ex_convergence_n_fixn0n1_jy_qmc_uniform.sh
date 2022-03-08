#!/bin/bash
# 方法: MCSampEn (with and without replacement),
# 目标: 固定 N_0 和 N_1, N 变化的误差
# 参数: N_1 = 150, N_0 = 2000, r = 0.15, m = 4
# 数据集: chfdb/chf01.txt, ltafdb/00.txt, ltstdb/s20011.txt,
#     mghdb/mgh001.txt, chbmit/chb07_01.txt,
#     mit-bih-long-term-ecg-database-1.0.0/14046.txt,
#     pink/pink-1m.txt, uniform/uniform-1m.txt, gaussian/gaussian-1m.txt
set -o noclobber

DoExperimentConvergenceN()
{
    local filename=$1
    local m=$2
    local r=$3
    local line_offset=$4
    local output_file=$5
    local sample_size=2000
    local sample_num=150
    n=2048
    date >> $output_file
    for i in `seq 0 9`; do
        ./build/bin/fast_sampen \
            --input $filename \
            --input-format multirecord \
            --input-type double \
            --line-offset $line_offset \
            --sample-size $sample_size \
            --sample-num $sample_num \
            -n $n -m $m -r $r \
            --swr -q --uniform --random \
            --variance --n-computation 50 \
            --output-level 1 >> $output_file
        n=$(python -c "print(int($n * 2))")
    done
}


output_level=1
m=3
r=0.15
line_offset=100000
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
subdir=convergence_n_fixedn0n1_m${m}_r${r}_n02kn1150_220308

if [ ! -e result/$subdir ]; then
    mkdir -p result/$subdir
fi
echo $comment > result/${subdir}/convergence.txt
for f in ${input_files[@]}; do
    input_file="$INPUT_DIR"/$f
    echo "input_file: $input_file"
    database=${input_file%/*}
    database=${input_file##*/}
    output_file=result/${subdir}/convergence_r${r}_m${m}_${database}_$(date +%Y-%m-%d).txt 
    DoExperimentConvergenceN $input_file $m $r $line_offset $output_file
done 
