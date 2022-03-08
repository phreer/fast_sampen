#!/bin/bash
# 方法: KDTreeMCSampEn
# 目标: 误差与 N_0 的变化关系, N_0 = [200, 400, 600, 800,..., 4000]
# 参数: N_1 = 1, r = 0.15, m = 4, n = 1000010
# 数据集: chfdb/chf01.txt, ltafdb/00.txt, ltstdb/s20011.txt,
#     mghdb/mgh001.txt, chbmit/chb07_01.txt,
#     mit-bih-long-term-ecg-database-1.0.0/14046.txt,
#     pink/pink-1m.txt, uniform/uniform-1m.txt, gaussian/gaussian-1m.txt
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
    local sample_size_array=$8
    date >> $output_file
    ./build/bin/experiment_n0n1 \
            --input $filename \
            --input-format multirecord \
            --input-type double \
            --line-offset $line_offset \
            -m $m -r $r -n $n \
            --sample-num-array $sample_num \
            --sample-size-array $sample_size_array \
            --kdtree-sample --random --variance \
            --output-level 1 >> $output_file
}

INPUT_DIR=./data.PhysioNet
input_files=(chfdb/chf01.txt\
             ltafdb/00.txt\
             ltstdb/s20011.txt\
             mghdb/mgh001.txt\
             chbmit/chb07_01.txt\
             mit-bih-long-term-ecg-database-1.0.0/14046.txt\
             uniform/uniform-1m.txt\
             pink/pink-1m.txt\
             gaussian/gaussian-1m.txt)

m=4
r=0.15
line_offset=100000
n=1000010
sample_num=1 # N1
sample_size_array=$(seq -s',' 2000 2000 100000)
output_level=1
subdir=convergence_n0_m${m}_r${r}_220308
mkdir -p result/$subdir 2>/dev/null
for f in ${input_files[@]}; do
    input_file=${INPUT_DIR}/$f
    database=${input_file%/*}
    database=${input_file##*/}
    output_file=result/${subdir}/n1${sample_num}_r${r}_m${m}_${database}_$(date +%Y-%m-%d).txt 
    DoExperimentConvergenceSampleSize $input_file $m $r $line_offset $output_file $n $sample_num $sample_size_array &
done 