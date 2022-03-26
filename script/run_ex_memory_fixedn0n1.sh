#!/usr/bin/env bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

DoExperimentMemoryFixedn0n1() {
  local input_filename=$1
  local r=$2
  local m=$3
  local output_dir=$4
  for i in `seq 4 20`; do
    local n=$(python -c "print((2 ** $i) * 1024)")
    local sample_size=2048
    local sample_num=150
    if (( $i < 11 )); then
      ${SCRIPT_DIR}/experiments/ex_memory_skd.sh $input_filename $n $r $m $sample_size $sample_num $output_dir
    fi
    ${SCRIPT_DIR}/experiments/ex_memory_mcwor.sh $input_filename $n $r $m $sample_size $sample_num $output_dir
  done
}

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

m=4
r=0.15
INPUT_DIR="$SCRIPT_DIR"/../data.PhysioNet
OUTPUT_DIR=result/memory_220326/fixedn0n1
for input_file in ${input_files[@]}; do
  input_filename=${INPUT_DIR}/$input_file
  DoExperimentMemoryFixedn0n1 $input_filename $r $m $OUTPUT_DIR &
done