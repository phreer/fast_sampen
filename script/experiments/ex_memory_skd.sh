#!/usr/bin/env bash

input_filename=$1
n=$2
r=$3
m=$4
n0=$5
n1=$6
output_dir=$7
record_name=${input_filename##*/}
record_name=${record_name%.*}

output_name=skd_"$record_name"_n"$n"_n0"$n0"_n1"$n1"
mkdir -p $output_dir 2>/dev/null
build/bin/fast_sampen --input $input_filename -n $n -m $m -r $r \
  --input-format multirecord \
  --sample-num $n1 -skd \
  --output-level 1 > ${output_dir}/${output_name}.txt