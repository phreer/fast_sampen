#!/usr/bin/env bash

PROG=$HOME/workspace/sampen_kdtree/build/bin/sampen2d
DoExperimentSampEn2DTime() {
  local filename=$1
  local m=$2
  local r=$3
  local sample_size=$4
  local sample_num=$5
  local output_filename=$6
  "$PROG" -m $m -r $r --sample-size $sample_size --sample-num $sample_num \
      --sampling1 \
      --input $filename > $output_filename
}

if [ $# -ne 6 ]; then
  echo "$#" parameters got.
  echo "Usage: ex_time_2d.sh <filename> <m> <r> <sample_size> <sample_num> <output_dir>" >&2
  exit 1
fi

filename=$1
m=$2
r=$3
sample_size=$4
sample_num=$5
outputdir=$6
outputname=$(basename "$filename")
outputname="$outputdir"/"${outputname%.*}".txt

DoExperimentSampEn2DTime $filename $m $r $sample_size $sample_num $outputname