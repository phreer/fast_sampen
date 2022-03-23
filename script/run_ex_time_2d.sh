#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
INPUT_DIR=$HOME/workspace/entropy/data/original_image

m=2
r=0.25
n0=2048
n1=20
outputdir=result/time_2d_original/m"$m"_r"$r"_n0"$n0"_n1"$n1"_220322
mkdir -p $outputdir 2>/dev/null
for f in `ls "$INPUT_DIR"/*`; do
  "$SCRIPT_DIR"/experiments/ex_time_2d.sh $f $m $r $n0 $n1 $outputdir
done