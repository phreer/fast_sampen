#!/bin/bash
build/bin/fast_sampen -m 4 -n 1048576 -r 0.15 \
  --input-format multirecord \
  --input data.PhysioNet/mit-bih-long-term-ecg-database-1.0.0/14046.txt  \
  --sample-size 2000 --sample-num 150 \
  --line-offset 100000 --random \
  --variance --n-computation 50 \
  --swr \
  --output-level 1
#  --kdtree-sample -rkd --simple-kdtree \
