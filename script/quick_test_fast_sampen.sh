#!/bin/bash
build/bin/fast_sampen -m 3 -n 100003 -r 0.15 \
  --input-format multirecord \
  --input data.PhysioNet/chfdb/chf01.txt  \
  --sample-size 2000 --sample-num 1 \
  --line-offset 100000 \
  --kdtree-sample -rkd --simple-kdtree
