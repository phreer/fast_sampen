#!/bin/bash
build/bin/fast_sampen -m 5 -n 2005 -r 0.15 --input-format multirecord \
  --input data.PhysioNet/chfdb/chf01.txt -rkd --line-offset 100000 \
  --kd-sample --sample-size 1000 --sample-num 1
