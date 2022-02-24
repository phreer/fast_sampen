#!/bin/bash
build/bin/fast_sampen -m 4 -n 200000 -r 0.15 --input-format multirecord \
  --input data.PhysioNet/chfdb/chf01.txt -rkd --line-offset 100000