#!/usr/bin/bash

r=0.2
m=2
INPUT_DIR=$HOME/workspace/entropy/2d/generation_dataset_kq_256
OUTPUT_DIR=$HOME/workspace/entropy/2d/kd_256_result_r"$r"_m"$m"

for folder in `ls "$INPUT_DIR"`; do
    echo "Processing folder $folder..."
    mkdir -p "$OUTPUT_DIR"/"$folder" 2>/dev/null
    for f in `ls "$INPUT_DIR"/"$folder"`; do
        echo "Processing file "$folder"/$f ..."
        build/bin/sampen2d -r $r -m $m --multiscale-depth 10\
            --multiscale-factor 0.8\
            --input "$INPUT_DIR"/"$folder"/"$f"\
            > "$OUTPUT_DIR"/"$folder"/"$f".txt
    done
done
