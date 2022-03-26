#!/usr/bin/env bash

inputdir=$HOME/workspace/entropy/additional_data/uurgeg
outputdir=$inputdir

for f in `ls "$inputdir"/*.txt`; do
  output_filename=${f%.*}.csv
  sed -n '32,$p' $f | sed 's/ //g' | sed -e 's/^.//' | sed -e '2d' > $output_filename
done