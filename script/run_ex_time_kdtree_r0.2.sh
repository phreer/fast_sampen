#!/usr/bin/env bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

if [ $# -ne 1 ]; then
  echo "Usage: $0 <tag>" >&2
  exit -1
fi

tag=$1

for m in `seq 2 4`; do
  ${SCRIPT_DIR}/experiments/ex_time_kdtree.sh 0.2 $m $tag &
done