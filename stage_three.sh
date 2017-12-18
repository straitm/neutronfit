#!/bin/bash

name=$1
region=$2
ndim=$3
shift 3

# This makes us return failure if grep can't find the requested files
set -o pipefail

if [ $name == compile ]; then
  root -n -l -b -q rhc_stage_three.C+'("'$name'", "'$region'")'
else
  grep -ihE "${name:0:2} $region" $@ | \
    awk '{print $4, $6, $8, $10, $12, $14, $16, $17}' | \
    root -n -l -b -q rhc_stage_three.C+'("'$name'", "'$region'")' | \
    tee stage_three.$name.$region.$ndim.out.txt
fi
