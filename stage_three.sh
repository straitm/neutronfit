#!/bin/bash

name=$1
region=$2
shift 2

# This makes us return failure if grep can't find the requested files
set -o pipefail

grep -ihE "${name:0:2} $region" $@ | \
  awk '{print $4, $6, $8, $10, $12, $14, $16}' | \
  root -n -l -b -q rhc_stage_three.C+'("'$name'", "'$region'")' | \
  tee stage_three.$name.$region.out.txt
