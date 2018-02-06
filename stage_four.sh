#!/bin/bash

# Make sure we fail if root fails at the beginning of a pipeline
set -o pipefail

if [ $1 == compile ]; then
  root -b -q rhc_stage_four.C++'(NULL, false)'
  exit
fi

name=$1
dim=$2
shift 2


cat $@ | \
grep Probdensity | \
awk '{print $2, $3}' > stage_four.$name.$dim.tmp
root -n -l -b -q rhc_stage_four.C+'("'stage_four.$name.$dim.'",'$(if [ $name == nm ]; then echo true; else echo false; fi)')' | \
  tee stage_four.$name.$dim.out.txt
