#!/bin/bash

name=$1
dim=$2
shift 2

cat $@ | \
grep Probdensity | \
awk '{print $2, $3}' > stage_four.$name.$dim.tmp
root -b -q rhc_stage_four.C+'("'stage_four.$name.$dim.tmp'")' | tee stage_four.$name.$dim.out.txt
