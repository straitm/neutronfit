#!/bin/bash

# Make sure we fail if root fails at the beginning of a pipeline
set -o pipefail

# Fragile. Must match enum in common.C
if [ "$5" == TWOD ]; then
  fiftharg=0
elif [ "$5" == THREED ]; then
  fiftharg=1
else
  echo $0: did not understand fifth argument
  exit 1
fi

if [ $1 -ge 0 ]; then
  out=fit_stage_one_mindist$1_nslc$2_$3_$4_$5.out.txt
  if ! root -n -b -q rhc_stage_one.C+O'("savedhists_mindist'$1'_nslc'$2'_'$3'_'$4'_'$5'.C", '$1','$2','$3',"'$4'", '$fiftharg')' | tee $out; then
    rm $out
    exit 1
  fi
else
  root -n -b -q rhc_stage_one.C+O'(NULL, '$1','$2','$3',"'$4'",'$fiftharg')'
fi
