#!/bin/bash

# Make sure we fail if root fails at the beginning of a pipeline
set -o pipefail

if [ $1 -ge 0 ]; then
  out=fit_stage_one_mindist$1_nslc$2_$3_$4.out.txt
  if ! root -b -q rhc_stage_one.C+O'("savedhists_mindist'$1'_nslc'$2'_'$3'_'$4'.C", '$1','$2','$3',"'$4'")' | tee $out; then
    rm $out
    exit 1
  fi
else
  root -b -q rhc_stage_one.C+O'(NULL, '$1','$2','$3',"'$4'")'
fi
