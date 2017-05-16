#!/bin/bash

if [ $1 -ge 0 ]; then
  root -b -q rhc_stage_one.C+O'("savedhists_mindist'$1'_nslc'$2'_'$3'.C", '$1','$2','$3')' \
    | tee fit_stage_one_mindist$1_nslc$2_$3.out.txt
else
  root -b -q rhc_stage_one.C+O'(NULL, '$1','$2','$3')'
fi
