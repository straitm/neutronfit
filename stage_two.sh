#!/bin/bash

if [ $1 -ge 0 ]; then
  root -l -b -q rhc_stage_two.C+'("for_stage_two_ready_mindist'$1'_nslc'$2'_'$3'.C",'$1','$2','$3')' | \
    tee fit_stage_two_mindist$1_nslc${2}_$3.out.txt
else
  root -l -b -q rhc_stage_two.C+'(NULL, '$1','$2','$3')'
fi
