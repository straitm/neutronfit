#!/bin/bash

if [ $1 -ge 0 ]; then
  root -b -q rhc_stage_one.C+O'("savedhists_mindist'$1'.C", '$1')' | tee fit_stage_one_mindist$1.out.txt
else
  root -b -q rhc_stage_one.C+O'(NULL, '$1')'
fi
