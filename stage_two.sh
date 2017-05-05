#!/bin/bash

if [ $1 -ge 0 ]; then
  root -b -q rhc_stage_two.C+'("for_stage_two_ready_mindist'$1'.C", '$1')' | \
    tee fit_stage_two_mindist$1.out.txt
else
  root -b -q rhc_stage_two.C+'(NULL, '$1')'
fi
