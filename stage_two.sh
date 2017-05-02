#!/bin/bash

root -b -q rhc_stage_two.C++'("for_stage_two_ready_mindist'$1'.C", '$1')' | \
  tee fit_stage_two_mindist$1.out.txt
