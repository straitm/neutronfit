#!/bin/bash

if ! root -l -b -q rhc_stage_zero.C+O'('$1','$2','$3',"'$4'")'; then
  rm savedhists_mindist${1}_nslc${2}_${3}_${4}.C
  exit 1
fi
