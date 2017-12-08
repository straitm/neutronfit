#!/bin/bash

# Fragile. Must match enum in common.C
if [ "$5" == TWOD ]; then
  fiftharg=0
elif [ "$5" == THREED ]; then
  fiftharg=1
else
  echo $0: did not understand fifth argument
  exit 1
fi

if ! root -n -l -b -q rhc_stage_zero.C+O'('$1','$2','$3',"'$4'", '$fiftharg')'; then
  rm savedhists_mindist${1}_nslc${2}_${3}_${4}_${5}.C
  exit 1
fi
