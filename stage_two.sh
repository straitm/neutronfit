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

if [ $1 -ge 0 ]; then
  root -n -l -b -q rhc_stage_two.C+'("for_stage_two_ready_mindist'$1'_nslc'$2'_'$3'_'$4'_'$5'.C",'$1','$2','$3',"'$4'",'$fiftharg')' | \
    tee fit_stage_two_mindist$1_nslc$2_$3_$4_$5.out.txt
else
  root -n -l -b -q rhc_stage_two.C+'(NULL, '$1','$2','$3',"'$4'",'$fiftharg')'
fi
