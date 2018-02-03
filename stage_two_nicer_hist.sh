#!/bin/bash

cat $1 | ./make_stage_two_nicer_hists_ready.awk > $1.ready.C

root -b -q -l -n rhc_stage_two_nicer_hist.C'("'$1.ready.C'")'
