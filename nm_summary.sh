#!/bin/bash

grep -ih "$1 scale" fit_stage_two_mindist?_nslc*.out.txt | awk '{print $4, $10, $12, $14}' | root -l -b -q nm_summary.C'("'$1'")';
