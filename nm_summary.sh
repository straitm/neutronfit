#!/bin/bash

grep -ih "${1:0:2} scale" $2 | \
  awk '{print $4, $6, $8, $10, $12, $14}' | \
  root -l -b -q nm_summary.C'("'$1'")';
