#!/bin/bash

grep -ihE "${1:0:2} (muoncatcher|main)" $3 | \
  awk '{print $4, $6, $8, $10, $12, $14}' | \
  root -l -b -q nm_summary.C'("'$1'", "'$2'")';
