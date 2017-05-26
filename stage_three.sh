#!/bin/bash

name=$1
region=$2
shift 2

grep -ihE "${name:0:2} $region" $@ | \
  awk '{print $4, $6, $8, $10, $12, $14}' | \
  root -l -b -q rhc_stage_three.C'("'$name'", "'$region'")';
