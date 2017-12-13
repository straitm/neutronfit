#!/bin/bash

d="results-`date -I`-$1"

# *Do* fail if the directory exists
if ! mkdir "$d"; then
  exit 1
fi

mv for_stage_two* savedhists* *.pdf *.out.txt "$d"
