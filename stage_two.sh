#!/bin/bash

cat for_stage_two_mindist$1.C | ./make_stage_two_ready.awk > for_stage_two_ready_mindist$1.C
