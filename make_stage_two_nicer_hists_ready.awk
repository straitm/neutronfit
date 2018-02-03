#!/usr/bin/awk -f

# Lordy. Fix up the macro so it can be run with gROOT::Macro(). For each
# TGraph, find its name and make all the subsequent calls use it instead
# of a meaningless variable.
/SetName/ {match($1, "\"(.*)\"", name)}

/SetPoint/ {match($1, ">(.*)", set); print "  " name[1] "-" set[0]}

!/SetName|SetPoint/ { print $0; }
