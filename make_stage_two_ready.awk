#!/usr/bin/awk -f

BEGIN { print "{"; }
/SetName/ {match($1, "\"(.*)\"", name)}
/SetPoint/ {match($1, ">(.*)", set); print "  " name[1] "-" set[0]}
END { print "}"; }

