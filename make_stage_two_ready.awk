#!/usr/bin/awk -f

BEGIN { print "{"; }
/SetName/ {match($1, "\"(.*)\"", name)}
/SetPoint/ {match($1, ">(.*)", set); print "  " name[1] "-" set[0]}
/error_matrix/ {print $0}
END { print "}"; }

