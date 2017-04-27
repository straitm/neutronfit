all: fit_stage_two.pdf fit.pdf

fit_stage_two.out.txt fit_stage_two.pdf: rhc_stage_two.C for_stage_two.C for_stage_two_ready.C common.C
	root -b -q rhc_stage_two.C++'("for_stage_two_ready.C")' | tee fit_stage_two.out.txt

for_stage_two_ready.C: for_stage_two.C make_stage_two_ready.awk
	cat for_stage_two.C | ./make_stage_two_ready.awk > for_stage_two_ready.C

fit_stage_one.out.txt fit.pdf for_stage_two.C: rhc.C savedhists.C common.C
	root -b -q rhc.C+O'("savedhists.C")' | tee fit_stage_one.out.txt

savedhists.C: rhc.C common.C
	root -b -q rhc.C+O
