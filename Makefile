all: fit_stage_two_mindist6.pdf fit_stage_two_mindist5.pdf \
     fit_stage_two_mindist4.pdf fit_stage_two_mindist3.pdf \
     fit_stage_two_mindist2.pdf fit_stage_two_mindist1.pdf \
     fit_stage_two_mindist0.pdf

fit_stage_two_mindist6.pdf: \
  rhc_stage_two.C for_stage_two.C for_stage_two_ready_mindist6.C common.C Makefile
	root -b -q rhc_stage_two.C++'("for_stage_two_ready_mindist6.C", 6)' | tee fit_stage_two_mindist6.out.txt

for_stage_two_ready_mindist6.C: for_stage_two_mindist6.C make_stage_two_ready.awk Makefile
	cat for_stage_two_mindist6.C | ./make_stage_two_ready.awk > for_stage_two_ready_mindist6.C

for_stage_two_mindist6.C: rhc.C savedhists_mindist6.C common.C Makefile
	root -b -q rhc.C+O'("savedhists_mindist6.C", 6)' | tee fit_stage_one_mindist6.out.txt

savedhists_mindist6.C: rhc_stage_zero.C common.C Makefile
	root -b -q rhc_stage_zero.C+O'(6)'

# Maybe I should learn how to write Makefile rules better?

fit_stage_two_mindist5.pdf: \
  rhc_stage_two.C for_stage_two.C for_stage_two_ready_mindist5.C common.C Makefile
	root -b -q rhc_stage_two.C++'("for_stage_two_ready_mindist5.C", 5)' | tee fit_stage_two_mindist5.out.txt

for_stage_two_ready_mindist5.C: for_stage_two_mindist5.C make_stage_two_ready.awk Makefile
	cat for_stage_two_mindist5.C | ./make_stage_two_ready.awk > for_stage_two_ready_mindist5.C

for_stage_two_mindist5.C: rhc.C savedhists_mindist5.C common.C Makefile
	root -b -q rhc.C+O'("savedhists_mindist5.C", 5)' | tee fit_stage_one_mindist5.out.txt

savedhists_mindist5.C: rhc_stage_zero.C common.C Makefile
	root -b -q rhc_stage_zero.C+O'(5)'



fit_stage_two_mindist4.pdf: \
  rhc_stage_two.C for_stage_two.C for_stage_two_ready_mindist4.C common.C Makefile
	root -b -q rhc_stage_two.C++'("for_stage_two_ready_mindist4.C", 4)' | tee fit_stage_two_mindist4.out.txt

for_stage_two_ready_mindist4.C: for_stage_two_mindist4.C make_stage_two_ready.awk Makefile
	cat for_stage_two_mindist4.C | ./make_stage_two_ready.awk > for_stage_two_ready_mindist4.C

for_stage_two_mindist4.C: rhc.C savedhists_mindist4.C common.C Makefile
	root -b -q rhc.C+O'("savedhists_mindist4.C", 4)' | tee fit_stage_one_mindist4.out.txt

savedhists_mindist4.C: rhc_stage_zero.C common.C Makefile
	root -b -q rhc_stage_zero.C+O'(4)'



fit_stage_two_mindist3.pdf: \
  rhc_stage_two.C for_stage_two.C for_stage_two_ready_mindist3.C common.C Makefile
	root -b -q rhc_stage_two.C++'("for_stage_two_ready_mindist3.C", 3)' | tee fit_stage_two_mindist3.out.txt

for_stage_two_ready_mindist3.C: for_stage_two_mindist3.C make_stage_two_ready.awk Makefile
	cat for_stage_two_mindist3.C | ./make_stage_two_ready.awk > for_stage_two_ready_mindist3.C

for_stage_two_mindist3.C: rhc.C savedhists_mindist3.C common.C Makefile
	root -b -q rhc.C+O'("savedhists_mindist3.C", 3)' | tee fit_stage_one_mindist3.out.txt

savedhists_mindist3.C: rhc_stage_zero.C common.C Makefile
	root -b -q rhc_stage_zero.C+O'(3)'



fit_stage_two_mindist2.pdf: \
  rhc_stage_two.C for_stage_two.C for_stage_two_ready_mindist2.C common.C Makefile
	root -b -q rhc_stage_two.C++'("for_stage_two_ready_mindist2.C", 2)' | tee fit_stage_two_mindist2.out.txt

for_stage_two_ready_mindist2.C: for_stage_two_mindist2.C make_stage_two_ready.awk Makefile
	cat for_stage_two_mindist2.C | ./make_stage_two_ready.awk > for_stage_two_ready_mindist2.C

for_stage_two_mindist2.C: rhc.C savedhists_mindist2.C common.C Makefile
	root -b -q rhc.C+O'("savedhists_mindist2.C", 2)' | tee fit_stage_one_mindist2.out.txt

savedhists_mindist2.C: rhc_stage_zero.C common.C Makefile
	root -b -q rhc_stage_zero.C+O'(2)'



fit_stage_two_mindist1.pdf: \
  rhc_stage_two.C for_stage_two.C for_stage_two_ready_mindist1.C common.C Makefile
	root -b -q rhc_stage_two.C++'("for_stage_two_ready_mindist1.C", 1)' | tee fit_stage_two_mindist1.out.txt

for_stage_two_ready_mindist1.C: for_stage_two_mindist1.C make_stage_two_ready.awk Makefile
	cat for_stage_two_mindist1.C | ./make_stage_two_ready.awk > for_stage_two_ready_mindist1.C

for_stage_two_mindist1.C: rhc.C savedhists_mindist1.C common.C Makefile
	root -b -q rhc.C+O'("savedhists_mindist1.C", 1)' | tee fit_stage_one_mindist1.out.txt

savedhists_mindist1.C: rhc_stage_zero.C common.C Makefile
	root -b -q rhc_stage_zero.C+O'(1)'


fit_stage_two_mindist0.pdf: \
  rhc_stage_two.C for_stage_two.C for_stage_two_ready_mindist0.C common.C Makefile
	root -b -q rhc_stage_two.C++'("for_stage_two_ready_mindist0.C", 0)' | tee fit_stage_two_mindist0.out.txt

for_stage_two_ready_mindist0.C: for_stage_two_mindist0.C make_stage_two_ready.awk Makefile
	cat for_stage_two_mindist0.C | ./make_stage_two_ready.awk > for_stage_two_ready_mindist0.C

for_stage_two_mindist0.C: rhc.C savedhists_mindist0.C common.C Makefile
	root -b -q rhc.C+O'("savedhists_mindist0.C", 0)' | tee fit_stage_one_mindist0.out.txt

savedhists_mindist0.C: rhc_stage_zero.C common.C Makefile
	root -b -q rhc_stage_zero.C+O'(0)'

