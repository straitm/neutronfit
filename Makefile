all: \
     fit_stage_two_mindist0.pdf \
     fit_stage_two_mindist1.pdf \
     fit_stage_two_mindist2.pdf \
     fit_stage_two_mindist3.pdf \
     fit_stage_two_mindist4.pdf \
     fit_stage_two_mindist5.pdf \
     fit_stage_two_mindist6.pdf

define mindist_rule
fit_stage_two_mindist$(1).out.txt fit_stage_two_mindist$(1).pdf: \
  rhc_stage_two_C.so for_stage_two.C for_stage_two_ready_mindist$(1).C \
  common.C stage_two.sh
	./stage_two.sh $(1)

fit_stage_one_mindist$(1).out.txt for_stage_two_mindist$(1).C: rhc_stage_one_C.so savedhists_mindist$(1).C common.C stage_one.sh
	./stage_one.sh $(1)

for_stage_two_ready_mindist$(1).C: for_stage_two_mindist$(1).C make_stage_two_ready.awk Makefile
	cat for_stage_two_mindist$(1).C | ./make_stage_two_ready.awk > for_stage_two_ready_mindist$(1).C

savedhists_mindist$(1).C: rhc_stage_zero_C.so common.C stage_zero.sh
	./stage_zero.sh $(1)	
endef

MINDISTS := 0 1 2 3 4 5 6
$(foreach mindist, $(MINDISTS), $(eval $(call mindist_rule,$(mindist))))

rhc_stage_zero_C.so: rhc_stage_zero.C
	./stage_zero.sh -1

rhc_stage_one_C.so: rhc_stage_one.C
	./stage_one.sh -1

rhc_stage_two_C.so: rhc_stage_two.C
	./stage_two.sh -1
