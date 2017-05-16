all: nm_summary.pdf nc_summary.pdf

define mindist_rule
fit_stage_two_mindist$(1)_nslc$(2)_$(3).out.txt fit_stage_two_mindist$(1)_nslc$(2)_$(3).pdf: \
  rhc_stage_two_C.so for_stage_two.C for_stage_two_ready_mindist$(1)_nslc$(2)_$(3).C \
  common.C stage_two.sh
	./stage_two.sh $(1) $(2) $(3)

fit_stage_one_mindist$(1)_nslc$(2)_$(3).out.txt for_stage_two_mindist$(1)_nslc$(2)_$(3).C: \
  rhc_stage_one_C.so savedhists_mindist$(1)_nslc$(2)_$(3).C common.C stage_one.sh
	./stage_one.sh $(1) $(2) $(3)

for_stage_two_ready_mindist$(1)_nslc$(2)_$(3).C: \
  for_stage_two_mindist$(1)_nslc$(2)_$(3).C make_stage_two_ready.awk
	cat for_stage_two_mindist$(1)_nslc$(2)_$(3).C | ./make_stage_two_ready.awk \
          > for_stage_two_ready_mindist$(1)_nslc$(2)_$(3).C

savedhists_mindist$(1)_nslc$(2)_$(3).C: rhc_stage_zero_C.so common.C stage_zero.sh
	./stage_zero.sh $(1) $(2) $(3)
endef

MINDISTS := 0 1 2 3 4 5 6
$(foreach mindist, $(MINDISTS), $(eval $(call mindist_rule,$(mindist),0,20)))

rhc_stage_zero_C.so: rhc_stage_zero.C
	./stage_zero.sh -1 0 0

rhc_stage_one_C.so: rhc_stage_one.C
	./stage_one.sh -1 0 0

rhc_stage_two_C.so: rhc_stage_two.C
	./stage_two.sh -1 0 0

nm_summary.pdf: nm_summary.C nm_summary.sh \
                fit_stage_two_mindist6_nslc0_20.out.txt \
		fit_stage_two_mindist5_nslc0_20.out.txt \
		fit_stage_two_mindist4_nslc0_20.out.txt \
		fit_stage_two_mindist3_nslc0_20.out.txt \
		fit_stage_two_mindist2_nslc0_20.out.txt \
		fit_stage_two_mindist1_nslc0_20.out.txt \
		fit_stage_two_mindist0_nslc0_20.out.txt
	./nm_summary.sh nm

nc_summary.pdf: nm_summary.C nm_summary.sh \
                fit_stage_two_mindist6_nslc0_20.out.txt \
		fit_stage_two_mindist5_nslc0_20.out.txt \
		fit_stage_two_mindist4_nslc0_20.out.txt \
		fit_stage_two_mindist3_nslc0_20.out.txt \
		fit_stage_two_mindist2_nslc0_20.out.txt \
		fit_stage_two_mindist1_nslc0_20.out.txt \
		fit_stage_two_mindist0_nslc0_20.out.txt
	./nm_summary.sh nc
