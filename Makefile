all: nm_mindist_summary_main.pdf        nc_mindist_summary_main.pdf \
     nm_mindist_summary_muoncatcher.pdf nc_mindist_summary_muoncatcher.pdf \
     nm_slc_summary_main.pdf            nc_slc_summary_main.pdf \
     nm_slc_summary_muoncatcher.pdf     nc_slc_summary_muoncatcher.pdf

define mindist_rule
fit_stage_two_mindist$(1)_nslc$(2)_$(3)_$(4).out.txt fit_stage_two_mindist$(1)_nslc$(2)_$(3)_$(4).pdf: \
  rhc_stage_two_C.so for_stage_two_ready_mindist$(1)_nslc$(2)_$(3)_$(4).C \
  common.C stage_two.sh
	./stage_two.sh $(1) $(2) $(3) $(4)

fit_stage_one_mindist$(1)_nslc$(2)_$(3)_$(4).pdf fit_stage_one_mindist$(1)_nslc$(2)_$(3)_$(4).out.txt for_stage_two_mindist$(1)_nslc$(2)_$(3)_$(4).C: \
  rhc_stage_one_C.so savedhists_mindist$(1)_nslc$(2)_$(3)_$(4).C common.C stage_one.sh
	./stage_one.sh $(1) $(2) $(3) $(4)

for_stage_two_ready_mindist$(1)_nslc$(2)_$(3)_$(4).C: \
  for_stage_two_mindist$(1)_nslc$(2)_$(3)_$(4).C make_stage_two_ready.awk
	cat for_stage_two_mindist$(1)_nslc$(2)_$(3)_$(4).C | ./make_stage_two_ready.awk \
          > for_stage_two_ready_mindist$(1)_nslc$(2)_$(3)_$(4).C

savedhists_mindist$(1)_nslc$(2)_$(3)_$(4).C: rhc_stage_zero_C.so common.C stage_zero.sh
	./stage_zero.sh $(1) $(2) $(3) $(4)
endef

REGIONS := main muoncatcher

MINDISTS := 0 1 2 3 4 5 6
MINSLCS := 0.0 1.0 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 16.0 17.0 18.0 19.0 20.0
MAXSLCS := 0.0 1.0 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 16.0 17.0 18.0 19.0 20.0
$(foreach region, $(REGIONS), \
  $(foreach minslc, $(MINSLCS), \
    $(foreach maxslc, $(MAXSLCS), \
      $(foreach mindist, $(MINDISTS), \
	$(eval $(call mindist_rule,$(mindist),$(minslc),$(maxslc),$(region))) \
       ) \
     ) \
   ) \
 )

define summary_rule
# 99.98% of primary contained tracks are in events with <= 20 slices
$(1)_slc_summary_$(2).pdf: common.C rhc_stage_three.C stage_three.sh \
                fit_stage_two_mindist3_nslc0.0_20.0_$(2).out.txt \
                fit_stage_two_mindist3_nslc0.0_2.8_$(2).out.txt \
                fit_stage_two_mindist3_nslc2.8_3.5_$(2).out.txt \
                fit_stage_two_mindist3_nslc3.5_4.5_$(2).out.txt \
                fit_stage_two_mindist3_nslc4.5_5.5_$(2).out.txt \
                fit_stage_two_mindist3_nslc5.5_7.0_$(2).out.txt
	./stage_three.sh $(1)_slc $(2) \
          fit_stage_two_mindist3_nslc{0.0_20.0,0.0_2.8,2.8_3.5,3.5_4.5,4.5_5.5,5.5_7.0}_$(2).out.txt
$(1)_mindist_summary_$(2).pdf: common.C rhc_stage_three.C stage_three.sh \
                fit_stage_two_mindist6_nslc0_10_$(2).out.txt \
		fit_stage_two_mindist5_nslc0_10_$(2).out.txt \
		fit_stage_two_mindist4_nslc0_10_$(2).out.txt \
		fit_stage_two_mindist3_nslc0_10_$(2).out.txt \
		fit_stage_two_mindist2_nslc0_10_$(2).out.txt \
		fit_stage_two_mindist1_nslc0_10_$(2).out.txt \
		fit_stage_two_mindist0_nslc0_10_$(2).out.txt
	./stage_three.sh $(1)_mindist $(2) \
          fit_stage_two_mindist?_nslc0_10_$(2).out.txt
endef

REACTIONS := nm nc
$(foreach region, $(REGIONS), \
  $(foreach reaction, $(REACTIONS), \
    $(eval $(call summary_rule,$(reaction),$(region))) \
   ) \
 )

rhc_stage_zero_C.so: rhc_stage_zero.C common.C
	./stage_zero.sh -1 0 0 main

rhc_stage_one_C.so: rhc_stage_one.C common.C
	./stage_one.sh -1 0 0 main

rhc_stage_two_C.so: rhc_stage_two.C common.C
	./stage_two.sh -1 0 0 main

output = fit_stage_*.out.txt \
         fit_stage_one_mindist*nslc*.pdf \
         fit_stage_two*pdf \
         for_stage_two*C \
         savedhists_*.C \
         for_stage_two*.C \
         n?_{mindist_,slc_}summary_{main,muoncatcher}.pdf \
         stage_three.*.*.out.txt

clean:
	rm -f $(output) *.so gmon.out *_C.d
