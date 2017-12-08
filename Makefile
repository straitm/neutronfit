all: nm_slc_summary_mindist1_main_TWOD.pdf          \
     nc_slc_summary_mindist1_main_TWOD.pdf          \
     nm_slc_summary_mindist1_muoncatcher_TWOD.pdf   \
     nc_slc_summary_mindist1_muoncatcher_TWOD.pdf   \
     nm_slc_summary_mindist1_main_THREED.pdf        \
     nc_slc_summary_mindist1_main_THREED.pdf        \
     nm_slc_summary_mindist1_muoncatcher_THREED.pdf \
     nc_slc_summary_mindist1_muoncatcher_THREED.pdf \
     nm_slc_summary_mindist3_main_TWOD.pdf          \
     nc_slc_summary_mindist3_main_TWOD.pdf          \
     nm_slc_summary_mindist3_muoncatcher_TWOD.pdf   \
     nc_slc_summary_mindist3_muoncatcher_TWOD.pdf   \
     nm_slc_summary_mindist3_main_THREED.pdf        \
     nc_slc_summary_mindist3_main_THREED.pdf        \
     nm_slc_summary_mindist3_muoncatcher_THREED.pdf \
     nc_slc_summary_mindist3_muoncatcher_THREED.pdf \
     nm_slc_summary_mindist6_main_TWOD.pdf          \
     nc_slc_summary_mindist6_main_TWOD.pdf          \
     nm_slc_summary_mindist6_muoncatcher_TWOD.pdf   \
     nc_slc_summary_mindist6_muoncatcher_TWOD.pdf   \
     nm_slc_summary_mindist6_main_THREED.pdf        \
     nc_slc_summary_mindist6_main_THREED.pdf        \
     nm_slc_summary_mindist6_muoncatcher_THREED.pdf \
     nc_slc_summary_mindist6_muoncatcher_THREED.pdf

define mindist_rule
fit_stage_two_mindist$(1)_nslc$(2)_$(3)_$(4)_$(5).out.txt fit_stage_two_mindist$(1)_nslc$(2)_$(3)_$(4)_$(5).pdf: \
  rhc_stage_two_C.so for_stage_two_ready_mindist$(1)_nslc$(2)_$(3)_$(4)_$(5).C \
  common.C stage_two.sh
	./stage_two.sh $(1) $(2) $(3) $(4) $(5)

fit_stage_one_mindist$(1)_nslc$(2)_$(3)_$(4)_$(5).pdf fit_stage_one_mindist$(1)_nslc$(2)_$(3)_$(4)_$(5).out.txt for_stage_two_mindist$(1)_nslc$(2)_$(3)_$(4)_$(5).C: \
  rhc_stage_one_C.so savedhists_mindist$(1)_nslc$(2)_$(3)_$(4)_$(5).C common.C stage_one.sh
	./stage_one.sh $(1) $(2) $(3) $(4) $(5)

for_stage_two_ready_mindist$(1)_nslc$(2)_$(3)_$(4)_$(5).C: \
  for_stage_two_mindist$(1)_nslc$(2)_$(3)_$(4)_$(5).C make_stage_two_ready.awk
	cat for_stage_two_mindist$(1)_nslc$(2)_$(3)_$(4)_$(5).C | ./make_stage_two_ready.awk \
          > for_stage_two_ready_mindist$(1)_nslc$(2)_$(3)_$(4)_$(5).C

savedhists_mindist$(1)_nslc$(2)_$(3)_$(4)_$(5).C: rhc_stage_zero_C.so common.C stage_zero.sh
	./stage_zero.sh $(1) $(2) $(3) $(4) $(5)
endef

REGIONS := main muoncatcher

CUTTYPES := TWOD THREED
MINDISTS := 0 1 2 3 4 5 6
# 99.98% of primary contained tracks are in events with <= 20 slices
# Must match the cuts in the summary_rule below
MINSLCS := 0.0 2.1 2.8 3.5 4.5 5.3 7.0 20.0
MAXSLCS := 0.0 2.1 2.8 3.5 4.5 5.3 7.0 20.0
$(foreach region, $(REGIONS), \
  $(foreach minslc, $(MINSLCS), \
    $(foreach maxslc, $(MAXSLCS), \
      $(foreach mindist, $(MINDISTS), \
        $(foreach cuttype, $(CUTTYPES), \
          $(eval $(call mindist_rule,$(mindist),$(minslc),$(maxslc),$(region),$(cuttype))) \
         ) \
       ) \
     ) \
   ) \
 )

define summary_rule
$(1)_slc_summary_mindist$(3)_$(2)_$(4).pdf: common.C util.C rhc_stage_three.C stage_three.sh \
                fit_stage_two_mindist$(3)_nslc0.0_20.0_$(2)_$(4).out.txt \
                fit_stage_two_mindist$(3)_nslc0.0_2.1_$(2)_$(4).out.txt \
                fit_stage_two_mindist$(3)_nslc2.1_2.8_$(2)_$(4).out.txt \
                fit_stage_two_mindist$(3)_nslc2.8_3.5_$(2)_$(4).out.txt \
                fit_stage_two_mindist$(3)_nslc3.5_4.5_$(2)_$(4).out.txt \
                fit_stage_two_mindist$(3)_nslc4.5_5.3_$(2)_$(4).out.txt \
                fit_stage_two_mindist$(3)_nslc5.3_7.0_$(2)_$(4).out.txt
	./stage_three.sh $(1)_slc $(2) \
          fit_stage_two_mindist$(3)_nslc{0.0_20.0,0.0_2.1,2.1_2.8,2.8_3.5,3.5_4.5,4.5_5.3,5.3_7.0}_$(2)_$(4).out.txt
endef

REACTIONS := nm nc
$(foreach region, $(REGIONS), \
    $(foreach reaction, $(REACTIONS), \
      $(foreach mindist, $(MINDISTS), \
        $(foreach cuttype, $(CUTTYPES), \
          $(eval $(call summary_rule,$(reaction),$(region),$(mindist),$(cuttype))) \
         ) \
       ) \
     ) \
 )

rhc_stage_zero_C.so: rhc_stage_zero.C common.C util.C
	./stage_zero.sh -1 0 0 main TWOD

rhc_stage_one_C.so: rhc_stage_one.C common.C util.C
	./stage_one.sh -1 0 0 main TWOD

rhc_stage_two_C.so: rhc_stage_two.C common.C util.C
	./stage_two.sh -1 0 0 main TWOD

output = fit_stage_*.out.txt \
         fit_stage_one_mindist*nslc*.pdf \
         fit_stage_two*pdf \
         for_stage_two*C \
         savedhists_*.C \
         for_stage_two*.C \
         n?_{mindist_,slc_}summary_mindist*_{main,muoncatcher}_*.pdf \
         *.out.txt

clean:
	rm -f $(output) *.so gmon.out *_C.d
