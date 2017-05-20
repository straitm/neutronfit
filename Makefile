all: nm_summary_main.pdf nc_summary_main.pdf nm_slc_summary_main.pdf nc_slc_summary_main.pdf \
     nm_summary_muoncatcher.pdf nc_summary_muoncatcher.pdf \
     nm_slc_summary_muoncatcher.pdf nc_slc_summary_muoncatcher.pdf

define mindist_rule
fit_stage_two_mindist$(1)_nslc$(2)_$(3)_$(4).out.txt fit_stage_two_mindist$(1)_nslc$(2)_$(3)_$(4).pdf: \
  rhc_stage_two_C.so for_stage_two.C for_stage_two_ready_mindist$(1)_nslc$(2)_$(3)_$(4).C \
  common.C stage_two.sh
	./stage_two.sh $(1) $(2) $(3) $(4)

fit_stage_one_mindist$(1)_nslc$(2)_$(3)_$(4).out.txt for_stage_two_mindist$(1)_nslc$(2)_$(3)_$(4).C: \
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
MINSLCS := 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
MAXSLCS := 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
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
$(1)_slc_summary_$(2).pdf: nm_summary.C nm_summary.sh \
                fit_stage_two_mindist4_nslc2_4_$(2).out.txt \
                fit_stage_two_mindist4_nslc5_5_$(2).out.txt \
                fit_stage_two_mindist4_nslc6_6_$(2).out.txt \
                fit_stage_two_mindist4_nslc7_7_$(2).out.txt \
                fit_stage_two_mindist4_nslc8_12_$(2).out.txt
	./nm_summary.sh $(1)_slc 'fit_stage_two_mindist4_nslc*_$(2).out.txt'
$(1)_summary_$(2).pdf: nm_summary.C nm_summary.sh \
                fit_stage_two_mindist6_nslc0_20_$(2).out.txt \
		fit_stage_two_mindist5_nslc0_20_$(2).out.txt \
		fit_stage_two_mindist4_nslc0_20_$(2).out.txt \
		fit_stage_two_mindist3_nslc0_20_$(2).out.txt \
		fit_stage_two_mindist2_nslc0_20_$(2).out.txt \
		fit_stage_two_mindist1_nslc0_20_$(2).out.txt \
		fit_stage_two_mindist0_nslc0_20_$(2).out.txt
	./nm_summary.sh $(1) 'fit_stage_two_mindist?_nslc0_20_$(2).out.txt'
endef

REACTIONS := nm nc
$(foreach region, $(REGIONS), \
  $(foreach reaction, $(REACTIONS), \
    $(eval $(call summary_rule,$(reaction),$(region))) \
   ) \
 )

rhc_stage_zero_C.so: rhc_stage_zero.C
	./stage_zero.sh -1 0 0 main

rhc_stage_one_C.so: rhc_stage_one.C
	./stage_one.sh -1 0 0 main

rhc_stage_two_C.so: rhc_stage_two.C
	./stage_two.sh -1 0 0 main

clean:
	rm -f fit_stage_two_mindist*nslc*.out.txt \
              fit_stage_one_mindist*nslc*.pdf \
              for_stage_two*.C \
              n?_summary.pdf n?_slc_summary.pdf
