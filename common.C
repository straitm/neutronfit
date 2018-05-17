// Whether we are using neutrons with 2D information alone (usually just
// one hit) or if we are requiring requiring hits in both x and y.
enum two_or_three_d { TWOD, THREED, MAX_TWO_OR_THREE_D };
const char * two_or_three_d_names[MAX_TWO_OR_THREE_D] = { "TWOD", "THREED" };

//#define BGSUB

//#define LOWINTENSITYSTUDY
//#define DOUBLERHCSTUDY

// Number of neutron pileup background samples per signal sample.  Set by my
// ntuple maker.  Set to zero to skip background subtraction.
const int bgmult =
#ifdef BGSUB
  1;
#else
  0;
#endif

// Everything is doubled because there is the signal and there is the
// off-space pileup background sample
const int SIG_AND_BG = bgmult?2:1;

// We aren't going to use anything between -1 and 2 microseconds because
// the detector conditions are just too awful.  And then don't use
// anything before holex_hi microseconds because I'm now using that region
// to define Michels that cut further interactions.
const double holex_lo = -1, holex_hi = 5;

const double nnegbins = 209;
const double maxrealtime = 269;
const double additional = 0;

// Restrict x and y to 100 and trkz_cutlo to 600
// as a blunt way to reduce pile-up subtraction problems.
const double trkx_cut = 170,
             trky_cut = 170,
             trkz_cutlo = 0,
             trkz_cuthi = 1250,
             mucatch_trky_cut = 45,
             mucatch_trkz_cutlo = 1310,
             mucatch_trkz_cuthi = 1560,
             mucatch_trkstartz_cut = 1100;
const double trklen_cut = 200;
const double remid_cut = 0.75;


// *** Parameters controlling how we deal with pileup neutrons from other ***
// *** neutrino interactions                                              ***

// How much weight to give the visible slices, as opposed to the intensity.
// i.e. What fraction of the pileup neutrons we think come from these
// rather than from interactions in the rock.
//
// From a Geant study using both FHC and RHC MC.  Fortunately, even though of
// course some difference is expected, it is very small.  For the main
// detector, FHC is 0.02 higher and for the muon catcher it is 0.03 higher.
// This is somewhat bigger than the statistical error (I only looked at one
// file each), but I think this can safely be considered within systematics
// given just how many poorly handled neutron processes there are that we are
// summing over.  Therefore, I've put the averages in.
//
// I used:
// neardet_genie_nonswap_genierw_fhc_v08_1162_r00012036_s20_c001_R17-03-01-
//   prod3reco.d_v1_20170322_204739_sim
// neardet_genie_nonswap_genierw_rhc_v08_1197_r00011654_s00_c001_R16-12-20-
//   prod3recopreview.b_v3_20161220_134502_sim
//
// The only reconstruction used in the study was slicing.  Otherwise, only
// truth information was used.
//
// The rock is a better neutron moderator than we give it credit for, since it
// certainly contains some hydrogen, but the model has none.  Therefore higher
// values of these numbers would probably be more accurate.  But since I don't
// know what order of magnitude this effect is, I am not attempting to adjust
// for it.
//
// (Note that the figure for 2D cut and the main detector is much lower than
// the others.  This is because the main detector is an excellent neutron
// moderator, so many captures are around its edges.  With only a 2D cut, many
// of these pileup neutrons from the rock are selected.)
//
// I have adjusted these numbers using the ratio of neutrons reaching the
// detector in a small MC run with the rock composition we have been using, and
// with a sensible composition that includes the water content documented in
// MINOS-doc-2777.  My MC is a little dubious because it does not simulate
// events very far away from the walls, and maybe not even in the ceiling, but
// the ratio should be more or less right.  To be updated.  The numbers without
// this adjustment were { 0.64, 0.80 }, { 0.76, 0.81 }.
const double npileup_sliceweight        [MAX_TWO_OR_THREE_D] = { 0.75, 0.86 };
const double npileup_sliceweight_mucatch[MAX_TWO_OR_THREE_D] = { 0.78, 0.83 };

// Number of slices per 10**12 POT in RHC and FHC.  The ratio between these
// sets the relative number of neutrons we think are coming out of the rock.
// The absolute values are unimportant since we can tune either them or
// npileup_sliceweight.
const double slc_per_twp_rhc = 0.075;
const double slc_per_twp_fhc = 0.172;

// Maximum fraction of energy in the event that can be from the muon.
// This is a crude way to get at resolution bins.  To allow any resolution
// bin, set to 1.0 (or to be safe, a bit more, since I use a really crude
// approximation for muon energy).  To only allow roughly the worst bin,
// use ~0.55.
const double max_frac_e_mu = 1.5;

#define SEPARATED_PERIODS

#ifdef SEPARATED_PERIODS
  const int nperiodrhc = 3; // 4, 6, 7
  const int nperiodfhc = 4; // 1, 2, 3, 5
#else
  const int nperiodrhc = 1;
  const int nperiodfhc = 1;
#endif
const int nperiod    = nperiodrhc + nperiodfhc;

const char * const inputfiles[nperiod] = {
#ifdef SEPARATED_PERIODS
  "../ndcosmic_data/201712-period4goodbad.root",
  "../ndcosmic_data/201712-period6goodbad.root",
  "../ndcosmic_data/201805-period7abcgoodbad.root",

  "../ndcosmic_data/201712-period1goodbad.root",
  "../ndcosmic_data/201712-period2goodbad.root",
  "../ndcosmic_data/201712-period3goodbad.root",
  "../ndcosmic_data/201712-period5goodbad.root"
#else
  "../ndcosmic_data/201712-period46goodbad_rhc.root",
  "../ndcosmic_data/201712-period1235goodbad_fhc.root"
#endif
};

const char * const Speriodnames[SIG_AND_BG*nperiod] = {
#ifdef SEPARATED_PERIODS
      "P4",   "P6",  "P7",  "P1",   "P2",   "P3",   "P5",
  #ifdef BGSUB
      "P4BG", "P6BG", "P7BG", "P1BG", "P2BG", "P3BG", "P5BG"
  #endif
#else
  "PR",
  "PF",
  #ifdef BGSUB
  "PRBG",
  "PFBG",
  #endif
#endif
};

const char * const Lperiodnames[SIG_AND_BG*nperiod] = {
#ifdef SEPARATED_PERIODS
     "Period 4 (RHC)",
     "Period 6 (RHC)",
     "Period 7 (RHC)",
     "Period 1 (FHC)",
     "Period 2 (FHC)",
     "Period 3 (FHC)",
     "Period 5 (FHC)",
  #ifdef BGSUB
     "Period 4 (RHC) Pileup",
     "Period 6 (RHC) Pileup",
     "Period 7 (RHC) Pileup",
     "Period 1 (FHC) Pileup",
     "Period 2 (FHC) Pileup",
     "Period 3 (FHC) Pileup",
     "Period 5 (FHC) Pileup",
  #endif
#else
     "RHC",
     "FHC",
  #ifdef BGSUB
     "RHC Pileup",
     "FHC Pileup",
  #endif
#endif
};

const int nbeam = 2; // not really generalizable as it stands

const int nbins_e[MAX_TWO_OR_THREE_D] = { 6, 3 };
const double bins_e[MAX_TWO_OR_THREE_D][6+1] = { // janky...
  { 0.5, 1.375, 2.250, 3.125, 4.0, 5.0, 6.0 },
  { 0.5, 1.5  , 3.0  , 6.0  ,  -1,  -1,  -1 }
};
// Other possibilities
// 5: 0.5, 1.6, 2.7, 3.8, 4.9, 6.0
// 8: 0.5, 1, 1.5, 2, 2.5, 3, 4.0, 5.0, 6.0

// From a loose cut running of stage one
const double nominal_neutron_lifetime_main        = 53.6; // us
const double nominal_neutron_lifetime_muoncatcher = 59.4;
