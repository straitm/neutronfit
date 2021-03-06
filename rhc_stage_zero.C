#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TROOT.h"
#include "TGraphAsymmErrors.h"
#include <fstream>
using std::string;

#include "TMinuit.h"
TMinuit * mn = NULL; // dumb, because of common.C

#include "common.C"
#include "util.C"

bool muoncatcher = true;// set at entry point

struct data{
  int run;
  float slce;
  int type;

  int i;
  int primary;
  float timeleft;
  float timeback;
  float remid;
  float trklen;
  int contained;
  float trkx;
  float trky;
  float trkz;
  float trkstartz;
  int nslc;
  float pot;

  int event;
  int trk;
  float t;
  float mindist;
  float e;
  float pe;

  int nhitx;
  int nhity;

  int subrun;
  int nhit;
  float maxdist;
  float cosx;
  float cosy;
};

enum branchstat { FOR2D, FOR1D };

static void setbranchaddresses(data * dat, TTree * t,
                               const branchstat how,
                               const two_or_three_d cut_dimensions)
{
  t->SetBranchStatus("*", 0);
  setbranchaddress("run", &dat->run, t);
  setbranchaddress("slce", &dat->slce, t);
  setbranchaddress("type", &dat->type, t);

  setbranchaddress("i", &dat->i, t);
  setbranchaddress("primary", &dat->primary, t);
  setbranchaddress("timeleft", &dat->timeleft, t);
  setbranchaddress("timeback", &dat->timeback, t);
  setbranchaddress("remid", &dat->remid, t);
  setbranchaddress("trklen", &dat->trklen, t);
  setbranchaddress("contained", &dat->contained, t);
  setbranchaddress("trkx", &dat->trkx, t);
  setbranchaddress("trky", &dat->trky, t);
  setbranchaddress("trkz", &dat->trkz, t);
  setbranchaddress("trkstartz", &dat->trkstartz, t);
  setbranchaddress("nslc", &dat->nslc, t);
  setbranchaddress("pot", &dat->pot, t);

  if(how == FOR2D){
    setbranchaddress("event", &dat->event, t);
    setbranchaddress("trk", &dat->trk, t);
    setbranchaddress("t", &dat->t, t);
    setbranchaddress("mindist", &dat->mindist, t);
    setbranchaddress("e", &dat->e, t);
    setbranchaddress("pe", &dat->pe, t);
    if(cut_dimensions == THREED){
      setbranchaddress("nhitx", &dat->nhitx, t);
      setbranchaddress("nhity", &dat->nhity, t);
    }
  }

  //setbranchaddress("subrun", &dat->subrun, t);
  //setbranchaddress("nhit", &dat->nhit, t);
  //setbranchaddress("maxdist", &dat->maxdist, t);
  //setbranchaddress("cosx", &dat->cosx, t);
  //setbranchaddress("cosy", &dat->cosy, t);
}

// Continue calling the cuts minslc and maxslc, but actually translate
// them into an intensity in 10**12 POTs ("twp") so that we are taking
// into account neutrons from invisible rock events too.
static bool pass_intensity(data * dat, const float minslc,
                           const float maxslc, const bool rhc,
                           const two_or_three_d cut_dimensions)
{
  // This may not agree with other conversions, but since what I actually
  // want is neutrons produced per POT --- something no one knows ---
  // it is necessarily approximate.
  const float slc_per_twp = rhc? slc_per_twp_rhc: slc_per_twp_fhc;

  // Less one for the noise slice and less another for the slice
  // we're considering
  const int other_physics_nslc = dat->nslc - 2;

  const double slicew = muoncatcher? npileup_sliceweight_mucatch[cut_dimensions]:
                                     npileup_sliceweight[cut_dimensions];

  const float eff_slc = dat->pot * slc_per_twp * (1-slicew)
                          + other_physics_nslc *    slicew;

  return eff_slc >= minslc && eff_slc < maxslc;
}

// true if it passes the cut
static bool track_itself_cut(data * dat, const float minslc,
                             const float maxslc, const bool rhc,
                             const two_or_three_d cut_dimensions)
{
  return dat->primary &&
    3 == dat->type
      #ifdef BGSUB
        %10
      #endif
    && dat->timeleft > maxrealtime && dat->timeback > -nnegbins
    && dat->remid > remid_cut
    && dat->trklen > trklen_cut
    && dat->contained

    // This cut can more-or-less restrict us to the worst
    // resolution bin in the 2017 numu analysis.
    && (dat->trklen * 4.5/2000)/dat->slce < max_frac_e_mu

    // Sufficient to catch all neutrons within 6 cell widths. Maybe not
    // conservative enough, since neutrons that spill out into the air
    // probably don't ever come back? Or do they?  If they do, that might
    // explain the long time constant which gets fit as "Boron-12".
    && fabs(dat->trkx) < trkx_cut
    && ((!muoncatcher && fabs(dat->trky) < trky_cut) ||
        ( muoncatcher && dat->trky < mucatch_trky_cut && dat->trky > -trky_cut))
    && ((!muoncatcher && dat->trkz < trkz_cuthi && dat->trkz > trkz_cutlo) ||
        ( muoncatcher && dat->trkz > mucatch_trkz_cutlo
       && dat->trkz < mucatch_trkz_cuthi))

    // Try to patch up a deficiency in my ntuples: because calibration
    // is/was broken for the muon catcher in my input files, I don't
    // have the CAF variables that limit the amount of hadronic energy
    // in the muon catcher. So instead make sure that the track starts
    // solidly in the main detector, which should do much of the same
    // job. I put this cut in because I saw a large excess of neutrons
    // with events reconstructed over 3 GeV with the primary track
    // ending in the muon catcher. The data was much higher than the MC
    // no matter how the NC fraction was floated. Let's see if avoiding
    // muon catcher mismodeling/miscalibration of hadronic energy
    // reduces that discrepancy.
    //
    // Update: Pretty sure that's pileup.
    && (!muoncatcher || dat->trkstartz < mucatch_trkstartz_cut)

    && pass_intensity(dat, minslc, maxslc, rhc, cut_dimensions)
    ;
}

// True if the track passes
static bool track_followers_cut(const vector<data> & dats)
{
  bool has_pion_gamma = false, has_xray = false;

  for(unsigned int i = 0; i < dats.size(); i++){
    const data * const dat = &dats[i];

    // reject the event if there is something that looks like a Michel
    // closely following the track.  Don't look too close in time, or
    // we may cut lots of events for late track light.  Don't look after
    // holex_hi, which is where the time series fit starts.
    //
    // I checked that this supresses the background more than the signal,
    // although the effect is not magical.
    if(dat->t > 0.75 && dat->t < holex_hi &&
       dat->mindist <= 3 &&
       dat->e > 10. && dat->e < 70.) return false;

    // In principle you can select pions here by looking for a 100MeV
    // shower coincident with the track, either from charge exchange (a
    // few percent) or radiative capture (a few percent). This does seem
    // to work: I get consistently higher fractions of neutrons. Can
    // this be made interesting?
    /*
    if(dat->t < 0.25 && dat->t > -0.25 &&
       dat->mindist > 3 &&
       dat->e > 50 &&
       dat->remid > 0.5 &&
       dat->type == 3 &&
       dat->trkz < 1250 &&
       ((dat->cosx > -1 && fabs(dat->cosx) < 0.6) ||
        (dat->cosy > -1 && fabs(dat->cosy) < 0.6))
      ){
       has_pion_gamma = true;
       printf("Pion-gamma: %d %d %d %6.1f %6.1f %4.2f %5.2f %5.2f %6.3f\n",
         dat->run, dat->subrun, dat->event,
         dat->trkz, dat->e, dat->mindist, dat->cosx, dat->cosy, dat->remid);
    }

    // Or try to preferentially select muon captures on Ti and Sn using
    // their x-rays.
    if(dat->t < 0.75 && dat->t > -0.75 &&
       dat->mindist > 2 &&
       dat->pe > 35 && dat->e > 1 && dat->e/0.7 < 8){
       has_xray = true;
       printf("Ti-xray: %d %d %d %6.1f %6.1f %4.2f %5.2f %5.2f %6.3f\n",
         dat->run, dat->subrun, dat->event,
         dat->trkz, dat->e, dat->mindist, dat->cosx, dat->cosy, dat->remid);
    }
    */
  }

  return true;
}

static bool clustercut(data * dat, const int mindist,
                       const two_or_three_d cut_dimensions)
{
  return !(dat->t >= holex_lo && dat->t < holex_hi) &&
    dat->t > -nnegbins && dat->t < maxrealtime &&
    (cut_dimensions == TWOD || (dat->nhitx >= 1 && dat->nhity >= 1)) &&
    dat->mindist <= mindist &&
    dat->pe > 35 && dat->e < 20;
}

static void
fill_2dhist(TH1D ** trackcounts, TH2D ** h, data * dat, TTree * t,
            const int mindist, const float minslc, const float maxslc,
            const bool rhc, const two_or_three_d cut_dimensions)
{
  int lastrun = 0, lastevent = 0, lasttrk = 0;
  vector<data> dats;
  const int progint = t->GetEntries()/10;
  int progtarg = progint;
  for(int i = 0; i < t->GetEntries(); i++){
    if(i > progtarg){
      printf("%.0f%% ", (100.*i)/t->GetEntries());
      fflush(stdout);
      progtarg += progint;
    }
    t->GetEntry(i);

    // read in all of the clusters for this track
    if(i != t->GetEntries()-1 &&
       lastrun == dat->run && lastevent == dat->event && lasttrk == dat->trk){
      dats.push_back(*dat);
    }

    // analyze them.  If we like the track, write out its clusters.
    else{
      if(track_followers_cut(dats)){
        for(unsigned int j = 0; j < dats.size(); j++){
          if(track_itself_cut(&dats[j], minslc, maxslc, rhc, cut_dimensions) &&
             clustercut(&dats[j], mindist, cut_dimensions))
            h[dats[j].type > 10]->Fill(dats[j].t, dats[j].slce);
        }
      }

      dats.clear();
      dats.push_back(*dat);
    }
    lastrun = dat->run;
    lastevent = dat->event;
    lasttrk = dat->trk;
  }
}

// In this scheme, all tracks go into the denominator
static void fill_1dhist(TH1D ** h, data * dat, TTree * t,
                        const float minslc, const float maxslc,
                        const bool rhc,
                        const two_or_three_d cut_dimensions)
{
  const int progint = t->GetEntries()/10;
  int progtarg = progint;
  TBranch * ibranch = t->GetBranch("i");
  for(int i = 0; i < t->GetEntries(); i++){
    if(i > progtarg){
      printf("%.0f%% ", (100.*i)/t->GetEntries());
      fflush(stdout);
      progtarg += progint;
    }
    ibranch->GetEntry(i);
    if(dat->i != 0) continue;

    t->GetEntry(i);
    // Don't test for dat->i==0 with type>10.  See tag
    // stagehotelcurlknotgoat in ntuple generation.
    if(track_itself_cut(dat, minslc, maxslc, rhc, cut_dimensions)){
      h[0]->Fill(dat->slce); // The tracks are the same for signal
      #ifdef BGSUB
        h[1]->Fill(dat->slce); // and off-space background
      #endif
    }
  }
}

int rhc_stage_zero(const int mindist, const float minslc,
                   const float maxslc, const string region,
                   const /* two_or_three_d */ int cut_dimensions_)
{
  if(mindist < 0) return 0; // used to compile only

  const two_or_three_d cut_dimensions = (two_or_three_d) cut_dimensions_;

  muoncatcher = region == "muoncatcher";

  data dat;

  std::ofstream o(
    Form("savedhists_mindist%d_nslc%.1f_%.1f_%s_%s.C",
         mindist, minslc, maxslc, region.c_str(),
         two_or_three_d_names[cut_dimensions]));
  o << "{\n";
  for(int i = 0; i < nperiod; i++){
    TFile * inputTFiles = NULL;
    TTree * trees = NULL;

    inputTFiles = new TFile(
      Form("%s", inputfiles[i]), "read");
    if(!inputTFiles || inputTFiles->IsZombie()){
      fprintf(stderr, "Couldn't read a file.  See above.\n");
      return 1;
    }
    trees = dynamic_cast<TTree *>(inputTFiles->Get("t"));
    if(!trees){
      fprintf(stderr, "Couldn't read a tree.  See above.\n");
      return 1;
    }

    TH1D * all_tcounts[SIG_AND_BG] = {
      new TH1D(Form("%s_tcounts", Speriodnames[i]), "",
                    nbins_e[cut_dimensions], bins_e[cut_dimensions]),
      #ifdef BGSUB
        new TH1D(Form("%sBG_tcounts",  Speriodnames[i]), "",
                      nbins_e[cut_dimensions], bins_e[cut_dimensions])
      #endif
    };
    TH2D * fithist[SIG_AND_BG] = {
      new TH2D(Form("%s", Speriodnames[i]), "",
        nnegbins + maxrealtime + additional,
        -nnegbins, maxrealtime + additional,
        nbins_e[cut_dimensions], bins_e[cut_dimensions]),
      #ifdef BGSUB
        new TH2D(Form("%sBG",  Speriodnames[i]), "",
          nnegbins + maxrealtime + additional,
          -nnegbins, maxrealtime + additional,
          nbins_e[cut_dimensions], bins_e[cut_dimensions])
      #endif
    };

    const bool rhc = i < nperiodrhc;

    setbranchaddresses(&dat, trees, FOR2D, cut_dimensions);
    fill_2dhist(all_tcounts, fithist, &dat, trees, mindist, minslc, maxslc,
                rhc, cut_dimensions);
    printf("\nGot %s 2D\n", Speriodnames[i]);

    setbranchaddresses(&dat, trees, FOR1D, cut_dimensions);
    fill_1dhist(all_tcounts, &dat, trees, minslc, maxslc, rhc, cut_dimensions);
    printf("Got %s 1D\n", Speriodnames[i]);
    fflush(stdout);

    for(int super = 0; super < SIG_AND_BG; super++){
      fithist[super]->SavePrimitive(o);
      all_tcounts[super]->SavePrimitive(o);
    }
    inputTFiles->Close();
  }
  o << "}\n";
  o.close();
  return 0;
}
