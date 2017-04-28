#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TROOT.h"
#include "TGraphAsymmErrors.h"
#include <fstream>

#include "TMinuit.h"
TMinuit * mn = NULL; // dumb, because of common.C

#include "common.C"

struct data{
  int run;
  int event;
  int trk;
  int primary;
  int type;
  int contained;
  int nslc;
  int nhitx;
  int nhity;
  int nhit;
  int i;
  float timeleft;
  float timeback;
  float remid;
  float trklen;
  float trkx;
  float trky;
  float trkz;
  float mindist;
  float pe;
  float e;
  float t;
  float slce;
};

void setbranchaddress(const char * const name, float * d, TTree * t)
{
  t->SetBranchStatus(name, 1);
  t->SetBranchAddress(name, d);
}

void setbranchaddress(const char * const name, int * d, TTree * t)
{
  t->SetBranchStatus(name, 1);
  t->SetBranchAddress(name, d);
}

void setbranchaddresses(data * dat, TTree * t)
{
  t->SetBranchStatus("*", 0);
  setbranchaddress("i", &dat->i, t);
  setbranchaddress("run", &dat->run, t);
  setbranchaddress("event", &dat->event, t);
  setbranchaddress("trk", &dat->trk, t);
  setbranchaddress("primary", &dat->primary, t);
  setbranchaddress("type", &dat->type, t);
  setbranchaddress("contained", &dat->contained, t);
  setbranchaddress("nslc", &dat->nslc, t);
  setbranchaddress("nhitx", &dat->nhitx, t);
  setbranchaddress("nhity", &dat->nhity, t);
  setbranchaddress("nhit", &dat->nhit, t);
  setbranchaddress("timeleft", &dat->timeleft, t);
  setbranchaddress("timeback", &dat->timeback, t);
  setbranchaddress("remid", &dat->remid, t);
  setbranchaddress("trklen", &dat->trklen, t);
  setbranchaddress("trkx", &dat->trkx, t);
  setbranchaddress("trky", &dat->trky, t);
  setbranchaddress("trkz", &dat->trkz, t);
  setbranchaddress("mindist", &dat->mindist, t);
  setbranchaddress("pe", &dat->pe, t);
  setbranchaddress("e", &dat->e, t);
  setbranchaddress("t", &dat->t, t);
  setbranchaddress("slce", &dat->slce, t);
}

// true if it passes the cut
bool basecut(data * dat)
{
  return dat->run != 11601 // noise at t = 106 in this run.
                   // Will go away with a better run list.
    && dat->primary && dat->type == 3
    && dat->timeleft > maxrealtime && dat->timeback > -nnegbins 
    && dat->remid > 0.75 
    && dat->trklen > 200  // standard
    //&& dat->trklen > 600  // longer -- drops nearly all NC background
    // standard numu ND containment cuts, except for some muon-catcher
    // track checks which are irrelevant because I'm about to cut
    // those in the next line anyhow. I need the track ends to be more
    // contained than the standard cuts.
    && dat->contained

    // Sufficient to catch all neutrons within 6 cell widths. Maybe not
    // conservative enough, since neutrons that spill out into the air
    // probably don't ever come back? Or do they?
    && fabs(dat->trkx) < 170 && fabs(dat->trky) < 170 && dat->trkz < 1250
    //&& dat->nslc <= 6 // reduce pileup
    ;
}

// True if the track passes
bool trackcut(const vector<data> & dats)
{
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
  }

  return true;
}

// Other possibilities:
// 
// Attempt to agressively reduce neutrons while still getting B-12
/*
  "nhitx >= 1 && nhity >= 1" // maximum range is about 5.7cm
  "&& nhitx <= 2 && nhity <= 2"
  "&& nhit <= 3"
  "&& mindist < 2.8" // allow up to 1 plane and 2 cells off
  "&& pe > 35 && e > 8*0.62 && e < 25*0.62";
*/

// a pretty strict, reasonable cut
/*
   "nhitx >= 1 && nhity >= 1 && mindist <= 6"
   "&& pe > 70 && e < 20";
*/
bool clustercut(data * dat, const int mindist)
{
  return !(dat->t >= -1 && dat->t < 2) &&
    dat->t > -nnegbins && dat->t < maxrealtime &&
    dat->nhit >= 1 && dat->mindist <= mindist
    && dat->pe > 35 && dat->e < 20;
}

void fill_2dhist(TH2D * h, data * dat, TTree * t, const int mindist)
{
  int lastrun = 0, lastevent = 0, lasttrk = 0;
  vector<data> dats;
  for(int i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);

    // read in all of the clusters for this track
    if(i != t->GetEntries()-1 &&
       lastrun == dat->run && lastevent == dat->event && lasttrk == dat->trk){
      dats.push_back(*dat);
    }

    // analyze them.  If we like the track, write out its clusters.
    else{
      if(trackcut(dats))
        for(unsigned int j = 0; j < dats.size(); j++)
          if(basecut(&dats[j]) && clustercut(&dats[j], mindist))
            h->Fill(dats[j].t, dats[j].slce);

      dats.clear();
      dats.push_back(*dat);
    }
    lastrun = dat->run;
    lastevent = dat->event;
    lasttrk = dat->trk;
  }
}

void fill_1dhist(TH1D * h, data * dat, TTree * t)
{
  for(int i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);
    if(basecut(dat) && dat->i == 0) h->Fill(dat->slce);
  }
}

void rhc_stage_zero(const int mindist)
{
  const char * const inputfiles[nperiod] = {
  "prod_pid_S16-12-07_nd_period6_keepup/1278-type3.root",
  "prod_pid_R16-12-20-prod3recopreview.b_nd_numi_rhc_epoch4a_v1_goodruns/all-type3.root",

  "prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period1_v1_goodruns/all-type3.root",
  "prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period2_v1_goodruns/all-type3.root",
  "prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period3_v1_goodruns/all-type3.root",
  "prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period5_v1_goodruns/all-type3.root"
  };

  data dat;

  std::ofstream o(Form("savedhists_mindist%d.C", mindist));
  o << "{\n";
  for(int i = 0; i < nperiod; i++){
    TFile * inputTFiles = NULL;
    TTree * trees = NULL;

    inputTFiles = new TFile(
      Form("/nova/ana/users/mstrait/ndcosmic/%s", inputfiles[i]), "read");
    if(!inputTFiles || inputTFiles->IsZombie()){
      fprintf(stderr, "Couldn't read a file.  See above.\n");
      return;
    }
    trees = dynamic_cast<TTree *>(inputTFiles->Get("t"));
    if(!trees){
      fprintf(stderr, "Couldn't read a tree.  See above.\n");
      return;
    }
    setbranchaddresses(&dat, trees);

    TH1D * all_tcounts = new
      TH1D(Form("tcounts_%s", Speriodnames[i]), "", nbins_e, bins_e);
    TH2D * fithist = new TH2D(Form("%s_s", Speriodnames[i]), "",
      nnegbins + maxrealtime + additional,
      -nnegbins, maxrealtime + additional, nbins_e, bins_e);

    fill_2dhist(fithist, &dat, trees, mindist);
    fill_1dhist(all_tcounts, &dat, trees);

    printf("Got %s_s\n", Speriodnames[i]);
    fflush(stdout);
    fithist->SavePrimitive(o);
    all_tcounts->SavePrimitive(o);
    inputTFiles->Close();
  }
  o << "}\n";
  o.close();
}
