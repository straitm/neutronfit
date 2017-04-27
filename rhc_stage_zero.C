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

/*
 * TODO: exclude neutron and B-12 candidates following a Michel.
 */

void rhc_stage_zero(const int mindist)
{
  const char * const basecut = Form(
    "run != 11601" // noise at t = 106 in this run.
                   // Will go away with a better run list.
    "&& primary && type == 3 && timeleft > %f && timeback > %f "
    "&& remid > 0.75 "
    "&& trklen > 200 " // standard
    //"&& trklen > 600 " // longer -- drops nearly all NC background
    // standard numu ND containment cuts, except for some muon-catcher
    // track checks which are irrelevant because I'm about to cut
    // those in the next line anyhow. I need the track ends to be more
    // contained than the standard cuts.
    "&& contained"

    // Sufficient to catch all neutrons within 6 cell widths. Maybe not
    // conservative enough, since neutrons that spill out into the air
    // probably don't ever come back? Or do they?
    "&& abs(trkx) < 170 && abs(trky) < 170 && trkz < 1250"
    //"&& nslc <= 6" // reduce pileup
    , maxrealtime, -nnegbins);


  const char * const clustercut = 
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

    Form("nhit >= 1 && mindist <= %d"
    "&& pe > 35 && e < 20", mindist);

  // a very loose cut
  // "1";

  const std::string tcut = Form("i == 0 && %s", basecut);

  const std::string cut = Form(
   "%s && %s && t > %f && t < %f && !(t >= -1 && t < 2)",
    basecut, clustercut, -nnegbins, maxrealtime);

  TFile * inputTFiles[nperiod] = { NULL };
  TTree * trees[nperiod] = { NULL };

  const char * const inputfiles[nperiod] = {
  "/nova/ana/users/mstrait/ndcosmic/prod_pid_"
    "S16-12-07_nd_period6_keepup/1278-type3.root",
  "/nova/ana/users/mstrait/ndcosmic/prod_pid_"
    "R16-12-20-prod3recopreview.b_nd_numi_rhc_epoch4a_v1_goodruns/all-type3.root",

  "/nova/ana/users/mstrait/ndcosmic/prod_pid_"
    "R17-03-01-prod3reco.b_nd_numi_fhc_period1_v1_goodruns/all-type3.root",
  "/nova/ana/users/mstrait/ndcosmic/prod_pid_"
    "R17-03-01-prod3reco.b_nd_numi_fhc_period2_v1_goodruns/all-type3.root",
  "/nova/ana/users/mstrait/ndcosmic/prod_pid_"
    "R17-03-01-prod3reco.b_nd_numi_fhc_period3_v1_goodruns/all-type3.root",
  "/nova/ana/users/mstrait/ndcosmic/prod_pid_"
    "R17-03-01-prod3reco.b_nd_numi_fhc_period5_v1_goodruns/all-type3.root"
  };


  for(int i = 0; i < nperiod; i++){
    inputTFiles[i] = new TFile(inputfiles[i], "read");
    if(!inputTFiles[i] || inputTFiles[i]->IsZombie()){
      fprintf(stderr, "Couldn't read a file.  See above.\n");
      return;
    }
    trees[i] = dynamic_cast<TTree *>(inputTFiles[i]->Get("t"));
    if(!trees[i]){
      fprintf(stderr, "Couldn't read a tree.  See above.\n");
      return;
    }
  }

  std::ofstream o(Form("savedhists_mindist%d.C", mindist));
  o << "{\n";
  for(int i = 0; i < nperiod; i++){
    all_tcounts[i] = new
      TH1D(Form("tcounts_%s", Speriodnames[i]), "", nbins_e, bins_e);
    fithist[i] = new TH2D(Form("%s_s", Speriodnames[i]), "",
      nnegbins + maxrealtime + additional,
      -nnegbins, maxrealtime + additional, nbins_e, bins_e);
    trees[i]->Draw(Form("slce:t >> %s_s"    , Speriodnames[i]),  cut.c_str());
    trees[i]->Draw(Form("slce >> tcounts_%s", Speriodnames[i]), tcut.c_str());
    printf("Got %s_s\n", Speriodnames[i]);
    fflush(stdout);
    fithist[i]->SavePrimitive(o);
    all_tcounts[i]->SavePrimitive(o);
    //inputTFiles[i]->Close(); // XXX why does this seg fault?
  }
  o << "}\n";
  o.close();
}
