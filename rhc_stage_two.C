#include "TLine.h"
#include "TMarker.h"
#include "TMinuit.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TError.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TMath.h"
#include <fstream>
#include <iostream>

static TMinuit * mn = NULL;

#include "common.C"

// filled by Macro call at start of running
double hessian[nbins_e*nbeam*2][nbins_e*nbeam*2];

enum conttype { oned68, twod68, twod90 };

const conttype contour_type = twod90;

static const double tsize = 0.04;

const int npar = 4;

const double mum_nyield_nominal = 0.181;
const double mum_nyield_error   = 0.030;

const double neff_nominal = 0.5;
const double neff_error = 0.2;

const double b12eff_nominal = 0.5;
const double b12eff_error = 0.2;

// Ratio of the neutron yield from pions to that of muons, per stop.
const double npimu_nominal = 19.1;
double npimu_error = 3.7;

const double nm_nominal = 1;
const double nm_error = 0.1;

bool useb12 = true; // changed between fits

TH1D *fhc_reco_numubar = new TH1D("fhc_reco_numubar","",100,0,10);
TH1D *fhc_reco_numu = new TH1D("fhc_reco_numu","",100,0,10);
TH1D *rhc_reco_numu = new TH1D("rhc_reco_numu","",100,0,10);
TH1D *rhc_reco_numubar = new TH1D("rhc_reco_numubar","",100,0,10);

TH1D *fhc_reco_nc = new TH1D("fhc_reco_nc","",100,0,10);
TH1D *rhc_reco_nc = new TH1D("rhc_reco_nc","",100,0,10);

TH1D * fhc_tracks = NULL;
TH1D * rhc_tracks = NULL;

TH1D * doublerat_ncomplications = NULL;
TH1D * doublerat_b12complications = NULL;

TGraphAsymmErrors * n_result = new TGraphAsymmErrors;
TGraphAsymmErrors * b12_result = new TGraphAsymmErrors;

TGraphAsymmErrors * g_n_rhc = new TGraphAsymmErrors;
TGraphAsymmErrors * g_n_fhc = new TGraphAsymmErrors;

TGraphAsymmErrors * g_b12_rhc = new TGraphAsymmErrors;
TGraphAsymmErrors * g_b12_fhc = new TGraphAsymmErrors;

// Out here so we can draw them
TH1D * rhc_neutrons = (TH1D*)rhc_reco_numu->Clone("rhc_neutrons");
TH1D * fhc_neutrons = (TH1D*)fhc_reco_numu->Clone("fhc_neutrons");
TH1D * rhc_b12 = (TH1D*)rhc_reco_numu->Clone("rhc_b12");
TH1D * fhc_b12 = (TH1D*)fhc_reco_numu->Clone("fhc_b12");

// Exist just for drawing
TH1D * rhc_neutrons_nc = (TH1D*)rhc_reco_numu->Clone("rhc_neutrons_nc");
TH1D * rhc_neutrons_numu = (TH1D*)rhc_reco_numu->Clone("rhc_neutrons_numu");
TH1D * rhc_neutrons_numubar = (TH1D*)rhc_reco_numu->Clone("rhc_neutrons_numubar");
TH1D * fhc_neutrons_nc = (TH1D*)fhc_reco_numu->Clone("fhc_neutrons_nc");
TH1D * fhc_neutrons_numu = (TH1D*)fhc_reco_numu->Clone("fhc_neutrons_numu");
TH1D * fhc_neutrons_numubar = (TH1D*)fhc_reco_numu->Clone("fhc_neutrons_numubar");
TH1D * rhc_b12_nc = (TH1D*)rhc_reco_numu->Clone("rhc_b12_nc");
TH1D * rhc_b12_numu = (TH1D*)rhc_reco_numu->Clone("rhc_b12_numu");
TH1D * rhc_b12_numubar = (TH1D*)rhc_reco_numu->Clone("rhc_b12_numubar");
TH1D * fhc_b12_nc = (TH1D*)fhc_reco_numu->Clone("fhc_b12_nc");
TH1D * fhc_b12_numu = (TH1D*)fhc_reco_numu->Clone("fhc_b12_numu");
TH1D * fhc_b12_numubar = (TH1D*)fhc_reco_numu->Clone("fhc_b12_numubar");

static void reset_hists()
{
  rhc_neutrons->Reset();
  fhc_neutrons->Reset();
  rhc_b12->Reset();
  fhc_b12->Reset();
  rhc_neutrons_nc->Reset();
  fhc_neutrons_nc->Reset();
  rhc_b12_nc->Reset();
  fhc_b12_nc->Reset();
  rhc_neutrons_numu->Reset();
  fhc_neutrons_numu->Reset();
  rhc_b12_numu->Reset();
  fhc_b12_numu->Reset();
  rhc_neutrons_numubar->Reset();
  fhc_neutrons_numubar->Reset();
  rhc_b12_numubar->Reset();
  fhc_b12_numubar->Reset();
}

static void update_hists(const double mum_nyield, const double b12eff,
                         const double neff, const double npimu,
                         const double nmscale, const double ncscale)
{
  // Probability of getting a neutron from a mu+ via electrodisintigration, roughly
  const double mup_nyield = 1e-4;

  // this could be another nuisance parameter.  Don't want to double count
  // it with the ratio of pion to muon yields, though, so careful.
  const double muminus_capture_frac = 0.182;

  /* We're going to ignore the corner case of pi+ -> mu+ Michel decays
     making a neutron, since we've fudged lots of other bigger stuff. */

  reset_hists();

  rhc_neutrons_nc     ->Add(rhc_reco_nc, ncscale*npimu*mum_nyield*neff);
  rhc_neutrons_numu   ->Add(rhc_reco_numu,    mum_nyield * nmscale /* note */ *neff);
  rhc_neutrons_numubar->Add(rhc_reco_numubar, mup_nyield*neff);

  rhc_neutrons->Add(rhc_neutrons_nc);
  rhc_neutrons->Add(rhc_neutrons_numu);
  rhc_neutrons->Add(rhc_neutrons_numubar);

  rhc_neutrons->Divide(rhc_tracks);
  rhc_neutrons_nc->Divide(rhc_tracks);
  rhc_neutrons_numu->Divide(rhc_tracks);
  rhc_neutrons_numubar->Divide(rhc_tracks);

  // Estimate of number of neutrons from FHC per muon
  fhc_neutrons_nc->Add(fhc_reco_nc, ncscale*npimu*mum_nyield*neff);
  fhc_neutrons_numu->Add(fhc_reco_numu, mum_nyield*neff);
  fhc_neutrons_numubar->Add(fhc_reco_numubar, mup_nyield*neff);

  fhc_neutrons->Add(fhc_neutrons_nc);
  fhc_neutrons->Add(fhc_neutrons_numu);
  fhc_neutrons->Add(fhc_neutrons_numubar);

  fhc_neutrons->Divide(fhc_tracks);
  fhc_neutrons_nc->Divide(fhc_tracks);
  fhc_neutrons_numu->Divide(fhc_tracks);
  fhc_neutrons_numubar->Divide(fhc_tracks);

  // My toy mc atomic cap frac, nucl cap frac on C-12, and Double Chooz!
  // If it ever mattered, the error on this should be taken into account
  const double mum_b12yield = 0.82 * 0.077 * 0.177;

  // I think really zero, although of course the mu+ in flight (or the
  // e+ Michel) could make N-12, which looks the same (in fact, doubled
  // due to half the lifetime). I'd guess this is 1e-5 at most, and
  // maybe much lower.
  const double mup_b12yield = 1e-6;

  // H. Hilscher, W.-D. Krebs, G. Sepp, and V. Soergel. An experimental
  // test of the analogy between radiative pion absorption and muon
  // capture in 12C. Nuclear Physics A, 158(2):584-592, 1970
  // http://www.sciencedirect.com/science/article/pii/037594747090206X
  //
  // pi- have a lower B-12 yield despite always capturing!
  //
  // According to the same paper, pi- make a lot of Li-8. If you're not
  // careful, their figure 2 makes it look that only 20% of the activity
  // at t=0 is B-12, which would make the Li-8 yield almost 100%. But
  // the cycle time they used was much shorter than the Li-8 lifetime,
  // so this is activity produced over many cycles. My reanalysis, if I
  // assume zero uncorrelated background, says that the Li-8 yield is
  // about 4 times the B-12 yield, and with its half-life being 41 times
  // longer, it only adds about 9% to the activity near t=0. Since the
  // error on the B-12 yield is 22%, this can be neglected.
  const double pim_b12yield = mum_b12yield * 5.8e-3 / 1.36e-2;

  // Estimate of number of B-12 from FHC per track
  rhc_b12_nc->Add(rhc_reco_nc, ncscale*pim_b12yield*b12eff);
  rhc_b12_numu->Add(rhc_reco_numu,    mum_b12yield*b12eff);
  rhc_b12_numubar->Add(rhc_reco_numubar, mup_b12yield*b12eff);

  rhc_b12->Add(rhc_b12_nc);
  rhc_b12->Add(rhc_b12_numu);
  rhc_b12->Add(rhc_b12_numubar);

  rhc_b12->Divide(rhc_tracks);
  rhc_b12_nc->Divide(rhc_tracks);
  rhc_b12_numu->Divide(rhc_tracks);
  rhc_b12_numubar->Divide(rhc_tracks);

  // Estimate of number of b12 from FHC per track
  fhc_b12_nc->Add(fhc_reco_nc, ncscale*pim_b12yield*b12eff);
  fhc_b12_numu->Add(fhc_reco_numu,    mum_b12yield*b12eff);
  fhc_b12_numubar->Add(fhc_reco_numubar, mup_b12yield*b12eff);

  fhc_b12->Add(fhc_b12_nc);
  fhc_b12->Add(fhc_b12_numu);
  fhc_b12->Add(fhc_b12_numubar);

  fhc_b12->Divide(fhc_tracks);
  fhc_b12_nc->Divide(fhc_tracks);
  fhc_b12_numu->Divide(fhc_tracks);
  fhc_b12_numubar->Divide(fhc_tracks);

  if(!doublerat_ncomplications)
    doublerat_ncomplications = (TH1D*)rhc_neutrons->Clone("doublerat_ncomplications");

  if(!doublerat_b12complications)
    doublerat_b12complications = (TH1D*)rhc_b12->Clone("doublerat_b12complications");

  doublerat_ncomplications->Divide(rhc_neutrons, fhc_neutrons);
  doublerat_b12complications->Divide(rhc_b12, fhc_b12);
}

static double compare(double * dat, double * datup, double * datdn,
                      double * model)
{
  double chi2 = 0;
  for(int i = 0; i < nbins_e*nbeam*2 /* n and b12 */; i++){
    for(int j = 0; j < nbins_e*nbeam*2; j++){
      if(!useb12 && (i >= nbins_e*nbeam || j >= nbins_e*nbeam)) continue;

      const double r = dat[i];
      const double pred = model[i];
      const double rj = dat[j];
      const double predj = model[j];

      // for diagonal elements, try to incorporate MINOS information
      const double rup = datup[i];
      const double rdn = datdn[i];
      const double correction = i != j || rdn == 0 || rup == 0?  1
      : (r > pred?rup:rdn)/(rup + rdn)*2;

      chi2 += correction*hessian[i][j]*(pred - r)*(predj - rj);
    }
  }

  return chi2;
}

static void fcn(__attribute__((unused)) int & np,
  __attribute__((unused)) double * gin, double & chi2, double *par,
  __attribute__((unused)) int flag)  
{
  chi2 = 0;

  const double ncscale = par[0];
  const double nmscale = par[1];
  const double npimu   = par[2];
  const double neff    = par[3];
  const double b12eff  = par[4];
  const double mum_nyield = par[5];
  update_hists(mum_nyield, b12eff, neff, npimu, nmscale, ncscale);

  // penalty terms
  
  // Justified by external data
  chi2 += pow((npimu - npimu_nominal)/npimu_error, 2);

  // Justified by external data
  chi2 += pow((mum_nyield - mum_nyield_nominal)/mum_nyield_error, 2);

  // not really justified
  //chi2 += pow((neff - neff_nominal)/neff_error, 2);
  //chi2 += pow((b12eff - b12eff_nominal)/b12eff_error, 2);

  // Can use this if we really think we have an external handle
  // on the numu contamination in RHC, but better to leave it out
  // and see that the fit finds a reasonable value
  //chi2 += pow((nmscale - nm_nominal)/nm_error, 2);

  static double alldat[nbins_e*4],
                alldatup[nbins_e*4],
                alldatdn[nbins_e*4];
  static bool first = true;
  if(first){
    first = false;
    TGraphAsymmErrors * ingraphs[4] = { g_n_rhc, g_n_fhc, g_b12_rhc, g_b12_fhc };

    int i = 0;
    for(int g = 0; g < 4; g++){
      for(int j = 0; j < ingraphs[g]->GetN(); j++){
        alldat[i] = ingraphs[g]->GetY()[j];
        alldatup[i] = ingraphs[g]->GetErrorYhigh(j);
        alldatdn[i] = ingraphs[g]->GetErrorYlow (j);
        i++;
      }
    }
  }

  double allmodel[nbins_e*4];
  TH1D * inhists[4] = { rhc_neutrons, fhc_neutrons, rhc_b12, fhc_b12 };
  int i = 0;
  for(int g = 0; g < 4; g++)
    for(int b = 1; b <= inhists[g]->GetNbinsX(); b++)
       allmodel[i++] = inhists[g]->GetBinContent(b);

  chi2 += compare(alldat, alldatup, alldatdn, allmodel);
}

void fill_hists(const char * const file, TH1D * const numu,
                TH1D * const numubar, TH1D * const nc)
{
  TFile * f = new TFile(file, "read");

  if(!f || f->IsZombie()){
    fprintf(stderr, "Could not open the file %s\n", file);
    _exit(1);
  }

  TTree * t = dynamic_cast<TTree *>(f->Get("t"));

  if(!t){
    fprintf(stderr, "Could not get the tree from %s\n", file);
    _exit(1);
  }

  // For I-don't-know-why TTree::Draw isn't working for me here, so
  // do it the hard way
  int i, primary, true_pdg, true_nupdg, true_nucc, contained;
  float slce, trkx, trky, trkz, trklen, remid;
  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("slce", 1); t->SetBranchAddress("slce", &slce);
  t->SetBranchStatus("i", 1); t->SetBranchAddress("i", &i);
  t->SetBranchStatus("primary", 1); t->SetBranchAddress("primary", &primary);
  t->SetBranchStatus("true_pdg", 1); t->SetBranchAddress("true_pdg", &true_pdg);
  t->SetBranchStatus("true_nupdg", 1); t->SetBranchAddress("true_nupdg", &true_nupdg);
  t->SetBranchStatus("true_nucc", 1); t->SetBranchAddress("true_nucc", &true_nucc);
  t->SetBranchStatus("contained", 1); t->SetBranchAddress("contained", &contained);
  t->SetBranchStatus("trkx", 1); t->SetBranchAddress("trkx", &trkx);
  t->SetBranchStatus("trky", 1); t->SetBranchAddress("trky", &trky);
  t->SetBranchStatus("trkz", 1); t->SetBranchAddress("trkz", &trkz);
  t->SetBranchStatus("trklen", 1); t->SetBranchAddress("trklen", &trklen);
  t->SetBranchStatus("remid", 1); t->SetBranchAddress("remid", &remid);

  for(int e = 0; e < t->GetEntries(); e++){
    t->GetEntry(e);
    if(!(i == 0 && primary && contained)) continue;

    // MUST MATCH cut at top of rhc.C
    if(fabs(trkx) > 170) continue;
    if(fabs(trky) > 170) continue;
    if(     trkz  >1250) continue;
    if(trklen < 200) continue;
    if(remid  < 0.75) continue;

    // XXX What should I be doing about the case where there's a 
    // true numu CC event where a pion is selected as the primary
    // track?

    // based entirely on what the track actually is
    if(true_pdg == -321 /* take K- too, but they are very rare */
       || true_pdg == -211) nc->Fill(slce);
    else if(true_pdg == 13) numu->Fill(slce);

    // XXX use numubar as a stand-in for all other tracks.  Can't just
    // throw out tracks because everything is done per-track!
    // This will make the guesses at neutron yield for mu+ a bit wrong,
    // but they are only guesses anyhow
    else /* if(true_pdg == -13) */ numubar->Fill(slce);

    if(true_pdg == -321) printf("K-!\n");

    /*
    // based on the neutrino interaction
    if(true_nucc){
      if     (true_nupdg ==  14) numu->Fill(slce);
      else if(true_nupdg == -14) numubar->Fill(slce);
    }
    else{
      // Picking out true NC events where the primary track is a pi-
      if(true_pdg == -211 || true_pdg == -321) nc->Fill(slce);
    }
    */
  }

  // can I safely close the file here?! You'd think so, but sometimes
  // ROOT really surprises me!
}


void set_hists()
{   
  // I don't need to scale these by POT because everything is done in a ratio
  // between the number of primary tracks in the given beam type and the truth
  // of those tracks.  But it might be confusing if the individual histograms
  // are drawn without normalization.
  fill_hists("prod_pid_R17-03-01-prod3reco.d_nd_genie_nonswap_"
    "fhc_nova_v08_period5_v1/all-type3.root", fhc_reco_numu, fhc_reco_numubar,
    fhc_reco_nc);

  fill_hists("prod_pid_R16-12-20-prod3recopreview.b_nd_genie_nonswap_"
    "rhc_nova_v08_epoch4a_v3/all-type3.root", rhc_reco_numu, rhc_reco_numubar,
    rhc_reco_nc);

  rhc_reco_numubar=(TH1D*)rhc_reco_numubar->Rebin(nbins_e,"rhc_reco_numubar",bins_e);
  fhc_reco_numubar=(TH1D*)fhc_reco_numubar->Rebin(nbins_e,"fhc_reco_numubar",bins_e);

  rhc_reco_numu   =(TH1D*)rhc_reco_numu   ->Rebin(nbins_e,"rhc_reco_numu"   ,bins_e);
  fhc_reco_numu   =(TH1D*)fhc_reco_numu   ->Rebin(nbins_e,"fhc_reco_numu"   ,bins_e);

  fhc_reco_nc    = (TH1D*)fhc_reco_nc   ->Rebin(nbins_e, "fhc_reco_nc"   , bins_e);
  rhc_reco_nc    = (TH1D*)rhc_reco_nc   ->Rebin(nbins_e, "rhc_reco_nc"   , bins_e);

  rhc_neutrons     =
    (TH1D*)rhc_neutrons     ->Rebin(nbins_e,"rhc_neutrons"     , bins_e);
  rhc_neutrons_nc  =
    (TH1D*)rhc_neutrons_nc  ->Rebin(nbins_e,"rhc_neutrons_nc"  , bins_e);
  rhc_neutrons_numu=
    (TH1D*)rhc_neutrons_numu->Rebin(nbins_e,"rhc_neutrons_numu", bins_e);
  rhc_neutrons_numubar = 
    (TH1D*)rhc_neutrons_numubar->Rebin(nbins_e, "rhc_neutrons_numubar", bins_e);
  fhc_neutrons = 
    (TH1D*)fhc_neutrons->Rebin(nbins_e, "fhc_neutrons", bins_e);
  fhc_neutrons_nc = 
    (TH1D*)fhc_neutrons_nc->Rebin(nbins_e, "fhc_neutrons_nc", bins_e);
  fhc_neutrons_numu = 
    (TH1D*)fhc_neutrons_numu->Rebin(nbins_e, "fhc_neutrons_numu", bins_e);
  fhc_neutrons_numubar = 
    (TH1D*)fhc_neutrons_numubar->Rebin(nbins_e, "fhc_neutrons_numubar", bins_e);

  rhc_b12 = (TH1D*)rhc_b12->Rebin(nbins_e, "rhc_b12", bins_e);
  rhc_b12_nc = (TH1D*)rhc_b12_nc->Rebin(nbins_e, "rhc_b12_nc", bins_e);
  rhc_b12_numu = (TH1D*)rhc_b12_numu->Rebin(nbins_e, "rhc_b12_numu", bins_e);
  rhc_b12_numubar = (TH1D*)rhc_b12_numubar->Rebin(nbins_e, "rhc_b12_numubar", bins_e);
  fhc_b12 = (TH1D*)fhc_b12->Rebin(nbins_e, "fhc_b12", bins_e);
  fhc_b12_nc = (TH1D*)fhc_b12_nc->Rebin(nbins_e, "fhc_b12_nc", bins_e);
  fhc_b12_numu = (TH1D*)fhc_b12_numu->Rebin(nbins_e, "fhc_b12_numu", bins_e);
  fhc_b12_numubar = (TH1D*)fhc_b12_numubar->Rebin(nbins_e, "fhc_b12_numubar", bins_e);


  // Total numu+numubar in RHC
  rhc_tracks = (TH1D*)rhc_reco_numu->Clone("rhc_tracks");
  rhc_tracks->Add(rhc_reco_numubar);
  rhc_tracks->Add(rhc_reco_nc);

  // Total numu+numubar in FHC
  fhc_tracks = (TH1D*)fhc_reco_numu->Clone("fhc_tracks");
  fhc_tracks->Add(fhc_reco_numubar);
  fhc_tracks->Add(fhc_reco_nc);
}

void make_mn()
{
  if(mn) delete mn;
  mn = new TMinuit(npar);
  mn->fGraphicsMode = false;
  const int print = 0;

  // not sure which of these works
  mn->SetPrintLevel(print);
  mn->Command(Form("SET PRINT %d", print));

  // Observed to help get MINOS errors out
  mn->Command("SET STRATEGY 2");

  mn->SetFCN(fcn);
  mn->Command(Form("SET ERR %f",
    contour_type == oned68? 1.00:
    contour_type == twod68? 2.30:
                            4.61
  ));

  int mnparmerr = 0;
  mn->mnparm(0, "NCscale", 1, 0.1, 0, 0, mnparmerr);
  mn->mnparm(1, "NMscale", 1, 0.1, 0, 0, mnparmerr);
  mn->mnparm(2, "npimu", npimu_nominal, 0.5, 0, 0, mnparmerr);
  mn->mnparm(3, "neff", neff_nominal, 0.1, 0, 0, mnparmerr);
  mn->mnparm(4, "b12eff", b12eff_nominal, 0.1, 0, 0, mnparmerr);
  mn->mnparm(5, "mum_nyield", mum_nyield_nominal, 0.03, 0, 0, mnparmerr);
}

// I believe that TGraph::Integral does a different thing
static double getscale(TGraphAsymmErrors * g, TH1D * h)
{
  g->Fit("pol0", "q0", "");
  const double gscale = g->GetFunction("pol0")?
                        g->GetFunction("pol0")->GetParameter(0):
                        1;

  for(int i = 1; i <= h->GetNbinsX(); i++)
    h->SetBinError(i, g->GetErrorYhigh(i));
  h->Fit("pol0", "q0", "");
  const double hscale = h->GetFunction("pol0")?
                        h->GetFunction("pol0")->GetParameter(0):
                        1;
  for(int i = 1; i <= h->GetNbinsX(); i++)
    h->SetBinError(i, 0);
 
  return gscale/hscale;
}

static void styleleg(TLegend * leg)
{
  leg->SetTextFont(42);
  leg->SetBorderSize(1);
  leg->SetFillStyle(0);
  leg->SetTextSize(tsize);
}

void draw()
{
  //////////////////////////////////////////////////////////////////////
  TCanvas * c1 = new TCanvas("rhc1", "rhc1");
  c1->Print("fit_stage_two.pdf[");
  c1->Print("fit_stage_two.pdf");
  const double leftmargin = 0.12;
  const double topmargin  = 0.025;
  const double rightmargin= 0.03;
  const double bottommargin=0.14;
  c1->SetMargin(leftmargin, rightmargin, bottommargin, topmargin);
  const bool logy = false;
  c1->SetLogy(logy);
  TH2D * dum = new TH2D("dm", "", 100, 0, 10, 1000, logy?1e-6:0, 2);
  TH2D * dum2 = (TH2D*) dum->Clone("dm2");
  dum->GetXaxis()->SetRangeUser(bins_e[0], bins_e[nbins_e]);
  if(!logy) dum->GetYaxis()->SetRangeUser(0, 
    min(1.99, 1.05*max(gdrawmax(n_result), gdrawmax(b12_result))));
  dum->GetXaxis()->CenterTitle();
  dum->GetYaxis()->CenterTitle();
  dum->GetXaxis()->SetTitle("Reconstructed E_{#nu} (GeV)");
  dum->GetYaxis()->SetTitle("RHC/FHC");
  dum->Draw();

  n_result->SetMarkerStyle(kFullCircle);
  b12_result->SetMarkerStyle(kOpenCircle);

  const int gray = kGray; // or kGray+1

  b12_result->SetMarkerColor(gray);
  b12_result->SetLineColor(gray);

  doublerat_b12complications->SetMarkerColor(gray);
  doublerat_b12complications->SetLineColor(gray);

  doublerat_b12complications->SetLineStyle(kDashed);
  doublerat_ncomplications->SetLineStyle(kDashed);
  doublerat_ncomplications->SetLineWidth(2);

  b12_result->Draw("pz");
  doublerat_b12complications->Draw("same");
  n_result->Draw("pz");
  doublerat_ncomplications->Draw("same");

  TLegend * leg = new TLegend(0.73, 0.73, 1-rightmargin, 1-topmargin);
  leg->AddEntry(n_result, "Neutrons, data", "lpe");
  leg->AddEntry(doublerat_ncomplications, "Neutrons, fit", "l");
  leg->AddEntry(b12_result, "^{12}B, data", "lpe");
  leg->AddEntry(doublerat_b12complications, "^{12}B, fit", "l");
  styleleg(leg);
  leg->Draw();

  c1->Print("fit_stage_two.pdf");
  
  //////////////////////////////////////////////////////////////////////
  stylegraph(g_n_rhc, kRed+3, kSolid, kOpenSquare, 1, 0.0);
  stylegraph(g_n_fhc, kBlue+3, kSolid, kOpenCircle, 1, 0.0);

  const double rhcscale = getscale(g_n_rhc, rhc_neutrons);
  const double fhcscale = getscale(g_n_fhc, fhc_neutrons);

  rhc_neutrons->SetLineColor(kRed);
  rhc_neutrons->SetLineWidth(2);
  rhc_neutrons_nc->SetLineColor(kRed);
  rhc_neutrons_nc->SetLineStyle(kDotted);
  rhc_neutrons_nc->SetLineWidth(2);
  rhc_neutrons_numu->SetLineColor(kRed);
  rhc_neutrons_numu->SetLineStyle(kDashed);
  rhc_neutrons_numu->SetLineWidth(2);
  rhc_neutrons_numubar->SetLineColor(kRed);
  rhc_neutrons_numubar->SetLineStyle(9);
  rhc_neutrons_numubar->SetLineWidth(2);

  fhc_neutrons->SetLineWidth(2);
  fhc_neutrons_nc->SetLineWidth(2);
  fhc_neutrons_numu->SetLineWidth(2);
  fhc_neutrons_numubar->SetLineWidth(2);
  fhc_neutrons->SetLineColor(kBlue);
  fhc_neutrons_nc->SetLineColor(kBlue);
  fhc_neutrons_numu->SetLineColor(kBlue);
  fhc_neutrons_numubar->SetLineColor(kBlue);
  fhc_neutrons_nc->SetLineStyle(kDotted);
  fhc_neutrons_numu->SetLineStyle(kDashed);
  fhc_neutrons_numubar->SetLineStyle(9);

  rhc_neutrons->Scale(rhcscale);
  fhc_neutrons->Scale(fhcscale);
  rhc_neutrons_nc->Scale(rhcscale);
  fhc_neutrons_nc->Scale(fhcscale);
  rhc_neutrons_numu->Scale(rhcscale);
  fhc_neutrons_numu->Scale(fhcscale);
  rhc_neutrons_numubar->Scale(rhcscale);
  fhc_neutrons_numubar->Scale(fhcscale);

  dum2->GetXaxis()->SetRangeUser(bins_e[0], bins_e[nbins_e]);
  if(!logy) dum2->GetYaxis()->SetRangeUser(0, 
    min(2.5, 1.05*max(gdrawmax(g_n_rhc), gdrawmax(g_n_fhc))));
  dum2->GetYaxis()->SetTitle("Neutrons/track");
  dum2->GetXaxis()->SetTitle("Reconstructed E_{#nu} (GeV)");
  dum2->GetYaxis()->CenterTitle();
  dum2->GetXaxis()->CenterTitle();

  //
  TCanvas * c2r = new TCanvas("rhc2r", "rhc2r");
  c2r->SetLogy(logy);
  c2r->SetMargin(leftmargin, rightmargin, bottommargin, topmargin);
  dum2->Draw();

  g_n_rhc->Draw("pz");

  TH1D * c2hists[4] = {
    rhc_neutrons,
    rhc_neutrons_nc,
    rhc_neutrons_numu,
    rhc_neutrons_numubar,
  };
    
  for(int i = 1; i <= rhc_neutrons->GetNbinsX(); i++)
    for(int h = 0; h < 4; h++)
      c2hists[h]->Draw("histsame][");

  leg = new TLegend(leftmargin+0.1, 0.70, leftmargin+0.33, 1-topmargin);
  styleleg(leg);
  leg->SetMargin(0.4);
  leg->AddEntry(g_n_rhc, "RHC data", "lpe");
  leg->AddEntry(rhc_neutrons, "RHC fit", "l");
  leg->AddEntry(rhc_neutrons_nc, "RHC NC", "l");
  leg->AddEntry(rhc_neutrons_numu, "RHC #nu_{#mu}", "l");
  leg->AddEntry(rhc_neutrons_numubar, "RHC #bar{#nu}_{#mu}", "l");
  leg->Draw();

  c2r->Print("fit_stage_two.pdf");

  //
  TCanvas * c2f = new TCanvas("rhc2f", "rhc2f");
  c2f->SetLogy(logy);
  c2f->SetMargin(leftmargin, rightmargin, bottommargin, topmargin);
  dum2->Draw();

  g_n_fhc->Draw("pz");
    
  TH1D * c2histsf[4] = {
    fhc_neutrons,
    fhc_neutrons_nc,
    fhc_neutrons_numu,
    fhc_neutrons_numubar
  };
  for(int i = 1; i <= rhc_neutrons->GetNbinsX(); i++)
    for(int h = 0; h < 4; h++)
      c2histsf[h]->Draw("histsame][");

  TLegend * legf = new TLegend(leftmargin+0.1, 0.70, leftmargin+0.33, 1-topmargin);
  styleleg(legf);
  legf->SetMargin(0.4);
  legf->AddEntry(g_n_fhc, "FHC data", "lpe");
  legf->AddEntry(fhc_neutrons, "FHC fit", "l");
  legf->AddEntry(fhc_neutrons_nc, "FHC NC", "l");
  legf->AddEntry(fhc_neutrons_numu, "FHC #nu_{#mu}", "l");
  legf->AddEntry(fhc_neutrons_numubar, "FHC #bar{#nu}_{#mu}", "l");
  legf->Draw();

  c2f->Print("fit_stage_two.pdf");

  //////////////////////////////////////////////////////////////////////
  stylegraph(g_b12_rhc, kRed+3, kSolid, kOpenSquare, 1, 0.0);
  stylegraph(g_b12_fhc, kBlue+3, kSolid, kOpenCircle, 1, 0.0);

  const double brhcscale = getscale(g_b12_rhc, rhc_b12);
  const double bfhcscale = getscale(g_b12_fhc, fhc_b12);

  rhc_b12->SetLineColor(kRed);
  rhc_b12->SetLineWidth(2);
  rhc_b12_nc->SetLineColor(kRed);
  rhc_b12_nc->SetLineStyle(kDotted);
  rhc_b12_nc->SetLineWidth(2);
  rhc_b12_numu->SetLineColor(kRed);
  rhc_b12_numu->SetLineStyle(kDashed);
  rhc_b12_numu->SetLineWidth(2);
  rhc_b12_numubar->SetLineColor(kRed);
  rhc_b12_numubar->SetLineStyle(9);
  rhc_b12_numubar->SetLineWidth(2);

  fhc_b12->SetLineWidth(2);
  fhc_b12_nc->SetLineWidth(2);
  fhc_b12_numu->SetLineWidth(2);
  fhc_b12_numubar->SetLineWidth(2);
  fhc_b12->SetLineColor(kBlue);
  fhc_b12_nc->SetLineColor(kBlue);
  fhc_b12_numu->SetLineColor(kBlue);
  fhc_b12_numubar->SetLineColor(kBlue);
  fhc_b12_nc->SetLineStyle(kDotted);
  fhc_b12_numu->SetLineStyle(kDashed);
  fhc_b12_numubar->SetLineStyle(9);

  rhc_b12->Scale(brhcscale);
  fhc_b12->Scale(bfhcscale);
  rhc_b12_nc->Scale(brhcscale);
  fhc_b12_nc->Scale(bfhcscale);
  rhc_b12_numu->Scale(brhcscale);
  fhc_b12_numu->Scale(bfhcscale);
  rhc_b12_numubar->Scale(brhcscale);
  fhc_b12_numubar->Scale(bfhcscale);

  dum2->GetXaxis()->SetRangeUser(bins_e[0], bins_e[nbins_e]);
  if(!logy) dum2->GetYaxis()->SetRangeUser(0, 
    min(2.5, 1.03*max(gdrawmax(g_b12_rhc), gdrawmax(g_b12_fhc))));
  dum2->GetYaxis()->SetTitle("^{12}B/track");
  dum2->GetXaxis()->SetTitle("Reconstructed E_{#nu} (GeV)");
  dum2->GetYaxis()->CenterTitle();
  dum2->GetXaxis()->CenterTitle();

  //
  TCanvas * c2rb = new TCanvas("rhc2rb", "rhc2rb");
  c2rb->SetLogy(logy);
  c2rb->SetMargin(leftmargin, rightmargin, bottommargin, topmargin);
  dum2->Draw();

  g_b12_rhc->Draw("pz");

  TH1D * c2histsb[4] = {
    rhc_b12,
    rhc_b12_nc,
    rhc_b12_numu,
    rhc_b12_numubar,
  };
    
  for(int i = 1; i <= rhc_b12->GetNbinsX(); i++)
    for(int h = 0; h < 4; h++)
      c2histsb[h]->Draw("histsame][");

  leg = new TLegend(leftmargin+0.1, 0.70, leftmargin+0.33, 1-topmargin);
  styleleg(leg);
  leg->SetMargin(0.4);
  leg->AddEntry(g_n_rhc, "RHC data", "lpe");
  leg->AddEntry(rhc_b12, "RHC fit", "l");
  leg->AddEntry(rhc_b12_nc, "RHC NC", "l");
  leg->AddEntry(rhc_b12_numu, "RHC #nu_{#mu}", "l");
  leg->AddEntry(rhc_b12_numubar, "RHC #bar{#nu}_{#mu}", "l");
  leg->Draw();

  c2rb->Print("fit_stage_two.pdf");

  //
  TCanvas * c2fb = new TCanvas("rhc2fb", "rhc2fb");
  c2fb->SetLogy(logy);
  c2fb->SetMargin(leftmargin, rightmargin, bottommargin, topmargin);
  dum2->Draw();

  g_b12_fhc->Draw("pz");
    
  TH1D * c2histsfb[4] = {
    fhc_b12,
    fhc_b12_nc,
    fhc_b12_numu,
    fhc_b12_numubar
  };
  for(int i = 1; i <= rhc_b12->GetNbinsX(); i++)
    for(int h = 0; h < 4; h++)
      c2histsfb[h]->Draw("histsame][");

  legf = new TLegend(leftmargin+0.1, 0.70, leftmargin+0.33, 1-topmargin);
  styleleg(legf);
  legf->SetMargin(0.4);
  legf->AddEntry(g_n_fhc, "FHC data", "lpe");
  legf->AddEntry(fhc_b12, "FHC fit", "l");
  legf->AddEntry(fhc_b12_nc, "FHC NC", "l");
  legf->AddEntry(fhc_b12_numu, "FHC #nu_{#mu}", "l");
  legf->AddEntry(fhc_b12_numubar, "FHC #bar{#nu}_{#mu}", "l");
  legf->Draw();

  c2fb->Print("fit_stage_two.pdf");

  //////////////////////////////////////////////////////////////////////
  TCanvas * c3 = new TCanvas("rhc3", "rhc3");
  c3->SetMargin(leftmargin, rightmargin, bottommargin, topmargin);

  const double ncmin = 0, ncmax = 1.2,
    nmmin = 0.8, nmmax = 2.6;
  TH2D * dum3 = new TH2D("dm", "", 100, ncmin, ncmax, 1, nmmin, nmmax);
  dum3->GetXaxis()->SetTitle("NC scale");
  dum3->GetYaxis()->SetTitle("#nu_{#mu} scale");
  dum3->GetXaxis()->CenterTitle();
  dum3->GetYaxis()->CenterTitle();
  dum3->Draw();

  TLine * l1 = new TLine(1, nmmin, 1, nmmax);
  l1->Draw();
  TLine * l2 = new TLine(ncmin, 1, ncmax, 1);
  l2->Draw();

  mn->fGraphicsMode = true;

  mn->Command("MNCONT 1 2 50");
  TGraph * cont1 = dynamic_cast<TGraph *>(mn->GetPlot());

  useb12 = false;
  mn->Command("MIGRAD");
  mn->Command("MNCONT 1 2 50");
  TGraph * cont2 = dynamic_cast<TGraph *>(mn->GetPlot());

  useb12 = true;
  mn->Command("MIGRAD");
  const double proper_npimu_error = npimu_error;
  npimu_error = 1.1;
  mn->Command("MNCONT 1 2 50");
  TGraph * cont3 = dynamic_cast<TGraph *>(mn->GetPlot());

  mn->Command("FIX 3");
  mn->Command("MNCONT 1 2 50");
  TGraph * cont4 = dynamic_cast<TGraph *>(mn->GetPlot());

  // restore
  npimu_error = proper_npimu_error;
  mn->Command("REL 3");
  mn->Command("MIGRAD");

  if(cont4 != NULL){
    cont4->SetFillStyle(1001);
    cont4->SetFillColor(kViolet);
    cont4->SetLineColor(kViolet);
    cont4->Draw("f");
  }
  if(cont3 != NULL){
    cont3->SetFillStyle(1001);
    cont3->SetFillColor(kRed);
    cont3->SetLineColor(kRed);

    // get rid of the visual gap 
    cont3->SetPoint(cont3->GetN(), cont3->GetX()[0], cont3->GetY()[0]);
    cont3->Draw("l");
  }
  if(cont1 != NULL){
    cont1->SetFillStyle(1001);
    cont1->SetFillColor(kBlue-2);
    cont1->SetLineColor(kBlue-2);
    cont1->SetPoint(cont1->GetN(), cont1->GetX()[0], cont1->GetY()[0]);
    cont1->Draw("l");
  }
  if(cont2 != NULL){
    cont2->SetFillColor(kBlue);
    cont2->SetLineColor(kBlue);
    cont2->SetLineStyle(kDashed);
    cont2->Draw("l");
  }
  
  mn->fGraphicsMode = false;

  TMarker * bestfit = new TMarker(getpar(0), getpar(1), kFullCircle);
  bestfit->Draw();

  mn->Command("MINOS 10000 1");
  mn->Command("MINOS 10000 2");

  TGraphAsymmErrors * onederrs = new TGraphAsymmErrors;
  onederrs->SetPoint(0, getpar(0), getpar(1) - getminerrdn(1)*1.1);
  onederrs->SetPoint(1, getpar(0)-getminerrdn(0)*1.1, getpar(1));

  mn->Command("SET ERR 1"); // get one D 68% errors
  mn->Command("MIGRAD");
  mn->Command("MINOS 10000 1");
  mn->Command("MINOS 10000 2");
  mn->Command("MINOS 10000 4");
  mn->Command("MINOS 10000 5");

  onederrs->SetPointError(0, getminerrdn(0), getminerrup(0), 0, 0);
  onederrs->SetPointError(1, 0, 0, getminerrdn(1), getminerrup(1));
  onederrs->SetMarkerStyle(kFullCircle);
  onederrs->Draw("pz");

  mn->Command(Form("SET ERR %f",
    contour_type == oned68? 1.00:
    contour_type == twod68? 2.30:
                            4.61
  )); // put back

  leg = new TLegend(0.42, 0.75, 1-rightmargin, 1-topmargin);
  styleleg(leg);
  leg->SetFillStyle(1001);
  leg->SetMargin(0.1);
  leg->AddEntry(cont1,
    Form("%s, #pi/#mu neutron yield is %.1f #pm %.1f",
    contour_type == oned68?"1D 68%":
    contour_type == twod68?"2D 68%":
                           "2D 90%",
    npimu_nominal, npimu_error),
    "l");
  leg->AddEntry(onederrs, "1D 68% errors", "lp");
  leg->AddEntry(cont2, "without using ^{12}B", "l");
  leg->AddEntry(cont3, "perfectly known nuclear #pi/#mu neutron yield", "l");
  leg->AddEntry(cont4, "and perfectly known atomic capture ratios", "f");
  leg->Draw();
  c3->Print("fit_stage_two.pdf");

  //////////////////////////////////////////////////////////////////////
  TCanvas * c4 = new TCanvas("rhc4", "rhc4");
  c4->SetLogy();
  c4->SetMargin(leftmargin, rightmargin, bottommargin, topmargin);

  fhc_reco_numu->SetLineColor(kBlue);
  fhc_reco_numu->SetMarkerColor(kBlue);
  fhc_reco_numu->SetLineWidth(2);

  rhc_reco_numu->SetLineColor(kRed);
  rhc_reco_numu->SetMarkerColor(kRed);
  rhc_reco_numu->SetLineWidth(2);

  fhc_reco_numubar->SetLineColor(kBlue-9);
  fhc_reco_numubar->SetMarkerColor(kBlue-9);
  fhc_reco_numubar->SetLineWidth(1);

  rhc_reco_numubar->SetLineColor(kRed-9);
  rhc_reco_numubar->SetMarkerColor(kRed-9);
  rhc_reco_numubar->SetLineWidth(1);

  fhc_reco_nc->SetLineColor(kBlue);
  fhc_reco_nc->SetMarkerColor(kBlue);
  fhc_reco_nc->SetLineWidth(4);

  rhc_reco_nc->SetLineColor(kRed);
  rhc_reco_nc->SetMarkerColor(kRed);
  rhc_reco_nc->SetLineWidth(4);
 
  fhc_reco_numu->GetYaxis()->SetRangeUser(max(0.1, rhc_reco_numu->GetMinimum()*0.75),
                                          fhc_reco_numu->GetMaximum()*1.25);
  fhc_reco_numu->GetYaxis()->CenterTitle();
  fhc_reco_numu->GetXaxis()->CenterTitle();
  fhc_reco_numu->GetYaxis()->SetTitle("Selected tracks/bin");
  fhc_reco_numu->GetXaxis()->SetTitle("Reconstructed E_{#nu} (GeV)");

  // neutron producers
  fhc_reco_numu->Draw("ehist");
  rhc_reco_numu->Draw("ehistsame");
  fhc_reco_nc->Draw("ehistsame");
  rhc_reco_nc->Draw("ehistsame");

  // non-neutron producers
  fhc_reco_numubar->Draw("ehistsame");
  rhc_reco_numubar->Draw("ehistsame");
  
  

  TLegend * leg4 = new TLegend(0.73, 0.63, 1-rightmargin, 1-topmargin);
  leg4->AddEntry((TH1D*)NULL, "Unscaled", "");
  leg4->AddEntry(fhc_reco_numu, "#mu^{-} in FHC", "l");
  leg4->AddEntry(rhc_reco_numu, "#mu^{-} in RHC", "l");
  leg4->AddEntry(fhc_reco_numubar, "#mu^{+} in FHC", "l");
  leg4->AddEntry(rhc_reco_numubar, "#mu^{+} in RHC", "l");
  leg4->AddEntry(fhc_reco_nc, "#pi^{-} in FHC", "l");
  leg4->AddEntry(rhc_reco_nc, "#pi^{-} in RHC", "l");
  styleleg(leg4);
  leg4->Draw();

  c4->Print("fit_stage_two.pdf");


  c1->Print("fit_stage_two.pdf]"); // doesn't print anything, just closes
}


void rhc_stage_two(const char * const input)
{
  TH1::SetDefaultSumw2();

  set_hists();

  gROOT->Macro(input);

  make_mn();

  mn->Command("SET LIM 1 0.01 10"); // NC scale
  mn->Command("SET LIM 2 0.01 10");  // numubar scale
  mn->Command(Form("SET LIM 3 %f %f", 
    max(0, npimu_nominal - 5*npimu_error),
           npimu_nominal + 5*npimu_error)); // pion to muon neutron yield ratio
  mn->Command("SET LIM 4 0.01 1"); // neutron efficiency
  mn->Command("SET LIM 5 0.01 1.1"); // B-12 efficiency
  mn->Command("SET LIM 6 0.01 1"); // mu- neutron yield
  mn->Command("MIGRAD");
  mn->Command("IMPROVE");

  update_hists(getpar(5), getpar(4), getpar(3), getpar(2), getpar(1), getpar(0));

  draw();
}
