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

static const double tsize = 0.04;

const int npar = 3;

// Ratio of the neutron yield from pions to that of muons, per stop.
const double npimu_nominal = 19.1;
const double npimu_error = 3.7;

const double nm_nominal = 1;
const double nm_error = 0.5;

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

// For drawing
TGraphAsymmErrors * g_n_rhc = new TGraphAsymmErrors;
TGraphAsymmErrors * g_n_fhc = new TGraphAsymmErrors;

// For drawing later
//TGraphAsymmErrors * g_b12_rhc = new TGraphAsymmErrors;
//TGraphAsymmErrors * g_b12_fhc = new TGraphAsymmErrors;

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

static void update_hists(const double npimu, const double nmscale, const double ncscale)
{
  // Probability of getting a neutron from a mu- and a mu+
  const double mum_nyield = 0.181;
  const double mup_nyield = 1e-4;

  // this could be another nuisance parameter
  const double muminus_capture_frac = 0.182;

  /* We're going to ignore the corner case of pi+ -> mu+ Michel decays
     making a neutron, since we've fudged lots of other bigger stuff. */

  reset_hists();

  rhc_neutrons_nc     ->Add(rhc_reco_nc, ncscale*npimu*mum_nyield);
  rhc_neutrons_numu   ->Add(rhc_reco_numu,    mum_nyield * nmscale /* note */);
  rhc_neutrons_numubar->Add(rhc_reco_numubar, mup_nyield);

  rhc_neutrons->Add(rhc_neutrons_nc);
  rhc_neutrons->Add(rhc_neutrons_numu);
  rhc_neutrons->Add(rhc_neutrons_numubar);

  rhc_neutrons->Divide(rhc_tracks);
  rhc_neutrons_nc->Divide(rhc_tracks);
  rhc_neutrons_numu->Divide(rhc_tracks);
  rhc_neutrons_numubar->Divide(rhc_tracks);

  // Estimate of number of neutrons from FHC per muon
  fhc_neutrons_nc->Add(fhc_reco_nc, ncscale*npimu*mum_nyield);
  fhc_neutrons_numu->Add(fhc_reco_numu, mum_nyield);
  fhc_neutrons_numubar->Add(fhc_reco_numubar, mup_nyield);

  fhc_neutrons->Add(fhc_neutrons_nc);
  fhc_neutrons->Add(fhc_neutrons_numu);
  fhc_neutrons->Add(fhc_neutrons_numubar);

  fhc_neutrons->Divide(fhc_tracks);
  fhc_neutrons_nc->Divide(fhc_tracks);
  fhc_neutrons_numu->Divide(fhc_tracks);
  fhc_neutrons_numubar->Divide(fhc_tracks);

  // My toy mc atomic cap frac, nucl cap frac on C-12, and Double Chooz!
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

  // Estimate of number of B-12 from FHC per muon
  rhc_b12_nc->Add(rhc_reco_nc, ncscale*pim_b12yield);
  rhc_b12_numu->Add(rhc_reco_numu,    mum_b12yield);
  rhc_b12_numubar->Add(rhc_reco_numubar, mup_b12yield);

  rhc_b12->Add(rhc_b12_nc);
  rhc_b12->Add(rhc_b12_numu);
  rhc_b12->Add(rhc_b12_numubar);

  rhc_b12->Divide(rhc_tracks);
  rhc_b12_nc->Divide(rhc_tracks);
  rhc_b12_numu->Divide(rhc_tracks);
  rhc_b12_numubar->Divide(rhc_tracks);

  // Estimate of number of b12 from FHC per muon
  fhc_b12_nc->Add(fhc_reco_nc, ncscale*pim_b12yield);
  fhc_b12_numu->Add(fhc_reco_numu,    mum_b12yield);
  fhc_b12_numubar->Add(fhc_reco_numubar, mup_b12yield);

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

static double compare(TGraphAsymmErrors * datgraph, TH1D * predhist)
{
  double chi2 = 0;
  for(int i = 0; i < datgraph->GetN(); i++){
    const double e = datgraph->GetX()[i];
    const double r = datgraph->GetY()[i];
    double rup = datgraph->GetErrorYhigh(i);
    double rdn = datgraph->GetErrorYlow (i);

    if(rdn == 0) rdn = rup;  // protect against bad fits
    if(rup == 0) rup = rdn;  // of course, it'd be good if this didn't trigger

    const int bin = predhist->FindBin(e);

    const double pred = predhist->GetBinContent(bin);
    chi2 += pow((pred - r)/(r > pred?rup:rdn), 2);
  }
  return chi2;
}

static void fcn(__attribute__((unused)) int & np,
  __attribute__((unused)) double * gin, double & chi2, double *par,
  __attribute__((unused)) int flag)  
{
  chi2 = 0;

  // TODO: add covariances

  const double ncscale = par[0];
  const double nmscale = par[1];
  const double npimu = par[2];
  update_hists(npimu, nmscale, ncscale);

  // penalty terms
  chi2 += pow((npimu - npimu_nominal)/npimu_error, 2);

  // TODO: May or may not want this to be here!
  //chi2 += pow((nmscale - nm_nominal)/nm_error, 2);

  chi2 += compare(n_result, doublerat_ncomplications);
  if(useb12) chi2 += compare(b12_result, doublerat_b12complications);
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


void set_leo_hists()
{   
  // I don't need to scale these by POT because everything is done in a ratio
  // between the number of primary tracks in the given beam type and the truth
  // of those tracks.  But it might be confusing if the individual histograms
  // are drawn without normalization.
  fill_hists("prod_pid_R17-03-01-prod3reco.d_nd_genie_nonswap_"
    "fhc_nova_v08_period5_v1/all-type3.root", fhc_reco_numu, fhc_reco_numubar,
    fhc_reco_nc);

  fill_hists(//"prod_pid_R16-12-20-prod3recopreview.b_nd_genie_nonswap_"
    /*"rhc_nova_v08_epoch4a_v3/all-type3.root"*/ "/nova/app/users/mstrait/201704-muons-backport/4-type3.root", rhc_reco_numu, rhc_reco_numubar,
    rhc_reco_nc);

  rhc_reco_numubar=(TH1D*)rhc_reco_numubar->Rebin(nbins_e,"rhc_reco_numubar",bins_e);
  rhc_reco_numu   =(TH1D*)rhc_reco_numu   ->Rebin(nbins_e,"rhc_reco_numu"   ,bins_e);
  fhc_reco_numubar=(TH1D*)fhc_reco_numubar->Rebin(nbins_e,"fhc_reco_numubar",bins_e);
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
  mn->Command("SET ERR 1"); // chi^2

  int mnparmerr = 0;
  mn->mnparm(0, "NCscale", 1, 0.1, 0, 0, mnparmerr);
  mn->mnparm(1, "NMscale", 1, 0.1, 0, 0, mnparmerr);
  mn->mnparm(2, "npimu", npimu_nominal, 0.5, 0, 0, mnparmerr);
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
  stylegraph(g_n_rhc, kRed, kSolid, kOpenSquare, 1, 0.7);
  stylegraph(g_n_fhc, kBlue, kSolid, kOpenCircle, 1, 0.7);

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
    min(2.5, 1.5*max(gdrawmax(g_n_rhc), gdrawmax(g_n_fhc))));
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

  leg = new TLegend(leftmargin+0.1, 0.73, 0.28+0.1, 1-topmargin);
  styleleg(leg);
  leg->AddEntry(g_n_rhc, "RHC data", "lpe");
  leg->AddEntry(rhc_neutrons, "RHC Fit", "l");
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

  TLegend * legf = new TLegend(0.28, 0.73, 0.28 + (0.28-leftmargin), 1-topmargin);
  styleleg(legf);
  legf->AddEntry(g_n_fhc, "FHC data", "lpe");
  legf->AddEntry(fhc_neutrons, "FHC Fit", "l");
  legf->AddEntry(fhc_neutrons_nc, "FHC NC", "l");
  legf->AddEntry(fhc_neutrons_numu, "FHC #nu_{#mu}", "l");
  legf->AddEntry(fhc_neutrons_numubar, "FHC #bar{#nu}_{#mu}", "l");
  legf->Draw();

  c2f->Print("fit_stage_two.pdf");

  //////////////////////////////////////////////////////////////////////
  TCanvas * c3 = new TCanvas("rhc3", "rhc3");
  c3->SetMargin(leftmargin, rightmargin, bottommargin, topmargin);

  TH2D * dum3 = new TH2D("dm", "", 100, 0, 10, 1, 0, 2);
  dum3->GetXaxis()->SetRangeUser(0, 2.5);
  dum3->GetYaxis()->SetRangeUser(0, 1.9);
  dum3->GetXaxis()->SetTitle("NC scale");
  dum3->GetYaxis()->SetTitle("#nu_{#mu} scale");
  dum3->GetXaxis()->CenterTitle();
  dum3->GetYaxis()->CenterTitle();
  dum3->Draw();

  mn->fGraphicsMode = true;

  mn->Command("MNCONT 1 2 99");
  TGraph * cont1 = dynamic_cast<TGraph *>(mn->GetPlot());

  useb12 = false;
  mn->Command("MIGRAD");
  mn->Command("MNCONT 1 2 99");
  TGraph * cont2 = dynamic_cast<TGraph *>(mn->GetPlot());

  useb12 = true;
  mn->Command("MIGRAD");
  mn->Command("FIX 3");
  mn->Command("MNCONT 1 2 99");
  TGraph * cont3 = dynamic_cast<TGraph *>(mn->GetPlot());

  // restore
  mn->Command("REL 3");
  mn->Command("MIGRAD");

  if(cont3 != NULL){
    cont3->SetFillStyle(1001);
    cont3->SetFillColor(kRed);
    cont3->SetLineColor(kRed);
    cont3->Draw("f");
  }
  if(cont1 != NULL){
    cont1->SetFillStyle(1001);
    cont1->SetFillColor(kBlue-2);
    cont1->SetLineColor(kBlue-2);

    // get rid of the visual gap 
    cont1->SetPoint(cont1->GetN(), cont1->GetX()[0], cont1->GetY()[0]);
    cont1->Draw("l");
  }
  if(cont2 != NULL){
    cont2->SetFillStyle(1001);
    cont2->SetFillColor(kBlue);
    cont2->SetLineColor(kBlue);
    cont2->SetLineStyle(kDashed);
    cont2->Draw("l");
  }
  
  mn->fGraphicsMode = false;

  TMarker * bestfit = new TMarker(getpar(0), getpar(1), kFullCircle);
  bestfit->Draw();

  leg = new TLegend(0.3, 0.83, 1-rightmargin, 1-topmargin);
  styleleg(leg);
  leg->SetMargin(0.1);
  leg->AddEntry(cont1,
    Form("1D 68%%, #pi/#mu neutron yield is %.0f#pm%.0f", npimu_nominal, npimu_error),
    "l");
  leg->AddEntry(cont2, "1D 68%, without using ^{12}B", "l");
  leg->AddEntry(cont3, "1D 68%, perfectly known #pi/#mu neutron yield", "f");
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

  set_leo_hists();

  gROOT->Macro(input);

  make_mn();

  //mn->Command("SET LIM 1 0.01 10"); // NC scale
  //mn->Command("SET LIM 2 0.01 10");  // numubar scale
  /*mn->Command(Form("SET LIM 3 %f %f", 
    max(0, npimu_nominal - 5*npimu_error),
           npimu_nominal + 5*npimu_error)); // pion to muon neutron yield ratio */
  mn->Command("MIGRAD");
  mn->Command("IMPROVE");

  update_hists(getpar(2), getpar(1), getpar(0));

  draw();
}
