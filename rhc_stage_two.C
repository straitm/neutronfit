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

static TMinuit * mn = NULL;

#include "common.C"

static const double tsize = 0.04;

const int npar = 3;

// Ratio of the neutron yield from pions to that of muons, per stop.
const double npimu_nominal = 19.1;

const double npimu_error = 3.7;

const double nm_nominal = 1;
const double nm_error = 0.2;

bool useb12 = true; // changed between fits

// Leo hists
TH1D *fhc_reco_numubar = new TH1D("fhc_reco_numubar","",100,0,10);
TH1D *fhc_reco_numu = new TH1D("fhc_reco_numu","",100,0,10);
TH1D *rhc_reco_numu = new TH1D("rhc_reco_numu","",100,0,10);
TH1D *rhc_reco_numubar = new TH1D("rhc_reco_numubar","",100,0,10);

// fake hist
TH1D *reco_nc = new TH1D("reco_nc","",100,0,10);

TH1D * fhc = NULL;
TH1D * rhc = NULL;

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
  const double mum_nyield = 0.15;
  const double mup_nyield = 1e-4;

  // Assume the track in NC events is a pi- this fraction of the time
  // Rough estimate from 
  // prod_pid_R17-03-01-prod3reco.d_nd_genie_nonswap_fhc_nova_v08_period5_v1
  const double piminus_frac = 0.20;

  const double muminus_capture_frac = 0.182;
  const double piminus_relative_nyield = npimu/muminus_capture_frac;

  /* We're going to ignore the corner case of pi+ -> mu+ Michel decays
     making a neutron, since we've fudged lots of other bigger stuff. */

  reset_hists();

  rhc_neutrons_nc->Add(reco_nc, ncscale*piminus_frac*piminus_relative_nyield*mum_nyield);
  rhc_neutrons_numu->Add(rhc_reco_numu,    mum_nyield * nmscale /* note */);
  rhc_neutrons_numubar->Add(rhc_reco_numubar, mup_nyield);

  rhc_neutrons->Add(rhc_neutrons_nc);
  rhc_neutrons->Add(rhc_neutrons_numu);
  rhc_neutrons->Add(rhc_neutrons_numubar);

  rhc_neutrons->Divide(rhc);
  rhc_neutrons_nc->Divide(rhc);
  rhc_neutrons_numu->Divide(rhc);
  rhc_neutrons_numubar->Divide(rhc);

  // Estimate of number of neutrons from FHC per muon
  fhc_neutrons_nc->Add(reco_nc, ncscale*piminus_frac*piminus_relative_nyield*mum_nyield);
  fhc_neutrons_numu->Add(fhc_reco_numu,    mum_nyield);
  fhc_neutrons_numubar->Add(fhc_reco_numubar, mup_nyield);

  fhc_neutrons->Add(fhc_neutrons_nc);
  fhc_neutrons->Add(fhc_neutrons_numu);
  fhc_neutrons->Add(fhc_neutrons_numubar);

  fhc_neutrons->Divide(fhc);
  fhc_neutrons_nc->Divide(fhc);
  fhc_neutrons_numu->Divide(fhc);
  fhc_neutrons_numubar->Divide(fhc);

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
  rhc_b12_nc->Add(reco_nc, ncscale*piminus_frac*pim_b12yield);
  rhc_b12_numu->Add(rhc_reco_numu,    mum_b12yield);
  rhc_b12_numubar->Add(rhc_reco_numubar, mup_b12yield);

  rhc_b12->Add(rhc_b12_nc);
  rhc_b12->Add(rhc_b12_numu);
  rhc_b12->Add(rhc_b12_numubar);

  rhc_b12->Divide(rhc);
  rhc_b12_nc->Divide(rhc);
  rhc_b12_numu->Divide(rhc);
  rhc_b12_numubar->Divide(rhc);

  // Estimate of number of b12 from FHC per muon
  fhc_b12_nc->Add(reco_nc, ncscale*piminus_frac*pim_b12yield);
  fhc_b12_numu->Add(fhc_reco_numu,    mum_b12yield);
  fhc_b12_numubar->Add(fhc_reco_numubar, mup_b12yield);

  fhc_b12->Add(fhc_b12_nc);
  fhc_b12->Add(fhc_b12_numu);
  fhc_b12->Add(fhc_b12_numubar);

  fhc_b12->Divide(fhc);
  fhc_b12_nc->Divide(fhc);
  fhc_b12_numu->Divide(fhc);
  fhc_b12_numubar->Divide(fhc);

  doublerat_ncomplications = (TH1D*)rhc_neutrons->Clone("doublerat_ncomplications");
  doublerat_ncomplications->Divide(rhc_neutrons, fhc_neutrons);

  doublerat_b12complications = (TH1D*)rhc_b12->Clone("doublerat_b12complications");
  doublerat_b12complications->Divide(rhc_b12, fhc_b12);
}

static double compare(TGraphAsymmErrors * datgraph, TH1D * predhist)
{
  double chi2 = 0;
  for(int i = 0; i < datgraph->GetN(); i++){
    const double e = datgraph->GetX()[i];
    const double r = datgraph->GetY()[i];
    const double rup = datgraph->GetErrorYhigh(i);
    const double rdn = datgraph->GetErrorYlow (i);
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

  const double ncscale = par[0];
  const double nmscale = par[1];
  const double npimu = par[2];
  update_hists(npimu, nmscale, ncscale);

  // penalty terms
  chi2 += pow((npimu - npimu_nominal)/npimu_error, 2);

  // TODO: May or may not want this to be here!
  chi2 += pow((nmscale - nm_nominal)/nm_error, 2);

  chi2 += compare(n_result, doublerat_ncomplications);
  if(useb12) chi2 += compare(b12_result, doublerat_b12complications);
}


void set_leo_hists()
{   
  /* From Leonidas 2017-03-28 */

  fhc_reco_numubar->SetBinContent(4,33);
  fhc_reco_numubar->SetBinContent(5,236);
  fhc_reco_numubar->SetBinContent(6,570);
  fhc_reco_numubar->SetBinContent(7,786);
  fhc_reco_numubar->SetBinContent(8,1014);
  fhc_reco_numubar->SetBinContent(9,1148);
  fhc_reco_numubar->SetBinContent(10,1174);
  fhc_reco_numubar->SetBinContent(11,1209);
  fhc_reco_numubar->SetBinContent(12,1228);
  fhc_reco_numubar->SetBinContent(13,1273);
  fhc_reco_numubar->SetBinContent(14,1177);
  fhc_reco_numubar->SetBinContent(15,1193);
  fhc_reco_numubar->SetBinContent(16,1151);
  fhc_reco_numubar->SetBinContent(17,1016);
  fhc_reco_numubar->SetBinContent(18,998);
  fhc_reco_numubar->SetBinContent(19,909);
  fhc_reco_numubar->SetBinContent(20,813);
  fhc_reco_numubar->SetBinContent(21,813);
  fhc_reco_numubar->SetBinContent(22,740);
  fhc_reco_numubar->SetBinContent(23,669);
  fhc_reco_numubar->SetBinContent(24,614);
  fhc_reco_numubar->SetBinContent(25,585);
  fhc_reco_numubar->SetBinContent(26,526);
  fhc_reco_numubar->SetBinContent(27,463);
  fhc_reco_numubar->SetBinContent(28,416);
  fhc_reco_numubar->SetBinContent(29,411);
  fhc_reco_numubar->SetBinContent(30,353);
  fhc_reco_numubar->SetBinContent(31,314);
  fhc_reco_numubar->SetBinContent(32,246);
  fhc_reco_numubar->SetBinContent(33,222);
  fhc_reco_numubar->SetBinContent(34,206);
  fhc_reco_numubar->SetBinContent(35,158);
  fhc_reco_numubar->SetBinContent(36,106);
  fhc_reco_numubar->SetBinContent(37,100);
  fhc_reco_numubar->SetBinContent(38,77);
  fhc_reco_numubar->SetBinContent(39,63);
  fhc_reco_numubar->SetBinContent(40,50);
  fhc_reco_numubar->SetBinContent(41,38);
  fhc_reco_numubar->SetBinContent(42,45);
  fhc_reco_numubar->SetBinContent(43,35);
  fhc_reco_numubar->SetBinContent(44,31);
  fhc_reco_numubar->SetBinContent(45,18);
  fhc_reco_numubar->SetBinContent(46,16);
  fhc_reco_numubar->SetBinContent(47,8);
  fhc_reco_numubar->SetBinContent(48,10);
  fhc_reco_numubar->SetBinContent(49,5);
  fhc_reco_numubar->SetBinContent(50,2);
  fhc_reco_numubar->SetBinContent(51,1);
  fhc_reco_numubar->SetBinContent(52,1);

  fhc_reco_numu->SetBinContent(3,3);
  fhc_reco_numu->SetBinContent(4,212.541);
  fhc_reco_numu->SetBinContent(5,1284.79);
  fhc_reco_numu->SetBinContent(6,4248.33);
  fhc_reco_numu->SetBinContent(7,9328.47);
  fhc_reco_numu->SetBinContent(8,15396.5);
  fhc_reco_numu->SetBinContent(9,21537.5);
  fhc_reco_numu->SetBinContent(10,26678.5);
  fhc_reco_numu->SetBinContent(11,32098.6);
  fhc_reco_numu->SetBinContent(12,36620.6);
  fhc_reco_numu->SetBinContent(13,40362.4);
  fhc_reco_numu->SetBinContent(14,43287.9);
  fhc_reco_numu->SetBinContent(15,44585.5);
  fhc_reco_numu->SetBinContent(16,45251.6);
  fhc_reco_numu->SetBinContent(17,44704.7);
  fhc_reco_numu->SetBinContent(18,42626.8);
  fhc_reco_numu->SetBinContent(19,39552.2);
  fhc_reco_numu->SetBinContent(20,35662.1);
  fhc_reco_numu->SetBinContent(21,31291.4);
  fhc_reco_numu->SetBinContent(22,27101.4);
  fhc_reco_numu->SetBinContent(23,22832.5);
  fhc_reco_numu->SetBinContent(24,18872.8);
  fhc_reco_numu->SetBinContent(25,15295.8);
  fhc_reco_numu->SetBinContent(26,12327);
  fhc_reco_numu->SetBinContent(27,9607.72);
  fhc_reco_numu->SetBinContent(28,7526.39);
  fhc_reco_numu->SetBinContent(29,5752.59);
  fhc_reco_numu->SetBinContent(30,4508.04);
  fhc_reco_numu->SetBinContent(31,3414.41);
  fhc_reco_numu->SetBinContent(32,2729.95);
  fhc_reco_numu->SetBinContent(33,2101.8);
  fhc_reco_numu->SetBinContent(34,1573);
  fhc_reco_numu->SetBinContent(35,1225.34);
  fhc_reco_numu->SetBinContent(36,886.063);
  fhc_reco_numu->SetBinContent(37,695.44);
  fhc_reco_numu->SetBinContent(38,663.329);
  fhc_reco_numu->SetBinContent(39,531.704);
  fhc_reco_numu->SetBinContent(40,430.067);
  fhc_reco_numu->SetBinContent(41,400.083);
  fhc_reco_numu->SetBinContent(42,347.758);
  fhc_reco_numu->SetBinContent(43,300.585);
  fhc_reco_numu->SetBinContent(44,259.687);
  fhc_reco_numu->SetBinContent(45,219.636);
  fhc_reco_numu->SetBinContent(46,171.117);
  fhc_reco_numu->SetBinContent(47,92.8628);
  fhc_reco_numu->SetBinContent(48,67.5887);
  fhc_reco_numu->SetBinContent(49,41.9139);
  fhc_reco_numu->SetBinContent(50,29.2711);
  fhc_reco_numu->SetBinContent(51,11);
  fhc_reco_numu->SetBinContent(52,7);
  fhc_reco_numu->SetBinContent(53,3);
  fhc_reco_numu->SetBinContent(54,1);
  fhc_reco_numu->SetBinContent(55,1);
  fhc_reco_numu->SetBinContent(57,1);
  fhc_reco_numu->SetBinContent(77,1);
  fhc_reco_numu->SetBinContent(90,1);

  rhc_reco_numu->SetBinContent(4,48.2483);
  rhc_reco_numu->SetBinContent(5,239.09);
  rhc_reco_numu->SetBinContent(6,681.456);
  rhc_reco_numu->SetBinContent(7,1181.58);
  rhc_reco_numu->SetBinContent(8,1667.23);
  rhc_reco_numu->SetBinContent(9,2118.18);
  rhc_reco_numu->SetBinContent(10,2212.95);
  rhc_reco_numu->SetBinContent(11,2275.83);
  rhc_reco_numu->SetBinContent(12,2293.39);
  rhc_reco_numu->SetBinContent(13,2203.7);
  rhc_reco_numu->SetBinContent(14,2214.57);
  rhc_reco_numu->SetBinContent(15,2063.23);
  rhc_reco_numu->SetBinContent(16,1969.59);
  rhc_reco_numu->SetBinContent(17,1860.93);
  rhc_reco_numu->SetBinContent(18,1655.79);
  rhc_reco_numu->SetBinContent(19,1732.2);
  rhc_reco_numu->SetBinContent(20,1480.12);
  rhc_reco_numu->SetBinContent(21,1370.12);
  rhc_reco_numu->SetBinContent(22,1316.47);
  rhc_reco_numu->SetBinContent(23,1194.72);
  rhc_reco_numu->SetBinContent(24,1094.11);
  rhc_reco_numu->SetBinContent(25,988.928);
  rhc_reco_numu->SetBinContent(26,831.687);
  rhc_reco_numu->SetBinContent(27,813.377);
  rhc_reco_numu->SetBinContent(28,743.707);
  rhc_reco_numu->SetBinContent(29,727.005);
  rhc_reco_numu->SetBinContent(30,590.358);
  rhc_reco_numu->SetBinContent(31,539.509);
  rhc_reco_numu->SetBinContent(32,477.372);
  rhc_reco_numu->SetBinContent(33,419.067);
  rhc_reco_numu->SetBinContent(34,401.82);
  rhc_reco_numu->SetBinContent(35,351.079);
  rhc_reco_numu->SetBinContent(36,330.392);
  rhc_reco_numu->SetBinContent(37,289.085);
  rhc_reco_numu->SetBinContent(38,276.722);
  rhc_reco_numu->SetBinContent(39,196.527);
  rhc_reco_numu->SetBinContent(40,169.621);
  rhc_reco_numu->SetBinContent(41,146.6);
  rhc_reco_numu->SetBinContent(42,144.985);
  rhc_reco_numu->SetBinContent(43,97.0034);
  rhc_reco_numu->SetBinContent(44,87.65);
  rhc_reco_numu->SetBinContent(45,93.3);
  rhc_reco_numu->SetBinContent(46,63.3);
  rhc_reco_numu->SetBinContent(47,36);
  rhc_reco_numu->SetBinContent(48,28);
  rhc_reco_numu->SetBinContent(49,12);
  rhc_reco_numu->SetBinContent(50,10);
  rhc_reco_numu->SetBinContent(51,7);
  rhc_reco_numu->SetBinContent(53,1);
  rhc_reco_numu->SetBinContent(58,1);

  rhc_reco_numubar->SetBinContent(3,2);
  rhc_reco_numubar->SetBinContent(4,122);
  rhc_reco_numubar->SetBinContent(5,862);
  rhc_reco_numubar->SetBinContent(6,2807);
  rhc_reco_numubar->SetBinContent(7,5695);
  rhc_reco_numubar->SetBinContent(8,8336);
  rhc_reco_numubar->SetBinContent(9,10513);
  rhc_reco_numubar->SetBinContent(10,12524);
  rhc_reco_numubar->SetBinContent(11,14365);
  rhc_reco_numubar->SetBinContent(12,16360);
  rhc_reco_numubar->SetBinContent(13,18252);
  rhc_reco_numubar->SetBinContent(14,19786);
  rhc_reco_numubar->SetBinContent(15,20909);
  rhc_reco_numubar->SetBinContent(16,21326);
  rhc_reco_numubar->SetBinContent(17,21288);
  rhc_reco_numubar->SetBinContent(18,20851);
  rhc_reco_numubar->SetBinContent(19,20037);
  rhc_reco_numubar->SetBinContent(20,18176);
  rhc_reco_numubar->SetBinContent(21,16599);
  rhc_reco_numubar->SetBinContent(22,14346);
  rhc_reco_numubar->SetBinContent(23,12270);
  rhc_reco_numubar->SetBinContent(24,10134);
  rhc_reco_numubar->SetBinContent(25,8345);
  rhc_reco_numubar->SetBinContent(26,6884);
  rhc_reco_numubar->SetBinContent(27,5381);
  rhc_reco_numubar->SetBinContent(28,4103);
  rhc_reco_numubar->SetBinContent(29,3208);
  rhc_reco_numubar->SetBinContent(30,2528);
  rhc_reco_numubar->SetBinContent(31,1919);
  rhc_reco_numubar->SetBinContent(32,1506);
  rhc_reco_numubar->SetBinContent(33,1019);
  rhc_reco_numubar->SetBinContent(34,738);
  rhc_reco_numubar->SetBinContent(35,570);
  rhc_reco_numubar->SetBinContent(36,382);
  rhc_reco_numubar->SetBinContent(37,329);
  rhc_reco_numubar->SetBinContent(38,205);
  rhc_reco_numubar->SetBinContent(39,201);
  rhc_reco_numubar->SetBinContent(40,147);
  rhc_reco_numubar->SetBinContent(41,117);
  rhc_reco_numubar->SetBinContent(42,111);
  rhc_reco_numubar->SetBinContent(43,86);
  rhc_reco_numubar->SetBinContent(44,54);
  rhc_reco_numubar->SetBinContent(45,49);
  rhc_reco_numubar->SetBinContent(46,24);
  rhc_reco_numubar->SetBinContent(47,27);
  rhc_reco_numubar->SetBinContent(48,13);
  rhc_reco_numubar->SetBinContent(49,9);
  rhc_reco_numubar->SetBinContent(50,10);
  rhc_reco_numubar->SetBinContent(51,3);
  rhc_reco_numubar->SetBinContent(52,2);
  rhc_reco_numubar->SetBinContent(54,2);

  // From a small sample of MC.  I don't yet know how to count POTs.
  // For that matter, I don't know how many POTs are above.  Actually, if
  // I did everything on my own, I wouldn't need to since ratios would be good 
  // enough(?).
  reco_nc->SetBinContent(7,8);
  reco_nc->SetBinContent(8,16);
  reco_nc->SetBinContent(9,14);
  reco_nc->SetBinContent(10,15);
  reco_nc->SetBinContent(11,24);
  reco_nc->SetBinContent(12,34);
  reco_nc->SetBinContent(13,23);
  reco_nc->SetBinContent(14,21);
  reco_nc->SetBinContent(15,18);
  reco_nc->SetBinContent(16,28);
  reco_nc->SetBinContent(17,20);
  reco_nc->SetBinContent(18,23);
  reco_nc->SetBinContent(19,18);
  reco_nc->SetBinContent(20,12);
  reco_nc->SetBinContent(21,11);
  reco_nc->SetBinContent(22,16);
  reco_nc->SetBinContent(23,16);
  reco_nc->SetBinContent(24,14);
  reco_nc->SetBinContent(25,9);
  reco_nc->SetBinContent(26,9);
  reco_nc->SetBinContent(27,5);
  reco_nc->SetBinContent(28,6);
  reco_nc->SetBinContent(29,3);
  reco_nc->SetBinContent(30,1);
  reco_nc->SetBinContent(31,5);
  reco_nc->SetBinContent(32,5);
  reco_nc->SetBinContent(33,6);
  reco_nc->SetBinContent(34,3);
  reco_nc->SetBinContent(35,5);
  reco_nc->SetBinContent(36,5);
  reco_nc->SetBinContent(37,7);
  reco_nc->SetBinContent(38,2);
  reco_nc->SetBinContent(39,3);
  reco_nc->SetBinContent(40,2);
  reco_nc->SetBinContent(41,5);
  reco_nc->SetBinContent(42,4);
  reco_nc->SetBinContent(43,1);
  reco_nc->SetBinContent(44,1);
  reco_nc->SetBinContent(45,3);
  reco_nc->SetBinContent(46,3);
  reco_nc->SetBinContent(47,2);
  reco_nc->SetBinContent(48,2);
  reco_nc->SetBinContent(49,1);
  reco_nc->SetBinContent(50,5);
  reco_nc->SetBinContent(51,2);
  reco_nc->SetBinContent(52,2);
  reco_nc->SetBinContent(53,3);
  reco_nc->SetBinContent(54,3);
  reco_nc->SetBinContent(55,3);
  reco_nc->SetBinContent(56,2);
  reco_nc->SetBinContent(57,2);
  reco_nc->SetBinContent(58,1);
  reco_nc->SetBinContent(59,1);
  reco_nc->SetBinContent(60,1);
  reco_nc->SetBinContent(61,3);
  reco_nc->SetBinContent(63,1);
  reco_nc->SetBinContent(65,1);
  reco_nc->SetBinContent(69,1);
  reco_nc->SetBinContent(71,1);
  reco_nc->SetBinContent(72,2);
  reco_nc->SetBinContent(73,2);
  reco_nc->SetBinContent(77,1);
  reco_nc->SetBinContent(79,1);
  reco_nc->SetBinContent(80,1);
  reco_nc->SetBinContent(81,1);
  reco_nc->SetBinContent(84,2);
  reco_nc->SetBinContent(85,1);
  reco_nc->SetBinContent(87,1);
  reco_nc->SetBinContent(88,2);
  reco_nc->SetBinContent(91,3);
  reco_nc->SetBinContent(92,1);
  reco_nc->SetBinContent(93,1);
  reco_nc->SetBinContent(94,2);
  reco_nc->SetBinContent(95,2);
  reco_nc->SetBinContent(96,2);
  reco_nc->SetBinContent(98,2);

  // Eyeballed from PRL plot
  const double rescale_nc = fhc_reco_numu->GetMaximum()/reco_nc->GetMaximum() * 0.05;
  reco_nc->Scale(rescale_nc);

  rhc_reco_numubar=(TH1D*)rhc_reco_numubar->Rebin(nbins_e,"rhc_reco_numubar",bins_e);
  rhc_reco_numu   =(TH1D*)rhc_reco_numu   ->Rebin(nbins_e,"rhc_reco_numu"   ,bins_e);
  fhc_reco_numubar=(TH1D*)fhc_reco_numubar->Rebin(nbins_e,"fhc_reco_numubar",bins_e);
  fhc_reco_numu   =(TH1D*)fhc_reco_numu   ->Rebin(nbins_e,"fhc_reco_numu"   ,bins_e);

  reco_nc    = (TH1D*)reco_nc   ->Rebin(nbins_e, "reco_nc"   , bins_e);

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
  rhc = (TH1D*)rhc_reco_numu->Clone("rhc");
  rhc->Add(rhc_reco_numubar);

  // Total numu+numubar in FHC
  fhc = (TH1D*)fhc_reco_numu->Clone("fhc");
  fhc->Add(fhc_reco_numubar);
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
  mn->mnparm(0, "NCscale", 1, 0.010, 0, 0, mnparmerr);
  mn->mnparm(1, "NMscale", 1, 0.010, 0, 0, mnparmerr);
  mn->mnparm(2, "npimu", npimu_nominal, 0.025, 0, 0, mnparmerr);
}

// I believe that TGraph::Integral does a different thing
static double getscale(TGraphAsymmErrors * g, TH1D * h)
{
  g->Fit("pol0", "q0", "");
  const double gscale = g->GetFunction("pol0")->GetParameter(0);

  for(int i = 1; i <= h->GetNbinsX(); i++)
    h->SetBinError(i, g->GetErrorYhigh(i));
  h->Fit("pol0", "q0", "");
  const double hscale = h->GetFunction("pol0")->GetParameter(0);
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
  TCanvas * c1 = new TCanvas("rhc1", "rhc1");
  TH2D * dum = new TH2D("dm", "", 100, 0, 10, 10000, 0, 10);
  TH2D * dum2 = (TH2D*) dum->Clone("dm2");
  dum->GetXaxis()->SetRangeUser(bins_e[0], bins_e[nbins_e]);
  dum->GetYaxis()->SetRangeUser(0, 
    min(1.99, 1.05*max(gdrawmax(  n_result), gdrawmax(b12_result)))
  );
  dum->GetXaxis()->CenterTitle();
  dum->GetYaxis()->CenterTitle();
  dum->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  dum->GetYaxis()->SetTitle("RHC/FHC");
  dum->Draw();

  n_result->SetMarkerStyle(kFullCircle);
  b12_result->SetMarkerStyle(kOpenCircle);
  b12_result->SetMarkerColor(kGray+1);
  b12_result->SetLineColor(kGray+1);

  doublerat_b12complications->SetMarkerColor(kGray+1);
  doublerat_b12complications->SetLineColor(kGray+1);

  doublerat_b12complications->SetLineStyle(kDashed);
  doublerat_ncomplications->SetLineStyle(kDashed);

  b12_result->Draw("pz");
  doublerat_b12complications->Draw("same");
  n_result->Draw("pz");
  doublerat_ncomplications->Draw("same");

  TLegend * leg = new TLegend(0.7, 0.7, 0.95, 0.97);
  leg->AddEntry(n_result, "Neutrons, data", "lpe");
  leg->AddEntry(doublerat_ncomplications, "Neutrons, fit", "l");
  leg->AddEntry(b12_result, "^{12}B, data", "lpe");
  leg->AddEntry(doublerat_b12complications, "^{12}B, fit", "l");
  styleleg(leg);
  leg->Draw();

  c1->Print("fit_stage_two.pdf(");
  
  TCanvas * c2 = new TCanvas("rhc2", "rhc2");
  dum2->GetXaxis()->SetRangeUser(bins_e[0], bins_e[nbins_e]);
  dum2->GetYaxis()->SetRangeUser(0, 
    min(2.5, 1.5*max(gdrawmax(g_n_rhc), gdrawmax(g_n_fhc)))
  );
  dum2->Draw();
  dum2->GetYaxis()->SetTitle("Neutrons per track");
  dum2->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  dum2->GetYaxis()->CenterTitle();
  dum2->GetXaxis()->CenterTitle();

  g_n_rhc->SetLineColor(kRed);
  g_n_rhc->SetMarkerColor(kRed);
  g_n_rhc->SetMarkerStyle(kOpenSquare);

  g_n_fhc->SetLineColor(kBlue);
  g_n_fhc->SetMarkerColor(kBlue);
  g_n_fhc->SetMarkerStyle(kOpenCircle);

  g_n_rhc->Draw("pz");
  g_n_fhc->Draw("pz");

  const double rhcscale = getscale(g_n_rhc, rhc_neutrons);
  const double fhcscale = getscale(g_n_fhc, fhc_neutrons);

  rhc_neutrons->SetLineWidth(2);
  rhc_neutrons_nc->SetLineWidth(2);
  rhc_neutrons_numu->SetLineWidth(2);
  rhc_neutrons_numubar->SetLineWidth(2);
  rhc_neutrons->SetLineColor(kRed);
  rhc_neutrons_nc->SetLineColor(kRed);
  rhc_neutrons_numu->SetLineColor(kRed);
  rhc_neutrons_numubar->SetLineColor(kRed);
  rhc_neutrons_nc->SetLineStyle(kDotted);
  rhc_neutrons_numu->SetLineStyle(kDashed);
  rhc_neutrons_numubar->SetLineStyle(9);

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

  rhc_neutrons->Draw("histsame][");
  fhc_neutrons->Draw("histsame][");
  rhc_neutrons_nc->Draw("histsame][");
  fhc_neutrons_nc->Draw("histsame][");
  rhc_neutrons_numu->Draw("histsame][");
  fhc_neutrons_numu->Draw("histsame][");
  rhc_neutrons_numubar->Draw("histsame][");
  fhc_neutrons_numubar->Draw("histsame][");

  leg = new TLegend(0.14, 0.74, 0.3, 0.97);
  styleleg(leg);
  leg->AddEntry(g_n_rhc, "RHC data", "lpe");
  leg->AddEntry(rhc_neutrons, "RHC Fit", "l");
  leg->AddEntry(rhc_neutrons_nc, "RHC NC", "l");
  leg->AddEntry(rhc_neutrons_numu, "RHC #nu_{#mu}", "l");
  leg->AddEntry(rhc_neutrons_numubar, "RHC #bar{#nu}_{#mu}", "l");
  leg->Draw();
  
  TLegend * legf = new TLegend(0.3, 0.74, 0.46, 0.97);
  styleleg(legf);
  legf->AddEntry(g_n_fhc, "FHC data", "lpe");
  legf->AddEntry(fhc_neutrons, "FHC Fit", "l");
  legf->AddEntry(fhc_neutrons_nc, "FHC NC", "l");
  legf->AddEntry(fhc_neutrons_numu, "FHC #nu_{#mu}", "l");
  legf->AddEntry(fhc_neutrons_numubar, "FHC #bar{#nu}_{#mu}", "l");
  legf->Draw();

  c2->Print("fit_stage_two.pdf");

  TCanvas * c3 = new TCanvas("rhc3", "rhc3");

  TH2D * dum3 = (TH2D*) dum2->Clone("dm3");
  dum3->GetXaxis()->SetRangeUser(0, 1.5);
  dum3->GetYaxis()->SetRangeUser(0, 1.5);
  dum3->GetXaxis()->SetTitle("NC scale");
  dum3->GetYaxis()->SetTitle("#nu_{#mu} scale");
  dum3->Draw();

  mn->fGraphicsMode = true;

  mn->Command("MNCONT 1 2 99");
  TGraph * cont1 = (TGraph *)mn->GetPlot();
  cont1->SetFillStyle(1001);
  cont1->SetFillColor(kBlue-2);
  cont1->SetLineColor(kBlue-2);

  useb12 = false;
  mn->Command("MIGRAD");
  mn->Command("MNCONT 1 2 99");
  TGraph * cont2 = (TGraph *)mn->GetPlot();
  cont2->SetFillStyle(1001);
  cont2->SetFillColor(kBlue);
  cont2->SetLineColor(kBlue);
  cont2->SetLineStyle(kDashed);

  useb12 = true;
  mn->Command("MIGRAD");
  mn->Command("FIX 3");
  mn->Command("MNCONT 1 2 99");
  TGraph * cont3 = (TGraph *)mn->GetPlot();
  cont3->SetFillStyle(1001);
  cont3->SetFillColor(kRed);
  cont3->SetLineColor(kRed);

  cont3->Draw("f");
  cont1->Draw("f");
  cont2->Draw("f");
  
  
  mn->fGraphicsMode = false;

  TMarker * bestfit = new TMarker(getpar(0), getpar(1), kFullCircle);
  bestfit->Draw();

  leg = new TLegend(0.3, 0.85, 0.96, 0.98);
  styleleg(leg);
  leg->AddEntry(cont1,
    Form("1D 68\%, #pi/#mu neutron yield is %.0f#pm%.0f", npimu_nominal, npimu_error),
    "f");
  leg->AddEntry(cont2, "1D 68\%, perfectly known #pi/#mu neutron yield", "f");
  leg->Draw();

  c3->Print("fit_stage_two.pdf)");
}


void rhc_stage_two(const char * const input = "for_stage_two.C")
{
  set_leo_hists();

  gROOT->Macro(input);

  make_mn();

  mn->Command("SET LIM 1 0 5");
  mn->Command("SET LIM 2 0 5");
  mn->Command("SET LIM 3 0 50");
  mn->Command("MIGRAD");

  update_hists(getpar(2), getpar(1), getpar(0));

  draw();
}
