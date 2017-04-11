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
const int npar = 3;

const int nbins_e = 6;
const double bins_e[nbins_e+1] = {0.5, 1, 1.5, 2.0, 2.5, 3.0, 6.0 };

// Take the neutron yield from pions to be 2.5 times that of muons,
// and they all capture. (R. Madey, et al. Neutrons from nuclear
// capture of negative pions. Phys. Rev. C, 25:3050­306
// XXX to be refined by reading more papers from the 70's
const double npimu_nominal = 2.5;
const double npimu_error = 1.0;

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

static void update_hists(const double npimu, const double nmscale, const double ncscale)
{
  // Probability of getting a neutron from a mu- and a mu+
  const double mum_nyield = 0.15;
  const double mup_nyield = 1e-4;

  // Assume the track in NC events is a pi- this fraction of the time
  // XXX need to get a much better idea of this
  const double piminus_frac = 0.33;

  const double muminus_capture_frac = 0.182;
  const double piminus_relative_nyield = npimu/muminus_capture_frac;

  /* We're going to ignore the corner case of pi+ -> mu+ Michel decays
     making a neutron, since we've fudged lots of other bigger stuff. */

  TH1D * rhc_neutrons = (TH1D*)rhc_reco_numu->Clone("rhc_neutrons");
  rhc_neutrons->Reset();
  rhc_neutrons->Add(reco_nc, ncscale*piminus_frac*piminus_relative_nyield*mum_nyield);
  rhc_neutrons->Add(rhc_reco_numu,    mum_nyield * nmscale /* note */);
  rhc_neutrons->Add(rhc_reco_numubar, mup_nyield);
  rhc_neutrons->Divide(rhc);

  // Estimate of number of neutrons from FHC per muon
  TH1D * fhc_neutrons = (TH1D*)fhc_reco_numu->Clone("fhc_neutrons");
  fhc_neutrons->Reset();
  fhc_neutrons->Add(reco_nc, ncscale*piminus_frac*piminus_relative_nyield*mum_nyield);
  fhc_neutrons->Add(fhc_reco_numu,    mum_nyield);
  fhc_neutrons->Add(fhc_reco_numubar, mup_nyield);
  fhc_neutrons->Divide(fhc);

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
  const double piminus_relative_b12yield = 1/30./muminus_capture_frac;

  // Estimate of number of B-12 from FHC per muon
  TH1D * rhc_b12 = (TH1D*)rhc_reco_numu->Clone("rhc_b12");
  rhc_b12->Reset();
  rhc_b12->Add(reco_nc, ncscale*piminus_frac*piminus_relative_b12yield*mum_b12yield);
  rhc_b12->Add(rhc_reco_numu,    mum_b12yield);
  rhc_b12->Add(rhc_reco_numubar, mup_b12yield);
  rhc_b12->Divide(rhc);

  // Estimate of number of b12 from FHC per muon
  TH1D * fhc_b12 = (TH1D*)fhc_reco_numu->Clone("fhc_b12");
  fhc_b12->Reset();
  fhc_b12->Add(reco_nc, ncscale*piminus_frac*piminus_relative_b12yield*mum_b12yield);
  fhc_b12->Add(fhc_reco_numu,    mum_b12yield);
  fhc_b12->Add(fhc_reco_numubar, mup_b12yield);
  fhc_b12->Divide(fhc);

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

  // penalty term
  chi2 += pow((npimu - npimu_nominal)/npimu_error, 2);

  chi2 += compare(n_result, doublerat_ncomplications);
  chi2 += compare(b12_result, doublerat_b12complications);
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

  // very fake data
  // Eyeballed from numu PRL
  const double peaknc = fhc_reco_numu->GetMaximum() * 0.015;
  for(int i = 1; i < reco_nc->GetNbinsX(); i++)
    reco_nc->SetBinContent(i, peaknc*TMath::Gaus(reco_nc->GetBinCenter(i), 1.3, 1));

  rhc_reco_numubar = (TH1D*)rhc_reco_numubar->Rebin(nbins_e, "rhc_reco_numubar", bins_e);
  rhc_reco_numu    = (TH1D*)rhc_reco_numu   ->Rebin(nbins_e, "rhc_reco_numu"   , bins_e);
  fhc_reco_numubar = (TH1D*)fhc_reco_numubar->Rebin(nbins_e, "fhc_reco_numubar", bins_e);
  fhc_reco_numu    = (TH1D*)fhc_reco_numu   ->Rebin(nbins_e, "fhc_reco_numu"   , bins_e);

  reco_nc    = (TH1D*)reco_nc   ->Rebin(nbins_e, "reco_nc"   , bins_e);

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

void draw()
{
  TCanvas * c1 = new TCanvas("rhc1", "rhc1");
  n_result->Draw("apz");
  
  doublerat_ncomplications->Draw("same");
  
}

void mncommand()
{
  std::string command;
  while(true){
    printf("MINUIT> ");
    if(!getline(std::cin, command)) break;
    if(command == "exit") break;
    mn->Command(command.c_str());
  }
}

void rhc_stage_two(const char * const input = "for_stage_two.C")
{
  set_leo_hists();

  gROOT->Macro(input);

  make_mn();

  mn->Command("MIGRAD");

  draw();
}
