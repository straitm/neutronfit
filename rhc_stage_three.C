#include <algorithm>
#include <string>

#include "TBranch.h"
#include "TTree.h"
#include "TMinuit.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TError.h"
#include "TFile.h"
#include "TMarker.h"
#include "TGraphAsymmErrors.h"

TMinuit * mn;

#include "common.C"
#include "util.C"
#include "bayes.C"

// selfpileup subtracts off what you get from one event sending high
// energy neutrons from the vertex forward to its own muon's track
// end. This should be some non-zero number, but it's not obvious what
// the value is. It's also different for FHC and RHC *tears-hair-out*.
// Anyway, it can be found by comparing how often the wrong part of an
// event gives us a selected neutron divided by how often, per slice,
// other events do. I've done a low-statistics MC study and found that
// this number is mercifully very small, so we can safely ignore the
// RHC/FHC differences.
static double nominalselfpileup = 1/24.; // multiplied by sliceweight below
double selfpileup = nominalselfpileup;

TGraphAsymmErrors g;

static void fcn(int & np, double * gin, double & chi2, double *par, int flag)
{
  const double intercept = par[0] + selfpileup * par[1];
  const double slope     = par[1];
  chi2 = 0;

  for(int i = 0; i < g.GetN(); i++){
    const double x = g.GetX()[i];
    const double y = g.GetY()[i];
    const double theory = intercept + x*slope;
    const double e = theory > y? g.GetErrorYhigh(i):
                     g.GetErrorYlow(i) == 0? g.GetErrorYhigh(i):
                     fabs(g.GetErrorYlow(i)) > y? y:
                     g.GetErrorYlow(i);

    chi2 += pow((y-theory)/e, 2);
  }

  // penalize heavily for negative intercept. I don't think it is valid
  // to require the slope to be positive, though, because of the way the
  // components play against each other.
  if(intercept < 0) chi2 += pow(intercept/0.01, 2);
}

struct err_t{
  double up, dn;
};


// Find the Bayesian errors for the given central value at the given confidence
// level (0..1).  The ordering principle is greatest likelihood, i.e.  all
// points inside of the confidence band have a higher likelihood than all points
// outside.
//
// Apply a systematic error which is a two-sided gaussian of error systup on
// the top side and systdn on the bottom, with top and bottom each having 50%
// probability.
static err_t bayes(const double besticept, const double CL,
                   const bool verbose, const double systup,
                   const double systdn)
{
  const int N = 4000;
  const double inc = 20./N;
  vector<double> pbyp;
  mn->Command("MIGRAD");
  const double gmin = mn->fAmin;
  mn->Command("SET PRINT -1");
  for(int i = 0; i < N; i++){
    mn->Command("REL 1");
    mn->Command(Form("SET PAR 1 %f", inc*i));
    mn->Command("FIX 1");
    mn->Command("MIGRAD");
    pbyp.push_back(exp(0.5*(gmin - mn->fAmin)));
  }
  mn->Command("REL 1");

  // Apply systematic error and normalize
  convolve_syst(pbyp, inc, systup, systdn);

  vector<double> pbyi = pbyp;
  std::sort   (pbyp.begin(), pbyp.end());
  std::reverse(pbyp.begin(), pbyp.end());

  if(verbose)
    for(unsigned int i = 0; i < pbyi.size(); i++)
      printf("Probdensity %10f %17.14g\n", inc*i, pbyi[i]);

  double acc = 0;
  double cutoff = -1;
  for(unsigned int i = 0; i < pbyp.size(); i++){
    acc += pbyp[i];
    if(acc >= CL){
      cutoff = (pbyp[i] + pbyp[(i>0?i:1)-1])/2;
      break;
    }
  }

  if(cutoff == -1) printf("Failed to find cutoff\n");

  err_t ans;

  for(unsigned int i = 0; i < pbyi.size(); i++){
    if(pbyi[i] > cutoff){
      ans.dn = besticept - inc*(i-0.5);
      break;
    }
  }

  for(unsigned int i = int(besticept/inc); i < pbyi.size(); i++){
    if(pbyi[i] < cutoff){
      ans.up = inc*(i-0.5) - besticept;
      break;
    }
  }

  return ans;
}

static void styletext(TLatex * t, const double tsize)
{
  t->SetTextFont(42);
  t->SetTextSize(tsize*9/10.95); // footnotesize/normalsize for 11pt
}

static void stylearrow(TArrow * a)
{
  a->SetLineWidth(2);
}

void rhc_stage_three(const string name, const string region)
{
  if(name == "compile") return; // to compile only

  gErrorIgnoreLevel = kError;
  gStyle->SetOptStat(0);
  gStyle->SetFrameLineWidth(2);

  const bool mindistscan = name.size() == 2;
  if(mindistscan){
    printf("mindistscan not supported at the moment\n");
    exit(1);
  }

  const bool nm = name.substr(0,2) == "nm";

  const double tsize = 0.06;

  TCanvas * c1 = new TCanvas("c1", "c1");
  c1->SetTickx(1);
  c1->SetTicky(1);
  const double leftmargin = 0.125;
  const double topmargin  = tsize*1.33;
  const double rightmargin= 0.03;
  const double bottommargin=0.13;
  c1->SetMargin(leftmargin, rightmargin, bottommargin, topmargin);

  const double minx = -0.5,
               maxx = 9.0,
               miny = 0;

  TH2D * dum = new TH2D("dum", "", 1, minx, maxx, 1000, miny, 6);
  dum->GetYaxis()->SetTickSize(0.015);
  dum->GetXaxis()->SetTickSize(0.02);

  g.SetMarkerStyle(kFullCircle);
  g.SetLineWidth(2);
  TGraphAsymmErrors gall; // for any number of slices
  gall.SetLineWidth(3);
  gall.SetMarkerStyle(kOpenCircle);
  gall.SetLineColor(kGray+1);
  gall.SetMarkerColor(kGray+1);
  double mindist, y, yeup, yedn, minslc, maxslc, meanslc;
  string cut_dimension_s;

  double systval = 0, systup = 0, systdn = 0;

  while(cin >> mindist >> minslc >> maxslc >> y >> yeup >> yedn
            >> meanslc >> cut_dimension_s){
    TGraphAsymmErrors * G = NULL;

    // Very conservatively take the error on the combined sample
    // with "any" number of slices to be the systematic error.
    // Currently this is small compared to the fit error.
    if(minslc < 1 && maxslc >= 10){
      printf("Got systematic from combined sample\n");
      systval = y;
      systup = yeup;
      systdn = yedn;
      gall.SetPoint(gall.GetN(), meanslc, y);
      gall.SetPointError(gall.GetN()-1,
#if 0
                         meanslc-minslc, maxslc-meanslc,
#else
                         0, 0,
#endif
                         yedn, yeup);
    }
    else{
      g.SetPoint(g.GetN(), meanslc, y);
      g.SetPointError(g.GetN()-1, meanslc-minslc, maxslc-meanslc, yedn, yeup);
    }
  }

  two_or_three_d cut_dimensions =
    cut_dimension_s == two_or_three_d_names[TWOD]  ? TWOD:
    cut_dimension_s == two_or_three_d_names[THREED]? THREED:
    (printf("BAD CUT DIMENSIONS NAME\n"), TWOD);

  selfpileup *= region == "main"?npileup_sliceweight[cut_dimensions]:
                                 npileup_sliceweight_mucatch[cut_dimensions];

  const double legbottom = 1-topmargin-0.16;

  double maxy = dum->GetYaxis()->GetBinLowEdge(dum->GetNbinsY()+1);
  while(gdrawmax(&g) > 0 && gdrawmax(&g) < maxy*legbottom*0.98) maxy *= 0.99;
  maxy = int(maxy*10.)/10.;
  if(maxy < 3) maxy = 3;
  dum->GetYaxis()->SetRangeUser(miny, maxy);

  dum->GetYaxis()->SetTitle(Form("%s scale relative to MC",
                                 nm?"RHC #nu_{#mu}":"NC"));
  dum->GetXaxis()->SetTitle(
    "Effective pile-up slices per spill");
  dum->GetYaxis()->CenterTitle();
  dum->GetXaxis()->CenterTitle();
  dum->GetYaxis()->SetTitleSize(tsize);
  dum->GetXaxis()->SetTitleSize(tsize);
  dum->GetYaxis()->SetLabelSize(tsize);
  dum->GetXaxis()->SetLabelSize(tsize);

  dum->Draw();

  TLine l(minx + (maxx-minx)/500., 1, maxx - (maxx-minx)/500., 1);
  l.SetLineColor(kBlack);
  l.SetLineStyle(kDashed);
  l.Draw();

  g.Draw("pz");
  gall.Draw("pz");

  TLegend leg(leftmargin+0.05, legbottom, leftmargin+0.2, 1-topmargin-0.03);
  leg.SetMargin(0.01);
  leg.SetTextFont(42);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(tsize);

  leg.AddEntry((TH1D*)NULL, Form("%s: %s",
    nm?"#nu_{#mu} RHC":"Neutral current",
    region == "main"?"Main ND":"Muon Catcher"), "");
  leg.AddEntry((TH1D*)NULL,
    Form("%s clusters %.0f cell widths from track end",
         cut_dimensions == TWOD?"2D":"3D", mindist), "");

  const int ans_color = kRed-5;

  gStyle->SetEndErrorSize(4);

  int ierr = 0;
  mn = new TMinuit(2);
  mn->mnparm(0, "icept",   1, 0.01, 0, 0, ierr);
  mn->mnparm(1, "slope", 0.1, 0.01, 0, 0, ierr);
  mn->fGraphicsMode = false;
  mn->SetFCN(fcn);
  mn->Command("SET ERR 1");

#if 0
  mn->Command("MIGRAD");

  // Sensitivity study: What if we got a similar point at pile-up
  // of 0.25 ~ mean 1 slice per spill.
  g.SetPoint(g.GetN(), 0.25, getpar(0) + 0.25*getpar(1));
  g.SetPointError(g.GetN()-1, 0.25, 0.25, g.GetErrorYlow(0), g.GetErrorYhigh(0));
#endif

  const double upvals[3] = {9, 4, 1};
  const double CLs[3] = {0.997300203937, 0.954499736104, 0.682689492137};
  const double upcols[3] = {kRed, kViolet, kBlue};

  TGraphAsymmErrors * ideal = new TGraphAsymmErrors;
  const int Nband = 200;
  for(int u = 0; u < 3; u++){
    mn->Command(Form("SET ERR %f", upvals[u]));
    TGraph * band = new TGraph(Nband*2);
    band->SetFillColorAlpha(upcols[u], 0.2);

    for(int b = 0; b < Nband; b++){
      selfpileup = nominalselfpileup - 15./Nband*(b-10);

      mn->Command("MIGRAD");

      const double val = std::max(0., getpar(0));

      const err_t stat_berr = bayes(val, CLs[u], false, systup, systdn);

      band->SetPoint(b, -selfpileup, val+stat_berr.up);
      band->SetPoint(Nband*2-b-1, -selfpileup, val-stat_berr.dn);
    }

    band->Draw("f");
  }

  mn->Command("SET ERR 1");
  selfpileup = nominalselfpileup;

  mn->Command("MIGRAD");

  const double val = std::max(0., getpar(0));

  printf("Systematics: + %f - %f\n", systup, systdn);

  const err_t stat_berr = bayes(val, CLs[2], true, systup, systdn);

  printf("%s %11s : %.2f + %.2f - %.2f\n", name.c_str(), region.c_str(),
         val, stat_berr.up, stat_berr.dn);

  ideal->SetPoint(0, -selfpileup, val);
  ideal->SetPointError(0, 0, 0, stat_berr.dn, stat_berr.up);

  ideal->SetMarkerStyle(kFullCircle);
  ideal->SetMarkerColor(ans_color);
  ideal->SetLineColor(ans_color);
  ideal->SetLineWidth(2);
  ideal->Draw("p");

  leg.Draw();

#if 0
  // too close to the other arrow
  const float lowlab = maxy/40, highlab = lowlab + maxy/16.;
  TArrow * a = new TArrow(0, 0, 0, lowlab*0.7, 0.01, "<");
  stylearrow(a);
  const int zios_color = kGreen+2;
  a->SetLineColor(zios_color);
  a->Draw();

  TLatex * t = new TLatex(0, lowlab, "Zero intensity, 1 slice");
  styletext(t, tsize);
  t->SetTextColor(zios_color);
  t->Draw();

  // Visually noisy - explain in text
  a = new TArrow(-selfpileup, maxy*0.006, -selfpileup, highlab*0.85, 0.02, "<");
  stylearrow(a);
  a->SetLineColor(ans_color);
  a->Draw();

  TLatex * t = new TLatex(-selfpileup, highlab, "Self-pileup subtracted");
  styletext(t, tsize);
  t->SetTextColor(ans_color);
  t->Draw();
#endif

  TLatex * t = new TLatex(0, 0, "NOvA Preliminary");
  t->SetTextColor(kBlue);
  t->SetTextSize(tsize);
  t->SetTextFont(42);
  t->SetNDC();
  t->SetTextAlign(33);
  t->SetX(1-rightmargin-0.005);
  t->SetY(1-0.02);
  t->Draw();

  c1->Print(Form("stage_three.%s.mindist%d.%s.%s.pdf", name.substr(0, 2).c_str(),
    int(mindist), region.c_str(), two_or_three_d_names[cut_dimensions]));
}
