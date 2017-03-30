#include "TMinuit.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TError.h"
#include "TRandom.h"

const double maxfitt = 269;
const double additional = 5000;

TCanvas * c1 = new TCanvas;

bool onegoodminos(const int par, const bool no_low_ok)
{
  if((!no_low_ok && gMinuit->fErn[par] == 0) || gMinuit->fErp[par] == 0){
    printf("Incomplete MINOS errors for free par %d\n", par);
    return false;
  }
  return true;
}

void reset_ee(TF1 * ee)
{
  for(int i = 0; i < ee->GetNpar(); i++)
    ee->ReleaseParameter(i);
  ee->SetParameters(3, 1.5e5, 2.1, 2e4, 55, 300, 1.5e5, 300);
  ee->FixParameter(4, 55);
  ee->FixParameter(5, 300);
}

double max(double a, double b)
{
  return a > b? a: b;
}

struct fitanswers{ // TODO
  bool n_good, b12_good;
  double  n_mag,  n_mage_up,  n_mage_dn,
         b12mag, b12mage_up, b12mage_dn;
};

static fitanswers fit(TF1 * ee, TH1D * hist, const bool is_rhc,
                      const double ntrack)
{
  fitanswers ans;

  int status = 0;

  const int c_nneutron_par = 3;
  const int c_nb12_par = 6;

  int nneutron_par = 3;
  int nb12_par = 4;
  int nfree_par = 6;

  reset_ee(ee);

  hist->GetYaxis()->SetRangeUser(0.01, 2e3);

  // Start with the muon lifetime fixed so that it doesn't try to swap
  // with the neutron lifetime.
  nfree_par--;
  nb12_par--;
  nneutron_par--;
  ee->FixParameter(2, 2.1);
  for(int i = 0; i < 8; i++)
    if(0 == (status = hist->Fit("ee", "ql", "e", -210, maxfitt)))
      break;

  // Now that we're (hopefully) converged, let muon lifetime float
  ee->ReleaseParameter(2);
  nfree_par++;
  nb12_par++;
  nneutron_par++;

  // But check if there is almost no B-12, if so, fix it at the best fit
  // because otherwise we'll probably get problems running MINOS.
  const bool fix_b12 = fabs(ee->GetParameter(6))
                     < fabs(ee->GetParameter(0))/10.;
  if(fix_b12){
    ee->FixParameter(6, ee->GetParameter(6));
    nfree_par--;
  }

  for(int i = 0; i < 8; i++)
    if(0 == (status = hist->Fit("ee", "ql", "e", -210, maxfitt)))
      break;

  static bool first = true; 
  ee->Draw("same");
  c1->Print(Form("fit.pdf%s", first?"(":""));
  c1->Update(); c1->Modified();
  first = false;

  // MINOS errors are wrong if we used the abs trick, so limit
  ee->SetParameter(6, fabs(ee->GetParameter(6)));
  ee->SetParameter(3, fabs(ee->GetParameter(3)));
  ee->SetParLimits(6, 0, max(1e5, 10*ee->GetParameter(6)));
  ee->SetParLimits(3, 0, max(1e5, 10*ee->GetParameter(3)));

  if(!status)
    for(int i = 0; i < 2; i++)
      gMinuit->Command(Form("minos 30000 %d", c_nneutron_par+1));
  gMinuit->Command("show min");
  ans.n_good = onegoodminos(nneutron_par, false);


  ans.n_mag = fabs(ee->GetParameter(3)); // ROOT
  ans.n_mage_up =  gMinuit->fErp[nneutron_par]; // MINUIT
  ans.n_mage_dn = -gMinuit->fErn[nneutron_par];

  // release neutron lifetime, which we had held constant for a fair
  // comparison of neutrons, but is a serious nusiance parameter for
  // B-12. The better way to handle this would be a simulatanous fit
  // of all the histograms.
  ee->ReleaseParameter(4);
  ee->SetParLimits(4, 40, 70);
  nb12_par++;
  nfree_par++;

  ee->ReleaseParameter(6);
  ee->SetParLimits(6, 0, max(1e5, 10*ee->GetParameter(6)));
  for(int i = 0; i < 8; i++)
    if(0 == (status = hist->Fit("ee", "ql", "e", -210, maxfitt)))
      break;
  if(!status)
    for(int i = 0; i < 2; i++)
      gMinuit->Command(Form("minos 30000 %d", c_nb12_par+1));
  gMinuit->Command("show min");
  ans.b12_good = onegoodminos(nb12_par, true);

  ee->Draw("same");
  c1->Print(Form("fit.pdf%s", first?"(":""));
  c1->Update(); c1->Modified();

  /* Cheating for sensitivity study! */

  // Assume that we look at cosmic trigger data or some such to get
  // the noise level to high precision.
  ee->FixParameter(0, ee->GetParameter(0));
  nfree_par--;
  nb12_par--;
  nneutron_par--;

  // Put in how many B-12 there ought to be
  ee->SetParameter(6, ntrack * (is_rhc? 0.15: 0.97) * 0.82 * 0.077 * 0.177 * 0.5);
  printf("Generating fake data with B-12 = %f\n", ee->GetParameter(6));

  for(int i = 210 + maxfitt + 1; i <= hist->GetNbinsX(); i++)
    hist->SetBinContent(i, gRandom->Poisson(ee->Eval(hist->GetBinCenter(i))));

  for(int i = 0; i < 8; i++)
    if(0 == (status = hist->Fit("ee", "ql", "e", -210, maxfitt + additional)))
      break;
  if(!status)
    for(int i = 0; i < 2; i++)
      gMinuit->Command(Form("minos 30000 %d", c_nb12_par+1));
  gMinuit->Command("show min");
  ans.b12_good = onegoodminos(nb12_par, true);

  ee->Draw("same");
  c1->Print(Form("fit.pdf%s", first?"(":""));
  c1->Update(); c1->Modified();

  /* End cheating for sensitivity study! */

  ans.b12mag = fabs(ee->GetParameter(6)); // ROOT

  ans.b12mage_up =  gMinuit->fErp[nb12_par]; // MINUIT
  ans.b12mage_dn = (-gMinuit->fErn[nb12_par] != 0 && 
                 -gMinuit->fErn[nb12_par] != 54321.0)?
                 -gMinuit->fErn[nb12_par]: fabs(ee->GetParameter(6));
  printf("b12mag = %f + %f - %f\n", ans.b12mag, ans.b12mage_up, ans.b12mage_dn);
  printf("n_mag  = %f + %f - %f\n",  ans.n_mag,  ans.n_mage_up,  ans.n_mage_dn);

  return ans;
}

static double ratio_error(const double x, const double y,
                          const double xe, const double ye)
{
  return 1/y * sqrt(pow(xe, 2) + x*x*pow(ye/y, 2));
}

void rhc()
{
  TFile * fhcfile = new TFile(
  "/nova/ana/users/mstrait/ndcosmic/period235-type3.root", "Read");
  //"/nova/ana/users/mstrait/ndcosmic/prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period2_v1/all-type3.root", "Read");
  //"/nova/ana/users/mstrait/ndcosmic/prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period3_v1/all.root", "Read");
  //"/nova/ana/users/mstrait/ndcosmic/prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period5_v1_goodruns/all-type3.root", "Read");
  TFile * rhcfile = new TFile("/nova/ana/users/mstrait/ndcosmic/prod_pid_S16-12-07_nd_period6_keepup/770-type3.root"                          , "Read");

  TTree * t5  = (TTree *)fhcfile->Get("t");
  TTree * rt6 = (TTree *)rhcfile->Get("t");

  if(!t5 || !rt6 || t5->IsZombie() || rt6->IsZombie()){
    fprintf(stderr, "Couldn't read something.  See above.\n");
    return;
  }

  // Let the TFile errors go to the screen, then suppress the rest
  gErrorIgnoreLevel = kError;

  TF1 * ee = new TF1("ee",
    "(x <= -1 || x >= 2)*abs([0]) + "
    "(x >= 2)*("
     "abs([1])/[2]    * exp(-x/[2]   ) + "
     "abs([3])/[4]    * exp(-x/[4]   ) "
     //"*(TMath::Erf(sqrt([5]/x))-2/sqrt(TMath::Pi())*sqrt([5]/x)*exp(-[5]/x))"
     "+ abs([6])/29.1e3 * exp(-x/29.1e3)"
    ") + "
    "((x >= -10 && x <= -1) || (x >= 2 && x <= 10))*([7]*abs(abs(x)-10))",
    -210, maxfitt);

  ee->SetParName(0, "flat");
  ee->SetParName(1, "Nmich");
  ee->SetParName(2, "Tmich");
  ee->SetParName(3, "Nneut");
  ee->SetParName(4, "Tneut");
  ee->SetParName(5, "Aneut");
  ee->SetParName(6, "NB12");
  ee->SetParName(7, "pileup");

  reset_ee(ee);
  ee->SetNpx(460);
  ee->SetLineColor(kRed);

  const int nbins_e = 4;
  const double bins_e[nbins_e+1] = {0.5, 1.25, 2.0, 3.0, 6.0 };

  TH2D * dum = new TH2D("dum", "", 100, 0, bins_e[nbins_e], 1000, 0, 10);

  TH2D * fhc_s = new TH2D("fhc_s", "",
    210 + maxfitt + additional, -210, maxfitt + additional,
    nbins_e, bins_e);
  TH2D * rhc_s = (TH2D *)fhc_s->Clone("rhc_s");

  const char * const basecut =
    "primary && type == 3 && timeleft > 269 && timeback > 210 "
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
    "&& abs(trkx)      < 170 && abs(trky)      < 170 && trkz < 1250"
    "&& nclu < 20" // cut very noisy spills
    //"&& nslc <= 10" // reduce pileup
    ;

  // Attempt to agressively reduce neutrons while still getting B-12
  /*
  const std::string ccut = Form("%s"
    "&& t > -210 && t < 269"
    "&& !(t >= -1 && t < 2)"
    "&& nhitx == 1 && nhity == 1 && mindist < 1.99 && dist2 < 4"
    "&& pe > 70 && e < 20",
    basecut);
  */

  const std::string ccut = Form("%s"
     "&& t > -210 && t < 269"
     "&& !(t >= -1 && t < 2)"
     "&& nhitx >= 1 && nhity >= 1 && mindist <= 6"
     "&& pe > 70 && e < 20", basecut);

  t5 ->Draw("slce:t >> fhc_s", ccut.c_str(), "colz");
  rt6->Draw("slce:t >> rhc_s", ccut.c_str(), "colz");
  c1->SetLogy();

  TGraphAsymmErrors * n_result = new TGraphAsymmErrors;
  TGraphAsymmErrors * n_resultbad = new TGraphAsymmErrors;
  n_result->SetMarkerStyle(kOpenCircle);
  n_resultbad->SetMarkerStyle(kOpenCircle);
  n_resultbad->SetLineColor(kRed);
  n_resultbad->SetMarkerColor(kRed);
  n_result->SetName("n_result");
  n_resultbad->SetName("n_resultbad");
  n_result->SetMarkerSize(0.7);
  n_resultbad->SetMarkerSize(0.7);

  TGraphAsymmErrors * b12_result = new TGraphAsymmErrors;
  TGraphAsymmErrors * b12_resultbad = new TGraphAsymmErrors;
  b12_result->SetLineStyle(kDashed);
  b12_resultbad->SetLineStyle(kDashed);
  b12_result->SetMarkerStyle(kOpenSquare);
  b12_resultbad->SetMarkerStyle(kOpenSquare);
  b12_resultbad->SetLineColor(kRed);
  b12_resultbad->SetMarkerColor(kRed);
  b12_result->SetName("b12_result");
  b12_resultbad->SetName("b12_resultbad");
  b12_result->SetMarkerSize(0.7);
  b12_resultbad->SetMarkerSize(0.7);

  TH1D * tcounts_fhc = new TH1D("tcounts_fhc", "",
    fhc_s->GetNbinsY(), fhc_s->GetYaxis()->GetBinLowEdge(1),
    fhc_s->GetYaxis()->GetBinLowEdge(fhc_s->GetNbinsY()+1));
  TH1D * tcounts_rhc = (TH1D *)tcounts_fhc->Clone("tcounts_rhc");

  const std::string tcut = Form("i == 0 && %s", basecut);

  t5 ->Draw("slce >> tcounts_fhc", tcut.c_str());
  rt6->Draw("slce >> tcounts_rhc", tcut.c_str());

  for(int s = 1; s <= nbins_e; s++){
    const double loslce = fhc_s->GetYaxis()->GetBinLowEdge(s);
    const double hislce = fhc_s->GetYaxis()->GetBinLowEdge(s+1);
  
    TH1D * fhc = fhc_s->ProjectionX("fhc", s, s);
    TH1D * rhc = rhc_s->ProjectionX("rhc", s, s);

    const double rt6scale = tcounts_rhc->GetBinContent(s);
    const double  t5scale = tcounts_fhc->GetBinContent(s);
    const double scale = t5scale/rt6scale;

    printf("Scale %f/%f = %f\n", t5scale, rt6scale, scale);

    fitanswers rhc_ans = fit(ee, rhc, true, rt6scale);

    rhc_ans.n_mag      *= scale;
    rhc_ans.n_mage_up  *= scale;
    rhc_ans.n_mage_dn  *= scale;
    rhc_ans.b12mag     *= scale;
    rhc_ans.b12mage_up *= scale;
    rhc_ans.b12mage_dn *= scale;
    
    const fitanswers fhc_ans = fit(ee, fhc, false, t5scale);

    const double n_rat     = rhc_ans.n_mag/fhc_ans.n_mag;
    const double n_rat_err_up =
      ratio_error(rhc_ans.n_mag, fhc_ans.n_mag, rhc_ans.n_mage_up, fhc_ans.n_mage_dn);
    const double n_rat_err_dn =
      ratio_error(rhc_ans.n_mag, fhc_ans.n_mag, rhc_ans.n_mage_dn, fhc_ans.n_mage_up);

    const double b12_rat     = rhc_ans.b12mag/fhc_ans.b12mag;
    const double b12_rat_err_up =
      ratio_error(rhc_ans.b12mag, fhc_ans.b12mag, rhc_ans.b12mage_up, fhc_ans.b12mage_dn);
    const double b12_rat_err_dn =
      ratio_error(rhc_ans.b12mag, fhc_ans.b12mag, rhc_ans.b12mage_dn, fhc_ans.b12mage_up);

    const bool n_good = fhc_ans.n_good && rhc_ans.n_good;
    const bool b12_good = fhc_ans.b12_good && rhc_ans.b12_good;

    printf("%s/%s RHC/FHC neutron (%4.2f-%4.2f)GeV: %.3f + %.3f - %.3f\n",
      rhc_ans.n_good?"Good":"Bad", fhc_ans.n_good?"Good":"Bad", loslce, hislce,
      n_rat, n_rat_err_up, n_rat_err_dn);
    TGraphAsymmErrors * addto = n_good?n_result:n_resultbad;
    addto->SetPoint(addto->GetN(), (loslce+hislce)/2, n_rat);
    addto->SetPointError(addto->GetN()-1, (hislce-loslce)/2,
      (hislce-loslce)/2, n_rat_err_dn, n_rat_err_up);

    printf("%s/%s RHC/FHC B-12    (%4.2f-%4.2f)GeV: %.3f + %.3f - %.3f\n",
      rhc_ans.b12_good?"Good":"Bad", fhc_ans.b12_good?"Good":"Bad", loslce, hislce,
      b12_rat, b12_rat_err_up, b12_rat_err_dn);
    addto = b12_good?b12_result:b12_resultbad;
    addto->SetPoint(addto->GetN(),
      (loslce+hislce)/2 + (hislce-loslce)/20 /* visual shift */, b12_rat);
    addto->SetPointError(addto->GetN()-1, (hislce-loslce)/2,
      (hislce-loslce)/2, b12_rat_err_dn, b12_rat_err_up);
  }

  c1->SetLogy(0);

  dum->Draw();
  n_result->Draw("pz");
  n_resultbad->Draw("pz");
  b12_result->Draw("pz");
  b12_resultbad->Draw("pz");

  c1->Print("fit.pdf)");
}
