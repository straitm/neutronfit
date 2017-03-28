#include "TMinuit.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TError.h"

TCanvas * c1 = new TCanvas;

bool goodminos(const int nfreepar, const bool no_low_ok)
{
  for(int i = 0; i < nfreepar; i++)
    if((!no_low_ok && gMinuit->fErn[i] == 0) || gMinuit->fErp[i] == 0){
      printf("Incomplete MINOS errors for free parameter %d/%d\n", i, nfreepar);
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

void fit(TF1 * ee, TH1D * hist, bool & n_good, bool & b12_good, double & n_mag, double & n_mage_up, double & n_mage_dn,
                                             double & b12_mag, double & b12_mage_up, double & b12_mage_dn)
{
  int status = 0;

  int nneutron_par = 3;
  int nb12_par = 4;
  int nfree_par = 6;

  reset_ee(ee);

  hist->GetYaxis()->SetRangeUser(0.01, 2e3);

  // Start with the muon lifetime fixed so that it doesn't try to swap with
  // the neutron lifetime.
  nfree_par--;
  nb12_par--;
  nneutron_par--;
  ee->FixParameter(2, 2.1);
  for(int i = 0; i < 8; i++) if(0 == (status = hist->Fit("ee", "ql"))) break;

  // Now that we're (hopefully) reasonably converged, let muon lifetime float
  ee->ReleaseParameter(2);
  nfree_par++;
  nb12_par++;
  nneutron_par++;

  // But check if there is almost no B-12, if so, fix it at the best fit because
  // otherwise we'll probably get problems running MINOS.
  const bool fix_b12 = fabs(ee->GetParameter(6)) < fabs(ee->GetParameter(0))/10.;
  if(fix_b12){
    ee->FixParameter(6, ee->GetParameter(6));
    nfree_par--;
  }

  for(int i = 0; i < 8; i++) if(0 == (status = hist->Fit("ee", "ql"))) break;

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

  if(!status) for(int i = 0; i < 2; i++) gMinuit->Command("minos 30000");
  n_good = goodminos(nfree_par, false);
  if(!n_good){
    printf("here was the attempt to fit neutrons:\n");
    gMinuit->Command("show min");
  }

  n_mag = fabs(ee->GetParameter(3)); // ROOT
  n_mage_up =  gMinuit->fErp[nneutron_par]; // MINUIT
  n_mage_dn = -gMinuit->fErn[nneutron_par];


  // release neutron lifetime, which we had held constant for a fair
  // comparison of neutrons, but is a serious nusiance parameter for
  // B-12. The better way to handle this would be a simulatanous fit
  // of all the histograms.
  /*ee->ReleaseParameter(4);
  nb12_par++;
  nfree_par++;
  */

  ee->ReleaseParameter(6);
  ee->SetParLimits(6, 0, max(1e5, 10*ee->GetParameter(6)));
  for(int i = 0; i < 8; i++) if(0 == (status = hist->Fit("ee", "ql"))) break;
  if(!status) for(int i = 0; i < 2; i++) gMinuit->Command("minos 30000");
  b12_good = goodminos(nfree_par, true);
  if(!b12_good){
    printf("here was the attempt to fit B-12\n");
    gMinuit->Command("show min");
  }

  ee->Draw("same");
  c1->Print(Form("fit.pdf%s", first?"(":""));
  c1->Update(); c1->Modified();

  b12_mag = fabs(ee->GetParameter(6)); // ROOT

  gMinuit->Command("show min");
  

  b12_mage_up =  gMinuit->fErp[nb12_par]; // MINUIT
  b12_mage_dn = (-gMinuit->fErn[nb12_par] != 0 && 
                 -gMinuit->fErn[nb12_par] != 54321.0)? -gMinuit->fErn[nb12_par]: fabs(ee->GetParameter(6));
  printf("b12_mag = %f + %f - %f\n", b12_mag, b12_mage_up, b12_mage_dn);
  printf("n_mag   = %f + %f - %f\n",   n_mag,   n_mage_up,   n_mage_dn);
}

double ratio_error(const double x, const double y, const double xe, const double ye)
{
  return 1/y * sqrt(pow(xe, 2) + x*x*pow(ye/y, 2));
}

void rhc()
{
  gErrorIgnoreLevel = kError;

  TFile * f5 = new TFile("/nova/ana/users/mstrait/ndcosmic/prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period5_v1_goodruns/all.root-type3.root", "Read");
  //TFile * f5 = new TFile("/nova/ana/users/mstrait/ndcosmic/prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period2_v1/all.root", "Read");
  //TFile * f5 = new TFile("/nova/ana/users/mstrait/ndcosmic/prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period3_v1/all.root", "Read");
  TFile * f6 = new TFile("/nova/ana/users/mstrait/ndcosmic/prod_pid_S16-12-07_nd_period6_keepup/727.root-type3.root"                          , "Read");

  TTree * t5  = (TTree *)f5->Get("t");
  TTree * rt6 = (TTree *)f6->Get("t");

  if(!t5 || !rt6 || t5->IsZombie() || rt6->IsZombie()){
    fprintf(stderr, "Couldn't read something.  Helpful, eh?\n");
    return;
  }

  TF1 * ee = new TF1("ee", "(x <= -1 || x >= 2)*abs([0]) + "
               "(x >= 2)*("
                 "abs([1])/[2]    * exp(-x/[2]   ) + "
                 "abs([3])/[4]    * exp(-x/[4]   )*(TMath::Erf(sqrt([5]/x)) - 2/sqrt(TMath::Pi())*sqrt([5]/x)*exp(-[5]/x)) + "
                 "abs([6])/29.1e3 * exp(-x/29.1e3)"
               ") + "
               "((x >= -10 && x <= -1) || (x >= 2 && x <= 10))*([7]*abs(abs(x)-10))", -210, 270);

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

  const int nbins_e = 5;
  const double bins_e[nbins_e+1] = {0, 1.0, 2.0, 3.0, 4.0, 15 };

  TH2D * dum = new TH2D("dum", "", 100, 0, bins_e[nbins_e], 1000, 0, 10);

  TH2D * fhc_s = new TH2D("fhc_s", "", 480, -210, 270, nbins_e, bins_e);
  TH2D * rhc_s = (TH2D *)fhc_s->Clone("rhc_s");

  const char * const basecut = "primary && type == 3 && timeleft > 270 && timeback > 210 "
                               "&& remid > 0.9 && "
                               "trklen > 200 && " // standard
                               //"trklen > 300 && " // longer
                               "abs(trkstartx) < 180 && abs(trkstarty) < 180 && trkstartz > 50 && "
                               "abs(trkx)      < 170 && abs(trky)      < 170 && trkz < 1250";

  // Attempt to agressively reduce neutrons while still getting some B-12
  const std::string ccut = Form("%s && t > -210 && t < 270 && !(t >= -1 && t < 2) && nhit <= 3 && mindist <= 2 && dist2 < 4 && pe > 35 && e < 20", basecut);

  //const std::string ccut = Form("%s && t > -210 && t < 270 && !(t >= -1 && t < 2) && nhit <= 3 && mindist <= 2 && e < 20 && pe > 35", basecut);

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

  TH1D * tcounts_fhc = new TH1D("tcounts_fhc", "", fhc_s->GetNbinsY(), fhc_s->GetYaxis()->GetBinLowEdge(1),
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

    bool n_goodrhc, n_goodfhc;
    bool b12_goodrhc, b12_goodfhc;
    double rhc_nmag, rhc_nmage_up, rhc_nmage_dn;
    double fhc_nmag, fhc_nmage_up, fhc_nmage_dn;
    double rhc_b12mag, rhc_b12mage_up, rhc_b12mage_dn;
    double fhc_b12mag, fhc_b12mage_up, fhc_b12mage_dn;

    fit(ee, rhc, n_goodrhc, b12_goodrhc, rhc_nmag, rhc_nmage_up, rhc_nmage_dn,
                          rhc_b12mag, rhc_b12mage_up, rhc_b12mage_dn);
    rhc_nmag   *= scale; rhc_nmage_up   *= scale; rhc_nmage_dn   *= scale;
    rhc_b12mag *= scale; rhc_b12mage_up *= scale; rhc_b12mage_dn *= scale;
    
    fit(ee, fhc, n_goodfhc, b12_goodfhc, fhc_nmag, fhc_nmage_up, fhc_nmage_dn,
                          fhc_b12mag, fhc_b12mage_up, fhc_b12mage_dn);

    const double n_rat     = rhc_nmag/fhc_nmag;
    const double n_rat_err_up = ratio_error(rhc_nmag, fhc_nmag, rhc_nmage_up, fhc_nmage_dn);
    const double n_rat_err_dn = ratio_error(rhc_nmag, fhc_nmag, rhc_nmage_dn, fhc_nmage_up);

    const double b12_rat     = rhc_b12mag/fhc_b12mag;
    const double b12_rat_err_up = ratio_error(rhc_b12mag, fhc_b12mag, rhc_b12mage_up, fhc_b12mage_dn);
    const double b12_rat_err_dn = ratio_error(rhc_b12mag, fhc_b12mag, rhc_b12mage_dn, fhc_b12mage_up);

    const bool n_good = n_goodfhc && n_goodrhc;
    const bool b12_good = b12_goodfhc && b12_goodrhc;

    printf("%s/%s RHC/FHC neutron (%4.2f-%4.2f)GeV: %.3f + %.3f - %.3f\n", n_goodrhc?"Good":"Bad", n_goodfhc?"Good":"Bad", loslce, hislce, n_rat, n_rat_err_up, n_rat_err_dn);
    TGraphAsymmErrors * addto = n_good?n_result:n_resultbad;
    addto->SetPoint(addto->GetN(), (loslce+hislce)/2, n_rat);
    addto->SetPointError(addto->GetN()-1, (hislce-loslce)/2, (hislce-loslce)/2, n_rat_err_dn, n_rat_err_up);

    printf("%s/%s RHC/FHC B-12    (%4.2f-%4.2f)GeV: %.3f + %.3f - %.3f\n", b12_goodrhc?"Good":"Bad", b12_goodfhc?"Good":"Bad", loslce, hislce, b12_rat, b12_rat_err_up, b12_rat_err_dn);
    addto = b12_good?b12_result:b12_resultbad;
    addto->SetPoint(addto->GetN(), (loslce+hislce)/2 + (hislce-loslce)/20 /* visual shift */, b12_rat);
    addto->SetPointError(addto->GetN()-1, (hislce-loslce)/2, (hislce-loslce)/2, b12_rat_err_dn, b12_rat_err_up);
  }

  c1->SetLogy(0);

  dum->Draw();
  n_result->Draw("pz");
  n_resultbad->Draw("pz");
  b12_result->Draw("pz");
  b12_resultbad->Draw("pz");

  c1->Print("fit.pdf)");
}
