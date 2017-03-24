#include "TMinuit.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"

bool goodminos(const int nfreepar)
{
  for(int i = 0; i < nfreepar; i++)
    if(gMinuit->fErn[i] == 0 || gMinuit->fErp[i] == 0){
      printf("Incomplete MINOS errors for free parameter %d/%d\n", i, nfreepar);
      gMinuit->Command("show min");
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

void fit(TF1 * ee, TH1D * hist, bool & good, double & mag, double & mage_up, double & mage_dn)
{
  int status = 0;

  int nneutron_par = 3;
  int nfree_par = 6;

  reset_ee(ee);

  // Start with the muon lifetime fixed so that it doesn't try to swap with
  // the neutron lifetime.
  ee->FixParameter(2, 2.1);
  for(int i = 0; i < 8; i++) if(0 == (status = hist->Fit("ee", "ql"))) break;

  // Now that we're (hopefully) reasonably converged, let muon lifetime float
  ee->ReleaseParameter(2);

  // But check if there is almost no B-12, if so, fix it at the best fit because
  // otherwise we'll probably get problems running MINOS.
  if(fabs(ee->GetParameter(6)) < fabs(ee->GetParameter(0))/10.){
    ee->FixParameter(6, ee->GetParameter(6));
    nfree_par--;
  }


  for(int i = 0; i < 8; i++) if(0 == (status = hist->Fit("ee", "ql"))) break;
  if(!status) for(int i = 0; i < 2; i++) gMinuit->Command("minos 30000");
  good = goodminos(nfree_par);

  mag = fabs(ee->GetParameter(nneutron_par));
  mage_up =  gMinuit->fErp[nneutron_par];
  mage_dn = -gMinuit->fErn[nneutron_par];
}

void rhc()
{
  TFile * f5 = new TFile("/nova/ana/users/mstrait/ndcosmic/prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period5_v1_goodruns/all.root", "Read");
  //TFile * f5 = new TFile("/nova/ana/users/mstrait/ndcosmic/prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period2_v1/all.root", "Read");
  //TFile * f5 = new TFile("/nova/ana/users/mstrait/ndcosmic/prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period3_v1/all.root", "Read");
  TFile * f6 = new TFile("/nova/ana/users/mstrait/ndcosmic/prod_pid_S16-12-07_nd_period6_keepup/414.root"                          , "Read");

  TTree * t5  = (TTree *)f5->Get("t");
  TTree * rt6 = (TTree *)f6->Get("t");

  TF1 * ee = new TF1("ee", "(x <= -1 || x >= 2)*abs([0]) + "
               "(x >= 2 && x <= 270)*("
                 "abs([1])/[2]*exp(-x/[2]) + "
                 "abs([3])/[4]*exp(-x/[4])*(TMath::Erf(sqrt([5]/x)) - 2/sqrt(TMath::Pi())*sqrt([5]/x)*exp(-[5]/x)) + "
                 "abs([6])*log(2)/20.2e3 * exp(-x*log(2)/20.2e3)"
               ") + "
               "(x <= -1 || (x >= 2 && x <= 10))*([7]*abs(abs(x)-10))", -10, 470);

  reset_ee(ee);
  ee->SetNpx(500);

  const int nbins_e = 9;
  const double bins_e[nbins_e+1] = {0, 1.0, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0, 4, 15 };

  TH2D * dum = new TH2D("dum", "", 100, 0, bins_e[nbins_e], 100, 0, 3);

  TH2D * fhc_s = new TH2D("fhc_s", "", 480, -10, 470, nbins_e, bins_e);
  TH2D * rhc_s = (TH2D *)fhc_s->Clone("rhc_s");

  const char * const basecut = "primary && type == 3 && timeleft > 270 && timeback > 210 "
                               "&& remid > 0.75 && "
                               "trklen > 200 && " // standard
                               //"trklen > 600 && " // long
                               "abs(trkstartx) < 180 && abs(trkstarty) < 180 && trkstartz > 50 && "
                               "abs(trkx)      < 180 && abs(trky)      < 180 && trkz < 1250";

  const std::string ccut = Form("%s && t > -210 && t < 270 && !(t >= -1 && t < 2) && pe > 35 && mindist <= 2", basecut);

  t5 ->Draw("slce:(t < -10)*(t+480) + (t >= -10)*t >> fhc_s", ccut.c_str(), "colz");
  rt6->Draw("slce:(t < -10)*(t+480) + (t >= -10)*t >> rhc_s", ccut.c_str(), "colz");

  TGraphAsymmErrors * result = new TGraphAsymmErrors;
  TGraphAsymmErrors * resultbad = new TGraphAsymmErrors;
  resultbad->SetLineColor(kRed);
  resultbad->SetMarkerColor(kRed);
  result->SetName("result");
  resultbad->SetName("resultbad");

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

    bool goodrhc, goodfhc;
    double rhc_nmag, rhc_nmage_up, rhc_nmage_dn;
    double fhc_nmag, fhc_nmage_up, fhc_nmage_dn;

    fit(ee, rhc, goodrhc, rhc_nmag, rhc_nmage_up, rhc_nmage_dn);
    fit(ee, fhc, goodfhc, fhc_nmag, fhc_nmage_up, fhc_nmage_dn);

    const double rat     = rhc_nmag/fhc_nmag * t5scale/rt6scale;
    const double rat_err_up = sqrt(pow(rhc_nmage_up/fhc_nmag,2) + pow(rhc_nmag/fhc_nmag/fhc_nmag * fhc_nmage_dn,2));
    const double rat_err_dn = sqrt(pow(rhc_nmage_dn/fhc_nmag,2) + pow(rhc_nmag/fhc_nmag/fhc_nmag * fhc_nmage_up,2));

    const bool good = goodfhc && goodrhc;

    printf("%sRHC/FHC (%4.2f-%4.2f)GeV: %.3f + %.3f - %.3f\n", good?"":"Bad: ", loslce, hislce, rat, rat_err_up, rat_err_dn);
    TGraphAsymmErrors * addto = good?result:resultbad;
    addto->SetPoint(addto->GetN(), (loslce+hislce)/2, rat);
    addto->SetPointError(addto->GetN()-1, (hislce-loslce)/2, (hislce-loslce)/2, rat_err_dn, rat_err_up);
  }
 
  dum->Draw();
  result->Draw("pz*");
  resultbad->Draw("pz*");
}
