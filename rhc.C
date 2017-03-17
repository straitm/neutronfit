{
  TFile * f5 = new TFile("/nova/ana/users/mstrait/ndcosmic/prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period5_v1_goodruns/all.root", "Read");
  TFile * f6 = new TFile("/nova/ana/users/mstrait/ndcosmic/prod_pid_S16-12-07_nd_period6_keepup/all.root"                          , "Read");

  TTree * t5  = (TTree *)f5->Get("t");
  TTree * rt6 = (TTree *)f6->Get("t");

  TF1 ee("ee", "(x <= -1 || x >= 2)*abs([0]) + "
               "(x >= 2 && x <= 270)*("
                 "abs([1])/[2]*exp(-x/[2]) + "
                 "abs([3])/[4]*exp(-x/[4])*(TMath::Erf(sqrt([5]/x)) - 2/sqrt(TMath::Pi())*sqrt([5]/x)*exp(-[5]/x)) + "
                 "abs([6])*log(2)/20.2e3 * exp(-x*log(2)/20.2e3)"
               ") + "
               "(x <= -1 || (x >= 2 && x <= 10))*([7]*abs(abs(x)-10))", -10, 470);

  ee.SetParameters(30, 1.5e5, 2.1, 2e4, 50, 300, 1.5e5, 300);
  ee.FixParameter(4, 55);
  ee.FixParameter(5, 300);
  ee.SetNpx(500);

  const int nbins_e = 8;
  const double bins_e[nbins_e+1] = {0, 1.0, 1.5, 2.0, 2.5, 3.0, 4, 5, 15 };

  TH2D * fhc5s = new TH2D("fhc5s", "", 480, -10, 470, nbins_e, bins_e);
  TH2D * rhc6s = (TH2D *)fhc5s->Clone("rhc6s");

  const char * const basecut = "primary && type == 3 && timeleft > 270 && timeback > 210 "
                               "&& remid > 0.75 && trklen > 200 && "
                               "abs(trkstartx) < 180 && abs(trkstarty) < 180 && trkstartz > 50 && "
                               "abs(trkx)      < 180 && abs(trky)      < 180 && trkz < 1250";

  const std::string ccut = Form("          %s && t > -210 && t < 270 && !(t >= -1 && t < 2) && pe > 35", basecut);

  t5 .Draw("slce:(t < -10)*(t+480) + (t >= -10)*t >> fhc5s", ccut.c_str(), "colz");
  rt6.Draw("slce:(t < -10)*(t+480) + (t >= -10)*t >> rhc6s", ccut.c_str(), "colz");

  TGraphErrors * result = new TGraphErrors;

  TH1D * tcounts_fhc5 = new TH1D("tcounts_fhc5", "", fhc5s->GetNbinsY(), fhc5s->GetYaxis()->GetBinLowEdge(1),
                                                     fhc5s->GetYaxis()->GetBinLowEdge(fhc5s->GetNbinsY()+1));
  TH1D * tcounts_rhc6 = (TH1D *)tcounts_fhc5->Clone("tcounts_rhc6");

  const std::string tcut = Form("i == 0 && %s", basecut);

  t5 .Draw("slce >> tcounts_fhc5", tcut.c_str());
  rt6.Draw("slce >> tcounts_rhc6", tcut.c_str());

  for(int s = 1; s <= nbins_e; s++){
    const double loslce = fhc5s->GetYaxis()->GetBinLowEdge(s);
    const double hislce = fhc5s->GetYaxis()->GetBinLowEdge(s+1);
  
    TH1D * fhc5 = fhc5s->ProjectionX("fhc5", s, s);
    TH1D * rhc6 = rhc6s->ProjectionX("rhc6", s, s);

    const double rt6scale = tcounts_rhc6.GetBinContent(s);
    const double  t5scale = tcounts_fhc5.GetBinContent(s);

    for(int i = 0; i < 4; i++) rhc6->Fit("ee", "lq");
    if(0 != (int)rhc6->Fit("ee", "lq")){ printf("Failed to fit RHC\n"); continue; }

    const double rhc6nmag = fabs(ee->GetParameter(3));
    const double rhc6nmage = ee->GetParError(3);

    for(int i = 0; i < 4; i++) fhc5->Fit("ee", "lq");
    if(0 != (int)fhc5->Fit("ee", "lq")){ printf("Failed to fit FHC\n"); continue; }

    const double fhc5nmag = fabs(ee->GetParameter(3));
    const double fhc5nmage = ee->GetParError(3);

    const double rat     = rhc6nmag/fhc5nmag * t5scale/rt6scale;
    const double rat_err = sqrt(pow(rhc6nmage/fhc5nmag,2) + pow(rhc6nmag/fhc5nmag/fhc5nmag * fhc5nmage ,2));

    printf("RHC/FHC (%.1f-%.1f)GeV: %.3f +- %.3f\n", loslce, hislce, rat, rat_err);
    result.SetPoint(result.GetN(), (loslce+hislce)/2, rat);

    result.SetPointError(result.GetN()-1, (hislce-loslce)/2, rat_err);
  }
 
  result.Draw("apz*");
}
