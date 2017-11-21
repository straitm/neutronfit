#include "common.C"

static double mean_slice(const bool nm, const float minslc, const float maxslc)
{
  if(minslc == maxslc) return minslc;

  TChain rhc("t"), fhc("t");

  for(int i = 0; i < nperiodrhc; i++)
    rhc.Add(inputfiles[i]);
#if 0
  // Do not include FHC for the nu-mu study, since we're answering
  // a question specifically about RHC.  Does this make sense?
  // No, I don't think it does.  We are using a ratio of FHC to RHC
  // to answer a question, so both matter.
  //
  // TODO: What's the effect of the average being different for
  // FHC and RHC?
  if(!nm)
#endif
    for(int i = nperiodrhc; i < nperiod; i++)
      fhc.Add(inputfiles[i]);

  TH1D tmp_r("tmp_r", "", 5000, 0, 50);
  TH1D tmp_f("tmp_f", "", 5000, 0, 50);

  // *Within* the allowed range (minslc to maxslc), find the mean
  rhc.Draw(Form("(nslc-1)*%f + %f * pot * %f >> tmp_r",
      npileup_sliceweight, slc_per_twp_rhc, 1-npileup_sliceweight),
    Form("i == 0 && contained && primary && (nslc-1)*%f + %f * pot * %f >= %f "
                                        "&& (nslc-1)*%f + %f * pot * %f <  %f",
              npileup_sliceweight, slc_per_twp_rhc, 1-npileup_sliceweight, minslc,
              npileup_sliceweight, slc_per_twp_rhc, 1-npileup_sliceweight, maxslc));
  fhc.Draw(Form("(nslc-1)*%f + %f * pot * %f >> tmp_f",
      npileup_sliceweight, slc_per_twp_fhc, 1-npileup_sliceweight),
    Form("i == 0 && contained && primary && (nslc-1)*%f + %f * pot * %f >= %f "
                                        "&& (nslc-1)*%f + %f * pot * %f <  %f",
              npileup_sliceweight, slc_per_twp_fhc, 1-npileup_sliceweight, minslc,
              npileup_sliceweight, slc_per_twp_fhc, 1-npileup_sliceweight, maxslc));

  printf("Getting mean: RHC %d events, FHC %d\n", tmp_r.GetEntries(), tmp_f.GetEntries());

  TH1D & tmp = tmp_r;
  tmp.Add(&tmp_f);

  return tmp.GetMean();
}

void rhc_stage_three(const string name, const string region)
{
  gErrorIgnoreLevel = kError;

  const bool mindistscan = name.size() == 2;
  const bool nm = name.substr(0,2) == "nm";

  const double tsize = 0.06;

  TCanvas * c1 = new TCanvas("c1", "c1");
  const double leftmargin = 0.12;
  const double topmargin  = 0.05;
  const double rightmargin= 0.03;
  const double bottommargin=0.14;
  c1->SetMargin(leftmargin, rightmargin, bottommargin, topmargin);

  const double minx = mindistscan?-0.5:0, maxx = mindistscan?7:15;

  TH2D * dum = new TH2D("dum", "", 100, minx, maxx, 1000, 0, 6);

  TGraphAsymmErrors g;
  g.SetMarkerStyle(kFullCircle);
  double mindist, y, yeup, yedn, minslc, maxslc;

  // Do not want to set the xerrors before fitting because they get used 
  // inappropriately by TGraph::Fit().  Set them afterwards for display.
  vector<double> drawxerr_dn, drawxerr_up;

  double systval = 0, systup = 0, systdn = 0;

  while(cin >> mindist >> minslc >> maxslc >> y >> yeup >> yedn){
    if(minslc < npileup_sliceweight) minslc = npileup_sliceweight;

    // Very conservatively take the error on the combined sample
    // with "any" number of slices to be the systematic error
    if(minslc < 1 && maxslc >= 10){
      printf("Got systematic from combined sample\n");
      systval = y;
      systup = yeup;
      systdn = yedn;
      continue;
    }

    g.SetPoint(g.GetN(), mindistscan?mindist:(maxslc+minslc)/2, y);
    if(mindistscan){
      g.SetPointError(g.GetN()-1, 0, 0, yedn, yeup);
    }
    else{
      // Put the center point in the R(F)HC-weighted mean number of slices for
      // numu (NC)
      const double mid = mean_slice(nm, minslc, maxslc);
      g.SetPoint(g.GetN()-1, mid, y);
      drawxerr_dn.push_back(mid-minslc);
      drawxerr_up.push_back(maxslc-mid);
      g.SetPointError(g.GetN()-1, 0, 0, yedn, yeup);
    }
  }

  double maxy = dum->GetYaxis()->GetBinLowEdge(dum->GetNbinsY()+1);
  while(gdrawmax(&g) > 0 && gdrawmax(&g)*1.5 < maxy) maxy *= 0.8;
  dum->GetYaxis()->SetRangeUser(0, maxy);

  dum->GetYaxis()->SetTitle("Scale relative to MC");
  dum->GetXaxis()->SetTitle(mindistscan?
    "Number of cells widths around track end searched":
    "Effective physics slices per spill");
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

  TLegend leg(leftmargin+0.05, 1-topmargin-0.19, leftmargin+0.2, 1-topmargin-0.01);
  leg.SetMargin(0.01);
  leg.SetTextFont(42);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(tsize);

  leg.AddEntry((TH1D*)NULL, nm?"#nu_{#mu} RHC":"Neutral current", "");
  leg.AddEntry((TH1D*)NULL, mindistscan?"Any number of slices":
    Form("Searching %.0f cell widths from track end", mindist), "");
  leg.AddEntry((TH1D*)NULL, region == "main"?"Main ND":"Muon Catcher", "");

  if(!mindistscan){

    const int ans_color = kRed-5;

    // number of effective slices giving no pileup.
    const double minpileup = npileup_sliceweight;

    TF1 * f = new TF1("f", "[0] + [1]*(x-2)", minpileup, 20);
    f->SetParameters(1, 0.1);
    f->SetParLimits(0, 0, 10);

    TGraphAsymmErrors * ideal = new TGraphAsymmErrors;

    g.Fit("f", "em", "", minpileup, maxx);
    if(g.GetFunction("f") == NULL){
      fprintf(stderr, "Fit failed\n");
    }
    else{
      g.GetFunction("f")->SetLineColor(ans_color);
      g.GetFunction("f")->SetNpx(500);
      g.GetFunction("f")->SetLineWidth(2);

      const double val = g.GetFunction("f")->Eval(minpileup);

      const double eup = sqrt(pow(+MINUIT->fErp[0],2) +
                              pow(systup * val/systval, 2));
      const double edn = sqrt(pow(-MINUIT->fErn[0],2) +
                              pow(systdn * val/systval, 2));

      printf("%s %11s : %.2f + %.2f - %.2f\n", name.c_str(), region.c_str(),
             val, eup, edn);

      ideal->SetPoint(0, minpileup, val);
      ideal->SetPointError(0, 0, 0, edn, eup);

      ideal->SetMarkerStyle(kFullCircle);
      ideal->SetMarkerColor(ans_color);
      ideal->SetLineColor(ans_color);
      ideal->Draw("p");
    }
  }

  for(int i = 0; i < drawxerr_dn.size(); i++)
    g.SetPointError(i, drawxerr_dn[i], drawxerr_up[i],
                       g.GetErrorYlow(i), g.GetErrorYhigh(i));

  leg.Draw();
  
  c1->Print(Form("%s_summary_%s.pdf", name.c_str(), region.c_str()));
}
