#include "common.C"

double mean_slice(const bool nm, const int minslc, const int maxslc)
{
  if(minslc == maxslc) return minslc;

  TChain c("t");
  if(!nm)
    for(int i = nperiodrhc; i < nperiod; i++)
      c.Add(inputfiles[i]);
  for(int i = 0; i < nperiodrhc; i++)
    c.Add(inputfiles[i]);

  TH1D tmp("tmp", "", 100, 1, 101);
  c.Draw("nslc >> tmp",
    Form("i == 0 && contained && primary && nslc >= %d && nslc <= %d",
         minslc, maxslc));

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

  const double minx = mindistscan?-0.5:0.5, maxx = mindistscan?7:15;

  TH2D * dum = new TH2D("dum", "", 100, minx, maxx, 1000, 0, 6);

  TGraphAsymmErrors g;
  g.SetMarkerStyle(kFullCircle);
  double mindist, y, yeup, yedn, minslc, maxslc;

  // Do not want to set the xerrors before fitting because they get used 
  // in appropriately by TGraph::Fit().  Set them afterwards for display.
  vector<double> drawxerr_dn, drawxerr_up;

  double systval = 0, systup = 0, systdn = 0;

  while(cin >> mindist >> minslc >> maxslc >> y >> yeup >> yedn){
    if(!mindistscan && minslc == 0 && maxslc >= 20) continue;
    if(minslc < 1) minslc = 1;

    // Very conservatively take the error on the combined sample
    // with "any" number of slices to be the systematic error
    if(minslc == 2 && maxslc >= 10){
      systval = y;
      systup = yeup;
      systdn = yedn;
      continue;
    }

    g.SetPoint(g.GetN(), mindistscan?mindist:(maxslc+minslc)/2, y);
    g.SetPointError(g.GetN()-1, mindistscan?0:(maxslc-minslc)/2+0.25,
                    mindistscan?0:(maxslc-minslc)/2+0.25, yedn, yeup);

    // Put the center point in the R(F)HC-weighted mean number of slices for
    // numu (NC)
    if(!mindistscan) {
      const double mid = mean_slice(nm, minslc, maxslc);
      g.SetPoint(g.GetN()-1, mid, y);
      drawxerr_dn.push_back(mid-minslc+0.25);
      drawxerr_up.push_back(maxslc+0.25-mid);
      g.SetPointError(g.GetN()-1, 0, 0, yedn, yeup);
    }
  }

  double maxy = dum->GetYaxis()->GetBinLowEdge(dum->GetNbinsY()+1);
  while(gdrawmax(&g) > 0 && gdrawmax(&g)*1.5 < maxy) maxy *= 0.8;
  dum->GetYaxis()->SetRangeUser(0, maxy);

  dum->GetYaxis()->SetTitle("Scale relative to MC");
  dum->GetXaxis()->SetTitle(mindistscan?
    "Number of cells widths around track end searched":
    "Slices per spill");
  dum->GetYaxis()->CenterTitle();
  dum->GetXaxis()->CenterTitle();
  dum->GetYaxis()->SetTitleSize(tsize);
  dum->GetXaxis()->SetTitleSize(tsize);
  dum->GetYaxis()->SetLabelSize(tsize);
  dum->GetXaxis()->SetLabelSize(tsize);

  dum->Draw();

  TLine l(minx + (maxx-minx)/500., 1, maxx - (maxx-minx)/500., 1);
  l.SetLineColor(kGray);
  l.Draw();

  g.Draw("pz");

  TLegend leg(leftmargin+0.05, 1-topmargin-0.19, leftmargin+0.2, 1-topmargin-0.01);
  leg.SetMargin(0.01);
  leg.SetTextFont(42);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(tsize);

  leg.AddEntry((TH1D*)NULL, nm?"#nu_{#mu} RHC":"Neutral current", "");
  leg.AddEntry((TH1D*)NULL, mindistscan?"Any number of slices":Form("Searching %.0f cell widths from track end", mindist), "");
  leg.AddEntry((TH1D*)NULL, region == "main"?"Main":"Muon Catcher", "");
  leg.Draw();

  if(!mindistscan){

    TF1 * f = new TF1("f", "[0] + [1]*(x-2)", 2, 20);
    f->SetParameters(1, 0);

    TGraphAsymmErrors * ideal = new TGraphAsymmErrors;

    g.Fit("f", "em", "", 2, 13);
    if(g.GetFunction("f") == NULL){
      fprintf(stderr, "Fit failed\n");
    }
    else{
      g.GetFunction("f")->SetLineColor(kGray);

      const double val = g.GetFunction("f")->Eval(2);

      // Crudely factor out the systematic error by making the 
      // reduced chi2 = 1 if it is smaller than 1.
      const double reduced_chi2 = g.GetFunction("f")->GetChisquare() /
                                  g.GetFunction("f")->GetNDF();
      double eup = +MINUIT->fErp[0]*
        (reduced_chi2 > 1? sqrt(reduced_chi2): 1);
      double edn = -MINUIT->fErn[0]*
        (reduced_chi2 > 1? sqrt(reduced_chi2): 1);

      eup = sqrt(pow(eup,2) + pow(systup * val/systval, 2));
      edn = sqrt(pow(edn,2) + pow(systdn * val/systval, 2));

      printf("%s %11s : %.2f + %.2f - %.2f\n", name.c_str(), region.c_str(),
             val, eup, edn);

      ideal->SetPoint(0, 2, val);
      ideal->SetPointError(0, 0, 0, edn, eup);

      ideal->SetMarkerStyle(kFullCircle);
      ideal->SetMarkerColor(kGray);
      ideal->SetLineColor(kGray);
      ideal->Draw("pz");
    }
  }

  for(int i = 0; i < drawxerr_dn.size(); i++)
    g.SetPointError(i, drawxerr_dn[i], drawxerr_up[i], g.GetErrorYlow(i), g.GetErrorYhigh(i));
  
  c1->Print(Form("%s_summary_%s.pdf", name.c_str(), region.c_str()));
}
