void nm_summary(const string name)
{
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

  TH2D * dum = new TH2D("dum", "", 100, minx, maxx, 100, 0, 2.5);

  TGraphAsymmErrors g;
  g.SetMarkerStyle(kFullCircle);
  double mindist, y, yeup, yedn, minslc, maxslc;
  while(cin >> mindist >> minslc >> maxslc >> y >> yeup >> yedn){
    printf("%f %f %f\n", y, minslc, maxslc);
    if(!mindistscan && minslc == 0 && maxslc >= 20) continue;
    if(minslc < 1) minslc = 1;

    g.SetPoint(g.GetN(), mindistscan?mindist:(maxslc+minslc)/2, y);
    g.SetPointError(g.GetN()-1, mindistscan?0:(maxslc-minslc)/2+0.25, mindistscan?0:(maxslc-minslc)/2+0.25, yedn, yeup);

    // Put the center point in the RHC-weighted mean number of slices
    if(!mindistscan && minslc == 2 && maxslc == 4) {
      const double mid = 3.49;
      g.SetPoint(g.GetN()-1, mid, y);
      g.SetPointError(g.GetN()-1, mid-1.75, 4.25-mid, yedn, yeup);
    }
    if(!mindistscan && minslc == 8 && maxslc == 12) {
      const double mid = 8.81;
      g.SetPoint(g.GetN()-1, mid, y);
      g.SetPointError(g.GetN()-1, mid-7.75, 12.25-mid, yedn, yeup);
    }

  }

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

  TLegend leg(leftmargin+0.05, 1-topmargin-0.16, leftmargin+0.2, 1-topmargin);
  leg.SetMargin(0.01);
  leg.SetTextFont(42);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(tsize);

  leg.AddEntry((TH1D*)NULL, nm?"#nu_{#mu} RHC":"Neutral Current", "");
  leg.AddEntry((TH1D*)NULL, mindistscan?"Any number of slices":Form("Searching %.0f cell widths from track end", mindist), "");
  leg.Draw();

  if(!mindistscan){

    TF1 * f = new TF1("f", "[0] + [1]*(x-2)", 2, 20);

    g.Fit("f", "", "", 2, 12);
    g.GetFunction("f")->SetLineColor(kGray);

    TGraphAsymmErrors * ideal = new TGraphAsymmErrors;

    const double val = g.GetFunction("f")->Eval(2);

    // Crudely factor out the systematic error by making the 
    // reduced chi2 = 1
    const double fiterror = g.GetFunction("f")->GetParError(0)
     * g.GetFunction("f")->GetChisquare() / g.GetFunction("f")->GetNDF();

    const double edn = sqrt(pow(g.GetErrorYlow(0),  2) + pow(fiterror, 2)),
                 eup = sqrt(pow(g.GetErrorYhigh(0), 2) + pow(fiterror, 2));

    printf("%.2f + %.2f - %.2f (%.2f)\n", val, eup, edn, fiterror);

    ideal->SetPoint(0, 2, val);
    ideal->SetPointError(0, 0, 0, edn, eup);
    ideal->SetMarkerStyle(kFullCircle);
    ideal->SetMarkerColor(kGray);
    ideal->SetLineColor(kGray);
    ideal->Draw("pz");
  }
  
  c1->Print(Form("%s_summary.pdf", name.c_str()));
}
