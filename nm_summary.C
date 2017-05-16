void nm_summary(const char * const name)
{
  const double tsize = 0.06;

  TCanvas * c1 = new TCanvas("c1", "c1");
  const double leftmargin = 0.12;
  const double topmargin  = 0.05;
  const double rightmargin= 0.03;
  const double bottommargin=0.14;
  c1->SetMargin(leftmargin, rightmargin, bottommargin, topmargin);

  const double minx = -1, maxx = 7;

  TH2D * dum = new TH2D("dum", "", 100, minx, maxx, 100, 0, 2);

  TGraphAsymmErrors g;
  g.SetMarkerStyle(kFullCircle);
  double x, y, yeup, yedn;
  while(cin >> x >> y >> yeup >> yedn){
    g.SetPoint(g.GetN(), x, y);
    g.SetPointError(g.GetN()-1, 0, 0, yedn, yeup);
  }

  dum->GetYaxis()->SetTitle("Scale relative to MC");
  dum->GetXaxis()->SetTitle("Number of cells widths around track end searched");
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

  TLegend leg(0.5, 0.85, 1-rightmargin, 1-topmargin);
  leg.SetMargin(0.01);
  leg.SetTextFont(42);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(tsize);

  const char * const longname = !strcmp(name, "nm")?"#nu_{#mu} RHC":"Neutral Current";

  leg.AddEntry((TH1D*)NULL, longname, "");
  leg.Draw();
  
  c1->Print(Form("%s_summary.pdf", name));
}
