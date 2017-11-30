#include "common.C"

// selfpileup subtracts off what you get from one event sending high
// energy neutrons from the vertex forward to its own muon's track
// end. This should be some non-zero number, but it's not obvious what
// the value is. It's also different for FHC and RHC *tears-hair-out*.
// Anyway, it can be found by comparing how often the wrong part of an
// event gives us a selected neutron divided by how often, per slice,
// other events do. I've done a low-statistics MC study and found that
// this number is mercifully very small, so we can safely ignore the
// RHC/FHC differences.
const double nominalselfpileup = 1/24. * npileup_sliceweight;
double selfpileup = nominalselfpileup;

TGraphAsymmErrors g;

TMinuit * mn;

static void fcn(int & np, double * gin, double & chi2, double *par, int flag)
{
  const double intercept = par[0] + selfpileup * par[1];
  const double slope     = par[1];
  chi2 = 0;

  for(int i = 0; i < g.GetN(); i++){
    const double x = g.GetX()[i];
    const double y = g.GetY()[i];
    const double theory = intercept + x*slope;
    const double e = theory > y? g.GetErrorYhigh(i): g.GetErrorYlow(i);
    chi2 += pow((y-theory)/e, 2);
  }

  // penalize heavily for negative slope.  Doing this rather than setting
  // a hard limit allows getting MINOS errors.
  if(slope < 0) chi2 += pow(slope/0.001, 2);
}

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

  TH1D slc_r("slc_r", "", 5000, 0, 50);
  TH1D slc_f("slc_f", "", 5000, 0, 50);

  // *Within* the allowed range (minslc to maxslc), find the mean.  Definition
  // must stay in sync with that in pass_intensity() in rhc_stage_one.C.
  rhc.Draw(Form("(nslc-2)*%f + %f * pot >> slc_r",
      npileup_sliceweight, slc_per_twp_rhc*(1-npileup_sliceweight)),
    Form("i == 0 && contained && primary && "
         "abs((nslc-2)*%f + %f * pot  -  %f) < %f",
         npileup_sliceweight, slc_per_twp_rhc*(1-npileup_sliceweight),
         (minslc+maxslc)/2, (maxslc-minslc)/2));

  fhc.Draw(Form("(nslc-2)*%f + %f * pot >> slc_f",
      npileup_sliceweight, slc_per_twp_fhc*(1-npileup_sliceweight)),
    Form("i == 0 && contained && primary && "
         "abs((nslc-2)*%f + %f * pot  -  %f) < %f",
         npileup_sliceweight, slc_per_twp_fhc*(1-npileup_sliceweight),
         (minslc+maxslc)/2, (maxslc-minslc)/2));

  printf("Getting mean for %4.1f-%4.1f slices: RHC %7d events, FHC %7d\n",
         minslc, maxslc, slc_r.GetEntries(), slc_f.GetEntries());

  return (slc_r.GetMean() + slc_f.GetMean())/2;
}

static double fixerr(const double e)
{
  if(e < 0 || e == 54321){
    printf("BAD ERROR\n");
    return 0;
  }
  return e;
}

static static void styletext(TLatex * t, const double tsize)
{
  t->SetTextFont(42);
  t->SetTextSize(tsize*9/10.95); // footnotesize/normalsize for 11pt
}

static static void stylearrow(TArrow * a)
{
  a->SetLineWidth(2);
}

void rhc_stage_three(const string name, const string region)
{
  gErrorIgnoreLevel = kError;
  gStyle->SetOptStat(0);
  gStyle->SetFrameLineWidth(2);

  const bool mindistscan = name.size() == 2;
  const bool nm = name.substr(0,2) == "nm";

  const double tsize = 0.06;

  TCanvas * c1 = new TCanvas("c1", "c1");
  const double leftmargin = 0.125;
  const double topmargin  = tsize*1.33;
  const double rightmargin= 0.03;
  const double bottommargin=0.13;
  c1->SetMargin(leftmargin, rightmargin, bottommargin, topmargin);

  const double minx = mindistscan?-0.5:-selfpileup - 0.5,
               maxx = mindistscan?7:9;

  TH2D * dum = new TH2D("dum", "", 1, minx, maxx, 1000, 0, 6);
  dum->GetYaxis()->SetTickSize(0.015);
  dum->GetXaxis()->SetTickSize(0.02);

  g.SetMarkerStyle(kFullCircle);
  TGraphAsymmErrors gall; // for any number of slices
  gall.SetMarkerStyle(kOpenCircle);
  gall.SetLineColor(kGray);
  gall.SetMarkerColor(kGray);
  double mindist, y, yeup, yedn, minslc, maxslc;

  double systval = 0, systup = 0, systdn = 0;

  while(cin >> mindist >> minslc >> maxslc >> y >> yeup >> yedn){
    TGraphAsymmErrors * G = NULL;

    if(mindistscan){
      g.SetPoint(g.GetN(), mindist, y);
      g.SetPointError(g.GetN()-1, 0, 0, yedn, yeup);
    }
    else{
      // Put the center point in the R(F)HC-weighted mean number of slices for
      // numu (NC)
      const double mid = mean_slice(nm, minslc, maxslc);

      // Very conservatively take the error on the combined sample
      // with "any" number of slices to be the systematic error.
      // Currently this is small compared to the fit error.
      if(minslc < 1 && maxslc >= 10){
        printf("Got systematic from combined sample\n");
        systval = y;
        systup = yeup;
        systdn = yedn;
        gall.SetPoint(gall.GetN(), mid, y);
        gall.SetPointError(gall.GetN()-1, mid-minslc, maxslc-mid, yedn, yeup);
      }
      else{
        g.SetPoint(g.GetN(), mid, y);
        g.SetPointError(g.GetN()-1, mid-minslc, maxslc-mid, yedn, yeup);
      }
    }
  }

  const double legbottom = 1-topmargin-0.15;

  double maxy = dum->GetYaxis()->GetBinLowEdge(dum->GetNbinsY()+1);
  while(gdrawmax(&g) > 0 && gdrawmax(&g) < maxy*legbottom*0.98) maxy *= 0.99;
  maxy = int(maxy*10.)/10.;
  dum->GetYaxis()->SetRangeUser(0, maxy);

  dum->GetYaxis()->SetTitle(Form("%s scale relative to MC",
                                 nm?"RHC #nu_{#mu}":"NC"));
  dum->GetXaxis()->SetTitle(mindistscan?
    "Number of cells widths around track end searched":
    "Effective other physics slices per spill");
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

  TLegend leg(leftmargin+0.05, legbottom, leftmargin+0.2, 1-topmargin-0.01);
  leg.SetMargin(0.01);
  leg.SetTextFont(42);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(tsize);

  leg.AddEntry((TH1D*)NULL, Form("%s: %s",
    nm?"#nu_{#mu} RHC":"Neutral current",
    region == "main"?"Main ND":"Muon Catcher"), "");
  leg.AddEntry((TH1D*)NULL, mindistscan?"Any number of slices":
    Form("Searching %.0f cell widths from track end", mindist), "");

  const int ans_color = kRed-5;

  gStyle->SetEndErrorSize(4);

  if(!mindistscan){
    int ierr = 0;
    mn = new TMinuit(2);
    mn->mnparm(0, "icept",   1, 0.01, 0, 0, ierr);
    mn->mnparm(1, "slope", 0.1, 0.01, 0, 0, ierr);
    mn->Command("SET STRATEGY 2");
    mn->fGraphicsMode = false;
    mn->SetFCN(fcn);
    mn->Command("SET ERR 1");

    double upvals[3] = {9, 4, 1};
    double upcols[3] = {kRed, kViolet, kBlue};

    TGraphAsymmErrors * ideal = new TGraphAsymmErrors;
    const int Nband = 100;
    for(int u = 0; u < 3; u++){
      mn->Command(Form("SET ERR %f", upvals[u]));
      TGraph * band = new TGraph(Nband*2);
      band->SetFillColorAlpha(upcols[u], 0.2);

      for(int b = 0; b < Nband; b++){
        selfpileup = nominalselfpileup - 0.15*(b-10);

        mn->Command("MIGRAD");
        mn->Command("MINOS");

        const double val = getpar(0);
        const double stat_eup = fixerr(+mn->fErp[0]);
        const double stat_edn = fixerr(-mn->fErn[0]);

        const double eup = sqrt(pow(stat_eup            , 2) +
                                pow(systup * val/systval, 2));
        const double edn = sqrt(pow(stat_edn            , 2) +
                                pow(systdn * val/systval, 2));

        if(stat_eup != 0)
          band->SetPoint(b, -selfpileup, val+eup);
        else // MINOS failure - use previous
          band->SetPoint(b, -selfpileup, band->GetY()[b-1]);

        if(stat_edn != 0)
          band->SetPoint(Nband*2-b-1, -selfpileup, val-edn);
        else
          band->SetPoint(Nband*2-b-1, -selfpileup, band->GetY()[Nband*2-(b-1)-1]);
      }

      band->Draw("f");
    }

    mn->Command("SET ERR 1");
    selfpileup = nominalselfpileup;

    mn->Command("MIGRAD");
    mn->Command("MINOS");

    TF1 f("f", Form("[0] + [1]*(x - %f)", -selfpileup), minx, maxx);
    f.SetParameters(getpar(0), getpar(1));
    f.SetLineColor(ans_color);
    f.SetNpx(500);
    f.SetLineWidth(2);

    const double val = f.Eval(-selfpileup);

    const double eup = sqrt(pow(fixerr(+mn->fErp[0]),2) +
                            pow(systup * val/systval, 2));
    const double edn = sqrt(pow(fixerr(-mn->fErn[0]),2) +
                              pow(systdn * val/systval, 2));
    mn->Command("MNCONT 1 2");
    mn->Command("show min");

    printf("%s %11s : %.2f + %.2f - %.2f\n", name.c_str(), region.c_str(),
           val, eup, edn);

    ideal->SetPoint(0, -selfpileup, val);
    ideal->SetPointError(0, 0, 0, edn, eup);

    ideal->SetMarkerStyle(kFullCircle);
    ideal->SetMarkerColor(ans_color);
    ideal->SetLineColor(ans_color);
    ideal->SetLineWidth(2);
    ideal->Draw("p");
    f.Draw("same");
  }

  leg.Draw();

  const float lowlab = maxy/40, highlab = lowlab + maxy/16.;

#if 0 // too close to the other arrow
  TArrow * a = new TArrow(0, 0, 0, lowlab*0.7, 0.01, "<");
  stylearrow(a);
  const int zios_color = kGreen+2;
  a->SetLineColor(zios_color);
  a->Draw();

  TLatex * t = new TLatex(0, lowlab, "Zero intensity, 1 slice");
  styletext(t, tsize);
  t->SetTextColor(zios_color);
  t->Draw();
#endif

  a = new TArrow(-selfpileup, maxy*0.006, -selfpileup, highlab*0.85, 0.02, "<");
  stylearrow(a);
  a->SetLineColor(ans_color);
  a->Draw();

  TLatex * t = new TLatex(-selfpileup, highlab, "Self-pileup subtracted");
  styletext(t, tsize);
  t->SetTextColor(ans_color);
  t->Draw();

  t = new TLatex(0, 0, "NOvA Preliminary");
  t.SetTextSize(tsize);
  t.SetTextFont(42);
  t.SetNDC();
  t.SetTextAlign(33);
  t.SetX(1-rightmargin-0.005);
  t.SetY(1-0.02);
  t.Draw();

  c1->Print(Form("%s_summary_%s.pdf", name.c_str(), region.c_str()));
}
