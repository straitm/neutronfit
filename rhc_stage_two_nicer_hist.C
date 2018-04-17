TMinuit * mn = NULL; // dummy so util.C compiles

#include "util.C"

// These are the same...
TGraphAsymmErrors * foo_0 = new TGraphAsymmErrors;
TGraphAsymmErrors * foo_1 = new TGraphAsymmErrors;
TGraphAsymmErrors * g_n_fhc = new TGraphAsymmErrors;
TGraphAsymmErrors * g_n_rhc = new TGraphAsymmErrors;

const double leftmargin = 0.12;

void stylepad(TPad * pad)
{
  pad->SetBorderMode(0);
  pad->SetBorderSize(2);
  pad->SetTickx(1);
  pad->SetTicky(1);
  pad->SetFrameLineWidth(2);
  pad->SetFrameBorderMode(0);
  pad->SetLeftMargin(leftmargin);
  pad->SetRightMargin(0.03);
}

void rhc_stage_two_nicer_hist(const string fmacroin)
{
  // Unforgivable ROOT magic in use
  gROOT->Macro(fmacroin.c_str());

  const bool isfhc =
    fmacroin.find("fhc") != string::npos;

  TCanvas *can = new TCanvas("can", "can",0,0,700,700);

  const double divheight = 0.37;

  TPad raw("raw", "raw", 0, divheight, 1, 1);
  raw.Draw();
  TPad rat("raw", "raw", 0, 0, 1, divheight);
  rat.Draw();
  stylepad(&raw);
  stylepad(&rat);
  raw.cd();

  gStyle->SetOptStat(0);

  const double topmargin = 0.09;

  raw.SetTopMargin(topmargin);
  raw.SetBottomMargin(0);

  rat.SetTopMargin(0);
  rat.SetBottomMargin(0.23);

  const double tsize = dm2__1->GetXaxis()->GetLabelSize();

  dm2__1->Draw();

  // Dirty trick to get rid of the half-drawn zero
  TPave coverzero;
  coverzero.SetX1NDC(leftmargin/2);
  coverzero.SetY1NDC(0);
  coverzero.SetX2NDC(leftmargin-0.005);
  coverzero.SetY2NDC(0.03);
  coverzero.SetFillStyle(1001);
  coverzero.SetFillColor(kWhite);
  coverzero.SetBorderSize(0);
  coverzero.Draw();

  TH1D * tot_neut         = isfhc? tot_fhc_neut:      tot_rhc_neut;
  TH1D * neut_numu        = isfhc?     fhc_neut_numu:     rhc_neut_numu;
  TH1D * neut_piflight    = isfhc?     fhc_neut_piflight: rhc_neut_piflight;
  TGraphAsymmErrors * g_n = isfhc? g_n_fhc:           g_n_rhc;


  dm2__1->GetYaxis()->SetRangeUser(0,
    1.1*max(gdrawmax(g_n), tot_neut->GetMaximum()));

  tot_neut->Draw("histsame");
  neut_numu->Draw("histsame");
  neut_piflight->Draw("histsame"); // acutally all NC if NCCOMBINE
  stylegraph(g_n, (isfhc?kBlue:kRed)+3, kSolid,
    isfhc?kOpenCircle:kOpenSquare, 2, 1.0);

  g_n->Draw("pz");

  leg->SetY2(1-topmargin-0.04);
  leg->Draw();
  tex->Draw();
  /********************************************************************/

  rat.cd();

  // Thanks ROOT for making this hard
  const double textratio = (1-divheight)/divheight;

  TH2D * dum = dm2__1->Clone("dum");

  dum->GetYaxis()->SetTitle("Data/Fit");
  dum->GetYaxis()->SetNdivisions(308);

  const double ytitleoff = 1.1;

  dum->GetYaxis()->SetTitleOffset(ytitleoff/textratio);
  dm2__1->GetYaxis()->SetTitleOffset(ytitleoff);

  const double halfwidth = 0.49;

  dum->GetYaxis()->SetRangeUser(1-halfwidth, 1+halfwidth);

  dum->Draw();

  dum->GetXaxis()->SetTitleSize(tsize*textratio);
  dum->GetYaxis()->SetTitleSize(tsize*textratio);
  dum->GetXaxis()->SetLabelSize(tsize*textratio);
  dum->GetYaxis()->SetLabelSize(tsize*textratio);

  TH1D * one = tot_neut->Clone("one");
  one->Divide(one); // :-/
  one->Draw("samehist][");

  TGraphAsymmErrors * datarat = g_n->Clone("datarat");

  for(int i = 0; i < datarat->GetN(); i++){
    const double rawy = datarat->GetY()[i];
    const double rawyeup = datarat->GetErrorYhigh(i);
    const double rawyedn = datarat->GetErrorYlow(i);
    const double dem = tot_neut->GetBinContent(i+1);

    datarat->GetY()[i] = rawy/dem;
    datarat->SetPointError(i,
      datarat->GetErrorXlow(i), datarat->GetErrorXhigh(i),
      rawyedn/dem, rawyeup/dem);
  }

  datarat->Draw("pz");

  const string pdfout =
    fmacroin.substr(0, fmacroin.size()-sizeof("C.ready.C")+1) + "pdf";
  can->SaveAs(pdfout.c_str());
}
