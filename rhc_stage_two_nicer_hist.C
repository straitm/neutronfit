#include "util.C"

// These are the same...
TGraphAsymmErrors * foo_1 = new TGraphAsymmErrors;
TGraphAsymmErrors * g_n_fhc = new TGraphAsymmErrors;
TGraphAsymmErrors * g_n_rhc = new TGraphAsymmErrors;

void rhc_stage_two_nicer_hist(const string fmacroin)
{
  // Unforgivable ROOT magic in use
  gROOT->Macro(fmacroin.c_str());

  TCanvas *can = new TCanvas("can", "can",0,0,700,500);
  gStyle->SetOptStat(0);
  can->SetBorderMode(0);
  can->SetBorderSize(2);
  can->SetTickx(1);
  can->SetTicky(1);
  can->SetLeftMargin(0.135);
  can->SetRightMargin(0.03);
  can->SetTopMargin(0.07315);
  can->SetBottomMargin(0.13);
  can->SetFrameLineWidth(2);
  can->SetFrameBorderMode(0);


  dm2->Draw();

  tot_fhc_neut->Draw("histsame");
  fhc_neut_numu->Draw("histsame");
  fhc_neut_piflight->Draw("histsame"); // acutally all NC if NCCOMBINE
  g_n_fhc->Draw("pz");
  stylegraph(g_n_fhc, kBlue+3, kSolid, kOpenCircle, 2, 1.0);
  ci = TColor::GetColor("#000066");
  g_n_fhc->SetMarkerColor(ci);
  g_n_fhc->SetMarkerStyle(24);
  leg->Draw();
  tex->Draw();

  const string pdfout = 
    fmacroin.substr(0, fmacroin.size()-sizeof("C.ready.C")+1) + "pdf";
  can->SaveAs(pdfout.c_str());
}
