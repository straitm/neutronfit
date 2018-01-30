#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include <algorithm>
#include <fstream>
#include <stdio.h>
#include <vector>
using std::vector;

#include "bayes.C"


int rhc_stage_four(const char * const inbase, const bool nm)
{
  if(inbase == NULL) return 0;

  double scale, prob;
  vector<double> probs[2];
  double oldscale = -1e100;

  int probi = 0;

  double interval = 0;

  std::ifstream finfile(Form("%stmp", inbase));
  if(!finfile.is_open()){
    fprintf(stderr, "File not found\n");
    return 1;
  }

  while(finfile >> scale >> prob){
    if(scale < oldscale){
      interval = oldscale/(probs[0].size()-1);
      printf("Interval = %f\n", interval);

      probi++;
      if(probi > 1){
        fprintf(stderr, "Too many series\n");
        return 1;
      }
    }
    oldscale = scale;
    probs[probi].push_back(prob);
  }

  if(probi < 1){
    fprintf(stderr, "Not enough series, got size1=%d size2=%d\n",
            probs[0].size(), probs[1].size());
    return 1;
  }

  if(probs[0].size() != probs[1].size()){
    fprintf(stderr, "Series don't have the same size! size1=%d size2=%d\n",
            probs[0].size(), probs[1].size());
    return 1;
  }

  vector<double> combined;

  for(unsigned int i = 0; i < probs[0].size(); i++)
    combined.push_back(probs[0][i]*probs[1][i]);

  normalize(combined);
  normalize(probs[0]);
  normalize(probs[1]);

  const ve vemain        = valerr(probs[0], interval);
  const ve vemuoncatcher = valerr(probs[1], interval);
  const ve vecombined    = valerr(combined, interval);

  /*******************************************************************/

  const double tsize = 0.06;
  TCanvas c;
  c.SetCanvasSize(500, 400);
  const double leftmargin = 0.135;
  const double topmargin  = tsize*1.33;
  const double rightmargin= 0.01;
  const double bottommargin=0.10;
  c.SetMargin(leftmargin, rightmargin, bottommargin, topmargin);

  c.SetFrameLineWidth(3);

  gStyle->SetOptStat(0);

  const double ymax = nm?1.9:3.79;

  TH2D dum("dum", "", 3, 0, 1, 1, 1e-10, ymax);

  dum.GetXaxis()->SetTickLength(0);
  dum.GetXaxis()->SetLabelSize(tsize*1.5);
  dum.GetYaxis()->SetLabelSize(tsize);
  dum.GetYaxis()->SetTitleSize(tsize);
  dum.GetYaxis()->SetTitleOffset(1.05);
  dum.GetYaxis()->SetTitle(Form("%s scale relative to MC",
    nm?"RHC #nu_{#mu}":"NC"));
  dum.GetYaxis()->CenterTitle();

  dum.GetXaxis()->SetBinLabel(1, "Main");
  dum.GetXaxis()->SetBinLabel(2, "#mu Catcher");
  dum.GetXaxis()->SetBinLabel(3, "Combined");

  dum.Draw();

  TGraph nominal;
  const double nerr = nm?
    sqrt(pow((0.15+0.20)/2, 2) + pow((0.03+0.08)/2, 2)): // doc-26205
    0.29; // Jeff on Slack, docs 22562 and 23138
  nominal.SetPoint(nominal.GetN(), 0, 1-nerr);
  nominal.SetPoint(nominal.GetN(), 1, 1-nerr);
  nominal.SetPoint(nominal.GetN(), 1, 1+nerr);
  nominal.SetPoint(nominal.GetN(), 0, 1+nerr);
  nominal.SetPoint(nominal.GetN(), 0, 1-nerr);
  nominal.SetFillColor(kGray);
  nominal.Draw("f");

  TLine divider(0.66, 0, 0.66, ymax);
  divider.SetLineWidth(2);
  divider.Draw();

  TGraphAsymmErrors g;
  g.SetPoint(0, 0.18, vemain.val);
  g.SetPoint(1, 0.50, vemuoncatcher.val);
  g.SetPoint(2, 0.84, vecombined.val);
  g.SetPointError(0, 0, 0, vemain.dn, vemain.up);
  g.SetPointError(1, 0, 0, vemuoncatcher.dn, vemuoncatcher.up);
  g.SetPointError(2, 0, 0, vecombined.dn, vecombined.up);
  g.SetLineWidth(3);
  g.SetMarkerStyle(kOpenCircle);
  g.Draw("pz");

  c.SetTicky();

  TLatex t(0, 0, "NOvA Preliminary");
  t.SetTextColor(kBlue);
  t.SetTextSize(tsize);
  t.SetTextFont(42);
  t.SetNDC();
  t.SetTextAlign(33);
  t.SetX(1-rightmargin-0.005);
  t.SetY(1-0.02);
  t.Draw();

  c.RedrawAxis();
  c.Print(Form("%spdf", inbase));

  return 0;
}
