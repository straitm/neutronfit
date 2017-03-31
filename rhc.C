#include "TMinuit.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TError.h"
#include "TRandom.h"

/*
 * TODO: exclude neutron and B-12 candidates following a Michel.
 */


struct fitanswers{
  bool n_good, b12_good;
  double  n_mag,  n_mage_up,  n_mage_dn,
         b12mag, b12mage_up, b12mage_dn;
};

const double nnegbins = 209;
const double maxrealtime = 269;
const double additional = 0;

double maxfitt = maxrealtime; // varies

const int nbins_e = 4;
const double bins_e[nbins_e+1] = {0.5, 1.25, 2.0, 3.0, 6.0 };

static TMinuit * mn = NULL;

const int npar = 8;
const char * const parnames[npar] = {
  "flat", "NMich", "Tmich", "Nneut", "Tneut", "Aneut", "NB12", "pileup" };
const double parinit[npar] = {
  3, 1.5e5, 2.1, 2e4, 55, 300, 1.5e5, 300 };
const int flat_nc   = 0, // parameter numbers, C numbering
          nmich_nc  = 1, // for use with TMinuit functions
          tmich_nc  = 2,
          nneut_nc  = 3,
          tneut_nc  = 4,
          aneut_nc  = 5,
          nb12_nc   = 6,
          pileup_nc = 7;
const int flat_nf   = flat_nc  +1, // and FORTRAN numbering,
          nmich_nf  = nmich_nc +1, // for use with native MINUIT
          tmich_nf  = tmich_nc +1, // commands through TMinuit::Command()
          nneut_nf  = nneut_nc +1,
          tneut_nf  = tneut_nc +1,
          aneut_nf  = aneut_nc +1,
          nb12_nf   = nb12_nc  +1,
          pileup_nf = pileup_nc+1;

static TH2D * fhc_s = NULL;
static TH2D * rhc_s = NULL;

static TH1D * fithist = NULL;

static TF1 * ee = new TF1("ee",
  "(x <= -1 || x >= 2)*abs([0]) + "
  "(x >= 2)*("
   "abs([1])/[2]    * exp(-x/[2]   ) + "
   "abs([3])/[4]    * exp(-x/[4]   ) "
   "*(TMath::Erf(sqrt([5]/x))-2/sqrt(TMath::Pi())*sqrt([5]/x)*exp(-[5]/x))"
   "+ abs([6])/29.1e3 * exp(-x/29.1e3)"
  ") + "
  "((x >= -10 && x <= -1) || (x >= 2 && x <= 10))*(abs([7])*abs(abs(x)-10))",
  -nnegbins, maxrealtime+additional);

static TCanvas * c1 = new TCanvas;

double getminerrup(int i) // 0-indexed!
{
  // Figure out how many fixed parameters there are below the one
  // desired. Fixed parameters are, of course, FORTRAN-indexed, but the
  // array fErp is C-indexed.

  // If there are no fixed parameters, it is just the C-index
  int wherearewe = i;

  for(int fixpar = 0; fixpar < mn->fNpfix; fixpar++)
    if(mn->fIpfix[fixpar] < i+1)
      wherearewe--;

  return mn->fErp[wherearewe];
}

double getminerrdn(int i) // 0-indexed!
{
  int wherearewe = i;
  for(int fixpar = 0; fixpar < mn->fNpfix; fixpar++)
    if(mn->fIpfix[fixpar] < i+1)
      wherearewe--;
  return -mn->fErn[wherearewe];
}

double getpar(int i) // 0-indexed!
{
  double answer, dum;
  mn->GetParameter(i, answer, dum);
  return answer;
}

double geterr(int i) // 0-indexed!
{
  double val, err;
  mn->GetParameter(i, val, err);
  return err;
}

static double getlimlo(int i) // 0-indexed!
{
  double answer, dum;
  int idum;
  TString sdum;
  mn->mnpout(i, sdum, dum, dum, answer, dum, idum);
  return answer;
}

static double getlimup(int i) // 0-indexed!
{
  double answer, dum;
  int idum;
  TString sdum;
  mn->mnpout(i, sdum, dum, dum, dum, answer, idum);
  return answer;
}

static void fixat(int i, float v) // 1-indexed!
{
  mn->Command(Form("REL %d", i));
  if(getlimup(i-1)) mn->Command(Form("SET LIM %d", i));
  mn->Command(Form("SET PAR %d %g", i, v));
  mn->Command(Form("FIX %d", i));
}

static void fixatzero(int i) // 1-indexed!
{
  mn->Command(Form("REL %d", i));
  mn->Command(Form("SET LIM %d", i));
  mn->Command(Form("SET PAR %d 0", i));
  mn->Command(Form("FIX %d", i));
}


static bool onegoodminos(const int par, const bool no_low_ok)
{
  if((!no_low_ok && getminerrdn(par) == 0) || getminerrup(par) == 0){
    printf("Incomplete MINOS errors for free par %d\n", par);
    return false;
  }
  return true;
}

static double max(const double a, const double b)
{
  return a > b? a: b;
}

void fcn(int & np, double * gin, double & like, double *par, int flag)
{
  like = 0;
  const double flat   = par[0],
               Nmich  = par[1],
               Tmich  = par[2],
               Nneut  = par[3],
               Tneut  = par[4],
               Aneut  = par[5],
               NB12   = par[6],
               pileup = par[7];

  ee->SetParameters(par);

  /*printf("%9f %9f %9f %9f %9f %9f %9f %9f\n",
   par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7]);*/

  /*
  printf(".");  fflush(stdout);
  ee->Draw("same");
  c1->Update();
  c1->Modified();
  */

  for(int i = 1; i <= fithist->GetNbinsX() && i <= maxfitt+nnegbins; i++){
    const double model = ee->Eval(fithist->GetBinCenter(i));
    const double data = fithist->GetBinContent(i);
    
    like += model - data;
    if(model > 0 && data > 0) like += data * log(data/model);
  }
}

void make_mn()
{
  if(mn) delete mn;
  mn = new TMinuit(npar);
  mn->fGraphicsMode = false;
  mn->SetFCN(fcn);
  mn->Command("SET ERR 0.5"); // log likelihood
  mn->SetPrintLevel(-1);
  mn->Command("SET PRINT -1");

  for(int i = 0; i < npar; i++){
    int mnparmerr = 0;
    mn->mnparm(i, parnames[i], parinit[i], parinit[i]/100., 0., 0., mnparmerr);
    ee->SetParName(i, parnames[i]);
  }
}


static void setee_to_mn()
{
  for(int i = 0; i < npar; i++) ee->SetParameter(i, getpar(i));
}

static fitanswers dothefit(TH1D * hist, const bool is_rhc,
                           const double ntrack)
{
  // Start each fit with a clean MINUIT slate.  This is crucial!
  make_mn();

  fitanswers ans;
  fithist = hist;
  hist->GetXaxis()->SetTitle(is_rhc?"RHC":"FHC");
  hist->GetYaxis()->SetRangeUser(ntrack*1e-7, ntrack*1e-1);
  hist->Draw("e");
  int status = 0;

  // May as well set this to near the right value
  mn->Command(Form("SET PAR %d %f", flat_nf, 
    hist->Integral(0, nnegbins-10)/(nnegbins - 10)));

  // And this is a reasonable starting point for nmich
  mn->Command(Form("SET PAR %d %f", nmich_nf,
    ntrack*0.8));

  // Ditto nneut
  mn->Command(Form("SET PAR %d %f", nneut_nf,
    ntrack*0.1));

  // Again, something reasonable based ont the data.
  mn->Command(Form("SET PAR %d %f", pileup_nf,
    hist->GetBinContent(nnegbins - 2)/8.));


  // Fix neutron lifetime and diffusion parameters
  fixat(tneut_nf, 55);
  fixat(aneut_nf, 300);

  // Start with the muon lifetime fixed so that it doesn't try to swap
  // with the neutron lifetime.
  mn->Command(Form("FIX %d", tmich_nf));
  for(int i = 0; i < 8; i++)
    if(0 == (status = mn->Command("MIGRAD")))
      break;

  // Now that we're (hopefully) converged, let muon lifetime float
  mn->Command(Form("REL %d", tmich_nf));

  // Hold the Michel lifetime to something reasonable. Among other
  // concerns, this prevents it from swapping with the neutron lifetime,
  // supposing we let that float.
  mn->Command(Form("SET LIM %d 1.6 2.6", tmich_nf));

  // But check if there is almost no B-12, if so, fix it at the best fit
  // because otherwise we'll probably get problems running MINOS.
  if(fabs(getpar(nb12_nc)) < fabs(getpar(flat_nc))/10.)
    mn->Command(Form("FIX %d", nb12_nf));

  for(int i = 0; i < 8; i++)
    if(0 == (status = mn->Command("MIGRAD")))
      break;

  static bool first = true; 
  setee_to_mn();
  ee->Draw("same");
  c1->Print(Form("fit.pdf%s", first?"(":""));
  c1->Update(); c1->Modified();
  first = false;

  // MINOS errors are wrong if we used the abs trick, so limit
  mn->Command(Form("SET PAR %d %f",  nb12_nf, getpar( nb12_nc)));
  mn->Command(Form("SET PAR %d %f", nneut_nf, getpar(nneut_nc)));
  mn->Command(Form("SET LIM %d 0 %f",  nb12_nf, max(1e6, 10*getpar( nb12_nc))));
  mn->Command(Form("SET LIM %d 0 %f", nneut_nf, max(1e6, 10*getpar(nneut_nc))));

  for(int i = 0; i < 8; i++)
    if(0 == (status = mn->Command("MIGRAD")))
      break;

  if(!status)
    for(int i = 0; i < 2; i++)
      gMinuit->Command(Form("minos 30000 %d", nneut_nf));
  gMinuit->Command("show min");
  ans.n_good = onegoodminos(nneut_nc, false);

  ans.n_mag = getpar(nneut_nc);
  ans.n_mage_up = getminerrup(nneut_nc);
  ans.n_mage_dn = getminerrdn(nneut_nc);

  // release neutron lifetime, which we had held constant for a fair
  // comparison of neutrons, but is a serious nusiance parameter for
  // B-12. The better way to handle this would be a simulatanous fit
  // of all the histograms.
  mn->Command(Form("REL %d", tneut_nf));
  mn->Command(Form("SET LIM %d 40 70", tneut_nf));

  mn->Command(Form("REL %d", nb12_nf));

  for(int i = 0; i < 8; i++)
    if(0 == (status = mn->Command("MIGRAD")))
      break;
  if(!status)
    for(int i = 0; i < 2; i++)
      gMinuit->Command(Form("minos 30000 %d", nb12_nf));
  gMinuit->Command("show min");
  ans.b12_good = onegoodminos(nb12_nc, true);

  setee_to_mn();
  ee->Draw("same");
  c1->Print(Form("fit.pdf%s", first?"(":""));
  c1->Update(); c1->Modified();

  /* Cheating for sensitivity study! */
#if 0

  // Assume that we look at cosmic trigger data or some such to get
  // the noise level to high precision.
  mn->Command(Form("FIX %d", flat_nf));

  // Put in how many B-12 there ought to be
  ee->SetParameter(6, ntrack * (is_rhc? 0.15: 0.97) * 0.82 * 0.077 * 0.177 * 0.5);
  printf("Generating fake data with B-12 = %f\n", ee->GetParameter(6));

  for(int i = nnegbins + maxrealtime + 1; i <= hist->GetNbinsX(); i++)
    hist->SetBinContent(i, gRandom->Poisson(ee->Eval(hist->GetBinCenter(i))));

  maxfitt = maxrealtime + additional;

  for(int i = 0; i < 8; i++)
    if(0 == (status = mn->Command("MIGRAD")))
      break;
  if(!status)
    for(int i = 0; i < 2; i++)
      gMinuit->Command(Form("minos 30000 %d", nb12_nf));
  gMinuit->Command("show min");
  ans.b12_good = onegoodminos(nb12_nc, true);

  setee_to_mn();
  ee->Draw("same");
  c1->Print(Form("fit.pdf%s", first?"(":""));
  c1->Update(); c1->Modified();
#endif
  /* End cheating for sensitivity study! */

  ans.b12mag = getpar(nb12_nc);
  ans.b12mage_up = getminerrup(nb12_nc);
  ans.b12mage_dn = (getminerrdn(nb12_nc) != 0 && 
                 getminerrdn(nb12_nc) != 54321.0)?
                 getminerrdn(nb12_nc): getpar(nb12_nc);
  printf("b12mag = %f + %f - %f\n", ans.b12mag, ans.b12mage_up, ans.b12mage_dn);
  printf("n_mag  = %f + %f - %f\n",  ans.n_mag,  ans.n_mage_up,  ans.n_mage_dn);

  return ans;
}

static double ratio_error(const double x, const double y,
                          const double xe, const double ye)
{
  return 1/y * sqrt(pow(xe, 2) + x*x*pow(ye/y, 2));
}

void mncommand()
{
  std::string command;
  while(true){
    printf("MINUIT> ");
    if(!getline(cin, command)) break;
    if(command == "exit") break;
    mn->Command(command.c_str());
  }
}

void rhc()
{
  TFile * fhcfile = new TFile(
  "/nova/ana/users/mstrait/ndcosmic/period235-type3.root", "Read");
  //"/nova/ana/users/mstrait/ndcosmic/prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period2_v1/all-type3.root", "Read");
  //"/nova/ana/users/mstrait/ndcosmic/prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period3_v1/all.root", "Read");
  //"/nova/ana/users/mstrait/ndcosmic/prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period5_v1_goodruns/all-type3.root", "Read");
  TFile * rhcfile = new TFile("/nova/ana/users/mstrait/ndcosmic/prod_pid_S16-12-07_nd_period6_keepup/770-type3.root"                          , "Read");

  TTree * fhc_tree = (TTree *)fhcfile->Get("t");
  TTree * rhc_tree = (TTree *)rhcfile->Get("t");

  if(!fhc_tree || !rhc_tree || fhc_tree->IsZombie() || rhc_tree->IsZombie()){
    fprintf(stderr, "Couldn't read something.  See above.\n");
    return;
  }

  // Let the TFile errors go to the screen, then suppress the rest
  gErrorIgnoreLevel = kError;

  ee->SetNpx(460);
  ee->SetLineColor(kRed);
  ee->SetLineWidth(1);

  TH2D * dum = new TH2D("dum", "", 100, 0, bins_e[nbins_e], 1000, 0, 5);
  TH2D * dum2 = (TH2D*) dum->Clone("dum2");
  TH2D * dum3 = (TH2D*) dum->Clone("dum3");

  const char * const basecut = Form(
    "run != 11601" // noise at t = 106 in this run.
                   // I have no idea how that's possible.
    "&& primary && type == 3 && timeleft > %f && timeback > %f "
    "&& remid > 0.75 "
    "&& trklen > 200 " // standard
    //"&& trklen > 600 " // longer -- drops nearly all NC background
    // standard numu ND containment cuts, except for some muon-catcher
    // track checks which are irrelevant because I'm about to cut
    // those in the next line anyhow. I need the track ends to be more
    // contained than the standard cuts.
    "&& contained"

    // Sufficient to catch all neutrons within 6 cell widths. Maybe not
    // conservative enough, since neutrons that spill out into the air
    // probably don't ever come back? Or do they?
    "&& abs(trkx) < 170 && abs(trky) < 170 && trkz < 1250"
    "&& nclu < 12" // cut very noisy spills
    //"&& nslc <= 10" // reduce pileup
    , maxrealtime, -nnegbins);

  // Attempt to agressively reduce neutrons while still getting B-12
  /*
  const std::string ccut = Form("%s"
    "&& t > %f && t < %f"
    "&& !(t >= -1 && t < 2)"
    "&& nhitx == 1 && nhity == 1 && mindist < 1.99 && dist2 < 4"
    "&& pe > 70 && e < 20", -nnegbins, maxrealtime,
    basecut);
  */

  // a pretty strict, reasonable cut 
  /*
  const std::string ccut = Form("%s"
     "&& t > %f && t < %f"
     "&& !(t >= -1 && t < 2)"
     "&& nhitx >= 1 && nhity >= 1 && mindist <= 6"
     "&& pe > 70 && e < 20", basecut, -nnegbins, maxrealtime);
  */

  // a fairly loose, reasonable cut 
  const std::string ccut = Form("%s"
     "&& t > %f && t < %f"
     "&& !(t >= -1 && t < 2)"
     "&& nhit >= 2 && mindist <= 6"
     "&& pe > 70 && e < 20", basecut, -nnegbins, maxrealtime);

  // a loose cut
  /*
  const std::string ccut = Form("%s"
     "&& t > %f && t < %f"
     "&& !(t >= -1 && t < 2)", basecut, -nnegbins, maxrealtime);
  */

  fhc_s = new TH2D("fhc_s", "",
  nnegbins + maxrealtime + additional, -nnegbins, maxrealtime + additional,
  nbins_e, bins_e);
  rhc_s = (TH2D *)fhc_s->Clone("rhc_s");

  rhc_tree->Draw("slce:t >> rhc_s", ccut.c_str());
  fhc_tree->Draw("slce:t >> fhc_s", ccut.c_str());

  c1->SetLogy();

  TGraphAsymmErrors * g_n_rhc = new TGraphAsymmErrors;
  TGraphAsymmErrors * g_n_fhc = new TGraphAsymmErrors;
  TGraphAsymmErrors * g_b12_rhc = new TGraphAsymmErrors;
  TGraphAsymmErrors * g_b12_fhc = new TGraphAsymmErrors;
  g_n_rhc->SetMarkerStyle(kOpenSquare);
  g_n_fhc->SetMarkerStyle(kOpenCircle);

  g_b12_rhc->SetMarkerStyle(kOpenSquare);
  g_b12_fhc->SetMarkerStyle(kOpenCircle);

  g_b12_rhc->SetLineStyle(kDashed);
  g_b12_fhc->SetLineStyle(kDashed);

  g_n_rhc->SetLineWidth(1);
  g_n_fhc->SetLineWidth(2);

  g_b12_rhc->SetLineWidth(1);
  g_b12_fhc->SetLineWidth(2);

  TGraphAsymmErrors * g_n_rhc_bad = new TGraphAsymmErrors;
  TGraphAsymmErrors * g_n_fhc_bad = new TGraphAsymmErrors;
  TGraphAsymmErrors * g_b12_rhc_bad = new TGraphAsymmErrors;
  TGraphAsymmErrors * g_b12_fhc_bad = new TGraphAsymmErrors;

  g_n_rhc_bad->SetMarkerStyle(kOpenSquare);
  g_n_fhc_bad->SetMarkerStyle(kOpenCircle);

  g_b12_rhc_bad->SetMarkerStyle(kOpenSquare);
  g_b12_fhc_bad->SetMarkerStyle(kOpenCircle);

  g_b12_rhc_bad->SetLineStyle(kDashed);
  g_b12_fhc_bad->SetLineStyle(kDashed);

  g_n_rhc_bad->SetLineWidth(1);
  g_n_fhc_bad->SetLineWidth(2);

  g_b12_rhc_bad->SetLineWidth(1);
  g_b12_fhc_bad->SetLineWidth(2);

  g_n_rhc_bad->SetLineColor(kRed);
  g_n_fhc_bad->SetLineColor(kRed);
  g_b12_rhc_bad->SetLineColor(kRed);
  g_b12_fhc_bad->SetLineColor(kRed);

  g_n_rhc_bad->SetMarkerColor(kRed);
  g_n_fhc_bad->SetMarkerColor(kRed);
  g_b12_rhc_bad->SetMarkerColor(kRed);
  g_b12_fhc_bad->SetMarkerColor(kRed);



  TGraphAsymmErrors * n_result = new TGraphAsymmErrors;
  TGraphAsymmErrors * n_resultbad = new TGraphAsymmErrors;
  n_result->SetMarkerStyle(kOpenCircle);
  n_resultbad->SetMarkerStyle(kOpenCircle);
  n_resultbad->SetLineColor(kRed);
  n_resultbad->SetMarkerColor(kRed);
  n_result->SetMarkerSize(0.7);
  n_resultbad->SetMarkerSize(0.7);

  TGraphAsymmErrors * b12_result = new TGraphAsymmErrors;
  TGraphAsymmErrors * b12_resultbad = new TGraphAsymmErrors;
  b12_result->SetLineStyle(kDashed);
  b12_resultbad->SetLineStyle(kDashed);
  b12_result->SetMarkerStyle(kOpenCircle);
  b12_resultbad->SetMarkerStyle(kOpenCircle);
  b12_resultbad->SetLineColor(kRed);
  b12_resultbad->SetMarkerColor(kRed);
  b12_result->SetName("b12_result");
  b12_resultbad->SetName("b12_resultbad");
  b12_result->SetMarkerSize(0.7);
  b12_resultbad->SetMarkerSize(0.7);

  TH1D * tcounts_fhc = new TH1D("tcounts_fhc", "",
    fhc_s->GetNbinsY(), fhc_s->GetYaxis()->GetBinLowEdge(1),
    fhc_s->GetYaxis()->GetBinLowEdge(fhc_s->GetNbinsY()+1));
  TH1D * tcounts_rhc = (TH1D *)tcounts_fhc->Clone("tcounts_rhc");

  const std::string tcut = Form("i == 0 && %s", basecut);

  fhc_tree->Draw("slce >> tcounts_fhc", tcut.c_str());
  rhc_tree->Draw("slce >> tcounts_rhc", tcut.c_str());

  for(int s = 1; s <= nbins_e; s++){
    const double loslce = fhc_s->GetYaxis()->GetBinLowEdge(s);
    const double hislce = fhc_s->GetYaxis()->GetBinLowEdge(s+1);
  
    TH1D * fhc = fhc_s->ProjectionX("fhc", s, s);
    TH1D * rhc = rhc_s->ProjectionX("rhc", s, s);

    const double rhc_scale = tcounts_rhc->GetBinContent(s);
    const double fhc_scale = tcounts_fhc->GetBinContent(s);
    const double scale = fhc_scale/rhc_scale;
    printf("Number of tracks %f/%f = %f\n", fhc_scale, rhc_scale, scale);

    if(fhc->Integral() < 10 || rhc->Integral() < 10){
      printf("Only %f, %f events, skipping\n", rhc->Integral(), fhc->Integral());
      continue;
    }

    const fitanswers rhc_ans = dothefit(rhc, true, rhc_scale);
    const fitanswers fhc_ans = dothefit(fhc, false, fhc_scale);

    TGraphAsymmErrors * g_n_r = rhc_ans.n_good? g_n_rhc: g_n_rhc_bad;
    TGraphAsymmErrors * g_n_f = fhc_ans.n_good? g_n_fhc: g_n_fhc_bad;
    TGraphAsymmErrors * g_b12_r = rhc_ans.n_good? g_b12_rhc: g_b12_rhc_bad;
    TGraphAsymmErrors * g_b12_f = fhc_ans.n_good? g_b12_fhc: g_b12_fhc_bad;

    const double graph_x = (loslce+hislce)/2;
    const double graph_xe = (hislce-loslce)/2;

    g_n_r->SetPoint(g_n_r->GetN(), graph_x + graph_xe/10 /* visual offset */,
      rhc_ans.n_mag/rhc_scale);
    g_n_r->SetPointError(g_n_r->GetN()-1, graph_xe, graph_xe,
      rhc_ans.n_mage_dn/rhc_scale, rhc_ans.n_mage_up/rhc_scale);

    g_n_f->SetPoint(g_n_f->GetN(), graph_x, fhc_ans.n_mag/fhc_scale);
    g_n_f->SetPointError(g_n_f->GetN()-1, graph_xe, graph_xe,
      fhc_ans.n_mage_dn/fhc_scale, fhc_ans.n_mage_up/fhc_scale);

    g_b12_r->SetPoint(g_b12_r->GetN(), graph_x + graph_xe/10,
      rhc_ans.b12mag/rhc_scale);
    g_b12_r->SetPointError(g_b12_r->GetN()-1, graph_xe, graph_xe,
      rhc_ans.b12mage_dn/rhc_scale, rhc_ans.b12mage_up/rhc_scale);

    g_b12_f->SetPoint(g_b12_f->GetN(), graph_x, fhc_ans.b12mag/fhc_scale);
    g_b12_f->SetPointError(g_b12_f->GetN()-1, graph_xe, graph_xe,
      fhc_ans.b12mage_dn/fhc_scale, fhc_ans.b12mage_up/fhc_scale);

    const bool n_good = fhc_ans.n_good && rhc_ans.n_good;
    const bool b12_good = fhc_ans.b12_good && rhc_ans.b12_good;

    const double n_rat     = scale*rhc_ans.n_mag/fhc_ans.n_mag;
    const double n_rat_err_up =
      ratio_error(scale*rhc_ans.n_mag, fhc_ans.n_mag,
                  scale*rhc_ans.n_mage_up, fhc_ans.n_mage_dn);
    const double n_rat_err_dn =
      ratio_error(scale*rhc_ans.n_mag, fhc_ans.n_mag,
                  scale*rhc_ans.n_mage_dn, fhc_ans.n_mage_up);

    const double b12_rat     = scale*rhc_ans.b12mag/fhc_ans.b12mag;
    const double b12_rat_err_up =
      ratio_error(scale*rhc_ans.b12mag, fhc_ans.b12mag,
                  scale*rhc_ans.b12mage_up, fhc_ans.b12mage_dn);
    const double b12_rat_err_dn =
      ratio_error(scale*rhc_ans.b12mag, fhc_ans.b12mag,
                  scale*rhc_ans.b12mage_dn, fhc_ans.b12mage_up);

    printf("%s/%s RHC/FHC neutron (%4.2f-%4.2f)GeV: %.3f + %.3f - %.3f\n",
      rhc_ans.n_good?"Good":"Bad", fhc_ans.n_good?"Good":"Bad", loslce, hislce,
      n_rat, n_rat_err_up, n_rat_err_dn);
    TGraphAsymmErrors * addto = n_good?n_result:n_resultbad;
    addto->SetPoint(addto->GetN(), graph_x, n_rat);
    addto->SetPointError(addto->GetN()-1, graph_xe, graph_xe,
      n_rat_err_dn, n_rat_err_up);

    printf("%s/%s RHC/FHC B-12    (%4.2f-%4.2f)GeV: %.3f + %.3f - %.3f\n",
      rhc_ans.b12_good?"Good":"Bad", fhc_ans.b12_good?"Good":"Bad", loslce, hislce,
      b12_rat, b12_rat_err_up, b12_rat_err_dn);
    addto = b12_good?b12_result:b12_resultbad;
    addto->SetPoint(addto->GetN(), graph_x + graph_xe/10, b12_rat);
    addto->SetPointError(addto->GetN()-1, graph_xe, graph_xe,
      b12_rat_err_dn, b12_rat_err_up);
  }

  TCanvas * c2 = new TCanvas;
  dum2->GetYaxis()->SetTitle("Neutrons per track");
  dum2->GetYaxis()->SetRangeUser(0, 1.5);
  dum2->Draw();
  g_n_rhc->Draw("pz");
  g_n_rhc_bad->Draw("pz");
  g_n_fhc->Draw("pz");
  g_n_rhc_bad->Draw("pz");
  c2->Print("fit.pdf");
  
  TCanvas * c3 = new TCanvas;
  dum3->GetYaxis()->SetTitle("B-12 per track");
  dum3->Draw();
  g_b12_rhc->Draw("pz");
  g_b12_rhc_bad->Draw("pz");
  g_b12_fhc->Draw("pz");
  g_b12_rhc_bad->Draw("pz");
  c3->Print("fit.pdf");

  c1->cd();
  c1->SetLogy(0);

  dum->GetYaxis()->SetTitle("Ratios");
  dum->Draw();
  n_result->Draw("pz");
  n_resultbad->Draw("pz");
  b12_result->Draw("pz");
  b12_resultbad->Draw("pz");

  c1->Print("fit.pdf)");
}
