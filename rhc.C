#include "TMinuit.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TError.h"
#include "TRandom.h"
#include "TROOT.h"
#include <fstream>

/*
 * TODO: exclude neutron and B-12 candidates following a Michel.
 *
 * TODO: joint fit that directly extracts the ratios of interest
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

const double n_lifetime_nominal = 50.;
const double n_lifetime_priorerr = 5.;

/*
   Based on my MC, 400us is reasonable for a mindist of 6, but more like
   25us for a mindist of 2. However, it's murky because we don't measure
   the capture point, but rather the position where the gammas compton,
   and that usually only in 2D. Assuming 2D, I get more like 80us. And
   of course the functional form isn't quite right for 2D...
*/
const double n_diffusion_nominal = 400.;
const double n_diffusion_priorerr = n_diffusion_nominal*0.5;

const double markersize = 0.7;

struct PAR {
  char * name;
  int nc;
  int nf;
  double start;
};

static PAR makepar(const char * const name_, const int n, const double start_)
{
  PAR p;

  if(strlen(name_) > 30){
    fprintf(stderr, "Your parameter name is out of control\n");
    exit(1);
  }

  p.name = (char *)malloc(strlen(name_)+1);
  strcpy(p.name, name_);
  
  p.nc = n;
  p.nf = n+1;

  p.start = start_;

  return p;
}

const int ncommonpar = 2;
const int nperbinpar = 6;

const int npar    = ncommonpar + nperbinpar*nbins_e;
const int npar_ee = ncommonpar + nperbinpar + 1;
const char * const ee_parnames[npar_ee] = {
"Tneut",
"Aneut",

"Nneut",
"NB12",

"Tmich",
"flat",
"NMich",
"pileup",

"scale"
};

const int ee_scale_par = npar_ee-1;

const double tmich_start = 2.1;
const double nneut_start = 1.5e5;
const double nb12_start = 2e4;
const double flat_start = 3;
const double nmich_start = 1.5e5;
const double pileup_start = 300;

// parameter numbers, C numbering, for use with TMinuit functions
const int
  // Same for all histograms
  tneut_nc  = 0,
  aneut_nc  = 1,

  // Parameters of interest, energy dependent
  nneut_nc  = aneut_nc + 1,
  nb12_nc   = nneut_nc + nbins_e,

  // Nuisance parameters, energy dependent
  tmich_nc  = nb12_nc  + nbins_e,
  flat_nc   = tmich_nc + nbins_e,
  nmich_nc  = flat_nc  + nbins_e,
  pileup_nc = nmich_nc + nbins_e;

const int flat_nf   = flat_nc  +1, // and FORTRAN numbering,
          nmich_nf  = nmich_nc +1, // for use with native MINUIT
          tmich_nf  = tmich_nc +1, // commands through TMinuit::Command()
          nneut_nf  = nneut_nc +1,
          tneut_nf  = tneut_nc +1,
          aneut_nf  = aneut_nc +1,
          nb12_nf   = nb12_nc  +1,
          pileup_nf = pileup_nc+1;

static std::vector<PAR> makeparameters()
{
  std::vector<PAR> p;
  p.push_back(makepar("Tneut", tneut_nc, n_lifetime_nominal));
  p.push_back(makepar("Aneut", aneut_nc, n_diffusion_nominal));

  for(int i = 0; i < nbins_e; i++)
    p.push_back(makepar(Form("Nneut%d", i), nneut_nc+i, nneut_start));
  for(int i = 0; i < nbins_e; i++)
    p.push_back(makepar(Form("NB12_%d", i), nb12_nc+i, nb12_start));

  for(int i = 0; i < nbins_e; i++)
    p.push_back(makepar(Form("TMich%d", i), tmich_nc, tmich_start));
  for(int i = 0; i < nbins_e; i++)
    p.push_back(makepar(Form("flat%d", i), flat_nc+i, flat_start));
  for(int i = 0; i < nbins_e; i++)
    p.push_back(makepar(Form("NMich%d", i), nmich_nc+i, nmich_start));
  for(int i = 0; i < nbins_e; i++)
    p.push_back(makepar(Form("pileup%d", i), pileup_nc+i, pileup_start));

  return p;
}

static std::vector<PAR> PARAMETERS = makeparameters();

// TF1::Draw() has a lot of trouble with the discontinuity at zero, 
// so split into the positive and negative parts
static TF1 * ee_pos = new TF1("ee_pos",
  "[8]*(abs([5]) + "
  "(x >= 0)*("
   "abs([6])/[4] * exp(-x/[4]) + "
   "abs([2])/[0] * exp(-x/[0]) "
   "*(TMath::Erf(sqrt([1]/x))-2/sqrt(TMath::Pi())*sqrt([1]/x)*exp(-[1]/x))"
   "+ abs([3])/29.14e3 * exp(-x/29.14e3)"
  ") + "
  "((x >= -10 && x <= 10))*(abs([7])*abs(abs(x)-10)))",
  0, maxrealtime+additional);

static TF1 * ee_neg = new TF1("ee_neg",
  "[8]*(abs([5]) + "
  "((x >= -10 && x <= 10))*(abs([7])*abs(abs(x)-10)))",
  -nnegbins, 0);

static TF1 * ee_flat =
  new TF1("ee_flat", "[8]*abs([5])", -nnegbins, maxrealtime+additional);

static TF1 * ee_mich =
  new TF1("ee_mich", "[8]*(abs([6])/[4] * exp(-x/[4]))", 0, maxrealtime+additional);

static TF1 * ee_neut =
  new TF1("ee_mich", 
   "[8]*(abs([2])/[0] * exp(-x/[0]) "
   "*(TMath::Erf(sqrt([1]/x))-2/sqrt(TMath::Pi())*sqrt([1]/x)*exp(-[1]/x)))"
   , 0, maxrealtime+additional);

static TF1 * ee_b12 =
  new TF1("ee_b12", "[8]*(abs([3])/29.14e3 * exp(-x/29.14e3))", 0, maxrealtime+additional);

static TF1 * ee_pileup =
  new TF1("ee_b12", "[8]*((x >= -10 && x <= 10))*(abs([7])*abs(abs(x)-10))",
         -nnegbins, maxrealtime+additional);

const int ntf1s = 7;
static TF1 * ees[ntf1s] =
  { ee_neg, ee_pos, ee_flat, ee_mich, ee_neut, ee_b12, ee_pileup };

static TH2D * fhc_s = NULL;
static TH2D * rhc_s = NULL;

static TH1D * tcounts_fhc = NULL;
static TH1D * tcounts_rhc = NULL;

static TH2D * fithist = NULL;

static TCanvas * c1 = new TCanvas("rhc1", "rhc1");
static TCanvas * c2 = new TCanvas("rhc2", "rhc2");
static TCanvas * c3 = new TCanvas("rhc3", "rhc3");
static TCanvas * c4 = new TCanvas("rhc4", "rhc4");

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

static double min(const double a, const double b)
{
  return a < b? a: b;
}

static double max(const double a, const double b)
{
  return a > b? a: b;
}

void fcn(int & np, double * gin, double & like, double *par, int flag)
{
  like = 0;

  const double b12life = 29.14e3;

  for(int eb = 0; eb < fithist->GetNbinsY(); eb++){
    for(int tb = 1; tb <= fithist->GetNbinsX() && tb <= maxfitt+nnegbins; tb++){
      const double x = fithist->GetXaxis()->GetBinCenter(tb);
      double model = 0;
      if(x <= -1 || x >= 2) model += abs(par[flat_nc+eb]);

      if(x >= 2)
        model += fabs(par[nmich_nc+eb])/par[tmich_nc] * exp(-x/par[tmich_nc]) +
                 fabs(par[nneut_nc+eb])/par[tneut_nc] * exp(-x/par[tneut_nc])
           *(erf(sqrt(par[aneut_nc]/x))
             -2/sqrt(M_PI)*sqrt(par[aneut_nc]/x)*exp(-par[aneut_nc]/x))
               + fabs(par[nb12_nc+eb])/b12life * exp(-x/b12life);

      if((x >= -10 && x <= -1) || (x >= 2 && x <= 10))
        model += fabs(par[pileup_nc+eb])*fabs(fabs(x)-10);

      const double data = fithist->GetBinContent(tb, eb+1 /* 0->1 index */);

      like += model - data;
      if(model > 0 && data > 0) like += data * log(data/model);
    }
  }

  const double tneut_penalty =
    0.5 * pow((par[tneut_nc] - n_lifetime_nominal)/n_lifetime_priorerr, 2);

  const double aneut_penalty =
    0.5 * pow((par[aneut_nc] - n_diffusion_nominal)/n_diffusion_priorerr, 2);

  like += tneut_penalty + aneut_penalty;
}

void make_mn()
{
  if(mn) delete mn;
  mn = new TMinuit(npar);
  mn->fGraphicsMode = false;
  mn->SetPrintLevel(-1);
  mn->Command("SET PRINT -1");
  mn->SetFCN(fcn);
  mn->Command("SET ERR 0.5"); // log likelihood

  int mnparmerr = 0;
  for(int i = 0; i < npar; i++)
    mn->mnparm(i, PARAMETERS[i].name,
                  PARAMETERS[i].start,
                  PARAMETERS[i].start/100., 0., 0., mnparmerr);
}

static void set_ee_to_mn(const int bin) // 0-indexed
{
  for(int i = 0; i < ncommonpar; i++){
    ee_neg->SetParameter(i, getpar(i));
  }

  for(int i = ncommonpar; i < npar_ee; i++){
    ee_neg->SetParameter(i, getpar(ncommonpar + (i-ncommonpar)*nbins_e + bin));
  }

  ee_neg->SetParameter(ee_scale_par, 1); // scale

  for(int f = 0; f < ntf1s; f++)
    for(int i = 0; i < npar_ee; i++)
      ees[f]->SetParameter(i, ee_neg->GetParameter(i));
}

// Have to get rid of any bins that have partial contents.
// Just for drawing.
static void sanitize_after_rebinning(TH1D * x)
{
  // if the bin contains the range [-1,2], zero it
  for(int b = 0; b <= x->GetNbinsX(); b++){
    if(-1.0 > x->GetBinLowEdge(b) &&
       -1.0 < x->GetBinLowEdge(b+1)) x->SetBinContent(b, 0);
    if(2 > x->GetBinLowEdge(b) &&
       2 < x->GetBinLowEdge(b+1)) x->SetBinContent(b, 0);
    // does not consider the case that the bin is entirely
    // inside of the range, because then it will be zero anyway
  }
}

static void draw_ee(const int bin, const double ntrack) // 0-indexed
{
  const int rebin = 2;

  c1->cd();
  static bool first = true;
  set_ee_to_mn(bin);
  TH1D * x = fithist->ProjectionX("x", bin+1, bin+1 /* 0->1 */);
  x->Rebin(rebin);

  sanitize_after_rebinning(x);

  if(rebin == 1)
    x->GetYaxis()->SetTitle("Delayed clusters/#mus");
  else
    x->GetYaxis()->SetTitle(Form("Delayed clusters/%d#mus", rebin));
  x->GetXaxis()->SetTitle("Time since muon stop (#mus)");
  x->Draw("e");

  x->GetYaxis()->SetRangeUser(min(0.5, ntrack*1e-5)*rebin, ntrack*rebin*0.1);

  for(int i = 0; i < ntf1s; i++){
    ees[i]->SetParameter(ee_scale_par, rebin);
    ees[i]->Draw("same");
  }

  c1->SetLogy();
  c1->Print(Form("fit.pdf%s", first?"(":""));
  c1->Update(); c1->Modified();
  first = false;
}

  /* Cheating for sensitivity study! */
#if 0
  // Assume that we look at cosmic trigger data or some such to get
  // the noise level to high precision.
  mn->Command(Form("FIX %d", flat_nf));

  set_ee_to_mn(0 /* XXX */);

  // Keep however many neutrons the fit above said.
  // And put in how many B-12 there ought to be.
  ee_pos->SetParameter(nb12_nc,
    ntrack * (is_rhc? 0.2: 0.9) // mu- fraction
           * 0.82  // atomic capture
           * 0.077 // nuclear capture
           * 0.177 // B-12 yield
           * 0.4); // efficiency
  printf("Generating fake data with B-12 = %f\n", ee_pos->GetParameter(nb12_nc));

  for(int i = nnegbins + maxrealtime + 1; i <= hist->GetNbinsX(); i++)
    hist->SetBinContent(i, gRandom->Poisson(ee->Eval(hist->GetBinCenter(i))));

  maxfitt = maxrealtime + additional;

  for(int i = 0; i < migrad_tries; i++)
    if(0 == (status = mn->Command("MIGRAD")))
      break;
  if(!status)
    for(int i = 0; i < 2; i++){
      gMinuit->Command(Form("MINOS 30000 %d", nb12_nf));
      gMinuit->Command(Form("MINOS 30000 %d", nneut_nf));
    }

  gMinuit->Command("show min");
  for(int i = 0; i < nbins_e; i++) draw_ee(i, ntrack[i]);
#endif
  /* End cheating for sensitivity study! */

static std::vector<fitanswers> dothefit(TH2D * hist, const bool is_rhc,
                                        const std::vector<double> & ntrack)
{
  maxfitt = maxrealtime; // modified later for sensitivity study

  // Start each fit with a clean MINUIT slate.  This is crucial!
  make_mn();

  fithist = hist;

  int status = 0;

  // May as well set this to near the right value
  mn->Command(Form("SET PAR %d %f", flat_nf,
    hist->Integral(0, nnegbins-10)/(nnegbins - 10)));

  for(int bin = 0; bin < nbins_e; bin++){
    // And this is a reasonable starting point for nmich
    mn->Command(Form("SET PAR %d %f", nmich_nf+bin, ntrack[bin]*0.8));

    // Ditto nneut
    mn->Command(Form("SET PAR %d %f", nneut_nf+bin, ntrack[bin]*0.1));

    // Again, something reasonable based ont the data.
    mn->Command(Form("SET PAR %d %f", pileup_nf+bin,
      hist->GetBinContent(nnegbins - 2, bin+1 /* 0->1 */)/8.));

    // Start with the muon lifetime fixed so that it doesn't try to swap
    // with the neutron lifetime.
    mn->Command(Form("FIX %d", tmich_nf+bin));
  }

  // Constraining neutron lifetime with a pull term, but also
  // set hard limits to hold it to something reasonable while
  // the rest of the fit settles.
  mn->Command(Form("SET LIM %d 30 80", tneut_nf));

  // Letting the diffusion parameter float makes the fit not converge,
  // so don't do that.
  mn->Command(Form("FIX %d", aneut_nf));

  const int migrad_tries = 8;

  for(int i = 0; i < migrad_tries; i++)
    if(0 == (status = mn->Command("MIGRAD")))
      break;
  status = mn->Command("HESSE");

  // Now that we're (hopefully) converged, let muon lifetime float
  for(int bin = 0; bin < nbins_e; bin++)
    mn->Command(Form("REL %d", tmich_nf+bin));

  // And the diffusion parameter
  mn->Command(Form("REL %d", aneut_nf));

  // This becomes badly behaved if allowed to wander too far, so limit
  mn->Command(Form("SET LIM %d %f %f", aneut_nf,
    n_diffusion_nominal*0.5, n_diffusion_nominal*1.5));

  // Hold the Michel lifetime to something reasonable. Among other
  // concerns, this prevents it from swapping with the neutron lifetime,
  // supposing we let that float.
  for(int bin = 0; bin < nbins_e; bin++)
    mn->Command(Form("SET LIM %d 1.6 2.6", tmich_nf));

  // MINOS errors are wrong if we used the abs trick, so limit
  for(int bin = 0; bin < nbins_e; bin++){
    mn->Command(Form("SET PAR %d %f",    nb12_nf+bin, getpar( nb12_nc+bin)));
    mn->Command(Form("SET PAR %d %f",   nneut_nf+bin, getpar(nneut_nc+bin)));
    mn->Command(Form("SET LIM %d 0 %f", 
      nb12_nf+bin, max(1e6, 10*getpar( nb12_nc+bin))));
    mn->Command(Form("SET LIM %d 0 %f",
      nneut_nf+bin, max(1e6, 10*getpar(nneut_nc+bin))));
  }

  for(int i = 0; i < migrad_tries; i++)
    if(0 == (status = mn->Command("MIGRAD")))
      break;
  status = mn->Command("HESSE");

  if(!status)
    for(int bin = 0; bin < nbins_e; bin++){
      for(int i = 0; i < 2; i++){
        gMinuit->Command(Form("MINOS 30000 %d", nneut_nf+bin));
        gMinuit->Command(Form("MINOS 30000 %d", nb12_nf +bin));
      }
    }

  gMinuit->Command("show min");
  for(int i = 0; i < nbins_e; i++) draw_ee(i, ntrack[i]);


  std::vector<fitanswers> anses;

  for(int bin = 0; bin < nbins_e; bin++){
    fitanswers ans;
    ans.n_good =   onegoodminos(nneut_nc+bin, false);
    ans.n_mag =          getpar(nneut_nc+bin);
    ans.n_mage_up = getminerrup(nneut_nc+bin);
    ans.n_mage_dn = getminerrdn(nneut_nc+bin);

    ans.b12_good =   onegoodminos(nb12_nc+bin, true);
    ans.b12mag =           getpar(nb12_nc+bin);
    ans.b12mage_up =  getminerrup(nb12_nc+bin);
    ans.b12mage_dn =
      (getminerrdn(nb12_nc+bin) != 0 && getminerrdn(nb12_nc+bin) != 54321.0)?
       getminerrdn(nb12_nc+bin): getpar(nb12_nc+bin);

    printf("b12mag = %f + %f - %f\n", ans.b12mag, ans.b12mage_up, ans.b12mage_dn);
    printf("n_mag  = %f + %f - %f\n",  ans.n_mag,  ans.n_mage_up,  ans.n_mage_dn);

    anses.push_back(ans);
  }

  return anses;
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
    if(!getline(std::cin, command)) break;
    if(command == "exit") break;
    mn->Command(command.c_str());
  }
}

static TGraphAsymmErrors *
newgraph(const int color, const int linestyle, const int marker,
         const int linewidth)
{
  TGraphAsymmErrors * g = new TGraphAsymmErrors;
  g->SetMarkerSize(markersize);

  g->SetMarkerStyle(marker);
  g->SetLineStyle(linestyle);
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetLineWidth(linewidth);
  return g;
}

static void addpoint(TGraphAsymmErrors * g, const double x,
                     const double y, const double xe,
                     const double ye_down, const double ye_up)
{
  g->SetPoint(g->GetN(), x, y);
  g->SetPointError(g->GetN()-1, xe, xe, ye_down, ye_up);
}

static void init_ee()
{
  for(int i = 0; i < ntf1s; i++){
    TF1 * e = ees[i];
    switch(i){
      case 0: case 1:
        e->SetLineColor(kRed);
        e->SetLineWidth(2);
        e->SetNpx(400);
        break;
      case 2:
        e->SetLineStyle(kDashed);
        e->SetLineColor(kBlack);
        e->SetLineWidth(1);
        break;
      case 3:
        e->SetLineStyle(kDashed);
        e->SetLineColor(kBlue);
        e->SetLineWidth(1);
        break;
      case 4:
        e->SetLineStyle(kDashed);
        e->SetLineColor(kViolet);
        e->SetLineWidth(2);
        break;
      case 5:
        e->SetLineStyle(kDashed);
        e->SetLineColor(kGreen+2);
        e->SetLineWidth(2);
        break;
      case 6:
        e->SetLineStyle(kDashed);
        e->SetNpx(400);
        e->SetLineColor(kOrange+2);
        e->SetLineWidth(1);
        break;
    }
    for(int p = 0; p < npar_ee; p++) e->SetParName(p, ee_parnames[p]);
  }
}

void rhc(const char * const savedhistfile = NULL)
{
  TFile * fhcfile = new TFile(
  "/nova/ana/users/mstrait/ndcosmic/period235-type3.root", "Read");
  //"/nova/ana/users/mstrait/ndcosmic/prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period2_v1/all-type3.root", "Read");
  //"/nova/ana/users/mstrait/ndcosmic/prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period3_v1/all-type3.root", "Read");
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

  init_ee();

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
    //"&& nslc <= 10" // reduce pileup
    , maxrealtime, -nnegbins);

  // Attempt to agressively reduce neutrons while still getting B-12
  /*
  const std::string ccut = Form("%s"
    "&& t > %f && t < %f"
    "&& !(t >= -1 && t < 2)"
    "&& nhitx >= 1 && nhity >= 1" // maximum range is about 5.7cm
    "&& nhitx <= 2 && nhity <= 2"
    "&& nhit <= 3"
    "&& mindist < 2.8" // allow up to 1 plane and 2 cells off
    "&& pe > 35 && e > 8*0.62 && e < 25*0.62", -nnegbins, maxrealtime,
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
     "&& nhit >= 1 && mindist <= 0"
     "&& pe > 35 && e < 20", basecut, -nnegbins, maxrealtime);

  // a loose cut
  /*
  const std::string ccut = Form("%s"
     "&& t > %f && t < %f"
     "&& !(t >= -1 && t < 2)", basecut, -nnegbins, maxrealtime);
  */

  fhc_s = new TH2D("fhc_s", "", nnegbins + maxrealtime + additional,
                                -nnegbins, maxrealtime + additional,
                                nbins_e, bins_e);
  rhc_s = (TH2D *)fhc_s->Clone("rhc_s");

  // Square means RHC, solid means neutron, black means good fit
  TGraphAsymmErrors
    * g_n_rhc       = newgraph(kBlack, kSolid, kOpenSquare, 1),
    * g_n_fhc       = newgraph(kBlack, kSolid, kOpenCircle, 1),
    * g_b12_rhc     = newgraph(kBlack, kDashed,kOpenSquare, 1),
    * g_b12_fhc     = newgraph(kBlack, kDashed,kOpenCircle, 1),
    * g_n_rhc_bad   = newgraph(kRed,   kSolid, kOpenSquare, 1),
    * g_n_fhc_bad   = newgraph(kRed,   kSolid, kOpenCircle, 1),
    * g_b12_rhc_bad = newgraph(kRed,   kDashed,kOpenSquare, 1),
    * g_b12_fhc_bad = newgraph(kRed,   kDashed,kOpenCircle, 1),
    * n_result      = newgraph(kBlack, kSolid, kOpenCircle, 1),
    * n_resultbad   = newgraph(kRed,   kSolid, kOpenCircle, 1),
    * b12_result    = newgraph(kBlack, kDashed,kOpenCircle, 1),
    * b12_resultbad = newgraph(kRed,   kDashed,kOpenCircle, 1);

  b12_result->SetLineStyle(kDashed);
  tcounts_fhc = new TH1D("tcounts_fhc", "", nbins_e, bins_e);
  tcounts_rhc = (TH1D *)tcounts_fhc->Clone("tcounts_rhc");

  const std::string tcut = Form("i == 0 && %s", basecut);

  if(savedhistfile == NULL){
    rhc_tree->Draw("slce:t >> rhc_s", ccut.c_str());
    fhc_tree->Draw("slce:t >> fhc_s", ccut.c_str());

    fhc_tree->Draw("slce >> tcounts_fhc", tcut.c_str());
    rhc_tree->Draw("slce >> tcounts_rhc", tcut.c_str());
  }
  else{
    gROOT->Macro(savedhistfile);
  }

  std::vector<double> rhc_scales, fhc_scales;

  for(int s = 1; s <= nbins_e; s++){
    rhc_scales.push_back(tcounts_rhc->GetBinContent(s));
    fhc_scales.push_back(tcounts_fhc->GetBinContent(s));
  }

  const std::vector<fitanswers> rhc_anses = dothefit(rhc_s, true , rhc_scales);
  const std::vector<fitanswers> fhc_anses = dothefit(fhc_s, false, fhc_scales);

  for(int s = 0; s < nbins_e; s++){
    const double rhc_scale = rhc_scales[s];
    const double fhc_scale = fhc_scales[s];
    const fitanswers rhc_ans = rhc_anses[s];
    const fitanswers fhc_ans = fhc_anses[s];

    const double scale = fhc_scale/rhc_scale;
    printf("N tracks %.0f, %.0f. Scale = %f\n", fhc_scale, rhc_scale, scale);

    const double loslce = fhc_s->GetYaxis()->GetBinLowEdge(s+1);
    const double hislce = fhc_s->GetYaxis()->GetBinLowEdge(s+2);

    const double graph_x  = (loslce+hislce)/2;
    const double graph_xe = (hislce-loslce)/2;

    addpoint(rhc_ans.n_good? g_n_rhc: g_n_rhc_bad,
      graph_x + graph_xe/10 /* visual offset */,
      rhc_ans.n_mag/rhc_scale, graph_xe,
      rhc_ans.n_mage_dn/rhc_scale, rhc_ans.n_mage_up/rhc_scale);

    addpoint(fhc_ans.n_good? g_n_fhc: g_n_fhc_bad,
      graph_x, fhc_ans.n_mag/fhc_scale, graph_xe,
      fhc_ans.n_mage_dn/fhc_scale, fhc_ans.n_mage_up/fhc_scale);

    addpoint(rhc_ans.n_good? g_b12_rhc: g_b12_rhc_bad,
      graph_x + graph_xe/10, rhc_ans.b12mag/rhc_scale,
      graph_xe, rhc_ans.b12mage_dn/rhc_scale, rhc_ans.b12mage_up/rhc_scale);

    addpoint(fhc_ans.n_good? g_b12_fhc: g_b12_fhc_bad,
      graph_x, fhc_ans.b12mag/fhc_scale, graph_xe,
      fhc_ans.b12mage_dn/fhc_scale, fhc_ans.b12mage_up/fhc_scale);

    const bool n_good   = fhc_ans.  n_good && rhc_ans.  n_good;
    const bool b12_good = fhc_ans.b12_good && rhc_ans.b12_good;

    const double n_rat = scale*rhc_ans.n_mag/fhc_ans.n_mag;
    const double n_rat_err_up =
      ratio_error(scale*rhc_ans.n_mag, fhc_ans.n_mag,
                  scale*rhc_ans.n_mage_up, fhc_ans.n_mage_dn);
    const double n_rat_err_dn =
      ratio_error(scale*rhc_ans.n_mag, fhc_ans.n_mag,
                  scale*rhc_ans.n_mage_dn, fhc_ans.n_mage_up);

    const double b12_rat = scale*rhc_ans.b12mag/fhc_ans.b12mag;
    const double b12_rat_err_up =
      ratio_error(scale*rhc_ans.b12mag, fhc_ans.b12mag,
                  scale*rhc_ans.b12mage_up, fhc_ans.b12mage_dn);
    const double b12_rat_err_dn =
      ratio_error(scale*rhc_ans.b12mag, fhc_ans.b12mag,
                  scale*rhc_ans.b12mage_dn, fhc_ans.b12mage_up);

    printf("%s/%s RHC/FHC neutron (%4.2f-%4.2f)GeV: %.3f + %.3f - %.3f\n",
      rhc_ans.n_good?"Good":"Bad", fhc_ans.n_good?"Good":"Bad", loslce, hislce,
      n_rat, n_rat_err_up, n_rat_err_dn);

    addpoint(n_good?n_result:n_resultbad,
             graph_x, n_rat, graph_xe, n_rat_err_dn, n_rat_err_up);

    printf("%s/%s RHC/FHC B-12    (%4.2f-%4.2f)GeV: %.3f + %.3f - %.3f\n",
      rhc_ans.b12_good?"Good":"Bad", fhc_ans.b12_good?"Good":"Bad", loslce, hislce,
      b12_rat, b12_rat_err_up, b12_rat_err_dn);

    addpoint(b12_good?b12_result:b12_resultbad,
             graph_x + graph_xe/10, b12_rat,
             graph_xe, b12_rat_err_dn, b12_rat_err_up);
  }

  TH2D * dum = new TH2D("dum", "", 100, 0, bins_e[nbins_e], 10000, 0, 10);
  TH2D * dum2 = (TH2D*) dum->Clone("dum2");
  TH2D * dum3 = (TH2D*) dum->Clone("dum3");

  c2->cd();
  dum2->GetYaxis()->SetTitle("Neutrons per track");
  dum2->GetYaxis()->SetRangeUser(0, 0.29);
  dum2->Draw();
  g_n_rhc->Draw("pz");
  g_n_rhc_bad->Draw("pz");
  g_n_fhc->Draw("pz");
  g_n_rhc_bad->Draw("pz");
  c2->Print("fit.pdf");

  c3->cd();
  dum3->GetYaxis()->SetTitle("B-12 per track");
  c3->SetLogy();
  dum3->GetYaxis()->SetRangeUser(0.001, 6);
  dum3->Draw();
  g_b12_rhc->Draw("pz");
  g_b12_rhc_bad->Draw("pz");
  g_b12_fhc->Draw("pz");
  g_b12_rhc_bad->Draw("pz");
  c3->Print("fit.pdf");

  c4->cd();
  dum->GetYaxis()->SetTitle("Ratios");
  dum->GetYaxis()->SetRangeUser(0, 1.5);
  dum->Draw();
  n_result->Draw("pz");
  n_resultbad->Draw("pz");
  b12_result->Draw("pz");
  b12_resultbad->Draw("pz");

  c4->Print("fit.pdf)");

  if(savedhistfile == NULL){
    std::ofstream o("savedhists.C");

    rhc_s->SavePrimitive(o);
    fhc_s->SavePrimitive(o);
    tcounts_fhc->SavePrimitive(o);
    tcounts_rhc->SavePrimitive(o);
  }
}
