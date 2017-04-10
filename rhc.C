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
#include "TLegend.h"
#include <fstream>

/*
 * TODO: exclude neutron and B-12 candidates following a Michel.
 *
 * TODO: Directly extract the NC and numu fractions
 */

const double nnegbins = 209;
const double maxrealtime = 269;
const double additional = 0;

double maxfitt = maxrealtime; // varies

const char * const basecut = Form(
  "run != 11601" // noise at t = 106 in this run.
                 // Will go away with a better run list.
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
  "&& nslc <= 6" // reduce pileup
  , maxrealtime, -nnegbins);


const char * const clustercut = 
// Attempt to agressively reduce neutrons while still getting B-12
/*
  "nhitx >= 1 && nhity >= 1" // maximum range is about 5.7cm
  "&& nhitx <= 2 && nhity <= 2"
  "&& nhit <= 3"
  "&& mindist < 2.8" // allow up to 1 plane and 2 cells off
  "&& pe > 35 && e > 8*0.62 && e < 25*0.62";
*/

// a pretty strict, reasonable cut
/*
   "nhitx >= 1 && nhity >= 1 && mindist <= 6"
   "&& pe > 70 && e < 20";
*/

// a fairly loose, reasonable cut
  "nhit >= 1 && mindist <= 6"
  "&& pe > 35 && e < 20";

// a very loose cut
// "";

struct fitanswers{
  bool n_good, b12_good;
  double  n_mag,  n_mage_up,  n_mage_dn,
         b12mag, b12mage_up, b12mage_dn;
};

const int nbins_e = 6;
const double bins_e[nbins_e+1] = {0.5, 1, 1.5, 2.0, 2.5, 3.0, 6.0 };
//const double bins_e[nbins_e+1] = {0.5, 3, 6.0 };

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
const double n_diffusion_nominal = 200.;
const double n_diffusion_priorerr = n_diffusion_nominal*0.5;

// I am modeling the time of neutron captures as the convolution of
// two exponentials. This comes from two effects. (1) Neutrons are not
// produced until the muon captures, so if the muon lifetime and neutron
// capture time were simple exponential, you'd already get the sum of
// two exponentials. In fact, the captures that give neutrons are mostly
// on chlorine, but with lots of of titanium and carbon too, so it is
// the sum of several convolutions of exponentials. (2) The neutrons
// have to thermalize before capturing. It is roughly correct to model
// this as another convolution of two exponentials. Art gives 0.65us as
// the rising exponential. All of this put together looks close to just
// the convolution of two exponentials, and since it is all hidden under
// Michel decays, we have no power to fit for the details. This constant
// gives the effective approximate rising exponential. For just muons it
// is 1.6us.
//
// Note, though, that if half of the neutrons are from *pion* capture,
// this is more like 0.9us. If that number is used, you get a 1.4%
// increase in the number of neutrons inferred, i.e. for a fixed number
// here, the "efficiency" for pion-produced neutrons is about 1.4% lower
// than muon-produced mpared to muons.
//
// I'll try to get the shape right, and just note this efficiency
// difference.
static const double n_start_conv_time = (1.60 + 0.9)/2;

static const double markersize = 0.7;

struct PAR {
  char * name;
  double start;
};

static PAR makepar(const char * const name_, const double start_)
{
  PAR p;

  if(strlen(name_) > 30){
    fprintf(stderr, "Your parameter name %s of length %d is out of control\n",
           name_, (int)strlen(name_));
    exit(1);
  }

  p.name = (char *)malloc(strlen(name_)+1);
  strcpy(p.name, name_);

  p.start = start_;

  return p;
}

const int ncommonpar = 2;
const int nbeam      = 2; // not really generalizable as it stands
const int nperbinpar = 4; // number per energy bin parameters
const int nperperpar = 2; // number per period parameters

const int nperiodrhc = 2; // 4, 6
const int nperiodfhc = 4; // 1, 2, 3, 5
const int nperiod    = nperiodrhc + nperiodfhc;

const char * const Speriodnames[nperiod] =
    { "P6", "P4", "P1", "P2", "P3", "P5" };

const char * const inputfiles[nperiod] = {
"/nova/ana/users/mstrait/ndcosmic/prod_pid_"
  "S16-12-07_nd_period6_keepup/985-type3.root",
"/nova/ana/users/mstrait/ndcosmic/prod_pid_"
  "R16-12-20-prod3recopreview.b_nd_numi_rhc_epoch4a_v1_goodruns/all-type3.root",

"/nova/ana/users/mstrait/ndcosmic/prod_pid_"
  "R17-03-01-prod3reco.b_nd_numi_fhc_period1_v1_goodruns/all-type3.root",
"/nova/ana/users/mstrait/ndcosmic/prod_pid_"
  "R17-03-01-prod3reco.b_nd_numi_fhc_period2_v1_goodruns/all-type3.root",
"/nova/ana/users/mstrait/ndcosmic/prod_pid_"
  "R17-03-01-prod3reco.b_nd_numi_fhc_period3_v1_goodruns/all-type3.root",
"/nova/ana/users/mstrait/ndcosmic/prod_pid_"
  "R17-03-01-prod3reco.b_nd_numi_fhc_period5_v1_goodruns/all-type3.root"
};

const char * const Lperiodnames[nperiod] = {
     "Period 6 (RHC)",
     "Period 4 (RHC)",
     "Period 1 (FHC)",
     "Period 2 (FHC)",
     "Period 3 (FHC)",
     "Period 5 (FHC)" };

const int npar = ncommonpar
               + nbeam  *nperbinpar*nbins_e
               + nperiod*nperperpar*nbins_e;
const int npar_ee = ncommonpar + nperbinpar + nperperpar + 1;

// Woe betide you if you rearrange these, since the positions are
// hardcoded in various places.
const char * const ee_parnames[npar_ee] = {
"Tneut",
"Aneut",

"Nneut",
"NB12",

"Tmich",
"NMich",

"flat",
"pileup",

"scale"
};

const int ee_scale_par = npar_ee-1;

const double nmich_start = 0.8;
const double nneut_start = 0.05;
const double nb12_start  = 0.001;
const double tmich_start = 2.1;
const double flat_start = 3;
const double pileup_start = 300;

// parameter numbers, C numbering, for use with TMinuit functions
const int
  // Same for all histograms
  tneut_nc  = 0,
  aneut_nc  = 1,

  // Parameters of interest, energy dependent
  nneut_nc  = aneut_nc + 1,
  nb12_nc   = nneut_nc + nbins_e*nbeam,

  // Nuisance parameters, energy dependent
  tmich_nc  = nb12_nc  + nbins_e*nbeam,
  nmich_nc  = tmich_nc + nbins_e*nbeam,

  // Nuisance parameters, energy and period dependent
  flat_nc   = nmich_nc + nbins_e*nbeam,
  pileup_nc = flat_nc  + nbins_e*nperiod;

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
  p.push_back(makepar("Tneut", n_lifetime_nominal));
  p.push_back(makepar("Aneut", n_diffusion_nominal));

  const char * const bname[nbeam] = { "R_", "F_" };

  // first set means the number of neutrons or B-12 seen per track in
  // RHC. The second set means the ratio of that with the same for FHC.
  const char * const ratname[nbeam] = { "RoF_", "F_" };

  for(int b = 0; b < nbeam; b++)
    for(int i = 0; i < nbins_e; i++)
      p.push_back(makepar(Form("%sNneut%d", ratname[b], i), nneut_start));
  for(int b = 0; b < nbeam; b++)
    for(int i = 0; i < nbins_e; i++)
      p.push_back(makepar(Form("%sNB12_%d", ratname[b], i), nb12_start));

  for(int b = 0; b < nbeam; b++)
    for(int i = 0; i < nbins_e; i++)
      p.push_back(makepar(Form("%sTMich%d", bname[b], i), tmich_start));
  for(int b = 0; b < nbeam; b++)
    for(int i = 0; i < nbins_e; i++)
      p.push_back(makepar(Form("%sNMich%d", bname[b], i), nmich_start));

  for(int ip = 0; ip < nperiod; ip++)
    for(int i = 0; i < nbins_e; i++)
      p.push_back(makepar(Form("%sflat%d", Speriodnames[ip], i), flat_start));
  for(int ip = 0; ip < nperiod; ip++)
    for(int i = 0; i < nbins_e; i++)
      p.push_back(makepar(Form("%spileup%d", Speriodnames[ip], i), pileup_start));

  return p;
}

static std::vector<PAR> PARAMETERS = makeparameters();

static TF1 * ee_flat =
  new TF1("ee_flat", "[8]*abs([6])", -nnegbins, maxrealtime+additional);

static TF1 * ee_mich =
  new TF1("ee_mich", "[8]*(abs([5])/[4] * exp(-x/[4]))",
          2, maxrealtime+additional);

static TF1 * ee_neut =
  new TF1("ee_mich",
   Form(
   // Convolution of muon lifetime and neutron lifetime
   "[8]*(abs([2])/([0]-%f) * (exp(-x/[0]) - exp(-x/%f))"

   // Neutron diffusion for a spherical cut.  Approximate for an 
   // intersection-of-two-cylinders cut
   "*(TMath::Erf(sqrt([1]/x))-2/sqrt(TMath::Pi())*sqrt([1]/x)*exp(-[1]/x)))",
  n_start_conv_time, n_start_conv_time),
   2, maxrealtime+additional);

static TF1 * ee_b12 =
  new TF1("ee_b12", "[8]*(abs([3])/29.14e3 * exp(-x/29.14e3))",
          2, maxrealtime+additional);

static TF1 * ee_pileup = new TF1("ee_b12",
  "[8]*((x >= -10 && x <= 10))*(abs([7])*abs(abs(x)-10))",
  -nnegbins, maxrealtime+additional);

// The sum of all the above. TF1::Draw() has a lot of trouble with the
// discontinuity at zero, so split into the positive and negative parts.
static TF1 * ee_pos = new TF1("ee_pos",
  Form("[8]*(abs([6]) + "
  "(x >= 0)*("
    "abs([5])/[4] * exp(-x/[4]) + "
    "abs([2])/([0]-%f) * (exp(-x/[0]) - exp(-x/%f))"
    "*(TMath::Erf(sqrt([1]/x))-2/sqrt(TMath::Pi())*sqrt([1]/x)*exp(-[1]/x))"
    "+ abs([3])/29.14e3 * exp(-x/29.14e3)"
  ") + "
  "((x >= -10 && x <= 10))*(abs([7])*abs(abs(x)-10)))",
  n_start_conv_time, n_start_conv_time),
  2, maxrealtime+additional);

static TF1 * ee_neg = new TF1("ee_neg",
  "[8]*(abs([6]) + "
  "((x >= -10 && x <= 10))*(abs([7])*abs(abs(x)-10)))",
  -nnegbins, -1);

const int ntf1s = 7;
static TF1 * ees[ntf1s] =
  { ee_neg, ee_pos, ee_flat, ee_mich, ee_neut, ee_b12, ee_pileup };

static const char * const ees_description[ntf1s] =
  { "Full fit", "Full fit", "Uncorrelated", "Michels", "Neutrons",
    "^{12}B", "pileup" };

static std::vector< std::vector<double> > scales;

static TH2D ** fithist     = (TH2D**)malloc(nperiod*sizeof(TH2D*));
static TH1D ** all_tcounts = (TH1D**)malloc(nperiod*sizeof(TH1D*));

static TCanvas * c1 = new TCanvas("rhc1", "rhc1"); // hist
static TCanvas * c2 = new TCanvas("rhc2", "rhc2"); // neutrons
static TCanvas * c3 = new TCanvas("rhc3", "rhc3"); // B-12
static TCanvas * c4 = new TCanvas("rhc4", "rhc4"); // ratios

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
  return -mn->fErn[wherearewe]; // define all errors as positive
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

__attribute__((unused)) static double getlimlo(int i) // 0-indexed!
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

__attribute__((unused)) static void fixat(int i, float v) // 1-indexed!
{
  mn->Command(Form("REL %d", i));
  if(getlimup(i-1)) mn->Command(Form("SET LIM %d", i));
  mn->Command(Form("SET PAR %d %g", i, v));
  mn->Command(Form("FIX %d", i));
}

__attribute__((unused)) static void fixatzero(int i) // 1-indexed!
{
  mn->Command(Form("REL %d", i));
  mn->Command(Form("SET LIM %d", i));
  mn->Command(Form("SET PAR %d 0", i));
  mn->Command(Form("FIX %d", i));
}

static bool onegoodminosup(const int par) // 0-indexed
{
  return getminerrup(par) > 0 && getminerrup(par) != 54321.0;
}

static bool onegoodminosdn(const int par) // 0-indexed
{ // I define all errors as positive
  return getminerrdn(par) > 0 && getminerrdn(par) != 54321.0;
}

static bool onegoodminos(const int par, const bool no_low_ok) // 0-indexed
{
  return (no_low_ok || onegoodminosdn(par)) && onegoodminosup(par);
}

// Get the MINOS error if present (regardless of whether, otherwise the
// MIGRAD/HESSE error
static double getbesterrup(const int par) // 0-indexed
{
  return onegoodminosup(par)? getminerrup(par): geterr(par);
}

static double getbesterrdn(const int par) // 0-indexed
{
  return onegoodminosdn(par)? getminerrdn(par): geterr(par);
}

static double min(const double a, const double b)
{
  return a < b? a: b;
}

__attribute__((unused)) static double max(const double a, const double b)
{
  return a > b? a: b;
}

static void fcn(__attribute__((unused)) int & np,
  __attribute__((unused)) double * gin, double & like, double *par,
  __attribute__((unused)) int flag)  
{
  like = 0;

  const double b12life = 29.14e3;

  // Assumption: all histograms have same dimensions and it don't change
  const int nxbins = fithist[0]->GetNbinsX();
  const int nybins = fithist[0]->GetNbinsY();
  const int maxtb = min(nxbins, maxfitt+nnegbins);

  // Cache all the bin contents. Yes, this is measured to make it faster
  // (~15%) as compared to calling GetBinContents in the inner loop.
  static std::vector< std::vector< std::vector<double> > > alldata;
  if(alldata.empty()){
    for(int period = 0; period < nperiod; period++){
      std::vector< std::vector<double> > vv;
      std::vector< double > emptyv;
      vv.push_back(emptyv); // 0->1
      for(int tb = 1; tb <= maxtb; tb++){
        std::vector<double> v;
        for(int eb = 0; eb < nybins; eb++)
          v.push_back(fithist[period]->GetBinContent(tb, eb+1));
        vv.push_back(v);
      }
      alldata.push_back(vv);
    }
  }

  static std::vector<double> xs;
  if(xs.empty()){
    xs.push_back(0); // 0->1
    for(int tb = 1; tb <= maxtb; tb++)
      xs.push_back(fithist[0]->GetXaxis()->GetBinCenter(tb));
  }

  // nesting these loops in this order is the fastest
  // Found to reduce time by 15% over other orders
  for(int tb = 1; tb <= maxtb; tb++){

    // loop invariants.  Who knows if ACLIC can optimize them?
    // Experimentally it would seem that it doesn't.
    const double x = xs[tb];
    const double sqrt_par_aneut_nc_x = sqrt(par[aneut_nc]/x);
    const double nstuff =
      // Neutron lifetime convolved with muon lifetime. Various mu-
      // capture times, plus thermalization, messed up by detector
      // effects
      1/(par[tneut_nc]-n_start_conv_time)
      * (exp(-x/par[tneut_nc]) - exp(-x/n_start_conv_time))
      // Diffusion
      *(erf(sqrt_par_aneut_nc_x)
          -M_2_SQRTPI*sqrt_par_aneut_nc_x*exp(-par[aneut_nc]/x));
    const double exp_x_b12life = exp(-x/b12life)/b12life;

    for(int period = 0; period < nperiod; period++){
      const int beam = period < nperiodrhc? 0: 1; // resist making cleverer

      for(int eb = 0; eb < nybins; eb++){
        const int off_beam   = beam  *nbins_e + eb;
        const int off_period = period*nbins_e + eb;

        double model = 0;

        const double sca = scales[period][eb];

        if(x <= -1 || x >= 2) model += fabs(par[flat_nc+off_period]);

        if(x >= 2){
          const double norm_n   =
            beam == 0 /* RHC */? 
            par[nneut_nc + nbins_e + eb]*par[nneut_nc+eb] /* R = F*(R/F) */:
            par[nneut_nc + nbins_e + eb]; /*F*/ 
          const double norm_b12 =
            beam == 0? 
            par[nb12_nc + nbins_e + eb]*par[nb12_nc+eb]:
            par[nb12_nc + nbins_e + eb];

          // Michel lifetime - ~average of mu+ and mu-, messed up by
          // detector effects
          model += fabs(par[nmich_nc+off_beam]*sca)/par[tmich_nc+off_beam]
                    * exp(-x/par[tmich_nc+off_beam]) +

                   fabs(norm_n  *sca) * nstuff
                + fabs(norm_b12*sca) * exp_x_b12life;
        }

        if((x >= -10 && x <= -1) || (x >= 2 && x <= 10))
          model += fabs(par[pileup_nc+off_period])*(10-fabs(x));

        const double data = alldata[period][tb][eb];

        like += model - data;
        if(model > 0 && data > 0) like += data * log(data/model);
      }
    }
  }

  const double tneut_penalty = 0.5 *
    pow((par[tneut_nc] - n_lifetime_nominal)/n_lifetime_priorerr, 2);

  const double aneut_penalty = 0.5 *
    pow((par[aneut_nc] - n_diffusion_nominal)/n_diffusion_priorerr, 2);

  like += tneut_penalty + aneut_penalty;
}

void make_mn()
{
  if(mn) delete mn;
  mn = new TMinuit(npar);
  mn->fGraphicsMode = false;
  const int print = 0;

  // not sure which of these works
  mn->SetPrintLevel(print);
  mn->Command(Form("SET PRINT %d", print));

  // Observed to help get MINOS errors out
  mn->Command("SET STRATEGY 2");

  mn->SetFCN(fcn);
  mn->Command("SET ERR 0.5"); // log likelihood

  int mnparmerr = 0;
  for(int i = 0; i < npar; i++)
    mn->mnparm(i, PARAMETERS[i].name,
                  PARAMETERS[i].start,
                  PARAMETERS[i].start/100., 0., 0., mnparmerr);
}

static void set_ee_to_mn(const int periodi, const int bin) // 0-indexed
{
  int beam;
  if(periodi < nperiodrhc) beam = 0;
  else                     beam = 1;

  for(int i = 0; i < ncommonpar; i++){
    ee_neg->SetParameter(i, getpar(i));
  }

  // walk forward through the blocks of parameters to the right place.
  for(int i = ncommonpar; i < ncommonpar+nperbinpar; i++){
    ee_neg->SetParameter(i,

      (i==2||i==3||i==5? scales[periodi][bin]:1)* // per track -> total

      ((i == 2 || i == 3) && beam == 0? // Check if this is a ratio
       getpar(ncommonpar +
              nbeam*nbins_e*(i-ncommonpar) +
              nbins_e +
              bin)
       : 1) *
      getpar(ncommonpar +
             nbeam*nbins_e*(i-ncommonpar) +
             beam * nbins_e +
             bin
      )
    );
  }

  for(int i = ncommonpar+nperbinpar; i < npar_ee-1; i++){
    ee_neg->SetParameter(i,
      getpar(ncommonpar +
             nbeam*nbins_e*nperbinpar + // end of per beam block
             nperiod*nbins_e * (i-ncommonpar-nperbinpar) + // par block
             nbins_e*periodi +
             bin));
  }

  ee_neg->SetParameter(ee_scale_par, 1); // scale

  for(int f = 0; f < ntf1s; f++)
    for(int i = 0; i < npar_ee; i++)
      ees[f]->SetParameter(i, ee_neg->GetParameter(i));
}

// Set the drawing functions for the sum of all periods in a beam type
static void set_ee_to_mn_beam(const int beam, const int bin)
{
  // Oh so hacky!
  double pars[npar_ee];
  for(int i = 0; i < (beam == 0?nperiodrhc:nperiodfhc); i++ ){
    const int idx = i + (beam == 1)*nperiodrhc;
    set_ee_to_mn(idx, bin);
    if(i == 0)
      for(int p = 0; p < npar_ee; p++)
        pars[p] = ee_pos->GetParameter(p);
    else
      for(int p = 0; p < npar_ee; p++)
        // Only the additive parameters
        if(p == 2 || p == 3 || p == 5 || p == 6 || p == 7 || p == 8)
          pars[p] += fabs(ee_pos->GetParameter(p));
  }

  for(int f = 0; f < ntf1s; f++)
    for(int i = 0; i < npar_ee; i++)
      ees[f]->SetParameter(i, pars[i]);
}

// Get rid of any bins that have partial contents.  Just for drawing.
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

static double beam_scale(const int beam, const int bin)
{
  double ans = 0;
  for(int i = 0; i < (beam == 0?nperiodrhc:nperiodfhc); i++ ){
    const int idx = i + (beam == 1)*nperiodrhc;
    ans += scales[i][bin];
  }
  return ans;
}

static void draw_ee_common(TH1D * x, const int rebin)
{
  x->Rebin(rebin);
  sanitize_after_rebinning(x);

  x->GetYaxis()->SetTitle(rebin== 1?"Delayed clusters/#mus":
                          Form("Delayed clusters/%d#mus", rebin));
  x->GetYaxis()->CenterTitle();
  x->GetXaxis()->CenterTitle();

  x->Draw("e");

  TLegend * leg = new TLegend(0.63, 0.7, 0.93, 0.98);
  for(int i = 0; i < ntf1s; i++){
    ees[i]->SetParameter(ee_scale_par, rebin);
    ees[i]->Draw("same");
    if(i != 0 /*ee_neg*/) leg->AddEntry(ees[i], ees_description[i], "l");
  }

  leg->SetTextFont(42);
  leg->SetBorderSize(1);
  leg->SetFillStyle(0);
  leg->Draw();


  c1->SetLogy();
  c1->Print("fit.pdf");
  c1->Update(); c1->Modified();
}

// Draw the sum of the fit histograms for a beam type and energy bin.
// 0-indexed
static void draw_ee_beam(const int beam, const int bin)
{
  c1->cd();

  const int rebin = 5;

  set_ee_to_mn_beam(beam, bin);

  // Get the form of the histogram and zero it
  TH1D * x = fithist[0]->ProjectionX("x", bin+1, bin+1 /* 0->1 */);
  x->Reset();

  for(int i = 0; i < (beam == 0?nperiodrhc:nperiodfhc); i++ ){
    const int idx = i + (beam == 1)*nperiodrhc;
    TH1D * tmp = fithist[idx]->ProjectionX("tmp", bin+1, bin+1);
    x->Add(tmp);
  }

  x->GetXaxis()->SetTitle(
    Form("All %s E_{#nu} bin %d: Time since muon stop (#mus)",
    beam == 0?"RHC":"FHC", bin));
  x->GetYaxis()->SetRangeUser(beam_scale(beam, bin)*1e-5*rebin,
                              beam_scale(beam, bin)* 0.1*rebin);

  draw_ee_common(x, rebin);
}

// Draw the fit histogram for a period and energy bin. 0-indexed bins
static void draw_ee(const int per, const int bin)
{
  c1->cd();

  const int rebin = 5;

  set_ee_to_mn(per, bin);

  TH1D * x = fithist[per]->ProjectionX("x", bin+1, bin+1 /* 0->1 */);

  x->GetXaxis()->SetTitle(
    Form("%s E_{#nu} bin %d: Time since muon stop (#mus)",
    Lperiodnames[per], bin));
  x->GetYaxis()->SetRangeUser(scales[per][bin]*1e-5*rebin,
                              scales[per][bin]* 0.1*rebin);
  draw_ee_common(x, rebin);
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
  for(int b = 0; b < nbeam; b++)
    for(int i = 0; i < nbins_e; i++)
      draw_ee(b, i, ntrack[beam][i]);
#endif
  /* End cheating for sensitivity study! */

static std::vector< std::vector<fitanswers> > dothefit()
{
  maxfitt = maxrealtime; // modified later for sensitivity study

  // Start each fit with a clean MINUIT slate.  This is crucial!
  make_mn();

  int status = 0;

  for(int period = 0; period < nperiod; period++){
    for(int bin = 0; bin < nbins_e; bin++){
      // May as well set this to near the right value
      mn->Command(Form("SET PAR %d %f", flat_nf+period*nbins_e + bin,
        fithist[period]->ProjectionX("x", bin+1, bin+1)
          ->Integral(0, nnegbins-10)/(nnegbins - 10)));
      // Again, something reasonable based on the data.
      mn->Command(Form("SET PAR %d %f", pileup_nf+period*nbins_e + bin,
        fithist[period]->ProjectionX("x", bin+1, bin+1)
          ->GetBinContent(nnegbins - 2)/8.));
    }
  }


  for(int beam = 0; beam < nbeam; beam++){
    for(int bin = 0; bin < nbins_e; bin++){
      // And this is a reasonable starting point for nmich
      mn->Command(Form("SET PAR %d %f", nmich_nf+beam*nbins_e + bin,
                       0.8));
      // Ditto nneut and B-12
      mn->Command(Form("SET PAR %d %f", nneut_nf+beam*nbins_e + bin,
                       beam == 0?0.2:0.1));

      mn->Command(Form("SET PAR %d %f", nb12_nf+beam*nbins_e + bin,
                       beam == 0?0.2:0.01));

      // Start with the muon lifetime fixed so that it doesn't try to
      // swap with the neutron lifetime.
      mn->Command(Form("FIX %d", tmich_nf+beam*nbins_e + bin));
    }
  }

  // Constraining neutron lifetime with a pull term, but also
  // set hard limits to hold it to something reasonable while
  // the rest of the fit settles.
  mn->Command(Form("SET LIM %d 30 80", tneut_nf));

  // Letting the neutron parameters float makes the fit not converge,
  // so don't do that.
  mn->Command(Form("FIX %d", tneut_nf));
  mn->Command(Form("FIX %d", aneut_nf));

  const int migrad_tries = 8;

  status = mn->Command("SIMPLEX 30000");
  for(int i = 0; i < migrad_tries; i++)
    if(0 == (status = mn->Command("MIGRAD 100000")))
      break;
  status = mn->Command("HESSE");

  // Now that we're (hopefully) converged, let muon lifetime float
  for(int beam = 0; beam < nbeam; beam++)
    for(int bin = 0; bin < nbins_e; bin++)
      mn->Command(Form("REL %d", tmich_nf+beam*nbins_e+bin));

  // And the neutron parameters
  mn->Command(Form("REL %d", tneut_nf));
  mn->Command(Form("REL %d", aneut_nf));

  // This becomes badly behaved if allowed to wander too far, so limit
  mn->Command(Form("SET LIM %d %f %f", aneut_nf,
    n_diffusion_nominal*0.2, n_diffusion_nominal*5));

  // Hold the Michel lifetime to something reasonable. Among other
  // concerns, this prevents it from swapping with the neutron lifetime,
  // supposing we let that float.
  for(int beam = 0; beam < nbeam; beam++)
    for(int bin = 0; bin < nbins_e; bin++)
      mn->Command(Form("SET LIM %d 1.6 2.6", tmich_nf+beam*nbins_e + bin));

  // MINOS errors are wrong if we used the abs trick, so limit
  for(int beam = 0; beam < nbeam; beam++){
    for(int bin = 0; bin < nbins_e; bin++){
      mn->Command(Form("SET PAR %d %f",
             nb12_nf+beam*nbins_e+bin, fabs(getpar( nb12_nc+beam*nbins_e+bin))));
      mn->Command(Form("SET PAR %d %f",
            nneut_nf+beam*nbins_e+bin, fabs(getpar(nneut_nc+beam*nbins_e+bin))));

      // Half of these are observed probabilities and the other half
      // are ratios. Same limits are reasonable at the moment.
      mn->Command(Form("SET LIM %d 0 10",  nb12_nf+beam*nbins_e+bin));
      mn->Command(Form("SET LIM %d 0 10", nneut_nf+beam*nbins_e+bin));
    }
  }

  status = mn->Command("SIMPLEX 30000");
  for(int i = 0; i < migrad_tries; i++)
    if(0 == (status = mn->Command("MIGRAD 100000")))
      break;
  status = mn->Command("HESSE");

  if(!status)
    // could skip MINOS errors for FHC and just do them for that ratio,
    // but then I can't make plots of the raw amount of each with MINOS errors
    for(int beam = 0; beam < nbeam; beam++)
      for(int bin = 0; bin < nbins_e; bin++){
        gMinuit->Command(Form("MINOS 50000 %d", nneut_nf+beam*nbins_e+bin));
        gMinuit->Command(Form("MINOS 50000 %d", nb12_nf +beam*nbins_e+bin));
      }

  gMinuit->Command("show min");

  std::vector< std::vector<fitanswers> > anses;

  for(int beam = 0; beam < nbeam; beam++){
    std::vector<fitanswers> anss;
    for(int bin = 0; bin < nbins_e; bin++){
      fitanswers ans;
      ans.n_good =   onegoodminos(nneut_nc+beam*nbins_e+bin, false);
      ans.n_mag =          getpar(nneut_nc+beam*nbins_e+bin);
      ans.n_mage_up = getbesterrup(nneut_nc+beam*nbins_e+bin);
      ans.n_mage_dn = getbesterrdn(nneut_nc+beam*nbins_e+bin);

      ans.b12_good =   onegoodminos(nb12_nc+beam*nbins_e+bin, true);
      ans.b12mag =           getpar(nb12_nc+beam*nbins_e+bin);
      ans.b12mage_up = getbesterrup(nb12_nc+beam*nbins_e+bin);
      ans.b12mage_dn = getbesterrdn(nb12_nc+beam*nbins_e+bin);

      printf("b12 = %f +%f -%f\n", ans.b12mag, ans.b12mage_up, ans.b12mage_dn);
      printf("n   = %f +%f -%f\n",  ans.n_mag,  ans.n_mage_up,  ans.n_mage_dn);

      anss.push_back(ans);
    }
    anses.push_back(anss);
  }

  return anses;
}

static double product_error(const double x, const double y,
                            const double xe, const double ye)
{
  return sqrt(y*y*pow(xe, 2) + x*x*pow(ye, 2));
}

__attribute__((unused))
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
        e->SetNpx(400);
        break;
      case 3:
        e->SetLineStyle(kDashed);
        e->SetLineColor(kBlue);
        e->SetLineWidth(1);
        e->SetNpx(400);
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

static double gdrawmax(const TGraphAsymmErrors * const g)
{
  double max = -1e300;
  for(int i = 0; i < g->GetN(); i++){
    const double y = g->GetY()[i] + g->GetErrorYhigh(i);
    if(y > max) max = y;
  }
  return max;
}

void rhc(const char * const savedhistfile = NULL)
{
  const std::string tcut = Form("i == 0 && %s", basecut);

  const std::string cut = Form(
   "%s && %s && t > %f && t < %f && !(t >= -1 && t < 2)",
    basecut, clustercut, -nnegbins, maxrealtime);

  init_ee();

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
    * n_result      = newgraph(kBlack, kSolid, kFullCircle, 1),
    * n_resultbad   = newgraph(kRed,   kSolid, kFullCircle, 1),
    * b12_result    = newgraph(kGray+1,kDashed,kOpenCircle, 1),
    * b12_resultbad = newgraph(kRed-10,kDashed,kOpenCircle, 1);

  if(savedhistfile == NULL){
    TFile * inputTFiles[nperiod] = { NULL };
    TTree * trees[nperiod] = { NULL };

    for(int i = 0; i < nperiod; i++){
      inputTFiles[i] = new TFile(inputfiles[i], "read");
      if(!inputTFiles[i] || inputTFiles[i]->IsZombie()){
        fprintf(stderr, "Couldn't read a file.  See above.\n");
        return;
      }
      trees[i] = dynamic_cast<TTree *>(inputTFiles[i]->Get("t"));
      if(!trees[i]){
        fprintf(stderr, "Couldn't read a tree.  See above.\n");
        return;
      }
    }

    std::ofstream o("savedhists.C");
    o << "{\n";
    for(int i = 0; i < nperiod; i++){
      all_tcounts[i] = new
        TH1D(Form("tcounts_%s", Speriodnames[i]), "", nbins_e, bins_e);
      fithist[i] = new TH2D(Form("%s_s", Speriodnames[i]), "",
        nnegbins + maxrealtime + additional,
        -nnegbins, maxrealtime + additional, nbins_e, bins_e);
      trees[i]->Draw(Form("slce:t >> %s_s", Speriodnames[i]), cut.c_str());
      trees[i]->Draw(Form("slce >> tcounts_%s", Speriodnames[i]), tcut.c_str());
      printf("Got %s_s\n", Speriodnames[i]);
      fflush(stdout);
      fithist[i]->SavePrimitive(o);
      all_tcounts[i]->SavePrimitive(o);
      //inputTFiles[i]->Close(); // XXX why does this seg fault?
    }
    o << "}\n";
    o.close();
  }
  else{
    gROOT->Macro(savedhistfile);
    for(int i = 0; i < nperiod; i++){
      if(NULL == (all_tcounts[i] = dynamic_cast<TH1D*>(
        gROOT->FindObject(Form("tcounts_%s", Speriodnames[i]))))){
        fprintf(stderr, "Could not read %s from %s\n",
               Form("tcounts_%s", Speriodnames[i]), savedhistfile);
        _exit(1);
      }
      if(NULL == (fithist[i] = dynamic_cast<TH2D*>(
        gROOT->FindObject(Form("%s_s", Speriodnames[i]))))){
        fprintf(stderr, "Could not read %s from %s\n",
               Form("%s_s", Speriodnames[i]), savedhistfile);
        _exit(1);
      }
    }
  }

  gErrorIgnoreLevel = kError; // after new TFile and Get above

  scales.resize(nperiod);
      
  for(int s = 1; s <= nbins_e; s++)
    for(int p = 0; p < nperiod; p++)
      scales[p].push_back(all_tcounts[p]->GetBinContent(s));

  const std::vector< std::vector<fitanswers> > anses = dothefit();

  for(int s = 0; s < nbins_e; s++){
    const fitanswers rof_ans = anses[0][s];
    const fitanswers fhc_ans = anses[1][s];

    const double loslce = fithist[0]->GetYaxis()->GetBinLowEdge(s+1);
    const double hislce = fithist[0]->GetYaxis()->GetBinLowEdge(s+2);

    const double graph_x  = (loslce+hislce)/2;
    const double graph_xe = (hislce-loslce)/2;
    const double graph_xoff = // visual offset
      graph_x + min(graph_xe/10, (bins_e[nbins_e]-bins_e[0])*0.02);

    // Directly fit for: easy
    addpoint(fhc_ans.n_good? g_n_fhc: g_n_fhc_bad, graph_x,
      fhc_ans.n_mag, graph_xe, fhc_ans.n_mage_dn, fhc_ans.n_mage_up);

    addpoint(fhc_ans.n_good? g_b12_fhc: g_b12_fhc_bad, graph_x,
      fhc_ans.b12mag, graph_xe, fhc_ans.b12mage_dn, fhc_ans.b12mage_up);

    printf("%s/%s RoF/FHC neutron (%4.2f-%4.2f)GeV: %.3f +%.3f -%.3f\n",
      rof_ans.n_good?"Good":"Bad",
      fhc_ans.n_good?"Good":"Bad", loslce, hislce,
      rof_ans.n_mag, rof_ans.n_mage_up, rof_ans.n_mage_dn);

    printf("%s/%s RoF/FHC B-12    (%4.2f-%4.2f)GeV: %.3f +%.3f -%.3f\n",
      rof_ans.b12_good?"Good":"Bad",
      fhc_ans.b12_good?"Good":"Bad", loslce, hislce,
      rof_ans.b12mag, rof_ans.b12mage_up, rof_ans.b12mage_dn);

    // Indirect, harder
    addpoint(rof_ans.n_good && fhc_ans.n_good? g_n_rhc: g_n_rhc_bad,
      graph_xoff,
      fhc_ans.n_mag * rof_ans.n_mag, graph_xe,
      product_error(fhc_ans.n_mag, rof_ans.n_mag,
                    fhc_ans.n_mage_dn, rof_ans.n_mage_dn),
      product_error(fhc_ans.n_mag, rof_ans.n_mag,
                    fhc_ans.n_mage_up, rof_ans.n_mage_up));

    addpoint(rof_ans.b12_good && fhc_ans.b12_good? g_b12_rhc: g_b12_rhc_bad,
      graph_xoff,
      fhc_ans.b12mag * rof_ans.b12mag, graph_xe,
      product_error(fhc_ans.b12mag, rof_ans.b12mag,
                    fhc_ans.b12mage_dn, rof_ans.b12mage_dn),
      product_error(fhc_ans.b12mag, rof_ans.b12mag,
                    fhc_ans.b12mage_up, rof_ans.b12mage_up));

    addpoint(rof_ans.n_good?n_result:n_resultbad, graph_x,
      rof_ans.n_mag, graph_xe, rof_ans.n_mage_dn, rof_ans.n_mage_up);

    addpoint(rof_ans.b12_good?b12_result:b12_resultbad,
      graph_xoff, rof_ans.b12mag,
      graph_xe, rof_ans.b12mage_dn, rof_ans.b12mage_up);
  }

  TH2D * dum = new TH2D("dm", "", 100, 0, bins_e[nbins_e], 10000, 0, 10);
  TH2D * dum2 = (TH2D*) dum->Clone("dm2");
  TH2D * dum3 = (TH2D*) dum->Clone("dm3");

  c4->cd();
  dum->GetYaxis()->SetTitle("RHC/FHC");
  dum->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  dum->GetYaxis()->CenterTitle();
  dum->GetXaxis()->CenterTitle();
  dum->GetYaxis()->SetRangeUser(0, 1.5);
  dum->Draw();
  dum->GetYaxis()->SetRangeUser(0, 
    min(2.5, 1.05*max(max(gdrawmax(  n_result), gdrawmax(  n_resultbad)),
             max(gdrawmax(b12_result), gdrawmax(b12_resultbad))))
  );
  n_result->Draw("pz");
  n_resultbad->Draw("pz");
  b12_result->Draw("pz");
  b12_resultbad->Draw("pz");
  TLegend * ratleg = new TLegend(0.6, 0.2, 0.85, 0.3);
  ratleg->SetTextFont(42);
  ratleg->AddEntry(n_result, "Neutrons", "lpe");
  ratleg->AddEntry(b12_result, "^{12}B", "lpe");
  ratleg->SetBorderSize(1);
  ratleg->SetFillStyle(0);
  ratleg->Draw();
  c4->Print("fit.pdf(");
 
  c2->cd();
  dum2->GetYaxis()->SetTitle("Neutrons per track");
  dum2->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  dum2->GetYaxis()->CenterTitle();
  dum2->GetXaxis()->CenterTitle();
  dum2->GetYaxis()->SetRangeUser(0, 
    min(2.5, 1.05*max(max(gdrawmax(g_n_rhc), gdrawmax(g_n_rhc_bad)),
                  max(gdrawmax(g_n_fhc), gdrawmax(g_n_rhc_bad))))
  );
  dum2->Draw();
  g_n_rhc->Draw("pz");
  g_n_rhc_bad->Draw("pz");
  g_n_fhc->Draw("pz");
  g_n_fhc_bad->Draw("pz");

  TLegend * nleg = new TLegend(0.6, 0.2, 0.9, 0.3);
  nleg->SetTextFont(42);
  nleg->AddEntry(g_n_fhc, "FHC", "lpe");
  nleg->AddEntry(g_n_rhc, "RHC", "lpe");
  nleg->SetBorderSize(1);
  nleg->SetFillStyle(0);
  nleg->Draw();

  c2->Print("fit.pdf");

  c3->cd();
  dum3->GetYaxis()->SetTitle("^{12}B per track");
  dum3->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  dum3->GetYaxis()->SetRangeUser(0, 
    min(2.5, 1.05*max(max(gdrawmax(g_b12_rhc), gdrawmax(g_b12_rhc_bad)),
                  max(gdrawmax(g_b12_fhc), gdrawmax(g_b12_rhc_bad))))
  );
  dum3->GetYaxis()->CenterTitle();
  dum3->GetXaxis()->CenterTitle();
  dum3->Draw();
  g_b12_rhc->Draw("pz");
  g_b12_rhc_bad->Draw("pz");
  g_b12_fhc->Draw("pz");
  g_b12_rhc_bad->Draw("pz");
  TLegend * b12leg = new TLegend(0.6, 0.2, 0.9, 0.3);
  b12leg->SetTextFont(42);
  b12leg->AddEntry(g_b12_fhc, "FHC", "lpe");
  b12leg->AddEntry(g_b12_rhc, "RHC", "lpe");
  b12leg->SetBorderSize(1);
  b12leg->SetFillStyle(0);
  b12leg->Draw();
  c3->Print("fit.pdf");

  for(int i = 0; i < nbins_e; i++)
    for(int beam = 0; beam < nbeam; beam++)
      draw_ee_beam(beam, i);

  for(int i = 0; i < nbins_e; i++)
    for(int period = 0; period < nperiod; period++)
      draw_ee(period, i);

  c4->Print("fit.pdf]");
}
