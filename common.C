// True if we are finding neutrons from 2D information alone, i.e not
// requiring hits in both x and y
#define TWO_D_CUT 1

//#define BGSUB

// Number of neutron pileup background samples per signal sample.  Set by my
// ntuple maker.  Set to zero to skip background subtraction.
const int bgmult =
#ifdef BGSUB
  1;
#else
  0;
#endif

// Everything is doubled because there is the signal and there is the
// off-space pileup background sample
const int SIG_AND_BG = bgmult?2:1;

// We aren't going to use anything between -1 and 2 microseconds because
// the detector conditions are just too awful.  And then don't use
// anything before holex_hi microseconds because I'm now using that region
// to define Michels that cut further interactions.
const double holex_lo = -1, holex_hi = 8;

const double nnegbins = 209;
const double maxrealtime = 269;
const double additional = 0;

// Restrict x and y to 100 and trkz_cutlo to 600
// as a blunt way to reduce pile-up subtraction problems.
const double trkx_cut = 170,
             trky_cut = 170,
             trkz_cutlo = 0,
             trkz_cuthi = 1250,
             mucatch_trky_cut = 45,
             mucatch_trkz_cutlo = 1310,
             mucatch_trkz_cuthi = 1560,
             mucatch_trkstartz_cut = 1100;
const double trklen_cut = 200;
const double remid_cut = 0.75;


// *** Parameters controlling how we deal with pileup neutrons from other ***
// *** neutrino interactions                                              ***

// How much weight to give the visible slices, as opposed to the intensity.
// i.e. What fraction of the pileup neutrons we think come from these
// rather than from interactions in the rock.
//
// From a Geant study using RHC MC.

const double npileup_sliceweight =
#if TWO_D_CUT == 1
  0.63;
#else
  0.79;
#endif

const double npileup_sliceweight_mucatch =
#if TWO_D_CUT == 1
  0.74;
#else
  0.80;
#endif

// Number of slices per 10**12 POT in RHC and FHC.  The ratio between these
// sets the relative number of neutrons we think are coming out of the rock.
// The absolute values are unimportant since we can tune either them or
// npileup_sliceweight.
const double slc_per_twp_rhc = 0.075;
const double slc_per_twp_fhc = 0.172;

// Maximum fraction of energy in the event that can be from the muon.
// This is a crude way to get at resolution bins.  To allow any resolution
// bin, set to 1.0 (or to be safe, a bit more, since I use a really crude
// approximation for muon energy).  To only allow roughly the worst bin,
// use ~0.55.
const double max_frac_e_mu = 1.5;

const int nperiodrhc = 1; // 4, 6
const int nperiodfhc = 1; // 1, 2, 3, 5
const int nperiod    = nperiodrhc + nperiodfhc;

const char * const inputfiles[nperiod] = {
  "201711-period46_rhc.root",
  "201711-period1235_fhc.root"
  //"periods46_rhc.root",
  //"periods1235_fhc.root"
/*  "prod_pid_S16-12-07_nd_period6_keepup/1947-type3.root",
  "prod_pid_R16-12-20-prod3recopreview.b_nd_numi_rhc_epoch4a_v1_goodruns/all-type3.root",

  "prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period1_v1_goodruns/all-type3.root",
  "prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period2_v1_goodruns/all-type3.root",
  "prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period3_v1_goodruns/all-type3.root",
  "prod_pid_R17-03-01-prod3reco.b_nd_numi_fhc_period5_v1_goodruns/all-type3.root"
*/
};

const char * const Speriodnames[SIG_AND_BG*nperiod] = {
"PR",
"PF",
#ifdef BGSUB
"PRBG",
"PFBG",
#endif
#if 0
    "P6",   "P4",   "P1",   "P2",   "P3",   "P5",
    "P6BG", "P4BG", "P1BG", "P2BG", "P3BG", "P5BG"
#endif
};

const char * const Lperiodnames[SIG_AND_BG*nperiod] = {
     "RHC",
     "FHC",
#ifdef BGSUB
     "RHC Pileup",
     "FHC Pileup",
#endif
#if 0
     "Period 6 (RHC)",
     "Period 4 (RHC)",
     "Period 1 (FHC)",
     "Period 2 (FHC)",
     "Period 3 (FHC)",
     "Period 5 (FHC)",
     "Period 6 (RHC) Pileup",
     "Period 4 (RHC) Pileup",
     "Period 1 (FHC) Pileup",
     "Period 2 (FHC) Pileup",
     "Period 3 (FHC) Pileup",
     "Period 5 (FHC) Pileup",
#endif
};

const int nbeam = 2; // not really generalizable as it stands

//const double bins_e[nbins_e+1] = {0.5, 1, 1.5, 2, 2.5, 3, 4.0, 5.0, 6.0 }; // 8

#if TWO_D_CUT == 1
  const int nbins_e = 6;
  const double bins_e[nbins_e+1] = { 0.5, 1.375, 2.250, 3.125, 4.0, 5.0, 6.0 };
#else
  const int nbins_e = 3;
  const double bins_e[nbins_e+1] = { 0.5, 1.5, 3.0, 6.0 };
#endif
#if 0
  const int nbins_e = 5;
  const double bins_e[nbins_e+1] = { 0.5, 1.6, 2.7, 3.8, 4.9, 6.0 };
#endif

double getpar(int i) // 0-indexed!
{
  double answer, dum;
  mn->GetParameter(i, answer, dum);
  return answer;
}

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

static double min(const double a, const double b)
{
  return a < b? a: b;
}

__attribute__((unused)) static double max(const double a, const double b)
{
  return a > b? a: b;
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


static void stylegraph(TGraphAsymmErrors * g, const int color,
                       const int linestyle, const int marker,
                       const int linewidth, const double markersize_)
{
  g->SetMarkerStyle(marker);
  g->SetLineStyle(linestyle);
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetLineWidth(linewidth);
  g->SetMarkerSize(markersize_);
}

static TGraphAsymmErrors *
newgraph(const int color, const int linestyle, const int marker,
         const int linewidth)
{
  TGraphAsymmErrors * g = new TGraphAsymmErrors;
  stylegraph(g, color, linestyle, marker, linewidth, 0.3);
  return g;
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

TBranch * setbranchaddress(const char * const name, float * d, TTree * t)
{
  t->SetBranchStatus(name, 1);
  t->SetBranchAddress(name, d);
  return t->GetBranch(name);
}

TBranch * setbranchaddress(const char * const name, int * d, TTree * t)
{
  t->SetBranchStatus(name, 1);
  t->SetBranchAddress(name, d);
  return t->GetBranch(name);
}
