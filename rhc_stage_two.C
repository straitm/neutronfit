#include "TLine.h"
#include "TMarker.h"
#include "TMinuit.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TMath.h"
#include "TLatex.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

// If defined, simplify output plots to combine the two
// categories of NC events.  No effect on the fit.
#define COMBINENC

using std::string;
using std::vector;

static const int npar = 9;
static TMinuit * mn = NULL;

#include "common.C"
#include "util.C"

static const int N_AND_B12 = 2;

static double ** residuals;

// filled by Macro call at start of running
// Apparently the way Macro() is implemented means this can't be static.
double ** hessian;

enum conttype { oned68, twod68, twod90 };

static const double tsize = 0.055;

static const double leftmargin = 0.135;
static const double topmargin  = tsize*1.33;
static const double rightmargin= 0.03;
static const double bottommargin=0.13;

static const double b12eff_nominal = 0.5;
static const double b12eff_error = 0.2;

/*********************************************************************/

static bool muoncatcher = true; // set at entry point
static two_or_three_d cut_dimensions = TWOD; // ditto
static int g_nbins_e = 0; // ditto

static double mum_nyield_nominal = 0;
static double mum_nyield_error   = 0;
static double neff_nominal = 0;
static double neff_error   = 0;
static double timing_eff_difference_for_pions = 0;

// Ratio of the neutron yield from stopping pions to that of muons, per
// stop, slightly adjusted by the greater fraction of neutrons from pions
// that undergo invisible interactions instead of being captured because of
// their higher energy.
//
// Multiplied by a pretty ad hoc correction for tracking difficulties, see
// below.
static double npimu_stop_nominal = 0;

// Not const because I want to adjust it to make different contours, and it
// is used by fcn().
//
// This is the error from physics added in quadrature with an error
// I made up for bad tracking which makes us miss the pions' neutrons
// more often than the muons'.
static double npimu_stop_error = 0;

/*********************************************************************/

// Ratio of the neutron yield from strongly interacting particles *in flight*
// to what Geant says.
static const double n_flight_nominal = 1;

// Somewhat arbitrary error on neutron production by things in
// flight (mostly pions) - factor of two on a log scale, i.e.
// 1 unit of chi2 for down 50% (half) or up 100% (double).
//const double n_flight_error = M_LN2;

// Put a 20% error on the actual yield. I justify this as follows:
// The real physics error on the *stopping* pion yield is 8%.  Geant
// is off by 22%, but I can explain 12% as it doing way too much
// capture on hydrogen, which is not a problem with in-flight
// interactions.  So without that mistake, Geant is off by about 1 sigma, i.e.
// we agree on the size of the error.  Let's suppose in-flight interactions
// are twice as hard.  Twice 8% plus a bit is 20%.
//
// Put a 20% error on the problems identifiying the neutron production
// point using the track end.  I justify this because it is double
// the error I put on the stopping pion case.
//
// Since this is a big error, I implement it as symmetric on a log scale.
//
// Silly extra work to make a error that averages to what I said on a linear scale
// but is asymmetric on that scale
static const double desired_average_n_flight_error = sqrt(pow(0.2, 2) + pow(0.2,2));
static const double so_this_is_the_error =
  (-2+2*desired_average_n_flight_error +
  sqrt(pow(2-2*desired_average_n_flight_error,2)
  + 8*desired_average_n_flight_error))/2;
static const double n_flight_error = log(1 + so_this_is_the_error);

// NOTE: Arbitrarily invent a small correlation coefficient for the neutron
// yield for in-flight interactions and stopped pions.
static const double corr_pi_stop_flight = 0.5;

static const double nm_nominal = 1;
static const double nm_error = 0.1;

static const double nc_nominal = 1;
static const double nc_error = 0.3;

static bool useb12 = false; // can be changed between fits

enum trktypes{ numubar, numu, stoppi, piflight, noneutrons,
               pileup /* not really a track type... */, MAXTRKTYPE };
static const char * const trktypenames[MAXTRKTYPE] =
  { "numubar", "numu", "stoppi", "piflight", "noneutrons", "pileup" };

static TH1D * makehist(const char * const name)
{
  return new TH1D(name, "", g_nbins_e, bins_e[cut_dimensions]);
}

static TH1D ** makehistset(const char * const basename)
{
  TH1D ** set = (TH1D **)malloc(MAXTRKTYPE * sizeof(TH1D*));
  for(int i = 0; i < MAXTRKTYPE; i++)
    set[i] = makehist(Form("%s_%s", basename, trktypenames[i]));
  return set;
}

static TH1D ** fhc_reco = NULL;
static TH1D ** rhc_reco = NULL;

static TH1D ** rhc_neut = NULL;
static TH1D ** rhc_b12  = NULL;
static TH1D ** fhc_neut = NULL;
static TH1D ** fhc_b12  = NULL;

// Total number of tracks in each beam's MC
static TH1D * fhc_tracks = NULL;
static TH1D * rhc_tracks = NULL;

// Makes either one or two TGraphAsymmErrors, depending on SIG_AND_BG
static TGraphAsymmErrors ** make_tgae_sigandbg_set()
{
  TGraphAsymmErrors ** set =
    (TGraphAsymmErrors **)malloc(SIG_AND_BG * sizeof(TGraphAsymmErrors*));
  for(int i = 0; i < SIG_AND_BG; i++) set[i] = new TGraphAsymmErrors;
  return set;
}

// Set via call to Macro()
TGraphAsymmErrors ** raw_g_n_rhc   = make_tgae_sigandbg_set();
TGraphAsymmErrors ** raw_g_n_fhc   = make_tgae_sigandbg_set();
TGraphAsymmErrors ** raw_g_b12_rhc = make_tgae_sigandbg_set();
TGraphAsymmErrors ** raw_g_b12_fhc = make_tgae_sigandbg_set();

static TGraphAsymmErrors * g_n_rhc = new TGraphAsymmErrors;
static TGraphAsymmErrors * g_n_fhc = new TGraphAsymmErrors;

static TGraphAsymmErrors * g_b12_rhc = new TGraphAsymmErrors;
static TGraphAsymmErrors * g_b12_fhc = new TGraphAsymmErrors;

// Out here so we can draw them
static TH1D * tot_rhc_neut = NULL;
static TH1D * tot_fhc_neut = NULL;
static TH1D * tot_rhc_b12  = NULL;
static TH1D * tot_fhc_b12  = NULL;

// All this has to happen after we get nbins_e from the entry point
static void init()
{
  residuals = (double **)malloc(sizeof(double *) * g_nbins_e*nbeam*N_AND_B12);
  hessian   = (double **)malloc(sizeof(double *) * g_nbins_e*nbeam*N_AND_B12);
  for(int i = 0; i < g_nbins_e*nbeam*N_AND_B12; i++){
    residuals[i] = (double *)malloc(sizeof(double) * g_nbins_e*nbeam*N_AND_B12);
    hessian  [i] = (double *)malloc(sizeof(double) * g_nbins_e*nbeam*N_AND_B12);
  }

  fhc_reco = makehistset("fhc_reco");
  rhc_reco = makehistset("rhc_reco");

  rhc_neut = makehistset("rhc_neut");
  rhc_b12  = makehistset("rhc_b12");
  fhc_neut = makehistset("fhc_neut");
  fhc_b12  = makehistset("fhc_b12");

  fhc_tracks = makehist("fhc_tracks");
  rhc_tracks = makehist("rhc_tracks");

  tot_rhc_neut = makehist("tot_rhc_neut");
  tot_fhc_neut = makehist("tot_fhc_neut");
  tot_rhc_b12  = makehist("tot_rhc_b12");
  tot_fhc_b12  = makehist("tot_fhc_b12");
}

static std::vector<double> getpars()
{
  std::vector<double> ans;
  for(int i = 0; i < npar; i++) ans.push_back(getpar(i));
  return ans;
}

static void reset_hists()
{
  tot_fhc_b12->Reset();
  tot_fhc_neut->Reset();
  tot_rhc_b12->Reset();
  tot_rhc_neut->Reset();

  for(int i = 0; i < MAXTRKTYPE; i++){
    fhc_b12[i]->Reset();
    fhc_neut[i]->Reset();
    rhc_b12[i]->Reset();
    rhc_neut[i]->Reset();
  }
}

static void update_hists(const double * const pars)
{
  const double ncscale         = pars[0];
  const double r_nmscale       = pars[1];
  const double real_npimu_stop = pars[2];
  const double neff            = pars[3];
  const double b12eff          = pars[4];
  const double mum_nyield      = pars[5];
  const double real_n_flight   = pars[6];
  const double pileup_rhc      = pars[7];
  const double pileup_fhc      = pars[8];

  const double npimu_stop = real_npimu_stop*timing_eff_difference_for_pions;
  const double n_flight   = real_n_flight  *timing_eff_difference_for_pions;

  // Probability of getting a neutron from a mu+ via
  // electrodisintigration, roughly.
  const double mup_nyield = 1e-4;

  // this could be another nuisance parameter.  Don't want to double count
  // it with the ratio of pion to muon yields, though, so careful.
  //const double muminus_capture_frac = 0.182;

  reset_hists();

  // Estimate of neutrons from RHC per muon, in total and for each truth
  rhc_neut[piflight]->Add(rhc_reco[piflight], ncscale*n_flight             *neff);
  rhc_neut[stoppi]  ->Add(rhc_reco[stoppi],   ncscale*npimu_stop*mum_nyield*neff);
  rhc_neut[numu]    ->Add(rhc_reco[numu],     r_nmscale         *mum_nyield*neff);
  rhc_neut[numubar] ->Add(rhc_reco[numubar],                     mup_nyield*neff);
  rhc_neut[pileup]  ->Add(rhc_reco[pileup],   pileup_rhc);

  // Ditto FHC
  fhc_neut[piflight]->Add(fhc_reco[piflight], ncscale*n_flight             *neff);
  fhc_neut[stoppi]  ->Add(fhc_reco[stoppi],   ncscale*npimu_stop*mum_nyield*neff);
  fhc_neut[numu]    ->Add(fhc_reco[numu],                        mum_nyield*neff);
  fhc_neut[numubar] ->Add(fhc_reco[numubar],                     mup_nyield*neff);
  fhc_neut[pileup]  ->Add(fhc_reco[pileup],   pileup_fhc);

  for(int i = 0; i < pileup; i++){
    tot_rhc_neut->Add(rhc_neut[i]);
    tot_fhc_neut->Add(fhc_neut[i]);
  }

  for(int i = 0; i < pileup; i++){
    rhc_neut[i]->Divide(rhc_tracks);
    fhc_neut[i]->Divide(fhc_tracks);
  }

  tot_rhc_neut->Divide(rhc_tracks);
  tot_fhc_neut->Divide(fhc_tracks);

  tot_rhc_neut->Add(rhc_neut[pileup]);
  tot_fhc_neut->Add(fhc_neut[pileup]);

  /*******************************************************************/

  const double carbon_acap_frac = 0.82; // my evaluation

  // nucl cap frac on C-12, and Double Chooz (including C-13)!
  // If it ever mattered, the error on this should be taken into account.
  const double mum_b12yield = carbon_acap_frac * (1-2028./2197.) * 0.177
  // Add in an ad-hoc estimate of how often capture results in a neutron
  // over the 12C(n,p)12B threshold (call it 5%) and then does that reaction
  // (call it 1% - cross section is ~30mb compared to few b elastic scattering).
  + carbon_acap_frac * (1-2028./2197.) * 0.05 * 0.01;

  // I think really zero, although of course the mu+ in flight (or the
  // e+ Michel) could make N-12, which looks the same (in fact, doubled
  // due to half the lifetime). I'd guess this is 1e-5 at most, and
  // maybe much lower.
  const double mup_b12yield = 1e-6;

  // H. Hilscher, W.-D. Krebs, G. Sepp, and V. Soergel. An experimental
  // test of the analogy between radiative pion absorption and muon
  // capture in 12C. Nuclear Physics A, 158(2):584-592, 1970
  // http://www.sciencedirect.com/science/article/pii/037594747090206X
  //
  // pi- have a lower B-12 yield despite always capturing!
  //
  // According to the same paper, pi- make a lot of Li-8. If you're not
  // careful, their figure 2 makes it look that only 20% of the activity
  // at t=0 is B-12, which would make the Li-8 yield almost 100%. But
  // the cycle time they used was much shorter than the Li-8 lifetime,
  // so this is activity produced over many cycles. My reanalysis, if I
  // assume zero uncorrelated background, says that the Li-8 yield is
  // about 4 times the B-12 yield, and with its half-life being 41 times
  // longer, it only adds about 9% to the activity near t=0. Add this in, but
  // note that since the error on the B-12 yield is 22%, it is negligible.
  const double hilscher_b12_per_pim_acap = 3.8e-3 * 1.09;  // +- 1.3e-3

  // Do not have to correct for the atomic capture ratio, because Hilscher's
  // figures are per stop, not per nuclear capture
  const double stop_pim_b12yield = carbon_acap_frac * hilscher_b12_per_pim_acap
  // Add in an ad-hoc estimate of how often capture results in a neutron
  // over the 12C(n,p)12B threshold (call it 10%) and then does that reaction
  // (call it 1% - cross section is ~30mb compared to few b elastic scattering).
  + carbon_acap_frac * 0.1 * 0.01;

  // I assume the probability is much lower for in flight interactions.
  // The factor of ten is only based on the number of fingers on my hands.
  const double flight_b12yield = stop_pim_b12yield / 10.
    // But add the same estimate for B-12 production via (n,p) as for stopped pions.
    + 0.1 * 0.01;

  // Estimate of number of B-12 from RHC per track
  rhc_b12[piflight]->Add(rhc_reco[piflight], ncscale*  flight_b12yield*b12eff);
  rhc_b12[stoppi]  ->Add(rhc_reco[stoppi],   ncscale*stop_pim_b12yield*b12eff);
  rhc_b12[numu]    ->Add(rhc_reco[numu],     r_nmscale*   mum_b12yield*b12eff);
  rhc_b12[numubar] ->Add(rhc_reco[numubar],               mup_b12yield*b12eff);

  // Ditto FHC
  fhc_b12[piflight]->Add(fhc_reco[piflight], ncscale*  flight_b12yield*b12eff);
  fhc_b12[stoppi]  ->Add(fhc_reco[stoppi],   ncscale*stop_pim_b12yield*b12eff);
  fhc_b12[numu]    ->Add(fhc_reco[numu],                  mum_b12yield*b12eff);
  fhc_b12[numubar] ->Add(fhc_reco[numubar],               mup_b12yield*b12eff);

  for(int i = 0; i < MAXTRKTYPE; i++){
    tot_rhc_b12->Add(rhc_b12[i]);
    tot_fhc_b12->Add(fhc_b12[i]);
  }
  for(int i = 0; i < MAXTRKTYPE; i++){
    rhc_b12[i]->Divide(rhc_tracks);
    fhc_b12[i]->Divide(fhc_tracks);
  }
  tot_rhc_b12->Divide(rhc_tracks);
  tot_fhc_b12->Divide(fhc_tracks);
}

static double compare(const double * const __restrict__ dat,
                      const double * const __restrict__ datup,
                      const double * const __restrict__ datdn,
                      const double * const __restrict__ model)
{
  double chi2 = 0;
  for(int i = 0; i < g_nbins_e*nbeam*N_AND_B12; i++){
    for(int j = i; j < g_nbins_e*nbeam*N_AND_B12; j++){
      if(!useb12 && (i >= g_nbins_e*nbeam || j >= g_nbins_e*nbeam)) continue;

      const double datai = dat[i];
      const double predi = model[i];
      const double dataj = dat[j];
      const double predj = model[j];

      // for diagonal elements, try to incorporate MINOS information
      const double rup = datup[i];
      const double rdn = datdn[i];
      const double correction = (i != j || rdn == 0 || rup == 0)?  1
      : (predi > datai?rup:rdn)/(rup/2 + rdn/2);

      // Need to count all off-diagonal elements twice since
      // we're only looking at the top half.
      const double offdiagonalmult = i == j? 1 : 2;

      // If the error and the hessian disagree sharply, it probably means that
      // there was a problem inverting the error matrix.  However, the error
      // itself can still be fine, and probably in this case the errors are
      // large and we can neglect the covariances. (Actually, I bet we could
      // always neglect the covariances...)
      if((rup*rup)*hessian[i][j] > 100 ||
         (i == j && (rup*rup)*hessian[i][j] < 1./100)){
        static int crazyprint = 0;
        if(i == j){
          chi2 += pow((predi - datai)/(rup/2 +rdn/2), 2);
          if(crazyprint++ < 100)
            fprintf(stderr, "Hessian for %d is crazy nuts: %g. "
                    "Using raw error instead.\n", i, hessian[i][i]);
        }
        else{
          if(crazyprint++ < 100)
            fprintf(stderr, "Hessian for %d %d is crazy nuts: %f.\n"
                   "Skipping.  Fit results dubious.\n", i, j, hessian[i][j]);
        }
        continue;
      }

      chi2 += offdiagonalmult * (residuals[i][j]
               = correction*hessian[i][j]*(predi - datai)*(predj - dataj));
    }
  }

  return chi2;
}

static double pi_penalty(const double npimu_stop, const double n_flight)
{
  const double delta_stop   = npimu_stop - npimu_stop_nominal;
  const double delta_flight = log(n_flight) - log(n_flight_nominal);

  const double cov = corr_pi_stop_flight * n_flight_error * npimu_stop_error;

  const int flight_and_stop = 2;

  const double det = pow(n_flight_error,2) * pow(npimu_stop_error,2) - cov*cov;
  const double phessian[flight_and_stop][flight_and_stop] = {
    { pow(n_flight_error, 2)/det, -cov/det                     },
    {  -cov/det                 , pow(npimu_stop_error, 2)/det }
  };

  double chi2 = 0;

  for(int i = 0; i < flight_and_stop; i++)
    for(int j = 0; j < flight_and_stop; j++)
      chi2 += phessian[i][j] * (i == 0?delta_stop: delta_flight)
                             * (j == 0?delta_stop: delta_flight);
  return chi2;
}

static void fcn(__attribute__((unused)) int & np,
  __attribute__((unused)) double * gin, double & chi2, double *par,
  __attribute__((unused)) int flag)
{
  chi2 = 0;

  update_hists(par);
  const double ncscale    = par[0];
  const double r_nmscale  = par[1];
  const double npimu_stop = par[2];
  const double neff       = par[3];
  const double b12eff     = par[4];
  const double mum_nyield = par[5];
  const double n_flight   = par[6];
  const double pileup_rhc = par[7];
  const double pileup_fhc = par[8];

  // penalty terms

  chi2 += pi_penalty(npimu_stop, n_flight);

  // Justified by external data
  chi2 += pow((mum_nyield - mum_nyield_nominal)/mum_nyield_error, 2);

  // not really justified
  //chi2 += pow((neff - neff_nominal)/neff_error, 2);
  //chi2 += pow((b12eff - b12eff_nominal)/b12eff_error, 2);

  // Can use this if we really think we have an external handle
  // on the numu contamination in RHC, but better to leave it out
  // and see that the fit finds a reasonable value
  //chi2 += pow((r_nmscale - nm_nominal)/nm_error, 2);

  // Can use this if we think we have an external handle
  // on the NC contamination, which we do, from GENIE
  //chi2 += pow((ncscale - nc_nominal)/nc_error, 2);

  static double * alldat   = (double *)malloc(sizeof(double)*g_nbins_e*nbeam*N_AND_B12),
                * alldatup = (double *)malloc(sizeof(double)*g_nbins_e*nbeam*N_AND_B12),
                * alldatdn = (double *)malloc(sizeof(double)*g_nbins_e*nbeam*N_AND_B12);
  static bool first = true;
  if(first){
    first = false;
    TGraphAsymmErrors * ingraphs[nbeam*N_AND_B12] = {
      g_n_rhc, g_n_fhc, g_b12_rhc, g_b12_fhc };

    int i = 0;
    for(int g = 0; g < nbeam*N_AND_B12; g++){
      for(int j = 0; j < ingraphs[g]->GetN(); j++){
        alldat[i] = ingraphs[g]->GetY()[j];
        alldatup[i] = ingraphs[g]->GetErrorYhigh(j);
        alldatdn[i] = ingraphs[g]->GetErrorYlow (j);
        i++;
      }
    }
  }

  double allmodel[g_nbins_e*nbeam*N_AND_B12];
  TH1D * inhists[nbeam*N_AND_B12] =
    { tot_rhc_neut, tot_fhc_neut, tot_rhc_b12, tot_fhc_b12 };
  int i = 0;
  for(int g = 0; g < nbeam*N_AND_B12; g++)
    for(int b = 1; b <= inhists[g]->GetNbinsX(); b++)
       allmodel[i++] = inhists[g]->GetBinContent(b);

  chi2 += compare(alldat, alldatup, alldatdn, allmodel);
}

/* Fill with the number of tracks in the MC for each category of truth */
void fill_hists(const char * const file, TH1D ** h, TH1D * tracks,
                const float minslc, const float maxslc)
{
  printf("Reading %s\n", file); fflush(stdout);

  TFile * f = new TFile(file, "read");
  if(!f || f->IsZombie()){
    fprintf(stderr, "Could not open the file %s\n", file);
    _exit(1);
  }

  TTree * t = dynamic_cast<TTree *>(f->Get("t"));
  if(!t){
    fprintf(stderr, "Could not get the tree from %s\n", file);
    _exit(1);
  }

  int i, primary, true_pdg, true_nupdg, true_nucc, contained,
      true_atom_cap, true_neutrons, nslc;
  float slce, trkx, trky, trkz, trkstartz, trklen, remid, timeleft, timeback, mcweight;
  t->SetBranchStatus("*", 0);

  setbranchaddress("slce", &slce, t);
  TBranch *       ibranch   = setbranchaddress("i", &i, t);
  TBranch *   primarybranch = setbranchaddress("primary", &primary, t);
  TBranch * containedbranch = setbranchaddress("contained", &contained, t);
  TBranch *      trkxbranch = setbranchaddress("trkx", &trkx, t);
  TBranch *      trkybranch = setbranchaddress("trky", &trky, t);
  TBranch *      trkzbranch = setbranchaddress("trkz", &trkz, t);
  TBranch * trkstartzbranch = setbranchaddress("trkstartz", &trkstartz, t);
  setbranchaddress("true_pdg", &true_pdg, t);
  setbranchaddress("true_nupdg", &true_nupdg, t);
  setbranchaddress("true_nucc", &true_nucc, t);
  setbranchaddress("trklen", &trklen, t);
  setbranchaddress("remid", &remid, t);
  setbranchaddress("timeleft", &timeleft, t);
  setbranchaddress("timeback", &timeback, t);
  setbranchaddress("true_atom_cap", &true_atom_cap, t);
  setbranchaddress("mcweight", &mcweight, t);
  setbranchaddress("true_neutrons", &true_neutrons, t);
  setbranchaddress("nslc", &nslc, t);

  const int MUMINUS =  13, MUPLUS  =  -13, // leptons are backwards
            EMINUS  =  11, EPLUS   =  -11, // from mesons!
            PIPLUS  = 211, PIMINUS = -211,
            KPLUS   = 321, KMINUS  = -321,
            PROTON = 2212,
            GAMMA = 22;
  const int MU = MUMINUS, E = EMINUS,
            PI = PIPLUS,  K = KPLUS;

  for(int e = 0; e < t->GetEntries(); e++){
    containedbranch->GetEntry(e);
    if(!contained) continue;

    primarybranch->GetEntry(e);
    if(!primary) continue;

    ibranch->GetEntry(e);
    if(i != 0) continue;

    trkzbranch->GetEntry(e);
    if(muoncatcher){
      if(trkz < mucatch_trkz_cutlo || trkz > mucatch_trkz_cuthi) continue;
      trkstartzbranch->GetEntry(e);
      if(trkstartz > mucatch_trkstartz_cut) continue;
    }
    else{
      if(trkz < trkz_cutlo || trkz > trkz_cuthi) continue;
    }

    trkxbranch->GetEntry(e);
    if(fabs(trkx) > trkx_cut) continue;

    trkybranch->GetEntry(e);
    if(muoncatcher){
      if(trky > mucatch_trky_cut || trky < -trky_cut) continue;
    }
    else{
      if(fabs(trky) > trky_cut) continue;
    }

    t->GetEntry(e);

    if(trklen < trklen_cut) continue;
    if(remid < remid_cut) continue;
    if(slce > bins_e[cut_dimensions][g_nbins_e] ||
       slce < bins_e[cut_dimensions][0]) continue;

    //if(nslc < minslc || nslc > maxslc) continue;

    // Should have no effect, unless there is cosmic overlay
    // In fact, MC times seem to be totally different, so never mind
    //if(timeleft < maxrealtime || timeback > -nnegbins) continue;

    // If it's not filled (perhaps for data or old MC), let it have no effect
    if(mcweight == 0) mcweight = 1;

    // based entirely on what the track actually is, NOT what the
    // neutrino interaction is. I think this makes more sense.

    // true_atom_cap is 1 if this is a mu-, pi- or K- that truly stops.
    // true_atom_cap is -1 if it is any particle that decays

    // Can I neglect the case that true mu- don't stop?
    // I think so, since Geant says that 99.90% of our mu- stop
    //
    // If the track is a pi- that decays in flight, it is now a mu- and
    // should be counted here.
    if(true_pdg == MUMINUS || (true_pdg == PIMINUS && true_atom_cap == -1))
      h[numu]->Fill(slce, mcweight);

    // This category is primarily for mu+ that stop and decay and
    // therefore mostly don't produce neutrons or B-12. But they do
    // produce neutrons just a little bit through electrodisintigration.
    // Put decaying pi+ into this category too, since they produce mu+
    // which then decay. Include both pi+ decaying at rest and pi+ and
    // pi- decaying in flight. The in-flight decays complicate things,
    // but this whole electrodisintigration thing is a very minor
    // contribution to the neutron yield, so don't worry about it. Also
    // lump decaying kaons in here as another very minor thing.
    else if(true_pdg == MUPLUS ||
            (true_atom_cap == -1 && (true_pdg == PIPLUS || abs(true_pdg) == K)))
      h[numubar]->Fill(slce, mcweight);

    // Lump K- into this category; they are very rare.
    else if(true_atom_cap == 1 && (true_pdg == KMINUS || true_pdg == PIMINUS))
       h[stoppi]->Fill(slce, mcweight);

    // Strongly interacting particles that interact in flight. Mostly
    // pions, with a substantial contribution from protons and a very
    // few kaons. Include nuclei too, even though they are never the
    // primary track for contained events (but do appear as tracks). For
    // all these, I'm taking the Geant estimate of the effective neutron
    // yield, which is questionable...
    else if(true_atom_cap == 0 &&
            (abs(true_pdg) == PI || abs(true_pdg) == K ||
             abs(true_pdg) == PROTON || abs(true_pdg) > 1000000000))
      h[piflight]->Fill(slce, mcweight * true_neutrons /* note */);

    else if(abs(true_pdg) == E || abs(true_pdg) == GAMMA)
      h[noneutrons]->Fill(slce, mcweight); // not used at this time

    else
      fprintf(stderr, "Unused track. It's ok, but: PDG = %d, true_atom_cap = %d\n",
              true_pdg, true_atom_cap);

    tracks->Fill(slce, mcweight);
  }

  for(int i = 1; i <= h[pileup]->GetNbinsX(); i++)
    h[pileup]->SetBinContent(i, 1);

  f->Close();
}

void make_mn()
{
  if(mn) delete mn;
  mn = new TMinuit(npar);
  mn->fGraphicsMode = false;

  // MUST be at least zero to get the contours out because MINUIT
  // decides whether you really want contours by the print level...!
  const int print = 0;

  // not sure which of these works
  mn->SetPrintLevel(print);
  mn->Command(Form("SET PRINT %d", print));

  // Observed to help get MINOS errors out
  mn->Command("SET STRATEGY 2");

  mn->SetFCN(fcn);

  int mnparmerr = 0;
  mn->mnparm(0, "NCscale", 1, 0.2, 0, 0, mnparmerr);
  mn->mnparm(1, "NMscale", 1, 0.1, 0, 0, mnparmerr);
  mn->mnparm(2, "npimu_stop", npimu_stop_nominal, npimu_stop_error,
             0, 0, mnparmerr);
  mn->mnparm(3, "neff", neff_nominal, neff_error, 0, 0, mnparmerr);
  mn->mnparm(4, "b12eff", b12eff_nominal, b12eff_error, 0, 0, mnparmerr);
  mn->mnparm(5, "mum_nyield", mum_nyield_nominal, mum_nyield_error,
                0, 0, mnparmerr);
  mn->mnparm(6, "n_flight", n_flight_nominal, 0.1,
             0, 0, mnparmerr);
  mn->mnparm(7, "pileup_rhc", rhc_tracks->GetMaximum()/10.,
                              rhc_tracks->GetMaximum()/100., 0, 0, mnparmerr);
  mn->mnparm(8, "pileup_fhc", fhc_tracks->GetMaximum()/10.,
                              fhc_tracks->GetMaximum()/100., 0, 0, mnparmerr);
}

static void styleleg(TLegend * leg)
{
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(tsize);
  leg->SetNColumns(1);
}

TGraph * wrest_contour(const char * const name)
{
  TGraph * cont = dynamic_cast<TGraph *>(mn->GetPlot());
  if(cont){
    cont = (TGraph *)cont->Clone(name);
    cont->SetLineWidth(2);
  }
  return cont;
}

void stylehist(TH1 * h, const int color, const int width,
               const int style = kSolid)
{
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetLineWidth(width);
  h->SetLineStyle(style);
  h->GetYaxis()->CenterTitle();
  h->GetXaxis()->CenterTitle();
}

void stylehistset(TH1D ** hh, const int color)
{
  stylehist(hh[piflight], color, 3, kDashed);
  stylehist(hh[stoppi],   color, 2, kDotted);
  stylehist(hh[numu],     color, 3, 9);
  stylehist(hh[numubar],
    color-9 /*probably bad except for blue and red */, 2, kDashDotted);
  stylehist(hh[pileup], color-9, 3, kDashDotted);
}

static void draw_with_visible_errors_and_leak_memory(TGraphAsymmErrors * g)
{
  g->Draw("pz");
  static int i = 0;
  TGraphAsymmErrors * foo = (TGraphAsymmErrors *)g->Clone(Form("foo-%d", i++));
  foo->SetMarkerSize(0);
  foo->Draw("pz");
}

static double gdrawmin(const TGraphAsymmErrors * const g)
{
  double min = 1e300;
  for(int i = 0; i < g->GetN(); i++){
    const double y = g->GetY()[i] - g->GetErrorYlow(i);
    if(y < min) min = y;
  }
  return min;
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

  // cd away so these Draws don't get printed to the PDF output
  TCanvas driveyarnspoon;

  // *Within* the allowed range (minslc to maxslc), find the mean.  Definition
  // must stay in sync with that in pass_intensity() in rhc_stage_one.C.
  rhc.Draw(Form("(nslc-2)*%f + %f * pot >> slc_r",
      npileup_sliceweight[cut_dimensions],
      slc_per_twp_rhc*(1-npileup_sliceweight[cut_dimensions])),
    Form("i == 0 && contained && primary && "
         "abs((nslc-2)*%f + %f * pot  -  %f) < %f",
         npileup_sliceweight[cut_dimensions],
         slc_per_twp_rhc*(1-npileup_sliceweight[cut_dimensions]),
         (minslc+maxslc)/2, (maxslc-minslc)/2));

  fhc.Draw(Form("(nslc-2)*%f + %f * pot >> slc_f",
      npileup_sliceweight[cut_dimensions],
      slc_per_twp_fhc*(1-npileup_sliceweight[cut_dimensions])),
    Form("i == 0 && contained && primary && "
         "abs((nslc-2)*%f + %f * pot  -  %f) < %f",
         npileup_sliceweight[cut_dimensions],
         slc_per_twp_fhc*(1-npileup_sliceweight[cut_dimensions]),
         (minslc+maxslc)/2, (maxslc-minslc)/2));

  printf("Mean for %4.1f-%4.1f slices: RHC %7d events: %4.2f, FHC %7d: %4.2f\n",
         minslc, maxslc,
         int(slc_r.GetEntries()), slc_r.GetMean(),
         int(slc_f.GetEntries()), slc_f.GetMean());

  return (slc_r.GetMean() + slc_f.GetMean())/2;
}

void draw(const int mindist, const float minslc, const float maxslc)
{
  gStyle->SetOptStat(0);
  gStyle->SetFrameLineWidth(2);

  TLatex t(0, 0, "NOvA Preliminary");
  t.SetTextColor(kBlue);
  t.SetTextSize(tsize);
  t.SetTextFont(42);
  t.SetNDC();
  t.SetTextAlign(33);
  t.SetX(1-rightmargin-0.005);
  t.SetY(1-0.02);

  //////////////////////////////////////////////////////////////////////
  const bool logy = false;
  TH2D * dum = new TH2D("dm", "", 100, 0, 10, 10000, logy?1e-4:-2.0, 2.0);
  dum->GetYaxis()->SetTitleSize(tsize);
  dum->GetXaxis()->SetTitleSize(tsize);
  dum->GetYaxis()->SetLabelSize(tsize);
  dum->GetXaxis()->SetLabelSize(tsize);
  TH2D * dum2 = (TH2D*) dum->Clone("dm2");

  //////////////////////////////////////////////////////////////////////
  stylegraph(g_n_rhc, kRed +3, kSolid, kOpenSquare, 2, 1.0);
  stylegraph(g_n_fhc, kBlue+3, kSolid, kOpenCircle, 2, 1.0);

  stylehist(tot_rhc_neut, kRed, 2);
  stylehistset( rhc_neut, kRed);

  stylehist(tot_fhc_neut, kBlue, 2);
  stylehistset( fhc_neut, kBlue);

  dum2->GetXaxis()->SetRangeUser(bins_e[cut_dimensions][0],
                                 bins_e[cut_dimensions][g_nbins_e]);
  dum2->GetYaxis()->SetRangeUser(
    min(0, 1.01*min(min(gdrawmin(g_n_rhc), tot_rhc_neut->GetMinimum()),
             min(gdrawmin(g_n_fhc), tot_fhc_neut->GetMinimum()))),
    1.1*max(max(gdrawmax(g_n_rhc), tot_rhc_neut->GetMaximum()),
            max(gdrawmax(g_n_fhc), tot_fhc_neut->GetMaximum()))
            );
  dum2->GetYaxis()->SetTitle("Visible neutrons/track");
  dum2->GetYaxis()->SetTitleOffset(1.25);
  dum2->GetXaxis()->SetTitle("Reconstructed E_{#nu} CC (GeV)");
  dum2->GetYaxis()->CenterTitle();
  dum2->GetXaxis()->CenterTitle();

  //
  TCanvas * c2r = new TCanvas("rhc2r", "rhc2r");
  const std::string outpdfname =
    Form("fit_stage_two_mindist%d_nslc%.1f_%.1f_%s_%s.pdf", mindist, minslc, maxslc,
         muoncatcher?"muoncatcher":"main", two_or_three_d_names[cut_dimensions]);
  c2r->Print((outpdfname + "[").c_str());
  c2r->SetLogy(logy);
  c2r->SetMargin(leftmargin, rightmargin, bottommargin, topmargin);
  dum2->Draw();

  const double one_entry_leg_y1  = 0.90;
  const double four_entry_leg_y1 = 1-topmargin-0.3
#ifdef COMBINENC
    +0.05
#endif
  ;
  const double leg_y2 = 1-topmargin-0.01;
  const double leg_x1 = leftmargin + 0.02;
  const double leg_x2 = leftmargin + 0.33;

  // true if you want a copy of each plot that doesn't yet have the fit
  // histograms drawn on it.
  const bool print_wo_fit = false;

  draw_with_visible_errors_and_leak_memory(g_n_rhc);
  TLegend * leg = new TLegend(leg_x1,
                              print_wo_fit?one_entry_leg_y1:four_entry_leg_y1,
                              leg_x2,
                              leg_y2);
  styleleg(leg);
  leg->AddEntry(g_n_rhc, "RHC data", "lpe");
  leg->Draw();

  if(print_wo_fit) c2r->Print(outpdfname.c_str());

  const int ndrawhists =
#ifdef COMBINENC
    3;
#else
    4;
#endif

  TH1D * c2hists[ndrawhists] = {
    tot_rhc_neut,
    rhc_neut[piflight],
    rhc_neut[numu],
#ifndef COMBINENC
    rhc_neut[stoppi],
#endif
  };

  for(int i = 1; i <= tot_rhc_neut->GetNbinsX(); i++) // wait, what?
    for(int h = 0; h < ndrawhists; h++)
      c2hists[h]->Draw("histsame][");

  leg->SetY1NDC(four_entry_leg_y1);
  leg->AddEntry(tot_rhc_neut, "RHC fit", "l");
#ifdef COMBINENC
  // Pretty sketchy, but turns out to be ok since they are all updated
  // during the fit anyway.
  rhc_neut[piflight]->Add(rhc_neut[stoppi]);
  leg->AddEntry(rhc_neut[piflight], "RHC neutral current", "l");
#else
  leg->AddEntry(rhc_neut[piflight], "RHC #pi^{#pm}, p in flight", "l");
  leg->AddEntry(rhc_neut[stoppi], "RHC stopped #pi^{#minus}", "l");
#endif
  leg->AddEntry(rhc_neut[numu], "RHC #nu_{#mu}", "l");

  t.Draw();
  c2r->Print(outpdfname.c_str());

  //
  TCanvas * c2f = new TCanvas("rhc2f", "rhc2f");
  c2f->SetLogy(logy);
  c2f->SetMargin(leftmargin, rightmargin, bottommargin, topmargin);
  dum2->Draw();

  draw_with_visible_errors_and_leak_memory(g_n_fhc);
  TLegend * legf = new TLegend(leg_x1,
                               print_wo_fit?one_entry_leg_y1:four_entry_leg_y1,
                               leg_x2,
                               leg_y2);
  styleleg(legf);
  legf->AddEntry(g_n_fhc, "FHC data", "lpe");
  legf->Draw();
  t.Draw();
  if(print_wo_fit) c2f->Print(outpdfname.c_str());

  TH1D * c2histsf[ndrawhists] = {
    tot_fhc_neut,
    fhc_neut[piflight],
    fhc_neut[numu],
#ifndef COMBINENC
    fhc_neut[stoppi],
#endif
  };
  for(int i = 1; i <= tot_rhc_neut->GetNbinsX(); i++)
    for(int h = 0; h < ndrawhists; h++)
      c2histsf[h]->Draw("histsame][");

  legf->SetY1NDC(four_entry_leg_y1);
  legf->AddEntry(tot_fhc_neut, "FHC fit", "l");
#ifdef COMBINENC
  fhc_neut[piflight]->Add(fhc_neut[stoppi]);
  legf->AddEntry(fhc_neut[piflight], "FHC neutral current", "l");
#else
  legf->AddEntry(fhc_neut[piflight], "FHC #pi^{#pm}, p in flight", "l");
  legf->AddEntry(fhc_neut[stoppi], "FHC stopped #pi^{#minus}", "l");
#endif
  legf->AddEntry(fhc_neut[numu], "FHC #nu_{#mu}", "l");

  c2f->Print(outpdfname.c_str());

  //////////////////////////////////////////////////////////////////////
  TCanvas * c3 = new TCanvas("rhc3", "rhc3");
  c3->SetMargin(leftmargin, rightmargin, bottommargin, topmargin);

  const double ncmin = 0, ncmax = 2.4,
    nmmin = 0.0, nmmax = 2.6;
  TH2D * dum3 = new TH2D("dm3", "", 100, ncmin, ncmax, 1, nmmin, nmmax);
  dum3->GetYaxis()->SetTitleSize(tsize);
  dum3->GetXaxis()->SetTitleSize(tsize);
  dum3->GetYaxis()->SetLabelSize(tsize);
  dum3->GetXaxis()->SetLabelSize(tsize);
  dum3->GetXaxis()->SetTitle("FHC & RHC NC scale relative to MC");
  dum3->GetYaxis()->SetTitle("RHC #nu_{#mu} scale relative to MC");
  dum3->GetXaxis()->CenterTitle();
  dum3->GetYaxis()->CenterTitle();
  dum3->Draw();

  TLine * l1 = new TLine(1, nmmin, 1, nmmax);
  l1->SetLineColor(kGray);
  l1->Draw();
  TLine * l2 = new TLine(ncmin, 1, ncmax, 1);
  l2->SetLineColor(kGray);
  l2->Draw();

  mn->fGraphicsMode = true;

  TMarker * bestfit = new TMarker(getpar(0), getpar(1), kFullCircle);

  mn->Command(Form("SET ERR %f", TMath::ChisquareQuantile(0.90, 2)));

  mn->Command("MINOS 10000 1");
  mn->Command("MINOS 10000 2");
  mn->Command("MNCONT 1 2 50");
  TGraph * cont_full = wrest_contour("cont_full");

  const double twodminerrdn = getminerrdn(1);
  const double twodminerrup = getminerrup(1);
  const double twodminerrleft  = getminerrdn(0);
  const double twodminerrright = getminerrup(0);

  /*useb12 = false;
  mn->Command("MIGRAD");
  mn->Command("MNCONT 1 2 50");
  TGraph * cont_nob12 = wrest_contour("cont_nob12");

  useb12 = true; */
  const double proper_npimu_stop_error = npimu_stop_error;
  npimu_stop_error = 1.1;
  mn->Command("MIGRAD");
  mn->Command("MNCONT 1 2 50");
  TGraph * cont_perfect_nuclear = wrest_contour("cont_perfect_nuclear");

  mn->Command(Form("SET PAR 3 %f", npimu_stop_nominal));
  mn->Command("FIX 3");
  mn->Command("MIGRAD");
  mn->Command("MNCONT 1 2 50");
  TGraph * cont_perfect_ratio = wrest_contour("cont_perfect_ratio");

  mn->Command(Form("SET PAR 7 %f", n_flight_nominal));
  mn->Command("FIX 7");
  mn->Command("MIGRAD");
  mn->Command("MNCONT 1 2 50");
  TGraph * cont_double_perfect_ratio = wrest_contour("cont_double_perfect_ratio");

  // restore
  npimu_stop_error = proper_npimu_stop_error;
  mn->Command("REL 3");
  mn->Command("REL 7");
  mn->Command("MIGRAD");

  if(cont_double_perfect_ratio != NULL){
    cont_double_perfect_ratio->SetLineColor(kOrange+2);
    cont_double_perfect_ratio->SetLineStyle(7);
    // get rid of the visual gap
    cont_double_perfect_ratio->SetPoint(cont_double_perfect_ratio->GetN(),
                                 cont_double_perfect_ratio->GetX()[0],
                                 cont_double_perfect_ratio->GetY()[0]);
    cont_double_perfect_ratio->Draw("l");
  }
  if(cont_perfect_ratio != NULL){
    cont_perfect_ratio->SetLineColor(kViolet);
    cont_perfect_ratio->SetLineStyle(kDashed);
    // get rid of the visual gap
    cont_perfect_ratio->SetPoint(cont_perfect_ratio->GetN(),
                                 cont_perfect_ratio->GetX()[0],
                                 cont_perfect_ratio->GetY()[0]);
    cont_perfect_ratio->Draw("l");
  }
  if(cont_perfect_nuclear != NULL){
    cont_perfect_nuclear->SetLineColor(kRed);
    cont_perfect_nuclear->SetLineStyle(kDotted);
    cont_perfect_nuclear->SetPoint(cont_perfect_nuclear->GetN(),
                                   cont_perfect_nuclear->GetX()[0],
                                   cont_perfect_nuclear->GetY()[0]);
    cont_perfect_nuclear->Draw("l");
  }
  if(cont_full != NULL){
    cont_full->SetLineColor(kBlue-2);
    cont_full->SetLineStyle(kSolid);
    cont_full->SetPoint(cont_full->GetN(),
                        cont_full->GetX()[0],
                        cont_full->GetY()[0]);
    cont_full->Draw("l");
  }
  /*if(cont_nob12 != NULL){
    cont_nob12->SetLineColor(kBlue);
    cont_nob12->SetLineStyle(7);
    cont_nob12->Draw("l");
  } */

  mn->fGraphicsMode = false;

  bestfit->Draw();

  mn->Command("MINOS 10000 1");
  mn->Command("MINOS 10000 2");

  TGraphAsymmErrors * onederrs = new TGraphAsymmErrors;
  onederrs->SetLineWidth(2);
  onederrs->SetPoint(0, getpar(0), getpar(1) - twodminerrdn*1.2);
  onederrs->SetPoint(1, getpar(0)-twodminerrleft*1.2, getpar(1));

  mn->Command("SET ERR 1"); // get one D 68% errors
  mn->Command("MIGRAD");
  mn->Command("MINOS 10000 1");
  mn->Command("MINOS 10000 2");
  mn->Command("MINOS 10000 4");
  if(useb12) mn->Command("MINOS 10000 5");

  printf("NC %s mindist %d slcrange %.1f - %.1f : %f + %f - %f meanslc %f %s\n",
    muoncatcher?"muoncatcher":"main",
    mindist, minslc, maxslc, getpar(0), getbesterrup(0), getbesterrdn(0),
    mean_slice(false, minslc, maxslc),
    two_or_three_d_names[cut_dimensions]);

  printf("NM %s mindist %d slcrange %.1f - %.1f : %f + %f - %f meanslc %f %s\n",
    muoncatcher?"muoncatcher":"main",
    mindist, minslc, maxslc, getpar(1), getbesterrup(1), getbesterrdn(1),
    mean_slice(true , minslc, maxslc),
    two_or_three_d_names[cut_dimensions]);
  c3->cd(); // mean_slice cd()s away

  onederrs->SetPointError(0, getminerrdn(0), getminerrup(0), 0, 0);
  onederrs->SetPointError(1, 0, 0, getminerrdn(1), getminerrup(1));
  onederrs->SetMarkerStyle(kFullCircle);
  onederrs->Draw("pz");

  mn->Command(Form("SET ERR %f", TMath::ChisquareQuantile(0.90, 2))); // put back

  leg = new TLegend(leftmargin, 0.67, 1-rightmargin, 1-topmargin);
  leg->SetTextFont(42);
  leg->SetBorderSize(1);
  leg->SetTextSize(tsize*9./11.);
  leg->SetFillStyle(1001);
  leg->SetMargin(0.1);
  if(cont_full != NULL){
    leg->AddEntry((TH1D*)NULL, "", "");
    leg->AddEntry(cont_full,
#if 0
      Form("90%%, effective stopped #pi^{#minus}/#mu^{#minus} n yield %.1f#pm%.1f; "
           "#pi,p in-flight %.2f^{+%.2f}_{-%.2f} #times Geant",
      npimu_stop_nominal, npimu_stop_error,
      n_flight_nominal, n_flight_nominal*exp(n_flight_error) - 1,
      1 - n_flight_nominal/exp(n_flight_error))
#else
    "90% CL"
#endif
     , "l");
    leg->AddEntry(onederrs, "1D errors", "lp");
  }
  if(cont_perfect_nuclear != NULL)
    leg->AddEntry(cont_perfect_nuclear,
                  "perfectly known nuclear stopping #pi/#mu neutron yield", "l");
  if(cont_perfect_ratio != NULL)
    leg->AddEntry(cont_perfect_ratio,
                  "and perfectly known atomic capture ratios", "l");
  if(cont_double_perfect_ratio != NULL)
    leg->AddEntry(cont_double_perfect_ratio,
                  "and perfectly known in-flight neutron yields", "l");
  leg->Draw();
  t.Draw();
  c3->Print(outpdfname.c_str());

  //////////////////////////////////////////////////////////////////////
  TCanvas * c4 = new TCanvas("rhc4", "rhc4");
  c4->SetMargin(leftmargin, rightmargin, bottommargin, topmargin);

  stylehistset(fhc_reco, kBlue);
  stylehistset(rhc_reco, kRed);

  rhc_reco[numu]->GetYaxis()->SetTitle("Probability/bin");
  rhc_reco[numu]->GetXaxis()->SetTitle("Reconstructed E_{#nu} CC (GeV)");

  // neutron producers
  TH1 * norm1 = fhc_reco[numu]->DrawNormalized("hist");
  TH1 * norm2 = rhc_reco[numu]->DrawNormalized("histsame");
  TH1 * norm3 = fhc_reco[stoppi]->DrawNormalized("histsame");
  TH1 * norm4 = rhc_reco[stoppi]->DrawNormalized("histsame");
  TH1 * norm5 = fhc_reco[piflight]->DrawNormalized("ehistsame");
  TH1 * norm6 = rhc_reco[piflight]->DrawNormalized("histsame");

  if(norm1 && norm2 && norm3 && norm4 && norm5 && norm6){
    double maxy = max(norm1->GetMaximum(), norm2->GetMaximum());
    maxy = max(maxy, norm3->GetMaximum());
    maxy = max(maxy, norm4->GetMaximum());
    maxy = max(maxy, norm5->GetMaximum());
    maxy = max(maxy, norm6->GetMaximum());
    norm1->GetYaxis()->SetRangeUser(0, maxy*1.1);


    TLegend * leg4 = new TLegend(leg_x1, 1-topmargin, leg_x2, leg_y2);
    leg4->AddEntry(fhc_reco[numu], "#mu^{#minus} in FHC", "l");
    leg4->AddEntry(rhc_reco[numu], "#mu^{#minus} in RHC", "l");
    leg4->AddEntry(fhc_reco[stoppi], "Stopped #pi^{#minus} in FHC", "l");
    leg4->AddEntry(rhc_reco[stoppi], "Stopped #pi^{#minus} in RHC", "l");
    leg4->AddEntry(fhc_reco[piflight], "#pi^{#pm}, p in flight in FHC", "l");
    leg4->AddEntry(rhc_reco[piflight], "#pi^{#pm}, p in flight in RHC", "l");
    styleleg(leg4);
    leg4->Draw();

    t.Draw();
    c4->Print(outpdfname.c_str());
  }


  TH2D * hresid = new TH2D("hresid", "",
    N_AND_B12*nbeam*g_nbins_e, 0, N_AND_B12*nbeam*g_nbins_e,
    N_AND_B12*nbeam*g_nbins_e, 0, N_AND_B12*nbeam*g_nbins_e);
  mn->Command("CALL");

  // If an input variable is up against a limit and the fit results
  // are garbage, check if one of these is huge.  If so, probably
  // the hessian is messed up.

  for(int i = 0; i < g_nbins_e*nbeam*N_AND_B12; i++)
    for(int j = 0; j < g_nbins_e*nbeam*N_AND_B12; j++)
      hresid->SetBinContent(i+1, j+1, fabs(residuals[i][j]));

  c4->SetLogz();
  c4->SetRightMargin(0.18);
  hresid->Draw("colz");
  t.Draw();
  c4->Print(outpdfname.c_str());

  c4->Print((outpdfname + "]").c_str()); // doesn't print anything, just closes
}

// Subtract pileup background histograms from sig+bg histograms.  If not using
// pileup background histograms (bgmult == 0), no subtraction occurs.
static void do_background_subtraction()
{
  TGraphAsymmErrors ** raw_gs[nbeam*N_AND_B12]
    = { raw_g_n_rhc, raw_g_n_fhc, raw_g_b12_rhc, raw_g_b12_fhc };
  TGraphAsymmErrors *      gs[nbeam*N_AND_B12]
    = {     g_n_rhc,     g_n_fhc,     g_b12_rhc,     g_b12_fhc };

  for(int g = 0; g < nbeam*N_AND_B12; g++){
    for(int i = 0; i < raw_gs[g][0]->GetN(); i++){
      gs[g]->SetPoint(i, raw_gs[g][0]->GetX()[i],
          raw_gs[g][0]->GetY()[i]
        - (bgmult==0?0:raw_gs[g][1]->GetY()[i]/bgmult));

      gs[g]->SetPointError(i, raw_gs[g][0]->GetErrorXlow (i),
                              raw_gs[g][0]->GetErrorXhigh(i),
        // Down error is signal down error (+) bg up error
        sqrt(pow(raw_gs[g][0]->GetErrorYlow (i), 2) +
             (bgmult==0?0:pow(raw_gs[g][1]->GetErrorYhigh(i)/bgmult, 2))),
        // Up error is signal up error (+) bg down error
        sqrt(pow(raw_gs[g][0]->GetErrorYhigh(i), 2) +
             (bgmult==0?0:pow(raw_gs[g][1]->GetErrorYlow (i)/bgmult, 2))));
    }
  }
}

void rhc_stage_two(const char * const input, const int mindist,
                   const float minslc, const float maxslc, const string region,
                   const two_or_three_d cut_dimensions_)
{
  if(mindist < 0) return; // to compile only

  muoncatcher = region == "muoncatcher";
  cut_dimensions = cut_dimensions_;
  g_nbins_e = nbins_e[cut_dimensions];

  init();

  mum_nyield_nominal = muoncatcher?0.855:0.1810;
  mum_nyield_error   = muoncatcher?0.036:0.0305;

  neff_nominal = 0.5  * (muoncatcher?0.5:1);
  neff_error   = 0.05 * (muoncatcher?0.5:1);

  const double rough_neut_lifetime = muoncatcher?60:52.; // us
  const double mean_time_of_muon_capture_weighted_by_neut_yield =
  muoncatcher?0.21:1.07; // us
  timing_eff_difference_for_pions =
  exp(-mean_time_of_muon_capture_weighted_by_neut_yield/rough_neut_lifetime);

  // Ratio of the neutron yield from stopping pions to that of muons, per
  // stop, slightly adjusted by the greater fraction of neutrons from pions
  // that undergo invisible interactions instead of being captured because of
  // their higher energy.
  //
  // Multiplied by a pretty ad hoc correction for tracking difficulties, see
  // below.
  npimu_stop_nominal = (muoncatcher? 3.98 : 14.29) *
                       (muoncatcher? 0.9  : 0.75);

  // Not const because I want to adjust it to make different contours, and it
  // is used by fcn().
  //
  // This is the error from physics added in quadrature with an error
  // I made up for bad tracking which makes us miss the pions' neutrons
  // more often than the muons'.
  npimu_stop_error = sqrt(pow(muoncatcher?0.6:2.52  , 2)
                        + pow(npimu_stop_nominal*0.1, 2));


  TH1::SetDefaultSumw2();

  // I don't need to scale these by POT because everything is done in a ratio
  // between the number of primary tracks in the given beam type and the truth
  // of those tracks.  But it might be confusing if the individual histograms
  // are drawn without normalization.
  fill_hists("fhc_mc/all-type3.root", fhc_reco, fhc_tracks, minslc, maxslc);
  fill_hists("rhc_mc/all-type3.root", rhc_reco, rhc_tracks, minslc, maxslc);

  gROOT->Macro(input);

  do_background_subtraction();

  make_mn();

  mn->Command("SET LIM 1 0 50"); // NC scale
  mn->Command("SET LIM 2 0 50");  // numubar scale

  // stopped pion to muon neutron yield ratio
  mn->Command(Form("SET LIM 3 %f %f",
    max(0, npimu_stop_nominal - 5*npimu_stop_error),
           npimu_stop_nominal + 5*npimu_stop_error));

  mn->Command("SET LIM 4 0.01 1.05"); // neutron efficiency
  mn->Command("SET LIM 5 0.01 1.05"); // B-12 efficiency
  mn->Command("SET LIM 6 0.01 1.5"); // mu- neutron yield

  // other-in-flight neutron yield
  mn->Command("SET LIM 7 0 20");

  // pile-up
  mn->Command(Form("SET PAR 8 %f", 0.));
  mn->Command(Form("SET PAR 9 %f", 0.));
  mn->Command(Form("FIX 8"));
  mn->Command(Form("FIX 9")); // Let's try later

  if(!useb12) mn->Command("FIX 5");
  mn->Command("MINIMIZE 10000");
  mn->Command("HESSE");
  mn->Command("CALL");
  update_hists(&(getpars()[0]));

  draw(mindist, minslc, maxslc);
}
