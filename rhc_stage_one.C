#include "TMatrixDSym.h"
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

TMinuit * mn = NULL;

#include "common.C"

bool muoncatcher = true;// set at entry point

double maxfitt = maxrealtime; // varies

struct fitanswers{
  bool n_good, b12_good;
  double  n_mag,  n_mage_up,  n_mage_dn,
         b12mag, b12mage_up, b12mage_dn;
};

double n_lifetime_nominal = 0; // set below

// From external considerations.  It doesn't substantively change the
// result of this program if a smaller number is put here, but it would
// be circular to do so in most cases, so don't.
double n_lifetime_priorerr = 0; // changed below

// This is the *effective* muon lifetime, with all detector effects
const double tmich_nominal = 2.1;
const double tmich_priorerr = 0.3;

double n_diffusion_nominal = 0; // set below using mindist

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
const int nperbinpar = 4; // number per energy bin parameters
const int nperperpar = 2; // number per period parameters

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

  for(int b = 0; b < nbeam; b++)
    for(int i = 0; i < nbins_e; i++)
      p.push_back(makepar(Form("%sNneut%d", bname[b], i), nneut_start));
  for(int b = 0; b < nbeam; b++)
    for(int i = 0; i < nbins_e; i++)
      p.push_back(makepar(Form("%sNB12_%d", bname[b], i), nb12_start));

  for(int b = 0; b < nbeam; b++)
    for(int i = 0; i < nbins_e; i++)
      p.push_back(makepar(Form("%sTMich%d", bname[b], i), tmich_nominal));
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

static TF1 * ee_flat =
  new TF1("ee_flat", "[8]*abs([6])", -nnegbins, maxrealtime+additional);

static TF1 * ee_mich =
  new TF1("ee_mich", "[8]*(abs([5])/[4] * exp(-x/[4]))",
          holex_hi, maxrealtime+additional);

static TF1 * ee_neut =
  new TF1("ee_mich",
   Form(
   // Convolution of muon lifetime and neutron lifetime
   "[8]*(abs([2])/([0]-%f) * (exp(-x/[0]) - exp(-x/%f))"

   // Neutron diffusion for a spherical cut.  Approximate for an 
   // intersection-of-two-cylinders cut
   "*(TMath::Erf(sqrt([1]/x))-2/sqrt(TMath::Pi())*sqrt([1]/x)*exp(-[1]/x)))",
  n_start_conv_time, n_start_conv_time),
   holex_hi, maxrealtime+additional);

static TF1 * ee_b12 =
  new TF1("ee_b12", "[8]*(abs([3])/29.14e3 * exp(-x/29.14e3))",
          holex_hi, maxrealtime+additional);

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
  holex_hi, maxrealtime+additional);

static TF1 * ee_neg = new TF1("ee_neg",
  "[8]*(abs([6]) + "
  "((x >= -10 && x <= 10))*(abs([7])*abs(abs(x)-10)))",
  -nnegbins, holex_lo);

const int ntf1s = 7;
static TF1 * ees[ntf1s] =
  { ee_neg, ee_pos, ee_flat, ee_mich, ee_neut, ee_b12, ee_pileup };

static const char * const ees_description[ntf1s] =
  { "Full fit", "Full fit", "Uncorrelated", "Michels", "Neutrons",
    "^{12}B", "pileup" };

static std::vector< std::vector<double> > scales;


static TCanvas * c1 = new TCanvas("rhc1", "rhc1"); // hist

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

        if(x <= holex_lo || x >= holex_hi) model += fabs(par[flat_nc+off_period]);

        if(x >= holex_hi){
          const double norm_n   =
            beam == 0 /* RHC */? 
            par[nneut_nc+eb] /* R */:
            par[nneut_nc + nbins_e + eb]; /*F*/ 
          const double norm_b12 =
            beam == 0? 
            par[nb12_nc+eb]:
            par[nb12_nc + nbins_e + eb];

          // Michel lifetime - ~average of mu+ and mu-, messed up by
          // detector effects
          model += fabs(par[nmich_nc+off_beam]*sca)/par[tmich_nc+off_beam]
                    * exp(-x/par[tmich_nc+off_beam]) +

                   fabs(norm_n  *sca) * nstuff
                + fabs(norm_b12*sca) * exp_x_b12life;
        }

        if((x >= -10 && x <= holex_lo) || (x >= holex_hi && x <= 10))
          model += fabs(par[pileup_nc+off_period])*(10-fabs(x));

        const double data = alldata[period][tb][eb];

        like += model - data;
        if(model > 0 && data > 0) like += data * log(data/model);
      }
    }
  }

  // penalty on Michel lifetimes.  Introduced when I started trying to 
  // implement cuts on later events based on Michels, which reduces
  // the lever arm of Michels considerably.
  for(int i = 0; i < nbeam*nbins_e; i++)
    like += 0.5 * pow((par[tmich_nc + i] - tmich_nominal)/tmich_priorerr, 2);

  // penalty on neutron lifetime
  like += 0.5 *
    pow((par[tneut_nc] - n_lifetime_nominal)/n_lifetime_priorerr, 2);

  // penalty on neutron diffusion constant - factor of 2 on a log scale
  like += 0.5 *
    pow((log(par[aneut_nc]) - log(n_diffusion_nominal))/M_LN2, 2);
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
  std::vector<PAR> PARAMETERS = makeparameters();
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
  for(int i = 0; i < (beam == 0?nperiodrhc:nperiodfhc); i++){
    const int idx = i + (beam == 1)*nperiodrhc;
    set_ee_to_mn(idx, bin);
    if(i == 0)
      for(int p = 0; p < npar_ee; p++)
        pars[p] = fabs(ee_pos->GetParameter(p));
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

static void draw_ee_common(TH1D * x, const int rebin, const char * const outname)
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

  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetLogy();
  c1->Print(outname);
  c1->Update(); c1->Modified();
}

// Draw the sum of the fit histograms for a beam type and energy bin.
// 0-indexed
static void draw_ee_beam(const int beam, const int bin, const char * const outname)
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
  x->GetYaxis()->SetRangeUser(beam_scale(beam, bin)*1e-6*rebin,
                              beam_scale(beam, bin)* 0.1*rebin);

  draw_ee_common(x, rebin, outname);
}

// Draw the fit histogram for a period and energy bin. 0-indexed bins
static void draw_ee(const int per, const int bin, const char * const outname)
{
  c1->cd();

  const int rebin = 5;

  set_ee_to_mn(per, bin);

  TH1D * x = fithist[per]->ProjectionX("x", bin+1, bin+1 /* 0->1 */);

  x->GetXaxis()->SetTitle(
    Form("%s E_{#nu} bin %d: Time since muon stop (#mus)",
    Lperiodnames[per], bin));
  x->GetYaxis()->SetRangeUser(scales[per][bin]*1e-6*rebin,
                              scales[per][bin]* 0.1*rebin);
  draw_ee_common(x, rebin, outname);
}

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

  status = mn->Command("SIMPLEX 100000");
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

      mn->Command(Form("SET LIM %d 0 10",  nb12_nf+beam*nbins_e+bin));
      mn->Command(Form("SET LIM %d 0 10", nneut_nf+beam*nbins_e+bin));
    }
  }

  status = mn->Command("SIMPLEX 100000");
  for(int i = 0; i < migrad_tries; i++)
    if(0 == (status = mn->Command("MIGRAD 100000")))
      break;
  status = mn->Command("HESSE");

  if(!status)
    for(int beam = 0; beam < nbeam; beam++)
      for(int bin = 0; bin < nbins_e; bin++){
        // XXX This is realy slow and only marginally useful
        //gMinuit->Command(Form("MINOS 50000 %d", nneut_nf+beam*nbins_e+bin));
        //gMinuit->Command(Form("MINOS 50000 %d", nb12_nf +beam*nbins_e+bin));
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

  // Just in case MINOS found a new minimum
  status = mn->Command("MIGRAD");
  status = mn->Command("HESSE");

  return anses;
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

static void save_for_stage_two(TGraphAsymmErrors * n_result,
                               TGraphAsymmErrors * b12_result,
                               TGraphAsymmErrors * g_n_rhc,
                               TGraphAsymmErrors * g_n_fhc,
                               TGraphAsymmErrors * g_b12_rhc,
                               TGraphAsymmErrors * g_b12_fhc,
                               const int mindist,
                               const int minslc, const int maxslc)
{
  ofstream for_stage_two(Form("for_stage_two_mindist%d_nslc%d_%d_%s.C",
                              mindist, minslc, maxslc, muoncatcher?"muoncatcher":"main"));
  for_stage_two << "{\n";
  g_n_rhc->SetName("g_n_rhc");
  g_n_fhc->SetName("g_n_fhc");
  g_n_rhc->SavePrimitive(for_stage_two);
  g_n_fhc->SavePrimitive(for_stage_two);
  g_b12_rhc->SetName("g_b12_rhc");
  g_b12_fhc->SetName("g_b12_fhc");
  g_b12_rhc->SavePrimitive(for_stage_two);
  g_b12_fhc->SavePrimitive(for_stage_two);

  // Following the example in TMinuitMinimizer::GetHessianMatrix()
  TMatrixDSym hessian(npar);
  mn->mnemat(hessian.GetMatrixArray(), npar);
  hessian.Invert();

  // Only care about the cross-terms for the 2*bins of the fit histograms
  for(int ii = 0; ii < nbins_e*nbeam*2 /* n and b12 */; ii++){
    for(int jj = ii; jj < nbins_e*nbeam*2; jj++){

      const int i = ii < nbins_e*nbeam? nneut_nc + ii:
               nb12_nc - nbins_e*nbeam + ii;
      const int j = jj < nbins_e*nbeam? nneut_nc + jj:
               nb12_nc - nbins_e*nbeam + jj;

      for_stage_two << "hessian[" << ii << "][" << jj << "] = "
        << hessian[i][j] << ";\n";
    }
    for_stage_two << "\n";
  }

  for_stage_two << "}\n";
}

void rhc_stage_one(const char * const savedhistfile, const int mindist,
                   const int minslc, const int maxslc, const string region)
{
  if(mindist < 0) return; // to compile only

  muoncatcher = region == "muoncatcher";

  /*
    Based on my MC, 400us is reasonable for a mindist of 6, but more like
    25us for a mindist of 2. However, it's murky because we don't measure
    the capture point, but rather the position where the gammas compton,
    and that usually only in 2D. Assuming 2D, I get more like 80us. And
    of course the functional form isn't quite right for 2D...
    
    Anyway, here's a stupid parameterization that skewers the 2D and 3D
    cases. Note the addition of 5/8 of a cellwidth because that's how
    much you get for mindist == 0 on average. In fcn() this is given a
    factor of 2 error. Mostly we just measure it.
  */
  n_diffusion_nominal = TWO_D_CUT?14.94 * pow(mindist + 0.625, 1.74)
                                 : 1.39 * pow(mindist + 0.625, 3.00);

  // From external Monte Carlo
  n_lifetime_priorerr = muoncatcher?10:5.;
  // Hackily provide information from the more data-rich fits to the 
  // data-poor fits so they don't spin out of control.
  if(mindist <= 2) n_lifetime_priorerr = 2.;

  // From a loose cut running of this program
  n_lifetime_nominal = muoncatcher?45:52.7;

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

  init_ee();

  // Square means RHC, solid means neutron, black means good fit
  TGraphAsymmErrors
    * g_n_rhc       = newgraph(kBlack, kSolid, kOpenSquare, 1),
    * g_n_fhc       = newgraph(kBlack, kSolid, kOpenCircle, 1),
    * g_b12_rhc     = newgraph(kBlack, kDashed,kOpenSquare, 1),
    * g_b12_fhc     = newgraph(kBlack, kDashed,kOpenCircle, 1),
    * n_result      = newgraph(kBlack, kSolid, kFullCircle, 1),
    * b12_result    = newgraph(kGray+1,kDashed,kOpenCircle, 1);

  gErrorIgnoreLevel = kError; // after new TFile and Get above

  scales.resize(nperiod);
      
  for(int s = 1; s <= nbins_e; s++)
    for(int p = 0; p < nperiod; p++)
      scales[p].push_back(all_tcounts[p]->GetBinContent(s));

  const std::vector< std::vector<fitanswers> > anses = dothefit();

  for(int s = 0; s < nbins_e; s++){
    const fitanswers rhc_ans = anses[0][s];
    const fitanswers fhc_ans = anses[1][s];

    const double loslce = fithist[0]->GetYaxis()->GetBinLowEdge(s+1);
    const double hislce = fithist[0]->GetYaxis()->GetBinLowEdge(s+2);

    const double graph_x  = (loslce+hislce)/2;
    const double graph_xe = (hislce-loslce)/2;
    const double graph_xoff = 
      min(graph_xe/30, (bins_e[nbins_e]-bins_e[0])*0.007);

    // Directly fit for: easy
    addpoint(g_n_rhc, graph_x,
      rhc_ans.n_mag, graph_xe, rhc_ans.n_mage_dn, rhc_ans.n_mage_up);

    addpoint(g_n_fhc, graph_x,
      fhc_ans.n_mag, graph_xe, fhc_ans.n_mage_dn, fhc_ans.n_mage_up);

    addpoint(g_b12_rhc, graph_x,
      rhc_ans.b12mag, graph_xe, rhc_ans.b12mage_dn, rhc_ans.b12mage_up);

    addpoint(g_b12_fhc, graph_x,
      fhc_ans.b12mag, graph_xe, fhc_ans.b12mage_dn, fhc_ans.b12mage_up);


    printf("%s RHC neutron (%5.3f-%5.3f)GeV: %.3f +%.3f -%.3f\n",
      rhc_ans.n_good?"Good":" Bad", loslce, hislce,
      rhc_ans.n_mag, rhc_ans.n_mage_up, rhc_ans.n_mage_dn);

    printf("%s RHC B-12    (%5.3f-%5.3f)GeV: %.3f +%.3f -%.3f\n",
      rhc_ans.b12_good?"Good":" Bad", loslce, hislce,
      rhc_ans.b12mag, rhc_ans.b12mage_up, rhc_ans.b12mage_dn);

    printf("%s FHC neutron (%5.3f-%5.3f)GeV: %.3f +%.3f -%.3f\n",
      fhc_ans.n_good?"Good":" Bad", loslce, hislce,
      fhc_ans.n_mag, fhc_ans.n_mage_up, fhc_ans.n_mage_dn);

    printf("%s FHC B-12    (%5.3f-%5.3f)GeV: %.3f +%.3f -%.3f\n",
      fhc_ans.b12_good?"Good":" Bad", loslce, hislce,
      fhc_ans.b12mag, fhc_ans.b12mage_up, fhc_ans.b12mage_dn);



    // Indirect, harder
    addpoint(n_result, graph_x,
      rhc_ans.n_mag / fhc_ans.n_mag, graph_xe,
      ratio_error(rhc_ans.n_mag, fhc_ans.n_mag,
                  rhc_ans.n_mage_dn, fhc_ans.n_mage_up),
      ratio_error(rhc_ans.n_mag, fhc_ans.n_mag,
                  rhc_ans.n_mage_up, fhc_ans.n_mage_dn));

    addpoint(b12_result, graph_x+graph_xoff,
      rhc_ans.b12mag / fhc_ans.b12mag, graph_xe,
      ratio_error(rhc_ans.b12mag, fhc_ans.b12mag,
                  rhc_ans.b12mage_dn, fhc_ans.b12mage_up),
      ratio_error(fhc_ans.b12mag, fhc_ans.b12mag,
                  rhc_ans.b12mage_up, fhc_ans.b12mage_dn));
  }

  const string filename = Form("fit_mindist%d_nslc%d_%d_%s.pdf", mindist, minslc, maxslc, region.c_str());

  c1->Print(Form("%s[", filename.c_str()));

  for(int i = 0; i < nbins_e; i++)
    for(int beam = 0; beam < nbeam; beam++)
      draw_ee_beam(beam, i, filename.c_str());

  for(int i = 0; i < nbins_e; i++)
    for(int period = 0; period < nperiod; period++)
      draw_ee(period, i, filename.c_str());

  c1->Print(Form("%s]", filename.c_str()));

  save_for_stage_two(n_result, b12_result, g_n_rhc, g_n_fhc,
                     g_b12_rhc, g_b12_fhc, mindist, minslc, maxslc);
}
