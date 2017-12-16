#include "TMatrixDSym.h"
#include "TMinuit.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TError.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TText.h"
#include "TGaxis.h"
#include <fstream>

using std::string;
using std::vector;

TMinuit * mn = NULL;

#include "common.C"
#include "util.C"

enum beam_type { RHC_BEAM, FHC_BEAM };

// It turns out, I think, that the main contributor to events with a
// longer lifetime than simple neutron capture are not B-12 decay but
// rather neutron capture from neutrons that have spent a long time
// bouncing around in the rock and/or in the air around the detector.
// This doesn't have a simple exponential form, but from the Monte
// Carlo, this lifetime is a rough approximation for the major relevant
// component. If you want to put in the B-12 lifetime instead, it is
// 29.14e3 microseconds;
const double b12life = 0.5e3;

// Control output plots
const float textsize = 0.05;
const int rebin = 5;

static TH2D ** fithist     = (TH2D**)malloc(nperiod*SIG_AND_BG*sizeof(TH2D*));
static TH1D ** all_tcounts = (TH1D**)malloc(nperiod*SIG_AND_BG*sizeof(TH1D*));

static bool muoncatcher = true;// set at entry point
static int g_nbins_e = 0; // ditto

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
// Note, though, that if half of the neutrons are from pions and/or protons,
// this is more like 0.9us. If that number is used, you get a 1.4%
// increase in the number of neutrons inferred, i.e. for a fixed number
// here, the "efficiency" for pion-produced neutrons is about 1.4% lower
// than muon-produced mpared to muons.
//
// I'll try to get the shape right, and just note this efficiency
// difference.
static const double n_start_conv_time = (1.60 + 0.9)/2;

// Length of a NuMI spill.  Controls prompt pile-up component of the fit.
static const double numi_len = 9.6; // us

// The meaning of the visual y axis in terms of number of selected
// cluster per track per bin.
static const double yminpertrack = 1e-6,
                    ymaxpertrack = 0.1;

struct PAR {
  char * name;
  double start;
};

static const double topmargin = textsize*1.33,
                 bottommargin = 0.11 * textsize/0.045,
                   leftmargin = 0.09 * textsize/0.045,
                  rightmargin = 0.10 * textsize/0.045;

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
const int nper_e_beam_sg = 4; // number per energy bin, beam, and sig/bg
const int nper_e_period = 2; // number per energy bin and period

enum PARnames { ee_Tneut_par, ee_Aneut_par, ee_Nneut_par, ee_NB12_par,
                ee_Tmich_par, ee_NMich_par, ee_flat_par,  ee_pileup_par,
                ee_scale_par, npar_ee };

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

const double nmich_start = 0.8;
const double nneut_start = 0.05;
const double nb12_start  = 0.001;
const double flat_start = 3;
const double pileup_start = 300;

// parameter numbers, C numbering, for use with TMinuit functions
// Set in makeparameters() once we know the number of energy bins.
static int tneut_nc,
           aneut_nc,
           nneut_nc,
           nb12_nc,
           tmich_nc,
           nmich_nc,
           flat_nc,
           pileup_nc,
           npar;

// and FORTRAN numbering, for use with native MINUIT commands through
// TMinuit::Command()
static int flat_nf,
           nmich_nf,
           tmich_nf,
           nneut_nf,
           tneut_nf,
           aneut_nf,
           nb12_nf,
           pileup_nf;

static std::vector<PAR> makeparameters()
{
  // Same for all histograms
  tneut_nc  = 0,
  aneut_nc  = 1,

  // Parameters of interest: energy, beam, and sig/bg dependent
  nneut_nc  = aneut_nc + 1,
  nb12_nc   = nneut_nc + nbeam*g_nbins_e*SIG_AND_BG,

  // Nuisance parameters: energy, beam, and sig/bg dependent
  tmich_nc  = nb12_nc  + nbeam*g_nbins_e*SIG_AND_BG,
  nmich_nc  = tmich_nc + nbeam*g_nbins_e*SIG_AND_BG,

  // Nuisance parameters: energy and period dependent
  // In principle, *not* dependent on sig vs. off-space bg
  flat_nc   = nmich_nc + nbeam *g_nbins_e*SIG_AND_BG,
  pileup_nc = flat_nc  + nperiod*g_nbins_e,

  npar      = pileup_nc+ nperiod*g_nbins_e;

  flat_nf   = flat_nc  +1,
  nmich_nf  = nmich_nc +1,
  tmich_nf  = tmich_nc +1,
  nneut_nf  = nneut_nc +1,
  tneut_nf  = tneut_nc +1,
  aneut_nf  = aneut_nc +1,
  nb12_nf   = nb12_nc  +1,
  pileup_nf = pileup_nc+1;

  std::vector<PAR> p;
  p.push_back(makepar("Tneut", n_lifetime_nominal));
  p.push_back(makepar("Aneut", n_diffusion_nominal));

  const char * const bname[nbeam] = { "R_", "F_" };

  for(int b = 0; b < nbeam; b++)
    for(int i = 0; i < g_nbins_e; i++)
      for(int sigorbg = 0; sigorbg < SIG_AND_BG; sigorbg++)
        p.push_back(makepar(Form("%sNn%d-%s",
          bname[b], i, sigorbg?"B":"S"), nneut_start));
  for(int b = 0; b < nbeam; b++)
    for(int i = 0; i < g_nbins_e; i++)
      for(int sigorbg = 0; sigorbg < SIG_AND_BG; sigorbg++)
        p.push_back(makepar(Form("%sNB12_%d-%s",
          bname[b], i, sigorbg?"B":"S"), nb12_start));

  for(int b = 0; b < nbeam; b++)
    for(int i = 0; i < g_nbins_e; i++)
      for(int sigorbg = 0; sigorbg < SIG_AND_BG; sigorbg++)
        p.push_back(makepar(Form("%sTMu%d-%s",
          bname[b], i, sigorbg?"B":"S"), tmich_nominal));
  for(int b = 0; b < nbeam; b++)
    for(int i = 0; i < g_nbins_e; i++)
      for(int sigorbg = 0; sigorbg < SIG_AND_BG; sigorbg++)
        p.push_back(makepar(Form("%sNMu%d-%s",
          bname[b], i, sigorbg?"B":"S"), nmich_start));

  for(int ip = 0; ip < nperiod; ip++)
    for(int i = 0; i < g_nbins_e; i++)
      p.push_back(makepar(Form("%sflat%d", Speriodnames[ip], i), flat_start));
  for(int ip = 0; ip < nperiod; ip++)
    for(int i = 0; i < g_nbins_e; i++)
      p.push_back(makepar(Form("%spile%d",Speriodnames[ip],i),pileup_start));

  return p;
}

static TF1 * ee_flat =
  new TF1("ee_flat", "[8]*abs([6])", -nnegbins, maxrealtime+additional);

static TF1 * ee_mich =
  new TF1("ee_mich", "[8]*(abs([5])/[4] * exp(-x/[4]))",
          0, maxrealtime+additional);

static TF1 * ee_neut =
  new TF1("ee_mich",
   Form(
   // Convolution of muon lifetime and neutron lifetime
   "[8]*(abs([2])/([0]-%f) * (exp(-x/[0]) - exp(-x/%f))"

   // Neutron diffusion for a spherical cut.  Approximate for an
   // intersection-of-two-cylinders cut
   "*(TMath::Erf(sqrt([1]/x))-2/sqrt(TMath::Pi())*sqrt([1]/x)*exp(-[1]/x)))",
  n_start_conv_time, n_start_conv_time),
   0, maxrealtime+additional);

static TF1 * ee_b12 =
  new TF1("ee_b12", Form("[8]*(abs([3])/%f * exp(-x/%f))", b12life, b12life),
          0, maxrealtime+additional);

static TF1 * ee_pileup = new TF1("ee_pileup",
  Form("[8]*((x >= -%f && x <= %f))*(abs([7])*abs(abs(x)-%f))",
       numi_len, numi_len, numi_len),
  -nnegbins, maxrealtime+additional);

// The sum of all the above. TF1::Draw() has a lot of trouble with the
// discontinuity at zero, so split into the positive and negative parts.
static TF1 * ee_pos = new TF1("ee_pos",
  Form("[8]*(abs([6]) + "
    "(x >= 0)*("
      "abs([5])/[4] * exp(-x/[4]) + "
      "abs([2])/([0]-%f) * (exp(-x/[0]) - exp(-x/%f))"
      "*(TMath::Erf(sqrt([1]/x))-2/sqrt(TMath::Pi())*sqrt([1]/x)*exp(-[1]/x))"
      "+ abs([3])/%f * exp(-x/%f)"
    ") + "
    "((x >= -%f && x <= %f))*(abs([7])*abs(abs(x)-%f)))",
    n_start_conv_time, n_start_conv_time, b12life, b12life,
    numi_len, numi_len, numi_len),
  0, maxrealtime+additional);

static TF1 * ee_neg = new TF1("ee_neg",
  Form("[8]*(abs([6]) + "
       "((x >= -%f && x <= %f))*(abs([7])*abs(abs(x)-%f)))",
       numi_len, numi_len, numi_len),
  -nnegbins, 0);

const int ntf1s = 7; // Number of TF1s that we are going to draw
static TF1 * ees[ntf1s] =
  { ee_neg, ee_pos, ee_flat, ee_mich, ee_neut, ee_b12, ee_pileup };

static const char * const ees_description[ntf1s] =
  { "Full fit", "Full fit", "Uncorrelated bkg",
    "Michel decays", "Neutron captures",
    "Air neutron pileup", "Prompt pileup" };

static std::vector< std::vector<double> > scales;


static TCanvas * c1 = new TCanvas("rhc1", "rhc1"); // hist

static int get_off_beam(const int beam, const int energybin, const int isbg)
{
  return (beam *g_nbins_e + energybin)*SIG_AND_BG + isbg;
}

static int get_off_period(const int period, int energybin)
{
  return period*g_nbins_e + energybin;
}

static void fcn(__attribute__((unused)) int & np,
  __attribute__((unused)) double * gin, double & like, double *par,
  __attribute__((unused)) int flag)
{
  like = 0;

  // Assumption: all histograms have same dimensions and they don't change
  const int ntbins = fithist[0]->GetNbinsX();
  const int nenergybins = fithist[0]->GetNbinsY();
  const int maxtb = min(ntbins, maxfitt+nnegbins);

  // Cache all the bin contents. Yes, this is measured to make it faster
  // (~15%) as compared to calling GetBinContents in the inner loop.
  static std::vector< std::vector< std::vector<double> > > alldata;
  if(alldata.empty()){
    for(int period = 0; period < SIG_AND_BG*nperiod; period++){
      std::vector< std::vector<double> > vv;
      std::vector< double > emptyv;
      vv.push_back(emptyv); // 0->1
      for(int tb = 1; tb <= maxtb; tb++){
        std::vector<double> v;
        for(int eb = 0; eb < nenergybins; eb++)
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

    for(int periodsg = 0; periodsg < SIG_AND_BG*nperiod; periodsg++){
      // resist making cleverer
      const int beam = (periodsg%nperiod) < nperiodrhc? 0: 1;

      const bool isbg = periodsg >= nperiod;
      const int period = periodsg%nperiod;

      for(int eb = 0; eb < nenergybins; eb++){
        const int off_beam   = get_off_beam(beam, eb, isbg);
        const int off_period = get_off_period(period, eb);

        double model = 0;

        const double sca = scales[periodsg][eb];

        if(x <= holex_lo || x >= holex_hi)
          model += fabs(par[flat_nc+off_period]) * (isbg?bgmult:1);

        if(x >= holex_hi){
          const double norm_n   =
            beam == RHC_BEAM?
            par[nneut_nc+eb*SIG_AND_BG+isbg] /* R */:
            par[nneut_nc + (g_nbins_e+eb)*SIG_AND_BG + isbg]; /*F*/
          const double norm_b12 =
            beam == RHC_BEAM?
            par[nb12_nc+eb*SIG_AND_BG+isbg]:
            par[nb12_nc + (g_nbins_e+eb)*SIG_AND_BG+isbg];

          // Michel lifetime - ~average of mu+ and mu-, messed up by
          // detector effects
          model += fabs(par[nmich_nc+off_beam]*sca)/par[tmich_nc+off_beam]
                    * exp(-x/par[tmich_nc+off_beam]) +

                   fabs(norm_n  *sca) * nstuff
                + fabs(norm_b12*sca) * exp_x_b12life;
        }

        if((x >= -numi_len && x <= holex_lo) || (x >= holex_hi && x <= numi_len))
          model += fabs(par[pileup_nc+off_period])*(numi_len-fabs(x))
                   * (isbg?bgmult:1);

        const double data = alldata[periodsg][tb][eb];

        if(model > 0){
          like += model - data;
          if(data > 0) like += data * log(data/model);
        }
      }
    }
  }

  // penalty on Michel lifetimes.  Introduced when I started trying to
  // implement cuts on later events based on Michels, which reduces
  // the lever arm of Michels considerably.
  for(int i = 0; i < nbeam*g_nbins_e; i++)
    like += 0.5 * pow((par[tmich_nc + i] - tmich_nominal)/tmich_priorerr, 2);

  // penalty on neutron lifetime
  like += 0.5 *
    pow((par[tneut_nc] - n_lifetime_nominal)/n_lifetime_priorerr, 2);

  // penalty on neutron diffusion constant - factor of 2 on a log scale
  // (offerbangheap)
  like += 0.5 *
    pow((log(par[aneut_nc]) - log(n_diffusion_nominal))/M_LN2, 2);
}

void make_mn()
{
  std::vector<PAR> PARAMETERS = makeparameters();
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

static void set_ee_to_mn(const int periodsg, const int energybin) // 0-indexed
{
  const int beam = periodsg%nperiod < nperiodrhc? 0: 1;
  const bool isbg = periodsg >= nperiod;

  for(int i = 0; i < ncommonpar; i++)
    ee_neg->SetParameter(i, getpar(i));

  // walk forward through the blocks of parameters to the right place.
  for(int i = ncommonpar; i < ncommonpar+nper_e_beam_sg; i++){
    ee_neg->SetParameter(i,

      (i==ee_Nneut_par || i==ee_NB12_par || i==ee_NMich_par?
       scales[periodsg][energybin]: 1)* // per track -> total

      getpar(ncommonpar +
             (nbeam*g_nbins_e*(i-ncommonpar) +
             beam * g_nbins_e +
             energybin)*SIG_AND_BG + isbg
      )
    );
  }

  for(int i = ncommonpar+nper_e_beam_sg; i < npar_ee-1; i++){
    ee_neg->SetParameter(i, (isbg?bgmult:1)*
      getpar(ncommonpar +
        SIG_AND_BG*nbeam  *g_nbins_e*      nper_e_beam_sg          +
                   nperiod*g_nbins_e*(i-ncommonpar-nper_e_beam_sg) +
        g_nbins_e*(periodsg%nperiod) +
        energybin));
  }

  ee_neg->SetParameter(ee_scale_par, 1); // scale

  for(int f = 0; f < ntf1s; f++)
    for(int i = 0; i < npar_ee; i++)
      ees[f]->SetParameter(i, ee_neg->GetParameter(i));
}

// Set the drawing functions for the sum of all periods in a beam type
static void set_ee_to_mn_beam(const int beam, const int bin,
                              const bool isbg)
{
  // Oh so hacky!
  double pars[npar_ee];
  for(int i = 0; i < (beam == RHC_BEAM?nperiodrhc:nperiodfhc); i++){

    // Get drawing parameters for each period in this beam and sig/bg type.
    // Add them together as appropriate.
    const int periodsg = i + (beam == FHC_BEAM)*nperiodrhc + isbg*nperiod;
    set_ee_to_mn(periodsg, bin);
    if(i == 0)
      for(int p = 0; p < npar_ee; p++)
        pars[p] = fabs(ee_pos->GetParameter(p));
    else
      for(int p = 0; p < npar_ee; p++)
        // Only the additive parameters
        if(p == ee_Nneut_par || p == ee_NB12_par   || p == ee_NMich_par ||
           p == ee_flat_par  || p == ee_pileup_par ||
           p == ee_scale_par /* the scale?  really? */)
          pars[p] += fabs(ee_pos->GetParameter(p));
  }

  for(int f = 0; f < ntf1s; f++)
    for(int i = 0; i < npar_ee; i++)
      ees[f]->SetParameter(i, pars[i]);
}

// Get rid of any bins that have partial contents.  Just for drawing.
static void sanitize_after_rebinning(TH1D * x)
{
  for(int b = 0; b <= x->GetNbinsX(); b++){
    if(holex_lo > x->GetBinLowEdge(b) &&
       holex_lo < x->GetBinLowEdge(b+1)) x->SetBinContent(b, 0);
    if(holex_hi > x->GetBinLowEdge(b) &&
       holex_hi < x->GetBinLowEdge(b+1)) x->SetBinContent(b, 0);
    // does not consider the case that the bin is entirely
    // inside of the range, because then it will be zero anyway
  }
}

// Get the drawing scale for this beam and energy bin.  Keep the
// scale the same for signal and background histograms.
static double beam_scale(const int beam, const int bin)
{
  double ans = 0;
  for(int i = 0; i < (beam == RHC_BEAM?nperiodrhc:nperiodfhc); i++ ){
    const int idx = i + (beam == FHC_BEAM)*nperiodrhc;
    ans += scales[idx][bin];
  }
  return ans;
}

static void draw_ee_common(TH1D * x, const int rebin,
                           const char * const header,
                           const char * const outname,
                           const double ymin, const double ymax)
{
  x->Rebin(rebin);
  sanitize_after_rebinning(x);

  x->SetLineColor(kBlack);
  x->SetMarkerColor(kBlack);
  x->GetYaxis()->SetTitle(rebin== 1?"Delayed clusters/#mus":
                          Form("Delayed clusters/%d#kern[-0.5]{ }#mus", rebin));
  x->Draw("e");

  TLegend * leg = new TLegend(1-0.47*textsize/0.045*0.9, // small
                              1-0.40*textsize/0.045*0.9, // small
                              1-rightmargin, 1-topmargin-0.03);
  leg->AddEntry((TH1D*)NULL, header, "");
  for(int i = 0; i < ntf1s; i++){
    ees[i]->SetParameter(ee_scale_par, rebin);
    ees[i]->Draw("same");
    if(i != 0 /*ee_neg*/) leg->AddEntry(ees[i], ees_description[i], "l");
  }

  leg->SetTextSize(textsize*0.9); // small
  leg->SetMargin(0.165);
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  TGaxis ya(x->GetXaxis()->GetBinLowEdge(x->GetNbinsX()+1), ymin,
            x->GetXaxis()->GetBinLowEdge(x->GetNbinsX()+1), ymax,
            yminpertrack, ymaxpertrack, 510, "G+S");
  ya.CenterTitle();
  ya.SetLabelSize(textsize);
  ya.SetTitleSize(textsize);
  ya.SetLabelFont(42);
  ya.SetTitleFont(42);
  ya.SetLabelOffset(0.05);
  ya.SetTitleOffset(1.05);
  ya.SetTickSize(0.015); // requires "S" option in constructor
  ya.SetTitle("Delayed clusters/5#kern[-0.5]{ }#mus/track");
  ya.Draw();


  TText t(0, 0, "NOvA Preliminary");
  t.SetTextColor(kBlue);
  t.SetTextSize(textsize);
  t.SetTextFont(42);
  t.SetNDC();
  t.SetTextAlign(33);
  t.SetX(1-rightmargin-0.01);
  t.SetY(1-0.01);
  t.Draw();


  TGraph box;
  box.SetPoint(0, holex_lo, ymin);
  box.SetPoint(1, holex_lo, ymax);
  box.SetPoint(2, holex_hi, ymax);
  box.SetPoint(3, holex_hi, ymin);
  box.SetPoint(4, holex_lo, ymin);
  box.SetFillColorAlpha(kGray, 0.3);
  box.Draw("f");

  c1->SetTickx(1);
  c1->SetTicky(0);
  x->GetYaxis()->SetTickSize(0.015);
  x->GetXaxis()->SetTickSize(0.02);
  c1->SetLogy();
  c1->Print(outname);
}

static void stylehist(TH1 * h)
{
  h->SetLineWidth(2);
  h->GetXaxis()->SetLabelSize(textsize);
  h->GetYaxis()->SetLabelSize(textsize);
  h->GetYaxis()->SetTitleOffset(0.95 * 0.05/textsize);
  h->GetXaxis()->SetTitleSize(textsize);
  h->GetYaxis()->SetTitleSize(textsize);
  h->GetYaxis()->CenterTitle();
  h->GetXaxis()->CenterTitle();
}

// Draw the sum of the fit histograms for a beam type and energy bin.
static void draw_ee_beam(const int beam, const int energybin,
                         const bool isbg, const char * const outname,
                         const two_or_three_d cut_dimensions)
{
  c1->cd();

  set_ee_to_mn_beam(beam, energybin, isbg);

  // Get the form of the histogram and zero it
  TH1D * x = fithist[0]->ProjectionX("x", energybin+1, energybin+1 /* 0->1 */);
  x->Reset();

  for(int i = 0; i < (beam == RHC_BEAM?nperiodrhc:nperiodfhc); i++ ){
    const int idx = i + (beam == FHC_BEAM)*nperiodrhc;
    TH1D * tmp = fithist[idx+nperiod*isbg]->ProjectionX(
      "tmp", energybin+1, energybin+1);
    x->Add(tmp);
  }

  stylehist(x);

  const char * const header = Form("%s%s E_{#nu} %.1f#minus%.1f GeV",
    isbg == 0?"":"Pileup: ",
    beam == RHC_BEAM?"RHC":"FHC", bins_e[cut_dimensions][energybin],
                                  bins_e[cut_dimensions][energybin+1]);
  x->GetXaxis()->SetTitle("Time since muon stop (#mus)");
  const double ymin = beam_scale(beam, energybin)*yminpertrack*rebin,
               ymax = beam_scale(beam, energybin)*ymaxpertrack*rebin;
  x->GetYaxis()->SetRangeUser(ymin, ymax);

  draw_ee_common(x, rebin, header, outname, ymin, ymax);
}

// Draw the fit histogram for a period and energy bin. 0-indexed bins
static void draw_ee(const int periodsg, const int energybin,
                    const char * const outname,
                    const two_or_three_d cut_dimensions)
{
  c1->cd();

  set_ee_to_mn(periodsg, energybin);

  TH1D * x = fithist[periodsg]->ProjectionX("x", energybin+1, energybin+1);

  stylehist(x);

  const char * const header = Form("%s E_{#nu} %.1f-%.1f GeV",
    Lperiodnames[periodsg], bins_e[cut_dimensions][energybin],
                            bins_e[cut_dimensions][energybin+1]);
  x->GetXaxis()->SetTitle("Time since muon stop (#mus)");
  const double ymin = scales[periodsg][energybin]*yminpertrack*rebin,
               ymax = scales[periodsg][energybin]*ymaxpertrack*rebin;
  x->GetYaxis()->SetRangeUser(ymin, ymax);
  draw_ee_common(x, rebin, header, outname, ymin, ymax);
}

static vector< vector< vector<fitanswers> > > dothefit()
{
  maxfitt = maxrealtime; // modified later for sensitivity study

  // Start each fit with a clean MINUIT slate.  This is crucial!
  make_mn();

  int status = 0;

  for(int period = 0; period < nperiod; period++){
    for(int bin = 0; bin < g_nbins_e; bin++){
      const int off_period = get_off_period(period, bin);
      // May as well set this to near the right value
      mn->Command(Form("SET PAR %d %f", flat_nf+off_period,
        fithist[period]->ProjectionX("x", bin+1, bin+1)
          ->Integral(0, nnegbins-10)/(nnegbins - 10)));

      // Again, something reasonable based on the data.
      mn->Command(Form("SET PAR %d %f", pileup_nf+off_period,
        fithist[period]->ProjectionX("x", bin+1, bin+1)
          ->GetBinContent(nnegbins - 2)/8.));
    }
  }

  for(int beam = 0; beam < nbeam; beam++){
    for(int bin = 0; bin < g_nbins_e; bin++){
      for(int sigorbg = 0; sigorbg < SIG_AND_BG; sigorbg++){
        const int off_beam = get_off_beam(beam, bin, sigorbg);
        // Reasonable starting points
        mn->Command(Form("SET PAR %d %f", nmich_nf+off_beam, 0.8));
        mn->Command(Form("SET PAR %d %f",nneut_nf+off_beam,beam == RHC_BEAM?0.2:0.1));
        mn->Command(Form("SET PAR %d %f",nb12_nf+off_beam, beam == RHC_BEAM?0.2:0.01));

        // Start with the muon lifetime fixed so that it doesn't try to
        // swap with the neutron lifetime.
        mn->Command(Form("FIX %d", tmich_nf+off_beam));
      }
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

  const int migrad_tries = 4;

  printf("+++++ SIMPLEX followed by MIGRAD with some things fixed ++++++\n");
  status = mn->Command("SIMPLEX 100000");
  for(int i = 0; i < migrad_tries; i++)
    if(0 == (status = mn->Command("MIGRAD 100000")))
      break;

  // Now that we're (hopefully) converged, let muon lifetime float
  for(int beam = 0; beam < nbeam; beam++)
    for(int bin = 0; bin < g_nbins_e; bin++)
      for(int sigorbg = 0; sigorbg < SIG_AND_BG; sigorbg++)
        mn->Command(Form("REL %d", tmich_nf+get_off_beam(beam, bin, sigorbg)));

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
    for(int bin = 0; bin < g_nbins_e; bin++)
      for(int sigorbg = 0; sigorbg < SIG_AND_BG; sigorbg++)
        mn->Command(Form("SET LIM %d 1.6 2.6",
          tmich_nf+get_off_beam(beam, bin, sigorbg)));

  // MINOS errors are wrong if we used the abs trick, so limit
  for(int beam = 0; beam < nbeam; beam++){
    for(int bin = 0; bin < g_nbins_e; bin++){
      for(int sigorbg = 0; sigorbg < SIG_AND_BG; sigorbg++){
        const int off_beam = get_off_beam(beam, bin, sigorbg);
        mn->Command(Form("SET PAR %d %f", nb12_nf+off_beam,
               fabs(getpar( nb12_nc+off_beam))));
        mn->Command(Form("SET PAR %d %f", nneut_nf+off_beam,
              fabs(getpar(nneut_nc+off_beam))));

        // Let "B-12" go quite high because really it is probably mostly air
        // neutrons, especially for the pileup background samples.
        mn->Command(Form("SET LIM %d 0 100",  nb12_nf+off_beam));
        mn->Command(Form("SET LIM %d 0 10", nneut_nf+off_beam));
      }
    }
  }

  printf("+++++ MIGRAD with things freed ++++++\n");
  mn->Command("MIGRAD 100000");

  // Fix any nuisance parameters (Michels, uncorrelated background, pileup)
  // that are up against zero so that we can get a sensible error matrix out.
  // This is a little dubious, but so is what happens if we don't do it.
   for(int beam = 0; beam < nbeam; beam++){
    for(int bin = 0; bin < g_nbins_e; bin++){
      for(int sigorbg = 0; sigorbg < SIG_AND_BG; sigorbg++){
        const int off_beam = get_off_beam(beam, bin, sigorbg);
        if(fabs(getpar(nmich_nc + off_beam)) <
                geterr(nmich_nc + off_beam)){
          mn->Command(Form("FIX %d", nmich_nf + off_beam));
          mn->Command(Form("FIX %d", tmich_nf + off_beam));
        }
      }
    }
  }

  for(int period = 0; period < nperiod; period++){
    for(int bin = 0; bin < g_nbins_e; bin++){
      const int off_period = get_off_period(period, bin);
      if(fabs(getpar(flat_nc + off_period) <
              geterr(flat_nc + off_period)))
        mn->Command(Form("FIX %d", flat_nf + off_period));
      if(fabs(getpar(pileup_nc + off_period) <
              geterr(pileup_nc + off_period)))
        mn->Command(Form("FIX %d", pileup_nf + off_period));
    }
  }

  printf("+++++ MIGRAD and HESSE with zero pars fixed again ++++++\n");
  for(int i = 0; i < migrad_tries; i++)
    if(0 == (status = mn->Command("MIGRAD 100000")))
      break;

  status = mn->Command("HESSE");

  const bool DOMINOS = true, // slow and only marginally useful
             EVENDOB12MINOS = false;

  if(!status)
    for(int beam = 0; beam < nbeam; beam++)
      for(int bin = 0; bin < g_nbins_e; bin++)
        for(int sigorbg = 0; sigorbg < SIG_AND_BG; sigorbg++){
          const int off_beam = get_off_beam(beam, bin, sigorbg);
          if(DOMINOS){
            printf("+++++ MINOS on parameters of interest ++++++\n");
            gMinuit->Command(Form("MINOS 50000 %d", nneut_nf+off_beam));
            if(EVENDOB12MINOS)
              gMinuit->Command(Form("MINOS 50000 %d", nb12_nf+off_beam));
          }
        }

  gMinuit->Command("show min");

  vector< vector< vector<fitanswers> > > anses;

  for(int beam = 0; beam < nbeam; beam++){
    vector< vector<fitanswers> > ans2;
    for(int bin = 0; bin < g_nbins_e; bin++){
      vector<fitanswers> ans1;
      for(int sigorbg = 0; sigorbg < SIG_AND_BG; sigorbg++){
        const int off_beam = get_off_beam(beam, bin, sigorbg);
        fitanswers ans;
        ans.n_good    = onegoodminos(nneut_nc+off_beam, false);
        ans.n_mag     =       getpar(nneut_nc+off_beam);
        ans.n_mage_up = getbesterrup(nneut_nc+off_beam);
        ans.n_mage_dn = std::min(ans.n_mag, getbesterrdn(nneut_nc+off_beam));

        ans.b12_good   = onegoodminos(nb12_nc+off_beam, true);
        ans.b12mag     =       getpar(nb12_nc+off_beam);
        ans.b12mage_up = getbesterrup(nb12_nc+off_beam);
        ans.b12mage_dn = std::min(ans.b12mag, getbesterrdn(nb12_nc+off_beam));

        ans1.push_back(ans);
      }
      ans2.push_back(ans1);
    }
    anses.push_back(ans2);
  }

  if(DOMINOS){
    printf("+++++ MIGRAD & HESSE in case MINOS found a new minimum ++++++\n");
    status = mn->Command("MIGRAD");
    status = mn->Command("HESSE");
  }

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
        e->SetLineColorAlpha(kRed, 0.75);
        e->SetLineWidth(3);
        e->SetNpx(400);
        break;
      case 2:
        e->SetLineStyle(kDashed);
        e->SetLineColor(kBlack);
        e->SetLineWidth(2);
        e->SetNpx(400);
        break;
      case 3:
        e->SetLineStyle(kDashed);
        e->SetLineColor(kBlue);
        e->SetLineWidth(2);
        e->SetNpx(400);
        break;
      case 4:
        e->SetLineStyle(kDashed);
        e->SetLineColor(kViolet);
        e->SetLineWidth(3);
        break;
      case 5:
        e->SetLineStyle(kDashed);
        e->SetLineColor(kGreen+2);
        e->SetLineWidth(3);
        break;
      case 6:
        e->SetLineStyle(kDashed);
        e->SetNpx(400);
        e->SetLineColor(kOrange+2);
        e->SetLineWidth(2);
        break;
    }
    for(int p = 0; p < npar_ee; p++) e->SetParName(p, ee_parnames[p]);
  }
}

static void save_for_stage_two(TGraphAsymmErrors ** g_n_rhc,
                               TGraphAsymmErrors ** g_n_fhc,
                               TGraphAsymmErrors ** g_b12_rhc,
                               TGraphAsymmErrors ** g_b12_fhc,
                               const int mindist,
                               const float minslc, const float maxslc,
                               const two_or_three_d cut_dimensions)
{
  ofstream for_stage_two(Form("for_stage_two_mindist%d_nslc%.1f_%.1f_%s_%s.C",
                              mindist, minslc, maxslc,
                              muoncatcher?"muoncatcher":"main",
                              two_or_three_d_names[cut_dimensions]));
  for_stage_two << "{\n";
  for(int sigorbg = 0; sigorbg < SIG_AND_BG; sigorbg++){
    g_n_rhc  [sigorbg]->SetName(Form("raw_g_n_rhc[%d]", sigorbg));
    g_n_fhc  [sigorbg]->SetName(Form("raw_g_n_fhc[%d]", sigorbg));
    g_n_rhc  [sigorbg]->SavePrimitive(for_stage_two);
    g_n_fhc  [sigorbg]->SavePrimitive(for_stage_two);
    g_b12_rhc[sigorbg]->SetName(Form("raw_g_b12_rhc[%d]", sigorbg));
    g_b12_fhc[sigorbg]->SetName(Form("raw_g_b12_fhc[%d]", sigorbg));
    g_b12_rhc[sigorbg]->SavePrimitive(for_stage_two);
    g_b12_fhc[sigorbg]->SavePrimitive(for_stage_two);
  }

  // Following the example in TMinuitMinimizer::GetHessianMatrix()
  // XXX can we just invert the block of interest?  Does that make sense
  // mathematically? Because parameters with values consistent with zero ruin
  // everything.
  TMatrixDSym hessian(npar);
  mn->mnemat(hessian.GetMatrixArray(), npar);
  hessian.Invert();

  // Only care about the cross-terms for the 2*bins of the fit histograms
  for(int ii = 0; ii < g_nbins_e*nbeam*2 /* n and b12 */; ii++){
    for(int jj = ii; jj < g_nbins_e*nbeam*2; jj++){

      const int i = ii < g_nbins_e*nbeam? nneut_nc + ii:
               nb12_nc - g_nbins_e*nbeam + ii;
      const int j = jj < g_nbins_e*nbeam? nneut_nc + jj:
               nb12_nc - g_nbins_e*nbeam + jj;

      for_stage_two << "hessian[" << ii << "][" << jj << "] = "
        << hessian[i][j] << ";\n";
    }
    for_stage_two << "\n";
  }

  for_stage_two << "}\n";
}

// TODO: This is getting to be a lot of arguments, suggesting that passing it
// in via a struct, or using a config file would be better.  Probably can't
// use a struct because there's no way to pass that in from the command line.
void rhc_stage_one(const char * const savedhistfile, const int mindist,
                   const float minslc, const float maxslc, const string region,
                   const two_or_three_d cut_dimensions)
{
  if(mindist < 0) return; // to compile only

  gStyle->SetOptStat(0);
  gStyle->SetFrameLineWidth(2);

  muoncatcher = region == "muoncatcher";

  // Needs to get into fcn(), so has to be global
  g_nbins_e = nbins_e[cut_dimensions];

  // From external Monte Carlo
  n_lifetime_priorerr = muoncatcher?10:5.;
  // XXX Hackily provide information from the more data-rich fits to the
  // data-poor fits so they don't spin out of control.
  if(cut_dimensions == THREED || mindist <= 2) n_lifetime_priorerr = 2.;

  n_lifetime_nominal = muoncatcher?nominal_neutron_lifetime_muoncatcher
                                  :nominal_neutron_lifetime_main;

  /*
    Based on my MC, 400us is reasonable for a mindist of 6, but more like
    25us for a mindist of 2. However, it's murky because we don't measure
    the capture point, but rather the position where the gammas compton,
    and that usually only in 2D. Assuming 2D, I get more like 80us. And
    of course the functional form isn't quite right for 2D...

    Anyway, here's a stupid parameterization that skewers the 2D and 3D
    cases. Note the addition of 5/8 of a cellwidth because that's how
    much you get for mindist == 0 on average. In fcn() this is given a
    factor of 2 error (offerbangheap). Mostly we just measure it.

    But don't start it smaller than the neutron lifetime, because
    then it competes with decay for explaining the shape of the data
    and the fit gets stuck in bad states.
  */
  n_diffusion_nominal = std::max(cut_dimensions == TWOD?
                        14.94 * pow(mindist + 0.625, 1.74):
                         1.39 * pow(mindist + 0.625, 3.00),
                        n_lifetime_nominal);

  gROOT->Macro(savedhistfile);
  for(int i = 0; i < nperiod * SIG_AND_BG; i++){
    const char * namt = Form("%s_tcounts", Speriodnames[i]);
    if(NULL == (all_tcounts[i] = dynamic_cast<TH1D*>(
      gROOT->FindObject(namt)))){
      fprintf(stderr, "Could not read %s from %s\n", namt, savedhistfile);
      _exit(1);
    }

    const char * nam2 = Form("%s", Speriodnames[i]);
    if(NULL == (fithist[i] = dynamic_cast<TH2D*>(
      gROOT->FindObject(nam2)))){
      fprintf(stderr, "Could not read %s from %s\n", nam2, savedhistfile);
      _exit(1);
    }
  }

  init_ee();

  TGraphAsymmErrors * g_n_rhc[SIG_AND_BG],
                    * g_n_fhc[SIG_AND_BG],
                    * g_b12_rhc[SIG_AND_BG],
                    * g_b12_fhc[SIG_AND_BG];

  for(int i = 0; i < SIG_AND_BG; i++){
    g_n_rhc[i] = new TGraphAsymmErrors;
    g_n_fhc[i] = new TGraphAsymmErrors;
    g_b12_rhc[i] = new TGraphAsymmErrors;
    g_b12_fhc[i] = new TGraphAsymmErrors;
  }

  gErrorIgnoreLevel = kError; // after new TFile and Get above

  scales.resize(nperiod*SIG_AND_BG);

  for(int s = 1; s <= g_nbins_e; s++)
    for(int p = 0; p < nperiod*SIG_AND_BG; p++)
      scales[p].push_back(all_tcounts[p]->GetBinContent(s));

  const vector< vector< vector<fitanswers> > > anses = dothefit();

  for(int s = 0; s < g_nbins_e; s++){
    for(int sigorbg = 0; sigorbg < SIG_AND_BG; sigorbg++){
      const fitanswers rhc_ans = anses[0][s][sigorbg];
      const fitanswers fhc_ans = anses[1][s][sigorbg];

      const double loslce = fithist[0]->GetYaxis()->GetBinLowEdge(s+1);
      const double hislce = fithist[0]->GetYaxis()->GetBinLowEdge(s+2);

      const double graph_x  = (loslce+hislce)/2;
      const double graph_xe = (hislce-loslce)/2;

      addpoint(g_n_rhc[sigorbg], graph_x,
        rhc_ans.n_mag, graph_xe, rhc_ans.n_mage_dn, rhc_ans.n_mage_up);

      addpoint(g_n_fhc[sigorbg], graph_x,
        fhc_ans.n_mag, graph_xe, fhc_ans.n_mage_dn, fhc_ans.n_mage_up);

      addpoint(g_b12_rhc[sigorbg], graph_x,
        rhc_ans.b12mag, graph_xe, rhc_ans.b12mage_dn, rhc_ans.b12mage_up);

      addpoint(g_b12_fhc[sigorbg], graph_x,
        fhc_ans.b12mag, graph_xe, fhc_ans.b12mage_dn, fhc_ans.b12mage_up);

      printf("%s RHC neutron (%5.3f-%5.3f)GeV %s: %.5f +%.5f -%.5f\n",
        rhc_ans.n_good?"MINOS ":"MIGRAD", loslce, hislce,
        sigorbg?"SIG":"BG ",
        rhc_ans.n_mag, rhc_ans.n_mage_up, rhc_ans.n_mage_dn);

      printf("%s FHC neutron (%5.3f-%5.3f)GeV %s: %.5f +%.5f -%.5f\n",
        fhc_ans.n_good?"MINOS ":"MIGRAD", loslce, hislce,
        sigorbg?"SIG":"BG ",
        fhc_ans.n_mag, fhc_ans.n_mage_up, fhc_ans.n_mage_dn);
    }
  }

  const string filename = Form("fit_stage_one_mindist%d_nslc%.1f_%.1f_%s_%s.pdf",
                               mindist, minslc, maxslc, region.c_str(),
                               two_or_three_d_names[cut_dimensions]);

  c1->SetTopMargin   (topmargin);
  c1->SetBottomMargin(bottommargin);
  c1->SetLeftMargin  (leftmargin);
  c1->SetRightMargin (rightmargin);

  c1->Print(Form("%s[", filename.c_str()));

  for(int beam = 0; beam < nbeam; beam++)
    for(int i = 0; i < g_nbins_e; i++)
      for(int sigorbg = 0; sigorbg < SIG_AND_BG; sigorbg++)
        draw_ee_beam(beam, i, sigorbg, filename.c_str(), cut_dimensions);

  // Print out a per-period set of histograms unless we have merged all RHC and
  // FHC periods together, in which case that would be redundant.
  if(nbeam < nperiod)
    for(int period = 0; period < nperiod; period++)
      for(int i = 0; i < g_nbins_e; i++)
        for(int sigorbg = 0; sigorbg < SIG_AND_BG; sigorbg++)
          draw_ee(period + sigorbg*nperiod, i, filename.c_str(), cut_dimensions);

  c1->Print(Form("%s]", filename.c_str()));

  save_for_stage_two(g_n_rhc, g_n_fhc, g_b12_rhc, g_b12_fhc,
                     mindist, minslc, maxslc, cut_dimensions);
}
