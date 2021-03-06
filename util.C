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
