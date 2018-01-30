double gaus(const double x, const double sigma)
{
  const double invsig = 1/sigma;
  // different sigmas are relatively normalized
  return exp(-0.5*x*x*invsig*invsig)*invsig;
}

void normalize(vector<double> & h)
{
  double totp = 0;
  for(unsigned int i = 0; i < h.size(); i++) totp += h[i];
  for(unsigned int i = 0; i < h.size(); i++) h[i] /= totp;
}

void convolve_syst(vector<double> & prob, const double interval,
                   const double up, const double dn)
{
  vector<double> result = prob;
  for(unsigned int i = 0; i < prob.size(); i++){
    double newp = 0;
    for(unsigned int j = 0; j < prob.size(); j++){
      const double delta = ((int)j-(int)i)*interval;
      if(prob[j] > 0)
        // Almost sure the sign is right
        newp += prob[j] * gaus(delta, delta > 0? dn: up);
    } 
    result[i] = newp;
  } 

  prob = result;

  normalize(prob);
}

void valerr(const vector<double> & combined, const double interval)
{
  const double CL = 0.682689492137;

  double best = 0;
  int besti = -1;
 
  for(unsigned int i = 0; i < combined.size(); i++){
    if(combined[i] > best){
      best = combined[i];
      besti = i;
    }
  }

  vector<double> combinedbyp = combined;
  std::sort   (combinedbyp.begin(), combinedbyp.end());
  std::reverse(combinedbyp.begin(), combinedbyp.end());

  double acc = 0;
  double cutoff = -1;
  for(unsigned int i = 0; i < combinedbyp.size(); i++){
    acc += combinedbyp[i];
    if(acc >= CL){
      cutoff = (combinedbyp[i] + combinedbyp[(i>0?i:1)-1])/2;
      break;
    }
  }

  if(cutoff == -1) fprintf(stderr, "Failed to find cutoff\n");

  double dn = 0, up = 0;

  for(unsigned int i = 0; i < combined.size(); i++){
    if(combined[i] > cutoff){
      dn = interval*(besti - (i-0.5));
      break;
    }
  }

  for(unsigned int i = besti; i < combined.size(); i++){
    if(combined[i] < cutoff){
      up = ((i-0.5) - besti)*interval;
      break;
    }
  }

  printf("%f + %f - %f\n", besti*interval, up, dn);
}

