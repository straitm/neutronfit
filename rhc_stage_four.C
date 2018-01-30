#include <algorithm>
#include <fstream>
#include <stdio.h>
#include <vector>
using std::vector;

#include "bayes.C"

int rhc_stage_four(const char * const infile)
{
  double scale, prob;
  vector<double> probs[2];
  double oldscale = -1e100;

  int probi = 0;

  double interval = 0;

  std::ifstream finfile(infile);
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
    fprintf(stderr, "Not enough series\n");
    return 1;
  }

  if(probs[0].size() != probs[1].size()){
    fprintf(stderr, "Series don't have the same size!\n");
    return 1;
  }

  vector<double> combined;

  for(unsigned int i = 0; i < probs[0].size(); i++)
    combined.push_back(probs[0][i]*probs[1][i]);

  normalize(combined);

  for(unsigned int i = 0; i < combined.size(); i++)
    printf("%10f %10f %10f -> %g\n", interval*i, probs[0][i], probs[1][i], combined[i]);

  valerr(combined, interval);

  return 0;
}
