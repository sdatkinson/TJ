
#include <fcntl.h>
#include <unistd.h>
#include <cstdlib>
#include <assert.h>
#include <gsl/gsl_rng.h>

#include "rng_class.h"

rng_class::rng_class() :
  RNGtype(gsl_rng_mt19937)
{
  gsl_rng_env_setup();
  RNG = gsl_rng_alloc(RNGtype);
}



rng_class::~rng_class()
{
  if (RNG) {
    gsl_rng_free(RNG);
    RNG=NULL;
  }
}

void rng_class::seed()
{
  int randomData = open("/dev/random", O_RDONLY);
  unsigned long int seed=0;
  int result = read(randomData, &seed, (sizeof seed));
  assert(result>=0);
  close(randomData);

  gsl_rng_set(RNG,seed);
}



double rng_class::rand_double()
{
  return gsl_rng_uniform(RNG);
}



euclidean_vector<DIM, double> rng_class::rand_sphere()
// Create a random vector within the unit sphere
// IMPORTANT: norm is not necessarily unity
{
  euclidean_vector<DIM,double> v;
  bool have_vector=false;
  while (!have_vector) {
    double v_norm_sq = 0.0;
    for (int i=0; i<DIM; i++) {
      const double val = 2.0*(rand_double()-1.0);
      v[i] = val;
      v_norm_sq += val*val;
    }
    have_vector = v_norm_sq<=1.0;
  }

  return v;
}
