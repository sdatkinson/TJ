#ifndef RNG_CLASS_H
#define RNG_CLASS_H

#include <gsl/gsl_rng.h>

#include "euclidean_vector.h"

class rng_class {
 public:
  rng_class();
  ~rng_class();

  void   seed();
  double rand_double();
  euclidean_vector<DIM,double> rand_sphere();

 private:
  const gsl_rng_type *RNGtype;
  gsl_rng *RNG;
};

#endif
