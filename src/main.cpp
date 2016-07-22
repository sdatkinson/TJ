//=======================================================================
//
// Torquato-Jiao Sequential Linear Programming Sphere Packing Algorithm
//
// Author:  Steven Atkinson
// contact: sdatkins@princeton.edu
//
//=======================================================================
//
//
// See README for details on the execution of this program
//
//
// Example execution:
// ./tj_3d iconfig.dat parameters.txt output

#include "macro.h"

#include "input_class.h"   // For parsing the arguments passed to the program
#include "packing_class.h" // Structure of the packing class
#include "slp_class.h"     // Class for the SLP algorithm
#include "rng_class.h"     // Wrapper for GSL's RNG



rng_class rng;



int main (int    argc,
	  char **argv)
{
  std::cout.precision(16);

  //Parse inputs:
  input_class input(argc,argv);

  if (!input.get_valid()) {
    input.print_usage();
    return EXIT_ERR_INVALID_INPUT;
  }

  //RNG initialization:
  rng.seed();

  //=======================================================================
  //  Read in the initial configuration

  packing_class packing(input.get_str(input_class::input_file));

  //=======================================================================
  //  Read in the running parameters

  parameters_class parameters(input.get_str(input_class::param_file));
  parameters.rescale_double_parameters(packing.get_diameter());

  //=======================================================================
  //  Compress the packing using the SLP algorithm

  slp_class slp(packing,parameters,input.get_str(input_class::output_file));

  //=======================================================================
  //  Write the results:

  packing.write(input.get_str(input_class::output_file)+".dat");

  return EXIT_SUCCESS;
}
