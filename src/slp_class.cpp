
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

#include <gsl/gsl_sf.h>

#include "const.h"
#include "packing_class.h"
#include "parameters_class.h"
#include "lp_class.h"

#include "slp_class.h"



slp_class::slp_class(packing_class&    packing,
		     parameters_class& parameters,
		     std::string       output_file) :
  vol_history_length(3),
  vol_history(vol_history_length,packing.get_vol())
{
  pack(packing,parameters,output_file);
}



void slp_class::pack(packing_class&    packing,
		     parameters_class& parameters,
		     std::string       output_file)
// Solve linear programs and incease the packing
// fraction until a termination criterion is met
{
  //Initialize the NNL, if requested
  if (parameters.bool_parameter[parameters_class::use_nnl])
    packing.update_nnls(parameters.double_parameter_scaled[parameters_class::influence_sphere],
			parameters.double_parameter_scaled[parameters_class::nnl_extra_distance]);

  //Resizing the initial configuration (if needed)
  packing.resize_if_overlap(parameters,false);

  //Find the total volume of the spheres:
  const double diameter   = packing.get_diameter();
  const double sphere_vol = double(packing.get_n()) * std::pow(PI,double(DIM)/2.0) / gsl_sf_gamma(DIM/2.0 +1.0) * std::pow(0.5*diameter,DIM);

  //The SLP iterations
  lp_iteration = 0;
  while (get_status(parameters)==slp_class::status_running) {
    lp_iteration++;

    // Run the LP routines:
    slp_step_class slp_step(packing,parameters);
    slp_step.run_once();

    vol_history[lp_iteration%vol_history_length] = packing.get_vol();

    const double phi = sphere_vol / packing.get_vol();
    std::cout << lp_iteration << "  " << phi << std::endl;

    if (parameters.int_parameter[parameters_class::print_every] > 0 && 
	lp_iteration % parameters.int_parameter[parameters_class::print_every] == 0) {
      // Print intermediate configurations
      std::stringstream ss_it; ss_it<<lp_iteration;
      std::string       out_name = output_file + "_" + ss_it.str() + ".dat";
      packing.write(out_name);
    }
  }//LP loop

  print_status(get_status(parameters));
}



int slp_class::get_status(parameters_class& parameters)
{
  // If we've exceeded the max iterations
  if (lp_iteration>=parameters.int_parameter[parameters_class::max_iterations])
    return slp_class::status_term_iterations;

  // If the packing volume has converged
  if (lp_iteration>=vol_history_length) {
    const int end_idx = lp_iteration  % vol_history_length;
    const int start_idx = (end_idx+1) % vol_history_length;
    const double vol_change = std::abs(vol_history[end_idx] - vol_history[start_idx]);
    if (vol_change <= parameters.double_parameter[parameters_class::term_tol])
      return slp_class::status_term_converged;
  }

  // No termination criteria met
  return slp_class::status_running;
}



void slp_class::print_status(int status)
{
  switch (status) {
  case slp_class::status_running :
    std::cout<<"SLP is still running\n";
    break;
  case slp_class::status_term_iterations :
    std::cout<<"SLP terminated due to max iterations exceeded\n";
    break;
  case slp_class::status_term_converged :
    std::cout<<"SLP terminated due to volume convergence\n";
    break;
  default :
    std::cout<<"Unknown SLP status\n";
    break;
  }
}
