
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cassert>

#include "parameters_class.h"


//Constructor for th
parameters_class::parameters_class()
{
  set_all_undefined();
}

parameters_class::parameters_class(std::string param_file)
{
  set_all_undefined();
  read_parameters_file(param_file);
}



void parameters_class::read_parameters_file(std::string param_file)
{
  // variables
  std::string   input;
  std::string   token;
  std::ifstream fin;

  // read input file directly, line by line
  fin.open(param_file.c_str(),std::ios::in);

  if (!fin.is_open()) {
    std::cerr << "ERROR: Failed to open parameters file at " << param_file << std::endl;
  }
  else {
    while (!fin.eof()) {//Go through the lines of the input file
      getline(fin,input);

      if (!input.empty()) {
	std::stringstream inputStream(input);

	if (input[0] == '#') {//This line specifies a parameter...
	  inputStream >> token;

	  //========================
	  // bool parameters

	  if (token == "#use_nnl") {
	    inputStream >> bool_parameter  [use_nnl];
	    bool_parameter_defined         [use_nnl] = true;
	  }
	  else if (token == "#random_moves") {
	    inputStream >> bool_parameter  [random_moves];
	    bool_parameter_defined         [random_moves] = true;
	  }

	  //========================
	  // int parameters

	  else if (token == "#max_iterations") {
	    inputStream >> int_parameter   [max_iterations];
	    int_parameter_defined          [max_iterations] = true;
	  }
	  else if (token == "#overboxes") {
	    inputStream >> int_parameter   [overboxes];
	    int_parameter_defined          [overboxes] = true;
	  }
	  else if (token == "#max_iterations") {
	    inputStream >> int_parameter   [max_iterations];
	    int_parameter_defined          [max_iterations] = true;
	  }
	  else if (token == "#print_every") {
	    inputStream >> int_parameter   [print_every];
	    int_parameter_defined          [print_every] = true;
	  }
	  else if (token == "#print_precision") {
	    inputStream >> int_parameter   [print_precision];
	    int_parameter_defined          [print_precision] = true;
	  }
	  else if (token == "#glpk_presolve_var") {
	    inputStream >> int_parameter   [glpk_presolve_var];
	    int_parameter_defined          [glpk_presolve_var] = true;
	  }
	  else if (token == "#glpk_basis_fact_type") {
	    inputStream >> int_parameter   [glpk_basis_fact_type];
	    int_parameter_defined          [glpk_basis_fact_type] = true;
	  }

	  //========================
	  // double parameters

	  else if (token == "#trans_max") {
	    inputStream >> double_parameter[trans_max];
	    double_parameter_defined       [trans_max] = true;
	  }
	  else if (token == "#comp_max") {
	    inputStream >> double_parameter[comp_max];
	    double_parameter_defined       [comp_max] = true;
	  }
	  else if (token == "#shear_max") {
	    inputStream >> double_parameter[shear_max];
	    double_parameter_defined       [shear_max] = true;
	  }
	  else if (token == "#influence_sphere") {
	    inputStream >> double_parameter[influence_sphere];
	    double_parameter_defined       [influence_sphere] = true;
	  }
	  else if (token == "#term_tol") {
	    inputStream >> double_parameter[term_tol];
	    double_parameter_defined       [term_tol] = true;
	  }
	  else if (token == "#resize_tol") {
	    inputStream >> double_parameter[resize_tol];
	    double_parameter_defined       [resize_tol] = true;
	  }
	  else if (token == "#resize_space") {
	    inputStream >> double_parameter[resize_space];
	    double_parameter_defined       [resize_space] = true;
	  }
	  else if (token == "#random_move_size") {
	    inputStream >> double_parameter[random_move_size];
	    double_parameter_defined       [random_move_size] = true;
	  }
	  else if (token == "#nnl_extra_distance") {
	    inputStream >> double_parameter[nnl_extra_distance];
	    double_parameter_defined       [nnl_extra_distance] = true;
	  }
	  else if (token == "#feasible_tol") {
	    inputStream >> double_parameter[feasible_tol];
	    double_parameter_defined       [feasible_tol] = true;
	  }
	  else if (token == "#glpk_pivot_tol") {
	    inputStream >> double_parameter[glpk_pivot_tol];
	    double_parameter_defined       [glpk_pivot_tol] = true;
	  }

	  //========================
	  // (end of parameter list)

	}
      } 
    }
  }

  assert(all_parameters_defined());

  if (int_parameter[overboxes] > 0 && bool_parameter[use_nnl]) {
    std::cerr << "You cannot use a nearest neighbor list (#useNNL = true) and have an overbox check (#overBoxes > 0). Setting use_nnl=false\n";
    bool_parameter[use_nnl]=false;
  }

  if (int_parameter[overboxes] == 0 && !bool_parameter[use_nnl]) {
    std::cerr << "Warning: checking only inside the unit cell for overlap. NOT checking periodic boundary conditions!\n";
  }

  fin.close();
}



void parameters_class::rescale_double_parameters(double scale)
// rescale all pertinent variables, e.g., delta, etc., by the 
// diameter.
//
// This also rescales a bunch of not-useful stuff (e.g. 
// feasible_tol, which has more to do with the precision of the 
// LP solver/machine precision), but it's currently up to the 
// programmer to use the right one (scaled or not).
{
  for (unsigned int i=0; i<num_double; i++)
    double_parameter_scaled[i] = scale*double_parameter[i];

  std::cout << "influence_sphere = " << double_parameter_scaled[parameters_class::influence_sphere]   << "\n"
	    << "nnl_extra_dist   = " << double_parameter_scaled[parameters_class::nnl_extra_distance] << "\n"
	    << "trans_max        = " << double_parameter_scaled[parameters_class::trans_max]          << "\n"
	    << "comp_max         = " << double_parameter_scaled[parameters_class::comp_max]           << "\n"
	    << "shear_max        = " << double_parameter_scaled[parameters_class::shear_max]          << "\n";
}



//====================================================================
// Private methods



void parameters_class::set_all_undefined()
{
  for (unsigned int i=0; i<num_bool  ; i++)
    bool_parameter_defined[i]   = false;

  for (unsigned int i=0; i<num_int   ; i++)
    int_parameter_defined[i]    = false;

  for (unsigned int i=0; i<num_double; i++)
    double_parameter_defined[i] = false;
}

bool parameters_class::all_parameters_defined()
{
  bool all_defined=true;
  for (unsigned int i=0; i<num_bool   && all_defined; i++) {
    if (!bool_parameter_defined[i]) {
      std::cerr<<"bool_parameter["<<i<<"] undefined\n";
      all_defined = false;
    }
  }

  for (unsigned int i=0; i<num_int    && all_defined; i++) {
    if (!int_parameter_defined[i]) {
      std::cerr<<"int_parameter["<<i<<"] undefined\n";
      all_defined = false;
    }
  }

  for (unsigned int i=0; i<num_double && all_defined; i++) {
    if (!double_parameter_defined[i]) {
      std::cerr<<"double_parameter["<<i<<"] undefined\n";
      all_defined = false;
    }
  }

  return all_defined;
}
