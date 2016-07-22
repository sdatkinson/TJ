#ifndef PARAMETERS_CLASS_H
#define PARAMETERS_CLASS_H

#include <string>

class parameters_class {
 public:
  parameters_class();
  parameters_class(std::string param_file);
  //~parameters_class();

  void read_parameters_file(std::string param_file);
  void rescale_double_parameters(double scale);

  enum e_bool {
    use_nnl=0,
    random_moves,
    num_bool
  };

  enum e_int {
    max_iterations=0,
    overboxes,
    print_every,
    print_precision,
    glpk_presolve_var,
    glpk_basis_fact_type,
    num_int
  };

  enum e_double {
    trans_max=0,
    comp_max,
    shear_max,
    influence_sphere,
    term_tol,
    resize_tol,
    resize_space,
    random_move_size,
    nnl_extra_distance,
    feasible_tol,
    glpk_pivot_tol,
    num_double
  };

  bool   bool_parameter          [num_bool];
  bool   bool_parameter_defined  [num_bool];
  int    int_parameter           [num_int];
  bool   int_parameter_defined   [num_int];
  double double_parameter        [num_double];
  double double_parameter_scaled [num_double];
  bool   double_parameter_defined[num_double];

 private:
  void set_all_undefined();
  bool all_parameters_defined();
};

#endif
