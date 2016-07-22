#ifndef GLPK_SOLVER_CLASS_H
#define GLPK_SOLVER_CLASS_H

#include <vector>
#include <glpk.h>
#include "lp_class.h"


// Implemented in glpk_solver_class.cpp:
class glpk_solver_class : public lp_solver_class {
 public:
  // Variables:
  glp_prob* lp;
  glp_smcp  solver_params;
  glp_bfcp  other_params;

  // Methods:
  glpk_solver_class(parameters_class* parameters);
  ~glpk_solver_class();

 protected:
  void                import_solver_parameters(parameters_class*                  parameters);
  void                import_problem(unsigned int                       num_variables_i,
				     std::vector<lp_constraint_struct>& lp_constraint,
				     objective_function_class&          objective_function);

  void                solve();
  std::vector<double> get_solution();

 private:
  // Private variables:
  objective_function_class objective_function;

  unsigned int num_constraints; // Number of rows in the constraint matrix
  unsigned int num_entries;     // Number of nonzero entries in LHS...
  unsigned int array_size_lhs;  // Allocation size for LHS
  unsigned int array_size_rhs;  // Allocation size for RHS

  int*    ia;
  int*    ja;
  double* ar;
  int*    ineq;
  double* bounds;

  //Solver numerics parameters:
  int    presolve_var;
  int    basis_fact_type;
  double feasible_tol;
  double pivot_tol;

  //Private methods:
  void set_lhs_size();
  void set_rhs_size();
  void import_constraints(std::vector<lp_constraint_struct>& lp_constraint);
};



#endif
