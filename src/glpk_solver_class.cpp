//Implementation of the GLPK solver class


#include <iostream>
#include <vector>
#include <iterator>
//#include <cassert>

#include <glpk.h>

#include "parameters_class.h"
#include "glpk_solver_class.h"



glpk_solver_class::glpk_solver_class(parameters_class* parameters) :
  lp_solver_class(parameters)
{
  lp     = NULL;
  ia     = NULL;
  ja     = NULL;
  ar     = NULL;
  ineq   = NULL;
  bounds = NULL;

  num_entries     = 0;
  num_constraints = 0;
  array_size_lhs  = 0;
  array_size_rhs  = 0;
  solved          = false;
  import_solver_parameters(parameters);
}



glpk_solver_class::~glpk_solver_class()
{
  if (lp) {
    glp_delete_prob(lp);
    glp_free_env();
    lp = NULL;
  }
  if (ia) {
    delete [] ia;
    ia = NULL;
  }
  if (ja) {
    delete [] ja;
    ja = NULL;
  }
  if (ar) {
    delete [] ar;
    ar = NULL;
  }
  if (ineq) {
    delete [] ineq;
    ineq = NULL;
  }
  if (bounds) {
    delete [] bounds;
    bounds = NULL;
  }
}



void glpk_solver_class::import_solver_parameters(parameters_class* parameters)
{
  presolve_var    = parameters->int_parameter   [parameters_class::glpk_presolve_var];
  basis_fact_type = parameters->int_parameter   [parameters_class::glpk_basis_fact_type];
  feasible_tol    = parameters->double_parameter[parameters_class::feasible_tol];
  pivot_tol       = parameters->double_parameter[parameters_class::glpk_pivot_tol];
}



void glpk_solver_class::import_problem(unsigned int                       num_variables_i,
				       std::vector<lp_constraint_struct>& lp_constraint,
				       objective_function_class&          objective_function_i)
{
  //Set number of variables:
  num_variables=num_variables_i;

  //Count the number of entries:
  num_entries=0;
  for (std::vector<lp_constraint_struct>::iterator it=lp_constraint.begin(); it!=lp_constraint.end(); ++it)
    num_entries += it->idx.size();

  num_constraints=lp_constraint.size();

  // Reallocate the LHS arrays if needed:
  set_lhs_size();

  //Reallocate the RHS arrays if needed:
  set_rhs_size();

  //Fill in the arrays:
  import_constraints(lp_constraint);

  //Copy the objective function since we need to wait until the "solve" step to use it.
  objective_function = objective_function_i;

  set_unsolved();
}



void glpk_solver_class::set_lhs_size()
{
  if (num_entries > array_size_lhs || !ia) {
    if (ia)
      delete [] ia;
    ia = new int[num_entries+1]; //"+1" because glpk doesn't use entry [0]
    if (ja)
      delete [] ja;
    ja = new int[num_entries+1];
    if (ar)
      delete [] ar;
    ar = new double[num_entries+1];
    array_size_lhs = num_entries;
  }
}



void glpk_solver_class::set_rhs_size()
// Sets variables and takes care of allocation
{
  if (num_constraints > array_size_rhs || !ineq) {
    if (ineq)
      delete [] ineq;
    ineq = new int[num_constraints+1]; //"+1" because glpk doesn't use entry [0]
    if (bounds)
      delete [] bounds;
    bounds = new double[num_constraints+1];
    array_size_rhs = num_constraints;
  }
}



void glpk_solver_class::import_constraints(std::vector<lp_constraint_struct>& lp_constraint)
{
  unsigned int lhs_index=0;
  for (unsigned int rhs_index=0; rhs_index < num_constraints; rhs_index++) {
    bounds[rhs_index+1] = lp_constraint[rhs_index].rhs;
    switch (lp_constraint[rhs_index].equality) {
    case lp_constraint_struct::less_than :
      ineq  [rhs_index+1] = GLP_UP;
      break;
    case lp_constraint_struct::greater_than :
      ineq  [rhs_index+1] = GLP_LO;
      break;
    case lp_constraint_struct::equal_to :
      ineq  [rhs_index+1] = GLP_FX;
      break;
    default :
      throw("Unexpected equality condition\n");
    }

    for (unsigned int entry_index=0; entry_index<lp_constraint[rhs_index].idx.size(); entry_index++) {
      // GLPK matrix indices start at 1, so we need to add 1 to ia and ja entries:
      ia[lhs_index    +1] = rhs_index+1;
      ja[lhs_index    +1] = lp_constraint[rhs_index].idx[entry_index] + 1;
      ar[(lhs_index++)+1] = lp_constraint[rhs_index].val[entry_index];
    }
  }
  //assert(lhs_index==num_entries);
  //std::cout<<"Imported "<<lp_constraint.size()<<" constraints and "<<num_entries<<" entries."<<std::endl;
}



void glpk_solver_class::solve()
{
  if (is_solved())
    return;

  // std::cout<<"glpk_solver.solve()"<<std::endl;

  if (!lp) {
    // create problem and read values in
    lp = glp_create_prob();
    glp_set_prob_name(lp,"SLP");
    glp_set_obj_dir(lp,GLP_MIN);

    // specify control parameters
    glp_init_smcp(&solver_params);
    solver_params.presolve = presolve_var;
    solver_params.tol_bnd  = feasible_tol;
    //solver_params.tol_piv  = pivot_tol;
    solver_params.msg_lev  = GLP_MSG_ERR; //Only print error and warning messages to the terminal

    glp_get_bfcp(lp,&other_params);
    other_params.type = basis_fact_type;
    glp_set_bfcp(lp,&other_params);
  }


  // first do the objective function and structural VARIABLES (i.e., position displacement and epsilons)
  glp_add_cols(lp,num_variables);
  // Set objective function
  // Remember that GLPK indices start at 1, not 0
  for (unsigned int i=0; i<objective_function.idx.size(); i++)
    glp_set_obj_coef(lp,objective_function.idx[i]+1,objective_function.val[i]);

  // Set variable bounds:
  {
    // Strains:
    unsigned int idx=0;
    for (int i=0; i<DIM; i++) {
      for (int j=i; j<DIM; j++) {
    	if (i==j)
    	  glp_set_col_bnds(lp,(idx++)+1,GLP_DB,-comp_max,comp_max);
    	else
    	  glp_set_col_bnds(lp,(idx++)+1,GLP_DB,-shear_max,shear_max);
      }
    }

    // Displacements:
    while (idx<num_variables)
      glp_set_col_bnds(lp,(idx++)+1,GLP_DB,-trans_max,trans_max);
  }

  // Constraint bounds
  glp_add_rows(lp,num_constraints);
  for (unsigned int i=0; i<num_constraints; i++) {
    switch (ineq[i+1]) {
    case GLP_LO :
      //glp_set_row_bnds(lp,i+1,GLP_LO,0.0,bounds[i+1]);
      glp_set_row_bnds(lp,i+1,GLP_LO,bounds[i+1],0.0);
      break;
    case GLP_UP :
      //glp_set_row_bnds(lp,i+1,GLP_UP,bounds[i+1],0.0);
      glp_set_row_bnds(lp,i+1,GLP_UP,0.0,bounds[i+1]);
      break;
    case GLP_FX :
      glp_set_row_bnds(lp,i+1,GLP_FX,bounds[i+1],bounds[i+1]);
      break;
    default :
      std::cout<<"Unexpected inequality in GLPK constraint."<<std::endl;
      throw;
    }
  }

  glp_load_matrix(lp,num_entries,ia,ja,ar);

  // Solve:
  const int solver_status = glp_simplex(lp,&solver_params);
  if (solver_status != 0) 
    std::cout<<"Warning: GLPK solver terminated with status "<<solver_status<<std::endl;

  set_solved();
}



std::vector<double> glpk_solver_class::get_solution()
{
  if (!is_solved())
    solve();

  //Allocate:
  std::vector<double> solution;
  solution.reserve(num_variables);

  //Export:
  for (unsigned int i=0; i<num_variables; i++)
    solution.push_back(glp_get_col_prim(lp,i+1)); //START AT 1!

  return solution;
}
