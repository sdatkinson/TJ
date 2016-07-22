//LPClass.cpp
//Implementation of the Linear Programming class
//Included from main.cpp

#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>

#include "const.h"
#include "euclidean_vector.h"
#include "euclidean_matrix.h"
#include "packing_class.h"
#include "sphere_class.h"
#include "parameters_class.h"
#include "nnl_class.h"
#include "glpk_solver_class.h"
#include "lp_class.h"



std::ostream& operator<<(std::ostream&               os,
			 const lp_constraint_struct& constr)
{
  //assert here that both vectors are the same length if you want to

  const int save_precision = os.precision();
  os.precision(3); //Could have the constraint have its own preferred precision...

  std::vector<int>::const_iterator    idx_iterator;
  std::vector<double>::const_iterator val_iterator;

  //LHS
  for (idx_iterator=constr.idx.begin(), val_iterator=constr.val.begin(); idx_iterator!=std::prev(constr.idx.end()) && val_iterator!=std::prev(constr.val.end()); ++idx_iterator, ++val_iterator)
    os<<(*val_iterator)<<" x"<<(*idx_iterator)<<" + ";
  os<<(*val_iterator)<<" x"<<(*idx_iterator);

  //Equality
  switch (constr.equality) {
  case lp_constraint_struct::less_than :
    os<<" < ";
    break;
  case lp_constraint_struct::equal_to :
    os<<" = ";
    break;
  case lp_constraint_struct::greater_than :
    os<<" > ";
    break;
  }

  //RHS
  os<<constr.rhs;
  os.precision(save_precision);

  return os;
}



std::ostream& operator<<(std::ostream&                   os,
			 const objective_function_class& obj)
{
  //assert here that both vectors are the same length if you want to

  std::vector<int>::const_iterator    idx_iterator;
  std::vector<double>::const_iterator val_iterator;

  for (idx_iterator=obj.idx.begin(), val_iterator=obj.val.begin(); idx_iterator!=std::prev(obj.idx.end()) && val_iterator!=std::prev(obj.val.end()); ++idx_iterator, ++val_iterator)
    os<<(*val_iterator)<<" x"<<(*idx_iterator)<<" + ";
  os<<(*val_iterator)<<" x"<<(*idx_iterator);

  return os;
}



lp_solver_class::lp_solver_class(parameters_class* parameters) :
 num_variables(0),
 solved(false)
{
  import_basic_parameters(parameters);
}



void lp_solver_class::import_basic_parameters(parameters_class* parameters)
{
  trans_max = parameters->double_parameter_scaled[parameters_class::trans_max];
  comp_max  = parameters->double_parameter_scaled[parameters_class::comp_max];
  shear_max = parameters->double_parameter_scaled[parameters_class::shear_max];
}



slp_step_class::slp_step_class(packing_class&    packing_i,
			       parameters_class& parameters_i) :
  num_epsilon(DIM*(DIM+1)/2)
{
  set_packing   (&packing_i);
  set_parameters(&parameters_i);

  get_objective_function();

  // Implement other solvers here
  lp_solver = new glpk_solver_class(parameters);
}



slp_step_class::~slp_step_class()
{
  if (lp_solver) {
    delete lp_solver;
    lp_solver=NULL;
  }
}



void slp_step_class::run_once()
{
  get_constraints();
  lp_solver->import_problem(num_epsilon+DIM*packing->get_n(),lp_constraint,objective_function);
  std::vector<double> solution = lp_solver->get_solution();
  apply(solution);
}



//======================================================================
// GENERATE THE LP



void slp_step_class::get_objective_function()
// Trace of the strain tensor
{
  int epsilon=0;
  for (int i=0; i<DIM; i++) {
    objective_function.idx.push_back(epsilon);
    objective_function.val.push_back(1.0);
    epsilon+=DIM-i;
  }

  //std::cout<<"Objective function is "<<objective_function<<std::endl;
}



void slp_step_class::get_constraints()
{
  //First, clear old constraints:
  lp_constraint.clear();
  
  //Now, get the new constraints:
  get_bulk_constraint();
  get_pair_constraints();
}



void slp_step_class::get_bulk_constraint()
// Tr(E) <= 0
{
  lp_constraint_struct bulk_constraint;
  //bulk_constraint.equality = lp_constraint_struct::less_than;
  bulk_constraint.equality = lp_constraint_struct::greater_than;
  bulk_constraint.rhs = 0.0;

  bulk_constraint.idx.reserve(DIM);
  bulk_constraint.val.reserve(DIM);

  unsigned int epsilon = 0;
  for (unsigned int i=0; i<DIM; i++) {
    bulk_constraint.idx.push_back(epsilon);
    //bulk_constraint.val.push_back(1.0);
    bulk_constraint.val.push_back(-1.0);
    epsilon += DIM-i;
  }

  lp_constraint.push_back(bulk_constraint);
}



void slp_step_class::get_pair_constraints()
// Look through the pairs of spheres and
// append constraints to the list as
// appropriate
//
// 2 different binary cases:
// 1: (A) NNL or (B) no NNL
// 2: (A) overboxes=0 (L/2 search) or (B) overboxes>1
//    Note: NNL was originally implemented with no overboxes,
//    and it still needs to be checked to ensure proper operation
{
  //=========
  // NNLs

  //Note: could create an object/class for making a list of particles to check against...

  const int overboxes = parameters->int_parameter[parameters_class::overboxes];

  //Compute the constraints between pairs of spheres:
  for (unsigned int i=0; i<packing->sphere.size(); i++) {
    sphere_class* sphere1=&(packing->sphere[i]);
    //for (std::vector<nnl_entry_periodic_class>::iterator nnl_entry=sphere1->nnl_influence.begin(); nnl_entry!=sphere1->nnl_influence.end(); ++nnl_entry)
    for (std::vector<nnl_entry_periodic_class>::iterator nnl_entry=std::prev(sphere1->nnl_influence.end()); nnl_entry!=std::prev(sphere1->nnl_influence.begin()); --nnl_entry)
      //for (nnl_entry_periodic_class& nnl_entry : sphere1->nnl_influence)
      check_for_constraint(sphere1,
			   nnl_entry->get_target(),
			   packing->lambda,
			   overboxes);
  }

  //========
  // no NNLs

  /*
  //Re-do using iterators.
  for (unsigned int i=0;     i<packing->get_n(); i++)
    for (unsigned int j=i+1; j<packing->get_n(); j++)
      check_for_constraint(&(packing->sphere[i]),
		   	   &(packing->sphere[j]),
			   packing->lambda,
			   overboxes);
  */

  //  std::cout<<"Generated "<<lp_constraint.size()<<" constraints."<<std::endl;
}



void slp_step_class::check_for_constraint(sphere_class*                 sphere_i,
					  sphere_class*                 sphere_j,
					  euclidean_matrix<DIM,double>& lambda,
					  int                           overboxes)
{
  if (sphere_i->get_index() > sphere_j->get_index())
    return;

  euclidean_vector<DIM,double> r_ij = sphere_j->get_local_position() - sphere_i->get_local_position();
  if (overboxes==0) {//Nearest image (L/2)
    r_ij.get_minimum_image();
    const double dist = euclidean_vector<DIM,double>::norm(lambda*r_ij);
    if (dist < packing->get_diameter() + parameters->double_parameter_scaled[parameters_class::influence_sphere])
      append_pair_constraint(r_ij,sphere_i->get_index(),sphere_j->get_index(),packing->get_diameter(),packing->lambda);
  }
  else {// overboxes
    throw("overboxes>0 not implemented right now.\n");
  }
}



void slp_step_class::append_pair_constraint(euclidean_vector<DIM,double>&       r_ij_local,
					    const unsigned int                  i,
					    const unsigned int                  j,
					    const double                        diameter,
					    const euclidean_matrix<DIM,double>& lambda)
{
  lp_constraint_struct pair_constraint;
  pair_constraint.idx.reserve(num_epsilon+2*DIM);
  pair_constraint.val.reserve(num_epsilon+2*DIM);

  euclidean_vector<DIM,double> r_ij_global = lambda*r_ij_local;

  // input the epsilon coefficients
  unsigned int epsilon_index = 0; //index of the strain entry
  for (unsigned int d1=0; d1<DIM; d1++) {//First index of the epsilon
    pair_constraint.idx.push_back(epsilon_index++);
    pair_constraint.val.push_back(pow(r_ij_global[d1],2));
    for (unsigned int d2=d1+1; d2<DIM; d2++) {//Second index of the epsilon
      pair_constraint.idx.push_back(epsilon_index++);
      pair_constraint.val.push_back(2.0*r_ij_global[d1]*r_ij_global[d2]);
    }
  }

  // calculate gram_local = G*r_ij_local (G is gram matrix = lambda^T*lambda, so we're computing Lambda^T*r_ij_global)
  euclidean_vector<DIM,double> gram_local = lambda.get_transpose()*r_ij_global; // Note the right multiplication by lambda

  // input the sphere translation coefficients (order is: x1,x2,y1,y2,z1,z2...)
  for (unsigned int d=0; d<DIM; d++) {
    //Sphere i:
    pair_constraint.idx.push_back(num_epsilon + i*DIM + d);
    pair_constraint.val.push_back(-gram_local[d]);
    //Sphere j:
    pair_constraint.idx.push_back(num_epsilon + j*DIM + d);
    pair_constraint.val.push_back( gram_local[d]);
  }

  // Get the right hand side and equality
  pair_constraint.equality = lp_constraint_struct::greater_than;
  pair_constraint.rhs = 0.5*(pow(diameter,2) - gram_local.dot(r_ij_local));

  lp_constraint.push_back(pair_constraint);
}



//======================================================================
// APPLY THE LP SOLVE



void slp_step_class::apply(std::vector<double>& solution)
//Update the packing with the results of the LP solve
{
  euclidean_matrix<DIM,double>               new_lambda   = get_new_lambda   (solution);
  std::vector<euclidean_vector<DIM,double> > displacement = get_displacements(solution);
  apply_changes(displacement,new_lambda); //Applies the solution and checks for overlap

  // Random moves
  if (parameters->bool_parameter[parameters_class::random_moves])
	packing->random_moves(parameters->double_parameter_scaled[parameters_class::random_move_size]);

  // Do an NNL update if needed:
  packing->update_nnls(parameters->double_parameter_scaled[parameters_class::influence_sphere  ],
  		  	  	  	   parameters->double_parameter_scaled[parameters_class::nnl_extra_distance]);

  // 5) now check for overlap
  packing->resize_if_overlap(*parameters,true); // use_space=true
}



euclidean_matrix<DIM,double> slp_step_class::get_new_lambda(std::vector<double>& solution)
{
  //Get the strain matrix:
  euclidean_matrix<DIM,double> strain;

  {
    unsigned int epsilon_idx = 0;
    for   (unsigned int i=0; i<DIM; i++) {
      for (unsigned int j=i; j<DIM; j++) {
    	const double val = solution[epsilon_idx++];
    	if (i==j)
    	  strain[i][i] = 1.0+val;
    	else {
    	  strain[i][j] = val;
    	  strain[j][i] = val;
    	}
      }
    }
  }

  return strain*packing->lambda;
}



std::vector<euclidean_vector<DIM,double> > slp_step_class::get_displacements(std::vector<double>& solution)
{
  std::vector<euclidean_vector<DIM,double> > displacement(packing->get_n());
  for (unsigned int i=0; i<displacement.size(); i++)
    for (unsigned int j=0; j<DIM; j++)
      displacement[i][j] = solution[num_epsilon + i*DIM + j];

  return displacement;
}



void slp_step_class::apply_changes(std::vector<euclidean_vector<DIM,double> >& displacement,
				   euclidean_matrix<DIM,double>&               new_lambda)
{
  // 1) New lambda:
  packing->lambda = new_lambda;
  euclidean_matrix<DIM,double> inv_lambda=packing->lambda;
  inv_lambda.invert();

  // 2) Displacements:
  for (unsigned int i=0; i<packing->get_n(); i++ ) {
    // Keep trcak of where we started
    euclidean_vector<DIM,double> pos_g0 = packing->sphere[i].get_global_position();
    // Move local position
    packing->sphere[i].displace_local_delay_nnl_update(displacement[i]); // Delay NNL distance update
    // Update global using new fundamental cell
    packing->sphere[i].set_global_position_from_local(packing->lambda);
    // Keep track of how far we moved for nnl update
    packing->sphere[i].update_global_displacement_since_last_nnl_update(euclidean_vector<DIM,double>::norm(packing->sphere[i].get_global_position()-pos_g0));
  }

  // 3) Keep inside fundamental cell!
  packing->keep_spheres_inside_fundamental_cell();
}
