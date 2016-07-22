//LPClass.h
//Structure of the Linear Programming class

#ifndef LP_CLASS_H
#define LP_CLASS_H

#include <iostream>
#include <vector>
#include <glpk.h>

#include "const.h"
#include "euclidean_vector.h"
#include "euclidean_matrix.h"
#include "parameters_class.h"
#include "packing_class.h"
#include "sphere_class.h"



//======================================================================
// Bits and pieces that make up an LP



// Constraint
struct lp_constraint_struct {
  enum e_equality {
    less_than=0,
    equal_to,
    greater_than
  };

  std::vector<int>    idx; // Variables entering the equation
  std::vector<double> val; // coefficients for those variables
  int equality;
  double rhs;
};
std::ostream& operator<<(std::ostream&               os,
			 const lp_constraint_struct& constr);



// Objective function expression
class objective_function_class {
 public:
  std::vector<int>    idx; // Index of a variable with a nonzero objective coefficient
  std::vector<double> val; // The coefficient

  objective_function_class& operator=(const objective_function_class& src) {
    idx = src.idx;
    val = src.val;
    return *this;
  };
};
std::ostream& operator<<(std::ostream&                   os,
			 const objective_function_class& obj);



//======================================================================
// LP solvers



// Abstract solver class
// For GLPK, see glpk_solver_class.h,cpp
class lp_solver_class {
 public:
  lp_solver_class(parameters_class* parameters);
  virtual ~lp_solver_class() {};

  // Virtual functions that outline the lp solver's use
  void                        import_basic_parameters(parameters_class* parameters);
  virtual void                import_problem(unsigned int                       num_variables_i,
					     std::vector<lp_constraint_struct>& lp_constraint,
					     objective_function_class&          objective_function) = 0;
  virtual void                solve()        = 0;
  virtual std::vector<double> get_solution() = 0;

 protected:
  unsigned int num_variables; // Number of columns in the constraint matrix
  bool         solved;        // Used so that we don't solve unless we need to

  double       trans_max;
  double       comp_max;
  double       shear_max;

  void         set_unsolved() {solved=false;};
  void         set_solved()   {solved=true;};
  bool         is_solved()    {return solved;};
};



//======================================================================
// Class for an SLP step



class slp_step_class {
 public:

  //===========================================
  // Public member functions

  slp_step_class(packing_class& packing,
		  	  	 parameters_class& parameters);
  ~slp_step_class();

  void set_packing   (packing_class*    packing_i)    {packing   =packing_i;};
  void set_parameters(parameters_class* parameters_i) {parameters=parameters_i;};

  //Calls the whole process of generating, solving, and applying
  void run_once();    

 private:

  //===========================================
  // Private members

  //Pointers to the packing and parameters that the slp step is working with
  packing_class*    packing;
  parameters_class* parameters;
  
  //The number of strain variables 
  const int num_epsilon;

  std::vector<lp_constraint_struct> lp_constraint;
  objective_function_class          objective_function;
  lp_solver_class* lp_solver;

  //===========================================
  // Private member functions

  //Generating constraints:
  void get_objective_function();
  void get_constraints();
  //I don't like the naming here: this adds the constraint too by calling append_pair_constraint if there's one to be had
  void check_for_constraint(sphere_class*                 sphere_i,
			    sphere_class*                 sphere_j,
			    euclidean_matrix<DIM,double>& lambda,
			    int                           overboxes);
  void append_pair_constraint(euclidean_vector<DIM,double>&       r_ij_local,
			      const unsigned int                  i,
			      const unsigned int                  j,
			      const double                        diameter,
			      const euclidean_matrix<DIM,double>& lambda);
  void get_bulk_constraint(); 
  void get_pair_constraints();

  //applying the solution to the LP:
  void                                       apply            (std::vector<double>& solution);
  euclidean_matrix<DIM,double>               get_new_lambda   (std::vector<double>& solution);
  std::vector<euclidean_vector<DIM,double> > get_displacements(std::vector<double>& solution);
  void                                       apply_changes    (std::vector<euclidean_vector<DIM,double> >& displacements,
							       euclidean_matrix<DIM,double>&               new_lambda);
  void                                       apply_strain     ();
  void                                       apply_movements  ();
};




#endif
