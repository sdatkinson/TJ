//slp_class.h
//Structure of the Sequential Linear Programming class

#ifndef SLP_CLASS_H
#define SLP_CLASS_H

#include <string>

#include "packing_class.h"
#include "parameters_class.h"

class slp_class {
 public:

  //===========================================
  // Public member classes

  slp_class(packing_class&    packing,
	        parameters_class& parameters,
	        std::string       output_file);
  //~slp_class();

  void pack(packing_class&    packing,
	        parameters_class& parameters,
	        std::string       ouptut_file);

 private:

  // Enumeration of the possible states of the SLP routine
  enum e_status {
    status_running,         // Still running
    status_term_iterations, // Terminated because the max number of iterations have been done
    status_term_converged   // Terminated because the packing's volume has converged
  };

  //===========================================
  // Private members

  int lp_iteration;
  const int  vol_history_length;
  std::vector<double> vol_history;

  //===========================================
  // Private member classes
   
  int    get_status(parameters_class& parameters);
  void   print_status(int status);
};
//static LPClass LP;

#endif
