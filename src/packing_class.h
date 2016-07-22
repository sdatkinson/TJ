#ifndef PACKING_CLASS_H
#define PACKING_CLASS_H

#include <string>
#include <vector>
#include <fstream>
#include <cmath>

#include "euclidean_vector.h"
#include "euclidean_matrix.h"
#include "sphere_class.h"
#include "parameters_class.h"

class packing_class {
 public:

  //=================================================================
  // Public members:

  // Note: make lambda into a class to make sure that things can be safely lazy-updated
  // currently relying on careful programming
  // Lattice generator matrix (columns span the fundamental cell)
  euclidean_matrix<DIM,double> lambda;
  // Its inverse:
  euclidean_matrix<DIM,double> inv_lambda;
  //bool inv_lambda_up_to_date;

  // The spheres in the packing:
  std::vector<sphere_class> sphere;


  //=================================================================
  // Public member functions:

  packing_class();
  packing_class(std::string input_file);
  ~packing_class();

  void read (std::string input_file);
  void write(std::string output_file);

  // Get simple properties:
  unsigned int get_n()       {return sphere.size();};
  double       get_diameter(){return diameter;};
  double       get_vol()     {return std::abs(lambda.determinant());};
  std::string  get_src()     {return src_file;};

  void update_nnls(const double       influence_sphere,
		   const double       nnl_extra_distance);

  void random_moves(const double move_size);

  void keep_spheres_inside_fundamental_cell();

  void resize_if_overlap(parameters_class& parameters,
			 bool              use_space);

 private:

  //=================================================================
  // Private members:


  std::string src_file; //The file from which the packing was read
  double      diameter;

  //=================================================================
  // Private member functions:

  void read_header          (std::ifstream& fin);
  void read_lattice         (std::ifstream& fin);
  void read_global_positions(std::ifstream& fin);
  void set_local_positions();

  void write_header          (std::ofstream& fout);
  void write_lattice         (std::ofstream& fout);
  void write_global_positions(std::ofstream& fout);

  bool nnls_need_update();

  // Adds entries to sphere i's NNL list 
  // (AND THE NNL LIST OF ANYTHING INTERACTING WITH IT!)
  void populate_nnl(const unsigned int i,
		    const double       influence_sphere,
		    const double       nnl_extra_distance);
};

#endif
