

#ifndef SPHERE_CLASS_H
#define SPHERE_CLASS_H

#include <vector>

#include "euclidean_vector.h"
#include "euclidean_matrix.h"

#include "nnl_class.h"

class sphere_class {
 public:
  sphere_class();
  ~sphere_class();

  std::vector<nnl_entry_periodic_class> nnl_influence; // Pairs that might cause SLP pair constraints (inside this's influence sphere)
  std::vector<nnl_entry_periodic_class> nnl_overlap;   // Pairs that might physically overlap
  double                                nnl_radius;    // This is a little clunky.  It'd be better to make a class for the nnls instead of using a vector of entries so that the radius can be part of the class...

  void         set_index(unsigned int index_i) {index=index_i;};
  unsigned int get_index() {return index;};

  void set_local_position (const euclidean_vector<DIM,double>& local_position_i);
  void set_global_position(const euclidean_vector<DIM,double>& global_position_i);

  euclidean_vector<DIM,double> get_local_position();
  euclidean_vector<DIM,double> get_global_position();

  void displace_local (const euclidean_vector<DIM,double>& disp,
		       const euclidean_matrix<DIM,double>& lambda);
  void displace_local_delay_nnl_update(const euclidean_vector<DIM,double>& disp);
  void displace_global(const euclidean_vector<DIM,double>& disp);

  void update_global_displacement_since_last_nnl_update(const double disp);

  void random_global_displacement(const double global_move_size);

  // These don't trip position_set_since_last_nnl_update to be true
  void set_local_position_from_global(const euclidean_matrix<DIM,double>& inv_lambda);
  void set_global_position_from_local(const euclidean_matrix<DIM,double>&     lambda);

  // make sure that the local position is inside [0,1)^DIM
  // Operates on local position only--need to use set_global_position_from_local() afterwards!
  void keep_inside_fundamental_cell();

  //void force_nnl_update(); //Not needed
  bool needs_nnl_update();
  void clear_nnls();
  void set_nnls(const std::vector<nnl_entry_periodic_class>& nnl_influence_i,
		const std::vector<nnl_entry_periodic_class>& nnl_overlap_i);

 private:
  euclidean_vector<DIM,double> local_position;
  euclidean_vector<DIM,double> global_position;
  //Add booleans for if the positions are up to date...
  unsigned int                 index;
  double                       global_displacement_since_last_nnl_update; // For Keeping track of if an nnl update is needed!
  bool                         position_set_since_last_nnl_update;
};


#endif
