

#include "const.h"
#include "euclidean_vector.h"
#include "euclidean_matrix.h"
#include "rng_class.h"
#include "sphere_class.h"



extern rng_class rng;



sphere_class::sphere_class()
{
  index=0;
  nnl_radius=0.0;
  global_displacement_since_last_nnl_update = 0.0;
  position_set_since_last_nnl_update = false;
}



sphere_class::~sphere_class()
{
}



void sphere_class::set_local_position(const euclidean_vector<DIM,double>& local_position_i)
{
  local_position = local_position_i;
  position_set_since_last_nnl_update = true;
}



void sphere_class::set_global_position(const euclidean_vector<DIM,double>& global_position_i)
{
  global_position = global_position_i;
  position_set_since_last_nnl_update = true;
}



euclidean_vector<DIM,double> sphere_class::get_local_position()
// Might benefit from something to check if it's up to date.
// But the current implementation was coded well and works correctly
{
  return local_position;
}




euclidean_vector<DIM,double> sphere_class::get_global_position()
// Might benefit from something to check if it's up to date.
// But the current implementation was coded well and works correctly
{
  return global_position;
}



void sphere_class::set_local_position_from_global(const euclidean_matrix<DIM,double>& invLambda)
{
  local_position = invLambda*global_position;
}



void sphere_class::set_global_position_from_local(const euclidean_matrix<DIM,double>&    lambda)
{
  global_position = lambda*local_position;
}



void sphere_class::keep_inside_fundamental_cell()
{
  local_position = (local_position + 1.0) % 1.0;
}



void sphere_class::displace_local (const euclidean_vector<DIM,double>& disp,
				   const euclidean_matrix<DIM,double>& lambda)
{
  displace_local_delay_nnl_update(disp);
  update_global_displacement_since_last_nnl_update(euclidean_vector<DIM,double>::norm(lambda*disp));
}



void sphere_class::displace_local_delay_nnl_update(const euclidean_vector<DIM,double>& disp)
{
  local_position += disp;
}



void sphere_class::displace_global(const euclidean_vector<DIM,double>& disp)
{
  global_position += disp;
  update_global_displacement_since_last_nnl_update(disp.norm());
}



void sphere_class::update_global_displacement_since_last_nnl_update(const double disp)
{
  global_displacement_since_last_nnl_update += disp;
}



void sphere_class::random_global_displacement(const double move_size)
{
  euclidean_vector<DIM,double> disp = rng.rand_sphere();
  disp*=move_size;
  displace_global(disp);
}



bool sphere_class::needs_nnl_update()
{
  //  if (force_nnl)
  //    return true;
  if (position_set_since_last_nnl_update)
    return true;
  if (global_displacement_since_last_nnl_update > nnl_radius)
    return true;

  return false;
}



void sphere_class::clear_nnls()
{
  nnl_influence.clear();
  nnl_overlap.clear();
}



void sphere_class::set_nnls(const std::vector<nnl_entry_periodic_class>& nnl_influence_i,
			    const std::vector<nnl_entry_periodic_class>& nnl_overlap_i)
{
  nnl_influence.clear();
  nnl_overlap.clear();

  nnl_influence = nnl_influence_i;
  nnl_overlap   = nnl_overlap_i;

  global_displacement_since_last_nnl_update = 0.0;
  position_set_since_last_nnl_update = false;
}




