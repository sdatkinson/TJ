
#include <string>
#include <vector>
#include <algorithm>

#include "macro.h"
#include "euclidean_vector.h"
#include "euclidean_matrix.h"
#include "packing_class.h"


packing_class::packing_class()
{
  diameter=0.0;
}



packing_class::packing_class(std::string input_file)
{
  diameter=0.0;
  this->read(input_file);
}



packing_class::~packing_class()
{
}



void packing_class::read(std::string input_file)
{
  src_file = input_file;
  std::ifstream fin;
  try {
    fin.open(input_file.c_str(),std::ios::in);
    if (!fin.is_open()) {
      std::cerr << "Can't open initial configuration at " << input_file << std::endl;
      throw;
    }
    read_header(fin);
    read_lattice(fin);
    read_global_positions(fin);
    fin.close();
    set_local_positions();

  }
  catch (...) {
    //If the read was a failure, then print an error and wipe the packing
    std::cerr<<"ERROR: failed to read the packing"<<std::endl;

    if (fin.is_open())
      fin.close();

    sphere.clear();
    lambda=0.0;
  }
}



void packing_class::write(std::string output_file)
{
  std::ofstream fout;
  try {
    fout.open(output_file.c_str(),std::ios::out);

    if (!fout.is_open()) {
      std::cerr << "ERROR: Failed to open output file at " << output_file << std::endl;
      throw;
    }
    fout.precision(FOUT_PRECISION);
    write_header(fout);
    write_lattice(fout);
    write_global_positions(fout);
    fout.close();
    std::cout<<"Wrote "<<output_file<<std::endl;
  }
  catch (...) {
    std::cerr<<"Some error while writing packing to file"<<std::endl;
    if (fout.is_open())
      fout.close();
  }
}



void packing_class::update_nnls(const double       influence_sphere,
				const double       nnl_extra_distance)
{
  if (!nnls_need_update())
    return;

  std::cout<<"update nnls...\n";

  // clear old nnls:
  std::for_each(sphere.begin(),sphere.end(),[](sphere_class& s){s.clear_nnls();});

  // calculate neighbor lists
  // Note: implement cell grid approach here for increased speed
  unsigned int influence_entries=0;
  unsigned int overlap_entries=0;
  for (unsigned int i=0; i<get_n(); i++) {
    populate_nnl(i,influence_sphere,nnl_extra_distance);
    influence_entries += sphere[i].nnl_influence.size();
    overlap_entries   += sphere[i].nnl_overlap.size();
  }
  std::cout<<"  influence_entries : "<<influence_entries
	   <<"  overlap_entries   : "<<overlap_entries<<"\n";
}



void packing_class::random_moves(const double move_size)
{
  inv_lambda=lambda;
  inv_lambda.invert();
  for (sphere_class& s : sphere) {
    s.random_global_displacement(move_size);
    s.set_local_position_from_global(inv_lambda);
  }

  keep_spheres_inside_fundamental_cell();
}



void packing_class::keep_spheres_inside_fundamental_cell()
{
  for (sphere_class& s : sphere) {
    s.keep_inside_fundamental_cell();
    s.set_global_position_from_local(lambda);
  }
}



void packing_class::resize_if_overlap(parameters_class& parameters,
				      bool              use_space)
{
  double max_resize  = 1.0; // 1.0 means keep the packing the same size; >1.0 means that lambda grows
  double max_overlap = 0.0; // 0.0 means contact, <0 means free space, >0 means there's overlap
  bool   first_comparison=true;

  const double resize_tol   = parameters.double_parameter[parameters_class::resize_tol];
  const double resize_space = parameters.double_parameter[parameters_class::resize_space];

  if (parameters.bool_parameter[parameters_class::use_nnl]) {
    // use the L/2 method if an NNL is in use
    for (unsigned int i=0; i<get_n(); i++) {
      for (unsigned int idx=0; idx<sphere[i].nnl_overlap.size(); idx++) {
	const unsigned int j = sphere[i].nnl_overlap[idx].get_target()->get_index();
	if (i<j) {
	  //get overlap
	  euclidean_vector<DIM,double> r_ij(sphere[j].get_local_position() - sphere[i].get_local_position());
	  r_ij.get_minimum_image();
	  const double dist = euclidean_vector<DIM,double>::norm(lambda*r_ij);
	  const double overlap = diameter - (dist + resize_tol);
	  if (first_comparison) {
		first_comparison=false;
	    const double resize = diameter/dist;
		max_resize  = resize;
		max_overlap = overlap;
	  }
	  if (overlap > 0) {
	    const double resize = diameter/dist;
	    if (resize > max_resize)
	      max_resize = resize;
	    if (overlap > max_overlap)
	      max_overlap = overlap;
	  }
	}
      }
    }
  }//end of "use nnls"
  else {//no NNL
    throw("No nnls not implemented right now.\n");
  }

  // rescale the lambdas matrix lattice vectors by maxResize so that no detected overlap is present
  if (max_overlap > 0.0) {
    // if no extra space is desired (e.g. last iteration), then simply rescale so that farthest overlap is now contacting
    if (use_space)
      max_resize *= (1.0 + resize_space);

    // rescale the lambdas matrix
    lambda *= max_resize;
    // Translate the spheres' global positions:
    for (sphere_class& s : sphere)
      s.set_global_position_from_local(lambda);
  }
}







//=======================================================================
//  Private methods:



void packing_class::read_header(std::ifstream& fin)
{
  unsigned int dim_test;
  fin >> dim_test;
  assert(dim_test==DIM);
  double numSpheres;
  fin >> numSpheres;
  sphere.resize(numSpheres);
  {
    unsigned int i=0;
    for (std::vector<sphere_class>::iterator it=sphere.begin(); it!=sphere.end(); ++it)
      it->set_index(i++);
  }
  fin >> diameter;
}



void packing_class::read_lattice(std::ifstream& fin)
{
  for (unsigned int i=0; i<DIM; i++)
    for (unsigned int j=0; j<DIM; j++)
      fin >> lambda[i][j];
}



void packing_class::read_global_positions(std::ifstream& fin)
{
  for (unsigned int i=0; i<sphere.size(); i++) {
    euclidean_vector<DIM,double> pos;
    for (int d=0; d<DIM; d++)
      fin >> pos[d];
    sphere[i].set_global_position(pos);
  }
}



void packing_class::set_local_positions()
{
  euclidean_matrix<DIM,double> inv_lambda = lambda;
  inv_lambda.invert();

  std::cout<<"lambda = \n"    <<lambda<<"\n"
	   <<"inv_lambda = \n"<<inv_lambda<<"\n";
  std::for_each(sphere.begin(),sphere.end(),[inv_lambda](sphere_class& s) {s.set_local_position_from_global(inv_lambda);});
}



void packing_class::write_header(std::ofstream& fout)
{
  fout << DIM << "\n";
  fout << get_n() << "\n";
  fout << diameter << "\n";
}



void packing_class::write_lattice(std::ofstream& fout)
{
  for (unsigned int i=0;i<DIM;i++) {
    for (unsigned int j=0; j<DIM; j++)
      fout << lambda[j][i] << "\t";
    fout << "\n";
  }
}



void packing_class::write_global_positions(std::ofstream& fout)
{
  for (unsigned int i=0; i<get_n(); i++) {
    euclidean_vector<DIM,double> pos_g = sphere[i].get_global_position();
    for (unsigned int d=0; d<DIM; d++) {
      fout << pos_g[d];
      if (d < DIM-1)
	fout << "\t";
    }
    if (i<get_n()-1)
      fout << "\n";
  }
}



bool packing_class::nnls_need_update()
{
  bool need_update=false;
  for (std::vector<sphere_class>::iterator it=sphere.begin(); it!=sphere.end() && !need_update; ++it)
    need_update = it->needs_nnl_update();
  return need_update;
}



void packing_class::populate_nnl(const unsigned int i,
				 const double       influence_sphere,
				 const double       nnl_extra_distance)
//uses nnl_entry_periodic_class entries because the SLP algorithm is implemented currently for periodic boundaries
{
  std::vector<nnl_entry_periodic_class> nnl_influence;
  std::vector<nnl_entry_periodic_class> nnl_overlap;

  for (unsigned int j=i+1; j<get_n(); j++) {
    euclidean_vector<DIM,double> r_ij(sphere[j].get_local_position() - sphere[i].get_local_position());
    euclidean_vector<DIM,int> image;
    r_ij.get_minimum_image(image);
    const double pair_distance = euclidean_vector<DIM,double>::norm(lambda*r_ij);
    if (pair_distance <= influence_sphere + diameter + nnl_extra_distance) {
      nnl_entry_periodic_class new_influence_entry(&sphere[j],image);
      new_influence_entry.set_image(image);
      nnl_influence.push_back(new_influence_entry);
      if (pair_distance <= diameter + nnl_extra_distance) {
    	nnl_entry_periodic_class new_overlap_entry(&sphere[j],image);
    	new_overlap_entry.set_image(image);
    	nnl_overlap.push_back(new_overlap_entry);
      }
    }
  }
  sphere[i].set_nnls(nnl_influence,nnl_overlap);
  sphere[i].nnl_radius = nnl_extra_distance;
}
