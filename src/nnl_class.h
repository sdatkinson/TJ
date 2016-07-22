#ifndef NNL_CLASS_H
#define NNL_CLASS_H

#include <vector>
#include <algorithm>

#include "const.h"
#include "euclidean_vector.h"

//Forward declaration of sphere_class:
class sphere_class;


//Basic NNL entry class
//Holds POINTERS to T
class nnl_entry_class {
 public:
  nnl_entry_class() {target=NULL;};
  nnl_entry_class(sphere_class* target_i) {set_target(target_i);};

  void          set_target(sphere_class* target_i) {target=target_i;};
  sphere_class* get_target() {return target;};
 private:
  sphere_class* target;
};



//NNL class for when neighbors might be in different periodic images
class nnl_entry_periodic_class : public nnl_entry_class {
 public:
  nnl_entry_periodic_class(sphere_class* target_i,
			   euclidean_vector<DIM,int>& image_i) {
    set_target(target_i);
    set_image(image_i);
  };

  void                      set_image(euclidean_vector<DIM,int>& image_i) {image=image_i;};
  euclidean_vector<DIM,int> get_image() {return image;};
 private:
  euclidean_vector<DIM,int> image;
};
  

#endif
