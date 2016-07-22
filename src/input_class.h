

#ifndef INPUT_CLASS_H
#define INPUT_CLASS_H

#include <string>

class input_class {
 public:

  //=================================================================
  //  Public members:

  enum e_str{
    input_file,
    param_file,
    output_file,
    num_str
  };

  //=================================================================
  //  Public member functions:

  input_class(int    argc,
	      char **argv);

  bool        get_valid()    {return valid;};
  std::string get_str(int i) {return str[i];};

  void        print_usage();

 private:

  //=================================================================
  //  Private members:

  bool valid;
  std::string str[num_str];
};

#endif
