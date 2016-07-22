
#include <iostream>

#include "macro.h"
#include "input_class.h"

input_class::input_class(int    argc,
			 char **argv) :
  valid(false)
{
  if (argc == NUM_INPUT_PARAMS) {
    str[input_file]  = argv[1];
    str[param_file]  = argv[2];
    str[output_file] = argv[3];
    valid=true;
  }  
}



void input_class::print_usage()
{
  std::cout<<"Usage: $ ./TJ {Initial Configuration} {Parameters} {Output}\n\n";
};
