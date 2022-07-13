#include <Rcpp.h>
#include <iostream>
#include <string>
#include "H5Cpp.h"

using namespace H5;
using namespace Rcpp;

void create_h5_dataset (const std::string& file_name_in, const std::string& dataset_name_in, const std::string& int_type, int length, int comp_level, int chunk_size) {
  
}
