#include <Rcpp.h>
#include <iostream>
#include <string>
#include <vector>
#include "H5Cpp.h"

using namespace H5;
using namespace Rcpp;


//' @title read intege vector from hdf5
//' @param file_name_in name of h5 file
//' @param dataset_name_in name of dataset to read
//' @param data_len length of integer vector
// [[Rcpp::export]]
IntegerVector read_integer_vector_hdf5(const std::string& file_name_in, const std::string& dataset_name_in, int data_len) {
  // set names
  const H5std_string file_name(&file_name_in[0]);
  const H5std_string dataset_name(&dataset_name_in[0]);
  // open file
  H5File* file = new H5File(file_name, H5F_ACC_RDONLY);
  // open dataset
  DataSet* dataset = new DataSet(file->openDataSet(dataset_name));
  // define out buffer
  IntegerVector buffer(data_len);
  // read data from disk into output buffer
  dataset->read(&buffer[0], PredType::NATIVE_INT);
  // close files and datasets
  delete dataset;
  delete file;
  // return integer vector
  return buffer;
}
