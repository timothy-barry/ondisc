#include <Rcpp.h>
#include <iostream>
#include <string>
#include "H5Cpp.h"

using namespace H5;
using namespace Rcpp;


//' @title write data h5
//' @param file_name_in name of h5 file
//' @param dataset_name_in name of the dataset within the h5 file
//' @param buffer the vector of integers
//' @param start_pos the position in the dataset in which to place the buffer
// [[Rcpp::export]]
void write_data_h5(const std::string& file_name_in, const std::string& dataset_name_in, IntegerVector buffer, int start_pos) {
  // set the file name
  const H5std_string file_name(&file_name_in[0]);
  // set the dataset within the file
  const H5std_string dataset_name(&dataset_name_in[0]);
  // open the file using the file_name
  H5File* file = new H5File(file_name, H5F_ACC_RDWR);
  // open the dataset
  DataSet* dataset = new DataSet(file->openDataSet(dataset_name));

  // define two dataspaces: on disk and in memory
  // 1. on disk (or file) data space
  DataSpace fspace = dataset->getSpace();
  //2 . in memory data space
  const hsize_t buffer_length = buffer.size();
  DataSpace mspace(1, &buffer_length);

  // set the hyperslab of both dataspaces
  // 1. on disk
  const hsize_t fstart = start_pos;
  const hsize_t fstride = 1;
  const hsize_t fcount = buffer_length;
  const hsize_t fblock = 1;
  fspace.selectHyperslab(H5S_SELECT_SET, &fcount, &fstart, &fstride, &fblock);
  // 2. in memory
  const hsize_t mstart = 0;
  const hsize_t mstride = 1;
  const hsize_t mcount = buffer_length;
  const hsize_t mblock = 1;
  mspace.selectHyperslab(H5S_SELECT_SET, &mcount, &mstart, &mstride, &mblock);

  // write the buffer data onto disk
  dataset->write(buffer.begin(), PredType::NATIVE_INT, mspace, fspace);

  // close the file and DataSet
  delete dataset;
  delete file;
}
