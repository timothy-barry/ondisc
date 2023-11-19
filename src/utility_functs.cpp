#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP init_ull_vect(int size) {
  std::vector<unsigned long long>* v = new std::vector<unsigned long long>(size, 0);
  Rcpp::XPtr<std::vector<unsigned long long>> ptr(v);
  return ptr;
}


// [[Rcpp::export]]
void print_ull_vect(SEXP ull_vect) {
  Rcpp::XPtr<std::vector<unsigned long long>> v(ull_vect);
  for (int i = 0; i < (*v).size(); i ++) {
    Rcout << (*v)[i] << " ";
  }
  return;
}

/*
// an ODM has the following fields.
// 1. a row pointer "p"
// length: number of features + 1.
// data type: H5T_NATIVE_ULLONG
// no chunking or compression
// 2. the column indices "j"
// length: number of nonzero data points
// data type:  H5T_NATIVE_INT
// chunking and compression (specified by chunk_size and compression_level)
// 3. the data "x"
// length: number of nonzero data points
// data type: H5T_NATIVE_INT
// chunking and compression (specified by chunk_size and compression_level)
// 4. the dimension "dimension"
// length: 2
// data type: H5T_NATIVE_INT
// no chunking or compression
// 5. the feature ids "feature_ids"
// length: number of features
// data type: H5T_VARIABLE
// no chunking or compression
//' Create an ODM
//'
//' @param file_name_in a string specifying the name of the odm file to create
//' @param n_nonzero_features an integer vector specifying the number of nonzero entries per feature
//' @param feature_ids a character vector specifying the feature IDs
//' @param n_cells an integer specifying the number of cells in the dataset
//' @param chunk_size the chunk size to set for the x and j vectors
//' @param compression_level the compression level to set for the x and j vectors
//' @return an ODM is created as a side-effect; NULL is returned
*/
