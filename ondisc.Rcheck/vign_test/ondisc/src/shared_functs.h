#include <Rcpp.h>
using namespace Rcpp;

#ifndef COMPUTE_ROW_PTR
#define COMPUTE_ROW_PTR
std::vector<unsigned long long> compute_row_ptr(IntegerVector n_nonzero_features);
#endif

#ifndef SPARSE_VECTOR
#define SPARSE_VECTOR
struct sparse_vector {
  std::vector<int> x, j;
};
#endif

#ifndef LOAD_SPARSE_ROW_LOW_LEVEL
#define LOAD_SPARSE_ROW_LOW_LEVEL
sparse_vector load_sparse_row_low_level(const std::string& file_name_in, SEXP f_row_ptr, int row_idx, bool load_x);
#endif
