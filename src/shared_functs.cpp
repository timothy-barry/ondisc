#include <Rcpp.h>
using namespace Rcpp;

std::vector<unsigned long long> compute_row_ptr(IntegerVector n_nonzero_features) {
  unsigned long long curr_sum = 0;
  std::vector<unsigned long long> row_ptr(n_nonzero_features.size() + 1);
  row_ptr[0] = 0;
  for (int i = 0; i < n_nonzero_features.size(); i ++) {
    curr_sum += (unsigned long long) n_nonzero_features[i];
    row_ptr[i + 1] = curr_sum;
  }
  return row_ptr;
}
