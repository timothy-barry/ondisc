#include <Rcpp.h>
using namespace Rcpp;

#ifndef COMPUTE_ROW_PTR
#define COMPUTE_ROW_PTR
std::vector<unsigned long long> compute_row_ptr(IntegerVector n_nonzero_features);
#endif
