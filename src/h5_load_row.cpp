#include <Rcpp.h>
#include "H5Cpp.h"
#include "shared_functs.h"
using namespace H5;
using namespace Rcpp;


// [[Rcpp::export]]
IntegerVector load_row_cpp(const std::string& file_name_in, SEXP f_row_ptr, int row_idx, int n_cells) {
  sparse_vector v = load_sparse_row_low_level(file_name_in, f_row_ptr, row_idx);
  std::vector<int> m_j = v.j, m_x = v.x;

  // 9. initialize the output vector
  IntegerVector out(n_cells, 0);
  for (int i = 0; i < m_x.size(); i ++) out[m_j[i]] = m_x[i];
  return out;
}


// [[Rcpp::export]]
IntegerVector threshold_count_matrix_cpp(const std::string& file_name_in, SEXP f_row_ptr, int row_idx, int threshold) {
  sparse_vector v = load_sparse_row_low_level(file_name_in, f_row_ptr, row_idx - 1);
  std::vector<int> j = v.j, x = v.x;
  IntegerVector out;
  int n_gte_threshold = 0, counter = 0;
  for (int k = 0; k < x.size(); k ++) if (x[k] >= threshold) n_gte_threshold ++;
  out = IntegerVector(n_gte_threshold);
  for (int k = 0; k < x.size(); k ++) if (x[k] >= threshold) out[counter ++] = j[k] + 1;
  return out;
}
