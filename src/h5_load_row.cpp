#include <Rcpp.h>
#include "H5Cpp.h"
#include "shared_functs.h"
using namespace H5;
using namespace Rcpp;


// [[Rcpp::export]]
IntegerVector load_row_cpp(const std::string& file_name_in, SEXP f_row_ptr, int row_idx, int n_cells) {
  sparse_vector v = load_sparse_row_low_level(file_name_in, f_row_ptr, row_idx, true);
  std::vector<int> m_j = v.j, m_x = v.x;

  // 9. initialize the output vector
  IntegerVector out(n_cells, 0);
  for (int i = 0; i < m_x.size(); i ++) out[m_j[i]] = m_x[i];
  return out;
}

