#include <Rcpp.h>
#include "H5Cpp.h"
using namespace H5;
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

// sparse vector struct
struct sparse_vector {
  std::vector<int> x, j;
};

// load sparse row low level
sparse_vector load_sparse_row_low_level(const std::string& file_name_in, SEXP f_row_ptr, int row_idx, bool load_x) {
  // 1. dereference f_row_ptr; determine the number of entries
  Rcpp::XPtr<std::vector<unsigned long long>> f_ptr(f_row_ptr);
  hsize_t start_pos = (*f_ptr)[row_idx];
  hsize_t n_entries = (*f_ptr)[row_idx + 1] - start_pos;

  // 2. open the file
  const H5std_string file_name(&file_name_in[0]);
  H5File file(file_name, H5F_ACC_RDONLY);

  // open dataset
  DataSet f_j_dataset = file.openDataSet("j");
  // initialize dataspace
  DataSpace f_j_dataspace = f_j_dataset.getSpace();
  // allocate output
  std::vector<int> m_j(n_entries);
  // initialize the dataspace for m_x and m_j
  DataSpace m_space(1, &n_entries);
  // set hyperslab
  f_j_dataspace.selectHyperslab(H5S_SELECT_SET, &n_entries, &start_pos);
  // read data
  f_j_dataset.read(&m_j[0], H5::PredType::NATIVE_INT, m_space, f_j_dataspace);
  // close resources
  f_j_dataset.close(); f_j_dataspace.close();

  std::vector<int> m_x;
  if (load_x) {
    // open dataset
    DataSet f_x_dataset = file.openDataSet("x");
    // initialize the dataspace
    DataSpace f_x_dataspace = f_x_dataset.getSpace();
    // allocate output
    m_x = std::vector<int>(n_entries);
    // set hyperslab
    f_x_dataspace.selectHyperslab(H5S_SELECT_SET, &n_entries, &start_pos);
    // read data
    f_x_dataset.read(&m_x[0], H5::PredType::NATIVE_INT, m_space, f_x_dataspace);
    // close resources
    f_x_dataset.close(); f_x_dataspace.close();
  }
  // close shared resources
  file.close(); m_space.close();

  // 9. return sparse vector struct
  sparse_vector out = {m_x, m_j};
  return out;
}
