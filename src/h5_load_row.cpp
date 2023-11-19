#include <Rcpp.h>
#include "H5Cpp.h"
using namespace H5;
using namespace Rcpp;


// [[Rcpp::export]]
IntegerVector load_row_cpp(const std::string& file_name_in, SEXP f_row_ptr, int row_idx, int n_cells) {
  // 1. dereference f_row_ptr; determine the number of entries
  Rcpp::XPtr<std::vector<unsigned long long>> f_ptr(f_row_ptr);
  hsize_t start_pos = (*f_ptr)[row_idx];
  hsize_t n_entries = (*f_ptr)[row_idx + 1] - start_pos;

  // 2. open the file
  const H5std_string file_name(&file_name_in[0]);
  H5File file(file_name, H5F_ACC_RDONLY);

  // 3. open the x and j datasets
  DataSet f_x_dataset = file.openDataSet("x");
  DataSet f_j_dataset = file.openDataSet("j");

  // 4. initialize the dataspace for x and j
  DataSpace f_x_dataspace = f_x_dataset.getSpace();
  DataSpace f_j_dataspace = f_j_dataset.getSpace();

  // 4. allocate the output m_x and m_j vectors
  std::vector<int> m_x(n_entries), m_j(n_entries);

  // 5. initialize the dataspace for m_x and m_j
  DataSpace m_space(1, &n_entries);

  // 6. set the hyperslab for f_x and f_j
  f_x_dataspace.selectHyperslab(H5S_SELECT_SET, &n_entries, &start_pos);
  f_j_dataspace.selectHyperslab(H5S_SELECT_SET, &n_entries, &start_pos);

  // 7. read data
  f_x_dataset.read(&m_x[0], H5::PredType::NATIVE_INT, m_space, f_x_dataspace);
  f_j_dataset.read(&m_j[0], H5::PredType::NATIVE_INT, m_space, f_j_dataspace);

  // 8. close resources
  f_j_dataset.close(); f_j_dataspace.close();
  f_x_dataset.close(); f_x_dataspace.close();
  file.close(); m_space.close();

  // 9. initialize the output vector
  IntegerVector out(n_cells, 0);
  for (int i = 0; i < m_x.size(); i ++) out[m_j[i]] = m_x[i];
  return out;
}