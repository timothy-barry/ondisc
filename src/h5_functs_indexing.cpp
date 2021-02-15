#include <Rcpp.h>
#include <iostream>
#include <string>
#include <vector>
#include "H5Cpp.h"

using namespace H5;
using namespace Rcpp;


// [[Rcpp::export]]
List index_h5_file(const std::string& file_name_in, const std::string& p_name_in, const std::string& idx_name_in, const std::string& umi_counts_name_in, IntegerVector contiguous_idx_range) {
  const H5std_string file_name(&file_name_in[0]);
  // open file
  H5File* file = new H5File(file_name, H5F_ACC_RDONLY);
  // set dataset names of (i) p, (ii) idxs, and (iii) umi_counts.
  const H5std_string p_name(&p_name_in[0]);
  const H5std_string idx_name(&idx_name_in[0]);
  const H5std_string umi_counts_name(&umi_counts_name_in[0]);

  // open the corresponding datasets
  DataSet* f_p = new DataSet(file->openDataSet(p_name));
  DataSet* f_idx = new DataSet(file->openDataSet(idx_name));
  DataSet* f_umi_counts = new DataSet(file->openDataSet(umi_counts_name));
  // set the dataspace for these datasets
  DataSpace f_p_space = f_p->getSpace();
  DataSpace f_idx_space = f_idx->getSpace();
  DataSpace f_umi_counts_space = f_umi_counts->getSpace();

  // define the variables to be used in the hyperslabs
  hsize_t major_idx_start, major_idx_end, major_idx_count;
  const hsize_t stride = 1, block = 1, zero = 0;
  major_idx_start = contiguous_idx_range[0];
  major_idx_end = contiguous_idx_range[1];
  major_idx_count = major_idx_end - major_idx_start + 2;

  // initialize m_p buffer
  std::vector<int> m_p(major_idx_count);
  DataSpace m_p_space(1, &major_idx_count);

  // set hyperslab for f_p and m_p; read data into buffer
  m_p_space.selectHyperslab(H5S_SELECT_SET, &major_idx_count, &zero, &stride, &block);
  f_p_space.selectHyperslab(H5S_SELECT_SET, &major_idx_count, &major_idx_start, &stride, &block);
  f_p->read(&m_p[0], PredType::NATIVE_INT, m_p_space, f_p_space);

  // convert m_p from integer vector to hsize_t vector (NOTE: think more carefully about data types later)
  std::vector<hsize_t> m_p_hsize(m_p.begin(), m_p.end());

  // center the pointer vector m_p
  hsize_t initial_p = m_p[0];
  for (int i = 0; i < major_idx_count; i ++) {
    m_p[i] -= initial_p;
  }
  hsize_t p_count = m_p[major_idx_count - 1];

  // define the idx and umi_counts buffers, along with the data spaces
  IntegerVector m_idx(p_count);
  DataSpace m_idx_space(1, &p_count);
  IntegerVector m_umi_counts(p_count);
  DataSpace m_umi_counts_space(1, &p_count);

  if (p_count > 0) { // if there are entries to extract, ...
    // set hyperslabs
    m_idx_space.selectHyperslab(H5S_SELECT_SET, &p_count, &zero, &stride, &block);
    m_umi_counts_space.selectHyperslab(H5S_SELECT_SET, &p_count, &zero, &stride, &block);
    f_idx_space.selectHyperslab(H5S_SELECT_SET, &p_count, &initial_p, &stride, &block);
    f_umi_counts_space.selectHyperslab(H5S_SELECT_SET, &p_count, &initial_p, &stride, &block);
    // read data
    f_idx->read(&m_idx[0], PredType::NATIVE_INT, m_idx_space, f_idx_space);
    f_umi_counts->read(&m_umi_counts[0], PredType::NATIVE_INT, m_umi_counts_space, f_umi_counts_space);
  }
  return List::create(m_p, m_idx, m_umi_counts);
}
