#include <Rcpp.h>
#include <iostream>
#include <string>
#include <vector>
#include "H5Cpp.h"

using namespace H5;
using namespace Rcpp;


//' @param file_name_in name of h5 file
//' @param p_name_in name of pointer
//' @param idx_name_in name of (minor) index
//' @param umi_counts_name_in name of umi counts
//' @param subset_vector the integer vector of indexes to extract
//' @param logical_mat boolean indicating whether the matrix is logical
// [[Rcpp::export]]
List index_h5_file(const std::string& file_name_in, const std::string& p_name_in, const std::string& idx_name_in, const std::string& umi_counts_name_in, IntegerVector subset_vector, bool logical_mat) {
  // open file
  const H5std_string file_name(&file_name_in[0]);
  H5File* file = new H5File(file_name, H5F_ACC_RDONLY);
  // set dataset names of (i) p, (ii) idxs, and (iii) umi_counts.
  const H5std_string p_name(&p_name_in[0]);
  const H5std_string idx_name(&idx_name_in[0]);

  // open the corresponding datasets
  DataSet* f_p = new DataSet(file->openDataSet(p_name));
  DataSet* f_idx = new DataSet(file->openDataSet(idx_name));

  // set the dataspace for these datasets
  DataSpace f_p_space = f_p->getSpace();
  DataSpace f_idx_space = f_idx->getSpace();

  // define vectors used in first round of computation, plus the data spaces
  int subset_vector_length = subset_vector.size();
  // ptr_storage
  std::vector<int> ptr_storage(2 * subset_vector_length);
  hsize_t ptr_storage_size = ptr_storage.size();
  DataSpace ptr_storage_space(1, &ptr_storage_size);
  // out_ptr
  IntegerVector out_ptr(subset_vector_length + 1);
  hsize_t out_ptr_size = out_ptr.size();
  DataSpace out_ptr_space(1, &out_ptr_size);
  // count_vect
  std::vector<hsize_t> count_vect(subset_vector_length);
  hsize_t count_vect_size = count_vect.size();
  DataSpace count_vect_space(1, &count_vect_size);

  // Next, define the variables to be used to set the hyperslabs.
  hsize_t m_start, f_start, count;
  const hsize_t stride = 1, block = 1, zero = 0;

  // Loop through subset vector, setting the values of ptr_storage
  count = 2;
  for (int i = 0; i < subset_vector_length; i++) {
    m_start = 2 * i;
    f_start = subset_vector[i];
    ptr_storage_space.selectHyperslab(H5S_SELECT_SET, &count, &m_start, &stride, &block);
    f_p_space.selectHyperslab(H5S_SELECT_SET, &count, &f_start, &stride, &block);
    f_p->read(&ptr_storage[0], PredType::NATIVE_INT, ptr_storage_space, f_p_space);
  }

  // populate count_vect and out_ptr
  for (int i = 0; i < subset_vector_length; i ++) {
    count_vect[i] = ptr_storage[(2 * i) + 1] - ptr_storage[2 * i];
  }

  out_ptr[0] = 0;
  for (int i = 1; i < subset_vector_length + 1; i ++) {
    out_ptr[i] += out_ptr[i - 1] + count_vect[i - 1];
  }

  // initialize idx
  hsize_t total_length = out_ptr[subset_vector_length];
  IntegerVector m_idx(total_length);
  DataSpace m_idx_space(1, &total_length);

  // initialize output
  List out;

  // take options on logical_mat
  if (!logical_mat) {
    // standard case -- integer matrix; define umi_count variables
    const H5std_string umi_counts_name(&umi_counts_name_in[0]);
    DataSet* f_umi_counts = new DataSet(file->openDataSet(umi_counts_name));
    DataSpace f_umi_counts_space = f_umi_counts->getSpace();
    IntegerVector m_umi_counts(total_length);
    DataSpace m_umi_counts_space(1, &total_length);

    // populate the vectors
    for (int i = 0; i < subset_vector_length; i ++) {
      m_start = (hsize_t) out_ptr[i];
      f_start = (hsize_t) ptr_storage[2 * i];
      count = count_vect[i];

      // first, the umis
      m_umi_counts_space.selectHyperslab(H5S_SELECT_SET, &count, &m_start, &stride, &block);
      f_umi_counts_space.selectHyperslab(H5S_SELECT_SET, &count, &f_start, &stride, &block);
      f_umi_counts->read(&m_umi_counts[0], PredType::NATIVE_INT, m_umi_counts_space, f_umi_counts_space);

      // next, the idxs
      m_idx_space.selectHyperslab(H5S_SELECT_SET, &count, &m_start, &stride, &block);
      f_idx_space.selectHyperslab(H5S_SELECT_SET, &count, &f_start, &stride, &block);
      f_idx->read(&m_idx[0], PredType::NATIVE_INT, m_idx_space, f_idx_space);
    }
    // set output
    out = List::create(out_ptr, m_idx, m_umi_counts);
    // close f_umi_counts file
    delete f_umi_counts;
  } else {
    // logical matrix case; no need to define umi vectors
    for (int i = 0; i < subset_vector_length; i ++) {
      m_start = (hsize_t) out_ptr[i];
      f_start = (hsize_t) ptr_storage[2 * i];
      count = count_vect[i];

      // the idxs only
      m_idx_space.selectHyperslab(H5S_SELECT_SET, &count, &m_start, &stride, &block);
      f_idx_space.selectHyperslab(H5S_SELECT_SET, &count, &f_start, &stride, &block);
      f_idx->read(&m_idx[0], PredType::NATIVE_INT, m_idx_space, f_idx_space);
    }
    // set output
    out = List::create(out_ptr, m_idx);
  }
  // close files and datasets
  delete f_p;
  delete f_idx;
  delete file;

  // return sparse matrix
  return out;
}
