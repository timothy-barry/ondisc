#include <Rcpp.h>
#include <iostream>
#include <string>
#include "H5Cpp.h"

using namespace H5;
using namespace Rcpp;


//' @title write data h5
//' @param file_name_in name of h5 file
//' @param dataset_name_in name of the dataset within the h5 file
//' @param buffer the vector of integers
//' @param start_pos the position in the dataset in which to place the buffer
// [[Rcpp::export]]
void write_data_h5(const std::string& file_name_in, const std::string& dataset_name_in, IntegerVector buffer, int start_pos) {
  // set the file name
  const H5std_string file_name(&file_name_in[0]);
  // set the dataset within the file
  const H5std_string dataset_name(&dataset_name_in[0]);
  // open the file using the file_name
  H5File* file = new H5File(file_name, H5F_ACC_RDWR);
  // open the dataset
  DataSet* dataset = new DataSet(file->openDataSet(dataset_name));

  // define two dataspaces: on disk and in memory
  // 1. on disk (or file) data space
  DataSpace fspace = dataset->getSpace();
  //2 . in memory data space
  const hsize_t buffer_length = buffer.size();
  DataSpace mspace(1, &buffer_length);

  // set the hyperslab of both dataspaces
  // 1. on disk
  const hsize_t fcount = buffer_length, fstart = start_pos, fstride = 1, fblock = 1;
  fspace.selectHyperslab(H5S_SELECT_SET, &fcount, &fstart, &fstride, &fblock);
  // 2. in memory
  const hsize_t mcount = buffer_length, mstart = 0, mstride = 1, mblock = 1;
  mspace.selectHyperslab(H5S_SELECT_SET, &mcount, &mstart, &mstride, &mblock);

  // write the buffer data onto disk
  dataset->write(buffer.begin(), PredType::NATIVE_INT, mspace, fspace);

  // close the file and DataSet
  delete dataset;
  delete file;
  return;
}


//' @title map memory to disk
//' @param file_name_in name of h5 file
//' @param m_cell_idxs the cell indexes integer vector
//' @param cell_idxs_name name of the cell_idxs dataset
//' @param m_umi_counts the m_umi_counts vector
//' @param umi_counts_name the name of the integer vector
//' @param n_features number of features
//' @param m_row_ptr the memory row pointer
//' @param f_row_ptr the disk row pointer
// [[Rcpp::export]]
void map_memory_to_disk(const std::string& file_name_in, IntegerVector m_cell_idxs, const std::string& cell_idxs_name, IntegerVector m_umi_counts, const std::string& umi_counts_name, int n_features, IntegerVector m_row_ptr, IntegerVector f_row_ptr) {
  // set the file name
  const H5std_string file_name(&file_name_in[0]);
  // open the file using the file_name
  H5File* file = new H5File(file_name, H5F_ACC_RDWR);

  // open the cell_idxs dataset; first, set name
  const H5std_string cell_idxs_dataset_name(&cell_idxs_name[0]);
  // open the dataset
  DataSet* f_cell_idxs = new DataSet(file->openDataSet(cell_idxs_dataset_name));
  // set the dataspace for f_cell_idxs
  DataSpace f_cell_idxs_space = f_cell_idxs->getSpace();
  // set the dataspace for m_cell_idxs
  const hsize_t m_cell_idxs_size = m_cell_idxs.size();
  DataSpace m_cell_idxs_space(1, &m_cell_idxs_size);


  // open the umi_counts dataset; first, set name
  const H5std_string umi_counts_dataset_name(&umi_counts_name[0]);
  // open the dataset
  DataSet* f_umi_counts = new DataSet(file->openDataSet(umi_counts_dataset_name));
  // set dataspace for f_umi_counts
  DataSpace f_umi_counts_space = f_umi_counts->getSpace();
  // set the dataspace for m_umi_counts
  const hsize_t m_umi_counts_size = m_umi_counts.size();
  DataSpace m_umi_counts_space(1, &m_umi_counts_size);

  // define the variables to be used in the hyperslabs
  hsize_t fstart = 0, mcount = 0, mstart = 0, mstride = 1, mblock = 1, mend = 0;

  // loop through the data, mapping memory to disk
  for (int i = 0; i < n_features; i++) {
    mstart = m_row_ptr[i];
    mend = m_row_ptr[i + 1] - 1;
    mcount = mend - mstart + 1;
    if (mcount >= 1) {
      fstart = f_row_ptr[i];
      // set the hyperslabs for cell_idxs; write the data
      m_cell_idxs_space.selectHyperslab(H5S_SELECT_SET, &mcount, &mstart, &mstride, &mblock);
      f_cell_idxs_space.selectHyperslab(H5S_SELECT_SET, &mcount, &fstart, &mstride, &mblock);
      f_cell_idxs->write(m_cell_idxs.begin(), PredType::NATIVE_INT, m_cell_idxs_space, f_cell_idxs_space);

      // repeat for umi_counts
      m_umi_counts_space.selectHyperslab(H5S_SELECT_SET, &mcount, &mstart, &mstride, &mblock);
      f_umi_counts_space.selectHyperslab(H5S_SELECT_SET, &mcount, &fstart, &mstride, &mblock);
      f_umi_counts->write(m_umi_counts.begin(), PredType::NATIVE_INT, m_umi_counts_space, f_umi_counts_space);
      }
    }

  // close the datasets and file
  delete f_cell_idxs;
  delete f_umi_counts;
  delete file;
  return;
}


//' @title map memory to disk
//' @param file_name_in name of h5 file
//' @param m_cell_idxs the cell indexes integer vector
//' @param cell_idxs_name name of the cell_idxs dataset
//' @param n_features number of features
//' @param m_row_ptr the memory row pointer
//' @param f_row_ptr the disk row pointer
// [[Rcpp::export]]
void map_memory_to_disk_logical_matrix(const std::string& file_name_in, IntegerVector m_cell_idxs, const std::string& cell_idxs_name, int n_features, IntegerVector m_row_ptr, IntegerVector f_row_ptr) {
  // set the file name
  const H5std_string file_name(&file_name_in[0]);
  // open the file using the file_name
  H5File* file = new H5File(file_name, H5F_ACC_RDWR);

  // open the cell_idxs dataset; first, set name
  const H5std_string cell_idxs_dataset_name(&cell_idxs_name[0]);
  // open the dataset
  DataSet* f_cell_idxs = new DataSet(file->openDataSet(cell_idxs_dataset_name));
  // set the dataspace for f_cell_idxs
  DataSpace f_cell_idxs_space = f_cell_idxs->getSpace();
  // set the dataspace for m_cell_idxs
  const hsize_t m_cell_idxs_size = m_cell_idxs.size();
  DataSpace m_cell_idxs_space(1, &m_cell_idxs_size);

  // define the variables to be used in the hyperslabs
  hsize_t fstart = 0, mcount = 0, mstart = 0, mstride = 1, mblock = 1, mend = 0;

  // loop through the data, mapping memory to disk
  for (int i = 0; i < n_features; i++) {
    mstart = m_row_ptr[i];
    mend = m_row_ptr[i + 1] - 1;
    mcount = mend - mstart + 1;
    if (mcount >= 1) {
      fstart = f_row_ptr[i];
      // set the hyperslabs for cell_idxs; write the data
      m_cell_idxs_space.selectHyperslab(H5S_SELECT_SET, &mcount, &mstart, &mstride, &mblock);
      f_cell_idxs_space.selectHyperslab(H5S_SELECT_SET, &mcount, &fstart, &mstride, &mblock);
      f_cell_idxs->write(m_cell_idxs.begin(), PredType::NATIVE_INT, m_cell_idxs_space, f_cell_idxs_space);
      }
    }

  // close the datasets and file
  delete f_cell_idxs;
  delete file;
  return;
}
