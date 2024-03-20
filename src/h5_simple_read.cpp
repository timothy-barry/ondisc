#include <Rcpp.h>
#include "H5Cpp.h"
#include <string>
using namespace H5;
using namespace Rcpp;


//' @title Read the dimension of the ODM into memory
//' @param file_name_in path to the odm file
//' @param dataset_name name of the dataset within the file
//' @param length dimension of the integer vector to read
// [[Rcpp::export]]
IntegerVector read_integer_vector(const std::string& file_name_in, const std::string& dataset_name, int length) {
   // set name
   const H5std_string file_name(&file_name_in[0]);
   // open file
   H5File file(file_name, H5F_ACC_RDONLY);
   // open dataset
   const H5std_string dataset_name_h5(&dataset_name[0]);
   DataSet dataset = file.openDataSet(dataset_name_h5);
   // define out buffer
   IntegerVector buffer(length);
   // read data from disk into output buffer
   dataset.read(&buffer[0], PredType::NATIVE_INT);
   // close files and datasets
   dataset.close(); file.close();
   // return integer vector
   return buffer;
}


//' @title Read the feature IDs of an ODM into memory
//' @param file_name_in path to the odm file
//' @param n_features the number of features in the dataset
// [[Rcpp::export]]
StringVector read_feature_ids(const std::string& file_name_in, int n_features) {
  // set names
  const H5std_string file_name(&file_name_in[0]);
  // open file
  H5File file(file_name, H5F_ACC_RDONLY);
  // open dataset
  DataSet dataset = file.openDataSet("feature_ids");
  // initialize string array
  char** c_str_array = new char*[n_features];
  // Define string datatype
  StrType datatype(PredType::C_S1, H5T_VARIABLE);
  // Read the data
  dataset.read(c_str_array, datatype);
  // convert the data to a string array; free c_str_array memory
  StringVector out(n_features);
  for (int i = 0; i < n_features; i ++) {
    out[i] = std::string(c_str_array[i]);
    free(c_str_array[i]);
  }
  // close file, dataset, and datatype
  delete[] c_str_array;
  dataset.close(); file.close(); datatype.close();
  // return string vector
  return out;
}


//' @title Read the row pointer into memory
//' @param file_name_in path to the odm file
//' @param n_features the number of features in the dataset
// [[Rcpp::export]]
SEXP read_row_ptr(const std::string& file_name_in, int n_features) {
  // set name
  const H5std_string file_name(&file_name_in[0]);
  // open file
  H5File file(file_name, H5F_ACC_RDONLY);
  // open dataset
  DataSet dataset = file.openDataSet("ptr");
  // define out buffer
  std::vector<unsigned long long>* buffer = new std::vector<unsigned long long>(n_features + 1);
  // read data from disk into output buffer
  dataset.read(&((*buffer)[0]), PredType::NATIVE_ULLONG);
  // close files and datasets
  dataset.close(); file.close();
  // wrap buffer into an external pointer ptr
  Rcpp::XPtr<std::vector<unsigned long long>> ptr(buffer);
  return(ptr);
}
