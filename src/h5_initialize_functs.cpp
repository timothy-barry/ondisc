#include <Rcpp.h>
#include "H5Cpp.h"
#include "shared_functs.h"
using namespace H5;
using namespace Rcpp;


char** convert_string_vector_to_c_array(StringVector x) {
  char** c_str_array = new char*[x.size()];
  char* temp_str;
  for (int i = 0; i < x.size(); i ++) {
    temp_str = new char[x[i].size() + 1];
    for (int k = 0; k < x[i].size(); k ++) temp_str[k] = x[i][k];
    temp_str[x[i].size()] = '\0';
    c_str_array[i] = temp_str;
  }
  return c_str_array;
}


void free_c_string_array(char** c_str_array, int size) {
  for (int i = 0; i < size; i++) delete[] c_str_array[i];
  delete[] c_str_array;
  return;
}


// [[Rcpp::export]]
void create_odm(const std::string& file_name_in, IntegerVector n_nonzero_features, StringVector feature_ids,
                int n_cells, int integer_id, int chunk_size, int compression_level) {
  // 0. define a few variables
  hsize_t n_features_ull = (hsize_t) n_nonzero_features.size(), chunk_size_ull = (hsize_t) chunk_size;
  int n_features_int = n_nonzero_features.size();

  // 1. compute the row pointer
  std::vector<unsigned long long> row_ptr = compute_row_ptr(n_nonzero_features);
  unsigned long long n_data_points = row_ptr[row_ptr.size() - 1];

  // 2. verify that n_data_points is less than h_size_t
  if (n_data_points > static_cast<unsigned long long>(std::numeric_limits<hsize_t>::max())) {
    Rcout << "The number of nonzero entries in the matrix exceeds the maximum integer size that HDF5 can handle; returning.";
    return;
  }

  // 3. create file
  const H5std_string file_name(&file_name_in[0]);
  H5File file(file_name, H5F_ACC_TRUNC);

  // 4. initialize the fields of the odm
  /****************
   * a. row pointer
   ****************/
  const H5std_string ptr_name("ptr");
  hsize_t ptr_dim[1] = {row_ptr.size()};
  DataSpace ptr_dataspace(1, ptr_dim);
  DataSet ptr_dataset = file.createDataSet(ptr_name, PredType::NATIVE_ULLONG, ptr_dataspace);
  // write the row ptr
  ptr_dataset.write(&row_ptr[0], PredType::NATIVE_ULLONG);
  // close dataset and dataspace
  ptr_dataset.close(); ptr_dataspace.close();

  /*******************
   * b. column indexes
   *******************/
  const H5std_string j_name("j");
  hsize_t j_dim[1] = {n_data_points};
  DataSpace j_dataspace(1, j_dim);
  DSetCreatPropList j_prop;
  j_prop.setChunk(1, &chunk_size_ull);
  j_prop.setDeflate(compression_level);
  DataSet j_dataset = file.createDataSet(j_name, PredType::NATIVE_INT, j_dataspace, j_prop);
  // close dataset and dataspace
  j_dataset.close(); j_dataspace.close(), j_prop.close();

  /**********
   * c. data
   *********/
  const H5std_string x_name("x");
  hsize_t x_dim[1] = {n_data_points};
  DataSpace x_dataspace(1, x_dim);
  DSetCreatPropList x_prop;
  x_prop.setChunk(1, &chunk_size_ull);
  x_prop.setDeflate(compression_level);
  DataSet x_dataset = file.createDataSet(x_name, PredType::NATIVE_INT, x_dataspace, x_prop);
  // close dataset and dataspace
  x_dataset.close(); x_dataspace.close(), x_prop.close();

  /**************
   * d. dimension
   **************/
  const H5std_string dimension_name("dimension");
  hsize_t dimension_dim[1] = {2};
  DataSpace dimension_dataspace(1, dimension_dim);
  DataSet dimension_dataset = file.createDataSet(dimension_name, PredType::NATIVE_INT, dimension_dataspace);
  // write the dimension
  std::vector<int> dimension {n_features_int, n_cells};
  dimension_dataset.write(&dimension[0], PredType::NATIVE_INT);
  // close dataset and dataspace
  dimension_dataset.close(); dimension_dataspace.close();

  /****************
   * e. feature IDs
   ****************/
  const H5std_string feature_ids_name("feature_ids");
  hsize_t feature_ids_dim[1] = {n_features_ull};
  DataSpace feature_ids_dataspace(1, feature_ids_dim);
  StrType feature_ids_datatype(PredType::C_S1, H5T_VARIABLE);
  DataSet feature_ids_dataset = file.createDataSet(feature_ids_name, feature_ids_datatype, feature_ids_dataspace);
  // convert the input string array to a char** array and write
  char** c_str_array = convert_string_vector_to_c_array(feature_ids);
  feature_ids_dataset.write(c_str_array, feature_ids_datatype);
  free_c_string_array(c_str_array, feature_ids.size());
  // close dataset, dataspace, and datatype
  feature_ids_dataset.close(); feature_ids_dataspace.close(); feature_ids_datatype.close();

  /***************
   * f. integer ID
   ***************/
  const H5std_string integer_id_name("integer_id");
  hsize_t integer_id_dim[1] = {1};
  DataSpace integer_id_dataspace(1, integer_id_dim);
  DataSet integer_id_dataset = file.createDataSet(integer_id_name, PredType::NATIVE_INT, integer_id_dataspace);
  // write the integer_id
  std::vector<int> integer_id_vect {integer_id};
  integer_id_dataset.write(&integer_id_vect[0], PredType::NATIVE_INT);
  // close dataset and dataspace
  integer_id_dataset.close(); integer_id_dataspace.close();

  // close file
  file.close();
  return;
}
