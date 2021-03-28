#include <Rcpp.h>
#include <fstream>
#include <sstream>

using Rcpp::CharacterVector;
using Rcpp::List;

//' Get mtx metadata
//'
//' @param mtx_fp path to the mtx file
//'
//' @return a list containing (i) n_genes, (ii) n_cells, (iii) the number of
//'     data points (i.e., fraction of entries that are zero),
//'     (iv) (TRUE/FALSE) matrix is logical,
//'     (v) number of rows to skip before reading the data
//' @noRd
// [[Rcpp::export]]
List get_mtx_metadata(CharacterVector mtx_fp)
{
  // Convert to C++ string
  std::string path = Rcpp::as<std::string>(mtx_fp);
  // Open file
  std::ifstream file(path);
  if(!file.is_open())
    Rcpp::stop("cannot open file '" + path + "'");

  // Find the first line that does not begin with "%"
  std::string line;
  int n_features = 0, n_cells = 0, n_data_points = 0, n_rows_to_skip = 0;
  // The return value of std::getline evaluates to be true
  // if a line is read successfully
  while(std::getline(file, line))
  {
    if(line[0] != '%')
    {
      // The first line that does not begin with "%" has the format
      //   rows columns num_of_nonzeros
      std::stringstream ss(line);
      ss >> n_features >> n_cells >> n_data_points;
      // Also skip this line
      n_rows_to_skip++;
      break;
    } else{
      // Skip this line if it is a comment
      n_rows_to_skip++;
    }
  }

  // Largest number of data points - 2147483647
  const int integer_max = ~(1 << 31);
  if(n_data_points > integer_max)
    Rcpp::stop("Numer of data points exceeds maximum value.");

  bool is_logical = false;
  // Read the first data line
  if(std::getline(file, line))
  {
    std::stringstream ss(line);
    int i, j, x;
    // Test whether there are three values
    bool status = bool(ss >> i >> j >> x);
    // If status is false, it means there are only two values,
    // indicating a logical matrix
    is_logical = !status;
  }

  return List::create(
    Rcpp::Named("n_features")     = n_features,
    Rcpp::Named("n_cells")        = n_cells,
    Rcpp::Named("n_data_points")  = n_data_points,
    Rcpp::Named("is_logical")     = is_logical,
    Rcpp::Named("n_rows_to_skip") = n_rows_to_skip
  );
}
