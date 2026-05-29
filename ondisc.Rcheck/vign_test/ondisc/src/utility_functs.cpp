#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP init_ull_vect(int size) {
  std::vector<unsigned long long>* v = new std::vector<unsigned long long>(size, 0);
  Rcpp::XPtr<std::vector<unsigned long long>> ptr(v);
  return ptr;
}


// [[Rcpp::export]]
void print_ull_vect(SEXP ull_vect) {
  Rcpp::XPtr<std::vector<unsigned long long>> v(ull_vect);
  for (int i = 0; i < (*v).size(); i ++) {
    Rcout << (*v)[i] << " ";
  }
  return;
}


// [[Rcpp::export]]
void update_dt_column(IntegerVector col, IntegerVector overwrite_vector, int start) {
  for (int i = start; i < start + overwrite_vector.size(); i ++) {
    col[i] = overwrite_vector[i - start];
  }
  return;
}
