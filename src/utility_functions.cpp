#include <Rcpp.h>
using namespace Rcpp;


//' @title decrement a vector of indexes
//' @param idxs an integer vector of indexes
//' @noRd
// [[Rcpp::export]]
void decrement_idxs(IntegerVector idxs) {
  for (int i = 0; i < idxs.size(); i ++) {
    idxs[i] --;
  }
}


//' @title increment a vector of indexes
//' @param idxs an integer vector of indexes
//' @noRd
// [[Rcpp::export]]
void increment_idxs(IntegerVector idxs) {
  for (int i = 0; i < idxs.size(); i ++) {
    idxs[i] ++;
  }
}


//' @title add vectors in-place
//' @param v1 the vector to be modified in-place
//' @param v2 the vector to add element-wise to vector v1; note: v1 may be longer than v2
//' @noRd
// [[Rcpp::export]]
void sum_in_place(IntegerVector v1, IntegerVector v2) {
  for (int i = 0; i < v2.size(); i ++) {
    v1[i] += v2[i];
  }
}
