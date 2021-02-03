#include <Rcpp.h>
using namespace Rcpp;


//' @title decrement a vector of indexes
//' @param idxs an integer vector of indexes
// [[Rcpp::export]]
void decrement_idxs(IntegerVector idxs) {
  for (int i = 0; i < idxs.size(); i ++) {
    idxs[i] --;
  }
}


//' @title add vectors in-place
//' @param v1 the vector to be modified in-place
//' @param v2 the vector to add element-wise to vector v1; note: v1 may be longer than v2
// [[Rcpp::export]]
void sum_in_place(IntegerVector v1, IntegerVector v2) {
  for (int i = 0; i < v2.size(); i ++) {
    v1[i] += v2[i];
  }
}


//' @title find contiguous subsequences
//' @param v an ORDERED integer vector
// [[Rcpp::export]]
List find_contig_subseqs(IntegerVector v) {
  IntegerVector starting_idxs;
  IntegerVector ending_idxs;
  starting_idxs.push_back(v[0]);
  for (int i = 1; i < v.size(); i ++) {
    if (v[i] != v[i - 1] + 1) {
      ending_idxs.push_back(v[i - 1]);
      starting_idxs.push_back(v[i]);
    }
  }
  ending_idxs.push_back(v[v.size() - 1]);
  return List::create(starting_idxs, ending_idxs);
}
