#include <Rcpp.h>
#include "H5Cpp.h"
#include "shared_functs.h"
using namespace H5;
using namespace Rcpp;


// [[Rcpp::export]]
IntegerVector threshold_count_matrix_cpp(const std::string& file_name_in, SEXP f_row_ptr, int row_idx, int threshold) {
  sparse_vector v = load_sparse_row_low_level(file_name_in, f_row_ptr, row_idx - 1, true);
  std::vector<int> j = v.j, x = v.x;
  IntegerVector out;
  int n_gte_threshold = 0, counter = 0;
  for (int k = 0; k < x.size(); k ++) if (x[k] >= threshold) n_gte_threshold ++;
  out = IntegerVector(n_gte_threshold);
  for (int k = 0; k < x.size(); k ++) if (x[k] >= threshold) out[counter ++] = j[k] + 1;
  return out;
}



void load_nonzero_posits_odm(const std::string& file_name_in, SEXP f_row_ptr, int row_idx, std::vector<bool>& y_orig,
                             std::vector<bool>& y_sub, std::vector<int>& cells_in_use_zero_idx) {
  sparse_vector v = load_sparse_row_low_level(file_name_in, f_row_ptr, row_idx, false);
  std::vector<int> j = v.j;
  for (int k = 0; k < y_orig.size(); k ++) y_orig[k] = false;
  for (int k = 0; k < j.size(); k ++) y_orig[j[k]] = true;
  for (int k = 0; k < y_sub.size(); k++) y_sub[k] = y_orig[cells_in_use_zero_idx[k]];
  return;
}


// [[Rcpp::export]]
List compute_nt_nonzero_matrix_and_n_ok_pairs_ondisc(const std::string& file_name_in, SEXP f_row_ptr, int n_genes, int n_cells_orig,
                                                     int n_cells_sub, List grna_group_idxs, List indiv_nt_grna_idxs, IntegerVector all_nt_idxs,
                                                     IntegerVector to_analyze_response_idxs, IntegerVector to_analyze_grna_idxs,
                                                     bool control_group_complement, IntegerVector cells_in_use) {
  // 0. initialize variables and objects
  int n_nt_grnas = indiv_nt_grna_idxs.size(), n_pairs = to_analyze_response_idxs.size(), pair_pointer = 0, n_nonzero = 0, curr_n_nonzero_cntrl = 0, curr_n_nonzero_trt = 0;
  IntegerVector curr_idxs, n_nonzero_tot(n_genes), n_nonzero_trt(n_pairs), n_nonzero_cntrl(n_pairs);
  IntegerMatrix M(n_nt_grnas, n_genes);

  std::vector<bool> y_sub(n_cells_sub), y_orig(n_cells_orig);
  std::vector<int> cells_in_use_zero_idx(n_cells_sub);
  for (int i = 0; i < cells_in_use_zero_idx.size(); i ++) cells_in_use_zero_idx[i] = cells_in_use[i] - 1;

  // 1. iterate over genes
  for (int row_idx = 0; row_idx < n_genes; row_idx++) {
    // load nonzero positions into the boolean vector y_sub
    load_nonzero_posits_odm(file_name_in, f_row_ptr, row_idx, y_orig, y_sub, cells_in_use_zero_idx);
    // 1.2 iterate over nt grnas, adding n nonzero trt to M matrix
    for (int grna_idx = 0; grna_idx < n_nt_grnas; grna_idx ++) {
      n_nonzero = 0;
      curr_idxs = indiv_nt_grna_idxs[grna_idx];
      for (int k = 0; k < curr_idxs.size(); k ++) {
        if (!control_group_complement) { // NT cells control group
          if (y_sub[all_nt_idxs[curr_idxs[k] - 1] - 1]) n_nonzero ++; // redirection
        } else { // discovery cells control group
          if (y_sub[curr_idxs[k] - 1]) n_nonzero ++; // no redirection
        }
      }
      M(grna_idx, row_idx) = n_nonzero;
    }

    // 1.2 update n_nonzero_tot for this gene
    n_nonzero_tot[row_idx] = 0;
    if (!control_group_complement) {
      for (int k = 0; k < n_nt_grnas; k ++) n_nonzero_tot[row_idx] += M(k, row_idx);
    } else {
      for (int k = 0; k < y_sub.size(); k++) if (y_sub[k]) n_nonzero_tot[row_idx] ++;
    }

    // iterate through the discovery pairs containing this gene
    while (to_analyze_response_idxs[pair_pointer] - 1 == row_idx) {
      curr_idxs = grna_group_idxs[to_analyze_grna_idxs[pair_pointer] - 1];
      curr_n_nonzero_trt = 0;
      // first, get n nonzero trt
      for (int k = 0; k < curr_idxs.size(); k ++) if (y_sub[curr_idxs[k] - 1]) curr_n_nonzero_trt ++;
      // second, get n nonzero cntrl
      if (!control_group_complement) { // NT cells
        curr_n_nonzero_cntrl = n_nonzero_tot[row_idx];
      } else { // complement set
        curr_n_nonzero_cntrl = n_nonzero_tot[row_idx] - curr_n_nonzero_trt;
      }
      n_nonzero_trt[pair_pointer] = curr_n_nonzero_trt;
      n_nonzero_cntrl[pair_pointer] = curr_n_nonzero_cntrl;
      pair_pointer ++;
    }
  }

  return List::create(Named("n_nonzero_mat") = M, Named("n_nonzero_tot") = n_nonzero_tot,
                      Named("n_nonzero_trt") = n_nonzero_trt, Named("n_nonzero_cntrl") = n_nonzero_cntrl);
}


// [[Rcpp::export]]
IntegerMatrix compute_n_trt_cells_matrix_ondisc(const std::string& file_name_in, SEXP f_row_ptr, int n_cells_orig,
                                                int n_cells_sub, int n_genes, List nt_grna_group_idxs, IntegerVector cells_in_use) {
  // define objects
  int n_nonzero;
  IntegerVector curr_idxs;
  IntegerMatrix M(nt_grna_group_idxs.size(), n_genes);
  std::vector<bool> y_sub(n_cells_sub), y_orig(n_cells_orig);
  std::vector<int> cells_in_use_zero_idx(n_cells_sub);
  for (int i = 0; i < cells_in_use_zero_idx.size(); i ++) cells_in_use_zero_idx[i] = cells_in_use[i] - 1;

  // decrement nt_grna_group_idxs
  for (int i = 0; i < nt_grna_group_idxs.size(); i ++) {
    curr_idxs = nt_grna_group_idxs[i];
    for (int k = 0; k < curr_idxs.size(); k ++) curr_idxs[k] --;
  }

  // iterate over genes
  for (int row_idx = 0; row_idx < n_genes; row_idx ++) {
    // load nonzero positions into the boolean vector y_sub
    load_nonzero_posits_odm(file_name_in, f_row_ptr, row_idx, y_orig, y_sub, cells_in_use_zero_idx);
    // iterate over nt grna groups
    for (int grna_idx = 0; grna_idx < nt_grna_group_idxs.size(); grna_idx ++) {
      n_nonzero = 0;
      curr_idxs = nt_grna_group_idxs[grna_idx];
      for (int k = 0; k < curr_idxs.size(); k ++) if (y_sub[curr_idxs[k]]) n_nonzero ++;
      M(grna_idx, row_idx) = n_nonzero;
    }
  }
  return M;
}


// [[Rcpp::export]]
List compute_n_ok_pairs_ondisc(const std::string& file_name_in, SEXP f_row_ptr, int n_genes, int n_cells_orig,
                               int n_cells_sub, List grna_group_idxs, IntegerVector all_nt_idxs,
                               IntegerVector to_analyze_response_idxs, IntegerVector to_analyze_grna_idxs,
                               bool control_group_complement, IntegerVector cells_in_use, IntegerVector unique_response_idxs) {
  // 0. initialize variables and objects
  int n_pairs = to_analyze_response_idxs.size(), pair_pointer = 0, curr_n_nonzero_cntrl = 0, curr_n_nonzero_trt = 0, row_idx = 0, n_nonzero_tot = 0;
  IntegerVector curr_idxs, n_nonzero_trt(n_pairs), n_nonzero_cntrl(n_pairs);
  std::vector<bool> y_sub(n_cells_sub), y_orig(n_cells_orig);
  std::vector<int> cells_in_use_zero_idx(n_cells_sub);
  for (int i = 0; i < cells_in_use_zero_idx.size(); i ++) cells_in_use_zero_idx[i] = cells_in_use[i] - 1;

  // 1. iterate over unique genes
  for (int i = 0; i < unique_response_idxs.size(); i++) {
    row_idx = unique_response_idxs[i] - 1;
    // load nonzero positions into the boolean vector y_sub
    load_nonzero_posits_odm(file_name_in, f_row_ptr, row_idx, y_orig, y_sub, cells_in_use_zero_idx);

    // update n_nonzero_tot for this gene
    n_nonzero_tot = 0;
    if (!control_group_complement) {
      for (int k = 0; k < all_nt_idxs.size(); k ++) if (y_sub[all_nt_idxs[k] - 1]) n_nonzero_tot ++;
    } else {
      for (int k = 0; k < y_sub.size(); k++) if (y_sub[k]) n_nonzero_tot ++;
    }

    // iterate through the discovery pairs containing this gene
    while (to_analyze_response_idxs[pair_pointer] - 1 == row_idx) {
      curr_idxs = grna_group_idxs[to_analyze_grna_idxs[pair_pointer] - 1];
      curr_n_nonzero_trt = 0;
      // first, get n nonzero trt
      for (int k = 0; k < curr_idxs.size(); k ++) if (y_sub[curr_idxs[k] - 1]) curr_n_nonzero_trt ++;
      // second, get n nonzero cntrl
      if (!control_group_complement) { // NT cells
        curr_n_nonzero_cntrl = n_nonzero_tot;
      } else { // complement set
        curr_n_nonzero_cntrl = n_nonzero_tot - curr_n_nonzero_trt;
      }
      n_nonzero_trt[pair_pointer] = curr_n_nonzero_trt;
      n_nonzero_cntrl[pair_pointer] = curr_n_nonzero_cntrl;
      pair_pointer ++;
    }
  }

  return List::create(Named("n_nonzero_trt") = n_nonzero_trt, Named("n_nonzero_cntrl") = n_nonzero_cntrl);
}
