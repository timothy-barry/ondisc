#include <Rcpp.h>
#include "H5Cpp.h"
#include "shared_functs.h"
#include <unordered_set>
using namespace H5;
using namespace Rcpp;


// [[Rcpp::export]]
void update_n_features_vector(IntegerVector feature_idx, IntegerVector n_nonzero_features_vector, int offset, int start_idx, int end_idx) {
  for (int i = start_idx; i <= end_idx; i++) {
    n_nonzero_features_vector[feature_idx[i] - offset] ++;
  }
  return;
}


// [[Rcpp::export]]
void decrement_vector(IntegerVector v) {
  for (unsigned long long i = 0; i < v.size(); i ++) v[i] --;
  return;
}


// [[Rcpp::export]]
void add_value_to_vector(IntegerVector v, int to_add) {
  for (unsigned long long i = 0; i < v.size(); i ++) v[i] += to_add;
  return;
}


// [[Rcpp::export]]
void compute_cellwise_covariates(IntegerVector n_umis, IntegerVector n_nonzero,
                                 NumericVector p_mito, IntegerVector feature_w_max_expression,
                                 NumericVector frac_umis_max_feature, IntegerVector feature_idx,
                                 IntegerVector j, IntegerVector x, IntegerVector mt_feature_idxs,
                                 int start_idx, int end_idx, int feature_offset, int cell_offset, int n_cells,
                                 bool compute_n_umis, bool compute_n_nonzero, bool compute_p_mito,
                                 bool compute_feature_w_max_expression, bool compute_frac_umis_max_feature) {
  int curr_feature_idx, curr_cell_idx, curr_x;
  bool mt_feature;
  std::vector<int> max_umi_count(compute_feature_w_max_expression ? n_cells : 0, 0);
  std::vector<int> n_mito_umis(compute_p_mito ? n_cells : 0, 0);
  // iterate through the entries of the mtx matrix
  for (int i = start_idx; i <= end_idx; i ++) {
    curr_feature_idx = feature_idx[i] - feature_offset;
    curr_cell_idx = j[i];
    curr_x = x[i];
    // n nonzero features
    if (compute_n_nonzero) n_nonzero[curr_cell_idx] ++;
    // n UMIs
    if (compute_n_umis) n_umis[curr_cell_idx] += curr_x;
    // max UMI count
    if (compute_feature_w_max_expression) {
      if (curr_x > max_umi_count[curr_cell_idx - cell_offset]) {
        max_umi_count[curr_cell_idx - cell_offset] = curr_x;
        feature_w_max_expression[curr_cell_idx] = curr_feature_idx;
      }
    }
    // n mito umis
    if (compute_p_mito) {
      mt_feature = std::find(mt_feature_idxs.begin(), mt_feature_idxs.end(), curr_feature_idx) != mt_feature_idxs.end();
      if (mt_feature) n_mito_umis[curr_cell_idx - cell_offset] += curr_x;
    }
  }

  if (compute_frac_umis_max_feature) {
    for (int cell_idx = cell_offset; cell_idx < cell_offset + n_cells; cell_idx ++) {
      frac_umis_max_feature[cell_idx] = (n_umis[cell_idx] == 0 ? 1 : (static_cast<double>(max_umi_count[cell_idx - cell_offset]))/(static_cast<double>(n_umis[cell_idx])));
    }
  }
  if (compute_p_mito) {
    for (int cell_idx = cell_offset; cell_idx < cell_offset + n_cells; cell_idx ++) {
      p_mito[cell_idx] = (n_umis[cell_idx] == 0 ? 0 : (static_cast<double>(n_mito_umis[cell_idx - cell_offset]))/(static_cast<double>(n_umis[cell_idx])));
    }
  }

  return;
}


// [[Rcpp::export]]
void compute_gene_expression_stats(
    NumericVector gene_sum,
    IntegerVector gene_count,
    IntegerVector feature_idx,
    IntegerVector x,
    int start_idx,
    int end_idx,
    int feature_offset
) {
  for (int i = start_idx; i <= end_idx; i++) {
    int curr_feature_idx = feature_idx[i] - feature_offset;
    int curr_x = x[i];
    
    gene_sum[curr_feature_idx] += curr_x;
    gene_count[curr_feature_idx] += 1;
  }
  return;
}


// [[Rcpp::export]]
void compute_normalized_gene_expression_stats(
    NumericVector gene_norm_sum,
    IntegerVector feature_idx,
    IntegerVector j,
    IntegerVector x,
    IntegerVector n_umis,
    int start_idx,
    int end_idx,
    int feature_offset,
    double scale_factor = 10000.0
) {
  // Only accumulate sum of normalized expression for non-zero entries
  // Zeros contribute 0 to the sum: log1p((0 / n_umis) * scale_factor) = 0
  // Mean will be computed in R by dividing by total number of cells
  for (int i = start_idx; i <= end_idx; i++) {
    int curr_feature_idx = feature_idx[i] - feature_offset;
    int curr_cell_idx = j[i];
    int curr_x = x[i];
    
    // Normalize expression: log1p((count / n_umis) * scale_factor)
    double normalized_expr = log1p((double(curr_x) / n_umis[curr_cell_idx]) * scale_factor);
    
    gene_norm_sum[curr_feature_idx] += normalized_expr;
  }
  return;
}


// [[Rcpp::export]]
void compute_cell_cycle_scores(
    NumericVector s_scores,
    NumericVector g2m_scores,
    IntegerVector s_gene_indices,
    IntegerVector g2m_gene_indices,
    IntegerVector s_control_indices,
    IntegerVector g2m_control_indices,
    IntegerVector feature_idx,
    IntegerVector j,
    IntegerVector x,
    IntegerVector n_umis,
    int start_idx,
    int end_idx,
    int feature_offset,
    int cell_offset,
    int n_cells,
    double scale_factor = 10000.0
) {
  // Initialize accumulators for each cell
  std::vector<double> s_target_sum(n_cells, 0.0);
  std::vector<double> g2m_target_sum(n_cells, 0.0);
  std::vector<double> s_control_sum(n_cells, 0.0);
  std::vector<double> g2m_control_sum(n_cells, 0.0);
  
  // Create lookup sets for fast membership testing
  std::unordered_set<int> s_genes(s_gene_indices.begin(), s_gene_indices.end());
  std::unordered_set<int> g2m_genes(g2m_gene_indices.begin(), g2m_gene_indices.end());
  std::unordered_set<int> s_controls(s_control_indices.begin(), s_control_indices.end());
  std::unordered_set<int> g2m_controls(g2m_control_indices.begin(), g2m_control_indices.end());
  
  // Process all non-zero entries
  for (int i = start_idx; i <= end_idx; i++) {
    int curr_feature_idx = feature_idx[i] - feature_offset;
    int curr_cell_idx = j[i] - cell_offset;
    int curr_x = x[i];
    
    // Skip if cell index is out of range
    if (curr_cell_idx < 0 || curr_cell_idx >= n_cells) continue;
    
    // Normalize expression: log1p((count / n_umis) * scale_factor)
    // This exactly matches Seurat's LogNormalize
    // j[i] is the global cell index after adjustment, use it directly for n_umis
    double normalized_expr = log1p((double(curr_x) / n_umis[j[i]]) * scale_factor);
    
    // Accumulate for S genes and controls
    if (s_genes.count(curr_feature_idx)) {
      s_target_sum[curr_cell_idx] += normalized_expr;
    }
    if (s_controls.count(curr_feature_idx)) {
      s_control_sum[curr_cell_idx] += normalized_expr;
    }
    
    // Accumulate for G2M genes and controls
    if (g2m_genes.count(curr_feature_idx)) {
      g2m_target_sum[curr_cell_idx] += normalized_expr;
    }
    if (g2m_controls.count(curr_feature_idx)) {
      g2m_control_sum[curr_cell_idx] += normalized_expr;
    }
  }
  
  // Compute final scores: mean(target) - mean(control)
  // CRITICAL: Compute means over ALL genes in module (including zeros), not just non-zero genes
  for (int cell_idx = 0; cell_idx < n_cells; cell_idx++) {
    // S genes: divide by total number of S genes, not just expressed ones
    double s_target_mean = s_target_sum[cell_idx] / s_gene_indices.size();
    
    // S controls: divide by total number of S control genes, not just expressed ones  
    double s_control_mean = s_control_sum[cell_idx] / s_control_indices.size();
    
    // G2M genes: divide by total number of G2M genes, not just expressed ones
    double g2m_target_mean = g2m_target_sum[cell_idx] / g2m_gene_indices.size();
    
    // G2M controls: divide by total number of G2M control genes, not just expressed ones
    double g2m_control_mean = g2m_control_sum[cell_idx] / g2m_control_indices.size();
    
    s_scores[cell_idx + cell_offset] = s_target_mean - s_control_mean;
    g2m_scores[cell_idx + cell_offset] = g2m_target_mean - g2m_control_mean;
  }
  
  return;
}


// [[Rcpp::export]]
void write_to_csr(const std::string& file_name_in, int start_idx, int end_idx, int feature_offset, int cell_offset,
                  int n_features, SEXP f_row_ptr, IntegerVector feature_idx, IntegerVector m_j, IntegerVector m_x) {
  // 1. dereference f_row_ptr
  Rcpp::XPtr<std::vector<unsigned long long>> f_ptr(f_row_ptr);

  // 2. compute the number of nonzero entries per feature and the memory row pointer (m_ptr)
  IntegerVector n_nonzero_features(n_features);
  update_n_features_vector(feature_idx, n_nonzero_features, feature_offset, start_idx, end_idx);
  std::vector<unsigned long long> m_ptr = compute_row_ptr(n_nonzero_features);

  // 3. open the h5 file with write access
  const H5std_string file_name(&file_name_in[0]);
  H5File file(file_name, H5F_ACC_RDWR);

  // 4. open the j dataset and get the dataspace
  const H5std_string f_j_name("j");
  DataSet f_j_dataset = file.openDataSet(f_j_name);
  DataSpace f_j_dataspace = f_j_dataset.getSpace();

  // 5. open the x dataset and get the dataspace
  const H5std_string f_x_name("x");
  DataSet f_x_dataset = file.openDataSet(f_x_name);
  DataSpace f_x_dataspace = f_x_dataset.getSpace();

  // 6. obtain the dataspace for m_j and m_x
  const hsize_t m_size = m_x.size();
  DataSpace m_j_dataspace(1, &m_size);
  DataSpace m_x_dataspace(1, &m_size);

  // 7. define the variables to be used in the hyperslabs
  hsize_t fstart = 0, mcount = 0, mstart = 0, mend = 0;

  // 8. loop through the data, mapping memory to disk
  for (int i = 0; i < n_features; i++) {
    mstart = m_ptr[i] + start_idx;
    mend = m_ptr[i + 1] - 1 + start_idx;
    mcount = mend - mstart + 1;
    if (mcount >= 1) {
      fstart = (*f_ptr)[i];

      // select the hyperslabs and write the data for j
      m_j_dataspace.selectHyperslab(H5S_SELECT_SET, &mcount, &mstart);
      f_j_dataspace.selectHyperslab(H5S_SELECT_SET, &mcount, &fstart);
      f_j_dataset.write(&m_j[0], PredType::NATIVE_INT, m_j_dataspace, f_j_dataspace);

      // repeat for x
      m_x_dataspace.selectHyperslab(H5S_SELECT_SET, &mcount, &mstart);
      f_x_dataspace.selectHyperslab(H5S_SELECT_SET, &mcount, &fstart);
      f_x_dataset.write(&m_x[0], PredType::NATIVE_INT, m_x_dataspace, f_x_dataspace);
    }
  }

  // 9. update f_ptr
  for (int i = 0; i < n_nonzero_features.size(); i ++) {
    (*f_ptr)[i] += static_cast<unsigned long long>(n_nonzero_features[i]);
  }

  // 10. close datasets, dataspaces, and file
  m_j_dataspace.close(); f_j_dataset.close(); f_j_dataspace.close();
  m_x_dataspace.close(); f_x_dataset.close(); f_x_dataspace.close();
  file.close();
  return;
}
