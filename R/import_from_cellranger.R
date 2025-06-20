#' Create `odm` object from Cell Ranger
#'
#' `create_odm_from_cellranger()` initializes an `odm` object, taking as input the output of one or more calls to Cell Ranger count. The number of `odm` objects returned corresponds to the number of modalities in the input data. Additionally, the cell-wise covariate data frame is computed and returned. `create_odm_from_cellranger()` supports the Cell Ranger modalities "Gene Expression", "CRISPR Guide Capture", and "Antibody Capture".
#'
#' @note The `grna_target_data_frame` is relevant only for CRISPR screen data (i.e., data for which the "CRISPR Guide Capture" modality is present). In single-cell CRISPR screens, gRNAs are delivered to cells via a viral vector. Some recent single-cell CRISPR screens involve a special design in which each viral vector harbors multiple gRNAs. For example, [Replogle 2022](https://pubmed.ncbi.nlm.nih.gov/35688146/) conducted a screen in which each viral vector contained two gRNAs, each targeting the same site. In such screens, users may wish to "collapse" the gRNA count matrix by summing over the UMI counts of gRNAs contained on the same vector. To do so, users can pass the argument `grna_target_data_frame`, which is a data frame containing two columns: `grna_id` and `vector_id`. `grna_id` should coincide with the gRNA IDs as contained within the `features.tsv` file, and `vector_id` should be a string indicating the vector to which a given gRNA ID belongs. The expression vectors of gRNAs contained within the same vector are summed.
#' @note The arguments `chunk_size` and `compression_level` control the extent to which the backing `.odm` files are compressed, with higher values corresponding to smaller file sizes (albeit possibly longer read and write times).
#'
#' @param directories_to_load a character vector specifying the locations of the directories to load. Each directory should contain the files "matrix.mtx.gz" and "features.tsv.gz" (and optionally "barcodes.tsv.gz", which is ignored).
#' @param directory_to_write a string indicating the directory in which to write the backing .odm files.
#' @param write_cellwise_covariates (optional; default `TRUE`) a logical value indicating whether to write the cellwise covariate data frame to disk (`TRUE`) in addition to returning it from `create_odm_from_cellranger()`.
#' @param chunk_size (optional; default `1000L`) a positive integer specifying the chunk size to use to store the data in the backing HDF5 file.
#' @param compression_level (optional; default `3L`) an integer in the inveral \[0, 9\] specifying the compression level to use to store the data in the backing HDF5 file.
#' @param grna_target_data_frame (optional) a data frame mapping each gRNA ID to its target. Relevant only if the CRISPR modality is present within the data. (See note.)
#' @param compute_cell_cycle (optional; default `FALSE`) a logical value indicating whether to compute cell cycle scores during import. Uses Seurat's default parameters and gene sets.
#'
#' @return a list containing the odm object(s) and cellwise covariates. If `compute_cell_cycle = TRUE`, the cellwise covariates will include columns `gene_s_score`, `gene_g2m_score`, and `gene_phase`.
#' @export
#'
#' @examples
#' library(sceptredata)
#' directories_to_load <- paste0(
#'  system.file("extdata", package = "sceptredata"),
#'  "/highmoi_example/gem_group_", c(1, 2)
#' )
#' directory_to_write <- tempdir()
#' out_list <- create_odm_from_cellranger(
#'   directories_to_load = directories_to_load,
#'   directory_to_write = directory_to_write,
#' )
create_odm_from_cellranger <- function(directories_to_load, directory_to_write, write_cellwise_covariates = TRUE,
                                       chunk_size = 1000L, compression_level = 3L, grna_target_data_frame = NULL,
                                       compute_cell_cycle = FALSE) {
  
  # Cell cycle scoring parameters (Seurat defaults)
  s_genes <- NULL
  g2m_genes <- NULL
  cc_ctrl_genes <- 100L  # Note: may need reduction for small datasets (100 * n_target_genes control genes needed)
  cc_nbin <- 24L
  cc_scale_factor <- 10000
  cc_seed <- 1L
  # 0. check that directory to write is valid; create it if it does not yet exist
  directory_to_write <- expand_tilde(directory_to_write)
  if (is.null(directory_to_write)) {
    stop("`directory_to_write` cannot be `NULL`.")
  }
  if (!dir.exists(directory_to_write)) dir.create(path = directory_to_write, recursive = TRUE)

  # 1. check that directories exist
  for (curr_directory in directories_to_load) {
    if (!dir.exists(curr_directory)) stop("The directory ", curr_directory, " does not exist.")
  }

  # 2. create the list of features and matrix files (exclude cell barcodes)
  input_files <- lapply(X = directories_to_load, function(curr_directory) {
    fs <- list.files(curr_directory)
    grep_strs <- c("*features.tsv($|.gz)", "*matrix.mtx($|.gz)")
    out <- vapply(grep_strs, function(grep_str) {
      file_names <- grep(pattern = grep_str, x = fs, value = TRUE)
      if (length(file_names) >= 2L) {
        stop("There are multiple ", grep_str, " files within the directory ", curr_directory, ".")
      }
      if (length(file_names) == 0L) {
        stop("The directory ", curr_directory, " contains zero ", grep_str, " files.")
      }
      return(paste0(curr_directory, "/", file_names))
    }, FUN.VALUE = character(1)) |> stats::setNames(c("features", "matrix"))
  }) |> stats::setNames(NULL)
  matrix_fps <- vapply(X = input_files, FUN = function(i) i[["matrix"]], FUN.VALUE = character(1))

  # 3. obtain the feature data frame
  feature_df <- data.table::fread(file = input_files[[1]][["features"]],
                                  colClasses = c("character", "character", "character"),
                                  col.names = c("feature_id", "feature_name", "modality"), header = FALSE)
  modality_names <- unique(feature_df$modality)
  if (!(all(modality_names %in% c("Gene Expression", "CRISPR Guide Capture", "Antibody Capture")))) {
    stop("The modality names must be 'Gene Expression' or 'CRISPR Guide Capture'.")
  }

  # 4. determine the start idx of each modality in the feature df
  modality_start_idx_features <- vapply(modality_names, function(modality_name) {
    feature_df[list(modality_name), which = TRUE, mult = "first", on = "modality"] - 1L
  }, FUN.VALUE = integer(1))
  modality_start_idx_features <- c(modality_start_idx_features, nrow(feature_df))

  # 5. obtain the feature ids of each modality
  feature_ids <- feature_df$feature_id
  modality_feature_ids <- lapply(seq_along(modality_names), function(k) {
    curr_modality_features <- feature_ids[seq(modality_start_idx_features[k], modality_start_idx_features[k + 1L] - 1L) + 1L]
  }) |> stats::setNames(modality_names)

  # 6. obtain the feature names of each modality
  feature_names <- feature_df$feature_name
  modality_feature_names <- lapply(seq_along(modality_names), function(k) {
    curr_modality_features <- feature_names[seq(modality_start_idx_features[k], modality_start_idx_features[k + 1L] - 1L) + 1L]
  }) |> stats::setNames(modality_names)

  # 7. determine idxs of MT- features
  mt_feature_idxs <- lapply(seq_along(modality_names), function(k) {
    curr_modality_features <- feature_names[seq(modality_start_idx_features[k], modality_start_idx_features[k + 1L] - 1L) + 1L]
    grep(pattern = "^MT-", x = curr_modality_features, ignore.case = TRUE) - 1L # - modality_start_idx_features[k]
  }) |> stats::setNames(modality_names)

  # 7.5. Compute n_features_per_modality first (needed for cell cycle setup)
  n_features_per_modality <- diff(modality_start_idx_features)

  # 7.6. prepare cell cycle scoring parameters if requested
  cc_gene_indices <- NULL
  gene_expression_sum <- NULL
  gene_expression_count <- NULL
  
  if (compute_cell_cycle) {
    # Check if Gene Expression modality is present
    if (!"Gene Expression" %in% modality_names) {
      stop("Cell cycle scoring requires Gene Expression modality data.")
    }
    
    # Use default gene sets if not provided
    if (is.null(s_genes) || is.null(g2m_genes)) {
      if (!exists("cc_genes_updated_2019")) {
        data("cc_genes_updated_2019", package = "ondisc", envir = environment())
      }
      if (is.null(s_genes)) s_genes <- cc_genes_updated_2019$s_genes
      if (is.null(g2m_genes)) g2m_genes <- cc_genes_updated_2019$g2m_genes
    }
    
    # Validate gene lists
    if (length(s_genes) < 3 || length(g2m_genes) < 3) {
      stop("Need at least 3 genes in each cell cycle gene set")
    }
    
    # Prepare gene indices for Gene Expression modality only
    ge_idx <- which(modality_names == "Gene Expression")
    ge_feature_ids <- modality_feature_ids[[ge_idx]]
    ge_feature_names <- modality_feature_names[[ge_idx]]
    
    cc_gene_indices <- prepare_cell_cycle_gene_indices(
      s_genes = s_genes,
      g2m_genes = g2m_genes,
      feature_ids = ge_feature_ids,
      feature_names = ge_feature_names,
      verbose = TRUE
    )
    
    # Check if we found enough genes
    if (length(cc_gene_indices$s_gene_indices) < 3 || length(cc_gene_indices$g2m_gene_indices) < 3) {
      stop("Insufficient cell cycle genes found in dataset. Need at least 3 genes per phase.")
    }
    
    # Initialize vectors for collecting gene expression statistics during Round 1
    # FIXED: Only allocate space for Gene Expression features, not all features
    ge_idx <- which(modality_names == "Gene Expression")
    n_ge_features <- n_features_per_modality[ge_idx]
    gene_expression_sum <- rep(0.0, n_ge_features)
    gene_expression_count <- rep(0L, n_ge_features)
    
    # Also collect normalized expression statistics for proper control gene binning
    gene_norm_sum <- rep(0.0, n_ge_features)
  }

  # 8 prepare feature_idx to vector_idx map
  feature_idx_to_vector_idx_map <- NULL
  if (!is.null(grna_target_data_frame)) {
    feature_idx_to_vector_idx_map <- generate_grna_idx_to_vector_idx_map(grna_target_data_frame = grna_target_data_frame,
                                                                         modality_start_idx_features = modality_start_idx_features,
                                                                         ordered_grna_ids = modality_feature_names[["CRISPR Guide Capture"]])
    modality_feature_ids[["CRISPR Guide Capture"]] <- unique(grna_target_data_frame$vector_id)
    n_features_per_modality[modality_names == "CRISPR Guide Capture"] <- length(unique(grna_target_data_frame$vector_id))
  }

  # 9. round 1
  round_1_out <- process_input_files_round_1(matrix_fps = matrix_fps,
                                             modality_names = modality_names,
                                             modality_start_idx_features = modality_start_idx_features,
                                             feature_idx_to_vector_idx_map = feature_idx_to_vector_idx_map,
                                             n_features_per_modality = n_features_per_modality,
                                             compute_cell_cycle = compute_cell_cycle,
                                             gene_expression_sum = gene_expression_sum,
                                             gene_expression_count = gene_expression_count,
                                             gene_norm_sum = gene_norm_sum,
                                             cc_scale_factor = cc_scale_factor)

  # 9.5. compute cell cycle control genes if requested
  cc_control_genes <- NULL
  if (compute_cell_cycle) {
    # Use normalized gene means for proper binning (matching Seurat)
    # Divide by total number of cells (including zeros), not just non-zero cells
    normalized_gene_means <- gene_norm_sum / round_1_out$n_cells
    
    ge_gene_means <- normalized_gene_means
    
    # Adjust indices to be relative to Gene Expression modality (0-based)
    cc_gene_indices_adjusted <- list(
      s_gene_indices = cc_gene_indices$s_gene_indices,
      g2m_gene_indices = cc_gene_indices$g2m_gene_indices
    )
    
    cc_control_genes <- compute_cell_cycle_control_genes(
      gene_means = ge_gene_means,
      s_gene_indices = cc_gene_indices_adjusted$s_gene_indices,
      g2m_gene_indices = cc_gene_indices_adjusted$g2m_gene_indices,
      nbin = cc_nbin,
      ctrl = cc_ctrl_genes,
      seed = cc_seed
    )
    
    cat("Cell cycle control gene selection:\n")
    cat("  S control genes:", length(cc_control_genes$s_control_indices), "\n")
    cat("  G2M control genes:", length(cc_control_genes$g2m_control_indices), "\n")
  }

  # 10. initialize cellwise covariates
  cellwise_covariates <- initialize_cellwise_covariates(modality_names = modality_names,
                                                        n_cells = round_1_out$n_cells,
                                                        compute_cell_cycle = compute_cell_cycle)

  # 11. obtain file paths to odms
  new_modality_names <- update_modality_names(modality_names)
  file_names_in <- vapply(new_modality_names, function(new_modality_name) {
    if (!dir.exists(directory_to_write)) dir.create(directory_to_write, recursive = TRUE)
    paste0(directory_to_write, "/", new_modality_name, ".odm")
  }, FUN.VALUE = character(1))

  # 12. sample integer id
  integer_id <- sample(x = seq(0L, .Machine$integer.max), size = 1L)

  # 13. initialize odms
  row_ptr_list <- initialize_odms(modality_names = modality_names,
                                  file_names_in = file_names_in,
                                  n_nonzero_features_vector_list = round_1_out$n_nonzero_features_vector_list,
                                  modality_feature_ids = modality_feature_ids,
                                  n_cells = round_1_out$n_cells,
                                  integer_id = integer_id,
                                  chunk_size = chunk_size,
                                  compression_level = compression_level)

  # 14. round 2
  process_input_files_round_2(matrix_fps = matrix_fps,
                              file_names_in = file_names_in,
                              modality_names = modality_names,
                              modality_start_idx_features = modality_start_idx_features,
                              row_ptr_list = row_ptr_list,
                              modality_start_idx_mtx_list = round_1_out$modality_start_idx_mtx_list,
                              mt_feature_idxs = mt_feature_idxs,
                              cellwise_covariates = cellwise_covariates,
                              feature_idx_to_vector_idx_map = feature_idx_to_vector_idx_map,
                              compute_cell_cycle = compute_cell_cycle,
                              cc_gene_indices = cc_gene_indices,
                              cc_control_genes = cc_control_genes,
                              cc_scale_factor = cc_scale_factor)
  gc() |> invisible()

  # 14.5. assign cell cycle phases if cell cycle scoring was performed
  if (compute_cell_cycle && "Gene Expression" %in% modality_names) {
    ge_idx <- which(modality_names == "Gene Expression")
    s_scores <- cellwise_covariates[[ge_idx]]$covariate_list$s_score
    g2m_scores <- cellwise_covariates[[ge_idx]]$covariate_list$g2m_score
    phases <- assign_cell_cycle_phase(s_scores, g2m_scores)
    # Add phase as a new covariate
    cellwise_covariates[[ge_idx]]$covariate_list$phase <- phases
    cellwise_covariates[[ge_idx]]$bool_vect[["phase"]] <- TRUE
  }

  # 15. save the covariates
  dt <- preprare_output_covariate_dt(cellwise_covariates = cellwise_covariates,
                                     new_modality_names = new_modality_names,
                                     n_cells_per_batch = round_1_out$n_cells_per_batch,
                                     modality_feature_ids = modality_feature_ids)
  if (write_cellwise_covariates) saveRDS(dt, file = paste0(directory_to_write, "/cellwise_covariates.rds"))

  # 16. return the odms and cell covariates
  l <- lapply(seq_along(modality_names), function(i) {
    initialize_odm_from_backing_file(file_names_in[i])
  }) |> stats::setNames(new_modality_names)
  l[["cellwise_covariates"]] <- dt

  return(l)
}
