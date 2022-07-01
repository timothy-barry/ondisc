#' Create an `ondisc` matrix from an R matrix.
#'
#' Initializes a `covariate_ondisc_matrix` object from an R matrix. Returns a `covariate_ondisc_matrix` including an `ondisc_matrix` alongside cell-specific and feature-specific covariate matrices.
#'
#' This function computes the following cell-specific and feature-specific covariates:
#' - cell-specific: (i) total number of features expressed in cell (n_nonzero_cell), (ii) total UMI count (n_umis_cell), and (iii) percentage of UMIs that map to mitochondrial genes (p_mito_cell).
#' - feature-specific: (i) total number of cells in which feature is expressed (n_nonzero_feature), (ii) mean expression of feature across cells (mean_expression_feature), (iii) coefficient of variation of feature expression across cells (coef_of_variation_feature).
#'
#' @param r_matrix an R matrix. The matrix can be represented as (i) a standard (dense) R matrix (integer or logical), (ii) a sparse integer matrix of class "dgTMatrix," "dgRMatrix," or "dgCMatrix", or (iii) a sparse logical matrix of class "lgTMatrix."
#' @param barcodes a character vector giving the cell barcodes.
#' @param features_df a data frame giving the names of the features. The first column (required) contains the feature IDs (e.g., ENSG00000186092), and the second column (optional) contains the human-readable feature names (e.g., OR4F5). Subsequent columns are discarded. Gene names starting with "MT-" are assumed to be mitochondrial genes and will be used to compute the p_mito covariate.
#' @param odm_fp location to write the backing .odm file.
#' @param metadata_fp (optional; default NULL) location to write the metadata .RDS file. By default, a file called "metadata.rds" stored in the same directory as the backing .odm file.
#' @param feature_access_only boolean value; if TRUE, fast acess to features; if FALSE (default), fast access to features AND cells
#'
#' @return A `covariate_ondisc_matrix` object.
#' @export
#'
#' @examples
#' # Please install the `ondiscdata` package before running the examples.
#' # install.packages("devtools")
#' # devtools::install_github("Katsevich-Lab/ondiscdata")
#'
#' ###########################
#' # EXAMPLE 1: integer counts
#' ###########################
#' file_locs <- system.file("extdata/r_matrix/gene", package = "ondiscdata",
#' c("matrix.rds", "features.rds", "barcodes.rds"))
#' r_matrix <- readRDS(file = file_locs[1])
#' features_df <- readRDS(file = file_locs[2])
#' barcodes <- readRDS(file = file_locs[3])
#' odm_fp <- paste0(create_new_directory(), "/integer_odm")
#' odm_integer <- create_ondisc_matrix_from_R_matrix(r_matrix, barcodes,
#' features_df, odm_fp)
#'
#' ####################
#' # EXAMPLE 2: logical
#' ####################
#' file_locs <- system.file("extdata/r_matrix/gRNA", package = "ondiscdata",
#' c("matrix.rds", "features.rds", "barcodes.rds"))
#' r_matrix_2 <- readRDS(file = file_locs[1])
#' features_df_2 <- readRDS(file = file_locs[2])
#' barcodes <- readRDS(file = file_locs[3])
#' odm_fp <- paste0(create_new_directory(), "/logical_odm")
#' odm_logical <- create_ondisc_matrix_from_R_matrix(r_matrix_2, barcodes,
#' features_df_2, odm_fp)
create_ondisc_matrix_from_R_matrix <- function(r_matrix, barcodes, features_df, odm_fp, metadata_fp = NULL, feature_access_only = FALSE) {
  # generate random ODM ID
  odm_id <- sample(seq(0L, .Machine$integer.max), size = 1)

  # set odm_fp, check if exists
  odm_fp <- append_file_extension(odm_fp, "odm")
  if (file.exists(odm_fp)) stop(paste0("File ", odm_fp, " already exists. Ending function."))

  ### STEP1: compute the cell- and feature- specific covariate matrices
  # Define "bag_of_variables" environment for storing args
  bag_of_variables <- new.env()

  # Extract features and expression metadata
  list2env(get_features_metadata(features_df, bag_of_variables, FALSE), bag_of_variables)
  list2env(get_expression_metadata_from_r_matrix(r_matrix), bag_of_variables)

  # Determine which covariates to compute
  covariates <- map_inputs_to_covariates(bag_of_variables)

  # Define a list of functions to compute the feature-specific covariates
  feature_specific_func_list <- list(n_nonzero_feature = function(x) as.integer(rowSums(x != 0)),
                                     mean_expression_feature = function(x) rowMeans(x),
                                     coef_of_variation_feature = function(x) {
                                       n_cells <- ncol(x)
                                       my_vars <- rowSums(x^2)/n_cells - (rowSums(x)/n_cells)^2
                                       my_means <- rowMeans(x)
                                       return (sqrt(my_vars)/my_means)
                                     })
  feature_covariates <- sapply(X = feature_specific_func_list[covariates$feature_covariates], FUN = function(f) f(r_matrix))
  feature_covariates <- as.data.frame(feature_covariates)
  # rename the column names to be consistent with mtx output
  colnames(feature_covariates) <- gsub('_feature', "", colnames(feature_covariates))

  # Define a list of functions to compute the cell-specific covariates
  cell_specific_func_list <- list(n_nonzero_cell = function(x) as.integer(colSums(x != 0)),
                                  n_umis_cell = function(x) as.integer(colSums(x)),
                                  p_mito_cell = function(x) {
                                    gene_names <- dplyr::pull(features_df, 2)
                                    mt_gene_index <- grep(pattern = "^MT-", x = gene_names)
                                    n_mito_cell = colSums(r_matrix[mt_gene_index,])
                                    return (n_mito_cell/colSums(x))
                                  })
  cell_covariates <- sapply(X = cell_specific_func_list[covariates$cell_covariates], FUN = function(f) f(r_matrix))
  cell_covariates <- as.data.frame(cell_covariates)
  # rename the column names to be consistent with mtx output
  colnames(cell_covariates) <- gsub('_cell', "", colnames(cell_covariates))

  ### STEP2: Generate ondisc_matrix and write the R matrix to disk in CSC and CSR format.
  # Create in-memory CSC and CSR representations of r_matrix
  if (is(r_matrix, "dgTMatrix")) { # sparse triplet form case
    csc_r_matrix <- Matrix::sparseMatrix(i = r_matrix@i,
                                         j = r_matrix@j,
                                         dims = r_matrix@Dim,
                                         repr = "C",
                                         index1 = FALSE,
                                         x = r_matrix@x)
    csr_r_matrix <- Matrix::sparseMatrix(i = r_matrix@i,
                                         j = r_matrix@j,
                                         dims = r_matrix@Dim,
                                         repr = "R",
                                         index1 = FALSE,
                                         x = r_matrix@x)
  } else if (is(r_matrix, "dgRMatrix")) { # CSR format
    r_matrix_t <- Matrix::t(r_matrix)
    csc_r_matrix <- Matrix::sparseMatrix(i = r_matrix_t@j,
                                         p = r_matrix_t@p,
                                         dims = r_matrix@Dim,
                                         repr = "C",
                                         index1 = FALSE,
                                         x = r_matrix_t@x)
    csr_r_matrix <- r_matrix
  } else if (is(r_matrix, "dgCMatrix")) { # CSC format
    csc_r_matrix <- r_matrix
    r_matrix_t <- Matrix::t(r_matrix)
    csr_r_matrix <- Matrix::sparseMatrix(j = r_matrix_t@i,
                                         p = r_matrix_t@p,
                                         dims = r_matrix@Dim,
                                         repr = "R",
                                         index1 = FALSE,
                                         x = r_matrix_t@x)
  } else if (is(r_matrix, "matrix")) { # dense case
    csc_r_matrix <- as(r_matrix, "dgCMatrix")
    csr_r_matrix <- as(r_matrix, "dgRMatrix")
  } else if (is(r_matrix, "lgTMatrix", "")) { # sparse matrix
    # remove all false entries
    true_posits <- r_matrix@x
    r_matrix@i <- r_matrix@i[true_posits]
    r_matrix@j <- r_matrix@j[true_posits]
    r_matrix@x <- r_matrix@x[true_posits]
    # convert to csc and csr formats
    csc_r_matrix <- Matrix::sparseMatrix(i = r_matrix@i,
                                         j = r_matrix@j,
                                         dims = r_matrix@Dim,
                                         repr = "C",
                                         index1 = FALSE,
                                         x = r_matrix@x)
    csr_r_matrix <- Matrix::sparseMatrix(i = r_matrix@i,
                                         j = r_matrix@j,
                                         dims = r_matrix@Dim,
                                         repr = "R",
                                         index1 = FALSE,
                                         x = r_matrix@x)
  } else { #invalid input
    stop("Input matrix must be of class matrix, dgTMatrix, dgCMatrix, dgRMatrix, or lgTMatrix.")
  }

  # initialize the ODM
  initialize_h5_file_on_disk(odm_fp, bag_of_variables, odm_id)

  # Write in memory matrix to the .h5 file on-disk (side-effect)
  write_matrix_to_h5(odm_fp, expression_metadata = bag_of_variables, csc_r_matrix = csc_r_matrix, csr_r_matrix = csr_r_matrix, feature_access_only = feature_access_only)

  ### STEP3: Prepare output
  odm <- ondisc_matrix(h5_file = odm_fp,
                       logical_mat = bag_of_variables$is_logical,
                       underlying_dimension = c(bag_of_variables$n_features, bag_of_variables$n_cells),
                       feature_ids = dplyr::pull(features_df, 1),
                       feature_names = if (bag_of_variables$feature_names) dplyr::pull(features_df, 2) else NA_character_,
                       cell_barcodes = barcodes,
                       odm_id = odm_id,
                       feature_access_only = feature_access_only)

  # initialize the metadata odm
  metadata_odm <- covariate_ondisc_matrix(ondisc_matrix = odm,
                                          cell_covariates = cell_covariates,
                                          feature_covariates = feature_covariates)

  # save the metadata
  if (is.null(metadata_fp)) metadata_fp <- paste0(dirname(odm_fp), "/metadata.rds")
  save_odm(metadata_odm, metadata_fp)

  # return initialized object
  return(metadata_odm)
}


#' Get metadata for r_matrix,an R matrix. The matrix can be either integer or logical.
#'
#' @param r_matrix an R matrix. The matrix can be either integer or logical
#'
#' @return a list with the following entries: (i) n_genes, (ii) n_cells, (iii) the number of data points (i.e., number of entries that are zero), (iv) (TRUE/FALSE) matrix is logical.
#' @noRd
get_expression_metadata_from_r_matrix <- function(r_matrix) {
  n_features <- nrow(r_matrix)
  n_cells <- ncol(r_matrix)
  if (is.logical(r_matrix) || is(r_matrix, "lgTMatrix")) {
    is_logical <- TRUE
    n_data_points <- sum(r_matrix)
  } else {
    is_logical <- FALSE
    n_data_points <- sum(r_matrix != 0)
  }
  return(list(n_features = n_features, n_cells = n_cells, n_data_points = n_data_points,
              is_logical = is_logical))
}

#' Initialize h5 file on-disk and write in memory r matrix to the file
#'
#' Initialize the on-disk portion on an ondisc_matrix, and write matrix to the h5 file.
#'
#' @param odm_fp file path to the .h5 file to be initialized
#' @param expression_metadata metadata of the r_matrix
#' @param bag_of_variables metadata of the features_df
#' @param barcodes a character vector giving the cell barcodes.
#' @param features_df a data frame giving the names of the features. The first column (required) contains the feature IDs (e.g., ENSG00000186092), and the second column (optional) contains the human-readable feature names (e.g., OR4F5). Subsequent columns are discarded. Gene names starting with "MT-" are assumed to be mitochondrial genes and will be used to compute the p_mito covariate.
#' @param csc_r_matrix a Matrix csc representation of the r matrix
#' @param csr_r_matrix a Matrix csr representation of the r matrix
#' @param feature_access_only a boolean value, TRUE if only allow gene-wise access; FALSE if allow both gene-wise and cell-wise access
#'
#' @return NULL
#' @noRd
write_matrix_to_h5 <- function(odm_fp, expression_metadata, csc_r_matrix, csr_r_matrix, feature_access_only) {
  # Write CSC, cell-wise access
  if (!feature_access_only) {
    rhdf5::h5write(csc_r_matrix@p, file = odm_fp, name="cell_ptr")
    rhdf5::h5write(csc_r_matrix@i, file = odm_fp, name="feature_idxs")
    if (!expression_metadata$is_logical) {
      rhdf5::h5write(csc_r_matrix@x, file = odm_fp, name="data_csc")
    }
  }

  # Write CSR, gene-wise access
  rhdf5::h5write(csr_r_matrix@p, file = odm_fp, name = "feature_ptr")
  rhdf5::h5write(csr_r_matrix@j, file = odm_fp, name = "cell_idxs")
  if (!expression_metadata$is_logical) {
    rhdf5::h5write(csr_r_matrix@x, file = odm_fp, name = "data_csr")
  }

  invisible(NULL)
}
