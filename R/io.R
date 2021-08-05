#' Save ODM
#'
#' Saves an ondisc matrix or metadata ondisc matrix.
#'
#' @param odm a metadata_ondisc_matrix object or ondisc_matrix object
#' @param name name to label save
#'
#' @return NULL
#' @export
#'
#' @examples
#' tempfile <- create_new_directory()
#' odm_fp <- paste0(tempfile, "/expression_odm")
#' file_locs <- system.file("extdata",package = "ondisc",
#' c("gene_expression.mtx", "genes.tsv", "cell_barcodes.tsv"))
#' names(file_locs) <- c("expressions", "features", "barcodes")
#' expression_data <- create_ondisc_matrix_from_mtx(mtx_fp = file_locs[["expressions"]],
#' barcodes_fp = file_locs[["barcodes"]],
#' features_fp = file_locs[["features"]],
#' odm_fp = odm_fp,
#' return_metadata_ondisc_matrix = TRUE)
#' expression_data <- expression_data[,1:10]
#' save_odm(expression_data, "filtered")
#' expression_data_loaded_unfiltered <- read_odm(odm_fp)
#' expression_data_loaded_filtered <- read_odm(odm_fp, "filtered")
save_odm <- function(odm, name) {
  aux_info <- list()
  if (is(odm, "multimodal_ondisc_matrix")) {
    stop("Save function not yet implemented for metadata_ondisc_matrix class.
         Save components of class individually.")
  }
  if (is(odm, "metadata_ondisc_matrix")) {
    aux_info$cell_covariates <- odm@cell_covariates
    aux_info$feature_covariates <- odm@feature_covariates
    odm <- odm@ondisc_matrix
  }
  # save aux info
  aux_info$cell_subset <- odm@cell_subset
  aux_info$feature_subset <- odm@feature_subset

  # determine file path to save
  append_rds <- !grepl(pattern = "*.rds$", x = name)
  aux_info_fp <- paste0(odm@odm_fp, "/", name, if (append_rds) ".rds" else "")

  # save
  saveRDS(object = aux_info, file = aux_info_fp)
}


#' Read ODM
#' @rdname save_odm
#' @export
#' @param odm_fp file path to an ondisc matrix (or metadata ondisc matrix)
read_odm <- function(odm_fp, name = NULL) {
  h5_file <- paste0(odm_fp, "/data.h5")
  # construct the ondisc_matrix object
  out <- ondisc_matrix(h5_file)
  # if name is non-null, modify the output futher
  if (!is.null(name)) {
    append_rds <- !grepl(pattern = "*.rds$", x = name)
    aux_info_fp <- paste0(odm_fp, "/", name, if (append_rds) ".rds" else "")
    # load aux info rds; update cell and feature idxs
    aux_info <- readRDS(aux_info_fp)
    out@cell_subset <- aux_info$cell_subset
    out@feature_subset <- aux_info$feature_subset
    # if available, load cell-specific and feature-specific covariates; return metadata_odm
    if ("cell_covariates" %in% names(aux_info) && "feature_covariates" %in% names(aux_info)) {
      out <- metadata_ondisc_matrix(ondisc_matrix = out,
                                    cell_covariates = aux_info$cell_covariates,
                                    feature_covariates = aux_info$feature_covariates)
    }
  }
  return(out)
}
