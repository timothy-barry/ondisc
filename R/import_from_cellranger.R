#' Create ODM from cellranger
#'
#' @param directories_to_load directories containing the expression data
#' @param directory_to_write directory in which to write the initialized `.odm` files and `cellwise_covariates.rds` file
#' @param write_cellwise_covariates boolean indicating whether to write the cellwise covariates to disk in `.rds` format
#' @param chunk_size chunk size to use in the backing HDF5 file
#' @param compression_level compression level to use in the backing HDF5 file
#'
#' @return list containing the odm object(s) and cellwise covariates
#' @export
#'
#' @examples
#' \dontrun{
#' base_dir <- "/Users/tib163/research_offsite/external/replogle-2022/raw/rd7/rpe1_other"
#' directories_to_load <- list.files(base_dir, full.names = TRUE)[1:3]
#' directory_to_write <- tempdir()
#' grna_target_data_frame <- readRDS("/Users/tib163/research_offsite/external/replogle-2022/raw/rd7/grna_table.rds")
#' output <- create_odm_from_cellranger(directories_to_load, directory_to_write)
#' }
create_odm_from_cellranger <- function(directories_to_load, directory_to_write, write_cellwise_covariates = TRUE,
                                       chunk_size = 1000L, compression_level = 3L, grna_target_data_frame = NULL) {
  # 0. check that directory to write is valid; create it if it does not yet exist
  directory_to_write <- expand_tilde(directory_to_write)
  if (is.null(directory_to_write)) {
    stop("`directory_to_write` cannot be `NULL`.")
  }
  if (!dir.exists(directory_to_write)) dir.create(path = directory_to_write, recursive = TRUE)

  # 1. check that directories exist
  for (curr_directory in directories_to_load) {
    if (!dir.exists(curr_directory)) stop(paste0("The directory ", curr_directory, " does not exist."))
  }

  # 2. create the list of features and matrix files (exclude cell barcodes)
  input_files <- sapply(X = directories_to_load, function(curr_directory) {
    fs <- list.files(curr_directory)
    grep_strs <- c("*features.tsv($|.gz)", "*matrix.mtx($|.gz)")
    out <- sapply(grep_strs, function(grep_str) {
      file_names <- grep(pattern = grep_str, x = fs, value = TRUE)
      if (length(file_names) >= 2L) {
        stop(paste0("There are multiple ", grep_str, " files within the directory ", curr_directory, "."))
      }
      if (length(file_names) == 0L) {
        stop(paste0("The directory ", curr_directory, " contains zero ", grep_str, " files."))
      }
      return(paste0(curr_directory, "/", file_names))
    }) |> stats::setNames(c("features", "matrix"))
  }, simplify = FALSE) |> stats::setNames(NULL)
  matrix_fps <- sapply(X = input_files, FUN = function(i) i[["matrix"]])

  # 3. obtain the feature data frame
  feature_df <- data.table::fread(file = input_files[[1]][["features"]],
                                  colClasses = c("character", "character", "character"),
                                  col.names = c("feature_id", "feature_name", "modality"), header = FALSE)
  modality_names <- unique(feature_df$modality)

  # 4. determine the start idx of each modality in the feature df
  modality_start_idx_features <- sapply(modality_names, function(modality_name) {
    feature_df[list(modality_name), which = TRUE, mult = "first", on = "modality"] - 1L
  })
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
    grep(pattern = "^MT-", x = curr_modality_features, ignore.case = TRUE) - 1L - modality_start_idx_features[k]
  }) |> stats::setNames(modality_names)

  # 8 prepare feature_idx to vector_idx map
  feature_idx_to_vector_idx_map <- NULL
  n_features_per_modality <- diff(modality_start_idx_features)
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
                                             n_features_per_modality = n_features_per_modality)

  # 10. initialize cellwise covariates
  cellwise_covariates <- initialize_cellwise_covariates(modality_names = modality_names,
                                                        n_cells = round_1_out$n_cells)

  # 11. obtain file paths to odms
  new_modality_names <- update_modality_names(modality_names)
  file_names_in <- sapply(new_modality_names, function(new_modality_name) {
    if (!dir.exists(directory_to_write)) dir.create(directory_to_write, recursive = TRUE)
    paste0(directory_to_write, "/", new_modality_name, ".odm")
  })

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
                              feature_idx_to_vector_idx_map = feature_idx_to_vector_idx_map)
  gc() |> invisible()

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
