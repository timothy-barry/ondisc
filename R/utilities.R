update_modality_names <- function(modality_names) {
  modality_rename <- c("gene" = "Gene Expression", "grna" = "CRISPR Guide Capture")
  match_idx <- match(x = modality_names, table = modality_rename)
  new_modality_names <- sapply(seq_along(modality_names), function(k) {
    if (is.na(match_idx[k])) gsub("\\s+", "_", modality_names[k]) else names(modality_rename)[match_idx[k]]
  })
  return(new_modality_names)
}

expand_tilde <- function(fp) {
  gsub(pattern = "~", fixed = TRUE, replacement = path.expand("~"), x = fp)
}

collapse_grna_matrix <- function(dt_sub, grna_ids, grna_target_data_frame) {
  # dt_sub <- dt[seq(modality_start_mtx[[2]] + 1, modality_start_mtx[[3]]),]
  # grna_ids <- modality_feature_ids$`CRISPR Guide Capture`
  # modality_start_idx_features

  # current id-to-idx map
  curr_id_to_idx_map <- data.frame(grna_id =  modality_feature_ids[["CRISPR Guide Capture"]],
                                   feature_idx = seq(modality_start_idx_features[[2]], modality_start_idx_features[[3]] - 1L))
  # id to vector map
  id_to_vector_map <- grna_target_data_frame |> dplyr::select(grna_id, vector_id)
  # vector to new idx map
  unique_vectors <- unique(grna_target_data_frame$vector_id)
  vector_to_new_id_map <- data.frame(vector_id = unique(grna_target_data_frame$vector_id),
                                     new_idx = seq(from = modality_start_idx_features[[2]], length.out = length(unique_vectors)))
  # id-to-idx-to-vector map
  map <- dplyr::left_join(x = dplyr::left_join(x = id_to_vector_map,
                                               y = curr_id_to_idx_map,
                                               by = "grna_id"),
                          y = vector_to_new_id_map, by = "vector_id") |>
    dplyr::select(feature_idx, new_idx)

  # create dt_sub_new, the collapsed version of dt_sub
  dt_sub_new <- dplyr::left_join(dt_sub, map, by = "feature_idx") |>
    dplyr::select(-feature_idx) |>
    dplyr::rename(feature_idx = new_idx) |>
    dplyr::relocate(feature_idx, j, x)
  data.table::setorderv(dt_sub_new, cols = c("feature_idx", "j"))
  dt_sub_new_summarized <- dplyr::group_by(dt_sub_new, feature_idx, j) |>
    dplyr::summarize(x = if (length(x) == 1) x else sum(x), .groups = "drop")

  # overwrite the relevant section of dt (write in C++!!!)
  dt[seq(modality_start_mtx[[2]] + 1, modality_start_mtx[[3]]),] <- dt_sub_new

  # update other stuff; modality end must be updated
  #

}
