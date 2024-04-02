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


generate_grna_idx_to_vector_idx_map <- function(grna_target_data_frame, modality_start_idx_features, ordered_grna_ids) {
  grna_target_data_frame <- grna_target_data_frame |> dplyr::select(grna_id, vector_id)
  unique_vectors <- unique(grna_target_data_frame$vector_id)
  vector_to_vector_idx <- data.frame(vector_id = unique_vectors,
                                     vector_idx = seq(from = modality_start_idx_features[["CRISPR Guide Capture"]],
                                                      length.out = length(unique_vectors)))
  grna_to_grna_idx <- data.frame(grna_id = ordered_grna_ids,
                                 feature_idx = seq(modality_start_idx_features[["CRISPR Guide Capture"]],
                                                   length.out = length(ordered_grna_ids)))
  grna_idx_to_vector_idx <- dplyr::left_join(x = dplyr::left_join(x = grna_to_grna_idx,
                                                                  y = grna_target_data_frame,
                                                                  by = "grna_id"),
                                             y = vector_to_vector_idx, by = "vector_id") |>
    dplyr::select(feature_idx, vector_idx)
  return(grna_idx_to_vector_idx)
}


collapse_grna_counts <- function(dt, feature_idx_to_vector_idx_map, modality_start_mtx, round_1) {
  # extract portion of data table corresponding to grna
  start <- modality_start_mtx[["CRISPR Guide Capture"]][1]
  end <- modality_start_mtx[["CRISPR Guide Capture"]][2]
  dt_sub <- dt[seq(start, end) + 1L,]
  # map grna id to vector id
  dt_sub_updated <- dplyr::left_join(dt_sub, feature_idx_to_vector_idx_map, by = "feature_idx") |>
    dplyr::select(feature_idx = vector_idx, j, if (!round_1) "x" else NULL)
  # sort on feature_idx and j
  data.table::setorderv(dt_sub_updated, c("feature_idx", "j"))
  if (round_1) {
    # if round 1, collapse feature_idx
    collapsed_dt_sub <- dt_sub_updated |>
      dplyr::group_by(feature_idx) |>
      dplyr::summarize(unique_count = length(unique(j)))
    feature_idx_new <- rep(collapsed_dt_sub$feature_idx, times = collapsed_dt_sub$unique_count)
    update_dt_column(dt$feature_idx, overwrite_vector = feature_idx_new, start = start)
  } else {
    # if round 2, collapse j, x
    collapsed_dt_sub <- dt_sub_updated |>
      dplyr::group_by(feature_idx, j) |>
      dplyr::summarize(s = if (dplyr::n() == 1) x else sum(x), .groups = "drop")
    feature_idx_new <- collapsed_dt_sub$feature_idx
    update_dt_column(dt$feature_idx, overwrite_vector = feature_idx_new, start = start)
    update_dt_column(dt$j, overwrite_vector = collapsed_dt_sub$j, start = start)
    update_dt_column(dt$x, overwrite_vector = collapsed_dt_sub$s, start = start)
  }
  # return updated modality_start_mtx
  new_end <- start + length(feature_idx_new) - 1L
  modality_start_mtx[["CRISPR Guide Capture"]] <- c(start, new_end)
  return(modality_start_mtx)
}
