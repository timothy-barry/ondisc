update_modality_names <- function(modality_names) {
  modality_rename <- c("gene" = "Gene Expression", "grna" = "CRISPR Guide Capture")
  match_idx <- match(x = modality_names, table = modality_rename)
  new_modality_names <- sapply(seq_along(modality_names), function(k) {
    if (is.na(match_idx[k])) gsub("\\s+", "_", modality_names[k]) else names(modality_rename)[match_idx[k]]
  })
  return(new_modality_names)
}
