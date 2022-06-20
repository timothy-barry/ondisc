test_that("convert assign list to sparse odm", {
  n_cells <- 50
  n_gRNAs <- 40
  cell_barcodes <- paste0("cell_", seq(1, n_cells))
  gRNA_ids <- paste0("gRNA_", seq(1, n_gRNAs))
  gRNA_assignment_list <- replicate(n_cells,
                                    sample(x = gRNA_ids, size = sample(x = seq(0, 3), size = 1), replace = FALSE),
                                    FALSE) |> lapply(FUN = function(x) if (length(x) == 0) "" else x)
  names(gRNA_assignment_list) <- cell_barcodes
  features_metadata_df <- data.frame(target_type = sample(x = c("gene", "non-targeting"), size = n_gRNAs, replace = TRUE),
                                     target = paste0("target", seq(1, n_gRNAs)))
  odm_fp <- paste0(tempdir(), "/matrix.odm")
  metadata_fp <- paste0(tempdir(), "/metadata.rds")
  if (file.exists(odm_fp)) file.remove(odm_fp)
  if (file.exists(metadata_fp)) file.remove(metadata_fp)
  odm <- convert_assign_list_to_sparse_odm(cell_barcodes, gRNA_ids, gRNA_assignment_list, odm_fp, metadata_fp, features_metadata_df)
  for (cell_barcode in sample(cell_barcodes, 4)) {
    for (gRNA_id in sample(gRNA_ids, 4)) {
      expect_equal(gRNA_id %in% gRNA_assignment_list[[cell_barcode]],
                   as.logical(odm[[gRNA_id, cell_barcode]]))
    }
  }
})
