##########################
if (FALSE) {
load_all(helpers = FALSE)
temp_dir <- tempdir()
expression_mtx <- system.file("extdata", "gene_expression.mtx", package = "ondisc")
genes_tsv <- system.file("extdata", "genes.tsv", package = "ondisc")
guides_tsv <- system.file("extdata", "guides.tsv", package = "ondisc")
barcodes <- system.file("extdata", "cell_barcodes.tsv", package = "ondisc")
perturbation_mtx <- system.file("extdata", "perturbation.mtx", package = "ondisc")
exp_mat <- create_ondisc_matrix_from_mtx(mtx_fp = expression_mtx,
                                         barcodes_fp = barcodes,
                                         features_fp = genes_tsv,
                                         return_metadata_ondisc_matrix = TRUE,
                                         on_disc_dir = temp_dir)
pert_mat <- create_ondisc_matrix_from_mtx(mtx_fp = perturbation_mtx,
                                          barcodes_fp = barcodes,
                                          features_fp = guides_tsv,
                                          return_metadata_ondisc_matrix = TRUE,
                                          on_disc_dir = temp_dir)
covariate_ondisc_matrix_list <- list(expressions = exp_mat, perturbations = pert_mat)
multimodal_mat <- multimodal_ondisc_matrix(covariate_ondisc_matrix_list)


raw_data_dir <- "/Users/timbarry/Box/onDisc_all/onDisc_offsite/raw_data/filtered_feature_bc_matrix"
mtx_fp <- paste0(raw_data_dir, "/matrix.mtx")
barcodes_fp <- paste0(raw_data_dir, "/barcodes.tsv")
features_fp <- paste0(raw_data_dir, "/features.tsv")

}
