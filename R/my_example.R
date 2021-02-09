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
                                         return_covariate_ondisc_matrix = TRUE)
pert_mat <- create_ondisc_matrix_from_mtx(mtx_fp = perturbation_mtx,
                                          barcodes_fp = barcodes,
                                          features_fp = guides_tsv,
                                          return_covariate_ondisc_matrix = TRUE,
                                          on_disc_dir = temp_dir)
}
