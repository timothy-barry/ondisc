if (FALSE) {

  # PBMC data
  raw_data_dir <- "/Users/timbarry/Box/onDisc_all/onDisc_offsite/raw_data/filtered_feature_bc_matrix"
  mtx_fp <- paste0(raw_data_dir, "/matrix.mtx")
  barcodes_fp <- paste0(raw_data_dir, "/barcodes.tsv")
  # three cols
  features_fp <- paste0(raw_data_dir, "/features.tsv")
  exp_mat <- ondisc_matrix(h5_file = "/Users/timbarry/Box/onDisc_all/onDisc_offsite/raw_data/filtered_feature_bc_matrix/ondisc_matrix_1.h5")

  # one col
  features_fp <- paste0(raw_data_dir, "/features_small.tsv")
  # two cols
  features_fp <- paste0(raw_data_dir, "/features_medium.tsv")

  # Logical gRNA data
  raw_data_dir <- "/Volumes/tims_new_drive/research/sceptre_example/raw_data"
  mtx_fp <- paste0(raw_data_dir, "/perturbations.mtx")
  features_fp <- paste0(raw_data_dir, "/gRNAs.tsv")
  barcodes_fp <- paste0(raw_data_dir, "/cell_barcodes.tsv")

  # Snare-seq data
  raw_data_dir <- "/Users/timbarry/Box/onDisc_all/onDisc_offsite/snare-seq"
  mtx_fp <- paste0(raw_data_dir, "/GSE126074_AdBrainCortex_SNAREseq_chromatin.counts.mtx")
  features_fp <- paste0(raw_data_dir, "/GSE126074_AdBrainCortex_SNAREseq_chromatin.peaks.tsv")
  barcodes_fp <- paste0(raw_data_dir, "/GSE126074_AdBrainCortex_SNAREseq_chromatin.barcodes.tsv")

  # ex code
  e2 <- exp_mat[1:10,]
  e3 <- e2[c("ENSG00000237613", "ENSG00000286448", "ENSG00000243485"),]
}
