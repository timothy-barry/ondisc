# Basic unit tests for cell cycle scoring functionality

test_that("Cell cycle gene dataset is loaded correctly", {
  data("cc_genes_updated_2019", package = "ondisc")

  expect_equal(length(cc_genes_updated_2019$s_genes), 43)
  expect_equal(length(cc_genes_updated_2019$g2m_genes), 54)
  expect_true("MCM5" %in% cc_genes_updated_2019$s_genes)
  expect_true("CDK1" %in% cc_genes_updated_2019$g2m_genes)
})

test_that("Cell cycle phase assignment works correctly", {
  s_scores <- c(-0.5, 0.8, -0.2, 0.1)
  g2m_scores <- c(-0.3, 0.2, 0.7, 0.1)

  phases <- ondisc:::assign_cell_cycle_phase(s_scores, g2m_scores)

  # Check that result is a factor
  expect_true(is.factor(phases))
  expect_equal(levels(phases), c("G1", "S", "G2M", "Undecided"))
  
  expect_equal(as.character(phases[1]), "G1")        # Both negative
  expect_equal(as.character(phases[2]), "S")         # S > G2M
  expect_equal(as.character(phases[3]), "G2M")       # G2M > S
  expect_equal(as.character(phases[4]), "Undecided") # S ≈ G2M
})

test_that("Basic cell cycle scoring functionality", {
  # Test that cell cycle scoring can be enabled in import
  set.seed(123)

  temp_input_dir <- tempfile()
  temp_output_dir <- tempfile()

  # Load cell cycle genes to ensure we include some
  data("cc_genes_updated_2019", package = "ondisc")
  s_genes_subset <- head(cc_genes_updated_2019$s_genes, 10)
  g2m_genes_subset <- head(cc_genes_updated_2019$g2m_genes, 10)
  other_genes <- paste0("GENE_", 1:30)
  all_genes <- c(s_genes_subset, g2m_genes_subset, other_genes)

  # Create synthetic dataset with actual cell cycle genes
  out <- write_example_cellranger_dataset(
    n_features = 50,
    n_cells = 20,
    n_batch = 1,
    modalities = "gene",
    directory_to_write = temp_input_dir,
    p_zero = 0.7
  )

  # Modify features file to include our cell cycle genes
  batch_dir <- list.files(temp_input_dir, pattern = "batch_", full.names = TRUE)[1]
  features_file <- file.path(batch_dir, "features.tsv.gz")
  features_df <- readr::read_tsv(features_file, col_names = c("id", "name", "modality"),
                                 show_col_types = FALSE)
  features_df$name <- all_genes

  temp_features_file <- file.path(batch_dir, "features_temp.tsv")
  readr::write_tsv(features_df, temp_features_file, col_names = FALSE)
  file.remove(features_file)
  R.utils::gzip(temp_features_file, destname = features_file)

  # Test without cell cycle
  result_no_cc <- create_odm_from_cellranger(
    directories_to_load = batch_dir,
    directory_to_write = paste0(temp_output_dir, "_no_cc"),
    compute_cell_cycle = FALSE,
    chunk_size = 10L
  )

  # Test with cell cycle
  result_with_cc <- create_odm_from_cellranger(
    directories_to_load = batch_dir,
    directory_to_write = paste0(temp_output_dir, "_with_cc"),
    compute_cell_cycle = TRUE,
    chunk_size = 10L
  )

  # Check that cell cycle columns are added
  no_cc_cols <- colnames(result_no_cc$cellwise_covariates)
  with_cc_cols <- colnames(result_with_cc$cellwise_covariates)

  expect_false("gene_s_score" %in% no_cc_cols)
  expect_true("gene_s_score" %in% with_cc_cols)
  expect_true("gene_g2m_score" %in% with_cc_cols)
  expect_true("gene_phase" %in% with_cc_cols)

  # Check that scores are reasonable
  s_scores <- result_with_cc$cellwise_covariates$gene_s_score
  g2m_scores <- result_with_cc$cellwise_covariates$gene_g2m_score
  phases <- result_with_cc$cellwise_covariates$gene_phase

  expect_true(all(is.finite(s_scores)))
  expect_true(all(is.finite(g2m_scores)))
  expect_true(all(phases %in% c("G1", "S", "G2M", "Undecided")))

  # Clean up
  unlink(temp_input_dir, recursive = TRUE)
  unlink(paste0(temp_output_dir, "_no_cc"), recursive = TRUE)
  unlink(paste0(temp_output_dir, "_with_cc"), recursive = TRUE)
})

test_that("C++ cell cycle scoring logic", {
  # Test C++ scoring with known values
  # 3 cells, 5 genes (2 S genes, 2 controls, 1 other)

  feature_idx <- as.integer(c(0, 2, 4, 1, 3, 0, 1, 2, 3, 4)) # Gene indices
  j <- as.integer(c(0, 0, 0, 1, 1, 2, 2, 2, 2, 2))         # Cell indices
  x <- as.integer(c(5, 3, 1, 4, 2, 3, 3, 1, 1, 2))         # Counts
  n_umis <- as.integer(c(9, 6, 10))                         # Total UMIs per cell

  s_gene_indices <- as.integer(c(0, 1))      # S genes
  g2m_gene_indices <- as.integer(c())        # No G2M genes for simplicity
  s_control_indices <- as.integer(c(2, 3))   # Control genes
  g2m_control_indices <- as.integer(c())

  n_cells <- 3L
  s_scores <- numeric(n_cells)
  g2m_scores <- numeric(n_cells)

  # Call C++ function directly
  ondisc:::compute_cell_cycle_scores(
    s_scores = s_scores,
    g2m_scores = g2m_scores,
    s_gene_indices = s_gene_indices,
    g2m_gene_indices = g2m_gene_indices,
    s_control_indices = s_control_indices,
    g2m_control_indices = g2m_control_indices,
    feature_idx = feature_idx,
    j = j,
    x = x,
    n_umis = n_umis,
    start_idx = 0L,
    end_idx = length(x) - 1L,
    feature_offset = 0L,
    cell_offset = 0L,
    n_cells = n_cells,
    scale_factor = 10000.0
  )

  # Calculate expected scores (mean over ALL genes including zeros)
  cell0_s_mean <- (log1p(5/9*10000) + log1p(0/9*10000)) / 2
  cell0_ctrl_mean <- (log1p(3/9*10000) + log1p(0/9*10000)) / 2
  cell0_score <- cell0_s_mean - cell0_ctrl_mean

  cell1_s_mean <- (log1p(0/6*10000) + log1p(4/6*10000)) / 2
  cell1_ctrl_mean <- (log1p(0/6*10000) + log1p(2/6*10000)) / 2
  cell1_score <- cell1_s_mean - cell1_ctrl_mean

  cell2_s_mean <- (log1p(3/10*10000) + log1p(3/10*10000)) / 2
  cell2_ctrl_mean <- (log1p(1/10*10000) + log1p(1/10*10000)) / 2
  cell2_score <- cell2_s_mean - cell2_ctrl_mean

  expected_scores <- c(cell0_score, cell1_score, cell2_score)

  # Verify C++ produces expected results
  expect_equal(s_scores, expected_scores, tolerance = 1e-4)
})

test_that("Cell cycle scores match Seurat", {
  skip_if_not_installed("Seurat")
  library(Seurat)

  set.seed(456)

  # Create temporary directories
  temp_input_dir <- tempfile()
  temp_output_dir <- tempfile()

  # Load cell cycle genes - use full default sets for realistic comparison
  data("cc_genes_updated_2019", package = "ondisc")
  s_genes_full <- cc_genes_updated_2019$s_genes
  g2m_genes_full <- cc_genes_updated_2019$g2m_genes

  # Create large enough dataset: 97 target genes + enough controls
  # Conservative estimate: 24 bins × 150 genes per bin = 3600 genes
  # This ensures each bin has plenty of non-target genes for sampling
  n_total_genes <- 4000
  other_genes <- paste0("GENE-", 1:(n_total_genes - length(s_genes_full) - length(g2m_genes_full)))
  all_genes <- c(s_genes_full, g2m_genes_full, other_genes)

  # Create synthetic 10x dataset - large enough for Seurat's default ctrl=100
  out <- write_example_cellranger_dataset(
    n_features = n_total_genes,
    n_cells = 120,
    n_batch = 1,
    modalities = "gene",
    directory_to_write = temp_input_dir,
    p_zero = 0.75,
    p_set_col_zero = 0.02,
    p_set_row_zero = 0.02
  )

  # Get the count matrix for Seurat
  original_counts <- out$matrix_list$gene

  # Modify features file to use our cell cycle genes
  batch_dir <- list.files(temp_input_dir, pattern = "batch_", full.names = TRUE)[1]
  features_file <- file.path(batch_dir, "features.tsv.gz")
  features_df <- readr::read_tsv(features_file, col_names = c("id", "name", "modality"),
                                 show_col_types = FALSE)
  features_df$name <- all_genes

  temp_features_file <- file.path(batch_dir, "features_temp.tsv")
  readr::write_tsv(features_df, temp_features_file, col_names = FALSE)
  file.remove(features_file)
  R.utils::gzip(temp_features_file, destname = features_file)

  # Run ondisc cell cycle scoring
  ondisc_result <- create_odm_from_cellranger(
    directories_to_load = batch_dir,
    directory_to_write = temp_output_dir,
    compute_cell_cycle = TRUE
  )

  # Extract ondisc results
  ondisc_s_scores <- ondisc_result$cellwise_covariates$gene_s_score
  ondisc_g2m_scores <- ondisc_result$cellwise_covariates$gene_g2m_score
  ondisc_phases <- ondisc_result$cellwise_covariates$gene_phase

  # Process same data with Seurat
  rownames(original_counts) <- all_genes
  colnames(original_counts) <- paste0("Cell_", 1:ncol(original_counts))

  seurat_obj <- CreateSeuratObject(counts = original_counts)
  seurat_obj <- NormalizeData(seurat_obj)

  seurat_obj <- CellCycleScoring(seurat_obj,
                                s.features = s_genes_full,
                                g2m.features = g2m_genes_full,
                                set.ident = FALSE,
                                random.seed = 1)

  seurat_s_scores <- seurat_obj$S.Score
  seurat_g2m_scores <- seurat_obj$G2M.Score
  seurat_phases <- seurat_obj$Phase

  # Compare ondisc and Seurat results
  s_correlation <- cor(ondisc_s_scores, seurat_s_scores, use = "complete.obs")
  g2m_correlation <- cor(ondisc_g2m_scores, seurat_g2m_scores, use = "complete.obs")
  phase_agreement <- mean(ondisc_phases == seurat_phases, na.rm = TRUE)

  # Clean up
  unlink(temp_input_dir, recursive = TRUE)
  unlink(temp_output_dir, recursive = TRUE)

  # Verify high correlation with Seurat
  expect_true(s_correlation > 0.9,
              info = paste("S correlation should be > 0.9, got:", round(s_correlation, 4)))
  expect_true(g2m_correlation > 0.9,
              info = paste("G2M correlation should be > 0.9, got:", round(g2m_correlation, 4)))
  expect_true(phase_agreement > 0.8,
              info = paste("Phase agreement should be > 0.8, got:", round(phase_agreement, 4)))
})
