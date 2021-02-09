## Example, small, synthetic expression matrix in .mtx format
## Note: it is assumed that the PBMC data have been downloaded from 10X; see the "Using ondisc" vignette for code to do this.

library(Matrix)
library(tidyverse)
data_dir <- "/Users/timbarry/Box/onDisc_all/onDisc_offsite/raw_data" # change this file path to reproduce!
package_dir <- "/Users/timbarry/Box/onDisc_all/ondisc/" # change this file path as well!

# Load package via load_all (to access non-exported functions)
load_all(path = package_dir, helpers = FALSE)
ext_data_dir <- paste0(package_dir, "inst/extdata")

# Load barcodes and gene ids, names
sub_data_dir <- paste0(data_dir, "/filtered_feature_bc_matrix")
all_barcodes <- read_tsv(file = paste0(sub_data_dir, "/barcodes.tsv"), col_types = "c", col_names = FALSE) %>% pull()
all_gene_ids_and_names <- read_tsv(file = paste0(sub_data_dir, "/features.tsv"), col_types = "ccc", col_names = c("gene_id", "gene_name", "feature_type"))

# Determine where the MT genes are
mt_locs <- grep(pattern = "^MT-", x = all_gene_ids_and_names$gene_name)
non_mt_locs <- (1:nrow(all_gene_ids_and_names))[-mt_locs]

# sample the genes and cells
n_genes <- 300
n_cells <- 900
non_mt_idxs <- sample(x = non_mt_locs, size = n_genes - length(mt_locs))
final_idxs <- sample(c(non_mt_idxs, mt_locs))
my_cell_barcodes <- sample(x = all_barcodes, size = n_cells, replace = FALSE)
my_gene_names <- all_gene_ids_and_names$gene_name[final_idxs]
my_gene_ids <- all_gene_ids_and_names$gene_id[final_idxs]

# generate the random matrix
m <- create_random_matrix(n_row = 300, n_col = 900, p_zero = 0.9, matrix_values = 1:11)

# save the gene expression matrix
writeMM(obj = m, file = paste0(ext_data_dir, "/gene_expression.mtx"))
write_tsv(x = tibble(my_cell_barcodes), file = paste0(ext_data_dir, "/cell_barcodes.tsv"), col_names = FALSE)
write_tsv(x = tibble(my_gene_ids, my_gene_names), file = paste0(ext_data_dir, "/genes.tsv"), col_names = FALSE)

# Next, generate a binary matrix of perturbations. Also, generate guide RNA IDs
n_guides <- 250
m <- create_random_matrix(n_row = n_guides, n_col = 900, p_zero = 0.95, matrix_values = 1)
guide_ids <- paste0("guide_", 1:n_guides)
writeMM(obj = m, file = paste0(ext_data_dir, "/perturbation.mtx"))
write_tsv(x = tibble(guide_ids), file = paste0(ext_data_dir, "/guides.tsv"), col_names = FALSE)
