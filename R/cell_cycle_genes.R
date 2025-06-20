#' Cell cycle genes: 2019 update
#'
#' A list of genes used in cell-cycle regression, updated with 2019 symbols.
#' This matches Seurat's cc.genes.updated.2019 dataset exactly.
#'
#' @format A list with two components:
#' \describe{
#'   \item{s_genes}{43 genes associated with S-phase}
#'   \item{g2m_genes}{54 genes associated with G2M-phase}
#' }
#'
#' @section Updated symbols:
#' The following symbols were updated from the original gene sets:
#' \describe{
#'   \item{s.genes}{
#'     \itemize{
#'       \item \emph{MCM2}: \emph{MCM7}
#'       \item \emph{MLF1IP}: \emph{CENPU}
#'       \item \emph{RPA2}: \emph{POLR1B}
#'       \item \emph{BRIP1}: \emph{MRPL36}
#'     }
#'   }
#'   \item{g2m.genes}{
#'     \itemize{
#'       \item \emph{FAM64A}: \emph{PIMREG}
#'       \item \emph{HN1}: \emph{JPT1}
#'     }
#'   }
#' }
#'
#' @source \url{https://www.science.org/doi/abs/10.1126/science.aad0501}
#' @concept data
#' @examples
#' # Load the cell cycle gene sets
#' data("cc_genes_updated_2019", package = "ondisc")
#' length(cc_genes_updated_2019$s_genes)   # 43 S-phase genes
#' length(cc_genes_updated_2019$g2m_genes) # 54 G2M-phase genes
#'
#' @export
"cc_genes_updated_2019"

# Create the actual data object
cc_genes_updated_2019 <- list(
  s_genes = c("MCM5", "PCNA", "TYMS", "FEN1", "MCM7", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "CENPU", "HELLS", "RFC2", "POLR1B", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "MRPL36", "E2F8"),
  g2m_genes = c("HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "PIMREG", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "JPT1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA")
)

#' Validate and prepare cell cycle gene sets
#'
#' Internal function to validate user-provided gene sets and map them to
#' feature indices for efficient computation during import.
#'
#' @param s_genes character vector of S phase genes
#' @param g2m_genes character vector of G2/M phase genes
#' @param feature_ids character vector of all feature IDs in the dataset
#' @param feature_names character vector of all feature names in the dataset
#' @param verbose logical, whether to print information about gene mapping
#'
#' @return A list containing:
#'   \describe{
#'     \item{s_gene_indices}{Integer vector of S gene indices (0-based)}
#'     \item{g2m_gene_indices}{Integer vector of G2/M gene indices (0-based)}
#'     \item{n_s_genes_found}{Number of S genes found in dataset}
#'     \item{n_g2m_genes_found}{Number of G2/M genes found in dataset}
#'     \item{missing_s_genes}{S genes not found in dataset}
#'     \item{missing_g2m_genes}{G2/M genes not found in dataset}
#'   }
#'
#' @keywords internal
prepare_cell_cycle_gene_indices <- function(s_genes, g2m_genes, feature_ids, feature_names, verbose = TRUE) {
  # Try to match genes by both feature_ids and feature_names
  find_gene_indices <- function(genes, feature_ids, feature_names) {
    # First try feature_names (more common for cell cycle genes)
    name_matches <- match(genes, feature_names)
    # Then try feature_ids as backup
    id_matches <- match(genes, feature_ids)

    # Combine matches, preferring feature_name matches
    indices <- ifelse(!is.na(name_matches), name_matches, id_matches)
    found_genes <- genes[!is.na(indices)]
    missing_genes <- genes[is.na(indices)]

    # Convert to 0-based indexing for C++
    indices_0based <- indices[!is.na(indices)] - 1L

    return(list(
      indices = indices_0based,
      found_genes = found_genes,
      missing_genes = missing_genes
    ))
  }

  # Process S phase genes
  s_result <- find_gene_indices(s_genes, feature_ids, feature_names)

  # Process G2/M phase genes
  g2m_result <- find_gene_indices(g2m_genes, feature_ids, feature_names)

  if (verbose) {
    cat("Cell cycle gene mapping:\n")
    cat("  S phase genes: ", length(s_result$found_genes), "/", length(s_genes), " found\n")
    cat("  G2/M phase genes: ", length(g2m_result$found_genes), "/", length(g2m_genes), " found\n")

    if (length(s_result$missing_genes) > 0) {
      cat("  Missing S genes: ", paste(utils::head(s_result$missing_genes, 5), collapse = ", "))
      if (length(s_result$missing_genes) > 5) cat(", ...")
      cat("\n")
    }

    if (length(g2m_result$missing_genes) > 0) {
      cat("  Missing G2/M genes: ", paste(utils::head(g2m_result$missing_genes, 5), collapse = ", "))
      if (length(g2m_result$missing_genes) > 5) cat(", ...")
      cat("\n")
    }
  }

  # Check if we have enough genes for scoring
  if (length(s_result$indices) < 5) {
    warning("Only ", length(s_result$indices), " S phase genes found. Cell cycle scoring may be unreliable.")
  }
  if (length(g2m_result$indices) < 5) {
    warning("Only ", length(g2m_result$indices), " G2/M phase genes found. Cell cycle scoring may be unreliable.")
  }

  return(list(
    s_gene_indices = s_result$indices,
    g2m_gene_indices = g2m_result$indices,
    n_s_genes_found = length(s_result$found_genes),
    n_g2m_genes_found = length(g2m_result$found_genes),
    missing_s_genes = s_result$missing_genes,
    missing_g2m_genes = g2m_result$missing_genes
  ))
}

#' Compute control gene selection for cell cycle scoring
#'
#' @param gene_means Numeric vector of mean expression per gene
#' @param s_gene_indices Integer vector of S gene indices (0-based)
#' @param g2m_gene_indices Integer vector of G2M gene indices (0-based)
#' @param nbin Number of expression bins (default 24)
#' @param ctrl Number of control genes per target gene (default 100)
#' @param seed Random seed for reproducibility
#'
#' @return List with s_control_indices and g2m_control_indices (0-based)
#' @keywords internal
compute_cell_cycle_control_genes <- function(
  gene_means,
  s_gene_indices,
  g2m_gene_indices,
  nbin = 24L,
  ctrl = 100L,
  seed = 1L
) {
  set.seed(seed)

  # Match Seurat's exact approach
  # 1. Order gene means (Seurat: data.avg <- data.avg[order(data.avg)])
  data.avg <- gene_means[order(gene_means)]

  # 2. Use cut_number with noise exactly like Seurat
  # Seurat: cut_number(x = data.avg + rnorm(n = length(data.avg))/1e30, n = nbin, labels = FALSE, right = FALSE)
  data.cut <- ggplot2::cut_number(x = data.avg + stats::rnorm(n = length(data.avg))/1e30, n = nbin, labels = FALSE, right = FALSE)
  names(data.cut) <- names(data.avg)

  # 3. Convert back to original gene indices
  # Create mapping from ordered indices back to original 0-based indices
  order_to_original <- order(gene_means) - 1  # Convert to 0-based

  # Helper function matching Seurat's exact logic
  sample_control_genes_seurat <- function(target_indices) {
    ctrl.use <- c()
    all_target_indices <- c(s_gene_indices, g2m_gene_indices)  # All target genes to exclude

    for (target_idx in target_indices) {
      # Find this gene in the ordered list
      target_gene_ordered_pos <- which(order_to_original == target_idx)
      if (length(target_gene_ordered_pos) > 0) {
        target_bin <- data.cut[target_gene_ordered_pos]
        # Find all genes in the same bin
        same_bin_genes_ordered <- which(data.cut == target_bin)
        # Convert back to original 0-based indices
        same_bin_genes_original <- order_to_original[same_bin_genes_ordered]

        # EXCLUDE target genes from potential controls (like Seurat does!)
        potential_controls <- setdiff(same_bin_genes_original, all_target_indices)

        # Sample ctrl genes from remaining candidates
        if (length(potential_controls) >= ctrl) {
          sampled <- sample(potential_controls, ctrl, replace = FALSE)
        } else {
          sampled <- potential_controls
        }
        ctrl.use <- c(ctrl.use, sampled)
      }
    }
    return(unique(ctrl.use))
  }

  list(
    s_control_indices = sample_control_genes_seurat(s_gene_indices),
    g2m_control_indices = sample_control_genes_seurat(g2m_gene_indices)
  )
}

#' Assign cell cycle phases based on scores
#'
#' @param s_scores Numeric vector of S phase scores
#' @param g2m_scores Numeric vector of G2M phase scores
#'
#' @return Factor vector of phase assignments with levels c("G1", "S", "G2M", "Undecided")
#' @keywords internal
assign_cell_cycle_phase <- function(s_scores, g2m_scores) {
  phases <- character(length(s_scores))

  for (i in seq_along(s_scores)) {
    s <- s_scores[i]
    g2m <- g2m_scores[i]

    if (s < 0 && g2m < 0) {
      phases[i] <- "G1"
    } else if (abs(s - g2m) < 1e-10) {  # Handle floating point equality
      phases[i] <- "Undecided"
    } else if (s > g2m) {
      phases[i] <- "S"
    } else {
      phases[i] <- "G2M"
    }
  }

  # Convert to factor for memory efficiency
  return(factor(phases, levels = c("G1", "S", "G2M", "Undecided")))
}
