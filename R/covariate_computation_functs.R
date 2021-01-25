#' Covariate Enum
#' Function to return a simple "enum" of all symbols (terminal and nonterminal) used in the grammar; helps to reduce number of string literals in the code.
#' @return an environment
symbols_enum <- function() {
  covariates <- c("mean_expression_feature",
                  "coef_of_variation_feature",
                  "n_nonzero_feature",
                  "n_nonzero_cell",
                  "n_umis_cell",
                  "p_mito_cell",
                  "mean_sq_expression_feature",
                  "n_mito_cell",
                  "sd_expression_feature")
  list2env(setNames(as.list(covariates), covariates))
}


#' Initialize grammar
#'
#' Initializes the context-free grammar used to compute the covariate matrices.
#'
#' @return an environment representing the context-free grammar.
initialize_grammar <- function() {
  sym_enum <- symbols_enum()
  e <- new.env()
  terminal_entry <- list(terminal = TRUE, value = NULL)

  # Define terminal symbols
  e[[sym_enum$mean_expression_feature]] <- terminal_entry
  e[[sym_enum$mean_sq_expression_feature]] <- terminal_entry
  e[[sym_enum$n_nonzero_feature]] <- terminal_entry
  e[[sym_enum$n_nonzero_cell]] <- terminal_entry
  e[[sym_enum$n_umis_cell]] <- terminal_entry
  e[[sym_enum$n_mito_cell]] <- terminal_entry

  # Define nonterminal symbols, along with their production rules
  e[[sym_enum$p_mito_cell]] <- list(terminal = FALSE,
                                    f = function(p1, p2) p1 / p2,
                                    symbols = c(sym_enum$n_mito_cell,
                                                sym_enum$n_nonzero_cell))
  e[[sym_enum$sd_expression_feature]] <- list(terminal = FALSE,
                                              f = function(p1, p2) p1 - p2^2,
                                              symbols = c(sym_enum$mean_sq_expression_feature,
                                                          sym_enum$mean_expression_feature))
  e[[sym_enum$coef_of_variation_feature]] <- list(terminal = FALSE,
                                                  f = function(p1, p2) p1 / p2,
                                                  symbols = c(sym_enum$sd_expression_feature,
                                                              sym_enum$mean_expression_feature))

  # Finally, store the symbol enum for convenience
  e[["sym_enum"]] <- sym_enum
  return(e)
}


#' Map inputs to covariates
#'
#' A function that maps the input metadata to the covariates to compute.
#'
#' @param mtx_metadata the mtx metadata
#' @param features_metadata the features.tsv metadata
#'
#' @return a list of two entries: feature_covariates and cell_covariates; each entry lists the covariates to compute
map_inputs_to_covariates <- function(mtx_metadata, features_metadata) {
  # Obtain enum
  sym_enum <- symbol_enum()
  # First, get feature covariates
  feature_covariates <- if (!mtx_metadata$is_logical) {
    # integer matrix
    c(sym_enum$mean_expression_feature,
      sym_enum$coef_of_variation_feature,
      sym_enum$n_nonzero_feature)
  } else {
    # logical matrix
    sym_enum$n_nonzero_feature
  }
  # Next, get cell covariates
  cell_covariates <- if (!mtx_metadata$is_logical) {
    if (features_metadata$mt_genes_present) {
      # integer matrix, mt genes present
      c(sym_enum$n_nonzero_cell,
        sym_enum$n_umis_cell,
        sym_enum$p_mito_cell)
    } else {
      # integer matrix, mt genes absent
      c(sym_enum$n_nonzero_cell,
        sym_enum$n_umis_cell)
    }
  } else {
    # logical matrix
    sym_enum$n_nonzero_cell
  }
  return(list(feature_covariates = feature_covariates, cell_covariates = cell_covariates))
}


#' Get nonterminals for covariate
#'
#' Obtains the nonterminals needed to compute a given covariate.
#'
#' @param covariate name of covariate
#' @param grammar a grammar, as initialized by initialize_grammar.
#'
#' @return character vector of nonterminals
get_nonterminals_for_covariate <- function(covariate, grammar) {
  if (grammar[[covariate]]$terminal) {
    ret <- covariate
  } else {
    ret <- lapply(grammar[[covariate]]$symbols,
                  function(s) get_nonterminals_for_covariate(s, grammar))
    ret <- unique(unlist(ret))
  }
  return(ret)
}


