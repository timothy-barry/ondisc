#' Covariate Enum
#' Function to return a simple "enum" of all symbols (terminal and nonterminal) used in the grammar; helps to reduce number of string literals in the code.
#' @return an environment
symbols_enum <- function() {
  symbols <- c("mean_expression_feature",
                "coef_of_variation_feature",
                "n_nonzero_feature",
                "n_nonzero_cell",
                "n_umis_cell",
                "p_mito_cell",
                "mean_sq_expression_feature",
                "n_mito_cell",
                "sd_expression_feature")
  list2env(setNames(as.list(symbols), symbols))
}


#' Accumulator function Enum
#' Function to return a simple enum of the accumulator functions
#' @return
accumulator_functs_enum <- function() {
  accumulator_functs <- c("inc_mean_count",
                          "inc_n_entries",
                          "inc_mean_sq_count",
                          "inc_count",
                          "inc_cell_count_if_feature_condition")
  list2env(setNames(as.list(accumulator_functs), accumulator_functs))
}


#' Arguments Enum
#' Function to return a simple enum of the possible arguments to the accumulator functions
#' @return
arguments_enum <- function() {
  arguments <- c("feature_idxs",
                 "cell_idxs",
                 "umi_counts",
                 "mt_gene_bool",
                 "n_cells",
                 "n_features",
                 "acc_vect")
  list2env(setNames(as.list(arguments), arguments))
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
  sym_enum <- symbols_enum()
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
get_terminals_for_covariate <- function(covariate, grammar) {
  if (grammar[[covariate]]$terminal) {
    ret <- covariate
  } else {
    ret <- lapply(grammar[[covariate]]$symbols,
                  get_terminals_for_covariate, grammar = grammar)
    ret <- unique(unlist(ret))
  }
  return(ret)
}


#' Get terminal accumulator and args
#'
#' Returns the accumulator function and ordered argument set for a given terminal symbol.
#'
#' @param terminal_symbol a terminal symbol
#'
#' @return a list containing (i) the name of the function, and (ii) the ordered argument names.
get_terminal_acc_and_args <- function(terminal_symbol) {
  sym_enum <- symbols_enum()
  acc_enum <- accumulator_functs_enum()
  arg_enum <- arguments_enum()
  # Mean UMI count of feature
  if (terminal_symbol == sym_enum$mean_expression_feature) {
    acc_funct <- acc_enum$inc_mean_count
    acc_args <- c(arg_enum$feature_idxs, arg_enum$umi_counts, arg_enum$n_cells)
    acc_constructor <- numeric
    acc_length <- arg_enum$n_features
  }
  # Mean squared UMI count of feature
  else if (terminal_symbol == sym_enum$mean_sq_expression_feature) {
    acc_funct <- acc_enum$inc_mean_sq_count
    acc_args <- c(arg_enum$feature_idxs, arg_enum$umi_counts, arg_enum$n_cells)
    acc_constructor <- numeric
    acc_length <- arg_enum$n_features
  }
  # N nonzero entries feature
  else if (terminal_symbol == sym_enum$n_nonzero_feature) {
    acc_funct <- acc_enum$inc_n_entries
    acc_args <- c(arg_enum$feature_idxs)
    acc_constructor <- integer
    acc_length <- arg_enum$n_features
  }
  # N nonzero entries cell
  else if (terminal_symbol == sym_enum$n_nonzero_cell) {
    acc_funct <- acc_enum$inc_n_entries
    acc_args <- c(arg_enum$cell_idxs)
    acc_constructor <- integer
    acc_length <- arg_enum$n_cells
  }
  # UMI count cell
  else if (terminal_symbol == sym_enum$n_umis_cell) {
    acc_funct <- acc_enum$inc_count
    acc_args <- c(arg_enum$cell_idxs, arg_enum$umi_counts)
    acc_constructor <- integer
    acc_length <- arg_enum$n_cells
  }
  else if (terminal_symbol == sym_enum$n_mito_cell) {
    acc_funct <- acc_enum$inc_cell_count_if_feature_condition
    acc_args <- c(arg_enum$feature_idxs, arg_enum$cell_idxs, arg_enum$umi_counts, arg_enum$mt_gene_bool)
    acc_constructor <- integer
    acc_length <- arg_enum$n_cells
  }
  return(list(acc_funct = acc_funct, acc_args = acc_args, acc_constructor = acc_constructor, acc_length = acc_length))
}
