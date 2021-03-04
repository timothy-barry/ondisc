#' setNames
#'
#' My personal setNames function
#'
#' @param object the object to name
#' @param nm the names
#'
#' @return the named object
#' @noRd
setNames <- function(object, nm) {
  names(object) <- nm
  object
}


#' Combine multimodal dataframes
#'
#' @param df_list list of dfs
#' @param modality_names names of the modalities
#'
#' @return a combined dataframe containing all the cell-specific covariates
#' @noRd
combine_multimodal_dataframes <- function(df_list, modality_names) {
  new_names <- lapply(seq(1, length(modality_names)), function(i) {
    paste0(modality_names[i], "_", colnames(df_list[[i]]))
  }) %>% unlist()
  out <- do.call(what = cbind, args = df_list)
  colnames(out) <- new_names
  return(out)
}
