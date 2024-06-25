#' Mixin Utility class
#'
#' @description
#' class that holds a few util methods used across all classes in the this R package
MixinUtilities <- R6::R6Class("MixinUtilities", public = list(

  # // @formatter:off
  #' @export
  # // @formatter:on
  create_empty_dataframe = function(col_names) {
    df <- data.frame(matrix(nrow = 0, ncol = length(col_names)))
    colnames(df) <- col_names
    return(df)
  }
))