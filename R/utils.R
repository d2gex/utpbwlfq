#' @title Mixin Utility class
#'
#' @description
#' It holds a few util methods used across all classes in the this R package
MixinUtilities <- R6::R6Class("MixinUtilities", public = list( # nolint

  # // @formatter:off
  #' @description
  #' Create an empty dataframe
  #'
  #' @param col_names vector with the column names of the empty dataframe
  #' @export
  # // @formatter:on
  create_empty_dataframe = function(col_names) {
    df <- data.frame(matrix(nrow = 0, ncol = length(col_names)))
    colnames(df) <- col_names
    return(df)
  }
))
