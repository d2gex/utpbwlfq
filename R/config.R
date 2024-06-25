#' @title DataCompositionContainer Class
#'
#' @description
#' A class container that holds the long and wide catch-at-length and mean-weight-at-length dataframes
#' @export
DataCompositionContainer <- R6Class("DataCompositionContainer", public = list( # nolint
  # (1)  Length Composition (lengths), ignoring the values of weight column
  catch_long_t = NULL,
  catch_wide_t = NULL,
  # (2)  Length Composition (lengths); only those rows that also have a value for weight are considered
  catch_long_wt = NULL,
  catch_wide_wt = NULL,
  # (3) Weight composition (weight)
  mwl = NULL,
  mww = NULL,
  initialize = function() {
  }
))
#' @title SpeciesGearDataContext Class
#'
#' @description
#' A class that provide context about
#' @export
DataCompositionContext <- R6Class("DataCompositionContext", public = list( # nolint
  bindwidth = 1, # cm
  time_col = "year",
  size_col = "TALLA",
  weight_col = "PESO",
  mean_weight_col = "mean_weight",
  interval_col = "interval",
  midpoint_col = "MeanLength",
  catch_col = "catch",
  col_prefix = "X",
  linf = NULL,
  initialize = function() {
  }
))
