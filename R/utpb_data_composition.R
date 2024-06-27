#' @title UtpbDataComposition  Class
#'
#' @description
#' Custom class to build the catches and mean-weight matrices for the UTPB database
UtpbDataComposition <- R6::R6Class("UtpbDataComposition", public = list( # nolint
  #' @field data a long dataframe with the length and weight data
  data = NULL,
  #' @field time_col name of the column that holds dates
  time_col = NULL,
  #' @field size_col name of the column that holds lengths
  size_col = NULL,
  #' @field weight_col name of the column that holds weights
  weight_col = NULL,
  #' @field mean_weight_col name of the column that will hold the mean weights
  mean_weight_col = NULL,
  #' @field interval_col name of the column that will hold the size interval (as opposed to midpoints)
  interval_col = NULL,
  #' @field midpoint_col name of the column that will hold the length classes
  midpoint_col = NULL,
  #' @field freq_col name of the column that will hold the the length totals per class
  freq_col = NULL,
  #' @field linf asymptotic length
  linf = NULL,

  # // @formatter:off
  #' @description
  #' Initialise the SpeciesDataComposition
  #'
  #' @param data a long dataframe with the length and weight data about a single species
  #' @param sp_arte_context a DataCompositionContext class
  #' @export
  # // @formatter:on
  initialize = function(data, sp_arte_context) {
    self$data <- data
    self$time_col <- sp_arte_context$time_col
    self$size_col <- sp_arte_context$size_col
    self$weight_col <- sp_arte_context$weight_col
    self$mean_weight_col <- sp_arte_context$mean_weight_col
    self$interval_col <- sp_arte_context$interval_col
    self$midpoint_col <- sp_arte_context$midpoint_col
    self$freq_col <- sp_arte_context$freq_col
    self$linf <- sp_arte_context$linf
  },
  # // @formatter:off
  #' @description
  #' Build a year-basis catch-at-length composition matrix when the weight column is
  #' not debugged, therefore it could be NA
  #'
  #' @param bindwidth size of each length class
  #' @param col_prefix prefix added to the each column representing a year
  #' @param min_padding padding that is substracted to the minimum length class
  #' @param up_to_linf boolean flag indicating whether the length classed should be capped
  #' @returns a composition result object with both long and wide catch-at-length
  #' @export
  # // @formatter:on
  build_catch_at_length = function(bindwidth, col_prefix, min_padding, up_to_linf = TRUE) {
    composition_results <- DataCompositionContainer$new()
    catch_compositor <- WLFeqComposition$new(
      self$data,
      self$size_col,
      self$weight_col,
      self$time_col,
      self$interval_col,
      self$midpoint_col,
      self$freq_col,
      self$linf
    )
    summary_catch <- catch_compositor$
      generate_catch_at_length(bindwidth, min_padding, up_to_linf)
    summary_catch_long_wide <-
      private$build_long_wide_variable_composition_matrix(
        summary_catch,
        col_prefix,
        self$freq_col
      )
    composition_results$catch_long_t <- summary_catch_long_wide$long
    composition_results$catch_wide_t <- summary_catch_long_wide$wide
    return(composition_results)
  },

  # // @formatter:off
  #' @description
  #' Build a year-basis catch-at-length and mean-weight composition matrices both in long and wide format
  #' (4 matrices, 2 for catch_at_length and 2 for mean-weight)
  #'
  #' @param bindwidth size of each length class
  #' @param col_prefix prefix added to the each column representing a year
  #' @param min_padding padding that is substracted to the minimum length class
  #' @param weight_na_as value to replace NA weights
  #' @param up_to_linf boolean flag indicating whether the length classed should be capped
  #' @returns a composition result object with both long and wide dataframes for both catch and mean-weight at length
  #' @export
  # // @formatter:on
  build_catch_and_mean_weight_at_length = function(bindwidth,
                                                   col_prefix,
                                                   min_padding = 0,
                                                   weight_na_as = 0,
                                                   up_to_linf = TRUE) {
    data <- self$data %>%
      dplyr::filter_at(.vars = self$weight_col, not_na)

    # Catch at length and mean weigh composition generator
    composition_results <- DataCompositionContainer$new()
    catch_weight_compositor <- WLFeqComposition$new(
      data,
      self$size_col,
      self$weight_col,
      self$time_col,
      self$interval_col,
      self$midpoint_col,
      self$freq_col,
      self$linf
    )
    # (1) Generate catch and mean-weight
    catch_mean_weight_at_length <-
      catch_weight_compositor$generate_catch_and_mean_weight_at_length(
        bindwidth,
        weight_na_as,
        min_padding,
        up_to_linf
      )

    # (2) Build long and wide catch-at-length dataframe
    catch_details <- catch_mean_weight_at_length %>% dplyr::select(-!!self$mean_weight_col)
    summary_catch <- private$build_long_wide_variable_composition_matrix(catch_details, col_prefix, self$freq_col)
    composition_results$catch_long_wt <- summary_catch$long
    composition_results$catch_wide_wt <- summary_catch$wide

    # (3) Build long and wide mean-weight-at-length dataframe
    mean_weight_details <- catch_mean_weight_at_length %>% dplyr::select(-!!self$freq_col)
    to_kg <- function(x) {
      round(x / 100, 2)
    }
    summary_mean_weight <- private$build_long_wide_variable_composition_matrix(mean_weight_details,
      col_prefix,
      self$mean_weight_col,
      converter_func = to_kg
    )
    composition_results$mwl <- summary_mean_weight$long
    composition_results$mww <- summary_mean_weight$wide
    return(composition_results)
  }
), private = list(
  build_long_wide_variable_composition_matrix = function(data,
                                                         col_prefix,
                                                         variable,
                                                         converter_func = NULL) {
    # return long frame
    summary_long <- data %>%
      dplyr::arrange_at(.vars = c(self$time_col, self$interval_col))
    if (!is.null(converter_func)) {
      summary_long <- summary_long %>%
        dplyr::mutate_at(.vars = variable, .funs = converter_func)
    }

    # return wide frame
    summary_wide <- summary_long %>%
      tidyr::pivot_wider(id_cols = -!!self$interval_col, names_from = !!self$time_col, values_from = !!variable)
    summary_wide <- summary_wide %>%
      dplyr::mutate_at(dplyr::vars(-self$midpoint_col), tidyr::replace_na, 0) %>%
      dplyr::rename_with(~ stringr::str_c(col_prefix, .x), tidyr::matches("^\\d{4}$"))

    summary_wide %>% assertr::assert(assertr::not_na, colnames(.), success_fun = assertr::success_logical)
    return(list(long = summary_long, wide = summary_wide))
  }
))
