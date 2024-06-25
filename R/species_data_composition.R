#' Review this script as now the species-arte pair should most likely be done through the arte column instead of
#' arte_grupo
library("R6")
library("dplyr")
library("tidyr")
library("stringr")
source("utils.R")
source("species_lwrel_lc_wc/catch_weight_composition.R")
source("species_lwrel_lc_wc/class_utils.R")

SpeciesDataComposition <- R6Class("SpeciesDataComposition", public = list(
  species_data = NULL,
  species = NULL,
  gears = NULL,
  time_col = NULL,
  size_col = NULL,
  weight_col = NULL,
  mean_weight_col = NULL,
  interval_col = NULL,
  midpoint_col = NULL,
  catch_col = NULL,
  data_holder = NULL,
  linf = NULL,

  initialize = function(db_data, sp_arte_context) {
    self$species <- sp_arte_context$species
    self$gears <- sp_arte_context$gears
    self$time_col <- sp_arte_context$time_col
    self$size_col <- sp_arte_context$size_col
    self$weight_col <- sp_arte_context$weight_col
    self$mean_weight_col <- sp_arte_context$mean_weight_col
    self$interval_col <- sp_arte_context$interval_col
    self$midpoint_col <- sp_arte_context$midpoint_col
    self$catch_col <- sp_arte_context$catch_col
    self$linf <- sp_arte_context$linf
    self$data_holder <- LCLWHolder$new()

    # TODO: Take this out the class. This is a bad example of unnecessary coupling
    self$species_data <- db_data %>%
      filter(ESPECIE == self$species) %>%
      filter(id_arte_grupo %in% self$gears) %>%
      select_at(.vars = c(self$time_col, self$size_col, self$weight_col))

  },
  build_talla_only_composition_matrix = function(bindwidth, col_prefix, min_padding, up_to_linf = TRUE) {
    # // @formatter:off
    #' Build a year-basis catch-at-length composition composition matrix when the weight column is
    #' not debugged, therefore it could be NA
    # // @formatter:on

    # Catch at length and mean weigh composition generator
    composition <- CatchWeightComposition$new(
      self$species_data,
      self$size_col,
      self$weight_col,
      self$time_col,
      self$interval_col,
      self$midpoint_col,
      self$catch_col,
      self$linf
    )
    summary_catch <- composition$
      generate_catch_at_length_composition(bindwidth, min_padding, up_to_linf)
    summary_catch_long_wide <-
      private$build_long_wide_variable_composition_matrix(summary_catch,
                                                          col_prefix,
                                                          self$catch_col)
    self$data_holder$catch_long_t <- summary_catch_long_wide$long
    self$data_holder$catch_wide_t <- summary_catch_long_wide$wide
  },

  build_talla_and_weight_composition_matrices = function(bindwidth, col_prefix, min_padding, up_to_linf = TRUE) {
    # // @formatter:off
    #' Build a year-basis catch-at-length and mean-weight composition matrices both in long and wide format
    #' (4 matrices, 2 for catch_at_length and 2 for mean-weight)
    # // @formatter:on

    species_data <- self$species_data %>%
      filter_at(.vars = self$weight_col, not_na)
    # Catch at length and mean weigh composition generator
    composition <- CatchWeightComposition$new(
      species_data,
      self$size_col,
      self$weight_col,
      self$time_col,
      self$interval_col,
      self$midpoint_col,
      self$catch_col,
      self$linf
    )
    # (1) Generate catch and mean-weight
    catch_mean_weight_at_length <-
      composition$generate_catch_and_m.weight_at_length_composition(bindwidth, min_padding = 1, up_to_linf)

    # (2) Build long and wide catch-at-length dataframe
    catch_details <- catch_mean_weight_at_length %>% select(-!!self$mean_weight_col)
    summary_catch <- private$build_long_wide_variable_composition_matrix(catch_details, col_prefix, self$catch_col)
    self$data_holder$catch_long_wt <- summary_catch$long
    self$data_holder$catch_wide_wt <- summary_catch$wide

    # (3) Build long and wide mean-weight-at-length dataframe
    mean_weight_details <- catch_mean_weight_at_length %>% select(-!!self$catch_col)
    to_kg <- function(x) { round(x / 100, 2) }
    summary_mean_weight <- private$build_long_wide_variable_composition_matrix(mean_weight_details,
                                                                               col_prefix,
                                                                               self$mean_weight_col,
                                                                               converter_func = to_kg)
    self$data_holder$mwl <- summary_mean_weight$long
    self$data_holder$mww <- summary_mean_weight$wide
  }

), private = list(
  build_long_wide_variable_composition_matrix = function(data,
                                                         col_prefix,
                                                         variable,
                                                         converter_func = NULL) {

    # return long frame
    summary_long <- data %>%
      arrange_at(.vars = c(self$time_col, self$interval_col))
    if (!is.null(converter_func)) {
      summary_long <- summary_long %>%
        mutate_at(.vars = variable, .funs = converter_func)
    }

    # return wide frame
    summary_wide <- summary_long %>%
      pivot_wider(id_cols = -!!self$interval_col, names_from = !!self$time_col, values_from = !!variable)
    summary_wide <- summary_wide %>%
      mutate_at(vars(-self$midpoint_col), replace_na, 0) %>%
      rename_with(~str_c(col_prefix, .x), matches("^\\d{4}$"))

    mute <- summary_wide %>% assertr::assert(assertr::not_na, colnames(.))
    return(list(long = summary_long, wide = summary_wide))
  }

))


