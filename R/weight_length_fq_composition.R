#' @title WLFeqComposition Class
#'
#' @description
#' Main class to build frequency matrices out of weight and lengths
WLFeqComposition <- R6::R6Class("WLFeqComposition", # nolint
  inherit = MixinUtilities,
  public = list(
    #' @field data a long dataframe with the length and weight data
    data = NULL,
    #' @field size_col name of the column that holds lengths
    size_col = NULL,
    #' @field weight_col name of the column that holds weights
    weight_col = NULL,
    #' @field time_col name of the column that holds dates
    time_col = NULL,
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
    #' Initialise the WLFeqComposition
    #'
    #' @param data a long dataframe with the length and weight data
    #' @param size_col name of the column that holds lengths
    #' @param weight_col name of the column that holds weights
    #' @param time_col name of the column that holds dates
    #' @param interval_col name of the column that will hold the size interval
    #' (as opposed to midpoints)
    #' @param midpoint_col name of the column that will hold the length classes
    #' @param freq_col name of the column that will hold the the length totals per class
    #' @param linf asymptotic length
    #' @export
    # // @formatter:on
    initialize = function(data, size_col, weight_col, time_col,
                          interval_col, midpoint_col, freq_col, linf) {
      self$data <- data
      self$size_col <- size_col
      self$weight_col <- weight_col
      self$time_col <- time_col
      self$interval_col <- interval_col
      self$midpoint_col <- midpoint_col
      self$freq_col <- freq_col
      self$linf <- linf
    },
    # // @formatter:off
    #' @description
    #' Generate the number of catches at length dataframe, starting from the minimum
    #' length class - min_padding and by capping the maximum length class to linf, if
    #' up_to_inf is set to TRUE
    #'
    #' @param bindwidth size of each length class
    #' @param min_padding padding that is substracted to the minimum length class
    #' @param up_to_linf boolean flag indicating whether the length classed should be capped
    #' @returns length frequency dataframe with dates as columns and length classes as rows
    # // @formatter:on
    generate_catch_at_length = function(bindwidth,
                                        min_padding,
                                        up_to_linf = TRUE) {
      interval_midpoints <- private$generate_interval_and_midpoint_sequences(
        bindwidth,
        min_padding,
        up_to_linf
      )
      catch_at_length <-
        private$generate_catch_at_length_df(
          interval_midpoints$size_intervals,
          interval_midpoints$mid_points
        )
      return(catch_at_length)
    },
    # // @formatter:off
    #' @description
    #' Generate the number of catches and mean weight at length dataframe,
    #' starting from the minimum length class - min_padding and by capping the maximum
    #' length class to linf, if up_to_inf is set to TRUE
    #'
    #' @param bindwidth size of each length class
    #' @param weight_na_as value to replace NA for weights
    #' @param min_padding padding that is substracted to the minimum length class
    #' @param up_to_linf boolean flag indicating whether the length classed should be capped
    #' @returns dataframe that keeps the LHT values per taxon
    # // @formatter:on
    generate_catch_and_mean_weight_at_length = function(bindwidth,
                                                        weight_na_as = 0,
                                                        min_padding,
                                                        up_to_linf = TRUE) {
      interval_midpoints <- private$generate_interval_and_midpoint_sequences(
        bindwidth,
        min_padding,
        up_to_linf
      )
      catch_at_length <-
        private$generate_catch_at_length_df(
          interval_midpoints$size_intervals,
          interval_midpoints$mid_points
        )
      size_intervas_weight_df <-
        private$generate_at_length_dataframe(interval_midpoints$size_intervals)
      mean_weight_at_length <-
        private$generate_mean_weight_at_length_df(size_intervas_weight_df)
      catch_mean_weight_length <- merge(catch_at_length,
        mean_weight_at_length,
        by = c(self$time_col, self$interval_col),
        all = TRUE
      )
      catch_mean_weight_length <- catch_mean_weight_length %>%
        tidyr::replace_na(list(mean_weight = weight_na_as))
      return(catch_mean_weight_length)
    }
  ),
  private = list(
    # // @formatter:off
    #' Generate interval and midpoint sequences for a given bidnwidth
    #'
    #' @returns a named list with the generated size intervals and midpoints
    # // @formatter:on
    generate_interval_and_midpoint_sequences = function(bindwidth,
                                                        min_padding = 0,
                                                        up_to_linf = TRUE) {
      if (up_to_linf & is.null(self$linf)) {
        stop("You need to provide a finite value to 'linf' if the flag 'up_to_linf' is set to TRUE.")
      }
      min_size <- floor(min(self$data[, self$size_col])) - min_padding
      max_size <- ceiling(max(self$data[, self$size_col]))
      if (up_to_linf) {
        max_size <- ceiling(max(max_size, self$linf))
      }
      half_bindwidth <- bindwidth / 2
      unique_size_intervals <- seq(min_size, max_size, bindwidth)
      mid_points <- seq(min_size + half_bindwidth, max_size - half_bindwidth, bindwidth)
      lengt_test <- length(unique_size_intervals) == length(mid_points) + 1
      testit::assert(deparse(lengt_test), lengt_test)
      return(list(size_intervals = unique_size_intervals, mid_points = mid_points))
    },
    # // @formatter:off
    #' Generate a dataframe organised by size intervals of width 'bindwidth' and
    #' its midpoints for each time period in the original dataframe. 'size_col' and
    #' 'weight_col' are as well added as columns.
    #'
    #' @returns a primitive dataframe with the initial frequency structure
    # // @formatter:on
    generate_at_length_dataframe = function(unique_size_intervals) {
      unique_time_periods <- unique(self$data[[self$time_col]])
      columns <- c(
        self$time_col, self$interval_col, self$midpoint_col,
        self$size_col, self$weight_col
      )
      size_weight_time_df <- self$create_empty_dataframe(columns)
      for (time_period in unique_time_periods) {
        yearly_data <- self$data %>%
          dplyr::filter_at(.vars = self$time_col, ~ .x == time_period) %>%
          dplyr::select_at(.vars = c(self$time_col, self$size_col, self$weight_col))
        intervals <- cut(yearly_data[[self$size_col]], unique_size_intervals)
        yearly_data <- yearly_data %>%
          dplyr::mutate(!!columns[2] := intervals)
        size_weight_time_df <- rbind(size_weight_time_df, yearly_data)
      }
      return(size_weight_time_df)
    },
    # // @formatter:off
    #' Generate the length frequency dataframe that contains the number of catches per
    #' length class
    #'
    #' @returns A dataframe with the number of catches at length
    # // @formatter:on
    generate_catch_at_length_df = function(unique_size_intervals, mid_points) {
      unique_time_periods <- unique(self$data[[self$time_col]])
      columns <- c(self$time_col, self$interval_col, self$midpoint_col, self$freq_col)
      catch_length_df <- self$create_empty_dataframe(columns)
      for (time_period in unique_time_periods) {
        yearly_data <- self$data %>%
          dplyr::filter_at(.vars = self$time_col, ~ .x == time_period)

        #  --> Build frequency table
        yearly_intervals <- as.data.frame(
          table(
            cut(yearly_data[[self$size_col]], unique_size_intervals),
            dnn = columns[2]
          ),
          responseName = columns[4]
        )
        # --> Add midpoints and year columns
        yearly_intervals <- yearly_intervals %>%
          dplyr::mutate(
            !!columns[1] := time_period,
            !!columns[3] := mid_points
          ) %>%
          dplyr::select_at(.vars = columns)
        # --> concat year intervals together
        catch_length_df <- rbind(catch_length_df, yearly_intervals)
      }
      return(catch_length_df)
    },
    # // @formatter:off
    #' Generate the weight frequency dataframe that contains the mean weight of catches per
    #' length class
    #'
    #' size_interval_weight_df primitive dataframe with the initial length frequency structure
    #' @returns A dataframe with the mean weight of catches at length
    # // @formatter:on
    generate_mean_weight_at_length_df = function(size_interval_weight_df) {
      average_weight <- size_interval_weight_df %>%
        dplyr::group_by_at(.vars = c(self$time_col, self$interval_col)) %>%
        dplyr::summarise(mean_weight = mean(.data[[self$weight_col]], na.rm = TRUE))
      return(as.data.frame(average_weight))
    }
  )
)
