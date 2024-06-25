library("R6")
library("tidyr")
source("utils.R")

CatchWeightComposition <- R6Class("CatchWeightComposition",
                                  public = list(
                                    data = NULL,
                                    size_col = NULL,
                                    weight_col = NULL,
                                    time_col = NULL,
                                    interval_col = NULL,
                                    midpoint_col = NULL,
                                    freq_col = NULL,
                                    linf = NULL,
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
                                    generate_catch_at_length_composition = function(bindwidth,
                                                                                    min_padding,
                                                                                    up_to_linf = TRUE) {
                                      interval_midpoints <- private$generate_interval_and_midpoint_sequences(bindwidth,
                                                                                                             min_padding,
                                                                                                             up_to_linf)
                                      catch_at_length <-
                                        private$generate_catch_at_length_frequency(interval_midpoints$size_intervals,
                                                                                   interval_midpoints$mid_points)
                                      return(catch_at_length)
                                    },
                                    generate_catch_and_m.weight_at_length_composition = function(bindwidth,
                                                                                                 weight_na_as = 0,
                                                                                                 min_padding,
                                                                                                 up_to_linf = TRUE) {
                                      interval_midpoints <- private$generate_interval_and_midpoint_sequences(bindwidth,
                                                                                                             min_padding,
                                                                                                             up_to_linf)
                                      catch_at_length <-
                                        private$generate_catch_at_length_frequency(interval_midpoints$size_intervals,
                                                                                   interval_midpoints$mid_points)
                                      size_intervas_weight_df <-
                                        private$generate_length_intervals(interval_midpoints$size_intervals)
                                      mean.weight_at_length <-
                                        private$generate_mean_weight_at_length(size_intervas_weight_df)
                                      catch_m.weight_length <- merge(catch_at_length,
                                                                     mean.weight_at_length,
                                                                     by = c(self$time_col, self$interval_col),
                                                                     all = TRUE)
                                      catch_m.weight_length <- catch_m.weight_length %>%
                                        replace_na(list(mean_weight = weight_na_as))
                                      return(catch_m.weight_length)
                                    }
                                  ),
                                  private = list(
                                    generate_interval_and_midpoint_sequences = function(bindwidth,
                                                                                        min_padding = 0,
                                                                                        up_to_linf = TRUE) {
                                      # // @formatter:off
                                      #' Generate interval and midpoint sequences for a given bidnwidth
                                      # // @formatter:on
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

                                    generate_length_intervals = function(unique_size_intervals) {
                                      # // @formatter:off
                                      #' Generate a dataframe organised by size intervals of width 'bindwidth' and
                                      #' its midpoints for each time period in the original dataframe. 'size_col' and
                                      #' 'weight_col' are as well added as columns.
                                      # // @formatter:off

                                      unique_time_periods <- unique(self$data[[self$time_col]])
                                      columns <- c(self$time_col, self$interval_col, self$midpoint_col,
                                                   self$size_col, self$weight_col)
                                      size_weight_time_df <- create_empty_dataframe(columns)
                                      for (time_period in unique_time_periods) {
                                        yearly_data <- self$data %>%
                                          filter_at(.vars = self$time_col, ~.x == time_period) %>%
                                          select_at(.vars = c(self$time_col, self$size_col, self$weight_col))
                                        intervals <- cut(yearly_data[[self$size_col]], unique_size_intervals)
                                        yearly_data <- yearly_data %>%
                                          mutate(
                                            !!columns[2] := intervals
                                          )
                                        size_weight_time_df <- rbind(size_weight_time_df, yearly_data)
                                      }
                                      return(size_weight_time_df)
                                    },
                                    generate_catch_at_length_frequency = function (unique_size_intervals, mid_points) {
                                      # // @formatter:off
                                      #' Generate catch at length frequency table (contingency table)
                                      # // @formatter:on

                                      unique_time_periods <- unique(self$data[[self$time_col]])
                                      columns <- c(self$time_col, self$interval_col, self$midpoint_col, self$freq_col)
                                      catch_length_df <- create_empty_dataframe(columns)
                                      for (time_period in unique_time_periods) {
                                        yearly_data <- self$data %>%
                                          filter_at(.vars = self$time_col, ~.x == time_period)

                                        #  --> Build frequency table
                                        yearly_intervals <- as.data.frame(
                                          table(
                                            cut(yearly_data[[self$size_col]], unique_size_intervals),
                                            dnn = columns[2]),
                                          responseName = columns[4])
                                        # --> Add midpoints and year columns
                                        yearly_intervals <- yearly_intervals %>%
                                          mutate(
                                            !!columns[1] := time_period,
                                            !!columns[3] := mid_points
                                          ) %>%
                                          select_at(.vars = columns)
                                        # --> concat year intervals together
                                        catch_length_df <- rbind(catch_length_df, yearly_intervals)
                                      }
                                      return(catch_length_df)
                                    },

                                    generate_mean_weight_at_length = function(size_interval_weight_df) {
                                      average_weight <- size_interval_weight_df %>%
                                        group_by_at(.vars = c(self$time_col, self$interval_col)) %>%
                                        summarise(mean_weight = mean(.data[[self$weight_col]], na.rm = TRUE))
                                      return(as.data.frame(average_weight))
                                    }
                                  ))