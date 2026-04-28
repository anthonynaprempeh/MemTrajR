# utils.R — Internal helper functions (not exported)

#' Generate Log10 Minor Tick Positions
#'
#' Internal helper used by plotting functions to place minor ticks on
#' logarithmic axes.
#'
#' @param major_ticks Numeric vector of major tick positions (powers of 10).
#' @return Numeric vector of minor tick positions (2x through 9x each decade).
#' @keywords internal
.log_minor_ticks <- function(major_ticks) {
  minor <- c()
  for (i in seq_len(length(major_ticks) - 1)) {
    decade <- major_ticks[i]
    minor  <- c(minor, decade * 2:9)
  }
  minor
}
