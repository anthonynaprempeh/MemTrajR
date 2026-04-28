# separation.R — Geodesic separation, interaction potential, MSSep, and plots
# Derived from: 3d_full_analysis.R


# ============================================================
#  Geodesic separation
# ============================================================

#' Compute Frame-by-Frame Geodesic Separation Between Two Domains
#'
#' Given 3D coordinates of two domains for each frame, computes the
#' geodesic (great-circle arc) distance between them on a sphere of
#' radius `R`.
#'
#' @param x1,y1,z1 Numeric vectors. 3D coordinates of domain 1 (micrometres).
#' @param x2,y2,z2 Numeric vectors. 3D coordinates of domain 2 (micrometres).
#' @param R Positive numeric. Vesicle radius in micrometres.
#'
#' @return Numeric vector of geodesic separations (micrometres), one value
#'   per frame.
#'
#' @export
#' @examples
#' R  <- 8
#' n  <- 50
#' x1 <- rnorm(n, sd = 2); y1 <- rnorm(n, sd = 2); z1 <- sqrt(R^2 - x1^2 - y1^2)
#' x2 <- rnorm(n, sd = 2); y2 <- rnorm(n, sd = 2); z2 <- sqrt(R^2 - x2^2 - y2^2)
#' d  <- compute_geodesic_separation(x1, y1, z1, x2, y2, z2, R = R)
#' summary(d)
compute_geodesic_separation <- function(x1, y1, z1, x2, y2, z2, R) {
  dot12 <- x1 * x2 + y1 * y2 + z1 * z2
  cosg  <- pmin(1, pmax(-1, dot12 / (R^2)))
  R * acos(cosg)
}


# ============================================================
#  Separation histogram and bin assignments
# ============================================================

#' Bin Geodesic Separations and Return Frame-Level Data
#'
#' Assigns each frame's geodesic separation to a histogram bin and
#' returns a data frame suitable for downstream analyses.
#'
#' @param d_geo Numeric vector of geodesic separations (micrometres),
#'   as returned by [compute_geodesic_separation()].
#' @param n_bins Positive integer. Number of equal-width bins spanning
#'   the range of `d_geo`. Default `20`.
#' @param frame_interval Positive numeric. Time between frames (seconds).
#'
#' @return A list with:
#'   \describe{
#'     \item{frame_bins}{Data frame with columns `frame`, `time_s`,
#'       `separation_um`, and `bin`.}
#'     \item{bin_counts}{Data frame with columns `bin` and
#'       `count_frames`.}
#'     \item{edges}{Numeric vector of bin-edge values.}
#'     \item{bin_midpoints}{Numeric vector of bin midpoint values.}
#'   }
#'
#' @export
#' @examples
#' d_geo <- abs(rnorm(200, mean = 5, sd = 1))
#' result <- bin_separations(d_geo, n_bins = 20, frame_interval = 0.1)
#' head(result$frame_bins)
#' result$bin_counts
bin_separations <- function(d_geo, n_bins = 20, frame_interval) {
  n_frames    <- length(d_geo)
  frame_index <- seq_len(n_frames)
  time_s      <- (frame_index - 1) * frame_interval

  start <- min(d_geo, na.rm = TRUE)
  end   <- max(d_geo, na.rm = TRUE)
  edges <- seq(start, end, length.out = n_bins + 1)
  bin_midpoints <- (utils::head(edges, -1) + utils::tail(edges, -1)) / 2

  bin_label  <- cut(d_geo, breaks = edges, include.lowest = TRUE,
                    right = FALSE)
  frame_bins <- data.frame(
    frame         = frame_index,
    time_s        = time_s,
    separation_um = d_geo,
    bin           = as.character(bin_label),
    stringsAsFactors = FALSE
  )

  counts     <- as.data.frame(table(bin_label))
  bin_counts <- data.frame(
    bin          = as.character(counts$bin_label),
    count_frames = counts$Freq,
    stringsAsFactors = FALSE
  )

  list(
    frame_bins    = frame_bins,
    bin_counts    = bin_counts,
    edges         = edges,
    bin_midpoints = bin_midpoints
  )
}


# ============================================================
#  Interaction potential
# ============================================================

#' Estimate Effective Interaction Potential via Boltzmann Analysis
#'
#' Converts a histogram of geodesic separations into an effective
#' interaction potential \eqn{U(r) / k_B T = -\ln[P(r) / P(r_{\max})]},
#' where \eqn{P(r_{\max})} (the last bin) is used as the reference state.
#'
#' @param bin_counts Data frame with columns `bin` and `count_frames`,
#'   as returned by [bin_separations()]\code{$bin_counts}.
#' @param bin_midpoints Numeric vector of bin midpoint values (micrometres),
#'   as returned by [bin_separations()]\code{$bin_midpoints}.
#'
#' @return A data frame with columns `bin`, `bin_midpoint_um`,
#'   `count_frames`, `P_r`, and `U_over_kBT`.
#'
#' @export
#' @examples
#' d_geo  <- abs(rnorm(200, mean = 5, sd = 1))
#' bins   <- bin_separations(d_geo, n_bins = 15, frame_interval = 0.1)
#' potential <- estimate_interaction_potential(
#'   bins$bin_counts, bins$bin_midpoints)
#' head(potential)
estimate_interaction_potential <- function(bin_counts, bin_midpoints) {
  counts         <- bin_counts$count_frames
  last_bin_count <- utils::tail(counts, 1)

  if (last_bin_count == 0)
    stop(paste(
      "The last bin has zero counts;",
      "cannot normalise P(r) using the last bin as reference."
    ))

  P_r        <- counts / last_bin_count
  U_over_kBT <- ifelse(P_r > 0, -log(P_r), NA_real_)

  data.frame(
    bin             = bin_counts$bin,
    bin_midpoint_um = bin_midpoints,
    count_frames    = counts,
    P_r             = P_r,
    U_over_kBT      = U_over_kBT,
    stringsAsFactors = FALSE
  )
}


# ============================================================
#  Mean Square Separation Increment (MSSep)
# ============================================================

#' Compute Mean Square Separation Increment (MSSep)
#'
#' Computes \eqn{\langle [\Delta r(\tau)]^2 \rangle} as a function of
#' lag time \eqn{\tau} for a geodesic separation time series.
#'
#' @param d_geo Numeric vector of geodesic separations (micrometres).
#' @param frame_interval Positive numeric. Time between frames (seconds).
#'
#' @return A data frame with columns `dt` (lag time in seconds) and
#'   `ms_inc_um2` (mean square separation increment in micrometres squared).
#'
#' @export
#' @examples
#' d_geo <- abs(rnorm(200, mean = 5, sd = 0.5))
#' ms    <- compute_mssep(d_geo, frame_interval = 0.1)
#' head(ms)
compute_mssep <- function(d_geo, frame_interval) {
  sep  <- as.numeric(d_geo)
  N    <- length(sep)
  lags <- seq_len(N - 1)
  dts  <- lags * frame_interval

  ms_inc <- sapply(lags, function(lag) {
    ds <- sep[(lag + 1):N] - sep[1:(N - lag)]
    mean(ds^2, na.rm = TRUE)
  })

  data.frame(dt = dts, ms_inc_um2 = ms_inc)
}


#' Fit MSSep Curve (Linear and Power-Law)
#'
#' Fits the MSSep vs lag-time data with a linear model and a log-log
#' power-law model over an early-time window, and estimates a separation
#' diffusion coefficient \eqn{D_{sep} = \text{slope} / 4}.
#'
#' @param mssep_data Data frame with columns `dt` and `ms_inc_um2`, as
#'   returned by [compute_mssep()].
#' @param fit_until Positive numeric. Upper time limit (seconds) for the
#'   fitting window. Default `1`.
#'
#' @return A named list with elements `fit_summary` (data frame),
#'   `D_sep` (numeric), `linear_slope`, `linear_intercept`,
#'   `power_prefactor`, `power_exponent`, and `log_intercept`.
#'
#' @export
#' @examples
#' d_geo <- abs(rnorm(300, mean = 5, sd = 0.5))
#' ms    <- compute_mssep(d_geo, frame_interval = 0.1)
#' fit   <- fit_mssep(ms, fit_until = 1)
#' fit$fit_summary
fit_mssep <- function(mssep_data, fit_until = 1) {
  result_pos <- mssep_data[mssep_data$dt > 0 & mssep_data$ms_inc_um2 > 0, ]
  fit_data   <- result_pos[result_pos$dt <= fit_until, ]

  if (nrow(fit_data) < 2)
    stop("Not enough points within the selected fitting range.")

  fit_linear <- stats::lm(ms_inc_um2 ~ dt, data = fit_data)
  c0         <- stats::coef(fit_linear)[1]
  m_linear   <- stats::coef(fit_linear)[2]
  D_sep      <- m_linear / 4

  fit_log    <- stats::lm(log10(ms_inc_um2) ~ log10(dt), data = fit_data)
  b          <- stats::coef(fit_log)[1]
  n_exp      <- stats::coef(fit_log)[2]
  A          <- 10^b

  fit_summary <- data.frame(
    parameter = c(
      "linear_equation", "power_law_equation",
      "linear_slope", "linear_intercept",
      "power_prefactor", "power_exponent",
      "diffusion_coefficient_um2_s"
    ),
    value = c(
      sprintf("y = %.4fx %+.4f", m_linear, c0),
      sprintf("y = %.4fx^%.5f",  A, n_exp),
      as.character(m_linear),
      as.character(c0),
      as.character(A),
      as.character(n_exp),
      as.character(D_sep)
    ),
    stringsAsFactors = FALSE
  )

  list(
    fit_summary      = fit_summary,
    D_sep            = D_sep,
    linear_slope     = m_linear,
    linear_intercept = c0,
    power_prefactor  = A,
    power_exponent   = n_exp,
    log_intercept    = b
  )
}


# ============================================================
#  Plotting helpers
# ============================================================

#' Plot Separation Histogram
#'
#' Saves a PNG histogram of geodesic separations.
#'
#' @param d_geo Numeric vector of geodesic separations (micrometres).
#' @param edges Numeric vector of bin edges, as returned by
#'   [bin_separations()]\code{$edges}.
#' @param sample_name Character. Used as plot title prefix and output
#'   file prefix.
#'
#' @return Invisibly `NULL`. A PNG file is written as a side effect.
#'
#' @export
#' @examples
#' \dontrun{
#' d_geo <- abs(rnorm(200, mean = 5, sd = 1))
#' bins  <- bin_separations(d_geo, n_bins = 20, frame_interval = 0.1)
#' plot_separation_histogram(d_geo, bins$edges, sample_name = "sample01")
#' }
plot_separation_histogram <- function(d_geo, edges, sample_name) {
  output_hist <- paste0(sample_name, "_separation_hist.png")
  grDevices::png(filename = output_hist, width = 1600, height = 1200, res = 220)
  graphics::par(mgp = c(2, 0.7, 0))
  graphics::hist(d_geo, breaks = edges, col = "red", border = "gray20",
    main = paste0(sample_name, " \u2014 Geodesic separation"),
    xlab = expression(bold("Geodesic separation (" * mu * "m)")),
    ylab = "Frequency", cex.lab = 1.2, font.lab = 2,
    cex.main = 1.2, font.main = 2)
  grDevices::dev.off()
  invisible(NULL)
}


#' Plot Effective Interaction Potential
#'
#' Saves a PNG of the interaction potential \eqn{U(r)/k_BT} vs. centre-to-
#' centre separation, with a reference point marking \eqn{(d_1 + d_2)/2}.
#'
#' @param potential_data Data frame as returned by
#'   [estimate_interaction_potential()].
#' @param min_separation Numeric. Contact distance \eqn{(d_1 + d_2)/2}
#'   in micrometres.
#' @param sample_name Character. Used as plot title prefix and output file
#'   prefix.
#'
#' @return Invisibly `NULL`. A PNG file is written as a side effect.
#'
#' @export
#' @examples
#' \dontrun{
#' d_geo <- abs(rnorm(200, mean = 5, sd = 1))
#' bins  <- bin_separations(d_geo, n_bins = 15, frame_interval = 0.1)
#' pot   <- estimate_interaction_potential(bins$bin_counts, bins$bin_midpoints)
#' plot_interaction_potential(pot, min_separation = 3.14, sample_name = "s01")
#' }
plot_interaction_potential <- function(potential_data, min_separation,
                                        sample_name) {
  output_file <- paste0(sample_name, "_potential_plot.png")
  x_min_plot  <- min(min_separation,
                     min(potential_data$bin_midpoint_um, na.rm = TRUE)) * 0.95
  x_max_plot  <- max(potential_data$bin_midpoint_um, na.rm = TRUE)

  grDevices::png(filename = output_file, width = 1800, height = 1400, res = 220)
  graphics::par(mgp = c(2, 0.7, 0))
  graphics::plot(
    potential_data$bin_midpoint_um, potential_data$U_over_kBT,
    type = "l", lwd = 3,
    xlim = c(x_min_plot, x_max_plot),
    xlab = expression("Centre to Centre Separations (" * mu * "m)"),
    ylab = expression("Interaction Potential (k"[B] * "T)"),
    main = paste0(sample_name, " \u2014 Effective interaction potential"),
    bty  = "l"
  )
  graphics::points(min_separation, 0, pch = 16, col = "red", cex = 1.2)
  graphics::text(min_separation, 0,
    labels = expression(frac(d[1] + d[2], 2)),
    col = "red", pos = 4, cex = 0.9)
  grDevices::dev.off()
  invisible(NULL)
}


#' Plot Separation vs Frame or vs Time
#'
#' Saves a PNG scatter plot of geodesic separation over time or frames,
#' with horizontal reference lines for the vesicle diameter and the contact
#' distance.
#'
#' @param x_vals Numeric vector. Frame indices or time values for the x-axis.
#' @param d_geo Numeric vector of geodesic separations (micrometres).
#' @param x_label Character. x-axis label.
#' @param output_file Character. Output PNG file path.
#' @param plot_title Character. Plot title.
#' @param Dv Positive numeric. Vesicle diameter in micrometres (reference
#'   maximum separation).
#' @param min_separation Positive numeric. Contact distance
#'   \eqn{(d_1 + d_2)/2} in micrometres.
#'
#' @return Invisibly `NULL`. A PNG file is written as a side effect.
#'
#' @export
#' @examples
#' \dontrun{
#' d_geo <- abs(rnorm(200, mean = 5, sd = 0.5))
#' plot_separation_vs_x(
#'   x_vals = seq_along(d_geo), d_geo = d_geo,
#'   x_label = "Frames", output_file = "sep_vs_frame.png",
#'   plot_title = "Sample 01 - Separation vs Frame",
#'   Dv = 16, min_separation = 3.14)
#' }
plot_separation_vs_x <- function(x_vals, d_geo, x_label, output_file,
                                  plot_title, Dv, min_separation) {
  grDevices::png(filename = output_file, width = 1800, height = 1400, res = 220)
  graphics::par(mgp = c(2.2, 0.8, 0), mar = c(5, 5, 2, 2))
  graphics::plot(x_vals, d_geo, type = "p", pch = 16, cex = 0.9,
    col  = "#2B6C8E",
    xlab = x_label,
    ylab = expression("Geodesic separation (" * mu * "m)"),
    main = plot_title, bty = "l", xaxs = "i", yaxs = "i",
    xlim = c(0, max(x_vals)),
    ylim = c(0, max(Dv, d_geo, na.rm = TRUE) + 1))
  graphics::abline(h = Dv,             col = "black", lwd = 2.5, lty = 1)
  graphics::abline(h = min_separation, col = "black", lwd = 2.0, lty = 2)
  graphics::legend("topright",
    legend = c("(d1+d2)/2", "Dv"),
    col = c("black", "black"), lwd = c(2.0, 2.5), lty = c(2, 1), bty = "n")
  grDevices::dev.off()
  invisible(NULL)
}


#' Plot r/R for Each Domain vs Frame or vs Time
#'
#' Saves a PNG scatter plot showing the radial coordinate ratio r/R for
#' two domains over time or frames.
#'
#' @param x_vals Numeric vector. Frame indices or time values for the x-axis.
#' @param r1_R,r2_R Numeric vectors. Radial coordinate ratios for domains
#'   1 and 2.
#' @param x_label Character. x-axis label.
#' @param output_file Character. Output PNG file path.
#' @param plot_title Character. Plot title.
#'
#' @return Invisibly `NULL`. A PNG file is written as a side effect.
#'
#' @export
#' @examples
#' \dontrun{
#' r1_R <- runif(200, 0.8, 1)
#' r2_R <- runif(200, 0.7, 1)
#' plot_rR_vs_x(
#'   x_vals = seq_along(r1_R), r1_R = r1_R, r2_R = r2_R,
#'   x_label = "Frames", output_file = "rR_vs_frame.png",
#'   plot_title = "Sample 01 - r/R vs Frame")
#' }
plot_rR_vs_x <- function(x_vals, r1_R, r2_R, x_label, output_file,
                          plot_title) {
  grDevices::png(filename = output_file, width = 1800, height = 1400, res = 220)
  graphics::par(mgp = c(2.2, 0.8, 0), mar = c(5, 5, 2, 2))
  graphics::plot(x_vals, r1_R, type = "p", pch = 16, cex = 0.6,
    col = "darkorange", xlab = x_label, ylab = "r / R",
    main = plot_title, bty = "l", xaxs = "i", yaxs = "i",
    xlim = c(0, max(x_vals)), ylim = c(0, 1))
  graphics::points(x_vals, r2_R, pch = 16, cex = 0.6, col = "darkgreen")
  graphics::legend("topright",
    legend = c("d1", "d2"), pch = 16,
    col = c("darkorange", "darkgreen"), bty = "n")
  grDevices::dev.off()
  invisible(NULL)
}


#' Plot MSSep on Log-Log Axes
#'
#' Saves a PNG log-log plot of the mean square separation increment vs
#' lag time.
#'
#' @param mssep_data Data frame with columns `dt` and `ms_inc_um2`, as
#'   returned by [compute_mssep()].
#' @param sample_name Character. Used as plot title prefix and output file
#'   prefix.
#' @param y_min Positive numeric. Lower y-axis limit. Default `0.01`.
#' @param y_max Positive numeric. Upper y-axis limit. Default `10`.
#'
#' @return Invisibly `NULL`. A PNG file is written as a side effect.
#'
#' @export
#' @examples
#' \dontrun{
#' d_geo <- abs(rnorm(300, mean = 5, sd = 0.5))
#' ms    <- compute_mssep(d_geo, frame_interval = 0.1)
#' plot_mssep(ms, sample_name = "sample01")
#' }
plot_mssep <- function(mssep_data, sample_name, y_min = 0.01, y_max = 10) {
  result_pos  <- mssep_data[mssep_data$dt > 0 & mssep_data$ms_inc_um2 > 0, ]
  x_min_data  <- min(result_pos$dt, na.rm = TRUE)
  x_max_data  <- max(result_pos$dt, na.rm = TRUE)
  x_min_mssep <- 10^floor(log10(x_min_data))
  x_max_mssep <- 10^ceiling(log10(x_max_data))
  y_max_auto  <- 10^ceiling(log10(max(result_pos$ms_inc_um2, na.rm = TRUE)))

  x_major <- 10^seq(log10(x_min_mssep), log10(x_max_mssep))
  y_major <- 10^seq(log10(y_min), log10(max(y_max, y_max_auto)))
  x_minor <- .log_minor_ticks(x_major)
  y_minor <- .log_minor_ticks(y_major)
  x_minor <- x_minor[x_minor > x_min_mssep & x_minor < x_max_mssep]
  y_minor <- y_minor[y_minor > y_min]

  output_file <- paste0(sample_name, "_mssep_plot.png")
  grDevices::png(filename = output_file, width = 1800, height = 1400, res = 220)
  graphics::par(mgp = c(2, 0.7, 0))
  graphics::plot(
    result_pos$dt, result_pos$ms_inc_um2,
    log = "xy", pch = 16, col = "black",
    xlim = c(x_min_mssep, x_max_mssep),
    ylim = c(y_min, max(y_max, y_max_auto)),
    xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", bty = "l",
    cex.lab = 1.5, font.lab = 2, cex.main = 1.6, font.main = 2,
    main = paste0(sample_name, " \u2014 MSSep"),
    xlab = "dt (s)",
    ylab = expression(bold(MSSep ~ "(" * mu * "m"^2 * ")"))
  )
  graphics::axis(1, at = x_major, labels = x_major, tcl = -0.5)
  graphics::axis(2, at = y_major, labels = y_major, tcl = -0.5)
  graphics::axis(1, at = x_minor, labels = FALSE, tcl = -0.25)
  graphics::axis(2, at = y_minor, labels = FALSE, tcl = -0.25)
  grDevices::dev.off()
  invisible(NULL)
}


# ============================================================
#  Main wrapper: analyze_separation
# ============================================================

#' Full Separation and Interaction Potential Analysis
#'
#' High-level wrapper that reads 3D domain coordinates, computes geodesic
#' separations, bins them, estimates the interaction potential, computes
#' MSSep, and writes all outputs (CSVs and PNGs) to the working directory.
#'
#' @param input_file Character. Path to a CSV file with columns `x1`, `y1`,
#'   `z1`, `x2`, `y2`, `z2` (3D coordinates per frame).
#' @param sample_name Character. Prefix for all output file names.
#' @param R Positive numeric. Vesicle radius in micrometres.
#' @param d1,d2 Positive numeric. Diameters of the two domains
#'   (micrometres). Used to compute the contact distance
#'   \eqn{(d_1 + d_2)/2}.
#' @param frame_interval Positive numeric. Time between frames (seconds).
#' @param fit_until Positive numeric. Upper time limit (seconds) for the
#'   MSSep fitting window. Default `1`.
#' @param n_bins Positive integer. Number of histogram bins. Default `20`.
#'
#' @return Invisibly, a named list with elements `d_geo`, `frame_bins`,
#'   `bin_counts`, `potential_data`, `mssep_data`, and `mssep_fit`.
#'   All output files are written as a side effect.
#'
#' @export
#' @examples
#' \dontrun{
#' results <- analyze_separation(
#'   input_file     = "0004_3d.csv",
#'   sample_name    = "0004_3d",
#'   R              = 8.02,
#'   d1             = 3.14,
#'   d2             = 3.14,
#'   frame_interval = 0.1,
#'   fit_until      = 1,
#'   n_bins         = 20
#' )
#' }
analyze_separation <- function(input_file, sample_name,
                                R, d1, d2,
                                frame_interval,
                                fit_until = 1,
                                n_bins    = 20) {

  dat <- utils::read.csv(input_file)
  dat <- dat[, !grepl("^X(\\.\\d+)?$", names(dat))]

  required_cols <- c("x1", "y1", "z1", "x2", "y2", "z2")
  missing_cols  <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0)
    stop(paste(
      "Missing required columns:", paste(missing_cols, collapse = ", "),
      "\nPlease name your CSV columns: x1, y1, z1, x2, y2, z2."
    ))

  x1 <- dat$x1; y1 <- dat$y1; z1 <- dat$z1
  x2 <- dat$x2; y2 <- dat$y2; z2 <- dat$z2

  n_frames    <- nrow(dat)
  frame_index <- seq_len(n_frames)
  time_s      <- (frame_index - 1) * frame_interval
  Dv          <- 2 * R
  min_sep     <- (d1 + d2) / 2

  # Geodesic separation
  d_geo  <- compute_geodesic_separation(x1, y1, z1, x2, y2, z2, R = R)
  bins   <- bin_separations(d_geo, n_bins = n_bins,
                             frame_interval = frame_interval)

  # Add raw coordinates to frame_bins
  bins$frame_bins <- cbind(
    bins$frame_bins[, c("frame", "time_s")],
    data.frame(x1 = x1, y1 = y1, z1 = z1, x2 = x2, y2 = y2, z2 = z2),
    bins$frame_bins[, c("separation_um", "bin")]
  )

  utils::write.csv(bins$frame_bins,
    paste0(sample_name, "_frame_bins.csv"), row.names = FALSE)
  utils::write.csv(bins$bin_counts,
    paste0(sample_name, "_bin_counts.csv"), row.names = FALSE)

  # Interaction potential
  potential_data <- estimate_interaction_potential(
    bins$bin_counts, bins$bin_midpoints)
  utils::write.csv(potential_data,
    paste0(sample_name, "_potential.csv"), row.names = FALSE)

  # Plots
  plot_separation_histogram(d_geo, bins$edges, sample_name)
  plot_interaction_potential(potential_data, min_sep, sample_name)
  plot_separation_vs_x(frame_index, d_geo, "Frames",
    paste0(sample_name, "_separation_vs_frame.png"),
    paste0(sample_name, " \u2014 Separation vs Frame"),
    Dv, min_sep)
  plot_separation_vs_x(time_s, d_geo, "Time (s)",
    paste0(sample_name, "_separation_vs_time.png"),
    paste0(sample_name, " \u2014 Separation vs Time"),
    Dv, min_sep)

  # r/R
  r1_R <- sqrt(x1^2 + y1^2) / R
  r2_R <- sqrt(x2^2 + y2^2) / R
  utils::write.csv(
    data.frame(frame = frame_index, time_s = time_s,
               r1_R = r1_R, r2_R = r2_R),
    paste0(sample_name, "_rR.csv"), row.names = FALSE)
  plot_rR_vs_x(frame_index, r1_R, r2_R, "Frames",
    paste0(sample_name, "_rR_vs_frame.png"),
    paste0(sample_name, " \u2014 r/R vs Frame"))
  plot_rR_vs_x(time_s, r1_R, r2_R, "Time (s)",
    paste0(sample_name, "_rR_vs_time.png"),
    paste0(sample_name, " \u2014 r/R vs Time"))

  # MSSep
  mssep_data <- compute_mssep(d_geo, frame_interval)
  utils::write.csv(mssep_data,
    paste0(sample_name, "_mssep.csv"), row.names = FALSE)
  mssep_fit <- fit_mssep(mssep_data, fit_until = fit_until)
  utils::write.csv(mssep_fit$fit_summary,
    paste0(sample_name, "_mssep_fit.csv"), row.names = FALSE)
  plot_mssep(mssep_data, sample_name)

  message("\n==============================")
  message("Summary: ", sample_name)
  message("==============================")
  message(sprintf("  Valid frames    : %d",       n_frames))
  message(sprintf("  Min separation  : %.4f um",  min(d_geo, na.rm = TRUE)))
  message(sprintf("  Max separation  : %.4f um",  max(d_geo, na.rm = TRUE)))
  message(sprintf("  Mean separation : %.4f um",  mean(d_geo, na.rm = TRUE)))
  message(sprintf("  D_sep           : %.6f um^2/s", mssep_fit$D_sep))

  invisible(list(
    d_geo          = d_geo,
    frame_bins     = bins$frame_bins,
    bin_counts     = bins$bin_counts,
    potential_data = potential_data,
    mssep_data     = mssep_data,
    mssep_fit      = mssep_fit
  ))
}
