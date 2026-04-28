# msd.R — MSD computation, fitting, and plotting functions
# Derived from: 2_3D_MSD_with_data_and_plots___xyz_displacements.R


# ============================================================
#  MSD computation
# ============================================================

#' Compute 2D Euclidean MSD
#'
#' Computes the mean squared displacement using flat Euclidean distances
#' for a single domain trajectory.
#'
#' @param relx Numeric vector of relative x coordinates (micrometres).
#' @param rely Numeric vector of relative y coordinates (micrometres).
#' @param max_lag Integer. Maximum lag to compute. Defaults to
#'   `length(relx) - 1`.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{lag}{Lag index (integer frames).}
#'     \item{msd}{Mean squared displacement (micrometres squared).}
#'   }
#'
#' @export
#' @examples
#' relx <- cumsum(rnorm(100, sd = 0.1))
#' rely <- cumsum(rnorm(100, sd = 0.1))
#' msd  <- compute_msd_2d(relx, rely, max_lag = 20)
#' head(msd)
compute_msd_2d <- function(relx, rely, max_lag = NULL) {
  n_points <- length(relx)

  if (n_points < 2)
    stop("At least two observations are required to compute MSD.")

  if (is.null(max_lag)) max_lag <- n_points - 1

  if (max_lag >= n_points)
    stop("'max_lag' must be smaller than the number of observations.")

  lags <- seq_len(max_lag)

  msd_values <- sapply(lags, function(tau) {
    dx <- relx[(tau + 1):n_points] - relx[1:(n_points - tau)]
    dy <- rely[(tau + 1):n_points] - rely[1:(n_points - tau)]
    mean(dx^2 + dy^2, na.rm = TRUE)
  })

  data.frame(lag = lags, msd = msd_values)
}


#' Compute 3D Geodesic MSD
#'
#' Computes the mean squared displacement using geodesic (great-circle)
#' distances on a sphere of radius `r`.
#'
#' @param relx Numeric vector of relative x coordinates (micrometres).
#' @param rely Numeric vector of relative y coordinates (micrometres).
#' @param r Positive numeric. Sphere (vesicle) radius in micrometres.
#' @param max_lag Integer. Maximum lag to compute. Defaults to
#'   `length(relx) - 1`.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{lag}{Lag index (integer frames).}
#'     \item{msd}{Mean squared geodesic displacement (micrometres squared).}
#'   }
#'
#' @export
#' @examples
#' r    <- 10
#' relx <- rnorm(100, sd = 0.5)
#' rely <- rnorm(100, sd = 0.5)
#' # keep points inside sphere
#' keep <- (relx^2 + rely^2) < r^2
#' msd  <- compute_msd_3d(relx[keep], rely[keep], r = r, max_lag = 20)
#' head(msd)
compute_msd_3d <- function(relx, rely, r, max_lag = NULL) {
  z_sq <- r^2 - relx^2 - rely^2
  if (any(z_sq < 0, na.rm = TRUE))
    stop(paste(
      "Some (relx, rely) points lie outside the sphere of radius r.",
      "Check your coordinates or radius."
    ))

  z        <- sqrt(z_sq)
  n_points <- length(relx)

  if (n_points < 2)
    stop("At least two observations are required to compute MSD.")

  if (is.null(max_lag)) max_lag <- n_points - 1

  if (max_lag >= n_points)
    stop("'max_lag' must be smaller than the number of observations.")

  lags <- seq_len(max_lag)

  msd_values <- sapply(lags, function(tau) {
    i_from <- 1:(n_points - tau)
    i_to   <- (tau + 1):n_points

    dot12 <- relx[i_from] * relx[i_to] +
              rely[i_from] * rely[i_to] +
              z[i_from]    * z[i_to]

    cosg  <- pmin(pmax(dot12 / (r^2), -1), 1)
    geo_dist <- r * acos(cosg)
    mean(geo_dist^2, na.rm = TRUE)
  })

  data.frame(lag = lags, msd = msd_values)
}


# ============================================================
#  MSD fitting
# ============================================================

#' Fit MSD Curve (Linear and Power-Law)
#'
#' Fits a raw MSD table (from [compute_msd_2d()] or [compute_msd_3d()])
#' with both a linear model and a log-log power-law model over an
#' early-time fitting window.
#'
#' @param msd_data Data frame with columns `lag` and `msd`, as returned
#'   by [compute_msd_2d()] or [compute_msd_3d()].
#' @param frame_interval Positive numeric. Time between frames in seconds.
#' @param fit_until Positive numeric. Upper time limit (seconds) for the
#'   fitting window. Default `1`.
#'
#' @return A named list containing:
#'   \describe{
#'     \item{results}{Data frame with added columns `dt` and `msd_um2`.}
#'     \item{fit_summary}{Data frame summarising fit parameters.}
#'     \item{diffusion_coefficient}{Estimated diffusion coefficient
#'       (micrometres squared per second), derived as slope / 4.}
#'     \item{power_prefactor}{Power-law prefactor A.}
#'     \item{power_exponent}{Power-law exponent n.}
#'     \item{linear_slope}{Slope of the linear fit.}
#'     \item{linear_intercept}{Intercept of the linear fit.}
#'     \item{fit_linear}{The `lm` object for the linear fit.}
#'     \item{fit_log}{The `lm` object for the log-log fit.}
#'   }
#'
#' @export
#' @examples
#' relx <- cumsum(rnorm(200, sd = 0.1))
#' rely <- cumsum(rnorm(200, sd = 0.1))
#' msd  <- compute_msd_2d(relx, rely, max_lag = 50)
#' fit  <- fit_msd(msd, frame_interval = 0.1, fit_until = 1)
#' fit$fit_summary
fit_msd <- function(msd_data, frame_interval, fit_until = 1) {
  results         <- msd_data
  results$dt      <- results$lag * frame_interval
  results$msd_um2 <- results$msd

  fit_data <- results[results$dt <= fit_until, ]

  if (nrow(fit_data) < 2)
    stop("Not enough points within the selected fitting range.")

  fit_linear            <- stats::lm(msd_um2 ~ dt, data = fit_data)
  linear_intercept      <- stats::coef(fit_linear)[1]
  linear_slope          <- stats::coef(fit_linear)[2]
  diffusion_coefficient <- linear_slope / 4

  log_fit_data <- fit_data[fit_data$dt > 0 & fit_data$msd_um2 > 0, ]

  if (nrow(log_fit_data) < 2)
    stop("Not enough positive points for log-log fitting.")

  fit_log         <- stats::lm(log10(msd_um2) ~ log10(dt), data = log_fit_data)
  log_intercept   <- stats::coef(fit_log)[1]
  power_exponent  <- stats::coef(fit_log)[2]
  power_prefactor <- 10^log_intercept

  fit_summary <- data.frame(
    parameter = c(
      "linear_equation",
      "power_law_equation",
      "linear_slope",
      "linear_intercept",
      "power_prefactor",
      "power_exponent",
      "diffusion_coefficient_um2_s"
    ),
    value = c(
      sprintf("y = %.4f x %+.4f", linear_slope, linear_intercept),
      sprintf("y = %.4f x^%.5f",  power_prefactor, power_exponent),
      as.character(linear_slope),
      as.character(linear_intercept),
      as.character(power_prefactor),
      as.character(power_exponent),
      as.character(diffusion_coefficient)
    ),
    stringsAsFactors = FALSE
  )

  list(
    results               = results,
    fit_data              = fit_data,
    log_fit_data          = log_fit_data,
    fit_linear            = fit_linear,
    fit_log               = fit_log,
    fit_summary           = fit_summary,
    diffusion_coefficient = diffusion_coefficient,
    power_prefactor       = power_prefactor,
    power_exponent        = power_exponent,
    linear_slope          = linear_slope,
    linear_intercept      = linear_intercept,
    log_intercept         = log_intercept
  )
}


# ============================================================
#  Displacement histograms
# ============================================================

#' Plot 1D Displacement Histograms (Delta x, Delta y, Delta z)
#'
#' For a given lag, computes frame-to-frame displacements in x, y, and z,
#' saves a three-panel PNG histogram, and writes the displacements to CSV.
#'
#' @param relx Numeric vector of relative x coordinates (micrometres).
#' @param rely Numeric vector of relative y coordinates (micrometres).
#' @param z Numeric vector of z coordinates on the sphere (micrometres).
#' @param lag Positive integer. Lag in frames.
#' @param n_bins Positive integer. Number of histogram bins.
#' @param frame_interval Positive numeric. Time between frames (seconds).
#' @param domain_name Character. Prefix for output file names.
#'
#' @return Invisibly, a list with elements `dx`, `dy`, `dz` (numeric
#'   displacement vectors). Output files are written as a side effect.
#'
#' @export
#' @examples
#' \dontrun{
#' r    <- 10
#' relx <- rnorm(200, sd = 0.3)
#' rely <- rnorm(200, sd = 0.3)
#' z    <- sqrt(r^2 - relx^2 - rely^2)
#' plot_displacement_histograms(relx, rely, z,
#'   lag = 2, n_bins = 20,
#'   frame_interval = 0.1, domain_name = "domain1")
#' }
plot_displacement_histograms <- function(relx, rely, z, lag, n_bins,
                                         frame_interval, domain_name) {
  n_points <- length(relx)

  if (lag >= n_points)
    stop(sprintf(
      "'hist_lag' (%d) must be smaller than the number of observations (%d).",
      lag, n_points
    ))

  i_from <- 1:(n_points - lag)
  i_to   <- (lag + 1):n_points

  dx <- relx[i_to] - relx[i_from]
  dy <- rely[i_to]  - rely[i_from]
  dz <- z[i_to]     - z[i_from]

  n_pairs  <- n_points - lag
  dt_label <- sprintf("dt = %.2f s (lag = %d, n = %d, bins = %d)",
                      lag * frame_interval, lag, n_pairs, n_bins)

  output_plot <- paste0(domain_name, "_displacement_hist_lag", lag, ".png")

  grDevices::png(filename = output_plot, width = 2600, height = 1000, res = 220)
  graphics::par(mfrow = c(1, 3), mgp = c(2, 0.7, 0), mar = c(4, 4, 3, 1))

  breaks_x <- seq(min(dx), max(dx), length.out = n_bins + 1)
  breaks_y <- seq(min(dy), max(dy), length.out = n_bins + 1)
  breaks_z <- seq(min(dz), max(dz), length.out = n_bins + 1)

  graphics::hist(dx, breaks = breaks_x, col = "steelblue", border = "white",
    main = bquote(bold(Delta * x ~ .(dt_label))),
    xlab = expression(bold(Delta * x ~ "(" * mu * "m)")),
    ylab = "Frequency", cex.lab = 1.2, font.lab = 2,
    cex.main = 1.1, font.main = 2)

  graphics::hist(dy, breaks = breaks_y, col = "firebrick", border = "white",
    main = bquote(bold(Delta * y ~ .(dt_label))),
    xlab = expression(bold(Delta * y ~ "(" * mu * "m)")),
    ylab = "Frequency", cex.lab = 1.2, font.lab = 2,
    cex.main = 1.1, font.main = 2)

  graphics::hist(dz, breaks = breaks_z, col = "forestgreen", border = "white",
    main = bquote(bold(Delta * z ~ .(dt_label))),
    xlab = expression(bold(Delta * z ~ "(" * mu * "m)")),
    ylab = "Frequency", cex.lab = 1.2, font.lab = 2,
    cex.main = 1.1, font.main = 2)

  grDevices::dev.off()

  output_csv <- paste0(domain_name, "_displacement_lag", lag, ".csv")
  utils::write.csv(data.frame(dx = dx, dy = dy, dz = dz),
                   output_csv, row.names = FALSE)

  message(sprintf("  Displacement histogram : %s", output_plot))
  message(sprintf("  Displacement CSV       : %s", output_csv))

  invisible(list(dx = dx, dy = dy, dz = dz))
}


# ============================================================
#  Combined MSD plot
# ============================================================

#' Plot Combined 2D and 3D MSD on Log-Log Axes
#'
#' Produces a single PNG showing both the Euclidean (2D) and geodesic (3D)
#' MSD curves on logarithmic axes, with optional vertical dashed lines
#' marking histogram lag times.
#'
#' @param results_2d Data frame (from [fit_msd()]\code{$results}) with
#'   columns `dt` and `msd_um2`.
#' @param results_3d Data frame (from [fit_msd()]\code{$results}) with
#'   columns `dt` and `msd_um2`.
#' @param plot_title Character. Title printed above the plot.
#' @param output_plot Character. File path for the output PNG.
#' @param y_min Positive numeric. Lower y-axis limit. Default `0.01`.
#' @param y_max Positive numeric. Upper y-axis limit. Default `100`.
#' @param hist_lag_seconds Numeric vector of lag times (seconds) to mark
#'   with vertical dashed lines, or `NULL` (default) for none.
#'
#' @return Invisibly `NULL`. A PNG file is written as a side effect.
#'
#' @export
#' @examples
#' \dontrun{
#' relx <- cumsum(rnorm(200, sd = 0.1))
#' rely <- cumsum(rnorm(200, sd = 0.1))
#' msd2 <- compute_msd_2d(relx, rely, max_lag = 50)
#' msd3 <- compute_msd_3d(relx, rely, r = 10, max_lag = 50)
#' f2   <- fit_msd(msd2, frame_interval = 0.1)
#' f3   <- fit_msd(msd3, frame_interval = 0.1)
#' plot_msd_combined(f2$results, f3$results,
#'   plot_title  = "Domain 1",
#'   output_plot = "domain1_msd.png")
#' }
plot_msd_combined <- function(results_2d, results_3d, plot_title, output_plot,
                               y_min = 0.01, y_max = 100,
                               hist_lag_seconds = NULL) {
  x_min_data <- min(c(results_2d$dt, results_3d$dt), na.rm = TRUE)
  x_max_data <- max(c(results_2d$dt, results_3d$dt), na.rm = TRUE)
  x_min_plot <- 10^floor(log10(x_min_data))
  x_max_plot <- 10^ceiling(log10(x_max_data))

  y_min_data <- min(c(results_2d$msd_um2, results_3d$msd_um2), na.rm = TRUE)
  y_max_data <- max(c(results_2d$msd_um2, results_3d$msd_um2), na.rm = TRUE)
  y_min_plot <- min(y_min, 10^floor(log10(y_min_data)))
  y_max_plot <- max(y_max, 10^ceiling(log10(y_max_data)))

  x_major <- 10^seq(log10(x_min_plot), log10(x_max_plot))
  y_major <- 10^seq(log10(y_min_plot), log10(y_max_plot))
  x_minor <- .log_minor_ticks(x_major)
  y_minor <- .log_minor_ticks(y_major)
  x_minor <- x_minor[x_minor > x_min_plot & x_minor < x_max_plot]
  y_minor <- y_minor[y_minor > y_min_plot & y_minor < y_max_plot]

  grDevices::png(filename = output_plot, width = 1800, height = 1400, res = 220)
  graphics::par(mgp = c(2, 0.7, 0))

  graphics::plot(
    results_2d$dt, results_2d$msd_um2,
    log = "xy", pch = 16, col = "steelblue",
    xlim = c(x_min_plot, x_max_plot),
    ylim = c(y_min_plot, y_max_plot),
    xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", bty = "l",
    cex.lab = 1.5, font.lab = 2, cex.main = 1.6, font.main = 2,
    main = plot_title,
    xlab = "dt (s)",
    ylab = expression(bold(MSD ~ "(" * mu * "m"^2 * ")"))
  )

  graphics::points(results_3d$dt, results_3d$msd_um2, pch = 17, col = "firebrick")

  if (!is.null(hist_lag_seconds)) {
    valid_lags <- hist_lag_seconds[
      hist_lag_seconds >= x_min_plot & hist_lag_seconds <= x_max_plot]
    if (length(valid_lags) > 0) {
      graphics::abline(v = valid_lags, lty = 2, col = "gray40", lwd = 1.2)
      graphics::axis(side = 3, at = valid_lags,
                     labels = paste0(valid_lags, "s"),
                     tcl = -0.3, cex.axis = 0.8,
                     col.axis = "gray40", font = 2, las = 2)
    }
  }

  graphics::legend("topleft",
    legend = c("2D (Euclidean)", "3D (Geodesic)"),
    pch = c(16, 17), col = c("steelblue", "firebrick"),
    bty = "n", cex = 1.1)

  graphics::axis(1, at = x_major, labels = x_major, tcl = -0.5)
  graphics::axis(2, at = y_major, labels = y_major, tcl = -0.5)
  graphics::axis(1, at = x_minor, labels = FALSE,   tcl = -0.25)
  graphics::axis(2, at = y_minor, labels = FALSE,   tcl = -0.25)

  grDevices::dev.off()
  invisible(NULL)
}


# ============================================================
#  Main wrapper: analyze_domain
# ============================================================

#' Analyse a Single Domain Trajectory (2D and 3D MSD)
#'
#' High-level wrapper that reads relative domain coordinates, computes both
#' 2D Euclidean and 3D geodesic MSD, fits each curve, generates combined
#' log-log plots, and optionally plots per-lag displacement histograms.
#' All outputs are written to the current working directory.
#'
#' @param df Data frame with columns `Dx`, `Dy`, `Vx`, `Vy` (domain and
#'   vesicle centroid coordinates in micrometres, image frame).
#' @param domain_name Character. Prefix used for all output file names.
#' @param frame_interval Positive numeric. Time between frames (seconds).
#' @param fit_until Positive numeric. Upper time limit (seconds) for the
#'   MSD fitting window. Default `1`.
#' @param r Positive numeric. Vesicle radius in micrometres.
#' @param domain_diameter Positive numeric. Domain diameter in micrometres,
#'   used to compute the characteristic diffusion time
#'   \eqn{\tau = d^2 / D_{3D}}.
#' @param hist_lags Data frame with columns `lag_seconds` (numeric) and
#'   `n_bins` (integer), specifying the lag times and bin counts for
#'   displacement histograms.
#' @param y_min Positive numeric. Lower y-axis limit for MSD plots.
#'   Default `0.01`.
#' @param y_max Positive numeric. Upper y-axis limit for MSD plots.
#'   Default `100`.
#'
#' @return Invisibly, a named list with elements `domain_name`, `relx`,
#'   `rely`, `msd_2d`, `fit_2d`, `msd_3d`, `fit_3d`, and `tau`. Output
#'   CSVs and PNGs are written as side effects.
#'
#' @export
#' @examples
#' \dontrun{
#' df <- read.csv("0022_d1.csv")
#' hist_lags <- data.frame(
#'   lag_seconds = c(0.2, 1, 10),
#'   n_bins      = c(20, 20, 15)
#' )
#' results <- analyze_domain(
#'   df             = df,
#'   domain_name    = "0022_d1",
#'   frame_interval = 0.1,
#'   fit_until      = 1,
#'   r              = 13.01,
#'   domain_diameter = 4.38,
#'   hist_lags      = hist_lags
#' )
#' }
analyze_domain <- function(df, domain_name,
                            frame_interval, fit_until = 1,
                            r, domain_diameter,
                            hist_lags,
                            y_min = 0.01, y_max = 100) {

  required_cols <- c("Dx", "Dy", "Vx", "Vy")
  missing_cols  <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0)
    stop(paste(
      "Missing required columns:", paste(missing_cols, collapse = ", "),
      "\nPlease name your CSV columns: Dx, Dy, Vx, Vy."
    ))

  if (!is.numeric(r) || r <= 0)
    stop("'r' must be a positive number (sphere radius in micrometres).")

  if (!is.numeric(domain_diameter) || domain_diameter <= 0)
    stop("'domain_diameter' must be a positive number (in micrometres).")

  relx <- df$Dx - df$Vx
  rely <- df$Vy - df$Dy
  z    <- sqrt(r^2 - relx^2 - rely^2)

  output_rel <- paste0(domain_name, "_rel_coords.csv")
  utils::write.csv(data.frame(relx = relx, rely = rely, z = z),
                   output_rel, row.names = FALSE)

  # 2D
  msd_2d <- compute_msd_2d(relx, rely)
  fit_2d <- fit_msd(msd_2d, frame_interval = frame_interval,
                    fit_until = fit_until)
  utils::write.csv(msd_2d,             paste0(domain_name, "_2d_msd.csv"),
                   row.names = FALSE)
  utils::write.csv(fit_2d$fit_summary, paste0(domain_name, "_2d_fit.csv"),
                   row.names = FALSE)

  # 3D
  msd_3d <- compute_msd_3d(relx, rely, r = r)
  fit_3d <- fit_msd(msd_3d, frame_interval = frame_interval,
                    fit_until = fit_until)
  tau    <- domain_diameter^2 / fit_3d$diffusion_coefficient

  fit_3d_summary_with_tau <- rbind(
    fit_3d$fit_summary,
    data.frame(
      parameter = c("domain_diameter_um",
                    "tau_characteristic_diffusion_time_s"),
      value     = c(as.character(domain_diameter), as.character(tau)),
      stringsAsFactors = FALSE
    )
  )
  utils::write.csv(msd_3d,                  paste0(domain_name, "_3d_msd.csv"),
                   row.names = FALSE)
  utils::write.csv(fit_3d_summary_with_tau, paste0(domain_name, "_3d_fit.csv"),
                   row.names = FALSE)

  # Plots
  plot_msd_combined(fit_2d$results, fit_3d$results,
    plot_title  = domain_name,
    output_plot = paste0(domain_name, "_msd_plot.png"),
    y_min = y_min, y_max = y_max, hist_lag_seconds = NULL)

  plot_msd_combined(fit_2d$results, fit_3d$results,
    plot_title  = domain_name,
    output_plot = paste0(domain_name, "_msd_plot_marked.png"),
    y_min = y_min, y_max = y_max,
    hist_lag_seconds = hist_lags$lag_seconds)

  # Displacement histograms
  hist_lags$lag_frames <- round(hist_lags$lag_seconds / frame_interval)
  invalid              <- hist_lags$lag_frames < 1
  if (any(invalid))
    warning(sprintf(
      "Some lag times are smaller than one frame interval (%.2f s) and will be skipped.",
      frame_interval))
  hist_lags <- hist_lags[!invalid, ]
  hist_lags <- hist_lags[!duplicated(hist_lags$lag_frames), ]

  message("\nDisplacement histograms:")
  lapply(seq_len(nrow(hist_lags)), function(i) {
    plot_displacement_histograms(
      relx = relx, rely = rely, z = z,
      lag            = hist_lags$lag_frames[i],
      n_bins         = hist_lags$n_bins[i],
      frame_interval = frame_interval,
      domain_name    = domain_name
    )
  })

  message("\n==============================")
  message("Results for ", domain_name)
  message("==============================")
  message(sprintf("  2D  D : %.6f um^2/s", fit_2d$diffusion_coefficient))
  message(sprintf("  3D  D : %.6f um^2/s", fit_3d$diffusion_coefficient))
  message(sprintf("  tau   : %.4f s",       tau))

  invisible(list(
    domain_name = domain_name,
    relx = relx, rely = rely,
    msd_2d = msd_2d, fit_2d = fit_2d,
    msd_3d = msd_3d, fit_3d = fit_3d,
    tau = tau
  ))
}
