# MemTrajR: MSD analysis and initial slope fitting for one or more 2D trajectories
#
# This script:
# 1. Reads 2D trajectory data from a CSV file
# 2. Computes mean squared displacement (MSD) for each selected domain
# 3. Fits the early-time MSD region with:
#    - a linear model
#    - a power-law model on log-log scales
# 4. Estimates the diffusion coefficient for each domain
# 5. Saves MSD results, fit summaries, and plots for each domain
#
# Each domain is analyzed separately in the same run.

compute_msd <- function(data, max_lag = NULL) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input 'data' must be a data frame or matrix.")
  }
  
  if (ncol(data) < 2) {
    stop("Input 'data' must contain at least two columns for x and y coordinates.")
  }
  
  x <- as.numeric(data[, 1])
  y <- as.numeric(data[, 2])
  
  if (length(x) != length(y)) {
    stop("x and y coordinate vectors must have the same length.")
  }
  
  n_points <- length(x)
  
  if (n_points < 2) {
    stop("At least two observations are required to compute MSD.")
  }
  
  if (is.null(max_lag)) {
    max_lag <- n_points - 1
  }
  
  if (max_lag >= n_points) {
    stop("'max_lag' must be smaller than the number of observations.")
  }
  
  lags <- seq_len(max_lag)
  
  msd_values <- sapply(lags, function(tau) {
    dx <- x[(tau + 1):n_points] - x[1:(n_points - tau)]
    dy <- y[(tau + 1):n_points] - y[1:(n_points - tau)]
    mean(dx^2 + dy^2, na.rm = TRUE)
  })
  
  data.frame(
    lag = lags,
    msd = msd_values
  )
}

fit_msd <- function(msd_data, frame_interval, fit_until = 1) {
  if (!all(c("lag", "msd") %in% names(msd_data))) {
    stop("msd_data must contain columns named 'lag' and 'msd'.")
  }
  
  if (!is.numeric(frame_interval) || length(frame_interval) != 1 || frame_interval <= 0) {
    stop("'frame_interval' must be a single positive number.")
  }
  
  if (!is.numeric(fit_until) || length(fit_until) != 1 || fit_until <= 0) {
    stop("'fit_until' must be a single positive number.")
  }
  
  results <- msd_data
  results$dt <- results$lag * frame_interval
  results$msd_um2 <- results$msd
  
  fit_data <- results[results$dt <= fit_until, ]
  
  if (nrow(fit_data) < 2) {
    stop("Not enough points within the selected fitting range.")
  }
  
  fit_linear <- lm(msd_um2 ~ dt, data = fit_data)
  linear_intercept <- coef(fit_linear)[1]
  linear_slope <- coef(fit_linear)[2]
  diffusion_coefficient <- linear_slope / 4
  
  log_fit_data <- fit_data[fit_data$dt > 0 & fit_data$msd_um2 > 0, ]
  
  if (nrow(log_fit_data) < 2) {
    stop("Not enough positive points for log-log fitting. Ensure dt and MSD are greater than zero.")
  }
  
  fit_log <- lm(log10(msd_um2) ~ log10(dt), data = log_fit_data)
  log_intercept <- coef(fit_log)[1]
  power_exponent <- coef(fit_log)[2]
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
      sprintf("y = %.4f x %+ .4f", linear_slope, linear_intercept),
      sprintf("y = %.4f x^%.5f", power_prefactor, power_exponent),
      as.character(linear_slope),
      as.character(linear_intercept),
      as.character(power_prefactor),
      as.character(power_exponent),
      as.character(diffusion_coefficient)
    )
  )
  
  list(
    results = results,
    fit_data = fit_data,
    log_fit_data = log_fit_data,
    fit_linear = fit_linear,
    fit_log = fit_log,
    fit_summary = fit_summary,
    diffusion_coefficient = diffusion_coefficient,
    power_prefactor = power_prefactor,
    power_exponent = power_exponent,
    linear_slope = linear_slope,
    linear_intercept = linear_intercept,
    log_intercept = log_intercept
  )
}

plot_msd <- function(results, plot_title, output_plot,
                     x_min = 0.1, y_min = 0.01, y_max = 100) {
  x_max <- max(results$dt, na.rm = TRUE)
  point_colors <- rep("black", nrow(results))
  
  png(filename = output_plot, width = 1800, height = 1400, res = 220)
  
  par(mgp = c(2, 0.7, 0))
  
  plot(
    results$dt, results$msd_um2,
    log = "xy",
    pch = 16,
    col = point_colors,
    xlim = c(x_min, x_max),
    ylim = c(y_min, y_max),
    xaxt = "n",
    yaxt = "n",
    xaxs = "i",
    yaxs = "i",
    bty = "l",
    cex.lab = 1.5,
    font.lab = 2,
    cex.main = 1.6,
    font.main = 2,
    main = plot_title,
    xlab = "dt (s)",
    ylab = expression(bold(MSD ~ "(" * mu * "m"^2 * ")"))
  )
  
  x_ticks <- c(0.1, 1, 10, 100, 1000)
  x_ticks <- x_ticks[x_ticks >= x_min & x_ticks <= x_max]
  
  y_ticks <- c(0.01, 0.1, 1, 10, 100)
  y_ticks <- y_ticks[y_ticks >= y_min & y_ticks <= y_max]
  
  axis(1, at = x_ticks, labels = x_ticks)
  axis(2, at = y_ticks, labels = y_ticks)
  
  dev.off()
}

analyze_domain <- function(df, coord_cols, domain_name,
                           frame_interval, fit_until = 1,
                           x_min = 0.1, y_min = 0.01, y_max = 100) {
  if (length(coord_cols) != 2) {
    stop(paste("Domain", domain_name, "must have exactly two coordinate columns."))
  }
  
  coords <- df[, coord_cols]
  msd_data <- compute_msd(coords)
  fit_results <- fit_msd(msd_data, frame_interval = frame_interval, fit_until = fit_until)
  
  output_msd <- paste0(domain_name, "_msd.csv")
  output_fit <- paste0(domain_name, "_fit.csv")
  output_plot <- paste0(domain_name, "_msd_plot.png")
  
  write.csv(msd_data, output_msd, row.names = FALSE)
  write.csv(fit_results$fit_summary, output_fit, row.names = FALSE)
  
  plot_msd(
    results = fit_results$results,
    plot_title = domain_name,
    output_plot = output_plot,
    x_min = x_min,
    y_min = y_min,
    y_max = y_max
  )
  
  cat("\n==============================\n")
  cat("Results for", domain_name, "\n")
  cat("==============================\n")
  
  cat("Linear equation:\n")
  cat(sprintf("y = %.4f x %+ .4f\n", fit_results$linear_slope, fit_results$linear_intercept))
  
  cat("\nPower-law equation:\n")
  cat(sprintf("y = %.4f x^%.5f\n", fit_results$power_prefactor, fit_results$power_exponent))
  
  cat("\nLog-log equation:\n")
  cat(sprintf(
    "log10(y) = %.4f + %.4f log10(dt)\n",
    fit_results$log_intercept,
    fit_results$power_exponent
  ))
  
  cat("\nDiffusion coefficient:\n")
  cat(sprintf("D = %.6f um^2/s\n", fit_results$diffusion_coefficient))
  
  cat("\nFiles written:\n")
  cat(sprintf("- MSD results: %s\n", output_msd))
  cat(sprintf("- Fit summary: %s\n", output_fit))
  cat(sprintf("- Plot: %s\n", output_plot))
  
  list(
    domain_name = domain_name,
    msd_data = msd_data,
    fit_results = fit_results
  )
}

# ==============================
# User input section
# ==============================

# Enter the name of your input CSV file
input_file <- "0012.csv"

# Set the frame interval in seconds
frame_interval <- 0.1

# Set the maximum time range to use for fitting, in seconds
fit_until <- 1

# Optional axis limits for the log-log plots
x_min <- 0.1
y_min <- 0.01
y_max <- 100

# Define each domain as a pair of x and y columns
# Example:
# domain_list <- list(
#   domain_1 = c(7, 8),
#   domain_2 = c(9, 10)
# )

domain_list <- list(
  domain_1 = c(7, 8),
  domain_2 = c(10, 11)
)

# ==============================
# Run analysis
# ==============================

df <- read.csv(input_file)

all_results <- lapply(names(domain_list), function(domain_name) {
  analyze_domain(
    df = df,
    coord_cols = domain_list[[domain_name]],
    domain_name = domain_name,
    frame_interval = frame_interval,
    fit_until = fit_until,
    x_min = x_min,
    y_min = y_min,
    y_max = y_max
  )
})

names(all_results) <- names(domain_list)