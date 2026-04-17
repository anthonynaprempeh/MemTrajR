# MemTrajR: frame-by-frame separation and interaction-potential analysis for two 2D domains
#
# This script:
# 1. Reads 2D trajectory data for two domains from a CSV file
# 2. Computes the separation distance between the domains for each frame
# 3. Assigns each frame to a separation bin
# 4. Exports frame-by-frame separation data and histogram bin counts
# 5. Estimates the interaction potential using:
#       U(r) / (k_B T) = -ln[P(r)]
#    where P(r) is defined relative to the last bin count
# 6. Saves a separation histogram and an interaction-potential plot
#
# Before running:
# - Place this script in your project folder, or run it from the folder containing your data
# - Set the input/output file names below
# - Specify the column names for x and y coordinates of both domains
# - Choose either a fixed bin width or a fixed number of bins

# ==============================
# User input section
# ==============================

# Enter the name of your input CSV file
input_file <- "0004.csv"

# Enter names for output files
output_frame_data <- "0004_2d_frame_bins.csv"
output_bin_counts <- "0004_2d_bin_counts.csv"
output_histogram <- "0004_2d_separation_histogram.png"
output_potential_data <- "0004_2d_interaction_potential.csv"
output_potential_plot <- "0004_2d_interaction_potential.png"

# Set plot titles
histogram_title <- "0004 2D separation histogram"
potential_title <- "0004 2D interaction potential"

# Specify the column names for the two domains
x1_col <- "x1"
y1_col <- "y1"
x2_col <- "x2"
y2_col <- "y2"

# Choose binning mode:
# "fixed_width" for a fixed bin width
# "fixed_count" for a fixed number of bins
bin_mode <- "fixed_count"

# Option A: fixed bin width
bin_width <- 0.2

# Option B: fixed number of bins
n_bins <- 20

# ==============================
# Read data
# ==============================

dat <- read.csv(input_file)

required_cols <- c(x1_col, y1_col, x2_col, y2_col)

if (!all(required_cols %in% names(dat))) {
  stop("One or more specified coordinate columns were not found in the input file.")
}

# Keep only frames where both domains have valid coordinates
valid <- complete.cases(dat[, required_cols])

if (sum(valid) == 0) {
  stop("No frames contain complete coordinate data for both domains.")
}

frame_index <- which(valid)

x1 <- dat[valid, x1_col]
y1 <- dat[valid, y1_col]
x2 <- dat[valid, x2_col]
y2 <- dat[valid, y2_col]

# ==============================
# Compute separation per frame
# ==============================

separation_um <- sqrt((x1 - x2)^2 + (y1 - y2)^2)

# ==============================
# Define histogram bins
# ==============================

if (bin_mode == "fixed_width") {
  start <- floor(min(separation_um, na.rm = TRUE) / bin_width) * bin_width
  end <- ceiling(max(separation_um, na.rm = TRUE) / bin_width) * bin_width
  edges <- seq(start, end, by = bin_width)
} else if (bin_mode == "fixed_count") {
  start <- min(separation_um, na.rm = TRUE)
  end <- max(separation_um, na.rm = TRUE)
  edges <- seq(start, end, length.out = n_bins + 1)
} else {
  stop("bin_mode must be either 'fixed_width' or 'fixed_count'.")
}

if (length(edges) < 2) {
  stop("Could not generate valid histogram bin edges.")
}

bin_label <- cut(separation_um, breaks = edges, include.lowest = TRUE, right = FALSE)

# ==============================
# Export frame-by-frame separation data
# ==============================

frame_bins <- data.frame(
  frame = frame_index,
  x1 = x1,
  y1 = y1,
  x2 = x2,
  y2 = y2,
  separation_um = separation_um,
  bin = as.character(bin_label)
)

write.csv(frame_bins, output_frame_data, row.names = FALSE)

# ==============================
# Export bin counts
# ==============================

bin_counts <- as.integer(table(bin_label))
bin_names <- levels(bin_label)

counts <- data.frame(
  bin = bin_names,
  count_frames = bin_counts
)

write.csv(counts, output_bin_counts, row.names = FALSE)

# ==============================
# Compute bin midpoints and interaction potential
# ==============================

# Bin midpoints are computed from the histogram edges
bin_midpoints <- (head(edges, -1) + tail(edges, -1)) / 2

# Reference state:
# P(r) is defined relative to the last bin count, so the last bin is the
# zero-potential reference state.
last_bin_count <- tail(bin_counts, 1)

if (last_bin_count == 0) {
  stop("The last bin has zero counts, so P(r) cannot be computed using the last bin as the reference.")
}

P_r <- bin_counts / last_bin_count
U_over_kBT <- ifelse(P_r > 0, -log(P_r), NA)

potential_data <- data.frame(
  bin = bin_names,
  bin_midpoint_um = bin_midpoints,
  count_frames = bin_counts,
  P_r = P_r,
  U_over_kBT = U_over_kBT
)

write.csv(potential_data, output_potential_data, row.names = FALSE)

# ==============================
# Save separation histogram
# ==============================

png(filename = output_potential_plot, width = 1800, height = 1400, res = 220)

par(
  bg = "gray92",
  fg = "gray35",
  col.axis = "gray25",
  col.lab = "gray25",
  mgp = c(2.2, 0.8, 0),
  mar = c(5, 5, 2, 2)
)

plot(
  potential_data$bin_midpoint_um,
  potential_data$U_over_kBT,
  type = "l",
  lwd = 4,
  col = "black",
  xaxt = "n",
  yaxt = "n",
  xlab = "",
  ylab = expression("Interaction Potential (k"[B] * "T)"),
  bty = "l",
  ylim = c(-3, 0),
  xlim = range(potential_data$bin_midpoint_um, na.rm = TRUE)
)

axis(3, at = pretty(potential_data$bin_midpoint_um), labels = TRUE, lwd = 1.2)
axis(2, at = seq(-3, 0, by = 0.5), las = 1, lwd = 1.2)

text(
  x = mean(range(potential_data$bin_midpoint_um, na.rm = TRUE)),
  y = -0.45,
  labels = expression("Centre to Centre Separations (" * mu * "m)"),
  font = 2,
  cex = 1.2,
  col = "gray35"
)

dev.off()

# ==============================
# Save interaction-potential plot
# ==============================

png(filename = output_potential_plot, width = 1800, height = 1400, res = 220)

par(mgp = c(2, 0.7, 0))

plot(
  potential_data$bin_midpoint_um,
  potential_data$U_over_kBT,
  type = "l",
  lwd = 3,
  pch = 16,
  xlab = expression("Centre to Centre Separations (" * mu * "m)"),
  ylab = expression("Interaction Potential (k"[B] * "T)"),
  main = potential_title,
  bty = "l"
)

dev.off()

# ==============================
# Print summary
# ==============================

cat("Separation and interaction-potential analysis complete.\n\n")
cat(sprintf("Number of valid frames: %d\n", length(separation_um)))
cat(sprintf("Minimum separation: %.4f um\n", min(separation_um, na.rm = TRUE)))
cat(sprintf("Maximum separation: %.4f um\n", max(separation_um, na.rm = TRUE)))
cat(sprintf("Mean separation: %.4f um\n", mean(separation_um, na.rm = TRUE)))

cat("\nFiles written:\n")
cat(sprintf("- Frame-by-frame separation data: %s\n", output_frame_data))
cat(sprintf("- Histogram bin counts: %s\n", output_bin_counts))
cat(sprintf("- Separation histogram: %s\n", output_histogram))
cat(sprintf("- Interaction potential data: %s\n", output_potential_data))
cat(sprintf("- Interaction potential plot: %s\n", output_potential_plot))