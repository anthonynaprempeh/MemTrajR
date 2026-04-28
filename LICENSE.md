# MemTrajR

**MemTrajR** is an R package for quantifying micron-scale rigid domain dynamics and interactions on spherical vesicles from 2D and 3D single-particle trajectory data.

## Installation

Once published on CRAN:

```r
install.packages("MemTrajR")
```

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("anthonynaprempeh/MemTrajR")
```

## Functions

### MSD analysis (`msd.R`)

| Function | Description |
|---|---|
| `compute_msd_2d()` | Euclidean MSD from 2D relative coordinates |
| `compute_msd_3d()` | Geodesic MSD on a sphere of radius *r* |
| `fit_msd()` | Linear + power-law fitting with diffusion coefficient |
| `plot_msd_combined()` | Combined log-log plot of 2D and 3D MSD |
| `plot_displacement_histograms()` | Δx, Δy, Δz histograms at a given lag |
| `analyze_domain()` | All-in-one wrapper for a single domain |

### Separation & interaction potential (`separation.R`)

| Function | Description |
|---|---|
| `compute_geodesic_separation()` | Frame-by-frame geodesic separation |
| `bin_separations()` | Histogram binning of separations |
| `estimate_interaction_potential()` | Boltzmann analysis → U(r)/kBT |
| `compute_mssep()` | Mean square separation increment vs lag |
| `fit_mssep()` | Linear + power-law fit of MSSep |
| `plot_separation_histogram()` | Separation histogram PNG |
| `plot_interaction_potential()` | U(r)/kBT plot PNG |
| `plot_separation_vs_x()` | Separation vs frame or time PNG |
| `plot_rR_vs_x()` | r/R for both domains vs frame or time |
| `plot_mssep()` | MSSep log-log plot PNG |
| `analyze_separation()` | All-in-one wrapper for two-domain analysis |

## Input format

**`analyze_domain()`** expects a CSV with columns: `Dx`, `Dy`, `Vx`, `Vy`
(domain and vesicle centroid coordinates in micrometres, image frame).

**`analyze_separation()`** expects a CSV with columns: `x1`, `y1`, `z1`, `x2`, `y2`, `z2`
(3D coordinates per frame for two domains).

## Quick start

```r
library(MemTrajR)

# Single-domain MSD analysis
df <- read.csv("my_domain.csv")  # columns: Dx, Dy, Vx, Vy

results <- analyze_domain(
  df              = df,
  domain_name     = "domain1",
  frame_interval  = 0.1,
  fit_until       = 1,
  r               = 13.01,
  domain_diameter = 4.38,
  hist_lags       = data.frame(lag_seconds = c(0.2, 1, 10),
                               n_bins      = c(20, 20, 15))
)

# Two-domain separation + interaction potential
sep_results <- analyze_separation(
  input_file     = "two_domains.csv",  # columns: x1,y1,z1,x2,y2,z2
  sample_name    = "sample01",
  R              = 8.02,
  d1             = 3.14,
  d2             = 3.14,
  frame_interval = 0.1,
  fit_until      = 1,
  n_bins         = 20
)
```

## Author

Anthony Prempeh

## License

MIT
