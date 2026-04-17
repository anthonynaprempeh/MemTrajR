# MemTrajR

MemTrajR is an R-based toolkit for membrane domain analysis from 2D trajectory data.

## Features

- Mean squared displacement (MSD) analysis
- Early-time linear and power-law fitting
- Diffusion coefficient estimation
- Frame-by-frame domain separation analysis
- Separation histograms
- Interaction potential estimation using Boltzmann analysis

## Scripts

### `MemTrajR_MSD.R`
Computes MSD curves for one or more 2D trajectories, performs early-time fitting, and estimates diffusion coefficients.

### `MemTrajR_Separation_&_Interaction potential.R`
Computes frame-by-frame separations between two domains, bins the separation distances, and estimates interaction potentials.

## Input
CSV files containing x and y coordinate data for tracked membrane domains.

## Output
Depending on the script, outputs may include:
- MSD tables
- fit summaries
- separation tables
- histogram bin counts
- interaction potential tables
- PNG plots

## Requirements
R base functions only.

## Author
Anthony Prempeh
