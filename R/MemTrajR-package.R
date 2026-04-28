#' MemTrajR: Quantify Micron-Scale Rigid Domain Dynamics and Interactions
#'
#' MemTrajR provides tools for analysing 2D and 3D membrane domain
#' trajectories from single-particle tracking data on spherical vesicles.
#'
#' @section MSD analysis:
#' Use [compute_msd_2d()] and [compute_msd_3d()] to compute mean squared
#' displacement curves, [fit_msd()] to fit them, and [analyze_domain()]
#' as a convenient all-in-one wrapper.
#'
#' @section Separation and interaction potential:
#' Use [compute_geodesic_separation()] to compute frame-by-frame separations,
#' [bin_separations()] to histogram them, [estimate_interaction_potential()]
#' for Boltzmann analysis, [compute_mssep()] and [fit_mssep()] for the mean
#' square separation increment, and [analyze_separation()] as an all-in-one
#' wrapper.
#'
#' @docType package
#' @name MemTrajR-package
#' @aliases MemTrajR
"_PACKAGE"
