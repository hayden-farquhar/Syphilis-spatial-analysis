# ==============================================================================
# 00_setup.R — Package loading, paths, and global options
# Syphilis Surge Spatial Analysis
# ==============================================================================

# ---- Install missing packages ------------------------------------------------
required_pkgs <- c(

  # Spatial
  "sf", "spdep", "spatialreg", "spgwr", "tmap",
  # Data wrangling
  "tidyverse", "readxl", "janitor", "lubridate",
  # Visualisation
  "patchwork", "viridis", "RColorBrewer",
  # Tables
  "gt", "gtsummary"
)

installed <- rownames(installed.packages())
to_install <- setdiff(required_pkgs, installed)
if (length(to_install) > 0) {
  message("Installing: ", paste(to_install, collapse = ", "))
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

# ---- Load packages -----------------------------------------------------------
library(tidyverse)
library(sf)
library(readxl)
library(janitor)
library(lubridate)
library(spdep)
library(spatialreg)
library(tmap)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(gt)

# ---- Project paths -----------------------------------------------------------
# Determine project root: try here::here(), then script location, then getwd()
proj_root <- tryCatch(here::here(), error = function(e) NULL)
if (is.null(proj_root) || !file.exists(file.path(proj_root, "CLAUDE.md"))) {
  # Try to find from script location (works when sourced)
  script_dir <- tryCatch({
    # When sourced from another script
    if (sys.nframe() > 0) dirname(sys.frame(1)$ofile)
    else NULL
  }, error = function(e) NULL)

  if (!is.null(script_dir)) {
    proj_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
  }

  # Final fallback: assume working directory is project root
  if (is.null(proj_root) || !file.exists(file.path(proj_root, "CLAUDE.md"))) {
    proj_root <- getwd()
  }
}

path_raw       <- file.path(proj_root, "data", "raw")
path_processed <- file.path(proj_root, "data", "processed")
path_figures   <- file.path(proj_root, "output", "figures")
path_tables    <- file.path(proj_root, "output", "tables")
path_maps      <- file.path(proj_root, "output", "maps")

# Ensure output directories exist
walk(c(path_processed, path_figures, path_tables, path_maps), dir.create,
     showWarnings = FALSE, recursive = TRUE)

# Raw data sub-paths
path_notifications <- file.path(path_raw, "notifications")
path_services      <- file.path(path_raw, "services")
path_demographics  <- file.path(path_raw, "demographics")
path_boundaries    <- file.path(path_raw, "boundaries")
path_congenital    <- file.path(path_raw, "congenital")

# ---- CRS --------------------------------------------------------------------
crs_gda2020 <- 7844  # EPSG:7844 GDA2020 geographic

# ---- Plotting defaults -------------------------------------------------------
theme_set(theme_minimal(base_size = 12))
tmap_mode("plot")

# Consistent colour palette for states
state_colours <- c(

  "NSW" = "#1b9e77", "VIC" = "#d95f02", "QLD" = "#7570b3",
  "SA"  = "#e7298a", "WA"  = "#66a61e", "TAS" = "#e6ab02",
  "NT"  = "#a6761d", "ACT" = "#666666"
)

# ---- Small-area suppression threshold ----------------------------------------
suppress_threshold <- 5  # Suppress counts < 5 in output tables

# ---- Helper functions --------------------------------------------------------

#' Suppress small counts for privacy
#' @param x numeric vector of counts
#' @param threshold minimum count to display (default 5)
#' @return character vector with small counts replaced by "<5"
suppress_small <- function(x, threshold = suppress_threshold) {
  ifelse(!is.na(x) & x > 0 & x < threshold, "<5", as.character(x))
}

#' Save a ggplot as both PDF and PNG
#' @param p ggplot object
#' @param filename base filename (no extension)
#' @param width width in inches
#' @param height height in inches
#' @param dpi resolution for PNG
save_plot <- function(p, filename, width = 8, height = 6, dpi = 300) {
  ggsave(file.path(path_figures, paste0(filename, ".pdf")),
         plot = p, width = width, height = height)
  ggsave(file.path(path_figures, paste0(filename, ".png")),
         plot = p, width = width, height = height, dpi = dpi)
}

#' Save a tmap object as both PDF and PNG
#' @param m tmap object
#' @param filename base filename (no extension)
#' @param width width in inches
#' @param height height in inches
#' @param dpi resolution for PNG
save_map <- function(m, filename, width = 8, height = 6, dpi = 300) {
  tmap_save(m, file.path(path_maps, paste0(filename, ".pdf")),
            width = width, height = height)
  tmap_save(m, file.path(path_maps, paste0(filename, ".png")),
            width = width, height = height, dpi = dpi)
}

message("Setup complete. Project root: ", proj_root)
