# Install (once) and load the packages used in this project.

cran_pkgs <- c("sf", "terra", "dplyr", "tidyr", "landscapemetrics",
               "lconnect", "spdep", "rgeoda", "ggplot2", "tmap", "remotes")
missing <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(missing)) install.packages(missing)

# Not on CRAN
if (!requireNamespace("Makurhini", quietly = TRUE)) {
  remotes::install_github("connectscape/Makurhini", dependencies = TRUE)
}
if (!requireNamespace("urbioconnect", quietly = TRUE)) {
  install.packages("urbioconnect",
                   repos = c("https://urbio-ecology.r-universe.dev", "https://cloud.r-project.org"))
}

# Shared settings used by every script
NLCD_PATH      <- "data/raw/Annual_NLCD_LndCov_2025_CU_C1V2.tif"
PILOT_COUNTY   <- "Kent"
PILOT_FIPS     <- "26081"
PILOT_MONTH    <- 8
HABITAT_CODES  <- c(52, 71, 81)   # shrub/scrub, grassland, pasture/hay
CELL_SIZE      <- c(5000, 3000)   # grid cell width x height in meters
BUFFER_M       <- 3000            # landscape buffer around the county
DIST_THRESH    <- c(1000, 3000)   # foraging-distance thresholds in meters
