# Connectivity and fragmentation metrics per grid cell.
#   % habitat (amount)         landscapemetrics  lsm_c_pland
#   patch cohesion (structural) landscapemetrics lsm_c_cohesion
#   effective mesh size         landscapemetrics lsm_c_mesh (roads burned in)
#   distance to habitat         terra distance + zonal mean
#   IIC (graph, binary)         lconnect
#   PC  (graph, probabilistic)  Makurhini

source("R/00_setup.R")
library(sf); library(terra); library(dplyr); library(tidyr); library(landscapemetrics)

grid            <- st_read("data/processed/grid_pilot.gpkg", quiet = TRUE)
habitat         <- rast("data/processed/habitat.tif")
habitat_barrier <- rast("data/processed/habitat_barrier.tif")
patches_sf      <- st_read("data/processed/habitat_patches.gpkg", quiet = TRUE)

# --- 1. Amount and structure (class 1 = habitat) ------------------------------
lsm <- sample_lsm(habitat, y = grid, plot_id = grid$cell_id,
                  what = c("lsm_c_pland", "lsm_c_cohesion")) |>
  filter(class == 1) |>
  select(cell_id = plot_id, metric, value) |>
  pivot_wider(names_from = metric, values_from = value)

# --- 2. Fragmentation: effective mesh size with roads as barriers ------------
mesh <- sample_lsm(habitat_barrier, y = grid, plot_id = grid$cell_id,
                   what = "lsm_c_mesh") |>
  filter(class == 1) |>
  select(cell_id = plot_id, mesh_ha = value)
# Note: landscapemetrics uses the cutting-out procedure at cell edges.

# --- 3. Distance to nearest habitat (simple functional) ----------------------
dist_r <- distance(ifel(habitat == 1, 1, NA))          # meters to nearest habitat pixel
dist_cell <- terra::extract(dist_r, vect(grid), fun = mean, na.rm = TRUE)
dist_cell <- data.frame(cell_id = grid$cell_id, dist_mean_m = dist_cell[, 2])
# HCI alternative: compute at regular sample points per cell, NOT at bee locations,
# so the landscape measure stays independent of the observations.

# --- 4. Graph-based connectivity per cell ------------------------------------
iic_for_cell <- function(cell, max_dist) {
  p <- st_intersection(patches_sf, st_geometry(cell))
  if (nrow(p) < 2) return(NA_real_)
  tmp <- tempfile(fileext = ".shp")
  st_write(p["class"], tmp, quiet = TRUE)
  land <- lconnect::upload_land(tmp, habitat = 1, max_dist = max_dist)
  # TODO: confirm which attribute field upload_land reads as the class
  as.numeric(lconnect::con_metric(land, metric = "IIC"))
}

graph <- grid |>
  st_drop_geometry() |>
  select(cell_id) |>
  mutate(
    iic_1km = sapply(seq_len(nrow(grid)), \(i) iic_for_cell(grid[i, ], DIST_THRESH[1])),
    iic_3km = sapply(seq_len(nrow(grid)), \(i) iic_for_cell(grid[i, ], DIST_THRESH[2]))
  )

# PC with Makurhini (per-cell landscape)
# TODO: wrap like iic_for_cell; see ?Makurhini::MK_dPCIIC
# Makurhini::MK_dPCIIC(nodes = p, attribute = "area_m2", distance = list(type = "edge"),
#                      metric = "PC", probability = 0.5,
#                      distance_thresholds = DIST_THRESH)

# --- Combine -----------------------------------------------------------------
metrics <- grid |>
  left_join(lsm,       by = "cell_id") |>
  left_join(mesh,      by = "cell_id") |>
  left_join(dist_cell, by = "cell_id") |>
  left_join(graph,     by = "cell_id")

st_write(metrics, "data/processed/grid_metrics.gpkg", delete_dsn = TRUE)

# Quick agreement check between metrics (Spearman)
metrics |>
  st_drop_geometry() |>
  select(pland, cohesion, mesh_ha, dist_mean_m, starts_with("iic")) |>
  cor(method = "spearman", use = "pairwise.complete.obs") |>
  round(2) |>
  print()
