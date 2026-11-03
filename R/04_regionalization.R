# PCA on landscape metrics, spatially constrained regionalization,
# then compare bee observation rates across regions.

source("R/00_setup.R")
library(sf); library(dplyr); library(rgeoda)

metrics <- st_read("data/processed/grid_metrics.gpkg", quiet = TRUE)

vars <- c("pland", "cohesion", "mesh_ha", "dist_mean_m", "iic_1km", "iic_3km")
X <- metrics |>
  st_drop_geometry() |>
  select(all_of(vars)) |>
  mutate(dist_mean_m = -dist_mean_m) |>      # flip so higher = more connected
  mutate(across(everything(), ~ ifelse(is.na(.x), 0, .x))) |>  # TODO: justify NA handling
  scale() |>
  as.data.frame()

# --- PCA ---------------------------------------------------------------------
pca <- prcomp(X)
summary(pca)
print(round(pca$rotation[, 1:3], 2))   # loadings: which metrics drive each gradient
metrics$PC1 <- pca$x[, 1]

# --- Regionalization ---------------------------------------------------------
w <- queen_weights(metrics)
# TODO: choose k with the indexes in Chen et al. (2024) / Mendler et al. (2025)
k <- 5
reg <- skater(k, w, X)
metrics$region <- factor(reg$Clusters)
print(reg$`The ratio of between to total sum of squares`)

# --- Compare bee observation rates across regions -----------------------------
bee_cols <- intersect(c("Apis mellifera", "Bombus impatiens", "Xylocopa virginica"),
                      names(metrics))
metrics$bees_total <- rowSums(st_drop_geometry(metrics)[, bee_cols, drop = FALSE])
metrics$bee_rate   <- metrics$bees_total / metrics$area_km2

kruskal.test(bee_rate ~ region, data = metrics)

st_write(metrics, "data/processed/grid_regions.gpkg", delete_dsn = TRUE)
