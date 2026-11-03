# Clip NLCD to the pilot county (+ buffer), build the binary habitat raster,
# burn in road barriers, and delineate habitat patches.

source("R/00_setup.R")
library(sf); library(terra); library(dplyr)

pilot <- st_read("data/processed/pilot_county.gpkg", quiet = TRUE)
nlcd  <- rast(NLCD_PATH)

# --- Clip --------------------------------------------------------------------
extent_v <- vect(st_buffer(pilot, BUFFER_M))
nlcd_p   <- crop(nlcd, extent_v, mask = TRUE)

# --- Habitat raster: 1 = habitat, 0 = not --------------------------------------
habitat <- ifel(nlcd_p %in% HABITAT_CODES, 1, 0)
names(habitat) <- "habitat"

# --- Roads as barriers (for effective mesh size) -----------------------------
roads <- st_read("data/raw/tl_2025_26081_roads/tl_2025_26081_roads.shp", quiet = TRUE) |>
  st_transform(crs(habitat))
road_r <- rasterize(vect(roads), habitat, field = 1, background = 0, touches = TRUE)
habitat_barrier <- ifel(road_r == 1, 0, habitat)
# TODO: decide whether to burn in all roads or only primary/secondary (MTFCC S1100, S1200)

# --- Habitat patches (8-neighbour) -------------------------------------------
patch_r <- patches(ifel(habitat == 1, 1, NA), directions = 8)
patches_sf <- as.polygons(patch_r) |>
  st_as_sf() |>
  rename(patch_id = 1) |>
  mutate(area_m2 = as.numeric(st_area(geometry)),
         class = 1)

writeRaster(habitat,         "data/processed/habitat.tif", overwrite = TRUE)
writeRaster(habitat_barrier, "data/processed/habitat_barrier.tif", overwrite = TRUE)
st_write(patches_sf,         "data/processed/habitat_patches.gpkg", delete_dsn = TRUE)
