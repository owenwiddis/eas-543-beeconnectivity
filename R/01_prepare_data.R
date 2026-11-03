# Clean observations, select the pilot county, build the 3 x 5 km grid,
# and count observations per grid cell.

source("R/00_setup.R")
library(sf); library(terra); library(dplyr); library(tidyr)

nlcd_crs <- crs(rast(NLCD_PATH))   # Albers Equal Area, WGS84 datum

# --- Observations ----------------------------------------------------------
obs <- read.csv("data/raw/inat_bees_MI_2025.csv") |>
  filter(usable_point) |>                       # not obscured, accuracy <= 1 km
  distinct(id, .keep_all = TRUE) |>
  mutate(observed_on = as.Date(observed_on),
         month = as.integer(format(observed_on, "%m")),
         species = sub("^(\\S+ \\S+).*", "\\1", taxon_name)) |>  # subspecies -> species
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(nlcd_crs)

# --- Counties --------------------------------------------------------------
counties <- st_read("data/raw/tl_2025_us_county/tl_2025_us_county.shp", quiet = TRUE) |>
  filter(STATEFP == "26") |>
  st_transform(nlcd_crs)

pilot <- filter(counties, GEOID == PILOT_FIPS)

# --- Pilot observations (August, inside the county) ------------------------
obs_pilot <- obs |>
  filter(month == PILOT_MONTH) |>
  st_filter(pilot)

# --- Grid ------------------------------------------------------------------
grid <- st_make_grid(pilot, cellsize = CELL_SIZE, square = TRUE) |>
  st_as_sf() |>
  st_intersection(st_geometry(pilot)) |>
  st_collection_extract("POLYGON") |>
  mutate(cell_id = row_number(),
         area_km2 = as.numeric(st_area(geometry)) / 1e6)
# TODO: decide what to do with small edge slivers (e.g. drop cells < 5 km2)

# --- Counts per cell and species --------------------------------------------
counts <- st_join(obs_pilot, grid["cell_id"]) |>
  st_drop_geometry() |>
  count(cell_id, species) |>
  pivot_wider(names_from = species, values_from = n, values_fill = 0)

grid <- grid |>
  left_join(counts, by = "cell_id") |>
  mutate(across(-c(cell_id, area_km2, geometry), ~ replace_na(.x, 0)))

st_write(obs_pilot, "data/processed/obs_pilot.gpkg", delete_dsn = TRUE)
st_write(grid,      "data/processed/grid_pilot.gpkg", delete_dsn = TRUE)
st_write(pilot,     "data/processed/pilot_county.gpkg", delete_dsn = TRUE)
