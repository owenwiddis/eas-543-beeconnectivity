# Kent County: August 2025, 3 bee species -> Bee Count grids + exploratory LISA
# ------------------------------------------------------------------------------
# Run from the repository root: source("R/pilots/Kent_BeeCount_Grid_LISA.R")
# Input CSVs and county shapefile belong in data/raw/ (ignored by Git).
# Outputs are saved under outputs/figures/kent_2025_08/.
# IMPORTANT: Counts represent iNaturalist reports, NOT biological abundance.

YEAR <- 2025L
MONTH <- 8L
MAIN_HEX_M <- 6000   # Hexagon distance between opposite sides (metres)
HEX_SENSITIVITY_M <- c(4000, 6000, 8000)
NSIM <- 999L

needed <- c("sf", "ggplot2", "spdep")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

out_dir <- "outputs/figures/kent_2025_08"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

species_files <- c(
  "Apis mellifera" = "data/pilot_inputs/Apis_mellifera.csv",
  "Bombus impatiens" = "data/pilot_inputs/Bombus_impatiens.csv",
  "Xylocopa virginica" = "data/pilot_inputs/Xylocopa_virginica.csv"
)

parse_observed_date <- function(x) {
  x <- trimws(as.character(x))
  x <- substr(x, 1L, 10L)
  x <- gsub("-", "/", x, fixed = TRUE)
  as.Date(x, format = "%Y/%m/%d")
}

read_species <- function(species_name, file_path) {
  if (!file.exists(file_path)) {
    stop("Missing ", file_path, ". Please put this CSV in the Bee project folder.")
  }
  x <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("uuid", "observed_on", "longitude", "latitude")
  if (!all(required %in% names(x))) {
    stop("CSV ", file_path, " is missing columns: ",
         paste(setdiff(required, names(x)), collapse = ", "))
  }
  county_text <- if ("place_county_name" %in% names(x)) {
    as.character(x$place_county_name)
  } else {
    rep(NA_character_, nrow(x))
  }
  data.frame(
    uuid = as.character(x$uuid),
    date = parse_observed_date(x$observed_on),
    longitude = suppressWarnings(as.numeric(x$longitude)),
    latitude = suppressWarnings(as.numeric(x$latitude)),
    species = species_name,
    county_label = county_text,
    stringsAsFactors = FALSE
  )
}

all_bees <- do.call(
  rbind,
  Map(read_species, names(species_files), unname(species_files))
)
all_bees <- all_bees[!is.na(all_bees$uuid) & nzchar(all_bees$uuid), ]
all_bees <- all_bees[!duplicated(all_bees$uuid), ]

selected <- all_bees[
  !is.na(all_bees$date) &
    as.integer(format(all_bees$date, "%Y")) == YEAR &
    as.integer(format(all_bees$date, "%m")) == MONTH &
    is.finite(all_bees$longitude) & is.finite(all_bees$latitude) &
    abs(all_bees$latitude) <= 90 & abs(all_bees$longitude) <= 180,
]
if (nrow(selected) == 0L) {
  stop("No valid records found for selected year/month.")
}

cat("\n--- Statewide CSV county labels for selected month (top 10) ---\n")
print(head(sort(table(selected$county_label), decreasing = TRUE), 10L))

# Use the boundary previously downloaded in the Kent NLCD pilot.
county_file <- "data/pilot_inputs/kent_2024.geojson"
if (!file.exists(county_file)) {
  stop("County boundary file not found: ", county_file,
       ". Copy the county_boundary_2024 folder from the Bee project into data/raw/.")
}
county_us <- sf::st_read(county_file, quiet = TRUE)
if (!"GEOID" %in% names(county_us)) {
  stop("County boundary is missing GEOID.")
}
kent <- county_us[as.character(county_us$GEOID) == "26081", ]
stopifnot(nrow(kent) == 1L)
# EPSG:5070 uses metres and is suitable for comparable-area grids.
kent <- sf::st_transform(sf::st_make_valid(kent), 5070)
kent <- sf::st_sf(name = "Kent County", geometry = sf::st_union(sf::st_geometry(kent)))

bees_sf <- sf::st_as_sf(
  selected, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE
)
bees_sf <- sf::st_transform(bees_sf, 5070)
is_in_county <- lengths(sf::st_intersects(bees_sf, kent)) > 0L
bees_kent <- bees_sf[is_in_county, ]
if (nrow(bees_kent) == 0L) stop("No selected points fall within Kent County geometry.")

cat("\n--- Kent County: verified points inside the county boundary ---\n")
cat("Month:", sprintf("%04d-%02d", YEAR, MONTH), "\n")
cat("All three species combined:", nrow(bees_kent), "\n")
print(table(bees_kent$species))
cat("Labels 'Kent' but outside boundary:",
    sum(trimws(selected$county_label) == "Kent", na.rm = TRUE) -
      sum(trimws(selected$county_label[is_in_county]) == "Kent", na.rm = TRUE), "\n")

# Construct fixed grid, retain all cells overlapping the county,
# and count EACH report exactly once by a unique cell assignment.
count_on_grid <- function(bees, county, grid_geom) {
  grid <- sf::st_sf(cell_id = seq_along(grid_geom), geometry = grid_geom)
  keep <- lengths(sf::st_intersects(grid, county)) > 0L
  grid <- grid[keep, , drop = FALSE]
  grid$cell_id <- seq_len(nrow(grid))

  matching_cells <- sf::st_intersects(bees, grid)
  assignment <- vapply(matching_cells, function(v) {
    if (length(v) == 0L) NA_integer_ else as.integer(v[1L])
  }, integer(1))
  if (anyNA(assignment)) {
    stop(sum(is.na(assignment)), " Kent bee point(s) not covered by the grid.")
  }
  for (nm in names(species_files)) {
    field <- gsub(" ", "_", nm, fixed = TRUE)
    grid[[field]] <- tabulate(
      assignment[as.character(bees$species) == nm], nbins = nrow(grid)
    )
  }
  fields <- gsub(" ", "_", names(species_files), fixed = TRUE)
  grid$bee_count <- rowSums(sf::st_drop_geometry(grid)[, fields, drop = FALSE])
  stopifnot(sum(grid$bee_count) == nrow(bees))

  # Use clipped shapes for the display, but original full cells for adjacency.
  display_grid <- suppressWarnings(sf::st_intersection(grid, sf::st_geometry(county)))
  list(grid = grid, display_grid = display_grid)
}

# Common map style. Note: part of an edge cell may be outside Kent.
plot_count <- function(result, title_text, subtitle_text, with_points = FALSE) {
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = result$display_grid,
      ggplot2::aes(fill = bee_count),
      color = "white", linewidth = 0.28
    ) +
    ggplot2::geom_sf(data = kent, fill = NA, color = "#263f4c", linewidth = 0.65) +
    ggplot2::scale_fill_gradient(
      low = "#ecf2eb", high = "#24734f", name = "Reported\nrecords"
    ) +
    ggplot2::labs(
      title = title_text, subtitle = subtitle_text,
      caption = "Three species combined | iNaturalist reports, not abundance | zero cells retained"
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(size = 10),
      legend.position = "right", plot.caption = ggplot2::element_text(size = 9)
    )
  if (with_points) {
    p <- p + ggplot2::geom_sf(
      data = bees_kent, color = "#c75a38", size = 0.40, alpha = 0.45
    )
  }
  p
}

# 1) Simple partitions: 2, 4, and 9 windows. These are for visualization,
# NOT meaningful standalone LISA tests.
partitions <- list(`2_cells` = c(2L, 1L), `4_cells` = c(2L, 2L),
                   `9_cells` = c(3L, 3L))
for (nm in names(partitions)) {
  gr <- sf::st_make_grid(kent, n = partitions[[nm]], square = TRUE)
  result <- count_on_grid(bees_kent, kent, gr)
  p <- plot_count(
    result,
    paste0("Kent: ", gsub("_", " ", nm), " (illustrative partition)"),
    sprintf("%04d-%02d | Three bee species", YEAR, MONTH)
  )
  ggplot2::ggsave(
    file.path(out_dir, paste0("Kent_", nm, "_bee_count.png")),
    plot = p, width = 8.2, height = 6.2, dpi = 240, bg = "white"
  )
  utils::write.csv(
    sf::st_drop_geometry(result$grid),
    file.path(out_dir, paste0("Kent_", nm, "_counts.csv")), row.names = FALSE
  )
}

# 2) Comparable fixed-area hexagons: 4, 6, 8 km.
scale_results <- list()
scale_summary <- list()
for (metres in sort(unique(c(HEX_SENSITIVITY_M, MAIN_HEX_M)))) {
  gr <- sf::st_make_grid(kent, cellsize = metres, square = FALSE)
  result <- count_on_grid(bees_kent, kent, gr)
  key <- paste0(metres / 1000, "km")
  scale_results[[key]] <- result
  values <- result$grid$bee_count
  scale_summary[[key]] <- data.frame(
    hex_width_m = metres,
    n_cells = length(values),
    total_count = sum(values),
    average_all_cells = mean(values),
    median_all_cells = stats::median(values),
    median_nonzero_cells = if (any(values > 0)) stats::median(values[values > 0]) else NA_real_,
    zero_cells = sum(values == 0),
    cells_with_at_least_10 = sum(values >= 10),
    max_cell_count = max(values)
  )
  p <- plot_count(
    result,
    sprintf("Kent bee records | %d km hexagons", metres / 1000),
    sprintf("%04d-%02d | Total N = %d | Grid cells = %d", YEAR, MONTH,
            sum(values), length(values)),
    with_points = metres == MAIN_HEX_M
  )
  ggplot2::ggsave(
    file.path(out_dir, sprintf("Kent_hex_%dkm_counts.png", metres / 1000)),
    plot = p, width = 9, height = 7.5, dpi = 260, bg = "white"
  )
  utils::write.csv(
    sf::st_drop_geometry(result$grid),
    file.path(out_dir, sprintf("Kent_hex_%dkm_counts.csv", metres / 1000)),
    row.names = FALSE
  )
}
scale_summary <- do.call(rbind, scale_summary)
row.names(scale_summary) <- NULL
utils::write.csv(
  scale_summary, file.path(out_dir, "Kent_hex_scale_diagnostics.csv"),
  row.names = FALSE
)
cat("\n--- Fixed hexagon scale diagnostics ---\n")
print(scale_summary, row.names = FALSE)

# 3) Exploratory LISA on the PREDEFINED 6 km hexagons.
# Queen-contiguity adjacency of fixed hex cells (full cell geometry).
main <- scale_results[[paste0(MAIN_HEX_M / 1000, "km")]]
nb <- spdep::poly2nb(main$grid, queen = TRUE)
w <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
x <- as.numeric(main$grid$bee_count)
main$grid$local_I <- NA_real_
main$grid$p_perm <- NA_real_
main$grid$p_BH <- NA_real_
main$grid$LISA <- "Not significant"

if (length(x) >= 20L && length(unique(x)) > 1L && any(spdep::card(nb) > 0L)) {
  set.seed(2026)
  lm <- spdep::localmoran_perm(
    x, w, nsim = NSIM, zero.policy = TRUE, alternative = "two.sided",
    iseed = 2026
  )
  main$grid$local_I <- as.numeric(lm[, "Ii"])
  # The simulation p-value has a version-dependent column title.
  sim_cols <- colnames(lm)[grepl("Sim$", colnames(lm)) &
                             !grepl("folded", colnames(lm), ignore.case = TRUE)]
  if (length(sim_cols) < 1L) stop("Cannot find permutation p-value in localmoran_perm output")
  main$grid$p_perm <- as.numeric(lm[, sim_cols[1L]])
  main$grid$p_BH <- stats::p.adjust(main$grid$p_perm, method = "BH")
  z <- x - mean(x)
  lag_z <- as.numeric(spdep::lag.listw(w, z, zero.policy = TRUE))
  quadrant <- ifelse(z > 0 & lag_z > 0, "High-High",
              ifelse(z < 0 & lag_z < 0, "Low-Low",
              ifelse(z > 0 & lag_z < 0, "High-Low",
              ifelse(z < 0 & lag_z > 0, "Low-High", "Not significant"))))
  significant <- !is.na(main$grid$p_BH) & main$grid$p_BH < 0.05
  main$grid$LISA[significant] <- quadrant[significant]
  main$grid$LISA[spdep::card(nb) == 0] <- "No neighbor"
} else {
  warning("Insufficient variation or too few cells for an exploratory LISA test")
}

# Replace clipped display geometry with current LISA/count attributes.
main$display_grid <- suppressWarnings(sf::st_intersection(main$grid, sf::st_geometry(kent)))
main$display_grid$LISA <- factor(
  main$display_grid$LISA,
  levels = c("High-High", "Low-Low", "High-Low", "Low-High", "Not significant", "No neighbor")
)

p_lisa <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = main$display_grid, ggplot2::aes(fill = LISA),
    color = "white", linewidth = 0.28
  ) +
  ggplot2::geom_sf(data = kent, fill = NA, color = "#263f4c", linewidth = 0.65) +
  ggplot2::scale_fill_manual(
    values = c("High-High" = "#b3474b", "Low-Low" = "#4676a8",
               "High-Low" = "#d69a59", "Low-High" = "#7caf93",
               "Not significant" = "#e0e3e3", "No neighbor" = "#b3b9bd"),
    drop = FALSE, name = "Local Moran"
  ) +
  ggplot2::labs(
    title = sprintf("Kent: exploratory LISA | %d km hexagons", MAIN_HEX_M / 1000),
    subtitle = sprintf("Three-species reported counts | %04d-%02d | 999 permutations | BH-adjusted p < 0.05", YEAR, MONTH),
    caption = paste(
      "No forced minimum records per cell; zeros mean no report, not confirmed absence.",
      "Boundary hexagons are partially inside the county; interpret with caution.", sep = "\n"
    )
  ) +
  ggplot2::theme_void(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 15),
    plot.subtitle = ggplot2::element_text(size = 10),
    plot.caption = ggplot2::element_text(size = 9)
  )
ggplot2::ggsave(
  file.path(out_dir, sprintf("Kent_hex_%dkm_LISA.png", MAIN_HEX_M / 1000)),
  plot = p_lisa, width = 9, height = 7.5, dpi = 260, bg = "white"
)

utils::write.csv(
  sf::st_drop_geometry(main$grid),
  file.path(out_dir, sprintf("Kent_hex_%dkm_LISA.csv", MAIN_HEX_M / 1000)),
  row.names = FALSE
)
sf::st_write(
  main$display_grid,
  dsn = file.path(out_dir, sprintf("Kent_hex_%dkm_LISA.gpkg", MAIN_HEX_M / 1000)),
  quiet = TRUE, delete_dsn = TRUE
)
cat("\n--- LISA groups (BH-adjusted, exploratory) ---\n")
print(table(main$grid$LISA, useNA = "ifany"))
cat("\nALL DONE. Output folder: ", normalizePath(out_dir), "\n", sep = "")
cat("Main maps: Kent_hex_6km_counts.png and Kent_hex_6km_LISA.png\n")
