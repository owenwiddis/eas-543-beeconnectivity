# Kent County: three-species reported counts and local Moran's I, August 2023–2025
# Run from the repository root; inputs are in data/pilot_inputs/ (not committed).
# Year comparisons use identical 6 km hexagons, assignment, and queen-neighbor rules.
# Full-county LISA remains the verified baseline. A smaller FIXED 58-hex study
# area is supplied by an auditable ID mask, not regenerated using LISA p-values.
# Its selection (Developed >= 30%, largest connected component, one queen ring)
# was informed by observed reports, so it is exploratory, not an independent test.

if (!file.exists("data/pilot_inputs/kent_2024.geojson")) {
  stop("Set your RStudio working directory to the cloned GitHub repository root first.")
}

YEARS <- 2023:2025
MONTH <- 8L
HEX_WIDTH_M <- 6000   # sf::st_make_grid hex distance across opposite edges
NSIM <- 999L
FOCUS_N_HEXES <- 58L     # verified one-ring study area, held constant across years
DEVELOPED_THRESHOLD <- 30 # percent NLCD 2025 classes 21-24 used for AOI selection
# Focus identifiers and geometry are included in the restricted team data ZIP.
OUT_DIR <- "outputs/figures/kent_2023_2025"

pkgs <- c("sf", "ggplot2", "spdep", "terra")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Date parsing mirrors the original Kent pilot for YYYY/M/D or YYYY-MM-DD.
parse_observed_date <- function(x) {
  x <- substr(trimws(as.character(x)), 1L, 10L)
  as.Date(gsub("-", "/", x, fixed = TRUE), format = "%Y/%m/%d")
}

read_year <- function(year) {
  # This canonical 5-column Kent-only August CSV is all teammates need.
  # The original three 2025 source files stay untouched on the author's machine.
  filename <- sprintf("data/pilot_inputs/kent_bees_%d_%02d.csv", year, MONTH)
  if (!file.exists(filename)) {
    stop("Missing ", filename,
         "; run python/download_kent_bees_2023_2024.py or unzip the team data ZIP into data/.")
  }
  d <- read.csv(filename, stringsAsFactors = FALSE)
  required <- c("uuid", "observed_on", "longitude", "latitude", "species")
  if (!all(required %in% names(d))) stop("Missing required columns for ", year)
  d <- d[!is.na(d$uuid) & nzchar(as.character(d$uuid)), , drop = FALSE]
  d <- d[!duplicated(as.character(d$uuid)), , drop = FALSE]
  d$date <- parse_observed_date(d$observed_on)
  d$longitude <- suppressWarnings(as.numeric(d$longitude))
  d$latitude <- suppressWarnings(as.numeric(d$latitude))
  valid <- !is.na(d$date) & format(d$date, "%Y") == as.character(year) &
    format(d$date, "%m") == sprintf("%02d", MONTH) &
    is.finite(d$longitude) & is.finite(d$latitude) &
    abs(d$longitude) <= 180 & abs(d$latitude) <= 90
  d[valid, , drop = FALSE]
}

county_file <- "data/pilot_inputs/kent_2024.geojson"
if (!file.exists(county_file)) stop("Missing Kent boundary: ", county_file)
kent <- sf::st_read(county_file, quiet = TRUE)
if ("GEOID" %in% names(kent)) kent <- kent[as.character(kent$GEOID) == "26081", ]
if (nrow(kent) != 1L) stop("Expected exactly one Kent polygon in GeoJSON")
kent <- sf::st_transform(sf::st_make_valid(kent), 5070)
kent <- sf::st_sf(name = "Kent County", geometry = sf::st_union(sf::st_geometry(kent)))

# Build ONE grid and ONE neighborhood graph for all three years.
grid_geom <- sf::st_make_grid(kent, cellsize = HEX_WIDTH_M, square = FALSE)
hex_all <- sf::st_sf(cell_id = seq_along(grid_geom), geometry = grid_geom)
hex_all <- hex_all[lengths(sf::st_intersects(hex_all, kent)) > 0L, , drop = FALSE]
hex_all$cell_id <- seq_len(nrow(hex_all))
nb_full <- spdep::poly2nb(hex_all, queen = TRUE)
stopifnot(nrow(hex_all) > 20L)

# Focus selection must run AFTER counts are read. No county boundary or grid is
# ever overwritten: the small study-area polygon is a new derived GeoJSON.

# Counts are assigned to the full cells before clipping for display.
counts_by_year <- list()
summary_by_year <- list()
for (year in YEARS) {
  obs <- read_year(year)
  if (!nrow(obs)) warning("No source observations for ", year)
  if (nrow(obs)) {
    points <- sf::st_as_sf(obs, coords = c("longitude", "latitude"), crs = 4326)
    points <- sf::st_transform(points, 5070)
    inside <- lengths(sf::st_intersects(points, kent)) > 0L
    points <- points[inside, , drop = FALSE]
    matches <- sf::st_intersects(points, hex_all)
    assignment <- vapply(matches, function(idx) if (length(idx)) idx[1L] else NA_integer_, integer(1))
    if (anyNA(assignment)) stop("Some Kent observations are not covered by the grid for ", year)
    v <- tabulate(assignment, nbins = nrow(hex_all))
    per_species <- table(points$species)
  } else {
    v <- integer(nrow(hex_all))
    per_species <- integer(0)
  }
  n_county <- if (nrow(obs) > 0L) nrow(points) else 0L
  if (sum(v) != n_county) {
    stop("Grid totals did not match county observation count for ", year)
  }
  counts_by_year[[as.character(year)]] <- v
  summary_by_year[[as.character(year)]] <- data.frame(
    year = year, observed_reports = sum(v), hex_cells = length(v),
    zero_cells = sum(v == 0), max_cell_count = max(v),
    average_all_cells = mean(v),
    nonzero_cells = sum(v > 0)
  )
  cat(sprintf("%d: %d county records, %d hexes, %d zero-report hexes\n",
              year, sum(v), length(v), sum(v == 0)))
  print(per_species)
}
# Compare the 2025 counts to the already-tested single-year pilot, if present.
old_counts_file <- "outputs/figures/kent_2025_08/Kent_hex_6km_counts.csv"
if (file.exists(old_counts_file)) {
  old_counts <- read.csv(old_counts_file)
  if (!all(c("cell_id", "bee_count") %in% names(old_counts)) ||
      nrow(old_counts) != nrow(hex_all) ||
      !identical(as.integer(old_counts$cell_id), as.integer(hex_all$cell_id)) ||
      !identical(as.integer(old_counts$bee_count), as.integer(counts_by_year[["2025"]]))) {
    stop("2025 6-km grid does NOT match the validated Kent pilot. Check inputs/CRS/grid.")
  }
  cat("2025 per-cell grid counts match the validated single-year Kent pilot.\n")
} else {
  warning("Original 2025 counts CSV is unavailable for independent grid consistency check")
}

# Freeze the previously verified Kent 2025 baseline (do not silently alter it).
if (sum(counts_by_year[["2025"]]) != 208L || length(counts_by_year[["2025"]]) != 95L) {
  stop("2025 baseline mismatch (expected 208 reports / 95 hexes). Check snapshots and grid.")
}

# ---- FIXED, CONTIGUOUS 58-HEX URBAN LANDSCAPE STUDY AREA ------------------
# Selection was completed separately using 2025 NLCD Developed (classes 21-24):
# >=30% -> largest queen-connected component (33 cells) -> one-ring neighbors
# in the original 95-hex graph (58 cells). The choice was informed by observed
# bee reports, so this smaller-area LISA is EXPLORATORY, not confirmatory.
# A frozen cell-ID mask means all years and all collaborators use the SAME AOI,
# and no one needs the large NLCD raster to reproduce the LISA statistics.

focus_id_file <- "data/pilot_inputs/kent_urban_focus_58hex_ids.csv"
focus_geojson_file <- "data/pilot_inputs/kent_urban_focus_58hex.geojson"
if (!file.exists(focus_id_file) || !file.exists(focus_geojson_file)) {
  stop("Missing urban58 input(s): ", focus_id_file, " and/or ", focus_geojson_file,
       ". Download the restricted team ZIP and unzip into data/.")
}
focus_id_data <- utils::read.csv(focus_id_file)
if (!identical(names(focus_id_data), "cell_id")) {
  stop("Urban focus IDs CSV must have exactly one column: cell_id")
}
focus_ids <- suppressWarnings(as.integer(focus_id_data$cell_id))
if (anyNA(focus_ids) || length(focus_ids) != FOCUS_N_HEXES ||
    anyDuplicated(focus_ids) || !all(focus_ids %in% hex_all$cell_id)) {
  stop("Urban focus IDs must contain 58 unique cell IDs from the frozen 95-cell grid")
}
focus_idx <- sort(match(focus_ids, hex_all$cell_id))
focus_hex <- hex_all[focus_idx, , drop = FALSE]
focus_nb <- spdep::poly2nb(focus_hex, queen = TRUE)
if (spdep::n.comp.nb(focus_nb)$nc != 1L) {
  stop("Urban58 study area is not one connected queen-neighbor component")
}
# The shared polygons are provided for GIS, but the LISA geometry is rebuilt
# from the original grid and frozen cell IDs to prevent alignment differences.
saved_focus_hex <- sf::st_read(focus_geojson_file, quiet = TRUE)
if (nrow(saved_focus_hex) != FOCUS_N_HEXES ||
    !"cell_id" %in% names(saved_focus_hex) ||
    !setequal(as.integer(saved_focus_hex$cell_id), focus_ids)) {
  stop("Urban58 GeoJSON must have the same 58 cell_id values as the CSV mask")
}
focus_geom <- suppressWarnings(sf::st_intersection(
  sf::st_union(sf::st_geometry(focus_hex)), sf::st_geometry(kent)))
if (length(focus_geom) != 1L || any(sf::st_is_empty(focus_geom))) {
  stop("Could not construct the clipped urban58 study-area polygon")
}
focus_area <- sf::st_sf(
  area_name = "Kent urban core + one ring (exploratory)",
  selected_hexes = FOCUS_N_HEXES,
  developed_threshold_percent = DEVELOPED_THRESHOLD,
  geometry = focus_geom
)
# Create a ready-to-use, dissolved smaller county study-area polygon for GIS.
sf::st_write(sf::st_transform(focus_area, 4326),
             "data/pilot_inputs/kent_urban_focus_area.geojson",
             delete_dsn = TRUE, quiet = TRUE)
utils::write.csv(
  data.frame(cell_id = hex_all$cell_id,
             in_urban_focus = hex_all$cell_id %in% focus_ids),
  file.path(OUT_DIR, "urban_focus_58hex_mask.csv"), row.names = FALSE)
cat(sprintf("Urban-focused study area: %d/%d connected hexes; one frozen region for 2023-25.\n",
            length(focus_idx), nrow(hex_all)))
for (y in YEARS) {
  vals <- counts_by_year[[as.character(y)]]
  cat(sprintf("  Urban focus August %d: %d/%d reports (%0.1f%%), %d/%d cells with reports.\n",
              y, sum(vals[focus_idx]), sum(vals), 100 * sum(vals[focus_idx])/sum(vals),
              sum(vals[focus_idx] > 0L), length(focus_idx)))
}
# Delete stale graphic/summary output names from the earlier exploratory
# interior-core and shrub-only versions; leave data/pilot_inputs untouched.
stale_files <- c(
  unlist(lapply(YEARS, function(y) file.path(OUT_DIR,
    paste0("Kent_6km_interior_core_LISA_", y, c(".png", ".csv"))))),
  unlist(lapply(YEARS, function(y) file.path(OUT_DIR,
    paste0("Kent_6km_shrub_observed_focus_LISA_", y, c(".png", ".csv"))))),
  unlist(lapply(YEARS, function(y) file.path(OUT_DIR,
    paste0("Kent_6km_shrub_focus_observed_counts_", y, ".png")))),
  file.path(OUT_DIR, c("grid_cell_ids_and_core_mask.csv",
                       "kent_shrub_focus_cells.csv", "kent_shrub_focus.geojson",
                       "Kent_6km_shrub_fraction_and_focus.png"))
)
# These are derived outputs only, and can be regenerated from prior branches.
unlink(stale_files[file.exists(stale_files)])
# ---------------------------------------------------------------------------

summary_table <- do.call(rbind, summary_by_year)
row.names(summary_table) <- NULL
summary_table$urban_focus_reports <- vapply(YEARS, function(y)
  sum(counts_by_year[[as.character(y)]][focus_idx]), integer(1))
summary_table$urban_focus_nonzero_cells <- vapply(YEARS, function(y)
  sum(counts_by_year[[as.character(y)]][focus_idx] > 0L), integer(1))
utils::write.csv(summary_table, file.path(OUT_DIR, "year_summary.csv"), row.names = FALSE)

wide <- data.frame(cell_id = hex_all$cell_id,
                   in_urban_focus = hex_all$cell_id %in% focus_ids)
for (year in YEARS) wide[[paste0("count_", year)]] <- counts_by_year[[as.character(year)]]
wide$delta_2024_minus_2023 <- wide$count_2024 - wide$count_2023
wide$delta_2025_minus_2024 <- wide$count_2025 - wide$count_2024
utils::write.csv(wide, file.path(OUT_DIR, "per_hex_counts_and_changes.csv"), row.names = FALSE)

lisa_results <- list()
focus_lisa_results <- list()
calc_lisa <- function(grid, x, year, label) {
  nb <- spdep::poly2nb(grid, queen = TRUE)
  out <- grid
  out$bee_count <- as.integer(x)
  out$local_I <- NA_real_
  out$p_perm <- NA_real_
  out$p_BH <- NA_real_
  out$LISA <- "Not significant"
  if (length(x) >= 20L && length(unique(x)) > 1L && any(spdep::card(nb) > 0L)) {
    weights <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
    set.seed(2026)
    test <- spdep::localmoran_perm(x, weights, nsim = NSIM, zero.policy = TRUE,
                                  alternative = "two.sided", iseed = 2026)
    perm_col <- colnames(test)[grepl("Sim$", colnames(test)) &
                                   !grepl("folded", colnames(test), ignore.case = TRUE)]
    if (!length(perm_col)) stop("No simulation p-value in localmoran_perm output")
    out$local_I <- as.numeric(test[, "Ii"])
    out$p_perm <- as.numeric(test[, perm_col[1L]])
    out$p_BH <- p.adjust(out$p_perm, method = "BH")
    z <- x - mean(x)
    lag_z <- as.numeric(spdep::lag.listw(weights, z, zero.policy = TRUE))
    quadrant <- ifelse(z > 0 & lag_z > 0, "High-High",
                ifelse(z < 0 & lag_z < 0, "Low-Low",
                ifelse(z > 0 & lag_z < 0, "High-Low",
                ifelse(z < 0 & lag_z > 0, "Low-High", "Not significant"))))
    sig <- !is.na(out$p_BH) & out$p_BH < 0.05
    out$LISA[sig] <- quadrant[sig]
    out$LISA[spdep::card(nb) == 0L] <- "No neighbor"
  }
  out$year <- year
  out$analysis <- label
  out
}

plot_geometry <- function(data, study_area = kent) {
  suppressWarnings(sf::st_intersection(data, sf::st_geometry(study_area)))
}
colors <- c("High-High" = "#b3474b", "Low-Low" = "#4676a8",
            "High-Low" = "#d69a59", "Low-High" = "#7caf93",
            "Not significant" = "#e0e3e3", "No neighbor" = "#b3b9bd")

save_lisa <- function(result, year, label, study_area = kent) {
  stem <- paste0("Kent_6km_", label, "_", year)
  utils::write.csv(sf::st_drop_geometry(result), file.path(OUT_DIR, paste0(stem, ".csv")), row.names = FALSE)
  display <- plot_geometry(result, study_area)
  display$LISA <- factor(display$LISA, levels = names(colors))
  plot <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = display, ggplot2::aes(fill = LISA), color = "white", linewidth = 0.3) +
    ggplot2::geom_sf(data = study_area, fill = NA, color = "#263f4c", linewidth = 0.6) +
    ggplot2::scale_fill_manual(values = colors, drop = FALSE, name = "Local Moran") +
    ggplot2::labs(title = sprintf("Kent | %s | August %d | 6 km hexagons", label, year),
                  subtitle = "999 permutations | BH-adjusted p < 0.05 | 3 species",
                  caption = if (label == "urban_focus_58hex_LISA") {
                    "EXPLORATORY: Developed >=30%, largest contiguous core + one queen ring.\nStudy area was informed by observations; all zero-report cells are RETAINED."
                  } else {
                    "Full county including zero-report hexes. Observation records are NOT abundance."
                  }) +
    ggplot2::theme_void(base_size = 12)
  ggplot2::ggsave(file.path(OUT_DIR, paste0(stem, ".png")), plot = plot,
                  width = 9, height = 7.5, dpi = 240, bg = "white")
  cat(sprintf("%s %d: ", label, year))
  print(table(result$LISA))
}

# Common color scale is fixed across the three years for count-map comparison.
common_max <- max(unlist(counts_by_year))
for (year in YEARS) {
  key <- as.character(year)
  value <- counts_by_year[[key]]
  full <- hex_all
  full$bee_count <- value
  display_counts <- plot_geometry(full)
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = display_counts, ggplot2::aes(fill = bee_count), color = "white", linewidth = 0.25) +
    ggplot2::geom_sf(data = kent, fill = NA, color = "#263f4c", linewidth = 0.6) +
    ggplot2::scale_fill_gradient(low = "#ecf2eb", high = "#24734f", limits = c(0, max(1, common_max)),
                                 name = "Reports") +
    ggplot2::labs(title = sprintf("Kent bee reports | August %d | 6 km hexagons", year),
                  subtitle = sprintf("Three study species | N=%d | identical grid across years", sum(value)),
                  caption = "No report is not confirmed absence. Same color scale for 2023–2025.") +
    ggplot2::theme_void(base_size = 12)
  ggplot2::ggsave(file.path(OUT_DIR, sprintf("Kent_6km_counts_%d.png", year)),
                  plot = p, width = 9, height = 7.5, dpi = 240, bg = "white")

  lisa <- calc_lisa(hex_all, value, year, "full_county")
  save_lisa(lisa, year, "full_county_LISA")
  lisa_results[[key]] <- lisa
  # Same fixed 58 hexes for each year. Zero-report cells are INCLUDED.
  focus_value <- value[focus_idx]
  focus_lisa <- calc_lisa(focus_hex, focus_value, year,
                          "urban_focus_58hex_exploratory")
  save_lisa(focus_lisa, year, "urban_focus_58hex_LISA", focus_area)
  focus_lisa_results[[key]] <- focus_lisa
  focus_display <- focus_hex
  focus_display$bee_count <- as.integer(focus_value)
  p_focus <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = focus_display, ggplot2::aes(fill = bee_count),
                     color = "white", linewidth = 0.3) +
    ggplot2::geom_sf(data = focus_area, fill = NA, color = "#263f4c",
                     linewidth = 0.7) +
    ggplot2::scale_fill_gradient(low = "#ecf2eb", high = "#24734f",
                                 limits = c(0, max(1, common_max)), name = "Reports") +
    ggplot2::labs(title = sprintf("Kent urban focus | August %d | 6 km hexagons", year),
                  subtitle = sprintf("%d reports in 58 connected cells | fixed area across all years",
                                     sum(focus_value)),
                  caption = "EXPLORATORY: 2025 NLCD developed >=30%, largest core + one ring.\nZero-report cells retained; reports are not confirmed bee abundance.") +
    ggplot2::theme_void(base_size = 12)
  ggplot2::ggsave(file.path(OUT_DIR,
                  sprintf("Kent_6km_urban_focus_counts_%d.png", year)),
                  p_focus, width = 9, height = 7.5, dpi = 240, bg = "white")
}

# Report overlap of the fixed-cell HH patterns, not individual bee movements.
hh <- data.frame(cell_id = hex_all$cell_id)
for (year in YEARS) hh[[paste0("HH_", year)]] <- lisa_results[[as.character(year)]]$LISA == "High-High"
utils::write.csv(hh, file.path(OUT_DIR, "high_high_by_year.csv"), row.names = FALSE)
hh_focus <- data.frame(cell_id = focus_hex$cell_id)
for (year in YEARS) hh_focus[[paste0("HH_", year)]] <-
  focus_lisa_results[[as.character(year)]]$LISA == "High-High"
utils::write.csv(hh_focus, file.path(OUT_DIR, "urban_focus_high_high_by_year.csv"), row.names = FALSE)

cat("\nDONE. Outputs saved to: ", normalizePath(OUT_DIR), "\n", sep = "")
cat("PRIMARY: full_county_LISA (95 cells); EXPLORATORY: urban_focus_58hex_LISA (58 cells).\n")
cat("Frozen Urban58 cell IDs and GeoJSON are in data/pilot_inputs/ for restricted team sharing.\n")
cat("Data-selected focus is exploratory, NOT a p-value improvement strategy.\n")
cat("Do not interpret changes in reported counts as bee migration or abundance changes.\n")
