# Kent County: Urban58 pooled August 2023-2025 exploratory Local Moran's I.
# Reproducible standalone follow-up to Kent_BeeCount_LISA_2023_2025.R.
# Run from repository root AFTER running the multiyear script.
# The focus area was selected after inspecting observations. This is NOT an
# independent confirmatory hotspot test. Never present unadjusted p-values as
# surviving multiple-comparison correction.

pkgs <- c("sf", "ggplot2", "spdep")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install missing R packages first: ", paste(missing, collapse = ", "))
}

OUT_DIR <- "outputs/figures/kent_2023_2025"
hex_file <- "data/pilot_inputs/kent_urban_focus_58hex.geojson"
count_file <- file.path(OUT_DIR, "per_hex_counts_and_changes.csv")
if (!file.exists(hex_file) || !file.exists(count_file)) {
  stop("Missing Urban58 GeoJSON or per-hex counts CSV. Run the multiyear script first.")
}

hex58 <- sf::st_read(hex_file, quiet = TRUE)
counts <- utils::read.csv(count_file)
required_cols <- c("cell_id", "count_2023", "count_2024", "count_2025")
if (!all(required_cols %in% names(counts)) || !"cell_id" %in% names(hex58)) {
  stop("Missing cell_id or yearly count columns")
}
idx <- match(hex58$cell_id, counts$cell_id)
stopifnot(nrow(hex58) == 58L, !anyNA(idx), !anyDuplicated(hex58$cell_id))

x <- rowSums(counts[idx, c("count_2023", "count_2024", "count_2025")])
stopifnot(sum(x) == 412L, all(is.finite(x)))

nb <- spdep::poly2nb(sf::st_transform(hex58, 5070), queen = TRUE)
stopifnot(spdep::n.comp.nb(nb)$nc == 1L, all(spdep::card(nb) > 0L))
lw <- spdep::nb2listw(nb, style = "W")

# Keep both unadjusted and BH-adjusted permutation p-values.
set.seed(543)
test <- spdep::localmoran_perm(
  x, lw, nsim = 999, alternative = "two.sided", iseed = 543
)
p_cols <- colnames(test)[grepl("Sim$", colnames(test)) &
                             !grepl("folded", colnames(test), ignore.case = TRUE)]
if (length(p_cols) != 1L) stop("Cannot uniquely find permutation p-value column")
p_perm <- as.numeric(test[, p_cols])
p_BH <- stats::p.adjust(p_perm, method = "BH")

z <- x - mean(x)
lag_z <- as.numeric(spdep::lag.listw(lw, z))
quad <- ifelse(z > 0 & lag_z > 0, "High-High",
        ifelse(z < 0 & lag_z < 0, "Low-Low",
        ifelse(z > 0 & lag_z < 0, "High-Low",
        ifelse(z < 0 & lag_z > 0, "Low-High", "Not significant"))))

raw_LISA <- rep("Not significant", length(x))
raw_LISA[!is.na(p_perm) & p_perm < 0.05] <- quad[!is.na(p_perm) & p_perm < 0.05]
BH_LISA <- rep("Not significant", length(x))
BH_LISA[!is.na(p_BH) & p_BH < 0.05] <- quad[!is.na(p_BH) & p_BH < 0.05]

cat("Urban58 pooled August reports:", sum(x), "\n")
cat("Unadjusted Local Moran types:\n")
print(table(raw_LISA))
cat("BH-significant cells:", sum(p_BH < 0.05, na.rm = TRUE), "\n")

# CSV contains BOTH classifications to prevent confusion.
result <- data.frame(
  cell_id = hex58$cell_id,
  bee_count = x,
  local_I = as.numeric(test[, "Ii"]),
  p_perm = p_perm,
  p_BH = p_BH,
  LISA_unadjusted = raw_LISA,
  LISA_BH = BH_LISA
)
utils::write.csv(result,
  file.path(OUT_DIR, "Kent_Urban58_LISA_Pooled_2023_2025.csv"),
  row.names = FALSE
)

hex58$pooled_count <- x
hex58$raw_LISA <- factor(
  raw_LISA,
  levels = c("High-High", "Low-Low", "High-Low", "Low-High", "Not significant")
)
sig_hex <- hex58[!is.na(hex58$raw_LISA) &
                   hex58$raw_LISA != "Not significant", , drop = FALSE]

# Retain exactly the well-received count + outline style; shorten the caption.
p <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = hex58, ggplot2::aes(fill = pooled_count),
    color = "grey75", linewidth = 0.35
  ) +
  ggplot2::geom_sf(
    data = sig_hex, ggplot2::aes(color = raw_LISA),
    fill = NA, linewidth = 1.3, show.legend = TRUE
  ) +
  ggplot2::scale_fill_gradient(
    low = "#EDF3ED", high = "#176B3A", transform = "sqrt",
    breaks = c(0, 5, 10, 20, 40, 80), name = "Bee reports"
  ) +
  ggplot2::scale_color_manual(
    values = c("High-High" = "#C23B42", "Low-Low" = "#2C7FB8",
               "High-Low" = "#F28E2B", "Low-High" = "#7A3E9D",
               "Not significant" = NA),
    drop = FALSE, name = "LISA type"
  ) +
  ggplot2::coord_sf(datum = NA, expand = FALSE) +
  ggplot2::labs(
    title = "Kent County | Urban58 | Aug 2023–2025",
    subtitle = "Pooled bee reports and exploratory LISA categories",
    caption = paste(
      "Outlines: exploratory LISA (unadjusted p < 0.05).",
      "No significant clusters after BH correction.",
      sep = "\n"
    )
  ) +
  ggplot2::theme_void(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 17, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 11),
    plot.caption = ggplot2::element_text(size = 10, color = "grey40", hjust = 0),
    plot.margin = ggplot2::margin(20, 25, 20, 20),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 11),
    legend.text = ggplot2::element_text(size = 10)
  )

ggplot2::ggsave(
  file.path(OUT_DIR, "Kent_Urban58_Pooled_AllLISA.png"),
  plot = p, width = 11, height = 9, dpi = 200, bg = "white"
)
cat("Saved pooled CSV and final figure under: ", OUT_DIR, "\n", sep = "")
