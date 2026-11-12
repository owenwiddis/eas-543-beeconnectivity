# Run from repository root AFTER Kent_BeeCount_Grid_LISA.R.
out_dir <- "outputs/figures/kent_2025_08"
counts <- read.csv(file.path(out_dir, "Kent_hex_6km_counts.csv"))
lisa <- read.csv(file.path(out_dir, "Kent_hex_6km_LISA.csv"))
actual <- c(
  all_records = sum(counts$bee_count),
  hex_cells = nrow(counts),
  high_high_cells = sum(lisa$LISA == "High-High", na.rm = TRUE),
  zero_cells = sum(counts$bee_count == 0)
)
expected <- c(all_records = 208L, hex_cells = 95L, high_high_cells = 5L, zero_cells = 59L)
print(actual)
if (identical(as.integer(actual), as.integer(expected))) {
  cat("KENT PILOT CHECK PASSED\n")
} else {
  warning("Kent pilot values differ from the previous snapshot. Check source data, grid, and packages before pushing.")
}
