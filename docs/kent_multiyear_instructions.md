# Kent County LISA — August 2023–2025 | Urban Focus (58 hexagons)

## What changed

This is the replacement for the unsuccessful 30-hex Shrub Focus analysis. **Do not use the old shrub-focused outputs as the current result.** The original 95-cell Full County LISA is preserved as the primary comparison. The **Urban Focus is exploratory** because it was defined after looking at bee observation distributions.

**Study-area method (fixed across all three years):** Using the NLCD 2025 Land Cover raster and the original identical 6 km hex grid, calculate % Developed (NLCD classes 21–24). Retain hexes with Developed >= **30%** (36 hexes); identify the largest queen-connected component (**33 hexes**); and add every immediate queen-neighbor in the original 95-cell grid (**58 hexes**). Note that this threshold was explored after observing where bees were reported, so reduced-area LISA must not be interpreted as independent confirmation. Do not tune the boundary to obtain significant clusters.

The frozen IDs and geometry are shared as small files. Teammates **do not need NLCD 2025 GeoTIFF** to reproduce counts or LISA for the fixed region (but the original raster is needed to independently recreate the selection itself).

## Required inputs (under `data/pilot_inputs/`)

```
kent_bees_2023_08.csv           # 95 Kent-only, August, 3 species
kent_bees_2024_08.csv           # 116
kent_bees_2025_08.csv           # 208, validated single-year baseline
kent_2024.geojson              # original county polygon, GEOID=26081; unchanged
kent_urban_focus_58hex_ids.csv  # fixed 58 cell_id values
kent_urban_focus_58hex.geojson  # fixed 58 hexagon features with cell_id
```

The observation CSVs contain only `uuid, observed_on, longitude, latitude, species`, are de-duplicated by UUID, filtered by Kent polygon and August. They contain reported sightings, not population estimates or validated bee absences. The GeoJSON and mask derive from the identical 95-cell grid; the author already validated that the 58 cells form one connected component. A separate `kent_urban_focus_area.geojson` is generated automatically by R as a dissolved, clipped map area; it need not be downloaded.

## How to run

From the repository root in RStudio Console:

```r
source("R/pilots/Kent_BeeCount_LISA_2023_2025.R")
```

The script checks the fixed 58 IDs and geometry, verifies the 2025 baseline against the original single-year 95-cell pilot (if available), and computes Full County (95 cells) and Urban Focus (58 cells) LISA using 999 permutations, queen adjacency, and BH-corrected p-values. **It does not ask for an NLCD raster**.

Expected input totals and Urban Focus totals (previous interactive audit):

| Year | Full Kent | Urban58 | Retained |
|---|---:|---:|---:|
| 2023 | 95 | 94 | 98.9% |
| 2024 | 116 | 115 | 99.1% |
| 2025 | 208 | 203 | 97.6% |

## Outputs

All generated analyses are under `outputs/figures/kent_2023_2025/`:

- `Kent_6km_counts_<YEAR>.png` — same full county hexes and scale across years.
- `Kent_6km_full_county_LISA_<YEAR>.png/.csv` — primary, unchanged 95-cell study.
- `Kent_6km_urban_focus_counts_<YEAR>.png` — 58-cell count maps, zero-report cells included.
- `Kent_6km_urban_focus_58hex_LISA_<YEAR>.png/.csv` — exploratory fixed 58-cell local Moran's I.
- `urban_focus_58hex_mask.csv` — full 95-cell inclusion mask.
- `urban_focus_high_high_by_year.csv`, `per_hex_counts_and_changes.csv`, `year_summary.csv` — comparative tables.
- `data/pilot_inputs/kent_urban_focus_area.geojson` — dissolved map boundary (generated locally).

Earlier Shrub Focus and Interior Core artifacts in the output directory are retired by the new script. The original `data/pilot_inputs/kent_2024.geojson` and all source CSV files remain intact.

## Exploratory pooled LISA (the final HH / LL overlay)

After running `R/pilots/Kent_BeeCount_LISA_2023_2025.R`, run this standalone script from the **repository root**:

```r
source("R/pilots/Kent_Urban58_Pooled_LISA.R")
```

This pools August reports across 2023–2025 for the same **58 hexagons** (expected: 412 reports). It runs permutation-based Local Moran's I using 999 permutations, queen neighbors, and a fixed random seed. It saves both **unadjusted** (`p_perm < 0.05`) and **BH-adjusted** (`p_BH < 0.05`) classifications in:

- `Kent_Urban58_LISA_Pooled_2023_2025.csv`
- `Kent_Urban58_Pooled_AllLISA.png`

In the author's exploratory run, **6 High–High and 4 Low–Low** cells met the *unadjusted* cutoff, but **zero** cells survived BH correction. The six HH hexagons were connected. The colored outlines in the final figure refer **only to unadjusted, exploratory categories**; they are not confirmatory hotspots. The plot caption explicitly discloses this. Do not select study boundaries or p-value procedures to maximize significance. The pooled analysis concerns observations in August across three years, not a full-year bee abundance trend.

## Sharing the minimal inputs with the team

After the author has saved the 58 IDs and GeoJSON, run in PowerShell from the repository root:

```powershell
python python/download_kent_bees_2023_2024.py --zip-only
```

This **does not redownload any bees**. It checks and creates `data/kent_pilot_inputs_2023_2025_aug.zip` containing the six small files above. Share via a **restricted** Google Drive folder, with observation coordinate sensitivity and iNaturalist data-license/redistribution checks. Team members unzip the archive under their repo's `data/` directory.

The `data/pilot_inputs/` directory remains excluded by `.gitignore`, so these files do not appear on GitHub even after push or PR. Commit code and selected summary maps/tables only; do not `git add -f` sensitive inputs.

## Interpretation

An exploration-selected focus changes the spatial support, mean counts, and neighborhood graph. It is valid for descriptive or exploratory comparisons but not for claiming greater statistical evidence from higher significance. The 2025 NLCD map is a fixed landscape descriptor, not reconstructed land cover for 2023 or 2024. iNaturalist counts confound underlying occurrence and observer effort. Do not interpret LISA patterns as abundance, absence, or bee movement.
