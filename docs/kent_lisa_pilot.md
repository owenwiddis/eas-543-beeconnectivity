# Kent County 6 km hexagon bee records and exploratory LISA

This is a separately tested pilot script and **does not replace** the repository's `R/01_...` through `R/04_...` skeleton pipeline. Counts are iNaturalist reports, not abundance.

## Inputs (local only, do not push)

From the earlier `Bee` folder, copy into this repository:

- `data/pilot_inputs/Apis_mellifera.csv`
- `data/pilot_inputs/Bombus_impatiens.csv`
- `data/pilot_inputs/Xylocopa_virginica.csv`
- `data/pilot_inputs/kent_2024.geojson` **and all accompanying shapefile files** (.shx, .dbf, .prj, etc.).

`data/raw/` is ignored by Git. This pilot script **does not need the NLCD raster**.

## Reproduce in RStudio

Set working directory to the root of the cloned Git repository (the folder containing `README.md`). Then run:

```r
source("R/pilots/Kent_BeeCount_Grid_LISA.R")
source("R/pilots/check_kent_lisa_outputs.R")
```

The script installs missing `sf`, `ggplot2`, `spdep` packages and writes figures, per-cell tables, and diagnostics to `outputs/figures/kent_2025_08/`. A GeoPackage file is created locally but ignored by Git.

## Expected August 2025 pilot snapshot

- Three species combined inside Kent County: **208** (Apis 36, Bombus 159, Xylocopa 13).
- 6 km hexagon cells: **95**, including **59** with zero reports.
- 6 km exploratory LISA: **5** BH-adjusted p < 0.05 Highâ€“High cells, **90** nonsignificant.
- 4 and 8 km grids are count diagnostics only; no LISA analysis was run at those scales in this script.

The script uses 2024 Census cartographic boundary files and EPSG:5070. Results can vary if input CSV snapshots, packages, or boundary files differ. The hexagon `cellsize` argument should be described consistently with `sf` documentation; scale selection is exploratory. Zeros mean no *reported observations*, not demonstrated absence.

## What to commit

Commit the two `R/pilots/*.R` scripts, this README, and approved `.png` figures. Do not commit CSV observations, full NLCD TIFFs, shapefiles, or temporary GeoPackages.

