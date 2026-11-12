"""Prepare the SMALL, shareable Kent County / August / 3-species inputs.

Run from the Git repository root:
    python -m pip install requests shapely
    python python/download_kent_bees_2023_2024.py

2023/2024: query iNaturalist (research grade, photos) then clip to the
Kent County GeoJSON BEFORE saving; 2025: clip the existing verified local
August snapshots. Never commit the location CSVs to a public repository.

The archive contains three county-month CSVs, original Kent GeoJSON,
and the frozen 58-hex Urban Focus IDs + GeoJSON. It requires no NLCD TIFF.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from datetime import date
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

try:
    import requests
    from shapely.geometry import Point, shape
except ImportError as exc:
    raise SystemExit("Missing package. Run: python -m pip install requests shapely") from exc

API = "https://api.inaturalist.org/v1/observations"
MONTH = 8
YEARS = (2023, 2024, 2025)
SPECIES = {
    47219: "Apis mellifera",
    118970: "Bombus impatiens",
    51110: "Xylocopa virginica",
}
DATA_DIR = Path("data/pilot_inputs")
COUNTY_FILE = DATA_DIR / "kent_2024.geojson"
ARCHIVE = Path("data/kent_pilot_inputs_2023_2025_aug.zip")
CSV_FIELDS = ("uuid", "observed_on", "longitude", "latitude", "species")
BASELINE_2025_COUNT = 208


def date_in_month(value: str, year: int) -> bool:
    """Accept YYYY-MM-DD, YYYY/M/D, and timestamps starting with a date."""
    try:
        first = str(value).strip().replace("/", "-").split("T", 1)[0].split(" ", 1)[0]
        parts = first.split("-")
        if len(parts) != 3:
            return False
        day = date(*(int(p) for p in parts))
        return day.year == year and day.month == MONTH
    except (ValueError, TypeError):
        return False


def kent_geometry():
    if not COUNTY_FILE.is_file():
        raise FileNotFoundError(f"Missing county boundary: {COUNTY_FILE}")
    with COUNTY_FILE.open(encoding="utf-8") as src:
        doc = json.load(src)
    if doc.get("type") == "FeatureCollection":
        features = doc.get("features", [])
    elif doc.get("type") == "Feature":
        features = [doc]
    else:
        raise ValueError("Expected a GeoJSON Feature or FeatureCollection")
    chosen = [feat for feat in features if str(feat.get("properties", {}).get("GEOID")) == "26081"]
    if len(chosen) != 1:
        raise ValueError("Kent GeoJSON must contain exactly one feature with GEOID=26081")
    geometry = shape(chosen[0]["geometry"])
    if geometry.is_empty or not geometry.is_valid:
        raise ValueError("Invalid or empty Kent polygon")
    return geometry


def in_county(value_lon, value_lat, county) -> bool:
    try:
        lon, lat = float(value_lon), float(value_lat)
    except (ValueError, TypeError):
        return False
    if not (math.isfinite(lon) and math.isfinite(lat) and abs(lon) <= 180 and abs(lat) <= 90):
        return False
    return county.covers(Point(lon, lat))


def read_original_2025(county):
    """Use, but NEVER modify, the three original 2025 input snapshots."""
    rows = []
    for species in SPECIES.values():
        filename = DATA_DIR / (species.replace(" ", "_") + ".csv")
        if not filename.is_file():
            raise FileNotFoundError(
                f"Missing original 2025 snapshot {filename}. "
                "It is only needed once to create kent_bees_2025_08.csv."
            )
        with filename.open(newline="", encoding="utf-8-sig") as src:
            reader = csv.DictReader(src)
            missing = set(CSV_FIELDS[:4]) - set(reader.fieldnames or [])
            if missing:
                raise ValueError(f"{filename}: missing columns {sorted(missing)}")
            for entry in reader:
                if not date_in_month(entry.get("observed_on"), 2025):
                    continue
                if not in_county(entry.get("longitude"), entry.get("latitude"), county):
                    continue
                rows.append({
                    "uuid": entry["uuid"],
                    "observed_on": entry["observed_on"],
                    "longitude": entry["longitude"],
                    "latitude": entry["latitude"],
                    "species": species,
                })
    return deduplicate(rows)


def download_species(session, county, year, taxon_id, species):
    """Query only county BBOX; save only polygon-inside observations."""
    minlon, minlat, maxlon, maxlat = county.bounds
    results = []
    page = 1
    while True:
        params = {
            "taxon_id": taxon_id,
            "d1": f"{year}-{MONTH:02d}-01",
            "d2": f"{year}-{MONTH:02d}-31",
            "swlat": minlat,
            "swlng": minlon,
            "nelat": maxlat,
            "nelng": maxlon,
            "quality_grade": "research",
            "has[]": "photos",
            "per_page": 200,
            "page": page,
            "order_by": "id",
            "order": "asc",
        }
        response = session.get(API, params=params, timeout=60)
        response.raise_for_status()  # Do not quietly save partial data on HTTP errors.
        payload = response.json()
        total = int(payload.get("total_results", 0))
        if total > 10000:
            raise RuntimeError(f"{year} {species}: >10,000 API matches; split date queries")
        observations = payload.get("results", [])
        if not isinstance(observations, list):
            raise ValueError("Unexpected iNaturalist API response")
        for obs in observations:
            coords = (obs.get("geojson") or {}).get("coordinates")
            if not coords or len(coords) != 2:
                continue
            lon, lat = coords
            if not in_county(lon, lat, county):
                continue
            if not date_in_month(obs.get("observed_on"), year):
                continue
            uuid = obs.get("uuid")
            if not uuid:
                # Stable fallback only for genuinely missing API UUIDs.
                uuid = f"inat-{obs['id']}"
            results.append({
                "uuid": uuid,
                "observed_on": obs["observed_on"],
                "longitude": lon,
                "latitude": lat,
                "species": species,
            })
        print(f"  {year} {species}: page {page}, {len(results)} inside Kent")
        if page * 200 >= total or not observations:
            break
        page += 1
        time.sleep(1.1)  # Follow iNaturalist's recommended ~1 request/s.
    return results


def deduplicate(rows):
    unique = {}
    for entry in rows:
        uid = str(entry["uuid"]).strip()
        if uid and uid not in unique:
            unique[uid] = entry
    return list(unique.values())


def output_path(year):
    return DATA_DIR / f"kent_bees_{year}_{MONTH:02d}.csv"


def write_year(year, rows):
    path = output_path(year)
    if not rows:
        raise RuntimeError(f"No Kent County August observations found for {year}; refusing to save an empty snapshot")
    if year == 2025 and len(rows) != BASELINE_2025_COUNT:
        raise RuntimeError(
            f"2025 verified baseline expected {BASELINE_2025_COUNT} Kent records, "
            f"but got {len(rows)}. Check inputs and polygon before saving."
        )
    if any(not date_in_month(r["observed_on"], year) for r in rows):
        raise ValueError("Year/month audit failed")
    with path.open("w", newline="", encoding="utf-8") as dst:
        writer = csv.DictWriter(dst, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved {len(rows)} Kent-only August observations: {path}")


def zip_small_share():
    wanted = [COUNTY_FILE] + [output_path(y) for y in YEARS]
    missing = [str(path) for path in wanted if not path.is_file()]
    if missing:
        raise FileNotFoundError("Cannot package missing files: " + ", ".join(missing))
    county = kent_geometry()
    for year in YEARS:
        with output_path(year).open(newline="", encoding="utf-8-sig") as source:
            reader = csv.DictReader(source)
            if tuple(reader.fieldnames or []) != CSV_FIELDS:
                raise ValueError(f"{output_path(year)} must have exactly the five required columns")
            seen = set()
            n = 0
            for row in reader:
                uid = str(row.get("uuid", "")).strip()
                if not uid or uid in seen:
                    raise ValueError(f"Blank/duplicate UUID in {output_path(year)}")
                seen.add(uid)
                if not date_in_month(row["observed_on"], year):
                    raise ValueError(f"Record outside August {year}: {output_path(year)}")
                if not in_county(row["longitude"], row["latitude"], county):
                    raise ValueError(f"Record outside Kent County: {output_path(year)}")
                if row["species"] not in SPECIES.values():
                    raise ValueError(f"Unexpected species in {output_path(year)}")
                n += 1
        if year == 2025 and n != BASELINE_2025_COUNT:
            raise ValueError(f"2025 audit expected {BASELINE_2025_COUNT} rows; got {n}")
        print(f"Share audit passed: {year} / August / Kent / {n} rows")
    ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
    focus_ids = DATA_DIR / "kent_urban_focus_58hex_ids.csv"
    focus_geometry = DATA_DIR / "kent_urban_focus_58hex.geojson"
    for path in (focus_ids, focus_geometry):
        if not path.is_file():
            raise FileNotFoundError(f"Missing fixed Urban58 input: {path}")
    with focus_ids.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames != ["cell_id"]:
            raise ValueError("Urban58 IDs CSV must contain exactly cell_id")
        focus_rows = list(reader)
    ids = [int(r["cell_id"]) for r in focus_rows]
    if len(ids) != 58 or len(set(ids)) != 58 or not all(1 <= i <= 95 for i in ids):
        raise ValueError("Urban58 mask must contain 58 unique IDs in 1..95")
    with focus_geometry.open(encoding="utf-8") as f:
        aoi = json.load(f)
    feats = aoi.get("features", [])
    if aoi.get("type") != "FeatureCollection" or len(feats) != 58:
        raise ValueError("Urban58 GeoJSON must contain 58 features")
    geo_ids = [int(feat.get("properties", {}).get("cell_id")) for feat in feats]
    if set(geo_ids) != set(ids) or len(set(geo_ids)) != 58:
        raise ValueError("Urban58 GeoJSON cell_ids do not match the IDs CSV")
    if any(shape(feat["geometry"]).is_empty or not shape(feat["geometry"]).is_valid
           for feat in feats):
        raise ValueError("Urban58 GeoJSON contains invalid/empty polygon")
    wanted.extend([focus_ids, focus_geometry])
    print("Including frozen Urban Focus (58 hexes, matching IDs and GeoJSON)")
    with ZipFile(ARCHIVE, "w", ZIP_DEFLATED, compresslevel=8) as archive:
        for path in wanted:
            archive.write(path, arcname=f"pilot_inputs/{path.name}")
    size_kb = ARCHIVE.stat().st_size / 1024
    print(f"Small team-share archive: {ARCHIVE} ({size_kb:.1f} KB)")
    print("Upload this ZIP to restricted Google Drive; unzip under the repo's data/ folder.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true",
                        help="Re-download 2023/2024; does not overwrite the frozen 2025 snapshot")
    parser.add_argument("--zip-only", action="store_true",
                        help="Package existing small inputs without API calls")
    args = parser.parse_args()
    if args.zip_only:
        zip_small_share()
        return

    county = kent_geometry()
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    print("Scope: Kent County polygon | August 2023, 2024, 2025 | three bee species")
    with requests.Session() as session:
        session.headers.update({
            "User-Agent": "EAS543 Kent Bee County-Month Snapshot / educational research",
            "Accept": "application/json",
        })
        for year in YEARS[:2]:
            path = output_path(year)
            if path.exists() and not args.refresh:
                print(f"Keeping previously downloaded frozen snapshot: {path}")
                continue
            rows = []
            for taxon_id, species in SPECIES.items():
                rows.extend(download_species(session, county, year, taxon_id, species))
                time.sleep(1.1)
            write_year(year, deduplicate(rows))

    if output_path(2025).exists():
        print(f"Keeping verified 2025 snapshot: {output_path(2025)}")
    else:
        write_year(2025, read_original_2025(county))
    zip_small_share()
    print("Done. No full-state or full-year input CSV was written.")


if __name__ == "__main__":
    main()
