"""Pull iNaturalist research-grade observations into a pandas DataFrame.

Michigan (place_id=29), March-November 2025, three study species.
Pages with id_above so it isn't capped at 10,000 results.
"""
import time
import requests
import pandas as pd

BASE = "https://api.inaturalist.org/v1/observations"
SPECIES = {
    47219: "Apis mellifera",
    118970: "Bombus impatiens",
    51110: "Xylocopa virginica",
}
PARAMS = {
    "place_id": 29,              # Michigan
    "d1": "2025-03-01",
    "d2": "2025-11-30",
    "quality_grade": "research",
    "per_page": 200,             # API max
    "order_by": "id",
    "order": "asc",
}


def fetch_observations(extra_params, pause=1.0):
    """Return a list of raw observation dicts, paging with id_above."""
    rows, last_id = [], 0
    with requests.Session() as s:
        s.headers["User-Agent"] = "EAS543-bee-project (University of Michigan)"
        while True:
            params = {**PARAMS, **extra_params, "id_above": last_id}
            r = s.get(BASE, params=params, timeout=60)
            r.raise_for_status()
            results = r.json()["results"]
            if not results:
                break
            rows.extend(results)
            last_id = results[-1]["id"]
            print(f"  {len(rows)} records...", end="\r")
            time.sleep(pause)    # stay under iNat's ~1 request/second guideline
    print()
    return rows


def flatten(obs):
    """Keep the fields we need from one raw observation."""
    coords = (obs.get("geojson") or {}).get("coordinates") or [None, None]
    taxon = obs.get("taxon") or {}
    return {
        "id": obs["id"],
        "observed_on": obs.get("observed_on"),
        "taxon_id": taxon.get("id"),
        "taxon_name": taxon.get("name"),
        "lon": coords[0],
        "lat": coords[1],
        "positional_accuracy": obs.get("positional_accuracy"),
        "obscured": obs.get("obscured"),
        "geoprivacy": obs.get("geoprivacy"),
        "taxon_geoprivacy": obs.get("taxon_geoprivacy"),
        "place_guess": obs.get("place_guess"),
        "user_login": (obs.get("user") or {}).get("login"),
        "quality_grade": obs.get("quality_grade"),
        "url": obs.get("uri"),
    }


def to_dataframe(raw):
    df = pd.DataFrame([flatten(o) for o in raw])
    df["observed_on"] = pd.to_datetime(df["observed_on"], errors="coerce")
    df["month"] = df["observed_on"].dt.month
    # Flags for the point-level (moving window) analysis
    df["usable_point"] = (
        ~df["obscured"].fillna(False).astype(bool)
        & df["positional_accuracy"].notna()
        & (df["positional_accuracy"] <= 1000)
    )
    return df


if __name__ == "__main__":
    frames = []
    for tid, name in SPECIES.items():
        print(f"Fetching {name} ({tid})")
        frames.append(to_dataframe(fetch_observations({"taxon_id": tid})))
    bees = pd.concat(frames, ignore_index=True)

    print(bees.groupby("taxon_name").agg(
        records=("id", "size"), usable_points=("usable_point", "sum")))
    bees.to_csv("data/raw/inat_bees_MI_2025.csv", index=False)

    # Optional: target-group background (all insects, ~89k records, ~8 min).
    # bg = to_dataframe(fetch_observations({"iconic_taxa": "Insecta"}))
    # bg = bg[~bg["id"].isin(bees["id"])]          # don't overlap with bee records
    # bg.to_csv("data/raw/inat_insects_MI_2025.csv", index=False)

    # Next step: drop the Upper Peninsula and add county names with a spatial join.
    # import geopandas as gpd
    # gdf = gpd.GeoDataFrame(bees, geometry=gpd.points_from_xy(bees.lon, bees.lat), crs=4326)
    # counties = gpd.read_file("mi_counties.shp").to_crs(4326)
    # gdf = gpd.sjoin(gdf, counties[["NAME", "geometry"]], how="inner", predicate="within")
