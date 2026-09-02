"""Fetch hourly 2025 model snowfall at the four DMI stations -- the model leg.

Source     open-meteo Historical Weather API (keyless), https://archive-api.open-meteo.com/v1/archive
Model      era5  ECMWF ERA5 0.25 deg (~25 km) reanalysis, hourly. The standard reference
                 reanalysis and the forcing of every regional Greenland model.
                 (The ECMWF IFS 9 km archive was dropped from this project on 2026-08-28;
                 re-adding it is one entry in MODELS, but every downstream notebook now
                 assumes a single model series.)
Variables  snowfall_water_equivalent (mm w.e.), precipitation (mm), rain (mm), temperature_2m (degC)
Stamps     every value is the sum over the PRECEDING hour, stamped at the hour end -- the same
           convention as DMI parameter 601 (covering interval (stamp-1h, stamp]), so the model and
           the gauge join on `ts` exactly, with no shifting.
Precision  values are stored at 0.1 mm, i.e. the same floor as the DMI gauge; "model snowing"
           (> 0) is therefore identical to (>= 0.1 mm/h).
Cell       cell_selection="land": open-meteo picks a LAND grid cell with an elevation similar to the
           station's (90 m DEM), 2.5-23.5 km from the station for ERA5. Neighbouring cells differ by
           up to ~27 % in annual snow, so the cell actually used is written to
           openmeteo_cells_2025.csv.
Checked    2026-08-28 -- Aasiaat annual snow w.e.: ERA5 107.9 mm (regression value).

Outputs    data/model/openmeteo_hourly_2025.csv   site, ts, model, snowfall_mm, precip_mm, rain_mm, t2m_C
           data/model/openmeteo_cells_2025.csv    site, model, cell_lat, cell_lon, cell_elev_m, offset_km,
                                                  station_lat, station_lon
Run from anywhere: paths are resolved relative to this file."""

from pathlib import Path

import numpy as np
import pandas as pd
import requests

PROJECT_DIR = Path(__file__).resolve().parent.parent
OUT_DIR = PROJECT_DIR / "data" / "model"
API_URL = "https://archive-api.open-meteo.com/v1/archive"
YEAR = 2025
MODELS = ["era5"]
# open-meteo variable name -> column name in the CSV
VARIABLES = {"snowfall_water_equivalent": "snowfall_mm", "precipitation": "precip_mm",
             "rain": "rain_mm", "temperature_2m": "t2m_C"}
# same keys and coordinates as SITES in compare_earthcare_dmi_stations.ipynb
SITES = {
    "aasiaat": (68.7081, -52.8517),
    "nuuk": (64.1833, -51.7308),
    "tasiilaq": (65.6111, -37.6367),
    "ittoqqortoormiit": (70.4844, -21.9511),
}


def great_circle_km(lat1, lon1, lat2, lon2):
    """Great-circle distance in km between two points given in degrees."""
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = (np.sin((lat2 - lat1) / 2) ** 2
         + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2)
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def fetch_model(model_id):
    """One GET for all four stations and the whole year for one model; returns the
    list of per-location JSON blocks (one call per model so each block carries
    that model's own grid-cell coordinates and elevation)."""
    response = requests.get(API_URL, params=dict(
        latitude=",".join(f"{lat:.4f}" for lat, _ in SITES.values()),
        longitude=",".join(f"{lon:.4f}" for _, lon in SITES.values()),
        start_date=f"{YEAR}-01-01", end_date=f"{YEAR}-12-31",
        hourly=",".join(VARIABLES), models=model_id,
        timezone="UTC", cell_selection="land",
    ), timeout=180)
    response.raise_for_status()
    blocks = response.json()
    return blocks if isinstance(blocks, list) else [blocks]


def site_of(block):
    """Name of the station nearest to a returned grid-cell centre -- the response
    order is not trusted (stations are hundreds of km apart, so this is unambiguous)."""
    distances = {name: great_circle_km(lat, lon, block["latitude"], block["longitude"])
                 for name, (lat, lon) in SITES.items()}
    return min(distances, key=distances.get)


# --- fetch, one call per model ---
hourly_frames, cell_rows = [], []
for model_id in MODELS:
    blocks = fetch_model(model_id)
    matched = [site_of(block) for block in blocks]
    assert sorted(matched) == sorted(SITES), f"{model_id}: sites matched {matched}"

    for site_name, block in zip(matched, blocks):
        hourly = block["hourly"]
        frame = pd.DataFrame({"site": site_name,
                              "ts": pd.to_datetime(hourly["time"]),
                              "model": model_id})
        for api_name, column in VARIABLES.items():
            frame[column] = pd.to_numeric(pd.Series(hourly[api_name]), errors="coerce")
        hourly_frames.append(frame)

        # --- the cell actually used, and its offset from the station ---
        station_lat, station_lon = SITES[site_name]
        cell_rows.append(dict(site=site_name, model=model_id,
                              cell_lat=block["latitude"], cell_lon=block["longitude"],
                              cell_elev_m=block.get("elevation"),
                              offset_km=round(float(great_circle_km(
                                  station_lat, station_lon,
                                  block["latitude"], block["longitude"])), 1),
                              station_lat=station_lat, station_lon=station_lon))

model_hourly = pd.concat(hourly_frames, ignore_index=True)
model_cells = pd.DataFrame(cell_rows)

# --- guards: complete hourly year, no gaps in snowfall, 0.1 mm floor ---
for (site_name, model_id), rows in model_hourly.groupby(["site", "model"]):
    assert len(rows) == 8760, (site_name, model_id, len(rows))
    assert rows["snowfall_mm"].notna().all(), (site_name, model_id, "NaN in snowfall_mm")
    smallest_positive = rows.loc[rows["snowfall_mm"] > 0, "snowfall_mm"].min()
    assert np.isclose(smallest_positive, 0.1), (site_name, model_id, smallest_positive)

# --- report ---
print(f"{'site':<17} {'model':<10} {'snow w.e. mm':>13} {'precip mm':>10} "
      f"{'snow hours':>11} {'cell offset':>12} {'cell elev':>10}")
for (site_name, model_id), rows in model_hourly.groupby(["site", "model"]):
    cell = model_cells[(model_cells.site == site_name) & (model_cells.model == model_id)].iloc[0]
    print(f"{site_name:<17} {model_id:<10} {rows['snowfall_mm'].sum():>13.1f} "
          f"{rows['precip_mm'].sum():>10.1f} {int((rows['snowfall_mm'] > 0).sum()):>11} "
          f"{cell['offset_km']:>9.1f} km {cell['cell_elev_m']:>8.0f} m")
print("NaN counts:", model_hourly[list(VARIABLES.values())].isna().sum().to_dict())

# --- write ---
OUT_DIR.mkdir(parents=True, exist_ok=True)
hourly_path = OUT_DIR / f"openmeteo_hourly_{YEAR}.csv"
cells_path = OUT_DIR / f"openmeteo_cells_{YEAR}.csv"
model_hourly.to_csv(hourly_path, index=False)
model_cells.to_csv(cells_path, index=False)
print(f"wrote {hourly_path.relative_to(PROJECT_DIR)}: {len(model_hourly):,} rows")
print(f"wrote {cells_path.relative_to(PROJECT_DIR)}: {len(model_cells)} rows")
