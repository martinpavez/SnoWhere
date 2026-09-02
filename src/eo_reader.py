"""Reader for EarthCARE ACM-CAP Level 2B granules (data/eo/*.h5).

ACM-CAP is a nadir "curtain": each granule is a strip of ~5000 profiles spaced
~1 km along the ground track, each profile holding 242 fixed height levels
(top-down, ~40 km down to -500 m). It is NOT a map -- there is no swath, so a
granule only tells you about the ~1 km wide line the satellite flew over.

Sentinel conventions that bite if ignored:
  9.96921e36  _FillValue on every float field
  -100 dBZ    CPR_reflectivity_factor: "no detectable echo", not a measurement
  0.0         ice_mass_flux / water contents in clear air -- real zeros, not gaps
  height < elevation   gates below the terrain; fixed grid, not masked for you
  ~0.7 km AGL CPR surface-clutter blind zone: no usable radar echo below it
"""

import glob
import os
import re
from datetime import datetime, timedelta

import h5py
import numpy as np
import pandas as pd

EO_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "data", "eo")

FILL = 1e30
NO_ECHO_DBZ = -100.0
EPOCH = datetime(2000, 1, 1)   # ACM-CAP time units: seconds since 2000-01-01

# CPR near-surface blind zone. Kim et al. 2025 (ACP 25, 15389) filter everything
# below 600 m; the Orbital-Radar simulator paper (Pfitzenmaier et al. 2025,
# GMD 18, 101) notes the 3.3 us pulse puts surface echo up to 500 m above ground.
# Measured directly in data/eo: lowest real echo sits at 540-720 m AGL. 700 m is
# the conservative choice -- drop to 600 to match Kim et al. exactly.
CLUTTER_AGL_M = 700.0

# Collocation radius. The three published CPR ground-validation studies split into
# two camps: 100 km with a +/-30 min window (Kim et al. 2025; Pfitzenmaier et al.
# 2026, EGUsphere 2026-4025), or 200 km with no temporal window at all, matching
# statistically instead (Feuillard et al. 2026, EGUsphere 2026-925).
# Note what widening actually buys: for Disko Bay it adds profiles within the same
# overpass (60 -> 378 going 50 -> 200 km), not extra overpasses. Orbit spacing,
# not the radius, sets how many overpasses you get.
DEFAULT_RADIUS_KM = 100.0
DEFAULT_WINDOW_MIN = 30.0

# Sampling snowfall at a FIXED height above ground, rather than at whatever the
# lowest valid gate happens to be, keeps profiles comparable to each other and to
# CloudSat's near-surface-bin convention (~1200 m over land). Souverijns et al.
# 2018 (TC 12, 3775) measured a 25% underestimate comparing CloudSat at 1200 m
# against an MRR at 300 m -- that gap is the blind zone, not retrieval error.
REFERENCE_AGL_M = 1000.0

# 1 kg m-2 s-1 of mass flux == 1 mm s-1 of liquid water equivalent
KGM2S_TO_MMHR = 3600.0

GRANULE_RE = re.compile(r"_(\d{5})([A-H])\.h5$")


def granule_id(path):
    """'09053C' from a full ACM-CAP filename."""
    m = GRANULE_RE.search(os.path.basename(path))
    return f"{m.group(1)}{m.group(2)}" if m else os.path.basename(path)


def list_granules(eo_dir=EO_DIR):
    return sorted(glob.glob(os.path.join(eo_dir, "ECA_*_ACM_CAP_2B_*.h5")))


def read_granule(path, variables=None):
    """Load an ACM-CAP granule into a dict of masked numpy arrays.

    Floats get _FillValue -> NaN. Reflectivity's -100 dBZ no-echo sentinel is
    turned into NaN as well, so `nanmean` over it means "mean of real echoes".
    ice_mass_flux keeps its clear-air zeros: they are physically meaningful
    (no snow there), so masking them here would bias any average upward.
    """
    default = [
        "latitude", "longitude", "time", "height", "elevation", "temperature",
        "CPR_reflectivity_factor", "CPR_doppler_velocity",
        "ice_mass_flux", "ice_water_content", "ice_water_path",
        "ice_effective_radius", "ice_riming_factor",
        "rain_rate", "rain_classification", "liquid_classification",
        "ice_mass_flux_error",
        "quality_status", "synergy_status", "convergence_status",
    ]
    variables = variables or default
    out = {}
    with h5py.File(path, "r") as f:
        sd = f["ScienceData"]
        for v in variables:
            if v not in sd:
                continue
            a = sd[v][:]
            if np.issubdtype(a.dtype, np.floating):
                a = a.astype("float64")
                a[a > FILL] = np.nan
            out[v] = a
    if "CPR_reflectivity_factor" in out:
        z = out["CPR_reflectivity_factor"]
        out["CPR_reflectivity_factor"] = np.where(z <= NO_ECHO_DBZ, np.nan, z)
    if "time" in out:
        out["datetime"] = np.array(
            [EPOCH + timedelta(seconds=float(s)) for s in out["time"]])
    if "height" in out and "elevation" in out:
        out["height_agl"] = out["height"] - out["elevation"][:, None]
    out["granule"] = granule_id(path)
    return out


def great_circle_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    return 2 * 6371.0 * np.arcsin(np.sqrt(
        np.sin((lat2 - lat1) / 2) ** 2
        + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2))


def surface_snowfall_rate(g, reference_agl_m=REFERENCE_AGL_M,
                          clutter_agl_m=CLUTTER_AGL_M, mode="reference"):
    """Snowfall rate (mm/h water equivalent) standing in for the surface rate.

    ACM-CAP has no surface-snowfall variable, so you have to pick a height and
    declare it. Two modes:

      "reference" (default) -- sample at the level nearest `reference_agl_m`
          above ground. Every profile is then measured at the same altitude,
          which is what makes profiles comparable to each other and to
          CloudSat's near-surface-bin convention.
      "lowest"    -- take the lowest gate above the clutter floor. Gets you
          closer to the ground but at a height that varies profile to profile.

    Why not simply take the bottom gate: the retrieval's along-track smoother
    fills in ice_mass_flux BELOW the lowest real radar echo -- 864 of 4707
    profiles in granule 09053C. Those numbers are not fill-valued and look
    perfectly valid. The product's own error estimate is what gives them away:
    median ln-error runs 0.09 above 1 km AGL but 0.37 in the lowest 300 m.
    """
    imf = g["ice_mass_flux"]
    agl = g["height_agl"]
    rate = np.full(imf.shape[0], np.nan)
    gate_h = np.full(imf.shape[0], np.nan)

    if mode == "reference":
        for p in range(imf.shape[0]):
            ok = np.where(np.isfinite(imf[p]) & (agl[p] >= clutter_agl_m))[0]
            if not ok.size:
                continue
            k = ok[int(np.argmin(np.abs(agl[p, ok] - reference_agl_m)))]
            rate[p] = imf[p, k] * KGM2S_TO_MMHR
            gate_h[p] = agl[p, k]
    elif mode == "lowest":
        usable = np.isfinite(imf) & (agl >= clutter_agl_m)
        for p in range(imf.shape[0]):
            ok = np.where(usable[p])[0]
            if ok.size:
                k = ok[-1]                  # top-down ordering: last = lowest
                rate[p] = imf[p, k] * KGM2S_TO_MMHR
                gate_h[p] = agl[p, k]
    else:
        raise ValueError("mode must be 'reference' or 'lowest'")
    return rate, gate_h


def collocate(path, station_lat, station_lon, radius_km=DEFAULT_RADIUS_KM,
              clutter_agl_m=CLUTTER_AGL_M):
    """Profiles of one granule falling within radius_km of a station.

    Returns a DataFrame (one row per profile) or None if the track never comes
    close enough. The default 100 km follows Kim et al. 2025 and Pfitzenmaier
    et al. 2026; widening it buys more profiles from the same overpass at the
    cost of comparing the satellite to weather the station never saw.
    """
    g = read_granule(path)
    dist = great_circle_km(station_lat, station_lon, g["latitude"], g["longitude"])
    sel = np.where(dist <= radius_km)[0]
    if sel.size == 0:
        return None

    rate, gate_h = surface_snowfall_rate(g, clutter_agl_m=clutter_agl_m)
    z = g["CPR_reflectivity_factor"]
    agl = g["height_agl"]
    near = np.isfinite(z) & (agl >= clutter_agl_m) & (agl < 3000)

    z_low = np.full(z.shape[0], np.nan)
    for p in sel:
        ok = np.where(near[p])[0]
        if ok.size:
            z_low[p] = z[p, ok[-1]]

    return pd.DataFrame({
        "granule": g["granule"],
        "datetime": g["datetime"][sel],
        "latitude": g["latitude"][sel],
        "longitude": g["longitude"][sel],
        "distance_km": dist[sel],
        "elevation_m": g["elevation"][sel],
        "snowfall_mmhr": rate[sel],
        "gate_agl_m": gate_h[sel],
        "reflectivity_dbz": z_low[sel],
        "ice_water_path": g["ice_water_path"][sel],
        "quality_status": g["quality_status"][sel],
        "synergy_status": g["synergy_status"][sel],
        "convergence_status": g["convergence_status"][sel],
    })


def converged(df):
    """Keep only profiles where the retrieval actually converged.

    Not optional: in granule 09053C only 74% of profiles converged (16% "empty
    state", 9% "failed", 1% "error"). The non-converged ones still carry
    plausible-looking numbers in every field.
    """
    return df[df["convergence_status"] == 0]


def overpass_summary(station_lat, station_lon, radius_km=DEFAULT_RADIUS_KM,
                     eo_dir=EO_DIR):
    """One row per granule: does its ground track reach the station, and what
    did it see if so."""
    rows = []
    for path in list_granules(eo_dir):
        g = read_granule(path, ["latitude", "longitude", "time"])
        dist = great_circle_km(station_lat, station_lon,
                               g["latitude"], g["longitude"])
        i = int(np.argmin(dist))
        row = {
            "granule": g["granule"],
            "start_utc": g["datetime"][0],
            "closest_km": dist[i],
            "closest_lat": g["latitude"][i],
            "closest_lon": g["longitude"][i],
            "within_radius": bool(dist.min() <= radius_km),
        }
        if row["within_radius"]:
            df = collocate(path, station_lat, station_lon, radius_km)
            row["n_profiles"] = len(df)
            row["snowing_profiles"] = int((df["snowfall_mmhr"] > 1e-4).sum())
            row["median_snow_mmhr"] = float(df["snowfall_mmhr"].median())
            row["max_snow_mmhr"] = float(df["snowfall_mmhr"].max())
        rows.append(row)
    return pd.DataFrame(rows)


def open_subset(subset_dir=None, region=None):
    """Open a consolidated subset (made by eo_subset.py) as one xarray
    Dataset, with the same masking conventions read_granule applies to full
    granules: _FillValue and the -100 dBZ no-echo sentinel become NaN, while
    ice_mass_flux keeps its physically-real clear-air zeros.

    region=None opens the original (unprefixed) Disko files; a region name
    (e.g. "nuuk") opens that region's eo_subset_<region>_*.nc files —
    typically in data/eo_subset_stations. Adds height_agl; time comes back
    CF-decoded as datetime64.
    """
    import xarray as xr

    subset_dir = subset_dir or os.path.join(
        os.path.dirname(EO_DIR),
        "eo_subset" if region is None else "eo_subset_stations")
    pattern = (f"eo_subset_{region}_*.nc" if region
               else "eo_subset_[0-9]*.nc")   # digits: exclude region files
    paths = sorted(glob.glob(os.path.join(subset_dir, pattern)))
    if not paths:
        raise FileNotFoundError(f"no {pattern} in {subset_dir} -- "
                                "run eo_subset.py first")
    # eager load: the whole subset is < 2 GB by design, so no dask needed
    parts = [xr.load_dataset(p) for p in paths]
    ds = (parts[0] if len(parts) == 1
          else xr.concat(parts, dim="profile", combine_attrs="drop_conflicts"))
    if "CPR_reflectivity_factor" in ds:
        z = ds["CPR_reflectivity_factor"]
        ds["CPR_reflectivity_factor"] = z.where(z > NO_ECHO_DBZ)
    # Non-converged profiles can carry ZERO-FILLED measured columns -- exactly
    # 0.0 dBZ at every gate, which renders as a plausible solid echo block
    # (seen as 40 consecutive status-6 profiles in granule 04525C). Mask the
    # measured fields for profiles matching that signature; genuine echoes in
    # partially-failed profiles are left alone.
    if "convergence_status" in ds and "CPR_reflectivity_factor" in ds:
        z = ds["CPR_reflectivity_factor"]
        zero_filled = ((ds["convergence_status"] != 0)
                       & (np.abs(z).max("level") == 0))
        for v in ("CPR_reflectivity_factor", "CPR_doppler_velocity"):
            if v in ds:
                ds[v] = ds[v].where(~zero_filled)
    ds["height_agl"] = ds["height"] - ds["elevation"]
    return ds


# The four active DMI stations (coords from results/data_analysis.md).
# GEM AWS2 is deliberately absent: its record ends 2024-12-31, so it cannot be
# compared against 2025 EO data at all. Of these four, only 422000 Aasiaat
# measures precipitation -- the SAVS airport stations essentially do not, so
# they can corroborate temperature and cloud but not snowfall itself.
DMI_STATIONS = {
    "421900 Qeqertarsuaq Heli.": (69.2514, -53.5147),   # fails QC in both periods
    "422000 Aasiaat":            (68.7081, -52.8517),   # only precipitation record
    "422100 Mitt. Ilulissat":    (69.2403, -51.0661),
    "422400 Mitt. Aasiaat":      (68.7219, -52.7847),   # record ends 2025-08-07
}


if __name__ == "__main__":
    print(f"{len(list_granules())} granules in {EO_DIR}\n")

    for name, (lat, lon) in DMI_STATIONS.items():
        summary = overpass_summary(lat, lon, radius_km=DEFAULT_RADIUS_KM)
        hits = summary[summary["within_radius"]]
        print(f"=== {name} ===")
        print(f"  closest track overall : {summary['closest_km'].min():.0f} km")
        print(f"  granules within {DEFAULT_RADIUS_KM:.0f} km : "
              f"{len(hits)}  {', '.join(hits['granule']) if len(hits) else '(none)'}")

        for gid in hits["granule"]:
            path = [p for p in list_granules() if granule_id(p) == gid][0]
            df = collocate(path, lat, lon)
            ok = converged(df)
            snow = ok["snowfall_mmhr"]
            print(f"    {gid} at {df['datetime'].iloc[0]:%H:%M} UTC — "
                  f"{len(df)} profiles, {len(ok)} converged")
            print(f"      snowfall at {REFERENCE_AGL_M:.0f} m AGL: "
                  f"max {snow.max():.5f} mm/h, "
                  f"{int((snow > 1e-4).sum())} profiles above 1e-4")
        print()
