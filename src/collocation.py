"""The satellite side of the study and its pairing with a station: the
profile usability screen, the per-profile arrays every notebook works from,
and the standard one-row-per-overpass table.

An "overpass" is one granule with at least one usable profile within the
collocation radius of the station; its timestamp is the mean profile time of
those profiles. Note what widening the radius buys: more profiles from the
same overpass, never more overpasses -- orbit spacing sets those.
"""

import numpy as np
import pandas as pd

import eo_reader as eo
from gauges import load_gauge_snow, snow_window_stats
from study_config import (COLLOCATION_RADIUS_KM, OUTLIER_MMHR,
                          SNOW_THRESHOLD_MMHR, THETA_FSS)


def usable_profiles(ds):
    """Boolean mask of profiles fit for analysis: converged retrieval, finite
    snowfall at the 1000 m reference level, and below the unphysical-rate
    screen. Not optional -- non-converged profiles carry plausible numbers."""
    rate = ds["snowfall_mmhr_1000m"].values
    return ((ds["convergence_status"].values == 0)
            & np.isfinite(rate) & (rate < OUTLIER_MMHR))


def profile_arrays(ds, dist_col):
    """The per-profile numpy arrays the analyses index into, pulled out of the
    xarray dataset once (keys: profile_time, profile_granule,
    profile_distance_km, profile_snowfall_mmhr, profile_usable, profile_lat,
    profile_lon)."""
    return dict(
        profile_time=pd.to_datetime(ds["time"].values),
        profile_granule=ds["granule"].values,
        profile_distance_km=ds[dist_col].values,
        profile_snowfall_mmhr=ds["snowfall_mmhr_1000m"].values,
        profile_usable=usable_profiles(ds),
        profile_lat=ds["latitude"].values,
        profile_lon=ds["longitude"].values,
    )


def build_overpass_table(profiles, station_hourly, gauge_snow_mmhr,
                         radius_km=COLLOCATION_RADIUS_KM):
    """One row per overpass (granule with usable profiles within radius_km):
    its mean time, closest approach, satellite snowing flags at 100/30 km, the
    fraction of profiles above the harmonized threshold, and the gauge's
    snow-phase precip in the +/-30 min window around the overpass."""
    rows = []
    usable = profiles["profile_usable"]
    granule = profiles["profile_granule"]
    distance = profiles["profile_distance_km"]
    for granule_id in sorted(set(granule[usable & (distance <= radius_km)])):
        in_overpass = usable & (granule == granule_id) & (distance <= radius_km)
        if not in_overpass.any():
            continue
        overpass_time = (pd.Series(profiles["profile_time"][in_overpass])
                         .mean().round("min"))
        pass_rates = profiles["profile_snowfall_mmhr"][in_overpass]
        pass_distances = distance[in_overpass]
        snow_mm_30, snow_frac_30, n_stamps_30 = snow_window_stats(
            station_hourly, gauge_snow_mmhr, overpass_time,
            pd.Timedelta("30min"))
        rows.append(dict(
            granule=granule_id, time=overpass_time,
            closest_km=float(pass_distances.min()),
            n_profiles=int(in_overpass.sum()),
            mean_rate=float(pass_rates.mean()),      # zeros included (Palerme)
            sat_any_100km=bool((pass_rates > SNOW_THRESHOLD_MMHR).any()),
            sat_any_30km=bool((pass_rates[pass_distances <= 30]
                               > SNOW_THRESHOLD_MMHR).any())
                         if (pass_distances <= 30).any() else False,
            sat_frac_theta=float((pass_rates > THETA_FSS).mean()),
            gauge_snow_pm30_mm=snow_mm_30, gauge_snow_frac_pm30=snow_frac_30,
            gauge_n_stamps_pm30=n_stamps_30,
        ))
    return pd.DataFrame(rows).sort_values("time").reset_index(drop=True)


def load_site(station, keep_dataset=False):
    """Everything one station's analyses need, in one dict: the gauge frame
    and snow series, the per-profile satellite arrays, the standard overpass
    table, and the station metadata (lat, lon, station_id, mask_before).
    keep_dataset=True also returns the open xarray dataset under "eo_subset"
    (needed for curtain plots); otherwise the dataset is closed."""
    eo_subset = eo.open_subset(region=station["region"])
    station_hourly, gauge_snow_mmhr = load_gauge_snow(station)
    profiles = profile_arrays(eo_subset, station["dist_col"])
    site = dict(
        station_hourly=station_hourly, gauge_snow_mmhr=gauge_snow_mmhr,
        overpasses=build_overpass_table(profiles, station_hourly,
                                        gauge_snow_mmhr),
        lat=station["lat"], lon=station["lon"],
        station_id=station["station_id"], mask_before=station["mask_before"],
        **profiles,
    )
    if keep_dataset:
        site["eo_subset"] = eo_subset
    else:
        eo_subset.close()
    return site
