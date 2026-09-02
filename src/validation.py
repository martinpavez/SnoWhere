"""Scores used across the comparison notebooks: 2x2 contingency counts and
the Fractions Skill Score family (Roberts & Lean 2008; point-observation
flavour after Mittermaier 2014), plus the small accumulation helpers.

The FSS convention of this project: per overpass, p_test is the fraction of
usable satellite profiles within the neighbourhood radius whose rate exceeds
theta, and p_ref is the fraction of hour-end stamps in the time window whose
snow is > 0 (gauge or model, read by the same function). FSS is then
1 - sum((p_ref - p_test)^2) / (sum(p_test^2) + sum(p_ref^2)) over overpasses,
read against the useful-skill yardstick 0.5 + f0/2.
"""

import numpy as np
import pandas as pd

from gauges import snow_window_stats
from study_config import WIND_FLOOR_MS


def contingency(test_snowing, reference_snowing):
    """2x2 contingency counts for two boolean arrays (test says snow,
    reference says snow) plus POD = hit/(hit+miss) and FAR = FA/(hit+FA)."""
    hit = int((test_snowing & reference_snowing).sum())
    miss = int((~test_snowing & reference_snowing).sum())
    false_alarm = int((test_snowing & ~reference_snowing).sum())
    correct_neg = int((~test_snowing & ~reference_snowing).sum())
    pod = hit / (hit + miss) if hit + miss else np.nan
    far = false_alarm / (hit + false_alarm) if hit + false_alarm else np.nan
    return dict(hit=hit, miss=miss, false_alarm=false_alarm,
                correct_neg=correct_neg,
                POD=round(pod, 2) if pod == pod else None,
                FAR=round(far, 2) if far == far else None)


def fss_score(pairs):
    """Fractions Skill Score of an (n, 2) array of (p_test, p_reference)
    pairs: 1 - sum((p_ref - p_test)^2) / (sum(p_test^2) + sum(p_ref^2));
    NaN when there are no pairs or both sides are all-zero."""
    if len(pairs) == 0:
        return np.nan
    denom = (pairs[:, 0] ** 2).sum() + (pairs[:, 1] ** 2).sum()
    return (1 - ((pairs[:, 1] - pairs[:, 0]) ** 2).sum() / denom
            if denom else np.nan)


def point_fss(test_snowing, reference_snowing):
    """FSS of paired 0/1 flags -- the point-scale limit, where there is no
    neighbourhood to grow. With binary fractions it reduces to
    2*hit / (2*hit + miss + false_alarm), the Dice/F1 score."""
    return fss_score(np.column_stack([np.asarray(test_snowing, dtype=float),
                                      np.asarray(reference_snowing,
                                                 dtype=float)]))


def useful_fss(base_rate):
    """The useful-skill yardstick 0.5 + f0/2 for a reference whose base rate
    (mean p_ref over the scored pairs) is base_rate."""
    return 0.5 + base_rate / 2


def advective_window(wind_series, radius_km, speed_floor=WIND_FLOOR_MS):
    """A half-window callable for fss_pairs implementing T = L/v: at each
    overpass, look up the wind speed (m/s, hour-stamped series) at the
    overpass hour, floor it at speed_floor, and return the Timedelta the
    neighbourhood takes to advect past the station."""
    def half_window(overpass_time):
        speed = wind_series.get(overpass_time.floor("h"), np.nan)
        speed = (max(float(speed), speed_floor)
                 if np.isfinite(speed) else speed_floor)
        return pd.Timedelta(hours=radius_km / (speed * 3.6))
    return half_window


def fss_pairs(site, radius_km, theta, half_window=pd.Timedelta("1h"),
              reference_frame=None, reference_snow=None):
    """One (p_test, p_ref) pair per overpass that has usable profiles within
    radius_km. p_test is the fraction of those profiles above theta; p_ref is
    the reference's snowing fraction over overpass +/- half_window, where
    half_window is a Timedelta or a callable(overpass_time) -> Timedelta
    (see advective_window). The reference defaults to the site's gauge; pass
    reference_frame/reference_snow (e.g. a model series) to swap it -- only
    the reference changes, never the overpasses, profiles or windows."""
    if reference_frame is None:
        reference_frame = site["station_hourly"]
        reference_snow = site["gauge_snow_mmhr"]
    pairs = []
    for _, overpass in site["overpasses"].iterrows():
        in_neighbourhood = (site["profile_usable"]
                            & (site["profile_granule"] == overpass["granule"])
                            & (site["profile_distance_km"] <= radius_km))
        if not in_neighbourhood.any():
            continue
        window = (half_window(overpass["time"]) if callable(half_window)
                  else half_window)
        p_ref = snow_window_stats(reference_frame, reference_snow,
                                  overpass["time"], window)[1]
        if np.isfinite(p_ref):
            p_test = float((site["profile_snowfall_mmhr"][in_neighbourhood]
                            > theta).mean())
            pairs.append((p_test, p_ref))
    return np.array(pairs)


def hours_in_month(month_str):
    """Number of hours in a 'YYYY-MM' month (the integration time T of the
    Palerme accumulation estimator)."""
    return pd.Period(month_str, "M").days_in_month * 24


def bootstrap_ci(rates, hours, rng, n_boot=10_000, percentiles=(10, 90)):
    """Percentile band of the accumulation (mean rate x hours) when the
    per-overpass rates are resampled with replacement n_boot times. Pass the
    notebook's own RNG so results are reproducible."""
    rates = np.asarray(rates)
    if len(rates) == 0:
        return np.nan, np.nan
    boot_idx = rng.integers(0, len(rates), size=(n_boot, len(rates)))
    boot_totals = rates[boot_idx].mean(axis=1) * hours
    lo, hi = np.percentile(boot_totals, list(percentiles))
    return lo, hi
