"""The gauge side of the study: loading the cleaned DMI hourly records,
splitting precipitation by phase, and reading them over time windows.

Stamp convention used throughout: an hourly value stamped at time T covers the
interval (T - 1 h, T]. DMI parameter 601 and the open-meteo model series both
follow it, which is what lets gauge and model join on the timestamp exactly.
"""

import os

import numpy as np
import pandas as pd

from study_config import SNOW_PHASE_MAX_TEMP_C

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(PROJECT_DIR, "results")


def load_gauge_snow(station, results_dir=RESULTS_DIR):
    """Load one station's cleaned hourly record and derive its snow-phase
    precipitation. Returns (station_hourly, gauge_snow_mmhr): the full frame
    (columns are DMI parameter ids as strings, e.g. "601" precip, "101" temp,
    "301" wind) and the snow-phase series -- precip clipped at 0, hours before
    the station's QC mask dropped, and warm-phase hours (T > 1.5 degC) zeroed.
    """
    station_hourly = (pd.read_csv(os.path.join(results_dir, station["clean_csv"]),
                                  parse_dates=["ts"])
                      .set_index("ts"))
    gauge_precip_mmhr = station_hourly["601"].clip(lower=0)
    if station["mask_before"]:
        gauge_precip_mmhr.loc[:station["mask_before"]] = np.nan
    gauge_snow_mmhr = gauge_precip_mmhr.where(
        station_hourly["101"] <= SNOW_PHASE_MAX_TEMP_C, 0.0)
    return station_hourly, gauge_snow_mmhr


def snow_window_stats(hourly_frame, snow_mmhr, overpass_time, half_window):
    """Snow over the hour-end stamps whose covering interval (stamp-1h, stamp]
    overlaps overpass_time +/- half_window. Works for any hour-end-stamped
    series (gauge or model), so both references are read the same way.
    Returns (total mm, fraction of stamps with snow > 0, number of stamps)."""
    lo, hi = overpass_time - half_window, overpass_time + half_window
    stamps = hourly_frame.index[(hourly_frame.index > lo)
                                & (hourly_frame.index - pd.Timedelta("1h") < hi)]
    vals = snow_mmhr.loc[stamps].dropna()
    return float(vals.sum()), float((vals > 0).mean()), len(vals)


def daily_snow(hourly_snow_mmhr, min_valid_hours=None):
    """Daily snow totals (mm); when min_valid_hours is given, days with fewer
    valid hours become NaN instead of a misleading partial sum."""
    daily_sum = hourly_snow_mmhr.resample("D").sum(min_count=1)
    if min_valid_hours is not None:
        valid_hours = hourly_snow_mmhr.resample("D").count()
        daily_sum = daily_sum.where(valid_hours >= min_valid_hours)
    return daily_sum
