"""Shared configuration of the SnoWhere validation study: the analysis
constants used by every comparison notebook, and the station registry.

Provenance of the choices (details in the notebooks that established them):

- SNOW_THRESHOLD_MMHR 0.01   the radar-side detection floor ("any snow at all")
- THETA_FSS 0.1              harmonized to the DMI gauge's 0.1 mm resolution
                             floor -- below it, satellite and gauge are
                             incommensurate binaries (events notebook, §5)
- THETA_ADJUSTED 0.3         sublimation / high-altitude adjustment: a rate at
                             1000 m AGL must survive the ~700 m radar blind
                             zone to reach the gauge; the accumulation ratios
                             (x2-3.7) independently predict 0.2-0.4 and the
                             theta sweep peaks at 0.3 (stations notebook, §3d).
                             Selected on this dataset -- not validated out of
                             sample, and it does not survive swapping the
                             reference to ERA5 (model notebook, §3).
- OUTLIER_MMHR 20            unphysical-rate screen (23 profiles > 20 mm/h in
                             the Disko subset, all in non-converged patches)
- SNOW_PHASE_MAX_TEMP_C 1.5  gauge precip counts as snow only when the hourly
                             air temperature is at or below this (the phase
                             proxy used since cross_source_consistency)
- MONTH_COMPLETENESS_MIN     gauge months with fewer valid hours are dropped
                             from accumulation totals
- WIND_FLOOR_MS 2            floor on the advective wind speed in T = L/v
                             windows, so calm hours do not produce day-long
                             windows
- FSS_RADII                  the neighbourhood radii swept by every FSS curve
"""

SNOW_THRESHOLD_MMHR = 0.01     # satellite "snowing" threshold at 1 km AGL (mm/h)
THETA_FSS = 0.1                # satellite threshold harmonized to the gauge floor
THETA_ADJUSTED = 0.3           # sublimation-adjusted threshold (stations §3d/3e)
OUTLIER_MMHR = 20.0            # unphysical-rate screen (mm/h)
SNOW_PHASE_MAX_TEMP_C = 1.5    # gauge precip counts as snow only at T <= this
MONTH_COMPLETENESS_MIN = 0.80  # min fraction of valid hours for a gauge month
WIND_FLOOR_MS = 2.0            # floor for advective wind speeds in T = L/v
FSS_RADII = list(range(10, 101, 10))     # neighbourhood radii L (km)
COLLOCATION_RADIUS_KM = 100.0  # overpasses are counted within this distance

# --- the four validated stations: one entry per site with everything the
#     notebooks need to load and collocate it.
#       station_id   5-digit WMO-style id (DMI station lists, EO subset files)
#       dmi_id       6-digit DMI id (cleaned gauge CSVs in results/)
#       region       EO-subset region name for eo_reader.open_subset
#                    (None = the original Disko Bay subset in data/eo_subset)
#       dist_col     per-profile distance column in the EO subset files
#       mask_before  gauge hours before this date are dropped (instrument
#                    commissioning artefacts; see dmi_sanity_checks)
STATIONS = {
    "aasiaat": dict(
        name="Aasiaat", station_id="04220", dmi_id="422000",
        lat=68.7081, lon=-52.8517, region=None,
        clean_csv="clean_422000_B.csv", dist_col="distance_km_422000",
        mask_before=None),
    "nuuk": dict(
        name="Nuuk", station_id="04250", dmi_id="425000",
        lat=64.1833, lon=-51.7308, region="nuuk",
        clean_csv="clean_425000_B.csv", dist_col="distance_km_04250",
        mask_before=None),
    "tasiilaq": dict(
        name="Tasiilaq", station_id="04360", dmi_id="436000",
        lat=65.6111, lon=-37.6367, region="tasiilaq",
        clean_csv="clean_436000_B.csv", dist_col="distance_km_04360",
        mask_before=None),
    "ittoqqortoormiit": dict(
        name="Ittoqqortoormiit", station_id="04339", dmi_id="433900",
        lat=70.4844, lon=-21.9511, region="ittoqqortoormiit",
        clean_csv="clean_433900_B.csv", dist_col="distance_km_04339",
        mask_before="2025-05-01"),
}
