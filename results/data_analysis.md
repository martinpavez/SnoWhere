# DMI Greenland Stations — Data Analysis

Analysis of the six DMI station CSVs in `data/`, based on the files themselves, the accompanying report **DMI Report 25-08** ("Weather Observations from Greenland 1958–2024"), and `Station_details.xlsx`. Goal: understand the data before building sanity-check scripts for a validation period.

---

## 1. Headline findings

1. **Snow depth is NOT in this dataset.** DMI Report 25-08 states explicitly that snow depth was dropped from this report series — "go to DMI Report 21-08 for snow data." These files provide the *meteorological reference variables* (temperature, precipitation, wind, humidity, pressure, radiation, cloud) against which snow variables can be validated. The snow observations themselves must come from DMI Report 21-08 or the DMI open-data API.
2. **All six stations are in the Disko Bay region** (Ilulissat, Aasiaat, Qeqertarsuaq), forming three near-neighbour pairs — very useful for cross-station consistency checks.
3. **DMI's own QC missed at least one gross failure**: station 421900 reports air temperatures of **+44 to +92 °C** during 8–24 April 2023. This alone justifies the sanity scripts.
4. **422000 Aasiaat is the anchor station**: the only one with a continuous 1958–2024 record and the only one with usable precipitation and radiation.
5. **Recommended sanity-check period: 2014–2024** (hourly cadence everywhere, hourly extremes + gusts + hourly precip exist only from 2014), with the caveats listed in §5.

## 2. Station inventory

All files are semicolon-separated CSV; timestamps are **UTC**; missing values are empty fields; values have 1-decimal precision except 201, 365, 371, 550, 801 (integers).

| File | Station (WMO) | Type | Operation | Lat | Lon | Elev (m) | Rows | Cadence |
|---|---|---|---|---|---|---|---|---|
| 421600 | Ilulissat (04216) | DMI synoptic | 1961-01 → 1992-08 | 69.2167 | −51.0500 | 39 | 49,999 | 3 h → 6 h (degrades after ~1980) |
| 421800 | Qeqertarsuaq (04218) | DMI synoptic | 1962-01 → 1980-06 | 69.2333 | −53.5167 | 24 | 49,761 | 3 h |
| 421900 | Qeqertarsuaq Heli. (04219) | MIT/SAVS airport | 2010-07 → active | 69.2514 | −53.5147 | 2.7 | 114,456 | 1 h |
| 422000 | Aasiaat (04220) | DMI/V98 automatic | 1958-01 → active | 68.7081 | −52.8517 | 43 | 355,373 | 3 h → 1 h (1997) |
| 422100 | Mitt. Ilulissat (04221) | MIT/SAVS airport | 1991-08 → active | 69.2403 | −51.0661 | 29 | 251,645 | 1 h (from 1997) |
| 422400 | Mitt. Aasiaat (04224) | MIT/SAVS airport | 2000-11 → active | 68.7219 | −52.7847 | 22.6 | 181,720 | 1 h |

**Near-neighbour pairs (for cross-validation):**

- Aasiaat 422000 ↔ Mitt. Aasiaat 422400: ~3 km apart, overlap 2000-11 → 2021-01 and 2024.
- Ilulissat 421600 ↔ Mitt. Ilulissat 422100: ~2.6 km apart, overlap only 1991-08 → 1992-08.
- Qeqertarsuaq 421800 ↔ Qeqertarsuaq Heli. 421900: ~2 km apart, **no overlap** (1980 vs 2010).

## 3. Parameter availability

Columns: `Station;Year;Month;Day;Hour(utc);101;112;113;122;123;201;301;305;365;371;401;504;550;601;603;609;801`

| Param | Meaning | 421600 | 421800 | 421900 | 422000 | 422100 | 422400 |
|---|---|---|---|---|---|---|---|
| 101 | Mean air temp (°C, 2 m) | 99.9% | 100% | 67% | 95% | 99.6% | 100% |
| 112/122 | Tmax/Tmin last hour | — | — | 2014+ | 2014+ | 2014+ | 2014+ |
| 113/123 | Tmax/Tmin last 12 h (06/18 UTC) | ✔ | ✔ | sparse | ✔ | ✔ | ✔ |
| 201 | Relative humidity (%) | 99.9% | 100% | 75% | 99.7% | 98.8% | 98.3% |
| 301 | Mean wind speed (m/s) | 98.7% | 99.9% | 85% | 99.2% | 99.6% | 99.9% |
| 305 | 3-s gust (m/s) | — | — | 2014+ | 2014+ | 2014+ | 2014+ |
| 365/371 | Wind direction (°) | ✔ | ✔ | 85% | ✔ | 94.6% | 98.7% |
| 401 | MSL pressure (hPa) | 97.2% | 100% | 82.5% | 99.6% | 99.8% | 99.6% |
| 504 | Sunshine | — (excluded at all 6; pyranometer problems) | — | — | — | — | — |
| 550 | Global radiation (W/m²) | — | — | — | **2014+** | — | — |
| 601 | Precip last 1 h (mm) | — | — | — | **2014+** | — | scraps (2024) |
| 603 | Precip last 6/12 h (mm) | 1961–1991 | 1962–1980 | — | 1958–2024 | — | scraps |
| 609 | Precip last 24 h (mm) | — | — | — | 2014+ (sparse) | — | — |
| 801 | Cloud cover (%) | 99.4% | 99.9% | — | 89% | 7.7% (mostly early years) | 2024 only |

Key takeaways:

- **Precipitation exists only at Aasiaat (422000)** for any recent period. The three SAVS airport stations have essentially none (the report notes automatic rain gauges have problems and SAVS airports don't measure it). For rain/snow partitioning or accumulation checks, Aasiaat is the only reference.
- **Radiation (550) only at Aasiaat, only from 2014** (relevant for snowmelt energy checks).
- **Cloud cover at the SAVS stations is permanently excluded** by DMI (ceilometers can't report clear sky): 422100 and 422400 (report lists Mitt. Ilulissat, Mitt. Aasiaat among the nine affected).
- Hourly extremes (112/122) and gusts (305) **exist only from 2014 onward** — a hard constraint on the validation period.

## 4. Data quality findings (verified in the files)

### 4.1 Confirmed bad data — must be caught by sanity scripts

- **421900, 8–24 April 2023**: 14 rows with T = +44…+92 °C across params 101/112/113/122/123 (e.g. 2023-04-08 20:00 → 101=71.2, 112=92.3). Physically impossible; slipped through DMI QC. A plausibility bound (T > +25 °C in Disko Bay) catches all of them.
- **421900**: 11 rows where gust (305) < mean wind (301).
- **422000**: 150 rows with wind direction = 0 (calm code) but speed > 2 m/s — direction should be treated as suspect/missing there.
- Small numbers of 12-h extreme inconsistencies (Tmax12h < T at ob time, or Tmin12h > T): ≤ 7 per station, 16 Tmin violations at 422000. Mostly rounding/window edge effects but worth flagging.

### 4.2 Encodings and sentinels the scripts must handle

- **Trace precipitation**: `-0.1` in 603 means "> 0 but < 0.1 mm" (957 occurrences at 422000, 2014–2024). Map to 0 (or an epsilon), never treat as negative rain.
- **Wind direction 0 = calm** (co-occurs with speed 0 in >99% of cases). 990 = variable direction does **not** appear (already cleaned).
- **Cloud cover 801 is already in %** in *all* eras of these files (octas were converted): pre-2014 values ∈ {0,10,25,40,50,60,75,90,100}; 2014+ values ∈ {0,10,40,60,90,100}. The octas-obscured code 9 does not appear.
- **Precip window semantics** (603): obs at 06/18 UTC cover 12 h; obs at 00/12 UTC cover 6 h and are *included in* the following 06/18 value — never simply sum the four synop values. Daily accumulation runs 06:01 → 06:00 UTC next day (report §5.2.1).
- Missing = empty field. No NaN codes like −99 were found.

### 4.3 Timestamps and gaps

Timestamps are clean: valid calendar dates, strictly sorted, **no duplicates** at any station. Each file ends at 2025-01-01 00:00 (the ob covering the last hour of 2024) or earlier.

Major gaps (> 30 days):

| Station | Gap |
|---|---|
| 422400 | **2021-01-02 → 2024-01-01 (3 years!)** |
| 421900 | 2017-06-09 → 2017-07-26 (47 d); 2018-07-29 → 2019-01-10 (165 d) |
| 422100 | 2006-04-06 → 2006-10-20 (197 d) |
| 421600 | Cadence degrades: 8/day (1960s) → 5/day (1970s) → 4/day (1980s) → 3/day (0,12,18 UTC only, 1990s); also odd observing hours (0,12,15,18) in some periods — matches the report's own example |

### 4.4 Observed physical ranges (sane bounds for the scripts)

Across all stations (excluding the 421900 April-2023 failure): T ∈ [−40.5, +22] °C, RH ∈ [10, 100] %, MSL pressure ∈ [943.9, 1057.6] hPa, wind ∈ [0, 38.6] m/s, gust ≤ 35.1 m/s, precip 12 h ≤ 76.7 mm, radiation ∈ [0, 1028] W/m², directions ∈ [0, 360].

## 5. Recommendation for the sanity-script period

**Primary window: 2014-01-01 → 2024-12-31.**

- All four active stations report hourly; hourly extremes (112/122), gusts (305), hourly precip (601) and radiation (550) exist only from 2014.
- Contains a *known-bad* episode (421900 April 2023) — a built-in positive control for the scripts.
- Caveats to encode: 422400 is silent 2021–2023; 421900 has the 2017/2018–19 outages; Aasiaat 601/550 have occasional excluded months (2016 radiation missing entirely).

**Long-record option:** 422000 Aasiaat alone supports 1958–2024 (3 h → 1 h) for climatological context; use 101/201/301/365/401/603/801 only.

**Cross-station checks:** Aasiaat pair (422000/422400, ~3 km) over 2014–2020 + 2024; Qeqertarsuaq pair has no overlap; Ilulissat pair overlaps only 1991–92.

## 6. Suggested sanity checks (blueprint for the scripts)

1. **Schema/parse**: semicolon separator, expected 22 columns, numeric coercion, station id matches filename.
2. **Timestamps**: valid calendar date/hour, strictly increasing, no duplicates, cadence conformity (expected 1 h post-2014), gap census with a report of gaps > N hours.
3. **Range (climatological plausibility)** per variable, e.g. T ∈ [−50, +25] °C, RH ∈ [5, 100], P ∈ [935, 1065] hPa, wind ∈ [0, 50], gust ∈ [0, 60], dir ∈ [0, 360], precip1h ∈ [0, 30], precip12h ∈ [0, 100], rad ∈ [0, 1100], cloud ∈ allowed value set.
4. **Sentinels**: map −0.1 precip → trace; dir 0 + speed 0 → calm; dir 0 + speed > 2 → flag.
5. **Internal consistency**: Tmin ≤ T ≤ Tmax (respecting 1 h vs 12 h windows), gust ≥ mean wind, spike test (|ΔT| > 10 °C between consecutive hourly obs), persistence test (value stuck > 24 h with variance 0), RH = 100 plateaus.
6. **Cross-station coherence**: |ΔT| between the Aasiaat pair beyond a seasonal threshold (e.g. > 8 °C) → flag; same for pressure (should agree within ~2 hPa after MSL reduction).
7. **Precip aggregation logic**: rebuild daily 06–06 UTC sums per report §5.2.1 and compare 601-derived vs 603-derived vs 609 where they coexist (2014+ at Aasiaat).
8. **Snow-relevance derived checks** (once snow data is added from Report 21-08): precipitation with T < ~+1.5 °C ⇒ snowfall expected ⇒ snow depth should not decrease; T ≫ 0 sustained ⇒ depth should not increase; wind > ~10 m/s ⇒ drifting caveat on depth changes.

---

*Generated 2026-08-24. Profiling scripts used for this analysis are reproducible from the CSVs (pandas); all numbers above were computed from the files in `data/`.*
