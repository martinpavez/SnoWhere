# SnoWhere

**At what spatial scale does EarthCARE's snowfall detection agree with the
ground?** A hackathon-week validation of EarthCARE ACM-CAP snowfall rates
against DMI's hourly gauges at four Greenland coastal stations (Aasiaat,
Nuuk, Tasiilaq, Ittoqqortoormiit) and against ERA5, over nine months of 2025
(Jan–Jun, Oct–Dec): 445 satellite overpasses, ~61,000 usable curtain
profiles, ~29,000 gauge hours.

## Headline results

- **Detection skill is scale-limited.** Scored with a Fractions Skill Score
  adapted to curtain-vs-point geometry, EarthCARE agrees usefully with a
  gauge to **~30 km**, extended to **~50 km** with a sublimation-adjusted
  threshold (θ = 0.3 mm/h) and an 850 hPa advective time window. Beyond
  ~60 km every configuration drops below useful skill.
- **A large part of that limit is the gauge, not the satellite.** Swapping
  the reference to ERA5's 25 km cell (identical overpasses and windows)
  keeps FSS at 0.67–0.70 flat out to 100 km — above its useful-skill line
  everywhere.
- **EarthCARE beats ERA5 at comparable scale.** ERA5's nearest-cell hourly
  matching against the same gauges scores FSS 0.408 — level with the
  satellite's *worst* configuration, below useful skill.
- **Accumulation: satellite ≈ 2–3.7× the gauge** over the common cold months
  (Palerme-style estimator), consistent with gauge undercatch plus the
  ~700 m radar blind zone read at 1000 m AGL.
- **No blind-zone miss found:** in every overpass where model and gauge
  snowed but the satellite (≤30 km) did not, the track was actually 33–89 km
  away — geometry, not a miss.

## Repository layout

| Path | What it is |
|---|---|
| `src/` | Shared Python modules (see below) |
| `notebooks/` | The analyses, one question per notebook |
| `data/in-situ/` | DMI hourly station CSVs + 850 hPa winds (committed) |
| `data/model/` | ERA5 hourly snowfall at the four stations (committed) |
| `data/eo_subset*/` | EarthCARE subsets (**not** committed — see Data) |
| `results/` | Cleaned gauge records, QC verdicts, derived FSS tables |
| `figures/` | All figures, incl. the `pres_*` talk slides |

### src/ modules

- `study_config.py` — every analysis constant and the station registry, with
  the provenance of each choice. Single source of truth.
- `eo_reader.py` — reads ACM-CAP granules and the consolidated subsets with
  all sentinel/masking conventions handled (`open_subset`). Read its
  docstring before touching raw EO data.
- `eo_subset.py` — subsets the ~300 GB archive on SharePoint down to the
  <1 GB station boxes (resumable CLI; run with `--help`).
- `fetch_model_snowfall.py` — ERA5 hourly snowfall via open-meteo (keyless).
- `gauges.py` — cleaned-gauge loading, snow-phase split (T ≤ 1.5 °C proxy),
  and the hour-end-stamp window logic shared by gauge and model.
- `collocation.py` — profile usability screen, per-profile arrays, and the
  standard one-row-per-overpass table (`load_site` gives a full station in
  one call).
- `validation.py` — contingency counts, the FSS family (fixed, advective and
  point-scale), and the accumulation helpers.

Notebooks import these with `sys.path.insert(0, "../src")` (they run with
`notebooks/` as the working directory).

### Notebooks

| Notebook | Question |
|---|---|
| `compare_earthcare_dmi_events` | Does EarthCARE see snow when Aasiaat sees snow? Case studies, contingency, first FSS |
| `compare_earthcare_dmi_stations` | Does it replicate at four stations? The **main-result** FSS ladder (§3e/3f) and accumulation |
| `compare_earthcare_dmi_accumulation` | From snapshots to monthly/seasonal totals (Palerme estimator, bootstrap) |
| `compare_model_earthcare_dmi` | The ERA5 leg: model vs gauge, model vs EarthCARE across scales, three-source agreement |
| `presentation` | The 5-minute-talk figures (problem transect + 5-step FSS ladder) |
| `visualize_earthcare` | Exploratory: what ACM-CAP data looks like |
| `dmi_sanity_checks`, `gem_sanity_checks`, `cross_source_consistency` | Station QC and cross-source checks that selected the usable gauges |
| `compare_earthcare_dmi_daily` | Parked: day-granularity comparison against manual gauges |

## Data

Committed: `data/in-situ/`, `data/model/`, `results/`, `figures/` (~90 MB).

**Not committed** (≈500 MB): the EarthCARE subsets `data/eo_subset/` (Disko
Bay) and `data/eo_subset_stations/` (Nuuk/Tasiilaq/Ittoqqortoormiit). Either
copy them from the shared drive, or regenerate from the archive:

```bash
python src/eo_subset.py --auth --link "<sharepoint link>"       # once
python src/eo_subset.py --source sharepoint --link "<link>" \
       --regions "nuuk,tasiilaq,ittoqqortoormiit,disko" --prefetch 6
python src/eo_subset.py --consolidate
```

(The sign-in token is cached locally and gitignored. Disko-region monthlies
are written region-prefixed and belong in `data/eo_subset/` unprefixed — see
the module docstring.)

## Environment

The environment is managed with [uv](https://docs.astral.sh/uv/):
`pyproject.toml` declares the dependencies and `uv.lock` pins the exact
versions, so everyone resolves the same environment.

```bash
# 1. install uv once (any of these)
winget install astral-sh.uv            # Windows
curl -LsSf https://astral.sh/uv/install.sh | sh   # macOS / Linux

# 2. from the repo root: create .venv and install everything
uv sync

# 3. run things inside it
uv run python src/eo_subset.py --help
```

For the notebooks, open the repo in VS Code / Jupyter and select `.venv` as
the kernel (`ipykernel` is installed in it). To add a package later:
`uv add <name>` — it updates `pyproject.toml` and the lock file together.

If you don't want to use virtual enviroments, then use Python ≥ 3.11; 
the main dependencies are `numpy pandas xarray h5py netCDF4
matplotlib cartopy requests` (+ `ipykernel`/`nbclient` for the notebooks and
`pyarrow` for the granule index) — see `pyproject.toml`. No API keys needed
for DMI or open-meteo; only the EO subsetter needs a Microsoft sign-in.

## Method references

Roberts & Lean 2008 (FSS); Mittermaier 2014 (point-obs neighbourhoods);
Kodamana & Fletcher 2021 (CloudSat point validation, 100 km / ±30 min);
Palerme et al. 2014 (accumulation estimator); Kochendorfer et al. 2017
(gauge undercatch). Full citations in the notebooks' References cells.
