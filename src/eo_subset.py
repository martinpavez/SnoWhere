"""Subset the full ACM-CAP archive (~5,259 granules, ~297 GB) down to the slice
relevant for Disko Bay surface-snowfall validation (< 2 GB).

Per granule: download (or read locally) -> keep only profiles inside the Disko
Bay box and only the 20 relevant variables, vertically cropped to below 6 km ASL
-> stage a small netCDF -> delete the original. Every granule leaves a manifest
row (id, per-station distances, profile count), so nothing is dropped silently
and an interrupted run resumes where it stopped.

Usage:
  python eo_subset.py --source local --path ../data/eo          # dry run on the 15 local files
  python eo_subset.py --probe "https://..."                     # check how the share link answers
  python eo_subset.py --auth --link "https://...sharepoint..."  # one-time device-code login
  python eo_subset.py --source sharepoint --link "https://..."  # the full run (resumable)
  python eo_subset.py --source onedrive --link "https://1drv.ms/..."   # anonymous consumer links
  python eo_subset.py --source rclone --remote "od:EC/ACM_CAP"  # last-resort fallback
  python eo_subset.py --consolidate                              # staged files -> monthly netCDF + parquet
  python eo_subset.py --source sharepoint --link "..." --limit 20      # small remote batch first

Regions: profiles are cropped to one or more named station boxes (--regions,
default nuuk,tasiilaq,ittoqqortoormiit — each a 100 km collocation radius
around its station). One download pass fills every requested region at once.
Each region set writes to its own --outdir (default data/eo_subset_stations)
because the resume manifest is keyed by granule id: reusing the Disko
manifest for different boxes would skip everything already processed there.
Monthly outputs are eo_subset_<region>_<YYYYMM>.nc; read them with
eo_reader.open_subset(dir, region=...).

The subset keeps the source's raw sentinel conventions (ice_mass_flux 0-filled
in clear air, -100 dBZ no-echo, 9.97e36 fill) so it is a faithful slice --
masking stays the job of eo_reader, same as for full granules.
"""

import argparse
import base64
import csv
import glob
import json
import os
import shutil
import sys
import tempfile
import time as _time
from datetime import timedelta

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import eo_reader as eo

# ---------------------------------------------------------------- configuration

# ---- regions: named station groups, each cropped to a box covering the
# 100 km collocation radius around its station(s) (+5 km margin).
COLLOCATION_RADIUS_KM = 100.0
BOX_MARGIN_KM = 5.0

def station_box(stations):
    """Bounding box (lat0, lat1, lon0, lon1) covering COLLOCATION_RADIUS_KM
    around every station in the dict."""
    import math
    lats = [c[0] for c in stations.values()]
    lons = [c[1] for c in stations.values()]
    r = COLLOCATION_RADIUS_KM + BOX_MARGIN_KM
    dlat = r / 111.32
    dlon = max(r / (111.32 * math.cos(math.radians(la))) for la in lats)
    return (min(lats) - dlat, max(lats) + dlat,
            min(lons) - dlon, max(lons) + dlon)

REGIONS = {
    # the three stations found to actually deliver hourly precip in 2025
    # (DMI API probe 2026-08-27); Station Nord has daily manual precip only
    "nuuk": {"stations": {"04250 Nuuk": (64.1833, -51.7308)}},
    "tasiilaq": {"stations": {"04360 Tasiilaq": (65.6111, -37.6367)}},
    "ittoqqortoormiit": {"stations": {"04339 Ittoqqortoormiit": (70.4844, -21.9511)}},
    # the original Disko Bay region; its generous hand-set box is kept verbatim
    # so re-runs stay byte-compatible with the existing data/eo_subset output
    "disko": {"stations": None,          # filled from eo.DMI_STATIONS below
              "box": (66.0, 72.0, -60.0, -44.0)},
}
REGIONS["disko"]["stations"] = eo.DMI_STATIONS
for _r in REGIONS.values():
    if "box" not in _r:
        _r["box"] = station_box(_r["stations"])

# default = the NEW station regions; the Disko box is already being filled by
# the running Jan-Jun sweep in data/eo_subset (its own manifest)
DEFAULT_REGIONS = ["nuuk", "tasiilaq", "ittoqqortoormiit"]
ACTIVE_REGIONS = list(DEFAULT_REGIONS)

# Vertical crop: keep the lowest N_LEVELS of the 242-level grid (top-down
# storage, so this is the slice [:, -N_LEVELS:]). 64 levels at ~103 m spacing
# reach ~6.1 km ASL -- the full snowfall column over the coastal stations with
# room for the melting layer and low-level cloud context. A fixed level COUNT,
# not a height threshold: the height grid drifts ~24 m between granules, so
# thresholding on height yields 63 levels in one file and 64 in another and
# the monthly concat fails.
N_LEVELS = 64

# 2-D variables to keep, cropped vertically. ice_extinction_kernel_sum is the
# obs-vs-prior QC proxy (ice_mass_flux is a derived variable and has no kernel).
VARS_2D = [
    "ice_mass_flux", "ice_mass_flux_error", "ice_water_content",
    "ice_riming_factor", "rain_rate", "rain_classification",
    "height", "temperature", "CPR_reflectivity_factor",
    "CPR_doppler_velocity", "liquid_water_content",
    "ice_extinction_kernel_sum",
]
# 1-D per-profile variables to keep.
VARS_1D = [
    "latitude", "longitude", "time", "elevation", "ice_water_path",
    "quality_status", "convergence_status", "synergy_status",
]

# Output dir. Each region-set gets its OWN dir (own manifest + staging):
# resume-done is keyed by granule id, so reusing data/eo_subset's manifest for
# a different box set would silently skip everything the Disko sweep already
# processed. The auth token alone is shared across dirs (TOKEN_CACHE below).
DEFAULT_OUT_DIR = os.path.join(os.path.dirname(eo.EO_DIR), "eo_subset_stations")
DISKO_OUT_DIR = os.path.join(os.path.dirname(eo.EO_DIR), "eo_subset")
OUT_DIR = DEFAULT_OUT_DIR
STAGING_DIR = os.path.join(OUT_DIR, "staging")
MANIFEST = os.path.join(OUT_DIR, "manifest.csv")


def set_outdir(path):
    global OUT_DIR, STAGING_DIR, MANIFEST
    OUT_DIR = path
    STAGING_DIR = os.path.join(OUT_DIR, "staging")
    MANIFEST = os.path.join(OUT_DIR, "manifest.csv")


def active_stations():
    """All stations of the active regions, merged."""
    st = {}
    for r in ACTIVE_REGIONS:
        st.update(REGIONS[r]["stations"])
    return st


def all_stations():
    """Every station of every known region, merged."""
    st = {}
    for r in REGIONS.values():
        st.update(r["stations"])
    return st


def manifest_fields():
    """REGION-INDEPENDENT superset schema: every row carries a count column
    for every known region (empty when that region was not active in the run
    that wrote the row) and the track's min distance to every station. This
    lets any --regions combination append to the same manifest — the earlier
    per-region schemas made combined single-download passes impossible."""
    return (["granule", "sensing_start", "status"]
            + [f"n_{r}" for r in REGIONS]
            + [f"min_km_{sid.split()[0]}" for sid in all_stations()])

# Frames that can reach 66N. A/E-H frames live south of 22.5N or in the SH.
NORTH_FRAMES = set("BCD")


# ------------------------------------------------------------------- manifest

def load_manifest():
    """Granules that need no re-processing. Error rows are deliberately left
    out so a resume retries them; the retry appends a second row and
    consolidate() keeps the last row per granule."""
    done = {}
    if os.path.exists(MANIFEST):
        with open(MANIFEST, newline="") as f:
            for row in csv.DictReader(f):
                if not row["status"].startswith("error"):
                    done[row["granule"]] = row["status"]
    return done


def append_manifest(row):
    new = not os.path.exists(MANIFEST)
    if not new:
        with open(MANIFEST, newline="") as f:
            header = next(csv.reader(f))
        if header != manifest_fields():
            raise SystemExit(
                f"manifest {MANIFEST} was written for a different region set "
                f"({header[3:]} vs {manifest_fields()[3:]}). Use a fresh "
                "--outdir for this region set instead of mixing manifests.")
    with open(MANIFEST, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=manifest_fields())
        if new:
            w.writeheader()
        w.writerow(row)


# ----------------------------------------------------------------- subsetting

def subset_granule(h5_path):
    """Subset one granule against every active region. Returns the manifest
    row; stages one netCDF per region the track crosses."""
    gid = eo.granule_id(h5_path)
    row = {"granule": gid, "status": "empty"}
    for r in REGIONS:
        row[f"n_{r}"] = 0 if r in ACTIVE_REGIONS else ""

    with h5py.File(h5_path, "r") as f:
        sd = f["ScienceData"]
        lat = sd["latitude"][:]
        lon = sd["longitude"][:]
        t0 = float(sd["time"][0])
        row["sensing_start"] = (
            eo.EPOCH + timedelta(seconds=t0)).strftime("%Y-%m-%d %H:%M:%S")

        for name, (slat, slon) in all_stations().items():
            d = eo.great_circle_km(slat, slon, lat, lon)
            row[f"min_km_{name.split()[0]}"] = round(float(d.min()), 1)

        sels = {}
        for r in ACTIVE_REGIONS:
            la0, la1, lo0, lo1 = REGIONS[r]["box"]
            sel = ((lat >= la0) & (lat <= la1)
                   & (lon >= lo0) & (lon <= lo1))
            n = int(sel.sum())
            row[f"n_{r}"] = n
            if n:
                sels[r] = sel
        if not sels:
            return row
        row["status"] = "ok"

        union = np.zeros(len(lat), dtype=bool)
        for sel in sels.values():
            union |= sel
        data_1d = {v: sd[v][:][union] for v in VARS_1D if v in sd}
        data_2d = {v: sd[v][:][union][:, -N_LEVELS:] for v in VARS_2D if v in sd}
        attrs = {}
        fills = {}
        for v in list(data_1d) + list(data_2d):
            a = sd[v].attrs
            attrs[v] = {k: a[k] for k in ("units", "long_name") if k in a}
            if "_FillValue" in a:      # must ride in the netCDF encoding, not
                fills[v] = a["_FillValue"].item()  # attrs, or xarray drops it

    # slice the union arrays back down to each region's profiles
    uidx = np.where(union)[0]
    for r, sel in sels.items():
        keep = np.isin(uidx, np.where(sel)[0])
        d1 = {v: a[keep] for v, a in data_1d.items()}
        d2 = {v: a[keep] for v, a in data_2d.items()}
        write_staged(gid, row["sensing_start"], r, REGIONS[r]["stations"],
                     d1, d2, attrs, fills, lat[sel], lon[sel])
    return row


def write_staged(gid, sensing_start, region, stations, data_1d, data_2d,
                 attrs, fills, lat, lon):
    """Write one granule's subset for one region as a small netCDF."""
    import xarray as xr

    n = len(lat)
    nlev = next(iter(data_2d.values())).shape[1]
    ds = xr.Dataset()
    for v, a in data_1d.items():
        ds[v] = ("profile", a)
    for v, a in data_2d.items():
        ds[v] = (("profile", "level"), a)
    for v in ds.data_vars:
        ds[v].attrs.update({k: (val.decode() if isinstance(val, bytes) else val)
                            for k, val in attrs.get(v, {}).items()})

    # derived convenience columns
    g = {"ice_mass_flux": mask_fill(data_2d["ice_mass_flux"]),
         "height_agl": (mask_fill(data_2d["height"])
                        - mask_fill(data_1d["elevation"])[:, None])}
    rate, gate = eo.surface_snowfall_rate(g)
    ds["snowfall_mmhr_1000m"] = ("profile", rate)
    ds["snowfall_mmhr_1000m"].attrs["long_name"] = (
        f"ice_mass_flux at ~{eo.REFERENCE_AGL_M:.0f} m AGL, mm/h water eq.")
    ds["snowfall_gate_agl_m"] = ("profile", gate)
    for name, (slat, slon) in stations.items():
        sid = name.split()[0]
        ds[f"distance_km_{sid}"] = (
            "profile", eo.great_circle_km(slat, slon, lat, lon))

    box = REGIONS[region]["box"]
    ds["granule"] = ("profile", np.full(n, gid, dtype="U6"))
    ds.attrs.update({
        "source_granule": gid,
        "sensing_start": sensing_start,
        "region": region,
        "box_lat": [box[0], box[1]], "box_lon": [box[2], box[3]],
        "n_levels": nlev,
        "note": ("Raw ACM-CAP conventions preserved: ice_mass_flux 0-filled in "
                 "clear air; CPR_reflectivity_factor -100 dBZ = no echo; "
                 "_FillValue 9.96921e36. Mask via eo_reader conventions."),
    })
    enc = {v: {"zlib": True, "complevel": 4, "shuffle": True}
           for v in ds.data_vars if ds[v].dtype.kind == "f"}
    for v, fv in fills.items():
        enc.setdefault(v, {})["_FillValue"] = fv
    os.makedirs(STAGING_DIR, exist_ok=True)
    out = os.path.join(STAGING_DIR, f"{region}_{gid}.nc")
    ds.to_netcdf(out, encoding=enc)
    ds.close()


def mask_fill(a):
    a = np.asarray(a, dtype="float64")
    a = a.copy()
    a[a > eo.FILL] = np.nan
    return a


# -------------------------------------------------------------- source layers

def iter_local(path, **_):
    for p in sorted(glob.glob(os.path.join(path, "*.h5"))):
        yield os.path.basename(p), lambda p=p: p, lambda: None


def _share_token(link):
    b = base64.urlsafe_b64encode(link.encode()).decode().rstrip("=")
    return "u!" + b


def _od_children(link):
    """Yield {name, size, downloadUrl} for every file behind a share link."""
    import requests
    url = (f"https://api.onedrive.com/v1.0/shares/{_share_token(link)}"
           f"/root/children?$top=200")
    while url:
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        j = r.json()
        for it in j.get("value", []):
            if "file" in it:
                yield {"name": it["name"], "size": it.get("size", 0),
                       "url": it.get("@content.downloadUrl")
                              or it.get("@microsoft.graph.downloadUrl")}
        url = j.get("@odata.nextLink")


def iter_onedrive(link, tmp_dir, **_):
    import requests

    def fetch(url, name):
        dest = os.path.join(tmp_dir, name)
        for attempt in range(5):
            try:
                with requests.get(url, stream=True, timeout=300) as r:
                    r.raise_for_status()
                    with open(dest, "wb") as f:
                        for chunk in r.iter_content(1 << 20):
                            f.write(chunk)
                return dest
            except Exception as e:
                if attempt == 4:
                    raise
                wait = 2 ** attempt * 5
                print(f"    retry {attempt+1} in {wait}s ({e})")
                _time.sleep(wait)

    for it in _od_children(link):
        name = it["name"]
        yield (name,
               lambda it=it: fetch(it["url"], it["name"]),
               lambda name=name: _rm(os.path.join(tmp_dir, name)))


# -------------------------------------------- SharePoint / OneDrive Business

# Microsoft Office's first-party public client: permits the device-code flow
# without registering an app of our own, is provisioned in the DTU tenant, and
# is preauthorized for Graph file access (Azure CLI's app is not -- AADSTS65002;
# the Graph CLI app is not provisioned there at all -- AADSTS700016).
# /.default resolves to the app's preauthorized delegated permissions.
GRAPH = "https://graph.microsoft.com/v1.0"
GRAPH_CLIENT_ID = "d3590ed6-52b3-4102-aeff-aad2292ab01c"
GRAPH_SCOPE = "https://graph.microsoft.com/.default offline_access openid"
# pinned to the original dir (not OUT_DIR): one sign-in serves every region set
TOKEN_CACHE = os.path.join(DISKO_OUT_DIR, ".graph_token.json")


def tenant_from_link(link):
    """'https://dtudk-my.sharepoint.com/...' -> 'dtudk.onmicrosoft.com'."""
    from urllib.parse import urlparse
    host = urlparse(link).hostname or ""
    stem = host.split(".")[0]
    for suf in ("-my", "-files"):
        if stem.endswith(suf):
            stem = stem[: -len(suf)]
    return f"{stem}.onmicrosoft.com"


class GraphAuth:
    """Device-code login against the share's tenant, with silent refresh.

    The access token lives ~1 h; the cached refresh token renews it without
    user interaction, so the multi-hour full run only needs the browser once.
    """

    def __init__(self, tenant):
        self.tenant = tenant
        self.tok = None
        if os.path.exists(TOKEN_CACHE):
            with open(TOKEN_CACHE) as f:
                cached = json.load(f)
            if cached.get("tenant") == tenant:
                self.tok = cached

    def _save(self, j):
        os.makedirs(OUT_DIR, exist_ok=True)
        self.tok = {"tenant": self.tenant,
                    "access_token": j["access_token"],
                    "refresh_token": j.get(
                        "refresh_token",
                        self.tok.get("refresh_token") if self.tok else None),
                    "expires_at": _time.time() + int(j.get("expires_in", 3600))}
        with open(TOKEN_CACHE, "w") as f:
            json.dump(self.tok, f)

    def device_login(self):
        import requests
        base = f"https://login.microsoftonline.com/{self.tenant}/oauth2/v2.0"
        r = requests.post(f"{base}/devicecode",
                          data={"client_id": GRAPH_CLIENT_ID,
                                "scope": GRAPH_SCOPE}, timeout=30)
        r.raise_for_status()
        dc = r.json()
        print("\n" + "=" * 62)
        print("  SIGN IN REQUIRED -- open this page in any browser:")
        print(f"    {dc['verification_uri']}")
        print(f"  and enter the code:  {dc['user_code']}")
        print("  (use the Microsoft account that can open the share link)")
        print("=" * 62 + "\n", flush=True)
        deadline = _time.time() + int(dc.get("expires_in", 900))
        while _time.time() < deadline:
            _time.sleep(int(dc.get("interval", 5)))
            r = requests.post(f"{base}/token", data={
                "grant_type": "urn:ietf:params:oauth:grant-type:device_code",
                "client_id": GRAPH_CLIENT_ID,
                "device_code": dc["device_code"]}, timeout=30)
            j = r.json()
            if "access_token" in j:
                self._save(j)
                print("signed in; token cached for later runs.")
                return
            if j.get("error") not in ("authorization_pending", "slow_down"):
                raise RuntimeError(f"device login failed: {j.get('error')}: "
                                   f"{j.get('error_description', '')[:200]}")
        raise RuntimeError("device login timed out (15 min)")

    def _refresh(self):
        import requests
        base = f"https://login.microsoftonline.com/{self.tenant}/oauth2/v2.0"
        r = requests.post(f"{base}/token", data={
            "grant_type": "refresh_token",
            "client_id": GRAPH_CLIENT_ID,
            "scope": GRAPH_SCOPE,
            "refresh_token": self.tok["refresh_token"]}, timeout=30)
        j = r.json()
        if "access_token" not in j:
            raise RuntimeError(f"token refresh failed: {j.get('error')} -- "
                               "run --auth again")
        self._save(j)

    def header(self):
        if self.tok is None:
            self.device_login()
        elif self.tok["expires_at"] - _time.time() < 300:
            self._refresh()
        return {"Authorization": f"Bearer {self.tok['access_token']}"}


def _graph_get(auth, url, stream=False):
    import requests
    for attempt in range(8):
        try:
            r = requests.get(url, headers=auth.header(), stream=stream,
                             timeout=300)
        except (requests.ConnectionError, requests.Timeout) as exc:
            # transient network loss (DNS failure, sleep/wake, wifi drop):
            # back off up to ~10 min total rather than killing a long run
            if attempt == 7:
                raise
            wait = min(2 ** attempt * 10, 300)
            print(f"    network error ({type(exc).__name__}), "
                  f"retrying in {wait}s", flush=True)
            _time.sleep(wait)
            continue
        if r.status_code == 401 and attempt == 0:
            auth._refresh()
            continue
        if r.status_code in (429, 503, 504):
            wait = int(r.headers.get("Retry-After", 2 ** attempt * 5))
            print(f"    throttled ({r.status_code}), waiting {wait}s",
                  flush=True)
            _time.sleep(wait)
            continue
        r.raise_for_status()
        return r
    raise RuntimeError(f"giving up on {url[:80]} after retries")


def _resolve_share(auth, link):
    r = _graph_get(auth, f"{GRAPH}/shares/{_share_token(link)}/driveItem")
    it = r.json()
    return it["parentReference"]["driveId"], it["id"], it["name"]


def iter_sharepoint(link, tmp_dir, tenant=None, **_):
    auth = GraphAuth(tenant or tenant_from_link(link))
    drive, item, name = _resolve_share(auth, link)
    print(f"share resolved: folder '{name}' on drive {drive[:20]}...")

    def fetch(item_id, fname):
        dest = os.path.join(tmp_dir, fname)
        # /content 302-redirects to a pre-authenticated download URL
        r = _graph_get(auth,
                       f"{GRAPH}/drives/{drive}/items/{item_id}/content",
                       stream=True)
        with open(dest, "wb") as f:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk)
        return dest

    url = (f"{GRAPH}/drives/{drive}/items/{item}/children"
           f"?$top=200&$select=id,name,size,file")
    while url:
        j = _graph_get(auth, url).json()
        for it in j.get("value", []):
            if "file" not in it:
                continue
            nm, iid = it["name"], it["id"]
            yield (nm,
                   lambda iid=iid, nm=nm: fetch(iid, nm),
                   lambda nm=nm: _rm(os.path.join(tmp_dir, nm)))
        url = j.get("@odata.nextLink")


def iter_rclone(remote, tmp_dir, **_):
    import json
    import subprocess

    out = subprocess.run(["rclone", "lsjson", remote, "--files-only"],
                         capture_output=True, text=True, check=True)
    items = json.loads(out.stdout)

    def fetch(name):
        subprocess.run(["rclone", "copy", f"{remote}/{name}", tmp_dir,
                        "--retries", "5"], check=True)
        return os.path.join(tmp_dir, name)

    for it in sorted(items, key=lambda x: x["Name"]):
        name = it["Name"]
        yield (name,
               lambda name=name: fetch(name),
               lambda name=name: _rm(os.path.join(tmp_dir, name)))


def _rm(p):
    try:
        os.remove(p)
    except OSError:
        pass


# ------------------------------------------------------------------ main loop

def relevant(name, month=None):
    if "_ACM_CAP_2B_" not in name or not name.endswith(".h5"):
        return False
    if month and f"_2B_{month}" not in name:
        return False   # sensing-start timestamp follows _2B_ in the filename
    m = eo.GRANULE_RE.search(name)
    return bool(m) and m.group(2) in NORTH_FRAMES


def run(source_iter, keep_original=False, limit=None, prefetch=3, month=None):
    """Process the source, downloading up to `prefetch` granules concurrently.

    The run is download-bound (~18 s/granule single-stream from SharePoint vs
    ~0.2 s to subset), so overlapping downloads is the whole ballgame. Results
    are processed in submission order to keep the manifest deterministic.
    """
    from collections import deque
    from concurrent.futures import ThreadPoolExecutor

    os.makedirs(OUT_DIR, exist_ok=True)
    done = load_manifest()
    n_proc = n_skip = n_hit = 0
    t0 = _time.time()
    pending = deque()

    def drain_one():
        nonlocal n_proc, n_hit
        gid, fut, cleanup = pending.popleft()
        n_proc += 1
        try:
            path = fut.result()
            row = subset_granule(path)
        except Exception as exc:
            row = {"granule": gid, "sensing_start": "",
                   "status": f"error:{type(exc).__name__}:{exc}"[:120]}
            for r in REGIONS:
                row[f"n_{r}"] = ""
            for sid in all_stations():
                row[f"min_km_{sid.split()[0]}"] = ""
        finally:
            if not keep_original:
                cleanup()
        append_manifest(row)
        if row["status"] == "ok":
            n_hit += 1
        if n_proc % 25 == 0 or row["status"] == "ok":
            el = _time.time() - t0
            counts = " ".join(f"{r}={row.get(f'n_{r}', 0)}"
                              for r in ACTIVE_REGIONS)
            print(f"  [{n_proc:>5}] {gid}  {row['status']:<8} {counts}  "
                  f"({el/max(n_proc,1):.1f} s/granule)", flush=True)

    with ThreadPoolExecutor(max_workers=prefetch) as ex:
        for name, fetch, cleanup in source_iter:
            if limit and n_proc + len(pending) >= limit:
                break
            if not relevant(name, month):
                continue
            gid = eo.granule_id(name)
            if gid in done:
                n_skip += 1
                continue
            pending.append((gid, ex.submit(fetch), cleanup))
            while len(pending) >= prefetch:
                drain_one()
        while pending:
            drain_one()

    print(f"\ndone: {n_proc} processed, {n_skip} already in manifest, "
          f"{n_hit} granules with box profiles", flush=True)


# --------------------------------------------------------------- consolidate

def consolidate():
    import pandas as pd
    import xarray as xr

    staged = sorted(glob.glob(os.path.join(STAGING_DIR, "*.nc")))
    print(f"consolidating {len(staged)} staged granules")
    groups = {}
    for p in staged:
        with xr.open_dataset(p) as ds:
            month = ds.attrs["sensing_start"][:7].replace("-", "")
            region = ds.attrs.get("region")   # absent in pre-multi-box files
        groups.setdefault((region, month), []).append(p)

    for (region, month), paths in sorted(groups.items(),
                                         key=lambda kv: (kv[0][0] or "",
                                                         kv[0][1])):
        parts = [xr.open_dataset(p) for p in sorted(paths)]
        merged = xr.concat(parts, dim="profile", combine_attrs="drop_conflicts")
        enc = {v: {"zlib": True, "complevel": 4, "shuffle": True}
               for v in merged.data_vars if merged[v].dtype.kind == "f"}
        stem = f"eo_subset_{region}_{month}" if region else f"eo_subset_{month}"
        out = os.path.join(OUT_DIR, f"{stem}.nc")
        merged.to_netcdf(out, encoding=enc)
        n = merged.sizes["profile"]
        for d in parts:
            d.close()
        merged.close()
        print(f"  {out}: {len(paths)} granules, {n} profiles, "
              f"{os.path.getsize(out)/1e6:.1f} MB")

    if os.path.exists(MANIFEST):
        df = pd.read_csv(MANIFEST)
        df = df.drop_duplicates("granule", keep="last")   # retried errors
        df.to_parquet(os.path.join(OUT_DIR, "granule_index.parquet"))
        n_err = int(df["status"].str.startswith("error").sum())
        print(f"  granule_index.parquet: {len(df)} granules, "
              f"{int((df['status']=='ok').sum())} with box profiles, "
              f"{n_err} still in error")


# --------------------------------------------------------------------- probe

def probe(link):
    print("probing share link ...")
    if "sharepoint.com" in link:
        print(f"  SharePoint / OneDrive-for-Business link "
              f"(tenant {tenant_from_link(link)}) -> needs Microsoft auth.")
        print("  1) --auth --link \"<link>\"      (one-time browser sign-in)")
        print("  2) --source sharepoint --link \"<link>\"")
        return
    try:
        n = size = 0
        first = None
        for it in _od_children(link):
            if first is None:
                first = it["name"]
            n += 1
            size += it["size"]
            if n >= 500:
                print("  (stopping enumeration preview at 500 files)")
                break
        if n == 0:
            raise RuntimeError("share answered but listed 0 files")
        print(f"  anonymous OneDrive shares API works: {n}+ files, "
              f"{size/1e9:.1f}+ GB, first: {first}")
        print("  -> use --source onedrive --link \"<link>\"")
    except Exception as e:
        print(f"  anonymous access failed ({type(e).__name__}: {e})")
        print("  -> try --auth + --source sharepoint (business links), or "
              "rclone as last resort")


# ----------------------------------------------------------------------- cli

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--source", choices=["local", "onedrive", "sharepoint",
                                         "rclone"])
    ap.add_argument("--path", help="local mode: folder containing .h5 granules")
    ap.add_argument("--link", help="onedrive/sharepoint mode: the share URL")
    ap.add_argument("--remote", help="rclone mode: remote:path")
    ap.add_argument("--tenant", help="override the tenant derived from --link")
    ap.add_argument("--auth", action="store_true",
                    help="run the one-time device-code sign-in and exit")
    ap.add_argument("--probe", metavar="LINK", help="test how a share link answers")
    ap.add_argument("--consolidate", action="store_true")
    ap.add_argument("--limit", type=int, help="stop after N new granules")
    ap.add_argument("--prefetch", type=int, default=3,
                    help="concurrent downloads (default 3)")
    ap.add_argument("--month", metavar="YYYYMM",
                    help="only granules whose sensing start is in this month")
    ap.add_argument("--regions", default=",".join(DEFAULT_REGIONS),
                    help="comma-separated region names to crop "
                         f"(known: {', '.join(REGIONS)}; "
                         f"default: {','.join(DEFAULT_REGIONS)})")
    ap.add_argument("--outdir", default=None,
                    help="output dir (manifest + staging + monthly files); "
                         "default data/eo_subset_stations, or data/eo_subset "
                         "when --regions is exactly 'disko'")
    args = ap.parse_args()

    global ACTIVE_REGIONS
    ACTIVE_REGIONS = [r.strip() for r in args.regions.split(",") if r.strip()]
    unknown = [r for r in ACTIVE_REGIONS if r not in REGIONS]
    if unknown:
        ap.error(f"unknown region(s): {unknown}; known: {list(REGIONS)}")
    set_outdir(args.outdir or (DISKO_OUT_DIR if ACTIVE_REGIONS == ["disko"]
                               else DEFAULT_OUT_DIR))

    if args.probe:
        probe(args.probe)
        return
    if args.auth:
        if not args.link and not args.tenant:
            ap.error("--auth needs --link (or --tenant)")
        GraphAuth(args.tenant or tenant_from_link(args.link)).device_login()
        return
    if args.consolidate:
        consolidate()
        return
    if not args.source:
        ap.error("need --source, --probe, --auth, or --consolidate")

    tmp = tempfile.mkdtemp(prefix="eo_subset_")
    try:
        if args.source == "local":
            if not args.path:
                ap.error("--source local needs --path")
            it = iter_local(args.path)
            run(it, keep_original=True, limit=args.limit, month=args.month)
        elif args.source == "onedrive":
            if not args.link:
                ap.error("--source onedrive needs --link")
            run(iter_onedrive(args.link, tmp), limit=args.limit)
        elif args.source == "sharepoint":
            if not args.link:
                ap.error("--source sharepoint needs --link")
            run(iter_sharepoint(args.link, tmp, tenant=args.tenant),
                limit=args.limit, prefetch=args.prefetch, month=args.month)
        else:
            if not args.remote:
                ap.error("--source rclone needs --remote")
            run(iter_rclone(args.remote, tmp), limit=args.limit)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
