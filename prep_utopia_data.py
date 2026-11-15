#!/usr/bin/env python3
"""
AresNav data prep
-----------------
Crops real Mars data to a small LOCAL METRIC GRID around a point in Utopia Planitia
and writes files the front-end can load.

Why a local metric grid?  The USGS DEM is on a lon/lat grid, so its pixels get narrower
(east-west) as you go north. At ~47 N a "200 m" pixel is only ~137 m wide. This script
resamples to square 200 m pixels (x = metres east of centre, y = metres north of centre),
so slope, distance, pan/zoom and route maths in the browser are all simple.

Sources
  DEM      USGS MOLA-HRSC blended DEM, 200 m, v2 (Fergason, Hare & Laura, 2018)
  Craters  Robbins & Hynek global crater database (>= 1 km) via the USGS OGC API
  Ice      OPTIONAL: SWIM ice-consistency GeoTIFFs you download yourself from
           https://swim.psi.edu/SWIM2Products.php  (one or several depth bands)

Usage
  pip install rasterio numpy
  python prep_utopia_data.py                       # reads the DEM straight from USGS (window only)
  python prep_utopia_data.py --dem ./Mars_HRSC_MOLA_BlendDEM_Global_200mp_v2.tif   # local copy
  python prep_utopia_data.py --lat 18.38 --lon 77.58 --out data_jezero --swim 0-1m=SWIM2_c0_1.tif 1-5m=SWIM2_c1_5.tif 5m+=SWIM2_c_5.tif
  (label=path pairs; keep paths short / free of spaces. Re-runs reuse the cached DEM window.)

Outputs (in ./data by default)
  aresnav_data.js   <- loads from a plain <script src> even when index.html is opened via file://
  terrain.json, elevation_i16.bin, swim_<label>_u8.bin, craters.json   <- same data, split up

Decoding in the browser
  const D = window.ARESNAV_DATA, m = D.meta;
  const b64 = s => Uint8Array.from(atob(s), c => c.charCodeAt(0));
  const elev = new Int16Array(b64(D.elevation).buffer);   // metres, row-major, row 0 = north, -32768 = no data
  const ice  = D.swim && Object.fromEntries(Object.entries(D.swim).map(([k, v]) => [k, b64(v)]));
  // each ice[label]: uint8 0..254 -> consistency -1..+1 (m.swim.vmin..vmax), 255 = no data. +1 = consistent with ice
  // pixel (col,row) -> metres from centre: x = (col-(m.width-1)/2)*m.res_m,  y = ((m.height-1)/2-row)*m.res_m
"""
import argparse
import base64
import json
import math
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import Window

DEM_URL = ("https://planetarymaps.usgs.gov/mosaic/Mars/HRSC_MOLA_Blend/"
           "Mars_HRSC_MOLA_BlendDEM_Global_200mp_v2.tif")
CRATER_URL = "https://astrogeology.usgs.gov/pygeoapi/collections/mars/robbinsv1/items"
R_MARS = 3396190.0          # sphere radius used by the USGS DEM (metres)
NODATA_I16 = -32768


# ----------------------------------------------------------------------------- geometry
def build_grid(lat0, lon0, width_km, height_km, res):
    """Square-pixel local grid. Row 0 is north, col 0 is west."""
    w = int(round(width_km * 1000 / res))
    h = int(round(height_km * 1000 / res))
    xs = (np.arange(w) - (w - 1) / 2.0) * res
    ys = ((h - 1) / 2.0 - np.arange(h)) * res
    X, Y = np.meshgrid(xs, ys)
    lat = lat0 + np.degrees(Y / R_MARS)
    lon = lon0 + np.degrees(X / (R_MARS * np.cos(np.radians(lat))))   # sinusoidal-style
    if np.abs(lat).max() > 80:
        sys.exit("Window is too close to a pole for this simple projection.")
    return w, h, lat, lon


def bilinear(arr, r, c):
    """Sample arr (NaN = no data) at fractional pixel-centre coords r, c. NaN outside / next to gaps."""
    h, w = arr.shape
    r0 = np.floor(r).astype(int)
    c0 = np.floor(c).astype(int)
    fr, fc = r - r0, c - c0
    inside = (r0 >= 0) & (c0 >= 0) & (r0 + 1 < h) & (c0 + 1 < w)
    r0 = np.clip(r0, 0, h - 2)
    c0 = np.clip(c0, 0, w - 2)
    v = (arr[r0, c0] * (1 - fr) * (1 - fc) + arr[r0, c0 + 1] * (1 - fr) * fc +
         arr[r0 + 1, c0] * fr * (1 - fc) + arr[r0 + 1, c0 + 1] * fr * fc)
    v[~inside] = np.nan
    return v


def infer_lon_left(src):
    """Longitude (deg) of the left edge of the DEM. USGS documents -180..180, but check."""
    b = src.bounds
    half = math.pi * R_MARS
    if src.crs is not None and src.crs.is_geographic:
        return float(b.left)
    if abs(b.left + half) < 5e4:
        return -180.0
    if abs(b.left) < 5e4:
        return 0.0
    print("  ! Could not confirm the DEM's longitude domain from its georeferencing; "
          "assuming -180..180 as documented by USGS.")
    return -180.0


# ----------------------------------------------------------------------------- DEM
def sample_dem(path, lat, lon):
    open_path = f"/vsicurl/{path}" if str(path).startswith("http") else path
    env = dict(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
               CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif,.tiff",
               GDAL_HTTP_MAX_RETRY="5", GDAL_HTTP_RETRY_DELAY="2")
    with rasterio.Env(**env), rasterio.open(open_path) as src:
        W, H = src.width, src.height
        print(f"DEM: {W} x {H} px, dtype={src.dtypes[0]}, crs={src.crs}, nodata={src.nodata}")
        print(f"     bounds={tuple(round(v, 3) for v in src.bounds)}")
        if (W, H) != (106694, 53347):
            print("  ! Size differs from the documented 106694 x 53347 grid. "
                  "Continuing, assuming it is still a global equirectangular raster.")
        lon_left = infer_lon_left(src)
        cpd, rpd = W / 360.0, H / 180.0
        lon_n = ((lon - lon_left) % 360.0) + lon_left
        pad = 3
        c_lo = max(0, int(math.floor((lon_n.min() - lon_left) * cpd)) - pad)
        c_hi = min(W, int(math.ceil((lon_n.max() - lon_left) * cpd)) + pad)
        r_lo = max(0, int(math.floor((90 - lat.max()) * rpd)) - pad)
        r_hi = min(H, int(math.ceil((90 - lat.min()) * rpd)) + pad)
        print(f"     reading window cols {c_lo}:{c_hi}, rows {r_lo}:{r_hi} "
              f"({(c_hi - c_lo) * (r_hi - r_lo) / 1e6:.2f} Mpx)")
        data = src.read(1, window=Window(c_lo, r_lo, c_hi - c_lo, r_hi - r_lo), masked=True)
    arr = data.astype("float64").filled(np.nan)
    arr[arr < -30000] = np.nan                                  # sentinel values
    r = (90 - lat) * rpd - 0.5 - r_lo
    c = (lon_n - lon_left) * cpd - 0.5 - c_lo
    return bilinear(arr, r, c)


# ----------------------------------------------------------------------------- SWIM / any lon-lat raster
def sample_swim_raster(path, lat, lon, valid=(-1.0, 1.0)):
    """Sample a raster in ANY georeferenced CRS (lon/lat degrees or a projected Mars CRS) onto our grid."""
    from rasterio.crs import CRS
    from rasterio.warp import transform
    with rasterio.open(path) as src:
        print(f"     raster {src.width} x {src.height} px, crs={src.crs}, "
              f"bounds={tuple(round(v, 3) for v in src.bounds)}, nodata={src.nodata}")
        if src.crs is None:
            sys.exit(f"{path} has no CRS. Run `gdalinfo {path}` and send me the output.")
        data = src.read(1, masked=True).astype("float64").filled(np.nan)
        inv = ~src.transform
        if src.crs.is_geographic:
            b = src.bounds
            lon_n = lon % 360.0 if (b.right > 180.5 and b.left >= -0.5) else ((lon + 180.0) % 360.0) - 180.0
            xs, ys = lon_n, lat
        else:
            mars = CRS.from_proj4(f"+proj=longlat +a={R_MARS} +b={R_MARS} +no_defs")
            lon_w = ((lon + 180.0) % 360.0) - 180.0
            xs, ys = transform(mars, src.crs, lon_w.ravel().tolist(), lat.ravel().tolist())
            xs, ys = np.array(xs).reshape(lat.shape), np.array(ys).reshape(lat.shape)
    c = inv.a * xs + inv.b * ys + inv.c - 0.5          # fractional pixel-centre coordinates
    r = inv.d * xs + inv.e * ys + inv.f - 0.5
    bad = (data < valid[0] - 0.01) | (data > valid[1] + 0.01)     # e.g. SWIM's -30 "no data" sentinel
    print(f"     {100 * bad.mean():.1f}% of raster values outside {valid} treated as no-data")
    data[bad] = np.nan
    return bilinear(data, r, c)


# ----------------------------------------------------------------------------- craters
def fetch_craters(lon0, lat0, half_w_km, half_h_km, lat_min, lat_max, lon_min, lon_max,
                  min_km, base_url=CRATER_URL, limit=1000, max_pages=60):
    q = urllib.parse.urlencode({"f": "json", "limit": limit,
                                "bbox": f"{lon_min:.5f},{lat_min:.5f},{lon_max:.5f},{lat_max:.5f}"})
    url, feats = f"{base_url}?{q}", []
    for _ in range(max_pages):
        req = urllib.request.Request(url, headers={"Accept": "application/geo+json, application/json",
                                                   "User-Agent": "aresnav-prep/0.1"})
        with urllib.request.urlopen(req, timeout=60) as resp:
            page = json.load(resp)
        feats += page.get("features", [])
        nxt = next((l.get("href") for l in page.get("links", []) if l.get("rel") == "next"), None)
        if not nxt:
            break
        url = nxt
    else:
        print("  ! Hit the page cap while fetching craters; the list may be truncated "
              "(the server may be ignoring the bbox filter).")
    out = []
    for f in feats:
        p = f.get("properties") or {}
        lon, la, d = p.get("lon_e"), p.get("lat"), p.get("diamkm")
        if lon is None or la is None:
            g = f.get("geometry") or {}
            if g.get("type") == "Point":
                lon, la = g["coordinates"][:2]
        if lon is None or la is None or d is None:
            continue
        lon, la, d = float(lon), float(la), float(d)
        if lon > 180:
            lon -= 360
        if d < min_km:
            continue
        x = R_MARS * math.radians(lon - lon0) * math.cos(math.radians(la))
        y = R_MARS * math.radians(la - lat0)
        rad = d * 500.0                               # radius in metres
        if abs(x) <= half_w_km * 1000 + rad and abs(y) <= half_h_km * 1000 + rad:
            out.append({"id": p.get("craterid", f.get("id")), "x": round(x, 1), "y": round(y, 1),
                        "diam_m": round(d * 1000, 1), "lat": round(la, 5), "lon": round(lon, 5)})
    return sorted(out, key=lambda c: -c["diam_m"])


# ----------------------------------------------------------------------------- misc
def slope_stats(elev, res):
    dzdy, dzdx = np.gradient(elev, res)
    s = np.degrees(np.arctan(np.hypot(dzdx, dzdy)))
    return {"mean_deg": float(np.nanmean(s)), "p95_deg": float(np.nanpercentile(s, 95)),
            "max_deg": float(np.nanmax(s))}


def b64(a):
    return base64.b64encode(np.ascontiguousarray(a).tobytes()).decode("ascii")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dem", default=DEM_URL, help="URL or local path of the 200 m DEM GeoTIFF")
    ap.add_argument("--lat", type=float, default=46.7, help="window centre latitude (deg N)")
    ap.add_argument("--lon", type=float, default=117.5, help="window centre longitude (deg E)")
    ap.add_argument("--width-km", type=float, default=200.0)
    ap.add_argument("--height-km", type=float, default=200.0)
    ap.add_argument("--res", type=float, default=200.0, help="output pixel size in metres")
    ap.add_argument("--swim", nargs="*", default=[], metavar="LABEL=PATH",
                    help="SWIM ice-consistency GeoTIFFs in lon/lat degrees, e.g. 0-1m=a.tif 1-5m=b.tif")
    ap.add_argument("--swim-range", default="-1,1", help=argparse.SUPPRESS)
    ap.add_argument("--refresh", action="store_true", help="ignore the cached DEM window")
    ap.add_argument("--no-craters", action="store_true")
    ap.add_argument("--crater-url", default=CRATER_URL, help=argparse.SUPPRESS)
    ap.add_argument("--min-crater-km", type=float, default=1.0)
    ap.add_argument("--out", default="data")
    a = ap.parse_args()

    lon_wrapped = ((a.lon + 180.0) % 360.0) - 180.0
    if lon_wrapped != a.lon:
        print(f"Note: --lon {a.lon} is outside -180..180; using {lon_wrapped} (same place, east longitude).")
        a.lon = lon_wrapped
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    w, h, lat, lon = build_grid(a.lat, a.lon, a.width_km, a.height_km, a.res)
    print(f"Grid: {w} x {h} px at {a.res:g} m  ({a.width_km:g} x {a.height_km:g} km) "
          f"centred on {a.lat}N, {a.lon}E")

    cache = out / f"_dem_cache_{a.lat}_{a.lon}_{a.width_km}_{a.height_km}_{a.res}.npy"
    if cache.exists() and not a.refresh:
        elev = np.load(cache)
        print(f"DEM window loaded from cache ({cache.name}); use --refresh to re-read it.")
    else:
        elev = sample_dem(a.dem, lat, lon)
        np.save(cache, elev)
    valid = np.isfinite(elev)
    if not valid.any():
        sys.exit("No valid elevation samples in the window - check the DEM path/georeferencing.")
    print(f"Elevation: {np.nanmin(elev):.0f} to {np.nanmax(elev):.0f} m, "
          f"{100 * (1 - valid.mean()):.2f}% no-data")
    elev_i16 = np.where(valid, np.round(elev), NODATA_I16).clip(-32768, 32767).astype("<i2")
    st = slope_stats(elev, a.res)
    print(f"Slope at {a.res:g} m baseline: mean {st['mean_deg']:.2f} deg, "
          f"95th pct {st['p95_deg']:.2f} deg, max {st['max_deg']:.2f} deg")

    sources = [{
        "layer": "Elevation",
        "name": "MOLA-HRSC blended DEM, 200 m, v2",
        "citation": "Fergason, R. L., Hare, T. M., & Laura, J. (2018). HRSC and MOLA Blended Digital "
                    "Elevation Model at 200m v2. Astrogeology PDS Annex, U.S. Geological Survey.",
        "url": "https://astrogeology.usgs.gov/search/map/mars_mgs_mola_mex_hrsc_blended_dem_global_200m",
        "license": "MOLA: CC0; HRSC: CC BY-SA 3.0 IGO"}]

    meta = {
        "region": "Utopia Planitia",
        "center": {"lat": a.lat, "lon": a.lon},
        "width": w, "height": h, "res_m": a.res,
        "extent_m": {"x": [-(w - 1) / 2 * a.res, (w - 1) / 2 * a.res],
                     "y": [-(h - 1) / 2 * a.res, (h - 1) / 2 * a.res]},
        "bounds_deg": {"lat_min": float(lat.min()), "lat_max": float(lat.max()),
                       "lon_min": float(lon.min()), "lon_max": float(lon.max())},
        "radius_m": R_MARS,
        "projection": "local sinusoidal: x = R*dLon*cos(lat), y = R*dLat (metres from centre; "
                      "inverse: lat = lat0 + y/R, lon = lon0 + x/(R*cos(lat)))",
        "layout": "row-major, row 0 = north, col 0 = west",
        "elevation": {"dtype": "int16-le", "units": "m", "nodata": NODATA_I16,
                      "min": float(np.nanmin(elev)), "max": float(np.nanmax(elev))},
        "slope_stats_at_res": st,
    }

    swim_layers = {}
    if a.swim:
        lo, hi = (float(v) for v in a.swim_range.split(","))
        for item in a.swim:
            label, _, path = item.partition("=")
            if not path:
                label, path = Path(item).stem, item
            label = "".join(ch if ch.isalnum() or ch in "-+_" else "_" for ch in label)
            print(f"SWIM layer '{label}':")
            sv = sample_swim_raster(path, lat, lon, (lo, hi))
            cover = 100 * np.isfinite(sv).mean()
            if cover == 0:
                print("  ! no valid values inside this window; skipping this layer.")
                continue
            q = np.where(np.isfinite(sv), np.round((np.clip(sv, lo, hi) - lo) / (hi - lo) * 254), 255)
            swim_layers[label] = q.astype("u1")
            print(f"     window coverage {cover:.1f}%, values {np.nanmin(sv):.2f} to {np.nanmax(sv):.2f}, "
                  f"mean {np.nanmean(sv):.2f}")
        if swim_layers:
            meta["swim"] = {"dtype": "uint8", "vmin": lo, "vmax": hi, "nodata": 255,
                            "layers": list(swim_layers),
                            "note": "0..254 maps linearly to vmin..vmax (ice consistency: +1 consistent with ice, "
                                    "0 no indication, -1 inconsistent)"}
            sources.append({"layer": "Subsurface water-ice", "name": "SWIM 2.0 ice-consistency maps",
                            "citation": "Mars Subsurface Water Ice Mapping (SWIM) project, Planetary Science "
                                        "Institute; see https://doi.org/10.1038/s41550-020-01290-z",
                            "url": "https://swim.psi.edu/SWIM2Products.php",
                            "license": "see SWIM project terms"})

    craters = None
    if not a.no_craters:
        try:
            craters = fetch_craters(a.lon, a.lat, a.width_km / 2, a.height_km / 2,
                                    float(lat.min()), float(lat.max()),
                                    float(lon.min()), float(lon.max()), a.min_crater_km,
                                    base_url=a.crater_url)
            print(f"Craters >= {a.min_crater_km:g} km in/overlapping window: {len(craters)}")
            sources.append({"layer": "Craters", "name": "Robbins & Hynek global Mars crater database",
                            "citation": "Robbins, S. J., & Hynek, B. M. (2012). A new global database of Mars "
                                        "impact craters >= 1 km. J. Geophys. Res. Planets, 117(E5). "
                                        "https://doi.org/10.1029/2011JE003966",
                            "url": "https://astrogeology.usgs.gov/pygeoapi/collections/mars/robbinsv1",
                            "license": "see USGS Astrogeology terms"})
        except Exception as e:                              # keep going: DEM is the important part
            print(f"  ! Crater download failed ({e}). Continuing without craters.")
    meta["sources"] = sources
    meta["crater_count"] = None if craters is None else len(craters)

    (out / "terrain.json").write_text(json.dumps(meta, indent=2))
    (out / "elevation_i16.bin").write_bytes(elev_i16.tobytes())
    for label, arr in swim_layers.items():
        (out / f"swim_{label}_u8.bin").write_bytes(arr.tobytes())
    if craters is not None:
        (out / "craters.json").write_text(json.dumps(craters))
    bundle = {"meta": meta, "elevation": b64(elev_i16),
              "swim": {k: b64(v) for k, v in swim_layers.items()} or None, "craters": craters}
    (out / "aresnav_data.js").write_text(
        "// generated by prep_utopia_data.py\nwindow.ARESNAV_DATA = " + json.dumps(bundle) + ";\n")
    print(f"\nDone. Wrote {', '.join(p.name for p in sorted(out.iterdir()))} to {out.resolve()}")


if __name__ == "__main__":
    main()
