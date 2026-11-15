"# nasa
# AresNav: Mars EVA traverse planner (prototype)

A Google Maps-style web app for planning Marswalks (EVAs) on Mars, built on real NASA/USGS science data. It helps mission planners compare candidate sites and traverses, find resources (especially buried water-ice), and see hazards along a planned route.

**Team:** [team name and members]  
**Event:** [NASA hackathon name, year]

## Quick start

1. Keep `index.html` in the same folder as the data folders (`data`, `data_zhurong`, `data_arcadia`, `data_jezero`).
2. Double-click `index.html` and open it in Chrome or Edge. No server or internet connection is needed.
3. Use the **Region** dropdown in the header to switch areas. This reloads the page and clears any planned route.

## What you can do

- **Explore:** pan, zoom, and click to inspect any location (position, elevation, slope, ice consistency, estimated temperature).
- **Layers:** elevation with shading, slope, hazard shading, craters (1 km and larger), and SWIM ice-consistency maps for three depth bands.
- **Plan marswalk:** place a Hab, add waypoints, and set a return point. The route is colored by hazard, and the sidebar reports distance, estimated EVA time, elevation gain and loss, steepest slope, hazardous distance, crater crossings, and an elevation profile.
- **Resource finder:** rank the best sites for water-ice, hydrated clays, hydrated silica, and basalt. Perchlorate is shown as a warning.
- **Temperature:** an estimated surface temperature for any inspected point, with a local-time slider.

## Rebuilding the data

The app reads a bundle (`aresnav_data.js`) produced by `prep_utopia_data.py`. It reads only the needed window of the global DEM, so you don't need to download the full 11 GB file.

```
pip install rasterio numpy
python prep_utopia_data.py --lat 46.7 --lon 117.5 --out data --swim 0-1m=SWIM2_c0_1.tif 1-5m=SWIM2_c1_5.tif 5m+=SWIM2_c_5.tif
```

Download the combined ice-consistency GeoTIFFs (0-1 m, 1-5 m, deeper than 5 m) from the SWIM 2.0 products page: https://swim.psi.edu/SWIM2Products.php. Change `--lat`, `--lon` and `--out` for other regions.

| Region folder | Center |
|---|---|
| `data` | Utopia Planitia, 46.7°N, 117.5°E |
| `data_zhurong` | Zhurong area, about 25.1°N, 109.9°E |
| `data_arcadia` | Arcadia Planitia, 47.2°N, 184.3°E |
| `data_jezero` | Jezero crater, 18.38°N, 77.58°E |

## Data sources and licenses

- **Elevation:** MOLA-HRSC blended DEM, 200 m, v2. Fergason, Hare & Laura (2018), USGS Astrogeology. MOLA is CC0; HRSC is CC BY-SA 3.0 IGO.
- **Water-ice:** SWIM 2.0 ice-consistency maps, Planetary Science Institute. https://doi.org/10.1038/s41550-020-01290-z
- **Craters:** Robbins & Hynek global crater database, served by USGS Astrogeology. https://doi.org/10.1029/2011JE003966
- **Rover detection:** hydrated sulfate/silica seen by the Zhurong rover, Liu et al. (2022), Science Advances. https://doi.org/10.1126/sciadv.abn8555 (marker position is approximate).
- **Temperature anchor:** NASA Viking lander measurements (daily variation of 35 to 50 K).

## What is real and what is estimated

| Layer | Status |
|---|---|
| Elevation, slope, craters | Real data (200 m resolution) |
| Ice consistency | Real SWIM data. Shows how consistent the data are with buried ice, not an ice amount |
| Water-ice site ranking | Derived from the SWIM data, ranked relative to each region |
| Hydrated clays, silica, basalt | **Estimates** from simple rules (crater ejecta, slope, ice). Not mineral maps |
| Perchlorate | **Warning only.** Not mapped from orbit |
| Surface temperature | **Estimate:** −63 °C baseline, −0.001 °C per meter from the region's mean elevation, and an assumed ±20 °C daily swing |
| EVA time | **Assumption:** adjustable walking speed with a slope penalty |

## Limitations

- At 200 m resolution, slopes and hazards are regional-scale and miss boulders and small scarps.
- Hazard thresholds (default caution 2°, danger 5°) are placeholders tuned for flat terrain, not official EVA limits.
- Switching region reloads the page and clears the route. Touch-screen pinch-zoom, dragging waypoints and route export are not implemented.
- This is a prototype and is not for operational use.

## Built with

Python 3, NumPy, rasterio/GDAL for data processing; HTML5 Canvas, vanilla JavaScript and CSS for the app. No frameworks or servers.

## AI assistance

Parts of the code and documentation were written with help from Claude (Anthropic). Data choices were checked against the sources above.

## License

[Choose a license for your own code. Note that the HRSC elevation data is CC BY-SA 3.0 IGO, which is share-alike.]" 
