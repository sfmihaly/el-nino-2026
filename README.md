# ONC slope and canyon monitor

Daily-updated view of temperature, salinity, oxygen and density (σθ) at Ocean Networks Canada
moorings off southwest Vancouver Island, compared with their pre-2026 record:

| Site | ONC location | Depth |
|---|---|---|
| Folger Pinnacle | FGPPN | 25 m |
| Upper Slope | NCBC | 396 m |
| Barkley Node | BACND | 643 m |
| Barkley Hydrates | BACHY | 871 m |
| Barkley Mid-East | BACME | 895 m |
| Barkley Canyon Axis | BACAX | 983 m |
| Clayoquot Slope (Bullseye) | NC89 | 1257 m |

Context: record subsurface warmth on the slope and in Barkley Canyon in autumn
2026, during a strong El Niño.

## What it produces

`index.html` (the published page) contains:

1. a station map on GEBCO bathymetry (`map.png`), each marker coloured by that day's
   temperature anomaly;
2. a summary table: for each site, the temperature and σθ anomaly against the
   pre-2026 normal for the date (in units and in sd), whether 2026 has set a new
   tide-removed temperature record, and whether σθ has been lighter than any pre-2026
   value; the full status table is folded underneath;
3. an interactive temperature plot for all sites (zoom to 10-minute detail;
   buttons for 3 days to the full year);
4. one figure per site (`site_<SITE>.png`) with temperature, salinity, oxygen and σθ:
   the full record on the left, the last 180 days against the seasonal cycle on
   the right.

Also `overview_year.png`, `overview_recent.png` (last 10 days at 10-minute resolution)
and `status.csv`.

## Method

- **Data**: `history/<SITE>_hourly_qc.parquet` is the quality-controlled hourly
  record to 6 Nov 2026 (frozen). Hourly averages from ONC Oceans 3.0, ONC good flags
  only, isolated salinity/oxygen spikes removed, and salinity artefact episodes
  (conductivity-cell fouling) masked by hand. Data after that are fetched live as
  10-minute averages (good flags only, isolated spikes removed automatically) and
  cached in `cache/`. Live data have **no** hand masking.
- **Map**: `assets/gebco_barkley.png` is a static GEBCO WMS image of the map area
  (`map.bbox` in `sites.json`); station positions are from the ONC locations service.
- **Tides**: temperature swings by up to ~0.5 °C over a tidal cycle at some sites,
  and a 24-h block mean leaks part of the M2 and O1 tides (aliased to ~14–15-day
  periods). All daily values here therefore use a Godin filter (24-, 24- and 25-h
  running means, gaps of up to 2 h interpolated) sampled at noon. The filtered series
  ends ~35 h before the latest data.
- **Baseline** (`baseline_doy.csv`): for each site, variable and calendar day, the
  mean and standard deviation of tide-removed daily values within ±7 days, from all
  years before 2026. **Records** (`sites.json`): pre-2026 tide-removed temperature
  maximum, salinity minimum and oxygen maximum.
- **Density**: σθ is the TEOS-10 potential density anomaly (reference pressure 0 dbar),
  computed from temperature and practical salinity at each site's nominal pressure
  (`pressure_dbar` in `sites.json`). Sensor depth varies by about ±2 dbar, which
  changes σθ by less than 0.001 kg m⁻³. A lower σθ than normal at a fixed sensor
  means lighter water, i.e. density surfaces pushed down. σθ needs valid salinity,
  so it is missing whenever salinity is (e.g. Clayoquot Slope since Dec 2024).
- **Anomalies in sd** (summary table and map) compare one tide-removed day with the
  spread of pre-2026 daily values within ±7 days of that date. They are relative
  indicators, not probabilities, and are not the same as the window-mean sd used in
  the report.
- Salinity, σθ and oxygen are not corrected for offsets between instruments; changes of a
  few hundredths in salinity or oxygen across an instrument change may be offsets.
  Clayoquot Slope salinity has been flagged bad by ONC since Dec 2024 (sensor drift).
- **Hand masks for live data**: live data get only automatic spike removal. A sensor
  fault that the spike filter misses (e.g. a salinity step from conductivity-cell
  fouling) can be blanked by adding an entry to that site's `masks` list in
  `sites.json`:
  `{"vars": ["salinity"], "start": "2026-11-11T10:30:00Z", "end": null, "reason": "..."}`.
  `end: null` keeps the mask open; set an end time (or delete the entry) when the sensor
  recovers. Masked salinity also removes σθ. Active masks are listed on the web page.
  Current mask: Hydrates salinity from 11 Nov 2026 10:30 UTC (step of about −0.12 with
  no change in temperature or oxygen).

## Run it on a Mac

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
echo 'YOUR-ONC-TOKEN' > ~/.onc_token && chmod 600 ~/.onc_token
.venv/bin/python monitor.py          # outputs in output/
```

Daily at 06:00 with launchd: edit the path in
`ca.uvic.onc-slope-temperature-monitor.plist`, copy it to `~/Library/LaunchAgents/`,
then `launchctl load ~/Library/LaunchAgents/ca.uvic.onc-slope-temperature-monitor.plist`.

## Run it on GitHub (daily web page)

The workflow `.github/workflows/monitor.yml` runs every day at 13:15 UTC and on
demand, and publishes `index.html` and the figures with GitHub Pages.

One-time setup in the repository on github.com:

1. Settings → Secrets and variables → Actions → New repository secret:
   name `ONC_API_TOKEN`, value your ONC token.
2. Settings → Pages → Build and deployment → Source: **GitHub Actions**.
3. Actions tab → Daily monitor → Run workflow (first run takes ~5–10 min; later
   runs ~1–2 min).

The data cache is kept between runs with the GitHub Actions cache. If it is lost,
the next run rebuilds it from ONC (slower, same result).

## Licence and attribution

Data: Ocean Networks Canada (CC-BY 4.0). Cite as "Ocean Networks Canada Society.
Oceans 3.0 data, accessed <date>." Code licence: to be decided.
