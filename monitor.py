#!/usr/bin/env python3
"""Near-real-time monitor of temperature, salinity and oxygen at ONC moorings
on the Vancouver Island shelf, continental slope and in Barkley Canyon.

Data
  history/<SITE>_hourly_qc.parquet  frozen, quality-controlled hourly record to 6 Nov 2026
  cache/<SITE>_10min.parquet        live 10-minute data from ONC since then (rebuilt if missing)

Tides
  Temperature at these sites swings by up to ~0.5 °C over a tidal cycle, and a
  24-h block mean does not remove the tide cleanly (aliasing). All "daily" values
  here are therefore Godin-filtered (24-h, 24-h, 25-h running means in sequence)
  and taken at noon each day. The filtered series ends ~35 h before the latest
  data; the most recent hours are shown unfiltered.

Density
  σθ (TEOS-10 potential density anomaly, reference pressure 0) from temperature and
  practical salinity at each site's nominal pressure (sites.json). The sensor depth
  varies by about ±2 dbar, which changes σθ by less than 0.001 kg m⁻³.

Outputs (output/)
  index.html                    overview page: status table, interactive temperature
                                plot (zoom to 10-min detail) and the per-site figures
  overview_year.png             temperature, this year, all sites (report Figure 3 layout)
  overview_recent.png           last N days at 10-minute resolution
  map.png                       station map, markers coloured by today's temperature anomaly
  site_<SITE>.png               temperature, salinity, oxygen, σθ: full record + last 180 days
  status.csv                    latest values and anomalies

Token: environment variable ONC_API_TOKEN or the file ~/.onc_token.
Usage: python monitor.py [--recent-days 10] [--hires-days 60]
"""
import argparse, json, os, sys, time
from datetime import datetime, timedelta, timezone
from pathlib import Path

import gsw
import numpy as np
import pandas as pd
import requests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ONC_API = "https://data.oceannetworks.ca/api"
HERE = Path(__file__).resolve().parent
RESAMPLE_S = 600
GOOD = {7}                         # ONC flag: valid resampled average
EV, BAND, CLIM, GREY = "#c0392b", "#cfd8e3", "#5b7a99", "0.6"
EN_COL, EN_TXT, BLOB_COL = "#f4a261", "#b5651d", "#7b6fb0"   # event shading
VARS = [("temperature", "T", "Temperature (°C)", "max"),
        ("salinity", "S", "Salinity (PSS-78)", "min"),
        ("O2", "O2", "Oxygen (ml/l)", "max"),
        ("sigma0", "SIG", "σθ (kg m⁻³)", "min")]
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "xtick.direction": "out", "ytick.direction": "out"})


# ---------------------------------------------------------------- data access
def token():
    t = os.environ.get("ONC_API_TOKEN")
    if not t:
        f = Path.home() / ".onc_token"
        t = f.read_text().strip() if f.exists() else None
    if not t:
        sys.exit("No ONC token: set ONC_API_TOKEN or create ~/.onc_token")
    return t


def onc_get(path, tok, **params):
    params["token"] = tok
    err = ""
    for attempt in range(3):
        try:
            r = requests.get(f"{ONC_API}/{path}", params=params, timeout=300)
            if r.status_code == 404:
                return None
            if r.status_code == 400:          # e.g. no device deployed in the window
                return None
            if r.status_code == 200:
                return r.json()
            err = f"HTTP {r.status_code}: {r.text[:200]}"
        except requests.RequestException as e:
            err = str(e)
        time.sleep(20 * (attempt + 1))
    raise RuntimeError(f"ONC {path} failed: {err}")


def fetch(code, cat, sensors, t0, t1, tok):
    """10-minute averages of the given sensors, ONC good flags only, columns renamed."""
    names = {"temperature": "T", "salinity": "S", "oxygen_corrected": "O2"}
    parts = []
    a = t0
    while a < t1:
        b = min(a + timedelta(days=30), t1)
        r = onc_get("scalardata/location", tok, locationCode=code, deviceCategoryCode=cat,
                    dateFrom=a.strftime("%Y-%m-%dT%H:%M:%S.000Z"), dateTo=b.strftime("%Y-%m-%dT%H:%M:%S.000Z"),
                    resamplePeriod=RESAMPLE_S, resampleType="avg", rowLimit=100000,
                    sensorCategoryCodes=",".join(sensors))
        for s in (r or {}).get("sensorData") or []:
            if s["sensorCategoryCode"] in sensors:
                d = s["data"]
                x = pd.Series(d["values"], index=pd.to_datetime(d["sampleTimes"], utc=True), dtype=float)
                x[~pd.Series(d["qaqcFlags"], index=x.index).isin(GOOD)] = np.nan
                parts.append(x.rename(names[s["sensorCategoryCode"]]))
        a = b
    if not parts:
        return pd.DataFrame(columns=list({names[s] for s in sensors}))
    df = pd.concat(parts, axis=1)
    df = df.T.groupby(level=0).first().T if df.columns.duplicated().any() else df
    return df[~df.index.duplicated()].sort_index()


def despike(s, floor, k=8):
    """Remove isolated spikes: |x - 5-h median| > k*MAD and > floor."""
    x = s.dropna()
    if len(x) < 50:
        return s
    med = x.rolling("5h", center=True, min_periods=3).median()
    res = x - med
    mad = 1.4826 * res.abs().rolling("30D", center=True, min_periods=48).median()
    bad = x.index[(res.abs() > k * mad) & (res.abs() > floor)]
    s = s.copy(); s.loc[bad] = np.nan
    return s


def update_cache(site, tok, hist_end, hires_days, cache_dir):
    f = cache_dir / f"{site['code']}_10min.parquet"
    now = datetime.now(timezone.utc)
    old = pd.read_parquet(f) if f.exists() else None
    if old is not None and len(old):
        start = old.index.max().to_pydatetime() - timedelta(days=1)
    else:   # rebuild: everything after the frozen history, and the high-resolution window
        start = min(hist_end.to_pydatetime(), now - timedelta(days=hires_days)) - timedelta(days=2)
    ctd = fetch(site["code"], "CTD", ["temperature", "salinity"], start, now, tok)
    oxy = fetch(site["code"], "OXYSENSOR", ["oxygen_corrected"], start, now, tok)
    new = pd.concat([ctd, oxy], axis=1) if len(oxy) else ctd
    for c in ["T", "S", "O2"]:
        if c not in new:
            new[c] = np.nan
    new = new[["T", "S", "O2"]]
    if len(new):
        new["S"] = despike(new["S"], site["despike_floor_salinity"])
        new["O2"] = despike(new["O2"], site["despike_floor_O2"])
    df = new if old is None else pd.concat([old[old.index < start], new])
    df = df[~df.index.duplicated(keep="last")].sort_index()
    df.to_parquet(f)
    return df


# ---------------------------------------------------------------- tidal filter
def godin(hourly, max_gap_h=2):
    """Godin filter (A24 A24 A25) of an hourly series; gaps <= max_gap_h are interpolated."""
    s = hourly.copy()
    s.index = s.index.floor("h")
    s = s[~s.index.duplicated()].asfreq("1h")
    si = s.interpolate(limit=max_gap_h, limit_area="inside")
    return (si.rolling(24, center=True, min_periods=24).mean()
              .rolling(24, center=True, min_periods=24).mean()
              .rolling(25, center=True, min_periods=25).mean())


def noon(g):
    d = g[g.index.hour == 12]
    d.index = d.index.normalize()
    return d.dropna()


def baseline(base, code, var, index):
    b = base[(base.code == code) & (base["var"] == var)].set_index("doy")
    doy = index.dayofyear
    return b["mean"].reindex(doy).to_numpy(), b["sd"].reindex(doy).to_numpy()


def apply_masks(h, s):
    """Blank hand-identified bad intervals listed under the site's "masks" in sites.json.
    end = null means the mask is still open (applies to all later data)."""
    for m in s.get("masks", []):
        t0 = pd.Timestamp(m["start"])
        t1 = pd.Timestamp(m["end"]) if m.get("end") else h.index.max()
        for v in m["vars"]:
            if v in h:
                h.loc[(h.index >= t0) & (h.index <= t1), v] = np.nan
    return h


def add_sigma0(h, s):
    """σθ (kg m⁻³) from temperature and practical salinity at the site's nominal pressure."""
    p = s["pressure_dbar"]
    SA = gsw.SA_from_SP(h["salinity"].to_numpy(float), p, s["lon"], s["lat"])
    CT = gsw.CT_from_t(SA, h["temperature"].to_numpy(float), p)
    h = h.copy(); h["sigma0"] = gsw.sigma0(SA, CT)
    return h


def draw_events(axs, events, end):
    """Shade El Niño (filled) and Blob (hatched) periods on every panel; label them above the top row.
    axs is the (rows, 2) array of site-figure axes; the x-limits are left as they are."""
    import matplotlib.transforms as mtransforms
    if not events:
        return
    ev = []
    for e in events:
        a = pd.Timestamp(e["start"], tz="UTC")
        b = pd.Timestamp(e["end"], tz="UTC") if e.get("end") else end
        ev.append((e["kind"], e["label"], a, b))
    lims = [ax.get_xlim() for ax in axs[0]]
    for row in axs:
        for ax in row:
            for kind, _, a, b in ev:
                if kind == "blob":
                    ax.axvspan(a, b, facecolor="none", edgecolor=BLOB_COL, hatch="////", lw=0, alpha=0.35, zorder=0)
                else:
                    ax.axvspan(a, b, color=EN_COL, alpha=0.22, lw=0, zorder=0)
    for col, ax in enumerate(axs[0]):
        x0, x1 = mdates.num2date(lims[col][0]), mdates.num2date(lims[col][1])
        tr = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
        for kind, lab, a, b in ev:
            lo, hi = max(a, x0), min(b, x1)
            if lo >= hi:
                continue
            mid = lo + (hi - lo) / 2
            if kind == "blob":
                ax.text(mid, 1.075, lab, transform=tr, ha="center", va="bottom", fontsize=7, color=BLOB_COL, clip_on=False)
            else:
                ax.text(mid, 1.015, lab if col == 0 else f"El Niño {lab}", transform=tr, ha="center", va="bottom",
                        fontsize=7, color=EN_TXT, clip_on=False)
        ax.set_title(ax.get_title(loc="left"), loc="left", fontsize=9, pad=24)
    for col in (0, 1):                       # spans must not widen the time axis
        for row in axs:
            row[col].set_xlim(lims[col])


# ---------------------------------------------------------------- figures
def site_figure(s, H, G, base, year, out, events=None):
    """Per-site figure: rows T, S, O2; left full record (tide-removed daily), right last 180 days."""
    fig, axs = plt.subplots(len(VARS), 2, figsize=(11, 10.2), sharex="col", sharey="row",
                            gridspec_kw={"width_ratios": [2.4, 1], "wspace": 0.12, "hspace": 0.32})
    end = H.index.max()
    r0 = end - pd.Timedelta(days=180)
    y0 = pd.Timestamp(f"{year}-01-01", tz="UTC")
    recs = s["pre2026_records_godin"]
    for i, (v, short, lab, stat) in enumerate(VARS):
        a, b = axs[i]
        name = lab.split(" (")[0]
        gd = noon(G[v]) if v in G else pd.Series(dtype=float)
        gdp = gd.asfreq("1D") if len(gd) else gd                      # NaN on missing days -> gaps in lines
        if len(gdp):
            a.plot(gdp[gdp.index < y0].index, gdp[gdp.index < y0], color=GREY, lw=0.5)
            cur = gdp[gdp.index >= y0 - pd.Timedelta("1D")]
            a.plot(cur.index, cur, color=EV, lw=1.0)
        hh = H[v].loc[r0:].dropna()
        hhp = hh.asfreq("1h") if len(hh) else hh
        days = pd.date_range(r0.normalize(), end.normalize() + pd.Timedelta("1D"), freq="1D")
        mu, sd = baseline(base, s["code"], v, days)
        b.fill_between(days, mu - sd, mu + sd, color=BAND, lw=0); b.plot(days, mu, color=CLIM, lw=0.9)
        b.plot(hhp.index, hhp, color=GREY, lw=0.3, alpha=0.8)
        if v in G:
            g = G[v].loc[r0:]; b.plot(g.index, g, color=EV, lw=1.3)
        # records: T max, S min, O2 max and min
        rlist = [(recs[v], stat)] if v in recs else []
        if v == "O2" and "O2_min" in recs:
            rlist.append((recs["O2_min"], "min"))
        if v == "sigma0" and "sigma0_max" in recs:
            rlist.append((recs["sigma0_max"], "max"))
        this = gd[gd.index >= y0]
        parts = []
        for rec, st_ in rlist:
            for ax in (a, b):
                ax.axhline(rec["value"], color="k", ls="--", lw=0.8)
            a.text(0.01, 0.97 if st_ == "max" else 0.03, f"pre-2026 tide-removed {st_} {rec['value']:.3f} ({rec['date']})",
                   transform=a.transAxes, fontsize=7, va="top" if st_ == "max" else "bottom")
            if len(this):
                x = this.max() if st_ == "max" else this.min()
                d = this.idxmax() if st_ == "max" else this.idxmin()
                new = (x > rec["value"]) if st_ == "max" else (x < rec["value"])
                parts.append(f"{year} {st_} {x:.3f} ({d:%d %b}) vs {rec['value']:.3f}" + (" — new record" if new else ""))
        if not len(this):
            title = f"{name}: no valid data in {year}" + (" (salinity sensor flagged bad by ONC)" if v in ("salinity", "sigma0") else "")
        else:
            title = f"{name}: " + "; ".join(parts)
        if v == "sigma0" and len(this):
            title += "  (lower = lighter water)"
        a.set_title(title, loc="left", fontsize=9)
        a.set_ylabel(lab)
        vals = pd.concat([gd, hh]).dropna()
        if len(vals):
            lo, hi = np.nanpercentile(vals, [0.2, 99.8])
            for rec, _ in rlist:
                lo, hi = min(lo, rec["value"]), max(hi, rec["value"])
            pad = 0.1 * (hi - lo) if hi > lo else 0.1
            a.set_ylim(lo - pad, hi + pad)
    axs[0, 1].set_title("Last 180 days: hourly (grey), tide-removed (red)", loc="left", fontsize=9)
    axs[-1, 0].xaxis.set_major_locator(mdates.YearLocator(2)); axs[-1, 0].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    axs[-1, 1].xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 3, 5, 7, 9, 11]))
    axs[-1, 1].xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    axs[-1, 1].set_xlim(r0, end + pd.Timedelta("2D"))
    draw_events(axs, events, end)
    if events:
        fig.legend([Patch(facecolor=EN_COL, alpha=0.4), Patch(facecolor="white", edgecolor=BLOB_COL, hatch="////")],
                   ["El Niño (NOAA ONI ≥ +0.5 °C, to nearest month)", "The Blob, NE Pacific marine heatwave"],
                   loc="lower left", bbox_to_anchor=(0.06, 0.035), fontsize=7, ncol=2, frameon=False)
    fig.legend([Line2D([], [], color=GREY, lw=1), Line2D([], [], color=EV, lw=1.3), Line2D([], [], color="k", ls="--", lw=0.8),
                Patch(color=BAND), Line2D([], [], color=CLIM, lw=0.9)],
               [f"tide-removed daily, before {year} (right: hourly)", f"tide-removed, {year}", "pre-2026 record (tide-removed)",
                "day-of-year mean ± 1 sd (pre-2026)", "day-of-year mean"],
               loc="lower left", bbox_to_anchor=(0.06, 0.01), fontsize=7, ncol=5, frameon=False)
    fig.suptitle(f"ONC {s['name']} ({s['code']}, {s['depth_m']} m) — data to {end:%Y-%m-%d %H:%M} UTC", x=0.06, ha="left", fontsize=10)
    fig.text(0.06, -0.02, "Tide-removed = Godin filter (24-24-25 h), noon value each day. σθ = TEOS-10 potential density anomaly at the "
             "nominal sensor pressure.\nSalinity, σθ and oxygen are not corrected for offsets between instruments; artefact masking "
             "was applied by hand to data before 6 Nov 2026; later data are masked only where listed on the monitor page.",
             fontsize=6.5, color="0.35")
    fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig)


def map_png(cfg, st, out_dir):
    """Station map on GEBCO bathymetry; markers coloured by the latest temperature anomaly (sd)."""
    import matplotlib.patheffects as pe
    from matplotlib.colors import TwoSlopeNorm
    m = cfg["map"]; x0, y0, x1, y1 = m["bbox"]
    img = plt.imread(HERE / m["image"])
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    ax.imshow(img, extent=(x0, x1, y0, y1), origin="upper", interpolation="bilinear")
    ax.set_aspect(1 / np.cos(np.radians((y0 + y1) / 2)))
    norm = TwoSlopeNorm(vmin=-5, vcenter=0, vmax=5); cmap = plt.get_cmap("RdBu_r")
    an = st.set_index("code")["anomaly_sd"] if "anomaly_sd" in st else pd.Series(dtype=float)
    halo = [pe.withStroke(linewidth=2.5, foreground="white")]
    canyon = [s for s in cfg["sites"] if s["code"] in ("BACHY", "BACME", "BACAX")]
    for s in cfg["sites"]:
        v = an.get(s["code"], np.nan)
        col = cmap(norm(np.clip(v, -5, 5))) if pd.notna(v) else "white"
        ax.plot(s["lon"], s["lat"], "o", color=col, mec="k", ms=9, mew=0.8, zorder=5)
        if s in canyon:
            continue
        txt = f"{s['name']} {s['depth_m']} m" + (f"\n{v:+.1f} sd" if pd.notna(v) else "")
        left = s["lon"] < -126.5
        ax.annotate(txt, (s["lon"], s["lat"]), xytext=(9 if left else -9, 0), textcoords="offset points", fontsize=7.5,
                    ha="left" if left else "right", va="center", path_effects=halo)
    if canyon:
        lab = "Barkley Canyon floor\n" + "\n".join(
            f"{s['name']} {s['depth_m']} m" + (f"  {an.get(s['code'], np.nan):+.1f} sd" if pd.notna(an.get(s["code"], np.nan)) else "")
            for s in canyon)
        c0 = canyon[-1]
        ax.annotate(lab, (c0["lon"], c0["lat"]), xytext=(14, 4), textcoords="offset points", fontsize=7.5, ha="left", va="center",
                    path_effects=halo)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_xlabel("Longitude (°E)"); ax.set_ylabel("Latitude (°N)")
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("Temperature anomaly today (sd)")
    ax.text(0.01, 0.01, m["attribution"], transform=ax.transAxes, fontsize=6.5, color="0.25")
    ax.text(0.73, 0.97, "Vancouver Island", transform=ax.transAxes, fontsize=8.5, style="italic", va="top", path_effects=halo)
    ax.text(0.10, 0.12, "Pacific Ocean", transform=ax.transAxes, fontsize=8.5, style="italic", color="0.2")
    fig.savefig(out_dir / "map.png", dpi=150, bbox_inches="tight"); plt.close(fig)


def overview_pngs(cfg, base, H, G, T10, year, recent_days, run_time, out_dir):
    sites = cfg["sites"]
    y0 = pd.Timestamp(f"{year}-01-01", tz="UTC"); y1 = pd.Timestamp(f"{year}-12-31", tz="UTC")
    days = pd.date_range(y0, y1, freq="1D")
    fig, axs = plt.subplots(len(sites), 1, figsize=(7.2, 9.4), sharex=True, gridspec_kw={"hspace": 0.18})
    for ax, s in zip(axs, sites):
        c = s["code"]; mu, sd = baseline(base, c, "temperature", days)
        ax.fill_between(days, mu - sd, mu + sd, color=BAND, lw=0); ax.plot(days, mu, color=CLIM, lw=0.9)
        h = H[c]["temperature"].loc[y0:].dropna().asfreq("1h"); ax.plot(h.index, h, color=GREY, lw=0.25, alpha=0.6)
        g = G[c]["temperature"].loc[y0:]; ax.plot(g.index, g, color=EV, lw=1.2)
        ax.axhline(s["pre2026_records_godin"]["temperature"]["value"], color="k", ls="--", lw=0.7)
        ax.text(1.005, 0.5, f"{s['name']}\n{s['depth_m']} m", transform=ax.transAxes, fontsize=7, va="center")
    axs[0].set_xlim(y0, y1)
    axs[-1].xaxis.set_major_locator(mdates.MonthLocator()); axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    fig.text(0.04, 0.5, "Temperature (°C)", rotation=90, va="center", fontsize=9)
    axs[0].legend([Line2D([], [], color=EV, lw=1.2), Line2D([], [], color=GREY, lw=0.6), Line2D([], [], color=CLIM, lw=0.9),
                   Patch(color=BAND), Line2D([], [], color="k", ls="--", lw=0.7)],
                  [f"{year}, tide-removed (Godin)", "hourly", "day-of-year mean (pre-2026)", "± 1 sd", "pre-2026 tide-removed maximum"],
                  ncol=3, fontsize=6.5, loc="lower left", bbox_to_anchor=(0, 1.02), frameon=False)
    fig.suptitle(f"Temperature at ONC slope, canyon and shelf sites — updated {run_time:%Y-%m-%d %H:%M} UTC", y=0.995, fontsize=9)
    fig.savefig(out_dir / "overview_year.png", dpi=170, bbox_inches="tight"); plt.close(fig)

    end = max((T10[s["code"]].index.max() if len(T10[s["code"]]) else H[s["code"]].index.max()) for s in sites)
    r0 = end - pd.Timedelta(days=recent_days)
    hrs = pd.date_range(r0.floor("h"), end.ceil("h"), freq="1h")
    fig, axs = plt.subplots(len(sites), 1, figsize=(7.2, 9.4), sharex=True, gridspec_kw={"hspace": 0.18})
    for ax, s in zip(axs, sites):
        c = s["code"]; mu, sd = baseline(base, c, "temperature", hrs)
        ax.fill_between(hrs, mu - sd, mu + sd, color=BAND, lw=0); ax.plot(hrs, mu, color=CLIM, lw=0.9)
        x = T10[c].loc[r0:] if len(T10[c]) else H[c]["temperature"].loc[r0:].dropna()
        ax.plot(x.index, x, color=EV, lw=0.6)
        g = G[c]["temperature"].loc[r0:]; ax.plot(g.index, g, color="k", lw=1.1)
        vals = np.r_[mu - sd, mu + sd, x.dropna().to_numpy()]
        lo, hi = np.nanmin(vals), np.nanmax(vals); rec = s["pre2026_records_godin"]["temperature"]["value"]
        if rec <= hi + 0.5 * (hi - lo):
            ax.axhline(rec, color="k", ls="--", lw=0.7); hi = max(hi, rec)
        else:
            ax.text(0.01, 0.95, f"pre-2026 tide-removed max {rec:.2f} °C (off scale)", transform=ax.transAxes, fontsize=6.5, va="top")
        pad = 0.06 * (hi - lo); ax.set_ylim(lo - pad, hi + pad)
        ax.text(1.005, 0.5, f"{s['name']}\n{s['depth_m']} m", transform=ax.transAxes, fontsize=7, va="center")
    axs[0].set_xlim(r0, end)
    axs[-1].xaxis.set_major_locator(mdates.DayLocator()); axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
    plt.setp(axs[-1].get_xticklabels(), rotation=30, ha="right")
    fig.text(0.04, 0.5, "Temperature (°C)", rotation=90, va="center", fontsize=9)
    axs[0].legend([Line2D([], [], color=EV, lw=0.8), Line2D([], [], color="k", lw=1.1), Line2D([], [], color=CLIM, lw=0.9),
                   Patch(color=BAND), Line2D([], [], color="k", ls="--", lw=0.7)],
                  ["10-min mean", "tide-removed (Godin; ends ~35 h before latest data)", "day-of-year mean (pre-2026)", "± 1 sd",
                   "pre-2026 tide-removed maximum"], ncol=2, fontsize=6.5, loc="lower left", bbox_to_anchor=(0, 1.02), frameon=False)
    fig.suptitle(f"Last {recent_days} days at 10-minute resolution — updated {run_time:%Y-%m-%d %H:%M} UTC", y=0.995, fontsize=9)
    fig.savefig(out_dir / "overview_recent.png", dpi=170, bbox_inches="tight"); plt.close(fig)


def mask_note(cfg):
    items = []
    for s in cfg["sites"]:
        for m in s.get("masks", []):
            t0 = pd.Timestamp(m["start"]).strftime("%d %b %Y %H:%M UTC")
            t1 = pd.Timestamp(m["end"]).strftime("%d %b %Y %H:%M UTC") if m.get("end") else "present"
            items.append(f"<li>{s['name']}: {', '.join(m['vars'])} (and σθ) removed {t0} – {t1}. {m.get('reason', '')}</li>")
    return f'<p class="note"><b>Data removed by hand:</b></p><ul class="note">{"".join(items)}</ul>' if items else ""


def write_html(cfg, base, H, G, T10, st, year, hires_days, run_time, out_dir, problems):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    sites = cfg["sites"]; n = len(sites)
    y0 = pd.Timestamp(f"{year}-01-01", tz="UTC")
    days = pd.date_range(y0, pd.Timestamp(f"{year}-12-31", tz="UTC"), freq="1D")
    fig = make_subplots(rows=n, cols=1, shared_xaxes=True, vertical_spacing=0.012,
                        subplot_titles=[f"{s['name']} ({s['depth_m']} m)" for s in sites])
    for i, s in enumerate(sites, start=1):
        c = s["code"]; show = i == 1
        mu, sd = baseline(base, c, "temperature", days)
        fig.add_trace(go.Scatter(x=days, y=mu + sd, line=dict(width=0), hoverinfo="skip", showlegend=False), i, 1)
        fig.add_trace(go.Scatter(x=days, y=mu - sd, fill="tonexty", fillcolor=BAND, line=dict(width=0), hoverinfo="skip",
                                 name="± 1 sd (pre-2026)", legendgroup="band", showlegend=show), i, 1)
        fig.add_trace(go.Scatter(x=days, y=mu, line=dict(color=CLIM, width=1), name="day-of-year mean (pre-2026)",
                                 legendgroup="clim", showlegend=show, hovertemplate="normal %{y:.3f} °C<extra></extra>"), i, 1)
        h = H[c]
        end = T10[c].index.max() if len(T10[c]) else h.index.max(); cut = end - pd.Timedelta(days=hires_days)
        hi_res = T10[c].loc[cut:] if len(T10[c]) else h["temperature"].loc[cut:]
        x = pd.concat([h["temperature"].loc[y0:cut], hi_res]).dropna()
        fig.add_trace(go.Scattergl(x=x.index, y=x.round(4), line=dict(color="rgba(120,120,120,0.6)", width=0.8),
                                   name=f"10-min (last {hires_days} d), hourly before", legendgroup="hi", showlegend=show,
                                   hovertemplate="%{x|%d %b %H:%M} UTC  %{y:.3f} °C<extra></extra>"), i, 1)
        g = G[c]["temperature"].loc[y0:].dropna()
        g = g[g.index.hour % 3 == 0]
        fig.add_trace(go.Scatter(x=g.index, y=g.round(4), line=dict(color=EV, width=2), name="tide-removed (Godin)",
                                 legendgroup="g", showlegend=show, hovertemplate="%{x|%d %b %H:%M}  filtered %{y:.3f} °C<extra></extra>"), i, 1)
        fig.add_hline(y=s["pre2026_records_godin"]["temperature"]["value"], line=dict(color="black", dash="dash", width=1), row=i, col=1)
        fig.update_yaxes(title_text="°C", row=i, col=1)
    last_end = max((T10[s["code"]].index.max() if len(T10[s["code"]]) else H[s["code"]].index.max()) for s in sites)
    fig.update_xaxes(range=[last_end - pd.Timedelta(days=30), last_end + pd.Timedelta(hours=6)])
    fig.update_xaxes(rangeselector=dict(buttons=[dict(count=3, label="3 d", step="day", stepmode="backward"),
                                                 dict(count=10, label="10 d", step="day", stepmode="backward"),
                                                 dict(count=1, label="1 mo", step="month", stepmode="backward"),
                                                 dict(count=3, label="3 mo", step="month", stepmode="backward"),
                                                 dict(step="all", label="year")]), row=1, col=1)
    fig.update_layout(height=230 * n, template="simple_white", hovermode="x unified",
                      legend=dict(orientation="h", y=1.03, x=0, font=dict(size=11)), margin=dict(l=60, r=20, t=110, b=40),
                      title=dict(text=f"Temperature, {year} — dashed line: pre-2026 tide-removed maximum", x=0, font=dict(size=14)))
    # ---- compact summary table (the one shown in the overview)
    def cell_anom(x, xs, unit, fmt):
        return "—" if pd.isna(x) else f"{x:{fmt}} {unit} ({xs:+.1f} sd)"
    srows = []
    for _, r in st.iterrows():
        rec = (f"yes, {r['max_filtered_this_year_C']:.2f} °C on {r['max_filtered_date']}"
               if r.get("new_record_this_year") is True else "no")
        light = "yes" if r.get("lightest_on_record_this_year") is True else ("—" if pd.isna(r.get("sigma0_anomaly")) else "no")
        srows.append(f"<tr><td class='l'>{r['site']}</td><td>{r['depth_m']} m</td>"
                     f"<td>{cell_anom(r.get('anomaly_C'), r.get('anomaly_sd'), '°C', '+.2f')}</td>"
                     f"<td>{cell_anom(r.get('sigma0_anomaly'), r.get('sigma0_anomaly_sd'), 'kg m⁻³', '+.3f')}</td>"
                     f"<td class='l'>{rec}</td><td class='l'>{light}</td></tr>")
    vt = st["filtered_valid_to_utc"].dropna()
    valid = pd.to_datetime(vt).min().strftime("%d %b %Y") if len(vt) else "—"
    summary = ("<table class='status'><tr><th class='l'>Site</th><th>Depth</th><th>Temperature anomaly</th><th>σθ anomaly</th>"
               f"<th class='l'>New {year} temperature record?</th><th class='l'>Lightest σθ on record in {year}?</th></tr>"
               + "".join(srows) + "</table>")
    cols = ["site", "depth_m", "last_data_utc", "last_25h_mean_C", "filtered_valid_to_utc", "filtered_C", "normal_C",
            "anomaly_C", "anomaly_sd", "max_filtered_this_year_C", "pre2026_record_C", "new_record_this_year",
            "sigma0_filtered", "sigma0_normal", "sigma0_anomaly", "sigma0_anomaly_sd", "sigma0_min_this_year",
            "pre2026_sigma0_min", "lightest_on_record_this_year"]
    detail = st.reindex(columns=cols).to_html(index=False, border=0, classes="status", na_rep="")
    warn = "".join(f"<p class='warn'>Data problem: {p}</p>" for p in problems)
    imgs = "".join(f"<h3 id='{s['code']}'>{s['name']} ({s['depth_m']} m)</h3><img src='site_{s['code']}.png' alt='{s['name']}'>" for s in sites)
    links = " · ".join(f"<a href='#{s['code']}'>{s['name']}</a>" for s in sites)
    page = f"""<!doctype html><html><head><meta charset="utf-8"><title>ONC slope and canyon monitor</title>
<style>body{{font-family:sans-serif;margin:16px;max-width:1500px}} table.status{{border-collapse:collapse;font-size:13px}}
table.status td,table.status th{{padding:3px 10px;border-bottom:1px solid #ddd;text-align:right}} .l{{text-align:left!important}}
img{{max-width:100%}} img.map{{max-width:900px}} p.warn{{color:#b00}} p.note{{color:#444;font-size:13px;max-width:1100px}}
details{{margin:8px 0 16px 0}} summary{{cursor:pointer;color:#333}}</style></head><body>
<h2>ONC slope, canyon and shelf moorings off southwest Vancouver Island</h2>
<p>Updated {run_time:%Y-%m-%d %H:%M} UTC. Site figures: {links}</p>
<img class="map" src="map.png" alt="Station map">
<h3>Today against the pre-2026 normal for the date (tide-removed values to {valid})</h3>
{summary}{warn}
<p class="note">Anomaly = latest tide-removed daily value minus the pre-2026 normal for that calendar date; "sd" is the spread of
pre-2026 daily values within ±7 days of the date. This is a day-against-day comparison and is not the window-mean sd used in the
report. σθ = TEOS-10 potential density anomaly; a negative σθ anomaly means lighter water than normal at that depth (density
surfaces pushed down). σθ needs valid salinity, so it is unavailable at Clayoquot Slope, whose salinity sensor has been flagged bad
by ONC since Dec 2024.</p>{mask_note(cfg)}
<details><summary>Full status table (also in <a href="status.csv">status.csv</a>)</summary>{detail}</details>
<p class="note">Tide-removed values use a Godin filter (24-, 24- and 25-h running means), so they lag the latest data by ~35 h;
'last 25 h mean' is the most recent near-tide-free value. Baseline: {cfg['baseline_note']} Data: Ocean Networks Canada Oceans 3.0
(10-minute averages, ONC good-quality flags only, isolated salinity/oxygen spikes removed). Contains data from Ocean Networks Canada
(CC-BY 4.0). {cfg['map']['attribution']}.</p>
{fig.to_html(include_plotlyjs="cdn", full_html=False)}
<h2>Site records: temperature, salinity, oxygen and σθ</h2>
<p class="note">{cfg.get("events_note", "")}</p>{imgs}</body></html>"""
    (out_dir / "index.html").write_text(page, encoding="utf-8")


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recent-days", type=int, default=10)
    ap.add_argument("--hires-days", type=int, default=60)
    ap.add_argument("--out", default=str(HERE / "output"))
    args = ap.parse_args()

    cfg = json.loads((HERE / "sites.json").read_text())
    base = pd.read_csv(HERE / "baseline_doy.csv")
    cache_dir, out_dir = HERE / "cache", Path(args.out)
    cache_dir.mkdir(exist_ok=True); out_dir.mkdir(parents=True, exist_ok=True)
    tok = token()
    run_time = datetime.now(timezone.utc)
    year = run_time.year

    H, G, T10, rows, problems = {}, {}, {}, [], []
    for s in cfg["sites"]:
        c = s["code"]
        hist = pd.read_parquet(HERE / "history" / f"{c}_hourly_qc.parquet").astype(float)
        hist_end = hist.index.max()
        try:
            live = update_cache(s, tok, hist_end, args.hires_days, cache_dir)
        except Exception as e:
            problems.append(f"{c}: {e}")
            f = cache_dir / f"{c}_10min.parquet"
            live = pd.read_parquet(f) if f.exists() else pd.DataFrame(columns=["T", "S", "O2"])
        lh = live.rename(columns={"T": "temperature", "S": "salinity"}).resample("1h").mean() if len(live) else live
        if len(lh):
            lh.index = lh.index + pd.Timedelta("30min")          # match history's hourly-centre stamps
            lh = lh[lh.index > hist_end]
        h = pd.concat([hist, lh[["temperature", "salinity", "O2"]]]) if len(lh) else hist
        h = h[~h.index.duplicated()].sort_index()
        h = apply_masks(h, s)
        h = add_sigma0(h, s)
        H[c] = h
        T10[c] = live["T"].dropna() if len(live) else pd.Series(dtype=float)
        G[c] = {v: godin(h[v].dropna()) for v, *_ in VARS if h[v].notna().sum() > 100}
        G[c] = pd.DataFrame(G[c])

        # status row
        t = h["temperature"].dropna(); g = G[c]["temperature"].dropna()
        rec = s["pre2026_records_godin"]["temperature"]
        row = dict(site=s["name"], code=c, depth_m=s["depth_m"])
        if len(t):
            src = T10[c] if len(T10[c]) else t
            last = src.index.max()
            row.update(last_data_utc=last.strftime("%Y-%m-%d %H:%M"),
                       last_25h_mean_C=round(src.loc[last - pd.Timedelta("25h"):].mean(), 3))
        if len(g):
            gt = g.index.max(); mu, sd = baseline(base, c, "temperature", pd.DatetimeIndex([gt]))
            gy = noon(G[c]["temperature"]); gy = gy[gy.index >= pd.Timestamp(f"{year}-01-01", tz="UTC")]
            row.update(filtered_valid_to_utc=gt.strftime("%Y-%m-%d %H:%M"), filtered_C=round(g.iloc[-1], 3),
                       normal_C=round(mu[0], 3), anomaly_C=round(g.iloc[-1] - mu[0], 3), anomaly_sd=round((g.iloc[-1] - mu[0]) / sd[0], 1),
                       max_filtered_this_year_C=round(gy.max(), 3) if len(gy) else np.nan,
                       pre2026_record_C=rec["value"], new_record_this_year=bool(len(gy) and gy.max() > rec["value"]))
        gs = G[c]["sigma0"].dropna() if "sigma0" in G[c] else pd.Series(dtype=float)
        if len(gs):
            gsy = noon(G[c]["sigma0"]).dropna(); gsy = gsy[gsy.index >= pd.Timestamp(f"{year}-01-01", tz="UTC")]
            srec = s["pre2026_records_godin"]["sigma0"]
            if len(gsy):                 # records use all valid σθ this year, even if salinity is not current
                row.update(sigma0_min_this_year=round(gsy.min(), 3), pre2026_sigma0_min=srec["value"],
                           lightest_on_record_this_year=bool(gsy.min() < srec["value"]))
        if len(gs) and len(g) and gs.index.max() < g.index.max() - pd.Timedelta("3D"):
            gs = gs.iloc[0:0]            # no recent valid salinity -> no current σθ anomaly
        if len(gs):
            mu, sd = baseline(base, c, "sigma0", pd.DatetimeIndex([gs.index.max()]))
            row.update(sigma0_filtered=round(gs.iloc[-1], 3), sigma0_normal=round(mu[0], 3),
                       sigma0_anomaly=round(gs.iloc[-1] - mu[0], 3), sigma0_anomaly_sd=round((gs.iloc[-1] - mu[0]) / sd[0], 1))
        if len(g):
            gy = noon(G[c]["temperature"]); gy = gy[gy.index >= pd.Timestamp(f"{year}-01-01", tz="UTC")]
            if len(gy):
                row.update(max_filtered_date=gy.idxmax().strftime("%d %b"))
        rows.append(row)
        site_figure(s, h, G[c], base, year, out_dir / f"site_{c}.png", events=cfg.get("events"))

    st = pd.DataFrame(rows); st.to_csv(out_dir / "status.csv", index=False)
    map_png(cfg, st, out_dir)
    overview_pngs(cfg, base, H, G, T10, year, args.recent_days, run_time, out_dir)
    write_html(cfg, base, H, G, T10, st, year, args.hires_days, run_time, out_dir, problems)
    print(f"Updated {run_time:%Y-%m-%d %H:%M} UTC")
    print(st.reindex(columns=["site", "last_data_utc", "filtered_valid_to_utc", "anomaly_C", "anomaly_sd",
                              "sigma0_anomaly", "sigma0_anomaly_sd", "new_record_this_year", "lightest_on_record_this_year"]).to_string(index=False))
    for p in problems:
        print("WARNING", p)
    return 1 if len(problems) == len(cfg["sites"]) else 0


if __name__ == "__main__":
    sys.exit(main())
