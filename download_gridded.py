#!/usr/bin/env python3
"""
OPTIONAL extension: North American composite maps from gridded reanalysis.

The analysis in the paper is deliberately self-contained and uses only the index
files supplied with the study.  This script adds the spatial picture over North
America -- 500 hPa height, surface air temperature and precipitation composites
for the ASO0 and FMA1 key periods -- by downloading four public gridded datasets.

It is NOT required to reproduce any figure in the paper, and it is kept separate
because it needs network access and several hundred megabytes of downloads.

Run it as:

    python download_gridded.py --download       # fetch the files (~700 MB)
    python download_gridded.py                  # make the maps from local files

Data sources (verify the current path at https://psl.noaa.gov/data/gridded/ if a
download fails; NOAA PSL occasionally reorganises its directories):

  500 hPa geopotential height   NCEP/NCAR Reanalysis 1, monthly
      https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/Monthlies/pressure/hgt.mon.mean.nc
  Surface air temperature       GHCN-CAMS land 0.5 deg, monthly
      https://downloads.psl.noaa.gov/Datasets/ghcncams/air.mon.mean.nc
  Precipitation                 GPCC / CPC global monthly precipitation
      https://downloads.psl.noaa.gov/Datasets/cpc_global_precip/precip.mon.mean.nc
  Sea surface temperature       NOAA ERSST v5, monthly
      https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc

Optional plotting dependency: cartopy (``pip install cartopy``).  Without it the
maps are drawn on a plain lon/lat grid, which is adequate for inspection.
"""
from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))

from enso_na.config import FIG_DIR, KEY_DEC_MONTHS, KEY_DEV_MONTHS
from enso_na.data import load_monthly
from enso_na.events import default_events

GRID_DIR = Path(__file__).resolve().parent / "data_gridded"

SOURCES = {
    "hgt": ("https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/Monthlies/"
            "pressure/hgt.mon.mean.nc", "hgt", "500 hPa height (m)"),
    "air": ("https://downloads.psl.noaa.gov/Datasets/ghcncams/air.mon.mean.nc",
            "air", "surface air temperature (K)"),
    "precip": ("https://downloads.psl.noaa.gov/Datasets/cpc_global_precip/"
               "precip.mon.mean.nc", "precip", "precipitation (mm day-1)"),
    "sst": ("https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc",
            "sst", "sea surface temperature (degC)"),
}

#: North American map domain
DOMAIN = dict(lon=(190, 310), lat=(15, 75))


def download():
    GRID_DIR.mkdir(exist_ok=True)
    for key, (url, _, _) in SOURCES.items():
        dest = GRID_DIR / f"{key}.nc"
        if dest.exists():
            print(f"{dest.name}: already present")
            continue
        print(f"downloading {url}")
        urllib.request.urlretrieve(url, dest)
        print(f"  -> {dest} ({dest.stat().st_size / 1e6:.0f} MB)")


def _anomaly(da, clim_period=("1991-01-01", "2020-12-01")):
    clim = da.sel(time=slice(*clim_period)).groupby("time.month").mean("time")
    return da.groupby("time.month") - clim


def composite_map(da, years, months):
    """Composite of a monthly DataArray over given life-cycle month indices."""
    import pandas as pd
    stack = []
    for y in years:
        for m in months:
            t = pd.Timestamp(f"{y}-01-01") + pd.DateOffset(months=int(m))
            try:
                stack.append(da.sel(time=t, method="nearest"))
            except Exception:
                pass
    return sum(stack) / len(stack)


def main():
    try:
        import xarray as xr
    except ImportError:
        raise SystemExit("xarray is required: pip install xarray netCDF4")

    missing = [k for k in SOURCES if not (GRID_DIR / f"{k}.nc").exists()]
    if missing:
        raise SystemExit(
            f"missing gridded files: {missing}\n"
            f"run `python {Path(__file__).name} --download` first, or place the "
            f"files listed in the module docstring in {GRID_DIR}")

    ev = default_events(load_monthly())
    fig, axes = plt.subplots(3, 2, figsize=(8.0, 7.5))

    for row, key in enumerate(["hgt", "air", "precip"]):
        url, var, label = SOURCES[key]
        ds = xr.open_dataset(GRID_DIR / f"{key}.nc")
        da = ds[var]
        if "level" in da.dims:
            da = da.sel(level=500)
        an = _anomaly(da)
        for col, (months, name) in enumerate(
                [(KEY_DEV_MONTHS, "ASO0"), (KEY_DEC_MONTHS, "FMA1")]):
            c = composite_map(an, ev, months)
            c = c.sortby("lat")
            sub = c.sel(lat=slice(*DOMAIN["lat"]))
            lons = sub.lon.values % 360
            order = np.argsort(lons)
            sub = sub.isel(lon=order)
            sub = sub.sel(lon=slice(*DOMAIN["lon"]))
            ax = axes[row, col]
            lim = float(np.nanpercentile(np.abs(sub.values), 98))
            m = ax.pcolormesh(sub.lon, sub.lat, sub.values, cmap="RdBu_r",
                              vmin=-lim, vmax=lim, shading="auto")
            ax.set_title(f"{label} -- {name}", fontsize=8)
            ax.tick_params(labelsize=6)
            fig.colorbar(m, ax=ax, fraction=0.035, pad=0.02).ax.tick_params(
                labelsize=6)
        ds.close()

    fig.suptitle("North American composite anomalies during El Nino "
                 "development and decay", fontsize=10)
    fig.tight_layout()
    out = FIG_DIR / "figS06_north_america_maps.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--download", action="store_true",
                    help="fetch the gridded files before plotting")
    args = ap.parse_args()
    if args.download:
        download()
    main()
