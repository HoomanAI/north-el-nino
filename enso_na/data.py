"""
Loading and homogenisation of the climate index files.

The 13 NetCDF files supplied with this study carry three different native
sampling rates (monthly, weekly, daily) and two different time-axis names
(``TSAXIS`` for the monthly CPC/PSL teleconnection indices, ``TIME`` for the
weekly OISST Nino indices and the daily OLR index).  Everything is converted
here to a single monthly-mean table on a month-start (``MS``) index.

Original public sources of the supplied files are documented in
``SOURCES`` below and reproduced in the README.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr

from .config import ANALYSIS_END, ANALYSIS_START, DATA_DIR

# ---------------------------------------------------------------------------
# File -> (variable name, native sampling) map
# ---------------------------------------------------------------------------
FILES = {
    "NINO12": ("nino12.nc", "NINO12", "weekly"),
    "NINO3": ("nino3.nc", "NINO3", "weekly"),
    "NINO34": ("nino34.nc", "NINO34", "weekly"),
    "NINO4": ("nino4.nc", "NINO4", "weekly"),
    "TNA": ("tna.nc", "TNA", "weekly"),
    "OLR": ("olren.nc", "OLREN", "daily"),
    "PNA": ("pna.nc", "PNA", "monthly"),
    "NAO": ("nao.nc", "NAO", "monthly"),
    "AO": ("ao.nc", "AO", "monthly"),
    "PDO": ("pdo.nc", "PDO", "monthly"),
    "SOI": ("soi.nc", "SOI", "monthly"),
    "AMO": ("amo.nc", "AMO", "monthly"),
}

#: Annual La Nina-like OLR day-count index (used only in the diversity section).
ANNUAL_FILES = {"OLRLN": ("olrln.nc", "OLRLN")}

SOURCES = {
    "NINO12/3/3.4/4": "NOAA CPC weekly ENSO indices (OISSTv2): "
    "https://www.cpc.ncep.noaa.gov/data/indices/wksst8110.for",
    "TNA": "NOAA PSL Tropical North Atlantic index: "
    "https://psl.noaa.gov/gcos_wgsp/Timeseries/TNA/",
    "OLR (OLREN)": "NOAA interpolated OLR anomaly averaged over the equatorial "
    "central Pacific: https://psl.noaa.gov/data/gridded/data.olrcdr.interp.html",
    "OLRLN": "Annual count of suppressed-convection (La Nina-like) OLR days, "
    "derived from the same daily OLR record.",
    "PNA/NAO/AO": "NOAA CPC monthly teleconnection indices: "
    "https://www.cpc.ncep.noaa.gov/data/teledoc/telecontents.shtml",
    "PDO": "NOAA NCEI/PSL Pacific Decadal Oscillation index: "
    "https://www.ncei.noaa.gov/access/monitoring/pdo/",
    "SOI": "NOAA CPC Southern Oscillation Index: "
    "https://www.cpc.ncep.noaa.gov/data/indices/soi",
    "AMO": "NOAA PSL Atlantic Multidecadal Oscillation index (unsmoothed): "
    "https://psl.noaa.gov/data/timeseries/AMO/",
}


def _read_series(fname: str, varname: str) -> pd.Series:
    """Read one index file and return a pandas Series on its native time axis."""
    ds = xr.open_dataset(DATA_DIR / fname, decode_times=True)
    tname = "TSAXIS" if "TSAXIS" in ds.dims else ("TIME" if "TIME" in ds.dims else "ANNUAL")
    s = pd.Series(
        np.asarray(ds[varname].values, dtype=float),
        index=pd.to_datetime(ds[tname].values),
        name=varname,
    )
    ds.close()
    return s.replace([-999.0, -99.9, -9.99e8], np.nan)


def _to_monthly(s: pd.Series, sampling: str) -> pd.Series:
    """Collapse a native-resolution series onto a month-start monthly mean."""
    if sampling == "monthly":
        out = s.copy()
        out.index = out.index.to_period("M").to_timestamp()
        return out.groupby(level=0).mean()
    # weekly / daily: require at least half the expected samples in a month
    grp = s.resample("MS")
    counts = grp.count()
    means = grp.mean()
    min_n = {"weekly": 2, "daily": 15}[sampling]
    return means.where(counts >= min_n)


def load_monthly() -> pd.DataFrame:
    """Return the homogenised monthly index table over the analysis period."""
    cols = {}
    for key, (fname, var, samp) in FILES.items():
        cols[key] = _to_monthly(_read_series(fname, var), samp)
    df = pd.DataFrame(cols).sort_index()
    df = df.loc[ANALYSIS_START:ANALYSIS_END]
    # linear interpolation of isolated gaps (<=2 months) only
    df = df.interpolate(limit=2, limit_area="inside")
    return df


def load_monthly_extended(start: str = "1950-01", end: str = "2026-06") -> pd.DataFrame:
    """Monthly table for the extended record.

    Only the indices that reach back to the pre-satellite era are available:
    PNA, NAO, AO, PDO, SOI and AMO.  The tropical convective (OLR) and tropical
    North Atlantic (TNA) indices begin in 1982 and 1981 respectively and are
    therefore absent.  SOI is inverted and standardised to serve as the ENSO
    forcing index.
    """
    keys = ["PNA", "NAO", "AO", "PDO", "SOI", "AMO"]
    cols = {}
    for key in keys:
        fname, var, samp = FILES[key]
        cols[key] = _to_monthly(_read_series(fname, var), samp)
    df = pd.DataFrame(cols).sort_index().loc[start:end]
    return df.interpolate(limit=2, limit_area="inside")


def load_annual() -> pd.DataFrame:
    cols = {}
    for key, (fname, var) in ANNUAL_FILES.items():
        s = _read_series(fname, var)
        s.index = s.index.year
        cols[key] = s
    return pd.DataFrame(cols)


def standardise(df: pd.DataFrame, cols=None, period=None) -> pd.DataFrame:
    """Remove the mean and divide by the standard deviation of each column.

    Indices supplied here are already anomalies, but they carry different units
    (degC, W m-2, non-dimensional).  Standardising over a common reference
    period puts every budget term in units of PNA standard deviations per month
    and makes the regression coefficients directly comparable.
    """
    out = df.copy()
    cols = cols if cols is not None else df.columns
    ref = out.loc[period[0]:period[1]] if period else out
    for c in cols:
        mu, sd = ref[c].mean(), ref[c].std(ddof=1)
        out[c] = (out[c] - mu) / sd
    return out


def running_mean(obj, window: int):
    """Centred running mean of ``window`` months.

    Odd windows use the plain centred mean.  Even windows use the standard
    centred even filter with half weights on the two end points
    (1/2, 1, ..., 1, 1/2)/window, so that the filtered value remains centred on
    a calendar month rather than falling between two months.
    """
    if window <= 1:
        return obj.copy()
    if window % 2 == 1:
        return obj.rolling(window=window, center=True, min_periods=window).mean()
    w = np.ones(window + 1)
    w[0] = w[-1] = 0.5
    w = w / w.sum()

    def _filt(a):
        out = np.full(a.shape, np.nan)
        h = window // 2
        for i in range(h, len(a) - h):
            seg = a[i - h:i + h + 1]
            if np.all(np.isfinite(seg)):
                out[i] = np.dot(seg, w)
        return out

    if isinstance(obj, pd.Series):
        return pd.Series(_filt(obj.values.astype(float)), index=obj.index,
                         name=obj.name)
    out = obj.copy()
    for c in out.columns:
        out[c] = _filt(out[c].values.astype(float))
    return out


def tendency(s: pd.Series) -> pd.Series:
    """Centred first difference (units: input units per month)."""
    return (s.shift(-1) - s.shift(1)) / 2.0
