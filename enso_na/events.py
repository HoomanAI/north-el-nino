"""
El Nino event identification, phase definition, life-cycle composites and the
point-symmetric inversion used to compare development against decay.

Conventions follow Song et al. (2026, npj Clim Atmos Sci):
  * month index 0  = January of the development year   ("Jan0")
  * month index 11 = December of the development year  ("Dec0", mature phase)
  * month index 23 = December of the decay year        ("Dec1")
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from .config import ENSO_THRESHOLD

LIFECYCLE_LEN = 24  # Jan0 .. Dec1
MONTH_LABELS = [
    "Jan0", "Feb0", "Mar0", "Apr0", "May0", "Jun0", "Jul0", "Aug0",
    "Sep0", "Oct0", "Nov0", "Dec0", "Jan1", "Feb1", "Mar1", "Apr1",
    "May1", "Jun1", "Jul1", "Aug1", "Sep1", "Oct1", "Nov1", "Dec1",
]


def identify_el_nino(nino: pd.Series, threshold: float = ENSO_THRESHOLD,
                     exclude: tuple = ()) -> list[int]:
    """Return development years of El Nino events.

    An event is declared when the DJF mean of the 3-month running-mean Nino3.4
    anomaly (Dec of year0 to Feb of year0+1) exceeds ``threshold``.  Multi-year
    events are attributed to their first qualifying year; a year whose
    predecessor already qualified is skipped so that each event contributes one
    development/decay pair.
    """
    r3 = nino.rolling(3, center=True, min_periods=3).mean()
    years = sorted({t.year for t in r3.index})
    djf = {}
    for y in years:
        sel = [pd.Timestamp(f"{y}-12-01"), pd.Timestamp(f"{y + 1}-01-01"),
               pd.Timestamp(f"{y + 1}-02-01")]
        vals = [r3.get(t, np.nan) for t in sel]
        if np.all(np.isfinite(vals)):
            djf[y] = float(np.mean(vals))
    ev, prev = [], -99
    for y in sorted(djf):
        if djf[y] >= threshold and y not in exclude:
            if y - 1 == prev:      # continuation of a multi-year warm state
                prev = y
                continue
            ev.append(y)
            prev = y
    return ev


def lifecycle_matrix(series: pd.Series, years: list[int]) -> np.ndarray:
    """(n_event, 24) matrix of a monthly series over the Jan0-Dec1 life cycle."""
    out = np.full((len(years), LIFECYCLE_LEN), np.nan)
    for i, y in enumerate(years):
        idx = pd.date_range(f"{y}-01-01", periods=LIFECYCLE_LEN, freq="MS")
        out[i] = series.reindex(idx).values
    return out


def composite(series: pd.Series, years: list[int]):
    """Ensemble mean and 95% confidence half-width of the life-cycle composite."""
    m = lifecycle_matrix(series, years)
    mean = np.nanmean(m, axis=0)
    n = np.sum(np.isfinite(m), axis=0)
    sd = np.nanstd(m, axis=0, ddof=1)
    ci = stats.t.ppf(0.975, np.maximum(n - 1, 1)) * sd / np.sqrt(np.maximum(n, 1))
    return mean, ci, m


def phase_speeds(series: pd.Series, years: list[int]):
    """Development and decay speed of an index for each event.

    Development speed  = mean monthly tendency from the minimum of the
    development year (or zero crossing) to the peak.
    Decay speed        = |mean monthly tendency| from the peak to the minimum of
    the decay year (or zero crossing).
    Returns a DataFrame indexed by development year.
    """
    m = lifecycle_matrix(series, years)
    rec = []
    for i, y in enumerate(years):
        v = m[i]
        if np.all(np.isnan(v)):
            continue
        # peak searched in Jul0-Mar1 (months 6-14)
        win = np.arange(6, 15)
        pk = win[np.nanargmax(v[win])]
        pre = v[: pk + 1]
        post = v[pk:]
        i0 = int(np.nanargmin(pre))
        i1 = int(np.nanargmin(post)) + pk
        dev = (v[pk] - v[i0]) / max(pk - i0, 1)
        dec = (v[pk] - v[i1]) / max(i1 - pk, 1)
        rec.append(dict(year=y, peak_month=pk, peak=v[pk],
                        dev_speed=dev, dec_speed=dec,
                        dev_len=pk - i0, dec_len=i1 - pk,
                        ratio=dec / dev if dev != 0 else np.nan))
    return pd.DataFrame(rec).set_index("year")


def point_symmetric(tend_composite: np.ndarray, pivot: int = 11) -> np.ndarray:
    """Point-symmetric inversion of a composite tendency about ``pivot``.

    The decay-phase tendency is rotated 180 degrees about the mature phase so
    that it can be overlaid on the development-phase tendency; where the two
    curves separate, development and decay proceed at different speeds.
    """
    n = len(tend_composite)
    out = np.full(n, np.nan)
    for k in range(n):
        j = 2 * pivot - k
        if 0 <= j < n:
            out[k] = -tend_composite[j]
    return out


def significant_months(tend_matrix: np.ndarray, pivot: int = 11, alpha: float = 0.05):
    """Months where development and (inverted) decay tendencies differ.

    A two-tailed paired Student's t-test is applied across events between the
    tendency at month k and minus the tendency at its point-symmetric partner.
    """
    n = tend_matrix.shape[1]
    p = np.full(n, np.nan)
    for k in range(n):
        j = 2 * pivot - k
        if not (0 <= j < n):
            continue
        a, b = tend_matrix[:, k], -tend_matrix[:, j]
        ok = np.isfinite(a) & np.isfinite(b)
        if ok.sum() > 3:
            p[k] = stats.ttest_rel(a[ok], b[ok]).pvalue
    return p, np.where(p < alpha)[0]


def ci95(sample: np.ndarray, axis: int = 0):
    """95% confidence half-width of a sample mean."""
    sample = np.asarray(sample, dtype=float)
    n = np.sum(np.isfinite(sample), axis=axis)
    sd = np.nanstd(sample, axis=axis, ddof=1)
    return stats.t.ppf(0.975, np.maximum(n - 1, 1)) * sd / np.sqrt(np.maximum(n, 1))


#: Events excluded because positive SSTAs persisted into the following event,
#: which makes their decay phase ill-defined (same criterion as Song et al. 2026,
#: who excluded 1986/1987 and 2014/2015).
EXCLUDED = (1986, 2014)


def default_events(df, nino_col: str = "NINO34") -> list[int]:
    """The El Nino event set used throughout the paper."""
    return identify_el_nino(df[nino_col], exclude=EXCLUDED)
