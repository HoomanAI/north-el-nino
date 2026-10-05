"""
Window-resolved tendency budget for the Pacific-North American (PNA)
teleconnection, and the three asymmetry diagnostics (AsyI, AAE, WTI).

This is the statistical analogue of the mixed-layer heat budget of Song et al.
(2026).  The tendency of the North American teleconnection index P is written
as the sum of an internal damping term, a set of remote forcing terms, and a
residual:

    dP/dt = -lambda(t) P                 internal damping
            + a_N(t)  N                  tropical Pacific SST forcing
            + a_C(t)  C                  tropical convective (OLR) forcing
            + a_S(t)  dN/dt              ENSO tendency forcing
            + a_D(t)  D                  North Pacific background (PDO)
            + a_T(t)  T                  tropical North Atlantic forcing
            + a_A(t)  A                  Arctic Oscillation
            + a_M(t)  M                  North Atlantic Oscillation
            + a_L(t)  NL                 nonlinear ENSO forcing
            + R                          residual

Two design choices matter and are documented in the Methods of the paper:

1. ``dP/dt`` is a **forward** finite difference.  For a first-order
   autoregressive process the centred difference is exactly orthogonal to the
   instantaneous state, so a centred tendency cannot identify the damping
   coefficient.

2. Every coefficient is allowed an **annual harmonic**,
   a(t) = a0 + a1 cos(2 pi m / 12) + a2 sin(2 pi m / 12), fitted by ordinary
   least squares over all months of the record.  Setting a1 = a2 = 0 gives the
   season-invariant *null model* used for the amplitude/seasonal decomposition.

All predictors and P are standardised, so every term carries units of PNA
standard deviations per month (sigma month^-1) and the coefficients are
directly comparable.

The free methodological parameter is the temporal integration window k
(1-7 months) over which every field is smoothed before the budget is formed.
It plays the role of the mixed-layer integration depth: it is arbitrary, it
must be chosen, and the attribution depends on it.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as _st

from .config import KEY_DEC_MONTHS, KEY_DEV_MONTHS, TERM_ORDER, WINDOWS
from .data import running_mean, standardise, tendency
from .events import LIFECYCLE_LEN, lifecycle_matrix

#: budget key -> design-matrix column
PREDICTORS = {
    "damp": "P", "nino": "N", "olr": "C", "dnino": "dN", "pdo": "D",
    "tna": "T", "ao": "A", "nao": "M", "nl": "NL",
}
PRED_ORDER = ["damp", "nino", "olr", "dnino", "pdo", "tna", "ao", "nao", "nl"]

BASE_COLS = ["PNA", "NAO", "AO", "PDO", "TNA", "OLR"]


def build_design(df: pd.DataFrame, window: int, nino_col: str = "NINO34",
                 flip_nino: bool = False, base_cols: list | None = None) -> pd.DataFrame:
    """Smoothed, standardised design matrix for one integration window."""
    cols = base_cols or BASE_COLS
    X = standardise(df[cols + [nino_col]].copy())
    X = X.rename(columns={"PNA": "P", "NAO": "M", "AO": "A", "PDO": "D",
                          "TNA": "T", "OLR": "C", nino_col: "N"})
    if flip_nino:                      # SOI inverted to an El Nino-positive sense
        X["N"] = -X["N"]
    X = running_mean(X, window)
    X["dN"] = X["N"].shift(-1) - X["N"]
    # nonlinear ENSO forcing, orthogonalised against the linear term
    nl = X["N"] * np.abs(X["N"])
    ok = np.isfinite(X["N"]) & np.isfinite(nl)
    b = np.polyfit(X["N"][ok], nl[ok], 1)
    X["NL"] = nl - (b[0] * X["N"] + b[1])
    X["dP"] = X["P"].shift(-1) - X["P"]
    X["tend_c"] = tendency(X["P"])     # centred tendency, display only
    return X


def _harmonics(months: np.ndarray):
    return np.cos(2 * np.pi * months / 12.0), np.sin(2 * np.pi * months / 12.0)


def fit_budget(df: pd.DataFrame, window: int, nino_col: str = "NINO34",
               seasonal: bool = True, flip_nino: bool = False,
               predictors: list | None = None, base_cols: list | None = None):
    """Fit the tendency budget; return (info, terms, coefficients, design).

    ``predictors`` restricts the budget to a subset of the standard terms, which
    is required for the extended (pre-satellite) record where the convective and
    tropical Atlantic indices do not exist.
    """
    keys = predictors or PRED_ORDER
    X = build_design(df, window, nino_col, flip_nino, base_cols)
    names = [PREDICTORS[k] for k in keys]
    sub = X[names + ["dP"]].dropna()
    m = sub.index.month.values
    cos_m, sin_m = _harmonics(m)

    blocks, labels = [], []
    for n in names:
        v = sub[n].values
        blocks.append(v)
        labels.append((n, "0"))
        if seasonal:
            blocks += [v * cos_m, v * sin_m]
            labels += [(n, "c"), (n, "s")]
    A = np.column_stack(blocks)
    y = sub["dP"].values
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    cf = dict(zip(labels, beta))

    mm = X.index.month.values
    c_all, s_all = _harmonics(mm)
    terms = pd.DataFrame(index=X.index)
    coefs = pd.DataFrame(index=X.index)
    for key in keys:
        n = PREDICTORS[key]
        a = np.full(len(X), cf[(n, "0")])
        if seasonal:
            a = a + cf[(n, "c")] * c_all + cf[(n, "s")] * s_all
        coefs[key] = a
        terms[key] = a * X[n].values
    terms["res"] = X["dP"] - terms[keys].sum(axis=1)
    for missing in [k for k in PRED_ORDER if k not in keys]:
        terms[missing] = 0.0
    terms["tend"] = X["dP"]
    for extra in ("P", "N", "C"):
        if extra in X:
            terms[extra] = X[extra]

    ok = np.isfinite(terms["res"]) & np.isfinite(X["dP"])
    info = {
        "coef_annual": {k: cf[(PREDICTORS[k], "0")] for k in keys},
        "lambda_annual": -cf[("P", "0")],
        "r2": 1.0 - np.nanvar(terms["res"][ok]) / np.nanvar(X["dP"][ok]),
        "n": int(ok.sum()), "seasonal": seasonal, "window": window,
    }
    return info, terms, coefs, X


# ---------------------------------------------------------------------------
# Phase integration
# ---------------------------------------------------------------------------

def phase_windows(scheme: str, peak_month: int = 11):
    """(development months, decay months) as life-cycle month indices."""
    if scheme == "key":
        return list(KEY_DEV_MONTHS), list(KEY_DEC_MONTHS)
    if scheme == "phase9":
        d0 = max(peak_month - 8, 0)
        return (list(range(d0, peak_month + 1)),
                list(range(peak_month + 1, min(peak_month + 10, LIFECYCLE_LEN))))
    if scheme == "aligned":
        d0 = max(peak_month - 8, 0)
        return (list(range(d0, peak_month + 1)),
                list(range(peak_month + 1, min(peak_month + 8, LIFECYCLE_LEN))))
    raise ValueError(scheme)


def integrate_terms(terms: pd.DataFrame, years: list[int], scheme: str = "key",
                    peak_months: dict | None = None):
    """Phase-integrated contribution of every term, per event (units: sigma)."""
    keys = TERM_ORDER + ["tend"]
    dev = pd.DataFrame(index=years, columns=keys, dtype=float)
    dec = pd.DataFrame(index=years, columns=keys, dtype=float)
    for k in keys:
        M = lifecycle_matrix(terms[k], years)
        for i, y in enumerate(years):
            pk = 11 if peak_months is None else int(peak_months.get(y, 11))
            di, ci = phase_windows(scheme, pk)
            dev.loc[y, k] = np.nansum(M[i, di]) if np.isfinite(M[i, di]).any() else np.nan
            dec.loc[y, k] = np.nansum(M[i, ci]) if np.isfinite(M[i, ci]).any() else np.nan
    return dev, dec


def asyi(dev: pd.DataFrame, dec: pd.DataFrame) -> pd.Series:
    """Asymmetry index: net contribution over the development and decay phases."""
    return dev.mean() + dec.mean()


def aae(dev: pd.DataFrame, dec: pd.DataFrame) -> pd.Series:
    """Asymmetric attribution efficiency, bounded in [-1, 1]."""
    a, b = dev.mean(), dec.mean()
    return (a + b) / (np.abs(a) + np.abs(b))


def window_scan(df: pd.DataFrame, years: list[int], nino_col: str = "NINO34",
                scheme: str = "key", windows=None, seasonal: bool = True,
                flip_nino: bool = False, peak_months=None):
    """Repeat the full diagnosis for every integration window."""
    windows = list(windows or WINDOWS)
    out = {"windows": windows, "dev": {}, "dec": {}, "info": {}, "terms": {},
           "coefs": {}}
    for k in windows:
        info, terms, coefs, _ = fit_budget(df, k, nino_col, seasonal, flip_nino)
        dev, dec = integrate_terms(terms, years, scheme, peak_months)
        out["dev"][k], out["dec"][k] = dev, dec
        out["info"][k], out["terms"][k], out["coefs"][k] = info, terms, coefs
    out["asyi_df"] = pd.DataFrame({k: asyi(out["dev"][k], out["dec"][k])
                                   for k in windows})
    out["aae_df"] = pd.DataFrame({k: aae(out["dev"][k], out["dec"][k])
                                  for k in windows})
    return out


def wti(scan: dict, phase: str = "dev") -> pd.Series:
    """Window trend index: change of a phase contribution across the window range.

    The phase integral of each term is regressed linearly on the integration
    window and the slope is rescaled by the full window range, exactly as the
    vertical trend index of Song et al. (2026) rescales the depth regression by
    the 30-100 m range.  Positive values indicate weakening with window length.
    """
    ws = np.array(scan["windows"], dtype=float)
    src = scan[phase]
    res = {}
    for t in TERM_ORDER + ["tend"]:
        y = np.array([src[k][t].mean() for k in scan["windows"]], dtype=float)
        ok = np.isfinite(y)
        res[t] = (np.polyfit(ws[ok], y[ok], 1)[0] * (ws.max() - ws.min())
                  if ok.sum() > 2 else np.nan)
    return pd.Series(res)


def closure_metric(scan: dict) -> pd.Series:
    """Normalised key-period closure error |I_res| / |I_tend| for every window.

    A scale-free criterion is required because the raw residual shrinks
    trivially as the smoothing window lengthens.
    """
    rec = {}
    for k in scan["windows"]:
        dev, dec = scan["dev"][k], scan["dec"][k]
        num = abs(dev["res"].mean()) + abs(dec["res"].mean())
        den = abs(dev["tend"].mean()) + abs(dec["tend"].mean())
        rec[k] = num / den
    return pd.Series(rec)


# ---------------------------------------------------------------------------
# Null-model decomposition of the coupling asymmetry
# ---------------------------------------------------------------------------

def coupling_decomposition(df: pd.DataFrame, years: list[int], window: int,
                           nino_col: str = "NINO34"):
    """Split the AsyI of each term into amplitude and seasonal components.

    Following Eq. (9) of Song et al. (2026), with the damping coefficient
    generalised to the full set of budget coefficients:

        AsyI[F_j] = a_j,bar AsyI[X_j]                                    (I)
                    + (a_j,dev - a_j,bar) sum(X_j,dev)
                    + (a_j,dec - a_j,bar) sum(X_j,dec)                  (II)
                    + epsilon                                           (III)

    Term I is the amplitude-driven component: what a season-invariant coupling
    would produce given the observed asymmetry of the forcing.  Term II
    isolates the seasonal (phase-dependent) modulation of the coupling.
    """
    info, terms, coefs, X = fit_budget(df, window, nino_col, seasonal=True)
    dev_i, dec_i = list(KEY_DEV_MONTHS), list(KEY_DEC_MONTHS)
    out = {}
    for key in PRED_ORDER:
        name = PREDICTORS[key]
        Xm = lifecycle_matrix(X[name], years)
        Am = lifecycle_matrix(coefs[key], years)
        Fm = lifecycle_matrix(terms[key], years)

        a_bar = np.nanmean(Am[:, dev_i + dec_i])
        a_dev = np.nanmean(Am[:, dev_i])
        a_dec = np.nanmean(Am[:, dec_i])
        sX_dev = np.nanmean(np.nansum(Xm[:, dev_i], axis=1))
        sX_dec = np.nanmean(np.nansum(Xm[:, dec_i], axis=1))
        asyi_X = sX_dev + sX_dec
        asyi_F = (np.nanmean(np.nansum(Fm[:, dev_i], axis=1))
                  + np.nanmean(np.nansum(Fm[:, dec_i], axis=1)))
        t1 = a_bar * asyi_X
        t2 = (a_dev - a_bar) * sX_dev + (a_dec - a_bar) * sX_dec
        out[key] = dict(a_bar=a_bar, a_dev=a_dev, a_dec=a_dec,
                        asyi_X=asyi_X, asyi_F=asyi_F,
                        amp=t1, seasonal=t2, residual=asyi_F - t1 - t2)
    return pd.DataFrame(out).T


def phase_regression(xs: pd.Series, ys: pd.Series, years: list[int]):
    """Event-by-event phase-mean regression of ys on xs for both key periods."""
    Xm = lifecycle_matrix(xs, years)
    Ym = lifecycle_matrix(ys, years)
    out = {}
    for name, idx in (("dev", list(KEY_DEV_MONTHS)), ("dec", list(KEY_DEC_MONTHS))):
        x = np.nanmean(Xm[:, idx], axis=1)
        y = np.nanmean(Ym[:, idx], axis=1)
        ok = np.isfinite(x) & np.isfinite(y)
        slope, icept = np.polyfit(x[ok], y[ok], 1)
        r = float(np.corrcoef(x[ok], y[ok])[0, 1])
        n = int(ok.sum())
        t = abs(r) * np.sqrt((n - 2) / max(1e-12, 1 - r ** 2))
        out[name] = dict(x=x, y=y, slope=slope, intercept=icept, r=r,
                         p=float(2 * (1 - _st.t.cdf(t, n - 2))), n=n)
    return out


def seasonal_coupling(df: pd.DataFrame, window: int, xkey: str, ykey: str,
                      nino_col: str = "NINO34"):
    """Calendar-month regression slope and correlation of ykey on xkey."""
    X = build_design(df, window, nino_col)
    rec = []
    for m in range(1, 13):
        d = X.loc[X.index.month == m, [xkey, ykey]].dropna()
        slope = np.polyfit(d[xkey], d[ykey], 1)[0]
        r = float(np.corrcoef(d[xkey], d[ykey])[0, 1])
        n = len(d)
        t = abs(r) * np.sqrt((n - 2) / max(1e-12, 1 - r ** 2))
        rec.append(dict(month=m, slope=slope, r=r, n=n,
                        p=float(2 * (1 - _st.t.cdf(t, n - 2)))))
    return pd.DataFrame(rec).set_index("month")


def select_reference_window(df: pd.DataFrame, years: list[int],
                            windows=None, max_window: int = 6):
    """Objective choice of the reference integration window.

    Two criteria are applied, in order:

    1. **Separation.** The smoothing windows centred on the two key periods must
       remain disjoint, which requires k <= 6 months for key periods whose
       centres (Sep0 and Mar1) are six months apart.
    2. **Closure.** Among the admissible windows for which the diagnosed
       tendency asymmetry is significant at the 95% level, the window minimising
       the normalised closure mismatch |AsyI_res| / |AsyI_tend| is selected.

    Returns (reference window, diagnostic table).
    """
    windows = list(windows or WINDOWS)
    scan = window_scan(df, years, windows=windows)
    rows = []
    for k in windows:
        tend = (scan["dev"][k]["tend"] + scan["dec"][k]["tend"]).dropna()
        res = (scan["dev"][k]["res"] + scan["dec"][k]["res"]).dropna()
        p = float(_st.ttest_1samp(tend, 0).pvalue)
        rows.append(dict(window=k, asyi_tend=tend.mean(), p_tend=p,
                         asyi_res=res.mean(),
                         mismatch=abs(res.mean()) / abs(tend.mean()),
                         admissible=(k <= max_window) and (p < 0.05)))
    tab = pd.DataFrame(rows).set_index("window")
    ok = tab[tab.admissible]
    ref = int(ok["mismatch"].idxmin()) if len(ok) else windows[0]
    return ref, tab
