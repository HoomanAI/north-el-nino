#!/usr/bin/env python3
"""
Table 7 - seasonal interaction test for the two links of the ENSO-North
America chain.

For each link, the regression of the response on the predictor is fitted
separately over the ASO0 and FMA1 calendar months, and the change in slope
between the two periods is estimated from the pooled sample with the
interaction model of Eq. (8):

    y = b0 + b1 x + b2 g + b3 (x g) + e,        g = 1 for FMA1, 0 for ASO0

so that b3 is the change in slope from ASO0 to FMA1. Regressions use every
year of the analysis period, not only El Nino years.

Note on inference: monthly samples within a three-month window are serially
dependent, so the reported p-values are nominal and are likely optimistic.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from enso_na.budget import build_design
from enso_na.config import OUT_DIR, REF_WINDOW
from enso_na.data import load_monthly

ASO0_MONTHS = (8, 9, 10)
FMA1_MONTHS = (2, 3, 4)

LINKS = [
    ("N", "C", "Tropical: Nino3.4 -> OLR"),
    ("C", "P", "Extratropical: OLR -> PNA"),
    ("N", "P", "Direct: Nino3.4 -> PNA"),
]


def _ols(x: np.ndarray, y: np.ndarray):
    """Slope, its standard error, and the two-sided p-value."""
    n = len(x)
    A = np.column_stack([np.ones(n), x])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ beta
    s2 = resid @ resid / (n - 2)
    se = np.sqrt(np.diag(s2 * np.linalg.inv(A.T @ A)))
    t = beta[1] / se[1]
    return beta[1], se[1], float(2 * (1 - stats.t.cdf(abs(t), n - 2))), n


def interaction(x: np.ndarray, y: np.ndarray, g: np.ndarray):
    """Eq. (8): return the interaction coefficient, its s.e., p and d.o.f."""
    A = np.column_stack([np.ones(len(x)), x, g, x * g])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ beta
    dof = len(x) - A.shape[1]
    s2 = resid @ resid / dof
    se = np.sqrt(np.diag(s2 * np.linalg.inv(A.T @ A)))
    t = beta[3] / se[3]
    return beta[3], se[3], float(2 * (1 - stats.t.cdf(abs(t), dof))), dof


def main(window: int = REF_WINDOW) -> pd.DataFrame:
    X = build_design(load_monthly(), window)
    rows = []
    for xk, yk, label in LINKS:
        d = X[[xk, yk]].dropna()
        aso = d[d.index.month.isin(ASO0_MONTHS)]
        fma = d[d.index.month.isin(FMA1_MONTHS)]

        b_aso, se_aso, p_aso, n_aso = _ols(aso[xk].values, aso[yk].values)
        b_fma, se_fma, p_fma, n_fma = _ols(fma[xk].values, fma[yk].values)

        x = np.concatenate([aso[xk].values, fma[xk].values])
        y = np.concatenate([aso[yk].values, fma[yk].values])
        g = np.concatenate([np.zeros(n_aso), np.ones(n_fma)])
        dbeta, se_d, p_d, dof = interaction(x, y, g)

        rows.append(dict(link=label, predictor=xk, response=yk,
                         slope_ASO0=round(b_aso, 3), se_ASO0=round(se_aso, 3),
                         p_ASO0=p_aso, n_ASO0=n_aso,
                         slope_FMA1=round(b_fma, 3), se_FMA1=round(se_fma, 3),
                         p_FMA1=p_fma, n_FMA1=n_fma,
                         delta_beta=round(dbeta, 3), se_delta=round(se_d, 3),
                         p_delta=p_d, n_pooled=len(x), dof=dof))

    tab = pd.DataFrame(rows)
    out = OUT_DIR / "table07_interaction.csv"
    tab.to_csv(out, index=False)
    print(tab.to_string(index=False))
    print(f"  wrote {out}")
    return tab


if __name__ == "__main__":
    main()
