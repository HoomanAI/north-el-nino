#!/usr/bin/env python3
"""
Regression test for the numbers quoted in the paper.

Recomputes the headline diagnostics from the shipped index files and asserts
that each matches the published value. Run it after any change to the
pipeline:

    python verify_numbers.py

Exit status is 0 if every claim matches, 1 otherwise.
"""
from __future__ import annotations

import sys

import numpy as np
from scipy import stats

from enso_na.budget import window_scan
from enso_na.config import KEY_DEC_MONTHS, KEY_DEV_MONTHS, REF_WINDOW
from enso_na.data import load_monthly, running_mean, standardise
from enso_na.events import ci95, default_events, lifecycle_matrix, phase_speeds

TOL = 0.005

CHECKS: list[tuple[str, float, float]] = []


def check(name: str, got: float, expected: float, tol: float = TOL) -> None:
    CHECKS.append((name, got, expected, tol))


def main() -> int:
    df = load_monthly()
    events = default_events(df)
    check("number of events", len(events), 10, 0)

    X = running_mean(standardise(df[["PNA", "NINO34"]]), REF_WINDOW)
    sp = phase_speeds(X["NINO34"], events)
    t = stats.ttest_rel(sp.dev_speed, sp.dec_speed)
    check("Nino3.4 development speed", sp.dev_speed.mean(), 0.22, 0.005)
    check("Nino3.4 decay speed", sp.dec_speed.mean(), 0.28, 0.005)
    check("paired speed p-value", t.pvalue, 0.033, 0.001)

    N = lifecycle_matrix(X["NINO34"], events)
    P = lifecycle_matrix(X["PNA"], events)
    dev, dec = list(KEY_DEV_MONTHS), list(KEY_DEC_MONTHS)
    check("Nino3.4 ASO0", np.nanmean(N[:, dev]), 1.25, 0.01)
    check("Nino3.4 FMA1", np.nanmean(N[:, dec]), 1.07, 0.01)
    check("PNA ASO0", np.nanmean(P[:, dev]), 0.05, 0.01)
    check("PNA FMA1", np.nanmean(P[:, dec]), 0.41, 0.01)

    scan = window_scan(df, events)
    a = (scan["dev"][REF_WINDOW]["tend"] + scan["dec"][REF_WINDOW]["tend"]).values
    _, p = stats.ttest_1samp(a, 0)
    check("AsyI tendency", a.mean(), -1.04, 0.01)
    check("AsyI 95% CI", ci95(a), 0.50, 0.01)
    check("AsyI p-value", p, 0.001, 0.0005)
    check("events with negative AsyI", int((a < 0).sum()), 9, 0)

    for term, expected in (("damp", -0.26), ("nino", 0.57), ("dnino", -0.46),
                           ("res", -0.53)):
        v = (scan["dev"][REF_WINDOW][term] + scan["dec"][REF_WINDOW][term]).mean()
        check(f"AsyI[{term}]", v, expected, 0.01)

    ok = True
    width = max(len(n) for n, *_ in CHECKS)
    for name, got, exp, tol in CHECKS:
        good = abs(got - exp) <= tol
        ok &= good
        print(f"{'PASS' if good else 'FAIL'}  {name:<{width}}  "
              f"computed {got:+.4f}   published {exp:+.4f}")
    print("\nAll checks passed." if ok else "\nSome checks FAILED.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
