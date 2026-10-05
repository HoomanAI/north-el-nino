"""
Figure S1. Sensitivity of the diagnosed asymmetry to the choice of ENSO index
and to the inclusion of the extreme events.

a  Composite PNA life cycle for event sets defined from Nino1+2, Nino3, Nino3.4,
   Nino4 and the inverted SOI.
b  Tendency asymmetry for each event set, with 95% confidence intervals.
c  Composite PNA with and without the three extreme events (1982, 1997, 2015).
d  Tendency asymmetry for the same two samples.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from enso_na.config import OUT_DIR, REF_WINDOW
from enso_na.data import load_monthly, running_mean, standardise
from enso_na.events import (EXCLUDED, ci95, composite, default_events,
                            identify_el_nino)
from enso_na.budget import fit_budget, integrate_terms
from enso_na.plotting import lifecycle_axis, panel_label, save, shade_key_periods

INDICES = [("NINO12", False), ("NINO3", False), ("NINO34", False),
           ("NINO4", False), ("SOI", True)]
COLS = ["#8c6d31", "#1b9e77", "#c1272d", "#1f4e9c", "#e08214"]
EXTREME = (1982, 1997, 2015)


def _asy(terms, years):
    dev, dec = integrate_terms(terms, years, "key")
    s = (dev["tend"] + dec["tend"]).dropna()
    return s.mean(), ci95(s.values), len(s)


def main():
    df = load_monthly()
    _, terms, _, _ = fit_budget(df, REF_WINDOW)
    Xs = running_mean(standardise(df[["PNA"]]), REF_WINDOW)
    x = np.arange(24)

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.0))
    fig.subplots_adjust(hspace=0.5, wspace=0.3)

    rows = []
    ax = axes[0, 0]
    shade_key_periods(ax)
    for (name, flip), col in zip(INDICES, COLS):
        s = -df[name] if flip else df[name]
        yrs = identify_el_nino(s if not flip else (s - s.mean()) / s.std(),
                               threshold=0.5, exclude=EXCLUDED)
        m, ci, _ = composite(Xs["PNA"], yrs)
        ax.plot(x, m, color=col, lw=1.2, label=f"{name} (n={len(yrs)})")
        a, e, n = _asy(terms, yrs)
        rows.append(dict(index=name, n=n, asyi=a, ci=e, years=yrs))
    lifecycle_axis(ax)
    ax.set_ylabel("PNA ($\\sigma$)")
    ax.set_title("Composite PNA by ENSO index", fontsize=8.5)
    ax.legend(ncol=2, fontsize=5.8, loc="upper left")
    panel_label(ax, "a")

    ax = axes[0, 1]
    tab = pd.DataFrame(rows)
    idx = np.arange(len(tab))
    ax.bar(idx, tab.asyi, 0.6, color=COLS, lw=0)
    ax.errorbar(idx, tab.asyi, yerr=tab.ci, fmt="none", ecolor="0.25",
                elinewidth=0.6, capsize=2)
    ax.axhline(0, color="0.4", lw=0.6)
    ax.set_xticks(idx); ax.set_xticklabels(tab["index"], rotation=45, fontsize=6.5)
    ax.set_ylabel("AsyI[$\\partial P/\\partial t$] ($\\sigma$)")
    ax.set_title("Tendency asymmetry by ENSO index", fontsize=8.5)
    panel_label(ax, "b")

    ev = default_events(df)
    noex = [y for y in ev if y not in EXTREME]
    ax = axes[1, 0]
    shade_key_periods(ax)
    for yrs, col, lab in ((ev, "#c1272d", f"all events (n={len(ev)})"),
                          (noex, "#1f4e9c", f"extremes removed (n={len(noex)})")):
        m, ci, _ = composite(Xs["PNA"], yrs)
        ax.plot(x, m, color=col, lw=1.3, label=lab)
        ax.fill_between(x, m - ci, m + ci, color=col, alpha=0.13, lw=0)
    lifecycle_axis(ax)
    ax.set_ylabel("PNA ($\\sigma$)")
    ax.set_title("Effect of the extreme events", fontsize=8.5)
    ax.legend(loc="upper left")
    panel_label(ax, "c")

    ax = axes[1, 1]
    vals = [_asy(terms, ev), _asy(terms, noex)]
    ax.bar([0, 1], [v[0] for v in vals], 0.5, color=["#c1272d", "#1f4e9c"], lw=0)
    ax.errorbar([0, 1], [v[0] for v in vals], yerr=[v[1] for v in vals],
                fmt="none", ecolor="0.25", elinewidth=0.6, capsize=2)
    ax.axhline(0, color="0.4", lw=0.6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["all events", "extremes removed"], fontsize=7)
    ax.set_ylabel("AsyI[$\\partial P/\\partial t$] ($\\sigma$)")
    ax.set_title("Tendency asymmetry", fontsize=8.5)
    panel_label(ax, "d")

    save(fig, "figS01")
    tab.round(3).to_csv(OUT_DIR / "figS01_index_sensitivity.csv", index=False)
    print(tab.round(3).to_string(index=False))
    print("all events %.3f +- %.3f | extremes removed %.3f +- %.3f"
          % (vals[0][0], vals[0][1], vals[1][0], vals[1][1]))


if __name__ == "__main__":
    main()
