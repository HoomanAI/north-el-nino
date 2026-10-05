"""
Figure S2. The extended 1950-2025 record.

The convective (OLR) and tropical North Atlantic indices start in the early
1980s, so the extended record supports only a reduced budget driven by the
inverted SOI, the PDO, the AO and the NAO.  Events are identified from the
standardised inverted SOI with the same DJF threshold.

a  Composite PNA life cycle, extended versus satellite-era samples.
b  Reduced-budget key-period integrals for the extended sample.
c  AsyI of the reduced budget across integration windows.
d  Tendency asymmetry of the individual extended-record events.
"""
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

from enso_na.config import OUT_DIR, REF_WINDOW, TERMS, WINDOWS
from enso_na.data import (load_monthly, load_monthly_extended, running_mean,
                          standardise)
from enso_na.events import (EXCLUDED, ci95, composite, default_events,
                            identify_el_nino)
from enso_na.budget import fit_budget, integrate_terms
from enso_na.plotting import lifecycle_axis, panel_label, save, shade_key_periods

RED_TERMS = ["damp", "nino", "dnino", "pdo", "ao", "nao", "nl"]
BASE = ["PNA", "NAO", "AO", "PDO"]


def _fit(df, k):
    return fit_budget(df, k, nino_col="SOI", flip_nino=True,
                      predictors=RED_TERMS, base_cols=BASE)


def main():
    dfe = load_monthly_extended()
    enso = -(dfe["SOI"] - dfe["SOI"].mean()) / dfe["SOI"].std()
    ev_e = identify_el_nino(enso, threshold=0.5, exclude=EXCLUDED)
    df = load_monthly()
    ev_s = default_events(df)

    info, terms, coefs, X = _fit(dfe, REF_WINDOW)
    dev, dec = integrate_terms(terms, ev_e, "key")
    asy_e = (dev["tend"] + dec["tend"]).dropna()

    Xs = running_mean(standardise(dfe[["PNA"]]), REF_WINDOW)
    x = np.arange(24)

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.0))
    fig.subplots_adjust(hspace=0.5, wspace=0.32)

    ax = axes[0, 0]
    shade_key_periods(ax)
    for yrs, col, lab in ((ev_e, "#1f4e9c", f"1950-2025, SOI events (n={len(ev_e)})"),
                          (ev_s, "#c1272d", f"1982-2025, Nino3.4 events (n={len(ev_s)})")):
        m, ci, _ = composite(Xs["PNA"], yrs)
        ax.plot(x, m, color=col, lw=1.3, label=lab)
        ax.fill_between(x, m - ci, m + ci, color=col, alpha=0.13, lw=0)
    lifecycle_axis(ax)
    ax.set_ylabel("PNA ($\\sigma$)")
    ax.set_title("Composite PNA, extended record", fontsize=8.5)
    ax.legend(loc="upper left", fontsize=5.8)
    panel_label(ax, "a")

    ax = axes[0, 1]
    show = RED_TERMS + ["res", "tend"]
    idx = np.arange(len(show))
    dv = [dev[t].mean() for t in show]
    dc = [dec[t].mean() for t in show]
    ax.bar(idx - 0.2, dv, 0.38, color="#4daf4a", lw=0, label="ASO0")
    ax.bar(idx + 0.2, dc, 0.38, color="#984ea3", lw=0, label="FMA1")
    ax.errorbar(idx - 0.2, dv, yerr=[ci95(dev[t].values) for t in show],
                fmt="none", ecolor="0.25", elinewidth=0.5, capsize=1.2)
    ax.errorbar(idx + 0.2, dc, yerr=[ci95(dec[t].values) for t in show],
                fmt="none", ecolor="0.25", elinewidth=0.5, capsize=1.2)
    ax.axhline(0, color="0.4", lw=0.6)
    ax.set_xticks(idx)
    ax.set_xticklabels([TERMS[t][0] if t != "tend" else "$\\partial P/\\partial t$"
                        for t in show], rotation=90, fontsize=6)
    ax.set_ylabel("phase integral ($\\sigma$)")
    ax.set_title("Reduced budget, key periods", fontsize=8.5)
    ax.legend(loc="lower left")
    panel_label(ax, "b")

    ax = axes[1, 0]
    rows = []
    for k in WINDOWS:
        i2, t2, _, _ = _fit(dfe, k)
        d1, d2 = integrate_terms(t2, ev_e, "key")
        a = (d1 + d2).mean()
        rows.append(dict(window=k, r2=i2["r2"],
                         **{t: a[t] for t in show}))
    scan = pd.DataFrame(rows).set_index("window")
    for t in RED_TERMS:
        ax.plot(WINDOWS, scan[t], color=TERMS[t][1], lw=1.2, marker="o", ms=3,
                label=TERMS[t][0])
    ax.plot(WINDOWS, scan["tend"], color="k", lw=1.4, marker="s", ms=3,
            label="$\\partial P/\\partial t$")
    ax.axhline(0, color="0.4", lw=0.6)
    ax.set_xlabel("integration window $k$ (months)")
    ax.set_ylabel("AsyI ($\\sigma$)")
    ax.set_title("Window dependence, reduced budget", fontsize=8.5)
    ax.legend(ncol=3, fontsize=5.8, loc="lower right")
    panel_label(ax, "c")

    ax = axes[1, 1]
    s = asy_e.sort_values()
    ax.bar(np.arange(len(s)), s.values, 0.66, color="#1f4e9c", lw=0)
    ax.axhline(0, color="0.4", lw=0.6)
    ax.set_xticks(np.arange(len(s)))
    ax.set_xticklabels(s.index, rotation=90, fontsize=5.5)
    ax.set_ylabel("AsyI[$\\partial P/\\partial t$] ($\\sigma$)")
    p = stats.ttest_1samp(s.values, 0).pvalue
    ax.set_title(f"Individual events (mean {s.mean():+.2f}, $p$ = {p:.3f})",
                 fontsize=8.5)
    panel_label(ax, "d")

    save(fig, "figS02")
    scan.round(3).to_csv(OUT_DIR / "figS02_extended_record.csv")
    print(scan.round(3).to_string())
    print(f"extended events ({len(ev_e)}): {ev_e}")
    print(f"AsyI_tend = {s.mean():+.3f} +- {ci95(s.values):.3f}, p = {p:.4f}")


if __name__ == "__main__":
    main()
