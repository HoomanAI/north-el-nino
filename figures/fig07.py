"""
Figure 7. Diversity of the asymmetry across El Nino events.

a  Tendency asymmetry of every event, ordered by amplitude.
b  Tendency asymmetry versus peak Nino3.4 amplitude.
c  Composite PNA life cycle for strongly and weakly asymmetric events.
d  Composite Nino3.4 life cycle for the same two groups.
"""
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

from enso_na.config import OUT_DIR, REF_WINDOW
from enso_na.data import load_monthly, running_mean, standardise
from enso_na.events import composite, default_events, lifecycle_matrix
from enso_na.budget import fit_budget, integrate_terms
from enso_na.plotting import lifecycle_axis, panel_label, save, shade_key_periods


def main():
    df = load_monthly()
    ev = default_events(df)
    _, terms, _, X = fit_budget(df, REF_WINDOW)
    dev, dec = integrate_terms(terms, ev, "key")
    asy = (dev["tend"] + dec["tend"])

    Xs = running_mean(standardise(df[["PNA", "NINO34"]]), REF_WINDOW)
    N = lifecycle_matrix(Xs["NINO34"], ev)
    amp = np.nanmax(N[:, 6:15], axis=1)
    tab = pd.DataFrame({"amp": amp, "asyi": asy.values}, index=ev)
    tab = tab.sort_values("amp")

    r, p = stats.pearsonr(tab.amp, tab.asyi)
    med = np.median(asy.values)
    strong = [y for y in ev if asy[y] <= med]
    weak = [y for y in ev if asy[y] > med]

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.0))
    fig.subplots_adjust(hspace=0.52, wspace=0.30)

    # ---- a ------------------------------------------------------------------
    ax = axes[0, 0]
    idx = np.arange(len(tab))
    sc = ax.bar(idx, tab.asyi.values, 0.66,
                color=plt.cm.OrRd((tab.amp.values - 1.0) / 2.2), lw=0.4,
                edgecolor="0.3")
    ax.axhline(0, color="0.4", lw=0.6)
    ax.set_xticks(idx); ax.set_xticklabels(tab.index, rotation=90, fontsize=6)
    ax.set_ylabel("AsyI[$\\partial P/\\partial t$] ($\\sigma$)")
    ax.set_title("Asymmetry of individual events", fontsize=8.5)
    n_neg = int((tab.asyi < 0).sum())
    ax.text(0.03, 0.06, f"{n_neg}/{len(tab)} events negative",
            transform=ax.transAxes, fontsize=6.5)
    panel_label(ax, "a")

    # ---- b ------------------------------------------------------------------
    ax = axes[0, 1]
    ax.scatter(tab.amp, tab.asyi, c=tab.amp, cmap="OrRd", vmin=1.0, vmax=3.2,
               s=38, edgecolor="0.3", lw=0.5, zorder=3)
    for y, row in tab.iterrows():
        ax.annotate(str(y), (row.amp, row.asyi), fontsize=5.5,
                    xytext=(3, 3), textcoords="offset points", color="0.3")
    xs = np.linspace(tab.amp.min(), tab.amp.max(), 10)
    b = np.polyfit(tab.amp, tab.asyi, 1)
    ax.plot(xs, np.polyval(b, xs), color="0.3", lw=1.0, ls="--")
    ax.axhline(0, color="0.6", lw=0.5)
    ax.set_xlabel("peak Nino3.4 amplitude ($\\sigma$)")
    ax.set_ylabel("AsyI[$\\partial P/\\partial t$] ($\\sigma$)")
    ax.set_title(f"Amplitude dependence ($r$ = {r:.2f}, $p$ = {p:.2f})",
                 fontsize=8.5)
    panel_label(ax, "b")

    # ---- c, d ---------------------------------------------------------------
    x = np.arange(24)
    for ax, var, name in ((axes[1, 0], "PNA", "PNA"),
                          (axes[1, 1], "NINO34", "Nino3.4")):
        shade_key_periods(ax)
        for grp, col, lab in ((strong, "#c1272d", "strong asymmetry"),
                              (weak, "#1f4e9c", "weak asymmetry")):
            m, ci, _ = composite(Xs[var], grp)
            ax.plot(x, m, color=col, lw=1.4, label=f"{lab} (n={len(grp)})")
            ax.fill_between(x, m - ci, m + ci, color=col, alpha=0.13, lw=0)
        lifecycle_axis(ax)
        ax.set_ylabel(f"{name} ($\\sigma$)")
        ax.set_title(f"{name} composites by asymmetry group", fontsize=8.5)
    axes[1, 0].legend(loc="upper left")
    panel_label(axes[1, 0], "c")
    panel_label(axes[1, 1], "d")

    save(fig, "fig07")

    tab["group"] = ["strong" if y in strong else "weak" for y in tab.index]
    tab.round(3).to_csv(OUT_DIR / "fig07_event_diversity.csv")
    print(tab.round(3).to_string())
    print(f"corr(amp, AsyI) = {r:.2f} (p = {p:.3f})")
    print("strong:", strong, "\nweak:", weak)


if __name__ == "__main__":
    main()
