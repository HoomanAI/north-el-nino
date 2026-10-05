"""
Figure 2. Tendency budget of the North American teleconnection during El Nino
development and decay.

a  Composite monthly evolution of every budget term over the El Nino life cycle
   at the reference integration window.
b  Key-period budgets (ASO0 versus FMA1) for integration windows of 1-7 months.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from enso_na.config import (OUT_DIR, REF_WINDOW, TERMS, TERM_ORDER, TERM_TEXT,
                            WINDOWS)
from enso_na.data import load_monthly
from enso_na.events import ci95, composite, default_events
from enso_na.budget import fit_budget, integrate_terms, window_scan
from enso_na.plotting import lifecycle_axis, panel_label, save, shade_key_periods

ROWS = ["tend"] + TERM_ORDER


def main():
    df = load_monthly()
    ev = default_events(df)
    info, terms, coefs, X = fit_budget(df, REF_WINDOW)
    scan = window_scan(df, ev, windows=WINDOWS)

    fig = plt.figure(figsize=(7.2, 8.6))
    gs = fig.add_gridspec(len(ROWS), 2, wspace=0.26, hspace=0.18,
                          width_ratios=[1.05, 1.0])
    x = np.arange(24)

    for i, key in enumerate(ROWS):
        lab = "$\\partial P/\\partial t$" if key == "tend" else TERMS[key][0]
        col = "k" if key == "tend" else TERMS[key][1]

        # --- left column: composite monthly evolution -----------------------
        ax = fig.add_subplot(gs[i, 0])
        m, ci, _ = composite(terms[key], ev)
        pos = np.where(m >= 0, m, 0)
        neg = np.where(m < 0, m, 0)
        ax.bar(x, pos, color="#c1272d", width=0.8, lw=0)
        ax.bar(x, neg, color="#1f4e9c", width=0.8, lw=0)
        ax.errorbar(x, m, yerr=ci, fmt="none", ecolor="0.3", elinewidth=0.5,
                    capsize=1.0)
        shade_key_periods(ax, alpha=0.10)
        lifecycle_axis(ax)
        lim_a = max(0.65, float(np.nanmax(np.abs(m) + np.nan_to_num(ci))) * 1.08)
        ax.set_ylim(-lim_a, lim_a)   # VERIFY (display only): no clipping
        ax.set_yticks([-0.5, 0, 0.5])
        ax.set_ylabel(lab, rotation=0, ha="right", va="center", labelpad=14)
        if i < len(ROWS) - 1:
            ax.set_xticklabels([])
        if i == 0:
            ax.set_title(f"Monthly budget, window $k$ = {REF_WINDOW} months",
                         fontsize=8.5)
            panel_label(ax, "a", dx=-0.22)

        # --- right column: key-period bars across windows -------------------
        ax = fig.add_subplot(gs[i, 1])
        w = np.arange(len(WINDOWS))
        dv = [scan["dev"][k][key].mean() for k in WINDOWS]
        dc = [scan["dec"][k][key].mean() for k in WINDOWS]
        ev_d = [ci95(scan["dev"][k][key].values) for k in WINDOWS]
        ev_c = [ci95(scan["dec"][k][key].values) for k in WINDOWS]
        ax.bar(w - 0.19, dv, 0.36, color="#4daf4a", lw=0, label="ASO0")
        ax.bar(w + 0.19, dc, 0.36, color="#984ea3", lw=0, label="FMA1")
        ax.errorbar(w - 0.19, dv, yerr=ev_d, fmt="none", ecolor="0.25",
                    elinewidth=0.5, capsize=1.0)
        ax.errorbar(w + 0.19, dc, yerr=ev_c, fmt="none", ecolor="0.25",
                    elinewidth=0.5, capsize=1.0)
        ax.axhline(0, color="0.4", lw=0.5)
        ax.set_xlim(-0.6, len(WINDOWS) - 0.4)
        vals = np.abs(np.r_[dv, dc]) + np.nan_to_num(np.r_[ev_d, ev_c])
        lim_b = max(1.0, float(np.nanmax(vals)) * 1.08)
        ax.set_ylim(-lim_b, lim_b)   # VERIFY (display only): no clipping
        ax.set_yticks([-1, 0, 1])
        ax.set_xticks(w)
        if i < len(ROWS) - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xticklabels([f"{k}m" for k in WINDOWS])
        ax.tick_params(labelsize=6.5)
        if i == 0:
            ax.set_title("Key-period budgets across windows", fontsize=8.5)
            ax.legend(loc="upper right", ncol=2, fontsize=6)
            panel_label(ax, "b", dx=-0.10)

    save(fig, "fig02")

    rows = []
    for k in WINDOWS:
        for key in ROWS:
            rows.append(dict(window=k, term=key, label=TERM_TEXT.get(key, "tendency"),
                             I_dev=scan["dev"][k][key].mean(),
                             I_dec=scan["dec"][k][key].mean()))
    pd.DataFrame(rows).to_csv(OUT_DIR / "fig02_phase_integrals.csv", index=False)
    print(pd.DataFrame(rows).query("window == @REF_WINDOW").round(3).to_string(index=False))


if __name__ == "__main__":
    main()
