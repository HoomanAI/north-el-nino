"""
Figure S5. Sensitivity of the attribution to the phase-integration scheme.

a  AsyI of every term under the three integration schemes at the reference window.
b-d  AsyI of the tendency and of the leading terms as a function of the
     integration window, for each of the three schemes.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from enso_na.config import OUT_DIR, REF_WINDOW, TERMS, TERM_ORDER, WINDOWS
from enso_na.data import load_monthly
from enso_na.events import ci95, default_events
from enso_na.budget import fit_budget, integrate_terms
from enso_na.plotting import panel_label, save

SCHEMES = [("key", "ASO0 / FMA1", "#1f4e9c"),
           ("phase9", "9-month phases", "#e08214"),
           ("aligned", "9/7-month phases", "#4daf4a")]
SHOW = ["damp", "nino", "olr", "dnino", "pdo", "tna", "nao", "nl", "res"]


def main():
    df = load_monthly()
    ev = default_events(df)

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.2))
    fig.subplots_adjust(hspace=0.55, wspace=0.32)
    axes = axes.ravel()

    # ---- a ------------------------------------------------------------------
    ax = axes[0]
    _, terms, _, _ = fit_budget(df, REF_WINDOW)
    idx = np.arange(len(SHOW) + 1)
    rows = []
    for j, (sch, lab, col) in enumerate(SCHEMES):
        dv, dc = integrate_terms(terms, ev, sch)
        a = (dv + dc).mean()
        keys = SHOW + ["tend"]
        ax.bar(idx + (j - 1) * 0.27, [a[t] for t in keys], 0.25, color=col,
               lw=0, label=lab)
        ax.errorbar(idx + (j - 1) * 0.27, [a[t] for t in keys],
                    yerr=[ci95((dv[t] + dc[t]).values) for t in keys],
                    fmt="none", ecolor="0.3", elinewidth=0.4, capsize=1.0)
        rows.append(dict(scheme=sch, **{t: a[t] for t in keys}))
    ax.axhline(0, color="0.4", lw=0.6)
    ax.set_xticks(idx)
    ax.set_xticklabels([TERMS[t][0] for t in SHOW] + ["$\\partial P/\\partial t$"],
                       rotation=90, fontsize=6)
    ax.set_ylabel("AsyI ($\\sigma$)")
    ax.set_title(f"AsyI by integration scheme, $k$ = {REF_WINDOW} months",
                 fontsize=8.5)
    ax.legend(loc="upper center", ncol=3, fontsize=5.5, handlelength=1.0, columnspacing=0.8)   # VERIFY (display only)
    yl = ax.get_ylim(); ax.set_ylim(yl[0], yl[1] + 0.22 * (yl[1] - yl[0]))
    panel_label(ax, "a")

    # ---- b, c, d ------------------------------------------------------------
    letters = "bcd"
    for j, (sch, lab, _) in enumerate(SCHEMES):
        ax = axes[j + 1]
        tab = {}
        for k in WINDOWS:
            _, tk, _, _ = fit_budget(df, k)
            dv, dc = integrate_terms(tk, ev, sch)
            tab[k] = (dv + dc).mean()
        T = pd.DataFrame(tab)
        for t in ["damp", "nino", "dnino", "pdo", "tna"]:
            ax.plot(WINDOWS, T.loc[t].values, color=TERMS[t][1], lw=1.1,
                    marker="o", ms=2.5, label=TERMS[t][0])
        ax.plot(WINDOWS, T.loc["tend"].values, color="k", lw=1.4, marker="s",
                ms=3, label="$\\partial P/\\partial t$")
        ax.axhline(0, color="0.4", lw=0.6)
        ax.set_xlabel("integration window $k$ (months)")
        ax.set_ylabel("AsyI ($\\sigma$)")
        ax.set_title(lab, fontsize=8.5)
        if j == 0:
            ax.legend(ncol=2, fontsize=5.8, loc="lower right")
        panel_label(ax, letters[j])

    save(fig, "figS05")
    pd.DataFrame(rows).round(3).to_csv(OUT_DIR / "figS05_schemes.csv", index=False)
    print(pd.DataFrame(rows).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
