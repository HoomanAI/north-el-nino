"""
Figure S4. Seasonal structure of the fitted budget coefficients.

a  Seasonal cycle of every budget coefficient at the reference window.
b  Coefficient of the damping term for each integration window.
c  Asymmetric attribution efficiency of every term across integration windows.
d  Amplitude and seasonal components of the AsyI, term by term.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from enso_na.config import OUT_DIR, REF_WINDOW, TERMS, WINDOWS
from enso_na.data import load_monthly
from enso_na.events import default_events
from enso_na.budget import (PRED_ORDER, coupling_decomposition, fit_budget,
                            window_scan)
from enso_na.plotting import panel_label, save

MON = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]


def main():
    df = load_monthly()
    ev = default_events(df)
    scan = window_scan(df, ev, windows=WINDOWS)

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.2))
    fig.subplots_adjust(hspace=0.52, wspace=0.34)

    # ---- a ------------------------------------------------------------------
    ax = axes[0, 0]
    _, _, coefs, _ = fit_budget(df, REF_WINDOW)
    seas = coefs.groupby(coefs.index.month).mean()
    for t in PRED_ORDER:
        ax.plot(range(1, 13), seas[t], color=TERMS[t][1], lw=1.1, marker="o",
                ms=2.2, label=TERMS[t][0])
    ax.axhline(0, color="0.4", lw=0.6)
    for mm, c in (([8, 9, 10], "#4daf4a"), ([2, 3, 4], "#984ea3")):
        ax.axvspan(min(mm) - 0.5, max(mm) + 0.5, color=c, alpha=0.10, lw=0)
    ax.set_xticks(range(1, 13)); ax.set_xticklabels(MON, fontsize=6.5)
    ax.set_ylabel("coefficient (month$^{-1}$)")
    ax.set_title(f"All coefficients, $k$ = {REF_WINDOW} months", fontsize=8.5)
    ax.legend(ncol=3, fontsize=5.2, loc="lower center")
    panel_label(ax, "a")

    # ---- b ------------------------------------------------------------------
    ax = axes[0, 1]
    for k in WINDOWS:
        c = scan["coefs"][k]
        s = c.groupby(c.index.month).mean()["damp"]
        ax.plot(range(1, 13), -s.values, lw=1.1, marker="o", ms=2.5,
                color=plt.cm.viridis((k - 1) / 6), label=f"$k$={k}")
    ax.axhline(0, color="0.4", lw=0.6)
    ax.set_xticks(range(1, 13)); ax.set_xticklabels(MON, fontsize=6.5)
    ax.set_ylabel("$\\lambda$ (month$^{-1}$)")
    ax.set_title("Seasonal cycle of the damping rate", fontsize=8.5)
    ax.legend(ncol=2, fontsize=6)
    panel_label(ax, "b")

    # ---- c ------------------------------------------------------------------
    ax = axes[1, 0]
    A = scan["aae_df"]
    for t in PRED_ORDER:
        ax.plot(WINDOWS, A.loc[t, WINDOWS].values.astype(float),
                color=TERMS[t][1], lw=1.1, marker="o", ms=2.5, label=TERMS[t][0])
    ax.axhline(0, color="0.4", lw=0.6)
    ax.set_ylim(-1.15, 1.15)
    ax.set_xlabel("integration window $k$ (months)")
    ax.set_ylabel("AAE")
    ax.set_title("Attribution efficiency across windows", fontsize=8.5)
    panel_label(ax, "c")

    # ---- d ------------------------------------------------------------------
    ax = axes[1, 1]
    dec = coupling_decomposition(df, ev, REF_WINDOW).astype(float)
    idx = np.arange(len(PRED_ORDER))
    ax.bar(idx - 0.2, dec.loc[PRED_ORDER, "amp"], 0.38, color="#c1272d", lw=0,
           label="amplitude (I)")
    ax.bar(idx + 0.2, dec.loc[PRED_ORDER, "seasonal"], 0.38, color="#1f4e9c",
           lw=0, label="seasonal (II)")
    ax.axhline(0, color="0.4", lw=0.6)
    ax.set_xticks(idx)
    ax.set_xticklabels([TERMS[t][0] for t in PRED_ORDER], rotation=90, fontsize=6)
    ax.set_ylabel("AsyI component ($\\sigma$)")
    ax.set_title("Amplitude versus seasonal, by term", fontsize=8.5)
    ax.legend(loc="lower left")
    panel_label(ax, "d")

    save(fig, "figS04")
    seas.round(3).to_csv(OUT_DIR / "figS04_seasonal_coefficients.csv")
    A.round(3).to_csv(OUT_DIR / "figS04_aae_windows.csv")
    print(seas.round(3).to_string())


if __name__ == "__main__":
    main()
