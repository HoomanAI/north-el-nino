"""
Figure 6. Window dependence of the diagnosis.

a  Phase integrals of the leading terms as a function of the integration window.
b  Objective selection of the reference window: normalised closure mismatch and
   the significance of the diagnosed tendency asymmetry.
c  Seasonal cycle of the fitted budget coefficients.
d  Variance of the tendency explained by the budget as a function of the window.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from enso_na.config import OUT_DIR, REF_WINDOW, TERMS, WINDOWS
from enso_na.data import load_monthly
from enso_na.events import default_events
from enso_na.budget import (fit_budget, select_reference_window, window_scan)
from enso_na.plotting import panel_label, save

MON = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]
SHOW = ["damp", "nino", "dnino", "pdo", "tna"]


def main():
    df = load_monthly()
    ev = default_events(df)
    scan = window_scan(df, ev, windows=WINDOWS)
    ref, tab = select_reference_window(df, ev)

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.2))
    fig.subplots_adjust(hspace=0.48, wspace=0.32)

    # ---- a -----------------------------------------------------------------
    ax = axes[0, 0]
    for t in SHOW + ["tend"]:
        col = "k" if t == "tend" else TERMS[t][1]
        lab = "$\\partial P/\\partial t$" if t == "tend" else TERMS[t][0]
        dv = [scan["dev"][k][t].mean() for k in WINDOWS]
        dc = [scan["dec"][k][t].mean() for k in WINDOWS]
        ax.plot(WINDOWS, dv, color=col, lw=1.2, marker="o", ms=3, label=lab)
        ax.plot(WINDOWS, dc, color=col, lw=1.2, marker="s", ms=3, ls="--")
    ax.axhline(0, color="0.4", lw=0.6)
    ax.axvline(ref, color="0.6", lw=0.8, ls=":")
    ax.set_xlabel("integration window $k$ (months)")
    ax.set_ylabel("phase integral ($\\sigma$)")
    ax.set_title("Phase integrals versus window\n(solid ASO0, dashed FMA1)",
                 fontsize=8)
    ax.legend(ncol=3, loc="lower right", fontsize=5.8, columnspacing=0.8)
    panel_label(ax, "a")

    # ---- b -----------------------------------------------------------------
    ax = axes[0, 1]
    adm = tab.admissible.values
    ax.bar(np.array(WINDOWS)[adm], tab.mismatch.values[adm], 0.6,
           color="#1f4e9c", lw=0, label="admissible")
    ax.bar(np.array(WINDOWS)[~adm], tab.mismatch.values[~adm], 0.6,
           color="0.75", lw=0, label="excluded")
    ax.axvline(ref, color="#c1272d", lw=1.0, ls="--")
    ax.set_xlabel("integration window $k$ (months)")
    ax.set_ylabel("$|$AsyI$_R| / |$AsyI$_{\\partial P/\\partial t}|$")
    ax.set_title("Reference-window selection", fontsize=8.5)
    ax.legend(loc="upper center")
    ax2 = ax.twinx()
    ax2.plot(WINDOWS, tab.p_tend.values, color="#e08214", marker="o", ms=3, lw=1.0)
    ax2.axhline(0.05, color="#e08214", lw=0.6, ls=":")
    ax2.set_ylabel("$p$ (AsyI tendency)", fontsize=7, color="#e08214")
    ax2.tick_params(labelsize=6, colors="#e08214")
    ax2.set_ylim(0, 1.3)
    panel_label(ax, "b")

    # ---- c -----------------------------------------------------------------
    ax = axes[1, 0]
    _, _, coefs, _ = fit_budget(df, ref)
    seas = coefs.groupby(coefs.index.month).mean()
    for t in SHOW:
        ax.plot(range(1, 13), seas[t].values, color=TERMS[t][1], lw=1.3,
                marker="o", ms=2.5, label=TERMS[t][0])
    ax.axhline(0, color="0.4", lw=0.6)
    for mm, c in (([8, 9, 10], "#4daf4a"), ([2, 3, 4], "#984ea3")):
        ax.axvspan(min(mm) - 0.5, max(mm) + 0.5, color=c, alpha=0.10, lw=0)
    ax.set_xticks(range(1, 13)); ax.set_xticklabels(MON, fontsize=6.5)
    ax.set_ylabel("coefficient (month$^{-1}$)")
    # VERIFY (display only): title shortened so it does not overlap the panel letter
    ax.set_title(f"Seasonal coefficients ($k$ = {ref} months)",
                 fontsize=8)
    ax.legend(ncol=2, loc="best", fontsize=6)
    panel_label(ax, "c")

    # ---- d -----------------------------------------------------------------
    ax = axes[1, 1]
    r2s = [scan["info"][k]["r2"] for k in WINDOWS]
    r2c = []
    for k in WINDOWS:
        info, _, _, _ = fit_budget(df, k, seasonal=False)
        r2c.append(info["r2"])
    ax.plot(WINDOWS, r2s, color="#1f4e9c", marker="o", ms=3.5, lw=1.3,
            label="seasonal coefficients")
    ax.plot(WINDOWS, r2c, color="#e08214", marker="s", ms=3.5, lw=1.3,
            label="constant coefficients")
    ax.axvline(ref, color="0.6", lw=0.8, ls=":")
    ax.set_xlabel("integration window $k$ (months)")
    ax.set_ylabel("explained variance of $\\partial P/\\partial t$")
    ax.set_title("Budget skill", fontsize=8.5)
    ax.legend(loc="upper right")
    panel_label(ax, "d")

    save(fig, "fig06")

    tab["r2_seasonal"] = r2s
    tab["r2_constant"] = r2c
    tab.round(3).to_csv(OUT_DIR / "fig06_window_selection.csv")
    print(tab.round(3).to_string())
    print(f"reference window = {ref}")


if __name__ == "__main__":
    main()
