"""
Figure 3. Asymmetric contributions of the budget terms.

a  Stacked asymmetry index (AsyI) of all terms for every integration window,
   with the summed AsyI and the AsyI of the tendency for reference.
b  AsyI of each individual term as a function of the integration window.
c  AsyI versus asymmetric attribution efficiency (AAE) at the reference window;
   open symbols show the five tropical-index configurations of the ensemble.
d  Window trend index (WTI) of every term, development versus decay phase.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from enso_na.config import OUT_DIR, REF_WINDOW, TERMS, TERM_ORDER, WINDOWS
from enso_na.data import load_monthly
from enso_na.events import ci95, default_events
from enso_na.budget import aae, asyi, window_scan, wti
from enso_na.plotting import panel_label, save

#: ensemble configurations: (tropical index column, invert sign?)
CONFIGS = [("NINO34", False), ("NINO3", False), ("NINO4", False),
           ("NINO12", False), ("SOI", True)]
MARKERS = ["o", "s", "v", "D", "^"]
BUDGET_TERMS = [t for t in TERM_ORDER if t != "res"]


def main():
    df = load_monthly()
    ev = default_events(df)
    scan = window_scan(df, ev, windows=WINDOWS)
    ens = {name: window_scan(df, ev, nino_col=name, flip_nino=flip,
                             windows=[REF_WINDOW])
           for name, flip in CONFIGS}

    fig = plt.figure(figsize=(7.2, 7.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 1.0],
                          hspace=0.62, wspace=0.30)

    # ---- a: stacked AsyI across windows ------------------------------------
    ax = fig.add_subplot(gs[0, :])
    w = np.arange(len(WINDOWS))
    A = scan["asyi_df"]
    bot_p = np.zeros(len(WINDOWS)); bot_n = np.zeros(len(WINDOWS))
    for t in BUDGET_TERMS:
        v = A.loc[t, WINDOWS].values.astype(float)
        pos = np.where(v > 0, v, 0); neg = np.where(v < 0, v, 0)
        ax.bar(w - 0.17, pos, 0.3, bottom=bot_p, color=TERMS[t][1], lw=0)
        ax.bar(w - 0.17, neg, 0.3, bottom=bot_n, color=TERMS[t][1], lw=0)
        bot_p += pos; bot_n += neg
    ssum = A.loc[BUDGET_TERMS, WINDOWS].sum().values
    tend = A.loc["tend", WINDOWS].values.astype(float)
    ax.bar(w + 0.17, ssum, 0.15, color="#1f4e9c", lw=0, label="$\\Sigma$AsyI")
    ax.bar(w + 0.33, tend, 0.15, color="0.6", lw=0, label="AsyI[$\\partial P/\\partial t$]")
    # VERIFY (display only): data-driven limits so stacks and labels are never clipped
    top = max(bot_p.max(), ssum.max(), tend.max(), 0); bot = min(bot_n.min(), ssum.min(), tend.min(), 0)
    span = top - bot; ytop = top + 0.06 * span
    for i in range(len(WINDOWS)):
        ax.text(w[i], ytop, f"$\\Delta$={ssum[i]-tend[i]:+.2f}", fontsize=5.5,
                color="#c1272d", ha="center", va="bottom")
    ax.set_ylim(bot - 0.08 * span, top + 0.22 * span)
    ax.axhline(0, color="0.4", lw=0.6)
    ax.set_xticks(w); ax.set_xticklabels([f"{k} m" for k in WINDOWS])
    ax.set_ylabel("AsyI ($\\sigma$)")
    ax.set_title("Asymmetry index summation versus tendency asymmetry", fontsize=8.5)
    ax.legend(loc="lower right", ncol=2)   # VERIFY (display only): legend off the bars
    panel_label(ax, "a", dx=-0.06)

    # ---- b: AsyI of each term versus window --------------------------------
    ax = fig.add_subplot(gs[1, :])
    for t in BUDGET_TERMS:
        y = A.loc[t, WINDOWS].values.astype(float)
        e = [ci95((scan["dev"][k][t] + scan["dec"][k][t]).values) for k in WINDOWS]
        ax.errorbar(WINDOWS, y, yerr=e, color=TERMS[t][1], marker="o", ms=3,
                    lw=1.2, elinewidth=0.5, capsize=1.5, label=TERMS[t][0])
    ax.axhline(0, color="0.4", lw=0.6)
    ax.set_xlabel("integration window $k$ (months)")
    ax.set_ylabel("AsyI ($\\sigma$)")
    ax.set_title("Window dependence of the term AsyI", fontsize=8.5)
    lo = min(float(np.nanmin(A.loc[BUDGET_TERMS, WINDOWS].values.astype(float))) - 0.6, -1.0)
    ax.set_ylim(lo, max(1.0, float(np.nanmax(A.loc[BUDGET_TERMS, WINDOWS].values.astype(float))) + 0.4))
    ax.legend(ncol=5, loc="lower center", fontsize=5.8, columnspacing=1.0,
              handletextpad=0.4)
    panel_label(ax, "b", dx=-0.06)

    # ---- c: AsyI versus AAE ------------------------------------------------
    ax = fig.add_subplot(gs[2, 0])
    for j, (name, _) in enumerate(CONFIGS):
        s = ens[name]
        av = asyi(s["dev"][REF_WINDOW], s["dec"][REF_WINDOW])
        ae = aae(s["dev"][REF_WINDOW], s["dec"][REF_WINDOW])
        for t in BUDGET_TERMS:
            ax.scatter(ae[t], av[t], marker=MARKERS[j], s=26 if j == 0 else 16,
                       facecolor=TERMS[t][1] if j == 0 else "none",
                       edgecolor=TERMS[t][1], linewidth=0.7, zorder=3,
                       alpha=1.0 if j == 0 else 0.75)
    ax.axhline(0, color="0.6", lw=0.6, ls="--")
    ax.axvline(0, color="0.6", lw=0.6, ls="--")
    ax.set_xlim(-1.15, 1.15)
    ax.set_xlabel("AAE"); ax.set_ylabel("AsyI ($\\sigma$)")
    ax.set_title(f"AsyI versus AAE, $k$ = {REF_WINDOW} months", fontsize=8.5)
    # VERIFY (display only): legend identifying the ENSO-index configurations
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], marker=MARKERS[j], ls="none", color="0.3",
                              markerfacecolor="0.3" if j == 0 else "none", ms=4,
                              label=("SOI (inv.)" if n == "SOI" else n.replace("NINO", "Nino")))
                       for j, (n, _) in enumerate(CONFIGS)],
              loc="upper left", fontsize=5.5, ncol=1, handletextpad=0.2, borderaxespad=0.2)
    panel_label(ax, "c")

    # ---- d: window trend index ---------------------------------------------
    ax = fig.add_subplot(gs[2, 1])
    wd, wc = wti(scan, "dev"), wti(scan, "dec")
    idx = np.arange(len(BUDGET_TERMS))
    ax.bar(idx - 0.2, [wd[t] for t in BUDGET_TERMS], 0.38,
           color=[TERMS[t][1] for t in BUDGET_TERMS], lw=0)
    ax.bar(idx + 0.2, [wc[t] for t in BUDGET_TERMS], 0.38,
           color="none", edgecolor=[TERMS[t][1] for t in BUDGET_TERMS],
           hatch="////", lw=0.8)
    ax.axhline(0, color="0.4", lw=0.6)
    ax.set_xticks(idx)
    ax.set_xticklabels([TERMS[t][0] for t in BUDGET_TERMS], rotation=90, fontsize=6)
    ax.set_ylabel("WTI ($\\sigma$)")
    ax.set_title("Window trend index", fontsize=8.5)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor="0.45", label="development"),
                       Patch(facecolor="none", edgecolor="0.45", hatch="////",
                             label="decay")], loc="upper right")   # VERIFY (display only): legend off the damping bar
    panel_label(ax, "d")

    save(fig, "fig03")

    tab = pd.DataFrame({
        "AsyI": A[REF_WINDOW], "AAE": scan["aae_df"][REF_WINDOW],
        "WTI_dev": wd, "WTI_dec": wc,
        "I_dev": scan["dev"][REF_WINDOW].mean(),
        "I_dec": scan["dec"][REF_WINDOW].mean(),
        "I_dev_ci": [ci95(scan["dev"][REF_WINDOW][t].values) for t in A.index],
        "I_dec_ci": [ci95(scan["dec"][REF_WINDOW][t].values) for t in A.index],
    })
    for name, _ in CONFIGS:
        s = ens[name]
        tab[f"AsyI_{name}"] = asyi(s["dev"][REF_WINDOW], s["dec"][REF_WINDOW])
    tab.round(3).to_csv(OUT_DIR / "fig03_asyi_aae_wti.csv")
    print(tab.round(3).to_string())


if __name__ == "__main__":
    main()
