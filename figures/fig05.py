"""
Figure 5. Seasonal gating of the two links of the ENSO-North America bridge.

a  Calendar-month regression of the convective index C on Nino3.4 (tropical link).
b  Calendar-month regression of the PNA index on C (extratropical link).
c  Calendar-month regression of the PNA index on Nino3.4 (net gain).
d  Product of the two link strengths, with the ASO0 and FMA1 means.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from enso_na.config import OUT_DIR, REF_WINDOW
from enso_na.data import load_monthly
from enso_na.events import default_events
from enso_na.budget import seasonal_coupling
from enso_na.plotting import panel_label, save

MON = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]
DEV_M = [8, 9, 10]     # Aug, Sep, Oct  (calendar months of ASO0)
DEC_M = [2, 3, 4]      # Feb, Mar, Apr  (calendar months of FMA1)


def _panel(ax, tab, title, ylab, col):
    m = tab.index.values
    sig = tab.p < 0.05
    ax.bar(m[~sig.values], tab.slope[~sig].values, 0.72, color=col, alpha=0.35,
           lw=0)
    ax.bar(m[sig.values], tab.slope[sig].values, 0.72, color=col, lw=0)
    ax.axhline(0, color="0.4", lw=0.6)
    ax.set_xticks(range(1, 13)); ax.set_xticklabels(MON, fontsize=6.5)
    ax.set_ylabel(ylab)
    ax.set_title(title, fontsize=8.5)
    ax2 = ax.twinx()
    ax2.plot(m, np.abs(tab.r.values), color="0.2", lw=1.0, marker="o", ms=2.5)
    ax2.set_ylim(0, 1.0)
    ax2.set_ylabel("|$r$|", fontsize=7, color="0.2")
    ax2.tick_params(labelsize=6, colors="0.2")
    for mm, c in ((DEV_M, "#4daf4a"), (DEC_M, "#984ea3")):
        ax.axvspan(min(mm) - 0.5, max(mm) + 0.5, color=c, alpha=0.10, lw=0)


def main():
    df = load_monthly()
    ev = default_events(df)
    k = REF_WINDOW
    t_nc = seasonal_coupling(df, k, "N", "C")
    t_cp = seasonal_coupling(df, k, "C", "P")
    t_np = seasonal_coupling(df, k, "N", "P")

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.8))
    fig.subplots_adjust(hspace=0.55, wspace=0.55)

    _panel(axes[0, 0], t_nc, "Tropical link: Nino3.4 $\\rightarrow$ OLR",
           "$\\partial C/\\partial N$", "#1b9e77")
    panel_label(axes[0, 0], "a")
    _panel(axes[0, 1], t_cp, "Extratropical link: OLR $\\rightarrow$ PNA",
           "$\\partial P/\\partial C$", "#7b3294")
    panel_label(axes[0, 1], "b")
    _panel(axes[1, 0], t_np, "Net gain: Nino3.4 $\\rightarrow$ PNA",
           "$\\partial P/\\partial N$", "#c1272d")
    panel_label(axes[1, 0], "c")

    ax = axes[1, 1]
    gain = t_nc.slope.values * t_cp.slope.values
    ax.bar(np.arange(1, 13), gain, 0.72, color="#1f4e9c", lw=0)
    ax.axhline(0, color="0.4", lw=0.6)
    g_dev = np.mean([gain[m - 1] for m in DEV_M])
    g_dec = np.mean([gain[m - 1] for m in DEC_M])
    # VERIFY (display only): key-period means in a legend (no overlap with axis/zero line)
    ax.axhline(g_dev, color="#4daf4a", lw=1.1, ls="--", label=f"ASO0 mean = {g_dev:+.2f}")
    ax.axhline(g_dec, color="#984ea3", lw=1.1, ls="--", label=f"FMA1 mean = {g_dec:+.2f}")
    ax.legend(loc="upper right", fontsize=6)
    yl = ax.get_ylim(); ax.set_ylim(yl[0], yl[1] + 0.25 * (yl[1] - yl[0]))
    ax.set_xticks(range(1, 13)); ax.set_xticklabels(MON, fontsize=6.5)
    ax.set_ylabel("$(\\partial C/\\partial N)(\\partial P/\\partial C)$")
    ax.set_title("Product of the two links", fontsize=8.5)
    panel_label(ax, "d")

    save(fig, "fig05")

    tab = pd.DataFrame({
        "slope_N_to_C": t_nc.slope, "r_N_to_C": t_nc.r, "p_N_to_C": t_nc.p,
        "slope_C_to_P": t_cp.slope, "r_C_to_P": t_cp.r, "p_C_to_P": t_cp.p,
        "slope_N_to_P": t_np.slope, "r_N_to_P": t_np.r, "p_N_to_P": t_np.p,
        "gain": gain})
    tab.round(3).to_csv(OUT_DIR / "fig05_seasonal_coupling.csv")
    print(tab.round(3).to_string())
    print(f"\nASO0 gain {g_dev:+.3f}   FMA1 gain {g_dec:+.3f}")
    print("tropical link ASO0 %.2f FMA1 %.2f" %
          (np.mean([t_nc.slope[m] for m in DEV_M]),
           np.mean([t_nc.slope[m] for m in DEC_M])))
    print("extratropical link ASO0 %.2f FMA1 %.2f" %
          (np.mean([t_cp.slope[m] for m in DEV_M]),
           np.mean([t_cp.slope[m] for m in DEC_M])))


if __name__ == "__main__":
    main()
