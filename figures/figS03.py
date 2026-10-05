"""
Figure S3. Overview of the index data used in the study.

a  Monthly time series of the six budget predictors and of the PNA index.
b  Record coverage of every supplied index file.
c  Cross-correlation matrix of the monthly indices over the common period.
d  Autocorrelation of the PNA index for the integration windows considered.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from enso_na.config import ANALYSIS_END, ANALYSIS_START, OUT_DIR, WINDOWS
from enso_na.data import FILES, _read_series, _to_monthly, load_monthly, running_mean, standardise
from enso_na.plotting import panel_label, save

SHOW = ["PNA", "NINO34", "OLR", "PDO", "TNA", "AO", "NAO"]


def main():
    df = load_monthly()
    Z = standardise(df[SHOW])

    fig = plt.figure(figsize=(7.2, 7.6))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.5, 1.0, 1.0], hspace=0.45,
                          wspace=0.32)

    # ---- a ------------------------------------------------------------------
    ax = fig.add_subplot(gs[0, :])
    for i, c in enumerate(SHOW):
        ax.plot(Z.index, Z[c] - 4 * i, lw=0.55, color=plt.cm.tab10(i % 10))
        ax.text(Z.index[0] - pd.Timedelta(days=150), -4 * i, c, fontsize=6.5, ha="right",
                va="center", color=plt.cm.tab10(i % 10))   # VERIFY (display only): labels outside the curves
    ax.set_yticks([])
    ax.set_xlim(Z.index[0] - pd.Timedelta(days=1500), Z.index[-1] + pd.Timedelta(days=120))
    ax.set_title(f"Standardised monthly indices, {ANALYSIS_START} to {ANALYSIS_END}",
                 fontsize=8.5)
    panel_label(ax, "a", dx=-0.05)

    # ---- b ------------------------------------------------------------------
    ax = fig.add_subplot(gs[1, :])
    rows = []
    for j, (key, (fname, var, samp)) in enumerate(FILES.items()):
        s = _to_monthly(_read_series(fname, var), samp).dropna()
        rows.append((key, s.index[0], s.index[-1], samp, len(s)))
        ax.plot([s.index[0], s.index[-1]], [j, j], lw=4,
                color={"monthly": "#1f4e9c", "weekly": "#1b9e77",
                       "daily": "#e08214"}[samp], solid_capstyle="butt")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=6.5)
    ax.set_xlim(pd.Timestamp("1945-01-01"), pd.Timestamp("2030-01-01"))
    ax.axvspan(pd.Timestamp(ANALYSIS_START), pd.Timestamp(ANALYSIS_END),
               color="0.8", alpha=0.4, lw=0)
    ax.set_title("Record coverage (shading: analysis period; PDO extends to 1854)",
                 fontsize=8.5)
    panel_label(ax, "b", dx=-0.05)

    # ---- c ------------------------------------------------------------------
    ax = fig.add_subplot(gs[2, 0])
    C = df[SHOW].corr()
    im = ax.imshow(C.values, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(SHOW))); ax.set_xticklabels(SHOW, rotation=90, fontsize=6)
    ax.set_yticks(range(len(SHOW))); ax.set_yticklabels(SHOW, fontsize=6)
    for i in range(len(SHOW)):
        for j in range(len(SHOW)):
            ax.text(j, i, f"{C.values[i, j]:.2f}", ha="center", va="center",
                    fontsize=4.8,
                    color="w" if abs(C.values[i, j]) > 0.55 else "0.1")
    fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03).ax.tick_params(labelsize=6)
    ax.set_title("Monthly cross-correlations", fontsize=8.5)
    panel_label(ax, "c", dx=-0.42)   # VERIFY (display only): no overlap with title

    # ---- d ------------------------------------------------------------------
    ax = fig.add_subplot(gs[2, 1])
    lags = np.arange(0, 13)
    for k in WINDOWS:
        p = running_mean(standardise(df[["PNA"]]), k)["PNA"].dropna()
        ac = [p.autocorr(l) for l in lags]
        ax.plot(lags, ac, lw=1.1, marker="o", ms=2.5,
                color=plt.cm.viridis((k - 1) / 6), label=f"$k$={k}")
    ax.axhline(0, color="0.4", lw=0.6)
    ax.set_xlabel("lag (months)"); ax.set_ylabel("autocorrelation")
    ax.set_title("PNA autocorrelation by window", fontsize=8.5)
    ax.legend(ncol=2, fontsize=6)
    panel_label(ax, "d")

    save(fig, "figS03")
    pd.DataFrame(rows, columns=["index", "start", "end", "sampling", "n_months"]
                 ).to_csv(OUT_DIR / "figS03_coverage.csv", index=False)
    C.round(2).to_csv(OUT_DIR / "figS03_correlations.csv")
    print(pd.DataFrame(rows, columns=["index", "start", "end", "sampling",
                                      "n_months"]).to_string(index=False))


if __name__ == "__main__":
    main()
