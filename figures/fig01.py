"""
Figure 1. Asymmetric development and decay of El Nino and of its North American
teleconnection.

a  Development versus decay speed of the Nino3.4 index for every event.
b  Composite life cycle of Nino3.4 and of the PNA index.
c  Composite PNA tendency, its point-symmetric counterpart, and the months at
   which development and decay speeds differ significantly.
"""
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

from enso_na.config import KEY_DEC_MONTHS, KEY_DEV_MONTHS, OUT_DIR, REF_WINDOW
from enso_na.data import load_monthly, running_mean, standardise
from enso_na.events import (MONTH_LABELS, composite, default_events,
                            lifecycle_matrix, phase_speeds, point_symmetric,
                            significant_months)
from enso_na.plotting import (lifecycle_axis, panel_label, save,
                              shade_key_periods)


def main():
    df = load_monthly()
    ev = default_events(df)
    X = running_mean(standardise(df[["PNA", "NINO34"]]), REF_WINDOW)

    sp = phase_speeds(X["NINO34"], ev)
    t = stats.ttest_rel(sp.dev_speed, sp.dec_speed)

    fig = plt.figure(figsize=(7.2, 5.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.45, wspace=0.42)

    # ---- a: development vs decay speed of Nino3.4 --------------------------
    ax = fig.add_subplot(gs[0, 0])
    sc = ax.scatter(sp.dev_speed, sp.dec_speed, c=sp.peak, cmap="OrRd",
                    vmin=1.0, vmax=3.2, s=42, edgecolor="0.25", linewidth=0.5,
                    zorder=3)
    # VERIFY (display only): label placement that avoids other labels AND data markers
    sc_pt = 268.0                                   # approx. points per data unit in this panel
    pts = [(r.dev_speed * sc_pt, r.dec_speed * sc_pt) for _, r in sp.iterrows()]
    ang = np.deg2rad(np.arange(0, 360, 30)); placed = []
    for y, r in sp.sort_values("dec_speed", ascending=False).iterrows():
        px, py = r.dev_speed * sc_pt, r.dec_speed * sc_pt
        best = None
        for rad in (9, 14, 20, 27):
            for a in ang:
                dx, dy = rad * np.cos(a), rad * np.sin(a)
                cx = px + dx + (8 if np.cos(a) >= 0 else -8); cy = py + dy   # label centre (w~16pt, h~6pt)
                ok_l = all(abs(cx - qx) > 17 or abs(cy - qy) > 7 for qx, qy in placed)
                ok_p = all(abs(cx - qx) > 12 or abs(cy - qy) > 7 for qx, qy in pts)
                if ok_l and ok_p:
                    best = (dx, dy, cx, cy); break
            if best: break
        dx, dy, cx, cy = best if best else (4, 4, px + 12, py + 4)
        placed.append((cx, cy))
        ax.annotate(str(y), (r.dev_speed, r.dec_speed), fontsize=5.5, xytext=(dx, dy),
                    textcoords="offset points", color="0.25", va="center",
                    ha="left" if dx >= 0 else "right",
                    arrowprops=dict(arrowstyle="-", lw=0.3, color="0.6"))
    lim = [0, 0.62]
    ax.plot(lim, lim, ls="--", color="0.4", lw=0.8)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("Development speed ($\\sigma$ month$^{-1}$)")
    ax.set_ylabel("Decay speed ($\\sigma$ month$^{-1}$)")
    ax.set_title("Nino3.4 phase speeds", fontsize=8.5)
    ax.text(0.96, 0.06, f"$p$ = {t.pvalue:.3f}", transform=ax.transAxes,
            ha="right", fontsize=7)
    cb = fig.colorbar(sc, ax=ax, pad=0.03, fraction=0.04)
    cb.set_label("peak amplitude ($\\sigma$)", fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    panel_label(ax, "a")

    # ---- b: composite life cycles ------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    m_n, ci_n, _ = composite(X["NINO34"], ev)
    m_p, ci_p, Mp = composite(X["PNA"], ev)
    x = np.arange(24)
    shade_key_periods(ax)
    ax.plot(x, m_n, color="k", lw=1.5, label="Nino3.4")
    ax.fill_between(x, m_n - ci_n, m_n + ci_n, color="k", alpha=0.12, lw=0)
    ax.plot(x, m_p, color="#c1272d", lw=1.5, label="PNA")
    ax.fill_between(x, m_p - ci_p, m_p + ci_p, color="#c1272d", alpha=0.15, lw=0)
    lifecycle_axis(ax)
    ax.set_ylabel("standardised anomaly ($\\sigma$)")
    ax.set_title("Composite life cycle", fontsize=8.5)
    ax.legend(loc="upper left")
    panel_label(ax, "b")

    # ---- c: PNA tendency and its point-symmetric counterpart ---------------
    ax = fig.add_subplot(gs[1, :])
    dP = X["PNA"].shift(-1) - X["PNA"]
    m_t, ci_t, Mt = composite(dP, ev)
    ps = point_symmetric(m_t, pivot=11)
    p, sig = significant_months(Mt, pivot=11)
    for k in sig:
        ax.axvspan(k - 0.4, k + 0.4, color="0.75", alpha=0.55, lw=0)
    ax.plot(x, m_t, color="#c1272d", lw=1.6, label="PNA tendency")
    ax.fill_between(x, m_t - ci_t, m_t + ci_t, color="#c1272d", alpha=0.15, lw=0)
    ax.plot(x, ps, color="#1b9e77", lw=1.4, ls="--",
            label="point-symmetric counterpart")
    lifecycle_axis(ax)
    ax.set_xticks(range(24))
    ax.set_xticklabels(MONTH_LABELS, rotation=90, fontsize=5.5)
    ax.set_ylabel("$\\partial P/\\partial t$ ($\\sigma$ month$^{-1}$)")
    ax.legend(loc="upper left", ncol=2)
    panel_label(ax, "c", dx=-0.05)

    save(fig, "fig01")

    # ---- numbers used in the text ------------------------------------------
    P = lifecycle_matrix(X["PNA"], ev)
    dev = np.nanmean(P[:, list(KEY_DEV_MONTHS)], 1)
    dec = np.nanmean(P[:, list(KEY_DEC_MONTHS)], 1)
    out = dict(
        n_events=len(ev), events=ev,
        nino_dev_speed=sp.dev_speed.mean(), nino_dec_speed=sp.dec_speed.mean(),
        nino_speed_p=t.pvalue,
        nino_dev_len=sp.dev_len.mean(), nino_dec_len=sp.dec_len.mean(),
        pna_peak_month=MONTH_LABELS[int(np.nanargmax(m_p))],
        pna_peak=float(np.nanmax(m_p)),
        pna_ASO=float(dev.mean()), pna_FMA=float(dec.mean()),
        nino_ASO=float(np.nanmean(lifecycle_matrix(X["NINO34"], ev)[:, list(KEY_DEV_MONTHS)])),
        nino_FMA=float(np.nanmean(lifecycle_matrix(X["NINO34"], ev)[:, list(KEY_DEC_MONTHS)])),
        sig_months=[MONTH_LABELS[i] for i in sig],
    )
    pd.Series(out).to_csv(OUT_DIR / "fig01_numbers.csv")
    print(pd.Series(out))


if __name__ == "__main__":
    main()
