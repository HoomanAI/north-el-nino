"""
Figure 4. Internal damping, remote forcing and the origin of the asymmetry.

a  Composite evolution of the damping term, of the grouped tropical and
   extratropical forcing, and of the residual, with Nino3.4 for reference.
b  Key-period integrals of the same quantities.
c  Event-by-event phase-mean regression of the PNA index on Nino3.4 for ASO0
   and FMA1.
d  AsyI of the leading terms under the three phase-integration schemes.
e  Decomposition of the AsyI into amplitude-driven, seasonally modulated and
   residual components.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from enso_na.config import (KEY_DEC_MONTHS, KEY_DEV_MONTHS, OUT_DIR, REF_WINDOW,
                            TERMS)
from enso_na.data import load_monthly
from enso_na.events import ci95, composite, default_events
from enso_na.budget import (build_design, coupling_decomposition, fit_budget,
                            integrate_terms, phase_regression)
from enso_na.plotting import lifecycle_axis, panel_label, save, shade_key_periods

TROPICAL = ["nino", "olr", "dnino", "nl"]
EXTRATROPICAL = ["pdo", "tna", "ao", "nao"]
GROUPS = [("$-\\lambda P$", ["damp"], "#c1272d"),
          ("tropical forcing", TROPICAL, "#1f4e9c"),
          ("extratropical forcing", EXTRATROPICAL, "#e08214"),
          ("$R$", ["res"], "0.55")]


def main():
    df = load_monthly()
    ev = default_events(df)
    info, terms, coefs, X = fit_budget(df, REF_WINDOW)
    grouped = pd.DataFrame({name: terms[keys].sum(axis=1)
                            for name, keys, _ in GROUPS})
    grouped["$\\partial P/\\partial t$"] = terms["tend"]

    fig = plt.figure(figsize=(7.2, 7.0))
    gs = fig.add_gridspec(3, 2, hspace=0.55, wspace=0.32,
                          height_ratios=[1.0, 1.0, 1.0])
    x = np.arange(24)

    # ---- a ------------------------------------------------------------------
    ax = fig.add_subplot(gs[0, :])
    shade_key_periods(ax)
    for name, keys, col in GROUPS:
        m, ci, _ = composite(terms[keys].sum(axis=1), ev)
        ax.plot(x, m, color=col, lw=1.3, label=name)
        ax.fill_between(x, m - ci, m + ci, color=col, alpha=0.12, lw=0)
    m, _, _ = composite(terms["tend"], ev)
    ax.plot(x, m, color="k", lw=1.0, ls=":", label="$\\partial P/\\partial t$")
    lifecycle_axis(ax)
    ax.set_ylabel("$\\sigma$ month$^{-1}$")
    ax2 = ax.twinx()
    mn, _, _ = composite(X["N"], ev)
    ax2.plot(x, mn, color="0.35", lw=1.0, ls="--")
    ax2.set_ylabel("Nino3.4 ($\\sigma$)", fontsize=7, color="0.35")
    ax2.tick_params(labelsize=6.5, colors="0.35")
    ax.legend(loc="lower left", ncol=3)
    ax.set_title("Grouped budget contributions", fontsize=8.5)
    panel_label(ax, "a", dx=-0.06)

    # ---- b ------------------------------------------------------------------
    ax = fig.add_subplot(gs[1, 0])
    names = [g[0] for g in GROUPS] + ["$\\partial P/\\partial t$"]
    dev, dec = integrate_terms(grouped.assign(**{
        k: terms[k] for k in ["res"]}), ev, "key") if False else (None, None)
    devs, decs, cid, cic = [], [], [], []
    for name in names:
        M_dev, M_dec = [], []
        s = grouped[name]
        from enso_na.events import lifecycle_matrix
        M = lifecycle_matrix(s, ev)
        d = np.nansum(M[:, list(KEY_DEV_MONTHS)], axis=1)
        c = np.nansum(M[:, list(KEY_DEC_MONTHS)], axis=1)
        devs.append(d.mean()); decs.append(c.mean())
        cid.append(ci95(d)); cic.append(ci95(c))
    idx = np.arange(len(names))
    ax.bar(idx - 0.2, devs, 0.38, color="#4daf4a", lw=0, label="ASO0")
    ax.bar(idx + 0.2, decs, 0.38, color="#984ea3", lw=0, label="FMA1")
    ax.errorbar(idx - 0.2, devs, yerr=cid, fmt="none", ecolor="0.25",
                elinewidth=0.5, capsize=1.2)
    ax.errorbar(idx + 0.2, decs, yerr=cic, fmt="none", ecolor="0.25",
                elinewidth=0.5, capsize=1.2)
    ax.axhline(0, color="0.4", lw=0.6)
    ax.set_xticks(idx); ax.set_xticklabels(names, rotation=25, ha="right", fontsize=6)
    ax.set_ylabel("phase integral ($\\sigma$)")
    ax.legend(loc="lower left")
    ax.set_title("Key-period integrals", fontsize=8.5)
    panel_label(ax, "b")

    # ---- c ------------------------------------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    pr = phase_regression(X["N"], X["P"], ev)
    for name, col, lab in (("dev", "#4daf4a", "ASO0"), ("dec", "#984ea3", "FMA1")):
        d = pr[name]
        ax.scatter(d["x"], d["y"], color=col, s=22, edgecolor="0.3", lw=0.4,
                   zorder=3)
        xs = np.linspace(min(d["x"]), max(d["x"]), 10)
        ax.plot(xs, d["slope"] * xs + d["intercept"], color=col, lw=1.2,
                label=f"{lab}: slope={d['slope']:.2f}, $r$={d['r']:.2f}")
    ax.axhline(0, color="0.6", lw=0.5)
    ax.set_xlabel("Nino3.4 ($\\sigma$)"); ax.set_ylabel("PNA ($\\sigma$)")
    ax.legend(loc="upper left")
    ax.set_title("Phase-mean ENSO-PNA relation", fontsize=8.5)
    panel_label(ax, "c")

    # ---- d ------------------------------------------------------------------
    ax = fig.add_subplot(gs[2, 0])
    schemes = [("key", "ASO0 / FMA1", "#1f4e9c"),
               ("phase9", "9-month phases", "#e08214"),
               ("aligned", "9/7-month phases", "#4daf4a")]
    show = ["damp", "nino", "dnino", "pdo", "tna", "nao", "res"]
    idx = np.arange(len(show))
    for j, (sch, lab, col) in enumerate(schemes):
        dv, dc = integrate_terms(terms, ev, sch)
        a = (dv + dc).mean()
        e = [ci95((dv[t] + dc[t]).values) for t in show]
        ax.bar(idx + (j - 1) * 0.27, [a[t] for t in show], 0.25, color=col,
               lw=0, label=lab)
        ax.errorbar(idx + (j - 1) * 0.27, [a[t] for t in show], yerr=e,
                    fmt="none", ecolor="0.3", elinewidth=0.4, capsize=1.0)
    ax.axhline(0, color="0.4", lw=0.6)
    ax.set_xticks(idx)
    ax.set_xticklabels([TERMS[t][0] for t in show], rotation=90, fontsize=6)
    ax.set_ylabel("AsyI ($\\sigma$)")
    ax.legend(loc="upper center", ncol=3, fontsize=5.5, handlelength=1.0, columnspacing=0.8)   # VERIFY: legend off bars
    yl = ax.get_ylim(); ax.set_ylim(yl[0], yl[1] + 0.22 * (yl[1] - yl[0]))
    ax.set_title("Sensitivity to the integration scheme", fontsize=8.5)
    panel_label(ax, "d")

    # ---- e ------------------------------------------------------------------
    ax = fig.add_subplot(gs[2, 1])
    dec_tab = coupling_decomposition(df, ev, REF_WINDOW).astype(float)
    labels = ["amplitude\n(I)", "seasonal\n(II)", "residual\n(III)", "total\nbudget"]
    vals = [dec_tab.amp.sum(), dec_tab.seasonal.sum(), dec_tab.residual.sum(),
            dec_tab.asyi_F.sum()]
    cols = ["#c1272d", "#1f4e9c", "#999999", "#4d4d4d"]
    ax.bar(np.arange(4), vals, 0.6, color=cols, lw=0)
    lo, hi = min(min(vals), 0), max(max(vals), 0); span = hi - lo    # VERIFY: labels inside the axes
    for i, v in enumerate(vals):
        ax.text(i, v + (0.03 * span if v >= 0 else -0.03 * span), f"{v:+.2f}", ha="center",
                fontsize=6.5, va="bottom" if v >= 0 else "top")
    ax.set_ylim(lo - 0.18 * span, hi + 0.18 * span)
    ax.axhline(0, color="0.4", lw=0.6)
    ax.set_xticks(np.arange(4)); ax.set_xticklabels(labels, fontsize=6.5)
    ax.set_ylabel("AsyI ($\\sigma$)")
    ax.set_title("Amplitude versus seasonal modulation", fontsize=8.5)
    panel_label(ax, "e")

    save(fig, "fig04")

    dec_tab.round(3).to_csv(OUT_DIR / "fig04_decomposition.csv")
    out = dict(
        lambda_annual=info["lambda_annual"], r2=info["r2"],
        slope_dev=pr["dev"]["slope"], r_dev=pr["dev"]["r"], p_dev=pr["dev"]["p"],
        slope_dec=pr["dec"]["slope"], r_dec=pr["dec"]["r"], p_dec=pr["dec"]["p"],
        amp_sum=dec_tab.amp.sum(), seasonal_sum=dec_tab.seasonal.sum(),
        resid_sum=dec_tab.residual.sum(), total=dec_tab.asyi_F.sum())
    pd.Series(out).round(4).to_csv(OUT_DIR / "fig04_numbers.csv")
    print(pd.Series(out).round(3))
    print(dec_tab[["a_dev", "a_dec", "amp", "seasonal", "residual"]].round(3))


if __name__ == "__main__":
    main()
