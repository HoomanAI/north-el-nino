"""Shared plotting style for all figures."""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .config import FIG_DIR, FIG_DPI
from .events import MONTH_LABELS

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.linewidth": 0.7,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6.5,
    "legend.frameon": False,
    "savefig.bbox": "tight",
    "figure.dpi": 110,
})

TICK_MONTHS = [0, 5, 11, 17, 23]
TICK_NAMES = [MONTH_LABELS[i] for i in TICK_MONTHS]


def lifecycle_axis(ax, mature=True):
    ax.set_xlim(-0.5, 23.5)
    ax.set_xticks(TICK_MONTHS)
    ax.set_xticklabels(TICK_NAMES)
    ax.axhline(0, color="0.4", lw=0.5)
    if mature:
        ax.axvline(11, color="0.3", lw=0.7, ls="--")


def shade_key_periods(ax, dev=(7, 9), dec=(13, 15), alpha=0.13):
    ax.axvspan(dev[0] - 0.5, dev[1] + 0.5, color="#2ca25f", alpha=alpha, lw=0)
    ax.axvspan(dec[0] - 0.5, dec[1] + 0.5, color="#7b3294", alpha=alpha, lw=0)


def panel_label(ax, text, dx=-0.10, dy=1.06):
    ax.text(dx, dy, text, transform=ax.transAxes, fontweight="bold",
            fontsize=9, va="bottom", ha="left")


def save(fig, name):
    """Write a figure as PNG and PDF.

    The PDF is vector and is what should be submitted to the journal; the PNG
    is for quick inspection. FIG_DPI only affects the PNG.
    """
    png = FIG_DIR / f"{name}.png"
    fig.savefig(png, dpi=FIG_DPI)
    fig.savefig(FIG_DIR / f"{name}.pdf")
    plt.close(fig)
    print(f"  wrote {png}")
    return png
