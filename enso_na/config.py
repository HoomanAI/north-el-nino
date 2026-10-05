"""
Configuration for the ENSO-North America phase-asymmetry analysis.

All paths can be overridden with environment variables so that the package runs
unchanged in any environment:

    ENSO_NA_DATA    directory holding the index NetCDF files (default: ../data)
    ENSO_NA_FIGDIR  directory for figure output          (default: ../figures)
    ENSO_NA_OUTDIR  directory for tables / numbers       (default: ../outputs)
"""
from __future__ import annotations

import os
from pathlib import Path

PKG_DIR = Path(__file__).resolve().parent
ROOT = PKG_DIR.parent

DATA_DIR = Path(os.environ.get("ENSO_NA_DATA", ROOT / "data"))
FIG_DIR = Path(os.environ.get("ENSO_NA_FIGDIR", ROOT / "figures"))
OUT_DIR = Path(os.environ.get("ENSO_NA_OUTDIR", ROOT / "outputs"))
for _d in (FIG_DIR, OUT_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ----------------------------------------------------------------------------
# Analysis constants
# ----------------------------------------------------------------------------

#: Period common to all monthly indices used in the tendency budget.
ANALYSIS_START = "1982-01"
ANALYSIS_END = "2026-07"

#: El Nino identification: DJF Nino3.4 (3-month running mean) threshold (degC).
ENSO_THRESHOLD = 0.5

#: Temporal integration (running-mean) windows, in months. This is the free
#: methodological parameter of the framework, analogous to the mixed-layer
#: integration depth of Song et al. (2026).
WINDOWS = [1, 2, 3, 4, 5, 6, 7]

#: Reference window used for the main diagnosis (selected objectively as the
#: window minimising the mean absolute budget residual; see budget.py).
REF_WINDOW = 2

#: Key development / decay periods (month indices relative to Jan of the
#: development year, 0-based: 0 = Jan0, 12 = Jan1).
#: Filled by events.significant_months(); the defaults below are the values
#: diagnosed from the observational record and are used for reproducibility.
KEY_DEV_MONTHS = (7, 8, 9)      # Aug0-Oct0
KEY_DEC_MONTHS = (13, 14, 15)   # Feb1-Apr1

#: Budget term keys, display labels (matplotlib mathtext) and plot colours.
TERMS = {
    "damp": (r"$-\lambda P$", "#c1272d"),
    "nino": (r"$a_N N$", "#0b6fa4"),
    "olr": (r"$a_C C$", "#1b9e77"),
    "dnino": (r"$a_{\dot N}\,\partial N/\partial t$", "#7b4ea3"),
    "pdo": (r"$a_D D$", "#e08214"),
    "tna": (r"$a_T T$", "#66a61e"),
    "ao": (r"$a_A A$", "#8c6d31"),
    "nao": (r"$a_M M$", "#3b9ab2"),
    "nl": (r"NL", "#d95f02"),
    "res": (r"$R$", "#808080"),
}
TERM_ORDER = ["damp", "nino", "olr", "dnino", "pdo", "tna", "ao", "nao", "nl", "res"]

#: Plain-text labels (for CSV export and for figures without mathtext).
TERM_TEXT = {
    "damp": "internal damping",
    "nino": "tropical Pacific SST forcing",
    "olr": "tropical convective forcing",
    "dnino": "ENSO tendency forcing",
    "pdo": "North Pacific background",
    "tna": "tropical North Atlantic forcing",
    "ao": "Arctic Oscillation",
    "nao": "North Atlantic Oscillation",
    "nl": "nonlinear ENSO forcing",
    "res": "residual",
}

#: Figure style.
FIG_DPI = 300
