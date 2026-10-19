#!/usr/bin/env python3
"""
Regenerate every figure and every table of the paper.

Usage
-----
    python run_all.py                  # all figures and tables
    python run_all.py fig01 fig05      # selected items only
    python verify_numbers.py           # assert the published values

Environment
-----------
    ENSO_NA_DATA    directory holding the index NetCDF files (default: ./data)
    ENSO_NA_FIGDIR  figure output directory                  (default: ./figures)
    ENSO_NA_OUTDIR  table output directory                   (default: ./outputs)
"""
from __future__ import annotations

import importlib
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "figures"))
sys.path.insert(0, str(ROOT / "tables"))

ALL = ([f"fig{i:02d}" for i in range(1, 8)]
       + [f"figS{i:02d}" for i in range(1, 6)]
       + ["table07_interaction"])


def main(names):
    names = names or ALL
    t0 = time.time()
    for name in names:
        print(f"\n=== {name} " + "=" * (60 - len(name)))
        mod = importlib.import_module(name)
        mod.main()
    print(f"\nDone in {time.time() - t0:.1f} s. "
          f"Figures in {ROOT / 'figures'}, tables in {ROOT / 'outputs'}.")


if __name__ == "__main__":
    main(sys.argv[1:])
