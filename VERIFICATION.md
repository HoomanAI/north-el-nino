# Verification notes

**Scope: verification only.** The question was whether the manuscript text,
captions and figures match the outputs of the authors' own code and data. No new
analysis was done and no method was changed.

## Result

| Outcome | Count |
|---|---|
| Quantitative statements checked | 58 |
| Match the code output | 32 |
| Did not match (corrected in the manuscript with tracked changes) | 24 |
| Cannot be checked from this code and data (Song et al. 2026, ref. 15) | 2 (still unverified) |
| Output tables identical to the original package | 15 of 15 |

The full list is in [`verification/verification_checks.csv`](verification/verification_checks.csv).
Reproduce it with:

```bash
python run_all.py
python verification/verify_claims.py
```

The "paper statement" column quotes the **original** manuscript. The 24 MISMATCH
rows are the statements corrected in the tracked manuscript
(`ENSO_North_America_teleconnection_paper_VERIFIED_tracked.docx`, not part of this repository).

## Code changes (display only)

Each changed line is marked `# VERIFY (display only)`. None of them changes a computed value.

| Figure | File | Change |
|---|---|---|
| All | `enso_na/plotting.py` | Output names get the suffix `_verified` (set `ENSO_NA_SUFFIX=` to disable); a one-line footnote says the figure is a verification copy |
| Fig. 1a | `figures/fig01.py` | Year labels placed so they no longer overlap each other or the markers |
| Fig. 2 | `figures/fig02.py` | Axis limits follow the data (the k = 1 bars were clipped) |
| Fig. 3a | `figures/fig03.py` | Axis limits follow the data; legend moved off the bars |
| Fig. 3c | `figures/fig03.py` | Legend added that labels the five existing marker types (no new data) |
| Fig. 3d | `figures/fig03.py` | Legend moved to the upper right, off the damping bar. Bars unchanged: the tendency and residual are **not** added |
| Fig. 4d | `figures/fig04.py` | Legend moved above the bars, with y-axis headroom |
| Fig. 4e | `figures/fig04.py` | Value labels kept inside the axes |
| Fig. 5d | `figures/fig05.py` | Key-period means shown in a legend instead of over the axis |
| Fig. 6c | `figures/fig06.py` | Panel title shortened so it does not overlap the panel letter |
| Fig. S3 | `figures/figS03.py` | Series labels and panel letter moved off the curves and title |
| Fig. S5a | `figures/figS05.py` | Legend moved above the bars |

Figs. 7, S1, S2 and S4 needed no change.

## README changes

Two statements in the original README did not match the data files and were aligned:

- The data have four sampling rates (weekly, daily, monthly, annual), not three.
- `olren.nc` runs to 4 September 2026, but September has fewer than the 15 daily
  values required, so the usable monthly OLR record ends in August 2026.

## Not verifiable here

Two manuscript statements come from Song et al. (2026), ref. 15:

1. about 91% of the oceanic asymmetry is amplitude-driven;
2. the same key periods emerge from the oceanic heat budget.

Only the abstract of that article could be accessed, and neither statement
appears there. Both remain **unverified** and are flagged with comments in the
tracked manuscript for the authors to check against the full article.
