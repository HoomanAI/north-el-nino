# Why El Niño's North American teleconnection engages late and collapses early (verified code)

[![reproduce](https://github.com/HoomanAI/north-el-nino/actions/workflows/reproduce.yml/badge.svg)](https://github.com/HoomanAI/north-el-nino/actions/workflows/reproduce.yml)

Reproduction package for the manuscript of the same name: a window-resolved
tendency-budget diagnosis of phase asymmetry in the Pacific–North American (PNA)
response to El Niño.

> **Verification copy.** This repository is the authors' original analysis package
> with **display-only** figure fixes. The data, the method and every number are
> unchanged. All 15 output tables are byte-identical to those of the original package.
> Every changed line is marked `# VERIFY (display only)`. See [`VERIFICATION.md`](VERIFICATION.md)
> for the list of changes and for how the manuscript text was checked against the code outputs.

Everything in the paper is produced from the thirteen index files in `data/`.
No network access is needed.

```bash
pip install -r requirements.txt
python run_all.py                         # ~25 s; writes figures/ and outputs/
python run_all.py fig03                   # one figure only
python verification/verify_claims.py      # re-check the 58 manuscript statements (after run_all.py)
```

Figures are written as `figXX_verified.png/.pdf`. Set `ENSO_NA_SUFFIX=` (empty) to
use the original file names. Output paths can be redirected with the environment
variables `ENSO_NA_DATA`, `ENSO_NA_FIGDIR` and `ENSO_NA_OUTDIR`.

## Layout

```
data/                     13 index NetCDF files (inputs)
enso_na/
  config.py               paths, analysis constants, term labels
  data.py                 NetCDF loading, monthly homogenisation, filtering
  events.py               El Niño identification, phases, composites, statistics
  budget.py               tendency budget, AsyI / AAE / WTI, decompositions
  plotting.py             shared figure style
figures/fig01..07.py      main figures  (also write the .csv tables they use)
figures/figS01..05.py     supplementary figures
run_all.py                regenerates everything
download_gridded.py       OPTIONAL: gridded North American composite maps (unchanged)
verification/
  verify_claims.py        checks every quantitative statement in the manuscript
  verification_checks.csv reference result (58 checks)
outputs/                  (generated) every number quoted in the paper, as CSV
```

Each figure script is standalone: `python figures/fig04.py` works if the package
root is on `PYTHONPATH` (`run_all.py` handles that automatically).

## Data

The supplied files carry four sampling rates (weekly, daily, monthly and annual)
and several time-axis names; `data.py` converts all of them to a single
monthly-mean table on a month-start index. Weekly and daily series require at
least 2 and 15 valid samples per month respectively, otherwise the month is set
missing.

| File | Variable | Sampling | Record | Public source |
|---|---|---|---|---|
| `nino12.nc`, `nino3.nc`, `nino34.nc`, `nino4.nc` | Niño SST anomalies | weekly | 1981-09 → 2026-08 | NOAA CPC weekly ENSO indices (OISSTv2), <https://www.cpc.ncep.noaa.gov/data/indices/wksst8110.for> |
| `tna.nc` | Tropical North Atlantic | weekly | 1981-09 → 2026-08 | NOAA PSL, <https://psl.noaa.gov/gcos_wgsp/Timeseries/TNA/> |
| `olren.nc` | equatorial central Pacific OLR anomaly | daily | 1982-01 → 2026-09-04 (usable monthly record ends 2026-08) | NOAA interpolated OLR, <https://psl.noaa.gov/data/gridded/data.olrcdr.interp.html> |
| `olrln.nc` | annual count of suppressed-convection days | annual | 1982 → 2026 | derived from the same daily OLR record |
| `pna.nc`, `nao.nc`, `ao.nc` | teleconnection indices | monthly | 1950-01 → 2026-08 | NOAA CPC, <https://www.cpc.ncep.noaa.gov/data/teledoc/telecontents.shtml> |
| `pdo.nc` | Pacific Decadal Oscillation | monthly | 1854-01 → 2026-07 | NOAA NCEI, <https://www.ncei.noaa.gov/access/monitoring/pdo/> |
| `soi.nc` | Southern Oscillation Index | monthly | 1951-01 → 2026-07 | NOAA CPC, <https://www.cpc.ncep.noaa.gov/data/indices/soi> |
| `amo.nc` | Atlantic Multidecadal Oscillation | monthly | 1948-01 → 2023-01 | NOAA PSL, <https://psl.noaa.gov/data/timeseries/AMO/> |

`olrln.nc` and `amo.nc` are not used in the main budget: the first is an annual
count and the second stops in January 2023, before the decay year of the
2023/24 event. `amo.nc` is used only in the extended-record loader.

### Optional gridded data

`download_gridded.py` adds North American composite maps of 500 hPa height,
surface air temperature and precipitation. It downloads roughly 700 MB from
NOAA PSL and is not needed for any figure in the paper. URLs are listed in the
module docstring; NOAA PSL occasionally reorganises its directories, so verify
the current path at <https://psl.noaa.gov/data/gridded/> if a download fails.

## Method in brief

**Events.** El Niño is declared when the DJF mean of the 3-month running-mean
Niño3.4 anomaly exceeds 0.5 °C. 1986/87 and 2014/15 are excluded because warm
anomalies persisted into the following event, leaving ten events over
1982/83–2023/24.

**Budget.** The tendency of the standardised PNA index is written as

```
dP/dt = -λ(t)·P + a_N(t)·N + a_C(t)·C + a_Ṅ(t)·dN/dt + a_D(t)·D
        + a_T(t)·T + a_A(t)·A + a_M(t)·M + a_L(t)·NL + R
```

with N = Niño3.4, C = OLR, D = PDO, T = TNA, A = AO, M = NAO and NL a nonlinear
ENSO term orthogonalised against N. Every coefficient carries an annual harmonic
and is fitted by ordinary least squares over all months of the record — no
information about the phase asymmetry enters the fit. All terms have units of
PNA standard deviations per month.

Two choices matter. `dP/dt` is a **forward** difference, because for a
first-order autoregressive process the centred difference is exactly orthogonal
to the instantaneous state and cannot identify a damping coefficient. And the
free methodological parameter is the **integration window** k (1–7 months) over
which every field is smoothed — the temporal analogue of the mixed-layer
integration depth in the ocean heat budgets this framework is modelled on.

**Diagnostics.** AsyI (net contribution over the key development and decay
periods), AAE (directional consistency, bounded in ±1) and WTI (sensitivity to
the integration window). The key periods, Aug–Oct of the development year (ASO0)
and Feb–Apr of the decay year (FMA1), are fixed in `enso_na/config.py`. The
reference window k = 2 is selected by a documented rule in
`budget.select_reference_window`, not by hand.

**Uncertainty.** 95% confidence intervals are Student-t intervals across the ten
events. Robustness to the ENSO index is assessed with five configurations
(Niño1+2, Niño3, Niño3.4, Niño4, inverted SOI), which play the role that the
multi-reanalysis ensemble plays in ocean heat-budget studies.

## Known limitations

- Ten events is a small sample. The composite tendency asymmetry is significant
  (p = 0.001) and holds for 9 of 10 individual events, but per-event speed
  comparisons are not significant, and several secondary budget terms change
  sign between ENSO-index configurations.
- The sign of the *total* tendency asymmetry depends on the phase-integration
  scheme (Fig. S5). Over a full phase the tendency integral telescopes to the
  difference between endpoint values, which is dominated by the La Niña-like
  negative PNA of the early development year. The key-period diagnosis compares
  two equal-length, seasonally matched windows and does not have this property.
- The regression budget is statistical. It quantifies which covariates account
  for the tendency asymmetry; it does not establish causality, and the residual
  is a substantial term during the decay phase.
- Two statements in the manuscript cite Song et al. (2026), ref. 15: the "91%"
  partition and the claim that the same key periods emerge from the oceanic heat
  budget. They cannot be checked from this code and data and remain unverified.

## License

The code is released under the [MIT License](LICENSE). The index files in `data/`
are derived from public NOAA products; please cite the original sources listed
above when reusing them.
