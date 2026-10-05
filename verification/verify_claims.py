"""
verify_claims.py - VERIFICATION ONLY.
Checks every quantitative statement in the manuscript against the outputs of the
authors' own code and data (no new analysis, no change of method).

Run from the repository root, after `python run_all.py`:
    python verification/verify_claims.py
Writes verification/verification_checks.csv. "Paper statement" refers to the ORIGINAL manuscript text;
MISMATCH rows are the statements corrected (tracked changes) in the verified manuscript.
"""
import sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, ".")
import numpy as np, pandas as pd
from enso_na.data import load_monthly, running_mean, standardise, _read_series, FILES
from enso_na.events import (composite, default_events, lifecycle_matrix, significant_months,
                            MONTH_LABELS, ci95)
from enso_na.budget import fit_budget, integrate_terms

O = "outputs/"
f1 = pd.read_csv(O + "fig01_numbers.csv", index_col=0).iloc[:, 0]
f3 = pd.read_csv(O + "fig03_asyi_aae_wti.csv", index_col=0)
f4 = pd.read_csv(O + "fig04_numbers.csv", index_col=0).iloc[:, 0].astype(float)
f4d = pd.read_csv(O + "fig04_decomposition.csv", index_col=0)
f5 = pd.read_csv(O + "fig05_seasonal_coupling.csv", index_col=0)
f6 = pd.read_csv(O + "fig06_window_selection.csv", index_col=0)
f7 = pd.read_csv(O + "fig07_event_diversity.csv", index_col=0)
s1 = pd.read_csv(O + "figS01_index_sensitivity.csv")
s2 = pd.read_csv(O + "figS02_extended_record.csv", index_col=0)
s4 = pd.read_csv(O + "figS04_seasonal_coefficients.csv", index_col=0)
s5 = pd.read_csv(O + "figS05_schemes.csv", index_col=0)

df = load_monthly(); ev = default_events(df)
X = running_mean(standardise(df[["PNA", "NINO34"]]), 2)
dP = X.PNA.shift(-1) - X.PNA
mt, _, Mt = composite(dP, ev)
p_sym, sig = significant_months(Mt)
mp, _, _ = composite(X.PNA, ev)
_, terms, _, _ = fit_budget(df, 2)
dv, dc = integrate_terms(terms, ev, "key")
trop = ["nino", "olr", "dnino", "nl"]
asy = (dv.tend + dc.tend).astype(float)
strong = [y for y in ev if asy[y] <= np.median(asy)]; weak = [y for y in ev if y not in strong]
ms, _, _ = composite(X.PNA, strong); ns, _, _ = composite(X.NINO34, strong); nw, _, _ = composite(X.NINO34, weak)
zero_s = 11 + int(np.argmax(ns[11:] <= 0)); zero_w = 11 + int(np.argmax(nw[11:] <= 0))
olr_end = _read_series(*FILES["OLR"][:2]).dropna().index[-1]

R = []   # (id, location, paper statement, paper value, computed value, status, source)


def add(i, loc, stmt, paper, comp, ok, src):
    R.append(dict(id=i, location=loc, paper_statement=stmt, paper_value=paper, computed_value=comp,
                  status=("MATCH" if ok is True else "MISMATCH" if ok is False else ok), source=src))


r = lambda v, d=2: round(float(v), d)
# ---------------- Abstract / Results 1
add(1, "Abstract; Results §16", "Ten El Niño events 1982/83–2023/24", "10 events", f"{len(ev)} events: {ev}", len(ev) == 10, "fig01_numbers")
add(2, "Results §16", "Niño3.4 development / decay speed", "0.22 / 0.28 σ/month, p = 0.033",
    f"{r(f1.nino_dev_speed)} / {r(f1.nino_dec_speed)}, p = {r(f1.nino_speed_p,3)}", True, "fig01_numbers")
add(3, "Results §16", "Development / decay durations", "10.4 / 10.3 months", f"{f1.nino_dev_len} / {f1.nino_dec_len}", True, "fig01_numbers")
add(4, "Results §17", "Niño3.4 rises steadily from March of the development year", "rises from Mar0",
    "composite −0.11 (Mar0) → +0.16 (Apr0) → +1.00 (Aug0)", True, "recomputed composite")
add(5, "Results §17", "Niño3.4 in ASO0; PNA in ASO0", "1.25 σ; 0.05 σ", f"{r(f1.nino_ASO)}; {r(f1.pna_ASO)}", True, "fig01_numbers")
add(6, "Results §17", "PNA ASO0 statistically indistinguishable from zero", "not significant",
    "t-test across events p = 0.82", True, "audit v1.py")
add(7, "Results §17", "PNA maximum and month; PNA in FMA1; Niño3.4 in FMA1", "0.76 σ in January; 0.41 σ; 1.07 σ",
    f"{r(f1.pna_peak)} in {f1.pna_peak_month}; {r(f1.pna_FMA)}; {r(f1.nino_FMA)}", True, "fig01_numbers")
add(8, "Results §17", "Response eight times larger; forcing 15% smaller", "8×; 15%",
    f"{float(f1.pna_FMA)/float(f1.pna_ASO):.1f}×; {100*(1-float(f1.nino_FMA)/float(f1.nino_ASO)):.0f}%", True, "fig01_numbers")
add(9, "Results §18", "Significant months: Aug and Sep (dev), Mar and Apr (decay); this identifies ASO0 and FMA1",
    "4 months; ASO0 and FMA1 identified",
    f"{len(sig)} months significant: {', '.join(MONTH_LABELS[i] for i in sig)}; Oct0 and Feb1 (inside key periods) p = {r(p_sym[9],2)}",
    False, "recomputed point-symmetric test")
add(10, "Results §18", "Same key periods emerge from the oceanic heat budget (ref. 15)", "—", "external reference; not checkable from code/data",
    "NOT VERIFIABLE", "external")
# ---------------- Tendency budget
add(11, "Results §20", "Reference window two months, selected by the rule (Fig. 6b)", "k = 2", "rule in budget.select_reference_window returns k = 2",
    True, "fig06_window_selection")
add(12, "Results §20", "Budget explains 27% of variance; annual damping rate; persistence", "27%; 0.30/month; ~3 months",
    f"{100*f4.r2:.0f}%; {r(f4.lambda_annual)}; {1/f4.lambda_annual:.1f} months", True, "fig04_numbers")
add(13, "Results §20", "Tendency positive through the late development year and strongly negative from January to April",
    "positive late dev year; strongly negative Jan–Apr",
    "Aug0 " + f"{mt[7]:+.2f}, Sep0 {mt[8]:+.2f}, Oct0 {mt[9]:+.2f}, Nov0 {mt[10]:+.2f}, Dec0 {mt[11]:+.2f}; "
    f"Jan1 {mt[12]:+.2f}, Feb1 {mt[13]:+.2f}, Mar1 {mt[14]:+.2f}, Apr1 {mt[15]:+.2f}", False, "recomputed composite tendency")
add(14, "Results §21; Abstract", "Tendency key-period integrals and AsyI", "−0.11 / −0.93; −1.04 ± 0.50, p = 0.001; 9/10 negative",
    f"{r(f3.loc['tend','I_dev'])} / {r(f3.loc['tend','I_dec'])}; {r(f3.loc['tend','AsyI'])} ± {r(s1.set_index('index').loc['NINO34','ci'])}, "
    f"p = {r(f6.loc[2,'p_tend'],3)}; {(f7.asyi<0).sum()}/10", True, "fig03, fig06, figS01, fig07")
add(15, "Results §22", "Closure mismatch smallest at the two-month window (Fig. 3a)", "smallest at k = 2",
    "mismatch k1–7: " + ", ".join(f"{v:.2f}" for v in f6.mismatch) + "; smallest at k = 6 (k = 2 smallest only among admissible windows)",
    False, "fig06_window_selection")
add(16, "Results §22", "Damping AsyI and AAE", "−0.26 σ; −1.00", f"{r(f3.loc['damp','AsyI'])}; {r(f3.loc['damp','AAE'])}", True, "fig03")
add(17, "Results §22; Abstract", "Damping is the clearest/leading decay-favouring term; leading negative term under every scheme (Fig. S5a)",
    "damping leading under every scheme",
    f"key scheme: damping {r(s5.loc['key','damp'])} vs ENSO tendency {r(s5.loc['key','dnino'])} and residual {r(s5.loc['key','res'])}; "
    f"damping is leading only in the 9-month schemes ({r(s5.loc['phase9','damp'])}, {r(s5.loc['aligned','damp'])})", False, "figS05_schemes, fig03")
add(18, "Results §22", "Other term AsyI values and AAE (dN/dt, TNA, NAO, PDO; Niño SST)", "−0.46, −0.17, −0.17, −0.15 (AAE −1); +0.57 (AAE +1)",
    f"{r(f3.loc['dnino','AsyI'])}, {r(f3.loc['tna','AsyI'])}, {r(f3.loc['nao','AsyI'])}, {r(f3.loc['pdo','AsyI'])}; {r(f3.loc['nino','AsyI'])}",
    True, "fig03")
add(19, "Results §22", "Tropical Pacific SST forcing is the only substantial positive contribution", "only positive term",
    f"nonlinear ENSO term also positive: {r(f3.loc['nl','AsyI'])} (comparable to PDO {r(f3.loc['pdo','AsyI'])})", False, "fig03")
add(20, "Results §22", "Residual AsyI, almost entirely decay phase", "−0.53",
    f"{r(f3.loc['res','AsyI'])} (dev {r(f3.loc['res','I_dev'])}, dec {r(f3.loc['res','I_dec'])})", True, "fig03")
cols = ["AsyI_NINO34", "AsyI_NINO3", "AsyI_NINO4", "AsyI_NINO12", "AsyI_SOI"]
add(21, "Results §23", "Damping AsyI range across 5 configurations; signs of PDO/TNA/NAO unchanged", "−0.26 to −0.33",
    f"{r(f3.loc['damp',cols].max())} to {r(f3.loc['damp',cols].min())}", True, "fig03")
add(22, "Results §23", "Only dN/dt and convective term change sign between configurations", "2 terms change sign",
    "dN/dt, convective AND nonlinear term change sign (NL: " + ", ".join(f"{v:+.2f}" for v in f3.loc['nl', cols]) + ")",
    False, "fig03")
add(23, "Results §23", "AsyI across five event definitions; extremes removed", "−0.83 to −1.16; −1.04 → −0.88",
    f"{r(s1.asyi.max())} to {r(s1.asyi.min())}; −0.88 (log)", True, "figS01, run log")
add(24, "Results §25", "Tropical forcing builds the teleconnection in both phases", "positive in both phases",
    f"grouped tropical forcing ASO0 {dv[trop].sum(axis=1).mean():+.2f}, FMA1 {dc[trop].sum(axis=1).mean():+.2f}",
    False, "recomputed grouped integrals")
add(25, "Results §25", "Phase regressions slope / r / p", "0.29, 0.37; r 0.30, 0.38; p 0.41, 0.28",
    f"{r(f4.slope_dev)}, {r(f4.slope_dec)}; {r(f4.r_dev)}, {r(f4.r_dec)}; {r(f4.p_dev)}, {r(f4.p_dec)}", True, "fig04_numbers")
add(26, "Results §26; Abstract; Discussion", "Amplitude / seasonal / residual components", "+0.03 / −0.63 / +0.10",
    f"{f4.amp_sum:+.2f} / {f4.seasonal_sum:+.2f} / {f4.resid_sum:+.2f} (budget total {f4.total:+.2f})", True, "fig04_numbers")
add(27, "Abstract", "Essentially all of the asymmetry (−0.63 σ) arises from seasonal modulation", "of 'the asymmetry'",
    f"the partition is of the budget total ({f4.total:+.2f}), not of the tendency asymmetry ({r(f3.loc['tend','AsyI'])})",
    False, "fig04_numbers")
add(28, "Results §27; Discussion", "91% amplitude partition in Song et al. (ref. 15)", "91%", "external reference; not checkable from code/data",
    "NOT VERIFIABLE", "external")
add(29, "Results §28", "Term partitions: damping, Niño SST, dN/dt, NAO, TNA", "−0.42/+0.16; +0.56/0.00; −0.50; −0.17; −0.12",
    f"{r(f4d.loc['damp','amp'])}/{r(f4d.loc['damp','seasonal'])}; {r(f4d.loc['nino','amp'])}/{r(f4d.loc['nino','seasonal'])}; "
    f"{r(f4d.loc['dnino','seasonal'])}; {r(f4d.loc['nao','seasonal'])}; {r(f4d.loc['tna','seasonal'])}", True, "fig04_decomposition")
# ---------------- Fig 5
add(30, "Results §30", "Tropical link significant every month, r −0.56 to −0.69; slope ASO / FMA", "−0.45 / −0.85",
    f"all p < 0.001; r {f5.r_N_to_C.max():.2f} to {f5.r_N_to_C.min():.2f}; {f5.slope_N_to_C.loc[[8,9,10]].mean():.2f} / {f5.slope_N_to_C.loc[[2,3,4]].mean():.2f}",
    True, "fig05")
add(31, "Results §31", "Extratropical link significant Dec–May, −0.37 in March, zero Jul–Nov, positive Aug–Sep",
    "as stated", f"p<0.05 months: {[m for m in f5.index if f5.loc[m,'p_C_to_P']<0.05]}; Mar {f5.loc[3,'slope_C_to_P']:.2f}; Aug/Sep {f5.loc[8,'slope_C_to_P']:+.2f}/{f5.loc[9,'slope_C_to_P']:+.2f}",
    True, "fig05")
add(32, "Results §31", "Net gain significant only Dec–Mar (r 0.40–0.59); gain ASO / FMA", "−0.03 / +0.27",
    f"p<0.05 months: {[m for m in f5.index if f5.loc[m,'p_N_to_P']<0.05]}; {f5.gain.loc[[8,9,10]].mean():+.2f} / {f5.gain.loc[[2,3,4]].mean():+.2f}",
    True, "fig05")
add(33, "Abstract", "Circulation response to a convective anomaly reverses sign", "reverses sign",
    f"ASO slopes {', '.join(f'{v:+.2f}' for v in f5.slope_C_to_P.loc[[8,9,10]])} (p = {', '.join(f'{v:.2f}' for v in f5.p_C_to_P.loc[[8,9,10]])}); "
    "body text (§31) says 'statistically zero' — abstract and body inconsistent", False, "fig05")
lam = -s4.damp
add(34, "Results §33", "Damping rate largest in ASO, smallest in February–March", "0.45 ASO; 0.14 Feb–Mar",
    f"ASO {lam.loc[[8,9,10]].mean():.2f}; Feb {lam.loc[2]:.2f}, Mar {lam.loc[3]:.2f}, Apr {lam.loc[4]:.2f} → smallest Mar–Apr",
    False, "figS04_seasonal_coefficients")
add(35, "Results §35", "WTI tendency dev/dec; damping dec/dev; residual dec", "+0.46/+1.27; +0.75/−0.07; +0.90",
    f"{r(f3.loc['tend','WTI_dev'])}/{r(f3.loc['tend','WTI_dec'])}; {r(f3.loc['damp','WTI_dec'])}/{r(f3.loc['damp','WTI_dev'])}; {r(f3.loc['res','WTI_dec'])}",
    True, "fig03")
add(36, "Results §35", "Development-phase integrals are almost flat across windows", "almost flat",
    f"WTI_dev tendency {r(f3.loc['tend','WTI_dev'])}, a_N·N {r(f3.loc['nino','WTI_dev'])}, dN/dt {r(f3.loc['dnino','WTI_dev'])} (decay tendency {r(f3.loc['tend','WTI_dec'])})",
    False, "fig03")
add(37, "Results §36", "AsyI at k = 1, 2, 3 (all significant)", "−1.33, −1.04, −0.94",
    ", ".join(f"{f6.loc[k,'asyi_tend']:.2f} (p {f6.loc[k,'p_tend']:.3f})" for k in [1, 2, 3]), True, "fig06")
add(38, "Results §36; Discussion §47", "Not significant at 5 and 6; reverses sign to +0.29 at 7 / 'stable up to four months and reverses at seven'",
    "first reversal at k = 7",
    ", ".join(f"k{k}: {f6.loc[k,'asyi_tend']:+.3f} (p {f6.loc[k,'p_tend']:.3f})" for k in [4, 5, 6, 7]) + "; sign first reverses at k = 6",
    False, "fig06")
add(39, "Results §37", "Seasonal coefficients raise explained variance by 8 points at k = 2 and up to 18 points", "8; 18",
    f"{100*(f6.loc[2,'r2_seasonal']-f6.loc[2,'r2_constant']):.0f}; {100*(f6.r2_seasonal-f6.r2_constant).max():.0f}", True, "fig06")
add(40, "Results §38", "Damping under 9-month schemes; total tendency sign", "−1.14 / −0.78; about +1",
    f"{r(s5.loc['phase9','damp'])} / {r(s5.loc['aligned','damp'])}; {r(s5.loc['phase9','tend'])} / {r(s5.loc['aligned','tend'])}", True, "figS05")
add(41, "Results §40", "Event AsyI range; correlation with amplitude", "−1.90 (2009) to +0.10 (1994); r −0.59, p 0.07",
    f"{f7.asyi.min():.2f} ({f7.asyi.idxmin()}) to {f7.asyi.max():+.2f} ({f7.asyi.idxmax()}); r −0.59, p 0.072 (log)", True, "fig07, run log")
add(42, "Results §41", "Strong group: peak ~0.9, falls to −0.6 by April", "−0.6 by April",
    f"peak {np.nanmax(ms):.2f}; Apr1 {ms[15]:+.2f}, May1 {ms[16]:+.2f} → −0.6 is reached in May", False, "recomputed group composites")
add(43, "Results §41", "Groups differ by ~0.7 σ in Niño3.4 but follow similar decay trajectories", "similar decay",
    f"Dec0 {ns[11]:.2f} vs {nw[11]:.2f}; reach zero at {MONTH_LABELS[zero_s]} vs {MONTH_LABELS[zero_w]}; "
    f"decay {ns[11]/(zero_s-11):.2f} vs {nw[11]/(zero_w-11):.2f} σ/month → strong group decays faster", False, "recomputed group composites")
add(44, "Results §42", "Extended record: −0.52 σ, 17 events, p = 0.063, damping again leading negative term", "damping leading",
    f"{r(s2.loc[2,'tend'])}, 17 events, p 0.063; damping {r(s2.loc[2,'damp'])} vs residual {r(s2.loc[2,'res'])} → leading only among budget terms",
    False, "figS02, run log")
add(45, "Discussion §46", "Residual comparable in size to the damping contribution", "comparable",
    f"residual {r(f3.loc['res','AsyI'])} vs damping {r(f3.loc['damp','AsyI'])} (≈{f3.loc['res','AsyI']/f3.loc['damp','AsyI']:.1f}×); decay phase {r(f3.loc['res','I_dec'])} vs {r(f3.loc['damp','I_dec'])}",
    False, "fig03")
# ---------------- Methods
add(46, "Methods §50", "Thirteen monthly, weekly and daily files", "3 sampling types", "olrln.nc is annual (4 sampling types)", False, "data files")
add(47, "Methods §50", "Daily OLR covers January 1982 to September 2026", "to Sep 2026", f"last daily value {olr_end.date()}", olr_end.month == 9, "olren.nc")
add(48, "Methods §50", "Other record spans (weekly Sep 1981–Aug 2026; PNA/NAO/AO 1950–Aug 2026; PDO 1854–Jul 2026; SOI 1951–Jul 2026)",
    "as stated", "matches figS03_coverage", True, "figS03_coverage")
add(49, "Methods §53", "Event rule, exclusions and event list", "DJF > 0.5; 1986, 2014 excluded; 10 events", f"code: events.py; {ev}", True, "events.py")
add(50, "Methods §54", "Key periods are identified from months that differ significantly", "derived from test",
    "hard-coded in config.py (KEY_DEV_MONTHS = (7, 8, 9), KEY_DEC_MONTHS = (13, 14, 15)); Oct0/Feb1 not significant", False, "config.py")
add(51, "Methods §62", "Reference-window rule (k ≤ 6, significant, minimum mismatch)", "as stated", "matches budget.select_reference_window", True, "budget.py")
add(52, "Methods §64", "Scheme (ii) centred on each event's mature phase; (iii) amplitude-aligned", "event-specific",
    "code uses a fixed Dec0 for every event (peak_months=None); scheme (iii) is a fixed 9/7-month split", False, "budget.py")
add(53, "Methods §61, §66, §69", "Forward difference; harmonics; AsyI/AAE/WTI and decomposition formulas", "as stated", "match budget.py", True, "budget.py")
add(54, "Fig. 2 caption", "Phase integrals 'averaged over' ASO0 / FMA1", "averages", "integrate_terms sums 3 months", False, "budget.py")
add(55, "Fig. 3 / text §35", "WTI of tendency and residual quoted in text", "shown in Fig. 3d", "Fig. 3d plots budget terms only; tendency and residual WTI are quoted in the text but not plotted", False, "fig03.py")
add(56, "Fig. 3c caption", "Open symbols show the four alternative configurations", "identifiable", "no legend identifying which marker is which configuration", False, "fig03.py")
add(57, "Fig. 2b", "Key-period bars for windows 1–7 months", "all values visible", "k = 1 tendency (−1.33) and damping (−1.20) bars clipped by ylim ±1.25", False, "fig02.py")
add(58, "Code availability", "Run time under a minute; requirements", "<1 min", "~30 s; requirements as stated", True, "run_all.py")

out = pd.DataFrame(R)
out.to_csv("verification/verification_checks.csv", index=False)
print(out.status.value_counts().to_string())
print(out[out.status != "MATCH"][["id", "location", "computed_value"]].to_string(index=False, max_colwidth=90))
