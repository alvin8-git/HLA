# Version History

## v2.4.0 — 2026-06-18

### Fixed — peer-review correctness pass (text-only, no recomputation)
- **C1 number drift**: removed EM-MLE values that leaked into prose where the
  bootstrap median is the reported estimate — §3.4 Chinese 42,871→42,847;
  Rec 4 Others 32,360→31,181; §3.5 per-group range ~32,000–45,000→~31,000–44,000.
  All in-text registry figures now agree with `registry_size_ci.csv` (Table 1/2).
- **C2 overclaim**: abstract and Recommendations 1–2 now flag Malay/Indian/Others
  N* as model projections pending validation; only Chinese is empirically validated
  (Spearman r=0.70, §3.6). "Mathematically necessary" softened in Rec 1 and §3.3.
- **C3 CI honesty**: bootstrap-CI lower-bound caveat (sampling variability only;
  excludes EM phasing and HWE model error) promoted from §2.4 into the abstract
  and Table 1/2 captions.
- **C5 EM citation/method**: replaced miscited Beatty [4] for EM phasing with
  Excoffier & Slatkin 1995 [20]; method now correctly described as full multi-locus
  phase-enumeration EM (matches `hwe_test._em_full_phase`), haplotypes retained ≥0.1%.
- `python-docx>=1.1` added to `analysis/requirements.txt` (report pipeline dep).

### Deferred (scope decision pending)
- C4 rare-haplotype cutoff sensitivity; broader patient validation for minority
  groups; Others cluster-stability bootstrap; external face-validity paragraph.

### Generated
- `HLA_Registry_Size_CMIO_v2.14.docx`

---

## v2.3.0 — 2026-04-29

### Added — quantified bias analyses (reviewer final polish)
- **`analysis/15_em_convergence.py`**: EM convergence test — reruns EM for Chinese
  at 500–45,018 donors; N* at 5k cap = 45,148 vs 41,727 at full sample (8.2%
  conservative overestimate). Figure S1 added to Supporting Analysis section.
- **`analysis/16_smoothing_sensitivity.py`**: Laplace pseudocount smoothing
  sensitivity — α=0.001 per haplotype; N* at 95% changes <3% for all groups
  (Chinese +0.9%, Malay +2.3%, Indian −3.1%, Others −1.9%); larger at 75%.
- **Tables 1 & 2**: Added "Signed-up target‡" row — N* ÷ 0.60 (40% attrition);
  shows range across CMIO groups per threshold.
- **§4.1 Limitations**: EM cap quantified (8.2% conservative bias); smoothing
  results cited; attrition adjustment formula stated.
- **Supporting Analysis** section with Figure S1 (EM convergence).

### Generated
- `HLA_Registry_Size_CMIO_v2.13.docx`

---

## v2.2.0 — 2026-04-29

### Changed
- **Figure 1 flowchart**: fixed arrow alignment — arrows now run between boxes
  (bottom of one box to top of the next) rather than inside boxes.
- **Figure numbering**: renumbered sequentially in document order —
  Fig 1 pipeline, Fig 2 CI bar chart (was unlabelled), Fig 3 10-locus partial match,
  Fig 4 8-locus partial match, Fig 5 sensitivity, Fig 6 validation scatter,
  Fig 7 Others PCA scatter (was Fig 1).

### Generated
- `HLA_Registry_Size_CMIO_v2.12.docx`

---

## v2.1.0 — 2026-04-29

### Changed — flow and comprehension improvements
- **Figure 0** (new): Methods pipeline flowchart (`analysis/14_pipeline_flowchart.py`)
  showing data → EM → HWE → C(N) → N* → bootstrap CI.
- **§2.3 Eq(4)**: Added verbal explanation paragraph — what C(N) means intuitively
  and why rare diplotypes make convergence slow.
- **§3 Results intro**: Added 7-line reading guide orienting the reader to which
  sections are primary, secondary, robustness, and exploratory.
- **§3 section order**: Reordered §3.4–3.7 — Partial Match (§3.4), Sensitivity (§3.5),
  Validation (§3.6), Others Exploratory (§3.7). Others clearly labelled as exploratory.
- **Table 1 footnote**: Weighted Average row explained as a mathematical convenience,
  not a policy target; per-group targets are the operative planning figures.
- **Glossary** (new table): 15 abbreviations defined (AFND, BMDP, CI, CMIO, EM,
  HSCT, HLA, HSA, HWE, LD, MLE, N*, PCA, RMSE, SCBB).
- All §3.x cross-references updated to match new numbering.

### Generated
- `HLA_Registry_Size_CMIO_v2.11.docx`

---

## v2.0.0 — 2026-04-29

### Changed — reviewer response (text-only, no recomputation)
- **§2.4 Bootstrap CI**: removed trivial "by construction" sentence; replaced with
  explicit scope statement — CIs capture donor-count sampling variability only, not
  EM phasing or HWE model uncertainty.
- **§3.7 Validation**: explicitly flagged Indian (1 shared haplotype = no validation)
  and stated that Malay/Indian/Others estimates are model-derived projections;
  Chinese is the primary validated result.
- **§4.1 Limitations**: expanded from 5 to 7 substantive points —
  (1) HWE bias direction uncertain; Indian/Others flagged as exploratory;
  (2) 5,000 cap binds materially only for Chinese — common haplotypes robust at 5k;
  (3) NEW: donor attrition — N* is biologically matched minimum; signed-up targets
      must exceed N* by 30–50% to account for real-world attrition;
  (4) NEW: N* is a lower bound — unobserved haplotypes assigned zero frequency,
      most material at 95% coverage threshold;
  (5) Others cluster: added note that silhouette reflects strong HLA–ancestry signal
      but cluster stability was not bootstrap-validated.

### Generated
- `HLA_Registry_Size_CMIO_v2.10.docx`

---

## v1.9.0 — 2026-04-29

### Changed
- **Table 4 & Table 5 cluster colours**: updated Cluster 2 and Cluster 3 cell backgrounds
  to match Figure 1 scatter colours exactly — light-blue `CCE0F5` (→ #377eb8) and
  light-green `CBF0CB` (→ #4daf4a). All three cluster colours now correspond to Figure 1.

### Generated
- `HLA_Registry_Size_CMIO_v2.9.docx`

---

## v1.8.0 — 2026-04-28

### Changed
- **Figure 1** (`analysis/11_others_stratification.py`): removed silhouette subplot;
  single 9×6 panel with ancestry labels in legend; cluster colours unchanged.
- **Table 5**: reduced to top-1 haplotype per cluster (was top-2); removed Rank column;
  Population association column widened to 7.5 cm.
- **Cluster 1 table colour**: changed from steel-blue `D6DCE4` → light-red `FFCCCC`
  to match Figure 1 red scatter colour for Cluster 1 (European/Eurasian).

### Generated
- `HLA_Registry_Size_CMIO_v2.8.docx`

---

## v1.7.0 — 2026-04-28

### Changed
- **§3.4 Others**: restored Figure 1 (PCA k-means scatter; knee plot omitted) and
  a condensed Table 5 showing top 2 haplotypes per cluster (ancestry validation).
  Removed the 5-row-per-cluster version and per-cluster narrative paragraphs.
- **§5 Conclusions**: softened Others references — from "requires particular attention"
  with alarming sub-cluster numbers to a brief parenthetical note (Table 4). Point (3)
  reworded to "ancestry sub-group data collection as registry grows".

### Generated
- `HLA_Registry_Size_CMIO_v2.7.docx`

---

## v1.6.0 — 2026-04-28

### Changed
- **§2.4 Bootstrap CI**: removed redundant left-skew technical paragraph (Jensen's
  inequality explanation); reasoning already stated in §2.4 para 1. Kept clinical
  planning sentence.
- **§3.4 "Others" subgroup**: condensed from full primary section (Table 4 + Table 5 +
  Figure 1 + Figure 2 + 3 cluster narratives) to a single brief paragraph + Table 4 only.
  Heading renamed to "Note on the 'Others' Subgroup". Cluster narratives and haplotype
  evidence table retained in code for reference.
- **Patient data source**: renamed "Actual patients" → "HSA Patient-Donor Data
  (Health Sciences Authority Singapore)" throughout Table 6, Figure 5, and §3.7.

### Generated
- `HLA_Registry_Size_CMIO_v2.6.docx`

---

## v1.5.0 — 2026-04-28

### Changed
- **Figures 3 & 4** (`analysis/06_partial_match_plots.py`): removed "Overall" panel;
  reformatted from 1×5 to 2×2 grid (Chinese+Malay top row, Indian+Others bottom row);
  larger figure size (14×10 in) for readability.
- **Tables 1 & 2**: renamed "Combined†" → "Weighted Average†" with explicit disclaimer
  that the weighted average does not guarantee equitable access for minority groups.
- **Cluster 3 name**: standardised to "Northeast Asian / Mixed" across Table 4 and
  all §3.4 narrative text (was "NE Asian / mixed" in table vs "Northeast Asian" in text).

### Generated
- `HLA_Registry_Size_CMIO_v2.5.docx`

---

## v1.4.0 — 2026-04-28

### Fixed
- **Bootstrap CI coverage** (`analysis/09_bootstrap_ci.py`): replaced EM MLE
  point estimate with bootstrap median (bias-corrected for Jensen's inequality
  on the concave N*(f) near saturation). All 32 point estimates now fall within
  their 95% CIs by construction.
- **n_eff**: changed from capped 5,000 to actual 5-locus donor counts
  (Chinese: 45,754; Malay: 5,868; Indian: 5,586; Others: 3,941).
- **B**: increased from 500 to 1,000 bootstrap resamples.
- EM MLE estimates preserved in CSV as `em_registry_size` for reference.

### Generated
- `HLA_Registry_Size_CMIO_v2.3.docx` — updated §2.4 methodology and §3.1
  narrative; all CI tables use bootstrap median values.

---

## v1.3.0 — 2026-04-20

### Added
- **Partial match coverage curves** (`analysis/06_partial_match_plots.py`)
  - 10-locus framework: 8/10, 9/10, 10/10 match levels for each CMIO group + Overall
  - 8-locus framework: 6/8, 7/8, 8/8 match levels for each CMIO group + Overall
  - Style matches WBMT (Aljurf et al. 2019) Figure 5
  - Model uses haplotype-pair enumeration capturing linkage disequilibrium (LD)
- `analysis/figures/partial_match_10locus.png`
- `analysis/figures/partial_match_8locus.png`

### Fixed
- Partial match model initially used per-locus allele frequencies (ignoring LD),
  giving Chinese 10/10 coverage of 2.9% at N=1M. Fixed to use EM haplotype-pair
  enumeration — Chinese 10/10 at N=1M now correctly reaches ~100%.
- x-axis adjusted to 1,000–10,000,000 so the S-curve body is visible in all panels.

---

## v1.2.0 — 2026-04-20

### Added
- `README.md` — badges, pipeline diagram, results tables, quick start, figure index
- `Documentation.md` — full technical documentation with HLA biology background,
  EM algorithm derivation, HWE theory, registry model math, figure interpretation
- `LICENSE` — MIT
- Figures embedded inline in `Documentation.md` (GitHub renders automatically)

---

## v1.1.0 — 2026-04-17

### Added
- `analysis/05_report.py` → `analysis/verification_summary.md` (580 lines, 6 sections)
- `analysis/run_all.sh` verified end-to-end; all 12 output files produced

### Pipeline outputs (complete)
| File | Description |
|------|-------------|
| `analysis/data/hla_clean.csv` | 305,745 rows, 61,149 samples |
| `analysis/data/allele_freq_comparison.csv` | 1,488 alleles, 0 flagged |
| `analysis/data/haplo_freqs_em.csv` | 251 haplotypes (≥0.1% freq) |
| `analysis/data/allele_freqs_per_locus.csv` | Per-locus allele frequencies |
| `analysis/data/hwe_results.csv` | 20 tests; 8 violations |
| `analysis/data/coverage_curves.csv` | 3,600 rows (200 N points × 18 scenarios) |
| `analysis/data/registry_size_targets.csv` | 72 rows (4 thresholds × 18 scenarios) |
| `analysis/figures/allele_freq_heatmap.png` | Discrepancy heatmap |
| `analysis/figures/coverage_curves_8of8.png` | Exact 8/8 coverage curves |
| `analysis/figures/coverage_curves_10of10.png` | Exact 10/10 coverage curves |

---

## v1.0.0 — 2026-04-17

### Initial pipeline implementation

**IMPL-1: Project setup**
- `analysis/requirements.txt`, `analysis/run_all.sh`

**IMPL-2: Data ingestion** (`analysis/01_ingest.py`, 13 tests)
- Reads BMDP + SCBB Excel sheets + HSA txt files
- Normalises to 2-field resolution, maps CMIO ethnicity codes
- Output: `analysis/data/hla_clean.csv` — 305,745 rows, 61,149 unique samples

**IMPL-3: Allele frequency verification** (`analysis/02_allele_freq.py`, 5 tests)
- Recomputes allele frequencies from BMDP+SCBB; compares to Gene[Rate] published values
- **Result: 0 flagged alleles; max discrepancy 0.27% (HLA-C) — fully reproducible**

**IMPL-4: EM haplotype estimation + HWE tests** (`analysis/03_hwe_test.py`, 5 tests)
- Product-approximation EM for 5-locus haplotype frequencies
- Chi-squared HWE test with Bonferroni correction (p < 0.0025)
- **Result: 8 violations — Indian (DQB1, HLA-B, HLA-C), Others (all 5 loci)**

**IMPL-5: Registry model core math** (`analysis/registry_model.py`, 6 tests)
- `get_diplotype_frequencies`: HWE expansion from haplotype freqs
- `compute_coverage`: Coverage(N) = Σ f_g · [1 − (1 − f_g)^N]
- `find_registry_size`: log-scale binary search for minimum N
- `get_combined_haplotype_freqs`: Singapore-weighted combined pool

**IMPL-6: Full registry model run** (`analysis/04_registry_model.py`)
- N sweep: log-spaced 1,000–10,000,000
- 18 scenarios: 2 match levels × 5 ethnicities × (1–2 model variants)
- 4 coverage thresholds: 75%, 85%, 90%, 95%

**IMPL-7: Report assembly** (`analysis/05_report.py`)
- `analysis/verification_summary.md` — 6-section narrative report

**IMPL-8: End-to-end verification**
- 29/29 tests passing; all outputs reproduced via `run_all.sh`

---

## Test Coverage

| Module | Tests | Status |
|--------|-------|--------|
| `01_ingest.py` | 13 | ✅ |
| `02_allele_freq.py` | 5 | ✅ |
| `03_hwe_test.py` (hwe_test.py) | 5 | ✅ |
| `04_registry_model.py` (registry_model.py) | 6 | ✅ |
| **Total** | **29** | **✅ All passing** |

---

## Key Findings

| Finding | Value |
|---------|-------|
| Allele frequency reproducibility | Max discrepancy 0.27%; **0 alleles flagged** |
| HWE violations | 8/20 tests (Indian: 3 loci; Others: 5 loci) |
| Chinese 10/10 registry (95% coverage) | **11,616 donors** (same-ethnicity) |
| Malay 10/10 registry (95% coverage) | **17,601 donors** (same-ethnicity) |
| Cross-ethnic 10/10 for Malay/Others | **Infeasible** (hits 10M ceiling) |
| Combined registry (95% coverage, 10/10) | **83,541 donors** |
