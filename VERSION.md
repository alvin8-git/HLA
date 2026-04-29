# Version History

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
