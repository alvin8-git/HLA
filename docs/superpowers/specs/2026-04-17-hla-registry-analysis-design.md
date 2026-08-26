# HLA Registry Analysis Design

> **Historical design spec — April 2026.** Kept as the original design record; see
> [README](../../../README.md) for the pipeline as it now stands.

**Date:** 2026-04-17
**Author:** Alvin Ng Yu-Jin
**Status:** Approved

---

## Overview

Two-part analysis of Singapore HLA typing data (BMDP + SCBB + HSA, n=59,186):

1. **Verification** — Reproduce and audit the 2022 Blood Cell Therapy paper's allele and haplotype frequency analysis, using independent tools to validate Gene[Rate] outputs and test methodological assumptions.

2. **Registry size modeling** — Estimate the optimal total number of donors needed in Singapore's bone marrow registry to achieve target match probabilities (75%, 85%, 90%, 95%) for 8/8 and 10/10 HLA match levels, broken down by CMIO ethnicity and matching model (same-ethnicity vs cross-ethnic).

---

## Data Sources

| File | Contents |
|------|----------|
| `HLA Data combined.xlsx` | Raw BMDP + SCBB data, alleles as A1/A2/B1/B2 etc. |
| `HLA Data.cleaned.xlsx` | Normalized to 2-field resolution |
| `BMDPnSCBB.results.xlsx` | Gene[Rate] allele + haplotype frequency outputs per CMIO group |
| `DonorPatient.txt` | HSA donor haplotypes (1,350 rows, 6 cols: ethnicity + 5 loci) |
| `Patient.txt` | HSA recipient haplotypes (~560 rows, same format) |
| `2022_HLA_BloodCellTherapy.pdf` | Published paper describing original methodology |

---

## Architecture

### Implementation language

**Hybrid Python + R:**
- Python: data ingestion, allele frequency recomputation, registry size modeling, figures
- R: independent EM haplotype frequency estimation (`haplo.stats`), formal HWE tests (`HardyWeinberg`)

### Project structure

```
/data/alvin/HLA/
├── analysis/
│   ├── 01_ingest.py          # Load & normalize Excel + txt data
│   ├── 02_allele_freq.py     # Recompute allele freqs, compare to paper
│   ├── 03_hwe_test.R         # Independent EM + HWE tests
│   ├── 04_registry_model.py  # Registry size coverage curve model
│   ├── 05_report.py          # Assemble tables, figures, verification summary
│   ├── run_all.sh            # Sequential driver script
│   ├── data/                 # Intermediate CSV/parquet outputs
│   └── figures/              # Output plots (PNG/PDF)
└── docs/
    └── superpowers/specs/    # This document
```

### Data flow

```
HLA Data combined.xlsx
    → 01_ingest.py
    → data/hla_clean.csv

data/hla_clean.csv
    → 02_allele_freq.py → data/allele_freq_comparison.csv
    → 03_hwe_test.R    → data/haplo_freqs_haplo_stats.csv
                       → data/hwe_results.csv

BMDPnSCBB.results.xlsx + data/haplo_freqs_haplo_stats.csv
    → 04_registry_model.py → data/coverage_curves.csv
                           → data/registry_size_targets.csv

All outputs → 05_report.py → figures/ + verification_summary.md
```

---

## Component Design

### 01_ingest.py — Data ingestion & normalization

- Reads BMDP and SCBB sheets from `HLA Data combined.xlsx`
- Reads `DonorPatient.txt` and `Patient.txt` (HSA)
- Reshapes from wide (A1, A2, B1, B2...) to tidy long format:

| sample_id | source | ethnicity | locus | allele1 | allele2 |
|-----------|--------|-----------|-------|---------|---------|
| ... | BMDP | Chinese | HLA-A | 11:01 | 24:02 |

- Ethnicity mapped to CMIO codes (C / M / I / O)
- Missing alleles ("-", blank) kept as `NaN` — not dropped
- Source tracked (BMDP / SCBB / HSA-Donor / HSA-Patient)
- Reports missingness rate per (source, ethnicity, locus)

Output: `data/hla_clean.csv`

---

### 02_allele_freq.py — Allele frequency verification

For each (ethnicity, locus):
- Compute observed allele frequency = count(allele) / total typed alleles (NaN excluded)
- Load published frequencies from `BMDPnSCBB.results.xlsx`
- Compute signed difference (observed − published) per allele
- Flag alleles where |difference| > 0.005 (0.5 percentage points)

Outputs:
- `data/allele_freq_comparison.csv`
- `figures/allele_freq_heatmap.png` — discrepancy magnitude across loci and CMIO groups

---

### 03_hwe_test.R — Independent verification

For each (ethnicity, locus):

1. **Independent EM haplotype estimation** via `haplo.stats::haplo.em()` on the allele pairs from `hla_clean.csv`
   - Produces per-ethnicity 5-locus haplotype frequency estimates
   - Compare top-20 haplotypes vs Gene[Rate] output: compute Pearson r and RMSE
   - Generate scatter plot: `haplo.stats` freq vs Gene[Rate] freq (should lie on y=x diagonal)

2. **HWE test** via `HardyWeinberg::HWExact()` per locus per ethnicity
   - Bonferroni-corrected significance threshold (p < 0.05 / 20 tests)
   - Report any loci deviating from HWE (violation of core Gene[Rate] assumption)

Outputs:
- `data/haplo_freqs_haplo_stats.csv`
- `data/hwe_results.csv`
- `figures/haplo_scatter_[ethnicity].png`

---

### 04_registry_model.py — Registry size model

**Mathematical framework:**

Haplotype frequencies h_i from `haplo_freqs_haplo_stats.csv` (independently estimated).

**Step 1 — Diplotype frequencies under HWE:**
```
f(h_i, h_i) = h_i²
f(h_i, h_j) = 2 · h_i · h_j   (i ≠ j)
```
Enumerate top K haplotypes per group, where K is the minimum number covering ≥99% cumulative frequency. Residual frequency is pooled as a single "other" type.

**Step 2 — Per-patient match probability:**

For a patient with genotype g at diplotype frequency f_g, given a registry of N donors:
```
P(≥1 match | N) = 1 − (1 − f_g)^N
```

**Step 3 — Population coverage:**
```
Coverage(N) = Σ_g  f_g · [1 − (1 − f_g)^N]
```
Sweep N from 1,000 to 5,000,000 on a log scale.

**Step 4 — Four model variants per ethnicity:**

| Variant | Donor pool | Patient pool |
|---------|-----------|--------------|
| Same-ethnicity | per-group haplotype freqs | per-group |
| Cross-ethnic | Singapore-weighted combined freqs | per-group |

Singapore population weights (approximate, from BMDP composition): Chinese 77%, Malay 8%, Indian 9%, Others 6%.

**Match levels:**
- **8/8**: HLA-A, B, C, DRB1 only (drop DQB1 from haplotype)
- **10/10**: all 5 loci (A, B, C, DRB1, DQB1)

**Target coverage thresholds:** 75%, 85%, 90%, 95%

Outputs:
- `data/coverage_curves.csv` — columns: N, coverage, match_level, ethnicity, model_variant
- `data/registry_size_targets.csv` — registry N needed per (threshold × match_level × ethnicity × variant)
- `figures/coverage_curves_8of8.png`
- `figures/coverage_curves_10of10.png`

Each figure: panel per ethnicity + combined panel; lines for same-ethnicity vs cross-ethnic; horizontal dashed lines at target thresholds.

---

### 05_report.py — Summary report assembly

Produces `verification_summary.md` containing:

1. **Allele frequency reproducibility** — was the paper's allele frequency calculation reproducible? Max discrepancy found, any systematic bias by source or ethnicity.

2. **HWE assessment** — list of any loci/groups deviating significantly from HWE. Implication: if HWE is violated, Gene[Rate]'s haplotype frequency estimates carry additional uncertainty.

3. **Haplotype frequency agreement** — Pearson r and RMSE between `haplo.stats` and Gene[Rate] estimates per group. Conclusion on whether the original Gene[Rate] analysis was methodologically sound.

4. **Missing data impact** — missingness rate per locus/group, discussion of whether exclusion of missing alleles biases frequency estimates.

5. **Suggested improvements** (written as recommendations, not code):
   - Higher resolution (4-field) typing where available
   - Bootstrap confidence intervals for haplotype frequencies
   - Linkage disequilibrium (D') reporting between loci pairs
   - Consideration of admixed individuals in "Others" group

6. **Registry size findings** — key table from `registry_size_targets.csv` with narrative interpretation.

---

## Key Assumptions

| Assumption | Where used | Justification |
|------------|-----------|---------------|
| Hardy-Weinberg equilibrium | Diplotype freq calculation, Gene[Rate] EM | Standard in HLA frequency analysis; tested explicitly in 03 |
| 2-field resolution sufficient | All frequency calculations | Consistent with original paper; higher resolution data not uniformly available |
| Haplotype frequencies stable (donor ≈ patient population) | Registry model donor pool | Validated by Donor vs Patient comparison in original paper |
| Random mating within ethnicity | HWE assumption | Reasonable approximation for large outbred population |
| Registry donors drawn at random from population | Registry size model | Idealization; real registries have demographic biases |

---

## Out of Scope

- Re-running Gene[Rate] online pipeline (used as reference only)
- DPB1 or other HLA loci beyond A, B, C, DRB1, DQB1
- Patient survival or transplant outcome modeling
- Cost modeling of registry recruitment
- Individual patient matching simulation (Monte Carlo)
