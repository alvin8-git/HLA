# HLA Registry Analysis — Singapore CMIO Population

[![Version](https://img.shields.io/badge/version-1.3.0-informational)](VERSION.md)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue?logo=python&logoColor=white)](https://www.python.org/)
[![Tests](https://img.shields.io/badge/tests-29%20passing-brightgreen?logo=pytest&logoColor=white)](tests/)
[![License](https://img.shields.io/badge/license-MIT-blue)](LICENSE)
[![Dataset](https://img.shields.io/badge/dataset-n%3D59%2C186%20donors-informational)](analysis/data/)
[![Paper](https://img.shields.io/badge/paper-Blood%20Cell%20Therapy%202022-red)](2022_HLA_BloodCellTherapy.pdf)
[![Loci](https://img.shields.io/badge/HLA%20loci-A%20%7C%20B%20%7C%20C%20%7C%20DRB1%20%7C%20DQB1-purple)](analysis/data/allele_freqs_per_locus.csv)
[![Status](https://img.shields.io/badge/pipeline-complete-success)](analysis/run_all.sh)

A reproducible Python pipeline for:

1. **Verification** — Independent re-analysis of HLA allele and haplotype frequencies from Singapore's BMDP + SCBB bone marrow donor registry, auditing the methodology of the [2022 *Blood Cell Therapy* publication](2022_HLA_BloodCellTherapy.pdf).
2. **Registry size modelling** — Predicting the minimum number of donors needed to achieve 75–95% patient match probability at 8/8 and 10/10 HLA match levels, stratified by CMIO ethnicity and matching model.

---

## Table of Contents

- [Background](#background)
- [Dataset](#dataset)
- [Pipeline Overview](#pipeline-overview)
- [Quick Start](#quick-start)
- [Results Summary](#results-summary)
- [Output Files](#output-files)
- [Figures](#figures)
- [Project Structure](#project-structure)
- [Dependencies](#dependencies)
- [Citation](#citation)

---

## Background

Human Leukocyte Antigen (HLA) matching between donor and recipient is critical for the success of haematopoietic stem cell transplantation (HSCT). The probability that a patient finds a suitably matched donor depends on:

- The **allele frequency distribution** in the donor registry
- The **haplotype structure** (linkage disequilibrium between HLA loci)
- The **size** of the donor registry
- Whether matching is within the same ethnic group (**same-ethnicity**) or across groups (**cross-ethnic**)

Singapore's multiethnic CMIO (Chinese, Malay, Indian, Others) population presents a unique challenge: minority groups have smaller effective donor pools, making ethnic-specific registry targets essential for equitable transplant access.

---

## Dataset

| Source | Description | n (samples) |
|--------|-------------|-------------|
| **BMDP** | Bone Marrow Donor Programme (Singapore) | ~44,400 Chinese, ~5,578 Malay, ~5,490 Indian, ~3,767 Others |
| **SCBB** | Singapore Cord Blood Bank | Included in above totals |
| **HSA-Donor** | Health Sciences Authority donor haplotypes | 1,350 |
| **HSA-Patient** | Health Sciences Authority recipient haplotypes | ~560 |

**Total BMDP+SCBB cohort: 59,186 donors**

HLA loci typed (2-field resolution): **A, B, C, DRB1, DQB1**

Published reference frequencies: `BMDPnSCBB.results.xlsx` (Gene[Rate] software outputs)

---

## Pipeline Overview

```
HLA Data.cleaned.xlsx (BMDP + SCBB)
    │
    ▼ 01_ingest.py
    analysis/data/hla_clean.csv (305,745 rows)
    │
    ├──▶ 02_allele_freq.py ──▶ allele_freq_comparison.csv
    │                          allele_freq_heatmap.png
    │
    ├──▶ 03_hwe_test.py ───▶ haplo_freqs_em.csv
    │                         allele_freqs_per_locus.csv
    │                         hwe_results.csv
    │
    ├──▶ 04_registry_model.py + plot_coverage.py
    │         ──▶ coverage_curves.csv
    │             registry_size_targets.csv
    │             coverage_curves_8of8.png
    │             coverage_curves_10of10.png
    │
    └──▶ 05_report.py ──▶ verification_summary.md
```

---

## Quick Start

### 1. Install dependencies

```bash
pip install -r analysis/requirements.txt
```

### 2. Run the full pipeline

```bash
bash analysis/run_all.sh
```

> **Note:** Step 1 (data ingestion) reads the Excel source files (~59K rows). If `hla_clean.csv` already exists, step 1 is skipped automatically.

### 3. Run tests

```bash
pytest tests/ -v
```

Expected: **29 tests passing**

### 4. Run individual steps

```bash
cd analysis/
python 02_allele_freq.py      # allele frequency verification
python 03_hwe_test.py         # EM haplotypes + HWE tests
python 04_registry_model.py   # registry size model
python plot_coverage.py       # generate coverage curve figures
python 05_report.py           # assemble verification_summary.md
```

---

## Results Summary

### Allele Frequency Verification

| Metric | Result |
|--------|--------|
| Total alleles compared | 1,488 |
| Flagged alleles (\|diff\| > 0.5%) | **0** |
| Max discrepancy (HLA-C) | 0.27% |
| Conclusion | ✅ Fully reproducible |

### Hardy–Weinberg Equilibrium

| Group | Violations (Bonferroni p < 0.0025) |
|-------|-------------------------------------|
| Chinese | 0/5 loci |
| Malay | 0/5 loci |
| Indian | 3/5 loci (DQB1, HLA-B, HLA-C — heterozygosity deficit) |
| Others | 5/5 loci (expected; heterogeneous population) |

### Registry Size Targets (same-ethnicity, 10/10 match)

| Ethnicity | 75% coverage | 85% coverage | 90% coverage | 95% coverage |
|-----------|-------------|-------------|-------------|-------------|
| Chinese | 3,883 | 6,008 | 7,926 | 11,616 |
| Malay | 5,840 | 9,096 | 12,017 | 17,601 |
| Indian | 1,459 | 2,121 | 2,692 | 3,759 |
| Others | 1,001 | 1,001 | 1,008 | 1,430 |
| **Combined** | **14,041** | **26,749** | **42,212** | **83,541** |

> Cross-ethnic matching is infeasible for Malay, Indian, and Others at high coverage targets (≥90%) — model reaches the 10M ceiling.

See [Documentation.md](Documentation.md) for full methodology and figure interpretation.

---

## Output Files

| File | Description |
|------|-------------|
| `analysis/data/hla_clean.csv` | Tidy HLA data (305,745 rows; sample_id, source, ethnicity, locus, allele1, allele2) |
| `analysis/data/allele_freq_comparison.csv` | Observed vs published allele frequencies with difference and flag columns |
| `analysis/data/allele_freqs_observed.csv` | Observed allele frequencies from BMDP+SCBB |
| `analysis/data/allele_freqs_per_locus.csv` | Per-locus allele frequencies from EM (all ethnicities) |
| `analysis/data/haplo_freqs_em.csv` | 5-locus haplotype frequencies from EM (freq ≥ 0.1%) |
| `analysis/data/hwe_results.csv` | HWE chi-squared test results (20 tests: 5 loci × 4 ethnicities) |
| `analysis/data/coverage_curves.csv` | Coverage(N) for N = 1,000–10,000,000 across all scenarios |
| `analysis/data/registry_size_targets.csv` | Minimum registry N per (match level × ethnicity × variant × threshold) |
| `analysis/verification_summary.md` | Full narrative verification report |

---

## Figures

| Figure | Description |
|--------|-------------|
| `analysis/figures/allele_freq_heatmap.png` | Heatmap of allele frequency discrepancies (observed − published) across all loci and CMIO groups |
| `analysis/figures/coverage_curves_8of8.png` | Registry coverage curves for 8/8 HLA match (A, B, C, DRB1) — exact match only |
| `analysis/figures/coverage_curves_10of10.png` | Registry coverage curves for 10/10 HLA match (A, B, C, DRB1, DQB1) — exact match only |
| `analysis/figures/partial_match_10locus.png` | Partial match coverage curves: 8/10, 9/10, 10/10 for each CMIO group + Overall (WBMT Fig 5 style) |
| `analysis/figures/partial_match_8locus.png` | Partial match coverage curves: 6/8, 7/8, 8/8 for each CMIO group + Overall |

See [Documentation.md](Documentation.md) for detailed figure interpretation.

---

## Project Structure

```
HLA/
├── analysis/
│   ├── 01_ingest.py           # Data ingestion & normalisation
│   ├── 02_allele_freq.py      # Allele frequency computation & comparison
│   ├── 03_hwe_test.py         # EM haplotype estimation + HWE tests
│   ├── 04_registry_model.py   # Registry size coverage model (runner)
│   ├── 05_report.py           # Summary report assembly
│   ├── hwe_test.py            # HWE library module
│   ├── registry_model.py      # Registry model library module
│   ├── plot_coverage.py       # Coverage curve figure generator
│   ├── run_all.sh             # Sequential pipeline driver
│   ├── requirements.txt       # Python dependencies
│   ├── data/                  # Intermediate & final CSV outputs
│   └── figures/               # Output PNG figures
├── tests/
│   ├── test_ingest.py         # 13 tests
│   ├── test_allele_freq.py    # 5 tests
│   ├── test_hwe_test.py       # 5 tests
│   └── test_registry_model.py # 6 tests
├── docs/
│   └── superpowers/           # Design spec & implementation plan
├── README.md
├── Documentation.md
└── VERSION.md
```

---

## Dependencies

```
pandas>=2.0
numpy>=1.24
scipy>=1.10
matplotlib>=3.7
seaborn>=0.12
openpyxl>=3.1
pytest>=7.4
```

Install with: `pip install -r analysis/requirements.txt`

---

## Citation

If you use this analysis or its outputs, please cite the original paper:

> Ng AYJ et al. (2022). HLA allele and haplotype frequencies of the Singapore bone marrow donor registry and cord blood bank. *Blood Cell Therapy*, 5(3), 86–95.

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

*Analysis by Alvin Ng Yu-Jin · Singapore, 2026*
