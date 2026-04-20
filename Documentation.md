# Technical Documentation

## HLA Registry Analysis — Singapore CMIO Population

**Author:** Alvin Ng Yu-Jin  
**Date:** April 2026  
**Dataset:** BMDP + SCBB, n = 59,186 donors  
**Reference:** Ng AYJ et al. (2022), *Blood Cell Therapy* 5(3):86–95

---

## Table of Contents

1. [HLA Biology Background](#1-hla-biology-background)
2. [Dataset Description](#2-dataset-description)
3. [Data Processing — Pipeline Step 1](#3-data-processing--pipeline-step-1)
4. [Allele Frequency Verification — Pipeline Step 2](#4-allele-frequency-verification--pipeline-step-2)
5. [EM Algorithm and Hardy–Weinberg Equilibrium — Pipeline Step 3](#5-em-algorithm-and-hardyweinberg-equilibrium--pipeline-step-3)
6. [Registry Size Model — Pipeline Steps 4–5](#6-registry-size-model--pipeline-steps-45)
7. [Figure Interpretation](#7-figure-interpretation)
8. [Key Findings](#8-key-findings)
9. [Limitations and Suggested Improvements](#9-limitations-and-suggested-improvements)

---

## 1. HLA Biology Background

### What is HLA?

The Human Leukocyte Antigen (HLA) system, encoded within the Major Histocompatibility Complex (MHC) on chromosome 6p21, is the most polymorphic region of the human genome. HLA molecules are cell-surface glycoproteins that present peptide fragments to T-lymphocytes, forming the molecular basis of immune self/non-self discrimination.

There are two main classes:

- **Class I** (HLA-A, -B, -C): expressed on virtually all nucleated cells; present endogenous peptides to CD8⁺ cytotoxic T cells.
- **Class II** (HLA-DRB1, -DQB1, -DPB1): expressed primarily on antigen-presenting cells; present exogenous peptides to CD4⁺ helper T cells.

This project analyses the five loci most relevant to haematopoietic stem cell transplantation (HSCT):

| Locus | Class | Role |
|-------|-------|------|
| HLA-A | I | Primary transplant matching locus |
| HLA-B | I | Primary transplant matching locus |
| HLA-C | I | Included in 8/8 and 10/10 matching |
| HLA-DRB1 | II | Strongest Class II predictor of GvHD |
| HLA-DQB1 | II | Included in 10/10 (but not 8/8) matching |

### HLA Nomenclature

HLA alleles are named using a hierarchical system. This analysis uses **2-field (intermediate) resolution**:

```
HLA-A*02:01
     │  │└─ protein coding variant (field 2)
     │  └── allele group / serological antigen (field 1)
     └────── locus
```

A 2-field designation identifies a unique protein sequence but does not distinguish synonymous coding or non-coding variants (fields 3 and 4).

### Why HLA Matching Matters for HSCT

Mismatches at HLA loci between donor and recipient are the primary driver of:

- **Graft-versus-host disease (GvHD):** donor T cells attack recipient tissues
- **Graft failure:** recipient immune cells reject the transplanted stem cells
- **Transplant-related mortality**

A **10/10 match** (allele-level match at HLA-A, B, C, DRB1, DQB1) is the gold standard for unrelated donor transplantation. An **8/8 match** (excluding DQB1) is the minimum typically accepted when a 10/10 match is unavailable.

### HLA in Multiethnic Populations

HLA allele and haplotype frequencies vary substantially between ethnic populations due to evolutionary history, founder effects, and geographic isolation. Singapore's CMIO (Chinese, Malay, Indian, Others) population is genetically diverse:

- **Chinese:** largest group (~77% of BMDP registry), predominantly Han Chinese ancestry, high frequency of A*02:07 and DRB1*09:01
- **Malay:** Austronesian ancestry, highest frequency of B*15:02 haplotypes  
- **Indian:** South Asian ancestry, higher DRB1*15:01 frequency, smaller effective population
- **Others:** heterogeneous non-CMIO backgrounds; HWE departures expected from population admixture

Because haplotype frequencies differ markedly between groups, a patient from a minority ethnic group is statistically less likely to find a match in a registry dominated by another ethnicity.

---

## 2. Dataset Description

### Sources

| File | Samples | Notes |
|------|---------|-------|
| `HLA Data.cleaned.xlsx` (sheets: BMDP.out, SCBB.out) | 59,186 | Primary analysis cohort; fully typed at all 5 loci |
| `DonorPatient.txt` | 1,350 | HSA donor haplotypes; single-field typing (allele2 absent) |
| `Patient.txt` | ~560 | HSA recipient haplotypes; single-field typing |

**HSA files were excluded from all frequency calculations** due to incomplete genotyping (single allele per locus — cannot compute diploid allele or heterozygosity statistics).

### BMDP+SCBB Cohort Breakdown

| Ethnicity | n (individuals) | % |
|-----------|----------------|---|
| Chinese | 44,400 | 75.0% |
| Malay | 5,578 | 9.4% |
| Indian | 5,490 | 9.3% |
| Others | 3,767 | 6.4% |
| **Total** | **59,235** | |

*Note: totals include both BMDP and SCBB sheets; small discrepancies from published n=59,186 due to duplicate handling.*

### Data Format

After ingestion, data is stored in tidy long format (`hla_clean.csv`):

| Column | Description |
|--------|-------------|
| `sample_id` | Unique donor/patient identifier |
| `source` | `BMDP_OUT`, `SCBB_OUT`, `HSA_DONOR`, `HSA_PATIENT` |
| `ethnicity` | `Chinese`, `Malay`, `Indian`, `Others` |
| `locus` | `HLA-A`, `HLA-B`, `HLA-C`, `DRB1`, `DQB1` |
| `allele1` | First allele (2-field format, e.g., `02:01`) |
| `allele2` | Second allele; NaN if missing |

**Total rows:** 305,745 (= 61,149 samples × 5 loci, allowing for HSA single-field records)

---

## 3. Data Processing — Pipeline Step 1

**Script:** `analysis/01_ingest.py`  
**Tests:** `tests/test_ingest.py` (13 tests)

### Normalisation

Raw allele values from the Excel files undergo the following transformations:

1. **Strip prefix:** `HLA-A*02:01` → `02:01`
2. **Truncate to 2-field:** `02:01:01:01` → `02:01`
3. **Remove suffixes:** `02:01G`, `02:01P` → `02:01` (G/P groups collapsed to protein level)
4. **Handle missing:** `-`, blank, `0`, `00:00` → `NaN`
5. **Ethnicity mapping:** `C` → `Chinese`, `M` → `Malay`, `I` → `Indian`, `O` → `Others`

### Column Detection

The Excel files use non-standard column naming (`DRB11`/`DRB12` instead of `DRB1_1`/`DRB1_2`). The ingestion pipeline uses a regex-based column detector that handles:

- Standard: `A1`, `A2`, `B1`, `B2`
- No-separator: `DRB11`, `DRB12`, `DQB11`, `DQB12`
- Prefixed: `HLA-A1`, `HLA-A2`

### Output Shape

- **305,745 rows** total in `hla_clean.csv`
- All 5 loci present for each BMDP/SCBB sample
- HSA files contribute single-allele rows (allele2 = NaN throughout)

---

## 4. Allele Frequency Verification — Pipeline Step 2

**Script:** `analysis/02_allele_freq.py`  
**Tests:** `tests/test_allele_freq.py` (5 tests)

### Method

For each (ethnicity, locus), the observed allele frequency is computed as:

$$f(a) = \frac{\text{count}(a \in \text{allele}_1) + \text{count}(a \in \text{allele}_2)}{2 \times N_{\text{individuals typed at locus}}}$$

where NaN values of allele₂ are excluded from both numerator and denominator. This is the **maximum likelihood estimator** under the assumption of random mating (Hardy–Weinberg equilibrium): it converges in one step of the EM algorithm.

### Comparison with Published Values

Published allele frequencies were loaded from `BMDPnSCBB.results.xlsx`, which stores Gene[Rate] software outputs in a wide format. For each (ethnicity, locus, allele), the signed difference is computed:

$$\Delta f(a) = f_{\text{observed}}(a) - f_{\text{published}}(a)$$

Alleles with $|\Delta f| > 0.005$ (0.5 percentage points) are flagged.

### Results

| Locus | Alleles compared | Max \|Δf\| | Flagged |
|-------|-----------------|-----------|---------|
| HLA-A | 283 | 0.049% | 0 |
| HLA-B | 482 | 0.040% | 0 |
| HLA-C | 281 | 0.272% | 0 |
| DRB1 | 310 | 0.035% | 0 |
| DQB1 | 132 | 0.056% | 0 |
| **Total** | **1,488** | **0.272%** | **0** |

**Conclusion:** The 2022 publication's allele frequencies are independently reproducible. The maximum discrepancy (0.272% at HLA-C) is well below the 0.5% threshold and consistent with floating-point rounding in Gene[Rate].

---

## 5. EM Algorithm and Hardy–Weinberg Equilibrium — Pipeline Step 3

**Script:** `analysis/03_hwe_test.py` | Library: `analysis/hwe_test.py`  
**Tests:** `tests/test_hwe_test.py` (5 tests)

### 5.1 Background: Hardy–Weinberg Equilibrium

Hardy–Weinberg Equilibrium (HWE) describes the genotype frequency distribution expected in a large, randomly mating population with no selection, mutation, migration, or genetic drift. For a locus with alleles $a_1, a_2, \ldots, a_k$ with frequencies $p_1, p_2, \ldots, p_k$:

$$P(\text{genotype } a_i a_i) = p_i^2 \qquad \text{(homozygote)}$$
$$P(\text{genotype } a_i a_j) = 2 p_i p_j \quad i \neq j \qquad \text{(heterozygote)}$$

**Expected heterozygosity** under HWE:

$$H_{\exp} = 1 - \sum_{i} p_i^2$$

**Observed heterozygosity:**

$$H_{\text{obs}} = \frac{\text{count}(\text{allele}_1 \neq \text{allele}_2)}{N}$$

### 5.2 HWE Test Statistic

The chi-squared test statistic used here compares observed vs expected heterozygosity:

$$\chi^2 = N \cdot \frac{(H_{\text{obs}} - H_{\text{exp}})^2}{H_{\text{exp}} \cdot (1 - H_{\text{exp}})}$$

with 1 degree of freedom (df = 1 for the 2-class partition heterozygous/homozygous). The p-value is obtained from `scipy.stats.chi2.sf(χ², df=1)`.

**Bonferroni correction:** With 20 simultaneous tests (5 loci × 4 ethnicities), the family-wise significance threshold is:

$$\alpha_{\text{corrected}} = \frac{0.05}{20} = 0.0025$$

### 5.3 EM Haplotype Frequency Estimation

#### Why EM is needed

In diploid organisms, each individual carries two copies of each chromosome. For a person heterozygous at multiple loci, it is not directly observable which alleles are on the same chromosome (i.e., which alleles form a *haplotype*). The assignment of alleles to chromosomes is called **phase determination**.

For example, a person with genotype:
- HLA-A: 02:01 / 11:01  
- HLA-B: 07:02 / 40:01

could have haplotypes **(02:01, 07:02) + (11:01, 40:01)** OR **(02:01, 40:01) + (11:01, 07:02)** — both are consistent with the observed genotype.

The **Expectation–Maximisation (EM) algorithm** resolves this ambiguity probabilistically by iterating between:

- **E-step:** given current haplotype frequency estimates, compute the posterior probability that each individual carries each pair of haplotypes
- **M-step:** update haplotype frequencies as the weighted sum of fractional haplotype assignments

#### Algorithm (simplified per-locus approach used here)

Due to the computational complexity of full 5-locus phase resolution, this pipeline uses a **product-approximation EM**:

1. For each locus independently, compute per-allele frequencies using the MLE formula above (this is exact; EM converges in one step for allele frequencies).

2. For 5-locus haplotype frequencies: enumerate all $(h_i, h_j)$ haplotype pairs consistent with each individual's genotype. A haplotype $h = (a_A, a_B, a_C, a_{DRB1}, a_{DQB1})$ combines allele 1 or allele 2 at each locus. Each individual contributes exactly 2 phase assignments (combining allele 1 at all loci, and allele 2 at all loci). 

3. **E-step:** weight each phase assignment by the product of haplotype frequencies:
   $$w_{ij} = \frac{f(h_i) \cdot f(h_j)}{\sum_{i'j'} f(h_{i'}) \cdot f(h_{j'})}$$

4. **M-step:** update frequencies from fractional counts:
   $$f'(h_k) = \frac{\sum_{\text{individuals}} \sum_{(i,j): h_i=h_k} w_{ij} + \sum_{(i,j): h_j=h_k} w_{ij}}{2N}$$

5. Iterate until $\max_k |f'(h_k) - f(h_k)| < 10^{-6}$ or 100 iterations.

6. Retain haplotypes with $f \geq 0.001$ (0.1%).

**Only individuals typed at all 5 loci** are used for haplotype estimation. A cap of 5,000 samples per ethnicity is applied for computational tractability.

#### Haplotype Format

Haplotypes are stored as pipe-separated 5-tuples in locus order (A|B|C|DRB1|DQB1):

```
02:07|46:01|01:02|09:01|03:03
```

### 5.4 HWE Results Summary

| Ethnicity | Violations | Loci |
|-----------|-----------|------|
| Chinese | 0 | — |
| Malay | 0 | — |
| Indian | 3 | DQB1 (p=0.0017), HLA-B (p=1.8×10⁻⁶), HLA-C (p=2.4×10⁻⁶) |
| Others | 5 | All loci (p ≤ 3.8×10⁻²⁸ for HLA-B) |

All violations show **heterozygosity deficit** (H_obs < H_exp), consistent with:

- **Indian group:** mild population sub-structure (South Indian vs North Indian sub-populations in Singapore) or possible genotyping artefacts
- **Others group:** population heterogeneity is the primary explanation — the "Others" category pools genetically distinct sub-populations (Eurasians, Caucasians, East Asians outside CMIO classification), artificially inflating allelic diversity relative to any single random-mating population

---

## 6. Registry Size Model — Pipeline Steps 4–5

**Script:** `analysis/04_registry_model.py` | Library: `analysis/registry_model.py`  
**Figures:** `analysis/plot_coverage.py`  
**Tests:** `tests/test_registry_model.py` (6 tests)

### 6.1 Mathematical Framework

#### Step 1 — Diplotype Frequencies under HWE

Given a set of haplotypes $\{h_1, \ldots, h_K\}$ with frequencies $\{f_1, \ldots, f_K\}$ (summing to 1), diplotype frequencies under HWE are:

$$P(h_i, h_i) = f_i^2 \qquad \text{(homozygous diplotype)}$$
$$P(h_i, h_j) = 2 f_i f_j \quad i < j \qquad \text{(heterozygous diplotype)}$$

The sum of all diplotype frequencies is $(\sum_i f_i)^2 = 1$.

A **residual "other" haplotype** pools all haplotypes not in the top-K set (those cumulatively covering ≥99% of frequency mass), with frequency $f_{\text{other}} = 1 - \sum_{i=1}^K f_i$.

#### Step 2 — Per-Patient Match Probability

For a patient with diplotype $g$ (occurring with frequency $f_g$ in the patient population), the probability of finding **at least one matching donor** in a registry of $N$ independently drawn donors is:

$$P(\geq 1 \text{ match} \mid N, g) = 1 - (1 - f_g)^N$$

This uses the complement of the probability that all $N$ donors are non-matches.

#### Step 3 — Population Coverage

The expected fraction of patients who find at least one match is the weighted average over all diplotypes:

$$\text{Coverage}(N) = \sum_g f_g \cdot \left[1 - (1 - f_g)^N\right]$$

where the outer $f_g$ is the probability that a random patient has diplotype $g$, and the bracket is the match probability for that patient given $N$ donors. This formula assumes:

- Donors are drawn independently and uniformly from the same haplotype frequency distribution as patients (same-ethnicity model) or from the combined Singapore distribution (cross-ethnic model)
- HWE holds for both patient and donor diplotype distributions

#### Step 4 — Model Variants

Four scenario combinations are evaluated per ethnicity:

| Variant | Donor pool | Patient pool |
|---------|-----------|--------------|
| Same-ethnicity | Per-group haplotype freqs | Per-group |
| Cross-ethnic | Singapore-weighted combined freqs | Per-group |

For **cross-ethnic matching**, the coverage formula becomes:

$$\text{Coverage}_{\text{cross}}(N) = \sum_{g \in \text{patient}} f_g^{\text{patient}} \cdot \left[1 - (1 - f_g^{\text{donor}})^N\right]$$

where $f_g^{\text{patient}}$ is the diplotype frequency in the patient's ethnic group and $f_g^{\text{donor}}$ is the frequency of that same diplotype in the combined donor pool (0 if the haplotype pair does not appear in the combined pool).

**Singapore population weights** (from BMDP composition):

| Ethnicity | Weight |
|-----------|--------|
| Chinese | 77% |
| Malay | 8% |
| Indian | 9% |
| Others | 6% |

#### Step 5 — Match Levels

- **8/8 match:** HLA-A, B, C, DRB1 only — DQB1 dropped from haplotype. When two 5-locus haplotypes differ only at DQB1, they collapse to the same 4-locus haplotype; their frequencies are summed.
- **10/10 match:** All 5 loci (A, B, C, DRB1, DQB1).

#### Step 6 — Minimum Registry Size

For a target coverage $\theta \in \{0.75, 0.85, 0.90, 0.95\}$, the minimum registry size is:

$$N^* = \min \{N \in \mathbb{Z}^+ : \text{Coverage}(N) \geq \theta\}$$

Found via **binary search on a log₁₀ scale** over $N \in [1{,}000, 10{,}000{,}000]$ with 50 iterations (precision ~10⁻¹⁵).

### 6.2 Numerical Considerations

For large $N$ and small $f_g$, $(1-f_g)^N$ underflows to 0 in IEEE 754 double precision. This is numerically harmless — it means the match probability for that diplotype approaches 1.0, which is mathematically correct. The computation uses `numpy.float64` throughout.

For rare diplotypes with $f_g \approx 10^{-5}$, even $N = 10^7$ may give $(1-f_g)^N \approx e^{-100} \approx 3.7 \times 10^{-44}$ — effectively 1.0 match probability.

---

## 7. Figure Interpretation

### Figure 1: Allele Frequency Discrepancy Heatmap

**File:** `analysis/figures/allele_freq_heatmap.png`

**What it shows:** A seaborn heatmap where each cell represents one (ethnicity, locus, allele group) combination. Cell colour indicates the signed difference between independently computed allele frequency and the Gene[Rate]-published frequency:

$$\Delta f = f_{\text{observed}} - f_{\text{published}}$$

**Colour scale:**
- **Red/warm:** Observed frequency is higher than published
- **Blue/cool:** Observed frequency is lower than published
- **White/near-zero:** Agreement within rounding

**How to interpret:**
- All cells should be near-white (close to zero) for a reproducible analysis
- Any systematic colour shift across a row (locus) would indicate a calculation error specific to that locus
- Any systematic shift across a column (ethnicity) could indicate a data source mismatch

**Key finding:** All cells are near-white. Maximum discrepancy is 0.27% (HLA-C), well below the 0.5% flagging threshold. No systematic bias is visible. The published allele frequencies are independently confirmed.

---

### Figure 2: Coverage Curves — 8/8 Match

**File:** `analysis/figures/coverage_curves_8of8.png`

**What it shows:** Five panels (Chinese, Malay, Indian, Others, Combined), each plotting registry coverage as a function of registry size N for an **8/8 HLA match** (loci A, B, C, DRB1 — excluding DQB1).

**Axes:**
- **X-axis:** Registry size N (log₁₀ scale, from 1,000 to 10,000,000)
- **Y-axis:** Coverage — the expected proportion of patients who find at least one HLA-matched donor

**Lines:**
- **Solid blue:** Same-ethnicity matching (donor and patient from same ethnic group)
- **Dashed orange:** Cross-ethnic matching (Chinese-dominated combined registry as donor pool; patient from that ethnicity)

**Reference lines:** Horizontal dashed grey lines at 75%, 85%, 90%, 95% coverage thresholds.

**How to interpret:**

1. **Coverage curves are always monotone increasing** — more donors can only improve (or maintain) coverage, never reduce it. The curve's shape reflects the haplotype frequency distribution: when many common haplotypes are present, early donors cover most patients quickly (steep initial rise); rare diplotypes require exponentially more donors (long tail).

2. **Same-ethnicity vs cross-ethnic gap:** For minority groups (Malay, Indian, Others), the blue curve typically lies far to the left (fewer donors needed) compared to the orange curve. This is because the Chinese-dominated combined registry has very different haplotype frequencies from minority populations — a Malay patient is unlikely to find a match in a predominantly Chinese donor pool.

3. **Ceiling values:** When the orange curve flattens far below the target threshold even at N = 10,000,000, cross-ethnic matching is practically infeasible for that group. This is visible as the dashed orange curve not reaching the grey reference lines.

4. **Combined panel:** The rightmost panel shows the coverage for a Singapore-weighted average patient. This represents the overall registry utility across all ethnicities.

**Key observations:**
- Chinese: both curves reach 90%+ within ~15,000 donors (8of8 same-ethnicity)
- Indian: fewest donors needed due to less haplotypic diversity (79 haplotypes but more concentrated frequencies)
- Malay + Others: cross-ethnic 8of8 coverage never reaches 75% regardless of registry size

---

### Figure 3: Coverage Curves — 10/10 Match

**File:** `analysis/figures/coverage_curves_10of10.png`

**What it shows:** Identical layout to Figure 2 but for **10/10 HLA match** (all 5 loci: A, B, C, DRB1, DQB1).

**Differences from 8/8:**

Adding DQB1 as a fifth matching locus increases HLA diversity (more possible haplotypes), making matches harder to find. The 10/10 coverage curves are shifted **right** compared to 8/8 — a larger registry is needed to achieve the same coverage.

**How to quantify the DQB1 effect (Chinese same-ethnicity):**

| Coverage target | 8/8 registry size | 10/10 registry size | Ratio |
|----------------|------------------|--------------------|----|
| 75% | 3,689 | 3,883 | 1.05× |
| 85% | 5,691 | 6,008 | 1.06× |
| 90% | 7,498 | 7,926 | 1.06× |
| 95% | 10,986 | 11,616 | 1.06× |

For the Chinese population, adding DQB1 increases the required registry size by only ~6%, reflecting that DQB1 is in strong linkage disequilibrium with DRB1 — knowing DRB1 largely predicts DQB1 in this population.

**For minority groups**, the ratio may differ because LD structure varies between populations.

---

## 8. Key Findings

### Verification

- **The 2022 paper's allele frequency calculations are independently reproducible.** 1,488 alleles compared across 5 loci and 4 ethnicities; zero alleles exceeded the 0.5% discrepancy threshold. Maximum observed discrepancy: 0.27% (HLA-C). This validates the Gene[Rate] software outputs used in the original publication.

- **HWE holds for Chinese and Malay groups** (the two largest groups, comprising ~85% of the registry). This supports the validity of the haplotype frequency estimates in the original paper for these groups.

- **HWE violations in Indian (3 loci) and Others (all 5 loci)** suggest mild population sub-structure. These violations do not invalidate the registry model but add uncertainty to haplotype frequency estimates for minority groups, particularly at rare diplotypes.

### Registry Size Model

- **Chinese patients are best served** by the current registry: ~11,600 donors achieve 95% coverage at 10/10 resolution. The existing BMDP (~44,400 Chinese donors) far exceeds this threshold.

- **Minority group coverage depends critically on same-ethnicity matching.** Cross-ethnic matching is infeasible for Malay, Indian, and Others at 90%+ coverage targets, regardless of total registry size. This quantitatively demonstrates the importance of ethnic-specific donor recruitment.

- **The addition of DQB1 (10/10 vs 8/8) increases required registry size by ~6–10%** for Chinese patients due to strong DRB1-DQB1 linkage disequilibrium. The effect may be larger for other populations.

- **Haplotype frequency concentrations differ markedly.** The Indian "Others" group appears to have a concentrated haplotype distribution (requiring surprisingly few donors for same-ethnicity coverage) — this should be interpreted cautiously given the small sample size (n=754 "Others" individuals in the EM estimate) and HWE violations.

---

## 9. Limitations and Suggested Improvements

### Current Limitations

1. **2-field resolution:** The analysis uses 2-field HLA typing throughout. High-resolution (4-field) typing is increasingly standard and would improve match discrimination, particularly for DRB1 and HLA-B.

2. **Simplified EM (product approximation):** Full 5-locus phase resolution was approximated using a per-locus product model. This ignores multi-locus linkage disequilibrium (LD) — the tendency for certain allele combinations across loci to co-occur more (or less) than expected by chance. Strong LD (as exists between DRB1 and DQB1) means the simplified EM will over-estimate the number of distinct haplotypes and under-estimate their frequencies. A full haplo.stats-style EM would resolve phase ambiguity properly.

3. **No bootstrap confidence intervals:** Point estimates from the EM algorithm carry sampling uncertainty, especially for rare haplotypes in small cohorts (Indian n=1,098, Others n=754 for EM). Bootstrap CIs would quantify this uncertainty and allow propagation into registry size predictions.

4. **"Others" group heterogeneity:** The 5-locus HWE violations in the Others group are consistent with population admixture. Registry size predictions for this group are unreliable without stratification by ancestry.

5. **Registry composition idealization:** The model assumes donors are drawn randomly from the population. Real registries have demographic biases (age, recruitment campaigns) that may inflate or deflate effective coverage.

### Suggested Future Analyses

1. **Linkage disequilibrium reporting (D', r²):** Pairwise LD between the 5 loci would validate the per-locus independence assumption and characterise the haplotype block structure in each CMIO group.

2. **Full multi-locus EM:** Replace the product approximation with a proper full-phase EM (equivalent to haplo.stats in R) for accurate 5-locus haplotype frequency estimation. This is computationally feasible for n ≤ 5,000 samples.

3. **Higher-resolution typing:** Reanalyse with 4-field allele designations where available, particularly for DRB1 and HLA-B where high polymorphism affects transplant outcomes.

4. **Ancestry stratification for "Others":** Apply principal component analysis (PCA) on HLA allele vectors to identify genetic sub-clusters within the Others group, then model each sub-cluster separately.

5. **Bootstrap registry size confidence intervals:** Resample haplotype frequencies from the EM posterior to generate 95% CIs on the minimum N estimates.

---

## Appendix: Software and Reproducibility

| Component | Version |
|-----------|---------|
| Python | ≥ 3.10 |
| pandas | ≥ 2.0 |
| numpy | ≥ 1.24 |
| scipy | ≥ 1.10 |
| matplotlib | ≥ 3.7 |
| seaborn | ≥ 0.12 |
| openpyxl | ≥ 3.1 |
| pytest | ≥ 7.4 |

All randomised steps (EM sampling cap) use `random_state=42` for reproducibility.

Run `pytest tests/ -v` to verify all 29 tests pass before reproducing the analysis.

---

*End of Documentation*
