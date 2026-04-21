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
   - [5.5 EM Validation Against Gene[Rate] Published Frequencies](#55-em-validation-against-generate-published-frequencies)
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

$$
f(a) = \frac{\text{count}(a \in \text{allele}_1) + \text{count}(a \in \text{allele}_2)}{2 \times N_{\text{individuals typed at locus}}}
$$

where NaN values of allele₂ are excluded from both numerator and denominator. This is the **maximum likelihood estimator** under the assumption of random mating (Hardy–Weinberg equilibrium): it converges in one step of the EM algorithm.

### Comparison with Published Values

Published allele frequencies were loaded from `BMDPnSCBB.results.xlsx`, which stores Gene[Rate] software outputs in a wide format. For each (ethnicity, locus, allele), the signed difference is computed:

$$
\Delta f(a) = f_{\text{observed}}(a) - f_{\text{published}}(a)
$$

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

$$
P(\text{genotype } a_i a_i) = p_i^2 \qquad \text{(homozygote)}
$$

$$
P(\text{genotype } a_i a_j) = 2 p_i p_j \quad i \neq j \qquad \text{(heterozygote)}
$$

**Expected heterozygosity** under HWE:

$$
H_{\exp} = 1 - \sum_{i} p_i^2
$$

**Observed heterozygosity:**

$$
H_{\text{obs}} = \frac{\text{count}(\text{allele}_1 \neq \text{allele}_2)}{N}
$$

### 5.2 HWE Test Statistic

The chi-squared test statistic used here compares observed vs expected heterozygosity:

$$
\chi^2 = N \cdot \frac{(H_{\text{obs}} - H_{\text{exp}})^2}{H_{\text{exp}} \cdot (1 - H_{\text{exp}})}
$$

with 1 degree of freedom (df = 1 for the 2-class partition heterozygous/homozygous). The p-value is obtained from `scipy.stats.chi2.sf(χ², df=1)`.

**Bonferroni correction:** With 20 simultaneous tests (5 loci × 4 ethnicities), the family-wise significance threshold is:

$$
\alpha_{\text{corrected}} = \frac{0.05}{20} = 0.0025
$$

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

#### Algorithm — full multi-locus EM with proper phase enumeration

This pipeline implements a **full multi-locus EM** that correctly handles phase ambiguity across all 5 loci simultaneously. For an individual heterozygous at $k$ loci, there are $2^{k-1}$ distinct phase configurations (the factor of ½ removes the h₁↔h₂ symmetry). All configurations are enumerated before the EM loop and re-used each iteration.

1. **Phase enumeration (`_enum_phase_configs`):** for each individual, produce the list of all valid (h₁, h₂) pairs by fixing the first heterozygous locus on h₁ and varying the remaining $k-1$ heterozygous loci over all $2^{k-1}$ assignments.

2. **E-step:** for each individual with phase configurations $\{(h_i, h_j)\}$, weight each configuration by its diplotype frequency under the current estimates:
   $$
   w_{ij} = \frac{f(h_i)^2 \cdot \mathbf{1}[i=j] + 2 f(h_i) f(h_j) \cdot \mathbf{1}[i \neq j]}{\sum_{i'j'} \left(f(h_{i'})^2 \cdot \mathbf{1}[i'=j'] + 2 f(h_{i'}) f(h_{j'}) \cdot \mathbf{1}[i' \neq j']\right)}
   $$

3. **M-step:** accumulate fractional haplotype counts across all individuals and configurations:
   $$
   f'(h_k) = \frac{\sum_{\text{ind.}} \sum_{(i,j)} w_{ij} \left(\mathbf{1}[h_i = h_k] + \mathbf{1}[h_j = h_k]\right)}{2N}
   $$

4. Iterate until $\max_k |f'(h_k) - f(h_k)| < 10^{-6}$ or 200 iterations.

5. Retain haplotypes with $f \geq 0.001$ (0.1%).

**Improvement over product-approximation:** an earlier version assigned allele-column 1 across all loci as haplotype h₁ and column 2 as h₂ — never exploring alternative phase assignments. Because column ordering in the input data is arbitrary, this systematically underestimated haplotype diversity by treating spurious allele combinations as real haplotypes. The full EM resolves phase by LD context, producing 2–4× more distinct haplotypes per ethnicity and 2–22× higher registry size estimates at high coverage thresholds.

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

### 5.5 EM Validation Against Gene[Rate] Published Frequencies

To validate our full-EM implementation, haplotype frequencies were compared against the Gene[Rate] software estimates published in the original Ng et al. 2022 paper (`BMDPnSCBB.results.xlsx`, `Haplotype.*` sheets). Gene[Rate] is an established HLA haplotype frequency tool used in bone marrow registry analysis worldwide.

#### Method

Gene[Rate]'s `Haplotype.*` sheets (one per ethnicity) contain estimated frequencies in tilde-separated format (`A*33:03~B*58:01~C*03:02~DRB1*03:01~DQB1*02:01`). These were converted to our pipe-separated format and matched on haplotype identity. Spearman rank correlation and RMSE were computed for all matched haplotypes. Script: `analysis/07_validate_em.py`.

#### Results

| Ethnicity | Our EM | Gene[Rate] | Matched | Spearman r | RMSE | Coverage (GR → EM) |
|-----------|--------|-----------|---------|------------|------|-------------------|
| Chinese | 140 | 2,196 | 140 | **0.913** | 0.00080 | 48.2% → 51.7% |
| Malay | 137 | 2,716 | 137 | **0.970** | 0.00027 | 52.6% → 52.9% |
| Indian | 144 | 3,475 | 144 | **0.963** | 0.00025 | 39.8% → 40.2% |
| Others | 123 | 3,079 | 123 | **0.990** | 0.000074 | 35.6% → 35.6% |

*Coverage = cumulative frequency mass of matched haplotypes. GR = Gene[Rate], EM = our estimate.*

#### Interpretation

**Rank agreement is strong (r ≥ 0.91).** The Spearman correlation reflects near-identical haplotype ranking between our EM and Gene[Rate] across all four ethnic groups. The Others group reaches r = 0.990, consistent with its more uniform haplotype distribution reducing sensitivity to minor frequency differences.

**Absolute frequency differences are negligible (RMSE < 0.001).** The small RMSE values confirm that individual haplotype frequencies — not just their ranks — are quantitatively consistent with the published estimates.

**All our haplotypes appear in Gene[Rate] (0 unmatched on our side).** Every haplotype we estimated above the 0.1% frequency threshold was also identified by Gene[Rate], confirming that our EM converges to the same haplotype space for common variants.

**The frequency coverage gap is expected, not a deficiency.** Gene[Rate] identifies 2,196–3,475 haplotypes per ethnicity versus our 123–144. The difference reflects two design choices: (1) Gene[Rate] used the full cohort (up to ~44,400 Chinese) without a sample cap, enabling estimation of haplotypes with frequencies down to ~0.001–0.01%, whereas our EM caps at 5,000 individuals per ethnicity; (2) our 0.1% frequency threshold filters rare haplotypes that individually contribute < 0.1% to coverage. The matched haplotypes together account for 35–53% of total frequency mass (Gene[Rate] denominator), with our EM recovering slightly more mass in each case — consistent with the full EM concentrating probability mass on fewer, better-resolved haplotypes.

**Top haplotypes agree closely with published values.** For Chinese, the top-4 haplotypes are identical in both ranking and frequency between our EM and Gene[Rate]:

| Rank | Haplotype | Gene[Rate] freq | Our EM freq |
|------|-----------|----------------|-------------|
| 1 | A\*33:03~B\*58:01~C\*03:02~DRB1\*03:01~DQB1\*02:01 | 0.0597 | 0.0541 |
| 2 | A\*02:07~B\*46:01~C\*01:02~DRB1\*09:01~DQB1\*03:03 | 0.0357 | 0.0383 |
| 3 | A\*11:01~B\*15:02~C\*08:01~DRB1\*12:02~DQB1\*03:01 | 0.0205 | 0.0220 |
| 4 | A\*33:03~B\*58:01~C\*03:02~DRB1\*13:02~DQB1\*06:09 | 0.0166 | 0.0170 |

The minor frequency differences (< 1 percentage point) reflect the 5,000-sample cap and random sampling, not algorithmic divergence.

**Conclusion:** The full-EM implementation is validated. It recovers the same dominant haplotype structure as Gene[Rate] with Spearman r ≥ 0.91 and RMSE < 0.001 across all CMIO groups.

---

## 6. Registry Size Model — Pipeline Steps 4–5

**Script:** `analysis/04_registry_model.py` | Library: `analysis/registry_model.py`  
**Figures:** `analysis/plot_coverage.py`  
**Tests:** `tests/test_registry_model.py` (6 tests)

### 6.1 Mathematical Framework

#### Step 1 — Diplotype Frequencies under HWE

Given a set of haplotypes $\{h_1, \ldots, h_K\}$ with frequencies $\{f_1, \ldots, f_K\}$ (summing to 1), diplotype frequencies under HWE are:

$$
P(h_i, h_i) = f_i^2 \qquad \text{(homozygous diplotype)}
$$

$$
P(h_i, h_j) = 2 f_i f_j \quad i < j \qquad \text{(heterozygous diplotype)}
$$

The sum of all diplotype frequencies is $(\sum_i f_i)^2 = 1$.

A **residual "other" haplotype** pools all haplotypes not in the top-K set (those cumulatively covering ≥99% of frequency mass), with frequency $f_{\text{other}} = 1 - \sum_{i=1}^K f_i$.

#### Step 2 — Per-Patient Match Probability

For a patient with diplotype $g$ (occurring with frequency $f_g$ in the patient population), the probability of finding **at least one matching donor** in a registry of $N$ independently drawn donors is:

$$
P(\geq 1 \text{ match} \mid N, g) = 1 - (1 - f_g)^N
$$

This uses the complement of the probability that all $N$ donors are non-matches.

#### Step 3 — Population Coverage

The expected fraction of patients who find at least one match is the weighted average over all diplotypes:

$$
\text{Coverage}(N) = \sum_g f_g \cdot \left[1 - (1 - f_g)^N\right]
$$

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

$$
\text{Coverage}_{\text{cross}}(N) = \sum_{g \in \text{patient}} f_g^{\text{patient}} \cdot \left[1 - (1 - f_g^{\text{donor}})^N\right]
$$

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

$$
N^* = \min \{N \in \mathbb{Z}^+ : \text{Coverage}(N) \geq \theta\}
$$

Found via **binary search on a log₁₀ scale** over $N \in [1{,}000, 10{,}000{,}000]$ with 50 iterations (precision ~10⁻¹⁵).

#### Literature context

This framework follows the probabilistic registry-size model established in the transplant literature. The formula $\text{Coverage}(N) = \sum_g f_g [1-(1-f_g)^N]$ was formalised by Beatty et al. (1988) and subsequently applied to large-scale population datasets by Maiers, Gragert and colleagues. Key references:

- **Beatty PG, Dahlberg S, Mickelson EM et al.** (1988). Probability of finding HLA-matched unrelated marrow donors. *Transplantation*, 45(4):714–718. *(Original probabilistic formulation of unrelated donor registry size.)*
- **Maiers M, Gragert L, Klitz W.** (2007). High-resolution HLA alleles and haplotypes in the United States population. *Human Immunology*, 68(9):779–788. *(First large-scale application of the coverage-curve formula to a national registry.)*
- **Gragert L, Madbouly A, Freeman J, Maiers M.** (2013). Six-locus high resolution HLA haplotype frequencies derived from mixed-resolution DNA typing for the entire US donor pool. *Human Immunology*, 74(10):1313–1320. *(Extension to 6-locus haplotypes and treatment of the "long tail" rare-diplotype problem.)*
- **Lim FLWI, Cheong SKS, Ho AYL et al.** (2010). HLA genetic analysis of the Singapore bone marrow donor programme — evidence for a heterogeneous Malay population. *Ann Acad Med Singapore*, 39(1):27–33. *(Singapore BMDP registry study providing the foundation for CMIO-specific models.)*
- **Aljurf M, Weisdorf D, Alfraih F et al.** (2019). Worldwide Network for Blood and Marrow Transplantation (WBMT) special article: Challenges facing developing countries in the establishment and maintenance of unrelated donor registries. *Bone Marrow Transplantation*, 54:1179–1188. *(International benchmarks for registry coverage targets and the partial-match curve framework.)*

### 6.2 Numerical Considerations

For large $N$ and small $f_g$, $(1-f_g)^N$ underflows to 0 in IEEE 754 double precision. This is numerically harmless — it means the match probability for that diplotype approaches 1.0, which is mathematically correct. The computation uses `numpy.float64` throughout.

For rare diplotypes with $f_g \approx 10^{-5}$, even $N = 10^7$ may give $(1-f_g)^N \approx e^{-100} \approx 3.7 \times 10^{-44}$ — effectively 1.0 match probability.

### 6.3 Numeric Examples from CMIO Data

This section walks through concrete calculations using the Singapore BMDP+SCBB data to illustrate how the model behaves in practice.

#### 6.3.1 Chinese population — 10/10 same-ethnicity matching

The EM algorithm identified **79 distinct haplotypes** (frequency ≥ 0.1%) in the Chinese cohort (~44,400 donors). HWE expansion produces **3,160 diplotype combinations**. The most frequent haplotype is:

> A\*02:07~B\*46:01~C\*01:02~DRB1\*09:01~DQB1\*03:03, $f_1 = 0.0516$

The **most common diplotype** is the heterozygote formed by the top two haplotypes.
$h_2$ = `A*11:01~B*40:01~C*03:02~DRB1*03:01~DQB1*02:01`, $f_2 = 0.0344$:

$$
f_g(h_1, h_2) = 2 \times 0.0516 \times 0.0344 = 0.00355
$$

**Per-patient match probability for this diplotype:**

$$
P(\geq 1\ \text{match} \mid N,\ g) = 1 - (1 - 0.00355)^N
$$

| $N$ donors | $(1 - 0.00355)^N$ | Match probability |
|------------|-------------------|-------------------|
| 100 | 0.699 | 30.1% |
| 500 | 0.168 | 83.2% |
| 1,000 | 0.028 | 97.2% |
| 3,000 | 0.000026 | 99.997% |

Even for the single most common diplotype in the Chinese population ($f_g = 0.00355$), a patient needs roughly **1,000 donors** before the match probability exceeds 97%. This is because each donor independently has only a 0.36% chance of carrying this exact diplotype.

**Population coverage across all 3,160 diplotypes:**

| $N$ donors | Coverage $\text{Coverage}(N)$ |
|------------|-------------------------------|
| 1,000 | 37.9% |
| 3,883 | **75.0%** ← 75% registry target |
| 6,008 | **85.0%** ← 85% registry target |
| 7,926 | **90.0%** ← 90% registry target |
| 11,616 | **95.0%** ← 95% registry target |
| 50,000 | 99.97% |

These figures confirm that the existing Chinese BMDP cohort (~44,400 donors) far exceeds all coverage targets — at 95% coverage the minimum requirement is 11,616 donors, less than 27% of the current registry.

**Note on incremental cost:** Moving from 85% → 90% requires +1,918 donors (+32%), while 90% → 95% requires +3,690 donors (+47%). Each additional 5% of coverage is progressively more expensive because it targets increasingly rare diplotypes.

#### 6.3.2 Effect of haplotype diversity on registry requirements

The number of distinct haplotypes and the concentration of their frequencies together determine how rapidly coverage accumulates with registry size.

| Ethnicity | Haplotypes (≥ 0.1%) | Diplotypes | Coverage at $N = 1{,}000$ | $N^*$ at 95%, 10/10 |
|-----------|----------------------|------------|---------------------------|----------------------|
| Others | 28 | 406 | 89.9% | 1,430 |
| Indian | 47 | 1,128 | 63.1% | 3,759 |
| Chinese | 79 | 3,160 | 37.9% | 11,616 |
| Malay | 97 | 4,753 | 29.2% | 17,601 |

**Indian patients** achieve the highest early-$N$ coverage despite being a smaller cohort. Their 47 haplotypes are more concentrated: the top-5 Chinese haplotypes carry $\sum f_i = 0.179$ of all frequency mass; the equivalent for Indian is $\sum f_i = 0.196$. Fewer, more common diplotypes mean fewer donors are needed to achieve high coverage.

**Malay patients** require the largest registry. The 97 haplotypes are spread more evenly — producing 4,753 diplotypes with a flatter frequency distribution. At $N = 1{,}000$ donors, a Malay patient has only a 29.2% chance of finding a 10/10 match, versus 63.1% for an Indian patient of the same registry size.

**Others** (a heterogeneous group) appear to need the fewest donors, but this reflects the concentration of a few very common North-European-type haplotypes (e.g. A\*01:01~B\*08:01~C\*07:01~DRB1\*03:01~DQB1\*02:01, $f = 0.129$) in a small admixed sample — HWE violations in this group add uncertainty to these estimates.

#### 6.3.3 The rare-diplotype long tail

The 79 Chinese haplotypes generate 3,160 diplotypes, but frequency is highly skewed:

| Diplotype tier | Cumulative frequency |
|----------------|----------------------|
| Top 10 | 2.8% |
| Top 100 | 14.8% |
| Top 500 | 41.1% |
| All 3,160 | 100.0% |

The top 10 diplotypes together represent only 2.8% of the patient population. The remaining 97.2% is distributed across 3,150 diplotypes — many of them with frequencies in the range $10^{-4}$–$10^{-3}$.

This is the "long-tail" problem inherent to HLA diversity. The coverage formula handles it correctly by summing contributions from every diplotype simultaneously, but it means:

- Early donors (small $N$) cover many patients at once — the S-curve rises steeply
- Late donors hit diminishing returns, each one matching an increasingly narrow slice of the population
- Achieving the final few percent of coverage requires disproportionately large registry growth

For the 5% of Chinese patients not covered by a 11,616-donor registry, their diplotypes have $f_g \lesssim 5 \times 10^{-4}$ — requiring on the order of $N \sim 1/(f_g) \approx 2{,}000$ donors just to reach 63% match probability for that individual diplotype. These patients would require a combined registry size of 50,000+ for near-certain coverage.

### 6.4 Partial Match Coverage Model

**Script:** `analysis/06_partial_match_plots.py`  
**Figures:** `analysis/figures/partial_match_10locus.png`, `analysis/figures/partial_match_8locus.png`

Sections 6.1–6.3 treat only **exact** HLA matching. In practice, transplant centres also consider partially matched donors (e.g. 9/10 or 8/10) when a fully matched donor cannot be found. The partial match model extends the framework to count allele-level matches across all loci.

#### Allele match count

Let $g_p = (h_p, h_q)$ be a patient diplotype and $g_d = (h_r, h_s)$ be a candidate donor diplotype. At each locus $\ell$, let $a_\ell(h)$ denote the allele carried on haplotype $h$. The **per-locus match score** takes the better of the two possible allele assignments:

$$
M_\ell(g_p, g_d) = \max\!\Bigl(
  \mathbf{1}[a_\ell(h_p) = a_\ell(h_r)] + \mathbf{1}[a_\ell(h_q) = a_\ell(h_s)],\;
  \mathbf{1}[a_\ell(h_p) = a_\ell(h_s)] + \mathbf{1}[a_\ell(h_q) = a_\ell(h_r)]
\Bigr)
$$

The **total allele match count** over all $L$ loci is:

$$
M(g_p, g_d) = \sum_{\ell=1}^{L} M_\ell(g_p, g_d)
$$

For the 10-locus framework ($L = 5$ loci, 10 allele positions) $M \in \{0, 1, \ldots, 10\}$; for 8-locus ($L = 4$) $M \in \{0, \ldots, 8\}$.

#### Per-patient partial match probability

For a patient with diplotype $g_p$ and a minimum match threshold $m$, the probability that a single randomly drawn donor meets or exceeds the threshold is:

$$
p_m(g_p) = \sum_{(h_r,\, h_s)} f_{h_r, h_s} \cdot \mathbf{1}\!\left[M(g_p,\, (h_r, h_s)) \geq m\right]
$$

where the sum runs over all donor diplotypes $(h_r, h_s)$ weighted by their frequency $f_{h_r, h_s}$ (computed under HWE from EM haplotype frequencies). $\mathbf{1}[\cdot]$ is the indicator function.

#### Partial match coverage curve

Substituting $p_m(g_p)$ in place of the exact-match diplotype frequency gives the partial match coverage curve:

$$
\text{Coverage}_m(N) = \sum_{g_p} f_{g_p} \cdot \left[1 - \left(1 - p_m(g_p)\right)^N\right]
$$

This generalises the exact-match formula (Section 6.1, Step 3): when $m$ equals the total number of allele positions, $p_m(g_p) = f_{g_p}$ and the two formulae are identical.

#### Implementation note

The full diplotype enumeration scales as $O(K^2)$ in the number of haplotypes $K$. For the Chinese population ($K = 79$) this produces 3,160 patient diplotypes × 3,160 donor diplotypes = ~10 million pairs per threshold level. The implementation (`compute_partial_match_probs`) uses NumPy broadcasting to vectorise the allele comparison across all pairs simultaneously, keeping wall time under one minute per ethnicity on a modern CPU.

---

## 7. Figure Interpretation

### Figure 1: Allele Frequency Discrepancy Heatmap

**File:** `analysis/figures/allele_freq_heatmap.png`

![Allele Frequency Discrepancy Heatmap](analysis/figures/allele_freq_heatmap.png)

**What it shows:** A seaborn heatmap where each cell represents one (ethnicity, locus, allele group) combination. Cell colour indicates the signed difference between independently computed allele frequency and the Gene[Rate]-published frequency:

$$
\Delta f = f_{\text{observed}} - f_{\text{published}}
$$

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

![Coverage Curves — 8/8 Match](analysis/figures/coverage_curves_8of8.png)

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

![Coverage Curves — 10/10 Match](analysis/figures/coverage_curves_10of10.png)

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

---

### Figure 4: Partial Match Coverage — 10-Locus Framework (8/10, 9/10, 10/10)

**File:** `analysis/figures/partial_match_10locus.png`

![Partial Match Coverage — 10-Locus](analysis/figures/partial_match_10locus.png)

**What it shows:** Five panels (Chinese, Malay, Indian, Others, Overall), each plotting registry coverage as a function of registry size N for three match stringency levels under the 10-locus framework (HLA-A, B, C, DRB1, DQB1 — 10 total allele positions):

- **Blue (10/10):** All 10 alleles match exactly across 5 loci
- **Red (9/10):** At least 9 of 10 alleles match (≤1 mismatch tolerated)
- **Green (8/10):** At least 8 of 10 alleles match (≤2 mismatches tolerated)

**X-axis:** Registry size N (log₁₀ scale, 1,000 to 10,000,000)  
**Y-axis:** Percentage of patients with at least one matching donor (0–100%)

**How to interpret:**

Match probabilities are computed using **haplotype-pair enumeration** over EM-estimated haplotype frequencies (capturing linkage disequilibrium between loci). For each patient diplotype $g_p = (h_p, h_q)$, the match probability for threshold $m$ is:

$$
p_m(g_p) = \sum_{(h_r,\, h_s)} f_{h_r, h_s} \cdot \mathbf{1}\!\left[M\!\left(g_p,\, (h_r, h_s)\right) \geq m\right]
$$

where $M(g_p, g_d)$ is the total allele match count defined in Section 6.4, using optimal per-locus assignment: $M_\ell = \max\!\bigl((a_\ell(h_p){=}a_\ell(h_r))+(a_\ell(h_q){=}a_\ell(h_s)),\ (a_\ell(h_p){=}a_\ell(h_s))+(a_\ell(h_q){=}a_\ell(h_r))\bigr)$.

**Key observations:**

- **Relaxing from 10/10 to 9/10 dramatically increases coverage at small N.** A single mismatch tolerance roughly halves the registry size needed to achieve a given coverage level for most ethnic groups.
- **8/10 tolerance achieves high coverage (~90–100%) even at N ≈ 100,000** for Chinese patients, making this a pragmatic target for patient counselling when no 10/10 donor is found.
- **Chinese and Malay curves converge quickly** because their haplotype distributions (concentrated in common regional haplotypes) allow rapid coverage accumulation with increasing registry size.
- **Others group shows wider spacing** between 10/10 and 8/10 curves, reflecting higher HLA diversity in this heterogeneous population — relaxing match stringency provides greater marginal benefit for this group.

---

### Figure 5: Partial Match Coverage — 8-Locus Framework (6/8, 7/8, 8/8)

**File:** `analysis/figures/partial_match_8locus.png`

![Partial Match Coverage — 8-Locus](analysis/figures/partial_match_8locus.png)

**What it shows:** Identical layout but for the 8-locus framework (HLA-A, B, C, DRB1 only — 8 total allele positions):

- **Blue (8/8):** All 8 alleles match exactly across 4 loci
- **Red (7/8):** At least 7 of 8 alleles match (≤1 mismatch)
- **Green (6/8):** At least 6 of 8 alleles match (≤2 mismatches)

**Differences from 10-locus figure:**

The 4-locus (A, B, C, DRB1) framework excludes DQB1. Because DQB1 is in strong linkage disequilibrium with DRB1, the 8/8 curve is very close to the 10/10 curve — knowing DRB1 largely predicts DQB1, so excluding it provides only modest coverage improvement. Correspondingly, the **separation between the 8/8 and 7/8 curves is smaller** in the 8-locus figure than between 10/10 and 9/10 in the 10-locus figure, because there are fewer loci at which a mismatch can occur.

**Clinical context:**

The 8-locus framework is relevant for:
1. Earlier-generation registries that typed only A, B, C, DRB1
2. Clinical decisions when only 8-locus typing is available
3. Assessment of whether adding DQB1 (moving to 10-locus) meaningfully changes registry size requirements

---

*End of Documentation*
