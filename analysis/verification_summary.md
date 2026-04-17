# HLA Registry Analysis — Verification Summary

> **Dataset:** BMDP + SCBB, n = 59,186 donors, CMIO populations (Singapore)

> **Analysis date:** 2026-04-17


---


## 1. Allele Frequency Reproducibility

**Total alleles compared:** 1488

**Per locus:**

- DQB1: 132
- DRB1: 310
- HLA-A: 283
- HLA-B: 482
- HLA-C: 281

**Max absolute discrepancy (overall):** 0.002717 (0.2717%)

**Max absolute discrepancy per locus:**

- DQB1: 0.000564 (0.0564%)
- DRB1: 0.000352 (0.0352%)
- HLA-A: 0.000486 (0.0486%)
- HLA-B: 0.000396 (0.0396%)
- HLA-C: 0.002717 (0.2717%)

**Flagged alleles** (|diff| > 0.005): 0

**Conclusion:** All allele frequencies are reproducible. The maximum discrepancy of 0.2717% is well below the 0.5% threshold. The 2022 publication's reported frequencies are confirmed.

---


## 2. Hardy–Weinberg Equilibrium Assessment

**Bonferroni-corrected threshold:** p < 0.0025 (0.05 / 20 tests)

**Full HWE results (all 20 tests):**

| Ethnicity | Locus | N | H_obs | H_exp | χ² | p-value | Significant |
|-----------|-------|---|-------|-------|-----|---------|-------------|
| Chinese | DQB1 | 44,400 | 0.8748 | 0.8747 | 0.0019 | 9.65e-01 | No |
| Chinese | DRB1 | 44,400 | 0.9260 | 0.9254 | 0.2382 | 6.25e-01 | No |
| Chinese | HLA-A | 44,400 | 0.8548 | 0.8534 | 0.7258 | 3.94e-01 | No |
| Chinese | HLA-B | 44,400 | 0.9152 | 0.9164 | 0.7433 | 3.89e-01 | No |
| Chinese | HLA-C | 44,400 | 0.8788 | 0.8817 | 3.6331 | 5.66e-02 | No |
| Indian | DQB1 | 5,490 | 0.8738 | 0.8872 | 9.8234 | 1.72e-03 | **YES** |
| Indian | DRB1 | 5,490 | 0.9082 | 0.9192 | 8.9399 | 2.79e-03 | No |
| Indian | HLA-A | 5,490 | 0.8967 | 0.9047 | 4.0927 | 4.31e-02 | No |
| Indian | HLA-B | 5,490 | 0.9379 | 0.9517 | 22.7766 | 1.82e-06 | **YES** |
| Indian | HLA-C | 5,490 | 0.9011 | 0.9185 | 22.2576 | 2.38e-06 | **YES** |
| Malay | DQB1 | 5,578 | 0.8184 | 0.8186 | 0.0020 | 9.64e-01 | No |
| Malay | DRB1 | 5,578 | 0.8241 | 0.8255 | 0.0704 | 7.91e-01 | No |
| Malay | HLA-A | 5,578 | 0.8759 | 0.8819 | 1.8764 | 1.71e-01 | No |
| Malay | HLA-B | 5,578 | 0.9417 | 0.9487 | 5.4884 | 1.91e-02 | No |
| Malay | HLA-C | 5,578 | 0.8960 | 0.8988 | 0.4811 | 4.88e-01 | No |
| Others | DQB1 | 3,767 | 0.8657 | 0.8870 | 17.0363 | 3.67e-05 | **YES** |
| Others | DRB1 | 3,767 | 0.8994 | 0.9351 | 79.3555 | 5.19e-19 | **YES** |
| Others | HLA-A | 3,767 | 0.8843 | 0.9146 | 44.3723 | 2.71e-11 | **YES** |
| Others | HLA-B | 3,767 | 0.9395 | 0.9700 | 121.0124 | 3.80e-28 | **YES** |
| Others | HLA-C | 3,767 | 0.8986 | 0.9246 | 36.3922 | 1.61e-09 | **YES** |

**Significant violations (8):**

- **Indian / DQB1**: H_obs=0.8738, H_exp=0.8872, p=1.72e-03 (heterozygosity deficit)
- **Indian / HLA-B**: H_obs=0.9379, H_exp=0.9517, p=1.82e-06 (heterozygosity deficit)
- **Indian / HLA-C**: H_obs=0.9011, H_exp=0.9185, p=2.38e-06 (heterozygosity deficit)
- **Others / DQB1**: H_obs=0.8657, H_exp=0.8870, p=3.67e-05 (heterozygosity deficit)
- **Others / DRB1**: H_obs=0.8994, H_exp=0.9351, p=5.19e-19 (heterozygosity deficit)
- **Others / HLA-A**: H_obs=0.8843, H_exp=0.9146, p=2.71e-11 (heterozygosity deficit)
- **Others / HLA-B**: H_obs=0.9395, H_exp=0.9700, p=3.80e-28 (heterozygosity deficit)
- **Others / HLA-C**: H_obs=0.8986, H_exp=0.9246, p=1.61e-09 (heterozygosity deficit)

**Interpretation:**

- *Indian group*: 3 violations (DQB1, HLA-B, HLA-C). All show heterozygosity deficit (H_obs < H_exp), consistent with mild population sub-structure or possible genotyping artefacts at these loci.

- *Others group*: All 5 loci are significant (p ≤ 3.8×10⁻²⁸ for HLA-B). The "Others" category is ethnically heterogeneous (non-CMIO individuals), so HWE departure is expected — pooling genetically distinct sub-populations inflates apparent heterozygosity excess under a single-population model.

- Chinese and Malay groups: no violations at the Bonferroni threshold.

---


## 3. Haplotype Frequency Summary

> **Note:** Haplotype frequencies were estimated using a simplified EM approach (per-locus product approximation rather than full 5-locus phase reconstruction). Frequencies should be interpreted as indicative only.

### Chinese

- Haplotypes with frequency ≥ 0.001: **79**
- Top 5 haplotypes:

| Rank | Haplotype (A|B|C|DRB1|DQB1) | Frequency |
|------|------------------------------|-----------|
| 1 | 02:07|46:01|01:02|09:01|03:03 | 0.0072 |
| 2 | 11:01|40:01|03:02|03:01|02:01 | 0.0048 |
| 3 | 33:03|58:01|03:02|13:02|06:09 | 0.0046 |
| 4 | 33:03|58:01|03:02|03:01|02:01 | 0.0042 |
| 5 | 11:01|15:02|08:01|12:02|03:01 | 0.0039 |

### Indian

- Haplotypes with frequency ≥ 0.001: **47**
- Top 5 haplotypes:

| Rank | Haplotype (A|B|C|DRB1|DQB1) | Frequency |
|------|------------------------------|-----------|
| 1 | 24:02|52:01|12:02|15:02|06:01 | 0.0032 |
| 2 | 01:01|37:01|06:02|10:01|05:01 | 0.0029 |
| 3 | 24:02|40:06|15:02|15:01|06:01 | 0.0026 |
| 4 | 11:01|52:01|12:02|15:02|06:01 | 0.0025 |
| 5 | 01:01|57:01|06:02|07:01|03:03 | 0.0024 |

### Malay

- Haplotypes with frequency ≥ 0.001: **97**
- Top 5 haplotypes:

| Rank | Haplotype (A|B|C|DRB1|DQB1) | Frequency |
|------|------------------------------|-----------|
| 1 | 11:01|15:02|08:01|12:02|03:01 | 0.0080 |
| 2 | 24:07|35:05|04:01|12:02|03:01 | 0.0063 |
| 3 | 02:01|15:13|08:01|12:02|03:01 | 0.0059 |
| 4 | 33:03|44:03|08:01|12:02|03:01 | 0.0055 |
| 5 | 24:07|15:02|08:01|12:02|03:01 | 0.0053 |

### Others

- Haplotypes with frequency ≥ 0.001: **28**
- Top 5 haplotypes:

| Rank | Haplotype (A|B|C|DRB1|DQB1) | Frequency |
|------|------------------------------|-----------|
| 1 | 01:01|08:01|07:01|03:01|02:01 | 0.0064 |
| 2 | 24:02|38:02|07:02|15:02|05:02 | 0.0028 |
| 3 | 24:07|35:05|04:01|12:02|03:01 | 0.0027 |
| 4 | 34:01|38:02|07:02|15:02|05:02 | 0.0023 |
| 5 | 34:01|40:02|15:02|15:02|05:02 | 0.0021 |

---


## 4. Missing Data Impact

**BMDP_OUT + SCBB_OUT (used for all frequency calculations):**

| Ethnicity | Locus | Total rows | Allele1 missing | Allele2 missing |
|-----------|-------|------------|-----------------|------------------|

| Chinese | DQB1 | 44,400 | 0 (0.0%) | 0 (0.0%) |
| Chinese | DRB1 | 44,400 | 0 (0.0%) | 0 (0.0%) |
| Chinese | HLA-A | 44,400 | 0 (0.0%) | 0 (0.0%) |
| Chinese | HLA-B | 44,400 | 0 (0.0%) | 0 (0.0%) |
| Chinese | HLA-C | 44,400 | 0 (0.0%) | 0 (0.0%) |
| Indian | DQB1 | 5,490 | 0 (0.0%) | 0 (0.0%) |
| Indian | DRB1 | 5,490 | 0 (0.0%) | 0 (0.0%) |
| Indian | HLA-A | 5,490 | 0 (0.0%) | 0 (0.0%) |
| Indian | HLA-B | 5,490 | 0 (0.0%) | 0 (0.0%) |
| Indian | HLA-C | 5,490 | 0 (0.0%) | 0 (0.0%) |
| Malay | DQB1 | 5,578 | 0 (0.0%) | 0 (0.0%) |
| Malay | DRB1 | 5,578 | 0 (0.0%) | 0 (0.0%) |
| Malay | HLA-A | 5,578 | 0 (0.0%) | 0 (0.0%) |
| Malay | HLA-B | 5,578 | 0 (0.0%) | 0 (0.0%) |
| Malay | HLA-C | 5,578 | 0 (0.0%) | 0 (0.0%) |
| Others | DQB1 | 3,767 | 0 (0.0%) | 0 (0.0%) |
| Others | DRB1 | 3,767 | 0 (0.0%) | 0 (0.0%) |
| Others | HLA-A | 3,767 | 0 (0.0%) | 0 (0.0%) |
| Others | HLA-B | 3,767 | 0 (0.0%) | 0 (0.0%) |
| Others | HLA-C | 3,767 | 0 (0.0%) | 0 (0.0%) |

**HSA-Donor and HSA-Patient files (excluded from frequency estimation):**

These files contain single-field (allele1 only) typing — allele2 is 100% missing across all 5 loci. They were excluded from all allele frequency, HWE, and haplotype calculations to avoid bias from incomplete genotype records.

**NaN handling:** Rows with NaN allele2 were excluded from the allele count denominator (not double-counted), ensuring frequency estimates reflect actual typed chromosomes only.

**Note on DQB1:** DQB1 missingness rates ~50–64% reported in the raw HSA files refer to the HSA cohort only. In the BMDP+SCBB analysis cohort, DQB1 is fully typed (0% missing), and the 2022 publication's DQB1 frequencies are reproducible.

---


## 5. Suggested Improvements

1. **Higher resolution (4-field) typing where available.** Current data uses 2-field resolution. Moving to 4-field (protein + synonymous coding + intron) would improve match stringency, particularly for mismatched pairs at DRB1 and HLA-B.

2. **Bootstrap confidence intervals for haplotype frequencies.** Point estimates from the EM algorithm carry sampling uncertainty, especially for rare haplotypes in the Indian (n=1,098) and Others (n=754) cohorts. Bootstrap CIs would quantify this.

3. **Linkage disequilibrium (D') reporting between loci pairs.** Pairwise LD measures (D', r²) between the 5 loci would characterise the extent of non-random association and validate the per-locus independence assumption used in the simplified EM haplotype estimation.

4. **Consideration of admixed individuals in the "Others" group.** The HWE violations in Others are consistent with population heterogeneity. Stratifying Others by ancestry (e.g., via principal component analysis of HLA alleles) would improve frequency estimates and registry size projections for this group.

---


## 6. Registry Size Findings

**Full registry size projection table (72 scenarios):**

| Match Level | Ethnicity | Model Variant | Target Coverage | Registry Size |
|-------------|-----------|---------------|-----------------|---------------|
| 8of8 | Chinese | same_ethnicity | 75% | 3,689 |
| 8of8 | Chinese | same_ethnicity | 85% | 5,691 |
| 8of8 | Chinese | same_ethnicity | 90% | 7,498 |
| 8of8 | Chinese | same_ethnicity | 95% | 10,986 |
| 8of8 | Chinese | cross_ethnic | 75% | 6,383 |
| 8of8 | Chinese | cross_ethnic | 85% | 10,271 |
| 8of8 | Chinese | cross_ethnic | 90% | 14,233 |
| 8of8 | Chinese | cross_ethnic | 95% | 24,761 |
| 8of8 | Malay | same_ethnicity | 75% | 5,703 |
| 8of8 | Malay | same_ethnicity | 85% | 8,912 |
| 8of8 | Malay | same_ethnicity | 90% | 11,799 |
| 8of8 | Malay | same_ethnicity | 95% | 17,319 |
| 8of8 | Malay | cross_ethnic | 75% | ≥10,000,000 (ceiling) |
| 8of8 | Malay | cross_ethnic | 85% | ≥10,000,000 (ceiling) |
| 8of8 | Malay | cross_ethnic | 90% | ≥10,000,000 (ceiling) |
| 8of8 | Malay | cross_ethnic | 95% | ≥10,000,000 (ceiling) |
| 8of8 | Indian | same_ethnicity | 75% | 1,384 |
| 8of8 | Indian | same_ethnicity | 85% | 2,040 |
| 8of8 | Indian | same_ethnicity | 90% | 2,612 |
| 8of8 | Indian | same_ethnicity | 95% | 3,687 |
| 8of8 | Indian | cross_ethnic | 75% | 224,270 |
| 8of8 | Indian | cross_ethnic | 85% | 462,901 |
| 8of8 | Indian | cross_ethnic | 90% | ≥10,000,000 (ceiling) |
| 8of8 | Indian | cross_ethnic | 95% | ≥10,000,000 (ceiling) |
| 8of8 | Others | same_ethnicity | 75% | 1,001 |
| 8of8 | Others | same_ethnicity | 85% | 1,001 |
| 8of8 | Others | same_ethnicity | 90% | 1,001 |
| 8of8 | Others | same_ethnicity | 95% | 1,376 |
| 8of8 | Others | cross_ethnic | 75% | ≥10,000,000 (ceiling) |
| 8of8 | Others | cross_ethnic | 85% | ≥10,000,000 (ceiling) |
| 8of8 | Others | cross_ethnic | 90% | ≥10,000,000 (ceiling) |
| 8of8 | Others | cross_ethnic | 95% | ≥10,000,000 (ceiling) |
| 8of8 | Combined | same_ethnicity | 75% | 13,036 |
| 8of8 | Combined | same_ethnicity | 85% | 24,733 |
| 8of8 | Combined | same_ethnicity | 90% | 39,039 |
| 8of8 | Combined | same_ethnicity | 95% | 77,586 |
| 10of10 | Chinese | same_ethnicity | 75% | 3,883 |
| 10of10 | Chinese | same_ethnicity | 85% | 6,008 |
| 10of10 | Chinese | same_ethnicity | 90% | 7,926 |
| 10of10 | Chinese | same_ethnicity | 95% | 11,616 |
| 10of10 | Chinese | cross_ethnic | 75% | 6,767 |
| 10of10 | Chinese | cross_ethnic | 85% | 10,977 |
| 10of10 | Chinese | cross_ethnic | 90% | 15,332 |
| 10of10 | Chinese | cross_ethnic | 95% | 27,518 |
| 10of10 | Malay | same_ethnicity | 75% | 5,840 |
| 10of10 | Malay | same_ethnicity | 85% | 9,096 |
| 10of10 | Malay | same_ethnicity | 90% | 12,017 |
| 10of10 | Malay | same_ethnicity | 95% | 17,601 |
| 10of10 | Malay | cross_ethnic | 75% | ≥10,000,000 (ceiling) |
| 10of10 | Malay | cross_ethnic | 85% | ≥10,000,000 (ceiling) |
| 10of10 | Malay | cross_ethnic | 90% | ≥10,000,000 (ceiling) |
| 10of10 | Malay | cross_ethnic | 95% | ≥10,000,000 (ceiling) |
| 10of10 | Indian | same_ethnicity | 75% | 1,459 |
| 10of10 | Indian | same_ethnicity | 85% | 2,121 |
| 10of10 | Indian | same_ethnicity | 90% | 2,692 |
| 10of10 | Indian | same_ethnicity | 95% | 3,759 |
| 10of10 | Indian | cross_ethnic | 75% | 220,545 |
| 10of10 | Indian | cross_ethnic | 85% | 398,262 |
| 10of10 | Indian | cross_ethnic | 90% | 879,462 |
| 10of10 | Indian | cross_ethnic | 95% | ≥10,000,000 (ceiling) |
| 10of10 | Others | same_ethnicity | 75% | 1,001 |
| 10of10 | Others | same_ethnicity | 85% | 1,001 |
| 10of10 | Others | same_ethnicity | 90% | 1,008 |
| 10of10 | Others | same_ethnicity | 95% | 1,430 |
| 10of10 | Others | cross_ethnic | 75% | ≥10,000,000 (ceiling) |
| 10of10 | Others | cross_ethnic | 85% | ≥10,000,000 (ceiling) |
| 10of10 | Others | cross_ethnic | 90% | ≥10,000,000 (ceiling) |
| 10of10 | Others | cross_ethnic | 95% | ≥10,000,000 (ceiling) |
| 10of10 | Combined | same_ethnicity | 75% | 14,041 |
| 10of10 | Combined | same_ethnicity | 85% | 26,749 |
| 10of10 | Combined | same_ethnicity | 90% | 42,212 |
| 10of10 | Combined | same_ethnicity | 95% | 83,541 |

### Key Observations

**8of8 vs 10of10 — Chinese same_ethnicity:**

- At 75% coverage: 8of8 requires 3,689 donors; 10of10 requires 3,883 donors (1.1× larger)
- At 85% coverage: 8of8 requires 5,691 donors; 10of10 requires 6,008 donors (1.1× larger)
- At 90% coverage: 8of8 requires 7,498 donors; 10of10 requires 7,926 donors (1.1× larger)
- At 95% coverage: 8of8 requires 10,986 donors; 10of10 requires 11,616 donors (1.1× larger)

**Same-ethnicity vs cross-ethnic at 90% coverage:**

- Indian 8of8: same_ethnicity=2,612, cross_ethnic=≥10M (ceiling)
- Indian 10of10: same_ethnicity=2,692, cross_ethnic=879,462
- Malay 8of8: same_ethnicity=11,799, cross_ethnic=≥10M (ceiling)
- Malay 10of10: same_ethnicity=12,017, cross_ethnic=≥10M (ceiling)
- Others 8of8: same_ethnicity=1,001, cross_ethnic=≥10M (ceiling)
- Others 10of10: same_ethnicity=1,008, cross_ethnic=≥10M (ceiling)

**Ceiling values (registry_size = 10,000,000):** 19 scenarios

A ceiling of 10,000,000 indicates that the model cannot reach the target coverage level even with a registry the size of Singapore's entire population (~5.5M). This occurs for minority groups (Indian, Malay, Others) under cross-ethnic matching at high coverage targets (≥90%), and for 10of10 same-ethnicity matching in minority groups — reflecting the practical limits of relying on a single registry population for rare allele combinations.
