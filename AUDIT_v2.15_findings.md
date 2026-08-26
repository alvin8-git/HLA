# Integrity audit of `HLA_Registry_Size_CMIO_v2.15.docx`


> **Historical record — describes manuscript v2.15 (superseded).** The figures here
> are the v2.15 ones and are left unchanged on purpose; several were corrected in
> v2.15b/v2.15c, and this file is the record of *why*. Current numbers are in
> [README](README.md) and [Documentation.md](Documentation.md). See
> [VERSION.md](VERSION.md) for the chain of corrections.

Prepared 2026-08-20 during preparation of the independent reanalysis
(`HLA_Registry_Size_CMIO_BoneMarrowTransplantation.md`). Every finding below was
verified against the pipeline's own CSV outputs or by re-running its own code; the
verification command or file is named in each entry.

Severity: **A** = invalidates a headline claim · **B** = wrong number or false statement
in the text · **C** = presentation/traceability.

---

## A1. Rare-haplotype truncation invalidates all absolute registry sizes

**Claim in v2.15:** "Same-ethnicity registries of approximately 40,000–45,000 donors are
required per CMIO group to achieve 95% coverage at 10/10 HLA matching" (Abstract, §3.1,
§5, Recommendation 1).

**Finding.** `hwe_test.run_em_haplotypes` retains only haplotypes with frequency
≥ 0.001, and `registry_model.extract_and_normalize` then rescales the survivors to sum
to 1.0. The retained set is 123–144 haplotypes per group and accounts for only
**36–53% of total haplotype frequency mass**:

| Group | Haplotypes retained | Frequency mass retained |
|---|---|---|
| Chinese | 140 | 51.7% |
| Malay | 137 | 52.9% |
| Indian | 144 | 40.2% |
| Others | 123 | 35.6% |

(`analysis/data/haplo_freqs_em.csv`, `sum_freq` per ethnicity; also
`em_validation_summary.csv` column `freq_coverage_em`.)

Renormalising deletes the rare tail and inflates the frequencies of the survivors. Since
`C(N) = Σ F_k[1−(1−F_k)^N]` approaches its limit slowly exactly when `F_k` is small,
the deleted tail is what determines high-coverage behaviour. Re-running the project's
own EM with the floor removed:

| Group | Floor | Haplotypes | Mass | Coverage at N=50,000 | N* at 95% |
|---|---|---|---|---|---|
| Chinese | 0.1% | 143 | 52.5% | 95.6% | 45,148 |
| Chinese | none | 57,934 | 100% | **36.8%** | **> 10,000,000** |
| Others | 0.1% | 126 | 35.7% | 97.4% | 34,110 |
| Others | none | 65,020 | 100% | **17.6%** | **> 10,000,000** |

The same manipulation applied to the GENE[RATE] reference distributions (the
manuscript's own designated validation standard, ref [3]) reproduces the effect:
N* understated by **120× to 480×**.

**External check the manuscript never performed.** The truncated model predicts 100.0%
10/10 coverage for every CMIO group from a 500,000-donor registry. No registry achieves
complete matching for any population; the US registry at 10.5 million donors delivers
~75% 8/8 for its best-served group (Gragert 2014, ref [2] — cited in the manuscript).
The untruncated model predicts 94.9% for Chinese at 10 million donors, which is the
right order. The headline numbers are therefore inconsistent with the manuscript's own
key reference.

**The damage is a threshold, and the pipeline sits on the wrong side of it.** Sweeping
the floor on a single unfloored EM estimate (`paper_BMT_workdir/floor_curve.py`) locates
it precisely for the Chinese group:

| Floor | Haplotypes | Mass | Coverage @ N=50,000 | N* at 90% |
|---|---|---|---|---|
| none | 57,934 | 100% | 36.8% | 6,339,849 |
| 10⁻⁵ | 2,697 | 99.8% | 37.0% | 6,225,401 |
| 10⁻⁴ | 2,309 | 96.9% | 39.5% | 4,741,908 |
| 5×10⁻⁴ | 326 | 65.1% | 80.3% | 120,838 |
| **10⁻³ (pipeline default)** | **143** | **52.5%** | **95.6%** | **24,868** |

A 10⁻⁴ floor would have been almost free (3.1% of mass, 2.7 coverage points). The
pipeline's 10⁻³ costs a further 44% of mass and understates N* at 90% coverage by
**255-fold** — a measured ratio, not a censored one.

**This is not an EM phase-ambiguity artefact.** The obvious objection is that EM over
unphased genotypes manufactures spurious rare haplotypes, so the tail is noise and
truncation is denoising. It is not: dropping all 55,237 haplotypes below 10⁻⁵ (95.3% of
the count) removes 0.19% of mass and moves coverage by 0.2 points. The bias is carried
by ~2,200 haplotypes in the 10⁻⁴–10⁻³ band, whose existence GENE[RATE] independently
corroborates (2,196 Chinese haplotypes ≥10⁻⁴, versus 2,309 in our own unfloored run).

**One-line fix with a large effect:** change `freq_threshold` from `0.001` to `0.0001` in
`hwe_test.run_em_haplotypes` and re-run. Everything downstream changes qualitatively.

**Consequence.** The 40,000–45,000 target cannot be defended. Reproduce with
`paper_BMT_workdir/floor_curve.py`, `em_notrunc.py` and `truncation_bias.csv`.

---

## A2. The 9/10 benefit is understated by half

**Claim in v2.15:** "Relaxing from 10/10 to 9/10 matching roughly halves the required
registry size… For Chinese patients, the 9/10 registry requirement at 95% coverage is
approximately 20,000–22,000 donors, compared with 42,847 at 10/10" (§3.4, Rec 3,
Abstract, §5).

**Finding.** Computed with the pipeline's own `06_partial_match_plots.py` functions
(`parse_haplotypes` → `compute_partial_match_probs` → `coverage_curve`), Chinese ≥9/10
reaches 95.1% coverage at **N = 11,000**, and N* at 95% is **10,740**. The reduction is
**~3.9-fold, not 2-fold**; the stated 20,000–22,000 is roughly double the true value.

| Group | 10/10 @95% | ≥9/10 @95% | Ratio | ≥8/10 @95% |
|---|---|---|---|---|
| Chinese | 41,974 | 10,740 | 3.9× | 1,794 |
| Malay | 41,234 | 8,351 | 4.9× | 1,651 |
| Indian | 44,371 | 10,350 | 4.3× | 2,381 |
| Others | 32,690 | 12,973 | 2.5× | 4,007 |

(`paper_BMT_workdir/partial_match_nstar.csv`.) The 10/10 column reproduces Table 1 to
within 2%, confirming the method agrees where the manuscript reports a value — the
9/10 figure appears never to have been computed, only estimated by eye from Figure 3.

**Consequence.** The paper's single most actionable recommendation is understated by
2×, in the direction that makes it look less attractive.

---

## A3. Table 1's "Weighted Average" row is not a weighted average

**Claim in v2.15:** Table 1 row "Weighted Average†" = 236,906 at 95%, footnoted as
"Singapore population weights (Chinese 77%, Malay 8%, Indian 9%, Others 6%)". §3.1 then
separately states the combined N* is "~42,332" (Table 6).

**Finding.** 236,906 is the `Combined` row of `registry_size_targets.csv`, produced by
`get_combined_haplotype_freqs` — a *pooled-haplotype* model in which both patient and
donor are drawn from a mixed pool. It is not an average of the per-group values. The
arithmetic weighted average is:

    0.77(42,871) + 0.08(41,779) + 0.09(44,863) + 0.06(32,360) = 42,332

which is exactly the Table 6 "SG population" figure. So the manuscript reports two
numbers that differ **5.6×**, labels the wrong one "Weighted Average", and never
explains the discrepancy to the reader.

**Fix.** Relabel the Table 1 row "Combined pooled-registry model" and state explicitly
that the weighted average of per-group targets is 42,332.

---

## B1. Bootstrap donor counts exceed the stated sample size

**Claim in v2.15 §2.4:** concentration parameters use "the actual 5-locus donor count per
ethnicity (Chinese: 45,754; Malay: 5,868; Indian: 5,586; Others: 3,941)".

**Finding.** These sum to **61,149**, which exceeds the stated total of 59,186 donors and
cord blood units. The counts reproducible from `hla_clean.csv` for individuals with all
five loci typed are Chinese 44,400 / Malay 5,578 / Indian 5,490 / Others 3,767
(= 59,235). The hardcoded values in `analysis/09_bootstrap_ci.py:51` match no
reproducible count. Inflated concentration parameters make the reported CIs
approximately 1.5% narrower than they should be — small, but the internal contradiction
(part > whole) is the kind a reviewer will find.

**Related:** three different values for the Others group circulate — 3,941 (bootstrap),
3,847 (clustering input), 3,767 (`hla_clean.csv` five-locus). §2.4 and §3.7 should state
which denominator each analysis used.

---

## B2. Linkage-disequilibrium claims are false for the Others group

**Claims in v2.15:** §2.2 "especially strong between DRB1 and DQB1 (D′ ≥ 0.94 in all CMIO
groups) and between B and C (D′ ≥ 0.95)"; §3.2 "D′ = 0.94–0.99"; Rec 5 "strong
DRB1–DQB1 LD".

**Finding** (`analysis/data/ld_report.csv`):

| Pair | Chinese | Malay | Indian | Others |
|---|---|---|---|---|
| DRB1–DQB1 | 0.9873 | 0.9416 | 0.9561 | **0.9339** |
| B–C | 0.9538 | 0.9755 | 0.9865 | **0.9488** |

Others fails both stated thresholds. The correct range for DRB1–DQB1 is **0.93–0.99**,
not 0.94–0.99. The qualitative conclusion (LD is strong; DQB1 is nearly free) is
unaffected.

---

## B3. The 8/8-versus-10/10 gap is misstated

**Claims in v2.15:** §3.2 "typically 600–1,200 fewer donors at 95% coverage"; Rec 5
"8/8 and 10/10 registry sizes differ by less than 5%".

**Finding** (Tables 1 and 2 of the manuscript itself):

| Group | 10/10 | 8/8 | Difference | % |
|---|---|---|---|---|
| Chinese | 42,847 | 42,115 | 732 | 1.7% |
| Malay | 40,032 | 36,176 | **3,856** | **9.6%** |
| Indian | 43,855 | 42,738 | 1,117 | 2.5% |
| Others | 31,181 | 30,525 | 656 | 2.1% |

Malay falls outside both stated bounds. Correct to "600–3,900 donors (1.7–9.6%)".

---

## B4. EM convergence is described as monotone; it is not

**Claim in v2.15 §4.1 / S1:** "N* stabilising near 42,000 donors above ~20,000 samples;
at 5,000 the estimate is 45,148 versus 41,727 at the full sample — a conservative 8.2%
overestimate… Convergence is near-complete above ~20,000 samples."

**Finding** (`analysis/data/em_convergence.csv`): the sequence is 482,681 → 111,647 →
59,324 → 50,685 → **45,148** (n=5,000) → 47,225 → 47,346 → 48,088 (n=15,000) → 44,326 →
41,787 → 41,727. The excursion **above** the 5,000-sample value between n = 7,500 and
n = 15,000 contradicts the "conservative overestimate" framing, which assumes the curve
descends monotonically from the cap.

**Additionally:** a fresh EM run at the same nominal 5,000 cap reproduces 45,148, but
the archived `haplo_freqs_em.csv` used for Table 1 yields 42,871 (143 vs 140 retained
haplotypes) — a 5.4% gap between nominally identical configurations, against a reported
CI of ±0.5%. The provenance of the archived table should be documented.

---

## B5. A validation output showing a 36–125× discrepancy is unreported

**Finding.** `analysis/data/match_rate_comparison.csv` records, for all four groups, a
model-predicted match probability of 0.0080–0.0275 against an observed match rate of
**1.00**. This appears nowhere in the manuscript.

The 1.00 is a selection artefact — the HSA pairs are transplanted pairs and therefore
matched by construction — so it is not evidence of model failure. But an unexplained
1.00 sitting in a validation output is precisely what a reviewer will seize on, and its
omission means the manuscript never states that this dataset **cannot** serve as an
outcome-based validation of the coverage model. Report it and explain it, or remove the
analysis.

---

## B6. Two different validations are both described as "validation"

§2.2 says haplotype frequencies were "validated against the HLA-net GENE[RATE] database
[3]" (Spearman ρ = 0.91–0.99, `em_validation_summary.csv`); §3.6 reports validation
against patient haplotypes (ρ = 0.70). These are different comparisons with very
different strengths, and a reader may carry the stronger figure into the weaker claim.
Name them distinctly.

---

## B7. Others fails HWE at all five loci; only Indian is disclosed

**Claim in v2.15 §4.1:** "HWE departures were detected in Indian and Others groups at
several loci."

**Finding** (`analysis/data/hwe_results.csv`): the Others group departs significantly at
**all five loci** — DQB1 p=3.7×10⁻⁵, DRB1 p=5.2×10⁻¹⁹, HLA-A p=2.7×10⁻¹¹,
HLA-B p=3.8×10⁻²⁸, HLA-C p=1.6×10⁻⁹ — with observed heterozygosity **below** expected at
every locus. Chinese and Malay depart at none; Indian at three.

Uniform heterozygote deficit across all loci is the textbook **Wahlund signature** of
pooling stratified subpopulations. It is therefore direct statistical corroboration of
the §3.7 three-cluster finding — and simultaneously invalidates the random-mating
assumption behind `F(hᵢ,hⱼ)=2fᵢfⱼ` for every pooled-Others figure in Tables 1–3. The
manuscript's Limitations paragraph singles out the *milder* Indian violation and treats
Others as an afterthought, which inverts the severity ordering. Pooled-Others estimates
should be withdrawn in favour of the cluster-level figures.

---

## C1. Blind spots — not errors, but unaddressed

1. **Cord blood units are pooled with adult donors** (CBUs are 2.3% Chinese, 1.7% Malay, 1.7% Indian, 6.2% Others — too few to bias frequencies, but) the output is
   framed as adult "donor recruitment targets". CBUs have different matching standards
   (the reviewer's own comment: minimum 4/6 or 6/8) and cannot be "recruited". The
   manuscript never separates them.
2. **No DPB1, hence no 12/12.** The reviewer states standard of care is a 12/12 panel
   with 10/10 minimum and 8/8 obsolete. v2.15 still presents 8/8 as a co-primary
   analysis (Table 2, Rec 5) rather than as a sensitivity analysis or historical
   comparator.
3. **The `Combined` / cross-ethnic model is never sanity-checked** against the fact that
   Singapore already operates a ~59,000-donor mixed registry with observable match
   rates. The observed rate is the obvious external validation and is not used.
4. **No comparison to published registry performance.** Gragert 2014 is cited for the
   method but never for its numbers, which is where the artefact would have surfaced.
5. **The Others silhouette (0.24)** is weak in absolute terms; the ancestry labels are
   carried by the haplotype signatures (Table 5), which is a much stronger argument than
   the clustering metric. The manuscript leans on the wrong evidence.

---

## Summary

| ID | Severity | One line |
|---|---|---|
| A1 | A | 0.1% haplotype floor discards 47–64% of frequency mass; headline N* understated 120–480× |
| A2 | A | 9/10 reduction is ~4×, not ~2×; Chinese 9/10 N* is 10,740, not 20,000–22,000 |
| A3 | A | Table 1 "Weighted Average" (236,906) is the pooled model; the real weighted average is 42,332 |
| B1 | B | Bootstrap donor counts sum to 61,149 > stated total 59,186 |
| B2 | B | D′ ≥ 0.94 / ≥ 0.95 claims false for Others (0.9339 / 0.9488) |
| B3 | B | 8/8-vs-10/10 gap is 600–3,900 donors (to 9.6%), not 600–1,200 (<5%) |
| B4 | B | EM convergence non-monotone; 5.4% gap between nominally identical runs vs ±0.5% CI |
| B5 | B | Unreported validation output: predicted 0.008–0.028 vs observed 1.00 |
| B6 | C | Two distinct "validations" share one name |
| B7 | B | Others fails HWE at all 5 loci (Wahlund); manuscript flags only the milder Indian case |
| C1 | C | Blind spots: cord blood pooling, no DPB1/12-12, no external benchmark, weak silhouette leaned on |

A1, A2 and A3 must be resolved before v2.15 is submitted anywhere. A1 is not
fixable by editing text — it requires recomputation without the frequency floor, and
the resulting conclusions differ qualitatively from those currently stated.
