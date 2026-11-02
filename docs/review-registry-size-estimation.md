# Independent review: registry-size estimation for the Singapore CMIO registry

*Reviewer: independent re-analysis, 2026-11-01. Scope: the registry-size model behind
`HLA_Registry_Size_CMIO_v2.15c` and the 1e-6 re-run in `analysis/data/`. No pipeline file was
changed. Every experiment script and its log are in `/data/alvin/tmp/hla_review/`.*

---

## 1. Verdict

**Neither current set of registry targets can be published.**

- **The v2.15c targets (1e-3 floor) are not credible.** For Chinese, the 0.1% model's 95% target is
  42,871. On the registry it was built from, it predicts that **95% of Chinese patients already find
  a 10/10 donor**. Independent cord-blood units find one **26% of the time**.
- **The 1e-6 targets rest on a model that is the only one consistent with independent data.**
  The targets themselves are 10 to 1,000 times beyond the registry's size, which is further than
  any method can be validated.

Two problems outrank the choice of floor:

1. **Duplicate records.** Between 4,370 and 12,058 of the 57,720 usable BMDP records appear to be
   duplicates. They inflate the "observed match rate" that was used to validate the 1e-3 floor.
2. **A mis-scaled validation.** The validation that favoured 1e-3 scaled conditional coverage by
   retained haplotype mass instead of retained *diplotype* mass, which is mass squared.

**What is publishable now:** coverage at today's size, validated against independent units; gains
from a two- to three-fold expansion; and a bounded statement that a 95% 10/10 target is out of reach
through same-ethnicity recruitment alone. That last statement agrees with the US registry
(Gragert 2014).

---

## 2. The question and the unit

- **Estimand.** C(N) is the probability that a patient drawn from a CMIO group has at least one
  allele-level 10/10 (or 8/8) genotype-identical donor in a same-group registry of N independent
  donors. N\* is the smallest N with C(N) ≥ target.
- **Model.** `registry_model.py` computes
  C(N) = Σ_d f_d·[1 − (1 − f_d)^N] over HWE diplotypes built from EM haplotype frequencies.
- **Independent unit.** The unit is a donor (a person), not a record. This is the assumption the
  data break (§4.1).

A correction to the brief behind this review: **0.5% is not a registry cut-off anywhere in the
repository.** It is the tolerance for the allele-frequency audit against the 2022 paper
(`02_allele_freq.py`, threshold 0.005). The manuscript's floor is 1e-3 (0.1%). The re-run used 1e-6.

---

## 3. Strengths

- The pipeline is reproducible, and its core maths is unit-tested.
- The allele-frequency audit is careful: zero alleles flagged, with a maximum discrepancy of 0.27%.
- The project has already done the hardest self-critique. `docs/explanation-frequency-floor.md`
  identifies the 1/(2n) boundary, shows the 2,098-fold swing, shows the ranking inversion, and
  separates the claims that survive a floor change from those that don't.
- `check_freshness.py` and the snapshot discipline make stale numbers detectable.
- The comparative conclusions remain sound: cross-ethnic penalty, 9/10 relaxation, and DRB1–DQB1
  linkage.

---

## 4. Concerns, by severity

### 4.1 Critical: duplicate donor records contaminate the benchmark and the sample size

**Observation.** Treat each allele-level 10/10 genotype as a type. Within BMDP, 41.9% of Chinese
records have a genotype twin. Twins should be scattered at random through the file. Instead they sit
next to each other (`twin_structure.log`):

| Group | Twin pairs | Gap ≤ 10 rows | Expected if random | Median gap (random) |
|---|---|---|---|---|
| Chinese | 12,277 | 27.9% | 0.02% | 2,101 (17,325) |
| Malay | 841 | 36.3% | 0.24% | 2,719 (17,375) |
| Indian | 575 | 64.5% | 0.00% | **1** (15,859) |
| Others | 345 | 77.4% | 0.29% | **1** (16,562) |

Every twin pair within 10 rows, **4,374 of 4,374**, is also character-identical in the original
high-resolution typing strings of `HLA Data combined.xlsx`. That includes the allele column order
and the G-group suffixes. Twin pairs further apart are identical only 77.7% of the time
(`dup_resolution.log`). Relatives typed independently would not reproduce the same report string
every time. These are duplicate entries of one donor.

The source sheet has no donor ID column, so the extent cannot be pinned down:

| Bound | Records |
|---|---|
| Lower: raw-identical within 10 rows | 4,370 (7.6%) |
| Upper: every raw-identical repeat | 12,058 (20.9%) of 57,720 usable BMDP records |

**Consequences.**

- The headline n=59,186 overstates the number of donors by 7.6% to 20.9%.
- EM haplotype frequencies are over-weighted toward duplicated genotypes.
- **Most importantly**, the "observed 10/10 twin rate" that `docs/explanation-frequency-floor.md`
  uses to show 1e-3 "predicts reality better" is inflated by records matching themselves.

**The clean benchmark.** The 1,450 SCBB cord units are typed and banked independently of BMDP. The
fraction of SCBB units with a 10/10 twin in BMDP is an uncontaminated estimate of C(N) at the
BMDP size. Deduplication cannot change it, because it depends only on which genotypes BMDP contains
(`dedup_coverage.log`):

| Group | Twin rate, as ingested | Nearby duplicates removed | All raw repeats removed | **SCBB probes** (95% CI) |
|---|---|---|---|---|
| Chinese | 41.9% | 31.0% | 9.8% | **26.1%** (23.4–28.9), n=1,000 |
| Malay | 25.5% | 16.3% | 7.7% | **7.8%** (3.2–15.4), n=90 |
| Indian | 20.0% | 7.1% | 2.4% | **4.4%** (1.2–10.9), n=91 |
| Others | 18.9% | 4.4% | 2.0% | **3.2%** (1.3–6.4), n=220 |

The as-ingested rate sits far above the independent probe rate in every group. The truth lies
between the two deduplication bounds.

*Certainty that duplicates exist: high. Certainty about their exact number: low until registry IDs
are obtained.*

### 4.2 Critical: the 1e-3 targets measure a different, much easier estimand

`04_registry_model.py` drops haplotypes below the floor and **renormalises** the survivors. At 1e-3
only 36–52% of haplotype mass survives:

| Group | Haplotype mass retained | Retained mass squared |
|---|---|---|
| Chinese | 0.4934 | 0.243 |
| Malay | 0.5235 | 0.274 |
| Indian | 0.4013 | 0.161 |
| Others | 0.3563 | 0.127 |

The resulting C(N) is coverage among patients whose two haplotypes are both common. The
manuscript's "95%" means 95% of that subgroup. For Chinese that is about 0.4934², or 24% of
patients.

Scored against the same registry at its actual size (`floor_sweep.log`):

| Floor | Chinese conditional C(n) | Chinese unconditional C(n) | MAE all groups, unconditional | MAE, conditional |
|---|---|---|---|---|
| 1e-6 | 0.295 | 0.295 | 12.6 pts | 12.6 pts |
| 1/(2n) per group | 0.321 | 0.294 | 12.7 | 7.8 |
| 1e-4 | 0.494 | 0.285 | 12.9 | 8.2 |
| 1e-3 | **0.949** | **0.195** | **16.2** | **46.9** |

*MAE is mean absolute error against the as-ingested twin rate, which is itself inflated (§4.1).*

The floor document reports that 1e-3 lands "within 7 points" of observation. That figure multiplies
conditional coverage by retained haplotype mass (0.949 × 0.49 ≈ 0.47). A patient needs both
haplotypes retained, so the correct factor is mass squared. With it, the 1e-3 model gives 19.5%,
not 49%. **The stated argument for 1e-3 ("the only version that has been validated") does not
hold.**

The external benchmark agrees. Gragert et al. (2014, Table 1) report a **41%** likelihood of an
available 8/8 donor for US Chinese patients, searching the whole multi-million-donor NMDP registry.
A model claiming 95% at 10/10 with 43,000 donors is irreconcilable with that.

*Certainty: high. This is arithmetic on the repository's own outputs.*

### 4.3 Important: once the benchmark is cleaned, the 1e-6 model is consistent with independent data

The 1e-6 haplotype–HWE model, unconditional and evaluated at the deduplicated BMDP sizes, against
the SCBB probes (`model_vs_probe.log`):

| Group | Model @ lower-bound dedup n | Model @ upper-bound dedup n | SCBB probe (95% CI) |
|---|---|---|---|
| Chinese | 28.5% | 26.6% | 26.1% (23.4–28.9) |
| Malay | 15.2% | 14.7% | 7.8% (3.2–15.4) |
| Indian | 5.7% | 5.6% | 4.4% (1.2–10.9) |
| Others | 2.9% | 2.9% | 3.2% (1.3–6.4) |

All four groups fall inside the probe interval.

- **Malay** sits at its upper edge. The model may be optimistic for Malay, so it may understate N\*.
- **Others** matches despite failing HWE at every locus. Per-locus inbreeding coefficients are
  F = 0.024–0.038 (Wahlund), which nudges genotype repetition up only slightly (`diagnose_gap.log`).
- **Chinese** shows no HWE departure at any locus (|F| ≤ 0.003). HWE is not the source of error.

The floor document reads 1e-6 as a "failed prediction" that comes from EM over-fragmentation. That
reading follows from the contaminated benchmark. Against clean probes the 1e-6 model is the best
calibrated option available. The remaining biases point toward the model being *optimistic*:
haplotypes never observed get zero mass, and EM was fitted on duplicated data. If anything, its N\*
values are lower bounds.

*Certainty: moderate. The probes are small for Malay, Indian and Others, and SCBB donors are a proxy
for patients.*

### 4.4 Important: below roughly 1e-4 the floor barely changes coverage today; it decides the 95% target

Coverage at current size hardly moves across floors from 1e-6 to 1e-4 (Chinese 0.295 → 0.285).
The far targets are a different matter.

- **An unconditional C(N) can never exceed the square of the retained mass.** At 1/(2n), Chinese
  keeps 96.7% of mass, so its ceiling is 93.5% and **a 95% target becomes unreachable**.
- **Reaching 95% requires keeping at least 97.5% of mass.** For every group that means keeping the
  sub-singleton tail, the haplotypes the sample has seen at most once.

So the 95% target is, by construction, a statement about haplotypes nobody has counted. That is why
it swings 2,098-fold. **No floor value resolves this. It is a limit of the data, not a tuning
choice.**

### 4.5 Important: high targets lie far beyond any validated extrapolation range

Floor-free genotype-level estimators avoid EM and phasing altogether
(`genotype_coverage.py`, `dedup_coverage.py`):

- **Rarefaction**, for N ≤ n. This is the exact expected coverage of a subsample, after Chao et
  al. 2014.
- **Chao–Jost coverage extrapolation**, Chao & Jost 2012. Its functional form is taken from
  the literature summary; the exact form was not retrieved.
- **Pitman–Yor posterior prediction**, a two-parameter heavy-tail model fitted by maximum likelihood.

**Held-out test.** Each estimator was fitted on half (or a quarter) of the donors and asked to
predict coverage at the full size (`genotype_coverage.log`, as-ingested data):

| Group | Chao–Jost error, 2× | Pitman–Yor error, 2× | Chao–Jost error, 4× | Pitman–Yor error, 4× |
|---|---|---|---|---|
| Chinese | +1.1 pts | −4.2 | +6.1 | −7.7 |
| Malay | +1.1 | −1.3 | +4.5 | −2.1 |
| Indian | +0.6 | −0.5 | −0.2 | −2.3 |
| Others | +0.2 | −0.7 | −0.1 | −1.6 |

- Chao–Jost is accurate to about 1 point at two-fold extrapolation and grows optimistic beyond it.
- Pitman–Yor is pessimistic.
- The literature agrees on the limits. Colwell et al. 2012 call extrapolation reliable to about
  2–3× the sample. Orlitsky et al. 2016 prove a limit near log n, which is about 8–11× here.

**After deduplication Pitman–Yor fails.** The minority groups become almost all singletons
(K/n ≈ 0.96 for Others), the tail parameter is unidentifiable, and the fitted model returns
infinite N\*. Converging on a fit did not mean the parameters were identified.

Extrapolation factors that the current targets demand:

| Group | n (lower-bound dedup) | 1e-6 model N\*95 | Factor | Chao–Jost N\*95 (dedup) | Factor |
|---|---|---|---|---|---|
| Chinese | 39,922 | 87,384,114 | 2,189× | 701,730 | 18× |
| Malay | 5,174 | 16,552,048 | 3,199× | 149,883 | 29× |
| Indian | 4,993 | 28,817,950 | 5,772× | 284,629 | 57× |
| Others | 3,261 | 26,762,600 | 8,207× | 276,059 | 85× |

The Chao–Jost figures are lower bounds. The method is optimistic at long range and residual
duplicates remain. **Every 90–95% target in this project lies outside the range where any estimator
has been, or can be, validated with these samples.**

### 4.6 Minor: pooling of the "other" haplotype

`get_diplotype_frequencies` and `diplotype_freq_vector` (`top_k_pct=0.99`) merge the last 1% of
haplotype mass into one pseudo-haplotype. Two patients carrying different rare haplotypes are then
scored as matched. The effect measured here is at most 1.3 points of coverage. It is the wrong
direction, though, and is easy to remove: use `top_k_pct=1.0` and the streamed computation below.
Scripts 09, 11, 13, 15 and 16 inherit it.

### 4.7 Minor: compute was never the binding constraint

One full-tail coverage evaluation at 1e-6 (Chinese, K=9,574, 45.8 million diplotypes) takes
**0.4 seconds** when streamed row by row over a float array (`model_vs_probe.py`).

- A binary search for N\* takes about 40 evaluations.
- The 9-hour, out-of-memory bootstrap and the 36-minute `04` run come from materialising
  string-labelled diplotype DataFrames, not from the mathematics.
- EM takes 59 seconds.

**A full non-parametric bootstrap that re-runs EM per replicate is feasible**: about 1–2 CPU-minutes
per replicate, so B=200 is a few hours across 8 cores.

---

## 5. Answers to the three questions

### 5.1 Is 0.1% the best trade-off between accuracy and compute?

**No, on both counts.**

- **Accuracy.** 1e-3 is the worst-performing option for the honest, unconditional estimand. It
  silently changes the estimand to "patients with two common haplotypes".
- **Compute.** Once the coverage sum is vectorised, compute does not separate the floors.
- **Where the floor should sit.** Any value at or below 1/(2n) gives almost identical coverage at
  today's size. Below about 1e-5 nothing changes. 1e-6 is a reasonable working value because it sits
  below every group's 1/(2n).
- **What a floor cannot do.** It cannot make 95% targets estimable. That needs data on haplotypes
  not yet observed.

### 5.2 Which methods can give a publishable, representative estimate?

| Method | Bias direction | Data needed (have?) | Compute | Publishable? |
|---|---|---|---|---|
| Haplotype EM + HWE, fixed floor 1e-3, renormalised | Strongly optimistic: wrong estimand | have | trivial | **N** |
| Haplotype EM + HWE, floor ≤ 1/(2n), unconditional, no pooling | Mildly optimistic: unseen haplotypes get zero mass | have; dedup needed | seconds per N\* | **Y, primary**, if deduplicated and probe-validated |
| Empirical calibration against independent probes (SCBB; patient genotypes) | None at N = n; validates only the current size | SCBB have, small; patient genotypes **lack** (Patient.txt is haplotype-level) | trivial | **Y**, as the validation table. Novel: the literature agents found no published predicted-versus-observed check |
| Genotype rarefaction, N ≤ n | Unbiased if records are independent | needs donor IDs to dedup | trivial | **Y**, after dedup |
| Chao–Jost extrapolation, genotype level | Optimistic beyond 2×; ~1 pt at 2× | have; dedup needed | trivial | **Y within 2–3×**, as a lower bound beyond |
| Good-Turing / Orlitsky unseen-mass | Valid to ~log n | have | trivial | **Y**, to quantify unseen mass, not N\*95 |
| Pitman–Yor / Ewens frequency spectrum | Pessimistic; unidentifiable at K/n ≈ 1 | have | seconds | **With caveat** for Chinese only; **N** for the minority groups |
| Power-law tail (Slater et al. 2015, NMDP) | Assumption-driven | have | moderate | **With caveat**, as a sensitivity bound |
| Dirichlet parametric bootstrap for CIs | Biased low for N\* on a long tail | have | heavy as written | **N**. Use donor-level non-parametric bootstrap with EM re-run |

### 5.3 Recommendation for the manuscript

**Primary analysis**

1. Obtain donor IDs from BMDP and SCBB and deduplicate. Report the true number of donors. Re-run
   EM, uncapped, on unique donors.
2. Use the haplotype–HWE model with a floor below every group's 1/(2n) (1e-6). Report the
   **unconditional** coverage curve C(N). Remove the 1% "other" pooling.
3. Validate at today's size against independent units: all 1,450 SCBB units and, if obtainable,
   genotype-level typing of real patients. Publish the predicted-versus-observed table.
4. Report the coverage **gain from realistic expansion**, two- and three-fold, where model and
   assumption-light estimator agree to within a few points (`publishable_range.log`). Each cell
   shows the 1e-6 model, then Chao–Jost:

| Group | Today (dedup) | 2× | 3× |
|---|---|---|---|
| Chinese | 0.285 / 0.310 | 0.352 / 0.411 | 0.394 / 0.497 |
| Malay | 0.152 / 0.163 | 0.206 / 0.244 | 0.243 / 0.316 |
| Indian | 0.057 / 0.071 | 0.087 / 0.118 | 0.110 / 0.163 |
| Others | 0.029 / 0.044 | 0.049 / 0.077 | 0.064 / 0.109 |

5. Report **high targets as bounds, not point estimates**:
   - lower bound from Chao–Jost on deduplicated genotypes;
   - model value from the 1e-6 haplotype model;
   - the extrapolation factor printed beside each.
   The message is: a 95% 10/10 match through same-ethnicity recruitment would need at least several
   hundred thousand to tens of millions of donors per group, far beyond any registry anywhere.

**Sensitivity analyses**, in a supplement:

- the floor sweep of §4.2 under the unconditional estimand;
- the deduplication bounds of §4.1;
- 8/8 and 9/10 alongside 10/10.

**Limitation text for the paper:**

> "Registry targets above roughly 50% coverage extrapolate one to three orders of magnitude beyond
> the observed donor sample. No estimator can be validated at that range, because the result depends
> on the frequency of HLA haplotypes not yet observed in Singapore. We therefore report these
> targets as bounds. We validated the model only at the current registry size, against independently
> banked cord-blood units."

**Reframe the policy message.** High-resolution 10/10 coverage of 95% is not attainable by domestic
recruitment for any CMIO group, consistent with the 41% 8/8 likelihood for Chinese patients in the
US registry (Gragert 2014). Realistic policy levers belong in the discussion: 9/10 and 7/8
acceptance, haploidentical transplantation, cord blood, and international search. The comparative
results that survive (§3) then carry the paper.

---

## 6. Overall assessment and remaining uncertainty

The comparative science is sound. The absolute registry targets in v2.15c are not, and they should
be withdrawn. The 1e-6 model is better supported than the project's own documentation claims, but
its 95% targets are model-dependent extrapolations, not estimates.

**Corrections the project should make**

- The section "Which floor predicts reality better" and "The real argument for 1e-3" in
  `docs/explanation-frequency-floor.md` should be corrected for the contaminated benchmark (§4.1)
  and the mass-versus-mass-squared scaling (§4.2).

**Limits of this review**

- The duplicate diagnosis rests on string identity and row adjacency, not donor IDs.
- The SCBB probes are small outside the Chinese group, and cord donors are a proxy for patients.
- EM was not re-run on deduplicated data.
- The Chao–Jost and Pitman–Yor formulas were implemented from the cited literature without the full
  texts in hand. The literature agent marked the exact Chao–Jost form [uncertain].
- The review's scripts live outside the repository and are not unit-tested.

---

## 7. Reproducing this review

Each script is run with `/home/alvin/miniconda3/bin/python -I`, from `/data/alvin/tmp/hla_review/`.

| Script | Purpose | Log |
|---|---|---|
| `genotype_coverage.py` | Floor-free coverage, extrapolation, held-out test | `run1.log` |
| `floor_sweep.py` | Floors scored, conditional versus unconditional, effect of pooling | `floor_sweep.log` |
| `diagnose_gap.py` | Per-locus HWE; twins by source | `diagnose_gap.log` |
| `twin_structure.py` | Twin adjacency against random expectation | `twin_structure.log` |
| `dup_resolution.py` | Raw-string identity of twins | `dup_resolution.log` |
| `dedup_coverage.py` | Dedup scenarios, SCBB probe validation | `dedup_coverage.log` |
| `model_vs_probe.py` | 1e-6 model against probes; timing | `model_vs_probe.log` |
| `publishable_range.py` | 2× and 3× coverage; model N\* tables | `publishable_range.log` |

---

## 8. References

Retrieved during this review unless marked otherwise.

- Gragert L, et al. HLA match likelihoods for hematopoietic stem-cell grafts in the U.S.
  registry. *N Engl J Med* 2014;371:339–348. doi:10.1056/NEJMsa1311707. Table 1 retrieved.
- Schmidt AH, et al. Toward an optimal global stem cell donor recruitment strategy. *PLoS ONE*
  2014;9:e86605. doi:10.1371/journal.pone.0086605.
- Kollman C, et al. *Transplantation* 2004;78:89–95. PMID 15257044. Abstract only.
- Slater N, Louzoun Y, et al. Power laws for heavy-tailed distributions: modeling allele and
  haplotype diversity for the National Marrow Donor Program. *PLoS Comput Biol* 2015.
  doi:10.1371/journal.pcbi.1004204.
- Good IJ. The population frequencies of species and the estimation of population parameters.
  *Biometrika* 1953;40:237–264. doi:10.1093/biomet/40.3-4.237.
- Chao A, Jost L. Coverage-based rarefaction and extrapolation. *Ecology* 2012;93:2533–2547.
  doi:10.1890/11-1952.1. Full text not read; formula [uncertain].
- Colwell RK, et al. Models and estimators linking individual-based and sample-based rarefaction,
  extrapolation and comparison of assemblages. *J Plant Ecol* 2012;5:3–21. doi:10.1093/jpe/rtr044.
- Chao A, et al. Rarefaction and extrapolation with Hill numbers. *Ecol Monogr* 2014;84:45–67.
  doi:10.1890/13-0133.1.
- Orlitsky A, Suresh AT, Wu Y. Optimal prediction of the number of unseen species. *PNAS*
  2016;113:13283–13288. doi:10.1073/pnas.1607774113.
- Schäfer C, et al. Hapl-o-Mat. *BMC Bioinformatics* 2017. doi:10.1186/s12859-017-1692-y.
- Eberhard HP, et al. *Tissue Antigens* 2013;82:93–105. Abstract only.
- Favaro S, Lijoi A, Mena RH, Prünster I. Bayesian non-parametric inference for species variety
  with a two-parameter Poisson–Dirichlet process prior. *JRSS-B* 2009. [uncertain]: not retrieved;
  this is the source of the posterior formula used.

---

## Addendum (2026-11-02): EM refit on deduplicated BMDP donors only

The live `haplo_freqs_em.csv` was fitted on BMDP and SCBB together, so the §4.3 check against
SCBB units was not strictly out of sample. `dedup_em.py` reruns the 5-locus EM at the 1e-6 floor
on BMDP donors only, after removing duplicates. It uses unconditional coverage with no pooling.
Each EM run takes about 50 s for all four groups (`dedup_em.log`, `dedup_em_results_point.csv`).

| Group | Dedup rule | n | SCBB probe | C(n) | C(2n) | C(3n) | N\*50 | N\*75 | N\*95 |
|---|---|---|---|---|---|---|---|---|---|
| Chinese | nearby | 39,943 | 26.1% | 28.5% | 35.3% | 39.5% | 317,067 | 3,723,956 | 109,111,095 |
| Malay | nearby | 5,179 | 7.8% | 15.2% | 20.7% | 24.3% | 147,020 | 1,291,789 | 19,088,074 |
| Indian | nearby | 5,026 | 4.4% | 5.7% | 8.7% | 10.9% | 605,921 | 3,769,158 | 32,918,769 |
| Others | nearby | 3,265 | 3.2% | 3.0% | 5.0% | 6.6% | 736,934 | 4,283,051 | 30,038,019 |
| Chinese | all repeats | 32,777 | 26.1% | 16.3% | 22.4% | 26.5% | 729,331 | 6,378,635 | 128,690,554 |
| Malay | all repeats | 4,846 | 7.8% | 11.8% | 16.7% | 20.1% | 194,677 | 1,525,935 | 19,217,142 |
| Indian | all repeats | 4,882 | 4.4% | 4.7% | 7.4% | 9.5% | 666,147 | 3,921,127 | 32,722,869 |
| Others | all repeats | 3,220 | 3.2% | 2.6% | 4.5% | 6.0% | 778,232 | 4,363,231 | 29,861,908 |

**Findings**

1. **The leakage did not matter.** The BMDP-only fit reproduces the §4.3 values to within
   0.1 point: Chinese 28.5%, Malay 15.2%, Indian 5.7%, Others 3.0%.
2. **The probes favour removing only nearby duplicates.**
   - Removing every raw repeat drops Chinese to 16.3%, far outside the probe interval (23.4–28.9%).
     That rule also deletes genuine, independent genotype twins.
   - Removing only nearby duplicates (≤10 rows apart) gives 28.5%, inside the interval. Use it as
     the primary rule until donor IDs are available.
3. **Today's size and modest expansion are robust; the targets are not.**
   - C(n) and C(2n) move by at most a few points between the two duplicate rules for Malay, Indian
     and Others.
   - N\*50 moves by up to 2.3× (Chinese), and already lies 8–226× beyond the sample size. Treat
     N\* values as model-based bounds, not estimates.
4. **Malay remains the weak cell.** The model gives 15.2% against a 7.8% probe (n = 90). The
   probe interval is too wide to separate a real bias from chance.

### Uncertainty

**The with-replacement bootstrap is biased.** It resamples donors with replacement, B = 100
(`dedup_em_bootstrap.csv`, `dedup_boot.log`). The full-sample estimate falls outside its own
95% percentile interval in every group: coverage is shifted up, and N\*50 is shifted down by
12–42%. A resample with replacement repeats donors, and the repeats concentrate the EM haplotype
distribution. It is the same artefact as the duplicate records, and it is why the Dirichlet
bootstrap was withdrawn. **Do not use percentile bootstrap intervals for this model.**

**Half-sampling without replacement** (`subsample_em.py`, B = 100, output in
`subsample_em_ci.csv`) draws half of each group's deduplicated donors and reruns EM. With m = n/2,
the half-sample variance approximates the full-sample variance, so each interval is the point
estimate plus the centred half-sample 2.5–97.5% quantiles:

| Group | C(n) | C(2n) | C(3n) | N\*50 | N\*75 |
|---|---|---|---|---|---|
| Chinese | 28.5% (28.2–28.9) | 35.3% (34.9–35.7) | 39.5% (39.1–39.9) | 317k (308–326k) | 3.72M (3.63–3.81M) |
| Malay | 15.2% (14.3–16.2) | 20.7% (19.6–21.8) | 24.3% (23.2–25.5) | 147k (138–157k) | 1.29M (1.23–1.36M) |
| Indian | 5.7% (5.2–6.2) | 8.7% (8.1–9.3) | 10.9% (10.2–11.6) | 606k (581–631k) | 3.77M (3.66–3.87M) |
| Others | 3.0% (2.6–3.4) | 5.0% (4.4–5.5) | 6.6% (5.9–7.2) | 737k (708–771k) | 4.28M (4.18–4.46M) |

**The sampling intervals are not the main uncertainty. Sample-size drift is.** Fitted on half the
donors, the model gives coverage 3–18% higher and N\* values 14–54% lower than when fitted on all
donors (`drift_pct` in `subsample_em_ci.csv`):

| Group | C(n) drift, half vs full | N\*50 drift | N\*75 drift |
|---|---|---|---|
| Chinese | +3.0% | −13.7% | −22.0% |
| Malay | +8.1% | −31.3% | −39.8% |
| Indian | +13.6% | −39.1% | −48.4% |
| Others | +18.3% | −45.8% | −54.3% |

Each doubling of the typed sample uncovers more rare haplotypes and raises N\*. Nothing shows this
has levelled off at the current sample size. It is the same unseen-haplotype bias described in
§4.3, now measured directly. Consequences for the manuscript:

1. **Coverage at today's size and at 2–3× is publishable with intervals**, especially for
   Chinese, where drift is small (3%). Even there, the drift is larger than the sampling interval,
   so say that values are mildly optimistic.
2. **N\* values are lower bounds** that grow with the typed sample. Their sampling intervals
   (±3–5%) understate the real uncertainty by an order of magnitude. Report them as "at least"
   values with the drift stated, never as point estimates with a CI.
3. The drift is largest where the samples are smallest (Others, Indian). Within this method, that
   ranks the cross-group comparison of N\* as the least reliable output.

## Addendum B (2026-11-02): sample-size projection and calibrated policy figures

**Purpose.** A registry-size figure used for policy must come from a projection whose error has
been measured. `sample_size_projection.py` refits the 1e-6 EM on 1/8, 1/4, 1/2 and all of the
near-deduplicated BMDP donors (12 reps per fraction, without replacement). Each fit is scored
against the SCBB probes. Fitting on n/k donors and predicting coverage at n tests a k-fold
extrapolation out of sample. That is the same step as projecting today's registry to k times its
size. Outputs: `sample_size_projection.{csv,log}`, `calibrated_projection.{py,csv,log}`.

The probe truth here comes from `hla_clean.csv` (Chinese 25.4%, n = 1,028). It differs slightly
from §4.1 (26.1%, n = 1,000), which used the regex-filtered Excel sheet.

### B.1 The model's optimism is systematic and grows with the extrapolation factor

Error is model minus probe truth at full BMDP size, in percentage points (median of 12 reps):

| Group | Fit at own size (1×) | 2× | 4× | 8× |
|---|---|---|---|---|
| Chinese | +3.1 | +4.0 | +5.4 | +7.8 |
| Malay | +7.8 | +8.9 | +11.2 | +15.7 |
| Indian | +1.4 | +2.2 | +3.9 | +6.1 |
| Others | +0.0 | +0.6 | +1.5 | +3.1 |

- The model is **always optimistic**, and the optimism grows by about 1–2 points per doubling of
  the extrapolation factor.
- Error at the model's own size is roughly scale-invariant. For Chinese it is +4.4, +3.8, +3.7 and
  +3.1 at n/8, n/4, n/2 and n. That supports correcting a k-fold projection by the error measured
  for the same k. The slight decline makes the correction mildly conservative.
- Malay carries the largest bias, but its probe is only 7 of 94 units.

### B.2 Calibrated coverage projections (publishable)

Corrected C(kn) = model C(kn) − measured error at factor k. The 3× error is interpolated in log k.
The 95% interval combines the probe Clopper–Pearson interval with the replicate spread:

| Group | Today (n) | 2× | 3× |
|---|---|---|---|
| Chinese | 25.4% (22.8–28.2) · n = 39,943 | **31.3%** (28.3–34.5) | **34.7%** (31.5–38.0) |
| Malay | 7.4% (3.0–14.7) · n = 5,179 | **11.7%** (6.4–19.9) | **14.0%** (8.7–22.3) |
| Indian | 4.3% (1.2–10.6) · n = 5,026 | **6.5%** (2.9–13.2) | **7.7%** (3.8–14.9) |
| Others | 3.0% (1.2–6.0) · n = 3,265 | **4.4%** (2.2–7.9) | **5.4%** (3.0–9.1) |

An independent check comes from scoring the probes against nested BMDP subsamples, with no model
at all. Observed Chinese coverage at n/8, n/4, n/2 and n is 9.8%, 14.4%, 19.3% and 25.4%, a gain of
about 5–6 points per doubling. One more doubling therefore gives about 31%, which matches the
calibrated 31.3%.

### B.3 Registry-size targets: only short-range targets can be calibrated

**N\* does not converge with sample size.** The elasticity d log N\* / d log m, between successive
doublings of the typed sample, measures how much the estimated target grows when twice as many
donors are typed:

| Group | N\*50 | N\*75 | N\*95 |
|---|---|---|---|
| Chinese | 0.46 → 0.32 → 0.21 | 0.66 → 0.49 → 0.36 | 1.05 → 0.87 → 0.68 |
| Malay | 0.85 → 0.69 → 0.53 | 1.02 → 0.91 → 0.73 | 1.07 → 1.25 → 1.15 |
| Indian | 0.92 → 0.91 → 0.73 | 1.13 → 1.10 → 0.97 | 0.75 → 1.24 → 1.30 |
| Others | 0.87 → 0.92 → 0.87 | 0.71 → 1.01 → 1.13 | 0.61 → 0.31 → 1.06 |

An elasticity near 1 means the estimated target roughly doubles every time the typed sample
doubles. For those cells the data contain no information about the target's size.
- Only Chinese N\*50 shows clear deceleration.
- Every N\*95, and every minority-group target, is still growing proportionally.
- **These targets cannot be estimated from these data by any method that relies on the observed
  haplotype spectrum**, because the answer depends on haplotypes not yet sampled.

**A calibrated target exists only where the required extrapolation is ≤8×**, the largest factor
validated in B.1. Solving model C(N) = target + err(N/n):
- **Chinese, 40% 10/10 coverage: N ≈ 236,000** (5.9× today; the error used is +6.8 points). The
  ±2.7-point probe interval and the observed gain of about 6 points per doubling widen this to
  roughly 170,000–320,000. The model-free doubling curve agrees: 25 → 31 → 37 → 43% at 40k → 80k
  → 160k → 320k.
- Chinese 50%: the calibrated N is about 930,000 (23×). That lies beyond the validated range,
  because the error has to be extrapolated. Report it only as "of the order of one million, not
  validated".
- No group reaches a calibrated 50% within 8×. At 8× today's size, the calibrated coverage is
  Chinese 42.3% (about 320,000 donors), Malay 18.9% (about 41,000), Indian 11.9% (about 40,000)
  and Others 8.8% (about 26,000).

### B.4 What this supports as policy guidance

1. **Publish calibrated coverage at today's size, 2× and 3×** (B.2), with intervals. These are
   validated projections.
2. **Publish one calibrated target: Chinese 40% at about 236,000 donors (170,000–320,000).**
   This is the highest coverage level any group can reach within the validated range.
3. **State the high targets qualitatively.** 75–95% 10/10 coverage cannot be reached by domestic
   recruitment at any size these data can validate. The models put it in the millions to hundreds
   of millions, and those figures grow with every doubling of the typed sample. This agrees with
   Gragert et al. (2014): 41% 8/8 for Chinese patients across the whole US registry.
4. **Marginal yield is the policy-relevant quantity.** Each doubling of the registry adds about
   6 points for Chinese, 2–4 points for Malay and Indian, and 1–2 points for Others. A minority
   donor recruited adds less absolute coverage than a Chinese donor. That argues for targeted
   recruitment together with 9/10, haploidentical, cord blood and international search, rather
   than a single headline donor count.

*Assumptions:*
- The extrapolation error depends on the extrapolation factor, not on absolute scale. This is
  checked for 1–8× only.
- SCBB cord units are representative of patients.
- The nearby-duplicate rule is close to the true duplicate set. Donor IDs would settle this.
