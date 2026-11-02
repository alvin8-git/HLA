# Why the frequency floor decides the answer

Every registry-size number this project reports depends on a preprocessing
parameter that never appears in the headline: the **rare-haplotype frequency
floor**. Haplotypes estimated below the floor are dropped before the coverage
model runs.

It looks like a housekeeping choice. It is not. Move the floor three orders of
magnitude and the Chinese 95% 10/10 target moves from roughly forty thousand
donors to roughly eighty-seven million, on the same data with the same code.
This document explains why, what each choice buys, and which conclusions
survive the change.

**Current status:** the submission manuscript (v2.15c) uses a floor of
**1×10⁻³**. A methodological re-run at **1×10⁻⁶** exists under the v2.17 line
and lives in `analysis/data/`. Neither supersedes the other — they answer
different questions, and this document is the map between them.

---

## The problem

The coverage model needs a frequency for every haplotype a patient might
carry. EM haplotype phasing produces one, but the tail is enormous: at full
resolution the Chinese sample yields hundreds of thousands of distinct
haplotypes, most of them estimated from a single ambiguous genotype.

Two things go wrong if you keep all of them:

1. **Cost.** The model expands K haplotypes into K(K+1)/2 diplotypes, so work
   grows with the square of the tail. See [Compute](#compute) below.
2. **Credibility.** Below a certain frequency, an "estimate" is not an
   observation. It is EM redistributing probability mass across phasings that
   the data cannot distinguish.

So a floor is necessary. The question is where to put it.

---

## The 1/(2n) rule

The natural boundary is the **singleton frequency**, 1/(2n), where n is the
number of individuals EM saw. A haplotype at that frequency corresponds to one
observed copy. Below it, nothing was counted.

```
                  1/(2n)
                    |
  EM noise          |          genuinely resolved haplotypes
  <-----------------+------------------------------------------>
  cutting here      |          cutting here destroys data
  is harmless       |
```

Per-group singleton frequencies, from the 5-locus donor counts:

| Group | n (5-locus) | 1/(2n) |
|---|---|---|
| Chinese | 44,400 | 1.13×10⁻⁵ |
| Malay | 5,578 | 8.96×10⁻⁵ |
| Indian | 5,490 | 9.11×10⁻⁵ |
| Others | 3,767 | 1.33×10⁻⁴ |

Both floors in use sit on the wrong side of that line, in opposite directions:

- **1×10⁻⁶** is *below* every group's singleton frequency. It keeps material
  the sample cannot resolve.
- **1×10⁻³** is *above* every group's singleton frequency, by 8× to 88×. It
  discards haplotypes the sample genuinely resolved.

Neither is the principled choice. The principled choice is a per-group floor at
roughly 1/(2n), which nothing in this repository currently computes.

---

## The two floors side by side

Haplotypes surviving the floor, and the resulting 95% 10/10 same-ethnicity
target:

| | K @ 1e-3 | N* @ 1e-3 | K @ 1e-6 | N* @ 1e-6 |
|---|---|---|---|---|
| Chinese | 140 | 42,871 | 9,574 | 87,384,114 |
| Malay | 137 | 41,779 | 3,134 | 16,552,048 |
| Indian | 144 | 44,863 | 4,024 | 28,817,950 |
| Others | 123 | 32,360 | 4,634 | 26,762,600 |

Read the K columns first, because they explain everything downstream.

At 1e-3 all four groups retain **123 to 144** haplotypes — nearly identical,
regardless of whether EM saw 3,767 people or 44,400. The floor is high enough
that it truncates every group to the same common core. That imposed equality is
what makes the four N* values land within 1.39× of each other.

At 1e-6 the floor stops binding and each group keeps as many haplotypes as its
sample can resolve. Chinese keeps 3× more than Malay, and the N* values spread
to **5.28×**.

**The ordering also inverts.** At 1e-3, Indian needs the largest registry and
Others the smallest. At 1e-6, Chinese needs the largest and Malay the smallest.
Any statement of the form "group X needs the biggest registry" is a statement
about the floor, not about the population.

---

## Which floor predicts reality better

> **Corrected 2026-11-02.** An earlier version of this section scored both floors against the
> in-registry twin rate and concluded that 1e-3 was closer to reality. That benchmark is
> inflated by duplicate donor records, and the 1e-3 figures multiplied by retained mass instead
> of mass squared. Details: `docs/review-registry-size-estimation.md` §4.1–4.3.

The quantity a model must reproduce is how often a patient finds an exact 10/10 genotype twin
among the donors the registry holds *today*. Two benchmarks exist.

**The in-registry twin rate is contaminated.** BMDP contains duplicate records: 4,370 records
repeat the raw typing string of a record at most 10 rows earlier, and every one of those nearby
pairs is identical character for character. A duplicate is its own twin, so the as-ingested
rates (Chinese 41.8%, Malay 25.3%, Indian 19.7%, Others 18.1%) overstate coverage.

**The SCBB cord units are clean.** They are typed and banked independently of BMDP. The share
of cord units with a 10/10 twin in BMDP is an uncontaminated estimate of coverage at BMDP's
size.

| | SCBB probe (95% CI) | 1e-3, unconditional | 1e-3, as reported (conditional) | 1e-6, BMDP-only fit, dedup |
|---|---|---|---|---|
| Chinese | **26.1%** (23.4–28.9) | 19.5% | 94.9% | 28.5% |
| Malay | **7.8%** (3.2–15.4) | 13.5% | 73.2% | 15.2% |
| Indian | **4.4%** (1.2–10.9) | 4.7% | 63.5% | 5.7% |
| Others | **3.2%** (1.3–6.4) | 2.6% | 61.1% | 3.0% |

The 1e-3 *unconditional* column is conditional coverage times retained mass **squared**. A
patient needs both haplotypes to clear the floor, so that column is what the 1e-3 model claims
about all patients. Both 1e-3 columns are evaluated at the as-ingested size, which flatters
them slightly. The 1e-6 column comes from a fit on deduplicated BMDP donors only, so the
cord units are out of sample.

How to read the table:

- **Today's coverage cannot tell the floors apart.** Both unconditional columns sit close to the
  probes.
- **The conditional figures the manuscript reports are far from reality.** For example, 94.9% for
  Chinese against a measured 26.1%.
- **The floors split when you extrapolate.** At 1e-3 unconditional coverage can never exceed
  retained mass squared, 13–27% depending on the group, so every target above that is
  unreachable by construction. The renormalised version reaches 95% only by dropping the
  patients who carry rare haplotypes.
- **1e-6 is the only floor that is both calibrated today and able to extrapolate.**

Scripts and logs: `/data/alvin/tmp/hla_review/` (`dedup_coverage.py`, `floor_sweep.py`,
`dedup_em.py`).

---

## Why bootstrap CIs stop meaning anything at 1e-6

The bootstrap is a Dirichlet parametric resample: concentration
`α = n_eff × f̂`, then recompute N* for each of 1,000 replicates. Two things
break when the floor drops.

**Most of the input stops being data.**

| | K | below 1/(2n) | α clipped to the 0.1 constant |
|---|---|---|---|
| Chinese | 9,574 | 1,571 (16%) | 130 (1.4%) |
| Malay | 3,134 | 602 (19%) | 356 (11%) |
| Indian | 4,024 | 728 (18%) | 496 (12%) |
| Others | 4,634 | 1,672 (36%) | 1,471 (32%) |

`alpha = np.maximum(freqs_norm * n_eff, 0.1)` in `analysis/09_bootstrap_ci.py`
replaces any concentration below 0.1 with that constant. For a third of the
Others haplotypes, the resampling spread is set by a hardcoded number rather
than by the sample. At 1e-3, nothing is clipped and the minimum α is 10.5.

**The bias is systematic, not random.** Because
`E[f²] = f² + f(1−f)/(n+1)`, and diplotype frequencies are built from `f²` and
`2fᵢfⱼ`, every replicate inflates them — overstating coverage, understating
N*. The inflation is roughly 525% for haplotypes in the 1e-6 to 1e-5 band and
134% in the 1e-5 to 1e-4 band, so a lower floor admits more of exactly the
material that drives the bias.

The fingerprint is visible in the output. In the 1e-6 run, Chinese 8/8 returned
a median of 52,266,536 against an EM estimate of 87,384,114, with a 95%
interval of 50.8M–53.9M and `pct_below = 1.000`. Every replicate fell below the
point estimate, and the interval never contains it. A confidence interval that
excludes the thing it brackets is not measuring sampling error.

**What the CI does and does not cover, at either floor.** The bootstrap
resamples EM's *output* and never re-runs EM, so it is blind to phase ambiguity
— the very reason EM exists. Ranked by how much they move N*:

| source | magnitude | in the CI? |
|---|---|---|
| frequency floor | 2,098× | no — reported separately |
| EM input cap | 8.2% at 1e-3, 264% at 1e-4 | no — reported separately |
| phase ambiguity | ~11 configurations per donor | **no, and not reported** |
| multinomial sampling | ±2–3% | yes — this is the CI |

The published interval is the smallest term on the list. Table 1 of the
manuscript says so: the CIs "reflect haplotype-frequency sampling variability
only and are a lower bound on total uncertainty."

---

## Compute

Measured wall-clock from the 2026-08-20 overnight chain at floor 1e-6
(`paper_BMT_workdir/overnight.log`):

| step | wall clock |
|---|---|
| `09_bootstrap_ci.py` | **9h 12m, then OOM-killed** (~1 of 8 units done) |
| `06_partial_match_plots.py` | 4h 16m |
| `04_registry_model.py` | 36m |
| `03` EM haplotype phasing | **59 seconds** |

EM — the algorithmically hard step — is the cheapest thing in the pipeline.
The cost is the O(K²) diplotype expansion downstream, repeated B=1,000 times by
the bootstrap:

| floor | Chinese K | diplotypes | bootstrap evaluations |
|---|---|---|---|
| 1e-3 | 140 | 9,870 | 3.7 × 10⁹ |
| 1e-6 | 9,574 | 45,835,525 | 7.0 × 10¹² |

**A corollary that matters for planning:** at a fixed floor, K *saturates* —
Chinese EM on 5,000 donors yields 143 haplotypes, on 45,018 it yields 136. So
growing the input dataset from 59k to 500k donors costs about seven extra
minutes in EM and changes nothing downstream. Lowering the floor on the same
59k donors costs roughly three days and more RAM than the machine has.

Compute is a genuine constraint, but it is a *consequence* of the floor choice,
never a justification for it. Do not defend 1e-3 on runtime grounds — a
reviewer will not accept it, and it is not the real argument.

---

## The real argument for 1e-3

Two things, neither of which is speed:

1. **It defines a stated estimand.** The manuscript reports coverage
   *conditional* on both of a patient's haplotypes clearing the floor, and says
   so. That is a well-defined quantity: the registry size needed to serve the
   ~52% of Chinese patients (and 36% of Others patients) who carry common
   haplotypes.
2. ~~It is the only version that has been validated.~~ *Withdrawn 2026-11-02.* The "within
   7 points" check used a duplicate-inflated benchmark and multiplied by retained mass instead
   of mass squared. Against independent SCBB cord units, the 1e-6 model is calibrated at
   today's size (see above). Spearman r = 0.70 tests the rank order of haplotypes, not the
   level of coverage.

---

## What survives a change of floor

**Survives — every claim that is a ratio or a direction:**

- Same-ethnicity matching beats cross-ethnic by a wide margin (1.7× for
  Chinese, 33–309× for the other three groups). More extreme at 1e-6, not less.
- Relaxing 10/10 to 9/10 roughly halves the requirement.
- 8/8 barely improves on 10/10, because DRB1 and DQB1 are in strong linkage
  disequilibrium (D′ = 0.94–0.99 at 1e-3, 0.90–0.96 at 1e-6). Malay is the
  exception at 9.6%.
- "Others" is not one population. It fails Hardy–Weinberg at all five loci,
  which is a property of population structure and independent of the floor.

**Does not survive:**

- **The 40,000–45,000 headline.** It moves 2,098× between floors.
- **The cross-group ranking.** It fully inverts, as shown above.
- **The demographic-robustness result.** The four patient-mix scenarios span
  3.4% at 1e-3 and **67.1%** at 1e-6. At 1e-3 the outlier group (Others) carries
  3.2% of the weight and cannot move the weighted mean; at 1e-6 the outlier is
  Chinese, carrying 74.3%, and dominates it. That flatness is a property of the
  floor, not of haplotype structure.

The safe way to read this project: its **comparative** conclusions are robust,
its **absolute targets and rankings** are conditional on a preprocessing choice
that should be stated wherever they are quoted.

---

## Reproducing either floor

The floor is `FREQ_THRESHOLD` in `analysis/hwe_test.py`. Two scripts,
`11_others_stratification.py` and `15_em_convergence.py`, have historically
held private copies of it — grep before changing it:

```bash
grep -rn "FREQ_THRESHOLD" analysis/
```

The frozen 1e-3 results the manuscript reports live in
`analysis/snapshot_1e-3/data/`. The 1e-6 re-run lives in the live
`analysis/data/`. Anything that reads the live directory silently gets the
1e-6 numbers, which is why every manuscript build points at the snapshot.

---

## Related

- [README](../README.md) — pipeline overview and how to run it
- [Technical Documentation](../Documentation.md) — full method reference,
  §5 (haplotype estimation) and §6 (registry model)
- [VERSION.md](../VERSION.md) — the change log, including commit `dbdd675`
  where the 1/(2n) analysis was first recorded
- `analysis/09_bootstrap_ci.py` — the Dirichlet bootstrap discussed above
- `analysis/15_em_convergence.py` — the EM saturation data behind the
  K-saturates argument
