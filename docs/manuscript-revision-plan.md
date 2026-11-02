# Manuscript revision plan: from "registry size" to "calibrated coverage"

Status: proposal, 2026-11-02. It follows from `docs/review-registry-size-estimation.md`
(Addenda A and B).

## Verdict: rewrite the registry half; keep the verification half

Patching v2.15c is not viable. Its central result is the table of N\* for 75–95% coverage. That
table changes estimand, changes by about 2,000×, and can no longer be stated as a point estimate.
A tracked-changes pass would turn almost every number in the abstract, results and discussion red,
and leftover literals are the documented failure mode of this repo's builders.

| Part of v2.15c | Fate |
|---|---|
| Allele and haplotype frequency verification of the 2022 paper | **Keep.** Unaffected by the floor or the duplicates. Re-check allele frequencies after dedup: duplicates over-weight some genotypes. |
| HWE, LD and Others stratification | Keep. Re-run on deduplicated donors. |
| Registry-size targets (N\* at 75–95%) and the 1e-3 snapshot | **Replace.** |
| Cross-group ranking of N\* | **Drop.** It inverts between floors, and minority-group N\* does not converge. |
| Bootstrap CIs (Dirichlet) | Replace with half-sampling intervals. |
| Policy message | Reframe from "recruit N donors" to marginal yield plus alternative-donor strategies. |

## New primary outcomes

1. **Coverage today, validated.** Predicted versus observed against independent SCBB cord units.
   This is the paper's calibration table and, as far as the literature agents found, a novel one.
2. **Calibrated coverage at 2× and 3×**, bias-corrected by the measured extrapolation error,
   with intervals (Addendum B.2).
3. **Marginal yield per doubling, by group.** Model-based, with a model-free check (B.2).
4. **One calibrated target:** Chinese 40% 10/10 at about 236,000 donors (170–320k).
5. **High targets stated as non-estimable**, together with the sample-size drift evidence (B.3).

## External triangulation to include

All of these are below 50%, which agrees with the BMDP experience that a ~150k-donor registry does
not find a match for half of Chinese patients:

| Source | Population | Registry | Match level | Likelihood |
|---|---|---|---|---|
| This work (calibrated) | Singapore Chinese | ~150k total (~100–120k Chinese) | 10/10 (A, B, C, DRB1, DQB1) | ~33–35% |
| BMDP national survey 2024 (press) | Singapore | current BMDP | "perfect match" | Chinese 40%, Malay 20%, Indian 2% |
| Hong Kong BMDR (Front Immunol 2021) | HK Chinese | 130,000 | 8/8 / 10/10 incl. DPB1 | 49% / 20% (27% at 250k) |
| Gragert et al. 2014 (NEJM) | US Chinese | whole NMDP | 8/8 | 41% |

Their definitions of a match differ, so they bound the result rather than test it. Note that 8/8
must exceed 10/10, so HK's 49% at 8/8 is consistent with roughly 35% at 10/10.

## WMDA position

No WMDA document found sets a numeric registry size. The WMDA and related literature frame
optimal size as population-specific, with diminishing returns. They prioritise donor age, male
donors, availability and diversity through targeted recruitment over raw size, and they treat
international search as the norm for small, diverse populations. (BMDP 2023: 2 in 3 patients
received cells from overseas donors.) The paper should cite this framing rather than a target
number.

## Blocking inputs before drafting

1. **BMDP donor IDs.** They resolve the duplicate count, the only assumption that moves C(n) by
   more than the sampling interval.
2. **What the ~59k typed extract is relative to the ~150k registry.** High-resolution subset? A
   time window? If the live search runs over 150k donors, today's coverage must be stated at
   150k, not 59k.
3. **BMDP's own match statistic and its definition** (10/10, 8/8, available, or searched),
   so it can be used as a second validation point rather than an anecdote.

## Build plan (when inputs arrive)

- Freeze a new snapshot, e.g. `analysis/snapshot_dedup_1e-6/`, from the review scripts promoted
  into `analysis/`. Keep `snapshot_1e-3/` frozen as the record of v2.15c.
- Write a new builder for v3.0. Every number is derived from snapshot CSVs, with no literals.
- Bump `VERSION.md`. The red-mark diff against v2.15 is meaningless for a rewrite, so mark
  against the plan instead.
