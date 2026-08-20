# Version History

## v2.6.0 — 2026-08-21 (v2.16 docx) — PIPELINE RE-RUN AT freq_threshold=1e-6

**Final configuration:** `freq_threshold=1e-6`, `cap=50000`, search ceiling 1e10,
EM maximum-likelihood point estimates, no confidence intervals. Retains 100% of
haplotype frequency mass in all four groups (3,134–9,574 haplotypes each). The
headline numbers are in "Final v2.16 figures" below.

This entry records the route as well as the destination, because two intermediate
configurations (1e-4 capped, then 1e-4 uncapped) produced numbers that circulated
during the work and are **superseded**. Sections marked SUPERSEDED are kept so the
reasoning is auditable, not because their figures stand.

### SUPERSEDED step 1 — haplotype retention floor 0.1% → 0.01% (`hwe_test.run_em_haplotypes`)
The 0.1% floor retained only 123–144 haplotypes per group, representing **36–53%
of total haplotype frequency mass**, which was then renormalised to 1.0. Because
`C(N)=Σ F·[1−(1−F)^N]` converges slowly precisely where F is small, the discarded
tail is what determines high-coverage behaviour.

At 0.01% the pipeline retains 2,310–3,035 haplotypes per group and **97.2–97.9%**
of frequency mass. Headline consequences (10/10, same-ethnicity):

| Group | N* at 95%, v2.15 (1e-3) | N* at 1e-4, capped (superseded) |
|---|---|---|
| Chinese | 42,847 | 12,001,379 |
| Malay | 40,032 | 11,591,997 |
| Indian | 43,855 | 20,877,121 |
| Others | 31,181 | 21,989,663 |

95% coverage is therefore **not attainable** by any national registry. The
manuscript is reframed around coverage attainable at feasible size: at 50,000
same-ethnicity donors, 10/10 coverage is 38.9% (Chinese), 40.7% (Malay), 23.1%
(Indian), 18.7% (Others).

Cross-ethnic matching from a Singapore-weighted pool, with the order fix below
applied, is now quantified rather than reported as ">10⁷": 18.1M for Chinese
(1.5× the same-ethnicity figure), 757M for Malay (65×), 1.88bn for Indian (90×),
and beyond the 10bn search ceiling for Others. The single-shared-registry
("Combined pooled") model requires 73.9M at 95%.

Bootstrap CIs are computed at B=1,000 (unchanged from prior releases) via
`run_overnight.sh`, which chains EM → downstream analyses → bootstrap → rebuild.

### SUPERSEDED step 2 — EM input cap 5,000 → 50,000 (must accompany the floor change)
`15_em_convergence.py`, re-run at the 1e-4 floor, showed the 5,000-individual cap
inflates the Chinese N* by **264%** (11,487,962 capped vs **3,153,571** at the
full 45,018 sample). A 5,000-individual EM cannot resolve phase in the rare tail
and retains spurious low-frequency haplotypes, which the coverage model reads as
real diversity. The curve is non-monotonic and peaks at the cap:

| EM sample size | haplotypes | N* at 95% |
|---|---|---|
| 500 | 566 | 596,303 |
| 5,000 (old cap) | 2,309 | 11,487,962 |
| 15,000 | 1,319 | 3,501,491 |
| 45,018 (full) | 1,253 | 3,153,571 |

At the 1e-3 floor the same cap cost only ~8% and in the conservative direction —
so **the floor and the input cap cannot be chosen independently**. Lowering the
floor without removing the cap converts a mild conservative bias into a large
anti-conservative one. The cap is now 50,000 and binds for no CMIO group.

**Consequence:** the 12.0M/11.6M/20.9M/22.0M figures produced by the capped 1e-4
run are themselves ~3.6× too high; the uncapped re-run supersedes them. The
qualitative conclusion (95% coverage unattainable domestically) is unaffected.

### WITHDRAWN — bootstrap confidence intervals
The Dirichlet parametric bootstrap is **biased downward for N\*** on a
long-tailed haplotype distribution, not merely imprecise. Under a Dirichlet draw
E[f²] = f² + f(1−f)/(n+1), so resampling inflates the squared and product terms
that form diplotype frequencies — by ~525% for haplotypes at 1e-6–1e-5, ~134% at
1e-5–1e-4, 12% at 1e-4–1e-3, 1% above. Since the rare tail governs high-coverage
behaviour, every replicate overstates coverage and understates N\*.

The symptom was present and misread in earlier releases: in v2.15, **18 of 32
rows had the EM estimate outside its own CI**, always above the upper bound,
with pct_below 0.98–1.00. Switching the reported point estimate to the bootstrap
median made the tables self-consistent without addressing the cause. At the 1e-6
floor the gap is unmistakable — Chinese 10/10 at 95%: EM 87,384,114 vs CI
50,772,223–53,870,680, all 1,000 replicates below.

v2.16 therefore reports **EM maximum-likelihood point estimates with no
intervals**, and §2.4 states plainly that uncertainty is real but unquantified.
Figure 2 is now the coverage curve rather than the CI forest plot. Restoring
intervals needs a resampling scheme that preserves rare-tail structure, or
analytic propagation through C(N); neither is attempted here. The bootstrap was
also impractical at this floor — 1 of 8 combinations in 9 hours (~72h projected),
killed (exit 137) on the second.

Post-mortem of the killed run also caught a bug in the parallelised
`_boot_one` (introduced 2026-08-20): it ignored `match_level`, so 8of8
replicates would have been computed on uncollapsed 5-locus haplotypes (log
signature: K=9574 for both levels). Fixed 2026-08-21 via a 4-locus
`collapse_idx` + `np.bincount`, verified identical to the serial path; no 8of8
CI was ever published from the broken path.

### Final v2.16 figures (10/10, same-ethnicity, EM point estimates)

| Group | 75% | 90% | 95% |
|---|---|---|---|
| Chinese | 3,148,792 | 26,222,315 | 87,384,114 |
| Malay | 1,121,822 | 6,716,756 | 16,552,048 |
| Indian | 3,319,587 | 13,826,361 | 28,817,950 |
| Others | 3,924,269 | 14,229,289 | 26,762,600 |

Cross-group comparison of these values remains unsafe (see below): they are
sampled to different depths.

### RESOLVED — the floor is safe below 1/(2n) and destructive above it
The full-sample floor sweep (`paper_BMT_workdir/floor_curve_full.py`) settles the
open question below. The harmful threshold is not a fixed frequency; it is the
frequency of a **singleton haplotype**, 1/(2n), and it moves with sample size.

Chinese (n=45,754, singleton 1.09e-5) — N* at 95%, and inflation vs unfloored:

| Floor | vs 1/(2n) | Haps | Mass | N* 95% | Inflation |
|---|---|---|---|---|---|
| none | — | 234,568 | 100% | 87,530,956 | 1.0× |
| 1e-6 | below | 9,595 | 100.0% | 86,971,552 | 1.0× |
| 1e-5 | below | 8,537 | 99.3% | 76,579,448 | 1.1× |
| 3e-5 | **above** | 3,198 | 91.3% | 16,405,166 | 5.3× |
| 1e-4 | **above** | 1,253 | 80.8% | 3,153,571 | 27.8× |
| 1e-3 | **above** | 136 | 49.0% | 41,727 | **2,097.7×** |

Others (n=3,941, singleton 1.27e-4) reproduces the rule with the break displaced
to a higher floor: flat to 1e-4 (1.2×), collapsing by 1e-3 (791.2×).

Two consequences:

1. **The sub-singleton tail is inert.** Dropping the 224,973 Chinese haplotypes
   below 1e-6 — 95.9% of all distinct haplotypes — costs 0.03% of mass and leaves
   N* within 0.6%. So the tail *is* partly EM phase-ambiguity noise, but that part
   does not matter. The damage above 1/(2n) is done to haplotypes the sample
   genuinely resolves.
2. **A single floor biases unequally-sampled groups unequally.** At the 1e-4 floor
   this release uses: Chinese sits 9.2× above its singleton (27.8× error), Malay
   1.2× above, Indian 1.1× above, Others 0.8× (i.e. below — 1.2× error). Registry
   sizes are therefore **not comparable across CMIO groups** at a common floor,
   and cross-group rankings in this analysis and its predecessors are unsafe.

**Adopted as the final configuration:** `freq_threshold = 1e-6` (below every group's
singleton; 4,654–9,595 haplotypes per group, computationally tractable), which
would give Chinese ~87.0M and Others ~27.0M at 95%. The present release's 1e-4
figures are lower bounds, most severely for Chinese.

### SUPERSEDED NOTE — the floor and the sample size are coupled
Retained frequency mass at a 1e-4 floor, before vs after removing the cap:

| Group | capped (n=5,000) | uncapped (full n) |
|---|---|---|
| Chinese | 2,356 haps, 97.2% | 1,257 haps, **80.8%** |
| Malay | 2,310 haps, 97.4% | 1,256 haps, 86.6% |
| Indian | 3,013 haps, 97.3% | 1,609 haps, 82.7% |
| Others | 3,035 haps, 97.9% | 3,035 haps, 97.9% (cap never bound) |

Removing the cap *lowers* retained mass. The reason is that a fixed frequency
floor means different things at different sample sizes: at n=5,000 (10,000
chromosomes) a haplotype seen once has frequency 1e-4 and survives the floor, so
the floor barely bites and retains sampling noise. At n=44,400 (88,800
chromosomes) a singleton sits at ~1.1e-5, well below the floor, so genuinely rare
haplotypes are now correctly distinguished from noise — but are also discarded.

A floor should therefore scale roughly as 1/(2n): ~1e-5 for the full Chinese
sample rather than 1e-4. The present release uses (full sample, 1e-4), which
retains 81–98% of mass — a large improvement on the 36–53% of v2.15 — but a
1e-5 floor is the logical next refinement and has **not** been run. Estimates
here should be read as lower bounds on N* for that reason.

A floor sweep shows the damage is a **threshold, not a gradient**: 0 → 1e-4
discards 96% of distinct haplotypes but only 3.1% of mass (coverage at 50k moves
36.8% → 39.5%), whereas 1e-4 → 1e-3 discards a further 44% of mass (39.5% →
95.6%). Removing everything below 1e-5 moves coverage 0.2 points, so the effect
is not EM phase-enumeration noise.

### Fixed — order-sensitive cross-ethnic merge (`registry_model.get_diplotype_frequencies`)
Diplotype pairs were labelled `(haplotype1, haplotype2)` in each population's own
frequency-rank order, so the same unordered pair was stored `(X,Y)` in one frame
and `(Y,X)` in another. `04_registry_model.compute_coverage_cross` merges on those
two columns, so mismatched orderings silently scored `donor_freq = 0`. Measured on
Malay-vs-combined: 62% of patient pairs unmatched, **30.4% of patient frequency
mass wrongly zeroed**; cross-ethnic coverage at N=1e6 rose 0.6353 → 0.7757 once
corrected. Pairs are now labelled in canonical lexicographic order.

### Fixed — Dirichlet bootstrap sample sizes exceeded the study total
`N_EFF` was `{45754, 5868, 5586, 3941}`, summing to 61,149 against a stated total
of 59,186. Corrected to the counts reproducible from `hla_clean.csv` for
five-locus-complete individuals: `{44400, 5578, 5490, 3767}` (= 59,235).

### Changed — search ceiling and sweep range
`find_registry_size` n_max 1e7 → 1e10 and `N_SWEEP` 1e3–1e7 → 1e3–1e9, because at
the corrected floor every 95% cell otherwise reported the censored ceiling value
rather than an estimate.

### Performance (no change to results; all verified identical)
- `get_diplotype_frequencies` vectorised over the upper triangle (was a Python
  double loop) — required at 2.3–4.6M pairs.
- `compute_coverage_cross` now aligns patient/donor vectors once and caches, rather
  than re-merging millions of rows inside every binary-search step.
- Added `diplotype_freq_vector` / `find_registry_size_vec` numeric fast paths, and
  parallelised the bootstrap replicate loop (identical rng draw order preserved).

### Known limitation of this release
`06_partial_match_plots.py` is O(N_diplo²) and cannot run at 1e-4 (≈8e12
comparisons). Partial-match results are recomputed by Monte-Carlo sampling of
patient diplotypes instead; Figures 3–4 are otherwise carried over and are flagged
in-text as computed at the previous floor.

## v2.5.0 — 2026-08-20

### Fixed — silhouette discrepancy (v2.15 docx)
- **Silhouette 0.97 was a hand-typed prose error** (introduced at v2.1), never a
  computed value. `11_others_stratification.py` computes s=0.24 at k=3 on the
  five-PC clustering space (0.43 in the PC1–PC2 projection); Figure 7's title
  was correct all along. §3.7, Figure 7 caption, and Limitations now cite 0.24
  and no longer claim "well-separated" on the strength of the phantom 0.97.
- `11_others_stratification.py` now writes `data/others_cluster_silhouette.csv`
  (silhouette per k) so the manuscript cites a traceable computed value, and the
  figure title reports `sil_dict[best_k]` explicitly.
- **§3.7/Figure 7 donor count**: 3,941 → 3,847. The clustering uses only Others
  donors with all five loci typed (cluster sizes 1,029+1,257+1,561 = 3,847);
  3,941 remains correct in §2.4 as the per-ethnicity bootstrap donor count.

## v2.4.0 — 2026-06-18

### Fixed — peer-review correctness pass (text-only, no recomputation)
- **C1 number drift**: removed EM-MLE values that leaked into prose where the
  bootstrap median is the reported estimate — §3.4 Chinese 42,871→42,847;
  Rec 4 Others 32,360→31,181; §3.5 per-group range ~32,000–45,000→~31,000–44,000.
  All in-text registry figures now agree with `registry_size_ci.csv` (Table 1/2).
- **C2 overclaim**: abstract and Recommendations 1–2 now flag Malay/Indian/Others
  N* as model projections pending validation; only Chinese is empirically validated
  (Spearman r=0.70, §3.6). "Mathematically necessary" softened in Rec 1 and §3.3.
- **C3 CI honesty**: bootstrap-CI lower-bound caveat (sampling variability only;
  excludes EM phasing and HWE model error) promoted from §2.4 into the abstract
  and Table 1/2 captions.
- **C5 EM citation/method**: replaced miscited Beatty [4] for EM phasing with
  Excoffier & Slatkin 1995 [20]; method now correctly described as full multi-locus
  phase-enumeration EM (matches `hwe_test._em_full_phase`), haplotypes retained ≥0.1%.
- `python-docx>=1.1` added to `analysis/requirements.txt` (report pipeline dep).

### Deferred (scope decision pending)
- C4 rare-haplotype cutoff sensitivity; broader patient validation for minority
  groups; Others cluster-stability bootstrap; external face-validity paragraph.

### Generated
- `HLA_Registry_Size_CMIO_v2.14.docx`

---

## v2.3.0 — 2026-04-29

### Added — quantified bias analyses (reviewer final polish)
- **`analysis/15_em_convergence.py`**: EM convergence test — reruns EM for Chinese
  at 500–45,018 donors; N* at 5k cap = 45,148 vs 41,727 at full sample (8.2%
  conservative overestimate). Figure S1 added to Supporting Analysis section.
- **`analysis/16_smoothing_sensitivity.py`**: Laplace pseudocount smoothing
  sensitivity — α=0.001 per haplotype; N* at 95% changes <3% for all groups
  (Chinese +0.9%, Malay +2.3%, Indian −3.1%, Others −1.9%); larger at 75%.
- **Tables 1 & 2**: Added "Signed-up target‡" row — N* ÷ 0.60 (40% attrition);
  shows range across CMIO groups per threshold.
- **§4.1 Limitations**: EM cap quantified (8.2% conservative bias); smoothing
  results cited; attrition adjustment formula stated.
- **Supporting Analysis** section with Figure S1 (EM convergence).

### Generated
- `HLA_Registry_Size_CMIO_v2.13.docx`

---

## v2.2.0 — 2026-04-29

### Changed
- **Figure 1 flowchart**: fixed arrow alignment — arrows now run between boxes
  (bottom of one box to top of the next) rather than inside boxes.
- **Figure numbering**: renumbered sequentially in document order —
  Fig 1 pipeline, Fig 2 CI bar chart (was unlabelled), Fig 3 10-locus partial match,
  Fig 4 8-locus partial match, Fig 5 sensitivity, Fig 6 validation scatter,
  Fig 7 Others PCA scatter (was Fig 1).

### Generated
- `HLA_Registry_Size_CMIO_v2.12.docx`

---

## v2.1.0 — 2026-04-29

### Changed — flow and comprehension improvements
- **Figure 0** (new): Methods pipeline flowchart (`analysis/14_pipeline_flowchart.py`)
  showing data → EM → HWE → C(N) → N* → bootstrap CI.
- **§2.3 Eq(4)**: Added verbal explanation paragraph — what C(N) means intuitively
  and why rare diplotypes make convergence slow.
- **§3 Results intro**: Added 7-line reading guide orienting the reader to which
  sections are primary, secondary, robustness, and exploratory.
- **§3 section order**: Reordered §3.4–3.7 — Partial Match (§3.4), Sensitivity (§3.5),
  Validation (§3.6), Others Exploratory (§3.7). Others clearly labelled as exploratory.
- **Table 1 footnote**: Weighted Average row explained as a mathematical convenience,
  not a policy target; per-group targets are the operative planning figures.
- **Glossary** (new table): 15 abbreviations defined (AFND, BMDP, CI, CMIO, EM,
  HSCT, HLA, HSA, HWE, LD, MLE, N*, PCA, RMSE, SCBB).
- All §3.x cross-references updated to match new numbering.

### Generated
- `HLA_Registry_Size_CMIO_v2.11.docx`

---

## v2.0.0 — 2026-04-29

### Changed — reviewer response (text-only, no recomputation)
- **§2.4 Bootstrap CI**: removed trivial "by construction" sentence; replaced with
  explicit scope statement — CIs capture donor-count sampling variability only, not
  EM phasing or HWE model uncertainty.
- **§3.7 Validation**: explicitly flagged Indian (1 shared haplotype = no validation)
  and stated that Malay/Indian/Others estimates are model-derived projections;
  Chinese is the primary validated result.
- **§4.1 Limitations**: expanded from 5 to 7 substantive points —
  (1) HWE bias direction uncertain; Indian/Others flagged as exploratory;
  (2) 5,000 cap binds materially only for Chinese — common haplotypes robust at 5k;
  (3) NEW: donor attrition — N* is biologically matched minimum; signed-up targets
      must exceed N* by 30–50% to account for real-world attrition;
  (4) NEW: N* is a lower bound — unobserved haplotypes assigned zero frequency,
      most material at 95% coverage threshold;
  (5) Others cluster: added note that silhouette reflects strong HLA–ancestry signal
      but cluster stability was not bootstrap-validated.

### Generated
- `HLA_Registry_Size_CMIO_v2.10.docx`

---

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
