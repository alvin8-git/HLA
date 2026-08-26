---
title: "A frequency floor above 1/(2n) inflates unrelated donor registry coverage estimates by up to three orders of magnitude"
running_title: "Rare-haplotype floors and registry coverage"
author: "Alvin Ng Yu-Jin"
affiliation: "National University Hospital / Singapore Cord Blood Bank, Singapore"
correspondence: "alvin1976sg@gmail.com"
target_journal: "Bone Marrow Transplantation (Nature Portfolio)"
article_type: "Article"
reference_style: "Vancouver (numbered), per BMT author guidelines"
version: "3.0 — full-sample re-derivation (EM cap removed); supersedes v2.0, which used a 5,000-individual EM cap"
date: "2026-08-20"
---

# A frequency floor above 1/(2n) inflates unrelated donor registry coverage estimates by up to three orders of magnitude

**Alvin Ng Yu-Jin**

National University Hospital / Singapore Cord Blood Bank, Singapore
Correspondence: alvin1976sg@gmail.com

---

## Abstract

Registry-size models estimate how many unrelated donors a country must recruit for a given HLA match probability, and are used to set national targets. Estimating the haplotype frequencies they consume requires discarding rare haplotypes below a frequency floor. We asked where that floor becomes harmful, and whether the discarded tail is real signal or expectation–maximisation (EM) phase-enumeration noise.

Five-locus (HLA-A, -B, -C, -DRB1, -DQB1) haplotype frequencies were estimated by EM with full phase enumeration from 59,186 Singapore Bone Marrow Donor Programme and Singapore Cord Blood Bank donors and cord blood units, using every five-locus-typed individual per group with no input cap. Population coverage C(N) = Σ F_k[1 − (1 − F_k)^N] was recomputed across floors from 0 to 10⁻³.

The floor is harmless below, and destructive above, the frequency of a singleton haplotype, 1/(2n). For Chinese donors (n = 45,754; singleton 1.09 × 10⁻⁵), discarding the 225,000 haplotypes below 10⁻⁶ — 96% of all distinct haplotypes — costs 0.03% of frequency mass and leaves the minimum registry for 95% coverage unchanged at 87 million. Raising the floor to 3 × 10⁻⁵ drops it to 16.4 million, to 10⁻⁴ to 3.2 million, and to 10⁻³ to 41,727: a 2,098-fold collapse, achieved entirely by deleting mass the model then renormalises away. The independently sampled Others group (n = 3,941; singleton 1.27 × 10⁻⁴) reproduces the rule with its threshold displaced accordingly, remaining flat to 10⁻⁴ and collapsing 791-fold by 10⁻³. Because the sub-singleton tail is provably inert, it cannot be dismissed as phase-ambiguity artefact; the damage is done to haplotypes the sample genuinely resolves.

A corollary follows that constrains all such analyses: since 1/(2n) differs between groups, a single floor applied across groups of unequal size biases them unequally, and minimum registry sizes are not comparable across groups whose rare tails are sampled to different depths. Floors should be set below 1/(2n) per group, retained frequency mass reported, and cross-group comparisons made only at matched sampling depth.

**Keywords:** unrelated donor registry; HLA haplotype frequency; match probability; expectation–maximisation; rare variants; haematopoietic cell transplantation; Singapore; health equity

---

## 1. Introduction

Allogeneic haematopoietic cell transplantation depends on locating an HLA-compatible donor. Most patients lack a matched sibling and rely on unrelated volunteer registries, where high-resolution matching is associated with survival [1]. Because HLA haplotype distributions differ markedly between ancestry groups, patients of non-European ancestry are systematically disadvantaged in registries dominated by donors of European descent [2,3]. For the emerging registries of Asia, where donor pools are small relative to population diversity, the resulting planning question is acute [14,15]: how many donors must be recruited, from which communities?

The standard framework, introduced by Beatty et al. [4] and developed for the US registry by Maiers et al. [5] and Gragert et al. [2], computes population coverage as a haplotype-frequency-weighted sum of per-patient match probabilities. Applied to national data it produces specific recruitment targets, which is precisely what makes it attractive to policymakers.

Estimating the haplotype frequencies that feed this calculation requires a decision about rare haplotypes. EM estimation over unphased multi-locus genotypes [6] yields a long-tailed distribution, and analysts routinely impose a frequency floor — commonly 10⁻³ or 10⁻⁴ — then renormalise the survivors to sum to unity.

**What is already known, and what is not.** That rare-haplotype handling matters is not new. Maiers et al. observed that their high-resolution US sample was "slightly tipped, favoring common haplotypes at the expense of rare haplotypes" [5]. The GENE[RATE] pipeline retains haplotypes to 10⁻⁴ rather than 10⁻³ [9], a design choice that itself reflects awareness of tail sensitivity. We therefore make no claim to having discovered that truncation biases coverage estimates.

What has not been established, and what this paper provides, is a rule for where the floor may safely be placed. We show that the boundary is the frequency of a singleton haplotype in the sample, 1/(2n): below it a floor removes only phase-enumeration noise and changes nothing, above it the estimate degrades rapidly and without any internal signal that it has done so. Four consequences follow: (i) the damage is not gradual but hinges on a sample-size-dependent threshold, so a floor that is safe in one study is unsafe in another with fewer donors; (ii) the discarded tail cannot be dismissed as EM artefact, because the artefactual portion is demonstrably inert; (iii) since 1/(2n) differs between groups, a single floor biases unequally-sampled groups unequally, and registry sizes are not comparable across groups sampled to different depths; and (iv) the effect is large enough — up to 2,098-fold here — to have carried an artefact undetected through a national planning analysis, including an earlier version of our own.

We also report which substantive conclusions about Singapore's CMIO populations [7,8] survive the artefact, since these are what a registry director must act on.

---

## 2. Materials and methods

### 2.1 Study population

HLA typing data were obtained for 59,186 volunteer donors and cord blood units accrued 2005–2020 by the Singapore Bone Marrow Donor Programme (BMDP) and Singapore Cord Blood Bank (SCBB); typing methodology and quality control are described in the source characterisation study [8]. Restricting to complete high-resolution typing at all five loci gave 44,400 Chinese, 5,578 Malay, 5,490 Indian and 3,767 Others individuals. An independent set of 564 patient–donor pairs was provided by the Health Sciences Authority (HSA), Singapore.

Cord blood units and adult volunteer donors were pooled for haplotype estimation, as in the source study. Because pooling two collection streams could in principle bias frequency estimates if their catchments differ, we quantified the composition: cord blood units constitute 2.3% of Chinese, 1.7% of Malay, 1.7% of Indian and 6.2% of Others individuals. At these proportions, cord blood cannot materially shift the estimated frequency distributions, though the higher Others fraction warrants note. Separate clinical implications of pooling are addressed in §4.5.

### 2.2 Haplotype frequency estimation

Five-locus haplotype frequencies were estimated per ethnicity by EM with full multi-locus phase enumeration in the Excoffier–Slatkin formulation [6], assuming Hardy–Weinberg equilibrium (HWE) to infer phase. **Every five-locus-typed individual in each group was used; no input cap was applied.** This matters: an earlier version of this analysis capped the EM at 5,000 individuals per group, which at a 10⁻⁴ floor inflated the Chinese minimum registry size by 264% (11,487,962 capped versus 3,153,571 uncapped), because a 5,000-individual EM cannot resolve phase in the rare tail and retains spurious low-frequency haplotypes. A cap and a floor therefore cannot be chosen independently, and results computed under a cap are not interpretable as tail behaviour.

Reference distributions for the same four populations were obtained from the HLA-net GENE[RATE] pipeline [9]. **GENE[RATE] applies its own 10⁻⁴ floor; we therefore label these distributions "10⁻⁴-floored", never "full", and treat N* values derived from them as lower bounds.**

### 2.3 Floor-sensitivity design

Our primary analysis re-estimates haplotype frequencies once with no floor and no cap, then applies floors of 0, 10⁻⁶, 10⁻⁵, 3×10⁻⁵, 10⁻⁴, 3×10⁻⁴ and 10⁻³ to that single distribution, renormalising the survivors at each step. Because every floor is applied to the same underlying estimate, differences isolate the floor's effect and are not confounded by re-estimation. The grid brackets the singleton frequency 1/(2n) rather than sitting above it, which the earlier capped design did not.

This was run for Chinese (n = 45,754, the largest and most deeply sampled group, singleton frequency 1.09 × 10⁻⁵) and Others (n = 3,941, the smallest, singleton frequency 1.27 × 10⁻⁴). The order-of-magnitude gap between their singleton frequencies is what allows the threshold rule to be tested rather than merely observed: if the boundary is set by sampling depth, the two curves should break at different floors, and in the same ratio as their sample sizes.

### 2.4 Coverage model

Diplotype frequencies were derived under HWE: F(h_i, h_i) = f_i² and F(h_i, h_j) = 2f_i f_j. For a patient of diplotype d_k at frequency F_k, the probability that at least one of N independently sampled donors carries it is 1 − (1 − F_k)^N, so

C(N) = Σ_{k=1}^{m} F_k · [1 − (1 − F_k)^N]

N* was located by logarithmic binary search, **capped at 10⁷ donors; values reaching the cap are reported as "> 10⁷" and any ratio computed against them is a lower bound, not a measurement.** Same-ancestry and cross-ethnic (Singapore-population-weighted pool [7]) variants were computed, as was partial-match coverage at ≥9/10 and ≥8/10.

### 2.5 External benchmark

We compared model predictions against published match likelihoods from the US registry, which at 10.5 million donors delivers approximately 75% 8/8 match likelihood for White patients of European descent [2]. **This benchmark is asymmetric and we use it only as such:** it has strong power to *reject* a model predicting complete matching at small registry size, but can only fail to reject a plausible model, never confirm it. It is also not like-for-like — the US figure is 8/8 from a multi-ancestry registry whose European-descent subset is smaller than the 10.5 million total, whereas our modelled figure is 10/10 from a hypothetical same-ancestry pool. Both differences make the comparison approximate; we draw only order-of-magnitude inferences from it.

### 2.6 Uncertainty, stratification, and software

Uncertainty quantification was attempted by Dirichlet parametric bootstrap (B = 1,000) [16] and subsequently withdrawn: §4.2 reports the failure, which is itself informative about how this class of model misbehaves on long-tailed haplotype distributions. Ancestry stratification of the Others group used principal component analysis on binary indicators of alleles at ≥1% frequency (n = 3,847 with complete typing), then k-means with k selected by silhouette coefficient; cluster ancestry was inferred from top-haplotype signatures against AFND [10] and published references [8,11]. Analyses used Python (NumPy, pandas, scikit-learn); code and derived data are at https://github.com/alvin8-git/HLA.

---

## 3. Results

### 3.1 The floor is harmless below 1/(2n) and destructive above it

Applying successive floors to a single uncapped, unfloored EM estimate localises the harm precisely, and it does not fall at a fixed frequency (Table 1). It falls at the frequency of a singleton haplotype in that group's sample.

**Table 1.** Effect of frequency floor, applied to one uncapped unfloored EM estimate per population. The singleton frequency 1/(2n) is 1.09 × 10⁻⁵ for Chinese and 1.27 × 10⁻⁴ for Others.

| Population | Floor | vs 1/(2n) | Haplotypes | Mass retained | Coverage at N = 50,000 | N* at 95% | Inflation vs unfloored |
|---|---|---|---|---|---|---|---|
| Chinese | none | — | 234,568 | 100% | 31.7% | 87,530,956 | 1.0× |
| Chinese | 10⁻⁶ | below | 9,595 | 100.0% | 31.7% | 86,971,552 | 1.0× |
| Chinese | 10⁻⁵ | below | 8,537 | 99.3% | 32.2% | 76,579,448 | 1.1× |
| Chinese | 3×10⁻⁵ | **above** | 3,198 | 91.3% | 39.7% | 16,405,166 | 5.3× |
| Chinese | 10⁻⁴ | **above** | 1,253 | 80.8% | 52.3% | 3,153,571 | 27.8× |
| Chinese | 3×10⁻⁴ | **above** | 488 | 67.9% | 70.9% | 528,501 | 165.6× |
| Chinese | 10⁻³ | **above** | 136 | 49.0% | 96.1% | 41,727 | **2,097.7×** |
| Others | none | — | 65,020 | 100% | 17.6% | 26,987,290 | 1.0× |
| Others | 10⁻⁶ | below | 4,654 | 100.0% | 17.6% | 26,987,264 | 1.0× |
| Others | 10⁻⁵ | below | 3,462 | 99.0% | 18.0% | 23,750,297 | 1.1× |
| Others | 3×10⁻⁵ | below | 3,182 | 98.6% | 18.3% | 23,056,443 | 1.2× |
| Others | 10⁻⁴ | below | 3,054 | 97.9% | 18.6% | 22,219,299 | 1.2× |
| Others | 3×10⁻⁴ | **above** | 571 | 58.7% | 59.4% | 732,600 | 36.8× |
| Others | 10⁻³ | **above** | 126 | 35.7% | 97.4% | 34,110 | **791.2×** |

The two curves are the same curve, displaced by sampling depth. In the Chinese group the estimate is stable to 10⁻⁵ and has already lost a factor of 5.3 by 3 × 10⁻⁵; in the Others group, whose singleton frequency is an order of magnitude higher, it remains stable all the way to 10⁻⁴ and breaks only at 3 × 10⁻⁴. The break tracks 1/(2n), not any absolute frequency.

This is why a 10⁻⁴ floor cannot be described as safe or unsafe in the abstract. Applied to the Others sample it sits below the singleton frequency and costs a factor of 1.2. Applied to the Chinese sample it sits 9.2-fold above it and costs a factor of 27.8.

### 3.2 The discarded tail is not EM phase-ambiguity noise

An EM run over unphased five-locus genotypes distributes fractional probability across many phase resolutions per ambiguous individual, and the uncapped Chinese run returned 234,568 distinct haplotypes from 45,754 individuals — more distinct haplotypes than there are chromosomes carrying any single one. It is therefore reasonable to ask whether the tail is largely estimation noise, in which case truncation would be legitimate denoising.

The floor curve answers this cleanly, and the answer separates the tail into two parts. Removing every haplotype below 10⁻⁶ eliminates 224,973 of the 234,568 Chinese haplotypes — 95.9% of all distinct haplotypes — and costs **0.03% of frequency mass**, leaving coverage at 50,000 donors unchanged at 31.7% and N* within 0.6%. That portion of the tail is indeed noise, and discarding it is free.

The damage is done above the singleton frequency, to haplotypes the sample genuinely resolves. Between 10⁻⁵ and 10⁻³ the Chinese estimate falls 1,835-fold while shedding half its frequency mass. These are not phase artefacts: they are haplotypes observed in multiple individuals, and they are exactly the ones that govern the approach to high coverage, since F_k[1 − (1 − F_k)^N] converges slowly precisely when F_k is small.

The practical rule is therefore specific and checkable before any coverage number is computed: **set the floor below 1/(2n) for the group being analysed, and report the retained frequency mass.** Retained mass, not haplotype count, is the diagnostic — 136 haplotypes sounds like a reasonable working set until one notices it represents half a population.

### 3.3 Replication in an independent distribution

Applying the same 10⁻³ truncation to the GENE[RATE] 10⁻⁴-floored distributions reproduces the effect in all four populations (Table 2). At 50,000 donors, estimated 10/10 coverage rises from 20.8–41.8% to 93.5–97.2% on truncation; the truncated model also predicts complete (100.0%) coverage for every population at 500,000 donors.

**Table 2.** Estimated 10/10 coverage (%) at fixed registry sizes, GENE[RATE] distributions, 10⁻⁴-floored versus additionally truncated at 10⁻³.

| Registry size N | Chinese 10⁻⁴ | Chinese 10⁻³ | Malay 10⁻⁴ | Malay 10⁻³ | Indian 10⁻⁴ | Indian 10⁻³ | Others 10⁻⁴ | Others 10⁻³ |
|---|---|---|---|---|---|---|---|---|
| 25,000 | 33.5 | 86.8 | 28.8 | 92.1 | 14.4 | 90.0 | 14.6 | 92.0 |
| 50,000 | 41.8 | 93.5 | 36.6 | 96.5 | 20.2 | 96.3 | 20.8 | 97.2 |
| 100,000 | 50.6 | 97.5 | 45.2 | 98.9 | 27.3 | 99.2 | 28.6 | 99.4 |
| 500,000 | 70.4 | 100.0 | 65.3 | 100.0 | 48.0 | 100.0 | 50.8 | 100.0 |
| 1,000,000 | 78.0 | 100.0 | 72.9 | 100.0 | 57.8 | 100.0 | 60.9 | 100.0 |
| 10,000,000 | 94.9 | 100.0 | 92.2 | 100.0 | 86.9 | 100.0 | 89.4 | 100.0 |

Agreement between two independent estimators rules out an implementation bug. It does not establish that either distribution's tail is correctly estimated, since both share the same rare-tail estimation problem (§4.5).

### 3.4 External benchmarking rejects the truncated model

The truncated model predicts 100.0% 10/10 coverage for every CMIO population from 500,000 donors. No registry achieves complete matching for any population; the US registry at 10.5 million donors delivers approximately 75% 8/8 for its best-served group [2]. The truncated prediction is not merely optimistic but structurally impossible, and the benchmark rejects it decisively.

The untruncated model predicts 94.9% for Chinese and 86.9–92.2% for other CMIO groups at 10 million same-ancestry donors — the right order of magnitude, and **not contradicted by** observed registry performance. We deliberately claim no more than that: as set out in §2.5, this benchmark cannot confirm a model, and the 8/8-versus-10/10 and whole-registry-versus-same-ancestry mismatches make it approximate. Its value is that it is cheap, uses only published data, and would have exposed the artefact immediately.

### 3.5 Conclusions robust across floors

Truncation acts in the same direction in every population, so comparative findings are far more robust than absolute sizes. The following held at every floor tested.

**Cross-ethnic matching cannot substitute for same-ancestry recruitment.** With a combined pool weighted to Singapore's composition (74.3% Chinese), Malay, Indian and Others patients failed to reach 75% coverage below the 10⁷ ceiling under every regime; Chinese patients were servable but required roughly twice the same-ancestry figure. Donors are not interchangeable across ancestry groups, consistent with other multiethnic registries [3].

**DQB1 matching remains comparatively cheap.** Recomputed on the uncapped distributions, DRB1–DQB1 composite D′ ranges from 0.90 (Others) to 0.96 (Chinese) and HLA-B–C from 0.85 (Others) to 0.93 (Indian). These are lower than the values reported at the 10⁻³ floor (0.93–0.99 and 0.95–0.99). We note the direction — a haplotype set restricted to common haplotypes yields higher apparent D′ — but did not isolate whether this reflects the floor, the removal of the sample cap, or both, so we do not attribute a mechanism. The qualitative conclusion holds: DRB1 and DQB1 remain strongly linked, so the incremental cost of requiring DQB1 identity is modest and routine 10-allele typing is justified.

**Single-allele relaxation is the largest available lever, and larger than previously estimated.** Recomputed on the uncapped distributions, relaxing to ≥9/10 reduces N* 17.9-fold for Chinese patients (2,998,168 to 167,194) and 5.2-fold for Others (20,751,252 to 3,991,544). At a fixed 50,000-donor registry the same relaxation lifts modelled coverage from 55.2% to 86.9% (Chinese) and 19.4% to 43.5% (Others) (Table 3).

**Table 3.** Effect of match stringency on the uncapped 10⁻⁴-floored distributions, by Monte-Carlo sampling of 3,000 patient diplotypes against the complete donor diplotype set (the exact algorithm is quadratic in diplotype count and is not computable at this scale). Absolute values inherit the floor caveat of §3.1 and are lower bounds; the ratios and the coverage gains are the quantities to use.

| Population | 10/10 at 95% | ≥9/10 at 95% | Reduction | Coverage at 50,000: 10/10 → ≥9/10 → ≥8/10ᵃ |
|---|---|---|---|---|
| Chinese | 2,998,168 | 167,194 | 17.9× | 55.2% → 86.9% → 98.6% |
| Others | 20,751,252 | 3,991,544 | 5.2× | 19.4% → 43.5% → 78.3% |

ᵃ **≥8/10 falls below currently accepted match stringency for most indications and must not be used as a recruitment target.** It is shown only because the continued steepening is mechanistically informative about linkage structure.

One caution applies to Table 3: the ≥9/10 model treats a mismatch at any of the five loci as equivalent, which the outcomes literature does not support — permissiveness is locus-specific, which is precisely why DPB1 mismatching is assessed through a T-cell-epitope framework rather than as a flat allele count [12,13]. The reduction factors should be read as upper bounds on the benefit of an undifferentiated single-mismatch policy. The Monte-Carlo 10/10 values agree with the deterministic pipeline to within 6% (2,998,168 versus 3,153,571 for Chinese; 20,751,252 versus 22,219,299 for Others), which is within the expected sampling error for 3,000 sampled patients.

**The "Others" category is genetically stratified.** Clustering selected k = 3, though silhouette coefficients were nearly flat across k = 2–5 (0.225, 0.242, 0.238, 0.240) and below the 0.25 threshold conventionally taken to indicate weak structure; the clustering metric alone would not support a strong claim. The decisive evidence is instead the top haplotype of each cluster, each an unambiguous population signature: A\*01:01~B\*08:01~C\*07:01~DRB1\*03:01~DQB1\*02:01 (12.2%; the 8.1 ancestral haplotype of Northern Europeans [11]); A\*24:07~B\*35:05~C\*04:01~DRB1\*12:02~DQB1\*03:01 (9.1%; Filipino/Southeast Asian [10]); and A\*02:07~B\*46:01~C\*01:02~DRB1\*09:01~DQB1\*03:03 (4.7%; Chinese-specific [8]).

Independent statistical corroboration comes from Hardy–Weinberg testing. The Others population shows significant departure at **all five loci** (p from 3.7×10⁻⁵ at DQB1 to 3.8×10⁻²⁸ at HLA-B), with observed heterozygosity below expected at every locus — the classical Wahlund signature of pooling stratified subpopulations. By contrast Chinese and Malay show no significant departure at any locus and Indian at three. Cluster-level requirements differ substantially (2,025,221–4,722,492 donors at 95%), and pooled-Others figures should be treated as inadmissible for planning: the HWE violation invalidates the random-mating assumption underlying F(h_i,h_j) = 2f_i f_j for that group specifically, independently of any truncation issue.

### 3.6 The available patient dataset cannot validate the coverage model

EM frequencies were compared against frequencies observed in the HSA patient–donor pairs. Only the Chinese comparison supported inference (33 shared haplotypes; Spearman ρ = 0.70, p < 0.001; RMSE 0.0094); Malay (11), Others (4) and Indian (1) had too few.

We additionally report an output omitted from our earlier analysis: model-predicted per-patient match probability was 0.008–0.028 against an observed match rate of 1.00 in all four groups. This is a selection artefact — the HSA pairs are transplanted, hence matched, by construction — and is not evidence of model failure. We report it to state plainly what follows: **this dataset cannot serve as an outcome-based validation of the coverage model**, and no claim in this paper rests on it.

---

## 4. Discussion

### 4.1 Principal finding

The harm done by a rare-haplotype frequency floor is a threshold phenomenon, and the threshold is the singleton frequency 1/(2n) rather than any fixed value. Below it, floors are nearly free: they remove the great majority of distinct haplotypes — 95.9% of the Chinese set — while costing 0.03% of frequency mass and leaving coverage estimates unchanged. Above it the estimate degrades rapidly, and a 10⁻³ floor understated the registry required for 95% Chinese coverage by 2,098-fold.

The mechanism is visible in C(N). The term F_k[1 − (1 − F_k)^N] approaches its limit slowly exactly when F_k is small, so the moderately-rare band governs high-coverage behaviour; truncation deletes that band while renormalisation inflates what remains. The reason the ultra-rare tail contributes little is the same equation read the other way: haplotypes at 10⁻⁶ contribute negligible weight regardless of N in the feasible range.

Practically, this means the guidance for the field is not "never truncate" — which would be computationally onerous and, on this evidence, unnecessary — but "truncate at 10⁻⁴, not 10⁻³, and report the mass you retained." Retained frequency mass, not haplotype count, is the diagnostic: 143 haplotypes sounds like a reasonable working set until one notices it represents half a population.

### 4.2 Why conventional uncertainty quantification did not detect it — and then failed outright

Bootstrap intervals computed at the 10⁻³ floor were narrow — approximately ±0.5% around the Chinese point estimate. They were correct as computed and misleading as a guide to total uncertainty, because they condition on the truncated distribution: resampling within a truncated support cannot recover deleted mass. A second variance source sat outside them entirely: at a 10⁻⁴ floor the Chinese N* was 11,487,962 from a 5,000-individual EM and 3,153,571 from all 45,754 — the target itself moving as the rare tail became resolvable.

At the corrected floor the bootstrap did not merely understate uncertainty; it failed in a way that forced its withdrawal. With K = 9,574 haplotypes and n = 44,400, the mean Dirichlet concentration per haplotype is ~4.6, and under a Dirichlet draw E[f²] = f² + f(1−f)/(n+1), so resampling systematically inflates the squared and product terms that form diplotype frequencies — by roughly five-fold for haplotypes near 10⁻⁶. Since precisely these terms govern high-coverage behaviour, every replicate overstates coverage: all 1,000 replicates fell below the point estimate, giving a "95% interval" of 50.8–53.9 million around an estimate of 87.4 million. An interval that excludes its own point estimate is not a confidence interval, and we report none. The symptom was already present at the old floor — in the superseded internal analysis, 18 of 32 rows had the estimate outside its interval, which had been accommodated by switching the reported point to the bootstrap median rather than by diagnosing the cause.

The general lessons: parametric bootstrap intervals are conditional on preprocessing, which here dominates; and on long-tailed frequency distributions the Dirichlet resampling scheme is biased for tail-driven functionals, not merely noisy. Restoring intervals would require a resampling scheme that preserves rare-tail structure, or analytic error propagation through C(N); we attempt neither, and state instead that uncertainty is real and unquantified. Headline registry sizes should be read to one significant figure.

### 4.3 Registry sizes are not comparable across unequally sampled groups

A consequence of the 1/(2n) rule constrains every cross-group statement in this literature, including our own earlier ones. Because the threshold moves with sample size, two groups analysed at the same floor are analysed at different effective resolutions. In this dataset the Chinese sample is 11.6 times larger than the Others sample, so at any common floor the Chinese rare tail is resolved and the Others rare tail is not.

The effect is large and acts in a counter-intuitive direction. At a 10⁻⁴ floor the Chinese estimate has already lost a factor of 27.8 relative to its unfloored value while the Others estimate has lost only 1.2, so the better-sampled group appears *easier* to serve than it is. Comparing the two at that floor would understate the Chinese requirement roughly twenty-fold relative to the Others requirement, purely as an artefact of who was sampled more deeply.

The same caution applies to the unfloored estimates, for the opposite reason: 26,987,290 for Others rests on 3,941 individuals, so most of that group's rare tail has never been observed and the figure is a loose lower bound. Chinese, at 45,754 individuals, is bounded more tightly.

We therefore make no quantitative claim ranking the four CMIO groups by registry requirement, and we regard such rankings in the prior literature — including versions of this analysis up to v2.15 — as unsafe unless sampling depth was matched or explicitly modelled. What survives is the within-group comparison, which shares a sample and therefore a resolution: relaxing match stringency, or drawing donors from a mismatched ancestry pool, changes the requirement for a given group by factors that are measured against that group's own baseline.

### 4.4 What a registry director should take from this

**Comparative conclusions should carry the planning weight.** That minority-ancestry patients cannot be served from a majority-ancestry pool; that Others is stratified and its pooled figure inadmissible; that single-allele relaxation buys 2.5–4.9-fold depending on population — these are ratio statements, robust to a bias acting in one direction, and directly actionable.

**Coverage at feasible size is the more honest headline, with a caveat.** Asking "how many donors for 95%?" invites an answer that may be unattainable in principle; asking "what does 50,000 donors deliver, and for whom?" yields a trackable number. Under the untruncated model a 50,000-donor same-ancestry registry delivers roughly 37% 10/10 coverage for Chinese patients and 18% for Others — a stark equity statement. We acknowledge the counter-argument: a registry director building a multi-year business case needs a numeric registrant goal to put before a ministry, and cost-effectiveness comparisons across strategies still require marginal reasoning about N. The honest resolution is to report both, with the coverage figure leading and N* explicitly flagged as preprocessing-dependent.

**When N* exceeds the national population, the answer is not domestic recruitment.** Our benchmark figure of 10 million same-ancestry donors for 94.9% Chinese coverage exceeds Singapore's entire ethnic-Chinese resident population several-fold. Reported without comment, such numbers invite misreading as aspirational domestic targets. They are better read as demonstrating that high same-ancestry coverage is reachable, if at all, only through international and diaspora registry linkage — partnership with regional registries serving the same ancestral populations [14,15] — and that domestic strategy must lean on the levers that do scale: formalised evidence-based single-mismatch protocols, expanded cord blood use with its more permissive matching standards, and haploidentical transplantation, which sidesteps the registry problem for patients with a family donor.

**Collect sub-ancestry data.** The three Others clusters cannot be targeted by outreach while the registration form records only "Others". Two questions at registration — parental birthplace, self-reported heritage — would convert a research finding into a recruitment instrument.

### 4.5 Limitations

*The unfloored distribution is not ground truth.* The unfloored EM tail remains an estimate. §3.2 shows the sub-singleton portion does not drive the result, and the band that does — between 1/(2n) and 10⁻³ — is now estimated from every available individual rather than a 5,000-person subsample. But no sample resolves haplotypes rarer than its own singleton frequency, so every estimate here is a lower bound on the true requirement, and the bound is looser for the smaller groups (§4.3). The GENE[RATE] distributions are floored at 10⁻⁴ and are bounded likewise. We claim the 10⁻³ estimates are wrong by two to three orders of magnitude and that estimates from floors below 1/(2n) are the right order; we do not claim they are precise.

*HWE.* Diplotype frequencies assume random mating. This fails significantly at all five loci in Others and three in Indian, so figures for those groups — the two already worst served — carry additional uncertainty. For Others the failure is structural (stratification), and cluster-level rather than pooled estimates should be used.

*Stratification compounds truncation for Others specifically.* A floor applied to a pooled, stratified population preferentially deletes haplotypes common within a sub-cluster but diluted in the pool, so the two biases act in the same direction for the group least able to absorb them. We have not disentangled them.

*DPB1 was not typed.* Contemporary standard of care in many centres is a 12/12 panel including DPB1, with permissive T-cell-epitope mismatching guiding selection [12,13]. Every coverage figure here is therefore an **upper bound** on true match probability under 12/12 practice: the DPB1 gap compounds the truncation problem rather than substituting for it. Our five-locus framework cannot address it.

*Cord blood.* Although cord blood units are only 1.7–6.2% of each population and cannot materially bias frequency estimates (§2.1), they carry a different clinical matching standard (approximately 4/6–6/8, lower resolution) and cannot be "recruited". A registry size expressed over a pooled population therefore conflates two products; adult-donor-only figures would be preferable and are a target for future work.

*Attrition and sampling.* The model assumes random donor sampling and biologically matched donors; real-world attrition is typically 30–50%, so recruitment targets must exceed any modelled N* accordingly.

*Scope of the field-wide claim.* We demonstrate that this artefact propagated undetected through one national planning analysis — our own. We have not audited published registry-size studies for the same error, and make no claim about its prevalence.

### 4.6 Relationship to our earlier internal analysis

This work supersedes an unpublished internal draft that reported per-population targets of approximately 40,000–45,000 donors for 95% 10/10 coverage. Those figures were computed under a 10⁻³ floor and are, on the present analysis, artefactual; we withdraw them. During this reanalysis we also identified several arithmetic and labelling errors in that draft, including an understatement of the ≥9/10 benefit by roughly an order of magnitude and a mislabelled summary row. A full itemised audit accompanies the deposited analysis code rather than this manuscript, where it would not belong; we note it here only to record that the headline figures of the earlier draft should not be cited.

---

## 5. Conclusions

Rare-haplotype frequency floors damage registry coverage estimates through a threshold set by sampling depth. A floor below the singleton frequency 1/(2n) is nearly free — discarding 95.9% of distinct Chinese haplotypes cost 0.03% of frequency mass and left the estimate unchanged — while a floor above it degrades the estimate rapidly, reaching 2,098-fold for Chinese and 791-fold for Others at 10⁻³. The discarded mass is not EM phase-enumeration artefact, because the artefactual sub-singleton portion is demonstrably inert. Only estimates from floors below 1/(2n) survive benchmarking against observed international registry performance.

We recommend that registry-size analyses use a floor no higher than 10⁻⁴; report retained frequency mass rather than haplotype count; benchmark predictions against observed registry match rates before publication; express findings as coverage attainable at feasible registry size; and treat bootstrap intervals as conditional on preprocessing. All estimates of this kind, including ours, remain upper bounds on match probability while DPB1 is untyped.

For Singapore, the equity conclusions are unchanged and sharpened: minority-ancestry patients cannot be served from a majority-ancestry pool at any registry size; the administrative "Others" category is genetically stratified and its pooled estimates inadmissible; and because the registry sizes implied by high coverage targets exceed the national population, the actionable levers are international registry linkage, evidence-based single-mismatch protocols, cord blood, and haploidentical transplantation rather than domestic recruitment alone.

---

## Data availability

Analysis code and derived data are available at https://github.com/alvin8-git/HLA. Individual-level HLA typing data are held by the Singapore Bone Marrow Donor Programme, the Singapore Cord Blood Bank and the Health Sciences Authority under their governance arrangements.

## Author contributions

A.N.Y.-J. designed the study, implemented the analyses, interpreted the results and wrote the manuscript.

## Competing interests

The author declares no competing interests.

## Funding

No external funding was received.

## Acknowledgements

The author thanks the Singapore Bone Marrow Donor Programme, the Singapore Cord Blood Bank and the Health Sciences Authority for data access, and Prof. Aloysius Ho Yew Leng, whose critical review of an earlier draft prompted the reconsideration of matching standards that led to this reanalysis.

---

## References

1. Lee SJ, Klein J, Haagenson M, Baxter-Lowe LA, Confer DL, Eapen M, et al. High-resolution donor-recipient HLA matching contributes to the success of unrelated donor marrow transplantation. Blood. 2007;110(13):4576–83.

2. Gragert L, Eapen M, Williams E, Freeman J, Spellman S, Baitty R, et al. HLA match likelihoods for hematopoietic stem-cell grafts in the U.S. registry. N Engl J Med. 2014;371(4):339–48.

3. Halagan M, Manor S, Shriki N, Yaniv I, Sela A, Maiers M, et al. East meets West: impact of ethnicity on donor match rates in the Ezer Mizion Bone Marrow Donor Registry. Biol Blood Marrow Transplant. 2017;23(8):1381–6.

4. Beatty PG, Mori M, Milford E. Impact of racial genetic polymorphism on the probability of finding an HLA-matched donor. Transplantation. 1995;60(8):778–83.

5. Maiers M, Gragert L, Klitz W. High-resolution HLA alleles and haplotypes in the United States population. Hum Immunol. 2007;68(9):779–88.

6. Excoffier L, Slatkin M. Maximum-likelihood estimation of molecular haplotype frequencies in a diploid population. Mol Biol Evol. 1995;12(5):921–7.

7. Singapore Department of Statistics. Census of Population 2020. Singapore: Department of Statistics; 2021.

8. Ng AYJ, Moshi GB, Prasath A, Teo D, Lim YA, Ang AL, et al. Human leukocyte antigen allele and haplotype frequencies in Singapore bone marrow donors and cord blood units. Blood Cell Ther. 2022;5(3):86–95.

9. Nunes JM, Buhler S, Roessli D, Sanchez-Mazas A; HLA-net 2013 collaboration. The HLA-net GENE[RATE] pipeline for effective HLA data analysis and its application to 145 population samples from Europe and neighbouring areas. Tissue Antigens. 2014;83(5):307–23.

10. Gonzalez-Galarza FF, McCabe A, Santos EJMD, Jones J, Takeshita L, Ortega-Rivera ND, et al. Allele frequency net database (AFND) 2020 update: gold-standard data classification, open access genotype data and new query tools. Nucleic Acids Res. 2020;48(D1):D783–8.

11. Price P, Witt C, Allcock R, Sayer D, Garlepp M, Kok CC, et al. The genetic basis for the association of the 8.1 ancestral haplotype (A1, B8, DR3) with multiple immunopathological diseases. Immunol Rev. 1999;167:257–74.

12. Fleischhauer K, Shaw BE, Gooley T, Malkki M, Bardy P, Bignon JD, et al. Effect of T-cell-epitope matching at HLA-DPB1 in recipients of unrelated-donor haematopoietic-cell transplantation: a retrospective study. Lancet Oncol. 2012;13(4):366–74.

13. Pidala J, Lee SJ, Ahn KW, Spellman S, Wang HL, Aljurf M, et al. Nonpermissive HLA-DPB1 mismatch increases mortality after myeloablative unrelated allogeneic hematopoietic cell transplantation. Blood. 2014;124(16):2596–606.

14. Aljurf M, Weisdorf D, Alfraih F, Szer J, Müller C, Confer D, et al. Worldwide Network for Blood & Marrow Transplantation (WBMT) special article: challenges facing emerging alternate donor registries. Bone Marrow Transplant. 2019;54(8):1179–88.

15. Lim YA, Teo D, Ang AL, Chan LL, Kuperan P. HLA allele and haplotype frequencies in unrelated bone marrow donor registries in Asia. Transpl Immunol. 2010;22(3–4):166–74.

16. Efron B, Tibshirani RJ. An introduction to the bootstrap. New York: Chapman & Hall; 1993.
