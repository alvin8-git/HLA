# Version History


> **This file is a change log, not a statement of current results.** Entries are
> dated and describe what changed on that date, so figures inside older entries are
> deliberately left as they were then — that trail is how a reader traces a
> correction. For current numbers see [README](README.md) and
> [Documentation.md](Documentation.md); the manuscript is **v2.15c** at a 1×10⁻³
> frequency floor. For why the floor sets the scale of every figure in this file,
> see [docs/explanation-frequency-floor.md](docs/explanation-frequency-floor.md).

## v2.15c — numeric consistency audit (2026-08-26)

Audited every comma-formatted figure in the built .docx against the snapshot
CSVs. Four defects, all from builder-side literals that survived a re-run of
the model:

  - Section 3.3 said Chinese cross-ethnic needs "93,348 donors ... roughly
    twice the same-ethnicity target" while Table 3, generated from the CSV,
    said 69,042 on the same page. 93,348 is the pre-census-weights value. Now
    derived; the multiple is 1.7x.
  - The deck was worse: its whole Table 3 was hardcoded and stale, showing
    Malay as ">10 million" at every threshold when the model gives 150,696 /
    360,396 / 634,813 / 1,391,093. Both deck tables now read from the CSVs.
  - "8/8 targets are typically 600-1,200 fewer donors" holds for Chinese (693),
    Indian (1,195) and Others (640) but not Malay, at 3,832 (9.6% of its 10/10
    target). Malay is named as the exception in both documents.
  - Section 3.7 described "3,847 fully five-locus-typed Others donors". The
    registry has 3,767; 11_others_stratification.py applies no MAIN_SOURCES
    filter, so its cohort includes 128 HSA donor and 46 HSA patient records.
    Relabelled, with a Limitations paragraph stating that Section 3.7 uses a
    wider cohort than every other result. NOT yet fixed at source — re-running
    script 11 with the filter would move Table 4, Table 5, Figure 7 and the
    63,856 planning target.
  - The abstract's "infeasible ... regardless of registry size" is false for
    Malay at 1,391,093; softened.

check_freshness.py does NOT audit the .docx or .pptx and would not have caught
any of these. The audit that did is not yet a committed script.

### Result 5 is floor-conditional — not yet reflected in the manuscript

Section 3.5 attributes the flat sensitivity bars to "a structural property of
CMIO haplotype diversity". It is a property of the 0.1% floor. At 1e-3 the
floor truncates every group to the same common core (140/144/137/123
haplotypes) so the four N* land within 1.39x and any weighted mean is flat. At
1e-6 the floor stops binding, retained haplotypes diverge (Chinese 9,574 vs
Malay 3,134), per-group N* spread to 5.28x, and the same four scenarios span
67.1% instead of 3.4%. The outlier also flips from Others (3.2% weight, cannot
move the mean) to Chinese (74.3% weight, dominates it).

Caveat: the live 1e-6 run appears to have been uncapped while the 1e-3 snapshot
capped the three large groups at 5,000, so part of the 5.28x spread is unequal
EM sample size. Same conclusion either way — the cap and the floor were both
equalising the groups, and Result 5's flatness depends on that equalisation.

## v2.15c — "combined" disambiguated; Table 6 labels; US registry citation (2026-08-25)

Three different quantities were all being called "combined", and two passages
described the wrong one:

  - per-group targets      41,183 / 39,831 / 43,785 / 31,129  (dedicated registry)
  - their sum                                        155,928  (all four groups)
  - pooled model (Table 1)                           243,849  (one merged pool)
  - weighted mean (Table 6)                           42,567  (expected per patient)

Table 1's row is the pooled model (registry_size_targets.csv, ethnicity
'Combined'); it is now labelled "Combined pooled registry" instead of
"Weighted Average", which read as an average and is not one. Section 3.1 and
the Table 1 footnote had described it as sizing "a single shared registry to
serve all groups proportionally" — correct for that row, but the same sentence
was also attached to Table 6's 42,567, which is a mean of per-group targets and
neither a pool nor a total. All three are now named wherever they can be
confused, and the 155,928 / 243,849 figures are computed from the CSVs rather
than hardcoded.

Pooling costs 1.6x more than four dedicated registries for the same patients
(243,849 vs 155,928) because a merged pool dilutes each patient's own
haplotypes. Stated explicitly for the first time.

Section 3.5 no longer claims the flat sensitivity bars show a structural
property of haplotype diversity: any weighted mean of four values in the
31,000-44,000 band is flat, and under dedicated registries the population
weights never enter the per-group targets at all. The substantive finding — the
four targets lie within 1.4x of each other — is now the claim.

Table 6's "Ethnic weights" column carried pre-census literals (77/8/9/6) that
contradicted the prose two paragraphs above once SG_WEIGHTS was corrected. All
four rows are now derived from the same CSV rows the sizes come from.

Section 3.3 cited only the Israeli Ezer Mizion registry [9] for a claim about
multiethnic registries plural; adds Gragert [2] for the U.S. registry and names
which is which.

Deck (build_slides_v215c.py) updated to match: slide 11 footnote carries the
243,849-vs-155,928 comparison, slide 19 states the bars are a weighted mean and
promotes the 1.4x spread to the finding, and the Result 3 spoken note's stale
"77 per cent Chinese" is corrected. Both slides rendered and checked for
overflow.

Red re-marked after every rebuild: 65/147 paragraphs and 81 table cells differ
from v2.15; verified by an independent differ (0 differing-but-black, 0
identical-but-red) with all 8 figures intact.

## v2.15c deck — acknowledgements slide added (2026-08-25)

New slide 26, "Acknowledgements", inserted before the closing slide (now 27):
Caitlyn Ng (Ngee Ann Polytechnic) contributed all of the coding for this
project — the ingestion and cleaning pipeline, the EM haplotype phasing, the
coverage and registry-size model, the bootstrap, and the figures. The slide
also lists the data providers (BMDP, SCBB, HSA) and the foundation dataset
(Ng et al. 2022).

The per-slide reader-notes dict is keyed by slide number, so the closing
slide's entry moved from 26 to 27; its text was also stale ("additions marked
red ... 5 paragraphs") and now states the current 65 paragraphs / 77 cells.
Deck is 27 slides, 21 images. Rendered and inspected.

Not added to the manuscript — its Declarations section has no Acknowledgements
subsection, and a contributorship statement is a decision for the authors.

## v2.15c — 2026-08-25 — POPULATION WEIGHTS CORRECTED; SUPERSEDES v2.15b

New label for the manuscript and deck carrying every correction from this
session: build_report_v215c.py -> HLA_Registry_Size_CMIO_v2.15c.docx,
build_slides_v215c.py -> HLA_Registry_Size_CMIO_v2.15c_slides.pptx.
Red means "differs from v2.15" and now covers EVERY difference, citation
renumbering included: 65 paragraphs + 77 table cells. The marker previously
masked against normalised text, which hid renumbered citations; normalisation
is now used for paragraph ALIGNMENT only (so the withdrawn entries read as
clean deletes) while the mask is computed against the raw v2.15 text.
Verified by an independent differ that imports none of the marker's logic:
0 differing-but-black and 0 identical-but-red, for both paragraphs and table
cells; 8 figures intact.

SG_WEIGHTS corrected from 0.77/0.08/0.09/0.06 to the Census of Population 2020
resident population: Chinese 3,006,770, Malays 545,500, Indians 362,270,
Others 129,670 of 4,044,210 = 74.3/13.5/9.0/3.2 per cent. The four counts sum
exactly to the total. Ng et al. quotes the same split (13.5/9/3.2; its "73.3%
Chinese" looks like a typo for 74.3%), so the old weights matched neither
source. The old Chinese 0.77 appears to have come from the census HOUSEHOLD
table (77.4% of household reference persons), not the population table.
The constant was duplicated in registry_model.py, 06 and 13 — all three fixed.

Re-ran 04, 06 and 13 against the snapshot. What moved:
  - Per-group same-ethnicity targets: UNCHANGED (they do not use weights), so
    the 40,000-45,000 headline and every bootstrap CI stand.
  - Combined pooled registry, 10/10 95%: 236,906 -> 243,849 (+2.9%).
  - Cross-ethnic Malay is no longer censored at the search ceiling:
    >10,000,000 -> 1,391,093 at 95% 10/10, because a correctly weighted pool
    contains proportionally more Malay haplotypes. Still infeasible, but now a
    real number rather than a ceiling hit. Indian and Chinese also improve.
  - Result 5: scenario spread 2.92% -> 3.50%, so the manuscript's "varies by
    less than 3%" claim is now FALSE and has been rewritten to state 3.5%
    (42,567 down to 41,129 at the 95% target). The qualitative conclusion —
    the target is robust to demographic reweighting — survives.

The re-run was done with n_max pinned to the snapshot's 1e7. Current code uses
1e10, and re-running unpinned silently rewrote ten ">10,000,000" ceiling values
as 1e10 — unrelated drift since the snapshot was frozen. Pinning isolates the
weights change and preserves the manuscript's ">10,000,000" convention. That
drift is still there for whoever next rebuilds the snapshot.

The Malay/Indian access claim is now cited to [1] Ng, which states it directly
("recipients from the Malay and Indian communities experience particular
challenges ... attributed to an underrepresentation of their haplotypes in
local donor pools"), with [8] Aljurf and [19] Kollman kept for the general
pattern. This is what the withdrawn Lim reference had been carrying.

Prior CSVs kept as *_PRE_WEIGHTS_FIX.csv.bak in the snapshot.

### v2.15c slide 23 — saturation argument added (2026-08-25)

Slide 23 stress-tested the 5,000-donor EM cap but never answered the question
it invites: would a much larger cohort change the targets? The data was already
in em_convergence.csv, in a column the slide ignored. Haplotypes clearing the
0.1% floor: 143 at n=5,000, 140 at n=20,000, 136 at n=45,018 — it saturates,
and sits at 123-144 across all four groups despite a 12-fold span in group size.
At a fixed floor extra donors cannot create common haplotypes, only measure
them better, which is why N* plateaus by n~5,000.

The corollary is the load-bearing half: retained mass (51.7/52.9/40.2/35.6%) is
population structure, not sample size, so no recruitment scale lifts a
below-floor patient above it — mismatch tolerance is the only lever for them.
The question is about the INPUT dataset, not the registry: had this analysis
started from 500,000 typed CMIO donors rather than 59,186, would the targets
differ? The convergence run answers it by re-estimating from subsamples of the
Chinese group. N* at 95%: 482,681 (n=500), 59,324 (n=2,000), 45,148 (n=5,000),
44,326 (n=20,000), 41,787 (n=30,000), 41,727 (n=45,018). The last step moves it
0.14% — the Chinese estimate is converged, and 500,000 would not change it.

The caveat is the other three groups. Malay (5,578), Indian (5,490) and Others
(3,767) sit where the Chinese curve was still 8-40% above its converged value,
so their reported targets are probably overestimates. That is the quantitative
content behind calling them exploratory: a larger CMIO dataset would sharpen
those three materially and the Chinese figure hardly at all.

Presented as a six-row table under Figure S1 with a "vs full sample" column,
a red bullet, and both notes layers. Slide rendered and inspected: figure
shrunk to 6.3in wide to make room, table clears the caption, no overlap.

Also corrected on that slide: the floor-sensitivity note quoted the Chinese 95%
target as 42,871, the EM point estimate, where the rest of the document reports
the bootstrap median 41,183.

### v2.15c addendum — census source pinned (2026-08-25)

The census PDF was fetched from the URL supplied by the user and is byte-identical
(md5 bc290a1c...) to the copy already in Publications/. It is Statistical Release 2
(Households, Geographic Distribution, Transport and Difficulty in Basic
Activities) — NOT Release 1, which is where ethnic composition would normally
sit. The CMIO figures are valid regardless: they come from its annex table
"Resident Population by Planning Area/Subzone of Residence, Ethnic Group and
Sex", Total row — Chinese 3,006,770, Malays 545,500, Indians 362,270, Others
129,670 of 4,044,210, which sums exactly.

Reference [10] now names Release 2 explicitly and carries the direct PDF URL;
Section 3.5 quotes the four counts so a reader can check the weights without
hunting. PDF filed as Publications/10_SingStat_Census2020_SR2.pdf. The old
cop2020-sr1 landing-page URL in DOWNLOAD_LINKS.txt was never verified and is
replaced. 17 of 19 references now have a PDF; only [04] Beatty and [16] Efron
(a book) are outstanding.

## v2.7.7 — 2026-08-25 (v2.15b) — CLAIM-LEVEL CITATION AUDIT

Every citation checked against the cited paper's full text (118 claim instances,
20 references, 17 PDFs). Full record: paper_BMT_workdir/citation_claim_audit_2026-08-25.md

CLINICAL CLAIMS CORRECTED — these were the serious ones:

- Recommendation 5 asserted "the clinical benefit of routine DQB1 typing is
  well-established [7]". Lee 2007 states the OPPOSITE: "Mismatching at HLA-DP or
  -DQ loci ... were not associated with survival" once 8/8 was achieved. The
  recommendation now rests on cost (DQB1 is nearly free given DRB1-DQB1 LD) and
  says plainly that withholding it is not known to cost lives.
- Recommendation 3 asserted "survival outcomes are comparable between 10/10 and
  permissive 9/10 matches [11,13]". Neither paper supports that. What
  Fleischhauer shows is narrower: a 10/10 with a NON-permissive DPB1 mismatch
  had outcomes not substantially different from a 9/10 with a PERMISSIVE one,
  while non-permissive mismatches raised mortality within both strata. Pidala
  found no significant DPB1 effect within 9/10 cases. Rewritten to state that a
  9/10 donor is not clinically equivalent, that permissiveness matters as much
  as the count, and that for patients with no 10/10 donor the comparison is
  against no transplant rather than against a full match.

REFERENCE WITHDRAWN — former [14] Anasetti 2012 (PBSC vs bone marrow) was cited
once, for a claim about centres accepting partial matches when time is critical.
It is a graft-source trial with no matching-level content. Citation removed and
the reference withdrawn; 15-21 shifted down to 13-19. PDF retained as
Publications/NOT_CITED_Anasetti_PBSCvsBM_NEJM2012.pdf.

CITATIONS RE-ATTRIBUTED to this manuscript's own results, where the cited paper
had no such content: the log-scale binary search and the same/cross-ethnic model
comparison (were [2,5]); the rare-diplotype tail (was [2]); the HWE departures in
Indian and Others (was [1] — Ng et al. reports no HWE testing); the 564 HSA
patients (was [1] — the count is ours, and is correct); the EM "more concentrated
distribution" claim (was [3,5], now cites [18] for EM and [3] for validation);
ancestry questions at registration (was [14] AFND, a database paper).

SCOPE NARROWED where a paper was stretched past what it covers: the
under-representation claim now states the generic finding Aljurf actually makes
and carries [19] Kollman, which says minority patients are less likely to find a
donor, with Singapore specifics moved to our own donor counts; "Asia-Pacific
region" softened to match Aljurf's worldwide framing; [6] Passweg moved off the
match-difficulty claim (it is an activity survey with no ethnicity content) onto
the transplant-volume statement it does support.

mark_v217_diffs.py: CITE_DROPPED/CITE_RENUMBER now compose both withdrawals
([8] and [14]) so renumbering still does not paint red.

Result: 19 references, all cited, none dangling. Red = 27 paragraphs + 60 cells.

## v2.7.6 — 2026-08-25 (v2.15b) — REFERENCE [08] WITHDRAWN; IDENTIFIER AUDIT

Reference [08] "Lim YA, Teo D, Ang AL, et al. HLA allele and haplotype
frequencies in unrelated bone marrow donor registries in Asia. Transpl
Immunol. 2010;22(3-4):166-174" could not be found in PubMed under any form of
its title, and Transpl Immunol 2010;22(3-4) runs ...165-71, 172-8, 179-83 —
there is no 166-174 slot. Treated as unverifiable and removed.

Its four citations rested on claims [09] Aljurf already supported, so each
now reads [08] Aljurf alone (or [1] alone for the haplotype-diversity
sentence). References 09-21 shifted down to 08-20.

Renumbering required two passes: cite() calls (74 rewritten) and 27 bracket
citations hardcoded in prose strings, which the first pass missed — the
abstract's "[1,11]" had silently come to mean Fleischhauer instead of the
Census. Verified afterwards that every surviving reference number still
resolves to the same paper it did in v2.15.

PMIDs — all 18 journal references checked against PubMed esummary (title,
journal, volume, pages must all match). 12 of the identifiers previously in
Publications/DOWNLOAD_LINKS.txt resolved to unrelated papers; the old Price
PMID returned a study of the avian MHC, the old Petersdorf PMID a paper on
dietary trans-fatty acids. Corrected: Nunes 24738646, Maiers 17869653,
Passweg 29540849, Lee 17785583, Aljurf 30778127, Halagan 28396163,
Fleischhauer 22340965, Pidala 25161269, Gonzalez-Galarza 31722398,
Price 10319267, Petersdorf 23878143.

[11] Fleischhauer title corrected to the published spelling,
"unrelated-donor haemopoietic-cell transplantation".

mark_v217_diffs.py — normalize_refs now also maps the old reference numbering
forward, so renumbering alone does not paint red; without it 51 of 148
paragraphs reddened. Normalisation moved to where the old text is read
(it is not idempotent — applying it twice would shift numbers twice), which
also lets the paragraph aligner see the withdrawn entry as a clean delete.
Red is back to 17 paragraphs + 60 cells, all substantive.

Publications/ renumbered to match (17 of 20 references now have a PDF;
[04] Beatty, [10] Census and [17] Efron outstanding). Nothing deleted.

## v2.7.5 — 2026-08-25 — REFERENCE IDENTIFIER AUDIT

Reference [13] was a mangled citation and is corrected in the manuscript.
Journal, volume, issue and pages (Tissue Antigens 2003;62(4):296-307) were
right, but the author list and subtitle belonged to no real paper:

  was:  Klitz W, Gragert L, Maiers M, Byard PJ. "New HLA haplotype frequency
        reference standards: five-locus haplotypes, and allele frequencies
        for 66 North American populations."
  now:  Klitz W, Maiers M, Spellman S, et al. "New HLA haplotype frequency
        reference standards: high-resolution and large sample typing of
        HLA DR-DQ haplotypes in a sample of European Americans."

Gragert L and Byard PJ are not authors of it; the real author list is Klitz,
Maiers, Spellman, Baxter-Lowe, Schmeckpeper, Williams, Fernandez-Vina
(PMID 12974796, verified via PubMed esummary; it is also cited as ref 3 of
Maiers 2007). The correction is apt: [13] is cited only for DRB1-DQB1 linkage
disequilibrium, and the real paper is specifically about DR-DQ haplotypes.

Publications/DOWNLOAD_LINKS.txt — identifiers audited against PubMed:
  [02] PMID 24989227 -> 25054717; DOI NEJMsa1311854 -> NEJMsa1311707
       (both were wrong and neither resolved); PMC5965695 added, so it is
       freely available after all.
  [04] PMID 7878796 -> 7482734 (7878796 is Burrows et al, Transplant Proc
       1995 — a different paper).
  [13] PMID 14617064 -> 12974796 (14617064 is a plant-biology paper on leaf
       senescence).
  [20] and [21] added — the file had no entries for them.

Publications/ PDFs renamed to NN_Author_ShortTitle_JournalYear.pdf, NN being
the manuscript reference number; each verified by reading its first page.
The file previously named 05_Gragret_HLA_Likelihood_NEJM2014.pdf is now
02_Gragert_MatchLikelihoods_NEJM2014.pdf and holds the real NEJM paper.
16 of 21 references now have a PDF; [04], [08], [11], [13], [18] outstanding.
Nothing was deleted.

## v2.7.4 — 2026-08-25 (v2.15b docx) — CITATION AUDIT AGAINST THE SOURCE PDFs

Verified reference [20] (Excoffier–Slatkin) against the primary literature:
Gragert et al. NEJM 2014 cites it directly as its own ref 20 ("haplotype
frequencies ... calculated ... with the use of the expectation-maximization
algorithm.15,16,20"), as does Maiers/Gragert/Klitz 2007 (our [5], its ref 6).
The method citation is the field standard, not an outlier. Note Gragert used a
STAGED EM (C-B and DRB1-DQB1 clusters, then a three-locus run) on four-locus
adult-donor haplotypes; this study does full five-locus phase enumeration.

Three miscitations corrected — Gragert 2014 is a four-locus 8/8/7/8 paper and
contains no "9/10" or "10/10" anywhere, so it cannot support a 9/10 claim:

- Abstract: "Relaxing to 9/10 ... halves the required registry size [14]" — [14]
  is Anasetti (PBSC vs bone marrow), which has zero content on matching level.
  Citation removed; the result is this study's own (Section 3.4).
- Discussion + Recommendation 3 + Conclusions: [2] / [2,19] re-pointed to
  Section 3.4 for the halving itself. Recommendation 3 now cites [2] correctly,
  as the analogous U.S. precedent for reporting 7/8 alongside 8/8. Conclusions
  retains [19] (Petersdorf) for the clinical meaning of a single mismatch —
  verified to use 9/10 terminology, though it says nothing about registry size.

Added reference [21] Kollman et al., Transplantation 2004;78(1):89–95
("Assessment of optimal size and composition of the U.S. National Registry of
hematopoietic stem cell donors") — cited by Gragert as their ref 18, and the
direct predecessor to this paper's question. Cited in Section 2.3 alongside
Beatty [4]. Appended as [21] so the existing twenty keep their numbers.

All 21 references now cited; no dangling citations and no uncited entries.
Red marking after rebuild: 15 paragraphs + 60 table cells.

Also: Publications/05_Gragret_HLA_Likelihood_NEJM2014.pdf previously held a
duplicate of Lee et al. Blood 2007; replaced with the actual NEJM paper.

## v2.7.3 — 2026-08-25 (v2.15b docx + slides) — BOOTSTRAP n_eff MATCHED TO THE EM SAMPLE

Four related defects, all inherited from the June `16bd246` snapshot, which was
frozen between two settings that the August pass later made consistent
(`hwe_test.py` EM cap 5,000 -> 50,000; `09_bootstrap_ci.py` N_EFF corrected).

1. **n_eff overstated.** The Dirichlet used the full 5-locus donor count while
   EM had seen only the capped sample: Chinese frequencies came from 5,000
   individuals but the interval asserted the precision of 44,400. Re-ran
   `09_bootstrap_ci.py` against the snapshot with `n_eff = min(cap, count)` =
   5,000 / 5,000 / 5,000 / 3,767. Chinese 95% CI widens ~4x and the Jensen
   correction deepens: 42,847 (42,649-43,058) -> **41,183 (40,184-42,153)**.
   Malay, Indian and Others move <1% (their caps barely bound).
2. **61,149 arithmetic.** Section 2.4 quoted n_eff values (45,754/5,868/5,586/
   3,941) summing to 61,149, above the 59,186 total stated in Section 2.1.
   True 5-locus counts, recomputed from the snapshot `hla_clean.csv`, are
   44,400/5,578/5,490/3,767 = 59,235.
3. **Section 3.1 causal claim.** Attributed the narrow Chinese CI to "the
   largest 5-locus sample (45,754)". Under the cap all three major groups fed
   EM the same 5,000 individuals, so the stated cause was contradicted by
   Section 2.3. Now reports comparable widths (+/-800-1,100) and explains why.
4. **Section 3.6 cap claim.** Said Malay/Indian/Others were uncapped ("sample
   sizes <= 5,868"). Malay (5,578) and Indian (5,490) both exceed the 5,000
   cap; only Others (3,767) does not.

No recommendation changes: 41,183 remains inside R1's 40,000-45,000 band and
every comparative result is unaffected. Figure 2 (`registry_ci_plot.png`)
regenerated. Previous CI table kept as
`analysis/snapshot_1e-3/data/registry_size_ci_PRE_NEFF_FIX.csv.bak`.

Red marking after rebuild: 12 paragraphs + 60 table cells (was 5 + 0) — the
corrected values genuinely differ from v2.15, so the marker paints them.
Deck: both results tables corrected and their changed cells marked red
(21 and 20 cells), slide-10 bullet and reader note rewritten.

## v2.7.2 — 2026-08-24 (v2.15b docx) — v2.15 MADE DEFENSIBLE BY ADDITION ONLY

`HLA_Registry_Size_CMIO_v2.15b.docx`, built by `build_report_v215b.py` — the
v2.15 script (commit 792fa1a) verbatim, reading the frozen 1e-3 snapshot, plus
the minimum text that makes it defensible. Tables, figures, recommendations,
and all v2.15 prose untouched; none of v2.17's extra table columns, Table 1b,
or §3.1.1/§4.2 sections.

Additions (exactly 3 paragraphs, marked red; all other text black):
1. §2.2 — one sentence: post-floor frequencies renormalised; coverage
   conditional on two haplotypes above the floor; see §4.1.
2. §4.1 Limitations — the conditional estimand: retained mass
   51.7/52.9/40.2/35.6%, unconditional ≈ one-third to one-half, the observed
   41.8% Chinese twin rate vs 49.1% conditional-equivalent, and rare-haplotype
   patients routed to partial-match (Recommendation 3).
3. §4.1 Limitations — floor sensitivity: 1e-6 floor → 8.7×10⁷; absolute
   targets are order-of-magnitude guidance, comparative findings robust.

Rationale (user decision): v2.15 is the simpler, more readable manuscript and
its conclusions are essentially those of v2.17; it need not be perfect, only
defensible on its own terms. The six recommendations survive with these
caveats (Recommendation 6's equity argument transfers to partial-match
tolerance — the rare-haplotype patients sit below the floor at any threshold).

Marking: `mark_v217_diffs.py` now takes optional OLD NEW arguments;
`python mark_v217_diffs.py HLA_Registry_Size_CMIO_v2.15.docx
HLA_Registry_Size_CMIO_v2.15b.docx` resets all runs to black (including
v2.15's legacy hand-set red) and paints only the diff. Verified: 3/148
paragraphs and 0 table cells red; body colours exactly {000000, C00000};
page 11 rendered and inspected.

**Same-day follow-ups.** (1) Cross-references now use "Section n" instead of
"§n" in both v2.15b and v2.17 (0 § remain in either docx); the diff marker
normalises §→Section before comparing, so pure restyling of inherited v2.15
text stays black. (2) Recommendation 6 restated additively: one red sentence
noting both coverage thresholds are conditional on the frequency floor and
that below-floor patients depend on the partial-match protocols of
Recommendation 3. v2.15b now carries exactly 4 red paragraphs (Section 2.2
sentence, Recommendation 6 sentence, two Limitations paragraphs); 0 red table
cells.

**v2.15b slides (same day).** `build_slides_v215b.py` →
`HLA_Registry_Size_CMIO_v2.15b_slides.pptx` (26 slides), adapted from the
v2.17 deck: version strings v2.15b; slide 12 reframed from
"correction to the previous version" to "the most important caveat"
(matching v2.15b, where this lives in Limitations); Recommendation 3 names
MMUD explicitly (a 9/10 donor is an MMUD; only route for below-floor
patients; safety from external trials — speaker note carries the full form);
Recommendation 6 restated per the manuscript (neither 90% nor 95% reaches
below-floor patients; their route is Recommendation 3); limitations slide
gains the frequency-floor bullet (conditional coverage; 1e-6 → 8.7×10⁷).

**Reader notes (same day).** Every v2.15b slide's notes now carry two
paragraphs: the speaking script, then an "Additional detail —" paragraph for
readers working through the deck without the talk (manuscript
section/table/figure cross-references, exact numbers, and scope caveats;
verified present on all 26 slides). Also fixed the last stale speaker note:
slide 25's take-home still promised "95 per cent of patients a fully matched
donor" — now states 95% of common-haplotype patients ≈ half of all patients,
41.8% observed.

**Deck figure-source bug found in passing:** both deck scripts read the live
`analysis/figures/`, which the 1e-6 re-run overwrote — slide 22/23's
em_convergence.png showed N* rising to 8.7e7 under a caption claiming
stability near 45k. Both `build_slides.py` (v2.17) and `build_slides_v215b.py`
now read `analysis/snapshot_1e-3/figures/`; all 9 deck figures exist in the
snapshot; both decks rebuilt and the limitations slide re-rendered and
inspected (curve now stabilises at ~45k, consistent with its caption).

**MMUD sentence (same day).** Recommendation 3 gains one red sentence drawing
the MMUD conclusion already implicit in v2.15b's own results: a 9/10 unrelated
donor is by definition an MMUD, so MMUD capability halves every group's
effective registry requirement and is the only access route for below-floor
patients (Section 4.1); safety rests on the external trial literature [12,15].
No new calculations. v2.15b now carries 5 red paragraphs, 0 red table cells,
8 figures.

**Marker bug fix (same day).** `mark_v217_diffs.py`'s picture guard used
`findall()` (direct children of w:r) and missed drawings nested at
w:r/w:drawing/wp:inline — figure paragraphs were stripped as empty text, so
both marked docx files silently lost all 8 embedded figures. Guard now uses
`iter()` over descendants and covers anchored drawings too; both documents
rebuilt and re-marked. Verified: 8 inline images in each, red counts
unchanged (v2.15b 4/148 + 0 cells; v2.17 37/162 + 62 cells), v2.15b page 7
rendered — Figure 2 present above the untouched v2.15-style Table 2.

Three manuscripts now stand: v2.15b (simple, additive caveats), v2.17 (full
conditional/unconditional reporting), v2.16 (1e-6 methodological).

## v2.7.1 — 2026-08-24 (v2.17 docx + v2.17 slides) — CONSISTENCY PASS + RED CHANGE-MARKING

Text and colour only. No analysis script re-run, no figure regenerated. Full
audit trail: `paper_BMT_workdir/v217_consistency_audit.md`.

**Cohort count — 59,186 retained, deliberately.** Three totals were in
circulation. Ground truth from `analysis/data/hla_clean.csv`:

| definition | value |
|---|---|
| Ng et al. 2022 published cohort, after inclusion criteria | **59,186** |
| BMDP_OUT 57,785 + SCBB_OUT 1,450 on re-ingest (all 5-locus complete) | 59,235 |
| + HSA-Donor 1,350 + HSA-Patient 564 | 61,149 |

An earlier pass here changed the headline to 59,235 on the reading that 59,186
was a three-source (BMDP+SCBB+HSA) total. That reading was wrong: 59,186 is
*below* BMDP+SCBB alone, so HSA cannot be inside it. 59,186 is BMDP+SCBB after
the published inclusion criteria, and the 49-record (0.08%) gap is what those
criteria remove. This study is a direct follow-up on the identical dataset, so
the published figure is what it reports. §2.1 now states the 59,235 re-ingest
count and the 49-record difference once, explicitly, and confirms the HSA
records enter only §3.5 and §3.6.

**Regression from the v2.17 branch point.** `build_report_v217.py` was branched
from `16bd246`, which predates `792fa1a`. That commit's fix — Others donor count
3,941 → 3,847, and removal of the phantom silhouette=0.97 — had been silently
undone in §3.7 and the Figure 7 caption. Its wording is grafted back verbatim.
(§4.1 Limitations already carried 0.24; the regression was confined to §3.7/Fig 7.)

**Also corrected**
- §3.1.1 caption: 42,871 labelled point estimate, alongside Table 1's 42,847
  bootstrap median. Both were always correct; neither was labelled.
- Table 1 footnote †: "Weighted Average" → "Combined pooled registry", matching
  the table body renamed at v2.17.
- §3.5: notes 77/8/9/6 is the resident-citizen breakdown, against the
  73.3/13.5/9.0/3.2 citizens-plus-PR figures cited by Ng et al.

**Verified consistent, unchanged:** Table 4 registry weights (75/9/9/6 vs actual
75.0/9.4/9.3/6.4), HSA scenario weights (72/15/5/8 vs actual 72.0/14.9/5.0/8.2),
§3.6's 564 HSA pairs, all Table 1 / Table 4 figures.

**Prose pass (same day).** The v2.17-new passages were rewritten into
manuscript register after review flagged them as reading machine-generated:
the abstract's labelled "Estimand:" block folded into Methods; ALL-CAPS
emphasis (CONDITIONAL/ALL/WITHIN) lowercased; editorial meta-commentary
removed ("the single most important check in this report", "companion
methodological analysis", "worth far more than any recruitment programme that
could plausibly be funded", "stark finding", code-file line citations);
87,384,114 now reported as 8.7×10⁷. Content, numbers, and claims unchanged —
register only. 36/162 paragraphs red after the pass.

**Red now means one thing: differs from v2.15.** `mark_v217_diffs.py` (new) is a
post-build pass that resets every run to black, then paints red the paragraph-
and word-level differences against `HLA_Registry_Size_CMIO_v2.15.docx`. It
replaces the hand-set red of the old `add_corrected_para` convention, which by
v2.17 marked v2.14-era corrections and no longer meant anything a reader could
name. Idempotent — colour is always recomputed from the two documents.

    python build_report_v217.py && python mark_v217_diffs.py

Result: 37/163 paragraphs and 62 table cells red; body text is exactly two
colours, 000000 and C00000. Table header rows stay white-on-navy (red on a dark
fill is unreadable, and the header change shows in the body cells beneath).
Tables are paired by whole-table text similarity, not by header — the 10/10 and
8/8 tables share a header, and header-matching floods the 8/8 table with false
differences.

**Slides — `HLA_Registry_Size_CMIO_v2.17_slides.pptx`, 26 slides.**
Rebuilt from `build_slides.py`; the v2.15 pptx is left untouched on disk. Slide 6
had the worst defect in either artifact: headline 59,186 over a table of
45,754/5,868/5,586/3,941 (per-locus, HSA included) under the header "Donors with
full 5-gene typing". Now the 5-locus counts.

The deck also still presented the conditional 95% as unconditional, which is what
v2.17 exists to correct. Brought into line:
- new slide 12 "Reality check" carrying Table 1b — observed 41.8% against
  modelled 49.1% — with speaker notes;
- slide 11 subtitle and headline bullet marked CONDITIONAL;
- recommendation 1 and the take-home bullet now state ~49% full-match coverage;
- Figure 7 caption and its bullet drop "cleanly separated" for s=0.24.

Slides 6 and 12, and docx pages 1–3, rendered and inspected: no overflow or
overlap; red marking lands where intended.

## v2.7.0 — 2026-08-21 (v2.17 docx) — v2.15 MODEL, CORRECTLY REPORTED

**There are now two manuscripts, answering two different questions. Neither
supersedes the other.**

| | v2.17 (`build_report_v217.py`) | v2.16 (`build_report.py`) |
|---|---|---|
| Floor | 1e-3 (v2.15 model, unchanged) | 1e-6 |
| Estimand | conditional on both haplotypes > floor | unconditional |
| Chinese 95% 10/10 | 42,847 | 87,384,114 |
| Intervals | bootstrap percentile ranges, retained | none (withdrawn) |
| Audience | clinical / policy | methodological |
| Data source | `analysis/snapshot_1e-3/` (frozen from 16bd246) | `analysis/data/` (live) |

v2.17 changes **nothing** about the v2.15 computation. Every registry size,
figure, and bootstrap range is the number v2.15 produced. What changed is that
the estimand is now declared, and the dropped denominator restored.

**The diagnosis.** `04_registry_model.py:139` renormalises the post-floor
haplotype frequencies to sum to 1. The v2.15 EM table sums to 0.5172 (Chinese),
0.5288 (Malay), 0.4024 (Indian), 0.3563 (Others) — so 47–64% of frequency mass
was discarded and the remainder rescaled. v2.15's coverage figures are coverage
*of the retained subpopulation*, reported as if unconditional. The arithmetic was
never wrong; the denominator was.

**The check that settles it.** `paper_BMT_workdir/empirical_match.py` counts, with
no EM, no HWE and no floor, how many donors have an exact 10/10 genotype twin
inside the registry. All five loci are 2-field typed, so this is allele-level.
Chinese: **41.8% of 44,400**. v2.15 claimed a ~42,871-donor registry gives 95% —
Singapore *has* that registry and gets 41.8%. Restoring the denominator gives
0.95 × 0.5172 = **49.1%**, which matches. The residual gap is in the expected
direction (HWE assumes random mating; substructure raises real match rates).

**Applied in v2.17** (all text-only; no script re-run, no figure regenerated):
1. Estimand declared conditional; retained mass per group in §2.2 and Table 1
2. N* retitled throughout to name the conditioning
3. Table 1 gains "Mass retained" and "Uncond. at 95%" columns
4. New §3.1.1 + Table 1b: observed registry match rate vs model
5. §2.4 retitled "Bootstrap Percentile Ranges" — kept, with the E[f²] bias
   disclosed. Small here (0.6–1.5% displacement) *because* the 1e-3 floor leaves
   no rare haplotypes; the same procedure at 1e-6 is uninterpretable
6. §4.1: cross-group ranking withdrawn (floor bites unequally, 0.36–0.53 mass)
7. §4.1: floor named as the dominant uncertainty, 1e-6 sensitivity cited
8. New §4.2 + Table 9: MMUD implications

**Also corrected in passing:** "Weighted Average" → "Combined pooled registry"
(it is the pooled model, not an average); Others silhouette 0.97 → 0.24;
bootstrap concentration parameters used per-locus typed counts (61,149) rather
than 5-locus-complete counts (59,235), now disclosed in §2.4.

**New headline (§4.2).** At a fixed 50,000-donor registry, coverage by stringency:

| Group | 10/10 | ≥9/10 | ≥8/10 |
|---|---|---|---|
| Chinese | 33.6% | 67.9% | 92.0% |
| Malay | 40.4% | 70.1% | 91.8% |
| Indian | 21.3% | 53.6% | 85.4% |
| Others | 17.3% | 41.0% | 76.0% |

One permitted mismatch is worth a 19–81× registry expansion, and is worth most to
the groups with the poorest same-ethnicity prospects. The paper's conclusion moves
from "recruit N donors" to "recruitment cannot close the gap; mismatch tolerance
can". PTCy/haplo/cord are cited as external literature — this analysis contains no
outcome data and says so explicitly.

**Robustness of the comparative claims** (`paper_BMT_workdir/robustness.py`):
cross-ethnic N* exceeds same-ethnicity N* in all four groups under *both* floors,
so that conclusion is floor-independent. Magnitude is not: three of four v2.15
ratios were censored at the old `n_max=1e7` ceiling and are lower bounds only.
Reported qualitatively with ≥ signs.

**Reproducing v2.17:** `python build_report_v217.py`. It reads
`analysis/snapshot_1e-3/`, a frozen copy of `analysis/data` + `analysis/figures`
taken from commit 16bd246, so it stays reproducible even though `analysis/data/`
has since been re-run at 1e-6. `check_freshness.py` governs the 1e-6 pipeline and
does not apply to the snapshot.

---

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
