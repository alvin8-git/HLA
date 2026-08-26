# v2.17 + v2.15 slides — sample-count consistency audit (2026-08-24)


> **Historical record — describes the v2.17 line at a 1×10⁻⁶ floor.** Not the
> submission manuscript, which is v2.15c at 1×10⁻³. Figures here are three orders of
> magnitude larger by design. Left unchanged as the record of that audit; see
> [../docs/explanation-frequency-floor.md](../docs/explanation-frequency-floor.md)
> for how the two floors relate.

Ground truth from `analysis/data/hla_clean.csv` (live; identical ingest to snapshot):

| source | unique samples |
|---|---|
| BMDP_OUT | 57,785 |
| SCBB_OUT | 1,450 |
| **BMDP+SCBB (MAIN_SOURCES)** | **59,235** |
| HSA-Donor | 1,350 |
| HSA-Patient | 564 |
| **all sources** | **61,149** |

All 59,235 MAIN samples are 5-locus complete. Per ethnicity (5-locus):
Chinese 44,400 · Malay 5,578 · Indian 5,490 · Others 3,767 → sums to 59,235.
Per-locus (any locus, all sources): 45,754 / 5,868 / 5,586 / 3,941 → 61,149.

Ng et al. 2022 (`2022_HLA_BloodCellTherapy.pdf`, lines 78–85, 142):
> "HLA typing data from **BMDP, SCBB, and HSA** were combined… High-resolution
> results were obtained from **59,186** individuals." (→ 118,372 alleles = 2 × 59,186)

So 59,186 is the 2022 **three-source** inclusion count. It is not the two-source
cohort this study models.

## Findings

| # | Where | Says | Should be | Severity |
|---|---|---|---|---|
| 1 | v2.17 §Abstract, §1, §2.1, §5 (4×) | "59,186 donors … from BMDP and SCBB" | 59,235; and 59,186 was BMDP+SCBB+**HSA** in Ng 2022 | **High** — headline N contradicts Table 1b, which sums to 59,235; also misattributes the source set |
| 2 | v2.17 §2.1 | "Recipients and donors processed by HSA were included for validation" | Correct as written, but sits in the same sentence as the 59,186 total, implying HSA is inside it | **High** (same fix as #1) |
| 3 | v2.17 §4.1, Fig. 6 caption (3×) | "3,941 Others donors" clustered | **3,847** (`others_cluster_assignments.csv` = 3,847 rows; slides already say 3,847) | **High** — wrong number, contradicts own data file |
| 4 | v2.17 §3.1.1 caption | Chinese "modelled 95% target (42,871)" | 42,871 is the point EM estimate; Table 1/Abstract use bootstrap median **42,847** | Medium — both real, but unlabelled |
| 5 | v2.17 Table 1 footnote † | "Weighted Average:" | Table body was renamed "Combined pooled registry†"; footnote not updated | Low |
| 6 | v2.15 slides, slide 6 | headline "59,186 … Sources: BMDP + SCBB" over a table of 45,754/5,868/5,586/3,941 | table sums to **61,149** (per-locus, incl. HSA) — three different cohorts on one slide | **High** |
| 7 | v2.15 slides, slide 11/14 | "Weighted average†" | same rename as #5 | Low |

## Confirmed consistent (no action)

- Table 1b donor counts 44,400 / 5,578 / 5,490 / 3,767 — exact.
- §2.4 disclosure of 61,149 vs 59,235 — exact, and correctly names both as per-locus vs 5-locus-complete.
- Table 4 "BMDP+SCBB donors 75/9/9/6" — actual 75.0/9.4/9.3/6.4. ✓
- Table 4 "HSA Patient-Donor Data 72/15/5/8" — actual HSA-Patient 72.0/14.9/5.0/8.2. ✓
- §3.6 "564 patient-donor pairs from HSA" — HSA-Patient = 564 exactly. ✓
- Slide 20 "3,847 fully-typed Others donors" — correct (v2.17 is the one that drifted).
- Slide 11 targets match Table 1 exactly (42,847 / 40,032 / 43,855 / 31,181 + CIs).
- Slide 18 scenario targets match Table 4 exactly.

## Note, not a defect

v2.17 §3.5 uses Singapore weights 77/8/9/6 citing [11]. Ng 2022 cites the 2020
census as 73.3 / 13.5 / 9 / 3.2 (citizens+PR vs citizens-only). Both are
defensible; a reviewer comparing the two papers will ask. One clause naming the
denominator would close it.

## Recommended minimal fix

Replace the 4 instances of "59,186" with 59,235 and add one clause to §2.1:

> "…from 59,235 donors and cord blood units from BMDP and SCBB. Ng et al. report
> 59,186 individuals for the same accrual; that figure additionally includes the
> 1,914 HSA donor/recipient records, which are used here only for the demographic
> scenario in §3.5 and the external validation in §3.6, not in frequency estimation."

Plus: 3,941 → 3,847 (3×), label 42,871 as the point estimate, fix the two
footnotes, and correct slide 6.

---

## RESOLVED — 2026-08-24

All seven findings fixed in `build_report_v217.py` and `build_slides.py`; both
artifacts rebuilt. Text-only: no analysis script re-run, no figure regenerated.

Root cause of #3: `build_report_v217.py` was branched from `16bd246`, which
predates `792fa1a` ("correct phantom silhouette 0.97 and Others donor count").
v2.17 had silently reintroduced both the 3,941 donor count and the phantom
silhouette=0.97 in §3.7 and the Figure 7 caption. The 792fa1a wording has been
grafted back in verbatim. (v2.17's §4.1 Limitations already carried 0.24 — the
regression was confined to §3.7/Fig 7.)

Verification after rebuild — all pass:

- docx: 59,235 ×4; 59,186 ×2 (both correctly attributed to Ng et al.); 3,847 ×2;
  3,941 ×1 (correct, §2.4 per-locus bootstrap count); "silhouette=0.97" ×0;
  "Weighted Average" ×0.
- pptx: 45,754/5,868/3,941 ×0; 44,400+5,578+5,490+3,767 = 59,235 ✓;
  "Weighted average" ×0; "cleanly separated" ×0; v2.15 strings ×0.
- Slides 6 and 12 rendered to PNG and inspected: no overflow or overlap.

Slides additionally brought in line with v2.17's central correction (they still
presented the conditional 95% as unconditional):

- new slide 12 "Reality check", carrying Table 1b (observed 41.8% vs modelled
  49.1%) with speaker notes;
- slide 11 subtitle and headline bullet now say CONDITIONAL;
- recommendation 1 and the take-home bullet now state ~49% full-match coverage,
  not 95%;
- output renamed `HLA_Registry_Size_CMIO_v2.17_slides.pptx` (26 slides).
  `HLA_Registry_Size_CMIO_v2.15_slides.pptx` is left untouched on disk.

---

## SUPERSEDED — finding #1 reversed, 2026-08-24

Finding #1 above concluded 59,186 was a three-source (BMDP+SCBB+HSA) total and
should become 59,235. **That conclusion was wrong**, on arithmetic already in
this file: 59,186 < 59,235, so HSA cannot be inside it. The 2022 sentence
"HLA typing data from BMDP, SCBB, and HSA were combined and checked for
discrepancies" describes the QC step, not the cohort feeding the 59,186.

Correct reading:

    BMDP 57,785 + SCBB 1,450          = 59,235  (re-ingest)
    less 49 records failing the 2022 inclusion criteria
                                      = 59,186  (published)
    HSA 1,350 + 564 = 1,914, held out for §3.5 / §3.6

This work is a direct follow-up on the identical dataset, so it reports the
published 59,186. §2.1 states the 59,235 re-ingest count and the 49-record
(0.08%) difference once, explicitly. Table 1b keeps the true per-group analysed
counts (44,400 / 5,578 / 5,490 / 3,767), which sum to 59,235 — the §2.1 sentence
is what reconciles that with the headline.

Findings #2–#7 stand as resolved.
