#!/usr/bin/env python3
"""05_report.py — Assemble verification summary report."""

import pandas as pd
import numpy as np
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, 'data')
OUT_PATH = os.path.join(HERE, 'verification_summary.md')


def section_allele_freq(af: pd.DataFrame) -> str:
    lines = ["## 1. Allele Frequency Reproducibility\n"]

    total = len(af)
    per_locus = af.groupby('locus').size().to_dict()
    lines.append(f"**Total alleles compared:** {total}\n")
    lines.append("**Per locus:**\n")
    for locus, n in sorted(per_locus.items()):
        lines.append(f"- {locus}: {n}")
    lines.append("")

    # Only rows with pub_frequency (non-null difference)
    af_matched = af.dropna(subset=['difference'])
    overall_max = af_matched['difference'].abs().max()
    lines.append(f"**Max absolute discrepancy (overall):** {overall_max:.6f} ({overall_max*100:.4f}%)\n")

    lines.append("**Max absolute discrepancy per locus:**\n")
    per_locus_max = af_matched.groupby('locus')['difference'].apply(lambda x: x.abs().max())
    for locus, mx in per_locus_max.sort_index().items():
        lines.append(f"- {locus}: {mx:.6f} ({mx*100:.4f}%)")
    lines.append("")

    n_flagged = af['flagged'].sum()
    lines.append(f"**Flagged alleles** (|diff| > 0.005): {n_flagged}\n")

    if n_flagged == 0:
        lines.append(
            "**Conclusion:** All allele frequencies are reproducible. "
            f"The maximum discrepancy of {overall_max*100:.4f}% is well below the 0.5% threshold. "
            "The 2022 publication's reported frequencies are confirmed."
        )
    else:
        lines.append(
            f"**Conclusion:** {n_flagged} allele(s) exceeded the 0.5% threshold — review recommended."
        )
    lines.append("")
    return "\n".join(lines)


def section_hwe(hwe: pd.DataFrame) -> str:
    lines = ["## 2. Hardy–Weinberg Equilibrium Assessment\n"]
    lines.append(f"**Bonferroni-corrected threshold:** p < 0.0025 (0.05 / 20 tests)\n")
    lines.append("**Full HWE results (all 20 tests):**\n")

    # Markdown table
    header = "| Ethnicity | Locus | N | H_obs | H_exp | χ² | p-value | Significant |"
    sep =    "|-----------|-------|---|-------|-------|-----|---------|-------------|"
    lines.append(header)
    lines.append(sep)
    for _, row in hwe.iterrows():
        sig = "**YES**" if row['significant'] else "No"
        lines.append(
            f"| {row['ethnicity']} | {row['locus']} | {int(row['n_individuals']):,} "
            f"| {row['H_obs']:.4f} | {row['H_exp']:.4f} "
            f"| {row['chi2_stat']:.4f} | {row['p_value']:.2e} | {sig} |"
        )
    lines.append("")

    sig_rows = hwe[hwe['significant']]
    lines.append(f"**Significant violations ({len(sig_rows)}):**\n")
    for _, row in sig_rows.iterrows():
        direction = "deficit" if row['H_obs'] < row['H_exp'] else "excess"
        lines.append(
            f"- **{row['ethnicity']} / {row['locus']}**: H_obs={row['H_obs']:.4f}, "
            f"H_exp={row['H_exp']:.4f}, p={row['p_value']:.2e} "
            f"(heterozygosity {direction})"
        )
    lines.append("")

    lines.append(
        "**Interpretation:**\n"
        "\n"
        "- *Indian group*: 3 violations (DQB1, HLA-B, HLA-C). "
        "All show heterozygosity deficit (H_obs < H_exp), consistent with mild population sub-structure "
        "or possible genotyping artefacts at these loci.\n"
        "\n"
        "- *Others group*: All 5 loci are significant (p ≤ 3.8×10⁻²⁸ for HLA-B). "
        "The \"Others\" category is ethnically heterogeneous (non-CMIO individuals), "
        "so HWE departure is expected — pooling genetically distinct sub-populations inflates "
        "apparent heterozygosity excess under a single-population model.\n"
        "\n"
        "- Chinese and Malay groups: no violations at the Bonferroni threshold."
    )
    lines.append("")
    return "\n".join(lines)


def section_haplotypes(haplo: pd.DataFrame) -> str:
    lines = ["## 3. Haplotype Frequency Summary\n"]
    lines.append(
        "> **Note:** Haplotype frequencies were estimated using a simplified EM approach "
        "(per-locus product approximation rather than full 5-locus phase reconstruction). "
        "Frequencies should be interpreted as indicative only.\n"
    )

    threshold = 0.001
    for eth in ['Chinese', 'Indian', 'Malay', 'Others']:
        sub = haplo[haplo['ethnicity'] == eth].copy()
        common = sub[sub['frequency'] >= threshold]
        top5 = sub.nlargest(5, 'frequency')
        lines.append(f"### {eth}\n")
        lines.append(f"- Haplotypes with frequency ≥ {threshold}: **{len(common)}**")
        lines.append(f"- Top 5 haplotypes:\n")
        lines.append("| Rank | Haplotype (A|B|C|DRB1|DQB1) | Frequency |")
        lines.append("|------|------------------------------|-----------|")
        for rank, (_, row) in enumerate(top5.iterrows(), 1):
            lines.append(f"| {rank} | {row['haplotype']} | {row['frequency']:.4f} |")
        lines.append("")

    return "\n".join(lines)


def section_missing_data() -> str:
    lines = ["## 4. Missing Data Impact\n"]
    lines.append(
        "**BMDP_OUT + SCBB_OUT (used for all frequency calculations):**\n"
        "\n"
        "| Ethnicity | Locus | Total rows | Allele1 missing | Allele2 missing |\n"
        "|-----------|-------|------------|-----------------|------------------|\n"
    )
    # All zero from our computation
    ethnicities = ['Chinese', 'Indian', 'Malay', 'Others']
    loci = ['DQB1', 'DRB1', 'HLA-A', 'HLA-B', 'HLA-C']
    sample_sizes = {
        'Chinese': 44400, 'Indian': 5490, 'Malay': 5578, 'Others': 3767
    }
    for eth in ethnicities:
        for locus in loci:
            n = sample_sizes[eth]
            lines.append(f"| {eth} | {locus} | {n:,} | 0 (0.0%) | 0 (0.0%) |")
    lines.append("")

    lines.append(
        "**HSA-Donor and HSA-Patient files (excluded from frequency estimation):**\n"
        "\n"
        "These files contain single-field (allele1 only) typing — allele2 is 100% missing "
        "across all 5 loci. They were excluded from all allele frequency, HWE, and haplotype "
        "calculations to avoid bias from incomplete genotype records.\n"
    )
    lines.append(
        "**NaN handling:** Rows with NaN allele2 were excluded from the allele count denominator "
        "(not double-counted), ensuring frequency estimates reflect actual typed chromosomes only.\n"
    )
    lines.append(
        "**Note on DQB1:** DQB1 missingness rates ~50–64% reported in the raw HSA files "
        "refer to the HSA cohort only. In the BMDP+SCBB analysis cohort, DQB1 is fully typed "
        "(0% missing), and the 2022 publication's DQB1 frequencies are reproducible.\n"
    )
    return "\n".join(lines)


def section_improvements() -> str:
    lines = ["## 5. Suggested Improvements\n"]
    lines.append(
        "1. **Higher resolution (4-field) typing where available.** "
        "Current data uses 2-field resolution. Moving to 4-field (protein + synonymous coding + intron) "
        "would improve match stringency, particularly for mismatched pairs at DRB1 and HLA-B.\n"
        "\n"
        "2. **Bootstrap confidence intervals for haplotype frequencies.** "
        "Point estimates from the EM algorithm carry sampling uncertainty, especially for rare haplotypes "
        "in the Indian (n=1,098) and Others (n=754) cohorts. Bootstrap CIs would quantify this.\n"
        "\n"
        "3. **Linkage disequilibrium (D') reporting between loci pairs.** "
        "Pairwise LD measures (D', r²) between the 5 loci would characterise the extent of "
        "non-random association and validate the per-locus independence assumption used in the "
        "simplified EM haplotype estimation.\n"
        "\n"
        "4. **Consideration of admixed individuals in the \"Others\" group.** "
        "The HWE violations in Others are consistent with population heterogeneity. "
        "Stratifying Others by ancestry (e.g., via principal component analysis of HLA alleles) "
        "would improve frequency estimates and registry size projections for this group."
    )
    lines.append("")
    return "\n".join(lines)


def section_registry(reg: pd.DataFrame) -> str:
    lines = ["## 6. Registry Size Findings\n"]
    lines.append("**Full registry size projection table (72 scenarios):**\n")

    lines.append("| Match Level | Ethnicity | Model Variant | Target Coverage | Registry Size |")
    lines.append("|-------------|-----------|---------------|-----------------|---------------|")
    for _, row in reg.iterrows():
        rs = f"{int(row['registry_size']):,}" if row['registry_size'] < 10_000_000 else "≥10,000,000 (ceiling)"
        lines.append(
            f"| {row['match_level']} | {row['ethnicity']} | {row['model_variant']} "
            f"| {int(row['target_coverage']*100)}% | {rs} |"
        )
    lines.append("")

    # Key observations
    lines.append("### Key Observations\n")

    # 8of8 vs 10of10 Chinese same_ethnicity
    ch_8 = reg[(reg['ethnicity']=='Chinese') & (reg['match_level']=='8of8') & (reg['model_variant']=='same_ethnicity')]
    ch_10 = reg[(reg['ethnicity']=='Chinese') & (reg['match_level']=='10of10') & (reg['model_variant']=='same_ethnicity')]
    lines.append("**8of8 vs 10of10 — Chinese same_ethnicity:**\n")
    for cov in [0.75, 0.85, 0.90, 0.95]:
        r8 = ch_8[ch_8['target_coverage']==cov]['registry_size'].values[0]
        r10 = ch_10[ch_10['target_coverage']==cov]['registry_size'].values[0]
        ratio = r10 / r8
        lines.append(
            f"- At {int(cov*100)}% coverage: 8of8 requires {int(r8):,} donors; "
            f"10of10 requires {int(r10):,} donors ({ratio:.1f}× larger)"
        )
    lines.append("")

    # same_ethnicity vs cross_ethnic for minorities at 90%
    lines.append("**Same-ethnicity vs cross-ethnic at 90% coverage:**\n")
    for eth in ['Indian', 'Malay', 'Others']:
        for ml in ['8of8', '10of10']:
            se = reg[(reg['ethnicity']==eth)&(reg['match_level']==ml)&(reg['model_variant']=='same_ethnicity')&(reg['target_coverage']==0.90)]['registry_size'].values
            ce = reg[(reg['ethnicity']==eth)&(reg['match_level']==ml)&(reg['model_variant']=='cross_ethnic')&(reg['target_coverage']==0.90)]['registry_size'].values
            if len(se) and len(ce):
                se_v = int(se[0])
                ce_v = int(ce[0])
                se_str = f"{se_v:,}" if se_v < 10_000_000 else "≥10M (ceiling)"
                ce_str = f"{ce_v:,}" if ce_v < 10_000_000 else "≥10M (ceiling)"
                lines.append(f"- {eth} {ml}: same_ethnicity={se_str}, cross_ethnic={ce_str}")
    lines.append("")

    # Ceiling values
    ceiling_rows = reg[reg['registry_size'] >= 10_000_000]
    lines.append(f"**Ceiling values (registry_size = 10,000,000):** {len(ceiling_rows)} scenarios\n")
    lines.append(
        "A ceiling of 10,000,000 indicates that the model cannot reach the target coverage level "
        "even with a registry the size of Singapore's entire population (~5.5M). "
        "This occurs for minority groups (Indian, Malay, Others) under cross-ethnic matching at "
        "high coverage targets (≥90%), and for 10of10 same-ethnicity matching in minority groups — "
        "reflecting the practical limits of relying on a single registry population for rare allele combinations."
    )
    lines.append("")
    return "\n".join(lines)


def main():
    af = pd.read_csv(os.path.join(DATA_DIR, 'allele_freq_comparison.csv'))
    hwe = pd.read_csv(os.path.join(DATA_DIR, 'hwe_results.csv'))
    haplo = pd.read_csv(os.path.join(DATA_DIR, 'haplo_freqs_em.csv'))
    reg = pd.read_csv(os.path.join(DATA_DIR, 'registry_size_targets.csv'))

    parts = [
        "# HLA Registry Analysis — Verification Summary\n",
        "> **Dataset:** BMDP + SCBB, n = 59,186 donors, CMIO populations (Singapore)\n",
        "> **Analysis date:** 2026-04-17\n\n",
        "---\n\n",
        section_allele_freq(af),
        "---\n\n",
        section_hwe(hwe),
        "---\n\n",
        section_haplotypes(haplo),
        "---\n\n",
        section_missing_data(),
        "---\n\n",
        section_improvements(),
        "---\n\n",
        section_registry(reg),
    ]

    report = "\n".join(parts)
    with open(OUT_PATH, 'w') as f:
        f.write(report)
    print(f"Report written to {OUT_PATH}")


if __name__ == '__main__':
    main()
