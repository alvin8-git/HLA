#!/usr/bin/env python3
"""Fail if any derived artifact is older than the haplotype table it derives from.

registry_size_ci.csv / registry_ci_plot.png are intentionally absent: the
Dirichlet bootstrap that produced them is biased downward for N* on a
long-tailed distribution, so CIs were withdrawn (see build_report.py 2.4).

Three separate incidents in the 2026-08-20 re-run were the same shape: an
artifact nobody scheduled for regeneration, silently carrying values from an
earlier pipeline configuration into the manuscript (10_ld_report's D' values,
11/15's private FREQ_THRESHOLD constants, and the 14/15 figures omitted from the
overnight chain). Reasoning about which outputs are current does not work; this
checks.

Root of the derivation chain is analysis/data/haplo_freqs_em.csv. Every artifact
below is produced by a script that reads it, or reads hla_clean.csv under
parameters that change with it. All must be at least as new.

Usage:  python check_freshness.py            # exit 1 if anything is stale
        python check_freshness.py --quiet    # only print failures
"""
import os
import sys
from datetime import datetime

ROOT = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(ROOT, 'analysis', 'data')
FIGS = os.path.join(ROOT, 'analysis', 'figures')
REFERENCE = os.path.join(DATA, 'haplo_freqs_em.csv')

# artifact -> producing script (for the failure message)
DERIVED = {
    # --- data, produced by scripts that read the haplotype table -------------
    'data/hwe_results.csv':             '03_hwe_test.py',
    'data/allele_freqs_per_locus.csv':  '03_hwe_test.py',
    'data/coverage_curves.csv':         '04_registry_model.py',
    'data/registry_size_targets.csv':   '04_registry_model.py',
    'data/em_validation.csv':           '07_validate_em.py',
    'data/em_validation_summary.csv':   '07_validate_em.py',
    'data/ld_report.csv':               '10_ld_report.py',
    'data/cross_ethnic_sensitivity.csv': '13_cross_ethnic_sensitivity.py',
    'data/smoothing_sensitivity.csv':   '16_smoothing_sensitivity.py',
    'data/match_validation.csv':        '12_match_validation.py',
    'data/match_rate_comparison.csv':   '12_match_validation.py',
    # own EM, but its FREQ_THRESHOLD tracks the pipeline floor
    'data/em_convergence.csv':          '15_em_convergence.py',
    'data/others_cluster_registry.csv': '11_others_stratification.py',
    'data/others_cluster_assignments.csv': '11_others_stratification.py',
    'data/others_cluster_haplotypes.csv':  '11_others_stratification.py',
    'data/others_cluster_silhouette.csv':  '11_others_stratification.py',
    # --- figures embedded in the manuscript ---------------------------------
    'figures/pipeline_flowchart.png':      '14_pipeline_flowchart.py',
    'figures/partial_match_10locus.png':   'paper_BMT_workdir/partial_match_figs.py',
    'figures/partial_match_8locus.png':    'paper_BMT_workdir/partial_match_figs.py',
    'figures/cross_ethnic_sensitivity.png': '13_cross_ethnic_sensitivity.py',
    'figures/match_validation_scatter.png': '12_match_validation.py',
    'figures/others_pca_scatter.png':       '11_others_stratification.py',
    'figures/em_convergence.png':           '15_em_convergence.py',
    # --- other figures ------------------------------------------------------
    'figures/coverage_curves_10of10.png':   '04_registry_model.py',
    'figures/coverage_curves_8of8.png':     '04_registry_model.py',
    'figures/diplotype_longtail.png':       '08_section_plots.py',
    'figures/registry_targets_bar.png':     '08_section_plots.py',
    'figures/ld_heatmap_dprime.png':        '10_ld_report.py',
    'figures/ld_heatmap_r2.png':            '10_ld_report.py',
    'figures/others_registry_by_cluster.png': '11_others_stratification.py',
}

# Figures embedded in the manuscript — a stale one here is a reader-visible error.
EMBEDDED = {k for k in DERIVED if k.startswith('figures/') and any(
    t in k for t in ('pipeline_flowchart', 'coverage_curves_10of10', 'partial_match',
                     'cross_ethnic_sensitivity', 'match_validation_scatter',
                     'others_pca_scatter', 'em_convergence'))}


def stamp(t):
    return datetime.fromtimestamp(t).strftime('%m-%d %H:%M')


def main():
    quiet = '--quiet' in sys.argv
    if not os.path.exists(REFERENCE):
        print(f"FAIL: reference {REFERENCE} missing"); return 1
    ref = os.path.getmtime(REFERENCE)
    print(f"reference: analysis/data/haplo_freqs_em.csv  ({stamp(ref)})\n")

    # Siblings written by the same script as the reference are emitted BEFORE it
    # within one invocation (03 writes HWE/allele tables, then spends ~1 min on
    # the EM), so they are legitimately a little older. Allow a window for those
    # only; everything else must be strictly newer.
    SIBLING_GRACE = 900  # seconds

    stale, missing, ok = [], [], []
    for rel, script in sorted(DERIVED.items()):
        path = os.path.join(ROOT, 'analysis', rel)
        if not os.path.exists(path):
            missing.append((rel, script)); continue
        mt = os.path.getmtime(path)
        threshold = ref - SIBLING_GRACE if script == '03_hwe_test.py' else ref
        (ok if mt >= threshold else stale).append((rel, script, mt))

    if not quiet:
        for rel, script, t in ok:
            print(f"  OK    {rel:<42} {stamp(t)}")
    for rel, script, t in stale:
        tag = 'STALE*' if rel in EMBEDDED else 'STALE '
        print(f"  {tag}{rel:<42} {stamp(t)}  <- rerun {script}")
    for rel, script in missing:
        print(f"  MISSING {rel:<40} <- run {script}")

    print(f"\n{len(ok)} current, {len(stale)} stale, {len(missing)} missing"
          f"   (* = embedded in the manuscript)")
    if stale or missing:
        print("RESULT: FAIL — do not build or commit until these are regenerated")
        return 1
    print("RESULT: PASS — every derived artifact is at least as new as the haplotype table")
    return 0


if __name__ == '__main__':
    sys.exit(main())
