"""
09_bootstrap_ci.py
Parametric bootstrap confidence intervals for registry size targets.

Method: Dirichlet resampling of haplotype frequencies (n_eff × f_hat as
concentration parameters), then compute N* at each coverage threshold from
the resampled diplotype distribution. B=1,000 resamples per ethnicity.

Point estimate: bootstrap MEDIAN (bias-corrected for the left-skew inherent in
N*(f) near saturation, which arises from Jensen's inequality on the concave
coverage function). The EM maximum-likelihood estimate is also stored for
reference as 'em_registry_size'.

CI: standard 2.5th–97.5th percentile interval. The bootstrap median is always
inside this interval by construction, satisfying the convention that point
estimates are reported within their confidence intervals.

n_eff: actual 5-locus donor counts per ethnicity (from hla_clean.csv), used as
Dirichlet concentration parameters to reflect true estimation precision.

Outputs:
  data/registry_size_ci.csv       — median point estimate + 95% CI per scenario
  figures/registry_ci_plot.png    — forest-style CI plot
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from registry_model import (get_diplotype_frequencies, find_registry_size,
                            diplotype_freq_vector, find_registry_size_vec)
from multiprocessing import Pool, cpu_count
from functools import partial

DATA_DIR = os.path.join(HERE, 'data')
FIG_DIR  = os.path.join(HERE, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

ETHNICITIES  = ['Chinese', 'Malay', 'Indian', 'Others']
THRESHOLDS   = [0.75, 0.85, 0.90, 0.95]
MATCH_LEVELS = ['10of10', '8of8']
B            = int(os.environ.get('HLA_BOOTSTRAP_B', 1000))
SEED         = 42

# Actual 5-locus donor counts from hla_clean.csv (used as Dirichlet n_eff).
# Using true sample sizes (not capped) to reflect genuine estimation precision
# and reduce the Jensen downward-bias in N* at high coverage thresholds.
# Corrected 2026-08-20: the previous values (45754/5868/5586/3941) summed to
# 61,149, exceeding the 59,186 stated study total. These are the counts
# reproducible from hla_clean.csv for individuals typed at all five loci.
N_EFF = {'Chinese': 44400, 'Malay': 5578, 'Indian': 5490, 'Others': 3767}

COLORS = {'Chinese': '#d62728', 'Malay': '#2ca02c',
          'Indian': '#9467bd',  'Others': '#DAA520'}


def _boot_one(sampled, thresholds):
    """Evaluate one bootstrap replicate. Module-level so it is picklable."""
    vec = diplotype_freq_vector(sampled)
    return {t: find_registry_size_vec(vec, t) for t in thresholds}


def bootstrap_registry_ci(freqs: np.ndarray, n_eff: int,
                           haplotypes: list, thresholds: list,
                           match_level: str, B: int, rng) -> dict:
    """
    Dirichlet parametric bootstrap for registry size CIs.

    freqs      : 1-D array of haplotype frequencies (sum ≤ 1; residual pooled)
    n_eff      : effective sample size (actual 5-locus donor count)
    haplotypes : list of haplotype strings (same order as freqs)
    thresholds : list of coverage targets
    match_level: '10of10' or '8of8'
    B          : number of bootstrap replicates
    rng        : np.random.Generator

    Returns dict: {threshold: {'point': int, 'em_point': int, 'lo': int,
                                'hi': int, 'pct_below': float}}
    where 'point' is the bootstrap median (bias-corrected point estimate)
    and 'em_point' is the EM maximum-likelihood estimate.
    """
    # EM maximum-likelihood point estimate
    haplo_df = _make_haplo_df('__eth__', haplotypes, freqs, match_level)
    diplo_df = get_diplotype_frequencies(haplo_df)
    em_point = {t: find_registry_size(diplo_df, t) for t in thresholds}

    # Normalise frequencies; Dirichlet concentration α_k = n_eff × f_k
    freqs_norm = freqs / freqs.sum()
    alpha = np.maximum(freqs_norm * n_eff, 0.1)

    # Draw every replicate in the parent so the rng sequence — and therefore
    # the result — is identical to the original serial implementation, then
    # evaluate them across cores. At a 1e-4 haplotype floor each replicate
    # expands to millions of diplotypes, so serial evaluation would take days.
    samples = [rng.dirichlet(alpha) for _ in range(B)]

    n_proc = max(1, min(cpu_count() - 2, 32))
    with Pool(n_proc) as pool:
        per_replicate = pool.map(
            partial(_boot_one, thresholds=thresholds), samples, chunksize=4)

    boot_results = {t: [r[t] for r in per_replicate] for t in thresholds}

    result = {}
    for t in thresholds:
        arr = np.array(boot_results[t])
        median = int(np.median(arr))
        lo = int(np.percentile(arr, 2.5))
        hi = int(np.percentile(arr, 97.5))
        pct_below = float(np.mean(arr < em_point[t]))
        result[t] = {
            'point':     median,       # bootstrap median (bias-corrected)
            'em_point':  em_point[t],  # EM MLE (for reference)
            'lo':        lo,
            'hi':        hi,
            'pct_below': round(pct_below, 4),
        }
    return result


def _make_haplo_df(eth, haplotypes, freqs, match_level):
    """Build a per-ethnicity haplotype DataFrame for get_diplotype_frequencies.
    Frequencies are normalised to sum to 1 so the diplotype distribution is proper."""
    total = freqs.sum()
    freqs_norm = freqs / total if total > 0 else freqs

    if match_level == '8of8':
        trimmed = ['|'.join(h.split('|')[:4]) for h in haplotypes]
        df = pd.DataFrame({'ethnicity': eth, 'haplotype': trimmed, 'frequency': freqs_norm})
        df = df.groupby(['ethnicity', 'haplotype'], as_index=False)['frequency'].sum()
        df['frequency'] = df['frequency'] / df['frequency'].sum()
    else:
        df = pd.DataFrame({'ethnicity': eth, 'haplotype': haplotypes, 'frequency': freqs_norm})
    return df


def make_ci_plot(ci_df: pd.DataFrame, out_path: str) -> None:
    """Forest-style CI plot: one panel per match level, one row per ethnicity × threshold."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), facecolor='white')
    threshold_labels = {0.75: '75%', 0.85: '85%', 0.90: '90%', 0.95: '95%'}
    y_gap = 1.2

    for ax, ml in zip(axes, MATCH_LEVELS):
        ax.set_facecolor('white')
        sub = ci_df[ci_df['match_level'] == ml].sort_values(
            ['target_coverage', 'ethnicity'], ascending=[True, True])

        y_pos = 0
        yticks, ylabels = [], []

        for t in sorted(THRESHOLDS):
            rows = sub[sub['target_coverage'] == t]
            for _, row in rows.iterrows():
                eth = row['ethnicity']
                color = COLORS[eth]
                ax.plot(row['registry_size'], y_pos, 'o', color=color, markersize=7, zorder=3)
                ax.hlines(y_pos, row['ci_lo'], row['ci_hi'], color=color,
                          linewidth=3, alpha=0.7, zorder=2)
                yticks.append(y_pos)
                ylabels.append(f"{threshold_labels[t]} – {eth}")
                y_pos += y_gap

            y_pos += y_gap * 0.5

        ax.set_yticks(yticks)
        ax.set_yticklabels(ylabels, fontsize=8)
        ax.set_xlabel('Registry size (donors)', fontsize=10)
        ax.set_title(f'{ml.replace("of", "/")} match\n(95% bootstrap CI, median estimate)',
                     fontsize=11, fontweight='bold')
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x):,}'))
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        ax.grid(axis='x', color='lightgrey', linestyle='--', linewidth=0.7, zorder=0)

    handles = [plt.Line2D([0], [0], color=COLORS[e], marker='o', linewidth=2,
                           markersize=6, label=e) for e in ETHNICITIES]
    fig.legend(handles=handles, loc='lower center', ncol=4, fontsize=9,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle('Registry size targets with 95% bootstrap confidence intervals\n'
                 '(Dirichlet resampling, B=1,000, median estimate, same-ethnicity matching)',
                 fontsize=12, fontweight='bold')
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {out_path}')


def main():
    haplo_df = pd.read_csv(os.path.join(DATA_DIR, 'haplo_freqs_em.csv'))
    rng = np.random.default_rng(SEED)

    rows = []
    for eth in ETHNICITIES:
        n_eff = N_EFF[eth]
        h = haplo_df[haplo_df['ethnicity'] == eth].sort_values(
            'frequency', ascending=False)
        haplotypes = list(h['haplotype'])
        freqs = h['frequency'].values

        for ml in MATCH_LEVELS:
            print(f'  Bootstrap {eth} {ml} (n_eff={n_eff}, K={len(haplotypes)}, B={B})...',
                  flush=True)
            ci = bootstrap_registry_ci(freqs, n_eff, haplotypes,
                                       THRESHOLDS, ml, B, rng)
            for t, vals in ci.items():
                rows.append({
                    'ethnicity':         eth,
                    'match_level':       ml,
                    'target_coverage':   t,
                    'registry_size':     vals['point'],    # bootstrap median
                    'em_registry_size':  vals['em_point'], # EM MLE (reference)
                    'ci_lo':             vals['lo'],
                    'ci_hi':             vals['hi'],
                    'ci_width':          vals['hi'] - vals['lo'],
                    'pct_below':         vals['pct_below'],
                })
            if ml == '10of10':
                print(f'    95% CI: {ci[0.95]["lo"]:,} – {ci[0.95]["hi"]:,}'
                      f'  (median={ci[0.95]["point"]:,}, EM={ci[0.95]["em_point"]:,},'
                      f' pct_below={ci[0.95]["pct_below"]:.3f})')

    ci_df = pd.DataFrame(rows)
    out_csv = os.path.join(DATA_DIR, 'registry_size_ci.csv')
    ci_df.to_csv(out_csv, index=False)
    print(f'\nSaved: {out_csv}')

    print('\n=== 95% CI summary (10/10, same-ethnicity, 95% coverage) ===')
    mask = (ci_df['match_level'] == '10of10') & (ci_df['target_coverage'] == 0.95)
    cols = ['ethnicity', 'registry_size', 'em_registry_size', 'ci_lo', 'ci_hi', 'pct_below']
    print(ci_df[mask][cols].to_string(index=False))

    # Sanity check: point estimate should be within CI
    violations = ci_df[(ci_df['registry_size'] < ci_df['ci_lo']) |
                        (ci_df['registry_size'] > ci_df['ci_hi'])]
    if violations.empty:
        print('\n✓ All point estimates are within their 95% CIs.')
    else:
        print(f'\n⚠ {len(violations)} point estimates outside CI:')
        print(violations[['ethnicity', 'match_level', 'target_coverage',
                           'registry_size', 'ci_lo', 'ci_hi']].to_string(index=False))

    make_ci_plot(ci_df, os.path.join(FIG_DIR, 'registry_ci_plot.png'))
    print('\nDone.')


if __name__ == '__main__':
    main()
