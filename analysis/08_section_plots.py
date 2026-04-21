"""
08_section_plots.py
Supplementary figures for Documentation.md sections 6.3.2 and 6.3.3.

Figures produced:
  figures/diplotype_longtail.png   — cumulative diplotype frequency coverage
                                     vs number of diplotypes (all ethnicities)
  figures/registry_targets_bar.png — grouped bar chart of minimum registry N
                                     at 75/85/90/95% coverage (10/10, same-eth)
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, 'data')
FIG_DIR  = os.path.join(HERE, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

ETHNICITIES = ['Chinese', 'Malay', 'Indian', 'Others']
COLORS = {'Chinese': '#d62728', 'Malay': '#2ca02c',
          'Indian': '#9467bd', 'Others': '#DAA520'}


# ---------------------------------------------------------------------------
# Figure A — Diplotype long-tail: cumulative frequency vs rank (6.3.3)
# ---------------------------------------------------------------------------

def make_longtail_figure(haplo_df: pd.DataFrame, out_path: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 5), facecolor='white')
    ax.set_facecolor('white')

    for eth in ETHNICITIES:
        h = haplo_df[haplo_df['ethnicity'] == eth].sort_values(
            'frequency', ascending=False).reset_index(drop=True)
        freqs = h['frequency'].values
        n_hap = len(freqs)

        # Build sorted diplotype frequencies
        dip_freqs = []
        for i in range(n_hap):
            for j in range(i, n_hap):
                f = freqs[i]**2 if i == j else 2 * freqs[i] * freqs[j]
                dip_freqs.append(f)
        dip_freqs = np.array(sorted(dip_freqs, reverse=True))
        total = dip_freqs.sum()

        cumfreq = np.cumsum(dip_freqs) / total * 100
        ranks = np.arange(1, len(dip_freqs) + 1)

        ax.plot(ranks, cumfreq, color=COLORS[eth], linewidth=2, label=eth)

    # Reference lines
    for y in [50, 75, 90, 95]:
        ax.axhline(y, color='lightgrey', linestyle='--', linewidth=0.8, zorder=0)
        ax.text(10, y + 0.8, f'{y}%', fontsize=8, color='grey')

    ax.set_xscale('log')
    ax.set_xlim(1, None)
    ax.set_ylim(0, 102)
    ax.set_xlabel('Number of diplotypes (ranked by frequency)', fontsize=11)
    ax.set_ylabel('Cumulative share of haplotype pool (%)', fontsize=11)
    ax.set_title('Diplotype long-tail: cumulative frequency coverage\nby top-K diplotypes (same-ethnicity EM pool)',
                 fontsize=11, fontweight='bold')
    ax.legend(fontsize=10, loc='lower right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {out_path}')


# ---------------------------------------------------------------------------
# Figure B — Registry size bar chart (6.3.2)
# ---------------------------------------------------------------------------

def make_registry_bar_figure(targets_df: pd.DataFrame, out_path: str) -> None:
    t = targets_df[
        (targets_df['match_level'] == '10of10') &
        (targets_df['model_variant'] == 'same_ethnicity') &
        (targets_df['ethnicity'].isin(ETHNICITIES))
    ].copy()

    thresholds = [0.75, 0.85, 0.90, 0.95]
    threshold_labels = ['75%', '85%', '90%', '95%']
    n_eth = len(ETHNICITIES)
    n_thr = len(thresholds)
    x = np.arange(n_thr)
    bar_width = 0.18

    fig, ax = plt.subplots(figsize=(9, 5), facecolor='white')
    ax.set_facecolor('white')

    for k, eth in enumerate(ETHNICITIES):
        eth_data = t[t['ethnicity'] == eth].set_index('target_coverage')['registry_size']
        vals = [eth_data.get(thr, np.nan) for thr in thresholds]
        offset = (k - (n_eth - 1) / 2) * bar_width
        bars = ax.bar(x + offset, vals, bar_width, label=eth,
                      color=COLORS[eth], alpha=0.85, edgecolor='white')
        # Value labels on top
        for bar, v in zip(bars, vals):
            if not np.isnan(v):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 300,
                        f'{int(v):,}', ha='center', va='bottom',
                        fontsize=7, rotation=90)

    ax.set_xticks(x)
    ax.set_xticklabels(threshold_labels, fontsize=11)
    ax.set_xlabel('Coverage target', fontsize=11)
    ax.set_ylabel('Minimum registry size (donors)', fontsize=11)
    ax.set_title('Minimum registry size for same-ethnicity 10/10 matching\nby CMIO group and coverage target (full multi-locus EM)',
                 fontsize=11, fontweight='bold')
    ax.legend(fontsize=10)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.15)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x):,}'))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for y in [10000, 20000, 30000, 40000]:
        ax.axhline(y, color='lightgrey', linestyle='--', linewidth=0.7, zorder=0)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {out_path}')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    haplo_df   = pd.read_csv(os.path.join(DATA_DIR, 'haplo_freqs_em.csv'))
    targets_df = pd.read_csv(os.path.join(DATA_DIR, 'registry_size_targets.csv'))

    print('Generating long-tail figure...')
    make_longtail_figure(haplo_df, os.path.join(FIG_DIR, 'diplotype_longtail.png'))

    print('Generating registry size bar chart...')
    make_registry_bar_figure(targets_df, os.path.join(FIG_DIR, 'registry_targets_bar.png'))

    print('Done.')


if __name__ == '__main__':
    main()
