"""
plot_coverage.py — Standalone script to generate coverage curve figures
from analysis/data/coverage_curves.csv -> analysis/figures/
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, 'data')
FIG_DIR = os.path.join(HERE, 'figures')

ETHNICITIES = ['Chinese', 'Malay', 'Indian', 'Others']
THRESHOLDS = [0.75, 0.85, 0.90, 0.95]

colors = {'same_ethnicity': 'steelblue', 'cross_ethnic': 'darkorange'}
styles = {'same_ethnicity': '-', 'cross_ethnic': '--'}
labels = {'same_ethnicity': 'Same ethnicity', 'cross_ethnic': 'Cross-ethnic (SG pool)'}


def make_figure(df_ml, match_level, fig_path):
    all_ethnicities = ETHNICITIES + ['Combined']
    fig, axes = plt.subplots(1, 5, figsize=(20, 4), sharey=True)

    for ax, eth in zip(axes, all_ethnicities):
        eth_df = df_ml[df_ml['ethnicity'] == eth]
        for mv, grp in eth_df.groupby('model_variant'):
            grp_sorted = grp.sort_values('N')
            ax.plot(grp_sorted['N'], grp_sorted['coverage'],
                    color=colors[mv], linestyle=styles[mv],
                    label=labels[mv], linewidth=1.5)

        for threshold in THRESHOLDS:
            ax.axhline(threshold, color='grey', linestyle='--', linewidth=0.8, alpha=0.6)

        ax.set_xscale('log')
        ax.set_ylim(0, 1)
        ax.set_title(eth)
        ax.set_xlabel('Registry size (N)')

    axes[0].set_ylabel('Coverage')

    # Legend in first panel only, deduplicated
    handles, lbls = axes[0].get_legend_handles_labels()
    seen = {}
    for h, l in zip(handles, lbls):
        if l not in seen:
            seen[l] = h
    axes[0].legend(seen.values(), seen.keys(), fontsize=8, loc='lower right')

    fig.suptitle(f'Registry Coverage — {match_level} match', fontsize=13, y=1.01)
    plt.tight_layout()
    os.makedirs(os.path.dirname(fig_path), exist_ok=True)
    fig.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {fig_path}")


def main():
    csv_path = os.path.join(DATA_DIR, 'coverage_curves.csv')
    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} rows from {csv_path}")

    os.makedirs(FIG_DIR, exist_ok=True)

    for match_level in ['8of8', '10of10']:
        df_ml = df[df['match_level'] == match_level]
        fig_path = os.path.join(FIG_DIR, f'coverage_curves_{match_level}.png')
        make_figure(df_ml, match_level, fig_path)

    print("Done.")


if __name__ == '__main__':
    main()
