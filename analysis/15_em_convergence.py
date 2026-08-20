"""
15_em_convergence.py
Tests whether N* stabilises at the 5,000-sample EM cap.
Subsamples Chinese donors at increasing sizes, runs EM, computes N* at 95%.
Outputs:
  data/em_convergence.csv
  figures/em_convergence.png
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from hwe_test import _pivot_to_individuals, _em_full_phase
from registry_model import get_diplotype_frequencies, find_registry_size

DATA_DIR = os.path.join(HERE, 'data')
FIG_DIR  = os.path.join(HERE, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

FREQ_THRESHOLD = 0.000001
SEED = 42
SAMPLE_SIZES = [500, 1000, 2000, 3000, 5000, 7500, 10000, 15000, 20000, 30000, None]
COVERAGE = 0.95


def run_for_size(pivoted_full: pd.DataFrame, n: int | None, seed: int) -> dict:
    if n is not None and len(pivoted_full) > n:
        pivoted = pivoted_full.sample(n=n, random_state=seed)
    else:
        pivoted = pivoted_full
    actual_n = len(pivoted)

    haplo_freqs = _em_full_phase(pivoted)
    rows = [{'haplotype': '|'.join(h), 'frequency': f}
            for h, f in haplo_freqs.items() if f >= FREQ_THRESHOLD]
    if not rows:
        return {'sample_size': actual_n, 'n_haplotypes': 0, 'registry_size': np.nan}

    haplo_df = pd.DataFrame(rows)
    haplo_df['frequency'] /= haplo_df['frequency'].sum()
    diplo = get_diplotype_frequencies(haplo_df)
    n_star = find_registry_size(diplo, COVERAGE)
    return {
        'sample_size':   actual_n,
        'n_haplotypes':  len(haplo_df),
        'registry_size': n_star,
    }


def main():
    print("Loading Chinese data...")
    hla_df = pd.read_csv(os.path.join(DATA_DIR, 'hla_clean.csv'))
    chn_df = hla_df[hla_df['ethnicity'] == 'Chinese'].copy()
    pivoted_full = _pivot_to_individuals(chn_df)
    print(f"  {len(pivoted_full):,} Chinese individuals with complete 5-locus typing")

    rows = []
    for n in SAMPLE_SIZES:
        label = f'{n:,}' if n else f'{len(pivoted_full):,} (full)'
        print(f"  Sample size {label}...", end=' ', flush=True)
        result = run_for_size(pivoted_full, n, SEED)
        print(f"N*={result['registry_size']:,}  ({result['n_haplotypes']} haplotypes)")
        rows.append(result)

    df = pd.DataFrame(rows)
    csv_path = os.path.join(DATA_DIR, 'em_convergence.csv')
    df.to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")

    # ── Figure ────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 4), facecolor='white')
    ax.set_facecolor('white')

    valid = df.dropna(subset=['registry_size'])
    ax.plot(valid['sample_size'], valid['registry_size'],
            'o-', color='#2c7bb6', linewidth=2, markersize=6, zorder=3)

    # Mark the 5,000 cap used in the main analysis
    ax.axvline(5000, color='#d7191c', linestyle='--', linewidth=1.5,
               label='Current EM cap (5,000)')
    n_star_5k = df[df['sample_size'] == 5000]['registry_size'].values
    if len(n_star_5k):
        ax.axhline(n_star_5k[0], color='#d7191c', linestyle=':', linewidth=1.0,
                   alpha=0.6)

    ax.set_xscale('log')
    ax.set_xlabel('EM input sample size (Chinese donors)', fontsize=11)
    ax.set_ylabel('N* at 95% coverage (10/10)', fontsize=11)
    ax.set_title('EM convergence: N* vs sample size (Chinese)',
                 fontsize=11, fontweight='bold')
    ax.legend(fontsize=9)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f'{int(v):,}'))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    fig_path = os.path.join(FIG_DIR, 'em_convergence.png')
    fig.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {fig_path}")

    # Report convergence
    n_star_full = df[df['sample_size'] == df['sample_size'].max()]['registry_size'].values[0]
    pct_diff = abs(n_star_5k[0] - n_star_full) / n_star_full * 100 if len(n_star_5k) else np.nan
    print(f"\nN* at cap (5k): {n_star_5k[0]:,}")
    print(f"N* at full ({len(pivoted_full):,}): {n_star_full:,}")
    print(f"Difference: {pct_diff:.2f}%")


if __name__ == '__main__':
    main()
