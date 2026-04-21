"""
12_match_validation.py
Donor-patient match rate validation using DonorPatient.txt and Patient.txt.

Approach:
  1. Parse Patient.txt (564 haplotype rows, cols: ethnicity A B C DRB1 DQB1).
  2. Compute observed haplotype frequencies per ethnicity from Patient.txt.
  3. Compare to EM-estimated frequencies: Spearman r and RMSE.
  4. For each patient haplotype, compute model-predicted match probability
     (P(donor carries same haplotype) from EM diplotype distribution).
  5. Compare predicted vs observed match rate (fraction of DonorPatient rows
     that share a haplotype with a Patient row of the same ethnicity).

Outputs:
  data/match_validation.csv         — per-ethnicity correlation stats
  data/match_rate_comparison.csv    — predicted vs observed match probability
  figures/match_validation_scatter.png — observed vs predicted freq scatter
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, 'data')
FIG_DIR  = os.path.join(HERE, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

ETH_MAP  = {'C': 'Chinese', 'M': 'Malay', 'I': 'Indian', 'O': 'Others'}
LOCI     = ['HLA-A', 'HLA-B', 'HLA-C', 'DRB1', 'DQB1']
COLORS   = {'Chinese': '#d62728', 'Malay': '#2ca02c',
            'Indian': '#9467bd', 'Others': '#DAA520'}


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_haplotype_file(path: str) -> pd.DataFrame:
    """Parse tab-separated haplotype file (no header, cols: eth A B C DRB1 DQB1)."""
    df = pd.read_csv(path, sep='\t', header=None,
                     names=['eth_code', 'A', 'B', 'C', 'DRB1', 'DQB1'])
    df['ethnicity'] = df['eth_code'].map(ETH_MAP)
    # Replace missing marker
    for col in ['A', 'B', 'C', 'DRB1', 'DQB1']:
        df[col] = df[col].replace('-', np.nan)
    return df


def make_haplotype_str(row, use_dqb1: bool = True) -> str | None:
    """Build pipe-separated haplotype string. Returns None if any allele missing."""
    alleles = [row['A'], row['B'], row['C'], row['DRB1']]
    if use_dqb1:
        alleles.append(row['DQB1'])
    if any(pd.isna(a) for a in alleles):
        return None
    return '|'.join(alleles)


# ---------------------------------------------------------------------------
# Frequency comparison
# ---------------------------------------------------------------------------

def observed_haplotype_freqs(df: pd.DataFrame, eth: str,
                              use_dqb1: bool = True) -> pd.Series:
    """Count haplotype occurrences in df for one ethnicity, normalise to freq."""
    sub = df[df['ethnicity'] == eth].copy()
    sub['haplotype'] = sub.apply(lambda r: make_haplotype_str(r, use_dqb1), axis=1)
    sub = sub.dropna(subset=['haplotype'])
    counts = sub['haplotype'].value_counts()
    return counts / counts.sum()


def compare_frequencies(obs: pd.Series, em: pd.Series,
                         eth: str) -> dict:
    """Spearman r and RMSE between observed and EM frequencies on shared haplotypes."""
    shared = obs.index.intersection(em.index)
    if len(shared) < 3:
        return {'ethnicity': eth, 'n_shared': len(shared),
                'spearman_r': np.nan, 'p_value': np.nan, 'rmse': np.nan}

    o = obs[shared].values
    e = em[shared].values
    r, p = spearmanr(o, e)
    rmse = np.sqrt(np.mean((o - e) ** 2))
    return {
        'ethnicity':  eth,
        'n_patient':  len(obs),
        'n_em':       len(em),
        'n_shared':   len(shared),
        'pct_covered': round(obs[shared].sum() * 100, 1),
        'spearman_r': round(r, 4),
        'p_value':    round(p, 6),
        'rmse':       round(rmse, 6),
    }


# ---------------------------------------------------------------------------
# Match rate validation
# ---------------------------------------------------------------------------

def compute_predicted_match_prob(em_freqs: pd.Series, patient_haps: list) -> float:
    """
    For each patient haplotype h, P(donor carries h) ≈ 2f(h)(1-f(h)) + f(h)^2 = 2f(h) - f(h)^2.
    Average across patient haplotypes.
    """
    probs = []
    for h in patient_haps:
        f = em_freqs.get(h, 0.0)
        probs.append(2 * f - f ** 2)
    return float(np.mean(probs)) if probs else 0.0


def compute_observed_match_rate(patient_df: pd.DataFrame,
                                 donor_patient_df: pd.DataFrame,
                                 eth: str) -> float:
    """
    Fraction of patient haplotypes (in Patient.txt) that also appear in
    DonorPatient.txt for the same ethnicity (proxy for donor availability).
    """
    pat = patient_df[patient_df['ethnicity'] == eth].copy()
    dp  = donor_patient_df[donor_patient_df['ethnicity'] == eth].copy()

    pat['haplotype'] = pat.apply(lambda r: make_haplotype_str(r, use_dqb1=False), axis=1)
    dp['haplotype']  = dp.apply(lambda r: make_haplotype_str(r, use_dqb1=False), axis=1)

    pat = pat.dropna(subset=['haplotype'])
    dp  = dp.dropna(subset=['haplotype'])

    if len(pat) == 0:
        return np.nan

    dp_haps = set(dp['haplotype'])
    matched = pat['haplotype'].isin(dp_haps).sum()
    return round(matched / len(pat), 4)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def make_scatter_plot(patient_df: pd.DataFrame, em_df: pd.DataFrame,
                      out_path: str) -> None:
    """Scatter: observed (patient) vs EM-predicted haplotype frequency per ethnicity."""
    ethnicities = [e for e in ['Chinese', 'Malay', 'Indian', 'Others']
                   if patient_df['ethnicity'].eq(e).any()]
    n = len(ethnicities)
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 4), facecolor='white')
    if n == 1:
        axes = [axes]

    for ax, eth in zip(axes, ethnicities):
        ax.set_facecolor('white')
        obs = observed_haplotype_freqs(patient_df, eth, use_dqb1=True)
        em_eth = em_df[em_df['ethnicity'] == eth].set_index('haplotype')['frequency']
        shared = obs.index.intersection(em_eth.index)

        if len(shared) > 0:
            x = em_eth[shared].values
            y = obs[shared].values
            ax.scatter(x, y, c=COLORS.get(eth, 'grey'), s=20, alpha=0.6)
            lim = max(x.max(), y.max()) * 1.1
            ax.plot([0, lim], [0, lim], 'k--', linewidth=0.8, alpha=0.5)
            r, _ = spearmanr(x, y)
            ax.text(0.05, 0.92, f'r={r:.3f}\nn={len(shared)}',
                    transform=ax.transAxes, fontsize=9)

        ax.set_xlabel('EM-estimated frequency', fontsize=9)
        ax.set_ylabel('Observed frequency (Patient.txt)', fontsize=9)
        ax.set_title(eth, fontsize=10, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.suptitle('Observed vs EM-predicted haplotype frequencies\n(Patient.txt validation)',
                 fontsize=11, fontweight='bold')
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {out_path}')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    root = os.path.dirname(HERE)

    print("Parsing haplotype files...")
    patient_df      = parse_haplotype_file(os.path.join(root, 'Patient.txt'))
    donor_pat_df    = parse_haplotype_file(os.path.join(root, 'DonorPatient.txt'))
    em_df           = pd.read_csv(os.path.join(DATA_DIR, 'haplo_freqs_em.csv'))

    print(f"  Patient.txt:      {len(patient_df)} haplotype rows")
    print(f"  DonorPatient.txt: {len(donor_pat_df)} haplotype rows")
    print(f"  EM haplotypes:    {len(em_df)}")

    # Ethnicity breakdown
    for eth_code, eth in ETH_MAP.items():
        n_pat = patient_df['ethnicity'].eq(eth).sum()
        n_dp  = donor_pat_df['ethnicity'].eq(eth).sum()
        print(f"  {eth}: Patient={n_pat}  DonorPatient={n_dp}")

    print("\nComparing haplotype frequencies (5-locus)...")
    val_rows = []
    for eth in ETH_MAP.values():
        obs_5l = observed_haplotype_freqs(patient_df, eth, use_dqb1=True)
        em_eth = em_df[em_df['ethnicity'] == eth].set_index('haplotype')['frequency']
        stats  = compare_frequencies(obs_5l, em_eth, eth)
        val_rows.append(stats)
        print(f"  {eth}: shared={stats['n_shared']}  r={stats['spearman_r']}  "
              f"RMSE={stats['rmse']}  coverage={stats.get('pct_covered','?')}%")

    val_df = pd.DataFrame(val_rows)
    val_df.to_csv(os.path.join(DATA_DIR, 'match_validation.csv'), index=False)

    print("\nComputing predicted vs observed match rates (4-locus, no DQB1)...")
    rate_rows = []
    for eth in ETH_MAP.values():
        # 4-locus EM frequencies (collapse DQB1)
        em_eth_5l = em_df[em_df['ethnicity'] == eth].copy()
        em_eth_5l['haplotype_4l'] = em_eth_5l['haplotype'].apply(
            lambda h: '|'.join(h.split('|')[:4]))
        em_4l = em_eth_5l.groupby('haplotype_4l')['frequency'].sum()
        em_4l = em_4l / em_4l.sum()

        pat_sub = patient_df[patient_df['ethnicity'] == eth]
        pat_haps_4l = [make_haplotype_str(r, use_dqb1=False)
                       for _, r in pat_sub.iterrows()]
        pat_haps_4l = [h for h in pat_haps_4l if h is not None]

        pred_match = compute_predicted_match_prob(em_4l, pat_haps_4l)
        obs_match  = compute_observed_match_rate(patient_df, donor_pat_df, eth)

        print(f"  {eth}: predicted={pred_match:.4f}  observed={obs_match:.4f}  "
              f"n_patient_haps={len(pat_haps_4l)}")
        rate_rows.append({
            'ethnicity':            eth,
            'n_patient_haplotypes': len(pat_haps_4l),
            'predicted_match_prob': round(pred_match, 4),
            'observed_match_rate':  obs_match,
        })

    rate_df = pd.DataFrame(rate_rows)
    rate_df.to_csv(os.path.join(DATA_DIR, 'match_rate_comparison.csv'), index=False)

    print("\nGenerating scatter plot...")
    make_scatter_plot(patient_df, em_df,
                      os.path.join(FIG_DIR, 'match_validation_scatter.png'))

    print('\n=== Summary ===')
    print(val_df[['ethnicity', 'n_shared', 'spearman_r', 'rmse']].to_string(index=False))
    print()
    print(rate_df.to_string(index=False))
    print('\nDone.')


if __name__ == '__main__':
    main()
